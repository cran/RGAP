# -------------------------------------------------------------------------------------------

#' Estimates the parameters and states of a two-dimensional state-space model by Bayesian
#' methods to obtain the tfp trend.
#'
#' @inheritParams fit.TFPmodel
#' @param FUN A function to be used to compute estimates from the posterior distribution.
#'   Possible options are \code{"mean"} and \code{"median"}. The default is \code{FUN = "mean"}.
#'   Only used if \code{method = "bayesian"}.
#'
#' @importFrom KFAS fitSSM KFS simulateSSM
#' @importFrom stats start end window ts lag frequency time Box.test coef
#' @importFrom zoo na.trim
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @keywords internal
.BayesFitTFP <- function(model, prior = initializePrior(model), R = 10000, burnin = ceiling(R / 10),
                         thin = 1, HPDIprob = 0.85, FUN = mean, MLEfit = NULL) {

  # check input
  .checkBayesInput(
    model = model, type = "tfp", prior = prior, R = R, burnin = burnin,
    thin = thin, HPDIprob = HPDIprob, FUN = FUN, MLEfit = MLEfit
  )

  # attributes
  trend <- attr(model, "trend")
  cycle <- attr(model, "cycle")
  cycleLag <- attr(model, "cubs")$cycleLag
  cubsAR <- attr(model, "cubs")$cubsAR
  errorARMA <- attr(model, "cubs")$errorARMA
  exoNames <- attr(model, "cubs")$exoNames
  freq <- ifelse(attr(model, "period")$frequency == "quarterly", 4, 1)

  # ----- initial parameters (from MLE)
  message("Inizializing parameters by MLE ...")
  if (!is.null(MLEfit) & !inherits(MLEfit, "TFPfit")) {
    MLEfit <- NULL
    message("The supplied object is not of class `TFPfit', starting MLE.")
  }
  if (is.null(MLEfit)) {
    parRestr <- initializeRestr(model = model, type = "hp")
    f <- suppressWarnings(fit(model = model, parRestr = parRestr))
  } else {
    .checkModelMLEfit(model = model, MLEfit = MLEfit)
    f <- MLEfit
  }
  pars <- f$parameters[sort(rownames(f$parameters)), 1]
  names(pars) <- sort(rownames(f$parameters))

  # get rid of trend variance if present
  if (trend != "RW1") {
    model <- .modifySSSystem(model, signalToNoise = NULL)
  }
  .checkModelPrior(model = model, prior = prior)

  # parameter names, locations and constraints
  loc <- model$loc
  nPar <- dim(loc)[1]

  # ----- prior

  # check whether prior is consistent with specified model
  .checkPrior(model = model, prior = prior)

  # convert prior mean and standard deviation to actual parameters
  # gamma and beta needs to be adjusted
  distrPar <- .priorMSd2Parameter(
    prior = cbind(prior$cycle, prior$cubs, prior$trend)[1:2, ],
    restr = cbind(prior$cycle, prior$cubs, prior$trend)[3:4, ],
    namesInvGammaDistr = loc$varName[loc$distribution == "invgamma"],
    namesNormalDistr = loc$varName[loc$distribution == "normal"],
    namesBetaDistr = loc$varName[loc$distribution == "beta"]
  )

  # model dimensions
  N <- dim(model$SSModel$y)[2]
  Tt <- dim(model$SSModel$y)[1]
  k <- dim(model$SSModel$T)[1]

  # update model
  SSModel <- .updateSSSystem(
    pars = pars, SSModel = model$SSModel, loc = loc, cycle = cycle,
    trend = trend, errorARMA = errorARMA, bayes = TRUE
  )

  # ----- assign gibbs functions and function input

  ####### order of names matters
  tmp <- .assignGibbsFUN(
    loc = loc, type = "tfp", trend = trend, cycle = cycle, cubsAR = cubsAR,
    cycleLag = cycleLag, errorARMA = errorARMA
  )
  FUNtrend <- tmp$FUN$trend
  FUNcycle <- tmp$FUN$cycle
  FUNcubs <- tmp$FUN$cubs
  names <- tmp$names

  # ----- Gibbs procedure

  # initialization
  state <- array(0, dim = c(Tt, k, R / thin))
  obsFitted <- array(0, dim = c(Tt, 2, R / thin))
  param <- array(0, c(R / thin, nPar))
  accept <- rep(0, R / thin)
  colnames(param) <- names(pars)
  count <- 0
  aProb <- 0

  # print details and progress
  message(
    paste0(
      "Bayesian Estimation of TFP model\n",
      "\tBurn-in period \t\t\t", burnin, "\n",
      "\tNumber of repititions \t\t", R, "\n",
      "\tSkipped draws (thinning) \t", thin - 1
    )
  )

  # loop
  message("Obtaining draws ...")
  pb <- utils::txtProgressBar(min = 0, max = R, style = 3)
  for (r in 1:R) {
    if (r %% 100 == 0) {
      utils::setTxtProgressBar(pb, r)
    }

    # ----- Step 1 ----- states
    # apply simulation smoother
    stateSmoothed <- ts(simulateSSM(SSModel, type = "states", nsim = 1)[,,1],
                        start = start(SSModel$y), frequency = freq)
    stateSmoothed_werror <- stateSmoothed
    stateSmoothed_werror[,grepl("E2error", colnames(stateSmoothed))] <- 0
    obs <- matmult3d(a = stateSmoothed_werror, b = SSModel$Z)
    
    # ----- Step 2 ----- trend
    Ytrend <- stateSmoothed[, "trend"]
    if (trend == "DT") {
      tmp <- FUNtrend(Y = Ytrend, par = pars, distr = distrPar, varName = names$trend)
      pars[names$trend] <- tmp$par[names$trend]
      aProb <- tmp$aProb
    } else if (trend == "RW1" | trend == "RW2") {
      n <- length(names$trend)
      tmp <- .postRW(
        Y = Ytrend,
        sigmaDistr = distrPar[, names$trend[n], drop = FALSE],
        sigmaLast = pars[names$trend[n]],
        muDistr = distrPar[, names$trend[n - 1], drop = FALSE]
      )
      pars[names$trend] <- tmp$par[names$trend]
    }

    # ----- Step 3 ----- cycle
    Ycycle <- stateSmoothed[, "cycle"]
    tmp <- FUNcycle(Y = Ycycle, parLast = pars, parDistr = distrPar, varNames = names$cycle)
    pars[names$cycle] <- tmp$par[names$cycle]

    # ----- Step 4 ----- cubs equation
    XYcubs <- .getXYcubs(
      stateSmoothed = stateSmoothed, model = model, cycleLag = cycleLag,
      cubsAR = cubsAR, names = names$cubs$betaNames, trim = TRUE
    )
    Yt <- XYcubs$Y
    X <- XYcubs$X
    if (cycle == "RAR2") {
      phiC <- .RAR2transform(par = pars[names$cycle[1:2]], phi2Atau = FALSE)
    } else {
      phiC <- pars[names$cycle][1:(length(names$cycle) - 1)]
    }
    tmp <- FUNcubs(
      Y = Yt, X = X,
      p = cubsAR,
      pk = cycleLag,
      pa = errorARMA[1],
      betaLast = pars[names$cubs$betaNames],
      sigmaLast = pars[names$cubs$varNames],
      betaDistr = distrPar[, names$cubs$betaNames, drop = FALSE],
      sigmaDistr = distrPar[, names$cubs$varNames, drop = FALSE],
      phiLast = pars[names$cubs$phiENames],
      phiDistr = distrPar[, names$cubs$phiENames, drop = FALSE],
      phiC = phiC,
      sigmaC = pars[loc$varName[loc$sysMatrix == "Q" & loc$variableRow == "cycle"]]
    )
    pars[unlist(names$cubs)] <- tmp$par[unlist(names$cubs)]

    # save draws and thin
    if (r %% thin == 0) {
      count <- count + 1
      state[, , count] <- stateSmoothed
      # state_f[, , count] <- stateFiltered
      param[count, ] <- pars
      accept[count] <- aProb
      obsFitted[, , count] <- obs
    }

    # update model parameters
    SSModel <- .updateSSSystem(
      pars = pars, SSModel = SSModel, loc = loc, cycle = cycle,
      trend = trend, errorARMA = errorARMA, bayes = TRUE
    )
  }
  close(pb)
  message("Done.")


  # get rid of burn in phase
  mcmc <- list(
    states = state,
    # states_f = state_f,
    parameters = param,
    fitted = obsFitted
  )
  state <- state[, , (burnin / thin + 1):(R / thin)]
  obsFitted <- obsFitted[, , (burnin / thin + 1):(R / thin)]
  param <- param[(burnin / thin + 1):(R / thin), ]
  accept <- accept[(burnin / thin + 1):(R / thin)]
  colnames(state) <- colnames(stateSmoothed)
  
  # ----- estimated parameters
  paramEstim <- apply(param, 2, FUN)
  dfRes <- mcmcSummary(x = param, HPDIprob = HPDIprob)
  SSModel <- .updateSSSystem(
    pars = paramEstim, SSModel = SSModel, loc = loc, cycle = cycle,
    trend = trend, errorARMA = errorARMA, bayes = TRUE
  )
  
  # ----- estimated states
  tslRes <- .SSresultsBayesian(model = model$SSModel, HPDIprob = HPDIprob, FUN = FUN,
                               state = state, obsFitted = obsFitted)
  start <- start(stateSmoothed)
  freq <- frequency(stateSmoothed)
  # tfp trend
  tsTmp <- exp(state[, "trend", ] / 100)
  tslRes$tfpTrendSummary <- ts(mcmcSummary(x = t(tsTmp), HPDIprob = HPDIprob), start = start, frequency = freq)
  tslRes$tfpTrend <- tslRes$tfpTrendSummary[, firstLetterUp(FUN)]
  # tfp trend growth rate
  tsTmp <- apply(exp(state[, "trend", ] / 100), 2, function(x) (x[2:length(x)] / x[1:(length(x) - 1)]) - 1)
  tslRes$tfpTrendGrowthSummary <- ts(mcmcSummary(x = t(tsTmp), HPDIprob = HPDIprob), start = start + c(0, 1), frequency = freq)
  tslRes$tfpTrendGrowth <- tslRes$tfpTrendGrowthSummary[, firstLetterUp(FUN)]
  # fitted cubs equation
  tsTmp <- obsFitted[, 2, ]
  tslRes$cubsFittedSummary <- ts(mcmcSummary(x = t(tsTmp), HPDIprob = HPDIprob), start = start, frequency = freq)
  tslRes$cubsFitted <- tslRes$cubsFittedSummary[, firstLetterUp(FUN)]

  # ----- convergence diagnostics: Geweke test
  .printGeweke(tsl = tslRes$stateSmoothedSummary, df = dfRes, alpha = 0.05)

  # ----- output
  
  # KFS out and fit output
  SSMout <- within(list(), {
    model <- model$SSModel
    dims <-  list(n = dim(state)[1])
    alphahat <- tslRes$stateSmoothed
    V <- array(apply(state, 3, var), c(rep(dim(state)[2], 2),dim(state)[3]))
    Ptt <- array(apply(state, 3, var), c(rep(dim(state)[2], 2),dim(state)[3]))
  })
  SSMfit <- list(model = SSModel)  

  # order parameters
  dfRes <- dfRes[order(rownames(dfRes)), ]

  # model fit
  info <- list()
  nameCoef <- firstLetterUp(FUN)
  info[["signal-to-noise"]] <- sum(dfRes[loc$varName[loc$equation == "trend" & loc$sysMatrix == "Q"], nameCoef], na.rm = TRUE) / dfRes["cSigma", nameCoef]
  info$R2 <- 1 - sum((tslRes$cubsFitted - SSModel$y[, 2])^2, na.rm = TRUE) / sum((SSModel$y[, 2] - mean(SSModel$y[, 2], na.rm = TRUE))^2, na.rm = TRUE) # 1 - sum((residuals(out)[,"cubs"])^2, na.rm = TRUE) / sum( (out$model$y[,"cubs"] - mean(out$model$y[,"cubs"], na.rm = TRUE) )^2, na.rm = TRUE)
  info$MRMSE <- mean(apply(apply(mcmc$fitted[, 2, , drop = TRUE], 
                                 2, function(x) x - SSModel$y[, 2]), 
                           2, function(x) sqrt(mean(x^2, na.rm = TRUE))))
  
  # ----- fitted bayesian tfp object

  TFPfit <- list(
    model = model,
    tsl = tslRes,
    SSMfit = SSMfit,
    SSMout = SSMout,
    parameters = dfRes,
    parametersSampled = mcmc$param,
    statesSampled = mcmc$state,
    mcmc = mcmc,
    prior = prior,
    fit = info,
    MLE = f
  )
  class(TFPfit) <- c("TFPfit", "fit")
  attr(TFPfit, "method") <- "bayesian"
  attr(TFPfit, "R") <- R
  attr(TFPfit, "burnin") <- burnin
  attr(TFPfit, "thin") <- thin
  attr(TFPfit, "HPDIprob") <- HPDIprob
  attr(TFPfit, "FUN") <- FUN

  invisible(TFPfit)
}
# -------------------------------------------------------------------------------------------

#' defines Y and X in the CUBS equation.
#' @param stateSmoothed The smoothed states.
#' @param names The names of the columns.
#' @param trim A logical indicating whether NAs should be trimmed.
#' @inheritParams TFPmodel
#' @inheritParams fit.TFPmodel
#' @keywords internal
.getXYcubs <- function(stateSmoothed, model, cycleLag, cubsAR, names, trim = FALSE) {
  loc <- model$loc
  Ycubs <- model$SSModel$y[, "cubs"]
  tmp <- cbind(Ycubs, rep(1, length(Ycubs)))
  if (cubsAR > 0) {
    tmp <- cbind(tmp, do.call(cbind, model$tsl[paste0("cubsAR", 1:cubsAR)]))
  }
  # cycle
  tmp <- cbind(tmp, stateSmoothed[, loc$variableRow[loc$equation == "cubs" & grepl("cycle", loc$variableRow)]])

  if (trim) {
    tmp <- na.trim(tmp)
  }
  Y <- tmp[, 1]
  X <- tmp[, -1]
  colnames(X) <- names

  res <- list(Y = Y, X = X)
  return(res)
}
