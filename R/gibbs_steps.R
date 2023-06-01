# -------------------------------------------------------------------------------------------

#' Draws from the posterior of the parameters of the cubs equation, conditional on the states.
#'
#' @param Y a \code{Tn x 1} vector.
#' @param X a \code{Tn x n} matrix, includes a constant, the contemporaneous cycle, cycle
#'   lags, and lags of cubs.
#' @param p integer, lag of cubs.
#' @param pk integer, contemporaneous cycle and cycle lags.
#' @param pa integer, autoregressive order of error term.
#' @param betaLast last draw from posterior of beta.
#' @param sigmaLast last draw from posterior of cubs innovation variance.
#' @param betaDistr prior distribution of beta.
#' @param sigmaDistr prior distribution of cubs innovation variance.
#' @param phiLast (optional) vector, last draw from posterior of the autoregressive
#'   parameter of the error term.
#' @param phiDistr (optional) prior of the autoregressive parameter of the error term.
#' @param phiC parameter vector of the cycle equation.
#' @param sigmaC innovation variance of the cycle equation
#'
#' @details The parameter vector beta and the innovation variance are Normal-inverse Gamma
#'   distributed. A draw from their posterior is obtained by conjugacy.
#' @details If there are additional lags of the cycle or cubs, conjugacy does not apply
#'   since the starting values are not given. In this case, a Metropolis-Hasting step is
#'   implemented.
#' @details If the error term is an AR(1) or AR(2) process, an additional Gibbs step draws
#'   from the posterior of the autoregressive parameter, given all other parameters.
#'
#' @importFrom stats runif
#' @keywords internal
.gibbsStep2Eq <- function(Y, X, p, pk, pa, betaLast, sigmaLast, betaDistr, sigmaDistr, phiLast = NULL, phiDistr = NULL, phiC, sigmaC) {
  # n0 shape prior
  # s0 mean of variance prior
  # prior
  # beta ~ N(beta0, Q0^-1)
  # h ~ gamma(s0, nu0)
  # posterior (NIG, conjugate)
  # h ~ gamma(s0^-2, nu0)
  # beta | h ~ N(beta0, Q0 %*% h^-1)

  # notation (see Bauwens et. al 1999)
  # prior: beta, sigma2 ~ NIG(beta0, M0, s0, nu0)
  # posterior: beta, sigma2 ~ NIG(beta*, M*, s*, nu*)
  # shape s = 2 * beta
  # degrees of freedom nu = 2 * alpha

  betaLast <- drop(betaLast)

  # dimension
  Tn <- length(Y)
  n <- dim(X)[2]
  pexo <- n - p - length(pk) - 1
  q <- max(pk) + p

  # Normal-inverse gamma posterior
  tmp <- .postNIG(Y = Y[(q + 1):Tn], X = X[(q + 1):Tn, ], betaLast = betaLast, betaDistr = betaDistr, sigmaDistr = sigmaDistr)
  betar <- tmp$beta
  sigmar <- tmp$sigma

  alpha <- NULL

  # which cycle lags
  pk_all <- 0:max(pk)
  pk_missing <- pk_all[!(pk_all %in% pk)]

  # additional cycle lags or lags of cubs
  # if (n - pexo > 2) {
  if (max(pk) > 0 | p > 0) {

    # metropolis hasting step
    mu <- betar[1]
    muL <- betaLast[1]
    beta <- rep(0, max(pk) + 1)
    betaL <- rep(0, max(pk) + 1)
    beta[pk + 1] <- betar[(2 + p + pexo):length(betar)]
    betaL[pk + 1] <- betaLast[(2 + p + pexo):length(betar)]
    Xcubs <- NULL
    phi <- NULL
    phiL <- NULL
    if (p > 0) {
      phi <- betar[2:(1 + p)]
      phiL <- betaLast[2:(1 + p)]
      Xcubs <- matrix(X[, 2:(1 + p)], Tn, p)
    }

    # acceptance probability
    Xc <- ts(matrix(0, nrow(X), max(pk) + 1), start = start(X), frequency = frequency(X))
    Xc[, pk + 1] <- X[, (2 + p + pexo):length(betar)]
    prod1 <- FUNcov(
      Xcubs = Xcubs, Xc = Xc, Y = Y,
      mu = mu, phi = phi, beta = beta, sigma = sigmar, phiC = phiC, sigmaC = sigmaC
    )
    prod2 <- FUNcov(
      Xcubs = Xcubs, Xc = Xc, Y = Y,
      mu = muL, phi = phiL, beta = betaL, sigma = sigmaLast, phiC = phiC, sigmaC = sigmaC
    )
    alpha <- min(1, prod1 / prod2)

    # accept/reject the current draw
    if (runif(1, 0, 1) > alpha) {
      betar <- betaLast
      sigmar <- sigmaLast
    }
  }

  # result
  res <- c(betar, sigmar)

  # 2nd equation error is AR(pa)
  if (pa == 1) {
    phiDistr <- as.matrix(phiDistr)
    eps <- Y - betar[1] - betar[2] * X[, 2]
    phir <- .postAR1(
      Y = eps, phi = phiLast, phi0 = phiDistr[1, ], Q0 = 1 / phiDistr[2, ],
      sigma = as.numeric(sigmar), lb = phiDistr[3, ], ub = phiDistr[4, ]
    )
    res <- c(betar, phir, sigmar)
  } else if (pa == 2) {
    phiDistr <- as.matrix(phiDistr)
    eps <- Y - betar[1] - betar[2] * X[, 2]
    tmp <- .postARp(
      Y = eps, phi = phiLast, phi0 = phiDistr[1, ], Q0 = solve(diag(phiDistr[2, ])),
      sigma = as.numeric(sigmar), lb = phiDistr[3, ], ub = phiDistr[4, ]
    )
    phir <- tmp$phi
    res <- c(betar, phir, sigmar)
  }

  names(res) <- c(names(betaLast), names(phiLast), names(sigmaLast))
  res <- list(
    "par" = res,
    "aProb" = alpha
  )
  return(res)
}


# -------------------------------------------------------------------------------------------

#' Computes mean and variance of the part of the posterior distribution that relies on
#' starting values. It then computes the density of the first p observations of Y.
#'
#' @param Xcubs matrix containing lags of cubs.
#' @param Xc matrix containing the contemporaneous cycle and lags thereof..
#' @param mu constant parameter.
#' @param phi cubs lag coeefficients.
#' @param beta cycle coefficients.
#' @param sigma innovation variance.
#' @param phiC cycle process parameter vector.
#' @param sigmaC cycle innovation variance
#' @inheritParams .gibbsStep2Eq
#' @keywords internal
FUNcov <- function(Xcubs, Xc, Y, mu, phi, beta, sigma, phiC, sigmaC) {

  # dimensions
  pk <- length(beta) - 1
  p <- length(phi)
  q <- max(pk, p)

  muhat <- mu
  if (p > 0) {
    muhat <- muhat / (1 - sum(phi))
  }

  if (p > 0 | pk > 0) {
    # unconditional cubs variance, unconditional cubs, cycle covariance
    covTmp <- .covCUBS(mu = mu, phi = phi, beta = beta, sigma = sigma, phiC = phiC, sigmaC = sigmaC)
    covCubs <- covTmp$covCubs
  }

  if (det(covCubs[1:q, 1:q, drop = FALSE]) > 0) {
    res <- det(covCubs[1:q, 1:q, drop = FALSE])^(-1 / 2) * exp(-1 / 2 * t(Y[1:q] - muhat) %*% solve(covCubs[1:q, 1:q, drop = FALSE]) %*% (Y[1:q] - muhat))
  } else {
    res <- 1e-6
    warning("Sigma of ind_1, ..., ind_q is singular. Rejecting draw.")
  }

  return(res)
}


# -------------------------------------------------------------------------------------------

#' Draws from the posterior of the parameters of the RAR2 cycle equation, conditional on the
#' states.
#'
#' @param Y a \code{Tn x 1} vector.
#' @param parLast A \code{3 x 1} vector containing the last draw for amplitude,
#'   mean cycle periodicity, innovation variance (in that order).
#' @param parDistr A \code{4 x 3} matrix with prior distribution and box constraints for
#'   the parameters of each variable (the order of the columns is as for parLast). In each
#'   column, the first two entries contain the prior hyperparameters and the last two entries
#'   the upper and lower bound.
#' @param varNames A vector with parameter names in the correct order, i.e., amplitude,
#'   mean cycle periodicity, variance.
#'
#' @details The parameters \eqn{A} and \eqn{\tau} are Beta distributed and the variance
#'   \eqn{\sigma^2} is Gamma distributed.
#' @details The posterior is not available in closed form, instead it is obtained via an ARMS
#'   step.
#'
#' @return A list with the draws.
#'
#' @importFrom dlm arms
#' @importFrom stats runif rgamma
#' @keywords internal
.gibbsStepRAR2 <- function(Y, parLast, parDistr, varNames) {
  # draws for A and tau from posterior are obtained via ARMS steps
  # draw for sigma from posterior is obtained through IG framework

  # varName must be in the correct order
  ALast <- parLast[varNames[1]]
  tauLast <- parLast[varNames[2]]
  sigmaLast <- parLast[varNames[3]]
  ADistr <- parDistr[, varNames[1]]
  tauDistr <- parDistr[, varNames[2]]
  sigmaDistr <- parDistr[, varNames[3]]

  # dimensions
  n <- length(Y)
  p <- 2

  # convert to phi
  phi1 <- 2 * ALast * cos(2 * pi / tauLast)
  phi2 <- -ALast^2
  phi <- c(phi1, phi2)

  # arms step
  # if target density is log-concave, then metropolis step will always accept -> no rejections
  FUNdensity <- function(x, Y, phi1, phi2, sigma, shape, lb, ub) {
    n <- length(Y)
    covC <- .covAR(k = 2, phi = c(phi1, phi2), sigma = sigma)
    prod1 <- 1 / ((2 * pi)^(n - 2) * sigma^(1 / 2)) * exp(-1 / 2 * 1 / sigma * sum((Y[3:n] - phi1 * Y[2:(n - 1)] - phi2 * Y[1:(n - 2)])^2))
    prod2 <- 1 / ((2 * pi)^2 * det(covC)^(1 / 2)) * exp(-1 / 2 * t(Y[1:2]) %*% solve(covC) %*% Y[1:2])
    # scaled beta distribution
    prod3 <- ((x - lb)^(shape[1] - 1) * (ub - x)^(shape[2] - 1)) / ((ub - lb)^(shape[1] + shape[2] - 1) * beta(shape[1], shape[2]))
    prod <- log(prod1 * prod2 * prod3)
    return(prod)
  }

  # ----- A
  lb <- ADistr[3]
  ub <- ADistr[4]
  shape0 <- ADistr[1:2]
  FUNrestriction <- function(x, Y, phi1, phi2, sigma, shape, lb, ub) (x > lb) * (x < ub)
  Ar <- dlm::arms(
    y.start = ALast, myldens = FUNdensity, indFunc = FUNrestriction, n.sample = 1,
    Y = Y, phi1 = phi1, phi2 = phi2, sigma = sigmaLast, shape = shape0, lb = lb, ub = ub
  )

  # update phi
  phi1 <- 2 * Ar * cos(2 * pi / tauLast)
  phi2 <- -Ar^2
  phi <- c(phi1, phi2)

  # ----- tau
  lb <- tauDistr[3]
  ub <- tauDistr[4]
  shape0 <- tauDistr[1:2]
  FUNrestriction <- function(x, Y, phi1, phi2, sigma, shape, lb, ub) (x > lb) * (x < ub)
  taur <- dlm::arms(
    y.start = tauLast, myldens = FUNdensity, indFunc = FUNrestriction, n.sample = 1,
    Y = Y, phi1 = phi1, phi2 = phi2, sigma = sigmaLast, shape = shape0, lb = lb, ub = ub
  )

  # update phi and variance of (c1, c2)
  phi1 <- 2 * Ar * cos(2 * pi / taur)
  phi2 <- -Ar^2
  phi <- c(phi1, phi2)
  covC <- .covAR(k = p, phi = phi, sigma = sigmaLast)

  # ----- sigma (IG conjugate framework)
  s0 <- sigmaDistr[1]
  nu0 <- sigmaDistr[2]
  lb <- max(0, sigmaDistr[3])
  ub <- sigmaDistr[4]
  RSS <- sum((Y[3:n] - phi1 * Y[2:(n - 1)] - phi2 * Y[1:(n - 2)])^2)
  nustar <- nu0 + n
  sstar <- s0 + t(Y[1:2]) %*% solve(covC / sigmaLast) %*% Y[1:2] + RSS
  # x in [lb, ub] <=> 1/x in [1/ub, 1/lb]
  tmp <- pgamma(1 / ub, shape = nustar / 2, rate = sstar / 2) +
    runif(1, 0, 1) * (pgamma(1 / lb, shape = nustar / 2, rate = sstar / 2) - pgamma(1 / ub, shape = nustar / 2, rate = sstar / 2))
  sigmar <- 1 / qgamma(tmp, shape = nustar / 2, rate = sstar / 2)

  # return
  parr <- c(Ar, taur, sigmar)
  names(parr) <- varNames
  res <- list("par" = parr)
  return(res)
}

# -------------------------------------------------------------------------------------------

#' Draws from the posterior of the parameters of the AR(p), \code{p = 1,2} cycle equation,
#' conditional on the states.
#'
#' @param Y a \code{Tn x 1} vector.
#' @param parLast A \code{(p + 1) x 1} vector containing the last draw for the
#'    autoregressive coefficients and the innovation variance (in that order).
#' @param parDistr A \code{4 x (p + 1)} matrix with prior distribution and box constraints for
#'   the parameters of each variable (the order of the columns is as for parLast). In each
#'   column, the first two entries contain the prior hyperparameters and the last two entries
#'   the upper and lower bound.
#' @param varNames A vector with parameter names in the correct order, i.e., autoregressive
#'   coefficients, variance.
#'
#' @details The autoregressive parameter and the innovation variance are drawn sequentially.
#' @details If the cycle is AR(1) process, the posterior is obtained by conjugacy. If it is
#'   an AR(2) process, a Metropolis-Hastings step is implemented.
#' @details Conditional on the autoregressive parameters, the innovation variance is drawn
#'   from the Inverse Gamma posterior which is obtained by conjugacy.
#'
#' @importFrom stats runif rgamma pgamma
#' @keywords internal
.gibbsStepAR <- function(Y, parLast, parDistr, varNames) {

  # dimension
  Tn <- length(Y)
  p <- length(varNames) - 1

  # ----- phi
  phiLast <- parLast[varNames[1:p]]
  sigmaLast <- parLast[varNames[p + 1]]
  phiDistr <- as.matrix(parDistr[, varNames[1:p]])
  sigmaDistr <- as.matrix(parDistr[, varNames[p + 1]])

  if (p == 1) {
    phir <- .postAR1(
      Y = Y, phi = phiLast, phi0 = phiDistr[1, ], Q0 = 1 / phiDistr[2, ],
      sigma = as.numeric(sigmaLast), lb = phiDistr[3, ], ub = phiDistr[4, ]
    )
    covCstar <- .covAR(k = 2, phi = phir, sigma = sigmaLast)
  } else if (p == 2) {
    tmp <- .postARp(
      Y = Y, phi = phiLast, phi0 = phiDistr[1, ], Q0 = solve(diag(phiDistr[2, ])),
      sigma = as.numeric(sigmaLast), lb = phiDistr[3, ], ub = phiDistr[4, ]
    )
    phir <- tmp$phi
    covCstar <- tmp$sigmaY12
  }

  # ----- sigma (IG conjugate framework)
  X <- sapply(1:p, function(x) {
    Y[seq((p + 1 - x), length(Y) - x)]
  })
  s0 <- sigmaDistr[1]
  nu0 <- sigmaDistr[2]
  lb <- max(0, sigmaDistr[3])
  ub <- sigmaDistr[4]
  RSS <- sum((Y[(p + 1):Tn] - X %*% phir)^2)
  nustar <- nu0 + Tn
  if (p == 1) {
    sstar <- s0 + t(Y[1]) %*% solve(covCstar[1, 1] / sigmaLast) %*% Y[1] + RSS
  } else if (p == 2) {
    sstar <- s0 + t(Y[1:2]) %*% solve(covCstar / sigmaLast) %*% Y[1:2] + RSS
  }
  # x in [lb, ub] <=> 1/x in [1/ub, 1/lb]
  tmp <- pgamma(1 / ub, shape = nustar / 2, rate = sstar / 2) +
    runif(1, 0, 1) * (pgamma(1 / lb, shape = nustar / 2, rate = sstar / 2) - pgamma(1 / ub, shape = nustar / 2, rate = sstar / 2))
  sigmar <- 1 / qgamma(tmp, shape = nustar / 2, rate = sstar / 2)

  # return
  par <- c(phir, sigmar)
  names(par) <- varNames
  res <- list("par" = par)
  return(res)
}

# -------------------------------------------------------------------------------------------

#' Draws from the posterior of the parameters of the damped trend equation, conditional on
#' the states.
#'
#' @param Y A \code{Tn x 1} vector.
#' @param par A \code{3 x 1} vector with parameters.
#' @param distr A \code{4 x 3} matrix with prior parameters and box constraints.
#' @param varName A code{3 x 1} vector with parameter names in the correct order, i.e.,
#'   mean reversion, autoregressive parameter, variance.
#'
#' @details The three parameters are drawn sequentially in a Gibbs procedure. (conditional
#'   on the two other parameters).
#' @details The parameter \eqn{\omega} is drawn from a normal posterior which is obtained
#'   by conjugancy.
#' @details The autoregressive parameter \eqn{\phi} is drawn via a Metropolis-Hastings step.
#' @details The innovation variance is drwan from the Inverse-Gamma distribution which
#'   obtained by conjugacy.
#'
#' @importFrom stats runif pgamma
#' @keywords internal
.gibbsStepDT <- function(Y, par, distr, varName) {

  # varName must be in the correct order
  phiLast <- par[varName[2]]
  sigmaLast <- par[varName[3]]

  omegaDistr <- distr[, varName[1]]
  phiDistr <- distr[, varName[2]]
  sigmaDistr <- distr[, varName[3]]

  # number of observations
  Ty <- length(Y)

  # compute difference process (stationary)
  pDiff <- Y[2:Ty] - Y[1:(Ty - 1)]
  Tp <- length(pDiff)

  # --- omega
  omega0 <- omegaDistr[1]
  omegaQ0 <- 1 / omegaDistr[2]
  lb <- omegaDistr[3]
  ub <- omegaDistr[4]
  omegaQstar <- solve(sigmaLast) * ((Ty - 2) * (1 - phiLast)^2 + 1 - phiLast^2) + omegaQ0
  omegastar <- solve(omegaQstar) * (solve(sigmaLast) * ((1 - phiLast) * (pDiff[1] + pDiff[Tp]) + (1 - phiLast)^2 * sum(pDiff[2:Tp])) +
    omega0 * omegaQ0)

  # Draw omega from density (procedure such that ub and lb hold)
  tmp <- pnorm(lb, mean = omegastar, sd = sqrt(solve(omegaQstar))) +
    runif(1, 0, 1) * (pnorm(ub, mean = omegastar, sd = sqrt(solve(omegaQstar))) -
      pnorm(lb, mean = omegastar, sd = sqrt(solve(omegaQstar))))
  omegar <- qnorm(tmp, mean = omegastar, sd = sqrt(solve(omegaQstar)))

  # --- phi (MH step)
  phi0 <- phiDistr[1]
  phiQ0 <- 1 / phiDistr[2]
  lb <- phiDistr[3]
  ub <- phiDistr[4]
  phiQstar <- sum((pDiff[1:(Tp - 1)] - omegar)^2) * solve(sigmaLast) + phiQ0
  phistar <- solve(phiQstar) * (sum((pDiff[2:Tp] - omegar) * (pDiff[1:(Tp - 1)] - omegar)) * solve(sigmaLast)
    + phi0 * phiQ0)

  # Draw theta from proposal density (procedure such that ub and lb hold)
  tmp <- pnorm(lb, mean = phistar, sd = sqrt(solve(phiQstar))) +
    runif(1, 0, 1) * (pnorm(ub, mean = phistar, sd = sqrt(solve(phiQstar))) -
      pnorm(lb, mean = phistar, sd = sqrt(solve(phiQstar))))
  phiprop <- qnorm(tmp, mean = phistar, sd = sqrt(solve(phiQstar)))

  # Compute acceptance probability
  alpha_ic <- min(1, exp((phiprop^2 - phiLast^2) * (pDiff[1] - omegar)^2 / (2 * sigmaLast)))

  # Accept-Reject the current draw
  if (runif(1, 0, 1) < alpha_ic) {
    phir <- phiprop
  } else {
    phir <- phiLast
  }

  # --- variance
  s0 <- sigmaDistr[1]
  nu0 <- sigmaDistr[2]
  lb <- max(0, sigmaDistr[3])
  ub <- sigmaDistr[4]
  nustar <- nu0 + Ty - 1
  sstar <- s0 + (1 - phir^2) * (pDiff[1] - omegar)^2 + sum((pDiff[2:Tp] - omegar * (1 - phir^2) - phir * pDiff[1:(Tp - 1)])^2)

  # Draw sigma from density (procedure such that ub and lb hold)
  # x in [lb, ub] <=> 1/x in [1/ub, 1/lb]
  tmp <- pgamma(1 / ub, shape = nustar / 2, rate = sstar / 2) +
    runif(1, 0, 1) * (pgamma(1 / lb, shape = nustar / 2, rate = sstar / 2) - pgamma(1 / ub, shape = nustar / 2, rate = sstar / 2))
  sigmar <- 1 / qgamma(tmp, shape = nustar / 2, rate = sstar / 2)

  parr <- c(omegar, phir, sigmar)
  names(parr) <- varName
  res <- list(
    "par" = parr,
    "aProb" = alpha_ic
  )
  return(res)
}

# -------------------------------------------------------------------------------------------

#' Assigns the appropriate function and its input variables for the Gibbs procedure.
#'
#' @param exoNames A character vector containing the names of the exogenous variables.
#' @inheritParams TFPmodel
#' @inheritParams NAWRUmodel
#' @inheritParams .updateSSSystem
#' @keywords internal
.assignGibbsFUN <- function(loc, type, trend, cycle, cubsAR, cycleLag, errorARMA, exoNames = NULL) {

  # initialize
  names <- list()
  FUN <- list()

  # trend
  names$trend <- loc$varName[loc$equation == "trend"]
  if (trend == "DT") {
    FUN$trend <- .gibbsStepDT
    names$trend <- c(
      loc$varName[loc$equation == "trend" & loc$variableRow == "const"],
      loc$varName[loc$equation == "trend" & loc$variableRow == "trendDrift" & loc$sysMatrix != "Q"],
      loc$varName[loc$equation == "trend" & loc$variableRow == "trendDrift" & loc$sysMatrix == "Q"]
    )
  } else if (trend == "RW1" | trend == "RW2") {
    FUN$trend <- .postRW
    names$trend <- c(
      loc$varName[loc$equation == "trend" & loc$sysMatrix != "Q"],
      loc$varName[loc$equation == "trend" & loc$sysMatrix == "Q"]
    )
  }

  # cycle
  names$cycle <- c(
    loc$varName[loc$equation == "cycle" & loc$variableRow == "cycle" & loc$sysMatrix != "Q"],
    loc$varName[loc$equation == "cycle" & loc$variableRow == "cycleLag1" & loc$sysMatrix != "Q"],
    loc$varName[loc$equation == "cycle" & loc$variableRow == "cycle" & loc$sysMatrix == "Q"]
  )
  if (cycle == "RAR2") {
    FUN$cycle <- .gibbsStepRAR2
  } else {
    FUN$cycle <- .gibbsStepAR
  }

  if (type == "tfp") {
    # cubs
    names$cubs$betaNames <- rep("", 2 + cubsAR + (length(cycleLag) - 1))
    names$cubs$phiENames <- rep("", errorARMA[1])
    names$cubs$varNames <- rep("", 1)
    names$cubs$betaNames[1] <- loc$varName[loc$equation == "cubs" & loc$variableRow == "const"]
    if (cubsAR > 0) {
      names$cubs$betaNames[2:(1 + cubsAR)] <- loc$varName[loc$equation == "cubs" & grepl("cubsAR", loc$variant)]
    }
    names$cubs$betaNames[(2 + cubsAR):(2 + cubsAR + length(cycleLag) - 1)] <- loc$varName[loc$equation == "cubs" & (grepl("cycleLag", loc$variableRow) | grepl("cycle", loc$variableRow))]
    names$cubs$phiENames <- loc$varName[loc$equation == "E2error" & grepl("errorAR", loc$variant)]
    names$cubs$varNames <- loc$varName[loc$equation == "E2error" & loc$sysMatrix == "Q"]
    FUN$cubs <- .gibbsStep2Eq
  } else if (type == "nawru") {
    # pcInd
    pcIndAR <- cubsAR
    nExo <- length(exoNames)
    names$pcInd$betaNames <- rep("", 2 + pcIndAR + (length(cycleLag) - 1))
    names$pcInd$phiENames <- rep("", errorARMA[1])
    names$pcInd$varNames <- rep("", 1)
    names$pcInd$betaNames[1] <- loc$varName[loc$equation == "pcInd" & loc$variableRow == "const"]
    if (nExo > 0) {
      names$pcInd$betaNames[2:(1 + nExo)] <- exoNames
    }
    if (pcIndAR > 0) {
      names$pcInd$betaNames[(2 + nExo):(1 + nExo + pcIndAR)] <- loc$varName[loc$equation == "pcInd" & grepl("pcIndAR", loc$variant)]
    }
    names$pcInd$betaNames[(2 + nExo + pcIndAR):(2 + nExo + pcIndAR + length(cycleLag) - 1)] <- loc$varName[loc$equation == "pcInd" & (grepl("cycleLag", loc$variableRow) | grepl("cycle", loc$variableRow))]
    names$pcInd$phiENames <- loc$varName[loc$equation == "E2error" & grepl("errorAR", loc$variant)]
    names$pcInd$varNames <- loc$varName[loc$equation == "E2error" & loc$sysMatrix == "Q"]
    FUN$pcInd <- .gibbsStep2Eq
  }

  # return
  result <- list(FUN = FUN, names = names)
  return(result)
}
