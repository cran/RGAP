
# -------------------------------------------------------------------------------------------

#' Estimation of a \code{TFPmodel}
#'
#' @description Estimates a two-dimensional state-space model and performs filtering and
#'   smoothing to obtain the TFP trend using either maximum likelihood estimation or
#'   bayesian methods.
#'
#' @param model An object of class TFPmodel.
#' @param parRestr A list of matrices containing the parameter restrictions for the cycle,
#'   trend, and the CUBS equation. Each matrix contains the lower and upper bound of the
#'   involved parameters. \code{NA} implies that no restriction is present. Autoregressive
#'   parameters are automatically restricted to the stationary region unless box constraints
#'   are specified. By default, \code{parRestr} is initialized by the function
#'   \code{initializeRestr(model)}. Only used if \code{method = "MLE"}.
#' @param signalToNoise (Optional) signal to noise ratio. Only used if \code{method = "MLE"}.
#' @param method The estimation method. Options are maximum likelihood estimation \code{"MLE"}
#'   and bayesian estimation \code{"bayesian"}. The default is \code{method = "MLE"}.
#' @param prior A list of matrices with parameters for the prior distribution and box
#'   constraints. By default, \code{prior} is initialized by \code{initializePrior(model)}.
#'   See details. Only used if \code{method = "bayesian"}.
#' @param R An integer specifying the number of MCMC draws. The default is \code{R = 10000}.
#'   Only used if \code{method = "bayesian"}.
#' @param burnin An integer specifying the burn-in phase of the MCMC chain. The default is
#'   \code{burnin = ceiling(R / 10)}. Only used if \code{method = "bayesian"}.
#' @param thin An integer specifying the thinning interval between consecutive draws. The
#'   default is \code{thin = 1}, implying that no draws are dopped. For \code{thin = 2},
#'   every second draw is dropped and so on. Only used if \code{method = "bayesian"}.
#' @param HPDIprob A numeric in the interval \code{(0,1)} specifying the target probability
#'   of the highest posterior density intervals. The default is \code{HPDIprob = 0.9}. Only
#'   used if \code{method = "bayesian"}.
#' @param pointEstimate Posterior distribution's statistic of central tendency. Possible
#'   options are \code{"mean"} and \code{"median"}. The default is \code{pointEstimate = "mean"}.
#'   Only used if \code{method = "bayesian"}.
#' @param MLEfit (Optional) An object of class \code{TFPfit} which is used for
#'   initialization. Only used if \code{method = "bayesian"}.
#' @param control (Optional) A list of control arguments to be passed on to \code{optim}.
#' @param ... additional arguments to be passed to the methods functions.
#'
#' @details The list object \code{prior} contains three list elements \code{cycle},
#'   \code{trend}, and \code{cubs}. Each list element is a \code{4 x n} matrix where \code{n}
#'   denotes the number of parameters involved in the respective equation. The upper two
#'   elements specify the distribution, the lower two parameters specify box constraints.
#'   \code{NA} denotes no constraints. Autoregressive parameters are automatically restricted
#'   to the stationary region unless box constraints are specified. For instance,
#'   \code{prior$cycle[, 1]} contains the mean, standard deviation, lower and upper bound for
#'   the first variable, in that respective order.
#' @details The respective prior distributions are defined through their mean and standard
#'   deviation.
#' @details The Gibbs sampling procedure is as follows. For each \eqn{r = 1, ..., R}
#'   \itemize{\item The states are sampled by running the Kalman filter and smoother
#'            conditional on the parameters of the previous step, \eqn{\theta_{r-1}}
#'            \item Trend equation parameters \eqn{\theta_{trend}}: Conditional on the
#'            states \eqn{\alpha_r}, a draw \eqn{\theta_{trend,k}} is obtained either
#'            by a sequential Gibbs step, a Metropolis Hasting step, or by conjugacy,
#'            depending on the trend model specification.
#'            \item Cycle equation parameters \eqn{\theta_{cycle}}: Conditional on the
#'            states \eqn{\alpha_r}, a draw \eqn{\theta_{cycle,k}} is obtained either
#'            by a sequential Gibbs step, a Metropolis Hasting step, or by conjugacy,
#'            depending on the cycle model specification.
#'            \item CUBS equation parameters \eqn{\theta_{cubs}}: Conditional on the
#'            states \eqn{\alpha_r}, a draw \eqn{\theta_{cubs,k}} is obtained either
#'            by a sequential Gibbs step, a Metropolis Hasting step, a combination
#'            thereof, or by conjugacy, depending on the CUBS equation specification.
#'   }
#'
#' @return For maximum likelihood estimation, an object of class \code{TFPit} containing
#'   the following components:
#'   \item{model}{The input object of class \code{TFPmodel}.}
#'   \item{SSMfit}{The estimation output from the funtcion \code{fitSSM} from \code{KFAS}.}
#'   \item{SSMout}{The filtering and smoothing output from the funtcion \code{KFS} from
#'         \code{KFAS}.}
#'   \item{parameters}{A data frame containing the estimated parameters, including standard
#'         errors, t-statistic, and p-values.}
#'   \item{parRestr}{A list of matrices containing the enforced parameter constraints.}
#'   \item{fit}{A list of model fit criteria (see below).}
#'   \item{call}{Original call to the function. }
#'   The list component \code{fit} contains the following model fit criteria:
#'   \item{loglik}{Log-likelihood function values.}
#'   \item{AIC}{Akaike information criterion.}
#'   \item{BIC}{Bayesian information criterion.}
#'   \item{AIC}{Hannan-Quinn information criterion.}
#'   \item{RMSE}{root mean squared error of the CUBS equation.}
#'   \item{R2}{R squared of the CUBS equation.}
#'   \item{LjungBox}{Ljung-Box test output of the CUBS equation.}
#'   \item{signal-to-noise}{Signal-to-noise ratio.}
#' For bayesian estimation, an object of class \code{TFPfit} containing the following components:
#'   \item{model}{The input object of class \code{TFPmodel}.}
#'   \item{tsl}{A list of time series containing the estimated states.}
#'   \item{parameters}{A data frame containing the estimated parameters, including standard
#'         errors, highest posterior density credible sets.}
#'   \item{prior}{A list of matrices containing the used prior distributions.}
#'   \item{fit}{A list of model fit criteria (see below).}
#'   \item{call}{Original call to the function. }
#'   The list component \code{fit} contains the following model fit criteria:
#'   \item{R2}{R squared of the CUBS equation.}
#'   \item{signal-to-noise}{Signal-to-noise ratio.}
#'   
#' @export
#' @family fitting methods
#' @examples
#' # load data for Italy
#' data("gap")
#' country <- "Italy"
#' tsList <- amecoData2input(gap[[country]])
#' # define tfp model
#' model <- TFPmodel(tsl = tsList, cycle = "RAR2")
#' # initialize parameter restrictions and estimate model
#' parRestr <- initializeRestr(model = model, type = "hp")
#' \donttest{
#' f <- fit(model = model, parRestr = parRestr)
#' }
#' # Bayesian estimation
#' prior <- initializePrior(model = model)
#' \donttest{
#' f <- fit(model = model, method = "bayesian", prior = prior, R = 5000, thin = 2)
#' }
fit.TFPmodel <- function(model, parRestr = initializeRestr(model = model), signalToNoise = NULL,
                   method = "MLE", control = NULL,
                   prior = initializePrior(model), R = 10000, burnin = ceiling(R / 10),
                   thin = 1, HPDIprob = 0.85, pointEstimate = "mean", MLEfit = NULL, ...) {
  # save call
  mc <- match.call(expand.dots = FALSE)

  if (method == "MLE") {
    f <- .MLEfitTFP(
      model = model, parRestr = parRestr, signalToNoise = signalToNoise,
      control = control
    )
  } else if (method == "bayesian") {
    if (!is.null(signalToNoise)) {
      warning("'signalToNoise' not applicable for method = 'bayesian'.")
    }
    f <- .BayesFitTFP(
      model = model, prior = prior, R = R, burnin = burnin,
      thin = thin, HPDIprob = HPDIprob, FUN = pointEstimate, MLEfit = MLEfit
    )
  } else {
    stop("Invalid estimation method: please use either \"MLE\" or \"bayesian\".")
  }

  f$call <- mc
  print(f)
  f
}

# -------------------------------------------------------------------------------------------

#' TFP trend model
#'
#' @description Creates a state space object object of class \code{TFPmodel} which can be 
#'   fitted using \code{fit}.
#'
#' @param tsl A list of time series objects, see details.
#' @param trend A character string specifying the trend model. \code{trend = "RW1"} denotes
#'   a first order random walk, \code{trend = "RW2"} a second order random walk (local linear
#'   trend) and \code{trend = "DT"} a damped trend model. The default is \code{trend = "DT"}.
#' @param cycle A character string specifying the cycle model. \code{cycle = "AR1"} denotes
#'   an AR(1) process, \code{cycle = "AR2"} an AR(2) process, \code{cycle = "RAR2"} a
#'   reparametrized AR(2) process. The default is \code{cycle = "AR2"}.
#' @param cycleLag A non-negative integer specifying the maximum cycle lag that is included
#'   in the CUBD equation. The default is \code{cycleLag = 0}, see details.
#' @param cubsAR A non-negative integer specifying the maximum CUBS lag that is included
#'   in the CUBS equation. The default is \code{cubsAR = 0}, see details.
#' @param cubsErrorARMA A vector with non-negative integers specifying the AR
#'   and MA degree of the error term in the CUBS equation. The default is
#'   \code{cubsErrorARMA = c(0, 0)}, see details.
#' @param start (Optional) Start vector for the estimation, e.g. \code{c(1980, 1)}.
#' @param end (Optional) End vector for the estimation, e.g. \code{c(2020, 1)}.
#' @param anchor (Optional) Snchor value for the log of the TFP trend.
#' @param anchor.h (Optional) Anchor horizon in the frequency of the given time series.
#'
#' @details The list of time series \code{tsl} needs to have the following components:
#' \describe{
#'     \item{tfp}{Total factor productivity.}
#'     \item{cubs}{Capacity utilization economic sentiment indicator.}
#'     }
#' @details A \code{cycleLag} equal to \code{0} implies that only the contemporaneous cycle
#'   is included in the CUBS equation.  A \code{cycleLag} equal to \code{0:1} implies that
#'   the contemporaneous as well as the lagged cycle are included.
#' @details A \code{cubsAR} equal to \code{0} implies that no autoregressive term is
#'   included in the CUBS equation.  \code{cubsAR = 1} implies that a lagged term is
#'   included, \code{cubsAR = 2} implies that a two lags are included, and so on.
#' @details A \code{cubsErrorARMA} equal to \code{c(0, 0)} implies that the error term in the
#'   CUBS equation is white noise. \code{cubsErrorARMA = c(1, 0)} implies that the error is
#'   an AR(1) process and for \code{cubsErrorARMA = c(1, 2)} the error follows an ARMA(1, 2)
#'   process.
#'
#' @return Object of class TFPmodel, which is a list with the following components:
#'   \item{tsl}{A list of used time series.}
#'   \item{SSModel}{An object of class SSModel specifying the state-space model.}
#'   \item{loc}{A data frame containing information on each involved parameter, for instance
#'         its corresponding system matrix, variable names, and parameter restrictions.}
#'   \item{call}{Original call to the function. }
#'   In addition, the object contains the following attributes:
#'   \item{cycle}{Cycle specification.}
#'   \item{trend}{Trend specification.}
#'   \item{cubs}{A list containing the components \code{cycleLag, cubsAR, errorARMA, exoVariables}.}
#'   \item{anchor}{A list containing the components \code{value, horizon}.}
#'   \item{period}{A list containing the components \code{start, end, frequency}.}
#'
#' @export
#' @importFrom KFAS SSModel SSMregression SSMcustom
#' @importFrom stats start end window ts lag frequency time
#' @importFrom zoo na.trim
#' @examples
#' # load data for Germany
#' data("gap")
#' data("indicator")
#' country <- "Germany"
#' tsList <- amecoData2input(gap[[country]], alpha = 0.65)
#'
#' # compute cubs indicator
#' namesCubs <- c("indu", "serv", "buil")
#' namesVACubs <- paste0("va", namesCubs)
#' tscubs <- cubs(
#'   tsCU = gap[[country]][, namesCubs],
#'   tsVA = gap[[country]][, namesVACubs]
#' )
#' tsList <- c(tsList, tscubs)
#'
#' # define tfp model
#' model <- TFPmodel(tsl = tsList, cycle = "RAR2", cubsErrorARMA = c(1,0))
TFPmodel <- function(tsl, trend = "DT", cycle = "AR2", cycleLag = 0, cubsAR = 0, cubsErrorARMA = c(0, 0),
                     start = NULL, end = NULL, anchor = NULL, anchor.h = NULL) {

  # local variables
  nPar <- nObs <- nState <- nStateV <- stateNames <- loc <- sys <- NULL

  # save call
  mc <- match.call(expand.dots = FALSE)

  # ----- check input for consistency
  errorARMA <- cubsErrorARMA
  if (length(errorARMA) == 1) {
    errorARMA <- c(errorARMA, 0)
  }
  list2env(.checkTfp(
    tsl = tsl,
    trend = trend,
    cycle = cycle,
    cycleLag = cycleLag,
    cubsAR = cubsAR,
    errorARMA = errorARMA,
    start = start,
    end = end,
    anchor = anchor,
    anchor.h = anchor.h
  ),
  envir = environment())

  # ----- preprocess data

  # assign end and start
  if (is.null(end)) {
    end <- end(tsl[[1]])
  }
  if (is.null(start)) {
    start <- start(tsl[[1]])
  }
  # collect used variables
  varUsed <- c("logtfp", "cubs")
  nExo <- 0
  exoNamesTmp <- NULL
  # demean cubs
  tsl$cubs <- 100 * (log(tsl$cubs) - mean(log(tsl$cubs), na.rm = TRUE))
  # log of tfp
  tsl$logtfp <- 100 * log(tsl$tfp)
  # cubs equation data
  if (cubsAR > 0) {
    for (ii in 1:cubsAR) {
      tsl[[paste0("cubsAR", ii)]] <- stats::lag(tsl$cubs, k = -ii)
      varUsed <- c(varUsed, paste0("cubsAR", ii))
      nExo <- nExo + 1
      exoNamesTmp <- c(exoNamesTmp, paste0("cubsAR", ii))
    }
  }
  # list of used time series
  tslUsed <- tsl[varUsed]
  tslUsed <- lapply(tslUsed, function(x) suppressWarnings(window(x, start = start, end = end, extend = TRUE))) # warns if start or end are not changed.
  tslUsed <- lapply(tslUsed, zoo::na.trim)

  # update start and end date if necessary
  start <- dateTsList(tslUsed, FUN1 = stats::start, FUN2 = base::max)
  end <- dateTsList(tslUsed, FUN1 = stats::end, FUN2 = base::max)
  tslUsed <- lapply(tslUsed, function(x) suppressWarnings(window(x, start = start, end = end, extend = TRUE))) # warns if start or end are not changed.
  tslUsed <- lapply(tslUsed, zoo::na.trim)

  # ----- system matrices pre processing
  list2env(.SSSystem(
    tsl = tslUsed,
    cycle = cycle,
    trend = trend,
    cycleLag = cycleLag,
    errorARMA = errorARMA
  ),
  envir = environment())
  tslUsed <- lapply(tslUsed, function(x) suppressWarnings(window(x, start = start, end = end, extend = TRUE))) # warns if start or end are not changed.

  # ----- state space model
  if (is.null(exoNamesTmp)) {
    modelSS <- SSModel(cbind(logtfp, cubs) ~ -1
    + SSMcustom(
        Z = sys$Zt, T = sys$Tt, R = sys$Rt, Q = sys$Qt,
        a1 = sys$a1, P1 = sys$P1, P1inf = sys$P1inf,
        state_names = stateNames
      ),
    H = sys$Ht,
    data = tslUsed
    )
  } else {
    tmp <- paste0("~ ", paste(exoNamesTmp, collapse = " + "))
    modelSS <- SSModel(cbind(logtfp, cubs) ~ -1
    + SSMregression(as.formula(tmp), index = 2, data = tslUsed, state_names = exoNamesTmp)
      + SSMcustom(
        Z = sys$Zt, T = sys$Tt, R = sys$Rt, Q = sys$Qt,
        a1 = sys$a1, P1 = sys$P1, P1inf = sys$P1inf,
        state_names = stateNames
      ),
    H = sys$Ht,
    data = tslUsed
    )
  }

  # -----  tfp model
  model <- list(
    tsl = tslUsed,
    SSModel = modelSS,
    call = mc
  )
  lcubs <- list(
    cycleLag = cycleLag,
    cubsAR = cubsAR,
    errorARMA = errorARMA,
    exoVariables = exoNamesTmp
  )
  lAnchor <- list(
    value = anchor,
    horizon = anchor.h
  )
  lPeriod <- list(
    start = paste0(start(tslUsed$logtfp)[1], ifelse(frequency(tslUsed$logtfp) == 4, paste0(" Q", start(tslUsed$logtfp)[2]), "")),
    end = paste0(end(tslUsed$logtfp)[1], ifelse(frequency(tslUsed$logtfp) == 4, paste0(" Q", end(tslUsed$logtfp)[2]), "")),
    frequency = ifelse(frequency(tslUsed$logtfp) == 4, "quarterly", "annual")
  )
  class(model) <- c("TFPmodel", "model")
  attr(model, "cycle") <- cycle
  attr(model, "trend") <- trend
  attr(model, "cubs") <- lcubs
  attr(model, "anchor") <- lAnchor
  attr(model, "period") <- lPeriod
  loc <- .initializeLoc(model)
  model$loc <- loc

  # return model
  invisible(model)
}

# -------------------------------------------------------------------------------------------

#' \code{TFPmodel} object check
#'
#' @description Tests whether the input object is a valid object of class \code{TFPmodel}.
#'
#' @param object An object to be tested.
#' @param return.logical If \code{return.logical = FALSE} (default), an error message is printed
#'   if the object is not of class \code{TFPmodel}. If \code{return.logical = TRUE}, a logical
#'   value is returned.
#'
#' @return A logical value or nothing, depending on the value of \code{return.logical}.
#'
#' @export
#' @importFrom KFAS is.SSModel
#' @examples
#' # load data for Germany
#' data("gap")
#' data("indicator")
#' country <- "Germany"
#' tsList <- amecoData2input(gap[[country]], alpha = 0.65)
#'
#' # compute cubs indicator
#' namesCubs <- c("indu", "serv", "buil")
#' namesVACubs <- paste0("va", namesCubs)
#' tscubs <- cubs(
#'   tsCU = gap[[country]][, namesCubs],
#'   tsVA = gap[[country]][, namesVACubs]
#' )
#' tsList <- c(tsList, tscubs)
#'
#' # define tfp model
#' model <- TFPmodel(
#'   tsl = tsList, trend = "DT", cycle = "RAR2",
#'   cycleLag = 2, cubsErrorARMA = c(1, 0)
#' )
#' is.TFPmodel(model, return.logical = TRUE)
#' attr(model, "cubs")$cycleLag <- 1
#' is.TFPmodel(model, return.logical = TRUE)
is.TFPmodel <- function(object, return.logical = FALSE) {
  cycle <- attr(object, "cycle")
  trend <- attr(object, "trend")
  cycleLag <- attr(object, "cubs")$cycleLag
  cubsAR <- attr(object, "cubs")$cubsAR
  errorARMA <- attr(object, "cubs")$errorARMA
  anchor <- attr(object, "anchor")$value
  anchor.h <- attr(object, "anchor")$horizon

  components <- c("tsl", "SSModel", "loc", "call")
  trendPossibilities <- c("RW1", "RW2", "DT")
  cyclePossibilities <- c("AR1", "AR2", "RAR2")
  cycleLagPossibilities <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  cubsARPossibilities <- c(0, 1, 2)
  errorARPossibilities <- c(0, 1, 2)
  errorMAPossibilities <- 0
  obsNames <- c("logtfp", "cubs")

  if (return.logical) {
    x <- inherits(object, "TFPmodel") &&
      all(components %in% names(object)) &&
      inherits(object$SSModel, "SSModel") &&
      is.SSModel(object$SSModel, return.logical = TRUE) &&
      all(trend %in% trendPossibilities) &&
      all(cycle %in% cyclePossibilities) &&
      all(cycleLag %in% cycleLagPossibilities) &&
      all(cubsAR %in% cubsARPossibilities) &&
      all(errorARMA[1] %in% errorARPossibilities) &&
      all(errorARMA[2] %in% errorMAPossibilities) &&
      !(!is.null(anchor) && is.null(anchor.h)) &&
      all(obsNames %in% names(object$tsl))
    x
  } else {
    if (!inherits(object, "TFPmodel")) {
      stop("Object is not of class 'TFPmodel'")
    }
    if (!all(components %in% names(object))) {
      stop(paste0(
        "Model is not a proper object of class 'TFPmodel'.",
        "The following components are missing: ", paste0(components[!(components %in% names(object))], collapse = ", ")
      ))
    }
    if (!inherits(object$SSModel, "SSModel")) {
      stop("Object component 'SSModel' is not of class 'SSModel'")
    }
    if (!is.SSModel(object$SSModel, return.logical = TRUE)) {
      stop("Object component 'SSModel' is not a proper object of class 'SSModel'.")
    }
    if (!all(trend %in% trendPossibilities)) {
      stop(paste0(
        "Invalid trend specification. ",
        "Valid trends: ", paste0(trendPossibilities, collapse = ", ")
      ))
    }
    if (!all(cycle %in% cyclePossibilities)) {
      stop(paste0(
        "Invalid cycle specification. ",
        "Valid cycles: ", paste0(cyclePossibilities, collapse = ", ")
      ))
    }
    if (!all(cycleLag %in% cycleLagPossibilities)) {
      stop(paste0(
        "Invalid cycleLag specification. ",
        "Valid cycleLags: ", paste0(cycleLagPossibilities, collapse = ",")
      ))
    }
    if (!all(cubsAR %in% cubsARPossibilities)) {
      stop(paste0(
        "Invalid cubsAR specification. ",
        "Valid cubsARs: ", paste0(cubsARPossibilities, collapse = ",")
      ))
    }
    if (!all(errorARMA[1] %in% errorARPossibilities)) {
      stop(paste0(
        "Invalid error AR specification. ",
        "Valid error AR specification: ", paste0(errorARPossibilities, collapse = ",")
      ))
    }
    if (!all(errorARMA[2] %in% errorMAPossibilities)) {
      stop(paste0(
        "Invalid error MA specification. ",
        "Valid error MA specification: ", paste0(errorMAPossibilities, collapse = ",")
      ))
    }
    if (!is.null(anchor) && is.null(anchor.h)) {
      stop("Anchor specified but no anchor horizon specified. ")
    }
    if (!all(obsNames %in% names(object$tsl))) {
      stop(paste0("Object component 'tsl' does not contain one or both of ", paste0(obsNames, collapse = ", ")))
    }
  }
}

# -------------------------------------------------------------------------------------------

#' \code{TFPfit} object check
#'
#' @description Tests whether the input object is a valid object of class \code{TFPfit}.
#'
#' @param object An object to be tested.
#' @param return.logical If \code{return.logical = FALSE} (default), an error message is printed
#'   if the object is not of class \code{TFPfit}. If \code{return.logical = TRUE}, a logical
#'   value is returned.
#'
#' @return A logical value or nothing, depending on the value of \code{return.logical}.
#'
#' @keywords internal
is.TFPfit <- function(object, return.logical = FALSE) {
  components <- list()
  components$MLE <- c("model", "tsl", "SSMfit", "SSMout", "parameters", "parRestr", "fit", "call")
  components$bayesian <- c("model", "tsl", "parameters", "parametersSampled",
                           "statesSampled", "mcmc", "prior", "fit", "MLE", "call")

  type <- attr(object, "method")

  if (return.logical) {
    x <- inherits(object, "TFPfit") &&
      all(components[[type]] %in% names(object))
    x
  } else {
    if (!inherits(object, "TFPfit")) {
      stop("Object is not of class 'TFPfit'")
    }
    if (!all(components[[type]] %in% names(object))) {
      stop(paste0(
        "Model is not a proper object of class 'TFPfit'.",
        "The following components are missing: ", paste0(components[!(components %in% names(object))], collapse = ", ")
      ))
    }
  }
}


# -------------------------------------------------------------------------------------------

#' Print \code{TFPmodel} object
#'
#' @description Prints the model specifications of an object of class \code{TFPmodel}.
#'
#' @param x An object of class \code{TFPmodel}.
#' @param call A logical. If \code{TRUE}, the call will be printed.
#' @param check A logical. If \code{TRUE}, the model class will be checked.
#' @param ... Ignored.
#' @return No return value, model information is printed.
#' @export
print.TFPmodel <- function(x, call = TRUE, check = TRUE, ...) {
  .printSSModel(x = x, call = call, check = check)
}

# -------------------------------------------------------------------------------------------

#' Print \code{TFPfit} object
#'
#' @description Prints the model specifications and the estimation results of an object of class \code{TFPfit}.
#'
#' @param x An object of class \code{TFPfit}.
#' @param ... Ignored.
#' @return No return value, results are printed.
#' @export
print.TFPfit <- function(x, ...) {
  .printSSModelFit(x = x, call = TRUE, check = FALSE, print.model = TRUE)
}

# -------------------------------------------------------------------------------------------

#' Plots for a \code{TFPfit} object
#'
#' @description Plots the TFP trend and the CUBS equation and gives diagnostic plots based on
#' standardized residuals for objects of class \code{TFPfit}.
#'
#' @param x An object of class \code{TFPfit}.
#' @param alpha The significance level for the TFP trend (\code{alpha in [0,1]}). Only used if
#'   \code{bounds = TRUE}.
#' @param bounds A logical indicating whether significance intervals should be plotted around
#'   tfp growth. The default is \code{bounds = TRUE}.
#' @param combine A logical indicating whether the diagnostic plots should be combined or not,
#'   the default is \code{TRUE}.
#' @param posterior A logical indicating whether posterior diagnostics should be plotted. The
#'   default is \code{FALSE}. Only applied in the case of bayesian estimation.
#' @inheritParams plot.gap
#'
#' @return No return value, plots are printed.
#' 
#' @export
plot.TFPfit <- function(x, alpha = 0.05, bounds = TRUE, path = NULL, combine = TRUE, prefix = NULL,
                        posterior = FALSE, device = "png", width = 10, height = 3, ...) {
  if (!is.TFPfit(x, return.logical = TRUE)) {
    stop("x is no valid object of class 'TFPfit'.")
  }
  if (!is.null(path)) {
    check <- dir.exists(paths = path)
    if (!check) {
      warning(paste0("The given file path '", path, "' does not exist. The plots will not be saved."))
      path <- NULL
    }
  }
  method <- attr(x, "method")

  # check whether prediction is present
  prediction <- "prediction" %in% names(attributes(x))
  
  if (posterior == FALSE) {

    if (!prediction) {
      
      # ----- estimation
      
      if (method == "MLE") {
        # confidence bounds
        tvalue <- -qnorm((alpha) / 2)
        boundName <- paste0(100 * (1 - alpha), "% CI")
        tslBounds <- list(
          ub = (x$tsl$tfpTrendGrowth + x$tsl$tfpTrendGrowthSE * tvalue),
          lb = (x$tsl$tfpTrendGrowth - x$tsl$tfpTrendGrowthSE * tvalue),
          ub2 = (x$tsl$obsFitted[, 2] + x$tsl$obsFittedSE[, 2] * tvalue),
          lb2 = (x$tsl$obsFitted[, 2] - x$tsl$obsFittedSE[, 2] * tvalue),
          ub3 = (x$tsl$tfpTrend + x$tsl$tfpTrendSE * tvalue),
          lb3 = (x$tsl$tfpTrend - x$tsl$tfpTrendSE * tvalue)
        )

        # residuals
        res <- x$tsl$obsResidualsRecursive[, "cubs"]
        
      } else { # Bayesian
        # confidence bounds
        HPDI <- attr(x, "HPDIprob")
        boundName <- paste0(100 * HPDI, "% HPDCI")
        tslBounds <- list(
          ub = x$tsl$tfpTrendGrowthSummary[, paste0(100 * HPDI, "% HPDI-UB")],
          lb = x$tsl$tfpTrendGrowthSummary[, paste0(100 * HPDI, "% HPDI-LB")],
          ub2 = x$tsl$cubsFittedSummary[, paste0(100 * HPDI, "% HPDI-UB")],
          lb2 = x$tsl$cubsFittedSummary[, paste0(100 * HPDI, "% HPDI-LB")],
          ub3 = x$tsl$tfpTrendSummary[, paste0(100 * HPDI, "% HPDI-UB")],
          lb3 = x$tsl$tfpTrendSummary[, paste0(100 * HPDI, "% HPDI-LB")]
        )

        # residuals
        res <- NULL
      }

      # --- data
      # tfp level
      tsl1 <- list(
        trend = x$tsl$tfpTrend,
        orig = exp(x$model$tsl$logtfp / 100),
        lb = tslBounds$lb3,
        ub = tslBounds$ub3
      )
      if (!is.null(x$tsl$tfpTrendAnchored)) {
        tsl1 <- c(tsl1, list(anchor = x$tsl$tfpTrendAnchored))
      }
      tsl1 <- do.call(cbind, tsl1)
  
      # cubs equation
      tsl2 <- do.call(cbind, list(
        fitted = x$tsl$obsFitted[, 2],
        E2 = x$model$tsl$cubs,
        lb = tslBounds$lb2,
        ub = tslBounds$ub2
      ))
      
      # tfp growth
      tsl3 <- list(
        trend = 100 * x$tsl$tfpTrendGrowth,
        orig = 100 * growth(exp(x$model$tsl$logtfp / 100)),
        lb = 100 * tslBounds$lb,
        ub = 100 * tslBounds$ub
      )
      if (!is.null(x$tsl$tfpTrendAnchored)) {
        tsl3 <- c(tsl3, list(anchor = 100 * growth(x$tsl$tfpTrendAnchored)))
      }
      tsl3 <- do.call(cbind, tsl3)
  
      # combine lists
      tsl <- list(tsl1, tsl2, tsl3)
  
      # --- legends and titles and print names
      legend <- list(
        c("trend tfp", "tfp", "anchored trend tfp"),
        c("fitted", "cubs"),
        c("trend tfp growth", "tfp growth", "anchored trend tfp growth")
      )
      title <- list(
        "Total factor productivity",
        "CUBS",
        "CUBS residuals",
        "Total factor productivity growth in %"
      )
      namesPrint <- c("tfp", "cubs", "tfp_growth")
      if (!is.null(prefix)) namesPrint <- paste(prefix, namesPrint , sep = "_")

      # plot
      plotSSresults(
        tsl = tsl, legend = legend, title = title,
        boundName = boundName, res = res, namesPrint = namesPrint,
        bounds = bounds, combine = combine, path = path, device = device,
        width = width, height = height
      )
      
    } else {
      
      # ----- prediction
      
      n.ahead <- attr(x, "prediction")$n.ahead
      
      if (method == "MLE") {
        # confidence bounds
        tvalue <- -qnorm((alpha) / 2)
        boundName <- paste0(100 * (1 - alpha), "% CI")
        tslBounds <- list(
          lb = (x$tsl$tfpTrend - x$tsl$tfpTrendSE * tvalue),
          ub = (x$tsl$tfpTrend + x$tsl$tfpTrendSE * tvalue),
          lb1 = (x$tsl$tfp - x$tsl$tfpSE * tvalue),
          ub1 = (x$tsl$tfp + x$tsl$tfpSE * tvalue),
          lb2 = (x$tsl$obs[, 2] - x$tsl$obsSE[, 2] * tvalue),
          ub2 = (x$tsl$obs[, 2] + x$tsl$obsSE[, 2] * tvalue),
          lb3 = ((x$tsl$stateSmoothed[, "cycle"] - x$tsl$stateSmoothedSE[, "cycle"] * tvalue)),
          ub3 = ((x$tsl$stateSmoothed[, "cycle"] + x$tsl$stateSmoothedSE[, "cycle"] * tvalue)),
          lb4 = 100 * (x$tsl$tfpTrendGrowth - x$tsl$tfpTrendGrowthSE * tvalue),
          ub4 = 100 * (x$tsl$tfpTrendGrowth + x$tsl$tfpTrendGrowthSE * tvalue),
          lb41 = 100 * (x$tsl$tfpGrowth - x$tsl$tfpGrowthSE * tvalue),
          ub41 = 100 * (x$tsl$tfpGrowth + x$tsl$tfpGrowthSE * tvalue)
        )
        
      } else { # Bayesian
        
        # confidence bounds
        HPDI <- attr(x, "HPDIprob")
        boundName <- paste0(100 * HPDI, "% HPDCI")
        tslBounds <- list(
          ub = x$tsl$tfpTrendSummary[, paste0(100 * HPDI, "% HPDI-UB")],
          lb = x$tsl$tfpTrendSummary[, paste0(100 * HPDI, "% HPDI-LB")],
          ub1 = x$tsl$tfpSummary[, paste0(100 * HPDI, "% HPDI-UB")],
          lb1 = x$tsl$tfpSummary[, paste0(100 * HPDI, "% HPDI-LB")],
          ub2 = x$tsl$obsSummary[, paste0("cubs.", 100 * HPDI, "% HPDI-UB")],
          lb2 = x$tsl$obsSummary[, paste0("cubs.", 100 * HPDI, "% HPDI-LB")],
          ub3 = x$tsl$stateSmoothedSummary[, paste0("cycle.", 100 * HPDI, "% HPDI-UB")],
          lb3 = x$tsl$stateSmoothedSummary[, paste0("cycle.", 100 * HPDI, "% HPDI-LB")],
          ub4 = 100 * x$tsl$tfpTrendGrowthSummary[, paste0(100 * HPDI, "% HPDI-UB")],
          lb4 = 100 * x$tsl$tfpTrendGrowthSummary[, paste0(100 * HPDI, "% HPDI-LB")],
          ub41 = 100 * x$tsl$tfpGrowthSummary[, paste0(100 * HPDI, "% HPDI-UB")],
          lb41 = 100 * x$tsl$tfpGrowthSummary[, paste0(100 * HPDI, "% HPDI-LB")]
          )
        
      }
      
      # tfp
      tsl1 <- do.call(cbind, list(
        trend = x$tsl$tfpTrend,
        orig = x$tsl$tfp,
        lb = tslBounds$lb,
        ub = tslBounds$ub,
        lb_fc = tslBounds$lb1,
        ub_fc = tslBounds$ub1
      ))

      # cubs
      tsl2 <- do.call(cbind, list(
        E2 = x$tsl$obs[, 2],
        lb_fc = tslBounds$lb2,
        ub_fc = tslBounds$ub2
      ))
      
      # cycle
      tsl3 <- do.call(cbind, list(
        cycle = x$tsl$stateSmoothed[, "cycle"],
        lb = tslBounds$lb3,
        ub = tslBounds$ub3
      ))
      
      # tfp growth
      tsl4 <- do.call(cbind, list(
        trend = 100 * x$tsl$tfpTrendGrowth,
        orig = 100 * x$tsl$tfpGrowth,
        lb = tslBounds$lb4,
        ub = tslBounds$ub4,
        lb_fc = tslBounds$lb41,
        ub_fc = tslBounds$ub41
      ))
      
      # combine lists
      tsl <- list(tsl1, tsl2, tsl3, tsl4)
      
      # --- legends and titles and print names
      legend <- list(
        c("trend tfp", "tfp", rep(paste0(boundName, " (tfp)"), 2)),
        c("cubs", rep(paste0(boundName, ""), 2)),
        c("tfp gap"),
        c("trend tfp growth", "tfp growth", rep(paste0(boundName, " (tfp growth)"), 2))
      )
      title <- list(
        "Total factor productivity",
        "CUBS",
        "Total factor productivity gap",
        "Total factor productivity growth in %"
      )
      namesPrint <- c("tfp", "cubs", "tfp_gap", "tfp_growth")
      if (!is.null(prefix)) namesPrint <- paste(prefix, namesPrint , sep = "_")
      
      # plot
      plotSSprediction(
        tsl = tsl, legend = legend, title = title, n.ahead = n.ahead,
        boundName = boundName, res = NULL, namesPrint = namesPrint,
        bounds = bounds, combine = combine, path = path, device = device,
        width = width, height = height
      )
      
    }

  # ----- bayesian estimation
  } else if (posterior == TRUE & method != "MLE") {

      R <- attr(x, "R")
      burnin <- attr(x, "burnin")
      thin <- attr(x, "thin")
      HPDIprob <- attr(x, "HPDIprob")
      FUN <- attr(x, "FUN")
      loc <- x$model$loc
      param <- x$parametersSampled
      parNames <- colnames(param)
      cycle <- attr(x$model, "cycle")
      prior <- x$prior

      # convert prior mean and standard deviation to actual parameters
      # gamma and beta needs to be adjusted
      distrPar <- .priorMSd2Parameter(
        prior = cbind(prior$cycle, prior$trend, prior[[3]])[1:2, ],
        restr = cbind(prior$cycle, prior$trend, prior[[3]])[3:4, ],
        namesInvGammaDistr = loc$varName[loc$distribution == "invgamma"],
        namesNormalDistr = loc$varName[loc$distribution == "normal"],
        namesBetaDistr = loc$varName[loc$distribution == "beta"]
      )

      for (jj in parNames) {
        distr <- loc$distribution[loc$varName == jj]

        inputl <- list(draws = param[, jj], burnin = burnin / thin, parName = gsub("E2", "cu", jj), lb = distrPar[3, jj], ub = distrPar[4, jj])
        switch(distr,
          "invgamma" = do.call(plot_gibbs_output, args = c(path, inputl,
            prefix = prefix, shape = distrPar[2, jj] / 2,
            scale = distrPar[1, jj] / 2, device = device,
            width = width, height = height
          )),
          "beta" = do.call(plot_gibbs_output, args = c(path, inputl,
            prefix = prefix, shape1 = distrPar[1, jj],
            shape2 = distrPar[2, jj], device = device,
            width = width, height = height
          )),
          "normal" = do.call(plot_gibbs_output, args = c(path, inputl,
            prefix = prefix, mu = distrPar[1, jj],
            prec = 1 / distrPar[2, jj], device = device,
            width = width, height = height
          ))
        )
      }
  }
}
