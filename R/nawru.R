
# -------------------------------------------------------------------------------------------

#' Estimation of a \code{NAWRUmodel}
#'
#' @description Estimates a two-dimensional state-space model and performs filtering and
#'   smoothing to obtain the NAWRU using either maximum likelihood estimation or bayesian
#'   methods.
#'
#' @param model An object of class NAWRUmodel.
#' @param parRestr A list of matrices containing the parameter restrictions for the cycle,
#'   trend, and the Phillip's curve. Each matrix contains the lower and upper bound of the
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
#' @param MLEfit (Optional) An object of class \code{NAWRUfit} which is used for
#'   initialization. Only used if \code{method = "bayesian"}.
#' @param control (Optional) A list of control arguments to be passed on to \code{optim}.
#' @param ... additional arguments to be passed to the methods functions.
#'
#' @details The list object \code{prior} contains three list elements \code{cycle},
#'   \code{trend}, and \code{pcInd}. Each list element is a \code{4 x n} matrix where \code{n}
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
#'            \item Phillip's curve equation parameters \eqn{\theta_{pcInd}}: Conditional on the
#'            states \eqn{\alpha_r}, a draw \eqn{\theta_{pcInd,k}} is obtained either
#'            by a sequential Gibbs step, a Metropolis Hasting step, a combination
#'            thereof, or by conjugacy, depending on the Phillip's curve equation specification.}
#'
#' @return For maximum likelihood estimation, an object of class \code{NAWRUfit} containing
#'   the following components:
#'   \item{model}{The input object of class \code{NAWRUmodel}.}
#'   \item{SSMfit}{The estimation output from the function \code{fitSSM} from \code{KFAS}.}
#'   \item{SSMout}{The filtering and smoothing output from the function \code{KFS} from
#'         \code{KFAS}.}
#'   \item{parameters}{A data frame containing the estimated parameters, including standard
#'         errors, t-statistics, and p-values.}
#'   \item{parRestr}{A list of matrices containing the enforced parameter constraints.}
#'   \item{fit}{A list of model fit criteria (see below).}
#'   \item{call}{Original call to the function. }
#'   The list component \code{fit} contains the following model fit criteria:
#'   \item{loglik}{Log-likelihood function values.}
#'   \item{AIC}{Akaike information criterion.}
#'   \item{BIC}{Bayesian information criterion.}
#'   \item{AIC}{Hannan-Quinn information criterion.}
#'   \item{RMSE}{root mean squared error of the Phillip's curve equation.}
#'   \item{R2}{R squared of the Phillip's curve equation.}
#'   \item{LjungBox}{Ljung-Box test output of the Phillip's curve equation.}
#'   \item{signal-to-noise}{Signal-to-noise ratio.}
#' For bayesian estimation, an object of class \code{NAWRUfit} containing the following
#'   components:
#'   \item{model}{The input object of class \code{NAWRUmodel}.}
#'   \item{tsl}{A list of time series containing the estimated states.}
#'   \item{parameters}{A data frame containing the estimated parameters, including standard
#'         errors, highest posterior density credible sets.}
#'   \item{prior}{A list of matrices containing the used prior distributions.}
#'   \item{fit}{A list of model fit criteria (see below).}
#'   \item{call}{Original call to the function. }
#'   The list component \code{fit} contains the following model fit criteria:
#'   \item{R2}{R squared of the phillips curve equation,}
#'   \item{signal-to-noise}{Signal-to-noise ratio.}
#'
#' @export
#' @family fitting methods
#' @examples
#' # define nawru model for France
#' data("gap")
#' country <- "France"
#' tsList <- amecoData2input(gap[[country]])
#' model <- NAWRUmodel(tsl = tsList)
#' # estimate nawru model via MLE
#' parRestr <- initializeRestr(model = model, type = "hp")
#' \donttest{
#' f <- fit(model = model, parRestr = parRestr)
#' }
#' # initialize priors and estimate model via Bayesian methods
#' prior <- initializePrior(model = model)
#' \donttest{
#' f <- fit(model = model, method = "bayesian", prior = prior, R = 5000, thin = 2)
#' }
fit.NAWRUmodel <- function(model, parRestr = initializeRestr(model = model), signalToNoise = NULL,
                     method = "MLE", control = NULL,
                     prior = initializePrior(model), R = 10000, burnin = ceiling(R / 10),
                     thin = 1, HPDIprob = 0.85, pointEstimate = "mean", MLEfit = NULL, ...) {
  # save call
  mc <- match.call(expand.dots = FALSE)

  if (method == "MLE") {
    f <- .MLEfitNAWRU(
      model = model, parRestr = parRestr, signalToNoise = signalToNoise,
      control = control
    )
  } else if (method == "bayesian") {
    if (!is.null(signalToNoise)) {
      warning("'signalToNoise' not applicable for method = 'bayesian'.")
    }
    f <- .BayesFitNAWRU(
      model = model, prior = prior, R = R, burnin = burnin,
      thin = thin, HPDIprob = HPDIprob, FUN = pointEstimate, MLEfit = MLEfit
    )
  } else {
    stop("Invalid estimation method: please use either 'MLE' or 'bayesian'.")
  }

  f$call <- mc
  print(f)
  f
}

# -------------------------------------------------------------------------------------------

#' NAWRU model
#'
#' @description Creates a state space object object of class \code{NAWRUmodel} which can be 
#'   fitted using \code{fit}.
#'
#' @param tsl A list of time series objects, see details.
#' @param trend A character string specifying the trend model. \code{trend = "RW1"} denotes
#'   a first order random walk, \code{trend = "RW2"} a second order random walk (local linear
#'   trend) and \code{trend = "DT"} a damped trend model. The default is \code{trend = "RW2"}.
#' @param cycle A character string specifying the cycle model. \code{cycle = "AR1"} denotes
#'   an AR(1) process, \code{cycle = "AR2"} an \code{AR(2)} process. The default is
#'   \code{cycle = "AR2"}.
#' @param type A character string specifying the type of the Phillip's curve.
#'   \code{type = "TKP"} denotes the traditional Keynesian Phillip's curve and
#'   \code{type = "NKP"} the New Keynesian Phillip's curve, see details. The default is
#'   \code{type = "TKP"}.
#' @param cycleLag A vector specifying the cycle lags that are included in the Phillip's
#'   curve. The default is \code{cycleLag = 0}, see details.
#' @param pcErrorARMA A \code{2 x 1} vector with non-negative integers specifying the AR
#'   and MA degree of the error term in the Phillip's curve equation. The default is
#'   \code{pcErrorARMA = c(0, 0)}, see details.
#' @param exoType An optional \code{n x m x 2} array specifying the possible difference
#'  and lag transformation for the variables. \code{exoType} can be initialized using the
#'  function \code{inizializeExo}. The column names give the variable names.
#'  \code{exoType[, , 1]} contains the difference transformations and \code{exoType[, , 2]}
#'  the subsequent lag transformations, see details.
#' @param start (Optional) Start vector for the estimation, e.g. \code{c(1980, 1)}.
#' @param end (Optional) End vector for the estimation, e.g. \code{c(2020, 1)}.
#' @param anchor (Optional) Anchor value for the unemployment rate.
#' @param anchor.h (Optional) Anchor horizon in the frequency of the given time series.
#'
#' @details The list of time series \code{tsl} needs to have the following components:
#' \describe{
#'     \item{ur}{Unemployment rate.}
#'     \item{nulc}{Nominal Unit labor costs, if \code{type = "TKP"}.}
#'     \item{rulc}{Real unit labor costs, if \code{type = "NKP"}.}
#'     }
#'   and optionally other variables included in \code{exoType}.
#' @details A \code{cycleLag} equal to \code{0} implies that only the contemporaneous cycle
#'   is included in the Phillip's curve.  A \code{cycleLag} equal to \code{0:1} implies that
#'   the contemporaneous as well as the lagged cycle are included.
#' @details A \code{pcErrorARMA} equal to \code{c(0, 0)} implies that the error term in the
#'   Phillip's curve is white noise. \code{pcErrorARMA = c(1, 0)} implies that the error is
#'   an AR(1) process and for \code{pcErrorARMA = c(1, 2)} the error follows an ARMA(1, 2)
#'   process.
#' @details For the New Keynesian Phillip's curve, the \code{cycleLag} cannot be chosen.
#'   \code{cycleLag} will be set to \code{0} if \code{cycle = "AR1"} and to \code{1} if
#'   \code{cycle = "AR2"}. In the latter case, the forward solution of the Phillip's curve
#'   implies parameter restrictions for the lagged cycle on the Phillip's curve. Moreover,
#'   exogenous variables will be ignored in the case of the New Keynesian Phillip's curve.
#' @details The array \code{exoType} consists of non-negative integers or \code{NA}s.
#'   \code{exoType[, , 1] = c(NA,1)} and \code{exoType[, , 2] = c(NA,2)} implies that
#'   the first variable is not included in the Phillip's curve whereas the second lag of
#'   the first difference of the second variable is included.
#'
#' @return Object of class \code{NAWRUmodel}, which is a list with the following components:
#'   \item{tsl}{A list of used time series.}
#'   \item{SSModel}{An object of class SSModel specifying the state-space model.}
#'   \item{loc}{A data frame containing information on each involved parameter, for instance
#'         its corresponding system matrix, variable names, and parameter restrictions.}
#'   \item{call}{Original call to the function.}
#'   In addition, the object contains the following attributes:
#'   \item{cycle}{Cycle specification.}
#'   \item{trend}{Trend specification.}
#'   \item{phillipsCurve}{A list containing the components \code{type, cycleLag, errorARMA, exoVariables}.}
#'   \item{anchor}{A list containing the components \code{value, horizon}.}
#'   \item{period}{A list containing the components \code{start, end, frequency}.}
#'
#' @export
#' @importFrom KFAS SSModel SSMregression SSMcustom
#' @importFrom stats start end window ts lag frequency time
#' @importFrom zoo na.trim
#' @examples
#' # load data for France
#' data("gap")
#' tsList <- amecoData2input(gap$France, alpha = 0.65)
#' # Traditional phillips curve
#' model <- NAWRUmodel(tsl = tsList, trend = "RW2", cycle = "AR2", type = "TKP", cycleLag = 0)
#'
#' # New-Keynesian Phillips curve
#' model <- NAWRUmodel(tsl = tsList, trend = "RW2", cycle = "AR2", type = "NKP", cycleLag = 0:1)
#'
#' # Traditional Phillips curve with 6 exogenous variables
#' # specify exogenous variable transformations
#' D <- matrix(c(2, 2, 2, 1, 1, 1), 2, 3, byrow = TRUE)
#' L <- matrix(c(0, 0, 0, 1, 1, 1), 2, 3, byrow = TRUE)
#' exoType <- initializeExo(varNames = c("tot", "prod","ws"), D = D, L = L)
#' model <- NAWRUmodel(tsl = tsList, cycleLag = 0:1, exoType = exoType)
NAWRUmodel <- function(tsl, trend = "RW2", cycle = "AR2", type = "TKP", cycleLag = 0, pcErrorARMA = c(0, 0),
                       exoType = NULL, start = NULL, end = NULL, anchor = NULL, anchor.h = NULL) {

  # local variables
  nPar <- nObs <- nState <- nStateV <- stateNames <- loc <- sys <- NULL

  # save call
  mc <- match.call(expand.dots = FALSE)

  # ----- check input for consistency
  exoNames <- colnames(exoType)
  errorARMA <- pcErrorARMA
  if (length(errorARMA) == 1) {
    errorARMA <- c(errorARMA, 0)
  }
  list2env(.checkNawru(
    tsl = tsl,
    trend = trend,
    cycle = cycle,
    type = type,
    cycleLag = cycleLag,
    errorARMA = errorARMA,
    exoNames = exoNames,
    exoType = exoType,
    start = start,
    end = end,
    anchor = anchor,
    anchor.h = anchor.h
  ),
  envir = environment())

  # ----- preprocess data

  # assign end and start
  if (is.null(end)) end <- end(tsl[[1]])
  if (is.null(start)) start <- start(tsl[[1]])
  # collect used variables
  varUsed <- c("ur")
  nExo <- 0
  exoNamesTmp <- NULL
  # merge observation equation data
  if (type == "TKP") {
    tsl$pcInd <- diff(tsl$winfl) * 100
    tsl$pcInd <- diff(diff(log(tsl$nulc))) * 100
    varUsed <- c(varUsed, "pcInd")
    obsNames <- c("ur", "ddlognulc")
  } else {
    tsl$pcInd <- diff(log(tsl$rulc)) * 100
    tsl$pcIndl <- stats::lag(tsl$pcInd, k = -1)
    varUsed <- c(varUsed, "pcInd", "pcIndl")
    obsNames <- c("ur", "dlogrulc")
    nExo <- 1
    exoNamesTmp <- c(exoNamesTmp, "pcIndl")
  }
  # assign and compute exogenous variables
  if (any(!is.na(exoType)) && type == "TKP") {
    tsl[exoNames[!(substr(exoNames, 1, 5) == "dummy")]] <- lapply(tsl[exoNames[!(substr(exoNames, 1, 5) == "dummy")]], log)
    for (jj in (1:dim(exoType)[2])) {
      name <- exoNames[jj]
      for (ii in (1:dim(exoType)[1])) {
        if (!is.na(exoType[ii, jj, 1])) {
          name_tmp <- paste0(
            "pc", paste(rep("d", exoType[ii, jj, 1]), collapse = ""), name,
            paste(rep("l", exoType[ii, jj, 2]), collapse = "")
          )
          tsl[[name_tmp]] <- tsl[[name]]
          try(tsl[[name_tmp]] <- diff(tsl[[name_tmp]], differences = exoType[ii, jj, 1]), silent = TRUE)
          tsl[[name_tmp]] <- stats::lag(tsl[[name_tmp]], k = -abs(exoType[ii, jj, 2]))
          exoNamesTmp <- c(exoNamesTmp, name_tmp)
        }
      }
    }
    nExo <- length(exoNamesTmp)
    varUsed <- c(varUsed, exoNamesTmp)
  }
  # list of used time series
  tslUsed <- tsl[varUsed]
  tslUsed <- lapply(tslUsed, function(x) suppressWarnings(window(x, start = start, end = end))) # warns if start or end are not changed.
  tslUsed <- lapply(tslUsed, zoo::na.trim)

  # update start and end date if necessary
  start <- dateTsList(tslUsed, FUN1 = stats::start, FUN2 = base::max)
  end <- dateTsList(tslUsed, FUN1 = stats::end, FUN2 = base::min)
  tslUsed <- lapply(tslUsed, function(x) suppressWarnings(window(x, start = start, end = end))) # warns if start or end are not changed.
  tslUsed <- lapply(tslUsed, zoo::na.trim)

  # ----- system matrices pre processing
  list2env(.SSSystem(
    tsl = tslUsed,
    cycle = cycle,
    trend = trend,
    cycleLag = cycleLag,
    type = type,
    errorARMA = errorARMA
  ),
  envir = environment())
  tslUsed <- lapply(tslUsed, function(x) suppressWarnings(window(x, start = start, end = end))) # warns if start or end are not changed.

  # ----- state space model

  if (is.null(exoNamesTmp)) {
    modelSS <- SSModel(cbind(ur, pcInd) ~ -1
    + SSMcustom(
        Z = sys$Zt, T = sys$Tt, R = sys$Rt, Q = sys$Qt,
        a1 = sys$a1, P1 = sys$P1, P1inf = sys$P1inf, state_names = stateNames
      ),
    H = sys$Ht,
    data = tslUsed
    )
  } else {
    tmp <- paste0("~ ", paste(exoNamesTmp, collapse = " + "))
    modelSS <- SSModel(cbind(ur, pcInd) ~ -1
    + SSMregression(as.formula(tmp), index = 2, data = tslUsed, state_names = exoNamesTmp)
      + SSMcustom(
        Z = sys$Zt, T = sys$Tt, R = sys$Rt, Q = sys$Qt,
        a1 = sys$a1, P1 = sys$P1, P1inf = sys$P1inf, state_names = stateNames
      ),
    H = sys$Ht,
    data = tslUsed
    )
  }

  # -----  nawru model
  model <- list(
    tsl = tslUsed,
    SSModel = modelSS,
    call = mc
  )
  lPc <- list(
    type = type,
    cycleLag = cycleLag,
    errorARMA = errorARMA,
    exoVariables = exoNamesTmp
  )
  lAnchor <- list(
    value = anchor,
    horizon = anchor.h
  )
  lPeriod <- list(
    start = paste0(start(tslUsed$ur)[1], ifelse(frequency(tslUsed$ur) == 4, paste0(" Q", start(tslUsed$ur)[2]), "")),
    end = paste0(end(tslUsed$ur)[1], ifelse(frequency(tslUsed$ur) == 4, paste0(" Q", end(tslUsed$ur)[2]), "")),
    frequency = ifelse(frequency(tslUsed$ur) == 4, "quarterly", "annual")
  )
  class(model) <- c("NAWRUmodel", "model")
  attr(model, "cycle") <- cycle
  attr(model, "trend") <- trend
  attr(model, "phillips curve") <- lPc
  attr(model, "anchor") <- lAnchor
  attr(model, "period") <- lPeriod
  loc <- .initializeLoc(model)
  if (type != "NKP") {
    loc$variableRow[loc$sysMatrix == "exo"] <- loc$varName[loc$sysMatrix == "exo"]
  }
  model$loc <- loc

  invisible(return(model))
}


# -------------------------------------------------------------------------------------------

#' \code{NAWRUodel} object check
#'
#' @description Tests whether the input object is a valid object of class \code{NAWRUmodel}.
#'
#' @param object An object to be tested.
#' @param return.logical If \code{return.logical = FALSE} (default), an error message is printed
#'   if the object is not of class \code{NAWRUmodel}. If \code{return.logical = TRUE}, a logical
#'   value is returned.
#'
#' @return A logical value or nothing, depending on the value of \code{return.logical}.
#'
#' @export
#' @importFrom KFAS is.SSModel
#' @examples
#' # load data for France
#' data("gap")
#' tsList <- amecoData2input(gap$France, alpha = 0.65)
#'
#' # Traditional phillips curve
#' model <- NAWRUmodel(tsl = tsList, trend = "RW2", cycle = "AR2", type = "NKP", cycleLag = 0:1)
#' is.NAWRUmodel(model, return.logical = TRUE)
#' attr(model, "phillips curve")$cycleLag <- 0
#' is.NAWRUmodel(model, return.logical = TRUE)
is.NAWRUmodel <- function(object, return.logical = FALSE) {
  cycle <- attr(object, "cycle")
  trend <- attr(object, "trend")
  type <- attr(object, "phillips curve")$type
  cycleLag <- attr(object, "phillips curve")$cycleLag
  errorARMA <- attr(object, "phillips curve")$errorARMA
  exoNames <- attr(object, "phillips curve")$exoVariables
  anchor <- attr(object, "anchor")$value
  anchor.h <- attr(object, "anchor")$horizon

  components <- c("tsl", "SSModel", "loc", "call")
  trendPossibilities <- c("RW1", "RW2", "DT")
  cyclePossibilities <- c("AR1", "AR2")
  typePossibilities <- c("TKP", "NKP")
  cycleLagPossibilities <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  errorARPossibilities <- c(0, 1, 2)
  errorMAPossibilities <- 0
  obsNames <- c("ur", "pcInd")
  exoNames <- exoNames[!grepl(obsNames[2], exoNames)]

  if (return.logical) {
    x <- inherits(object, "NAWRUmodel") &&
      all(components %in% names(object)) &&
      inherits(object$SSModel, "SSModel") &&
      is.SSModel(object$SSModel, return.logical = TRUE) &&
      all(trend %in% trendPossibilities) &&
      all(cycle %in% cyclePossibilities) &&
      all(type %in% typePossibilities) &&
      all(cycleLag %in% cycleLagPossibilities) &&
      all(errorARMA[1] %in% errorARPossibilities) &&
      all(errorARMA[2] %in% errorMAPossibilities) &&
      all(obsNames %in% names(object$tsl)) &&
      !(all("NKP" %in% type) && all("AR1" %in% cycle) && max(cycleLag) != 0) &&
      !(all("NKP" %in% type) && all("AR2" %in% cycle) && max(cycleLag) != 1) &&
      !(!is.null(anchor) && is.null(anchor.h)) &&
      !(all("NKP" %in% type) && !is.null(exoNames) && length(exoNames) != 0)
    x
  } else {
    if (!inherits(object, "NAWRUmodel")) {
      stop("Object is not of class 'NAWRUmodel'")
    }
    if (!all(components %in% names(object))) {
      stop(paste0(
        "Model is not a proper object of class 'NAWRUmodel'.",
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
    if (!all(type %in% typePossibilities)) {
      stop(paste0(
        "Invalid Phillips curve type specification. ",
        "Valid types: ", paste0(typePossibilities, collapse = ", ")
      ))
    }
    if (!all(cycleLag %in% cycleLagPossibilities)) {
      stop(paste0(
        "Invalid cycleLag specification. ",
        "Valid cycleLags: ", paste0(cycleLagPossibilities, collapse = ",")
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
    if (!all(obsNames %in% names(object$tsl))) {
      stop(paste0("Object component 'tsl' does not contain one or both of ", paste0(obsNames, collapse = ", ")))
    }
    if (all("NKP" %in% type) && all("AR1" %in% cycle) && max(cycleLag) != 0) {
      stop(paste0(
        "Invalid cycleLag for the New Keynesian Phillip's curve and an AR(1) cycle. ",
        "Maximum cycleLag must be set to 0."
      ))
    }
    if (all("NKP" %in% type) && all("AR2" %in% cycle) && max(cycleLag) != 1) {
      stop(paste0(
        "Invalid cycleLag for the New Keynesian Phillip's curve and an AR(2) cycle. ",
        "Maxmimum cycleLag must be set to 1."
      ))
    }
    if (!is.null(anchor) && is.null(anchor.h)) {
      stop("Nawru anchor specified but no anchor horizon specified. ")
    }
    if (all("NKP" %in% type) && !is.null(exoNames) && length(exoNames) != 0) {
      stop("Exogenous variables are not compatible with the NeW Keynesian Phillip's curve.")
    }
  }
}

# -------------------------------------------------------------------------------------------

#' \code{NAWRUfit} object check
#'
#' @description Tests whether the input object is a valid object of class \code{NAWRUfit}.
#'
#' @param object An object to be tested.
#' @param return.logical If \code{return.logical = FALSE} (default), an error message is printed
#'   if the object is not of class \code{NAWRUfit}. If \code{return.logical = TRUE}, a logical
#'   value is returned.
#'
#' @return A logical value or nothing, depending on the value of \code{return.logical}.
#'
#' @keywords internal
is.NAWRUfit <- function(object, return.logical = FALSE) {
  components <- list()
  components$MLE <- c("model", "tsl", "SSMfit", "SSMout", "parameters", "parRestr", "fit", "call")
  components$bayesian <- c("model", "tsl", "parameters", "parametersSampled",
                           "statesSampled", "mcmc", "prior", "fit", "MLE", "call")

  type <- attr(object, "method")

  if (return.logical) {
    x <- inherits(object, "NAWRUfit") &&
      all(components[[type]] %in% names(object))
    x
  } else {
    if (!inherits(object, "NAWRUfit")) {
      stop("Object is not of class 'NAWRUfit'")
    }
    if (!all(components[[type]] %in% names(object))) {
      stop(paste0(
        "Model is not a proper object of class 'NAWRUfit'.",
        "The following components are missing: ", paste0(components[!(components %in% names(object))], collapse = ", ")
      ))
    }
  }
}


# -------------------------------------------------------------------------------------------

#' Print \code{NAWRUmodel} object
#'
#' @description Prints the model specifications of an object of class \code{NAWRUmodel}.
#'
#' @param x An object of class \code{NAWRUmodel}.
#' @param call A logical. If \code{TRUE}, the call will be printed.
#' @param check A logical. If \code{TRUE}, the model class will be checked.
#' @param ... Ignored.
#' @return No return value, model information is printed.
#' @export
print.NAWRUmodel <- function(x, call = TRUE, check = TRUE, ...) {
  .printSSModel(x = x, call = call, check = check)
}

# -------------------------------------------------------------------------------------------

#' Print \code{NAWRUfit} object
#'
#' Prints the model specifications and the estimation results of an object of class \code{NAWRUfit}.
#'
#' @param x An object of class \code{NAWRUfit}.
#' @param ... Ignored.
#' @return No return value, results are printed.
#' @export
print.NAWRUfit <- function(x, ...) {
  .printSSModelFit(x = x, call = TRUE, check = FALSE, print.model = TRUE)
}


# -------------------------------------------------------------------------------------------

#' Plots for a \code{NAWRUfit} object
#'
#' @description Plots the NAWRU and the Phillip's curve and gives diagnostic plots based on
#' standardized residuals for objects of class \code{NAWRUfit}.
#'
#' @param x An object of class \code{NAWRUfit}.
#' @param alpha The significance level for the NAWRU (\code{alpha in [0,1]}). Only used if
#'   \code{bounds = TRUE}.
#' @param bounds A logical indicating whether significance intervals should be plotted around
#'   the nawru. The default is \code{bounds = TRUE}.
#' @param combine A logical indicating whether the diagnostic plots should be combined or not,
#'   the default is \code{TRUE}.
#' @param posterior A logical indicating whether posterior diagnostics should be plotted. The
#'   default is \code{FALSE}. Only applied in the case of bayesian estimation.
#' @inheritParams plot.gap
#' 
#' @return No return value, plots are printed.
#'
#' @export
plot.NAWRUfit <- function(x, alpha = 0.05, bounds = TRUE, path = NULL, combine = TRUE, prefix = NULL,
                          posterior = FALSE, device = "png", width = 10, height = 3, ...) {
  if (!is.NAWRUfit(x, return.logical = TRUE)) {
    stop("x is no valid object of class 'NAWRUfit'.")
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
          ub = (x$tsl$nawru + x$tsl$nawruSE * tvalue),
          lb = (x$tsl$nawru - x$tsl$nawruSE * tvalue),
          ub2 = (x$tsl$obsFitted[, 2] + x$tsl$obsFittedSE[, 2] * tvalue),
          lb2 = (x$tsl$obsFitted[, 2] - x$tsl$obsFittedSE[, 2] * tvalue)        )
        
        # residuals
        res <- x$tsl$obsResidualsRecursive[, "pcInd"]
        
      } else { # Bayesian
        # confidence bounds
        HPDI <- attr(x, "HPDIprob")
        boundName <- paste0(100 * HPDI, "% HPDCI")
        tslBounds <- list(
          ub = x$tsl$nawruSummary[, paste0(100 * HPDI, "% HPDI-UB")],
          lb = x$tsl$nawruSummary[, paste0(100 * HPDI, "% HPDI-LB")],
          ub2 = x$tsl$pcIndFittedSummary[, paste0(100 * HPDI, "% HPDI-UB")],
          lb2 = x$tsl$pcIndFittedSummary[, paste0(100 * HPDI, "% HPDI-LB")]
        )
        
        # residuals
        res <- NULL
      }
  
      # --- data
      # nawru
      tsl1 <- list(
        trend = x$tsl$nawru,
        # orig = x$model$tsl$ur,
        orig = x$tsl$obs[, 1],
        lb = tslBounds$lb,
        ub = tslBounds$ub
      )
      if (!is.null(x$tsl$nawruAnchored)) {
        tsl1 <- c(tsl1, list(anchor = x$tsl$nawruAnchored))
      }
      tsl1 <- do.call(cbind, tsl1)
  
      # phillips curve
      tsl2 <- do.call(cbind, list(
        fitted = x$tsl$obsFitted[, 2],
        # E2 = x$model$tsl$pcInd,
        E2 = x$tsl$obs[, 2],
        lb = tslBounds$lb2,
        ub = tslBounds$ub2
      ))
  
      # combine lists
      tsl <- list(tsl1, tsl2)
  
      # --- legends and titles and print names
      legend <- list(
        c("nawru", "unemployment rate", "anchored nawru"),
        c("fitted", "Phillips curve indicator")
      )
      title <- list(
        "Unemployment rate in %",
        "Phillips curve",
        "Phillips curve residuals"
      )
      namesPrint <- c("nawru", "phillips_curve")
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
          lb = (x$tsl$nawru - x$tsl$nawruSE * tvalue),
          ub = (x$tsl$nawru + x$tsl$nawruSE * tvalue),
          lb1 = (x$tsl$obs[, 1] - x$tsl$obsSE[, 1] * tvalue),
          ub1 = (x$tsl$obs[, 1] + x$tsl$obsSE[, 1] * tvalue),
          lb2 = (x$tsl$obs[, 2] - x$tsl$obsSE[, 2] * tvalue),
          ub2 = (x$tsl$obs[, 2] + x$tsl$obsSE[, 2] * tvalue),
          lb3 = (x$tsl$stateSmoothed[, "cycle"] - x$tsl$stateSmoothedSE[, "cycle"] * tvalue),
          ub3 = (x$tsl$stateSmoothed[, "cycle"] + x$tsl$stateSmoothedSE[, "cycle"] * tvalue)        )
        
      } else { # Bayesian
        
        # confidence bounds
        HPDI <- attr(x, "HPDIprob")
        boundName <- paste0(100 * HPDI, "% HPDCI")
        tslBounds <- list(
          ub = x$tsl$nawruSummary[, paste0(100 * HPDI, "% HPDI-UB")],
          lb = x$tsl$nawruSummary[, paste0(100 * HPDI, "% HPDI-LB")],
          ub1 = x$tsl$obsSummary[, paste0("ur.", 100 * HPDI, "% HPDI-UB")],
          lb1 = x$tsl$obsSummary[, paste0("ur.", 100 * HPDI, "% HPDI-LB")],
          ub2 = x$tsl$obsSummary[, paste0("pcInd.", 100 * HPDI, "% HPDI-UB")],
          lb2 = x$tsl$obsSummary[, paste0("pcInd.", 100 * HPDI, "% HPDI-LB")],
          ub3 = x$tsl$stateSmoothedSummary[, paste0("cycle.", 100 * HPDI, "% HPDI-UB")],
          lb3 = x$tsl$stateSmoothedSummary[, paste0("cycle.", 100 * HPDI, "% HPDI-LB")]       )
        
      }
      
      # nawru
      tsl1 <- list(
        trend = x$tsl$nawru,
        orig = x$tsl$obs[, 1],
        lb = tslBounds$lb,
        ub = tslBounds$ub,
        lb_fc = tslBounds$lb1,
        ub_fc = tslBounds$ub1
      )
      tsl1 <- do.call(cbind, tsl1)
      
      # phillips curve
      tsl2 <- do.call(cbind, list(
        E2 = x$tsl$obs[, 2],
        lb_fc = tslBounds$lb2,
        ub_fc = tslBounds$ub2
      ))
      
      # cycle
      tsl3 <- do.call(cbind, list(
        cycle = 100 * x$tsl$stateSmoothed[, "cycle"],
        lb = tslBounds$lb3,
        ub = tslBounds$ub3
      ))
      
      # combine lists
      tsl <- list(tsl1, tsl2, tsl3)
      
      # --- legends and titles and print names
      legend <- list(
        c("nawru", "unemployment rate", rep(paste0(boundName, " (ur)"), 2)),
        c("Phillips curve indicator", rep(paste0(boundName, ""), 2)),
        c("unemployment gap")
      )
      title <- list(
        "Unemployment rate in %",
        "Phillips curve",
        "Unemployment gap"
      )
      namesPrint <-  c("nawru", "phillips_curve", "u_gap")
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

      inputl <- list(
        draws = param[, jj], burnin = burnin / thin, parName = gsub("E2", "pc", jj),
        lb = distrPar[3, jj], ub = distrPar[4, jj]
      )
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
