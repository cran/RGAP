
# -------------------------------------------------------------------------------------------

#' Maximum likelihood estimation of a \code{KuttnerModel}
#'
#' @description Estimates a two-dimensional state-space model and performs filtering and
#'   smoothing to obtain the output gap.
#'
#' @param model An object of class KuttnerModel.
#' @param parRestr A list of matrices containing the parameter restrictions for the cycle,
#'   trend, and the inflation equation. Each matrix contains the lower and upper bound of the
#'   involved parameters. \code{NA} implies that no restriction is present. Autoregressive
#'   parameters are automatically restricted to the stationary region unless box constraints
#'   are specified. By default, \code{parRestr} is intitialized by the function
#'   \code{initializeRestr(model)}.
#' @param signalToNoise (Optional) signal to noise ratio.
#' @param control (Optional) A list of control arguments to be passed on to \code{optim}.
#' @param ... additional arguments to be passed to the methods functions.
#'
#' @return An object of class \code{KuttnerFit} containing the following components:
#'   \item{model}{The input object of class \code{KuttnerModel}.}
#'   \item{SSMfit}{The estimation output from the funtcion \code{fitSSM} from \code{KFAS}.}
#'   \item{SSMout}{The filtering and smoothing output from the funtcion \code{KFS} from
#'         \code{KFAS}.}
#'   \item{parameters}{A data frame containing the estimated parameters, including standard
#'         errors, t-statistics, and p-values.}
#'   \item{fit}{A list of model fit criteria (see below).}
#'   \item{call}{Original call to the function. }
#'   The list component \code{fit} contains the following model fit criteria:
#'   \item{loglik}{Log-likelihood function value.}
#'   \item{AIC}{Akaike information criterion.}
#'   \item{BIC}{Bayesian information criterion.}
#'   \item{AIC}{Hannan-Quinn information criterion.}
#'   \item{RMSE}{Root mean squared error of the inflation equation.}
#'   \item{R2}{R squared of the inflation equation.}
#'   \item{LjungBox}{Ljung-Box test output of the inflation equation.}
#'
#' @export
#' @importFrom KFAS fitSSM KFS
#' @importFrom stats start end window ts lag frequency time Box.test coef
#' @importFrom zoo na.trim
#' @family fitting methods
#' @examples
#' # load data for the Netherlands
#' data("gap")
#' country <- "Netherlands"
#' tsList <- as.list(gap[[country]][, c("cpih", "gdp")])
#' tsList$infl <- diff(tsList$cpih)
#' model <- KuttnerModel(tsl = tsList, trend = "RW2", cycleLag = 1, cycle = "AR2", start = 1980)
#' 
#' # estimate Kutter's model
#' parRestr <- initializeRestr(model = model, type = "hp")
#' \donttest{
#' gapKuttner <- fit(model, parRestr, signalToNoise = 1 / 10)
#' }
fit.KuttnerModel <- function(model, parRestr = initializeRestr(model), signalToNoise = NULL, control = NULL, ...) {

  # save call
  mc <- match.call(expand.dots = FALSE)

  # attributes
  trend <- attr(model, "trend")
  cycle <- attr(model, "cycle")
  exoNames <- attr(model, "inflation equation")$exoVariables
  errorARMA <- attr(model, "inflation equation")$errorARMA
  anchor <- attr(model, "anchor")$value
  anchor.h <- attr(model, "anchor")$horizon
  freq <- ifelse(attr(model, "period")$frequency == "quarterly", 4, 1)

  # time series list
  tslUsed <- model$tsl
  start <- start(tslUsed$loggdp)

  # update parameter constraints
  trendVarExist <- any(grepl("tSigma", colnames(parRestr$trend)))
  if (!is.null(parRestr)) {
    .checkParRestr(model = model, parRestr = parRestr)
  }
  model <- .updateParConstraints(model = model, parRestr = parRestr)
  if (trendVarExist) {
    # modify state system if necessary
    if (all(parRestr$trend[, "tSigma"] == 0) | !is.null(signalToNoise)) {
      parRestr$trend <- parRestr$trend[, colnames(parRestr$trend) != "tSigma", drop = FALSE]
      model <- .modifySSSystem(model, signalToNoise = signalToNoise)
    }
  }

  # parameter names, locations and constraints
  loc <- model$loc
  nPar <- dim(loc)[1]
  colnames(parRestr) <- cbind(loc$Q)[2, ]

  # ----- estimate system

  # initial parameters
  initPar <- rep(0, nPar)
  names(initPar) <- loc$varName

  # control settings for optim
  controlBase <- list(
    maxit = 2000,
    trace = 1,
    REPORT = 50,
    parscale = rep(0.1, nPar)
  )
  if (!is.null(control)) {
    control <- c(controlBase[!(names(controlBase) %in% names(control))], control)
  } else {
    control <- controlBase
  }
  controlList <- control

  # optimiziation
  message("Starting optimization ...")
  f <- fitSSM(
    model$SSModel,
    inits = initPar,
    method = "BFGS",
    updatefn = .updateSSSystem,
    update_args = list(loc, cycle, trend, errorARMA, signalToNoise),
    hessian = TRUE,
    control = controlList
  )
  # check covariance matrix
  rerun <- .checkCV(fit = f, loc = loc)
  loc <- .checkBoundaries(fit = f, loc = loc)
  if (rerun) {
    # check Fisher information
    message(
      "The covariance matrix contains negative values on its diagonal. \n",
      "Convergence tolerance parameter decreased."
    )
    controlList$reltol <- 0.01 * sqrt(.Machine$double.eps)
  } else if (any(loc$boundaries)) {
    rerun <- TRUE
    message(
      "Box constraints have been reached. \n",
      "Convergence tolerance parameter decreased."
    )
    controlList$reltol <- 0.01 * sqrt(.Machine$double.eps)
  }

  # rerun if necessary
  if (rerun) {
    message("Rerunning optimization ...")
    f <- fitSSM(
      model$SSModel,
      inits = initPar,
      method = "BFGS",
      updatefn = .updateSSSystem,
      update_args = list(loc, cycle, trend, errorARMA, signalToNoise),
      hessian = TRUE,
      control = controlList
    )
    rerun <- .checkCV(fit = f, loc = loc)
  }
  if (f$optim.out$convergence != 0) {
    message("The optimization could NOT be completed.")
  } else {
    message("The optimization was successfully completed.")
  }

  # check covariance matrix
  if (rerun) {
    message("The covariance matrix contains negative values on its diagonal.")
  }
  # check boundaries
  loc <- .checkBoundaries(fit = f, loc = loc)
  if (any(loc$boundaries)) {
    message("Box constraints have been reached. Consider modifying the constraints.")
  }

  # ----- filtering and smoothing
  out <- KFS(f$model, simplify = FALSE, filtering = c("state", "signal"), smoothing = c("state", "signal", "disturbance"))

  # non exogenous and non constant states
  indexTmp <- (!colnames(coef(out)) %in% c(exoNames, "const"))
  namesState <- colnames(coef(out))[indexTmp]
  namesObs <- colnames(f$model$y)

  # save filtered and smoothed time series and residuals
  tslRes <- .SSresults(out = out, model = model)

  # ----- inference
  dfRes <- inference(parOptim = f$optim.out$par, hessian = f$optim.out$hessian, loc = loc)

  # order parameters
  dfRes <- dfRes[order(rownames(dfRes)), ]

  # ----- model fit
  info <- .SSmodelfit(out = out, nPar = nPar)
  info[["signal-to-noise"]] <- signalToNoise
  if (is.null(signalToNoise)) {
    info[["signal-to-noise"]] <- sum(dfRes["tdSigma", "Coefficient"], dfRes["tSigma", "Coefficient"], na.rm = TRUE) / dfRes["cSigma", "Coefficient"]
  }

  # ----- fitted nawru object
  KuttnerFit <- list(
    model = model,
    tsl = tslRes,
    SSMfit = f,
    SSMout = out,
    parameters = dfRes,
    parRestr = parRestr,
    fit = info,
    call = mc
  )
  class(KuttnerFit) <- c("KuttnerFit", "fit")
  attr(KuttnerFit, "method") <- "MLE"

  # ----- anchor
  if (!is.null(anchor) & !is.null(anchor.h)) {
    KuttnerFit$tsl$trendAnchored <- trendAnchor(fit = KuttnerFit, returnFit = FALSE)
  }

  print(KuttnerFit)
  return(KuttnerFit)
}
