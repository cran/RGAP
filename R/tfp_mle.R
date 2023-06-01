# -------------------------------------------------------------------------------------------

#' Estimates a two-dimensional state-space model and performs filtering and smoothing
#' to obtain the tfp trend.
#'
#' @inheritParams fit.TFPmodel
#'
#' @importFrom KFAS fitSSM KFS
#' @importFrom stats start end window ts lag frequency time Box.test coef
#' @importFrom zoo na.trim
#' @keywords internal
.MLEfitTFP <- function(model, parRestr = initializeRestr(model), signalToNoise = NULL, control = NULL) {

  if (!is.TFPmodel(model, return.logical = TRUE)) {
    stop("Model object is not of class 'TFPmodel'.")
  }

  # attributes
  trend <- attr(model, "trend")
  cycle <- attr(model, "cycle")
  cycleLag <- attr(model, "cubs")$cycleLag
  cubsAR <- attr(model, "cubs")$cubsAR
  errorARMA <- attr(model, "cubs")$errorARMA
  exoNames <- attr(model, "cubs")$exoNames
  anchor <- attr(model, "anchor")$value
  anchor.h <- attr(model, "anchor")$horizon
  freq <- ifelse(attr(model, "period")$frequency == "quarterly", 4, 1)

  # time series list
  tslUsed <- model$tsl
  start <- start(tslUsed$logtfp)

  # update parameter constraints
  trendVarExist <- any(grepl("tSigma", colnames(parRestr$trend)))
  if (!is.null(parRestr)) {
    .checkParRestr(model = model, parRestr = parRestr)
  }
  model <- .updateParConstraints(model = model, parRestr = parRestr)
  if (trendVarExist) {
    if (all(parRestr$trend[, "tSigma"] == 0) | !is.null(signalToNoise)) {
      parRestr$trend <- parRestr$trend[, colnames(parRestr$trend) != "tSigma", drop = FALSE]
      model <- .modifySSSystem(model, signalToNoise = signalToNoise)
    }
  }

  # parameter names, locations and constraints
  loc <- model$loc
  nPar <- dim(loc)[1]

  # ----- estimate system

  # initial parameters
  initPar <- rep(0.1, nPar)
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

  # filtered and smoothed time series and residuals
  tslRes <- .SSresults(out = out, model = model)

  # ----- parameter inference
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
  TFPfit <- list(
    model = model,
    tsl = tslRes,
    SSMfit = f,
    SSMout = out,
    parameters = dfRes,
    parRestr = parRestr,
    fit = info
  )
  class(TFPfit) <- c("TFPfit", "fit")
  attr(TFPfit, "method") <- "MLE"

  # ----- anchor
  if (!is.null(anchor) & !is.null(anchor.h)) {
    TFPfit$tsl$tfpTrendAnchored <- trendAnchor(fit = TFPfit, returnFit = FALSE)
  }

  invisible(TFPfit)
}
