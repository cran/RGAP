
# -------------------------------------------------------------------------------------------

#' Estimates a two-dimensional state-space model and performs filtering and smoothing to
#' obtain the nawru.
#'
#' @inheritParams fit.NAWRUmodel
#'
#' @importFrom KFAS fitSSM KFS
#' @importFrom stats start end window ts lag frequency time Box.test coef
#' @importFrom zoo na.trim
#' @keywords internal
.MLEfitNAWRU <- function(model, parRestr = initializeRestr(model = model), signalToNoise = NULL, control = NULL) {

  if (!is.NAWRUmodel(model, return.logical = TRUE)) {
    stop("Model object is not of class 'NAWRUmodel'.")
  }

  # attributes
  trend <- attr(model, "trend")
  cycle <- attr(model, "cycle")
  type <- attr(model, "phillips curve")$type
  cycleLag <- attr(model, "phillips curve")$cycleLag
  errorARMA <- attr(model, "phillips curve")$errorARMA
  exoNames <- attr(model, "phillips curve")$exoVariables
  anchor <- attr(model, "anchor")$value
  anchor.h <- attr(model, "anchor")$horizon
  freq <- ifelse(attr(model, "period")$frequency == "quarterly", 4, 1)

  # time series list
  tslUsed <- model$tsl
  start <- start(tslUsed$ur)

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
    update_args = list(loc, cycle, trend, errorARMA, signalToNoise, type),
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
      update_args = list(loc, cycle, trend, errorARMA, signalToNoise, type),
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

  # add pre-specified parameter to output
  if (type == "NKP" && cycle == "AR2") dfRes <- rbind(pcC1 = c(0.99 * dfRes["pcC0", "Coefficient"] * dfRes["cPhi2", "Coefficient"], NA, NA, NA), dfRes)

  # order parameters
  dfRes <- dfRes[order(rownames(dfRes)), ]

  # ----- model fit
  info <- .SSmodelfit(out = out, nPar = nPar)
  info[["signal-to-noise"]] <- signalToNoise
  if (is.null(signalToNoise)) {
    info[["signal-to-noise"]] <- sum(dfRes["tdSigma", "Coefficient"], dfRes["tSigma", "Coefficient"], na.rm = TRUE) / dfRes["cSigma", "Coefficient"]
  }

  # ----- fitted nawru object
  NAWRUfit <- list(
    model = model,
    tsl = tslRes,
    SSMfit = f,
    SSMout = out,
    parameters = dfRes,
    parRestr = parRestr,
    fit = info
  )
  class(NAWRUfit) <- c("NAWRUfit", "fit")
  attr(NAWRUfit, "method") <- "MLE"

  # ----- anchor
  if (!is.null(anchor) & !is.null(anchor.h)) {
    NAWRUfit$tsl$nawruAnchored <- trendAnchor(fit = NAWRUfit, returnFit = FALSE)
  }

  invisible(NAWRUfit)
}
