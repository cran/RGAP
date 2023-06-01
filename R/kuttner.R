
# -------------------------------------------------------------------------------------------

#' Kuttner model
#'
#' @description Creates a state space object object of class \code{KuttnerModel} which can be 
#'   fitted using \code{fit}.
#'
#' @param tsl A list of time series objects, see details.
#' @param trend A character string specifying the trend model. \code{trend = "RW1"} denotes
#'   a first order random walk, \code{trend = "RW2"} a second order random walk (local linear
#'   trend) and \code{trend = "DT"} a damped trend model. The default is \code{trend = "RW1"}.
#' @param cycle A character string specifying the cycle model. \code{cycle = "AR1"} denotes
#'   an AR(1) process, \code{cycle = "AR2"} an AR(2) process. The default is
#'   \code{cycle = "AR2"}.
#' @param inflErrorARMA A \code{2 x 1} vector with non-negative integers specifying the AR
#'   and MA degree of the error term in the inflation equation. The default is
#'   \code{inflErrorARMA = c(0, 3)}, see details.
#' @param cycleLag A non-negative integer specifying the maximum cycle lag that is included
#'   in the inflation equation. The default is \code{cycleLag = 0}, see details.
#' @param start (Optional) Start vector for the estimation, e.g. \code{c(1980, 1)}.
#' @param end (Optional) End vector for the estimation, e.g. \code{c(2020, 1)}.
#' @param anchor (Optional) Anchor value for the logarithm of trend gdp.
#' @param anchor.h (Optional) Anchor horizon in the frequency of the given time series.
#'
#' @details The list of time series \code{tsl} needs to have the following components:
#' \describe{
#'   \item{gdp}{Real gross domestic product.}
#'   \item{infl}{Inflation.}
#'   }
#' @details A \code{cycleLag} equal to \code{0} implies that only the contemporaneous cycle
#'   is included in the inflation equation.  A \code{cycleLag} equal to \code{0:1} implies that
#'   the contemporaneous as well as the lagged cycle are included.
#' @details A \code{inflErrorARMA} equal to \code{c(0, 0)} implies that the error term in the
#'   inflation equation is white noise. \code{inflErrorARMA = c(1, 0)} implies that the error is
#'   an AR(1) process and for \code{inflErrorARMA = c(1, 2)} the error follows an ARMA(1, 2)
#'   process.
#'
#' @return Object of class \code{KuttnerModel}, which is a list with the following components:
#'   \item{tsl}{A list of used time series.}
#'   \item{SSModel}{An object of class SSModel specifying the state-space model.}
#'   \item{loc}{A data frame containing information on each involved parameter, for instance
#'         its corresponding system matrix, variable names, and parameter restrictions.}
#'   \item{call}{Original call to the function. }
#'   In addition, the object contains the following attributes:
#'   \item{cycle}{Cycle specification.}
#'   \item{trend}{Trend specification.}
#'   \item{inflation equation}{A list containing the components \code{cycleLag, errorARMA, exoVariables}.}
#'   \item{anchor}{A list containing the components \code{value, horizon}.}
#'   \item{period}{A list containing the components \code{start, end, frequency}.}
#'
#' @export
#' @importFrom KFAS SSModel SSMregression SSMcustom
#' @importFrom stats start end window ts lag frequency time
#' @importFrom zoo na.trim
#' @examples
#' # load data for the Netherlands
#' data("gap")
#' country <- "Netherlands"
#' tsList <- as.list(gap[[country]][, c("cpih", "gdp")])
#' tsList$infl <- diff(tsList$cpih)
#' model <- KuttnerModel(tsl = tsList, trend = "RW2", start = 1980)
KuttnerModel <- function(tsl, cycle = "AR2", cycleLag = 1, trend = "RW1", inflErrorARMA = c(0, 3),
                         start = NULL, end = NULL, anchor = NULL, anchor.h = NULL) {

  # local variables
  nPar <- nObs <- nState <- nStateV <- stateNames <- loc <- sys <- NULL

  # save call
  mc <- match.call(expand.dots = FALSE)

  # ----- check input for consistency
  errorARMA <- inflErrorARMA
  errorAR <- errorARMA[1]
  inflMA <- errorARMA[2]
  list2env(.checkKuttner(
    tsl = tsl,
    trend = trend,
    cycle = cycle,
    cycleLag = cycleLag,
    errorARMA = errorARMA,
    start = start,
    end = end,
    anchor = anchor,
    anchor.h = anchor.h
  ),
  envir = environment())

  # ----- preprocess data

  # assign end and start
  if (is.null(end)) end <- end(tsl$gdp)
  if (is.null(start)) start <- start(tsl$gdp)
  # collect used variables
  varUsed <- c("ur")
  nExo <- 0
  exoNamesTmp <- NULL
  # merge observation equation data
  tslUsed <- list()
  tslUsed$loggdp <- log(tsl$gdp) * 100
  tslUsed$dinfl <- diff(tsl$infl)
  tslUsed$gdpGL1 <- stats::lag(diff(log(tsl$gdp)), k = -1) * 100
  obsNames <- varUsed <- c("loggdp", "dinfl")
  nExo <- 1
  exoNamesTmp <- "gdpGL1"
  tslUsed <- tslUsed[c(obsNames, exoNamesTmp)]

  # list of used time series
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
    type = NULL,
    errorARMA = errorARMA
  ),
  envir = environment())

  tslUsed <- lapply(tslUsed, function(x) suppressWarnings(window(x, start = start, end = end))) # warns if start or end are not changed.

  # ----- state space model
  tmp <- paste0("~ ", paste(exoNamesTmp, collapse = " + "))
  modelSS <- SSModel(cbind(loggdp, dinfl) ~ -1
    + SSMregression(as.formula(tmp), index = 2, data = tslUsed, state_names = exoNamesTmp)
      + SSMcustom(
        Z = sys$Zt, T = sys$Tt, R = sys$Rt,
        Q = sys$Qt, a1 = sys$a1, P1 = sys$P1, P1inf = sys$P1inf, state_names = stateNames
      ),
  H = sys$Ht,
  data = tslUsed
  )


  # -----  Kuttner model
  model <- list(
    tsl = tslUsed,
    SSModel = modelSS,
    call = mc
  )
  lInfl <- list(
    cycleLag = cycleLag,
    errorARMA = errorARMA,
    exoVariables = exoNamesTmp
  )
  lAnchor <- list(
    value = anchor,
    horizon = anchor.h
  )
  lPeriod <- list(
    start = paste0(start(tslUsed$loggdp)[1], ifelse(frequency(tslUsed$loggdp) == 4, paste0(" Q", start(tslUsed$loggdp)[2]), "")),
    end = paste0(end(tslUsed$loggdp)[1], ifelse(frequency(tslUsed$loggdp) == 4, paste0(" Q", end(tslUsed$loggdp)[2]), "")),
    frequency = ifelse(frequency(tslUsed$loggdp) == 4, "quarterly", "annual")
  )
  class(model) <- c("KuttnerModel", "model")
  attr(model, "cycle") <- cycle
  attr(model, "trend") <- trend
  attr(model, "inflation equation") <- lInfl
  attr(model, "anchor") <- lAnchor
  attr(model, "period") <- lPeriod
  model$loc <- .initializeLoc(model)

  invisible(return(model))
}

# -------------------------------------------------------------------------------------------

#' \code{KuttnerModel} object check
#'
#' @description Tests whether the input object is a valid object of class \code{KuttnerModel}.
#'
#' @param object An object to be tested.
#' @param return.logical If \code{return.logical = FALSE} (default), an error message is printed
#'   if the object is not of class \code{KuttnerModel}. If \code{return.logical = TRUE}, a logical
#'   value is returned.
#'
#' @return A logical value or nothing, depending on the value of \code{return.logical}.
#'
#' @export
#' @importFrom KFAS is.SSModel
is.KuttnerModel <- function(object, return.logical = FALSE) {
  cycle <- attr(object, "cycle")
  trend <- attr(object, "trend")
  cycleLag <- attr(object, "inflation equation")$cycleLag
  errorARMA <- attr(object, "inflation equation")$errorARMA
  anchor <- attr(object, "anchor")$value
  anchor.h <- attr(object, "anchor")$horizon

  components <- c("tsl", "SSModel", "loc", "call")
  trendPossibilities <- c("RW1", "RW2", "DT")
  cyclePossibilities <- c("AR1", "AR2", "RAR2")
  cycleLagPossibilities <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  errorARPossibilities <- 0
  errorMAPossibilities <- 1:3
  obsNames <- c("loggdp", "dinfl")

  if (return.logical) {
    x <- inherits(object, "KuttnerModel") &&
      all(components %in% names(object)) &&
      inherits(object$SSModel, "SSModel") &&
      is.SSModel(object$SSModel, return.logical = TRUE) &&
      all(trend %in% trendPossibilities) &&
      all(cycle %in% cyclePossibilities) &&
      all(cycleLag %in% cycleLagPossibilities) &&
      all(errorARMA[1] %in% errorARPossibilities) &&
      all(errorARMA[2] %in% errorMAPossibilities) &&
      !(!is.null(anchor) && is.null(anchor.h)) &&
      all(obsNames %in% names(object$tsl))
    x
  } else {
    if (!inherits(object, "KuttnerModel")) {
      stop("Object is not of class 'KuttnerModel'")
    }
    if (!all(components %in% names(object))) {
      stop(paste0(
        "Model is not a proper object of class 'KuttnerModel'.",
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

#' \code{KuttnerFit} object check
#'
#' @description Tests whether the input object is a valid object of class \code{KuttnerFit}.
#'
#' @param object An object to be tested.
#' @param return.logical If \code{return.logical = FALSE} (default), an error message is printed
#'   if the object is not of class \code{KuttnerFit}. If \code{return.logical = TRUE}, a logical
#'   value is returned.
#'
#' @return A logical value or nothing, depending on the value of \code{return.logical}.
#'
#' @keywords internal
is.KuttnerFit <- function(object, return.logical = FALSE) {
  components <- list()
  components$MLE <- c("model", "tsl", "SSMfit", "SSMout", "parameters", "parRestr", "fit", "call")

  type <- attr(object, "method")

  if (return.logical) {
    x <- inherits(object, "KuttnerFit") &&
      all(components[[type]] %in% names(object))
    x
  } else {
    if (!inherits(object, "KuttnerFit")) {
      stop("Object is not of class 'KuttnerFit'")
    }
    if (!all(components[[type]] %in% names(object))) {
      stop(paste0(
        "Model is not a proper object of class 'KuttnerFit'.",
        "The following components are missing: ", paste0(components[!(components %in% names(object))], collapse = ", ")
      ))
    }
  }
}

# -------------------------------------------------------------------------------------------

#' Print \code{KuttnerModel} object
#'
#' @description Prints the model specifications of an object of class \code{KuttnerModel}.
#'
#' @param x An object of class \code{KuttnerModel}.
#' @param call A logical. If \code{TRUE}, the call will be printed.
#' @param check A logical. If \code{TRUE}, the model class will be checked.
#' @param ... Ignored.
#' @return No return value, model information is printed.
#' @export
print.KuttnerModel <- function(x, call = TRUE, check = TRUE, ...) {
  .printSSModel(x = x, call = call, check = check)
}

# -------------------------------------------------------------------------------------------

#' Print \code{KuttnerFit} object
#'
#' @description Prints the model specifications and the estimation results of an object of
#'   class \code{KuttnerFit}.
#'
#' @param x An object of class \code{KuttnerFit}.
#' @param ... Ignored.
#' @return No return value, results are printed.
#' @export
print.KuttnerFit <- function(x, ...) {
  .printSSModelFit(x = x, call = TRUE, check = FALSE, print.model = TRUE)
}

# -------------------------------------------------------------------------------------------

#' Plots for a \code{KuttnerFit} object
#'
#' @description Plots potential growth and the output gap and gives diagnostic plots based on
#' standardized residuals for objects of class \code{KuttnerFit}.
#'
#' @param x An object of class \code{KuttnerFit}.
#' @param alpha The significance level for the trend (\code{alpha in [0,1]}). Only used if
#'   \code{bounds = TRUE}.
#' @param bounds A logical indicating whether significance intervals should be plotted around
#'   gdp. The default is \code{bounds = TRUE}.
#' @param combine A logical indicating whether the diagnostic plots should be combined or not,
#'   the default is \code{TRUE}.
#' @inheritParams plot.gap
#'
#' @return No return value, plots are printed.
#' 
#' @export
plot.KuttnerFit <- function(x, alpha = 0.05, bounds = TRUE, path = NULL, combine = TRUE,
                            prefix = NULL, device = "png", width = 10, height = 3, ...) {
  # potential <- gap <- lb <- ub <- NULL

  if (!is.KuttnerFit(x, return.logical = TRUE)) {
    stop("x is no valid object of class 'KuttnerFit'.")
  }

  if (!is.null(path)) {
    check <- dir.exists(paths = path)
    if (!check) {
      warning(paste0("The given file path '", path, "' does not exist. The plots will not be saved."))
      path <- NULL
    }
  }
  
  # check whether prediction is present
  prediction <- "prediction" %in% names(attributes(x))
  
  # confidence bounds
  tvalue <- -qnorm((alpha) / 2)
  boundName <- paste0(100 * (1 - alpha), "% CI")
  
  
  if (!prediction) {
  
  
    # ----- SSModel plots
  
    # residuals
    res <- x$tsl$obsResidualsRecursive[, "dinfl"]
  
    # --- data
    # potential growth
    tsl1 <- list(
      trend = x$tsl$potential,
      orig = exp(x$tsl$obs[, 1] / 100),
      lb = (x$tsl$potential - x$tsl$potentialSE * tvalue),
      ub = (x$tsl$potential + x$tsl$potentialSE * tvalue)
    )
    if (!is.null(x$tsl$trendAnchored)) {
      tsl1 <- c(tsl1, list(anchor = exp(x$tsl$trendAnchored / 100)))
    }
    tsl1 <- do.call(cbind, tsl1)
  
    # cubs equation
    tsl2 <- do.call(cbind, list(
      fitted = x$tsl$obsFitted[, 2],
      E2 = x$tsl$obs[, 2],
      lb = (x$tsl$obsFitted[, 2] - x$tsl$obsFittedSE[, 2] * tvalue),
      ub = (x$tsl$obsFitted[, 2] + x$tsl$obsFittedSE[, 2] * tvalue)
    ))
  
    # combine lists
    tsl <- list(tsl1, tsl2)
  
    # --- legends and titles and print names
    legend <- list(
      c("potential", "gdp", "anchored potential"),
      c("fitted", "change in inflations")
    )
    title <- list(
      "Potential output",
      "Inflation",
      "Inflation residuals"
    )
    namesPrint <- paste(prefix, c("potential_growth", "inflation"), sep = "_")
  
    # plot
    plotSSresults(
      tsl = tsl, legend = legend, title = title,
      boundName = boundName, res = res, namesPrint = namesPrint,
      bounds = bounds, combine = combine, path = path, device = device,
      width = width, height = height
    )
  
    # ----- gap plots
  
    # --- data
    # potential and gdp growth
    tsl1 <- list(
      potential = 100 * x$tsl$potentialGrowth,
      gdp = 100 * growth(window(exp(x$model$tsl$loggdp / 100), start = start(x$tsl$potential))),
      lb = 100 * (x$tsl$potentialGrowth - x$tsl$potentialGrowthSE * tvalue),
      ub = 100 * (x$tsl$potentialGrowth + x$tsl$potentialGrowthSE * tvalue)
    )
    tsl1 <- do.call(cbind, tsl1)
  
    # gap
    tsl2 <- do.call(cbind, list(
      gap = x$tsl$gap,
      lb = (x$tsl$gap - x$tsl$gapSE * tvalue),
      ub = (x$tsl$gap + x$tsl$gapSE * tvalue)
    ))
  
    # combine lists
    tsl <- list(tsl1, tsl2)
  
    # --- legends and titles and print names
    legend <- list(
      c("potential", "gdp"),
      c("output gap")
    )
    title <- list(
      "GDP growth in %",
      "Output gap in %"
    )
    namesPrint <-  c("potential_growth", "gap")
    if (!is.null(prefix)) namesPrint <- paste(prefix, namesPrint , sep = "_")

    # plot
    plotGap(
      tsl = tsl, legend = legend, title = title, boundName = boundName,
      contribution = FALSE, res = res, namesPrint = namesPrint,
      bounds = bounds, combine = combine, path = path, device = device,
      width = width, height = height
    )
    
  } else { # plot predictions
    
    n.ahead <- attr(x, "prediction")$n.ahead
    
    # confidence bounds
    tvalue <- -qnorm((alpha) / 2)
    boundName <- paste0(100 * (1 - alpha), "% CI")

    # residuals
    res <- NULL
    
    # --- data
    # potential 
    tsl1 <- do.call(cbind, list(
      trend = x$tsl$potential,
      orig = x$tsl$gdp,
      lb = (x$tsl$potential - x$tsl$potentialSE * tvalue),
      ub = (x$tsl$potential + x$tsl$potentialSE * tvalue),
      lb_fc = (x$tsl$gdp - x$tsl$gdpSE * tvalue),
      ub_fc = (x$tsl$gdp + x$tsl$gdpSE * tvalue)
    ))
    
    # inflation equation
    tsl2 <- do.call(cbind, list(
      E2 = x$tsl$obs[, 2],
      lb_fc = (x$tsl$obs[, 2] - x$tsl$obsSE[, 2] * tvalue),
      ub_fc = (x$tsl$obs[, 2] + x$tsl$obsSE[, 2] * tvalue)
    ))
    
    # cycle
    tsl3 <- do.call(cbind, list(
      cycle = x$tsl$stateSmoothed[, "cycle"],
      lb = (x$tsl$stateSmoothed[, "cycle"] - x$tsl$stateSmoothedSE[, "cycle"] * tvalue),
      ub = (x$tsl$stateSmoothed[, "cycle"] + x$tsl$stateSmoothedSE[, "cycle"] * tvalue)
    ))
    
    # potential growth
    tsl4 <- NULL
    tsl4 <- do.call(cbind, list(
      trend = 100 * x$tsl$potentialGrowth,
      orig = 100 * x$tsl$gdpGrowth,
      lb = 100 * (x$tsl$potentialGrowth - x$tsl$potentialGrowthSE * tvalue),
      ub = 100 * (x$tsl$potentialGrowth + x$tsl$potentialGrowthSE * tvalue),
      lb_fc = 100 * (x$tsl$gdpGrowth - x$tsl$gdpGrowthSE * tvalue),
      ub_fc = 100 * (x$tsl$gdpGrowth + x$tsl$gdpGrowthSE * tvalue)
    ))
    
    # combine lists
    tsl <- list(tsl1, tsl2, tsl3, tsl4)
    
    # --- legends and titles and print names
    legend <- list(
      c("potential", "gdp", rep(paste0(boundName, " (gdp)"), 2)),
      c("change in inflation", rep(paste0(boundName, ""), 2)),
      c("output gap"),
      c("potential growth", "gdp growth", rep(paste0(boundName, " (gdp growth)"), 2))
    )
    title <- list(
      "Potential output",
      "Inflation",
      "Output gap",
      "Potential output growth"
    )
    namesPrint <-  c("potential", "inflation", "gap")
    if (!is.null(prefix)) namesPrint <- paste(prefix, namesPrint , sep = "_")
    
    # plot
    plotSSprediction(
      tsl = tsl, legend = legend, title = title, n.ahead = n.ahead,
      boundName = boundName, res = NULL, namesPrint = namesPrint,
      bounds = bounds, combine = combine, path = path, device = device,
      width = width, height = height
    )
    
  }
}
