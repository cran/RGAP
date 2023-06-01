# -------------------------------------------------------------------------------------------

#' Production function output gap
#'
#' @description Computes potential output and the output gap based on a production function
#'   methodology.
#'
#' @param tsl A list of time series objects, see details.
#' @param NAWRUfit An object of class \code{NAWRUfit} obtained via the function \code{fit}.
#' @param TFPfit An object of class \code{TFPfit} obtained via the function \code{fit}.
#' @param alpha A scalar between zero and one depicting the labor share. The default is
#'   \code{alpha = 0.65}.
#' @param lambda The smoothing parameter for the application of the HP filter (see details).
#'   If not supplied, \code{lambda = 6.25} for yearly data, \code{lambda = 1600} for quarterly
#'   data, and \code{lambda = 129600} for monthly data.
#' @param start (optional) A two element vector containing a year and a period specifying
#'   the start point for the estimation.
#' @param end (optional) A two element vector containing a year and a period specifying the
#'   end point for the estimation.
#'
#' @details The list of time series \code{tsl} needs to have the following components:
#' \describe{
#'   \item{lfnd}{Labor force non-domestic (unit: 1000 persons). (Set to zero if left unspecified).}
#'   \item{parts}{Participation rate.}
#'   \item{ahours}{Average hours worked (unit: hours).}
#'   \item{gdp}{Gross domestic product at constant prices (unit: bn National currency, code: OVGD).}
#'   \item{k}{Net capital stock at constant prices: total economy (unit: bn National currency, code: OKND).}
#'   \item{popw}{Population: 15 to 64 years (unit: 1000 persons, code: NPAN).}
#'   }
#'
#' @details The trend of the list components \code{parts, ahours and lfnd} (if available)
#'   is computed using the Hodrick-Prescott filter with the smoothing constant
#'   \code{lambda}, unless the supplied time series list \code{tsl} contains their trend 
#'   (for instance, denoted by \code{partsTrend}).
#'
#' @return Object of class \code{gap}, which is a list with the following components:
#'   \item{tsl}{List of time series including potential output \code{potential}, the
#'   output gap \code{gap}, all HP-filtered trend series, and all original series.}
#'   \item{NAWRUfit}{Provided \code{NAWRUfit} object.}
#'   \item{TFPfit}{Provided \code{TFPfit} object.}
#'   \item{call}{Original call to the function.}
#'
#' @export
#' @importFrom stats start end window ts lag frequency time
#' @examples
#' # compute the output gap given the previously obtained nawru and trend tfp

#' data("gap")
#' country <- "Belgium"
#' tsList <- amecoData2input(gap[[country]])
#' modelNAWRU <- NAWRUmodel(tsl = tsList)
#' modelTFP <- TFPmodel(tsl = tsList, cycle = "RAR2")
#' 
#' \donttest{
#' fittedNAWRU <- fit(model = modelNAWRU)
#' fittedTFP <- fit(model = modelTFP)
#' 
#' gapProd(tsl = tsList, NAWRUfit = fittedNAWRU, TFPfit = fittedTFP)
#' }
gapProd <- function(tsl, NAWRUfit, TFPfit, alpha = 0.65, start = NULL, end = NULL, lambda = NULL) {

  # save call
  mc <- match.call(expand.dots = FALSE)

  # smoothing constant
  if (is.null(lambda)) {
    freq <- frequency(NAWRUfit$model$tsl[[1]])
    lambda <- 1600 / ((4 / freq)^4)
  }

  # nawru
  tsl$nawru <- NAWRUfit$tsl$nawru
  if (!is.null(NAWRUfit$tsl$nawruAnchored)) tsl$nawru <- NAWRUfit$tsl$nawruAnchored

  # tfp
  tsl$tfpTrend <- TFPfit$tsl$tfpTrend

  # check for lfnd
  if (is.null(tsl$lfnd)) {
    tsl$lfnd <- ts(0, start = start(tsl$gdp), end = end(tsl$gdp))
  }

  # apply HP filter
  namesHP <- c("parts", "lfnd", "ahours")
  # tsl[paste0(namesHP, "Trend")] <- lapply(tsl[namesHP], hpfilter, lambda)
  for (k in namesHP) {
    k_trend <- paste0(k, "Trend")
    if (!(k_trend %in% names(tsl))) {
      tsl[[k_trend]] <- hpfilter(x = tsl[[k]], lambda = lambda)
    }
  }
  

  # trend labor
  tsl$lTrend <- ((tsl$popw * tsl$partsTrend) * (1 - tsl$nawru / 100) + tsl$lfndTrend) * tsl$ahoursTrend / 1000

  # potential gdp
  # tsl$potential <- 1000 * (tsl$lTrend)^alpha * (tsl$k / 1000)^(1 - alpha) * tsl$tfpTrend
  tsl$potential <- 1 / 1000 * (tsl$lTrend)^alpha * (tsl$k * 1000)^(1 - alpha) * tsl$tfpTrend

  # output gap
  tsl$gap <- (tsl$gdp / tsl$potential - 1) * 100

  # -----  gap object
  gap <- list(
    tsl = tsl,
    NAWRUfit = NAWRUfit,
    TFPfit = TFPfit,
    call = mc
  )
  class(gap) <- c("gap", "fit")
  attr(gap, "type") <- "prod"
  attr(gap, "alpha") <- alpha
  gap
}


# -------------------------------------------------------------------------------------------

#' \code{gap} object check
#'
#' @description Tests whether the input object is a valid object of class \code{gap}.
#'
#' @param object An object to be tested.
#' @param return.logical If \code{return.logical = FALSE} (default), an error message is printed
#'   if the object is not of class \code{gap}. If \code{return.logical = TRUE}, a logical
#'   value is returned.
#'
#' @return A logical value or nothing, depending on the value of \code{return.logical}.
#'
#' @export
is.gap <- function(object, return.logical = FALSE) {
  components <- list()
  components$prod <- c("tsl", "NAWRUfit", "TFPfit", "call")
  components$HP <- c("potential", "gap", "gdp")

  type <- attr(object, "type")

  if (return.logical) {
    x <- inherits(object, "gap") &&
      all(components[[type]] %in% names(object))
    x
  } else {
    if (!inherits(object, "gap")) {
      stop("Object is not of class 'gap'")
    }
    if (!all(components[[type]] %in% names(object))) {
      stop(paste0(
        "Model is not a proper object of class 'gap'.",
        "The following components are missing: ", paste0(components[!(components %in% names(object))], collapse = ", ")
      ))
    }
  }
}


# -------------------------------------------------------------------------------------------

#' Print \code{gap} object
#'
#' @description Prints the model specifications of an object of class \code{gap}.
#'
#' @param x An object of class \code{gap}.
#' @param ... Ignored.
#' 
#' @return No return value, results are printed.
#' @export
print.gap <- function(x, ...) {
  type <- attr(x, "type")
  if (type == "prod") {
    cat("Call:\n", paste0(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n")
    print(x$NAWRUfit)
    cat("\n\n")
    print(x$TFPfit)
  }
}


# ----------------------------------------------------------------------------------------------------------------

#' HP-filter output gap
#'
#' @description Computes a HP filtered output gap.
#'
#' @param x A time series object containing gdp.
#' @param lambda The smoothing parameter for the application of the HP filter. If not supplied,
#'   \code{lambda = 6.25} for yearly data, \code{lambda = 1600} for quarterly data, and \code{lambda = 129600}
#'   for monthly data.
#' @param start (optional) A two element vector containing a year and a period specifying the start point for the
#'   filter application.
#' @param end (optional) A two element vector containing a year and a period specifying the end point for the
#'   filter application.
#'
#' @return Object of class \code{gap}, which is a list containing the two elements \code{potential} and
#'   \code{gap} and additionally the original time series.
#'
#' @export
gapHP <- function(x, lambda = NULL, end = NULL, start = NULL) {
  if (is.null(lambda)) {
    freq <- frequency(x)
    lambda <- 1600 / ((4 / freq)^4)
  }
  # select window
  if (is.null(end)) end <- end(x)
  if (is.null(start)) start <- start(x)
  x <- window(x, start, end)

  # potential and gap
  potential <- hpfilter(x, lambda = lambda)
  gap <- 100 * (x / potential - 1)

  # return result
  res <- list(
    potential = potential,
    gap = gap,
    gdp = x
  )
  class(res) <- c("gap", "fit")
  attr(res, "type") <- "HP"
  res
}


# -------------------------------------------------------------------------------------------

#' Plots for a \code{gap} object
#'
#' @description Plots potential output growth and the output gap based on an objects of
#'   class \code{gap}.
#'
#' @param x An object of class \code{gap}.
#' @param contribution A boolean indicating whether the contributions to potential output
#'   growth and the output gap should be plotted (only applicable for production function
#'   type output gaps).
#' @param path An optional file path. If specified, the plots will be saved using the format
#'   in \code{device} under the given path.
#' @param combine A logical indicating whether the plots should be combined or not, the
#'   default is \code{TRUE}.
#' @param prefix An optional character string to be added to the names of the plots in case
#'   \code{path} is specified.
#' @param device Device passed on to \code{ggplot} for plot saving. Options are 
#'   \code{"eps", "ps", "tex"} (pictex), \code{"pdf", "jpeg", "tiff", "png", "bmp", "svg"} 
#'   or \code{"wmf"}.
#' @param width The plot width in case of printing.
#' @param height The plot height in case of printing.
#' @param ... Ignored.
#' 
#' @return No return value, plots are printed.
#'
#' @export
plot.gap <- function(x, contribution = FALSE, path = NULL, combine = TRUE,
                     prefix = NULL, device = "png", width = 10, height = 3, ...) {
  if (!is.gap(x, return.logical = TRUE)) {
    stop("x is no valid object of class 'gap'.")
  }
  type <- attr(x, "type")
  if (type == "HP") {
    x$tsl <- x
    if (contribution == TRUE) {
      stop("Output gap ad potential contributions not applicable for HP-filtered gaps.")
    }
  }
  if (!is.null(path)) {
    check <- dir.exists(paths = path)
    if (!check) {
      warning(paste0("The given file path '", path, "' does not exist. The plots will not be saved."))
      path <- NULL
    }
  }

  # confidence bounds do not exist
  tslBounds <- list(
    ub = NA,
    lb = NA
  )
  boundName <- ""
  tslBounds2 <- list(
    ub = NA,
    lb = NA
  )


  # --- data
  # potential and gdp growth
  tsl1 <- do.call(cbind, list(
    potential = 100 * growth(x$tsl$potential),
    gdp = 100 * growth(window(x$tsl$gdp, start = start(x$tsl$potential))),
    lb = tslBounds$lb,
    ub = tslBounds$ub
  ))

  # gap
  tsl2 <- do.call(cbind, list(
    gap = x$tsl$gap,
    lb = tslBounds2$lb,
    ub = tslBounds2$ub
  ))

  if (contribution) {
    lfTrend <- x$tsl$lTrend / x$tsl$ahoursTrend
    lf <- x$tsl$l / x$tsl$ahours
    alpha <- attr(x, "alpha")
    
    tsl1 <- na.trim(do.call(cbind, list(
      "average hours worked" = 100 * diff(x$tsl$ahoursTrend^alpha) / stats::lag(x$tsl$ahoursTrend, -1)^alpha,
      "working population" = 100 * diff(lfTrend^alpha) / stats::lag(lfTrend, -1)^alpha,
      "total factor productivity" = 100 * diff(x$tsl$tfpTrend) / stats::lag(x$tsl$tfpTrend, -1),
      "capital stock" = 100 * diff(x$tsl$k^(1 - alpha)) / stats::lag(x$tsl$k, -1)^(1 - alpha)
    )))
    tsl2 <- na.trim(do.call(cbind, list(
      "average hours worked" = 100 * alpha * (log(x$tsl$ahours) - log(x$tsl$ahoursTrend)),
      "working population" = 100 * alpha * (log(lf) - log(lfTrend)),
      "total factor productivity" = 100 * (log(x$tsl$tfp) - log(x$tsl$tfpTrend))
    )))
  }

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
  namesPrint <- paste(prefix, c("potential_growth", "gap"), sep = "_")

  if (contribution) {
    title <- list(
      "Potential GDP growth in %",
      "Output gap in %"
    )
  }

  # plot
  plotGap(
    tsl = tsl, legend = legend, title = title, boundName = boundName,
    contribution = contribution, res = NULL, namesPrint = namesPrint,
    bounds = FALSE, combine = combine, path = path, device = device,
    width = width, height = height
  )
}
