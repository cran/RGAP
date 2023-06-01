# -------------------------------------------------------------------------------------------

#' CUBS indicator
#'
#' @description Computes the capacity utilization economic sentiment (CUBS) indicator.
#'
#' @param tsCU A multiple time series containing three survey time series, the first
#'   element needs to be capacity utilization in industry, see details. Alternatively, a list
#'   of time series can be supplied.
#' @param tsVA A multiple time series containing three value added series that correspond to
#'   \code{tsCU}. Alternatively, a list of time series can be supplied.
#' @param frequency The frequency of the computed cubs indicator. Possible entries are
#'   \code{frequency = 1} (annual), \code{frequency = 4} (quarterly). The default is
#'   \code{frequency = 1}.
#' @param lambda The smoothing parameter for the application of the HP filter (see details).
#'   If not supplied, \code{lambda = 6.25} is used for yearly data and \code{lambda = 1600} 
#'   for quarterly data.
#'
#' @details The list \code{tslCU} contains capacity utilization in industry, and the relevant
#'   survey outcomes of the construction and service sector. The first list object needs to
#'   contain capacity utilization in industry.
#' @details The list \code{tslVA} contains the real value added series for the industry,
#'   construction and service sector in the same order as \code{tslCU}.
#' @details The computed CUBS indicator consists exclusively of capacity utilization in
#'   industry until both other series become available.
#'
#' @return A list containing the two time series capacity utilization in industry \code{cu}
#'  and the CUBS indicator \code{cubs}.
#'
#' @export
#' @importFrom stats start end window ts lag frequency time var
#' @importFrom zoo na.trim
#' @examples
#' # load data for Germany
#' data("gap")
#' country <- "Germany"
#'
#' # compute cubs indicator
#' namesCubs <- c("indu", "serv", "buil")
#' namesVACubs <- paste0("va", namesCubs)
#' tscubs <- cubs(
#'   tsCU = gap[[country]][, namesCubs],
#'   tsVA = gap[[country]][, namesVACubs]
#' )
cubs <- function(tsCU, tsVA, frequency = 1, lambda = NULL) {
  tslCU <- as.list(tsCU)
  tslVA <- as.list(tsVA)

  # adjust frequency
  if (is.null(lambda)) {
    lambda <- 1600 / ((4 / frequency)^4)
  }

  # check input variables
  list2env(.checkCubs(tsCU = tslCU, tsVA = tslVA, lambda = lambda, frequency = frequency),
           envir = environment())


  # adjust names such that they are not the same
  namesCU <- names(tslCU)
  namesVA <- names(tslVA)
  names(tslCU) <- paste0("cu", namesVA)

  # frequency adjustment
  tslCU <- suppressWarnings(lapply(tslCU, .cubsTa, conversion = "average", frequency = frequency))
  # warns when extending series
  tslVA <- lapply(tslVA, .cubsTa, conversion = "sum", frequency = frequency)
  tslCU <- lapply(tslCU, zoo::na.trim)
  tslVA <- lapply(tslVA, zoo::na.trim)

  # get last start date and first end date
  start_final <- dateTsList(x = c(tslCU[-1], tslVA[-1]), FUN1 = start, FUN2 = max)
  end_final <- dateTsList(x = c(tslCU[-1], tslVA[-1]), FUN1 = end, FUN2 = min)

  # adjust names
  names(tslCU) <- namesVA

  # weights and scaling factors
  tslVA <- lapply(tslVA, log)
  tslVAtrend <- list()
  tslVAtrend[namesVA] <- lapply(tslVA, hpfilter, lambda)
  tslVAtrend$sumExp <- Reduce("+", lapply(tslVAtrend, exp))

  tslVAcycle <- list()
  tslVAcycle[namesVA] <- operTsLists(tslVA[namesVA], tslVAtrend[namesVA], operator = "-")

  weights <- list()
  weights[namesVA] <- lapply(tslVAtrend[namesVA], function(x) {
    exp(x) / tslVAtrend$sumExp
  })

  scale <- list()
  scale[namesVA] <- lapply(tslVAcycle, function(x) {
    ts(sqrt(var(x, na.rm = TRUE)) / sqrt(var(tslVAcycle[[namesVA[1]]], na.rm = TRUE)),
      start = start(x), end = end(x), frequency = frequency(x)
    )
  })

  # STEP 1: normalize
  tslNorm <- lapply(tslCU, window, start = start_final, end = end_final)
  mCU <- mean(tslNorm[[1]], na.rm = TRUE)
  stdCU <- sqrt(var(tslNorm[[1]], na.rm = TRUE))
  tslNorm <- lapply(tslNorm, .normalize)

  # STEP 2: rescale
  tslScale <- operTsLists(tslNorm, scale, operator = "*")

  # STEP 3: weighted average
  tslwAve <- operTsLists(tslScale, weights, operator = "*")
  cubsComb <- Reduce("+", tslwAve)

  # STEP 4: normalize 2nd time
  cubsComb <- .normalize(cubsComb)

  # STEP 5: rescale 2nd time
  cubsComb <- (cubsComb * stdCU + mCU)

  # merge CU and CUBS
  tsl <- list()
  tsl$cu <- tsl$cubs <- tslCU[[1]]
  tsl$cubs <- window(tslCU[[1]], end = end(cubsComb))
  window(tsl$cubs, start = start(cubsComb), end = end(cubsComb)) <- cubsComb

  return(tsl)
}

# -------------------------------------------------------------------------------------------

#' Normalizes a time series / vector.
#'
#' @param x E.g. a time series.
#' @keywords internal
.normalize <- function(x) {
  m <- mean(x, na.rm = TRUE)
  std <- sqrt(var(x, na.rm = TRUE))
  y <- (x - m) / std
  return(y)
}

# -------------------------------------------------------------------------------------------

#' Adjusts the frequency of cubs input series.
#'
#' @param tsObj A time series object.
#' @param conversion An appropriate conversion method, i.e., \code{"average"} or \code{"sum"}.
#' @param frequency The frequency of the aggregated time series, i.e., \code{frequency = 4} for
#'     quarterly and \code{frequency = 1} for annual.
#'
#' @details If \code{frequency == 1}, then the values for the latest and the first year are
#'    averaged over the existing months/quarters. However, for the latest year, this is only
#'    done if the third quarter is available.
#'
#' @importFrom stats start end window ts lag frequency time
#' @keywords internal
.cubsTa <- function(tsObj, conversion, frequency) { # to annual
  freq <- frequency(tsObj)

  if (freq == frequency) {
    tsRes <- tsObj
  } else {
    if (frequency == 1) { # conversion to yearly frequency
      startYear <- start(tsObj)[1]
      endYear <- end(tsObj)[1]
      startValue <- mean(tsObj[floor(time(tsObj)) == startYear], na.rm = TRUE)
      endValue <- mean(tsObj[floor(time(tsObj)) == endYear], na.rm = TRUE)
      # temporal aggregation
      tsRes <- .aggregate(x = tsObj, freqLow = frequency, FUN = ifelse(conversion == "average", mean, sum))

      ## procedure according to EC:
      # add first and last year
      suppressWarnings(window(tsRes, start = startYear, end = startYear, extend = TRUE) <- ifelse(conversion == "average", startValue, freq * startValue))
      suppressWarnings(window(tsRes, start = endYear, end = endYear, extend = TRUE) <- ifelse(conversion == "average", endValue, freq * endValue))
      # delete last year if the date is before the third quarter
      if ((end(tsObj)[2] < 3 & freq == 4) | (end(tsObj)[2] < 9 & freq == 12)) {
        suppressWarnings(window(tsRes, start = endYear, end = endYear) <- NA)
      }
    } else { # conversion to quarterly frequency
      # temporal aggregation
      tsRes <- .aggregate(x = tsObj, freqLow = frequency, FUN = ifelse(conversion == "average", mean, sum))
    }
  }
  return(tsRes)
}

# -------------------------------------------------------------------------------------------

#' Temporally Aggregates time series.
#'
#' @param x A time series object.
#' @param freqLow Frequency of low frequency series.
#' @param FUN A function for aggregation.
#'
#' @importFrom stats start end window ts frequency time aggregate
#' @importFrom zoo na.trim
#' @keywords internal
.aggregate <- function(x, freqLow, FUN) {

  # frequency, start and end
  freq <- frequency(x)
  start <- start(x)
  end <- end(x)

  # remove NAs
  if (!all(is.na(x))) x <- zoo::na.trim(x)

  # years
  year <- floor(time(x))

  # fill first and last year/quarter to aggregate over correct period
  x <- suppressWarnings(window(x, start = c(year[1], 1), end = c(year[length(year)], freq), extend = TRUE))
  # warns due to extension

  # get year and quarter
  yearquarter <- floor(4 * time(x) + 0.001) / 4
  year <- floor(time(x))

  if (freqLow == 1) {
    sort <- year
  } else if (freqLow == 4) {
    sort <- yearquarter
  }

  # aggregate
  tmp <- aggregate(as.numeric(x), FUN = FUN, by = list(sort))

  # time series object
  xLow <- ts(tmp$x, start = tmp$Group.1[1], frequency = freqLow)
  if (!all(is.na(xLow))) xLow <- zoo::na.trim(xLow)
  xLow
}
