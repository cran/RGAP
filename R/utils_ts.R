# --------------------------------------------------------------------------------------------------- #

#' HP filter
#'
#' @description Applies the Hodrick Prescott Filter.
#'
#' @param x A univariate time series object.
#' @param lambda The smoothing parameter.
#'
#' @return A univariate time series object containing the trend of the original time series.
#'
#' @export
#' @importFrom stats start end window ts lag frequency time
#' @examples
#' # get data for France
#' data("gap")
#' country <- "France"
#' tsList <- amecoData2input(gap[[country]], alpha = 0.65)
#' hp <- hpfilter(x = tsList$gdp, lambda = 6.25)
hpfilter <- function(x, lambda) {
  n <- length(x[is.na(x) == FALSE])
  A <- 6 * diag(n)
  A[row(A) == (col(A) - 1)] <- -4
  A[row(A) == (col(A) + 1)] <- -4
  A[row(A) == (col(A) - 2)] <- 1
  A[row(A) == (col(A) + 2)] <- 1
  A[1:2, 1:2] <- matrix(c(1, -2, -2, 5), 2, 2)
  A[(n - 1):n, (n - 1):n] <- matrix(c(5, -2, -2, 1), 2, 2)

  trend <- ts(NA, start = start(x), end = end(x), frequency = frequency(x))
  trend[is.na(x) == FALSE] <- (solve(diag(n) + lambda * A)) %*% x[is.na(x) == FALSE]

  trend
}

# ---------------------------------------------------------------------------------------------------

#' Growth rate
#'
#' @description Computes growth rates for time series.
#'
#' @param ts A time series object.
#' @param k An integer specifying the lag size in the growth computation. The default is
#'   \code{k = 1}.
#'
#' @importFrom stats lag is.mts
#' @keywords internal
growth <- function(ts, k = 1) {
  g <- ts / lag(ts, -k) - 1
  if (is.mts(ts)) colnames(g) <- colnames(ts)
  g
}


# ---------------------------------------------------------------------------------------------------

#' Performs a mathematical operation to the ts elements of two lists with the same names
#' @param x A list.
#' @param y A list.
#' @param operator  mathematical operator in quotation marks, e.g., "*".
#' @keywords internal
operTsLists <- function(x, y, operator) {
  namesList <- names(x)
  resList <- list()
  for (ii in (1:length(namesList))) {
    resList[[namesList[ii]]] <- eval(parse(text = paste("x[[namesList[ii]]]", operator, "y[[namesList[ii]]]", sep = " ")))
  }
  resList
}

# ---------------------------------------------------------------------------------------------------

#' Finds first/last starting/end date in list of time series deepening on the input functions.
#'
#' @param x A list.
#' @param FUN1 A function, namely either \code{start} or \code{end} from \code{stats}.
#' @param FUN2 A function, namely either \code{max} or \code{min} from \code{base}.
#'
#' @importFrom stats start end window ts lag frequency time
#' @importFrom zoo na.trim
#' @keywords internal
dateTsList <- function(x, FUN1 = start, FUN2 = max) {
  # delete NAs
  x <- lapply(x, na.trim)
  # locate ts with appropriate year
  nameY <- names(which(unlist(lapply(lapply(x, FUN1), `[[`, 1)) == FUN2(unlist(lapply(lapply(x, FUN1), `[[`, 1)))))
  # get appropriate quarter in that year
  nameYQ <- names(which(unlist(lapply(lapply(x[nameY], FUN1), `[[`, 2)) == FUN2(unlist(lapply(lapply(x[nameY], FUN1), `[[`, 2)))))[1]
  # get appropriate date
  date <- FUN1(x[[nameYQ]])
  date
}
