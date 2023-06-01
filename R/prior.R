
# -------------------------------------------------------------------------------------------

#' Initialization of prior distributions
#'
#' @description Initializes the prior distributions for a model of class \code{TFPmodel} or
#'   \code{NAWRUmodel}.
#'
#' @param model An object of class \code{TFPmodel} or \code{NAWRUmodel}.
#' @param MLE (Optional) A logical indicating whether the MLE estimates should be used for
#'   the initialization. The default is \code{MLE = FALSE} if \code{MLEfit} is not provided
#'   and vice versa.
#' @param MLEfit (Optional) An object of class \code{TFPfit} or \code{NAWRUfit} which is
#'   used if \code{MLE = TRUE}.
#'
#' @return A list of three matrices with parameters for the prior distribution and box
#'   constraints. Each list item refers to an equation, namely the \code{cycle}, \code{trend},
#'   and second observation equation. Each list element is a \code{4 x n} matrix where \code{n}
#'   denotes the number of parameters involved in the respective equation. The upper two
#'   elements specify the distribution, the lower two parameters specify box constraints.
#'   \code{NA} denotes no constraints. Autoregressive parameters are automatically restricted
#'   to the stationary region unless box constraints are specified. The respective prior
#'   distributions are defined through their mean and standard deviation. For instance,
#'   \code{prior$cycle[, 1]} contains the mean, standard deviation, lower and upper bound for
#'   the first variable, in that respective order.
#' @export
#'
#' @importFrom utils capture.output
initializePrior <- function(model, MLE = !is.null(MLEfit), MLEfit = NULL) {

  # model class and attributes
  class <- class(model)
  trend <- attr(model, "trend")
  cycle <- attr(model, "cycle")
  attrib <- unlist(attributes(model), recursive = FALSE)
  errorARMA <- attrib[grepl("errorARMA", names(attrib))][[1]]

  # initialize list
  prior <- list()

  # load model data
  namesExtract <- c("varName", "mean", "std", "LB", "UB")

  prior <- .accessDfSystem(model = model)
  prior <- lapply(prior, function(x) {
    x[, namesExtract]
  })

  # rearrange data
  prior <- lapply(prior, function(x) {
    rownames(x) <- x[, "varName"]
    x <- x[, c("mean", "std", "LB", "UB")]
    x <- t(x)
    return(x)
  })

  # drop trend variance
  if (trend != "RW1") {
    prior$trend <- prior$trend[, !grepl("tSigma", colnames(prior$trend)), drop = FALSE]
  }

  # obtain variance restrictions
  varRestr <- .initializeVar(model = model, type = "hp", prior = TRUE, errorARMA = errorARMA)
  trendNames <- colnames(prior$trend)[grepl("Sigma", colnames(prior$trend))]
  # adjust variance priors
  prior$cycle[c("mean", "std"), "cSigma"] <- mean(varRestr[, "cSigma"])
  prior$trend[c("mean", "std"), trendNames] <- mean(varRestr[, 2:(dim(varRestr)[2] - 1)])
  prior[[3]][c("mean", "std"), grepl("Sigma", colnames(prior[[3]]))] <- mean(varRestr[, "E2Sigma"])

  # ----- priors from MLE
  if (MLE) {
    if (inherits(model, "TFPmodel")) {
      if (!is.null(MLEfit) & !inherits(MLEfit, "TFPfit")) {
        MLEfit <- NULL
        message("The supplied object is not of class `TFPfit', starting MLE.")
      }
      if (is.null(MLEfit)) {
        parRestr <- initializeRestr(model = model, type = "hp")
        invisible(utils::capture.output(f <- suppressWarnings(fit(model = model, parRestr = parRestr))))
      } else {
        f <- MLEfit
      }
      pars <- f$parameters[sort(rownames(f$parameters)), 1]
      names(pars) <- sort(rownames(f$parameters))
      name <- rownames(f$parameters)
      Tn <- length(f$tsl$tfpTrend)

      # in case of NaNs for standard error; set standard error equal to estimate
      f$parameters[is.na(f$parameters[, 2]), 2] <- f$parameters[is.na(f$parameters[, 2]), 1]

      prior2 <- list()
      prior2$cubs <- t(as.matrix(f$parameters[grepl("cu", name, fixed = TRUE), 1:2]))
      prior2$trend <- t(as.matrix(f$parameters[grepl("td", name, fixed = TRUE), 1:2]))
      prior2$cycle <- t(as.matrix(f$parameters[grepl("c", name, fixed = TRUE) & !grepl("cu", name, fixed = TRUE), 1:2]))

      prior$cubs[1:2, colnames(prior2$cubs)] <- matrix(c(prior2$cubs[1, ], prior2$cubs[2, ] * sqrt(Tn)), 2, dim(prior2$cubs)[2], byrow = TRUE)
      prior$trend[1:2, colnames(prior2$trend)] <- matrix(c(prior2$trend[1, ], prior2$trend[2, ] * sqrt(Tn)), 2, dim(prior2$trend)[2], byrow = TRUE)
      prior$cycle[1:2, colnames(prior2$cycle)] <- matrix(c(prior2$cycle[1, ], prior2$cycle[2, ] * sqrt(Tn)), 2, dim(prior2$cycle)[2], byrow = TRUE)

      prior$cubs[2, grepl("Sigma", colnames(prior$cubs))] <- prior$cubs[1, grepl("Sigma", colnames(prior$cubs))]
      prior$trend[2, grepl("Sigma", colnames(prior$trend))] <- prior$trend[1, grepl("Sigma", colnames(prior$trend))]
      prior$cycle[2, grepl("Sigma", colnames(prior$cycle))] <- prior$cycle[1, grepl("Sigma", colnames(prior$cycle))]
    } else if (inherits(model, "NAWRUmodel")) {
      if (!is.null(MLEfit) & !inherits(MLEfit, "NAWRUfit")) {
        MLEfit <- NULL
        message("The supplied object is not of class `NAWRUfit', starting MLE.")
      }
      if (is.null(MLEfit)) {
        parRestr <- initializeRestr(model = model, type = "hp")
        invisible(utils::capture.output(f <- suppressWarnings(.MLEfitNAWRU(model = model, parRestr = parRestr))))
      } else {
        f <- MLEfit
      }
      pars <- f$parameters[sort(rownames(f$parameters)), 1]
      names(pars) <- sort(rownames(f$parameters))
      name <- rownames(f$parameters)
      Tn <- length(f$tsl$nawru)

      # in case of NaNs for standard error; set standard error equal to estimate
      f$parameters[is.na(f$parameters[, 2]), 2] <- f$parameters[is.na(f$parameters[, 2]), 1]

      prior2 <- list()
      prior2$pcInd <- t(as.matrix(f$parameters[grepl("pc", name, fixed = TRUE), 1:2]))
      prior2$trend <- t(as.matrix(f$parameters[grepl("td", name, fixed = TRUE), 1:2]))
      prior2$cycle <- t(as.matrix(f$parameters[grepl("c", name, fixed = TRUE) & !grepl("pc", name, fixed = TRUE), 1:2]))

      prior$pcInd[1:2, colnames(prior2$pcInd)] <- matrix(c(prior2$pcInd[1, ], prior2$pcInd[2, ] * sqrt(Tn)), 2, dim(prior2$pcInd)[2], byrow = TRUE)
      prior$trend[1:2, colnames(prior2$trend)] <- matrix(c(prior2$trend[1, ], prior2$trend[2, ] * sqrt(Tn)), 2, dim(prior2$trend)[2], byrow = TRUE)
      prior$cycle[1:2, colnames(prior2$cycle)] <- matrix(c(prior2$cycle[1, ], prior2$cycle[2, ] * sqrt(Tn)), 2, dim(prior2$cycle)[2], byrow = TRUE)

      prior$pcInd[2, grepl("Sigma", colnames(prior$pcInd))] <- prior$pcInd[1, grepl("Sigma", colnames(prior$pcInd))]
      prior$trend[2, grepl("Sigma", colnames(prior$trend))] <- prior$trend[1, grepl("Sigma", colnames(prior$trend))]
      prior$cycle[2, grepl("Sigma", colnames(prior$cycle))] <- prior$cycle[1, grepl("Sigma", colnames(prior$cycle))]
    }
  }

  # ensure unimodal distribution for beta for RAR2 cycle (alpha, beta (shapes) > 1)
  if (cycle == "RAR2") {
    m1 <- (prior$cycle[1, 1] - prior$cycle[3, 1]) / (prior$cycle[4, 1] - prior$cycle[3, 1])
    m2 <- prior$cycle[2, 1]
    cond1 <- sqrt(m1^2 * (1 - m1) / (1 + m1)) * (prior$cycle[4, 1] - prior$cycle[3, 1])
    cond2 <- sqrt((m1^3 - 2 * m1^2 + m1) / (2 + m1)) * (prior$cycle[4, 1] - prior$cycle[3, 1])
    prior$cycle[2, 1] <- min(cond1, cond2, m2)
    m1 <- (prior$cycle[1, 2] - prior$cycle[3, 2]) / (prior$cycle[4, 2] - prior$cycle[3, 2])
    m2 <- prior$cycle[2, 2]
    cond1 <- sqrt(m1^2 * (1 - m1) / (1 + m1)) * (prior$cycle[4, 2] - prior$cycle[3, 2])
    cond2 <- sqrt((m1^3 - 2 * m1^2 + m1) / (2 + m1)) * (prior$cycle[4, 2] - prior$cycle[3, 2])
    prior$cycle[2, 2] <- min(cond1, cond2, m2)
  }

  # ensure that standard deviation is not too small (invertible)
  sd_min <- 1e-6 # such that var_min = 1e-12
  prior[[1]][2, ] <- sapply(prior[[1]][2, ], max, sd_min)
  prior[[2]][2, ] <- sapply(prior[[2]][2, ], max, sd_min)
  prior[[3]][2, ] <- sapply(prior[[3]][2, ], max, sd_min)

  # return
  return(prior)
}

# -------------------------------------------------------------------------------------------

#' Converts the mean and standard deviation of a Gamma-distributed variable into the
#' parameters \eqn{s} and \eqn{\nu}.
#'
#' @param m mean parameter
#' @param std standard deviation
#'
#' @details The parameters \code{s} and \code{nu} are related to the regular shape
#'   \eqn{\alpha} and rate \eqn{\beta} parametrization in the following way:
#'   \eqn{\alpha = \nu / 2}
#'   \eqn{\beta = s / 2}
#'
#' @return A vector with parameters \eqn{s} and \eqn{\nu}.
#' @keywords internal
.meanStd2GammasNu <- function(m, std) {
  nu <- 2 * m^2 / std^2 + 4
  s <- m * (nu - 2)
  result <- c(s, nu)
  return(result)
}

# -------------------------------------------------------------------------------------------

#' Transforms the prior distribution defined by mean and standard deviation to the
#' appropriate input parameters.
#'
#' @param prior A matrix containing the mean and standard deviation of the priors.
#' @param restr A matrix containing box constraints for the parameters.
#' @param namesInvGammaDistr A character vector specifying the names of all inverse Gamma
#'   distributed parameters.
#' @param namesNormalDistr A character vector specifying the names of all normally
#'   distributed parameters.
#' @param namesBetaDistr A character vector specifying the names of all Beta distributed
#'   parameters.
#' @keywords internal
.priorMSd2Parameter <- function(prior, restr, namesInvGammaDistr, namesNormalDistr, namesBetaDistr) {
  # transform mean and standard deviation to appropriate input parameters
  for (names in namesInvGammaDistr) {
    prior[, names] <- .meanStd2GammasNu(m = prior[1, names], std = prior[2, names])
  }
  for (names in namesBetaDistr) {
    prior[, names] <- .meanStd2Beta(
      m = prior[1, names], std = prior[2, names],
      lb = restr[1, names], ub = restr[2, names]
    )
  }
  prior[2, namesNormalDistr] <- prior[2, namesNormalDistr]^2
  # deal with no box constraints
  distrPar <- rbind(prior, restr)
  distrPar[3, is.na(distrPar[3, ])] <- -Inf
  distrPar[4, is.na(distrPar[4, ])] <- Inf
  # return
  return(distrPar)
}

# -------------------------------------------------------------------------------------------

#' Converts the mean and standard deviation of a (possibly scaled) Beta-distributed variable
#' into the two shape parameters \eqn{\alpha} and \eqn{\beta}.
#'
#' @param m mean parameter
#' @param std standard deviation
#' @param lb lower bound
#' @param ub upper bound
#'
#' @return A vector with shape parameters \eqn{\alpha} and \eqn{\beta}.
#' @keywords internal
.meanStd2Beta <- function(m, std, lb = 0, ub = 1) {
  # unscaled mean and standard deviation
  m <- (m - lb) / (ub - lb)
  std <- std / (ub - lb)
  # check bounds
  if (m < 0 | m > 1) {
    m <- (m < 0) * (0 + 1e-8) + (m > 1) * (1 - 1e-8)
  }
  if (std^2 > m * (1 - m)) {
    std <- sqrt(m * (1 - m) - 1e-8)
  }
  alpha <- m^2 / std^2 * (1 - m) - m
  beta <- alpha / m - alpha
  result <- c(alpha, beta)
  return(result)
}
