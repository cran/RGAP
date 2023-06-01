
# -------------------------------------------------------------------------------------------

#' Checks whether model and prior match.
#'
#' @inheritParams .BayesFitTFP
#' @keywords internal
.checkModelPrior <- function(model, prior) {
  namesModel <- model$loc$varName
  namesPrior <- colnames(do.call(cbind, prior))
  check <- all(namesModel %in% namesPrior) &&
    all(namesPrior %in% namesModel)
  if (!check) {
    stop("Specified model and prior do not match. Check model parameters.")
  }
}


# -------------------------------------------------------------------------------------------

#' Checks whether model, prior and MLE fit match.
#'
#' @inheritParams .BayesFitTFP
#' @keywords internal
.checkModelMLEfit <- function(model, MLEfit) {
  namesModel <- model$loc$varName
  namesFit <- MLEfit$model$loc$varName
  namesModel <- namesModel[namesModel != "tSigma"]
  namesFit <- namesFit[namesFit != "tSigma"]
  check <- all(namesModel %in% namesFit) &&
    all(namesFit %in% namesModel)
  if (!check) {
    stop("Specified model and MLEfit do not match. Check model parameters.")
  }
}

# -------------------------------------------------------------------------------------------

#' Checks the input variables for the procedure \code{cubs} for consistency and validity.
#'
#' @inheritParams cubs
#' @keywords internal
#' @keywords internal
.checkCubs <- function(tsCU, tsVA, lambda, frequency) {
  n <- length(tsCU)
  m <- length(tsVA)
  minFrequency <- min(sapply(c(tsCU, tsVA), frequency))
  maxFrequency <- max(sapply(c(tsCU, tsVA), frequency))
  frequencyPossibilities <- c(1, 4)


  if (any(unlist(lapply(lapply(tsCU, is.na), all))) | any(unlist(lapply(lapply(tsVA, is.na), all)))) {
    index <- (unlist(lapply(lapply(tsCU, is.na), all)) | unlist(lapply(lapply(tsVA, is.na), all)))
    if (index[1]) stop("Capacity utilization or value added in industry are not provided.")
    tsCU[index] <- NULL
    tsVA[index] <- NULL
  }
  if (n != m) {
    stop("The time series lists are not of the same length", call. = FALSE)
  }
  if (n != 3) {
    stop(paste0("The time series list 'tsCU' is of length ", n, " instead of 3."), call. = FALSE)
  }
  if (m != 3) {
    stop(paste0("The time series list 'tsVA' is of length ", m, " instead of 3."), call. = FALSE)
  }
  if (!(frequency %in% frequencyPossibilities)) {
    stop(paste0("The frequency is misspecified, please use one of the following: ", paste(frequencyPossibilities, collapse = ", ")), call. = FALSE)
  }
  if (maxFrequency < frequency) {
    stop("The frequency is larger than the provided time series, please respecify.", call. = FALSE)
  }
  if (lambda < 0) {
    stop("The smoothing constant is negative, please respecify.", call. = FALSE)
  }

  res <- list(
    tsCU = tsCU,
    tsVA = tsVA,
    lambda = lambda,
    frequency = frequency
  )
  return(res)
}

# -------------------------------------------------------------------------------------------

#' Checks the input parameters of \code{.BayesFitTFP} and \code{.BayesFitNAWRU} for
#' consistency.
#'
#' @param type A character specifying whether a "nawru" or "tfp" model should be checked.
#' @inheritParams .BayesFitTFP
#' @keywords internal
.checkBayesInput <- function(model, type, prior = NULL, R = NULL, burnin = NULL,
                             thin = NULL, HPDIprob = NULL, FUN = NULL, MLEfit = NULL) {
  if (type == "tfp") {
    if (!is.TFPmodel(model, return.logical = TRUE)) {
      stop("Model object is not of class 'TFPmodel'.")
    }
  } else if (type == "nawru") {
    if (!is.NAWRUmodel(model, return.logical = TRUE)) {
      stop("Model object is not of class 'NAWRUmodel'.")
    }
  }
  if (!is.null(R) && (!is.numeric(R) || R <= 0 || R %% 1)) {
    stop("R negative or not an integer. Please respecify.")
  }
  if (!is.null(burnin) && (!is.numeric(burnin) || burnin <= 0 || burnin %% 1)) {
    stop("burnin negative or not an integer. Please respecify.")
  }
  if (!is.null(thin) && (!is.numeric(thin) || thin <= 0 || thin %% 1)) {
    stop("Thinning parameter negative or not an integer. Please respecify.")
  }
  if (!is.null(HPDIprob) && (!is.numeric(HPDIprob) || HPDIprob <= 0 || HPDIprob >= 1)) {
    stop("Level for highest posterior denisty credible interval invalid. Please respecify a value between 0 and 1.")
  }
  if (!is.null(FUN) && (!is.function(FUN) && !identical("mean", FUN) & !identical("median", FUN))) {
    stop("Invalid specification for FUN. Possible options are 'mean' and 'median'.")
  }
}

# -------------------------------------------------------------------------------------------

#' Checks the given prior information for consistency and applicability.
#'
#' @param model An object of class \code{TFPmodel}.
#' @param prior A list of matrices with parameters for the prior distribution and box
#'   constraints.
#' @keywords internal
.checkPrior <- function(model, prior) {
  priorControl <- initializePrior(model)
  priorControl <- unlist(lapply(priorControl, colnames))
  if (any(!(priorControl %in% unlist(lapply(prior, colnames)))) | any(!(unlist(lapply(prior, colnames)) %in% priorControl))) {
    stop("Prior specification does not match the model. Please use the function 'initializePrior' to initialize your prior information.")
  }

  # beta distribution
  # ensure unimodal distribution for beta for RAR2 cycle (alpha, beta (shapes) > 1)
  cycle <- attr(model, "cycle")
  if (cycle == "RAR2") {
    ab1 <- .meanStd2Beta(
      m = prior$cycle[1, 1], std = prior$cycle[2, 1],
      lb = prior$cycle[3, 1], ub = prior$cycle[4, 1]
    )
    ab2 <- .meanStd2Beta(
      m = prior$cycle[1, 2], std = prior$cycle[2, 2],
      lb = prior$cycle[3, 2], ub = prior$cycle[4, 2]
    )
    if (ab1[1] < 1 | ab1[2] < 1 | ab2[1] < 1 | ab2[2] < 1) {
      warning("Prior distribution (beta-distribution) for A and tau is not unimodal.")
    }
  }

  # standard deviation large enough such that variance is invertible
  mat <- diag(c(prior[[1]][2, ], prior[[2]][2, ], prior[[3]][2, ])^2)
  invertible <- class(try(solve(mat), silent = T)) == "matrix"
  if (!invertible[1]) {
    stop("Prior standard deviation is too small (not invertible). Please respecify.")
  }
}

# -------------------------------------------------------------------------------------------

#' Checks the given variance restrictions for consistency.
#'
#' @param model An object of class \code{NAWRUmodel} or \code{TFPmodel}.
#' @param parRestr list of matrices containing the parameter restrictions for the cycle,
#'   trend, and the second observation equation (Phillips curve, CUBS equation). Each matrix
#'   contains the lower and upper bound of the involved parameters. \code{NA} implies that no
#'   restriction is present.
#' @keywords internal
.checkParRestr <- function(model, parRestr) {
  parRestrControl <- initializeRestr(model = model)

  # the list items of parRestr should be a matrix of dimension k x n
  k <- 2

  if (length(parRestr) != 3) {
    stop(paste0(
      "Misspecified parameter restriction.\n",
      "Please initialize parameter restrictions using the function 'initializeRestr(model)'.\n"
    ))
  }

  for (j in 1:3) {

    # dimensions
    k <- dim(parRestr[[j]])[1]
    n <- dim(parRestr[[j]])[2]
    kControl <- dim(parRestrControl[[j]])[1]
    nControl <- dim(parRestrControl[[j]])[2]

    # plausibility
    parRestr[[j]][1, is.na(parRestr[[j]][1, ])] <- -Inf
    parRestr[[j]][2, is.na(parRestr[[j]][2, ])] <- Inf
    index_eq <- (parRestr[[j]][2, ] - parRestr[[j]][1, ] == 0)
    index <- (parRestr[[j]][2, ] - parRestr[[j]][1, ] < 0)
    parNames <- colnames(parRestr[[j]])


    if (k != kControl | n != nControl) {
      stop(paste0(
        "Misspecified parameter restriction.\n",
        "Please initialize parameter restrictions using the function 'initializeRestr(model)'.\n"
      ))
    }
    if (any(index_eq) && parNames[index_eq] != model$loc$varName[model$loc$equation == "trend" & model$loc$variableRow == "trend" & model$loc$sysMatrix == "Q"]) {
      stop(paste0(
        "Misspecified parameter restriction.\n",
        "The upper bound of variable ", parNames[index], " is equal to the lower bound.\n",
        "Please respecify.\n"
      ))
    }
    if (any(index)) {
      stop(paste0(
        "Misspecified parameter restriction.\n",
        "The upper bound of variable ", parNames[index], " is smaller than the lower bound.\n",
        "Please respecify.\n"
      ))
    }
  }
}

# -------------------------------------------------------------------------------------------

#' Checks the input variables for the procedure \code{NAWRUmodel} for consistency and
#' validity.
#'
#' @param errorARMA A vector with non-negative integers specifying the AR
#'   and MA degree of the error term in the Phillip's curve equation.
#' @param exoNames A character vector containing the names of the exogenous variables.
#' @inheritParams NAWRUmodel
#' @keywords internal
.checkNawru <- function(tsl, trend, cycle, type, cycleLag, errorARMA, exoNames, exoType, start, end, anchor, anchor.h) {
  trendPossibilities <- c("RW1", "RW2", "DT")
  cyclePossibilities <- c("AR1", "AR2")
  typePossibilities <- c("TKP", "NKP")
  cycleLagPossibilities <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  errorARPossibilities <- c(0, 1, 2)
  errorMAPossibilities <- c(0)

  if (!is.list(tsl)) {
    tsl <- as.list(tsl)
  }
  if (!(trend %in% trendPossibilities)) {
    stop(paste0("The trend is misspecified, please use one of the following: ", paste(trendPossibilities, collapse = ", ")), call. = FALSE)
  }
  if (!(cycle %in% cyclePossibilities)) {
    stop(paste0("The cycle is misspecified, please use one of the following: ", paste(cyclePossibilities, collapse = ", ")), call. = FALSE)
  }
  if (!(type %in% typePossibilities)) {
    stop(paste0("The Phillips curve type is misspecified, please use one of the following: ", paste(typePossibilities, collapse = ", ")))
  }
  if (!is.null(exoNames) && !is.null(exoType)) {
    dimExo <- dim(exoType)
    dimExoNames <- length(exoNames)
    if (dimExo[2] != dimExoNames) {
      stop(paste("Array exoType or vector exoNames misspecified.",
        "\n exoType is a n x m x 2 array, where m is the number of exogenous variables for which n specifications can be included.",
        "\n exoNames is an m dimensional vector with corresponding variable names",
        sep = " "
      ))
    }
    if (length(unique(exoNames)) != dimExoNames) {
      stop("Variable names for exogenous variables are not unique. Please assign unique names.")
    }
    if (!all(exoNames %in% names(tsl))) {
      stop("Specified exogneous variables are not part of tsl. Please include all variables in tsl or adjust exoNames and exoType.")
    }
    if (type == "NKP") {
      warning("Exogenous variables are not compatible with the Ney Keynesian Phillip's curve. They will be dropped.",
        call. = FALSE
      )
      exoNames <- NULL
      exoType <- NULL
    }
  }
  if (type == "NKP") {
    if (cycle == "AR1" && any(cycleLag != 0)) {
      warning("Invalid cycleLag for the New Keynesian Phillip's curve and an AR(1) cycle. cycleLag is set to 0.",
        call. = FALSE
      )
      cycleLag <- 0
    }
    if (cycle == "AR2" && any(cycleLag != c(0, 1))) {
      warning(paste("Invalid cycleLag for the New Keynesian Phillip's curve and an AR(2) cycle. cycleLag is set to 0:1.",
        "The parameter on the lagged cycle is not estimated but instead implied by other model parameters.",
        sep = " "
      ),
      call. = FALSE
      )
      cycleLag <- 0:1
    }
  }
  if (!is.null(anchor)) {
    if (is.null(anchor.h)) {
      anchor.h <- frequency(tsl[[1]]) * 10
      message(paste("No anchor horizon specified. Horizon is set to", anchor.h, "periods.", sep = " "))
    }
  }
  if (!all(cycleLag %in% cycleLagPossibilities)) {
    stop(paste0(
      "The number of cycle terms in the Phillip's curve equation is misspecified, please use a combination or one of the following: ",
      paste(cycleLagPossibilities, collapse = ", ")
    ), call. = FALSE)
  }
  if (!(errorARMA[1] %in% errorARPossibilities)) {
    stop(paste0(
      "The AR part of the error term in the Phillip's curve equation is misspecified, please use one of the following: ",
      paste(errorARPossibilities, collapse = ", ")
    ), call. = FALSE)
  }
  if (!(errorARMA[2] %in% errorMAPossibilities)) {
    stop(paste0(
      "The MA part of the error term in the Phillip's curve equation is misspecified, please use one of the following: ",
      paste(errorMAPossibilities, collapse = ", ")
    ), call. = FALSE)
  }

  res <- list(
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
  )
  return(res)
}

# -------------------------------------------------------------------------------------------

#' Checks the input variables for the procedure \code{TFPmodel} for consistency and validity.
#'
#' @param errorARMA A vector with non-negative integers specifying the AR
#'   and MA degree of the error term in the CUBS equation.
#' @inheritParams TFPmodel
#' @keywords internal
.checkTfp <- function(tsl, trend, cycle, cycleLag, cubsAR, errorARMA, start, end, anchor, anchor.h) {
  trendPossibilities <- c("RW1", "RW2", "DT")
  cyclePossibilities <- c("AR1", "AR2", "RAR2")
  cycleLagPossibilities <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  cubsARPossibilities <- c(0, 1, 2)
  errorARPossibilities <- c(0, 1, 2)
  errorMAPossibilities <- 0

  if (!is.list(tsl)) {
    tsl <- as.list(tsl)
  }
  if (!(trend %in% trendPossibilities)) {
    stop(paste0("The trend is misspecified, please use one of the following: ", paste(trendPossibilities, collapse = ", ")), call. = FALSE)
  }
  if (!(cycle %in% cyclePossibilities)) {
    stop(paste0("The cycle is misspecified, please use one of the following: ", paste(cyclePossibilities, collapse = ", ")), call. = FALSE)
  }
  if (!all(cycleLag %in% cycleLagPossibilities)) {
    stop(paste0(
      "The number of cycle terms in the cubs equation is misspecified, please use a combination or one of the following: ",
      paste(cycleLagPossibilities, collapse = ", ")
    ), call. = FALSE)
  }
  if (!(cubsAR %in% cubsARPossibilities)) {
    stop(paste0(
      "The number of lagged cubs terms in the cubs equation is misspecified, please use one of the following: ",
      paste(cubsARPossibilities, collapse = ", ")
    ), call. = FALSE)
  }
  if (!(errorARMA[1] %in% errorARPossibilities)) {
    stop(paste0(
      "The AR part of the error term in the cubs equation is misspecified, please use one of the following: ",
      paste(errorARPossibilities, collapse = ", ")
    ), call. = FALSE)
  }
  if (!(errorARMA[2] %in% errorMAPossibilities)) {
    stop(paste0(
      "The MA part of the error term in the cubs equation is misspecified, please use one of the following: ",
      paste(errorMAPossibilities, collapse = ", ")
    ), call. = FALSE)
  }
  if (!is.null(anchor)) {
    if (is.null(anchor.h)) {
      anchor.h <- frequency(tsl[[1]]) * 10
      message(paste("No anchor horizon specified. Horizon is set to", anchor.h, "periods.", sep = " "))
    }
  }

  res <- list(
    tsl = tsl,
    trend = trend,
    cycle = cycle,
    cycleLag = cycleLag,
    cubsAR = cubsAR,
    start = start,
    end = end,
    anchor = anchor,
    anchor.h = anchor.h
  )
  return(res)
}

# -------------------------------------------------------------------------------------------

#' Checks the input variables for the procedure \code{KuttnerModel} for consistency and validity.
#'
#' @param errorARMA A vector with non-negative integers specifying the AR
#'   and MA degree of the error term in the inflation equation.
#' @inheritParams KuttnerModel
#' @keywords internal
.checkKuttner <- function(tsl, trend, cycle, cycleLag, errorARMA, start, end, anchor, anchor.h) {
  trendPossibilities <- c("RW1", "RW2", "DT")
  cyclePossibilities <- c("AR1", "AR2", "RAR2")
  cycleLagPossibilities <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  errorARPossibilities <- 0
  errorMAPossibilities <- 1:3

  if (!is.list(tsl)) {
    tsl <- as.list(tsl)
  }
  if (!(trend %in% trendPossibilities)) {
    stop(paste0(
      "The trend is misspecified, please use one of the following: ",
      paste(trendPossibilities, collapse = ", ")
    ), call. = FALSE)
  }
  if (!(cycle %in% cyclePossibilities)) {
    stop(paste0(
      "The cycle is misspecified, please use one of the following: ",
      paste(cyclePossibilities, collapse = ", ")
    ), call. = FALSE)
  }
  if (!all(cycleLag %in% cycleLagPossibilities)) {
    stop(paste0(
      "The number of cycle terms in the inflation equation is misspecified, please use a combination or one of the following: ",
      paste(cycleLagPossibilities, collapse = ", ")
    ), call. = FALSE)
  }
  if (!(errorARMA[1] %in% errorARPossibilities)) {
    stop(paste0(
      "The AR part of the error term in the inflation equation is misspecified, please use one of the following: ",
      paste(errorARPossibilities, collapse = ", ")
    ), call. = FALSE)
  }
  if (!(errorARMA[2] %in% errorMAPossibilities)) {
    stop(paste0(
      "The MA part of the error term in the inflation equation is misspecified, please use one of the following: ",
      paste(errorMAPossibilities, collapse = ", ")
    ), call. = FALSE)
  }
  if (!is.null(anchor)) {
    if (is.null(anchor.h)) {
      anchor.h <- frequency(tsl[[1]]) * 10
      message(paste("No anchor horizon specified. Horizon is set to", anchor.h, "periods.", sep = " "))
    }
  }

  res <- list(
    tsl = tsl,
    trend = trend,
    cycle = cycle,
    errorARMA = errorARMA,
    start = start,
    end = end,
    anchor = anchor,
    anchor.h = anchor.h
  )
  return(res)
}
