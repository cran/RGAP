# -------------------------------------------------------------------------------------------

#' Initialization of parameter restrictions
#'
#' @description Initializes parameter restrictions for objects of class \code{NAWRUmodel},
#' \code{TFPmodel}, or \code{KuttnerModel}.
#'
#' @param model An object of class \code{NAWRUmodel}, \code{TFPmodel}, or \code{KuttnerModel}.
#' @param type The variance restriction type. Possible options are \code{"basic"},
#'   \code{"hp"}, see details. The default is \code{type = "basic"}.
#' @param lambda The smoothing constant for the HP-filter if \code{type = "hp"}.
#' @param q Quantile for the Inverse Gamma distribution (only used if \code{type = "hp"}). The 
#'   default is \code{q = 0.01}.
#'
#' @details For \code{type = "hp"}, the HP filter is applied to the appropriately differences 
#'   first observation series to obtain its trend and cycle. Subsequently, the specified trend 
#'   and cycle models are fitted to obtain its innovation variance. Moreover, the second 
#'   observation series (according to its specification) is fitted to obtain its innovation 
#'   variance. Lastly, the obtained innovations variances are used to get lower and upper 
#'   bounds. To that end, the \code{q} and \code{1-q} quantiles of the inverse gamma 
#'   distribution are used, with mean and standard deviation set to the estimated variances.
#'
#' @return A list of three matrices containing the parameter restrictions for the cycle,
#'   trend, and the second observation equation. Each matrix contains the lower and upper
#'   bound of the involved parameters. \code{NA} implies that no
#'   restriction is present.
#'
#' @export
initializeRestr <- function(model, type = "basic", lambda = NULL, q = 0.01) {

  # model attributes
  trend <- attr(model, "trend")
  cycle <- attr(model, "cycle")
  attrib <- unlist(attributes(model), recursive = FALSE)
  errorARMA <- attrib[grepl("errorARMA", names(attrib))][[1]]

  # initialize
  restr <- list()

  # load model data
  namesExtract <- c("varName", "LB", "UB")
  restr <- .accessDfSystem(model = model)
  restr <- lapply(restr, function(x) {
    x[, namesExtract]
  })

  # rearrange data
  restr <- lapply(restr, function(x) {
    rownames(x) <- x[, "varName"]
    x <- x[, c("LB", "UB")]
    x <- t(x)
    return(x)
  })

  # variance restrictions
  varRestr <- .initializeVar(model = model, type = type, lambda = lambda, errorARMA = errorARMA, q = q)
  trendNames <- model$loc$varName[grepl("trend", model$loc$variableRow) & model$loc$sysMatrix == "Q"]
  restr$cycle[c("UB", "LB"), "cSigma"] <- varRestr[, "cSigma"]
  restr$trend[c("UB", "LB"), trendNames] <- varRestr[, 2:(dim(varRestr)[2] - 1)]
  restr[[3]][c("LB"), grepl("Sigma", colnames(restr[[3]]))] <- 1e-8
  restr[[3]][c("UB", "LB"), grepl("Sigma", colnames(restr[[3]]))][1:2] <- varRestr[, "E2Sigma"]

  # if tSigma is 0
  if (trend != "RW1" & dim(model$loc[model$loc$sysMatrix == "Q", ])[2] < 4) {
    restr$trend <- restr$trend[, !grepl("tSigma", colnames(restr$trend)), drop = FALSE]
  }

  # returns
  restr
}
# -------------------------------------------------------------------------------------------

#' Initializes variance restrictions.
#'
#' @param errorARMA The ARMA order of the second equation error process.
#' @param prior A logical indicating whether prior parameters should be returned.
#' @inheritParams initializeRestr
#' @inheritParams fit.NAWRUmodel
#' @inheritParams hpfilter
#'
#' @importFrom stats optim arima as.formula lm
#' @keywords internal
.initializeVar <- function(model, type = NULL, lambda = NULL, prior = FALSE, errorARMA = c(0, 0), q = 0.05) {

  # model attributes
  trend <- attr(model, "trend")
  cycle <- attr(model, "cycle")

  # set lambda if not supplied
  if (is.null(lambda)) {
    freq <- frequency(model$tsl[[1]])
    lambda <- 1600 / ((4 / freq)^4)
  }

  # 1) compute variance of second observation equation
  # 2) compute variance of HP-filtered trend of first observation equation
  if (inherits(model, "NAWRUmodel")) {
    tsE2 <- model$tsl$pcInd
    varE2 <- var(model$tsl$pcInd, na.rm = TRUE)
    tsE1 <- model$tsl$ur
  } else if (inherits(model, "TFPmodel")) {
    tsE2 <- model$tsl$cubs
    varE2 <- var(model$tsl$cubs, na.rm = TRUE)
    tsE1 <- model$tsl$logtfp
  } else if (inherits(model, "KuttnerModel")) {
    tsE2 <- stats::lag(model$tsl$dinfl, -2) # such that second equation is correct later ?????
    tsE2 <- stats::lag(model$tsl$dinfl, 0) #
    varE2 <- var(model$tsl$dinfl, na.rm = TRUE)
    tsE1 <- model$tsl$loggdp
  }

  diffOrder <- switch(trend,
    "DT" = 1,
    "RW2" = 2,
    "RW1" = 1
  ) # for stationarity of trend

  if (type == "basic") {
    eps <- 1e-8
    varRestrTrend <- varRestrCycle <- c(var(diff(tsE1, differences = diffOrder), na.rm = TRUE), eps)
    varRestrE2 <- c(varE2, eps)
  } else {
    tsE1Hp <- hpfilter(tsE1, lambda = lambda)

    # trend
    arOrder <- switch(trend,
      "DT" = 1,
      "RW2" = 0,
      "RW1" = 0
    )
    tsTrendStat <- diff(tsE1Hp, differences = diffOrder)
    varTrend <- tryCatch(
      {
        arima(x = tsTrendStat, order = c(arOrder, 0, 0))$sigma2
      },
      error = function(cont) {
        warning("The HP-filtered trend process does not fit the model specification in terms of stationarity,
                the proposed trend variance restrictions might be compromised. Consider to respecify.")
        return(arima(x = diff(tsTrendStat), order = c(arOrder, 0, 0))$sigma2)
      }
    )

    # cycle
    arOrder <- switch(cycle,
      "AR2" = 2,
      "RAR2" = 2,
      "AR1" = 1
    )
    tsCycle <- tsE1 - tsE1Hp
    varCycle <- arima(x = tsCycle, order = c(arOrder, 0, 0))$sigma2

    # 2nd equation
    nExo <- length(model$tsl) - 2
    cycleLag <- unlist(attributes(model), recursive = F)
    cycleLag <- unlist(cycleLag[grepl("cycleLag", names(cycleLag))])

    datal <- list(E2 = tsE2)
    for (k in cycleLag) {
      datal_tmp <- list(stats::lag(tsCycle, -k))
      names(datal_tmp) <- paste0("cycleLag", k)
      datal <- c(datal, datal_tmp)
    }
    if (nExo != 0) {
      datal <- c(datal, model$tsl[3:(2 + nExo)])
    }
    datal <- na.trim(do.call(cbind, datal))
    varE2 <- tryCatch(
      {
        arima(x = datal[, 1], order = c(errorARMA[1], 0, errorARMA[2]), xreg = datal[, -1])$sigma2
      },
      error = function(cont) {
        stop("The selected exogenous variables likely are linearly dependent. Consider to respecify.")
      }
    )
    
    
    # computes quantiles of gamma distribution with mean and standard deviation
    # set to the previously computed variances
    beta <- 1
    varRestrE2 <- rev(.intervalIGamma(varE2, beta * sqrt(varE2), qlower = q))
    varRestrTrend <- rev(.intervalIGamma(varTrend, beta * sqrt(varTrend), qlower = q))
    varRestrCycle <- rev(.intervalIGamma(varCycle, beta * sqrt(varCycle), qlower = q))
  }

  # for prior, use estimate for mean and standard deviation
  if (prior) {
    varRestrCycle <- rep(varCycle, 2)
    varRestrTrend <- rep(varTrend, 2)
    varRestrE2 <- rep(varE2, 2)
  }

  varRestr <- matrix(c(
    varRestrCycle,
    c(0, 0),
    varRestrTrend,
    varRestrE2
  ),
  2, 4,
  byrow = FALSE
  )
  colnames(varRestr) <- c("cSigma", "tSigma", "tdSigma", "E2Sigma")

  if (trend == "RW1" | dim(model$loc[model$loc$sysMatrix == "Q", ])[2] < 4) {
    varRestr <- varRestr[, -2]
  }
  # return
  varRestr
}

# -------------------------------------------------------------------------------------------

#' Updates the parameter constraints for on object of class \code{NAWRUmodel} or
#' \code{TFPmodel}.
#'
#' @param model A model of class \code{NAWRUmodel} or \code{TFPmodel}.
#' @inheritParams fit.NAWRUmodel
#'
#' @return The same model with updated list item \code{loc}.
#' @keywords internal
.updateParConstraints <- function(model, parRestr) {
  equ <- names(parRestr)
  sysMat <- unique(model$loc$sysMatrix)
  for (j in equ) {
    varName <- colnames(parRestr[[j]])
    for (k in varName) {
      ind <- which(model$loc$varName == k)
      lb <- parRestr[[j]][1, k]
      ub <- parRestr[[j]][2, k]
      if (!is.na(lb) & !is.na(ub)) {
        tmp <- paste0("I", lb, "_", ub)
        model$loc$LB[ind] <- lb
        model$loc$UB[ind] <- ub
      } else if (is.na(lb) & !is.na(ub)) {
        tmp <- paste0("I", "-Inf", "_", ub)
        model$loc$LB[ind] <- NA
        model$loc$UB[ind] <- ub
      } else if (!is.na(lb) & is.na(ub)) {
        tmp <- paste0("I", lb, "_", "Inf")
        model$loc$LB[ind] <- lb
        model$loc$UB[ind] <- NA
      } else {
        tmp <- "NA"
        model$loc$LB[ind] <- NA
        model$loc$UB[ind] <- NA
      }
      for (l in sysMat) {
        varname2 <- model$loc$varName[model$loc$sysMatrix == l]
        for (m in varname2) {
          if (m == k) {
            if ((tmp != "NA") | !grepl("AR", model$loc$restriction[model$loc$varName == m])) {
              model$loc$restriction[model$loc$varName == m] <- tmp
            }
          }
        }
      }
    }
  }
  model
}

# -------------------------------------------------------------------------------------------

#' Transforms the parameters of an AR(2) process to its re-parametrized version RAR(2) and
#' vice versa.
#'
#' @param par vector of parameters
#' @param phi2Atau logical indicating whether the transformation should be from AR(2) to
#'   RAR(2) or vice versa.
#'
#' @return A vector with the transformed parameters
#' @keywords internal
.RAR2transform <- function(par, phi2Atau = TRUE) {
  if (phi2Atau) {
    phi1 <- par[1]
    phi2 <- par[2]
    A <- sqrt(-phi2)
    tau <- 2 * pi / acos(phi1 / (2 * sqrt(-phi2)))
    result <- c(A, tau)
  } else {
    A <- par[1]
    tau <- par[2]
    phi1 <- 2 * A * cos(2 * pi / tau)
    phi2 <- -(A^2)
    result <- c(phi1, phi2)
  }
  names(result) <- names(par)
  result
}

# -------------------------------------------------------------------------------------------

#' Applies suitable contraining functions to parameters.
#'
#' @param par names vector with parameters
#' @param loc A data frame containing information on each involved parameter (list element
#'   of objects of class \code{NAWRUmodel}, \code{TFPmodel}, \code{KuttnerModel}).
#'
#' @keywords internal
assignConstraints <- function(par, loc) {
  namePar <- names(par)
  nPar <- length(par) - sum(loc$restriction == "AR2") / 2 - sum(loc$restriction == "RAR2") / 2
  parConstr <- NULL

  count <- 0
  for (ii in 1:nPar) {
    count <- count + 1
    constr <- loc$restriction[count]
    parTmp <- par[count]
    if (constr == "AR2" || constr == "RAR2") {
      parTmp <- par[count:(1 + count)]
      count <- count + 1
    }
    parConstrTmp <- constraint(parTmp, type = constr)
    parConstr <- c(parConstr, parConstrTmp)
  }

  names(parConstr) <- namePar
  parConstr
}

# -------------------------------------------------------------------------------------------

#' Applies contraints to parameters.
#'
#' @param type character with constraint type
#' @inheritParams assignConstraints
#'
#' @keywords internal
constraint <- function(par, type = NA) {
  parConstr <- switch(type,
    "NA"  = par,
    "AR2" = constrAR2(par),
    "AR1" = constrAR1(par),
    "neg" = constrInterval(par, -1000, -1e-8),
    "pos" = constrInterval(par, 1e-8, 1000)
  )

  if (substr(type, 1, 1) == "I") {
    a <- as.numeric(strsplit(substr(type, 2, 100), "_")[[1]][1])
    b <- as.numeric(strsplit(substr(type, 2, 100), "_")[[1]][2])
    if (a > -Inf & b < Inf) {
      parConstr <- constrInterval(par, a, b)
    } else if (a == -Inf) {
      parConstr <- constrInterval(par, -1000, b)
    } else if (b == Inf) {
      parConstr <- constrInterval(par, a, 1000)
    }
  }

  parConstr
}

# -------------------------------------------------------------------------------------------

#' Computes the covariance of the estimated parameters given restrictions.
#'
#' @param CV unconstrained variance covariance matrix
#' @inheritParams assignConstraints
#'
#' @keywords internal
computeCovar <- function(par, loc, CV) {
  namePar <- names(par)
  nPar <- length(par) - sum(loc$restriction == "AR2") / 2 - sum(loc$restriction == "RAR2") / 2

  count <- 0
  D <- matrix(0, length(par), length(par))
  for (ii in 1:nPar) {
    count <- count + 1
    parTmp <- par[count]
    constr <- loc$restriction[loc$varName == names(parTmp)]
    if (constr == "AR2" || constr == "RAR2") {
      parTmp <- par[count:(1 + count)]
    }
    DTmp <- Dconstraint(parTmp, type = constr)
    if (constr == "AR2" || constr == "RAR2") {
      D[count:(count + 1), count:(count + 1)] <- DTmp
      count <- count + 1
    } else {
      D[count, count] <- DTmp
    }
  }

  rownames(D) <- namePar
  colnames(D) <- namePar

  res <- D %*% CV %*% t(D)
  res
}

# -------------------------------------------------------------------------------------------

#' Extracts the derivative of the applied restriction function.
#'
#' @inheritParams constraint
#'
#' @keywords internal
Dconstraint <- function(par, type = NA) {
  D <- switch(type,
    "NA"  = 1,
    "AR2" = DconstrAR2(par),
    "AR1" = DconstrAR1(par),
    "neg" = DconstrInterval(par, -1000, -1e-8),
    "pos" = DconstrInterval(par, 1e-8, 1000)
  )

  if (substr(type, 1, 1) == "I") {
    a <- as.numeric(strsplit(substr(type, 2, 100), "_")[[1]][1])
    b <- as.numeric(strsplit(substr(type, 2, 100), "_")[[1]][2])
    D <- DconstrInterval(par, a, b)
  }

  D
}


# -------------------------------------------------------------------------------------------

constrInterval <- function(par, a, b) {
  a <- max(a, -1000)
  b <- min(b, 1000)
  parConstr <- a + (b - a) * (1 / (1 + exp(-par)))
  parConstr
}

# -------------------------------------------------------------------------------------------

DconstrInterval <- function(par, a, b) {
  a <- max(a, -1000)
  b <- min(b, 1000)
  logitInv <- (1 / (1 + exp(-par)))
  D <- (b - a) * logitInv * (1 - logitInv)
  D
}

# -------------------------------------------------------------------------------------------

constrAR1 <- function(par) {
  eps <- 0.01
  a <- -1 + eps
  b <- 1 - eps
  parConstr <- a + (b - a) * (1 / (1 + exp(-par)))
  parConstr
}

# -------------------------------------------------------------------------------------------

DconstrAR1 <- function(par) {
  eps <- 0.01
  a <- -1 + eps
  b <- 1 - eps
  logitInv <- (1 / (1 + exp(-par)))
  D <- (b - a) * logitInv * (1 - logitInv)
  D
}

# -------------------------------------------------------------------------------------------

constrAR2 <- function(par) {
  eps <- 0.01
  psi1 <- par[1]
  psi2 <- par[2]
  z1 <- psi1 / (1 + (1 + eps) * abs(psi1))
  z2 <- psi2 / (1 + (1 + eps) * abs(psi2))
  phi <- c(0, 0)
  phi[1] <- z1 + z2
  phi[2] <- -1 * z1 * z2
  phi
}

# -------------------------------------------------------------------------------------------

DconstrAR2 <- function(par) {
  eps <- 0.01
  psi1 <- par[1]
  psi2 <- par[2]
  a <- ((-abs(psi1)) / (1 + eps + abs(psi1))^2 + 1 / (1 + eps + abs(psi1)))
  b <- ((-abs(psi2)) / (1 + eps + abs(psi2))^2 + 1 / (1 + eps + abs(psi2)))
  c <- -psi2 / (1 + eps + abs(psi2)) + a
  d <- -psi1 / (1 + eps + abs(psi1)) + b
  D <- matrix(c(a, b, c, d), 2, 2, byrow = TRUE)
  D
}
