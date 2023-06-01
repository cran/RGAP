#' Fit Method
#'
#' @description Defines the fit method.
#'
#' @param model Some model.
#' @param ... Some stuff passed on to methods.
#' 
#' @return Depends on the model object, see documentation of specific methods.
#' 
#' @family fitting methods
#'
#' @export
fit <- function(model, ...) UseMethod("fit")

# -------------------------------------------------------------------------------------------

#' Prepares state space model system matrices to create an object of type \code{NAWRUmodel}
#' or \code{TFPmodel}.
#'
#' @param tsl A list of input time series computed during the procedure \code{NAWRUmodel}
#'   or \code{TFPmodel}.
#' @param errorARMA A vector with non-negative integers specifying the AR and MA degree of
#'   the error term in the second observation equation.
#' @inheritParams NAWRUmodel
#'
#' @importFrom stats start end window ts lag frequency time
#' @keywords internal
.SSSystem <- function(tsl, cycle, trend, cycleLag, type = NULL, errorARMA) {

  # ----- T, R, Q

  # cycle
  cycleLagMax <- max(cycleLag)
  T_cycle <- rbind(
    c(rep(NA, 2), rep(0, max(0, cycleLagMax - 1))),
    cbind(diag(1, max(1, cycleLagMax)), rep(0, max(1, cycleLagMax)))
  )
  R_cycle <- matrix(c(1, rep(0, max(1, cycleLagMax))), max(1, cycleLagMax) + 1, 1)
  names_cycle <- c("cycle", paste0("cycleLag", 1:(max(1, cycleLagMax))))
  if (cycle == "AR1") {
    T_cycle[1, 2] <- 0
  }

  # trend
  T_trend <- matrix(c(1, 1, 0, 1), 2, 2, byrow = TRUE)
  R_trend <- diag(1, 2)
  names_trend <- c("trend", "trendDrift")
  if (trend == "RW1") {
    T_trend[2, 2] <- 0
    R_trend <- matrix(c(1, 0), 2, 1)
  } else if (trend == "DT") {
    T_trend <- matrix(c(1, 1, 0, NA), 2, 2, byrow = TRUE)
    # constant needs to be added later
  }

  # 2nd equation error
  p <- errorARMA[1]
  q <- errorARMA[2]
  T_error11 <- rbind(
    rep((p > 0) * NA, max(1, p)),
    cbind(diag(1, max(0, p - 1)), rep(0, max(0, p - 1)))
  )
  if (p == 0) {
    T_error11[1, 1] <- 0
  }
  T_error12 <- rbind(
    rep(NA, q),
    matrix(0, max(0, p - 1), q)
  )
  T_error21 <- matrix(0, q, max(1, p))
  T_error22 <- rbind(
    rep(0, q),
    cbind(diag(1, max(0, q - 1)), rep(0, max(0, q - 1)))
  )
  if (q > 0) {
    T_error <- rbind(
      cbind(T_error11, T_error12),
      cbind(T_error21, T_error22)
    )
  } else {
    T_error <- T_error11
  }

  R_error1 <- matrix(c(1, rep(0, max(0, p - 1))), max(1, p), 1)
  R_error2 <- matrix(c(rep(1, min(1, q)), rep(0, max(0, q - 1))), q, 1)
  R_error <- rbind(R_error1, R_error2)
  names_error1 <- paste0("E2errorL", 0:p)[1:(1 + p - 1)]
  names_error2 <- paste0("E2innoL", 0:q)[1:(1 + q - 1)]
  names_error <- c(names_error1, names_error2[rep(q > 0, length(names_error2))])

  # constant
  T_constant <- matrix(1, 1, 1)
  names_const <- "const"

  # merge
  Tt <- rbind(
    cbind(T_cycle, matrix(0, nrow(T_cycle), ncol(T_trend) + ncol(T_error) + ncol(T_constant))),
    cbind(matrix(0, nrow(T_trend), ncol(T_cycle)), T_trend, matrix(0, nrow(T_trend), ncol(T_error) + ncol(T_constant))),
    cbind(matrix(0, nrow(T_error), ncol(T_cycle) + ncol(T_trend)), T_error, matrix(0, nrow(T_error), ncol(T_constant))),
    cbind(matrix(0, nrow(T_constant), ncol(T_cycle) + ncol(T_trend) + ncol(T_error)), T_constant)
  )
  Rt <- rbind(
    cbind(R_cycle, matrix(0, nrow(R_cycle), ncol(R_trend) + ncol(R_error))),
    cbind(matrix(0, nrow(R_trend), ncol(R_cycle)), R_trend, matrix(0, nrow(R_trend), ncol(R_error))),
    cbind(matrix(0, nrow(R_error), ncol(R_cycle) + ncol(R_trend)), R_error),
    cbind(matrix(0, 1, ncol(R_cycle) + ncol(R_trend) + ncol(R_error)))
  )
  Qt <- diag(1, ncol(Rt))
  diag(Qt) <- NA

  # name matrices
  obsNames <- names(tsl[1:2])
  stateNames <- c(names_cycle, names_trend, names_error, names_const)
  varNames <- stateNames[apply(Rt, 1, sum) > 0][1:ncol(Rt)]
  colnames(Tt) <- rownames(Tt) <- rownames(Rt) <- stateNames
  colnames(Rt) <- colnames(Qt) <- rownames(Qt) <- varNames

  # add constants
  if (trend != "RW2") {
    Tt[names_trend[2], names_const] <- NA
  }

  # ----- Z, H

  # initialize
  Zt <- matrix(0, 2, ncol(Tt))
  Ht <- diag(0, nrow(Zt))

  # name matrices
  rownames(Zt) <- colnames(Ht) <- rownames(Ht) <- obsNames
  colnames(Zt) <- stateNames

  # assign Z
  Zt[1, c(names_cycle[1], names_trend[1])] <- 1
  Zt[2, c(names_cycle[cycleLag + 1], names_const)] <- NA
  Zt[2, names_error[1]] <- 1

  # ----- P1, P1inf, a1
  # P1:     variance of stationary part should be assigned
  # P1inf:  diffuse parts should be set to 1
  # a1:     diffuse parts should be set to 0

  # initialize
  a1 <- matrix(0, nrow(Tt), 1)
  P1 <- P1inf <- diag(0, nrow(Tt))

  # name matrices
  colnames(P1) <- rownames(P1) <- colnames(P1inf) <- rownames(P1inf) <- rownames(a1) <- stateNames

  # assign
  a1[names_const, ] <- 1
  diag(P1[c(names_cycle, names_error[1]), c(names_cycle, names_error[1])]) <- NA
  diag(P1inf[names_trend, names_trend]) <- 1

  # ----- model dimensions
  nPar <- sum(is.na(Tt)) + sum(is.na(Zt)) + sum(is.na(Qt))
  nObs <- nrow(Zt)
  nState <- nrow(Tt)
  nStateV <- nrow(Qt)

  # return
  sys <- list(Zt = Zt, Tt = Tt, Qt = Qt, Rt = Rt, Ht = Ht, a1 = a1, P1 = P1, P1inf = P1inf)
  res <- list(
    tsl = tsl,
    nPar = nPar,
    nObs = nObs,
    nState = nState,
    nStateV = nStateV,
    stateNames = stateNames,
    sys = sys
  )
  return(res)
}

# -------------------------------------------------------------------------------------------

#' Updates the system matrices of an object of class \code{NAWRUmodel} or \code{TFPmodel}
#' during optimization or during a Bayesian Gibbs procedure.
#'
#' @param pars A vector of parameters.
#' @param SSModel An object of class \code{SSModel} specifying the state-space model.
#' @param loc A data frame containing information on each involved parameter, for instance
#'   its corresponding system matrix, variable names, and parameter restrictions.
#' @param errorARMA A vector with non-negative integers specifying the AR
#'   and MA degree of the error term in the second observation equation.
#' @param bayes A logical indicating whether the update is part of a Bayesian procedure, i.e.,
#'   the parameter constraints do not need to be enforced.
#' @inheritParams NAWRUmodel
#' @inheritParams fit.NAWRUmodel
#'
#' @importFrom stats start end window ts lag frequency time
#' @keywords internal
.updateSSSystem <- function(pars, SSModel, loc, cycle, trend, errorARMA, signalToNoise = NULL, type = NULL, bayes = FALSE) {

  # assign location information
  locT <- loc[loc$sysMatrix == "T", ]
  locZ <- loc[loc$sysMatrix == "Z", ]
  locQ <- loc[loc$sysMatrix == "Q", ]
  locTt <- loc[loc$sysMatrix == "Tt", ]
  locExo <- loc[loc$sysMatrix == "exo", ]

  # assign parameter constraints
  if (!bayes) {
    parConstr <- assignConstraints(pars, loc = loc)
  }
  else {
    parConstr <- pars
  }

  # assign stationary and non stationary variables
  indexStat <- c(
    colnames(SSModel$T)[grepl("cycle", colnames(SSModel$T))],
    colnames(SSModel$T)[grepl("error", colnames(SSModel$T))],
    colnames(SSModel$T)[grepl("inno", colnames(SSModel$T))]
  )
  indexRoot <- colnames(SSModel$T)[grepl("trend", colnames(SSModel$T))]

  # assign name of second observation equation
  nameE2 <- colnames(SSModel$y)[2]

  # ----- assign parameters

  # observation equation
  SSModel$Z[nameE2, locZ$variableRow, ] <- parConstr[locZ$varName]
  if (!is.null(type)) {
    if (type == "NKP" & cycle == "AR2") {
      SSModel$Z[nameE2, "cycleLag1", ] <- 0.99 * parConstr[locZ$varName[locZ$variableRow == "cycle"]] * parConstr[locT$varName[locT$variableRow == "cycleLag1"]]
    }
  }
  # variance
  if (!is.null(signalToNoise)) {
    k <- dim(SSModel["Q"])[1] - 1
    SSModel["Q"][-(2:k), -(2:k), ] <- diag(parConstr[locQ$varName])
    SSModel["Q"][k, k, ] <- parConstr["cSigma"] * signalToNoise
  } else {
    SSModel["Q"] <- diag(parConstr[locQ$varName])
    # SSModel$R[locQ$variableRow,,] <- diag(parConstr[locQ$varName]) ##### neu
  }
  # cycle
  if (cycle == "AR1") {
    SSModel$T["cycle", locT$variableRow, ] <- parConstr[locT$varName]
  } else if (cycle == "AR2") {
    SSModel$T["cycle", locT$variableRow[grepl("cycle", locT$variableRow)], ] <- parConstr[locT$varName[grepl("cycle", locT$variableRow)]]
  } else if (cycle == "RAR2") {
    A <- parConstr[locT$varName[locT$variableRow == "cycle"]]
    tau <- parConstr[locT$varName[locT$variableRow == "cycleLag1"]]
    SSModel$T["cycle", "cycle", ] <- 2 * A * cos(2 * pi / tau)
    SSModel$T["cycle", "cycleLag1", ] <- -(A^2)
  }
  # trend
  if (trend == "DT") {
    SSModel$T["trendDrift", "const", ] <- parConstr[locTt$varName[locTt$variableRow == "const"]] * (1 - parConstr[locTt$varName[locTt$variableRow == "trendDrift"]])
    SSModel$T["trendDrift", "trendDrift", ] <- parConstr[locTt$varName[locTt$variableRow == "trendDrift"]]
    indexRoot <- indexRoot[!(indexRoot %in% "trendDrift")]
    indexStat <- c(indexStat, "trendDrift")
  } else if (trend == "RW1") {
    SSModel$T["trendDrift", "const", ] <- parConstr[locT$varName[locT$variableRow == "trendDrift"]]
    indexRoot <- indexRoot[!(indexRoot %in% "trendDrift")]
    indexStat <- c(indexStat, "trendDrift")
  }
  # options
  if (errorARMA[1] > 0) {
    SSModel$T["E2errorL0", "E2errorL0", ] <- parConstr[locT$varName[locT$variableRow == "E2errorL0"]]
  }
  if (errorARMA[1] == 2) {
    SSModel$T["E2errorL0", "E2errorL1", ] <- parConstr[locT$varName[locT$variableRow == "E2errorL1"]]
  }
  if (errorARMA[2] > 0) {
    SSModel$T["E2errorL0", locT$variableRow[grepl("E2innoL", locT$variableRow)], ] <- parConstr[locT$varName[grepl("E2innoL", locT$variableRow)]]
  }

  # assign P1 and P1inf for diffuse inizialization
  nStat <- length(indexStat)
  nRoot <- length(indexRoot)
  tryCatch(
    {
      modR <- SSModel$R[, , 1] %*% chol(SSModel$Q[, , 1])
    },
    error = function(cont) {
      stop("The stationary part of the model is close to being non-stationary, please respecify.")
    }
  )
  modR <- modR[indexStat, ]
  modR <- modR[, colSums(modR) > 0]
  tryCatch(
    {
      SSModel$P1[indexStat, indexStat] <- matrix(solve(diag(nStat^2) - kronecker(SSModel$T[indexStat, indexStat, 1], SSModel$T[indexStat, indexStat, 1])) %*%
        c(modR %*% t(modR)), nStat, nStat)
    },
    error = function(cont) {
      stop("The stationary part of the model is close to being non-stationary, please respecify.")
    }
  )
  SSModel$P1inf[] <- 0
  SSModel$P1inf[indexRoot, indexRoot] <- diag(nRoot)

  # update unconditional mean
  SSModel$a1[indexStat, ] <- solve(diag(nStat) - SSModel$T[indexStat, indexStat, 1]) %*% SSModel$T[indexStat, "const", 1]
  SSModel$a1[locExo$variableRow, ] <- parConstr[locExo$varName]

  # return
  SSModel
}

# -------------------------------------------------------------------------------------------

#' Modifies a an object of type \code{NAWRUmodel} or \code{TFPmodel} in case the variance
#' constraint for the trend is set to zero or in case a signal-to-noise ratio is specified.
#'
#' @inheritParams fit.TFPmodel
#'
#' @importFrom KFAS SSModel SSMcustom
#' @keywords internal
.modifySSSystem <- function(model, signalToNoise) {

  # get trend
  trend <- attr(model, "trend")

  # obtain system variances
  Z <- model$SSModel$Z
  T <- model$SSModel$T
  R <- model$SSModel$R
  Q <- model$SSModel$Q
  a1 <- model$SSModel$a1
  P1 <- model$SSModel$P1
  P1inf <- model$SSModel$P1inf
  H <- model$SSModel$H
  stateNames <- colnames(model$SSModel$T)

  # get rid of trend variance
  if (trend != "RW1") {
    varNames <- model$loc$variableRow[model$loc$sysMatrix == "Q"]
    name_delete <- stateNames[grepl("trend", stateNames)][1]
    index_delete <- which(R[name_delete, , ] == 1)
    R <- R[, -index_delete, ]
    Q <- Q[-index_delete, -index_delete, ]
    model$loc <- model$loc[!(model$loc$sysMatrix == "Q" & model$loc$variableRow == "trend"), ]
  }

  # signal to noise ratio specified
  if (!is.null(signalToNoise)) {
    model$loc <- model$loc[!(model$loc$sysMatrix == "Q" & grepl("trend", model$loc$variableRow)), , drop = FALSE]
    if (trend != "RW1") {
      varNames <- model$loc$variableRow[model$loc$sysMatrix == "Q"]
      index <- varNames != "trend"
      Q <- Q[index, index]
      R <- R[, index]
    }
  }

  # ----- state space model
  if (inherits(model, "TFPmodel")) {
    modelSS <- SSModel(cbind(logtfp, cubs) ~ -1 + SSMcustom(Z = Z, T = T, R = R, Q = Q, a1 = a1, P1 = P1, P1inf = P1inf, state_names = stateNames),
      H = H,
      data = model$tsl
    )
  } else if (inherits(model, "NAWRUmodel")) {
    modelSS <- SSModel(cbind(ur, pcInd) ~ -1 + SSMcustom(Z = Z, T = T, R = R, Q = Q, a1 = a1, P1 = P1, P1inf = P1inf, state_names = stateNames),
      H = H,
      data = model$tsl
    )
  } else if (inherits(model, "KuttnerModel")) {
    modelSS <- SSModel(cbind(loggdp, dinfl) ~ -1 +
      +SSMcustom(Z = Z, T = T, R = R, Q = Q, a1 = a1, P1 = P1, P1inf = P1inf, state_names = stateNames),
    H = H,
    data = model$tsl
    )
  }

  # return
  model$SSModel <- modelSS
  model
}

# -------------------------------------------------------------------------------------------

#' Initialization of exogenous variables
#'
#' Initializes the transformations applied to exogenous variables.
#'
#' @param varNames A \code{(k x 1)} character vector containing the names of the exogenous
#'   variables.
#' @param D A \code{(n x k)} matrix containing the difference transformations, see details.
#' @param L A \code{(n x k)} matrix containing the lag transformations, see details.
#'
#' @details For the matrices \code{D} and \code{L}, the rows denote different transformations 
#'   to each of the variables in the columns. \code{NA} indicates no transformation.
#'
#' @return An array of size \code{(n, k, 2)}. The \code{[, , 1]} specifies
#'   the difference order and \code{[, , 2]} the lag order.
#'
#' @export
initializeExo <- function(varNames, D = NULL, L = NULL) {
  k <- length(varNames)
  # n <- max(maxLag, maxDiff) + 1
  n <- 2
  
  if (is.null(D) & is.null(L)) {
    D <- L <- matrix(NA, n, k)
  } else if (is.null(D) & !is.null(L)) {
    L <- as.matrix(L)
    D <- L
    D[] <- NA
    n <- dim(L)[1]
  } else if (is.null(D) & !is.null(L)) {
    D <- as.matrix(D)
    L <- D
    L[] <- NA
    n <- dim(D)[1]
  } else {
    L <- as.matrix(L)
    D <- as.matrix(D)
    n <- dim(D)[1]
    if (!all(dim(D) == dim(L))) {
      stop("If both D and L are supplied, they need to have the same dimension.")
    }
  }

  # check consistency of varNames
  if (k != dim(D)[2]) {
    stop("The dimension of D or L does not match the number or variables in varNames.")
  }
  
  # initialize
  exoType <- array(NA, dim = c(n, k, 2))
  dimnames(exoType)[[3]] <- c("difference", "lag")
  colnames(exoType) <- varNames
  
  exoType[, , 1] <- D
  exoType[, , 2] <- L

  exoType
}

# -------------------------------------------------------------------------------------------

#' Accesses the internal data frame dfSystem which contains data on the parameters to be
#' estimated.
#'
#' @param model A model of class \code{NAWRUmodel} or \code{TFPmodel}.
#'
#' @return A data frame containing information on each involved parameter, for instance
#'   its corresponding system matrix, variable names, and parameter restrictions.
#' @keywords internal
.accessDfSystem <- function(model) {

  # model attributes
  trend <- attr(model, "trend")
  cycle <- attr(model, "cycle")

  # initialize
  tmp <- list()
  type <- NULL

  # access cycle and trend data
  tmp$cycle <- dfSystem[dfSystem$equation == "cycle" & dfSystem$variant == cycle, ]
  tmp$trend <- dfSystem[dfSystem$equation == "trend" & dfSystem$variant == trend, ]

  if (inherits(model, "TFPmodel")) {
    nameE2 <- "cubs"
    # tfp attributes
    cycleLag <- attr(model, "cubs")$cycleLag
    cubsAR <- attr(model, "cubs")$cubsAR
    errorARMA <- attr(model, "cubs")$errorARMA
    errorAR <- errorARMA[1]
    errorMA <- errorARMA[2]
    tmp[[nameE2]] <- NULL
    if (cubsAR > 0) {
      tmp[[nameE2]] <- rbind(tmp[[nameE2]], dfSystem[dfSystem$equation == nameE2 & dfSystem$variant == paste0("cubsAR", cubsAR), ])
    }
  } else if (inherits(model, "NAWRUmodel")) {
    nameE2 <- "pcInd"
    # nawru attributes
    type <- attr(model, "phillips curve")$type
    cycleLag <- attr(model, "phillips curve")$cycleLag
    errorARMA <- attr(model, "phillips curve")$errorARMA
    errorAR <- errorARMA[1]
    errorMA <- errorARMA[2]
    exoNames <- attr(model, "phillips curve")$exoVariables

    tmp[[nameE2]] <- NULL
    if (length(exoNames) > 0) {
      for (name in exoNames) {
        if (grepl("pcInd", name)) {
          tmp[[nameE2]] <- rbind(tmp[[nameE2]], dfSystem[dfSystem$equation == "pcInd" & dfSystem$variant == "pcIndAR", ])
        } else {
          tmp[[nameE2]] <- rbind(tmp[[nameE2]], dfSystem[dfSystem$equation == "pcInd" & dfSystem$variant == "exo", ])
          tmp[[nameE2]]$varName[grepl("XXX", tmp[[nameE2]]$varName)] <- name
        }
      }
    }
  } else if (inherits(model, "KuttnerModel")) {
    nameE2 <- "infl"
    # kuttner attributes
    cycleLag <- attr(model, "inflation equation")$cycleLag
    errorARMA <- attr(model, "inflation equation")$errorARMA
    errorAR <- errorARMA[1]
    errorMA <- errorARMA[2]

    tmp[[nameE2]] <- NULL
  }

  tmp[[nameE2]] <- rbind(
    tmp[[nameE2]],
    dfSystem[dfSystem$equation == nameE2 & dfSystem$variant == "base", ],
    dfSystem[dfSystem$equation == "E2error" & dfSystem$variant == "base", ]
  )
  if (!(0 %in% cycleLag)) {
    tmp[[nameE2]] <- tmp[[nameE2]][!(tmp[[nameE2]]$sysMatrix == "Z" & tmp[[nameE2]]$variableRow == "cycle"), ]
  }
  if (max(cycleLag) > 0) {
    tmp[[nameE2]] <- rbind(tmp[[nameE2]], dfSystem[dfSystem$equation == nameE2 & dfSystem$variant == "cycleLag", ][cycleLag[cycleLag != 0], ])
  }
  if (errorAR > 0) {
    tmp[[nameE2]] <- rbind(tmp[[nameE2]], dfSystem[dfSystem$equation == "E2error" & dfSystem$variant == paste0("errorAR", errorAR), ])
  }
  if (errorMA > 0) {
    tmp[[nameE2]] <- rbind(tmp[[nameE2]], dfSystem[dfSystem$equation == "E2error" & dfSystem$variant == "errorMA", ][1:errorMA, ])
  }
  if (!is.null(type) && type == "NKP" && all(cycleLag == c(0, 1))) {
    tmp[[nameE2]] <- tmp[[nameE2]][!(tmp[[nameE2]]$sysMatrix == "Z" & tmp[[nameE2]]$variableRow == "cycleLag1"), ]
  }

  tmp
}

# -------------------------------------------------------------------------------------------

#' Initializes the location file containing default parameter constraints, among other things.
#'
#' @param model A model of class \code{NAWRUmodel} or \code{TFPmodel}.
#'
#' @return A data frame containing information on each involved parameter, for instance
#'   its corresponding system matrix, variable names, and parameter restrictions.
#' @keywords internal
.initializeLoc <- function(model) {

  # model attributes
  trend <- attr(model, "trend")
  cycle <- attr(model, "cycle")

  # initialize
  tmp <- list()

  # load model data
  namesExtract <- c("equation", "variant", "sysMatrix", "variableRow", "varName", "LB", "UB", "statRestr", "distribution")
  tmp <- .accessDfSystem(model = model)
  tmp <- lapply(tmp, function(x) {
    x[, namesExtract]
  })

  # merge equations
  tmpAll <- rbind(tmp$cycle, tmp$trend, tmp[[3]])

  # create column "restriction"
  loc <- cbind(tmpAll, data.frame(restriction = rep("NA", dim(tmpAll)[1]), stringsAsFactors = FALSE))

  for (k in tmpAll$varName) {
    lb <- tmpAll[tmpAll$varName == k, "LB"]
    ub <- tmpAll[tmpAll$varName == k, "UB"]
    stat <- tmpAll[tmpAll$varName == k, "statRestr"]
    if (!is.na(lb) & !is.na(ub)) {
      tmp <- paste0("I", lb, "_", ub)
    } else if (is.na(lb) & !is.na(ub)) {
      tmp <- paste0("I", "-Inf", "_", ub)
    } else if (!is.na(lb) & is.na(ub)) {
      tmp <- paste0("I", lb, "_", "Inf")
    } else if (stat != "") {
      tmp <- stat
    } else {
      tmp <- "NA"
    }
    loc$restriction[loc$varName == k] <- tmp
  }

  loc
}

# -------------------------------------------------------------------------------------------

#' Prints the model specifications.
#'
#' @param x An object of class \code{NAWRUmodel}, \code{TFPmodel}, or \code{KuttnerModel}.
#' @inheritParams print.NAWRUmodel
#' @keywords internal
.printSSModel <- function(x, call = TRUE, check = TRUE) {

  ARterms <- 0
  type <- NULL
  if (inherits(x, "NAWRUmodel")) {
    nameE2 <- "pcInd"
    nameE2long <- "phillips curve"
    checkFUN <- is.NAWRUmodel

    type <- attr(x, nameE2long)$type
  } else if (inherits(x, "TFPmodel")) {
    nameE2 <- "cubs"
    nameE2long <- "cubs"
    checkFUN <- is.TFPmodel

    ARterms <- attr(x, "cubs")$cubsAR
  } else if (inherits(x, "KuttnerModel")) {
    nameE2 <- "infl"
    nameE2long <- "inflation equation"
    checkFUN <- is.KuttnerModel
  }
  anchor <- attr(x, "anchor")$value
  anchor.h <- attr(x, "anchor")$horizon

  # attributes
  cycleLag <- attr(x, nameE2long)$cycleLag
  errorARMA <- attr(x, nameE2long)$errorARMA
  exoNames <- attr(x, nameE2long)$exoVariables
  trend <- attr(x, "trend")
  cycle <- attr(x, "cycle")
  freq <- attr(x, "period")$frequency
  start <- attr(x, "period")$start
  end <- attr(x, "period")$end
  n <- attr(x$SSModel, "n")

  digits <- 3
  if (call) {
    cat("Call:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n",
      sep = ""
    )
  }
  cat(paste0("\tState space model object of class ", class(x)[1], "\n\n"))
  cat(paste0("cycle ", "\t\t\t\t", cycle, "\n"))
  cat(paste0("trend ", "\t\t\t\t", trend, "\n"))
  cat(paste0(nameE2long, "\n"))
  if (!is.null(type)) {
    cat(paste0("  type ", "\t\t\t\t", type, "\n"))
  }
  cat(paste0("  cycle lags ", "\t\t\t", paste0(cycleLag, collapse = ","), "\n"))
  if (ARterms > 0) {
    cat(paste0("  AR terms ", "\t\t\t", paste0(1:ARterms, collapse = ","), "\n"))
  }
  cat(paste0("  error term", "\t\t\t", ifelse(any(errorARMA > 0),
    paste0("ARMA(", errorARMA[1], ",", errorARMA[2], ")"), "iid normal"
  ), "\n"))
  cat(paste0("  exogenous variables", "\t\t", ifelse(is.null(exoNames), "-", paste0(exoNames, collapse = ", ")), "\n"))
  if (!is.null(attr(x, "anchor"))) {
    cat(paste0("anchor\n"))
    cat(paste0("  value ", "\t\t\t", ifelse(is.null(anchor), "-", format(anchor, digits = digits)), "\n"))
    cat(paste0("  horizon ", "\t\t\t", ifelse(is.null(anchor.h), "-", anchor.h), "\n"))
  }
  cat("dimensions\n")
  cat(paste0("  number of observations", "\t", attr(x$SSModel, "n"), "\n"))
  cat(paste0("  period ", "\t\t\t", start, " - ", end, "\n"))
  cat(paste0("  frequency ", "\t\t\t", freq, "\n\n"))
  if (check) {
    if (checkFUN(x, return.logical = TRUE)) {
      cat(paste0("Object is a valid object of class ", class(x)[1], ".\n"))
    } else {
      checkFUN(x, return.logical = FALSE)
    }
  }
}

# -------------------------------------------------------------------------------------------

#' Prints the model fit and possibly specifications.
#'
#' @param x An object of class \code{NAWRUmodel}, \code{TFPmodel}, or \code{KuttnerModel}.
#' @param print.model A logical indicating whether the model specification should be printed.
#' @inheritParams print.NAWRUmodel
#' @keywords internal
.printSSModelFit <- function(x, call = TRUE, check = TRUE, print.model = TRUE) {

  if (attr(x, "method") == "bayesian") {
    bayes <- TRUE
    title <- "MCMC estimation results"
    dfIndex <- c(1:2, 4:5)
  } else {
    bayes <- FALSE
    title <- "Maximum likelihood estimation results"
    dfIndex <- 1:4
  }

  if (inherits(x$model, "NAWRUmodel")) {
    nameE2 <- "pcInd"
    nameE2long <- "phillips curve"
    nameE2short <- "pc"
  } else if (inherits(x$model, "TFPmodel")) {
    nameE2 <- "cubs"
    nameE2long <- "cubs"
    nameE2short <- "cu"
  } else if (inherits(x$model, "KuttnerModel")) {
    nameE2 <- "infl"
    nameE2long <- "inflation equation"
    nameE2short <- "infl"
  }

  digits <- 3
  if (call) {
    cat("Call:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  }
  # print nawru model
  .printSSModel(x = x$model, call = FALSE, check = FALSE)

  # print estimation results
  dfParBase <- x$parameters[, dfIndex]
  index <- list()
  index[c("cycle", "trend", nameE2, "E2error")] <- lapply(c("cycle", "trend", nameE2, "E2error"), function(y) {
    rownames(dfParBase) %in% x$model$loc$varName[x$model$loc$equation == y]
  })
  index[[nameE2]] <- index[[nameE2]] | index$E2error
  rownames(dfParBase) <- gsub("E2", nameE2short, rownames(dfParBase))
  rownames(dfParBase) <- paste0("  ", rownames(dfParBase))
  cat(paste0("\n\t", title, "\n\n"))
  # cycle
  cat("cycle\n")
  dfPar <- dfParBase[index$cycle, ]
  print(dfPar[order(rownames(dfPar)), ], digits = digits, quote = FALSE, right = TRUE, row.names = TRUE)
  # trend
  dfPar <- dfParBase[index$trend, ]
  if (dim(dfPar)[1] > 0) {
    cat("\ntrend\n")
    print(dfPar[order(rownames(dfPar)), ], digits = digits, quote = FALSE, right = TRUE, row.names = TRUE)
  }
  # E2
  cat(paste0("\n", nameE2long, "\n"))
  dfPar <- dfParBase[index[[nameE2]], ]
  print(dfPar[order(rownames(dfPar)), ], digits = digits, quote = FALSE, right = TRUE, row.names = TRUE)
  if (bayes) {
    print(unlist(x$fit[c("MRMSE","signal-to-noise")]), digits = 4)
  } else {
    cat(paste0("  RMSE: ", format(x$fit$RMSE, digits = digits), "\n"))
    cat(paste0("  R2: ", format(x$fit$R2, digits = digits), "\n"))
    cat(paste0(
      "  Box-Ljung test: X-squared = ", format(x$fit$LjungBox$statistic, digits = digits),
      ", df = ", format(x$fit$LjungBox$parameter, digits = digits),
      ", p-value = ", format(x$fit$LjungBox$p.value, digits = digits), "\n\n"
    ))
    print(unlist(x$fit[c("loglik", "AIC", "BIC", "HQC", "signal-to-noise")]), digits = 4)
  }
}



# -------------------------------------------------------------------------------------------

#' Computes additional results of the Kalman filter and smoother.
#'
#' @param out The return object of the function \code{KFS} from the package \code{KFAS}.
#' @param model An object of class \code{NAWRUmodel}, \code{TFPmodel}, or \code{KuttnerModel}.
#' @param prediction Logical indicating whether \code{out} includes predictions.
#' @importFrom KFAS mvInnovations
#' @importFrom stats coef ts start frequency
#' @keywords internal
.SSresults <- function(out, model, prediction = FALSE) {
  namesObs <- colnames(out$model$y)
  namesState <- colnames(out$model$T)

  start <- start(out$model$y)
  freq <- frequency(out$model$y)

  tsl <- list()
  tsl$obs <- out$model$y
  tsl$obsFitted <- out$m
  tsl$obsFittedSE <- ts(t(sqrt(apply(out$P_mu, 3, diag))), start = start, frequency = freq)
  tsl$stateSmoothed <- coef(out)
  tsl$stateFiltered <- out$att
  tsl$stateSmoothedSE <- ts(t(sqrt(apply(out$V, 3, diag))), start = start, frequency = freq)
  tsl$stateFilteredSE <- ts(t(sqrt(apply(out$Ptt, 3, diag))), start = start, frequency = freq)
  colnames(tsl$stateSmoothedSE) <- colnames(tsl$stateFilteredSE) <- namesState

  # standardized one step ahead residuals (=recursive residuals)
  tmp <- mvInnovations(out)
  FF <- tmp$F
  FF <- FF[, , !apply(FF, 3, function(x) any(is.na(x)))]
  v <- tmp$v
  B <- array(apply(FF, 3, function(x) chol(solve(x))), dim = dim(FF))

  vstd <- sapply(1:dim(FF)[3], function(x) B[, , x] %*% v[x, ])
  tsl$obsResidualsRecursive <- ts(t(vstd), start = start, frequency = freq)
  colnames(tsl$obsResidualsRecursive) <- namesObs[1:dim(v)[2]]
  
  # model specific series
  if (inherits(model, "NAWRUmodel")) {
    
    tsl$nawru <- with(tsl, stateSmoothed[, "trend"])
    tsl$nawruSE <- with(tsl, stateSmoothedSE[, "trend"])

  } else if (inherits(model, "TFPmodel")) {
    
    tsl$tfpTrend <- with(tsl, exp(stateSmoothed[, "trend"] / 100))
    tsl$tfpTrendGrowth <- with(tsl, growth(tfpTrend))
    tsl$tfp <- with(tsl, exp(obs[, 1] / 100))
    tsl$tfpGrowth <- with(tsl, growth(tfp))

    tsl[c("tfpTrendSE", "tfpTrendGrowthSE")] <- .deltaMethodState(out = out, nameState = "trend", constant = 100)[c("expStateSE", "diffStateSE")]
    if (prediction) {
      tsl[c("tfpSE", "tfpGrowthSE")] <- .deltaMethodObs(out = out, nameObs = "logtfp", model = model, constant = 100)[c("expObsSE", "diffObsSE")]
    }
    
  } else if (inherits(model, "KuttnerModel")) {
    
    tsl$potential <- with(tsl, exp(stateSmoothed[, "trend"] / 100))
    tsl$potentialGrowth <- with(tsl, growth(potential))
    tsl$gap <- with(tsl, stateSmoothed[, "cycle"])
    tsl$gdp <- with(tsl, exp(obs[, 1] / 100))
    tsl$gdpGrowth <- with(tsl, growth(gdp))

    tsl[c("potentialSE", "potentialGrowthSE")] <- .deltaMethodState(out = out, nameState = "trend", constant = 100)[c("expStateSE", "diffStateSE")]
    tsl[[c("gapSE")]] <- .deltaMethodState(out = out, nameState = "cycle")$StateSE
    if (prediction) {
      tsl[c("gdpSE", "gdpGrowthSE")] <- .deltaMethodObs(out = out, nameObs = "loggdp", model = model, constant = 100)[c("expObsSE", "diffObsSE")]
    }

  }
  

  tsl
}

# -------------------------------------------------------------------------------------------

#' Computes additional results of the Kalman filter and smoother for Bayesian output.
#'
#' @param model The return object of the function \code{fitSSM} from the package \code{KFAS}.
#' @param state An array with the smoothed state.
#' @param obsFitted An array with the fitted observables.
#' @inheritParams fit.TFPmodel
#' @importFrom KFAS mvInnovations
#' @importFrom stats coef ts start frequency
#' @keywords internal
.SSresultsBayesian <- function(model, HPDIprob, state, obsFitted, FUN) {
  start <- start(model$y)
  freq <- frequency(model$y)
  tsl <- list()
  
  # smoothed states
  tsl$stateSmoothedSummary <- lapply(setdiff(colnames(state), "const"), function(x) {
    ts(mcmcSummary(x = t(state[, x, ]), HPDIprob = HPDIprob),
       start = start, frequency = freq
    )
  })
  names(tsl$stateSmoothedSummary) <- setdiff(colnames(state), "const")
  tsl$stateSmoothedSummary <- do.call(cbind, tsl$stateSmoothedSummary)
  # fitted obs
  tsl$obsFittedSummary <- lapply(colnames(obsFitted), function(x) {
    ts(mcmcSummary(x = t(obsFitted[, x, ]), HPDIprob = HPDIprob),
       start = start, frequency = freq
    )
  })
  names(tsl$obsFittedSummary) <- colnames(obsFitted)
  tsl$obsFittedSummary <- do.call(cbind, tsl$obsFittedSummary)
  
  tsl$stateSmoothed <- ts(apply(state, c(1, 2), FUN), start = start, frequency = freq)
  tsl$obsFitted <- ts(apply(obsFitted, c(1, 2), FUN), start = start, frequency = freq)
  tsl$obs <- model$y
  tsl
}


# --------------------------------------------------------------------------------------------------- #

#' Computes standard errors, t-statistics, and p-values for the estimated state space parameters
#' using the delta method.
#'
#' @param parOptim The vector of optimized parameters, without transformations.
#' @param hessian The hessian from the optimization.
#' @param loc A \code{3 x n} array where \code{n} is the number of optimized parameters. The array
#'     contains information on the estimated parameters's name (\code{loc[1, ]}), its location
#'     (\code{loc[2, ]}) and possible parameter constraints (\code{loc[3, ]}).
#'
#' @importFrom stats pnorm
#' @keywords internal
inference <- function(parOptim, hessian, loc) {

  # set up data frame for optimized parameters
  dfRes <- data.frame(
    estim = assignConstraints(parOptim, loc = loc),
    se = NA,
    tstat = NA,
    pvalue = NA,
    row.names = names(parOptim)
  )

  # compute fisher information and covariance matrix
  fisherInfo <- hessian
  CV <- tryCatch(
    {
      solve(fisherInfo)
    },
    error = function(cont) {
      return(NA)
    }
  )
  # not invertible
  if (all(is.na(CV))) {
    ind <- !(loc$boundaries | apply(hessian == 0, 2, all))
    CV <- tryCatch(
      {
        solve(fisherInfo[ind, ind])
      },
      error = function(cont) {
        message("The fisher information matrix is not invertible.")
        return(matrix(NA))
      }
    )
  # invertible
  } else {
    ind <- !(loc$boundaries | apply(hessian == 0, 2, all))
    CV <- CV[ind, ind]
  }
  
  if (!all(is.na(CV))) {
    # compute covariance for constrained parameters
    Covar <- computeCovar(par = parOptim[ind], loc = loc[ind, ], CV = CV)
    # check diagonal values
    if (any(diag(Covar) < 0)) {
      ind[ind] <- ind[ind] & diag(Covar) > 0
      Covar <- Covar[diag(Covar) > 0, diag(Covar) > 0]
    }
  
    # compute standard error, t statistic and p values
    dfRes[loc$varName[ind], "se"] <- sqrt(diag(Covar))[loc$varName[ind]]
    dfRes[loc$varName[ind], "tstat"] <- dfRes[loc$varName[ind], "estim"] / dfRes[loc$varName[ind], "se"]
    dfRes[loc$varName[ind], "pvalue"] <- 2 * (1 - pnorm(abs(dfRes[loc$varName[ind], "tstat"])))
  }
  
  colnames(dfRes) <- c("Coefficient", "Standard Error", "t-statistic", "p-value")
  dfRes
}


# -------------------------------------------------------------------------------------------

#' Computes figures regarding the model fit of the maximum likelihood estimation.
#'
#' @param out The return object of the function \code{KFS} from the package \code{KFAS}.
#' @param nPar A scalar specifying the number of estimated parameters.
#' @keywords internal
.SSmodelfit <- function(out, nPar) {
  nTime <- out$dims$n

  # one step ahead residuals
  residuals <- mvInnovations(out)$v[, 2]

  info <- list()
  info$loglik <- out$logLik
  info$AIC <- 2 * nPar - 2 * out$logLik
  info$BIC <- log(nTime) * nPar - 2 * out$logLik
  info$HQC <- 2 * nPar * log(log(nTime)) - 2 * out$logLik
  info$RMSE <- sqrt(mean(residuals^2, na.rm = TRUE)) # sqrt(mean((residuals(out)[,"cubs"])^2, na.rm = TRUE)) # sqrt(mean(out$v[,2]^2))
  info$R2 <- 1 - sum((residuals)^2, na.rm = TRUE) / sum((out$model$y[, 2] - mean(out$model$y[, 2], na.rm = TRUE))^2, na.rm = TRUE) # 1 - sum((residuals(out)[,"cubs"])^2, na.rm = TRUE) / sum( (out$model$y[,"cubs"] - mean(out$model$y[,"cubs"], na.rm = TRUE) )^2, na.rm = TRUE)
  info$LjungBox <- stats::Box.test(residuals, lag = min(10, length(residuals) / 5), type = "Ljung-Box") # see https://robjhyndman.com/hyndsight/ljung-box-test/

  info
}
# -------------------------------------------------------------------------------------------

#' Computes standard errors of the state using the delta method.
#'
#' @param out The return object of the function \code{KFS} from the package \code{KFAS}.
#' @param nameState The name of the state as character.
#' @param constant A constant used in the transformation functions.
#' @keywords internal
.deltaMethodState <- function(out, nameState, constant = 1) {
  tsl <- list()
  nTime <- out$dims$n
  tsStateSmoothed <- out$alphahat

  start <- start(out$model$y)
  freq <- frequency(out$model$y)

  indexState <- (colnames(tsStateSmoothed) %in% nameState)

  # function g = identity
  Dg <- matrix(0, nTime, nTime)
  diag(Dg) <- 1 / constant #tsStateSmoothed[, nameState]
  # variance of trend
  varState <- diag(out$V[indexState, indexState, 1:nTime])
  tsl$StateSE <- ts(sqrt(diag(t(Dg) %*% varState %*% Dg)), start = start, frequency = freq)

  # function g = exponential
  Dg <- matrix(0, nTime, nTime)
  diag(Dg) <- exp(tsStateSmoothed[, nameState] / constant) / constant
  # variance of trend
  varState <- diag(out$V[indexState, indexState, 1:nTime])
  tsl$expStateSE <- ts(sqrt(diag(t(Dg) %*% varState %*% Dg)), start = start, frequency = freq)

  # function g = difference
  Dg <- matrix(0, nTime, nTime)
  diag(Dg) <- -1 / constant
  diag(Dg[-nrow(Dg), -1]) <- 1 / constant
  Dg <- Dg[1:(nrow(Dg) - 1), ]
  # variance of trend
  varState <- diag(out$V[indexState, indexState, 1:nTime])
  tsl$diffStateSE <- ts(sqrt(diag(Dg %*% varState %*% t(Dg))), start = start + c(0, 1), frequency = freq)

  return(tsl)
}

# -------------------------------------------------------------------------------------------

#' Computes standard errors of the observation equation using the delta method (for forecast).
#'
#' @param out The return object of the function \code{KFS} from the package \code{KFAS}.
#' @param nameObs The name of the observation equation as character.
#' @param constant A constant used in the transformation functions.
#' @inheritParams .SSresults
#' @keywords internal
.deltaMethodObs <- function(out, nameObs, model, constant = 1) {
  tsl <- list()
  nTime <- out$dims$n
  nData <- attr(model$SSModel, "n")
  nForecast <-  nTime - nData
  V <- out$V_mu[, ,(nTime - nForecast + 1):nTime]
  y <- out$model$y
  tsObs <- y[, nameObs]
  
  start <- start(y) + c(0, nData)
  freq <- frequency(y)
  
  indexObs <- (colnames(y) %in% nameObs)
  
  # function g = identity
  Dg <- matrix(0, nForecast, nForecast)
  diag(Dg) <- 1 / constant #tsObs[, nameObs]
  # variance of trend
  varObs <- diag(V[indexObs, indexObs, ])
  tsl$ObsSE <- ts(sqrt(diag(t(Dg) %*% varObs %*% Dg)), start = start, frequency = freq)
  
  # function g = exponential
  Dg <- matrix(0, nForecast, nForecast)
  diag(Dg) <- exp(y[(nTime - nForecast + 1):nTime, nameObs] / constant) / constant
  # variance of trend
  varObs <- diag(V[indexObs, indexObs, ])
  tsl$expObsSE <- ts(sqrt(diag(t(Dg) %*% varObs %*% Dg)), start = start, frequency = freq)
  
  # function g = difference
  Dg <- matrix(0, nForecast, nForecast)
  diag(Dg) <- -1 / constant
  diag(Dg[-nrow(Dg), -1]) <- 1 / constant
  Dg <- Dg[1:(nrow(Dg) - 1), ]
  # variance of trend
  varObs <- diag(V[indexObs, indexObs, ])
  tsl$diffObsSE <- ts(sqrt(diag(Dg %*% varObs %*% t(Dg))), start = start + c(0, 1), frequency = freq)
  
  tsl <- lapply(tsl, function(x) {
    window(x, start = start(y), extend = TRUE)
  })
  
  return(tsl)
}


# -------------------------------------------------------------------------------------------

#' Trend anchor
#'
#' @description Computes the anchored trend given a fitted object of class \code{NAWRUfit},
#'   \code{TFPfit}, or \code{KuttnerFit}.
#'
#' @param fit An object of class \code{NAWRUfit}, \code{TFPfit}, or \code{KuttnerFit}.
#' @param anchor A numeric specifying the anchor value. If unspecified, \code{anchor} is
#'   taken from the object \code{fit} (if specified).
#' @param h An integer specifying the anchor horizon in the frequency of the underlying model.
#'   If unspecified, \code{h} is taken from the object \code{fit} (if specified).
#' @param returnFit A logical. If \code{TRUE}, an object of the same class as \code{fit}
#'   including the anchored trend is returned. If \code{FALSE}, only the anchored trend time
#'   series is returned.
#'
#' @return A fitted object if \code{returnFit = TRUE} or a time series with the anchored 
#'   trend.
#' @export
#' @importFrom stats start end window ts lag frequency time window<-
#' @examples
#' # define nawru model for France
#' data("gap")
#' tsList <- amecoData2input(gap$France)
#' model <- NAWRUmodel(tsl = tsList)
#'
#' # estimate nawru model
#' \donttest{
#' f <- fit(model = model)
#' 
#' # compute anchored nawru
#' anchoredNawru <- trendAnchor(fit = f, anchor = 6.5, h = 10)
#' }
trendAnchor <- function(fit, anchor = NULL, h = NULL, returnFit = FALSE) {
  if (attr(fit, "method") == "bayesian") {
    stop("Anchor only implemented for MLE.")
  }
  if (is.null(anchor)) {
    anchor <- attr(fit$model, "anchor")$value
  }
  if (is.null(h)) {
    h <- attr(fit$model, "anchor")$horizon
  }
  if (is.null(anchor) | is.null(h)) {
    stop("Please specify 'anchor' and/or 'h'.")
  }

  # filtering object
  out <- fit$SSMout

  # nawru
  trend <- out$alphahat[, "trend"]

  # timing
  start <- start(trend)
  end <- end(trend)
  freq <- frequency(trend)
  timetrend <- time(trend)
  nTime <- length(trend)

  # selection vector
  nState <- length(colnames(out$model$T))
  nObs <- length(rownames(out$model$Z))
  s <- rep(0, nState)
  s[which(colnames(out$model$T) == "trend")] <- 1

  # initialize
  Tpower <- array(NA, c(nState, nState, h))
  Tpower[, , 1] <- out$model$T[, , 1]
  varTmp <- array(NA, c(nState, nState, h))
  varTmp[, , 1] <- Tpower[, , 1] %*% out$Ptt[, , nTime] %*% t(Tpower[, , 1]) + out$model$R[, , 1] %*% out$model$Q[, , 1] %*% t(out$model$R[, , 1])
  varMultiplier1 <- varMultiplier <- out$model$R[, , 1] %*% out$model$Q[, , 1] %*% t(out$model$R[, , 1])

  # matrix power and variance of predicted states
  for (ii in 1:(h - 1)) {
    varMultiplier <- varMultiplier1 + out$model$T[, , 1] %*% varMultiplier %*% t(out$model$T[, , 1])
    Tpower[, , ii + 1] <- out$model$T[, , 1] %*% Tpower[, , ii]
    varTmp[, , ii + 1] <- Tpower[, , ii + 1] %*% out$Ptt[, , nTime] %*% t(Tpower[, , ii + 1]) + varMultiplier
  }

  # correction factor
  correction <- (anchor - t(s) %*% Tpower[, , h] %*% out$alphahat[nTime, ])

  # location of constant (for inverse computation of P)
  locConst <- (1:nState)[rowSums(out$P[, , nTime]) == 0] # locate constants
  locConstInv <- (1:nState)[-locConst]

  # initialize
  covTmp <- array(NA, c(nState, nState, nTime))
  AA <- array(0, c(nState, nState, nTime))
  weight <- rep(1, nTime + h)
  anchoredtrend <- ts(NA, start, end + c(0, h), freq)

  # first step backward (at nTime)
  covTmp[, , nTime] <- out$Ptt[, , nTime]
  weight[nTime] <- (t(s) %*% covTmp[, , nTime] %*% t(Tpower[, , h]) %*% s) / (t(s) %*% varTmp[, , h] %*% s)
  window(anchoredtrend, start = timetrend[nTime], end = timetrend[nTime]) <- out$alphahat[nTime, "trend"] + weight[nTime] * correction

  # create temporary Pinf with zeros in all times > out$d
  PinfTmp <- array(0, dim = c(nState, nState, nTime))
  PinfTmp[, , (1:out$d)] <- out$Pinf

  # loop backward
  for (tt in ((nTime - 1):(1))) { # (out$d))) {

    ### De Jong& Mackinnon (1988)
    Pinv <- tryCatch(
      {
        solve(out$P[locConstInv, locConstInv, tt + 1] + PinfTmp[locConstInv, locConstInv, tt + 1])
      },
      error = function(cont) {
        warning("Error covariance matrix of predicted states is singular. The anchored trend could not be computed.")
        return(NA)
      }
    )
    if (all(is.na(Pinv))) { return(NA) }
    AA[locConstInv, locConstInv, tt] <- out$Ptt[locConstInv, locConstInv, tt] %*% t(out$model$T[locConstInv, locConstInv, 1]) %*% Pinv
    covTmp[, , tt] <- AA[, , tt] %*% covTmp[, , tt + 1]

    weight[tt] <- (t(s) %*% covTmp[, , tt] %*% t(Tpower[, , h]) %*% s) / (t(s) %*% varTmp[, , h] %*% s)

    # place value
    window(anchoredtrend, start = timetrend[tt], end = timetrend[tt]) <- out$alphahat[tt, "trend"] + weight[tt] * correction
  }

  # loop forward
  for (tt in 1:(h - 1)) {
    weight[nTime + tt] <- (t(s) %*% varTmp[, , tt] %*% t(Tpower[, , h - tt]) %*% s) / (t(s) %*% varTmp[, , h] %*% s)
    window(anchoredtrend, start = timetrend[nTime] + tt / freq, end = timetrend[nTime] + tt / freq) <- t(s) %*% Tpower[, , tt] %*% out$alphahat[nTime, ] + weight[nTime + tt] * correction
  }

  # last step in forward loop
  window(anchoredtrend, start = timetrend[nTime] + h / freq, end = timetrend[nTime] + h / freq) <- t(s) %*% Tpower[, , h] %*% out$alphahat[nTime, ] + correction

  # return
  if (returnFit) {
    if (inherits(fit, "NAWRUfit")) {
      fit$tsl$nawruAnchored <- anchoredtrend
    } else if (inherits(fit, "TFPfit")) {
      fit$tsl$tfpTrendAnchored <- exp(anchoredtrend / 100)
    } else {
      fit$tsl$trendAnchored <- anchoredtrend
    }
    attr(fit$model, "anchor")$value <- anchor
    attr(fit$model, "anchor")$horizon <- h
    return(fit)
  } else {
    return(anchoredtrend)
  }
}


# -------------------------------------------------------------------------------------------

#' Checks the covariance matrix for invertibility and negative entries on the diagonal.
#'
#' @param fit Output of \code{fitSSM} from \code{KFAS}.
#' @param loc A data frame containing information on each involved parameter (list element
#'   of objects of class \code{NAWRUmodel}, \code{TFPmodel}, \code{KuttnerModel}).
#'
#' @keywords internal
.checkCV <- function(fit, loc) {
  tryCatch(
    {
      CV <- solve(fit$optim.out$hessian)
      Covar <- computeCovar(par = fit$optim.out$par, loc = loc, CV = CV)
      any(diag(Covar) < 0)
    },
    error = function(cont) {
      return(TRUE)
    }
  )
}
# -------------------------------------------------------------------------------------------

#' Checks whether estimated parameters lie on boundaries.
#'
#' @inheritParams .checkCV
#'
#' @keywords internal
.checkBoundaries <- function(fit, loc) {
  ub <- loc$UB
  lb <- loc$LB
  par <- assignConstraints(par = fit$optim.out$par, loc = loc)
  par_index_lb <- round((lb - par) / (ub - lb), 4) >= 0
  par_index_ub <- round((ub - par) / (ub - lb), 4) <= 0
  par_index <- par_index_lb | par_index_ub
  par_index[is.na(par_index)] <- FALSE
  loc$boundaries <- par_index
  loc
}
