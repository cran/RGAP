
# -------------------------------------------------------------------------------------------

#' Fit best production function model
#'
#' @description Finds the most suitable model for the NAWRU and the TFP trend according to 
#'   the BIC or the RMSE. The function computes the output gap based on the chosen models.
#'
#' @param tsl A list of time series objects, see details.
#' @param type The variance restriction type. Possible options are \code{"basic"},
#'   \code{"hp"}, see \code{initializeRestr}. The default is \code{type = "hp"}.
#' @param q Quantile for the Inverse Gamma distribution (only used if \code{type = "hp"}), 
#'   see \code{initializeRestr}. The default is \code{q = 0.01}.
#' @param method The estimation method. Options are maximum likelihood estimation \code{"MLE"}
#'   and bayesian estimation \code{"bayesian"}. If \code{method = c("MLE", "bayesian")} the 
#'   NAWRU is fitted by MLE and the TFP trend by Bayesian methods. The default is 
#'   \code{method = "MLE"}.
#' @param criterion Model selection criterion. Options are the Bayesian information criterion 
#'   \code{"BIC"} and the root mean squared error \code{"RMSE"}, both computed for the second 
#'   observation equation. The default is \code{criterion = "BIC"}. For Bayesian estimation 
#'   and  \code{criterion = "RMSE"}, the mean RMSE is used.
#' @param fast Boolean, indicating whether a "fast" procedure should be used, see details.
#' @param nModels Integer, the maximum number of models for each unobserved component model.
#' @param nawruPoss List with possible model specifications for the NAWRU, see details.
#' @param tfpPoss List with possible model specifications for the NAWRU, see details.
#' @param auto If \code{auto = "NAWRU"} or \code{auto = "TFP"}, the function only 
#'   finds the most suitable NAWRU or TFP model, respectively. The default is 
#'   \code{auto = "gap"}.
#'   
#' @details For \code{fast = TRUE}, the function pre-selects suitable models by applying the 
#'   following procedure: A HP-filtered trend is computed based on which the best trend and 
#'   cycle models are chosen according to the BIC. Also based on the HP trend, a variety of 
#'   different specifications for the second observation equation are estimated in a 
#'   univariate regression and the best models are selected via the BIC. The \code{nModels} 
#'   best models are subsequently estimated in the usual bivariate unobserved component 
#'   model. For \code{fast = FALSE}, a variety of models is estimated in the usual bivariate 
#'   unobserved component framework.
#' @details The input component \code{nawruPoss} is a list containing a (sub-) set of the 
#'   following components: 
#'   \describe{
#'   \item{maxCycleLag}{Maximum cycle lag included in the second observation equation.}
#'   \item{trend}{Trend model specification.}
#'   \item{cycle}{Cycle model specification.}
#'   \item{errorARmax}{Maximum autoregressive order of the error term in the second 
#'   observation equation.}
#'   \item{errorMAmax}{Maximum moving average order of the error term in the second 
#'   observation equation.}
#'   \item{type}{Type of Phillip's curve.}
#'   \item{exoNames}{Names of the exogenous variables potentially included in the Phillip's 
#'   curve (need to be included in the list of time series \code{tsl}).}
#'   \item{signal-to-noise}{Signal-to-noise ratio.}
#'   }
#' @details The input component \code{tfpPoss} is a list containing a (sub-) set of the 
#'   following components: 
#'   \describe{
#'   \item{maxCycleLag}{Maximum cycle lag included in the second observation equation.}
#'   \item{trend}{Trend model specification.}
#'   \item{cycle}{Cycle model specification.}
#'   \item{cubsARmax}{Maximum CUBS autoregressive order.}
#'   \item{errorARmax}{Maximum autoregressive order of the error term in the second 
#'   observation equation.}
#'   \item{errorMAmax}{Maximum moving average order of the error term in the second 
#'   observation equation.}
#'   \item{signal-to-noise}{Signal-to-noise ratio.}
#'   }
#' @details The list of time series \code{tsl} needs to have the following components 
#'  (plus those series included in the list component \code{exoNames} in \code{nawruPoss}):
#' \describe{
#'   \item{ur}{Unemployment rate.}
#'   \item{nulc}{Nominal Unit labor costs, if \code{type = "TKP"}.}
#'   \item{rulc}{Real unit labor costs, if \code{type = "NKP"}.}
#'   \item{tfp}{Total factor productivity.}
#'   \item{cubs}{Capacity utilization economic sentiment indicator.}
#'   \item{lfnd}{Labor force non-domestic (unit: 1000 persons).}
#'   \item{parts}{Participation rate.}
#'   \item{ahours}{Average hours worked (unit: hours).}
#'   \item{gdp}{Gross domestic product at constant prices (unit: bn National currency, code: OVGD).}
#'   \item{k}{Net capital stock at constant prices: total economy (unit: bn National currency, code: OKND).}
#'   \item{popw}{Population: 15 to 64 years (unit: 1000 persons, code: NPAN).}
#'   }
#' @details The set of tested models is extensive but not exhaustive. The best model is 
#'   solely based on convergence and the chosen criterion (RMSE or BIC). A manual check 
#'   of the results is highly recommended. 
#' @details In some cases, more than \code{nModels} are checked. For instance, if a 
#'   re-parametrized and regular AR(2) process are options for the cycle.
#'     
#' @return A list containing three components: \code{gap} (the best model of class \code{"gap"}), 
#' \code{tfp} (a nested list of TFP models, fitted objects and model fit criteria), \code{nawru} 
#' (a nested list of NAWRU models, fitted objects and model fit criteria). The lists \code{nawru}
#' and \code{tfp} contain a list of models, a list of fitted objects and a dataframe \code{info},
#' which contains
#'   \item{loglik}{log-likelihood function at optimum}
#'   \item{AIC}{Akaike information criterion}
#'   \item{BIC}{Bayesian information criterion}
#'   \item{HQC}{Hannan-Quinn information criterion}
#'   \item{RMSE}{Root mean squared error}
#'   \item{R2}{Coefficient of determination (R squared)}
#'   \item{signal-to-noise}{Signal-to-noise ratio}
#'   \item{LjungBox}{p-value of Ljung-Box test for autocorrelation (H0 = no autocorrelation)}
#'   \item{convergence}{0 indicates convergence of the optimization}
#'   \item{rrange}{relative range of trend series w.r.t original series}
#'   \item{neg}{1 indicates that negative values are present in the trend series}
#'   \item{rev}{relative excess volatility w.r.t. original series (stationary series)}
#'   \item{rsd}{relative standard deviation w.r.t. original series (stationary series)}
#'   \item{cor}{correlation between trend and original series (stationary series)}
#'   \item{msdtg}{mean standardized deviation (stationary trend)}
#'   \item{magtg}{mean absolute growth of trend (stationaty trend)}
#'   \item{drop}{1 indicates the model should be dropped}
#' 
#' @export
autoGapProd <- function(tsl, 
                        type = "hp",
                        q = 0.01,
                        method = "MLE",
                        criterion = "BIC",
                        fast = TRUE,
                        nModels = 5,
                        nawruPoss = list(maxCycleLag = 2,
                                        trend = c("RW2", "DT"),
                                        cycle = c("AR1", "AR2"),
                                        errorARmax = 1,
                                        errorMAmax = 0,
                                        type = c("TKP", "NKP"),
                                        exoNames = c("ws", "prod", "tot"),
                                        signalToNoise = NULL),
                        tfpPoss = list(maxCycleLag = 2,
                                      trend = c("RW2", "DT"),
                                      cycle = c("AR1", "AR2", "RAR2"),
                                      cubsARmax = 0,
                                      errorARmax = 1,
                                      errorMAmax = 0,
                                      signalToNoise = NULL),
                        auto = "gap") {
  
  
  nawruPossBase <- list(maxCycleLag = 5,
                        trend = c("RW2", "DT"),
                        cycle = c("AR1", "AR2"),
                        errorARmax = 2,
                        errorMAmax = 0,
                        type = c("TKP", "NKP"),
                        exoNames = c("ws", "prod", "tot"),
                        signalToNoise = NULL)
  tfpPossBase <- list(maxCycleLag = 5,
                      trend = c("RW2", "DT"),
                      cycle = c("AR1", "AR2", "RAR2"),
                      cubsARmax = 2,
                      errorARmax = 2,
                      errorMAmax = 0,
                      signalToNoise = NULL)
  nawruPoss <- c(nawruPoss, nawruPossBase[!(names(nawruPossBase) %in% names(nawruPoss))])
  tfpPoss <- c(tfpPoss, tfpPossBase[!(names(tfpPossBase) %in% names(tfpPoss))])
  
  if (length(method) == 1) {
    method <- rep(method, 2)
  }
  
  result <- list()
  
  # NAWRU ------------------------------------
  
  if (auto == "gap" || auto == "NAWRU") {
    if (fast == TRUE) {
      
      comb <-  autoNAWRUmodel(tsl = tsl, poss = nawruPoss, nModels = nModels)
      comb <- lapply(comb, function(x) {
        y <- names(x)
        index <- which(names(x) == "BIC")
        x <- x[-index] 
        x
      })

    } else {
      
      # nawru possibilities
      possN <- list()
      possN$cycleLag <- lapply(0:nawruPoss$maxCycleLag, function(x) (0:10)[0:x + 1])
      possN$trend <- nawruPoss$trend
      possN$cycle <- nawruPoss$cycle
      possN$errorAR <- c(0:nawruPoss$errorARmax)
      possN$errorMA <- c(0:nawruPoss$errorMAmax)
      possN$type <- nawruPoss$type
      possN$exoType[[1]] <- NULL
      if (!is.null(nawruPoss$exoNames)) {
        possN$exoType[[2]] <- initializeExo(varNames = nawruPoss$exoNames)
        possN$exoType[[2]][1, , "difference"] <- 2
        possN$exoType[[2]][2, , "difference"] <- 1
        possN$exoType[[2]][1, , "lag"] <- 0
        possN$exoType[[2]][2, , "lag"] <- 1
        # if (length(nawruPoss$exoNames) > 1) {
        #   for (k in 1:length(nawruPoss$exoNames)) {
        #     possN$exoType[[2 + k]] <- initializeExo(varNames = nawruPoss$exoNames[k])
        #     possN$exoType[[2 + k]][1, , "difference"] <- 2
        #     possN$exoType[[2 + k]][2, , "difference"] <- 1
        #     possN$exoType[[2 + k]][1, , "lag"] <- 0
        #     possN$exoType[[2 + k]][2, , "lag"] <- 1
        #   }
        # }
      }
      
      # number of models 
      comb <- expand.grid(lapply(possN, function(x) 1:length(x)))
      
      comb <- lapply(1:nrow(comb), function(x) { 
        tmp <- lapply(1:ncol(comb), function(y) { 
          (possN[[y]][[comb[x,y]]] )
        })
        names(tmp) <- names(possN)
        tmp 
      })
      
    }
    
    # define and fit models
    res <- helper_model_fit(FUNmodel = NAWRUmodel, 
                            FUNfit = fit.NAWRUmodel, 
                            tsl = tsl,
                            comb = comb, 
                            poss = nawruPoss, 
                            method = method[1], 
                            type = type,
                            q = q,
                            modelName = "NAWRU")
    
    # compute criteria for model selection
    crit <- helper_fit_comparison(fit = res$fit, 
                                  E1name = "ur", 
                                  E1Trendname = "nawru",
                                  fitBayes = res$fitBayes)

    
    # eliminate models
    info <- cbind(helper_model_comparison(models = res$model), crit$info)
    info$drop <- 1 * (info$convergence != 0 
                      | info$neg == 1
                      | info$R2 < 0 
                      | info$LjungBox < 0.1)
    info <- info[order(info$drop, info[[criterion]]),]
  
    # order and save results
    index <-  as.numeric(rownames(info))
    nawru <- list()
    nawru$model <- res$model[index]
    nawru$fit <- res$fit[index]
    nawru$info <- info
    
    if (method[1] == "bayesian") {
      
      # eliminate models
      info <- cbind(helper_model_comparison(models = res$model), crit$infoBayes)
      info$drop <- 1 * (info$R2 < 0                   
                        | info$neg == 2)
      criterion_bayes <- grep(criterion, c("MRMSE", "R2"), value = TRUE)
      if (length(criterion_bayes) == 0) { criterion_bayes <- "MRMSE" }
      info <- info[order(info$drop, info[[criterion_bayes]]),]
      
      # order and save results
      index <-  as.numeric(rownames(info))
      nawru$modelBayes <- res$model[index]
      nawru$fitBayes <- res$fitBayes[index]
      nawru$infoBayes <- info
    }
    
    result$nawru <- nawru
    
  }
  
  # TFP ------------------------------------
  
  if (auto == "gap" || auto == "TFP") {
    if (fast == TRUE) {
      
      comb <- autoTFPmodel(tsl = tsl, poss = tfpPoss, nModels = nModels)
      comb <- lapply(comb, function(x) {
        y <- names(x)
        index <- which(names(x) == "BIC")
        x <- x[-index] 
        x
      })

    } else {
      
      possT <- list()
      possT$cycleLag <- lapply(0:tfpPoss$maxCycleLag, function(x) (0:10)[0:x + 1])
      possT$trend <- tfpPoss$trend
      possT$cycle <- tfpPoss$cycle
      possT$cubsAR <- tfpPoss$cubsARmax
      possT$errorAR <- c(0:tfpPoss$errorARmax)
      possT$errorMA <- c(0:tfpPoss$errorMAmax)
      
      comb <- expand.grid(lapply(possN, function(x) 1:length(x)))
      
      comb <- lapply(1:nrow(comb), function(x) { 
        tmp <- lapply(1:ncol(comb), function(y) { 
          (possN[[y]][[comb[x,y]]] )
        })
        names(tmp) <- names(possN)
        tmp })
      
    }
    
    # define and fit models
    res <- helper_model_fit(FUNmodel = TFPmodel, 
                            FUNfit = fit.TFPmodel, 
                            tsl = tsl,
                            comb = comb, 
                            poss = tfpPoss, 
                            method = method[2], 
                            type = type,
                            q = q,
                            modelName = "TFP")
    
    # compute criteria for model selection
    crit <- helper_fit_comparison(fit = res$fit, 
                                  E1name = "logtfp", 
                                  E1Trendname = "tfpTrend",
                                  E1trans = exp,
                                  fitBayes = res$fitBayes)
    
    # eliminate models
    info <- cbind(helper_model_comparison(models = res$model), crit$info)
    info$drop <- 1 * (info$convergence != 0 
                      | info$neg == 1
                      | info$R2 < 0 
                      | info$LjungBox < 0.1)
    info <- info[order(info$drop, info[[criterion]]),]

    # order and save results
    index <-  as.numeric(rownames(info))
    tfp <- list()
    tfp$model <- res$model[index]
    tfp$fit <- res$fit[index]
    tfp$info <- info
    
    if (method[2] == "bayesian") {
      
      # eliminate models
      info <- cbind(helper_model_comparison(models = res$model), crit$infoBayes)
      info$drop <- 1 * (info$R2 < 0 
                        | info$neg == 1)
      criterion_bayes <- grep(criterion, c("MRMSE", "R2"), value = TRUE)
      if (length(criterion_bayes) == 0) { criterion_bayes <- "MRMSE" }
      info <- info[order(info$drop, info[[criterion_bayes]]),]
      
      # order and save results
      index <-  as.numeric(rownames(info))
      tfp$modelBayes <- res$model[index]
      tfp$fitBayes <- res$fitBayes[index]
      tfp$infoBayes <- info
    }
    
    result$tfp <- tfp
  
  }
  
  # gap --------------------------------------
  
  if (auto == "gap") {
    
    if (nrow(tfp$info) >= 1 & nrow(nawru$info) >= 1) {
      if (method[1] == "MLE") {
        nawrufit <- nawru$fit[[1]]
      } else {
        nawrufit <- nawru$fitBayes[[1]]
      }
      if (method[2] == "MLE") {
        tfpfit <- tfp$fit[[1]]
      } else {
        tfpfit <- tfp$fitBayes[[1]]
      }
      bestgap <- gapProd(tsl = tsl, NAWRUfit = nawrufit, TFPfit = tfpfit)
    } else {
      warning("No valid model left.")
      bestgap <- NULL
    }
    result$gap <- bestgap
    
  }

  return(result)
}

# -------------------------------------------------------------------------------------------

#' NAWRU model suggestion
#'
#' @description Finds the most suitable NAWRU models.
#'
#' @param poss A list with the characteristics of possible models (see \code{autoGapProd)}
#' @inheritParams autoGapProd
#'
#' @return A nested list with one model specification per list entry.
#' @keywords internal
autoNAWRUmodel <- function(tsl, poss, nModels = 10) {
  
  # trend and cycle
  model <- NAWRUmodel(tsl = tsl)
  trend <- trendOptim(x = model$tsl$ur, opt = poss$trend)
  cycle <- cycleOptim(x = model$tsl$ur, opt = poss$cycle)
  
  # second observation equation
  typel <-list()
  for (k in poss$type) {
    model <- suppressWarnings(NAWRUmodel(tsl = tsl, type = k))
    xexo <- switch(k,
                   "TKP" = tsl[poss$exoNames],
                   "NKP" = NULL)
    typel[[k]] <- obs2Optim(x1 = model$tsl$ur, 
                            x2 = model$tsl$pcInd, 
                            xexo = xexo,
                            errorARmax = poss$errorARmax, 
                            errorMAmax = poss$errorMAmax, 
                            maxCycleLag = poss$maxCycleLag, 
                            maxAR = 0,
                            nModels = nModels)
  }
  modell <- c(typel[[1]],typel[[2]])
  
  # sort and discard models
  bic <- unlist(lapply(modell, "[[", "BIC"))
  index <- order(bic)[1:nModels]
  modell <- modell[index]
  
  for (k in 1:nModels) {
    # create exoType input
    if (!is.null(modell[[k]]$exo) && !is.na(modell[[k]]$exo)) {
      exoType <- initializeExo(varNames = unique(modell[[k]]$exo$var))
      for (k1 in unique(modell[[k]]$exo$var)) {
        if (k1 %in% modell[[k]]$exo$var) {
          exo_tmp <- modell[[k]]$exo[modell[[k]]$exo$var == k1, ]
          for (k2 in 1:nrow(exo_tmp)) {
            exoType[k2, k1, "difference"] <- exo_tmp$diff[k2]
            exoType[k2, k1, "lag"] <- exo_tmp$lag[k2]
          }
        }
      }
      modell[[k]]$exoType <- exoType
      modell[[k]]$exo <- NULL
      
    } 
    # assign type
    modell[[k]]$type <- rep(poss$type, each = nModels)[index][k]
  }
  # assign trend, cycle
  nModels <- length(modell)
  modell <- rep(modell, length(cycle))
  for (k in 1:length(cycle)) {
    for (j in 1:nModels) {
      modell[[j]]$cycle <- cycle[k]
    }
  }
  nModels <- length(modell)
  modell <- rep(modell, length(trend))
  for (k in 1:length(trend)) {
    for (j in 1:nModels) {
      modell[[j]]$trend <- trend[k]
    }
  }
  
  modell <- lapply(modell, function(x) {
    index <- which(names(x) == "errorARMA")
    names(x)[index] <- "pcErrorARMA"
    index <- which(names(x) == "ar")
    x <- x[-index]
    x
  })
  
  return(modell)
  
}

# -------------------------------------------------------------------------------------------

#' TFP model suggestion
#'
#' @description Finds the most suitable TFP models.
#'
#' @param poss A list with the characteristics of possible models (see \code{autoGapProd)}
#' @inheritParams autoGapProd
#'
#' @return A nested list with one model specification per list entry.
#' @keywords internal
autoTFPmodel <- function(tsl, poss, nModels = 10) {
  
  # trend and cycle
  model <- TFPmodel(tsl = tsl)
  trend <- trendOptim(x = model$tsl$logtfp, opt = poss$trend)
  cycle <- cycleOptim(x = model$tsl$logtfp, opt = poss$cycle)
  
  # second observation equation
  modell <- obs2Optim(x1 = model$tsl$logtfp, 
                      x2 = model$tsl$cubs, 
                      xexo = NULL,
                      errorARmax = poss$errorARmax, 
                      errorMAmax = poss$errorMAmax, 
                      maxCycleLag = poss$maxCycleLag, 
                      maxAR = poss$cubsARmax,
                      nModels = nModels)
  
  # sort by bic
  bic <- unlist(lapply(modell, "[[", "BIC"))
  
  # assign trend and cycle
  nModels <- length(modell)
  modell <- rep(modell, length(cycle))
  for (k in 1:length(cycle)) {
    for (j in 1:nModels) {
      modell[[j]]$cycle <- cycle[k]
    }
  }
  nModels <- length(modell)
  modell <- rep(modell, length(trend))
  for (k in 1:length(trend)) {
    for (j in 1:nModels) {
      modell[[j]]$trend <- trend[k]
    }
  }
  
  
  modell <- lapply(modell, function(x) {
    index <- which(names(x) == "errorARMA")
    names(x)[index] <- "cubsErrorARMA"
    index <- which(names(x) == "ar")
    names(x)[index] <- "cubsAR"
    index <- which(names(x) == "exo")
    x <- x[-index]
    x
  })
  
  
  return(modell)
  
}

# -------------------------------------------------------------------------------------------

#' model selection helper function
#'
#' @description Defines and fits models given the input parameters
#'
#' @param FUNmodel The model definition function.
#' @param FUNmodel The model fitting function.
#' @param poss A list with the characteristics of possible models (see \code{autoGapProd)}
#' @param comb A nested list with one model specification per list entry.
#' @param modelName Name of the model, i.e., NAWRU or TFP.
#' @inheritParams autoGapProd
#'
#' @return A nested list with the models and fitted objects.
#' @keywords internal
helper_model_fit <- function(FUNmodel, FUNfit, tsl, comb, poss, method, type, q, modelName) {
  n_poss <- length(comb)
  
  # define models
  model <- f <- fBayes <- list()
  valid <- rep(NA, n_poss)
  message(paste0("\nDefining ", n_poss," ", modelName, " models ..."))
  pb <- utils::txtProgressBar(min = 0, max = n_poss, style = 3)
  for (k in 1:n_poss) {
    if (k %% (n_poss/1) == 0) {
      utils::setTxtProgressBar(pb, k)
    }
    model_tmp <- tryCatch(
      {
        tmp <- suppressWarnings(do.call(FUNmodel, c(list(tsl = tsl), comb[[k]])))
        tmp
      },
      error=function(cond) {
        return(NULL)
      }
    )    
    model[[k]] <- model_tmp
    valid[k] <- !is.null(model_tmp)
  }
  model <- model[valid]
  comb <- comb[valid]
  n_poss <- length(comb)
  
  # find and delete duplicates
  valid <- rep(NA, n_poss)
  valid[n_poss] <- TRUE
  for (k in 1:(n_poss-1)) {
    valid_tmp <- NULL
    for (k2 in (k+1):n_poss) {
      valid_tmp[k2-k] <- tryCatch(
        {
          !all(unlist(attributes(model[[k]])) == unlist(attributes(model[[k2]])))
        },
        error=function(cond) {
          return(TRUE)
        },
        warning=function(cond) {
          return(TRUE)
        }
      )   
    }
    valid_tmp
    valid[k] <- all(valid_tmp)
  }
  model <- model[valid]
  comb <- comb[valid]
  n_poss <- length(comb)
  message(paste0("\n", n_poss," valid ", modelName, " models remain."))
  
  # fit models
  message(paste0("\nFitting ", n_poss," ", modelName, " models ..."))
  pb <- utils::txtProgressBar(min = 0, max = n_poss, style = 3)
  for (k in 1:n_poss) {
    if (k %% 1 == 0) {
      utils::setTxtProgressBar(pb, k)
      cat("\n")
    }
    parRestr <- initializeRestr(model = model[[k]], type = type, q = q)
    f[[k]] <- FUNfit(model = model[[k]], 
                       parRestr = parRestr,
                       signalToNoise = poss$signalToNoise,
                       control = list(maxit = 1000))
    if (method == "bayesian") {
      prior <- initializePrior(model = model[[k]])
      fBayes[[k]] <- tryCatch(
        {     FUNfit(model = model[[k]], 
                     parRestr = parRestr, 
                     prior = prior,
                     signalToNoise = poss$signalToNoise,
                     method = method,
                     R = 1000,
                     MLEfit = f[[k]])
        },
        error=function(cond) {
          return(NA)
        }
      )  
    }
  }
  
  res <- list(model = model,
              fit = f,
              fitBayes = fBayes)
  return(res)
  
}

# -------------------------------------------------------------------------------------------

#' model selection fit comparison helper function
#'
#' @description Defines and fits models given the input parameters
#'
#' @param fit Fitted object.
#' @param E1name Name of first observation equation.
#' @param E1Trendname Name of trend of first observation equation.
#' @param E1trans Transformation function for the first observation equation..
#' @param fitBayes Fitted Bayesian object.
#'
#' @return A data frame with information criteria, other goodness-of-fit measures, 
#'   convergence status, trend volatility measures.
#' @keywords internal
helper_fit_comparison <- function(fit, E1name, E1Trendname, E1trans = identity, fitBayes) {
  n_poss <- length(fit)
  
  
  # fit criteria
  crit <- names(fit[[1]]$fit)
  crit <- c(crit["LjungBox" != crit])
  info <- data.frame(matrix(NA, length(fit), length(crit)))
  colnames(info) <- crit
  for (k in 1:(length(crit))) {
    info[, k] <- unlist(lapply(lapply(fit, "[[", "fit"), "[[", crit[k]))
  }
  info$LjungBox <- unlist(lapply(lapply(lapply(fit, "[[", "fit"), "[[", "LjungBox"), "[[","p.value"))
  info$convergence <- unlist(lapply(lapply(lapply(fit, "[[", "SSMfit"), "[[", "optim.out"), "[[", "convergence"))
  
  # trend volatility measures
  df <- data.frame()
  for (k in 1:n_poss) {
    
    nDiff <- switch(attr(fit[[k]]$model, "trend"),
                    "RW1" = 1,
                    "RW2" = 2, 
                    "DT" = 1)
    tmp <- trendVolaMeasures(tsOriginal = E1trans(fit[[k]]$model$tsl[[E1name]]), 
                             tsTrend = fit[[k]]$tsl[[E1Trendname]], 
                             nDiff = nDiff)
    tmp <- as.data.frame(tmp)
    df <- rbind(df, tmp)
  }
  info <- cbind(info, df)
  info_mle <- info
  
  info <- NULL
  if (length(fitBayes)!=0) {
    
    nFit <- length(fitBayes)
    # index_na <- c(1:length(fitBayes))[unlist(lapply(fitBayes, function(x) all(is.na(x))))]
    index_na_invers <- c(1:length(fitBayes))[unlist(lapply(fitBayes, function(x) !all(is.na(x))))]
    fitBayes <- fitBayes[unlist(lapply(fitBayes, function(x) !all(is.na(x))))]
    # fit criteria
    crit <- names(fitBayes[[1]]$fit)
    info <- data.frame(matrix(NA, length(fitBayes), length(crit)))
    colnames(info) <- crit
    for (k in 1:(length(crit))) {
      info[, k] <- unlist(lapply(lapply(fitBayes, "[[", "fit"), "[[", crit[k]))
    }
    
    df <- data.frame()
    for (k in 1:length(fitBayes)) {
      
      nDiff <- switch(attr(fitBayes[[k]]$model, "trend"),
                      "RW1" = 1,
                      "RW2" = 2, 
                      "DT" = 1)
      tmp <- trendVolaMeasures(tsOriginal = fitBayes[[k]]$model$tsl[[E1name]], 
                               tsTrend = fitBayes[[k]]$tsl[[E1Trendname]],
                               nDiff = nDiff)
      tmp <- as.data.frame(tmp)
      df <- rbind(df, tmp)
    }
    info <- cbind(info, df)
    # add NA rows for error runs.
    info_na <- data.frame(matrix(NA, nFit, ncol(info)))
    info_na[index_na_invers, ] <- info
    colnames(info_na) <- colnames(info)
    info <- info_na
    
  }
  
  res <- list(info = info_mle,
              infoBayes = info)
  return(res)
}

# -------------------------------------------------------------------------------------------

#' model comparison helper function
#'
#' @description Gets model attributes
#'
#' @param models A list of models.
#'
#' @return A data frame with model attributes
#' @keywords internal
helper_model_comparison <- function(models) {
  
  model_info <- lapply(models, function(x) {
    att <- attributes(x)
    xclass <- att$class[1]
    E2name <- ifelse(xclass == "NAWRUmodel", "phillips curve", "cubs")
    df <- data.frame(cycle = att$cycle,
                     trend = att$trend,
                     E2cycleLag = paste0(att[[E2name]]$cycleLag, collapse = ","),
                     E2exoVariables = paste0(att[[E2name]]$exoVariables, collapse = ", "),
                     E2errorAR = att[[E2name]]$errorARMA[1],
                     E2errorMA = att[[E2name]]$errorARMA[2])
    
    if (xclass == "NAWRUmodel") {
      df$E2type <- att[[E2name]]$type
    }
    df
  })
  do.call(rbind, model_info)
  
}


# -------------------------------------------------------------------------------------------

#' Trend volatility measures
#'
#' @description Computes trend volatility measures.
#'
#' @param tsOriginal The original time series.
#' @param tsTrend The trend time series.
#' @param nDiff Integer indicating the order of differencing applied to the input series.
#'
#' @return A list containing the different measures.
#' @importFrom stats cor
#' @keywords internal
trendVolaMeasures <- function(tsOriginal, tsTrend, nDiff) {

  result <- list()

  # negative values
  result$neg <- 1 * any(tsTrend  < 0)
  # relative range
  result$rrange <- (max(tsTrend) - min(tsTrend)) / (max(tsOriginal) - min(tsOriginal))
  
  # transform to stationary series  
  tsOriginal <- diff(tsOriginal, differences = nDiff)
  tsTrend <- diff(tsTrend, differences = nDiff)
  
  
  # rev: relative excess volatility (result close of above 1 -> excess vola)
  result$rev <- (max(tsTrend) - min(tsTrend)) / (max(tsOriginal) - min(tsOriginal))
  
  # rstd: relative standard deviations
  result$rsd <- sqrt(var(tsTrend, na.rm = TRUE)) / sqrt(var(tsOriginal, na.rm = TRUE))
  
  # corr: correlation
  result$corr <- cor(tsTrend, tsOriginal)
  
  # normalized mean absolute deviation of trend growth (diff) (the smaller the smoother the trend)
  result$msdtg <-  mean(abs(tsTrend - mean(tsTrend))) / mean(abs(tsTrend))
  
  # magtg: mean absolute growth of trend growth (diff)
  result$magtg <- mean(abs(growth(tsTrend)))
  
  result
  
}

# -------------------------------------------------------------------------------------------

#' Find suitable trend specification
#'
#' @description Finds the most suitable trend model according to the BIC.
#'
#' @param x A time series.
#' @param opt A character vector with the trend models to be tested. The default is 
#'   \code{opt = c("RW1", "RW2", "DT")}.
#'
#' @return A character string with the chosen trend model.
#' @importFrom stats BIC arima frequency
#' @importFrom zoo na.trim
#' @keywords internal
trendOptim <- function(x, opt = c("RW1", "RW2", "DT")) {
  x <- na.trim(x)
  
  freq <- frequency(x)
  lambda <- 1600 / ((4 / freq)^4)
  trend <- hpfilter(x, lambda = lambda)

  bic <- NULL
  if ("RW1" %in% opt) { 
    bic["RW1"] <- BIC(arima(x = diff(trend), order = c(0, 0, 0))) 
  }
  if ("RW2" %in% opt) { 
    bic["RW2"] <- BIC(arima(x = diff(trend, differences = 2), order = c(0, 0, 0))) 
  }
  if ("DT" %in% opt) { 
    bic["DT"] <- BIC(arima(x = diff(trend), order = c(1, 0, 0))) 
  }
  
  res <- which(min(bic) == bic)
  return(names(res))  
  
}

# -------------------------------------------------------------------------------------------

#' Find suitable cycle specification
#'
#' @description Finds the most suitable cycle model according to the BIC.
#'
#' @param x A time series.
#' @param opt A character vector with the cycle models to be tested. The default is 
#'   \code{opt = c("AR1", "AR2", "RAR2")}.
#'
#' @return A character string with the chosen cycle model.
#' @importFrom stats BIC arima frequency
#' @importFrom zoo na.trim
#' @keywords internal
cycleOptim <- function(x, opt = c("AR1", "AR2", "RAR2")) {
  x <- na.trim(x)
  
  freq <- frequency(x)
  lambda <- 1600 / ((4 / freq)^4)
  trend <- hpfilter(x, lambda = lambda)
  cycle <- x - trend

  bic <- NULL
  if ("AR1" %in% opt) { 
    bic["AR1"] <- BIC(arima(cycle, order = c(1, 0, 0), method = "ML", include.mean = FALSE)) 
  }
  if ("AR2" %in% opt | "RAR2" %in% opt) { 
    bic["AR2"] <- BIC(arima(cycle, order = c(2, 0, 0), method = "ML", include.mean = FALSE)) 
  }
  
  res <- which(min(bic) == bic)
  if (names(res) == "AR2" & "RAR2" %in% opt) {
    res <- c("AR2", "RAR2")
  } else {
    res <- names(res)
  }
  return(res)  
  
}
# -------------------------------------------------------------------------------------------

#' Find suitable 2nd observation specification
#'
#' @description Finds the most suitable model for the second observation equation according 
#' to the BIC.
#'
#' @param x1 A time series, the first observation equation.
#' @param x1 A time series, the second observation equation.
#' @param xexo (Optional) A (multiple) time series with exogenous variables.
#' @param errorARmax Integer, maximal AR order of the error process of the 2nd observation 
#'   equation.
#' @param errorMAmax Integer, maximal MA order of the error process of the 2nd observation 
#'   equation.
#' @param maxCycleLag Integer, maximal cycle lag included in the 2nd observation 
#'   equation.
#' @param maxAR Integer, maximal AR order of the time series \code{x2} in the 2nd observation 
#'   equation. \code{0} means that no lag is included.
#' @param nModels Integer, maximum number of models chosen to be fitted.
#'   
#' @return A list containing the chosen parameters: \code{errorAR, errorMA, cycleLag, ar, exo}.
#' @keywords internal
#' @importFrom stats BIC predict arima frequency window lag
obs2Optim <- function(x1, x2, xexo = NULL, errorARmax = 2, errorMAmax = 2, maxCycleLag = 2, maxAR = 2, nModels = 1) {
  
  maxExo <- 10
  
  freq <- frequency(x1)
  lambda <- 1600 / ((4 / freq)^4)
  trend <- hpfilter(x1, lambda = lambda)
  cycle <- list()
  cycle[1:(maxCycleLag + 1)] <- lapply(0:maxCycleLag, function(x) {
    window(stats::lag(x1 - trend, -x), end = end(x1 - trend))
  })
  xregAR <- regExo <- list()
  if (maxAR > 0) {
    xregAR <- list()
    xregAR[1:(maxAR)] <- lapply(1:maxAR, function(x) {
      window(stats::lag(x2, -x), start = start(x2), end = end(x1 - trend))
    })
  }
  
  xreg_base  <- cycle[[1]]
  bic_base <- BIC(arima(x2, order = c(0, 0, 0), xreg = xreg_base, method="ML"))
  
  
  # xexo pre selection
  tmp <- NULL
  bic <- NULL
  par_max <- length(x2) - 2 - 1
  exo_lag_max <- min(ifelse(freq == 1, 2, 8), par_max)
  if (!is.null(xexo)) {
    count <- 0
    for (k in 1:length(xexo)) {
      for (p1 in 0:exo_lag_max) { # lag
        for (p2 in 1:2) { # diff
          count <- count + 1
          xregExo <- cbind(xreg_base, window(diff(stats::lag(xexo[[k]], -p1), differences = p2), start = start(x2), end = end(x2)))
          mod <- arima(x2, order = c(0, 0, 0), xreg = xregExo, method="ML")
          bic[count] <- BIC(mod)
        }
      }
    }
    tmp <- cbind(bic, expand.grid(lapply(list( 1:2, 0:exo_lag_max, 1:length(xexo)), function(x) 1:length(x))))
    tmp <- tmp[tmp[,1] < bic_base, ]
    tmp <- tmp[order(tmp[, 1]), ][1:min(maxExo, nrow(tmp)), 2:4]
    names(tmp) <- c("diff", "lag", "var")
    tmp$lag <- tmp$lag-1
    tmp$var <- names(xexo)[tmp$var]
    tmp
    # only 2 lags/diff combs, otherwise there can be linearly dependent combinations
    tmpl <- list()
    for (k in unique(tmp$var)) {
      tmpl[[k]] <- tmp[tmp$var==k, ]
      tmpl[[k]] <- tmpl[[k]][1:min(nrow(tmpl[[k]]), 2), ]
    }
    tmp <- do.call(rbind, tmpl)
    tmp
    regExo <- list()
    for (k in 1:nrow(tmp)) {
      regExo[[k]] <- window(diff(lag(xexo[[tmp$var[k]]], -tmp$lag[k]), differences = tmp$diff[k]), start = start(x2), end = end(x2))
    }
    
  }
  
  # ARMA errors
  par_cycle_max <- floor(length(x2) * 3 / 4) - nrow(tmp) - errorARmax - errorMAmax - maxAR
  poss <- list(cycle = 0:(length(cycle)-1),
               xregAR = 0,
               AR = 0:errorARmax, 
               MA = 0:errorMAmax)
  if (maxAR > 0) {
    poss$xregAR = 0:maxAR
  }
  if (!is.null(xexo)) {
    poss <- c(poss, lapply(1:length(regExo), function(x) 1:2))
    names(poss)[length(poss) - (length(regExo)-1):0] <- paste0("exo", 1:length(regExo))
  }
  combinations <- expand.grid(lapply(poss, function(x) 1:length(x)))
  comb <- as.data.frame(sapply(1:length(poss), function(x) poss[[x]][combinations[,x]]))
  colnames(comb) <- colnames(combinations)
  bic <- ljungp <- mae_ratio <- rep(NA, nrow(comb))
  h <- min(10, length(x2) / 5) # see https://robjhyndman.com/hyndsight/ljung-box-test/
  for (p in 1:nrow(comb)) {
    
    pAR <- comb$AR[p]
    pMA <- comb$MA[p]
    xreg <- cbind(do.call(cbind, cycle[1:(comb$cycle[p] + 1)]),
                  do.call(cbind, regExo[comb[p, grepl("exo", names(poss))] == 1]))
    if (comb$xregAR[p] > 0) {
      xreg <- cbind(xreg, do.call(cbind, xregAR[1:comb$xregAR[p]]))
    }
    mod <- arima(x2, order = c(pAR, 0, pMA), xreg = xreg, method="ML")
    bic[p] <- BIC(mod)
    ljungp[p] <- Box.test(mod$residuals, type = "Ljung-Box", lag = h)$p.value
    
    # out of sample performance
    n_par <- ncol(xreg) + comb$AR[p] + comb$MA[p] + 1
    index <- 1:max(floor(length(x2) * 3 / 4), n_par + 1)
    mae_ratio[p] <- tryCatch(
      {
        mod <- arima(x2[index], order = c(comb$AR[p], 0, comb$MA[p]), xreg = as.matrix(xreg)[index,], method="ML")
        pred <- stats::predict(mod, newxreg = as.matrix(xreg)[-index,])
        mae_in <- 100 / length(index) * sum( abs( mod$residuals ), na.rm = TRUE ) 
        mae_out <- 100 / length(pred$pred) * sum( abs( x2[-index] - pred$pred ), na.rm = TRUE ) 
        mae_out / mae_in
      },
      error=function(cond) {
        return(Inf)
      },
      warning=function(cond) {
        return(Inf)
      }
    )    
    
  }
  if (sum(ljungp < 0.1) < nrow(comb)) { bic[ljungp < 0.1] <- Inf }
  if (sum(mae_ratio > 1.5 | ljungp < 0.1) < nrow(comb)) { bic[mae_ratio > 1.5] <- Inf }
  nModels <- min(nModels, length(bic))
  index <- order(bic)[1:nModels]
  bic_min <- bic[index]
  errorAR <- comb$AR[index]
  errorMA <- comb$MA[index]
  cycleLag <- comb$cycle[index] 
  ar <- comb$xregAR[index]
  exo <- as.list(rep(NA, nModels))
  if (!is.null(xexo)) {
    exo <- list()
    for (k in 1:nModels)
      exo[[k]] <- tmp[comb[index[k], grepl("exo", names(poss))] == 1, ]
  }
  
  
  res <- list(
    errorARMA = as.list(as.data.frame(rbind(errorAR, errorMA))),
    cycleLag = lapply(cycleLag, function(x) 0:x),
    ar = as.list(ar),
    exo = exo,
    BIC = as.list(bic_min)
  )
  res <- lapply(1:nModels, function(x) lapply(res, "[[", x))
  return(res)  
  
}



