# -------------------------------------------------------------------------------------------

#' Draws from the posterior the autoregressive parameter of a stationary AR(1) process
#' without starting values.
#'
#' @param Y a \code{Tn x 1} vector.
#' @param phi a scalar containing the last draw of the autoregressive parameter \eqn{\phi}.
#' @param sigma a scalar containing the innovation variance.
#' @param phi0 a scalar containing the prior mean for \code{phi}.
#' @param Q0 a scalar containing the prior precision for \code{phi}.
#' @param lb (optional) lower bound for \code{phi}.
#' @param ub (optional) upper bound for \code{phi}.
#'
#' @details The corresponding model is given by \eqn{Y_t = \phi Y_{t-1} + e_t}, where
#'   \eqn{e_t ~ N(0, \sigma)} with prior distribution
#'   \eqn{p(\phi) = N(\phi_0, 1/Q_0 )}.
#' @details Conditional on the variance \eqn{\sigma}, the posterior is normal and known.
#' @details Stationarity and box constraints are enforced. If the stationarity constraint
#'   is not fulfilled, the last draw is returned.
#'
#' @importFrom stats runif
#' @keywords internal
.postAR1 <- function(Y, phi, phi0, Q0, sigma, lb = -Inf, ub = Inf) {
  # gibbs step for phi of AR(1) process: Yt = phi Yt-1 + et, et ~ N(0, sigma)
  # draws for phi | sigma, Y

  # dimension
  Tn <- length(Y)
  # posterior parameters
  Qstar <- (sum(Y[2:Tn]^2) - Y[1]^2) * solve(sigma) + Q0
  phistar <- solve(Qstar) * (Y[2:Tn] %*% Y[1:(Tn - 1)] * solve(sigma) + phi0 * Q0)
  # draw from posterior with box constraints
  tmp <- pnorm(lb, mean = phistar, sd = sqrt(solve(Qstar))) +
    runif(1, 0, 1) * (pnorm(ub, mean = phistar, sd = sqrt(solve(Qstar))) - pnorm(lb, mean = phistar, sd = sqrt(solve(Qstar))))
  phiProp <- qnorm(tmp, mean = phistar, sd = sqrt(solve(Qstar)))
  # stationarity
  if (!all(abs(polyroot(z = c(1, -phiProp))) > 1)) {
    phiProp <- phi
  }

  names(phiProp) <- "phi"
  return(phiProp)
}

# -------------------------------------------------------------------------------------------

#' Draws from the posterior of the autoregressive paramteres of a stationary AR(p),
#' \eqn{p > 1} process without starting values.
#'
#' @param Y A \code{Tn x 1} vector with the time series.
#' @param phi a \code{1 x p} vector containing the last draw of the autoregressive parameters
#'   \eqn{\phi}. \code{p} has to be larger than one.
#' @param sigma a scalar containing the innovation variance.
#' @param phi0 a \code{1 x p} vector containing the prior mean for \code{phi}.
#' @param Q0 a \code{p x p} matrix containing the prior precision for \code{phi}.
#' @param lb (optional) \code{1 x p} vector with lower bounds for \code{phi}.
#' @param ub (optional) \code{1 x p} vector with upper bounds for \code{phi}.
#'
#' @details The corresponding model is given by
#'   \eqn{Y_t = \phi_1 Y_{t-1} + ... + \phi_p Y_{t-p} + e_t}, where
#'   \eqn{e_t ~ N(0, \sigma)} with prior distribution \eqn{p(\phi) = N(\phi_0, 1/Q_0 )}.
#' @details The posterior draw is obtained via a Metropolis Hastings step with proposal density
#'   \eqn{ q = \prod_{t=p+1}^Tn p(Y_t, \phi, \sigma, Y_{t-1}, ..., Y_{t-p} )} which is
#'   known due to conjugacy. The acceptance probability is given by
#'   \eqn{ \alpha = \min{1, p(Y_1, ... Y_p | \phi_r, \sigma) / p(Y_1, ... Y_p | \phi_{r-1}, \sigma)}}
#'   where the subscript \eqn{r} denotes the r-th draw. \eqn{p(Y_1, ... Y_p | \phi_r, \sigma)}
#'   is itself normal.
#' @details Stationarity and box constraints are enforced. If the constraints
#'   are not fulfilled, the last draw is returned.
#'
#' @importFrom stats runif
#' @keywords internal
.postARp <- function(Y, phi, phi0, Q0, sigma, lb = -Inf, ub = Inf) {
  phi <- drop(phi)
  phi0 <- drop(phi0)

  # dimension
  Tn <- length(Y)
  p <- length(phi)

  X <- sapply(1:p, function(x) {
    Y[seq((p + 1 - x), Tn - x)]
  })
  y <- Y[(p + 1):Tn]
  Qstar <- t(X) %*% X / sigma + Q0
  phistar <- solve(Qstar) %*% (t(X) %*% y / sigma + Q0 %*% phi0)

  # Draw theta from proposal density (procedure such that ub and lb hold)
  phiProp <- .mvrnorm(mu = phistar, sigma = solve(Qstar))

  # variance covariance of c_1 and c_2
  covCLast <- .covAR(k = p, phi = phi, sigma = sigma)
  covCProp <- .covAR(k = p, phi = phiProp, sigma = sigma)

  # Compute acceptance probability
  alpha_ic <- min(1, det(covCLast)^(1 / 2) / det(covCProp)^(1 / 2) *
    exp(1 / 2 * Y[1:p] %*% (solve(covCLast) - solve(covCProp)) %*% Y[1:p]))

  # stationarity and box constraints for phi
  if (!all(abs(polyroot(z = c(1, -phiProp))) > 1) | !(all(phiProp > lb) & all(phiProp < ub))) {
    alpha_ic <- 0
  }

  # Accept-Reject the current draw
  if (runif(1, 0, 1) < alpha_ic) {
    phir <- phiProp
    covCstar <- covCProp
  } else {
    phir <- phi
    covCstar <- covCLast
  }

  res <- list(
    phi = phir,
    sigmaY12 = covCstar,
    aProb = alpha_ic
  )
  return(res)
}

# -------------------------------------------------------------------------------------------

#' Draws from the posterior of the variance parameter of a random walk or a random walk with
#' constant or stochastic drift.
#'
#' @param Y A \code{Tn x 1} vector with the dependent variable.
#' @param sigmaDistr A \code{1 x k} matrix with prior distribution and box constraints for
#'   the innovation variance. The first two entries contain the prior hyperparameters and
#'   the last two entries the upper and lower bound.
#' @param sigmaLast A scalar containing the last draw of the innovation variance.
#' @param muDistr A \code{k x 1} matrix with prior distribution and box constraints for
#'   the constant trend. The first two entries contain the prior hyperparameters and
#'   the last two entries the upper and lower bound.
#'
#' @details If the process follows a random walk with constant drift, the two parameters are
#'   drawn sequentially (conditional on the other parameter). The constant is drawn from a
#'   normal posterior given by conjugacy.
#' @details The innovation variance is drawn from a Inverse-Gamma posterior given by
#' conjugacy.
#'
#' @importFrom stats runif pgamma
#' @keywords internal
.postRW <- function(Y, sigmaDistr, sigmaLast = NULL, muDistr = NULL) {

  # RW1 or RW2
  parr <- NULL
  diff <- 1
  if (is.null(muDistr) | length(muDistr) == 0) { # RW2
    diff <- 2
    muProp <- 0
  }

  # dimension
  Tn <- length(Y)

  # difference series
  y <- diff(Y, lag = 1, differences = diff)

  if (diff == 1) {
    # ----- constant
    # prior parameters
    mu0 <- muDistr[1, ]
    V0 <- muDistr[2, ]
    lb <- muDistr[3, ]
    ub <- muDistr[4, ]
    # posterior parameters
    Vstar <- 1 / ((Tn - 1) / sigmaLast + 1 / V0)
    mustar <- Vstar * (sum(y) / sigmaLast + mu0 / V0)
    # draw from posterior with box constraints
    tmp <- pnorm(lb, mean = mustar, sd = sqrt(Vstar)) +
      runif(1, 0, 1) * (pnorm(ub, mean = mustar, sd = sqrt(Vstar)) - pnorm(lb, mean = mustar, sd = sqrt(Vstar)))
    muProp <- qnorm(tmp, mean = mustar, sd = sqrt(Vstar))
    parr <- muProp
    names(parr) <- colnames(muDistr)
  }

  # ----- sigma
  # prior parameters
  s0 <- sigmaDistr[1, ]
  nu0 <- sigmaDistr[2, ]
  lb <- max(0, sigmaDistr[3, ])
  ub <- sigmaDistr[4, ]
  # posterior parameters
  nustar <- nu0 + Tn - diff
  sstar <- s0 + sum((y - muProp)^2) # seems low
  # draw from posterior with box constraints
  # x in [lb, ub] <=> 1/x in [1/ub, 1/lb]
  tmp <- pgamma(1 / ub, shape = nustar / 2, rate = sstar / 2) +
    runif(1, 0, 1) * (pgamma(1 / lb, shape = nustar / 2, rate = sstar / 2) - pgamma(1 / ub, shape = nustar / 2, rate = sstar / 2))
  sigmar <- 1 / qgamma(tmp, shape = nustar / 2, rate = sstar / 2)
  names(sigmar) <- colnames(sigmaDistr)

  parr <- c(sigmar, parr)
  res <- list("par" = parr)
  return(res)
}

# -------------------------------------------------------------------------------------------

#' Draws from the posterior of parameters of a normal model with normal inverse gamma prior.
#'
#' @param Y A \code{Tn x 1} vector with the dependent variable.
#' @param X A \code{Tn x k} matrix with the \code{k} explanatory variables.
#' @param betaLast A \code{k x 1} vector containing the last draw.
#' @param betaDistr A \code{4 x k} matrix with prior distribution and box constraints for
#'   the parameters of each variable. In each column, the first two entries contain the
#'   prior hyperparameters and the last two entries the upper and lower bound.
#' @param sigmaDistr A \code{1 x k} matrix with prior distribution and box constraints for
#'   the innovation variance. The first two entries contain the prior hyperparameters and
#'   the last two entries the upper and lower bound.
#'
#' @details Draws from the posterior are obtained by conjugacy of the Normal-Inverse-Gamma
#' distribution.
#'
#' @return A list with a draws for beta and the innovation variance.
#'
#' @importFrom stats runif pgamma
#' @keywords internal
.postNIG <- function(Y, X, betaLast, betaDistr, sigmaDistr) {
  betaLast <- drop(betaLast)
  # dimension
  Tn <- length(Y)

  # ----- beta and gamma
  # prior parameters
  betaDistr <- as.matrix(betaDistr)
  beta0 <- betaDistr[1, ]
  Q0 <- solve(diag(betaDistr[2, ]))
  s0 <- sigmaDistr[1]
  nu0 <- sigmaDistr[2]

  # posterior parameters
  Qstar <- Q0 + t(X) %*% X
  betastar <- solve(Qstar) %*% (Q0 %*% beta0 + t(X) %*% Y)
  nustar <- nu0 + Tn
  sstar <- s0 + t(Y) %*% Y + t(beta0) %*% Q0 %*% beta0 - t(betastar) %*% Qstar %*% betastar
  # sstar    <- s0 + t(Y - X %*% beta0) %*%  solve( diag(Tn) + X %*% solve(Q0) %*% t(X) ) %*% (Y - X %*% beta0)
  # (is the same)

  # dram from inverse gamma enacting constraints if present
  lb <- max(0, sigmaDistr[3])
  ub <- sigmaDistr[4]
  # x in [lb, ub] <=> 1/x in [1/ub, 1/lb]
  tmp <- pgamma(1 / ub, shape = nustar / 2, rate = sstar / 2) +
    runif(1, 0, 1) * (pgamma(1 / lb, shape = nustar / 2, rate = sstar / 2) - pgamma(1 / ub, shape = nustar / 2, rate = sstar / 2))
  sigmar <- 1 / qgamma(tmp, shape = nustar / 2, rate = sstar / 2)

  # Get draw for mean from normal distribution conditional on draw for gamma
  lb <- betaDistr[3, ]
  ub <- betaDistr[4, ]
  betar <- .mvrnorm(mu = betastar, sigma = solve(Qstar) * sigmar)
  if (!(all(betar > lb) & all(betar < ub))) {
    betar <- betaLast
  }

  # result
  res <- list(
    "beta" = betar,
    "sigma" = sigmar
  )
  return(res)
}
