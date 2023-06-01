
# -------------------------------------------------------------------------------------------

#' Computes the covariance of an AR(q) process.
#'
#' @param k integer indicating the lag length of the covariance.
#' @param phi \code{q x 1} vector of parameters.
#' @param sigma the innovation variance of the AR(q) process.
#'
#' @return A \code{k x k} covariance matrix.
#'
#' @importFrom stats ARMAacf toeplitz
#' @keywords internal
.covAR <- function(k, phi, sigma) {
  phi <- drop(phi)
  # covariance
  cov <- sigma * toeplitz(ARMAacf(phi, lag.max = k - 1))

  return(cov)
}



# -------------------------------------------------------------------------------------------

#' computes the unconditional variance of the cubs equation with p lags of cubs and k
#' additional lags of the cycle. The cycle follows an AR process of order l.
#'
#' @param mu constant.
#' @param phi \code{p x 1} vector of parameters for cubs lags.
#' @param beta \code{k x 1} vector of parameters for contemporaneous cycle and cycle lags.
#' @param sigma the innovation variance of the cubs process.
#' @param phiC \code{l x 1} vector of cycle process parameters.
#' @param sigmaC the innovation variance of the cycle process.

#' @return A list with a two covariance matrices; the first one depicts the covariance
#' between the cycle and cubs and the second onf the covariance of cubs.
#'
#' @importFrom stats toeplitz
#' @keywords internal
.covCUBS <- function(mu, phi, beta, sigma, phiC, sigmaC) {
  # computes the unconditional variance of the cubs equation with p lags of cubs and k
  # additional lags of the cycle. The cycle follows an AR process of order l.

  # we must have k >= pc - 1
  # if not, we can extend Phi1tilde

  phi <- drop(phi)
  beta <- drop(beta)

  # dimensions
  p <- length(phi)
  kbase <- length(beta) - 1
  pc <- length(phiC)
  k <- max(pc - 1, kbase) # add extra equations such that all variances can be retrieved

  # covariance of cycle process
  covC <- .covAR(k = k + 1, phi = phiC, sigma = sigmaC)

  # set autoregressive parameters to zero if missing
  if (is.null(phi)) {
    p <- length(phiC)
    phi <- rep(0, p)
  }
  q <- max(p, k)
  
  # set up matrix for Yule Walker equations for cov(c, cubs)
  Phi1tilde <- matrix(0, k + 1, p + k + 1)
  Phi2tilde <- matrix(0, q, pc + q)
  for (j in 1:(k + 1)) {
    Phi1tilde[j, j:(p + j)] <- c(1, -phi)
  }
  for (j in 1:q) {
    Phi2tilde[j, (j):(j + pc)] <- rev(c(1, -phiC))
  }
  Phitilde <- rbind(
    cbind(Phi1tilde, matrix(0, k + 1, q - p)),
    cbind(matrix(0, q, k + 1 - pc), Phi2tilde)
  )
  ytilde <- rep(0, k + q + 1)
  gammac <- c(rev(covC[1, 2:(k + 1)]), covC[1, ]) # -k:k
  for (j in (0:k)) {
    ytilde[j + 1] <- t(beta) %*% gammac[(1 + j):(kbase + 1 + j)]
  }
  gammatilde <- as.numeric(solve(Phitilde) %*% ytilde) # -k:q

  # set up matrix for Yule Walker equations for cov(cubs, cubs)
  Y1 <- Y2 <- matrix(0, p + 1, p + 1)
  Y1[1, ] <- c(0, phi)
  for (j in 2:(p + 1)) {
    Y1[j, 1:(p - j + 2)] <- -phi[(j - 1):p]
  }
  if (p > 1) {
    for (j in 3:(p + 1)) {
      Y2[j, 2:(j - 1)] <- rev(-phi[1:(j - 2)])
    }
  }
  Phi1 <- diag(1, p + 1) + Y1 + Y2
  Phi2 <- matrix(0, q - p, q)
  if (q > p) {
    for (j in 1:(q - p)) {
      Phi2[j, j:(p + j)] <- rev(c(1, -phi))
    }
  }
  Phi <- rbind(
    cbind(Phi1, matrix(0, p + 1, q - p)),
    cbind(matrix(0, q - p, 1), Phi2)
  )
  y <- rep(0, q + 1)
  for (j in (0:q)) {
    y[j + 1] <- t(beta) %*% rev(gammatilde[(1 + j):(kbase + 1 + j)])
  }
  y[1] <- y[1] + sigma
  gamma <- as.numeric(solve(Phi) %*% y) # 0:q

  # covariance
  cov <- toeplitz(gamma)
  result <- list(
    covCubs = cov,
    covCubsC = gammatilde
  )
  return(result)
}

# -------------------------------------------------------------------------------------------

#' Computes the approximate highest posterior density interval (HPDI).
#'
#' @param x A \code{R x n} matrix with \code{R} draws of \code{n} variables.
#' @param prob The probability mass of the interval, a scalar between zero and one.
#' @keywords internal
HPDinterval <- function(x, prob = 0.95) {
  x <- as.matrix(x)
  # order values
  xOrder <- apply(x, 2, sort)
  if (!is.matrix(xOrder)) stop("x must have nsamp > 1")
  # number of samples and parameters
  R <- nrow(xOrder)
  n <- ncol(xOrder)
  # number or sampled values included in the HPD interval
  Rprob <- max(1, min(R - 1, round(R * prob)))
  # compute ranges of possible intervals
  inds <- apply(
    xOrder[(Rprob + 1):R, , drop = FALSE] - xOrder[1:(R - Rprob), , drop = FALSE],
    2, which.min
  )
  # choose intervals with shortest range
  result <- cbind(
    xOrder[cbind(inds, 1:n)],
    xOrder[cbind(inds + Rprob, 1:n)]
  )
  # rename
  dimnames(result) <- list(colnames(x), c("lower", "upper"))
  attr(result, "Probability") <- Rprob / R
  # return
  return(result)
}

# -------------------------------------------------------------------------------------------

#' Conducts a Geweke test for convergence of the draws.
#'
#' @param x A \code{R x n} matrix with \code{R} draws of \code{n} variables.
#' @param frac1 The probability mass of the first interval, a scalar between zero and one.
#' @param frac2 The probability mass of the second interval, a scalar between zero and one.
#' @param alpha The significance level used to compute the test decision, a scalar between
#'   zero and one.
#'
#' @details Under the H0 of convergence, the test statistic is standard normally distributed.
#' @details Naturally, \code{frac1 + frac2} is between zero and one.
#'
#' @return A list with the following items
#' \describe{
#'   \item{h}{Test decision.}
#'   \item{CD}{Convergence Diagnostic (test statistic)}
#'   \item{pvalue}{The p-value.}
#'   \item{alpha}{The applied signifcicance level.}
#'   \item{frac1}{The fraction of data contained in the first interval.}
#'   \item{frac2}{The fraction of data contained in the second interval.}
#' }
#' @importFrom stats window start end spec.ar pnorm
#' @keywords internal
gewekeTest <- function(x, frac1 = 0.1, frac2 = 0.5, alpha = 0.05) {
  if (frac1 + frac2 > 1 | any(c(frac1, frac2) < 0) | any(c(frac1, frac2) > 1)) {
    stop("The input parameters 'frac1' and/or 'frac2' are invalid.")
  }
  if (alpha > 1 | alpha < 0) {
    stop("The input parameters 'alpha' is invalid.")
  }

  # number of variables
  x <- as.matrix(x)
  n <- ncol(x)

  # initialize
  CD <- pvalue <- h <- rep(NA, n)

  for (j in 1:n) {
    if (var(x[, j]) != 0) {
      # start and end dates
      x1start <- start(x[, j])
      x2start <- floor(end(x[, j]) - frac2 * (end(x[, j]) - start(x[, j])))
      x1end <- ceiling(start(x[, j]) + frac1 * (end(x[, j]) - start(x[, j])))
      x2end <- end(x[, j])

      # two series
      x1 <- window(x[, j], start = x1start, end = x1end)
      x2 <- window(x[, j], start = x2start, end = x2end)

      if ((var(x1) != 0) & (var(x2) != 0)) {

        # means
        m1 <- mean(x1)
        m2 <- mean(x2)

        # spectral densities
        sd1 <- spec.ar(x = x1, plot = FALSE)$spec[1]
        sd2 <- spec.ar(x = x2, plot = FALSE)$spec[1]

        # convergence diagnostic
        CD[j] <- (m1 - m2) / sqrt(sd1 / length(x1) + sd2 / length(x2))

        # p-value and test decision
        pvalue[j] <- 2 * (1 - pnorm(abs(CD[j])))
        h[j] <- 0 + (pvalue[j] < alpha)
      }
    }
  }

  # prepare results
  names(h) <- names(CD) <- names(pvalue) <- colnames(x)
  result <- list(
    h = h,
    CD = CD,
    pvalue = pvalue,
    alpha = alpha,
    frac1 = frac1,
    frac2 = frac2
  )
  return(result)
}

# -------------------------------------------------------------------------------------------

#' Computes MCMC summary statistics.
#'
#' @param x A \code{R x n} matrix with \code{R} draws of \code{n} variables.
#' @param HPDIprob The probability mass of the HPDI, a scalar between zero and one.
#' @param frac1 The probability mass of the first interval used for the Geweke test, a scalar
#'   between zero and one.
#' @param frac2 The probability mass of the second interval used for the Geweke test, a scalar
#'   between zero and one.
#'
#' @details Naturally, \code{frac1 + frac2} is between zero and one.
#'
#' @return A data frame with the following columns
#' \describe{
#'   \item{Mean}{The posterior mean.}
#'   \item{Median}{The posterior median.}
#'   \item{SD}{Standard deviation.}
#'   \item{HPDI-LB}{Highest posterior density credible interval lower bound}
#'   \item{HPDI-UB}{Highest posterior density credible interval upper bound}
#'   \item{Naive SE}{Naive Standard error of the mean (ignoring chain autocorrelation.}
#'   \item{Time-series SE}{Time-series standard error (based on spectral density at 0).}
#'   \item{Geweke statistic}{The Geweke test statistic.}
#'   \item{frac1}{The fraction of data contained in the first interval.}
#'   \item{frac2}{The fraction of data contained in the second interval.}
#' }
#' @importFrom stats median spec.ar
#' @keywords internal
mcmcSummary <- function(x, HPDIprob, frac1 = 0.1, frac2 = 0.5) {

  # number of variables
  x <- as.matrix(x)
  x <- x[, !apply(x, 2, function(a) all(is.na(a)))]
  n <- ncol(x)

  # initialize
  m <- md <- sd <- seNaive <- seTs <- tGeweke <- rep(NA, n)
  hpd <- matrix(NA, n, 2)

  for (j in 1:n) {
    # mean
    m[j] <- mean(x[, j])
    # median
    md[j] <- median(x[, j])
    if (var(x[, j]) != 0) {
      # standard deviation
      sd[j] <- sqrt(var(x[, j]))
      # naive standard errors
      seNaive[j] <- sqrt(var(x[, j]) / length(x[, j]))
      # spectral densities standard errors
      seTs[j] <- sqrt(spec.ar(x = x[, j], plot = FALSE)$spec[1] / length(x[, j]))
      # HPDI
      hpd[j, ] <- HPDinterval(x[, j], prob = HPDIprob)
      # Geweke test
      tGeweke[j] <- gewekeTest(x[, j], frac1 = frac1, frac2 = frac2)$CD
    }
  }

  # prepare results
  result <- data.frame(m, md, sd, hpd, seNaive, seTs, tGeweke, row.names = colnames(x))
  names(result) <- c(
    "Mean", "Median", "SD",
    paste0(HPDIprob * 100, "% HPDI-", c("LB", "UB")),
    "Naive SE",
    "Time-series SE",
    "Geweke statistic"
  )
  # attributes
  attr(result, "HPD probability") <- HPDIprob
  attr(result, "frac1") <- frac1
  attr(result, "frac2") <- frac2
  return(result)
}

# -------------------------------------------------------------------------------------------

#' Prints the results of the Geweke test for the states and the parameters.
#'
#' @param tsl A multiple time series list containing the Geweke test statistic for the
#'   states.
#' @param df A data frame containing the Geweke test statistic for the parameters.
#' @param alpha The significance level used to reach a test decision, a scalar between zero
#'   and one.
#' @keywords internal
.printGeweke <- function(tsl, df, alpha = 0.05) {
  gewekeFail <- tsl[, grepl("Geweke", colnames(tsl))] > qnorm(1 - alpha)
  if (sum(gewekeFail, na.rm = TRUE) == 0) {
    message("Convergence of states satisfied.")
  } else {
    message(paste0(
      "Convergence of states in the following year(s) not satisfied:\n  ",
      paste0(time(gewekeFail)[apply(gewekeFail, 1, any, na.rm = TRUE) == 1], collapse = ", "),
      "\n  If possible, increase burnin phase and/or the number of draws R."
    ))
  }
  gewekeFail <- df[, "Geweke statistic"] > qnorm(1 - alpha) | is.na(df[, "Geweke statistic"])
  if (sum(gewekeFail) > 0) {
    message(paste0(
      "Convergence of the following parameter(s) not satisfied:\n  ",
      paste0(rownames(df)[gewekeFail], collapse = ", "),
      "\n  If possible, increase burnin phase and/or the number of draws R."
    ))
  } else {
    message("Convergence of parameters satisfied.")
  }
}


# -------------------------------------------------------------------------------------------

#' Draws from the multivariate normal distribution.
#'
#' @param mu A \code{n x 1} vector, the mean vector.
#' @param sigma A \code{n x n} matrix, the covariance matrix.
#'
#' @importFrom stats rnorm
#' @keywords internal
.mvrnorm <- function(mu, sigma) {
  # tolerance regarding positive-definiteness of the covariance matrix
  tol <- 1e-08
  # check dimension
  n <- length(drop(mu))
  if (!all(dim(sigma) == c(n, n))) {
    stop("mean vector and covariance matrix not compatible")
  }
  # check covariance for positive definiteness
  eig <- eigen(sigma, symmetric = TRUE)
  if (!all(eig$values >= -tol * abs(eig$values[1L]))) {
    stop("The covariance matrix is not positive definite")
  }
  # cholesky
  L <- chol(sigma)
  # draw
  Y <- drop(mu) + t(L) %*% rnorm(n)
  # names
  name <- names(mu)
  dimnames(Y) <- list(name, NULL)
  # return
  return(Y)
}
