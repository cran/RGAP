% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checks.R
\name{.checkBayesInput}
\alias{.checkBayesInput}
\title{Checks the input parameters of \code{.BayesFitTFP} and \code{.BayesFitNAWRU} for
consistency.}
\usage{
.checkBayesInput(
  model,
  type,
  prior = NULL,
  R = NULL,
  burnin = NULL,
  thin = NULL,
  HPDIprob = NULL,
  FUN = NULL,
  MLEfit = NULL
)
}
\arguments{
\item{model}{An object of class TFPmodel.}

\item{type}{A character specifying whether a "nawru" or "tfp" model should be checked.}

\item{prior}{A list of matrices with parameters for the prior distribution and box
constraints. By default, \code{prior} is initialized by \code{initializePrior(model)}.
See details. Only used if \code{method = "bayesian"}.}

\item{R}{An integer specifying the number of MCMC draws. The default is \code{R = 10000}.
Only used if \code{method = "bayesian"}.}

\item{burnin}{An integer specifying the burn-in phase of the MCMC chain. The default is
\code{burnin = ceiling(R / 10)}. Only used if \code{method = "bayesian"}.}

\item{thin}{An integer specifying the thinning interval between consecutive draws. The
default is \code{thin = 1}, implying that no draws are dopped. For \code{thin = 2},
every second draw is dropped and so on. Only used if \code{method = "bayesian"}.}

\item{HPDIprob}{A numeric in the interval \code{(0,1)} specifying the target probability
of the highest posterior density intervals. The default is \code{HPDIprob = 0.9}. Only
used if \code{method = "bayesian"}.}

\item{FUN}{A function to be used to compute estimates from the posterior distribution.
Possible options are \code{"mean"} and \code{"median"}. The default is \code{FUN = "mean"}.
Only used if \code{method = "bayesian"}.}

\item{MLEfit}{(Optional) An object of class \code{TFPfit} which is used for
initialization. Only used if \code{method = "bayesian"}.}
}
\description{
Checks the input parameters of \code{.BayesFitTFP} and \code{.BayesFitNAWRU} for
consistency.
}
\keyword{internal}
