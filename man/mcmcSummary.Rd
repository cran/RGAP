% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_bayes.R
\name{mcmcSummary}
\alias{mcmcSummary}
\title{Computes MCMC summary statistics.}
\usage{
mcmcSummary(x, HPDIprob, frac1 = 0.1, frac2 = 0.5)
}
\arguments{
\item{x}{A \code{R x n} matrix with \code{R} draws of \code{n} variables.}

\item{HPDIprob}{The probability mass of the HPDI, a scalar between zero and one.}

\item{frac1}{The probability mass of the first interval used for the Geweke test, a scalar
between zero and one.}

\item{frac2}{The probability mass of the second interval used for the Geweke test, a scalar
between zero and one.}
}
\value{
A data frame with the following columns
\describe{
  \item{Mean}{The posterior mean.}
  \item{Median}{The posterior median.}
  \item{SD}{Standard deviation.}
  \item{HPDI-LB}{Highest posterior density credible interval lower bound}
  \item{HPDI-UB}{Highest posterior density credible interval upper bound}
  \item{Naive SE}{Naive Standard error of the mean (ignoring chain autocorrelation.}
  \item{Time-series SE}{Time-series standard error (based on spectral density at 0).}
  \item{Geweke statistic}{The Geweke test statistic.}
  \item{frac1}{The fraction of data contained in the first interval.}
  \item{frac2}{The fraction of data contained in the second interval.}
}
}
\description{
Computes MCMC summary statistics.
}
\details{
Naturally, \code{frac1 + frac2} is between zero and one.
}
\keyword{internal}
