% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checks.R
\name{.checkModelPrior}
\alias{.checkModelPrior}
\title{Checks whether model and prior match.}
\usage{
.checkModelPrior(model, prior)
}
\arguments{
\item{model}{An object of class TFPmodel.}

\item{prior}{A list of matrices with parameters for the prior distribution and box
constraints. By default, \code{prior} is initialized by \code{initializePrior(model)}.
See details. Only used if \code{method = "bayesian"}.}
}
\description{
Checks whether model and prior match.
}
\keyword{internal}
