% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gap.R
\name{gapHP}
\alias{gapHP}
\title{HP-filter output gap}
\usage{
gapHP(x, lambda = NULL, end = NULL, start = NULL)
}
\arguments{
\item{x}{A time series object containing gdp.}

\item{lambda}{The smoothing parameter for the application of the HP filter. If not supplied,
\code{lambda = 6.25} for yearly data, \code{lambda = 1600} for quarterly data, and \code{lambda = 129600}
for monthly data.}

\item{end}{(optional) A two element vector containing a year and a period specifying the end point for the
filter application.}

\item{start}{(optional) A two element vector containing a year and a period specifying the start point for the
filter application.}
}
\value{
Object of class \code{gap}, which is a list containing the two elements \code{potential} and
  \code{gap} and additionally the original time series.
}
\description{
Computes a HP filtered output gap.
}
