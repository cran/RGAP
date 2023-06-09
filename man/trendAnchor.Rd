% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/state_space.R
\name{trendAnchor}
\alias{trendAnchor}
\title{Trend anchor}
\usage{
trendAnchor(fit, anchor = NULL, h = NULL, returnFit = FALSE)
}
\arguments{
\item{fit}{An object of class \code{NAWRUfit}, \code{TFPfit}, or \code{KuttnerFit}.}

\item{anchor}{A numeric specifying the anchor value. If unspecified, \code{anchor} is
taken from the object \code{fit} (if specified).}

\item{h}{An integer specifying the anchor horizon in the frequency of the underlying model.
If unspecified, \code{h} is taken from the object \code{fit} (if specified).}

\item{returnFit}{A logical. If \code{TRUE}, an object of the same class as \code{fit}
including the anchored trend is returned. If \code{FALSE}, only the anchored trend time
series is returned.}
}
\value{
A fitted object if \code{returnFit = TRUE} or a time series with the anchored 
  trend.
}
\description{
Computes the anchored trend given a fitted object of class \code{NAWRUfit},
  \code{TFPfit}, or \code{KuttnerFit}.
}
\examples{
# define nawru model for France
data("gap")
tsList <- amecoData2input(gap$France)
model <- NAWRUmodel(tsl = tsList)

# estimate nawru model
\donttest{
f <- fit(model = model)

# compute anchored nawru
anchoredNawru <- trendAnchor(fit = f, anchor = 6.5, h = 10)
}
}
