% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gap.R
\name{is.gap}
\alias{is.gap}
\title{\code{gap} object check}
\usage{
is.gap(object, return.logical = FALSE)
}
\arguments{
\item{object}{An object to be tested.}

\item{return.logical}{If \code{return.logical = FALSE} (default), an error message is printed
if the object is not of class \code{gap}. If \code{return.logical = TRUE}, a logical
value is returned.}
}
\value{
A logical value or nothing, depending on the value of \code{return.logical}.
}
\description{
Tests whether the input object is a valid object of class \code{gap}.
}
