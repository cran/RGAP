% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constraints.R
\name{.updateParConstraints}
\alias{.updateParConstraints}
\title{Updates the parameter constraints for on object of class \code{NAWRUmodel} or
\code{TFPmodel}.}
\usage{
.updateParConstraints(model, parRestr)
}
\arguments{
\item{model}{A model of class \code{NAWRUmodel} or \code{TFPmodel}.}

\item{parRestr}{A list of matrices containing the parameter restrictions for the cycle,
trend, and the Phillip's curve. Each matrix contains the lower and upper bound of the
involved parameters. \code{NA} implies that no restriction is present. Autoregressive
parameters are automatically restricted to the stationary region unless box constraints
are specified. By default, \code{parRestr} is initialized by the function
\code{initializeRestr(model)}. Only used if \code{method = "MLE"}.}
}
\value{
The same model with updated list item \code{loc}.
}
\description{
Updates the parameter constraints for on object of class \code{NAWRUmodel} or
\code{TFPmodel}.
}
\keyword{internal}
