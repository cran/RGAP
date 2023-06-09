% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kuttner_mle.R
\name{fit.KuttnerModel}
\alias{fit.KuttnerModel}
\title{Maximum likelihood estimation of a \code{KuttnerModel}}
\usage{
\method{fit}{KuttnerModel}(
  model,
  parRestr = initializeRestr(model),
  signalToNoise = NULL,
  control = NULL,
  ...
)
}
\arguments{
\item{model}{An object of class KuttnerModel.}

\item{parRestr}{A list of matrices containing the parameter restrictions for the cycle,
trend, and the inflation equation. Each matrix contains the lower and upper bound of the
involved parameters. \code{NA} implies that no restriction is present. Autoregressive
parameters are automatically restricted to the stationary region unless box constraints
are specified. By default, \code{parRestr} is intitialized by the function
\code{initializeRestr(model)}.}

\item{signalToNoise}{(Optional) signal to noise ratio.}

\item{control}{(Optional) A list of control arguments to be passed on to \code{optim}.}

\item{...}{additional arguments to be passed to the methods functions.}
}
\value{
An object of class \code{KuttnerFit} containing the following components:
  \item{model}{The input object of class \code{KuttnerModel}.}
  \item{SSMfit}{The estimation output from the funtcion \code{fitSSM} from \code{KFAS}.}
  \item{SSMout}{The filtering and smoothing output from the funtcion \code{KFS} from
        \code{KFAS}.}
  \item{parameters}{A data frame containing the estimated parameters, including standard
        errors, t-statistics, and p-values.}
  \item{fit}{A list of model fit criteria (see below).}
  \item{call}{Original call to the function. }
  The list component \code{fit} contains the following model fit criteria:
  \item{loglik}{Log-likelihood function value.}
  \item{AIC}{Akaike information criterion.}
  \item{BIC}{Bayesian information criterion.}
  \item{AIC}{Hannan-Quinn information criterion.}
  \item{RMSE}{Root mean squared error of the inflation equation.}
  \item{R2}{R squared of the inflation equation.}
  \item{LjungBox}{Ljung-Box test output of the inflation equation.}
}
\description{
Estimates a two-dimensional state-space model and performs filtering and
  smoothing to obtain the output gap.
}
\examples{
# load data for the Netherlands
data("gap")
country <- "Netherlands"
tsList <- as.list(gap[[country]][, c("cpih", "gdp")])
tsList$infl <- diff(tsList$cpih)
model <- KuttnerModel(tsl = tsList, trend = "RW2", cycleLag = 1, cycle = "AR2", start = 1980)

# estimate Kutter's model
parRestr <- initializeRestr(model = model, type = "hp")
\donttest{
gapKuttner <- fit(model, parRestr, signalToNoise = 1 / 10)
}
}
\seealso{
Other fitting methods: 
\code{\link{fit.NAWRUmodel}()},
\code{\link{fit.TFPmodel}()},
\code{\link{fit}()}
}
\concept{fitting methods}
