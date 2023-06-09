% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nawru.R
\name{NAWRUmodel}
\alias{NAWRUmodel}
\title{NAWRU model}
\usage{
NAWRUmodel(
  tsl,
  trend = "RW2",
  cycle = "AR2",
  type = "TKP",
  cycleLag = 0,
  pcErrorARMA = c(0, 0),
  exoType = NULL,
  start = NULL,
  end = NULL,
  anchor = NULL,
  anchor.h = NULL
)
}
\arguments{
\item{tsl}{A list of time series objects, see details.}

\item{trend}{A character string specifying the trend model. \code{trend = "RW1"} denotes
a first order random walk, \code{trend = "RW2"} a second order random walk (local linear
trend) and \code{trend = "DT"} a damped trend model. The default is \code{trend = "RW2"}.}

\item{cycle}{A character string specifying the cycle model. \code{cycle = "AR1"} denotes
an AR(1) process, \code{cycle = "AR2"} an \code{AR(2)} process. The default is
\code{cycle = "AR2"}.}

\item{type}{A character string specifying the type of the Phillip's curve.
\code{type = "TKP"} denotes the traditional Keynesian Phillip's curve and
\code{type = "NKP"} the New Keynesian Phillip's curve, see details. The default is
\code{type = "TKP"}.}

\item{cycleLag}{A vector specifying the cycle lags that are included in the Phillip's
curve. The default is \code{cycleLag = 0}, see details.}

\item{pcErrorARMA}{A \code{2 x 1} vector with non-negative integers specifying the AR
and MA degree of the error term in the Phillip's curve equation. The default is
\code{pcErrorARMA = c(0, 0)}, see details.}

\item{exoType}{An optional \code{n x m x 2} array specifying the possible difference
and lag transformation for the variables. \code{exoType} can be initialized using the
function \code{inizializeExo}. The column names give the variable names.
\code{exoType[, , 1]} contains the difference transformations and \code{exoType[, , 2]}
the subsequent lag transformations, see details.}

\item{start}{(Optional) Start vector for the estimation, e.g. \code{c(1980, 1)}.}

\item{end}{(Optional) End vector for the estimation, e.g. \code{c(2020, 1)}.}

\item{anchor}{(Optional) Anchor value for the unemployment rate.}

\item{anchor.h}{(Optional) Anchor horizon in the frequency of the given time series.}
}
\value{
Object of class \code{NAWRUmodel}, which is a list with the following components:
  \item{tsl}{A list of used time series.}
  \item{SSModel}{An object of class SSModel specifying the state-space model.}
  \item{loc}{A data frame containing information on each involved parameter, for instance
        its corresponding system matrix, variable names, and parameter restrictions.}
  \item{call}{Original call to the function.}
  In addition, the object contains the following attributes:
  \item{cycle}{Cycle specification.}
  \item{trend}{Trend specification.}
  \item{phillipsCurve}{A list containing the components \code{type, cycleLag, errorARMA, exoVariables}.}
  \item{anchor}{A list containing the components \code{value, horizon}.}
  \item{period}{A list containing the components \code{start, end, frequency}.}
}
\description{
Creates a state space object object of class \code{NAWRUmodel} which can be 
  fitted using \code{fit}.
}
\details{
The list of time series \code{tsl} needs to have the following components:
\describe{
    \item{ur}{Unemployment rate.}
    \item{nulc}{Nominal Unit labor costs, if \code{type = "TKP"}.}
    \item{rulc}{Real unit labor costs, if \code{type = "NKP"}.}
    }
  and optionally other variables included in \code{exoType}.

A \code{cycleLag} equal to \code{0} implies that only the contemporaneous cycle
  is included in the Phillip's curve.  A \code{cycleLag} equal to \code{0:1} implies that
  the contemporaneous as well as the lagged cycle are included.

A \code{pcErrorARMA} equal to \code{c(0, 0)} implies that the error term in the
  Phillip's curve is white noise. \code{pcErrorARMA = c(1, 0)} implies that the error is
  an AR(1) process and for \code{pcErrorARMA = c(1, 2)} the error follows an ARMA(1, 2)
  process.

For the New Keynesian Phillip's curve, the \code{cycleLag} cannot be chosen.
  \code{cycleLag} will be set to \code{0} if \code{cycle = "AR1"} and to \code{1} if
  \code{cycle = "AR2"}. In the latter case, the forward solution of the Phillip's curve
  implies parameter restrictions for the lagged cycle on the Phillip's curve. Moreover,
  exogenous variables will be ignored in the case of the New Keynesian Phillip's curve.

The array \code{exoType} consists of non-negative integers or \code{NA}s.
  \code{exoType[, , 1] = c(NA,1)} and \code{exoType[, , 2] = c(NA,2)} implies that
  the first variable is not included in the Phillip's curve whereas the second lag of
  the first difference of the second variable is included.
}
\examples{
# load data for France
data("gap")
tsList <- amecoData2input(gap$France, alpha = 0.65)
# Traditional phillips curve
model <- NAWRUmodel(tsl = tsList, trend = "RW2", cycle = "AR2", type = "TKP", cycleLag = 0)

# New-Keynesian Phillips curve
model <- NAWRUmodel(tsl = tsList, trend = "RW2", cycle = "AR2", type = "NKP", cycleLag = 0:1)

# Traditional Phillips curve with 6 exogenous variables
# specify exogenous variable transformations
D <- matrix(c(2, 2, 2, 1, 1, 1), 2, 3, byrow = TRUE)
L <- matrix(c(0, 0, 0, 1, 1, 1), 2, 3, byrow = TRUE)
exoType <- initializeExo(varNames = c("tot", "prod","ws"), D = D, L = L)
model <- NAWRUmodel(tsl = tsList, cycleLag = 0:1, exoType = exoType)
}
