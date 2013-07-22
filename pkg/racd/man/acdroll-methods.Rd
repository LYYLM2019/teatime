\name{acdroll-methods}
\docType{methods}
\alias{acdroll}
\alias{acdroll-methods}
\alias{acdroll,ANY-method}
\alias{acdroll,ACDspec-method}
\title{function: ACD Rolling Density Forecast and Backtesting}
\description{
Method for creating rolling density forecast from ACD models with option
for refitting every n periods with parallel functionality.
}
\usage{
acdroll(spec, data, n.ahead = 1, forecast.length = 500, 
n.start = NULL, refit.every = 25, refit.window = c("recursive", "moving"), 
window.size = NULL, solver = "ucminf", fit.control = list(), 
solver.control = list(), calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05), 
cluster = NULL, keep.coef = TRUE, fixARMA = TRUE, fixGARCH = TRUE, 
fixUBShape = TRUE, UBShapeAdd = 0, fixGHlambda = TRUE, 
compareGARCH=c("LL", "none"),...)
}
\arguments{
\item{spec}{
An ACD specification object.}
\item{data}{
A univariate xts dataset.
}
\item{n.ahead}{
The number of periods to forecast (only n.ahead=1 supported).
}
\item{forecast.length}{
The length of the total forecast for which out of sample data 
from the dataset will be used for testing.
}
\item{n.start}{
Instead of forecast.length, this determines the starting point in the dataset
from which to initialize the rolling forecast.
}
\item{refit.every}{
Determines every how many periods the model is re-estimated.
}
\item{refit.window}{
Whether the refit is done on an expanding window including all the previous data 
or a moving window where all previous data is used for the first estimation 
and then moved by a length equal to refit.every (unless the window.size option 
is used instead).
}
\item{window.size}{
If not NULL, determines the size of the moving window in the rolling estimation,
which also determines the first point used.
}
\item{solver}{
The solver to use.
}
\item{fit.control}{
Control parameters parameters passed to the fitting function.
}
\item{solver.control}{
Control parameters passed to the solver.
}
\item{calculate.VaR}{
Whether to calculate forecast Value at Risk during the estimation.
}
\item{VaR.alpha}{
The Value at Risk tail level to calculate.
}
\item{cluster}{
A cluster object created by calling \code{makeCluster} from the parallel 
package. If it is not NULL, then this will be used for parallel estimation
of the refits (remember to stop the cluster on completion).
}
\item{keep.coef}{
Whether to return the list of coefficients and their robust standard errors.
}
\item{fixARMA}{
}
\item{fixGARCH}{
}
\item{fixUBShape}{
}
\item{UBShapeAdd}{
}
\item{fixGHlambda}{
}
\item{compareGARCH}{
}
\item{...}{ . }
}
\details{
This is a wrapper function for creating rolling forecasts from the ACD model, 
and optionally calculating the Value at Risk at specified levels.
The argument \code{refit.every} determines every how many periods the model is
re-estimated. Given a dataset of length N, it is possible to choose either
how many periods from the end to use for out of sample forecasting (using the
forecast.length option), or the starting point for initializing the rolling
forecast (and using all the data after that for the out of sample forecast).
Only rolling 1-ahead forecasts are supported spanning the dataset, which should
be useful for backtesting models. Anything more complicated should be wrapped by 
the user by making use of the underlying functions in the package.
In case of no-convergence in some of all the windows, a new method called 
\code{\link{resume}} now allows to pass the returned (non-converged) object with
new solver and control parameters to be re-estimated (only the non-converged
windows are re-estimated). Non-convergence here implies both a failure of the
solver to converge to a solution (global failure) OR a failure to invert the 
resulting Hessian (local failure). The \code{convergence} method can be used
on an object (aside from the printed warning) to print out the number of the
non-converged estimation windows.\cr
Parallel functionality is now based entirely on the parallel package, and it is
up to the user to pass a cluster object, and then stop it once the routine is
completed.
}
\value{
An object of class \code{\linkS4class{ACDroll}}.
}
\author{Alexios Ghalanos}
\keyword{methods}