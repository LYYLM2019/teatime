\name{acdpath-methods}
\docType{methods}
\alias{acdpath}
\alias{acdpath-methods}
\alias{acdpath,ANY-method}
\alias{acdpath,ACDspec-method}
\title{function: ACD Path Simulation}
\description{
Method for simulating the path of an ACD model from a variety of univariate 
GARCH models. This is a convenience function which does not require a fitted 
object (see note below).
}
\usage{
acdpath(spec, n.sim = 1000, n.start = 0, m.sim = 1, presigma = NA, 
prereturns = NA, preresiduals = NA, preskew = NA, preshape = NA, rseed = NA,  
mexsimdata = NULL, vexsimdata = NULL, skxsimdata = NULL, shxsimdata = NULL, 
cluster = NULL, ...)
}
\arguments{
\item{spec}{
A univariate ACD spec object of class \code{\linkS4class{ACDspec}} with
the required parameters of the model supplied via the fixed.pars list argument 
or \code{\link{setfixed<-}} method.}
\item{n.sim}{
The simulation horizon.}
\item{n.start}{
The burn-in sample.}
\item{m.sim}{
The number of simulations.}
\item{presigma}{
Allows the starting sigma values to be provided by the user.}
\item{prereturns}{
Allows the starting return data to be provided by the user.}
\item{preresiduals}{
Allows the starting residuals to be provided by the user.}
\item{preskew}{
Allows the starting skew dynamic's parameter to be provided by the user.
}
\item{preshape}{
Allows the starting shape dynamic's parameter to be provided by the user.
}
\item{rseed}{
Optional seeding value(s) for the random number generator. For m.sim>1, it is 
possible to provide either a single seed to initialize all values, or one seed 
per separate simulation (i.e. m.sim seeds). However, in the latter case this may 
result in some slight overhead depending on how large m.sim is.}
\item{mexsimdata}{
List of matrices (size of list m.sim, with each matrix having n.sim rows) of 
simulated external regressor-in-mean data. If the fit object contains 
external regressors in the mean equation, this must be provided else will be
assumed zero.}
\item{vexsimdata}{
List of matrices (size of list m.sim, with each matrix having n.sim rows) of 
simulated external regressor-in-variance data. If the fit object contains 
external regressors in the mean equation, this must be provided else will be
assumed zero.}
\item{skxsimdata}{
List of matrices (size of list m.sim, with each matrix having n.sim rows) of 
simulated external regressor-in-skew data. If the fit object contains 
external regressors in the skew dynamics, this must be provided else will be
assumed zero.}
\item{shxsimdata}{
List of matrices (size of list m.sim, with each matrix having n.sim rows) of 
simulated external regressor-in-shape data. If the fit object contains 
external regressors in the shape dynamics, this must be provided else will be
assumed zero.}
\item{cluster}{
Pre-created cluster from the parallel package for parallel evaluation of 
independent (m.sim) streams.
}
\item{...}{ If the model is the \dQuote{csGARCH}, then \code{preq} can be used 
to denote the previous value of the permanent component of the variance model 
(q, e.g. tail(fit@fit$q,1)) so that the acdpath method with all pre-values 
included will evaluate to the same result as the acdsim method with method 
equal to \dQuote{sample} (assuming the same random seeding values are used).}
}
\value{
A \code{\linkS4class{ACDpath}} object containing details of the ACD path 
simulation.
}
\details{
This is a convenience method to allow path simulation of various ACD models 
without the need to supply a fit object as in the \code{\link{acdsim}} method. 
Instead, a ACD spec object is required with the fixed model parameters.
The mcsGARCH model is not supported for the path method-use \code{\link{acdsim}}
instead.
}
\author{Alexios Ghalanos}
\keyword{methods}