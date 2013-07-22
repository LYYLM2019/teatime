\name{acdsim-methods}
\docType{methods}
\alias{acdsim}
\alias{acdsim-methods}
\alias{acdsim,ANY-method}
\alias{acdsim,ACDfit-method}
\title{function: ACD Simulation}
\description{
Method for simulation from a variety of univariate ACD models.
}
\usage{
acdsim(fit, n.sim = 1000, n.start = 0, m.sim = 1, 
presigma = NA, prereturns = NA, preresiduals = NA, preskew = NA, 
preshape = NA, rseed = NA, mexsimdata = NULL, vexsimdata = NULL, 
skxsimdata = NULL, shxsimdata = NULL, ...)
}
\arguments{
\item{fit}{
A ACD fit object of class \code{\linkS4class{ACDfit}}.}
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
\item{...}{ For the multiplicative component sGARCH model (mcsGARCH), the
additional argument \sQuote{DailyVar} is required and should be an xts object 
of length floor(n.sim/increments-per-day) by m.sim of the the daily simulated 
variance to use with the intraday data.}
}
\value{
A \code{\linkS4class{ACDsim}} object containing details of the ACD simulation.
}
\details{
Unlike GARCH models, the conditional standardized innovations cannot be 
pre-generated but are dependent on the value of the conditional simulated skew 
and shape values at each point in time. This means that simulation is somewhat
slower and there are some additional issues discussed in the vignette,
}
\author{Alexios Ghalanos}
\keyword{methods}
