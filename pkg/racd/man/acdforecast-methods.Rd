\name{acdforecast-methods}
\docType{methods}
\alias{acdforecast-methods}
\alias{acdforecast}
\alias{acdforecast,ACDfit-method}
\alias{acdforecast,ACDspec-method}
\alias{acdforecast,ANY-method}
\title{ ACD Model Forecast }
\description{
Forecasting of ACD models.
}
\usage{
acdforecast(fitORspec, data = NULL, n.ahead = 10, n.roll = 0,
out.sample = 0, external.forecasts = list(mregfor = NULL, 
vregfor = NULL, skregfor = NULL, shregfor = NULL), 
m.sim = 1000, cluster = NULL, skew0 = NULL, shape0 = NULL, ...)
}
\arguments{
\item{fitORspec}{
A univariate ACD estimated object of class \code{\linkS4class{ACDfit}}, or
an specification object of class \code{\linkS4class{ACDspec}} with fixed 
parameters}
\item{data}{An xts dataset in the case when a specification is used.}
\item{n.ahead}{
The forecast horizon.}
\item{n.roll}{
The no. of rolling forecasts to create beyond the first one (see details).}
\item{out.sample}{
Used when a specification object with data is supplied.
}
\item{external.forecasts}{
A list with forecasts for the (optional) external regressors in the mean, 
variance, skew and shape dynamics equations if specified.}
\item{m.sim}{
For n.ahead>1, forecast of the higher moments is done via simulation, in which
case m.sim indicates the number of simulations per forecast to use for generating
the average forecast value.
}
\item{cluster}{
A pre-created cluster from the parallel package to use when simulation methods 
are called (for n.ahead>1).
}
\item{skew0}{Skew dynamics recursion initialization value in the case of 
a spec being supplied.}
\item{shape0}{Shape dynamics recursion initialization value in the case of 
a spec being supplied.}
\item{...}{For the multiplicative component sGARCH model (mcsGARCH), the
additional argument \sQuote{DailyVar} is required and should be an xts object of 
the daily forecasted variance for the period under consideration to be used 
with the intraday data.}
}
\value{
A \code{\linkS4class{ACDforecast}} object containing details of the ACD forecast.
}
\details{
The ability to roll the forecast 1 step at a time is implemented with the 
\code{n.roll} argument which controls how many times to roll the n.ahead 
forecast. The default argument of n.roll = 0 denotes no rolling and returns the 
standard n.ahead forecast. Critically, since n.roll depends on data being 
available from which to base the rolling forecast, the \code{\link{acdfit}} 
function needs to be called with the argument \code{out.sample} being at least 
as large as the n.roll argument, or in the case of a specification being used 
instead of a fit object, the \code{out.sample} argument directly in the forecast 
function.
}
\author{Alexios Ghalanos}
\keyword{methods}