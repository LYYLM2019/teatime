\name{acdfit-methods}
\docType{methods}
\alias{acdfit-methods}
\alias{acdfit}
\alias{acdfit,ACDspec-method}
\alias{acdfit,ANY-method}
\title{ ACD Model Estimation }
\description{
Estimation of ACD models.
}
\usage{
acdfit(spec, data, solver = "ucminf", out.sample = 0, solver.control = list(), 
fit.control = list(stationarity = 0, fixed.se = 0, scale = 0, n.sim = 2000), 
skew0 = NULL, shape0 = NULL, cluster = NULL, ...)
}
\arguments{
\item{data}{
A univariate xts data object (or one which can be coerced to such).}
\item{spec}{
A univariate ACD spec object of class \code{\linkS4class{ACDspec}}.}
\item{out.sample}{
A positive integer indicating the number of periods before the last to keep for 
out of sample forecasting (see details).}
\item{solver}{
One of either \dQuote{nlminb}, \dQuote{solnp}, \dQuote{optim},
\dQuote{ucminf} or \dQuote{cmaes} or their multistart equivalents (see details).
}
\item{solver.control}{
Control arguments list passed to optimizer.}
\item{fit.control}{
Control arguments passed to the fitting routine. Stationarity explicitly imposes 
the variance stationarity constraint during optimization. The fixed.se argument 
controls whether standard errors should be calculated for those parameters which 
were fixed (through the fixed.pars argument of the \code{\link{acdspec}}
function). The scale parameter controls whether the data should be scaled before 
being submitted to the optimizer. The n.sim option determines the number of
random parameter sets to generate and evaluate as starting values (see details).}
\item{skew0}{ Optional recursion starting parameter for the skew dynamics. If not
used, a restricted GARCH model is first estimated in order to obtain this. }
\item{shape0}{ Optional recursion starting parameter for the shape dynamics. If 
not used, a restricted GARCH model is first estimated in order to obtain this. }
\item{cluster}{
A pre-created cluster from the parallel package for use with the multistart
solvers.}
\item{...}{ In the case of the mcsGARCH model, a \sQuote{DailyVar}  xts object
needs to be passed of the forecast daily variance for the period under 
consideration.}
}
\value{
A \code{\linkS4class{ACDfit}} object containing details of the ACD fit.
}
\details{
For the optim and ucminf solvers (unconstrained optimization), the parameters
are transformed from the unbounded to the bounded domain (upper and lower bounds)
using a logistic transformation. This approach seems to provide better results for
most cases and enough data. For the optim solvers, the default used is \dQuote{BFGS},
but you can choose another option by passing the \sQuote{method} in solver.control
options in addition to any other solver specific options. It is suggested that 
the n.sim option is set high enough to allow for a good set of starting parameters
to be generated, where the Rsolnp package's \sQuote{startpars} function is used 
for this purpose. A sure indication that the solver has converged to only a local 
solution is a warning about failure to invert the hessian or that the likelihood
is lower than that of the restricted GARCH model. In that case, try with a 
different solver, higher n.sim etc. A visual inspection of the skewness
and kurtosis plots may also reveal problems. Setting tighter or wider bounds on
either the complete skew/shape dynamics or the individual parameters driving those
is also an option using the \code{setbounds<-} method on the specification.
The 4 main solvers have multistart versions (preceded by the letters \dQuote{ms}, e.g.
\dQuote{msoptim}) which take an extra parameter in the solver.control list of
\sQuote{restarts} denoting the number of independent optimizations to perform from
different starting parameters. When used with a cluster based object this leads
to a parallel evaluation of different starting parameters and a much better 
chance of obtaining a good solution. The \dQuote{cmaes} solver has a very large
and specific set of options which the user can consult.
}
\author{Alexios Ghalanos}
\examples{
\dontrun{
data(sp500ret)
spec = acdspec(variance.model=list(variance.targeting = TRUE),
mean.model=list(armaOrder=c(1,1)),distribution.model=list(model = "jsu",
skewOrder=c(1,1,0), shapeOrder=c(1,1,0)))
fit = acdfit(spec, sp500ret, solver="msoptim",solver.control=list(restarts=5))
convergence(fit)
head(residuals(fit))
head(fitted(fit))
head(sigma(fit))
head(skewness(fit))
head(kurtosis(fit))
head(quantile(fit, probs=c(0.01, 0.05)))
}
}
\keyword{methods}