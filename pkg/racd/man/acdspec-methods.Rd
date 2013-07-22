\name{acdspec-methods}
\docType{methods}
\alias{acdspec}
\alias{acdspec-methods}
\alias{acdspec,ANY-method}
\title{function: ACD Specification}
\description{
Method for creating a univariate ACD specification object prior to fitting.
}
\usage{
acdspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1), 
external.regressors = NULL, variance.targeting = FALSE), 
mean.model = list(armaOrder = c(1, 1), include.mean = TRUE, archm = FALSE, 
arfima = FALSE, external.regressors = NULL), 
distribution.model = list(model = "snorm", skewOrder = c(1, 1, 1), skewshock = 1, 
skewshocktype = 1, skewmodel = "quad", skew.regressors = NULL, 
shapeOrder = c(0, 1, 1), shapeshock = 1, shapeshocktype = 1, shapemodel = "quad", 
shape.regressors = NULL, exp.rate = 1), start.pars = list(), fixed.pars = list())
}
\arguments{
\item{variance.model}{
List containing the variance model specification:\cr
\code{model} Valid models (currently implemented) are \dQuote{sGARCH}, 
\dQuote{csGARCH}, \dQuote{mcsGARCH}.\cr
\code{garchOrder} The ARCH (q) and GARCH (p) orders.\cr
\code{external.regressors} A matrix object containing the external regressors to 
include in the variance equation with as many rows as will be included in the
data (which is passed in the fit function).
\code{variance.targeting} (Logical) If logical, indicates whether to use variance 
targeting for the conditional variance intercept \dQuote{omega}. 
}
\item{mean.model}{
List containing the mean model specification:\cr
\code{armaOrder} The autoregressive (ar) and moving average (ma) orders (if any).\cr
\code{include.mean} Whether to include the mean.\cr
\code{archm} Whether to include ARCH volatility in the mean regression.\cr
\code{arfima} Whether to fractional differencing in the ARMA regression.\cr
\code{external.regressors} A matrix object containing the external regressors to 
include in the mean equation with as many rows as will be included in the data 
(which is passed in the fit function).\cr
}
\item{distribution.model}{
List containing the conditional distribution dynamics specification:\cr
\code{model} The conditional density to use for the innovations. Valid choices 
are \dQuote{snorm} for the skew-normal distribution, \dQuote{std} for the 
student-t, \dQuote{sstd} for the skew-student, \dQuote{ged} for the generalized 
error distribution, \dQuote{sged} for the skew-generalized error distribution, 
\dQuote{nig} for the normal inverse gaussian distribution, \dQuote{ghyp} for the 
Generalized  Hyperbolic, and \dQuote{jsu} for Johnson's SU distribution.\cr
\code{skewOrder} The skew dynamics order (see vignette for details).\cr
\code{skewshock} Whether to use the residuals (2) or standardized residuals (1) 
to model the skew dynamics.\cr
\code{skewshocktype} Whether to model the skewshock dynamics using the squared (1)
or absolute (any other value) value function.\cr
\code{skewmodel} A choice of skew dynamic models including quadratic (\dQuote{quad}), 
piece-wise linear (\sQuote{pwl}), general threshold based (\sQuote{tar}) or 
excess shock over the conditional absolute sigma (\sQuote{xar}).\cr
\code{skew.regressors} Optional matrix of regressors to include in the conditional
skew dynamics.\cr
\code{shapeOrder} The shape dynamics order (see vignette for details).\cr
\code{shapeshock} Whether to use the residuals (2) or standardized residuals (1) 
to model the shape dynamics.\cr
\code{shapeshocktype} Whether to model the shapeshock dynamics using the squared (1)
or absolute (any other value) value function.\cr
\code{shapemodel} A choice of shape dynamic models including quadratic (\dQuote{quad}), 
piece-wise linear (\sQuote{pwl}), general threshold based (\sQuote{tar}) or 
excess shock over the conditional absolute sigma (\sQuote{xar}).\cr
\code{shape.regressors} Optional matrix of regressors to include in the conditional
shape dynamics.\cr
\code{exp.rate} The rate value for the exponential transformation used for the
shape dynamics (the skew dynamics uses the logistic transformation without 
additional parameters).
}
\item{start.pars}{
List of staring parameters for the optimization routine. These are not usually 
required unless the optimization has problems converging.}
\item{fixed.pars}{
List of parameters which are to be kept fixed during the optimization. It is 
possible that you designate all parameters as fixed so as to quickly recover 
just the results of some previous work or published work. The optional argument
\dQuote{fixed.se} in the \code{\link{acdfit}} function indicates whether to
calculate standard errors for those parameters fixed during the post 
optimization stage.}
}
\value{
A \code{\linkS4class{ACDspec}} object containing details of the ACD specification.
}
\details{
The ACD specification allows for a number of different parameterizations of the
conditional higher moment dynamics, with a host of specialized options. Not all
combinations of options and dynamics are available and the user should consult
the vignette where the models are more fully discussed.
}
\author{Alexios Ghalanos}
\keyword{methods}