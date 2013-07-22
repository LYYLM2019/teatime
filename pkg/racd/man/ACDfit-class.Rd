\name{ACDfit-class}
\docType{class}
\alias{ACDfit-class}
\alias{coef,ACDfit-method}
\alias{convergence,ACDfit-method}
\alias{fitted,ACDfit-method}
\alias{getspec,ACDfit-method}
\alias{infocriteria,ACDfit-method}
\alias{likelihood,ACDfit-method}
\alias{plot,ACDfit,missing-method}
\alias{quantile,ACDfit-method}
\alias{pit,ACDfit-method}
\alias{residuals,ACDfit-method}
\alias{shape,ACDfit-method}
\alias{shape,ANY-method}
\alias{shape}
\alias{kurtosis,ACDfit-method}
\alias{kurtosis,ANY-method}
\alias{kurtosis}
\alias{show,ACDfit-method}
\alias{sigma,ACDfit-method}
\alias{skew,ACDfit-method}
\alias{skew,ANY-method}
\alias{skew}
\alias{skewness,ACDfit-method}
\alias{skewness,ANY-method}
\alias{skewness}
\title{class: ACD Fit Class}
\description{
S4 class containing the estimated ACD object.
}
\section{Objects from the Class}{
Objects can be created by calling the \code{acdfit} function.
}
\section{Slots}{
  \describe{
    \item{\code{fit}:}{Object of class \code{"vector"} estimated model.}
    \item{\code{model}:}{Object of class \code{"vector"} model details. }
  }
}
\section{Methods}{
  \describe{
    \item{getspec}{\code{signature(object = "ACDfit")}: returns the specification
    used to estimate the model. }
    \item{coef}{\code{signature(object = "ACDfit")}: estimated coefficients. }
    \item{convergence}{\code{signature(object = "ACDfit")}: solver convergence 
    code (0==converged).}
    \item{infocriteria}{\code{signature(object = "ACDfit")}: the information 
    criteria for the ACD and restricted GARCH model.}
    \item{likelihood}{\code{signature(object = "ACDfit")}: log-likelihood of both
    the ACD and restricted GARCH model.}
    \item{show}{\code{signature(object = "ACDfit")}: summary method. }
    \item{plot}{\code{signature(x = "ACDfit", y = "missing")}: plots.}
    \item{fitted}{\code{signature(object = "ACDfit")}: conditional mean. }
    \item{residuals}{\code{signature(object = "ACDfit")}: the model residuals 
    (with logical option standardize set to FALSE as default).}
    \item{sigma}{\code{signature(object = "ACDfit")}: conditional sigma. }
    \item{skew}{\code{signature(object = "ACDfit")}: conditional skew. }
    \item{shape}{\code{signature(object = "ACDfit")}: conditional shape. }
    \item{skewness}{\code{signature(object = "ACDfit")}: conditional skewness. }    
    \item{kurtosis}{\code{signature(object = "ACDfit")}: conditional excess kurtosis.}
    \item{quantile}{\code{signature(x = "ACDfit")}: the conditional quantiles given
    a vector of probabilities. }
    \item{pit}{\code{signature(x = "ACDfit")}: the probability integral 
    transform of the data given the conditional density.}
	 }
}
\references{
Hansen, B. E. 1994, Autoregressive conditional density estimation, 
\emph{International Economic Review}, \bold{35(3)}, 705--730.\cr
Ghalanos, A., Rossi E., and Urga G. 2013, Independent Factor Autoregressive 
Conditional Density Model, \emph{Econometric Reviews}, forthcoming.
}
\author{Alexios Ghalanos}
\note{
If the \code{skew0} and \code{shape0} were both passed to the \code{\link{acdfit}} 
function, then no GARCH model will be estimated to obtain starting parameters, 
meaning that the \sQuote{infocriteria} and \sQuote{likelihood} methods will not 
have this information for comparison.\cr
The \sQuote{skew} and \sQuote{shape} methods return the time varying skew and 
shape parameters, and have an additional logical option \code{transformed} (set 
to TRUE) which returns the logistic transformed parameters.\cr
The \sQuote{quantile} method takes an additional vector \code{probs} as in the 
S3 stats package method.\cr
The \sQuote{residuals} method takes an additional option \code{standardize} 
(default is FALSE) which returns the standardized version of the same.
}
\examples{
\dontrun{
data(sp500ret)
spec = acdspec(variance.model=list(variance.targeting = TRUE),
mean.model=list(armaOrder=c(1,1)),distribution.model=list(model = "jsu",
skewOrder=c(1,1,0), shapeOrder=c(1,1,0)))
fit = acdfit(spec, sp500ret)
head(kurtosis(fit))
head(skewness(fit))
head(sigma(fit))
head(quantile(fit, probs=c(0.01)))
}
}
\keyword{classes}
