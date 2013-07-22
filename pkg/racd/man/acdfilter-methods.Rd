\name{acdfilter-methods}
\docType{methods}
\alias{acdfilter-methods}
\alias{acdfilter}
\alias{acdfilter,ACDspec-method}
\alias{acdfilter,ANY-method}
\title{ ACD Model Estimation }
\description{
Filtering of data with ACD dynamics.
}
\usage{
acdfilter(spec, data, out.sample = 0,  n.old = NULL, skew0 = NULL, 
shape0 = NULL, ...)
}
\arguments{
\item{data}{
A univariate xts data object (or one which can be coerced to such).}
\item{spec}{
A univariate ACD spec object of class \code{\linkS4class{ACDspec}}.}
\item{out.sample}{
A positive integer indicating the number of periods before the last to keep for 
out of sample forecasting.}
\item{n.old}{ For comparison with ACDfit models using the out.sample argument, 
this is the length of the original dataset.}
\item{skew0}{ Optional recursion starting parameter for the skew dynamics. If not
used, the transformed skew dynamics intercept value is used. }
\item{shape0}{ Optional recursion starting parameter for the shape dynamics. If not
used, the transformed shape dynamics intercept value is used. }
\item{...}{ . }
}
\value{
A \code{\linkS4class{ACDfilter}} object containing details of the ACD filter.
}
\author{Alexios Ghalanos}
\examples{
\dontrun{
# Check that fit/filter return the same values
library(racd)
library(rugarch)
data(sp500ret)
spec = acdspec(variance.model=list(variance.targeting = TRUE),
mean.model=list(armaOrder=c(1,1)),distribution.model=list(model = "jsu",
skewOrder=c(1,1,0), shapeOrder=c(1,1,0)))
fit = acdfit(spec, sp500ret)
# remove variance targeting:
spec = acdspec(variance.model=list(variance.targeting = FALSE),
mean.model=list(armaOrder=c(1,1)),distribution.model=list(model = "jsu",
skewOrder=c(1,1,0), shapeOrder=c(1,1,0)))
setfixed(spec)<-as.list(coef(fit))
# fit@fit$skhEst contains the untransformed recursion starting values for the skew and shape dynamics.
# fit@model$sbounds contain the: [skew.LB, skew.UB, shape.LB, shape.UB, shape.rate(for exponential transformation)]
# skew0 and shape0 take the transformed values (so use the functions logtransform and exptransform)
filt = acdfilter(spec, sp500ret, skew0 = logtransform(fit@fit$skhEst[1], fit@model$sbounds[1],fit@model$sbounds[2]), 
shape0 = exptransform(fit@fit$skhEst[2], fit@model$sbounds[3],fit@model$sbounds[4], rate = fit@model$sbounds[5]))
head(cbind(kurtosis(fit), kurtosis(filt)))
head(cbind(skewness(fit), skewness(filt)))
head(cbind(sigma(fit), sigma(filt)))
head(cbind(quantile(fit, probs=c(0.01)), quantile(filt, probs=c(0.01))))
}
}
\keyword{methods}