\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot}
\alias{plot,ANY,ANY-method}
\alias{plot,GPDFIT,missing-method}
\alias{plot,GPDTAILS,missing-method}
\title{Method: Plotting (S4) for implemented S4 classes}
\description{
Locally implemented and modified methods for plotting the fit of the 
\code{\linkS4class{GPDFIT}} object (taken from package \code{fExtremes}), and 
the overall fit of the \code{\linkS4class{GPDTAILS}} object.
}
\usage{
plot(x,y,...)
}
\arguments{
	\item{x}{ An object of class \code{\linkS4class{GPDFIT}} or \code{\linkS4class{GPDTAILS}}.}
	\item{y}{ missing }
	\item{...}{
	[which] - \cr
	Valid arguments are \dQuote{ask} for interactive plotting, \dQuote{all} to 
	plot all or a number from 1 to 6 indicating the desired charts to plot.}
}
\examples{\dontrun{
library(MASS)
x<-SP500/100
fit<-spdfit(x)
plot(fit,which=1)
# this in fact exctracts the GPDFIT object (from GPDTAILS) for which plot 
# methods exist.
plot(fit,which=3)
}}
\keyword{methods}