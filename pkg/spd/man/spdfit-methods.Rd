\name{spdfit-methods}
\docType{methods}
\alias{spdfit-methods}
\alias{spdfit,ANY-method}
\alias{spdfit}
\title{Method: Fitting the Semi-Parametric Distribution}
\description{
The semi-parametric distribution fitting method.}
\usage{
spdfit(data, upper = 0.9, lower = 0.1, tailfit="GPD", type = c("mle", "pwm"), 
kernelfit = c("normal", "box", "epanech", "biweight", "triweight"), 
information = c("observed", "expected"), title = NULL, description = NULL, ...)
}
\arguments{
  \item{data}{An object coercible to a \code{matrix}.}
  \item{upper}{Upper tail cutoff for fitting the generalized pareto or other distribution.}
  \item{lower}{Lower tail cutoff for fitting the generalized pareto or other distribution.}
  \item{tailfit}{Distribution to Use for fitting the tails.}
  \item{type}{A character string selecting the desired estimation method, either "mle" for the maximum likelihood method or
 		"pwm" for the probability weighted moment method. By default, the first will be selected.}
  \item{kernelfit}{Type of kernel to fit to the interior of the distribution.}
  \item{information}{Whether tail distribution standard errors should be calculated with "observed" or "expected" information. 
	This only applies to the maximum likelihood method; for the probability-weighted moments method "expected" information 
	is used if possible.}
  \item{title}{A character string which allows for a project title. }
  \item{description}{A character string which allows for a brief description. }
  \item{\dots}{Control parameters and plot parameters optionally passed to the optimization and/or plot function. 
	Parameters for the optimization function are passed to components of the control argument of optim.}
}
\value{
Returns an object of class \code{\linkS4class{SPD}}.
}
\examples{\dontrun{
library(MASS)
x<-SP500/100
fit<-spdfit(x)
show(fit)
#plot(fit,which="all")
}}
\keyword{methods}


