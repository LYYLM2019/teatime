\name{plot.dbm}
\alias{vcov.dbm}
\alias{plot.dbm}
\alias{fitted.dbm}
\alias{residuals.dbm}
\alias{coef.dbm}
\alias{likelihood.dbm}
\alias{score.dbm}
\alias{score}
\alias{likelihood.default}
\alias{likelihood}
\alias{deviance.dbm}
\alias{summary.dbm}
\alias{print.dbm}
\title{
dbm Methods
}
\description{
A number of methods for accessing and extracting information from an estimated
dbm object.
}
\usage{
\method{vcov}{dbm}(object, robust = FALSE, ...)
\method{plot}{dbm}(x, ...)
\method{fitted}{dbm}(object,...)
\method{residuals}{dbm}(object,...)
\method{coef}{dbm}(object,...)
\method{deviance}{dbm}(object, null=FALSE,...)
\method{likelihood}{dbm}(object, ...)
\method{likelihood}{default}(object, ...)
\method{summary}{dbm}(object,...)
\method{score}{dbm}(object, pars = NULL, analytic = TRUE, ... )
\method{print}{dbm}(x, digits = max(3L, getOption("digits") - 3L), ...)
}
\arguments{
\item{x}{
An object of class dbm.
}
\item{object}{
An object of class dbm.
}
\item{robust}{
[vcov] Extracts the robust covariance matrix.
}
\item{null}{
[deviance] Whether to return the null deviance.
}
\item{digits}{
No. of digits for priting the coefficient matrix.
}
\item{pars}{
[score] An (optional) named vector of coefficients to be used instead
of the object's esimated parameters.
}
\item{analytic}{
[score] Whether to calculate (if pars is not null, else no calculation is 
necessary) and return the analytic scores else the numerically estimated
values are returned (as calculated by the \code{jacobian} function of the 
numDeriv package).
}
\item{\dots}{
}
}
\details{
The methods provided allow for extracting or calculating useful values from
the estimated object of class \sQuote{dbm}.
}
\author{
Alexios Ghalanos
}