\name{ACDmultifilter-class}
\docType{class}
\alias{ACDmultifilter-class}
\alias{coef,ACDmultifilter-method}
\alias{fitted,ACDmultifilter-method}
\alias{likelihood,ACDmultifilter-method}
\alias{residuals,ACDmultifilter-method}
\alias{sigma,ACDmultifilter-method}
\alias{skew,ACDmultifilter-method}
\alias{shape,ACDmultifilter-method}
\title{class: Univariate ACD Multiple Filter Class}
\description{Class for the univariate ACD Multiple filter.}
\section{Methods}{
  \describe{
    \item{fitted}{\code{signature(object = "ACDmultifilter")}:
    Extracts the fitted values.}
    \item{residuals}{\code{signature(object = "ACDmultifilter")}: 
    Extracts the residuals. Optional logical argument \code{standardize} 
	(default is FALSE) allows to extract the standardized residuals.}
    \item{sigma}{\code{signature(object = "ACDmultifilter")}:
    Extracts the conditional sigma values.}
    \item{skew}{\code{signature(object = "ACDmultifilter")}:
    Extracts the conditional skew values. }
    \item{shape}{\code{signature(object = "ACDmultifilter")}:
    Extracts the conditional shape values. }
    \item{coef}{\code{signature(object = "ACDmultifilter")}:
    Extracts the coefficients.}
    \item{likelihood}{\code{signature(object = "ACDmultifilter")}:
    Extracts the likelihood. }
	}
}
\author{Alexios Ghalanos}
\seealso{
Classes \code{\linkS4class{ACDmultiforecast}}, \code{\linkS4class{ACDmultifit}} 
and \code{\linkS4class{ACDmultispec}}.
}
\keyword{classes}