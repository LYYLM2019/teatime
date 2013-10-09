\name{ACDmultifit-class}
\docType{class}
\alias{ACDmultifit-class}
\alias{coef,ACDmultifit-method}
\alias{fitted,ACDmultifit-method}
\alias{likelihood,ACDmultifit-method}
\alias{residuals,ACDmultifit-method}
\alias{sigma,ACDmultifit-method}
\alias{skew,ACDmultifit-method}
\alias{shape,ACDmultifit-method}
\title{class: Univariate ACD Multiple Filter Class}
\description{Class for the univariate ACD Multiple fit.}
\section{Methods}{
  \describe{
    \item{fitted}{\code{signature(object = "ACDmultifit")}:
    Extracts the fitted values.}
    \item{residuals}{\code{signature(object = "ACDmultifit")}: 
    Extracts the residuals. Optional logical argument \code{standardize} 
	(default is FALSE) allows to extract the standardized residuals.}
    \item{sigma}{\code{signature(object = "ACDmultifit")}:
    Extracts the conditional sigma values.}
    \item{skew}{\code{signature(object = "ACDmultifit")}:
    Extracts the conditional skew values. }
    \item{shape}{\code{signature(object = "ACDmultifit")}:
    Extracts the conditional shape values. }
    \item{coef}{\code{signature(object = "ACDmultifit")}:
    Extracts the coefficients.}
    \item{likelihood}{\code{signature(object = "ACDmultifit")}:
    Extracts the likelihood. }
	}
}
\author{Alexios Ghalanos}
\seealso{
Classes \code{\linkS4class{ACDmultiforecast}}, \code{\linkS4class{ACDmultifilter}} 
and \code{\linkS4class{ACDmultispec}}.
}
\keyword{classes}