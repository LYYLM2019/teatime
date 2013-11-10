\name{ACDmultifit-class}
\Rdversion{1.1}
\docType{class}
\alias{ACDmultifit-class}
\alias{coef,ACDmultifit-method}
\alias{convergence,ACDmultifit-method}
\alias{fitted,ACDmultifit-method}
\alias{likelihood,ACDmultifit-method}
\alias{multifilter,ACDmultifit-method}
\alias{residuals,ACDmultifit-method}
\alias{shape,ACDmultifit-method}
\alias{sigma,ACDmultifit-method}
\alias{skew,ACDmultifit-method}
\title{class: Univariate ACD Multiple Fit Class}
\description{
Class for the univariate ACD Multiple fit.
}
\section{Objects from the Class}{
Objects can be created by calls of the function \code{\link{multifit}}. 
}
\section{Slots}{
  \describe{
    \item{\code{fit}:}{Object of class \code{"vector"} }
    \item{\code{desc}:}{Object of class \code{"vector"} }
  }
}
\section{Methods}{
  \describe{
    \item{coef}{\code{signature(object = "ACDmultifit")}: 
    coefficient extraction.}
    \item{convergence}{\code{signature(object = "ACDmultifit")}: 
    estimation convergence indicator. }
    \item{fitted}{\code{signature(object = "ACDmultifit")}: 
    conditional mean dynamics xts matrix.}
    \item{likelihood}{\code{signature(object = "ACDmultifit")}: 
    conditional log-likelihood at the optimal (vector).}
    \item{multifilter}{\code{signature(multifitORspec = "ACDmultifit")}:\cr
    function:\cr
    \bold{multifilter(multifitORspec, data = NULL, out.sample = 0, 
		n.old = NULL, cluster = NULL, skew0 = NULL, shape0 = NULL, ...)}
    }
    \item{residuals}{\code{signature(object = "ACDmultifit")}: 
    residuals xts matrix. }
    \item{shape}{\code{signature(object = "ACDmultifit")}: 
    conditional shape dynamics xts matrix, with optional argument
    \sQuote{transformed} on whether to returned the transformed
    dynamics (default: TRUE).}
    \item{sigma}{\code{signature(object = "ACDmultifit")}: 
    conditional volatility xts matrix.}
    \item{skew}{\code{signature(object = "ACDmultifit")}:  
    conditional skew dynamics xts matrix, with optional argument
    \sQuote{transformed} on whether to returned the transformed
    dynamics (default: TRUE)}
	 }
}
\author{Alexios Ghalanos}
\keyword{classes}