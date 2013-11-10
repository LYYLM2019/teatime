\name{ACDmultifilter-class}
\Rdversion{1.1}
\docType{class}
\alias{ACDmultifilter-class}
\alias{coef,ACDmultifilter-method}
\alias{fitted,ACDmultifilter-method}
\alias{likelihood,ACDmultifilter-method}
\alias{residuals,ACDmultifilter-method}
\alias{shape,ACDmultifilter-method}
\alias{sigma,ACDmultifilter-method}
\alias{skew,ACDmultifilter-method}
\title{class: Univariate ACD Multiple Filter Class}
\description{
Class for the univariate ACD Multiple filter.
}
\section{Objects from the Class}{
Objects can be created by calls of the function \code{\link{multifilter}}. 
}
\section{Slots}{
  \describe{
    \item{\code{filter}:}{Object of class \code{"vector"}}
    \item{\code{desc}:}{Object of class \code{"vector"}}
  }
}
\section{Methods}{
  \describe{
     \item{coef}{\code{signature(object = "ACDmultifilter")}: 
    coefficient extraction.}
    \item{convergence}{\code{signature(object = "ACDmultifilter")}: 
    estimation convergence indicator. }
    \item{fitted}{\code{signature(object = "ACDmultifilter")}: 
    conditional mean dynamics xts matrix.}
    \item{likelihood}{\code{signature(object = "ACDmultifilter")}: 
    conditional log-likelihood at the optimal (vector).}
    \item{residuals}{\code{signature(object = "ACDmultifilter")}: 
    residuals xts matrix. }
    \item{shape}{\code{signature(object = "ACDmultifilter")}: 
    conditional shape dynamics xts matrix, with optional argument
    \sQuote{transformed} on whether to returned the transformed
    dynamics (default: TRUE).}
    \item{sigma}{\code{signature(object = "ACDmultifilter")}: 
    conditional volatility xts matrix.}
    \item{skew}{\code{signature(object = "ACDmultifilter")}:  
    conditional skew dynamics xts matrix, with optional argument
    \sQuote{transformed} on whether to returned the transformed
    dynamics (default: TRUE)}
	 }
}
\author{Alexios Ghalanos}
\keyword{classes}