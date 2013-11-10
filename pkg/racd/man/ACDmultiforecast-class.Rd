\name{ACDmultiforecast-class}
\Rdversion{1.1}
\docType{class}
\alias{ACDmultiforecast-class}
\alias{fitted,ACDmultiforecast-method}
\alias{shape,ACDmultiforecast-method}
\alias{sigma,ACDmultiforecast-method}
\alias{skew,ACDmultiforecast-method}
\title{class: Univariate ACD Multiple Forecast Class}
\description{
Class for the univariate ACD Multiple forecast.
}
\section{Objects from the Class}{
Objects can be created by calls of the function \code{\link{multiforecast}}. 
}
\section{Slots}{
  \describe{
    \item{\code{forecast}:}{Object of class \code{"vector"}}
    \item{\code{desc}:}{Object of class \code{"vector"}}
  }
}
\section{Methods}{
  \describe{
   \item{fitted}{\code{signature(object = "ACDmultiforecast")}: 
    conditional mean dynamics.}
    \item{residuals}{\code{signature(object = "ACDmultiforecast")}: 
    residuals.}
    \item{shape}{\code{signature(object = "ACDmultiforecast")}: 
    conditional shape dynamics, with optional argument
    \sQuote{transformed} on whether to returned the transformed
    dynamics (default: TRUE).}
    \item{sigma}{\code{signature(object = "ACDmultiforecast")}: 
    conditional volatility.}
    \item{skew}{\code{signature(object = "ACDmultiforecast")}:  
    conditional skew dynamics, with optional argument
    \sQuote{transformed} on whether to returned the transformed
    dynamics (default: TRUE)}
	 }
}
\author{Alexios Ghalanos}
\note{Depending on whether \sQuote{n.roll>0} was used 
(in which case \sQuote{n.ahead=1}), this will either return a matrix
or list of matrices.}
\keyword{classes}