\name{ACDmultiforecast-class}
\docType{class}
\alias{ACDmultiforecast-class}
\alias{sigma,ACDmultiforecast-method}
\alias{fitted,ACDmultiforecast-method}
\alias{skew,ACDmultiforecast-method}
\alias{shape,ACDmultiforecast-method}
\title{class: Univariate ACD Multiple Forecast Class}
\description{Class for the univariate ACD Multiple forecast.}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Methods}{
\describe{
  \item{sigma}{\code{signature(x = "ACDmultiforecast")}: 
   		extracts the n.ahead by (n.roll+1) by n.assets array of conditional sigma 
   		forecasts.}
 \item{fitted}{\code{signature(x = "ACDforecast")}: 
     extracts the n.ahead by (n.roll+1) by n.assets array of conditional mean
   		forecasts.}
  \item{skew}{\code{signature(x = "ACDforecast")}: 
     extracts the n.ahead by (n.roll+1) by n.assets array of conditional skew
   		forecasts.}
  \item{shape}{\code{signature(x = "ACDforecast")}:
     extracts the n.ahead by (n.roll+1) by n.assets array of conditional shape
   		forecasts.}
  }
}
\author{Alexios Ghalanos}
\seealso{
Classes \code{\linkS4class{ACDmultifilter}}, \code{\linkS4class{ACDmultifit}} and 
\code{\linkS4class{ACDmultispec}}.
}
\keyword{classes}