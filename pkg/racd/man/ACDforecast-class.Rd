\name{ACDforecast-class}
\docType{class}
\alias{ACDforecast-class}
\alias{fitted,ACDforecast-method}
\alias{shape,ACDforecast-method}
\alias{show,ACDforecast-method}
\alias{sigma,ACDforecast-method}
\alias{skew,ACDforecast-method}
\alias{skewness,ACDforecast-method}
\alias{kurtosis,ACDforecast-method}
\title{class: ACD Forecast Class}
\description{
S4 class containing the forecast ACD object.
}
\section{Objects from the Class}{
Objects can be created by calling the \code{acdforecast} function.
}
\section{Slots}{
  \describe{
    \item{\code{forecast}:}{Object of class \code{"vector"} forecast model.}
    \item{\code{model}:}{Object of class \code{"vector"} model details. }
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "ACDforecast")}: summary method. }
    \item{fitted}{\code{signature(object = "ACDforecast")}: conditional mean. }
    \item{sigma}{\code{signature(object = "ACDforecast")}: conditional sigma. }
    \item{skew}{\code{signature(object = "ACDforecast")}: conditional skew. }
    \item{shape}{\code{signature(object = "ACDforecast")}: conditional shape. }
    \item{skewness}{\code{signature(object = "ACDforecast")}: conditional skewness. }    
    \item{kurtosis}{\code{signature(object = "ACDforecast")}: conditional excess kurtosis.}
	 }
}
\author{Alexios Ghalanos}
\keyword{classes}