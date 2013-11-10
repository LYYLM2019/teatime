\name{ACDmultispec-class}
\Rdversion{1.1}
\docType{class}
\alias{ACDmultispec-class}
\alias{multifilter,ACDmultispec-method}
\alias{multifit,ACDmultispec-method}
\alias{multiforecast,ACDmultispec-method}

\title{class: Univariate ACD Multiple Specification Class}
\description{
Class for the univariate ACD Multiple specification.
}
\section{Objects from the Class}{
Created by calling the \code{\link{acdmultispec}} function.
}
\section{Slots}{
  \describe{
    \item{\code{spec}:}{Object of class \code{"vector"} containing the 
    details of the specification.}
    \item{\code{type}:}{Object of class \code{"character"} determining
    whether it is a common or diverse specification.}
  }
}
\section{Methods}{
  \describe{
    \item{multifilter}{\code{signature(multifitORspec = "ACDmultispec")}\cr
    function:\cr
    \bold{multifilter(multifitORspec, data = NULL, out.sample = 0, 
		n.old = NULL, cluster = NULL, skew0 = NULL, shape0 = NULL, ...)}
    }
    \item{multifit}{\code{signature(multispec = "ACDmultispec")}:\cr
    function:\cr
    \bold{multifit(multispec, data, out.sample = 0, solver = "solnp", 
		solver.control = list(), fit.control = list(stationarity = 1, 
		scale = 0, n.sim=2000), cluster = NULL, ...)}}
    \item{multiforecast}{\code{signature(multifitORspec = "ACDmultispec")}:\cr
    function:\cr
    \bold{multiforecast(multifitORspec, data = NULL, n.ahead = 1, 
		n.roll = 0, out.sample = 0,external.forecasts = list(mregfor = NULL, 
				vregfor = NULL, skregfor = NULL, shregfor = NULL), 
		m.sim = 1000, cluster = NULL, ...)}}
	 }
}
\note{See the univariate ACD functions corresponding to these multi-based functions
for a description of each argument.}
\author{Alexios Ghalanos}
\keyword{classes}