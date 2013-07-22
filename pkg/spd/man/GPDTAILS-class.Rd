\name{GPDTAILS-class}
\docType{class}
\alias{GPDTAILS-class}
\alias{GPDTAILS}
\title{Class: GPDTAILS}
\description{Class: SPD with Generalized Pareto Distribution (GPD) Tails }
\section{Objects from the Class}{
Objects can be created by calling \code{\link{spdfit}}. The main implemented 
class of the spd package holding the details of the fitted object with gpd tails.}
\section{Slots}{
	 \describe{
    \item{\code{call}:}{...}
    \item{\code{method}:}{the gpd fitting method}
    \item{\code{kernel}:}{the kernel type}
    \item{\code{data}:}{the original dataset}
    \item{\code{threshold}:}{the upper and lower thresholds fitted}
    \item{\code{ptails}:}{the upper and lower cutoff points.}
    \item{\code{fit}:}{the upper and lower gpd fit and the interior kernel fit objects.}
    \item{\code{title}:}{optional title of project}
    \item{\code{description}:}{optional description}
  }
}
\section{Extends}{
Class \code{\linkS4class{SPD}}, directly.
}
\references{Carmona, R. and J. Morrisson (2001). Heavy Tails and Copulas with 
Evanesce, ORFE Tech. Report, Princeton University\cr
Carmona, R. (2001). Statistical Analysis of Financial Data, with an 
implementation in Splus}
\author{ Alexios Ghalanos }
\keyword{classes}