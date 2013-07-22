\name{SPD-class}
\docType{class}
\alias{SPD-class}
\alias{SPD}
\alias{show,SPD-method}
\title{Class: Semi-Parametric Distribution }
\description{ Virtual Class for holding the tail fit and kernel interior objects. }
\section{Methods}{
  \describe{
    \item{dspd}{\code{signature(x = "numeric", fit = "SPD", linear = "logical")}: 
    density function}
    \item{pspd}{\code{signature(q = "numeric", fit = "SPD", linear = "logical")}: 
    distribution function}
    \item{qspd}{\code{signature(p = "numeric", fit = "SPD", linear = "logical")}: 
    quantile function }
    \item{rspd}{\code{signature(n = "numeric", fit = "SPD", linear = "logical")}: 
    random number generation function}
    \item{show}{\code{signature(object = "SPD")}: show method }
	 }
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\references{Carmona, R. and J. Morrisson (2001). Heavy Tails and Copulas with 
Evanesce, ORFE Tech. Report, Princeton University\cr
Carmona, R. (2001). Statistical Analysis of Financial Data, with an 
implementation in Splus}
\author{ Alexios Ghalanos }
\examples{
showClass("SPD")
}
\keyword{classes}
