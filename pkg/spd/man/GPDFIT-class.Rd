\name{GPDFIT-class}
\docType{class}
\alias{GPDFIT-class}
\alias{GPDFIT}
\alias{show,GPDFIT-method}

\title{Class: Generalized Pareto Distribution}
\description{ Locally implemented and slightly modified class for the generalized 
pareto distribution (gpd) fit, borrowed from package \code{fExtremes}. Created 
on modelling the tails of the data by \code{\link{spdfit}}}
\section{Objects from the Class}{
Objects of this class cannot be created by user as the methods are not exported.
}
\section{Slots}{
	 \describe{
    \item{\code{call}:}{Object of class \code{"call"} }
    \item{\code{method}:}{Object of class \code{"character"} }
    \item{\code{parameter}:}{Object of class \code{"list"} }
    \item{\code{data}:}{Object of class \code{"list"} }
    \item{\code{fit}:}{Object of class \code{"list"} }
    \item{\code{residuals}:}{Object of class \code{"numeric"} }
    \item{\code{title}:}{Object of class \code{"character"} }
    \item{\code{description}:}{Object of class \code{"character"} }
  }
}
\author{
Alec Stephenson for the functions from R\`s \code{"evd-package"}, 
Alec Stephenson for the functions from R\`s \code{"evir-package"}, 
Alexander McNeil for the EVIS functions underlying the \code{"evir-package"}, 
Diethelm Wuertz for the functions from R\`s \code{"fExtremes-package"},
M.P.Wand and M.C.Jones for the functions from R\`s \code{"KernSmooth-package"},
Alexios Ghalanos for this package.}
\references{
Carmona, R. and J. Morrisson (2001). Heavy Tails and Copulas with Evanesce, 
ORFE Tech. Report, Princeton University\cr
Carmona, R. (2001). Statistical Analysis of Financial Data, with an 
implementation in Splus\cr
Embrechts, P., Klueppelberg, C., Mikosch, T. (1997); \emph{Modelling Extremal 
Events}, Springer
}
\note{ S3 plot method exists which provides for visual inspection of the fit 
and is called by the higher level S3 plot method of the \code{\linkS4class{SPD}} 
class}
\keyword{classes}