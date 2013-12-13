\name{spd-package}
\alias{spd-package}
\alias{spd}
\docType{package}
\title{The Semi-Parametric Distribution (spd) package}
\description{The Semi-Parametric Distribution is a piecewise distribution 
constructed by parametrically modelling the tails of the distribution
using an appropriate distribution (e.g. generalized pareto) and the interior by 
kernel methods. The package implements fit, distribution, density,
quantile and random number generation. Currently, only the generalized pareto 
distribution is implemented for modelling the tails, but the package can easily 
be extended.
}
\details{
\tabular{ll}{
Package: \tab spd\cr
Type: \tab Package\cr
Version: \tab 2.0-0\cr
Date: \tab 2013-12-15\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
Depends: \tab methods\cr
}
The main functionality of the package is contained in the \code{\linkS4class{SPD}} 
class, created by calling \code{\link{spdfit}}. Methods for density 
\code{\link{dspd}}, distribution \code{\link{pspd}}, quantile \code{\link{qspd}} 
and random number generation \code{\link{rspd}} exist and take 2 main arguments, 
the input value and the fitted object.\cr
The spd package uses the \code{"bkde"} function from the package KernSmooth for 
the kernel interior fit, while for the tail fit borrows from the fExtremes 
package and implements a locally modified copy of the gpd functionality and 
methods.
}
\author{
Alec Stephenson for the functions from R\`s \code{"evd-package"},\cr
Alec Stephenson for the functions from R\`s \code{"evir-package"},\cr
Alexander McNeil for the EVIS functions underlying the \code{"evir-package"},\cr
Diethelm Wuertz for the functions from R\`s \code{"fExtremes-package"},\cr
M.P.Wand and M.C.Jones for the functions from R\`s \code{"KernSmooth-package"},\cr
Alexios Ghalanos for this package.}
\references{
Carmona, R. and J. Morrisson (2001). Heavy Tails and Copulas with Evanesce, 
ORFE Tech. Report, Princeton University\cr
Carmona, R. (2001). Statistical Analysis of Financial Data, with an 
implementation in Splus\cr
Embrechts, P., Klueppelberg, C., Mikosch, T. (1997); \emph{Modelling Extremal 
Events}, Springer
}
\examples{\dontrun{
library(MASS)
x<-SP500/100
fit<-spdfit(x)
show(fit)
}}

