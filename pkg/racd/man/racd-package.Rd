\name{racd-package}
\alias{racd-package}
\alias{racd}
\title{The racd package}
\description{
The racd package provides a set of methods for estimating, filtering, forecasting
and simulating Autoregressive Conditional Density (ACD) models. These are 
effectively extensions of the GARCH models but allow for time variation in the
conditional skew and shape parameters of the distribution. Because of the 
transformations required to restrict these parameters within their distribution
specific bounds, and the fact that skewness and kurtosis are driven by extreme
events making their identification with a particular law of motion difficult, 
the optimization of these models is particularly hard. The racd package provides
for a number of different dynamics for the skew and shape parameters, and a number
of different solvers and approaches to solving such problems and described in more
details in the package's vignette. The multivariate independent factor ACD model
(IFACD) of Ghalanos et al (2014) is also included along with methods similar to
those of the GO-GARCH model in the rmgarch package.
}
\details{
\tabular{ll}{
Package: \tab racd\cr
Type: \tab Package\cr
Version: \tab 1.0-3\cr
Date: \tab 2014-02-23\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
Depends: \tab methods, rugarch\cr
Imports: \tab xts, Rsolnp, ucminf, numDeriv, zoo, parma, Matrix\cr}
The racd package builds on the rugarch package and imports many of the methods
used there. In terms of GARCH models, only the simple GARCH (sGARCH) and 
component GARCH (csGARCH) and multiplicative component GARCH (mcsGARCH) are used, 
and this is because of the variation in the conditional skew and shape parameters 
of the distributions which makes the calculation of persistence and the use of 
variance targeting near impossible except for these models.
}
\section{How to cite this package}{
  Whenever using this package, please cite as\cr

 \preformatted{@Manual{Ghalanos_2013,
 author       = {Alexios Ghalanos},
 title        = {{racd}: Autoregressive Conditional Density Models.},
 year         = {2014},
 note 	      = {R package version 1.0-3.},}}
\emph{}
}
\section{License}{
  The releases of this package is licensed under GPL version 3.
}
\author{Alexios Ghalanos}
\references{
Hansen, B. E. 1994, Autoregressive conditional density estimation, 
\emph{International Economic Review}, \bold{35(3)}, 705--730.\cr
Ghalanos, A., Rossi E., and Urga G. 2013, Independent Factor Autoregressive 
Conditional Density Model, \emph{Econometric Reviews}, forthcoming.
}