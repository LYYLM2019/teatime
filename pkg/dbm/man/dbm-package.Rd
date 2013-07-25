\name{dbm-package}
\alias{dbm-package}
\title{The dbm package}
\description{
The dbm package implements the dynamic binary response model of Kauppi and 
Saikkonen (2008), including the error correction restriction of Nyberg (2011), 
with methods for maximum likelihood estimation and iterated multi-step ahead 
prediction. It allows for either probit, logit or skew-logit link functions 
and provides a set of methods for working with the estimated object.
}
\details{
\tabular{ll}{
Package: \tab dbm\cr
Type: \tab Package\cr
Version: \tab 1.0-0\cr
Date: \tab 2013-07-15\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
Imports: \tab Rsolnp, numDeriv, xts, zoo, chron, nloptr, Epi}
The type of data handled by the package is completely based on the xts package, 
and only data which can be coerced to such will be accepted by the package.}
\section{How to cite this package}{
  Whenever using this package, please cite as\cr
 \preformatted{@Manual{Ghalanos_2013,
 author       = {Alexios Ghalanos},
 title        = {{dbm}: Dynamic Binary Model.},
 year         = {2013},
 note 	      = {R package version 1.0-0.},}}
\emph{}
}
\section{License}{
  The releases of this package is licensed under GPL version 3.
}
\author{Alexios Ghalanos}
\references{
Kauppi, Heikki, and Pentti Saikkonen, 2008, Predicting US recessions with dynamic 
binary response models, \emph{The Review of Economics and Statistics}, 90(4), 
777--791.\cr
Nyberg, H, 2011, Forecasting the direction of the US stock market with dynamic 
binary probit models, \emph{International Journal of Forecasting}, 27(2), 
561--578.
}