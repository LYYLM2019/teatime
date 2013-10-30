\name{predict.dbm}
\alias{predict.dbm}
\title{
Iterated multi-step ahead predictions
}
\description{
Creates multi-step ahead iterated predictions from a dynamic binary model by
enumerating all possible future paths.
}
\usage{
\method{predict}{dbm}(object, newdata = NULL, n.ahead=1, ...)
}
\arguments{
 \item{object}{
An object of class dbm.
}
  \item{newdata}{
Not currently used by this package. Use the dbm function instead with fixed
parameters to get the filtered 1-ahead predictions.
}
\item{n.ahead}{
The number of forward forecasts to make (see details).
}
\item{\dots}{
...
}
}
\details{
The details of the iterated n-ahead prediction can be found in the paper by
Kauppi and Saikkonen (2008). Because of the complete enumeration of all possible 
paths, anything more than n.ahead=15 will likely result in out of memory problems.
}
\value{
A list with the following items:
\item{prediction}{The n.ahead forecasts.}
\item{path}{The matrix of all possible paths evaluated.}
}
\references{	
Kauppi, H., and Saikkonen,P. 2008, Predicting US recessions with dynamic 
binary response models, \emph{The Review of Economics and Statistics}, 90(4), 
777--791.\cr
Nyberg, H, 2011, Forecasting the direction of the US stock market with dynamic 
binary probit models, \emph{International Journal of Forecasting}, 27(2), 
561--578.
}
\author{
Alexios Ghalanos
}
\keyword{regression}
\keyword{classif}
