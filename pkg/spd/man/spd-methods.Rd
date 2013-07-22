\name{spd-methods}
\alias{pspd}
\alias{dspd}
\alias{qspd}
\alias{rspd}
\alias{dspd-methods}
\alias{pspd-methods}
\alias{qspd-methods}
\alias{rspd-methods}
\alias{rspd,ANY,ANY,ANY-method}
\alias{rspd,vector,SPD,logical-method}
\alias{rspd,vector,SPD,missing-method}
\alias{dspd,ANY,ANY,ANY-method}
\alias{dspd,vector,SPD,logical-method}
\alias{dspd,vector,SPD,missing-method}
\alias{pspd,ANY,ANY,ANY-method}
\alias{pspd,vector,SPD,logical-method}
\alias{pspd,vector,SPD,missing-method}
\alias{qspd,ANY,ANY,ANY-method}
\alias{qspd,vector,SPD,logical-method}
\alias{qspd,vector,SPD,missing-method}
\title{Method: Semi-Parametric Distribution}
\description{
Density, Distribution, Quantile and Random Number Generation methods for the 
Semi-Parametric Distribution.
}
\usage{
dspd(x, fit, linear)
pspd(q, fit, linear)
qspd(p, fit, linear)
rspd(n, fit, linear)
}
\arguments{
    
    \item{n}{
        [rspd] - \cr
        the number of random deviates to be generated from fitted distribution.
        } 
    \item{p}{
        a vector of probability levels, the desired probability for the 
        quantile estimate (e.g. 0.99 for the 99th percentile).
        } 
    \item{x,q}{
        [pspd,dspd] - \cr
        a numeric vector of quantiles.
        }
    \item{fit}{
        [all] - \cr
        the object of class SPDFIT created by calling the spdfit function.
        }
    \item{linear}{
	[all] - \cr
        logical, if TRUE interior smoothing function uses linear interpolation 
        rather than constant.
     }
}
\value{
    All values are numeric vectors: \cr
    \code{d*} returns the density (pdf), \cr
    \code{p*} returns the probability (cdf), \cr
    \code{q*} returns the quantiles (inverse cdf), and \cr
    \code{r*} generates random deviates.
}
\note{
The density is computed using the generalized pareto distribution in the tails, 
while for the middle, the density is computed by using a smooth gradient approach. 
Interpolation is used to splice together the ends with the middle portion,
providing for an approximate piecewise constant density function. As such, 
caution should be used when interpreting results obtained by use of this function.
}
\author{
    Alec Stephenson for the functions from R\`s \code{evd} package, \cr
    Alec Stephenson for the functions from R\`s \code{evir} package, \cr
    Alexander McNeil for the EVIS functions underlying the \code{evir} package, \cr
    Diethelm Wuetrz for the \code{fExtremes} Implementation of the gpd, \cr
    Alexios Ghalanos for the SPD Implementation, \cr
}
\references{
Embrechts, P., Klueppelberg, C., Mikosch, T. (1997); \emph{Modelling Extremal 
Events}, Springer.\cr
Carmona, R. (2004);\emph{Statistical Anlaysis of Financial Data in Splus}, 
Springer.
}
\examples{\dontrun{
library(MASS)
x = SP500/100
fit=spdfit(x, upper=0.9, lower=0.1)
## rspd  -
   par(mfrow = c(2, 2), cex = 0.7)
   r = rspd(n = 1000, fit)
   hist(r, n = 100, probability = TRUE, xlab = "r", 
   col = "steelblue", border = "white",main = "Density")
   box()
## dspd -
   # Plot empirical density and compare with true density:
   r = rspd(n = 1000, fit)
   hist(r, n = 100, probability = TRUE, xlab = "r", 
   col = "steelblue", border = "white",main = "Density")
   box()
   x = seq(-0.3, 0.3, length.out = 1000)
   lines(x, dspd(x, fit), col = "darkorange",lwd=2)
   
## pspd -
   # Plot df and compare with true df:
   plot(sort(r), (1:length(r)/length(r)), 
   ylim = c(0, 1), pch = 19, 
   cex = 0.5, ylab = "p", xlab = "q", main = "CDF")
   grid()
   q = seq(-0.3, 0.3, length.out = 1000)
   lines(q, pspd(q, fit), col = "darkorange",lwd=2)
}}
\keyword{distribution}