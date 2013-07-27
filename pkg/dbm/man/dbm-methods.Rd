\name{plot.dbm}
\alias{vcov.dbm}
\alias{plot.dbm}
\alias{fitted.dbm}
\alias{residuals.dbm}
\alias{coef.dbm}
\alias{likelihood.dbm}
\alias{score.dbm}
\alias{score}
\alias{likelihood.default}
\alias{likelihood}
\alias{deviance.dbm}
\alias{summary.dbm}
\alias{print.dbm}
\title{
dbm Methods
}
\description{
A number of methods for accessing and extracting information from an estimated
dbm object.
}
\usage{
\method{vcov}{dbm}(object, robust = FALSE, ...)
\method{plot}{dbm}(x, ...)
\method{fitted}{dbm}(object,...)
\method{residuals}{dbm}(object,...)
\method{coef}{dbm}(object,...)
\method{deviance}{dbm}(object, null=FALSE,...)
\method{likelihood}{dbm}(object, ...)
\method{likelihood}{default}(object, ...)
\method{summary}{dbm}(object,...)
\method{score}{dbm}(object, pars = NULL, analytic = TRUE, ... )
\method{print}{dbm}(x, digits = max(3L, getOption("digits") - 3L), ...)
}
\arguments{
\item{x}{
An object of class dbm.
}
\item{object}{
An object of class dbm.
}
\item{robust}{
[vcov] Extracts the robust covariance matrix.
}
\item{null}{
[deviance] Whether to return the null deviance.
}
\item{digits}{
No. of digits for priting the coefficient matrix.
}
\item{pars}{
[score] An (optional) named vector of coefficients to be used instead
of the object's esimated parameters.
}
\item{analytic}{
[score] Whether to calculate (if pars is not null, else no calculation is 
necessary) and return the analytic scores else the numerically estimated
values are returned (as calculated by the \code{jacobian} function of the 
numDeriv package).
}
\item{\dots}{
}
}
\details{
The methods provided allow for extracting or calculating useful values from
the estimated object of class \sQuote{dbm}.
}
\examples{
\dontrun{
require(xts)
require(TTR)
data(recession)
tmp = recession$T10Y - recession$TB3M
recession$spread =  tmp - EMA(tmp, 6)
recession$cpiadj = ROC(recession$CPI,n=12,type="discrete")-SMA(ROC(recession$CPI,n=12,type="discrete"), 6)

##############################
# EXAMPLE 1: Lagrange Multiplier Test (is alpha[1]=0?)
##############################
y = recession[,c("OECDR","spread","cpiadj")]
y = na.omit(y)
mod = dbm(y, x.vars = c("spread","cpiadj"), x.lags = 3, arp = 1, arq = 1,
fixed.pars = c("alpha[1]"=0))
# extract the coefficients and add back the alpha[1] parameter with the restricted
# value (zero) so that we can re-calculate the scores:
cf = c(coef(mod), "alpha[1]"=0)
S = score(mod, pars =  cf, analytic=TRUE)
# Use Outer-Product of Gradients Hessian
H = solve(crossprod(S))
G = colSums(S)
# stat = as.numeric(t(G) %*% H %*% G)
stat = as.numeric(t(G) \%*\% H \%*\% G)
pvalue = pchisq(stat, 1,lower.tail=FALSE)
if(pvalue<0.1) cat("\nReject H0 (restricted model) at the 10 pct significance level\n")
##############################
# EXAMPLE 2: Different Lags Same Explanatory Variable
##############################
# assume we want to include lags 1, 3, and 5 on the spread
y = recession[,c("OECDR","spread","spread","spread")]
colnames(y)[2:4]<-c("spread1","spread3","spread5")
y = na.omit(y)
mod = dbm(y, x.vars = c("spread1","spread3","spread5"), x.lags = c(1,3,5), 
arp = 1, arq = 1)
##############################
# EXAMPLE 3: Methods
##############################
y = recession[,c("OECDR","spread","cpiadj")]
y = na.omit(y)
mod = dbm(y, x.vars = c("spread","cpiadj"), x.lags = c(1,2), arp = 1, arq = 1)
likelihood(mod)
# null= TRUE returns the null deviance
c(deviance(mod), deviance(mod, null=TRUE))
# standard and robust covariance matrices
vcov(mod)
vcov(mod, robust=TRUE)
# check analytical and numeric scores:
head(score(mod, analytic=TRUE) - score(mod, analytic=FALSE))
}
}
\author{
Alexios Ghalanos
}