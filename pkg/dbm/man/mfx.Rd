\name{mfx.dbm}
\alias{mfx.dbm}
\alias{mfx}
\alias{mfxTrend}
\alias{mfxTrend.mfx}
\alias{summary.mfx}
\alias{print.mfx}
\alias{plot.mfxTrend}
\title{
Marginal Effects Diagnostics
}
\description{
Returns the marginal effects of the model.
}
\usage{
\method{mfx}{dbm}(model, x.mean = TRUE, rev.dum = TRUE, ...)
\method{mfxTrend}{mfx}(object, n = 300, cregressor, bregressor, ...)
\method{plot}{mfxTrend}(x,...)
\method{summary}{mfx}(object,...)
\method{print}{mfx}(x,...)
}
\arguments{
\item{model}{
[mfx] An object of class dbm.
}
\item{object}{
[mfxTrend, summary] An object of class mfx.
}
\item{x}{
[mfxTrend] An object of class mfxTrend.\cr
[print] An object of class mfx.\cr
}
\item{x.mean}{
[mfx] Whether to calculate marginal effects at the means of independent 
variables (default: TRUE), else marginal effects are calculated for each 
observation and then averaged.
}
\item{rev.dum}{
[mfx]  Whether to revise the estimates and standard errors for binary 
independent variables (default: TRUE), else derivatives are taken on binary 
independent variables as continuous variables.
}
\item{n}{
[mfxTrend] The number of points for calculating the probability, with larger 
numbers leading to a smoother curve.
}
\item{cregressor}{
[mfxTrend] Required name of a continuous explanatory variable in the regression
for which to calculate the change in probability.
}
\item{bregressor}{
[mfxTrend] An optional name of a binary independent variable used to stratify 
the probability.
}
\item{\dots}{
Additional arguments passed to the plot or print functions.
}
}
\value{
The mfx object returns a list with a number of value, the most important of
which is the \sQuote{out} which is a a data frame object of marginal effects, 
t-values, and p-values.
}
\details{
The functions are imported from the erer package of Changyou Sun, and adapted 
locally in this package to work with the dbm class. See the \code{maBina} and
\code{maTrend} for the original functions in the erer package.
}
\examples{
\dontrun{
require(xts)
require(TTR)
data(recession)
tmp = recession$T10Y - recession$TB3M
recession$spread =  tmp - EMA(tmp, 6)
recession$cpiadj = ROC(recession$CPI,n=12,type="discrete")-SMA(ROC(recession$CPI,n=12,type="discrete"), 6)
y = recession[,c("OECDR","spread","cpiadj")]
y = na.omit(y)
mod = dbm(y, x.vars = c("spread","cpiadj"), x.lags = 3, arp = 1, arq = 1)
marg = mfx(mod)
margt = mfxTrend(marg, cregressor = "cpiadj")
par(mfrow=c(2,1))
plot(recession$cpiadj)
plot(margt)
}
}
\author{
Alexios Ghalanos
}