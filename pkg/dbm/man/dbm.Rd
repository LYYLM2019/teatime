\name{dbm}
\alias{dbm}
\title{
Estimation of Dynamic Binary Model
}
\description{
Specifies and estimates a dynamic binary model using maximum likelihood, 
with a choice of link functions and making use of analytic gradient.
}
\usage{
dbm(y, x.vars = NULL, x.lags = 1, arp = 1, arq = 0, ecm = FALSE, constant = TRUE, 
link = "gaussian", fixed.pars = NULL, solver = "optim", control = list(), 
parsearch = TRUE, parsim = 5000, method = "Nelder-Mead", regularization = FALSE, 
reg.cost = 1, ...)
}
\arguments{
 \item{y}{
A multivariate xts object with named columns.
}
  \item{x.vars}{
An optional character vector of the column headings for the regressors. Whatever
is not a regressor is then taken to be the binary regressand variable.
}
  \item{x.lags}{
A vector of lags for each regressor. If scalar will assume that all regressors
have the same lags. If you want to include multiple lags for the same regressor,
then the y object should contain multiple instances of that regressor.
}
  \item{arp}{
The autoregressive lag for the dynamics.
}
  \item{arq}{
The autoregressive lags for the regressand variable.
}
  \item{ecm}{
(logical) Whether to enforce an error correction specification (see details).
}
  \item{constant}{
(logical) Whether to include a constant.
}
  \item{link}{
The link (cdf transformation) function (\dQuote{gaussian} for probit, 
\dQuote{logistic} for logit, and dQuote{glogistic} for skew-logit).
}
  \item{fixed.pars}{
An optional named vector of parameters. If all the parameters are fixed then
the data supplied is filtered (not estimated) with those parameters.
}
  \item{solver}{
A choice of either \sQuote{optim} or \sQuote{gosolnp}.
}
  \item{control}{
Solver control parameters.
}
  \item{parsearch}{
Whether to perform an initial parameter search to obtain starting values.
}
  \item{parsim}{
If parsearch is TRUE, the number of parameter vectors to simulate and evaluate.
}
  \item{method}{
The \sQuote{optim} sub-solver.
}
\item{regularization}{
Whether to use L2 regularization.
}
\item{reg.cost}{
The regularization penalty.
}
  \item{\dots}{
Additional parameters are passed to the \sQuote{gosolnp} solver.
}
}
\details{
The function estimates the dynamic binary model of Kauppi and Saikkonen (2008) 
with the option of enforcing an error correction restriction described in 
Nyberg (2011). This is best described and illustrated in the package vignette.
The likelihood and analytic gradient evaluation is coded in C for speed, and 
there is a choice of 3 solvers (only 2 of which make use of the analytic 
gradient), with the Nelder-Mead method of the optim solver being the default and
recommended way to solve this model. Because of the absence of upper and lower
bounds in the case of optim, the autoregressive parameter on the dynamics 
equation  (\dQuote{alpha}) is restricted between -1 and 1 using a logistic
transformation during the estimation phase, and similarly for the power parameter
when using the Generalized Logistic link function which is constrained to be 
positive. The bounding transformation is lifted once the parameters are estimated 
and transformed so that standard errors and related statistics are correctly 
calculated.
}
\value{
An S3 object of class dbm, as a list with the following items:
\item{fit}{Details of the estimated object.}
\item{model}{Details of the specification used.}
}
\references{	
Kauppi, Heikki, and Pentti Saikkonen, 2008, Predicting US recessions with dynamic 
binary response models, \emph{The Review of Economics and Statistics}, 90(4), 
777--791.\cr
Nyberg, H, 2011, Forecasting the direction of the US stock market with dynamic 
binary probit models, \emph{International Journal of Forecasting}, 27(2), 
561--578.
}
\author{
Alexios Ghalanos
}
\note{
It is possible that the routine will not converge or that the converged estimate
may not be globally optimal. It is therefore advisable to run the estimation a
few times (a random parameter search is conducted each time for starting values)
to make sure. Your mileage will vary depending on the degree of multicollinearity
in the dataset, non-stationarity, partial separation and heteroscedasticity.
}
\keyword{regression}
\keyword{classif}