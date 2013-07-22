\name{ACDroll-class}
\docType{class}
\alias{ACDroll-class}
\alias{as.data.frame,ACDroll-method}
\alias{show,ACDroll-method}
\alias{resume,ACDroll-method}
\alias{quantile,ACDroll-method}
\alias{pit,ACDroll-method}
\alias{convergence,ACDroll-method}
\alias{coef,ACDroll-method}
\alias{skew,ACDroll-method}
\alias{skewness,ACDroll-method}
\alias{shape,ACDroll-method}
\alias{kurtosis,ACDroll-method}
\alias{show,ACDroll-method}
\title{class: ACD Rolling Forecast Class}
\description{Class for the ACD rolling forecast.}
\section{Slots}{
  \describe{
    \item{\code{forecast}:}{Object of class \code{"vector"} }
    \item{\code{model}:}{Object of class \code{"vector"} }
  }
}
\section{Methods}{
  \describe{
    \item{as.data.frame}{\code{signature(x = "ACDroll")}: 
    Extracts various values from object (see note).}
	\item{resume}{\code{signature(object = "ACDroll")}:
	Resumes a rolling backtest which has non-converged windows using 
	alternative solver and control parameters.}
	\item{coef}{\code{signature(object = "ACDroll")}: 
    Extracts the list of coefficients for each estimated window in the 
    rolling backtest.}
    \item{show}{\code{signature(object = "ACDroll")}: 
    Summary.}
    \item{quantile}{\code{signature(x = "ACDroll")}: 
	Calculates and returns, given a vector of probabilities (additional argument 
	\dQuote{probs}), the conditional quantiles of the rolling object as an
	xts matrix.}
	\item{pit}{\code{signature(object = "ACDroll")}: 
	Calculates and returns the conditional probability integral transform given the 
	realized data and forecast density.}
	 \item{skew}{\code{signature(object = "ACDroll")}: conditional skew. }
    \item{shape}{\code{signature(object = "ACDroll")}: conditional shape. }
    \item{skewness}{\code{signature(object = "ACDroll")}: conditional skewness. }    
    \item{kurtosis}{\code{signature(object = "ACDroll")}: conditional excess kurtosis.}
    \item{convergence}{\code{signature(object = "ACDroll")}: 
    Returns the convergence code for the estimation windows, with 0 indicating
    that all have converged and 1 that there were non-converged windows. In 
    the latter case the \sQuote{nonconverged} attribute is also printed of those 
    windows which failed to converge.}
    \item{show}{\code{signature(object = "ACDroll")}:
	Summary.}
	 }
}
\note{
The \code{as.data.frame} extractor method allows the extraction of either the
conditional forecast density or the VaR. It takes additional argument
\code{which} with valid values either \dQuote{density} or \dQuote{VaR}.\cr
The \code{coef} method will return a list of the coefficients and their robust 
standard errors (assuming the keep.coef argument was set to TRUE in the 
acdroll function), and the ending date of each estimation window.\cr
}
\author{Alexios Ghalanos}
\keyword{classes}