\name{ACDfilter-class}
\docType{class}
\alias{ACDfilter-class}
\alias{coef,ACDfilter-method}
\alias{fitted,ACDfilter-method}
\alias{kurtosis,ACDfilter-method}
\alias{likelihood,ACDfilter-method}
\alias{plot,ACDfilter,missing-method}
\alias{quantile,ACDfilter-method}
\alias{pit,ACDfilter-method}
\alias{residuals,ACDfilter-method}
\alias{shape,ACDfilter-method}
\alias{show,ACDfilter-method}
\alias{sigma,ACDfilter-method}
\alias{skew,ACDfilter-method}
\alias{skewness,ACDfilter-method}
\title{class: ACD Filter Class}
\description{
Object holding the details of the ACD filter.
}
\section{Objects from the Class}{
Objects can be created by calling the \code{acdfilter} function.
}
\section{Slots}{
  \describe{
    \item{\code{filter}:}{Object of class \code{"vector"} filtered model.}
    \item{\code{model}:}{Object of class \code{"vector"} model details.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{rACD}"}, directly.
}
\section{Methods}{
  \describe{
    \item{coef}{\code{signature(object = "ACDfilter")}: estimated coefficients. }
    \item{likelihood}{\code{signature(object = "ACDfilter")}: log-likelihood of the 
    ACD model.}
    \item{show}{\code{signature(object = "ACDfilter")}: summary method. }
    \item{plot}{\code{signature(x = "ACDfilter", y = "missing")}: plots.}
    \item{fitted}{\code{signature(object = "ACDfilter")}: conditional mean. }
    \item{residuals}{\code{signature(object = "ACDfilter")}: model residuals 
    (with logical option standardize set to FALSE as default).}
    \item{sigma}{\code{signature(object = "ACDfilter")}: conditional sigma. }
    \item{skew}{\code{signature(object = "ACDfilter")}: conditional skew. Takes 
    optional argument \sQuote{transform} (default TRUE) which indicates whether
    to return the transformed value (logistic transformation).}
    \item{shape}{\code{signature(object = "ACDfilter")}: conditional shape. Takes 
    optional argument \sQuote{transform} (default TRUE) which indicates whether
    to return the transformed value (exponential transformation).}
    \item{skewness}{\code{signature(object = "ACDfilter")}: conditional skewness. }    
    \item{kurtosis}{\code{signature(object = "ACDfilter")}: conditional excess kurtosis.}
    \item{quantile}{\code{signature(x = "ACDfilter")}: the conditional quantiles given
    a vector of probabilities. }
    \item{pit}{\code{signature(x = "ACDfilter")}: the probability integral 
    transform of the data given the conditional density.}
	 }
}
\author{Alexios Ghalanos}
\keyword{classes}
