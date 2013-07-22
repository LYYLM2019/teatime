\name{ACDsim-class}
\docType{class}
\alias{ACDsim-class}
\alias{fitted,ACDsim-method}
\alias{kurtosis,ACDsim-method}
\alias{shape,ACDsim-method}
\alias{sigma,ACDsim-method}
\alias{skew,ACDsim-method}
\alias{skewness,ACDsim-method}
\title{class: ACD Simulation Class}
\description{Class for the ACD simulation.}
\section{Slots}{
  \describe{
    \item{\code{simulation}:}{Object of class \code{"vector"} The simulation. }
    \item{\code{model}:}{Object of class \code{"vector"} The model. }
    \item{\code{seed}:}{Object of class \code{"integer"} The random seeds. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{rACD}"}, directly.
}
\section{Methods}{
  \describe{
   \item{sigma}{\code{signature(object = "ACDsim")}: 
	Extracts the conditional sigma simulated values as a matrix of size n.sim x m.sim.}
	\item{fitted}{\code{signature(object = "ACDsim")}: 
	Extracts the conditional mean simulated values as a matrix of size n.sim x m.sim.}    
    \item{kurtosis}{\code{signature(object = "ACDsim")}: Calculates and 
    Returns the conditional excess simulated Kurtosis values as a matrix of 
    size n.sim x m.sim}
    \item{shape}{\code{signature(object = "ACDsim")}: Extracts the conditional 
    shape simulated values as a matrix of size n.sim x m.sim. Takes optional
    argument \sQuote{transform} (default TRUE) denoting whether to transform the
    value to the distribution domain bounds.}
    \item{skew}{\code{signature(object = "ACDsim")}: Extracts the conditional 
    skew simulated values as a matrix of size n.sim x m.sim. Takes optional
    argument \sQuote{transform} (default TRUE) denoting whether to transform the
    value to the distribution domain bounds. }
    \item{skewness}{\code{signature(object = "ACDsim")}: 
    Returns the conditional simulated Skewness values as a matrix of 
    size n.sim x m.sim}
    }
}
\author{Alexios Ghalanos}
\keyword{classes}
