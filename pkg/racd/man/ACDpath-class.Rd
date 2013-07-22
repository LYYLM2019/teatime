\name{ACDpath-class}
\docType{class}
\alias{ACDpath-class}
\alias{fitted,ACDpath-method}
\alias{kurtosis,ACDpath-method}
\alias{shape,ACDpath-method}
\alias{sigma,ACDpath-method}
\alias{skew,ACDpath-method}
\alias{skewness,ACDpath-method}
\title{class: ACD Path Simulation Class}
\description{Class for the ACD Path simulation.}
\section{Slots}{
  \describe{
    \item{\code{path}:}{Object of class \code{"vector"} The simulation. }
    \item{\code{model}:}{Object of class \code{"vector"} The model. }
    \item{\code{seed}:}{Object of class \code{"integer"} The random seeds. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{rACD}"}, directly.
}
\section{Methods}{
  \describe{
   \item{sigma}{\code{signature(object = "ACDpath")}: 
	Extracts the conditional sigma simulated values as a matrix of size n.sim x m.sim.}
	\item{fitted}{\code{signature(object = "ACDpath")}: 
	Extracts the conditional mean simulated values as a matrix of size n.sim x m.sim.}    
    \item{kurtosis}{\code{signature(object = "ACDpath")}: Calculates and 
    Returns the conditional excess simulated Kurtosis values as a matrix of 
    size n.sim x m.sim}
    \item{shape}{\code{signature(object = "ACDpath")}: Extracts the conditional 
    shape simulated values as a matrix of size n.sim x m.sim. Takes optional
    argument \sQuote{transform} (default TRUE) denoting whether to transform the
    value to the distribution domain bounds.}
    \item{skew}{\code{signature(object = "ACDpath")}: Extracts the conditional 
    skew simulated values as a matrix of size n.sim x m.sim. Takes optional
    argument \sQuote{transform} (default TRUE) denoting whether to transform the
    value to the distribution domain bounds. }
    \item{skewness}{\code{signature(object = "ACDpath")}: 
    Returns the conditional simulated Skewness values as a matrix of 
    size n.sim x m.sim}
    }
}
\author{Alexios Ghalanos}
\keyword{classes}