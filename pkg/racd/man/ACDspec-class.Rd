\name{ACDspec-class}
\docType{class}
\alias{ACDspec-class}
\alias{setbounds<-,ACDspec,vector-method}
\alias{setfixed<-,ACDspec,vector-method}
\alias{setstart<-,ACDspec,vector-method}
\alias{show,ACDspec-method}
\title{class: ACD Specification Class}
\description{Class for the ACD specification.}
\section{Slots}{
  \describe{
    \item{\code{model}:}{Object of class \code{"vector"}
     The model specification common to all objects.}
    }
}
\section{Methods}{
  \describe{
	\item{show}{\code{signature(object = "ACDspec")}: 
	Specification summary.}
	\item{setfixed<-}{\code{signature(object = "ACDspec", value = "vector")}: 
	Sets the fixed parameters (which must be supplied as a named list). }
	\item{setstart<-}{\code{signature(object = "ACDspec", value = "vector")}: 
	Sets the starting parameters (which must be supplied as a named list). }
	\item{setbounds<-}{\code{signature(object = "ACDspec", value = "vector")}: 
	Sets the parameters lower and upper bounds, which must be supplied as a named 
	list with each parameter being a numeric vector of length 2 i.e. 
	"alpha1"=c(0,1)). If the vector is of length 1, then this is assumed to be 
	the lower bound, and the upper bound will be set to its default value prior
	to estimation. Setting the upper and lower bounds on the shape or skew 
	parameter, when these are time varying enables to set the bounds on the
	transformed higher moment dynamics.}
	}
}
\author{Alexios Ghalanos}
\seealso{
Classes \code{\linkS4class{ACDfit}}, \code{\linkS4class{ACDsim}} and 
\code{\linkS4class{ACDforecast}}.
}
\keyword{classes}