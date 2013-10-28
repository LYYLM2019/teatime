\name{acdmultispec-methods}
\docType{methods}
\alias{acdmultispec}
\alias{acdmultispec-methods}
\alias{acdmultispec,ANY-method}
\alias{acdmultispec,vector-method}
\title{function: Univariate multiple ACD Specification}
\description{
Method for creating a univariate multiple ACD pecification object 
prior to fitting.
}
\usage{
acdmultispec( speclist )
}
\arguments{
\item{speclist}{
A list with as many univariate ACD specifications of class 
\code{\linkS4class{ACDspec}} as there will be columns in the data object passed 
to one of the other methods which uses a multiple specification object 
(fitting, filtering and forecasting).}
}
\value{
A \code{\linkS4class{ACDmultispec}} 
object containing details of the multiple ACD specifications.
}
\author{Alexios Ghalanos}
\keyword{methods}