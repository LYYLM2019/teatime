\name{logtransform}
\alias{logtransform}
\alias{exptransform}
\title{
  Bounding Transformations
}
\description{
The logistic and exponential transformation functions and their inverses used
in the ACD package.
}
\usage{
logtransform(x, lower, upper, inverse = FALSE)
exptransform(x, lower, upper, rate = 1, inverse = FALSE)
}
\arguments{
\item{x}{A numeric value to be transformed.}
\item{lower}{The lower bound.}
\item{upper}{The upper bound.}
\item{inverse}{(logical) Whether to return the inverse transformation (in which 
case x is the already transformed value).}
\item{rate}{For the exponential transformation, a scaling value.}
}
\details{
The logistic transformation is used for constraining the skew dynamics whilst the
exponential transformation is used for the shape dynamics. It may be of use when
passing the skew0 and shape0 parameters in the \code{acdfit} and \code{acdfilter}
methods.
}
\value{
A named numeric vector of the transformed values, the same length as the input
vector.
}
\author{Alexios Ghalanos}