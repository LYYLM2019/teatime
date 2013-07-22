\name{distbounds}
\alias{distbounds}
\title{
  Conditional Distribution Bounds
}
\description{
Returns the default distribution bounds used by the package.
}
\usage{
distbounds(distribution)
}
\arguments{
\item{distribution}{
A character string of a valid conditional distribution (i.e \dQuote{snorm}, 
\dQuote{std}, \dQuote{sstd}, \dQuote{sged}, \dQuote{ghyp}, \dQuote{nig},
\dQuote{ghst} and \dQuote{jsu}).}
}
\details{
This is a useful utility in helping to decide changes to the skew and shape 
bounds during the ACD estimation (by using the \code{setbounds<-} method.
}
\value{
A named numeric vector of the skew and shape parameters upper and lower bounds.
}
\author{Alexios Ghalanos}