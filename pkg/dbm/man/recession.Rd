\name{recession}
\docType{data}
\alias{recession}
\title{A mulivariate xts dataset for recession forecasting}
\description{
The monthly dataset contains the official NBER recession dates (\sQuote{NBER}) 
from the Federal Reserve Economic Data site (FRED), the OECD composite leading 
indicators recession indicator (\sQuote{OECDR}) based on point in time release 
(i.e. lagged by at least 3 months) and the OECD indicator without lagging 
corresponding to the actual recession occurrence (\sQuote{OECDT}) from the
OECD statistical data website, the 10 year Treasury (\sQuote{T10Y}) and 3 Month
T-Bill (\sQuote{TB3M}) from FRED with data appropriately lagged based on the
release H15 schedule, and the CPI point in time release value from the Archival
Federal Reserve Economic Data site (ALFRED) after adjusting for numerous historic
index re-bases.
}
\usage{data(recession)}
\format{An xts matrix with 618x6 observations.}
\source{FRED,ALFRED,OECD}
\keyword{datasets}
