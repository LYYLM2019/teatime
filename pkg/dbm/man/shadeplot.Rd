\name{shadeplot.xts}
\alias{shadeplot.xts}
\alias{shadeplot}
\title{
Signal Overlay Shade Plot
}
\description{
A simple plot function to provide shaded rectangular regions to binary 
signals.
}
\usage{
\method{shadeplot}{xts}(signal, series, signal.col = "WhiteSmoke", 
series.col = "steelblue", main = "", ylim = c(min(series), max(series)), 
cex.main = 1, cex.axis = 1, ...)
}
\arguments{
\item{signal}{
An xts signal object with the shading denoted by 1, and non-shading denoted by
0.
}
\item{series}{
An xts object on which shading overlay is to be performed given the signal.
}
\item{signal.col}{
The shaded region color.
}
\item{series.col}{
The color of the underlying series.
}
\item{main}{
Plot title.
}
\item{ylim}{
The upper and lower y-axis limits.
}
\item{cex.main}{
Title magnification.
}
\item{cex.axis}{
Axis magnification.
}
\item{\dots}{
Additional arguments passed to the \code{lines} function
}
}
\details{
The function takes an xts object and returns a plot with the original series 
shaded by rectangles for the required periods where the signal was equal to 1.
}
\author{
Alexios Ghalanos
}