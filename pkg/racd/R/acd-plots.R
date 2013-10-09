#################################################################################
##
##   R package racd by Alexios Ghalanos Copyright (C) 2012, 2013 
##   This file is part of the R package racd.
##
##   The R package racd is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package racd is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################

.plotacdfit = function(x, which="ask",...)
{
	dist = x@model$dmodel$model
	D = .DistributionBounds(dist)
	if(D$include.shape) sh = 1 else sh=0
	if(D$include.skew ) sk = 1 else sk=0
	if(sum(x@model$dmodel$skewOrder )>0) tsk = 1 else tsk = 0
	if(sum(x@model$dmodel$shapeOrder)>0) tsh = 1 else tsh = 0
	
	P = 4 + sh + sk + tsk + tsh
	choices = c(
			"Series with 1% Quantiles",
			"Conditional SD",
			if(tsk==1) "Conditional Skew",
			if(tsh==1) "Conditional Shape",
			if(sk==1) "Conditional Skewness",
			if(sh==1) "Conditional Kurtosis(ex)",
			"ACF of Standardized Residuals",
			"ACF of Squared Standardized Residuals")
	cp = c(1,2, if(tsk==1) 3 else NULL, if(tsh==1) 4 else NULL, if(sk==1) 5 else NULL, if(sh==1) 6 else NULL, 7:8)
	.interacdfitPlot(x, choices = choices, P = P, plotFUN = paste(".plot.acdfit", cp, sep = "."), which = which, ...)
	# Return Value:
	invisible(x)
}

.interacdfitPlot = function(x, choices, P, plotFUN, which, ...)
{
	if (is.numeric(which)) {
		if(which>length(choices)) stop("Not a valid choice.\n",call. = FALSE)
		FUN = match.fun(plotFUN[which])
		FUN(x)
	}
	if(is.character(which))
	{
		if(which!="all" && which!="ask") stop("Not a valid choice.\n",call. = FALSE)
		if(which[1] == "all") {
			old.par <- par(no.readonly = TRUE)
			#Which = rep(TRUE, times = length(choices))
			par(mfrow=c(3,3))
			for(i in 1:P){
				FUN = match.fun(plotFUN[i])
				FUN(x)
			}
			par(old.par)
		} else{
			.multacdfitPlot(x, choices, ...)
		}
	}
	invisible(x)
}

.multacdfitPlot = function(x, choices, ...)
{
	dist = x@model$dmodel$model
	D = .DistributionBounds(dist)
	if(D$include.shape) sh = 1 else sh=0
	if(D$include.skew ) sk = 1 else sk=0
	if(sum(x@model$dmodel$skewOrder )>0) tsk = 1 else tsk = 0
	if(sum(x@model$dmodel$shapeOrder)>0) tsh = 1 else tsh = 0
	pick = 1
	while (pick > 0) {
		pick = menu (
				choices = paste(" ", choices),
				title = "\nMake a plot selection (or 0 to exit):")
		switch (pick,
				.plot.acdfit.1(x, ...),  .plot.acdfit.2(x, ...),  if(tsk==1) .plot.acdfit.3(x, ...) else NULL,
				if(tsh==1) .plot.acdfit.4(x, ...) else NULL,  if(sk==1) .plot.acdfit.5(x, ...) else NULL,  
				if(sh==1) .plot.acdfit.6(x, ...) else NULL, .plot.acdfit.7(x, ...),  .plot.acdfit.8(x, ...))
	}
}

# Series with 2.5% VaR Limits
.plot.acdfit.1 = function(x, ...)
{
	vmodel  = x@model$vmodel$model
	T = x@model$modeldata$T
	insample = 1:T
	xseries = x@model$modeldata$data[insample]
	xdates  = x@model$modeldata$index[insample]
	distribution = x@model$dmodel$model
	xcmu 	= fitted(x)
	xsigma 	= sigma(x)
	tskew   = skew(x)
	tshape  = shape(x)
	if(distribution == "ghyp") ghlambda = coef(x)["ghlambda"] else ghlambda = 0
	cat("\nplease wait...calculating quantiles...\n")
	q025 	= xcmu + xsigma * qdist(distribution, 0.01, 0, 1, lambda = ghlambda, skew = tskew, shape = tshape)
	q975 	= xcmu + xsigma * qdist(distribution, 0.99, 0, 1, lambda = ghlambda, skew = tskew, shape = tshape)
	plot(xdates, xseries, type = "l", col = "steelblue", ylab = "Returns", xlab="Time", main = "Series with with 2.5% VaR Limits", cex.main = 0.8)
	lines(xdates, q025, col = "tomato1")
	lines(xdates, q975, col = "green")
	mtext(paste("racd"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.6)
	abline(h = 0, col = "grey", lty = 3)
	grid()
}

# Conditional SD
.plot.acdfit.2 = function(x, ...)
{
	vmodel  = x@model$vmodel$model
	T = x@model$modeldata$T
	insample = 1:T
	xdates  = x@model$modeldata$index[insample]
	y 	= sigma(x)
	plot(xdates, y, type = "l", col = "steelblue", ylab = "Volatility", xlab="Time", main = "Conditional SD", cex.main = 0.8)
	mtext(paste("racd"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.6)
	abline(h = 0, col = "grey", lty = 3)
	grid()
}

# Conditional Skew
.plot.acdfit.3 = function(x, ...)
{
	distribution = x@model$dmodel$model
	skmodel  = x@model$dmodel$skewmodel
	T = x@model$modeldata$T
	insample = 1:T
	xdates  = x@model$modeldata$index[insample]
	y 	= skew(x)
	plot(xdates, y, type = "l", col = "steelblue", ylab = "Skew", xlab="Time", cex.lab = 0.9, cex.axis = 0.9,
			main = paste("Conditional Skew [", toupper(distribution),"]", sep=""), cex.main = 0.8)
	mtext(paste("racd"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.6)
	abline(h = 0, col = "grey", lty = 3)
	grid()
}

# Conditional Shape
.plot.acdfit.4 = function(x, ...)
{
	distribution = x@model$dmodel$model
	shmodel  = x@model$dmodel$shapemodel
	T = x@model$modeldata$T
	insample = 1:T
	xdates  = x@model$modeldata$index[insample]
	y 	= shape(x)
	plot(xdates, y, type = "l", col = "steelblue", ylab = "Shape", xlab="Time", cex.lab = 0.9, cex.axis = 0.9,
			main = paste("Conditional Shape [", toupper(distribution),"]", sep=""), cex.main = 0.8)
	mtext(paste("racd"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.6)
	abline(h = 0, col = "grey", lty = 3)
	grid()
}

# Conditional Skewness
.plot.acdfit.5 = function(x, ...)
{
	distribution = x@model$dmodel$model
	T = x@model$modeldata$T
	insample = 1:T
	xdates  = x@model$modeldata$index[insample]
	sk = skew(x) 
	sh = shape(x)
	S = dskewness(distribution, skew = sk, shape = sh)
	plot(xdates, S, type = "l", col = "steelblue", ylab = "Skewness", xlab="Time", cex.lab = 0.9, cex.axis = 0.9,
			main = paste("Conditional Skewness [", toupper(distribution),"]", sep=""), cex.main = 0.8)
	abline(h = 0, col = "grey", lty = 3)
	mtext(paste("racd"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.6)
	grid()
}

# Conditional Kurtosis
.plot.acdfit.6 = function(x, ...)
{
	distribution = x@model$dmodel$model
	T = x@model$modeldata$T
	insample = 1:T
	xdates  = x@model$modeldata$index[insample]
	sk = skew(x) 
	sh = shape(x)
	K = dkurtosis(distribution, skew = sk, shape = sh)
	plot(xdates, K, type = "l", col = "steelblue", ylab = "Kurtosis (ex)", xlab="Time", cex.lab = 0.9, cex.axis = 0.9,
			main = paste("Conditional Kurtosis [", toupper(distribution),"]", sep=""), cex.main = 0.8)
	abline(h = 0, col = "grey", lty = 3)
	mtext(paste("racd"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.6)
	grid()
}


# ACF of standardized residuals
.plot.acdfit.7 = function(x, ...)
{
	zseries = residuals(x, standardize = TRUE)
	zseries[is.na(zseries)] = 0
	n 		= length(zseries)
	lag.max = as.integer(10*log10(n))
	acfx 	= acf(zseries, lag.max = lag.max, plot = FALSE)
	clim0	= qnorm((1 + 0.95)/2)/sqrt(acfx$n.used)
	ylim 	= range(c(-clim0, clim0, as.numeric(acfx$acf)[-1]))
	clx 	= vector(mode = "character", length = lag.max)
	clx[which(as.numeric(acfx$acf)[-1]>=0)] = "steelblue"
	clx[which(as.numeric(acfx$acf)[-1]<0)] = "orange"
	barplot(height = as.numeric(acfx$acf)[-1], names.arg = as.numeric(acfx$lag)[-1], ylim=1.2*ylim, col = clx,
			ylab = "ACF", xlab="lag", main = "ACF of Standardized Residuals", cex.main = 0.8)
	abline(h = c(clim0, -clim0), col = "tomato1", lty = 2)
	abline(h = 0, col = "black", lty = 1)
	box()
	mtext(paste("racd"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.6)
	grid()
}

# ACF of squared standardized residuals
.plot.acdfit.8 = function(x, ...)
{
	zseries = residuals(x, standardize = TRUE)
	zseries[is.na(zseries)] = 0
	n 		= length(zseries)
	lag.max = as.integer(10*log10(n))
	acfx 	= acf(zseries^2, lag.max = lag.max, plot = FALSE)
	clim0	= qnorm((1 + 0.95)/2)/sqrt(acfx$n.used)
	ylim 	= range(c(-clim0, clim0, as.numeric(acfx$acf)[-1]))
	clx 	= vector(mode="character",length=lag.max)
	clx[which(as.numeric(acfx$acf)[-1]>=0)]="steelblue"
	clx[which(as.numeric(acfx$acf)[-1]<0)]="orange"
	barplot(height = as.numeric(acfx$acf)[-1], names.arg = as.numeric(acfx$lag)[-1], ylim=1.2*ylim, col = clx,
			ylab = "ACF", xlab="lag", main = "ACF of Squared Standardized Residuals", cex.main = 0.8)
	abline(h = c(clim0, -clim0), col = "tomato1", lty = 2)
	abline(h = 0, col = "black", lty = 1)
	box()
	mtext(paste("racd"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.6)
	grid()
}