#################################################################################
##
##   R package dbm by Alexios Ghalanos Copyright (C) 2013.
##   This file is part of the R package dbm.
##
##   The R package dbm is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package dbm is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################

# Taken from the rugarch package
################################################################################
move = function(index, by=1){
	if(!is(index, "POSIXct") && !is(index, "Date") && !is(index, "numeric"))
		stop("\nunrecongnized time index")
	n = length(index)
	if(by>=n) stop("\nby must be less than length of index!")
	if(by==0) return(index)
	p = median(diff(index))
	newindex = c(index[-c(1:by)], generatefwd(index[n], length.out=1, by=p))
	return(newindex)
}

generatefwd = function(T0, length.out = 1, by = "days"){
	if(!is(T0, "POSIXct") && !is(T0, "Date") && !is(T0, "numeric"))
		stop("\nunrecongnized time index")
	Z = seq(T0, by = by, length.out=length.out*4)[-1]
	if(!is.numeric(T0)){
		W = weekdays(Z)
		idx = c(which(W=="Saturday"), which(W=="Sunday"))
		if(length(idx)>0) Z = Z[-idx]
	}
	Z = Z[1:length.out]
	return(Z)
}

# T0 = as.POSIXct("2001-01-01 16:00:00")
# interval = format(seq(as.POSIXct("2001-01-01 09:30:00"), as.POSIXct("2001-01-01 16:00:00"), by="min"), "%H:%M:%S")
# by = "mins"
# length.out=1000
ftseq = function(T0, length.out, by, interval, exclude.weekends = TRUE)
{
	start = T0
	# Just one check:
	if(!is(start, "POSIXct")) stop("\nstart must be a POSIXct object")
	U = format(start, "%H:%M:%S")
	if(is.na(match(U, interval))) stop("\nstart must match one of the supplied interval values.")
	zn = length.out-1
	k = 1
	while(zn<length.out){
		# function does not know what increment "by" is nor the interval.
		# start out with length.out*10 as an estimate and continue looping
		# until we have enough points to satisfy the requirements, including
		# exclusion of weekends
		z = seq(from = start, by = by, length.out = length.out*10*k)[-1]
		y = format(z, "%H:%M:%S")
		z = z[which(!is.na(match(y, interval)))]
		if(exclude.weekends){	
			z = z[-which(weekdays(z)=="Saturday")]
			z = z[-which(weekdays(z)=="Sunday")]
		}
		zn = length(z)
		k = k*2
	}
	return(z[1:length.out])
}

#R = ftseq(as.POSIXct("2001-01-01 16:00:00"), length.out=2000, by = 18000, interval = interval)
################################################################################
.extractdata = function(data, warn = FALSE)
{
	xdata = try(as.xts(data), silent = TRUE)
	if(inherits(xdata, "try-error")){
		if(warn) warning("\ndbm-->warning: data indexing not recognized by xts...coercing to Date with origin 1970-01-01.")
		if(is.data.frame(data) | is.matrix(data)) data = as.numeric(data[,1]) else data = as.numeric(data)
		data = unname(data)
		xdata = xts(data, as.POSIXct(as.Date(seq_along(data), origin="1970-01-01")))
	}
	obj = list()
	obj$data = as.numeric(coredata(xdata))
	obj$index = index(xdata)
	obj$period = median(diff(index(xdata)))
	return(obj)
}

.numeric2xts = function(data){
	data = as.numeric(data)
	return(xts(data, as.Date(1:NROW(data), origin="1950-01-01")))
}

.genxts = function(index0, length.out = 10, period = "days"){
	Z = seq(index0, by = period, length.out=length.out)
	return(Z)
}


logtransform = function(x, lower, upper, inverse = FALSE) 
{
	if(!inverse) {
		ans = lower + (upper - lower)/(1 + exp(-1 * x))
	}
	else {
		ans = -1 * log(-(upper - x)/(-x + lower))
	}
	return(ans)
}

minmaxtransform = function(x, inverse = FALSE, min.x = NULL, max.x = NULL)
{
	if(inverse){
		if(is.null(min.x)) stop("\ncannot inverse without min.x!")
		if(is.null(max.x)) stop("\ncannot inverse without max.x!")
		ans = (x * (max.x - min.x)) + min.x
	} else{
		min.x = min(x, na.rm = TRUE)
		max.x = max(x, na.rm = TRUE)
		ans = (x - min.x)/(max.x - min.x)
	}
	return(ans)
}

iqrtransform = function(x, type = 7, inverse = FALSE, median.x = FALSE, IQR.x = NULL)
{
	if(inverse){
		if(is.null(median.x)) stop("\ncannot inverse without median.x!")
		if(is.null(IQR.x)) stop("\ncannot inverse without IQR.x!")
		ans = (x * IQR.x) + median.x
	} else{
		ans = (x - median(x, na.rm = TRUE))/IQR(x, na.rm = TRUE, type = type)
	}
	return(ans)
}

.sqrtm = function (x) 
{
    tmp = svd(x)
    sqrtx = tmp$u %*% sqrt(diag(tmp$d)) %*% t(tmp$u)
    return(sqrtx)
}

shadeplot = function(signal, series, signal.col = "WhiteSmoke", 
		series.col = "steelblue", main = "", ylim = c(min(series), max(series)), 
		...)
{
	UseMethod("shadeplot")
}

shadeplot.xts = function(signal, series, signal.col = "WhiteSmoke", 
		series.col = "steelblue", main = "", ylim = c(min(series), max(series)), 
		...)
{
	ep <- axTicksByTime(index(series))
	par(mar = c(2.5, 2.5, 2, 1))
	plot(as.numeric(series[,1]), type="l", xaxt = "n", ylab="",xlab="", main = main, ylim = ylim, yaxs = "i")
    signal = signal[index(series)]
    start = which(diff(signal)==1)
    end = which(diff(signal)==(-1))-1
	if(length(start)==0 && length(end)>0) start=1	
    if(length(end)<length(start)) end = c(end, length(signal))
    if(length(end)>length(start)) start = c(1, start)
    miny = min(series, na.rm=TRUE)
    maxy = max(series, na.rm=TRUE)
    n = length(start)
    for(i in 1:n){
		rect(start[i], miny, end[i],  maxy, col= signal.col, border = FALSE)
    }
	lines(as.numeric(series), col = series.col, ...)
	axis(1, at = ep, labels = names(ep), tick = TRUE)
	box()
	par(mar = c(5, 4, 4, 2) + 0.1)
	return(invisible())
}