#################################################################################
##
##   R package rugarch by Alexios Ghalanos Copyright (C) 2008-2013.
##   This file is part of the R package rugarch.
##
##   The R package rugarch is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package rugarch is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################

.multispecall = function( speclist ){
	model = unlist( strsplit(class(speclist[[1]]), "spec") )
	if( model == "ACD" ){
		ans = .multispecacd( speclist )
	}
	return( ans )
}

.multispecacd = function( speclist )
{
	# first create a spec which goes through validation process
	tp = 1
	if( !all(unlist(lapply(speclist, FUN = function(x) is(x, "ACDspec"))) ) ){
		stop("\nNot a valid list of univariate ACD specs.")
	}
	# then check type
	n = length(speclist)
	for(i in 2:n){
		modelnames1 = rownames( speclist[[i]]@model$pars[speclist[[i]]@model$pars[,3]==1, ] )
		modelnames2 = rownames( speclist[[i-1]]@model$pars[speclist[[i-1]]@model$pars[,3]==1, ] )
		if(length(modelnames1) != length(modelnames2))
		{
			tp  = 0
			break()
		} else{
			if(!all(modelnames1 == modelnames2))
			{
				tp  = 0
				break()
			}
		}
	}
	if(tp) type = "equal" else type = "unequal"
	ans = new("ACDmultispec",
			spec = speclist,
			type = type)
	return(ans)
}

# a multifit function possible utilizing parallel execution returning a fitlist
# object
.multifitacd = function(multispec, data, out.sample = 0, solver = "solnp", 
		solver.control = list(), fit.control = list(stationarity = 1, scale = 0, n.sim=2000), 
		cluster = NULL, ...)
{
	n = length(multispec@spec)
	if(is.null(data)) stop("\nmultifit ACD-->error: multifit requires a data object", call. = FALSE)
	if(!is.matrix(data) & !is.data.frame(data) & !is.xts(data)) 
		stop("\nmultifit ACD-->error: multifit only supports matrix, data.frame or xts objects for the data", call. = FALSE)
	asset.names = colnames(data)
	if(dim(data)[2] != n)
		stop("\nmultifit ACD-->error: speclist length not equal to data length", call. = FALSE)
	fitlist = vector(mode = "list", length = n)
	if(length(out.sample) == 1 | length(out.sample) < n) out.sample = rep(out.sample, n)
	##################
	# Parallel Execution Prelim Checks
	if( !is.null(cluster) ){
		clusterEvalQ(cluster, library(rugarch))
		clusterExport(cluster, c("multispec", "data", "out.sample", "solver", 
							"solver.control", "fit.control"), envir = environment())
			fitlist = parLapply(cluster, as.list(1:n), fun = function(i){
						acdfit(spec = multispec@spec[[i]], data = data[, i, drop = FALSE], 
								out.sample = out.sample[i], solver = solver, 
								solver.control = solver.control, fit.control = fit.control)
					})
		} else{
			fitlist = lapply(as.list(1:n), FUN = function(i){
						acdfit(spec = multispec@spec[[i]], data = data[, i, drop = FALSE], 
								out.sample = out.sample[i], solver = solver, 
								solver.control = solver.control, fit.control = fit.control)
					})
	}
	# converged: print
	desc = list()
	desc$type = multispec@type
	desc$asset.names = asset.names
	ans = new("ACDmultifit",
			fit  = fitlist,
			desc = desc)
	return(ans)
}

.multifilteracd1 = function(multifitORspec, data = NULL, out.sample = 0, 
		n.old = NULL, cluster = NULL, ...)
{
	fitlist = multifitORspec
	n = length(fitlist@fit)
	if(is.null(data)) 
		stop("\nmultifilter ACD-->error: multifilter requires a data object", call. = FALSE)
	if(!is.matrix(data) & !is.data.frame(data) & !xts(data)) 
		stop("\nmultifilter ACD-->error: multifilter only supports matrix, data.frame and xts objects for the data", call. = FALSE)
	if(dim(data)[2] != n)
		stop("\nmultifilter ACD-->error: fitlist length not equal to data length", call. = FALSE)
	if(is.matrix(data)) data = as.data.frame(data)
	asset.names = colnames(data)
	filterlist = vector(mode = "list", length = n)
	if(length(out.sample) == 1 | length(out.sample) < n) out.sample = rep(out.sample, n)	
	specx = vector(mode = "list", length = n)
	for(i in 1:n){
		specx[[i]] = getspec(fitlist@fit[[i]])
		specx[[i]]@model$fixed.pars = as.list(coef(fitlist@fit[[i]]))
	}
	if( !is.null(cluster) ){
		clusterEvalQ(cluster, library(rugarch))
		clusterExport(cluster, c("specx", "data", "out.sample", "n.old"), envir = environment())
		filterlist = parLapply(cluster, as.list(1:n), fun = function(i){
					acdfilter(data = data[, i, drop = FALSE], spec = specx[[i]], 
							out.sample =  out.sample[i], n.old = n.old)
				})
	} else{
		filterlist = lapply(as.list(1:n), FUN = function(i){
					acdfilter(data = data[, i, drop = FALSE], spec = specx[[i]], 
							out.sample =  out.sample[i], n.old = n.old)
				})
	}
	desc = list()
	desc$type = "equal"
	desc$asset.names = asset.names
	
	ans = new("ACDmultifilter",
			filter = filterlist,
			desc = desc)
	return(ans)
}

.multifilteracd2 = function(multifitORspec, data = NULL, out.sample = 0, 
		n.old = NULL, cluster = NULL, ...)
{
	speclist = multifitORspec
	n = length(speclist@spec)
	if(is.null(data)) 
		stop("\nmultifilter ACD-->error: multifilter requires a data object", call. = FALSE)
	if(!is.matrix(data) & !is.data.frame(data) & !xts(data)) 
		stop("\nmultifilter ACD-->error: multifilter only supports matrix, data.frame and xts objects for the data", call. = FALSE)
	if(dim(data)[2] != n)
		stop("\nmultifilter ACD-->error: multispec length not equal to data length", call. = FALSE)
	if(is.matrix(data)) data = as.data.frame(data)
	asset.names = colnames(data)
	filterlist = vector(mode = "list", length = n)
	if(length(out.sample) == 1 | length(out.sample) < n) out.sample = rep(out.sample, n)
	
	if( !is.null(cluster) ){
		clusterEvalQ(cluster, library(rugarch))
		clusterExport(cluster, c("speclist", "data", "out.sample", "n.old"), envir = environment())
		filterlist = parLapply(cluster, as.list(1:n), fun = function(i){
					ugarchfilter(data = data[, i, drop = FALSE], spec = speclist@spec[[i]], 
							out.sample =  out.sample[i], n.old = n.old)
				})
	} else{
		filterlist = lapply(as.list(1:n), FUN = function(i){
					ugarchfilter(data = data[, i, drop = FALSE], spec = speclist@spec[[i]], 
							out.sample =  out.sample[i], n.old = n.old)
				})	
	}
	# converged: print
	desc = list()
	desc$type = speclist@type
	desc$asset.names = asset.names
	
	ans = new("ACDmultifilter",
			filter  = filterlist,
			desc = desc)
	return(ans)
}

.multiforecastacd1 = function(multifitORspec, data = NULL, n.ahead = 1, 
		n.roll = 0, out.sample = 0,external.forecasts = list(mregfor = NULL, 
				vregfor = NULL, skregfor = NULL, shregfor = NULL), 
		m.sim = 1000, cluster = NULL, ...)
{
	multifit = multifitORspec
	n = length(multifit@fit)
	asset.names = multifit@desc$asset.names
	forecastlist = vector(mode = "list", length = n)
	if( !is.null(cluster) ){
		clusterEvalQ(cluster, library(rugarch))
		clusterExport(cluster, c("multifit", "n.ahead", "n.roll", "external.forecasts", "m.sim"), envir = environment())
		forecastlist = parLapply(cluster, as.list(1:n), fun = function(i){
					acdforecast(fitORspec = multifit@fit[[i]], n.ahead = n.ahead, 
							n.roll = n.roll, external.forecasts = external.forecasts, m.sim = m.sim)
				})
	} else{
		forecastlist = lapply(as.list(1:n), FUN = function(i){
					acdforecast(fitORspec = multifit@fit[[i]], n.ahead = n.ahead, 
							n.roll = n.roll, external.forecasts = external.forecasts, m.sim = m.sim)
				})
	}
	desc = list()
	desc$type = "equal"
	desc$asset.names = asset.names
	ans = new("ACDmultiforecast",
			forecast  = forecastlist,
			desc = desc)
	return(ans)
}