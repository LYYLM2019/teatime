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
.mcsacdfit = function(spec,  data, solver = "ucminf", out.sample = 0, solver.control = list(), 
		fit.control = list(stationarity = 0, fixed.se = 0, scale = 0, n.sim = 2000), 
		skew0 = NULL, shape0 = NULL, cluster = NULL, DailyVar, ...)
{
	tic = Sys.time()
	vmodel = spec@model$vmodel$model
	if(is.null(solver.control$trace)) trace = 0 else trace = solver.control$trace
	# default for stationarity is off for ACD models
	if(is.null(fit.control$stationarity)) fit.control$stationarity = FALSE
	if(is.null(fit.control$fixed.se)) fit.control$fixed.se = FALSE
	if(is.null(fit.control$scale)){
		fit.control$scale = FALSE
	} else{
		if(fit.control$scale) stop("\nscaling not valid for mcsACD model.")
	}
	if(is.null(fit.control$n.sim)) fit.control$n.sim = 2000
	
	mm = match(names(fit.control), c("stationarity", "fixed.se", "scale", "n.sim"))
	if(any(is.na(mm))){
		idx = which(is.na(mm))
		enx = NULL
		for(i in 1:length(idx)) enx = c(enx, names(fit.control)[idx[i]])
		warning(paste(c("\nunidentified option(s) in fit.control:\n", enx), sep="", collapse=" "), call. = FALSE, domain = NULL)
	}
	# if there are fixed pars we do no allow scaling as there would be no way of mixing scaled
	# amd non scaled parameters	
	if(sum(spec@model$pars[,2]) > 0) fit.control$scale = FALSE
	# if we have arch, skewm or shapem in the model turn off scaling
	if(sum(spec@model$modelinc[5:7])>0) fit.control$scale = FALSE
	if(is.null(DailyVar)){
		stop("\nacdfit-->error: you must supply the daily forecast variance (DailyVar) for the msGARCH model\n")
	} else{
		if(!is(DailyVar, "xts")) stop("\nacdfit-->error: DailyVar must be an xts object\n")
		DailyVarIndex = format(index(DailyVar), format="%Y-%m-%d")
	}
	# we are not going to extract the data just yet
	UIndex = unique(format(index(data), format="%Y-%m-%d"))
	DIndex = format(index(data), format="%Y-%m-%d")
	RIndex = index(data)
	M = length(UIndex)
	matchD = all.equal(UIndex, DailyVarIndex)
	if(!matchD) stop("\nugarchfit-->error: DailyVar dates do not match the data dates (unique days).\n")
	Tb = lapply(1:M, function(i) RIndex[which(DIndex==UIndex[i])])
	DVar = lapply(1:M, function(i) rep(DailyVar[i], length(Tb[[i]])))
	# can't unlist a POSIXct object...need to manually concatentate (can't use 'c' with recursive option either)
	dTT = Tb[[1]]
	if(length(Tb)>1) for(i in 2:length(Tb)) dTT = c(dTT, Tb[[i]])
	DVar = xts(as.numeric(unlist(DVar)), dTT)
	xdata = rugarch:::.extractdata(data)
	if(!is.numeric(out.sample))  stop("\nacdfit-->error: out.sample must be numeric\n")
	if(as.numeric(out.sample)<0) stop("\nacdfit-->error: out.sample must be positive\n")
	n.start = round(out.sample,0)
	n = length(xdata$data)
	if((n-n.start)<100) stop("\nacdfit-->error: function requires at least 100 data\n points to run\n")
	data = xdata$data[1:(n-n.start)]
	index = xdata$index[1:(n-n.start)]
	origdata = xdata$data
	origindex = xdata$index
	period = xdata$period
	# for the mcsACD model we need to work with xts
	data = xts(data, index)
	itime = .unique_intraday(data)
	DV = DVar[1:(n-n.start)]
	idx1 = .unique_time(data)
	idx2 = .stime(data)
	# create a temporary environment to store values (deleted at end of function)
	garchenv = new.env(hash = TRUE)
	arglist = list()
	###################
	# needed for the startpars initialization (for the mcsGARCH model)
	arglist$DailyVar = DailyVar
	# needed for the msGARCH model
	arglist$idx1 = idx1
	arglist$idx2 = idx2
	arglist$DV = as.numeric(DV)
	###################
	arglist$garchenv <- garchenv
	arglist$sbounds = spec@model$sbounds
	arglist$pmode = 0
	model = spec@model
	modelinc = model$modelinc
	pidx = model$pidx
	# expand the spec object and assign spec lists
	if(modelinc[8] > 0){
		arglist$mexdata = model$modeldata$mexdata[1:(n-n.start), , drop = FALSE]
	} else{
		arglist$mexdata = 0
	}
	if(modelinc[17] > 0){
		arglist$vexdata = model$modeldata$vexdata[1:(n-n.start), ,drop = FALSE]
	} else{
		arglist$vexdata = 0
	}
	if(modelinc[25] > 0){
		arglist$skxdata = model$modeldata$skxdata[1:(n-n.start), ,drop = FALSE]
	} else{
		arglist$skxdata = 0
	}
	if(modelinc[31] > 0){
		arglist$shxdata = model$modeldata$shxdata[1:(n-n.start), ,drop = FALSE]
	} else{
		arglist$shxdata = 0
	}
	arglist$index = index
	arglist$trace = trace
	m =  model$maxOrder
	# store length of data for easy retrieval
	model$modeldata$T = T = length(as.numeric(data))
	dist = model$modeldesc$distribution
	# if(fit.control$scale) dscale = sd(data) else 
	dscale = 1
	zdata = data/dscale
	arglist$transform = FALSE
	arglist$fnscale = 1
	arglist$data = zdata
	arglist$dscale = dscale
	arglist$model = model
	arglist$shape0 = shape0
	arglist$skew0 = skew0
	ipars = model$pars
	# Optimization Starting Parameters Vector & Bounds
	tmp = .acdstart(ipars, arglist)
	arglist = tmp$arglist
	ipars = arglist$ipars = tmp$pars
	# arglist$skhEst
	arglist$tmph  = tmp$tmph
	arglist$model = model
	# we now split out any fixed parameters
	estidx = as.logical( ipars[,4] )
	arglist$estidx = estidx	
	arglist$fit.control = fit.control
	npars = sum(estidx)
	if(any(ipars[,2]==1)){
		if(npars == 0){
			if(fit.control$fixed.se==0) {
				# if all parameters are fixed an no standard erros are to
				# be calculated then we return a filter object
				cat("\nacdfit-->warning: all parameters fixed...returning ACDfilter object instead\n")
				return(acdfilter(data = xts(origdata, origindex), spec = spec, out.sample = out.sample, DailyVar = DailyVar))
			} else{
				# if all parameters are fixed but we require standard errors, we
				# skip the solver
				use.solver = 0
				ipars[ipars[,2]==1, 4] = 1
				ipars[ipars[,2]==1, 2] = 0
				arglist$ipars = ipars
				estidx = as.logical( ipars[,4] )
				arglist$estidx = estidx
			}
		} else{
			# with some parameters fixed we extract them (to be rejoined at end)
			# so that they do not enter the solver
			use.solver = 1
		}
	} else{
		use.solver = 1
	}
	# start counter
	assign("racd_llh", 1, envir = garchenv)
	arglist$fit.control = fit.control
	
	fun = racd:::.mcsacdLLH
	fname = "mcsACD"
	
	if(use.solver){
		parscale = rep(1, length = npars)
		names(parscale) = rownames(ipars[estidx,])
		if(modelinc[1] > 0) parscale["mu"] = abs(mean(zdata))
		if(modelinc[9] > 0) parscale["omega"] = var(zdata)
		arglist$returnType = "llh"
		solution = .acdsolver(solver, pars = ipars[estidx, 1], 
				fun = fun,	Ifn = NULL, ILB = NULL, IUB = NULL, gr = NULL, 
				hessian = NULL, parscale = parscale, control = solver.control, 
				LB = ipars[estidx, 5], UB = ipars[estidx, 6], cluster = cluster, 
				arglist = arglist)
		#-----------------------------------------------------------------------
		sol = solution$sol
		hess = solution$hess
		timer = Sys.time()-tic
		pars = solution$sol$pars
		if(!is.null(sol$par)){
			ipars[estidx, 1] = sol$par
			if(modelinc[9]==0){
				# call it once more to get omega
				tmpx = fun(sol$par, arglist)
				ipars[pidx["omega",1], 1] = get("omega", garchenv)
			}
			if(sum(ipars[,2]) == 0){
				if(modelinc[1] > 0) ipars[pidx["mu",1]:pidx["mu",2], 1] = ipars[pidx["mu",1]:pidx["mu",2], 1] * dscale
				if(modelinc[8] > 0){
					ipars[pidx["mxreg", 1]:pidx["mxreg", 2], 1] = ipars[pidx["mxreg", 1]:pidx["mxreg", 2], 1] * dscale
				}
				ipars[pidx["omega",1], 1] = ipars[pidx["omega",1],1] * dscale^2
			}
		} else{
			ipars[estidx, 1] = NA
		}
		arglist$ipars = ipars
		convergence = sol$convergence
		if(convergence != 0) warning("\nacdfit-->warning: solver failed to converge.")
	} else{
		solution = NULL
		hess = NULL
		timer = Sys.time()-tic
		convergence = 0
		sol = list()
		sol$message = "all parameters fixed"
	}
	fit = list()
	# check convergence else write message/return
	# create a copy of ipars in case we need to change it below to calculate standard errors
	# which we will need to reset later (because for example, infocriteria uses estimated
	# parameters, not fixed.
	ipars2 = ipars
	if(convergence == 0){
		arglist$dscale  = 1
		arglist$data = data
		if(sum(ipars[,2]) > 0 && fit.control$fixed.se == 1){
			ipars[ipars[,2]==1, 4] = 1
			ipars[ipars[,2]==1, 2] = 0
			arglist$ipars = ipars
			estidx = as.logical( ipars[,4] )
			arglist$estidx = estidx	
		}
		fit = .acdmakefitmodel(acdmodel = fname, f = fun, T = T, m = m, 
				timer = timer, convergence = convergence, message = sol$message, 
				hess, arglist = arglist)
		model$modelinc[7] = modelinc[7]
		model$modeldata$data = origdata
		model$modeldata$index = origindex
		model$modeldata$period = period
		model$pars[, 1] = fit$ipars[,1]
		model$pars[, 5:6] = ipars2[,5:6]
		fit$ipars[, 4] = ipars2[, 4]
		fit$ipars[, 2] = ipars2[, 2]
		fit$ipars[, 5:6] = ipars2[,5:6]
		# make sure omega is now included (for working with object post-estimation)
		fit$ipars["omega", 3] = 1
		model$pars["omega", 3] = 1
		########################################################################
		# V (daily forecast Variance) and S (diurnal Variance) are aligned to the
		# original time index
		# model$idx1 = .unique_time(xts(origdata, origindex))
		# model$idx2 = .stime(xts(origdata, origindex))
		# itime == unique intraday intervals on which diurnal vol is based and for
		# use with forecasting and simulation
		model$dtime = itime
		idx1 = .unique_time(xts(origdata[1:T], origindex[1:T]))
		idx2 = .stime(xts(origdata[1:T], origindex[1:T]))
		model$dvalues = .diurnal_series(fit$residuals, as.numeric(DVar)[1:T], idx1)
		# DailyVar will be of length T+out.sample (since it does not depend on any endogenous
		# variables and we can safely store the full values)
		model$DailyVar = DVar
		# DiurnalVar will be of length T (because it depends on residuals)		
		model$DiurnalVar = xts(.diurnal_series_aligned(xts(fit$residuals, index), DVar[1:T], idx1, idx2), origindex[1:T])
		# adjust sigma (q = component volatility i.e. on deasonalized data)
		fit$q = fit$sigma
		fit$sigma = fit$q * sqrt(DVar[1:T]*model$DiurnalVar[1:T])
		fit$z = fit$residuals/fit$sigma
		fit$skhEst = arglist$skhEst
	} else{
		fit$message = sol$message
		fit$convergence = 1
		fit$skhEst = arglist$skhEst
		model$modeldata$data = origdata
		model$modeldata$index = origindex
		model$modeldata$period = period
	}
	# make model list to return some usefule information which
	# will be called by other functions (show, plot, sim etc)
	model = model
	model$garchLL = get("garchLL", garchenv)
	model$n.start = n.start
	ans = new("ACDfit",
			fit = fit,
			model = model)
	rm(garchenv)
	return(ans)

}
.mcsacdLLH = function(pars, arglist)
{
	if(arglist$transform){ pars = logtransform(pars, arglist$LB, arglist$UB) }
	# prepare inputs
	eps = .Machine$double.eps
	data = arglist$data
	assign("x_pars", pars, envir = arglist$garchenv)
	if(!is.null(arglist$n.old)) Nx = arglist$n.old else Nx = length(data)
	model = arglist$model
	estidx = arglist$estidx
	idx = model$pidx
	ipars = arglist$ipars
	ipars[estidx, 1] = pars
	trace = arglist$trace
	T = length(data)
	fit.control = arglist$fit.control
	m = model$maxOrder
	N = c(m, T)
	distribution = model$modeldesc$distribution
	modelinc = model$modelinc
	# distribution number
	# 1 = norm, 2=snorm, 3=std, 4=sstd, 5=ged, 6=sged, 7=nig
	dist = model$modeldesc$distno
	hm = arglist$tmph
	rx = .arfimaxfilteracd(modelinc, ipars[,1], idx, mexdata = arglist$mexdata, h = hm, 
			tskew = 0, tshape = 0, data = data, N = N, arglist$garchenv)
	res = rx$res
	zrf = rx$zrf
	res[is.na(res) | !is.finite(res) | is.nan(res)] = 0
	# 1. Create the diurnal series (bins)
	# 2. Adjust res by denominator
	dseries = .diurnal_series_aligned(res, arglist$DV, arglist$idx1, arglist$idx2)
	eres = res/sqrt(arglist$DV*dseries)
	if( !is.null(arglist$n.old) ){
		rx = .arfimaxfilteracd(modelinc, ipars[,1], idx, mexdata = arglist$mexdata[1:Nx, , drop=FALSE], 
				h = hm, tskew = 0, tshape = 0, data = data[1:Nx], N = c(m, Nx), arglist$garchenv)
		res2 = rx$res
		res2[is.na(res2) | !is.finite(res2) | is.nan(res2)] = 0
		xdseries = .diurnal_series_aligned(res2, as.numeric(arglist$DV), arglist$idx1, arglist$idx2)
		xeres = res2/sqrt(as.numeric(arglist$DV[1:Nx])*xdseries[1:Nx])		
		mvar = mean(xeres*xeres)
	} else{
		mvar = mean(eres*eres)
	}

	# sgarch persistence value
	mexdata = as.double(as.vector(arglist$mexdata))
	vexdata = as.double(as.vector(arglist$vexdata))
	skxdata = as.double(as.vector(arglist$skxdata))
	shxdata = as.double(as.vector(arglist$shxdata))
	persist = (sum(ipars[idx["alpha",1]:idx["alpha",2],1]) + sum(ipars[idx["beta",1]:idx["beta",2],1]))
	if(modelinc[9]>0){
		ipars[idx["omega",1],1] = max(eps, ipars[idx["omega",1],1]) 
		hEst = mvar
	} else{
		if(modelinc[17]>0) {
			mv = sum(apply(matrix(arglist$vexdata, ncol = modelinc[17]), 2, "mean")*ipars[idx["vxreg",1]:idx["vxreg",2],1])
		} else{
			mv = 0
		}
		ipars[idx["omega",1],1] = mvar * (1 - persist) - mv
		hEst = mvar
		assign("omega", ipars[idx["omega",1],1], arglist$garchenv)
	}
	if(is.na(hEst) | !is.finite(hEst) | is.nan(hEst)) hEst = (1 - persist)
	assign("racd_ipars", ipars, envir = arglist$garchenv)
	if(fit.control$stationarity == 1 && modelinc[17] == 0){
		if(!is.na(persist) && persist >= 1) return(llh = get("racd_llh", arglist$garchenv) + 0.1*(abs(get("racd_llh", arglist$garchenv))))
	}	
	sbounds = model$sbounds
	skhEst = arglist$skhEst
	tempskew  = double(length = T)
	tempshape = double(length = T)
	tskew 	= double(length = T)
	tshape 	= double(length = T)
	h 		= double(length = T)
	z 		= double(length = T)
	constm 	= double(length = T)
	condm 	= double(length = T)
	llh 	= double(length = 1)
	LHT 	= double(length = T)
	
	ans = try(.C("mcsacdfilterC",
					model = as.integer(modelinc), 
					pars = as.double(ipars[,1]), 
					idx = as.integer(idx[,1]-1), 
					hEst = as.double(hEst), 
					res = as.double(res),
					e = as.double(eres*eres),
					eres = as.double(eres),
					s = as.double(dseries), 
					v = as.double(arglist$DV),
					vexdata = as.double(vexdata), 
					m = as.integer(m), 
					T = as.integer(T),
					h = double(T), 
					z = double(T), 
					tempskew = double(T), 
					tempshape = double(T),
					skhEst = as.double(skhEst),
					tskew = double(T),
					tshape = double(T),
					skxreg = as.double(skxdata),
					shxreg = as.double(shxdata),
					sbounds = as.double(sbounds),
					llh = double(1), 
					LHT = double(T),
					PACKAGE = "racd"), silent = TRUE )
	
	if( inherits(ans, "try-error") ){
		cat(paste("\nacdfit-->warning: ", ans,"\n", sep=""))
		return( llh = get("racd_llh", arglist$garchenv) + 0.1*abs( get("racd_llh", arglist$garchenv) ) )
	}
	
	z = ans$z
	h = ans$h
	res = ans$res
	llh  = ans$llh
	tskew  = ans$tskew
	tshape = ans$tshape
	tempskew = ans$tempskew
	tempshape = ans$tempshape
	
	if( is.finite(llh) && !is.na(llh) && !is.nan(llh) ){
		assign("racd_llh", llh, envir = arglist$garchenv) 
	} else {
		llh = (get("racd_llh", arglist$garchenv) + 100*(abs(get("racd_llh",arglist$garchenv))))
	}
	# LHT = raw scores
	LHT = -ans$LHT
	ans = switch(arglist$returnType,
			llh = arglist$fnscale*llh,
			LHT = LHT,
			all = list(llh = llh, h = h, res = res, z = z, kappa = 1, 
					tskew = tskew, tshape = tshape, tempshape = tempshape, 
					tempskew = tempskew, LHT = LHT, dseries = dseries))
	return( ans )
}



.mcsacdfilter = function(spec, data, out.sample = 0, n.old = NULL, skew0 = NULL, shape0 = NULL, DailyVar, ...)
{
	# n.old is optional and indicates the length of the original dataseries (in
	# cases when this represents a dataseries augmented by newer data). The reason
	# for using this is so that the old and new datasets agree since the original
	# recursion uses the sum of the residuals to start the recursion and therefore
	# is influenced by new data. For a small augmentation the values converge after
	# x periods, but it is sometimes preferable to have this option so that there is
	# no forward looking information contaminating the study.
	tic = Sys.time()
	if(missing(DailyVar)){
		stop("\nacdfilter-->error: you must supply the daily forecast variance (DailyVar) for the msGARCH model\n")
	} else{
		if(!is(DailyVar, "xts")) stop("\nacdfilter-->error: DailyVar must be an xts object\n")
		DailyVarIndex = format(index(DailyVar), format="%Y-%m-%d")
	}
	# we are not going to extract the data just yet
	UIndex = unique(format(index(data), format="%Y-%m-%d"))
	DIndex = format(index(data), format="%Y-%m-%d")
	RIndex = index(data)
	M = length(UIndex)
	matchD = all.equal(UIndex, DailyVarIndex)
	if(!matchD) stop("\nacdfilter-->error: DailyVar dates do not match the data dates (unique days).\n")
	Tb = lapply(1:M, function(i) RIndex[which(DIndex==UIndex[i])])
	DVar = lapply(1:M, function(i) rep(DailyVar[i], length(Tb[[i]])))
	# can't unlist a POSIXct object...need to manually concatentate (can't use 'c' with recursive option either)
	dTT = Tb[[1]]
	if(length(Tb)>1) for(i in 2:length(Tb)) dTT = c(dTT, Tb[[i]])
	DVar = xts(as.numeric(unlist(DVar)), dTT)
	vmodel = spec@model$vmodel$model
	xdata = rugarch:::.extractdata(data)
	data = xdata$data
	index = xdata$index
	period = xdata$period
	origdata = data
	origindex = xdata$index
	T = length(origdata)  - out.sample
	if(!is.null(n.old) && n.old>T) stop("\nn.old cannot be greater than length data - out.sample!")	
	if(!is.null(n.old)) Nx = n.old else Nx = length(data)
	
	data = origdata[1:T]
	index = origindex[1:T]
	
	data = xts(data, index)
	itime = .unique_intraday(data)
	DV = DVar[1:T]
	idx1 = .unique_time(data[1:Nx])
	idx2 = .stime(data)
	
	model = spec@model
	ipars = model$pars
	pars = unlist(model$fixed.pars)
	parnames = names(pars)
	modelnames = rugarch:::.checkallfixed(spec)
	if(is.na(all(match(modelnames, parnames), 1:length(modelnames)))) {
		cat("\nacdfilter-->error: parameters names do not match specification\n")
		cat("Expected Parameters are: ")
		cat(paste(modelnames))
		cat("\n")
		stop("Exiting", call. = FALSE)
	}
	# once more into the spec
	# NB Any changes made to the spec are not preserved once we apply set fixed
	setfixed(spec)<-as.list(pars)
	garchenv = new.env(hash = TRUE)
	racd_llh = 1
	assign("racd_llh", 1, envir = garchenv)
	arglist = list()
	###################
	# needed for the msGARCH model
	arglist$idx1 = idx1
	arglist$idx2 = idx2
	arglist$DV = as.numeric(DV)
	###################
	arglist$n.old = n.old
	arglist$transform = FALSE
	arglist$garchenv <- garchenv
	arglist$pmode = 0
	modelinc = model$modelinc
	pidx = model$pidx
	# expand the spec object and assign spec lists
	if(modelinc[8] > 0){
		arglist$mexdata = model$modeldata$mexdata[1:T, , drop = FALSE]
	} else{
		arglist$mexdata = 0
	}
	if(modelinc[17] > 0){
		arglist$vexdata = model$modeldata$vexdata[1:T, ,drop = FALSE]
	} else{
		arglist$vexdata = 0
	}
	if(modelinc[25] > 0){
		arglist$skxdata = model$modeldata$skxdata[1:T, ,drop = FALSE]
	} else{
		arglist$skxdata = 0
	}
	if(modelinc[31] > 0){
		arglist$shxdata = model$modeldata$shxdata[1:T, ,drop = FALSE]
	} else{
		arglist$shxdata = 0
	}
	arglist$index = index
	arglist$trace = 0
	# store length of data for easy retrieval
	arglist$data = data
	arglist$ipars  = ipars
	# starting parameter for skew+shape
	arglist$skhEst = rep(0,2)
	if(!is.null(shape0)){
		if(modelinc[27]==1) arglist$skhEst[2] = exptransform(shape0, model$sbounds[3], model$sbounds[4], rate = model$sbounds[5], inverse=TRUE)
	} else{
		if(modelinc[27]==1) arglist$skhEst[2] = pars["shcons"]
	}
	if(!is.null(skew0)){
		if(modelinc[21]==1) arglist$skhEst[1] = logtransform(skew0, model$sbounds[1], model$sbounds[2], inverse=TRUE)
	} else{
		if(modelinc[21]==1) arglist$skhEst[1] = pars["skcons"]
	}
	
	arglist$tmph  = 0
	arglist$model = model
	# we now split out any fixed parameters
	estidx = as.logical( ipars[,3] )
	arglist$estidx = estidx	
	arglist$returnType = "all"
	arglist$fit.control=list(stationarity = 0)
	ans  = .mcsacdLLH(pars, arglist)
	filter = list()
	filter$z = ans$z
	# V (daily forecast Variance) and S (diurnal Variance) are aligned to the
	# original time index
	# set idx1 to only the n.old value.
	model$idx1 = .unique_time(xts(origdata[1:Nx], origindex[1:Nx]))
	model$idx2 = .stime(xts(origdata, origindex))
	# itime == unique intraday intervals on which diurnal vol is based and for
	# use with forecasting and simulation
	itime = .unique_intraday(xts(origdata[1:Nx], origindex[1:Nx]))
	model$itime = itime
	model$DailyVar = DVar
	model$DiurnalVar = xts(.diurnal_series_aligned(ans$res, as.numeric(DVar), model$idx1, model$idx2), origindex)
	filter$q = sqrt(ans$h)
	filter$sigma = filter$q * sqrt(DVar[1:T]*model$DiurnalVar[1:T])
	filter$residuals = ans$res
	filter$z = filter$residuals/filter$sigma
	filter$fitted.values = data - ans$res
	filter$tskew = ans$tskew
	filter$tshape = ans$tshape
	filter$tempskew = ans$tempskew
	filter$tempshape = ans$tempshape
	filter$LLH = -ans$llh
	filter$log.likelihoods = ans$LHT
	filter$ipars = ipars
	filter$skhEst = arglist$skhEst
	model$modeldata$data = origdata
	model$modeldata$index = origindex
	model$modeldata$period = period
	model$modeldata$T = T
	model$n.start = out.sample
	filter$timer = Sys.time() - tic
	sol = new("ACDfilter",
			filter = filter,
			model = model)
	return(sol)
}
.mcsacdforecast = function(fit, n.ahead = 1, n.roll = 0, 
		external.forecasts = list(mregfor = NULL, vregfor = NULL, skxregfor = NULL,
				shxregfor = NULL), m.sim = 1000, cluster = NULL, DailyVar, ...)
{
	data = fit@model$modeldata$data
	Nor = length(as.numeric(data))
	index = fit@model$modeldata$index
	period = fit@model$modeldata$period
	if(n.ahead>1) stop("\nmcsACD model does not currently support multi-step ahead forecast (which is simulation based).")
	# check DailyVar forecast provided with what is available from the fitted object
	# (if out.sample was used it may not be needed).
	DiurnalVar = fit@model$DiurnalVar
	inDailyVar = fit@model$DailyVar
	lastDate = format(tail(index(inDailyVar), 1), "%Y-%m-%d")
	dtime = fit@model$dtime
	dvalues = fit@model$dvalues
	# prepare the diurnal, daily vols
	if(fit@model$n.start>0){
		if((n.ahead+n.roll)<=fit@model$n.start){
			# we don't require external DailyVaR
			# completely in-the-sample
			DVar = fit@model$DailyVar
			DiurnalVar = NULL
		} else{
			# mixed in and out sample
			needT = (n.ahead+n.roll) - fit@model$n.start
			outD = ftseq(T0 = as.POSIXct(tail(index, 1)), 
					length.out = needT, by = period, 
					interval = fit@model$dtime, 
					exclude.weekends = TRUE)
			Dmatch = match(format(outD, "%H:%M:%S"), dtime)
			D2 = xts(dvalues[Dmatch], outD)
			D1 = xts(dvalues[match(format(index, "%H:%M:%S"), dtime)], index)
			DiurnalVar = c(D1, D2)
			DVar = .intraday2daily(fit@model$DailyVar)
			# Check to see whether we need DailyVar Forecast
			U = unique(format(index(DiurnalVar), "%Y-%m-%d"))
			Y = unique(c(format(index(DVar), "%Y-%m-%d"), if(!missing(DailyVar)) format(index(DailyVar),  "%Y-%m-%d") else NULL))
			DV = c(as.numeric(DVar), if(!missing(DailyVar)) as.numeric(DailyVar) else NULL)
			test = match(U, Y)
			if(any(is.na(test))){
				idx = which(is.na(test))
				stop(paste(c("DailyVar requires forecasts for: ", U[idx],"...resubmit."), sep="",collapse=" "))
			} else{
				# create the DailyVar
				M = length(Y)
				RIndex = index(DiurnalVar)
				UIndex = unique(format(index(DiurnalVar), format="%Y-%m-%d"))
				DIndex = format(index(DiurnalVar), format="%Y-%m-%d")
				Tb = lapply(1:M, function(i) RIndex[which(DIndex==UIndex[i])])
				DVar = lapply(1:M, function(i) rep(DV[i], length(Tb[[i]])))
				# can't unlist a POSIXct object...need to manually concatentate (can't use 'c' with recursive option either)
				dTT = Tb[[1]]
				if(length(Tb)>1) for(i in 2:length(Tb)) dTT = c(dTT, Tb[[i]])
				DVar = xts(as.numeric(unlist(DVar)), dTT)
			}
		}
	} else{
		# completely out of the sample
		outD = ftseq(T0 = as.POSIXct(tail(index, 1)), 
				length.out = n.ahead+n.roll, by = period, 
				interval = fit@model$dtime, 
				exclude.weekends = TRUE)
		Dmatch = match(format(outD, "%H:%M:%S"), dtime)
		D2 = xts(dvalues[Dmatch], outD)
		D1 = xts(dvalues[match(format(index, "%H:%M:%S"), dtime)], index)
		DiurnalVar = c(D1, D2)
		DVar =  .intraday2daily(fit@model$DailyVar)
		# Check to see whether we need DailyVar Forecast
		U = unique(format(index(DiurnalVar), "%Y-%m-%d"))
		Y = unique(c(format(index(DVar), "%Y-%m-%d"), if(!missing(DailyVar)) format(index(DailyVar), "%Y-%m-%d") else NULL))
		DV = c(as.numeric(DVar), if(!missing(DailyVar)) as.numeric(DailyVar) else NULL)
		test = match(U, Y)
		if(any(is.na(test))){
			idx = which(is.na(test))
			stop(paste(c("DailyVar requires forecasts for: ", U[idx],"...resubmit."), sep="",collapse=" "))
		} else{
			# create the DailyVar
			M = length(Y)
			RIndex = index(DiurnalVar)
			UIndex = unique(format(index(DiurnalVar), format="%Y-%m-%d"))
			DIndex = format(index(DiurnalVar), format="%Y-%m-%d")
			Tb = lapply(1:M, function(i) RIndex[which(DIndex==UIndex[i])])
			DVar = lapply(1:M, function(i) rep(DV[i], length(Tb[[i]])))
			# can't unlist a POSIXct object...need to manually concatentate (can't use 'c' with recursive option either)
			dTT = Tb[[1]]
			if(length(Tb)>1) for(i in 2:length(Tb)) dTT = c(dTT, Tb[[i]])
			DVar = xts(as.numeric(unlist(DVar)), dTT)
		}
	}
	ns = fit@model$n.start
	N = Nor - ns
	model = fit@model
	ipars = fit@fit$ipars
	modelinc = model$modelinc
	idx = model$pidx
	if( n.roll > ns ) stop("\nacdforecast-->error: n.roll must not be greater than out.sample!")
	pars = fit@fit$coef
	ipars = fit@fit$ipars
	# check if necessary the external regressor forecasts provided first
	xreg = .acdforcregressors(model, external.forecasts$mregfor, 
			external.forecasts$vregfor, external.forecasts$skxregfor, 
			external.forecasts$shxregfor, n.ahead, Nor, out.sample = ns, n.roll)
	mxf = xreg$mxf
	vxf = xreg$vxf
	skxf = xreg$skxf
	shxf = xreg$shxf
	
	# filter data (check external regressor data - must equal length of origData)
	fcreq = ifelse(ns >= (n.ahead+n.roll), n.ahead+n.roll, ns)
	fspec = acdspec(variance.model = list(model = "mcsGARCH", 
					garchOrder = model$vmodel$garchOrder, 
					external.regressors = model$modeldata$vexdata, variance.targeting = FALSE), 
			mean.model = list(armaOrder = model$mmodel$armaOrder, 
					include.mean = model$mmodel$include.mean, archm = as.logical(modelinc[5]), 
					arfima = as.logical(modelinc[4]), external.regressors = model$modeldata$mexdata), 
			distribution.model = list(model = model$dmodel$model, 
					skewOrder = model$dmodel$skewOrder, skewshock = model$dmodel$skewshock, 
					skewmodel = model$dmodel$skewmodel,
					skew.regressors = model$modeldata$skxdata,
					shapeOrder = model$dmodel$shapeOrder, shapeshock = model$dmodel$shapeshock, 
					shapemodel = model$dmodel$shapemodel,
					shape.regressors = model$modeldata$shxdata, exp.rate=model$sbounds[5]))
	setfixed(fspec)<-as.list(fit@model$pars[fit@model$pars[,3]==1,1])
	setbounds(fspec)<-list(shape = fit@model$sbounds[3:4], skew = fit@model$sbounds[1:2])
	fspec@model$modeldata$mexdata = mxf
	fspec@model$modeldata$vexdata = vxf
	fspec@model$modeldata$skxdata = skxf
	fspec@model$modeldata$shxdata = shxf
	# Generate the 1 extra ahead forecast
	if((n.ahead+n.roll)<=fit@model$n.start){
		tmp =  xts(data[1:(N + fcreq)], index[1:(N + fcreq)])
		DailyV = .intraday2daily(DVar[1:(N + fcreq)])
		flt = .mcsacdfilter(spec = fspec, data = tmp, n.old = N, 
				skew0 = fit@fit$tskew[1], shape0 = fit@fit$tshape[1], DailyVar = DailyV)
		if(is.null(DiurnalVar)) DiurnalVar = flt@model$DiurnalVar
		sigmafilter = as.numeric(sigma(flt))
		qfilter = flt@filter$q
		resfilter = as.numeric(residuals(flt))
		zfilter = as.numeric(flt@filter$z)
		eresfilter = resfilter/sqrt(as.numeric(DVar[1:(N + fcreq)])*as.numeric(DiurnalVar[1:(N + fcreq)]))
		tskewfilter 	= flt@filter$tskew
		tshapefilter 	= flt@filter$tshape
		tempskewfilter 	= flt@filter$tempskew
		tempshapefilter = flt@filter$tempshape
		qfor = seriesfor = sigmafor = matrix(NA, ncol = n.roll+1, nrow = n.ahead)
		tskewfor = tshapefor = tempshafor = tempskewfor = matrix(NA, ncol = n.roll+1, nrow = n.ahead)
		# only 1-ahead for mcsACD model
		qfor[1,] = qfilter[(N+1):(N+n.roll+1)]
		seriesfor[1,] = as.numeric(fitted(flt)[(N+1):(N+n.roll+1)])
		sigmafor[1,]  =  as.numeric(sigmafilter[(N+1):(N+n.roll+1)])
		tskewfor[1,]  =  as.numeric(skew(flt)[(N+1):(N+n.roll+1)])
		tshapefor[1,] = as.numeric(shape(flt)[(N+1):(N+n.roll+1)])
		# n.roll x n.ahead (n.ahead=1 generted by model)
		# n.ahead>1 by simulation
		colnames(qfor) = colnames(seriesfor) = colnames(sigmafor) = as.character(index(fitted(flt))[(N+1):(N+n.roll+1)])
		colnames(tskewfor) = colnames(tshapefor) = as.character(index(fitted(flt))[(N+1):(N+n.roll+1)])
		rownames(qfor) = rownames(seriesfor) = rownames(sigmafor) = paste("n.ahead-", 1:n.ahead, sep="")
		rownames(tskewfor) = rownames(tshapefor) = paste("n.ahead-", 1:n.ahead, sep="")
	} else{
		tmp =  xts(c(data[1:(N + fcreq)], 0), c(index[1:(N + fcreq)], ftseq(index[N + fcreq], length.out=1, by=period, 
								interval = fit@model$dtime)))
		DailyV = .intraday2daily(DVar[1:(N + fcreq+1)])
		flt = .mcsacdfilter(spec = fspec, data = tmp, n.old = N, 
				skew0 = fit@fit$tskew[1], shape0 = fit@fit$tshape[1], DailyVar = DailyV)
		if(is.null(DiurnalVar)) DiurnalVar = flt@model$DiurnalVar
		sigmafilter = as.numeric(sigma(flt))
		qfilter = flt@filter$q
		resfilter = as.numeric(residuals(flt))
		zfilter = as.numeric(flt@filter$z)
		eresfilter = resfilter/sqrt(as.numeric(DVar[1:(N + fcreq+1)])*as.numeric(DiurnalVar[1:(N + fcreq+1)]))
		tskewfilter 	= flt@filter$tskew
		tshapefilter 	= flt@filter$tshape
		tempskewfilter 	= flt@filter$tempskew
		tempshapefilter = flt@filter$tempshape
		qfor = seriesfor = sigmafor = matrix(NA, ncol = n.roll+1, nrow = n.ahead)
		tskewfor = tshapefor = tempshafor = tempskewfor = matrix(NA, ncol = n.roll+1, nrow = n.ahead)
		# only 1-ahead for mcsACD model
		qfor[1,] = qfilter[(N+1):(N+n.roll+1)]
		seriesfor[1,] = as.numeric(fitted(flt)[(N+1):(N+n.roll+1)])
		sigmafor[1,]  =  as.numeric(sigmafilter[(N+1):(N+n.roll+1)])
		tskewfor[1,]  =  as.numeric(skew(flt)[(N+1):(N+n.roll+1)])
		tshapefor[1,] = as.numeric(shape(flt)[(N+1):(N+n.roll+1)])
		# n.roll x n.ahead (n.ahead=1 generted by model)
		# n.ahead>1 by simulation
		colnames(qfor) = colnames(seriesfor) = colnames(sigmafor) = as.character(index(fitted(flt))[(N+1):(N+n.roll+1)])
		colnames(tskewfor) = colnames(tshapefor) = as.character(index(fitted(flt))[(N+1):(N+n.roll+1)])
		rownames(qfor) = rownames(seriesfor) = rownames(sigmafor) = paste("n.ahead-", 1:n.ahead, sep="")
		rownames(tskewfor) = rownames(tshapefor) = paste("n.ahead-", 1:n.ahead, sep="")
	}
	mx = model$maxOrder
	fcst = list()
	fcst$n.ahead = n.ahead
	fcst$n.roll = n.roll
	fcst$N = N+ns
	fcst$n.start = ns
	fcst$seriesFor = seriesfor
	fcst$sigmaFor  = sigmafor
	fcst$qFor  = qfor
	fcst$tskewFor  = tskewfor
	fcst$tshapeFor = tshapefor
	model$modeldata$sigma = flt@filter$sigma
	model$modeldata$residuals = flt@filter$residuals
	model$modeldata$tskew  = flt@filter$tskew
	model$modeldata$tshape = flt@filter$tshape
	ans = new("ACDforecast",
			forecast = fcst,
			model = model)
	return(ans)
}

.mcsacdsim = function(fit, n.sim = 1000, n.start = 0, m.sim = 1, presigma = NA, 
		prereturns = NA, preresiduals = NA, preskew = NA, preshape = NA, rseed = NA,  
		mexsimdata = NULL, vexsimdata = NULL, skxsimdata = NULL, shxsimdata = NULL, DailyVar, ...)
{
	if(fit@model$modelinc[4]>0){
		if(n.start<fit@model$modelinc[3]){
			warning("\nugarchsim-->warning: n.start>=MA order for arfima model...automatically setting.")
			n.start = fit@model$modelinc[3]
		}
	}
	# 1. Create Diurnal Var (n.sim)
	T = fit@model$modeldata$T
	m = fit@model$maxOrder
	
	T0 = fit@model$modeldata$index[T-m]
	dtime = fit@model$dtime
	dvalues = fit@model$dvalues
	D = ftseq(T0, length.out = m+n.sim+n.start, by = fit@model$modeldata$period, interval = dtime)
	Dmatch = match(format(D, "%H:%M:%S"), dtime)
	DiurnalVar = xts(dvalues[Dmatch], D)
	# 2. Transform DailyVar into intraday (n.sim by m.sim)
	if(missing(DailyVar)) stop("\nDailyVar cannot be missing for the mcsGARCH model.")
	if(!is(DailyVar, "xts")) stop("\nDailyVar must be an xts object of daily variance simulations for the n.sim period")
	
	DailyVarOld = .intraday2daily(fit@model$DailyVar)
	Unique1 =  unique(format(index(DiurnalVar), "%Y-%m-%d"))
	Unique2 =  unique(format(index(DailyVarOld), "%Y-%m-%d"))
	# find the required dates
	ReqDates = setdiff(Unique1, Unique2)
	Unique3  = format(index(DailyVar), "%Y-%m-%d")
	if(any(is.na(match(ReqDates, Unique3)))){
		stop(paste(c("The required dates for the DailyVar are:", ReqDates), sep="", collapse=" "))
	}
	if(NCOL(DailyVar)!=m.sim){
		if(NCOL(DailyVar)==1){
			warning("\nDailyVar not equal to m.sim...will replicate to use the same for all independent simulations (m.sim)")
			DVar = matrix(NA, ncol = m.sim, nrow = m + n.sim+n.start)
			# need to align old Daily Var with new Daily Var
			DailyVarOld = .intraday2daily(fit@model$DailyVar)
			DailyVar = c(DailyVarOld, DailyVar)
			Dx = .daily2intraday(DiurnalVar, DailyVar)
			Dy = tail(coredata(Dx), m + n.sim+n.start)
			for(i in 1:m.sim) DVar[,i] = Dy
		} else{
			stop("\nNCOL(DailyVar) is greater than 1 and less than m.sim...resubmit something which makes more sense.")
		}
	} else{
		DVar = matrix(NA, ncol = m.sim, nrow = m + n.sim+n.start)
		DailyVarOld = .intraday2daily(fit@model$DailyVar)
		for(i in 1:m.sim){
			DailyVarx = c(DailyVarOld, DailyVar[,i])
			Dx = .daily2intraday(DiurnalVar, DailyVarx)
			DVar[,i] = tail(coredata(Dx), m + n.sim+n.start)
		}
	}
	if(is.na(rseed[1])){
		sseed = as.integer(runif(m.sim,0,as.integer(Sys.time())))
	} else{
		if(length(rseed) != m.sim) stop("\uacdsim-->error: rseed must be of length m.sim!\n")
		sseed = rseed[1:m.sim]
	}
	n = n.sim + n.start
	data = fit@model$modeldata$data
	N = length(as.numeric(data))
	data = data[1:(N - fit@model$n.start)]
	N = length(as.numeric(data))
	resids = fit@fit$residuals
	sigma = fit@fit$sigma
	
	model = fit@model
	modelinc = model$modelinc
	idx = model$pidx
	ipars = model$pars
	sbounds = model$sbounds
	N = 0
	m = model$maxOrder
	if(modelinc[8]>0) {
		mexdata = matrix(model$modeldata$mexdata, ncol = modelinc[8])
		N = dim(mexdata)[1]
	} else { mexdata = NULL }
	if(modelinc[17]>0) {
		vexdata = matrix(model$modeldata$vexdata, ncol = modelinc[17]) 
		N = dim(vexdata)[1]
	} else { vexdata = NULL }
	if(modelinc[25]>0) {
		skexdata = matrix(model$modeldata$skxdata, ncol = modelinc[25]) 
		N = dim(skexdata)[1]
	} else { skexdata = NULL }
	if(modelinc[31]>0) {
		shexdata = matrix(model$modeldata$shxdata, ncol = modelinc[31]) 
		N = dim(shexdata)[1]
	} else { shexdata = NULL }
	
	distribution = model$dmodel$model
	# check if necessary the external regressor forecasts provided first
	xreg = .acdsimregressors(model, mexsimdata, vexsimdata, skxsimdata, shxsimdata, N, n, m.sim, m)	
	mexsim  = xreg$mexsimlist
	vexsim  = xreg$vexsimlist
	skexsim = xreg$skexsimlist
	shexsim = xreg$shexsimlist
	
	if(!is.na(presigma[1])){
		presigma = as.vector(presigma)
		if(length(presigma)<m) stop(paste("\nacdpath-->error: presigma must be of length ", m, sep=""))
	} else{
		presigma = tail(as.numeric(sigma(fit)), m)
	}
	if(!is.na(prereturns[1])){
		prereturns = as.vector(prereturns)
		if(length(prereturns)<m) stop(paste("\nuacdsim-->error: prereturns must be of length ", m, sep=""))
	} else{
		prereturns = tail(model$modeldata$data[1:model$modeldata$T], m)
	}
	if(!is.na(preresiduals[1])){
		preresiduals = as.vector(preresiduals)
		if(length(preresiduals)<m) stop(paste("\nuacdsim-->error: preresiduals must be of length ", m, sep=""))
		preres = matrix(preresiduals, nrow = m, ncol = m.sim)
	} else{
		preres = matrix(tail(fit@fit$residuals, m), nrow = m, ncol = m.sim)
	}
	
	# Random Samples from the Distribution are calculated at every recursion in the
	# c-code as they depend on the actual time-varying skew & shape
	z = matrix(0, ncol = m.sim, nrow = n.sim + n.start)
	z = rbind(matrix(0, nrow = m, ncol = m.sim), z)
	# z = matrix(0, ncol = m.sim, nrow = n.sim+n.start)
	pretskew = pretempskew = rep(0, m)
	
	if(model$modelinc[21]>0)
	{
		if( is.na(preskew[1]) ){
			# The tempskew[1] is the transformed skew parameter of the 
			# non-time varying model from which we initiated the original fit.
			pretempskew = tail(fit@fit$tempskew, m)
			pretskew = tail(fit@fit$tskew, m)
		} else{
			# preskew is provided un-transformed
			pretempskew = logtransform(tail(preskew, m), sbounds[1], sbounds[2], inverse = TRUE)
			pretskew = tail(preskew, m)
		}
	}
	if(model$modelinc[18]>0){
		tskew = rep(ipars["skew", 1], n+m)
	} else{
		tskew = c(pretskew, rep(0, n))
	}
	
	pretshape = pretempshape = rep(0, m)
	if(model$modelinc[27]>0)
	{
		if( is.na(preshape[1]) ){
			# The tempshape[1] is the transformed shape parameter of the 
			# non-time varying model from which we initiated the original fit.
			pretempshape = tail(fit@fit$tempshape, m)
			pretshape = tail(fit@fit$tshape, m)
		} else{
			pretempshape = exptransform(tail(preshape, m), sbounds[3], sbounds[4], sbounds[5], inverse = TRUE)
			pretshape = tail(preshape, m)
		}
	}
	if(model$modelinc[19]>0){
		tshape = rep(ipars["shape", 1], n+m)
	} else{
		tshape = c(pretshape, rep(0, n))
	}
	if(!is.null(list(...)$preq)){
		preq = tail(list(...)$preq^2, m)
	} else{
		preq = tail(as.numeric(fit@fit$q)^2, m)
	}
	
	preeres = apply(preres, 2, function(x) x/tail(sqrt(as.numeric(fit@model$DiurnalVar[1:N])*as.numeric(fit@model$DailyVar[1:N])), m))
	preeres = matrix(preeres, nrow = m, ncol = m.sim)
	
	# input vectors/matrices
	h = c(preq, rep(0, n))
	x = c(prereturns, rep(0, n))
	tmpskew = c(pretempskew, rep(0, n))
	tmpshape = c(pretempshape, rep(0, n))
	constm = matrix(ipars[idx["mu",1]:idx["mu",2], 1], ncol = m.sim, nrow = n + m)
	
	# MATRIX
	if( !is.na(preresiduals) && !is.na(presigma) ){
		zz = preres[1:m]/presigma[1:m]
		for(j in 1:m.sim){
			z[1:m, j] = zz
		}
	} else{
		# ? Do we want the same for all m.sim? If yes, there is no uncertainty for
		# the n.sim = 1 for sigma, and higher moment (equal to their forecast value).
		# If no, then uncertainty is introduced in the n.sim=1 values.
		for(k in 1:m){
			z[k, ] = rugarch:::.makeSample(distribution, skew = tskew[k], shape = tshape[k], 
					lambda = ipars[idx["ghlambda",1],1], n = m.sim, seed = sseed[1]+k)
		}
	}
	#if(is.na(preresiduals[1])){
	#	preres = z[1:m, , drop = FALSE]*matrix(presigma, ncol = m.sim, nrow = m)
	#}
	res = rbind( preres, matrix(0, ncol = m.sim, nrow = n) )
	eres = rbind(preeres, matrix(0, ncol = m.sim, nrow = n) )
	# eres should be based on res!!!!

	# outpus matrices
	qSim = sigmaSim =  matrix(0, ncol = m.sim, nrow = n.sim)
	seriesSim = matrix(0, ncol = m.sim, nrow = n.sim)
	residSim =  matrix(0, ncol = m.sim, nrow = n.sim)
	skewSim 		= matrix(0, ncol = m.sim, nrow = n.sim)
	tempskewSim 	= matrix(0, ncol = m.sim, nrow = n.sim)
	shapeSim 		= matrix(0, ncol = m.sim, nrow = n.sim)
	tempshapeSim 	= matrix(0, ncol = m.sim, nrow = n.sim)
	zSim 			= matrix(0, ncol = m.sim, nrow = n.sim)
	
	for(i in 1:m.sim){
		set.seed(sseed[i])
		tmp = try(.C("mcsacdsimC", 
						model = as.integer(modelinc), 
						pars = as.double(ipars[,1]), 
						idx = as.integer(idx[,1]-1), 
						h = as.double(h), 
						s = as.double(sqrt(DiurnalVar)), 
						v = as.double(sqrt(DVar[,i])),
						z = as.double(z[,i]), 
						res = as.double(res[,i]),
						eres = as.double(eres[,i]),
						e = as.double(eres[,i]*eres[,i]), 
						tempskew = as.double(tmpskew), 
						tempshape = as.double(tmpshape), 
						tskew = as.double(tskew), 
						tshape = as.double(tshape), 
						sbounds = as.double(sbounds), 
						vexdata = as.double(vexsim[[i]]), 
						skxreg = as.double(skexsim[[i]]), 
						shxreg = as.double(shexsim[[i]]),
						T = as.integer(n+m), 
						m = as.integer(m),
						PACKAGE = "racd"), silent = TRUE)
		qSim[,i] = tmp$h[(n.start + m + 1):(n+m)]^(1/2)
		sigmaSim[,i] = qSim[,i]*sqrt(DVar[-c(1:m),i]*as.numeric(DiurnalVar[-c(1:m)]))
		res = c(preres[,i], tmp$eres[-c(1:m)]*sqrt(DVar[-c(1:m),i]*as.numeric(DiurnalVar[-c(1:m)])))
		residSim[,i] = res[(n.start + m + 1):(n+m)]
		skewSim[,i] 	 = tail(tmp$tskew, n.sim)
		shapeSim[,i] 	 = tail(tmp$tshape, n.sim)
		zSim[,i] 		 = tail(tmp$z, n.sim)
		if(modelinc[8]>0){
			mxreg = matrix( ipars[idx["mxreg",1]:idx["mxreg",2], 1], ncol = modelinc[8] )
			constm[,i] = constm[,i] + mxreg %*%t( matrix( mexsim[[i]], ncol = modelinc[8] ) )
		}
		if(modelinc[4]>0){
			fres = c(res[(m+1):(n+m)], if(modelinc[3]>0) rep(0, modelinc[3]) else NULL)
			ans2 = rugarch:::.arfimaxsim(modelinc[1:5], ipars, idx, constm[1:n, i], fres, T = n)
			seriesSim[,i] = head(ans2$series, n.sim)
		} else{
			ans2 = rugarch:::.armaxsim(modelinc[1:5], ipars = ipars, idx = idx, constm = constm[,i],  
					x = x, res = res, T = n + m, m)
			seriesSim[,i] = ans2$x[(n.start + m + 1):(n+m)]
		}
	}
	sim = list(qSim = qSim, sigmaSim = sigmaSim, seriesSim = seriesSim, residSim = residSim, skewSim = skewSim,
			shapeSim = shapeSim, zSim = zSim)
	sim$n.sim  = n.sim
	sim$m.sim  = m.sim
	sim$DiurnalVar = DiurnalVar[-c(1:m)]
	sim$DailyVar = DVar[-c(1:m),]
	model$modeldata$sigma = sigma
	sol = new("ACDsim",
			simulation = sim,
			model = model,
			seed = as.integer(sseed))
	return(sol)
}

################################################################################
# special purpose rolling method for the mcsGARCH model
################################################################################
.acdroll.mcs = function(spec, data, n.ahead = 1, forecast.length = 500, n.start = NULL, 
		refit.every = 25, refit.window = c("recursive", "moving"), 
		window.size = NULL, solver = "ucminf", fit.control = list(), 
		solver.control = list(), calculate.VaR = TRUE, VaR.alpha = c(0.01, 
				0.05), cluster = NULL, keep.coef = TRUE, fixARMA = TRUE, fixGARCH = TRUE, 
		fixUBShape = TRUE, UBShapeAdd = 0, fixGHlambda = TRUE, compareGARCH = c("LL", "none"),
		DailyVar, ...)
{
	tic = Sys.time()
	compareGARCH = compareGARCH[1]
	if(is.null(solver.control$trace)) trace = 0 else trace = solver.control$trace
	if(is.null(fit.control$stationarity)) fit.control$stationarity = TRUE
	if(is.null(fit.control$fixed.se)) fit.control$fixed.se = FALSE
	if(is.null(fit.control$scale)) fit.control$scale = FALSE
	if(is.null(fit.control$n.sim)) fit.control$n.sim = 2000
	mm = match(names(fit.control), c("stationarity", "fixed.se", "scale", "n.sim"))
	if(any(is.na(mm))){
		idx = which(is.na(mm))
		enx = NULL
		for(i in 1:length(idx)) enx = c(enx, names(fit.control)[idx[i]])
		warning(paste(c("unidentified option(s) in fit.control:\n", enx), sep="", collapse=" "), call. = FALSE, domain = NULL)
	}
	###########################################################################
	# Check DailyVar
	if(is.null(DailyVar)){
		stop("\nacdroll-->error: you must supply the daily forecast variance (DailyVar) for the msGARCH model\n")
	} else{
		if(!is(DailyVar, "xts")) stop("\nacdroll-->error: DailyVar must be an xts object\n")
		DailyVarIndex = format(index(DailyVar), format="%Y-%m-%d")
	}
	# we are not going to extract the data just yet
	UIndex = unique(format(index(data), format="%Y-%m-%d"))
	DIndex = format(index(data), format="%Y-%m-%d")
	RIndex = index(data)
	M = length(UIndex)
	matchD = all.equal(UIndex, DailyVarIndex)
	if(!matchD) stop("\nacdroll-->error: DailyVar dates do not match the data dates (unique days).\n")
	
	refit.window = refit.window[1]
	datanames = names(data)
	xdata = rugarch:::.extractdata(data)
	data = xdata$data
	index = xdata$index
	period = xdata$period
	T = NROW(data)
	modelinc = spec@model$modelinc
	if( modelinc[8]>0 ){
		mexdata = spec@model$modeldata$mexdata
		mex = TRUE
	} else{
		mex = FALSE
		mexdata = NULL
	}
	if( modelinc[17]>0 ){
		vexdata = spec@model$modeldata$vexdata
		vex = TRUE
	} else{
		vex = FALSE
		vexdata = NULL
	}
	if( modelinc[25]>0 ){
		skdata = spec@model$modeldata$skxdata
		skex = TRUE
	} else{
		skdata = NULL
		skex = FALSE
	}
	if( modelinc[31]>0 ){
		shdata = spec@model$modeldata$shxdata
		shex = TRUE
	} else{
		shdata = NULL
		shex = FALSE
	}
	
	
	if(n.ahead>1) stop("\nacdroll:--> n.ahead>1 not supported...try again.")
	if(is.null(n.start)){
		if(is.null(forecast.length)) stop("\nacdroll:--> forecast.length amd n.start are both NULL....try again.")
		n.start = T - forecast.length
	} else{
		forecast.length = T - n.start
	}
	if(T<=n.start) stop("\nacdroll:--> start cannot be greater than length of data")
	# the ending points of the estimation window
	s = seq(n.start+refit.every, T, by = refit.every)
	m = length(s)
	# the rolling forecast length
	out.sample = rep(refit.every, m)
	# adjustment to include all the datapoints from the end
	if(s[m]<T){
		s = c(s,T)
		m = length(s)
		out.sample = c(out.sample, s[m]-s[m-1])
	}
	if(refit.window == "recursive"){
		rollind = lapply(1:m, FUN = function(i) 1:s[i])
	} else{
		if(!is.null(window.size)){
			if(window.size<100) stop("\nacdroll:--> window size must be greater than 100.")
			rollind = lapply(1:m, FUN = function(i) max(1, (s[i]-window.size-out.sample[i])):s[i])
		} else{
			rollind = lapply(1:m, FUN = function(i) (1+(i-1)*refit.every):s[i])
		}
	}
	DailyV = lapply(rollind, FUN = function(x){
				dindex = unique(format(index[x], "%Y-%m-%d"))
				DailyVar[dindex]
			})
	# distribution
	distribution = spec@model$dmodel$model
	if(any(distribution==c("snorm","sstd","sged","jsu","nig","ghyp","ghst"))) skewed = TRUE else skewed = FALSE
	if(any(distribution==c("std","sstd","ged","sged","jsu","nig","ghyp","ghst"))) shaped = TRUE else shaped = FALSE
	if(any(distribution==c("ghyp"))) ghyp = TRUE else ghyp = FALSE
	
	gspec = .spec2GARCH(spec)
	if( !is.null(cluster) ){
		parallel::clusterEvalQ(cl = cluster, library(racd))
		parallel::clusterExport(cluster, c("data", "index", "s","refit.every", 
						"keep.coef", "shaped", "skewed", "ghyp", "gspec", "fixARMA",
						"fixGARCH", "fixUBShape", "UBShapeAdd", "fixGHlambda","compareGARCH",
						"rollind", "spec", "out.sample", "mex", "vex", "skex", "shex",
						"solver", "solver.control", "fit.control", "DailyV"), envir = environment())
		if(mex) parallel::clusterExport(cluster, c("mexdata"), envir = environment())
		if(vex)  parallel::clusterExport(cluster, c("vexdata"), envir = environment())
		if(skex)  parallel::clusterExport(cluster, c("skdata"), envir = environment())
		if(shex)  parallel::clusterExport(cluster, c("shdata"), envir = environment())
		tmp = parallel::parLapplyLB(cl = cluster, 1:m, fun = function(i){
					zspec = spec
					xspec = gspec
					if(mex){
						zspec@model$modeldata$mexdata = mexdata[rollind[[i]],,drop=FALSE]
						xspec@model$modeldata$mexdata = mexdata[rollind[[i]],,drop=FALSE]
					}
					if(vex){
						zspec@model$modeldata$vexdata = vexdata[rollind[[i]],,drop=FALSE]
						xspec@model$modeldata$vexdata = vexdata[rollind[[i]],,drop=FALSE]
					}
					if(skex) zspec@model$modeldata$skxdata = skdata[rollind[[i]],,drop=FALSE]
					if(shex) zspec@model$modeldata$shxdata = shdata[rollind[[i]],,drop=FALSE]
					gfit = ugarchfit(xspec, xts::xts(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i], 
							solver = "hybrid", DailyVar = DailyV[[i]])
					if(convergence(gfit)==0){
						if(fixGHlambda && zspec@model$dmodel$model=="ghyp"){
							setfixed(zspec)<-list(ghlambda=unname(coef(gfit)["ghlambda"]))
						}
						if(fixARMA && fixGARCH){
							setfixed(zspec)<-as.list(coef(gfit)[1:sum(gspec@model$modelinc[1:15])])
						} else if(fixARMA && !fixGARCH){
							setfixed(zspec)<-as.list(coef(gfit)[1:sum(gspec@model$modelinc[1:6])])
						} else if(!fixARMA && fixGARCH){
							setfixed(zspec)<-as.list(coef(gfit)[(sum(gspec@model$modelinc[1:6])+1):sum(gspec@model$modelinc[7:15])])
						} else{
							setstart(zspec)<-as.list(coef(gfit)[1:sum(gspec@model$modelinc[1:15])])
						}
						if(fixUBShape){
							sh = coef(gfit)["shape"]
							setbounds(zspec)<-list(shape = c(spec@model$sbounds[3], sh+UBShapeAdd))
						}
						if(xspec@model$modelinc[16]>0) skew0 = coef(gfit)["skew"] else skew0 = NULL
						if(xspec@model$modelinc[17]>0) shape0 = coef(gfit)["shape"] else shape0 = NULL
						glik = likelihood(gfit)
					} else{
						shape0 = NULL
						skew0 = NULL
						glik = NA
					}
					fit = try(acdfit(zspec, xts::xts(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i], 
									solver = solver, solver.control = solver.control, 
									fit.control = fit.control, shape0 = shape0, skew0 = skew0,
									DailyVar = DailyV[[i]]), silent=TRUE)
					if(inherits(fit, 'try-error') || convergence(fit)!=0 || is.null(fit@fit$cvar)){
						ans = list(y = NA, cf = NA, converge = FALSE, lik = c(NA, glik))
					} else{
						# compare GARCH likelihood with ACD model and reject if lik less than
						clik = likelihood(fit)[1]
						if(!is.na(glik) && clik[1]<glik[1] && compareGARCH=="LL"){
							ans = list(y = NA, cf = NA, q = NA, DiurnalVar = NA, DailyVar = NA, 
									converge = FALSE, lik = c(clik, glik))
						} else{
							if(mex) fmex = tail(mexdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fmex = NULL
							if(vex) fvex = tail(vexdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fvex = NULL
							if(skex) fskx = tail(skdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fskx = NULL
							if(shex) fshx = tail(shdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fshx = NULL
							# since the forecast is in the out.sample we don't need to supply DailyVar
							f = acdforecast(fit, n.ahead = 1, n.roll = out.sample[i]-1, 
									external.forecasts = list(mregfor = fmex, vregfor = fvex, skregfor = fskx,
											shregfor = fshx))
							sig = as.numeric(sigma(f))
							ret = as.numeric(fitted(f))
							if(shaped) shp = as.numeric(shape(f)) else shp = rep(0, out.sample[i])
							if(skewed) skw = as.numeric(skew(f))  else skw = rep(0, out.sample[i])
							if(ghyp) shpgig = rep(coef(fit)["ghlambda"], out.sample[i]) else shpgig = rep(0, out.sample[i])
							rlz = tail(data[rollind[[i]]], out.sample[i])
							# use xts for indexing the forecasts
							y = as.data.frame(cbind(ret, sig, skw, shp, shpgig, rlz))
							rownames(y) = tail(as.character(index[rollind[[i]]]), out.sample[i])				
							colnames(y) = c("Mu", "Sigma", "Skew", "Shape", "Shape(GIG)", "Realized")
							if(keep.coef) cf = fit@fit$robust.matcoef else cf = NA
							ans = list(y = y, cf = cf, q = f@forecast$qFor, DiurnalVar = f@model$DiurnalVar, 
									DailyVar = f@model$DailyVar, converge = TRUE, lik = c(likelihood(fit)[1], glik))
						}
					}
					return(ans)})
	} else{
		tmp = vector(mode = "list", length = m)
		for(i in 1:m){
			zspec = spec
			xspec = gspec
			if(mex){
				zspec@model$modeldata$mexdata = mexdata[rollind[[i]],,drop=FALSE]
				xspec@model$modeldata$mexdata = mexdata[rollind[[i]],,drop=FALSE]
			}
			if(vex){
				zspec@model$modeldata$vexdata = vexdata[rollind[[i]],,drop=FALSE]
				xspec@model$modeldata$vexdata = vexdata[rollind[[i]],,drop=FALSE]
			}
			if(skex) zspec@model$modeldata$skxdata = skdata[rollind[[i]],,drop=FALSE]
			if(shex) zspec@model$modeldata$shxdata = shdata[rollind[[i]],,drop=FALSE]
			gfit = ugarchfit(xspec, xts(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i], 
					solver = "hybrid", DailyVar = DailyV[[i]])
			if(convergence(gfit)==0){
				if(fixGHlambda && zspec@model$dmodel$model=="ghyp"){
					setfixed(zspec)<-list(ghlambda=unname(coef(gfit)["ghlambda"]))
				}
				if(fixARMA && fixGARCH){
					setfixed(zspec)<-as.list(coef(gfit)[1:sum(gspec@model$modelinc[1:15])])
				} else if(fixARMA && !fixGARCH){
					setfixed(zspec)<-as.list(coef(gfit)[1:sum(gspec@model$modelinc[1:6])])
				} else if(!fixARMA && fixGARCH){
					setfixed(zspec)<-as.list(coef(gfit)[(sum(gspec@model$modelinc[1:6])+1):sum(gspec@model$modelinc[7:15])])
				} else{
					setstart(zspec)<-as.list(coef(gfit)[1:sum(gspec@model$modelinc[1:15])])
				}
				if(fixUBShape){
					sh = coef(gfit)["shape"]
					setbounds(zspec)<-list(shape = c(spec@model$sbounds[3], sh+UBShapeAdd))
				}
				if(xspec@model$modelinc[16]>0) skew0 = coef(gfit)["skew"] else skew0 = NULL
				if(xspec@model$modelinc[17]>0) shape0 = coef(gfit)["shape"] else shape0 = NULL
				glik = likelihood(gfit)
			} else{
				shape0 = NULL
				skew0 = NULL
				glik = NA
			}
			fit = try(acdfit(zspec, xts::xts(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i], 
							solver = solver, solver.control = solver.control, 
							fit.control = fit.control, shape0 = shape0, skew0 = skew0, 
							DailyVar = DailyV[[i]]), silent=TRUE)
			if(inherits(fit, 'try-error') || convergence(fit)!=0 || is.null(fit@fit$cvar)){
				tmp[[i]] = list(y = NA, cf = NA, converge = FALSE, lik = c(NA, glik))
			} else{
				clik = likelihood(fit)[1]
				if(!is.na(glik) && clik[1]<glik[1] && compareGARCH=="LL"){
					ans = list(y = NA, cf = NA, q = NA, DiurnalVar = NA, DailyVar = NA, 
							converge = FALSE, lik = c(clik, glik))
				} else{
					if(mex) fmex = tail(mexdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fmex = NULL
					if(vex) fvex = tail(vexdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fvex = NULL
					if(skex) fskx = tail(skdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fskx = NULL
					if(shex) fshx = tail(shdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fshx = NULL
					
					f = acdforecast(fit, n.ahead = 1, n.roll = out.sample[i]-1, 
							external.forecasts = list(mregfor = fmex, vregfor = fvex, skregfor = fskx,
									shregfor = fshx))
					sig = as.numeric(sigma(f))
					ret = as.numeric(fitted(f))
					if(shaped) shp = as.numeric(shape(f)) else shp = rep(0, out.sample[i])
					if(skewed) skw = as.numeric(skew(f)) else skw = rep(0, out.sample[i])
					if(ghyp) shpgig = rep(coef(fit)["ghlambda"], out.sample[i]) else shpgig = rep(0, out.sample[i])
					rlz = tail(data[rollind[[i]]], out.sample[i])
					# use xts for indexing the forecasts
					y = as.data.frame(cbind(ret, sig, skw, shp, shpgig, rlz))
					rownames(y) = tail(as.character(index[rollind[[i]]]), out.sample[i])
					colnames(y) = c("Mu", "Sigma", "Skew", "Shape", "Shape(GIG)", "Realized")
					if(keep.coef) cf = fit@fit$robust.matcoef else cf = NA
					tmp[[i]] = list(y = y, cf = cf, q = f@forecast$qFor, DiurnalVar = f@model$DiurnalVar, 
							DailyVar = f@model$DailyVar, converge = TRUE, lik = c(likelihood(fit)[1], glik))
				}
			}
		}
	}
	conv = sapply(tmp, FUN = function(x) x$converge)
	if(any(!conv)){
		warning("\nnon-converged estimation windows present...resubsmit object with different solver parameters...")
		noncidx = which(!conv)
		model = list()
		model$fixARMA = fixARMA
		model$fixGARCH = fixGARCH
		model$fixUBShape = fixUBShape
		model$UBShapeAdd = UBShapeAdd
		model$fixGHlambda = fixGHlambda
		model$compareGARCH = compareGARCH
		model$spec = spec
		model$data = data
		model$index = index
		model$period = period
		model$datanames = datanames
		model$n.ahead = n.ahead
		model$forecast.length = forecast.length 
		model$n.start = n.start
		model$n.refits = m
		model$refit.every = refit.every
		model$refit.window = refit.window
		model$window.size = window.size
		model$calculate.VaR = calculate.VaR
		model$VaR.alpha = VaR.alpha
		model$keep.coef = keep.coef
		model$noncidx = noncidx
		model$rollind = rollind
		model$out.sample = out.sample
		model$DailyVar = DailyVar
		forecast = tmp
		toc = Sys.time()-tic
		model$elapsed = toc
		ans = new("ACDroll",
				model = model,
				forecast = forecast)
		return(ans)
	} else{
		noncidx = NULL
		forc = tmp[[1]]$y
		if(m>1){
			for(i in 2:m){
				forc = rbind(forc, tmp[[i]]$y)
			}
		}
		if(keep.coef){
			cf = vector(mode = "list", length = m)
			for(i in 1:m){
				cf[[i]]$index = index[tail(rollind[[i]],1) - out.sample[i]]
				cf[[i]]$coef = tmp[[i]]$cf
			}
		} else{
			cf = NULL
		}
		LL = vector(mode = "list", length = m)
		for(i in 1:m){
			LL[[i]]$index = index[tail(rollind[[i]],1) - out.sample[i]]
			LL[[i]]$log.likelihood = tmp[[i]]$lik
		}
		
		if(calculate.VaR){
			if(is.null(VaR.alpha)) VaR.alpha = c(0.01, 0.05)
			VaR.matrix = matrix(NA, ncol = length(VaR.alpha)+1, nrow = NROW(forc))
			for(i in 1:length(VaR.alpha)){
				VaR.matrix[,i] = qdist(p = VaR.alpha[i], mu = forc[,1], sigma = forc[,2], 
						skew = forc[,3], shape = forc[,4], lambda = forc[,5], 
						distribution = distribution)
			}
			VaR.matrix[,length(VaR.alpha)+1] = forc[,6]
			colnames(VaR.matrix) = c(paste("alpha(", round(VaR.alpha,2)*100, "%)",sep=""), "realized")
			VaR.matrix = as.data.frame(VaR.matrix)
			rownames(VaR.matrix) = rownames(forc)
		} else{
			VaR.matrix = NULL
		}
		model = list()
		model$spec = spec
		model$data = data
		model$index = index
		model$period = period
		model$n.ahead = n.ahead
		model$forecast.length = forecast.length 
		model$n.start = n.start
		model$refit.every = refit.every
		model$n.refits = m
		model$refit.window = refit.window
		model$window.size = window.size
		model$calculate.VaR = calculate.VaR
		model$VaR.alpha = VaR.alpha
		model$keep.coef = keep.coef
		model$noncidx = noncidx
		model$coef = cf
		model$LL = LL
		model$rollind = rollind
		model$out.sample = out.sample
		forecast = list(VaR = VaR.matrix, density = forc)
	}
	toc = Sys.time()-tic
	model$elapsed = toc
	model$fixARMA = fixARMA
	model$fixGARCH = fixGARCH
	model$fixUBShape = fixUBShape
	model$UBShapeAdd = UBShapeAdd
	model$fixGHlambda = fixGHlambda
	model$compareGARCH = compareGARCH
	ans = new("ACDroll",
			model = model,
			forecast = forecast)
	return( ans )
}

.acdresumeroll.mcs = function(object, spec = NULL, solver = "ucminf", fit.control = list(), 
		solver.control = list(), cluster = NULL, fixARMA = NULL, fixGARCH = NULL, 
		fixUBShape = NULL, UBShapeAdd = NULL, fixGHlambda = NULL, compareGARCH = NULL)
{
	if(!is.null(object@model$noncidx)){
		noncidx = object@model$noncidx
		tic = Sys.time()
		if(is.null(solver.control$trace)) trace = 0 else trace = solver.control$trace
		if(is.null(fit.control$stationarity)) fit.control$stationarity = TRUE
		if(is.null(fit.control$fixed.se)) fit.control$fixed.se = FALSE
		if(is.null(fit.control$scale)) fit.control$scale = FALSE
		if(is.null(fit.control$n.sim)) fit.control$n.sim = 2000
		mm = match(names(fit.control), c("stationarity", "fixed.se", "scale", "n.sim"))
		if(any(is.na(mm))){
			idx = which(is.na(mm))
			enx = NULL
			for(i in 1:length(idx)) enx = c(enx, names(fit.control)[idx[i]])
			warning(paste(c("unidentified option(s) in fit.control:\n", enx), sep="", collapse=" "), call. = FALSE, domain = NULL)
		}
		model = object@model
		if(is.null(fixARMA))  fixARMA = model$fixARMA
		if(is.null(fixGARCH)) fixGARCH = model$fixGARCH
		if(is.null(fixUBShape)) fixUBShape = model$fixUBShape
		if(is.null(UBShapeAdd)) UBShapeAdd = model$UBShapeAdd
		if(is.null(fixGHlambda)) fixGHlambda = model$fixGHlambda
		if(is.null(compareGARCH)) compareGARCH = model$compareGARCH
		
		keep.coef = model$keep.coef
		if(is.null(spec)) spec = model$spec
		gspec = .spec2GARCH(spec)
		datanames = model$datanames
		data = model$data
		index = model$index
		period = model$period
		T = NROW(data)
		modelinc = spec@model$modelinc
		calculate.VaR = model$calculate.VaR
		VaR.alpha = model$VaR.alpha
		if( modelinc[8]>0 ){
			mexdata = spec@model$modeldata$mexdata
			mex = TRUE
		} else{
			mex = FALSE
			mexdata = NULL
		}
		if( modelinc[17]>0 ){
			vexdata = spec@model$modeldata$vexdata
			vex = TRUE
		} else{
			vex = FALSE
			vexdata = NULL
		}
		if( modelinc[25]>0 ){
			skdata = spec@model$modeldata$skxdata
			skex = TRUE
		} else{
			skdata = NULL
			skex = FALSE
		}
		if( modelinc[31]>0 ){
			shdata = spec@model$modeldata$shxdata
			shex = TRUE
		} else{
			shdata = NULL
			shex = FALSE
		}
		n.ahead = model$n.ahead
		n.start = model$n.start
		forecast.length = model$forecast.length
		refit.every = model$refit.every
		refit.window = model$refit.window
		window.size = model$window.size
		if(n.ahead>1) stop("\nacdroll:--> n.ahead>1 not supported...try again.")
		if(is.null(n.start)){
			if(is.null(forecast.length)) stop("\nacdroll:--> forecast.length amd n.start are both NULL....try again.")
			n.start = T - forecast.length
		}
		if(T<=n.start) stop("\nacdroll:--> start cannot be greater than length of data")
		# the ending points of the estimation window
		s = seq(n.start+refit.every, T, by = refit.every)
		m = length(s)
		# the rolling forecast length
		out.sample = rep(refit.every, m)
		# adjustment to include all the datapoints from the end
		if(s[m]<T){
			s = c(s,T)
			m = length(s)
			out.sample = c(out.sample, s[m]-s[m-1])
		}
		if(refit.window == "recursive"){
			rollind = lapply(1:m, FUN = function(i) 1:s[i])
		} else{
			if(!is.null(window.size)){
				if(window.size<100) stop("\nacdroll:--> window size must be greater than 100.")
				rollind = lapply(1:m, FUN = function(i) max(1, (s[i]-window.size-out.sample[i])):s[i])
			} else{
				rollind = lapply(1:m, FUN = function(i) (1+(i-1)*refit.every):s[i])
			}
		}
		DailyVar = object@model$DailyVar
		DailyV = lapply(rollind, FUN = function(x){
					dindex = unique(format(index[x], "%Y-%m-%d"))
					DailyVar[dindex]
				})
		# distribution
		distribution = model$spec@model$dmodel$model
		if(any(distribution==c("snorm","sstd","sged","jsu","nig","ghyp","ghst"))) skewed = TRUE else skewed = FALSE
		if(any(distribution==c("std","sstd","ged","sged","jsu","nig","ghyp","ghst"))) shaped = TRUE else shaped = FALSE
		if(any(distribution==c("ghyp"))) ghyp = TRUE else ghyp = FALSE
		if( !is.null(cluster) ){
			parallel::clusterEvalQ(cl = cluster, library(racd))
			parallel::clusterExport(cluster, c("data", "index","s","refit.every",
							"keep.coef", "shaped", "skewed", "ghyp", "gspec", "fixARMA",
							"fixGARCH", "fixUBShape", "UBShapeAdd", "fixGHlambda","compareGARCH",
							"rollind", "spec", "out.sample", "mex", "vex", "skex", "shex",
							"noncidx", "solver", "solver.control", "fit.control", "DailyV"),
					envir = environment())
			if(mex) parallel::clusterExport(cluster,  c("mexdata"), envir = environment())
			if(vex)  parallel::clusterExport(cluster, c("vexdata"), envir = environment())
			if(skex)  parallel::clusterExport(cluster, c("skdata"), envir = environment())
			if(shex)  parallel::clusterExport(cluster, c("shdata"), envir = environment())
			tmp = parallel::parLapplyLB(cl = cluster, as.list(noncidx), fun = function(i){
						zspec = spec
						xspec = gspec
						if(mex){
							zspec@model$modeldata$mexdata = mexdata[rollind[[i]],,drop=FALSE]
							xspec@model$modeldata$mexdata = mexdata[rollind[[i]],,drop=FALSE]
						}
						if(vex){
							zspec@model$modeldata$vexdata = vexdata[rollind[[i]],,drop=FALSE]
							xspec@model$modeldata$vexdata = vexdata[rollind[[i]],,drop=FALSE]
						}
						if(skex) zspec@model$modeldata$skxdata = skdata[rollind[[i]],,drop=FALSE]
						if(shex) zspec@model$modeldata$shxdata = shdata[rollind[[i]],,drop=FALSE]
						gfit = ugarchfit(xspec, xts::xts(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i], 
								solver = "hybrid", DailyVar = DailyV[[i]])
						if(convergence(gfit)==0){
							if(fixGHlambda && zspec@model$dmodel$model=="ghyp"){
								setfixed(zspec)<-list(ghlambda=unname(coef(gfit)["ghlambda"]))
							}
							if(fixARMA && fixGARCH){
								setfixed(zspec)<-as.list(coef(gfit)[1:sum(gspec@model$modelinc[1:15])])
							} else if(fixARMA && !fixGARCH){
								setfixed(zspec)<-as.list(coef(gfit)[1:sum(gspec@model$modelinc[1:6])])
							} else if(!fixARMA && fixGARCH){
								setfixed(zspec)<-as.list(coef(gfit)[(sum(gspec@model$modelinc[1:6])+1):sum(gspec@model$modelinc[7:15])])
							} else{
								setstart(zspec)<-as.list(coef(gfit)[1:sum(gspec@model$modelinc[1:15])])
							}
							if(fixUBShape){
								sh = coef(gfit)["shape"]
								setbounds(zspec)<-list(shape = c(spec@model$sbounds[3], sh+UBShapeAdd))
							}
							if(xspec@model$modelinc[16]>0) skew0 = coef(gfit)["skew"] else skew0 = NULL
							if(xspec@model$modelinc[17]>0) shape0 = coef(gfit)["shape"] else shape0 = NULL
							glik = likelihood(gfit)
						} else{
							shape0 = NULL
							skew0 = NULL
							glik = NA
						}
						fit = try(acdfit(spec, xts::xts(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i], 
										solver = solver, solver.control = solver.control, 
										fit.control = fit.control, shape0 = shape0, 
										skew0 = skew0, DailyVar = DailyV[[i]]), silent=TRUE)
						if(inherits(fit, 'try-error') || convergence(fit)!=0 || is.null(fit@fit$cvar)){
							ans = list(y = NA, cf = NA, q = NA, DiurnalVar = NA, DailyVar = NA, converge = FALSE, lik = c(NA, glik))
						} else{
							clik = likelihood(fit)[1]
							if(!is.na(glik) && clik[1]<glik[1] && compareGARCH=="LL"){
								ans = list(y = NA, cf = NA, converge = FALSE, lik = c(clik, glik))
							} else{
								if(mex) fmex = tail(mexdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fmex = NULL
								if(vex) fvex = tail(vexdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fvex = NULL
								if(skex) fskx = tail(skdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fskx = NULL
								if(shex) fshx = tail(shdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fshx = NULL
								f = acdforecast(fit, n.ahead = 1, n.roll = out.sample[i]-1, 
										external.forecasts = list(mregfor = fmex, vregfor = fvex, skregfor = fskx,
												shregfor = fshx))
								sig = as.numeric(sigma(f))
								ret = as.numeric(fitted(f))
								if(shaped) shp = as.numeric(shape(f)) else shp = rep(0, out.sample[i])
								if(skewed) skw = as.numeric(skew(f))  else skw = rep(0, out.sample[i])
								if(ghyp) shpgig = rep(coef(fit)["ghlambda"], out.sample[i]) else shpgig = rep(0, out.sample[i])
								rlz = tail(data[rollind[[i]]], out.sample[i])
								# use xts for indexing the forecasts
								y = as.data.frame(cbind(ret, sig, skw, shp, shpgig, rlz))
								rownames(y) = tail(as.character(index[rollind[[i]]]), out.sample[i])
								colnames(y) = c("Mu", "Sigma", "Skew", "Shape", "Shape(GIG)", "Realized")
								if(keep.coef) cf = fit@fit$robust.matcoef else cf = NA
								ans = list(y = y, cf = cf, q = f@forecast$qFor, DiurnalVar = f@model$DiurnalVar, 
										DailyVar = f@model$DailyVar, converge = TRUE, lik = c(likelihood(fit)[1], glik))
							}
						}
						return(ans)})
		} else{
			tmp = lapply(as.list(noncidx), FUN = function(i){
						zspec = spec
						xspec = gspec
						if(mex){
							zspec@model$modeldata$mexdata = mexdata[rollind[[i]],,drop=FALSE]
							xspec@model$modeldata$mexdata = mexdata[rollind[[i]],,drop=FALSE]
						}
						if(vex){
							zspec@model$modeldata$vexdata = vexdata[rollind[[i]],,drop=FALSE]
							xspec@model$modeldata$vexdata = vexdata[rollind[[i]],,drop=FALSE]
						}
						if(skex) zspec@model$modeldata$skxdata = skdata[rollind[[i]],,drop=FALSE]
						if(shex) zspec@model$modeldata$shxdata = shdata[rollind[[i]],,drop=FALSE]
						gfit = ugarchfit(xspec, xts::xts(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i], 
								solver = "hybrid", DailyVar = DailyV[[i]])
						if(convergence(gfit)==0){
							if(fixGHlambda && zspec@model$dmodel$model=="ghyp"){
								setfixed(zspec)<-list(ghlambda=unname(coef(gfit)["ghlambda"]))
							}
							if(fixARMA && fixGARCH){
								setfixed(zspec)<-as.list(coef(gfit)[1:sum(gspec@model$modelinc[1:15])])
							} else if(fixARMA && !fixGARCH){
								setfixed(zspec)<-as.list(coef(gfit)[1:sum(gspec@model$modelinc[1:6])])
							} else if(!fixARMA && fixGARCH){
								setfixed(zspec)<-as.list(coef(gfit)[(sum(gspec@model$modelinc[1:6])+1):sum(gspec@model$modelinc[7:15])])
							} else{
								setstart(zspec)<-as.list(coef(gfit)[1:sum(gspec@model$modelinc[1:15])])
							}
							if(fixUBShape){
								sh = coef(gfit)["shape"]
								setbounds(zspec)<-list(shape = c(spec@model$sbounds[3], sh+UBShapeAdd))
							}
							if(xspec@model$modelinc[16]>0) skew0 = coef(gfit)["skew"] else skew0 = NULL
							if(xspec@model$modelinc[17]>0) shape0 = coef(gfit)["shape"] else shape0 = NULL
							glik = likelihood(gfit)
						} else{
							shape0 = NULL
							skew0 = NULL
							glik = NA
						}
						fit = try(acdfit(spec, xts::xts(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i], 
										solver = solver, solver.control = solver.control, 
										fit.control = fit.control, shape0 = shape0, 
										skew0 = skew0, DailyVar = DailyV[[i]]), silent=TRUE)
						if(inherits(fit, 'try-error') || convergence(fit)!=0 || is.null(fit@fit$cvar)){
							ans = list(y = NA, cf = NA, converge = FALSE, lik = c(NA, glik))
						} else{
							clik = likelihood(fit)[1]
							if(!is.na(glik) && clik[1]<glik[1] && compareGARCH=="LL"){
								ans = list(y = NA, cf = NA, q = NA, DiurnalVar = NA, DailyVar = NA, converge = FALSE, lik = c(clik, glik))
							} else{
								if(mex) fmex = tail(mexdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fmex = NULL
								if(vex) fvex = tail(vexdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fvex = NULL
								if(skex) fskx = tail(skdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fskx = NULL
								if(shex) fshx = tail(shdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fshx = NULL
								
								f = acdforecast(fit, n.ahead = 1, n.roll = out.sample[i]-1, 
										external.forecasts = list(mregfor = fmex, vregfor = fvex, skregfor = fskx,
												shregfor = fshx))
								sig = as.numeric(sigma(f))
								ret = as.numeric(fitted(f))
								if(shaped) shp = as.numeric(shape(f)) else shp = rep(0, out.sample[i])
								if(skewed) skw = as.numeric(skew(f)) else skw = rep(0, out.sample[i])
								if(ghyp) shpgig = rep(coef(fit)["ghlambda"], out.sample[i]) else shpgig = rep(0, out.sample[i])
								rlz = tail(data[rollind[[i]]], out.sample[i])
								# use xts for indexing the forecasts
								y = as.data.frame(cbind(ret, sig, skw, shp, shpgig, rlz))
								rownames(y) = tail(as.character(index[rollind[[i]]]), out.sample[i])
								colnames(y) = c("Mu", "Sigma", "Skew", "Shape", "Shape(GIG)", "Realized")
								if(keep.coef) cf = fit@fit$robust.matcoef else cf = NA
								ans = list(y = y, cf = cf, q = f@forecast$qFor, DiurnalVar = f@model$DiurnalVar, 
										DailyVar = f@model$DailyVar, converge = TRUE, lik = c(likelihood(fit)[1], glik))
							}
						}
						return(ans)})
		}
		forecast = object@forecast
		conv = sapply(tmp, FUN = function(x) x$converge)
		for(i in 1:length(noncidx)){
			if(conv[i]) forecast[[noncidx[i]]] = tmp[[i]]
		}
		if(any(!conv)){
			warning("\nnon-converged estimation windows present...resubsmit object with different solver parameters...")
			noncidx = object@model$noncidx[which(!conv)]
			model = list()
			model$fixARMA = fixARMA
			model$fixGARCH = fixGARCH
			model$fixUBShape = fixUBShape
			model$UBShapeAdd = UBShapeAdd
			model$fixGHlambda = fixGHlambda
			model$compareGARCH = compareGARCH
			model$spec = spec
			model$data = data
			model$index = index
			model$period = period
			model$n.ahead = n.ahead
			model$forecast.length = forecast.length 
			model$n.start = n.start
			model$refit.every = refit.every
			model$n.refits = m
			model$refit.window = refit.window
			model$window.size = window.size
			model$calculate.VaR = calculate.VaR
			model$VaR.alpha = VaR.alpha
			model$keep.coef = keep.coef
			model$noncidx = noncidx
			model$rollind = rollind
			model$out.sample = out.sample
			model$DailyVar = DailyVar
			forecast = forecast
			toc = Sys.time()-tic
			model$elapsed = toc
			ans = new("ACDroll",
					model = model,
					forecast = forecast)
			return( ans )			
		} else{
			noncidx = NULL
			forc = forecast[[1]]$y
			if(m>1){
				for(i in 2:m){
					forc = rbind(forc, forecast[[i]]$y)
				}
			}
			if(keep.coef){
				cf = vector(mode = "list", length = m)
				for(i in 1:m){
					cf[[i]]$index = index[tail(rollind[[i]],1) - out.sample[i]]
					cf[[i]]$coef = forecast[[i]]$cf
				}
			} else{
				cf = NULL
			}
			LL = vector(mode = "list", length = m)
			for(i in 1:m){
				LL[[i]]$index = index[tail(rollind[[i]],1) - out.sample[i]]
				LL[[i]]$log.likelihood = forecast[[i]]$lik
			}
			if(calculate.VaR){
				if(is.null(VaR.alpha)) VaR.alpha = c(0.01, 0.05)
				VaR.matrix = matrix(NA, ncol = length(VaR.alpha)+1, nrow = NROW(forc))
				for(i in 1:length(VaR.alpha)){
					VaR.matrix[,i] = qdist(p = VaR.alpha[i] , mu = forc[,1], sigma = forc[,2], 
							skew = forc[,3], shape = forc[,4], lambda = forc[,5], 
							distribution = distribution)
				}
				VaR.matrix[,length(VaR.alpha)+1] = forc[,6]
				colnames(VaR.matrix) = c(paste("alpha(", round(VaR.alpha,2)*100, "%)",sep=""), "realized")
				VaR.matrix = as.data.frame(VaR.matrix)
				rownames(VaR.matrix) = rownames(forc)
			} else{
				VaR.matrix = NULL
			}
			model = list()
			model$fixARMA = fixARMA
			model$fixGARCH = fixGARCH
			model$fixUBShape = fixUBShape
			model$UBShapeAdd = UBShapeAdd
			model$fixGHlambda = fixGHlambda
			model$compareGARCH = compareGARCH
			model$spec = spec
			model$data = data
			model$index = index
			model$period = period
			model$n.ahead = n.ahead
			model$forecast.length = forecast.length 
			model$n.start = n.start
			model$refit.every = refit.every
			model$n.refits = m
			model$refit.window = refit.window
			model$window.size = window.size
			model$calculate.VaR = calculate.VaR
			model$VaR.alpha = VaR.alpha
			model$keep.coef = keep.coef
			model$noncidx = noncidx
			model$rollind = rollind
			model$out.sample = out.sample
			model$coef = cf
			model$LL = LL
			forecast = list(VaR = VaR.matrix, density = forc)
		}
		toc = Sys.time()-tic
		model$elapsed = toc
		ans = new("ACDroll",
				model = model,
				forecast = forecast)
	} else{
		# do nothing...all converged
		ans = object
	}
	return( ans )
}
################################################################################
# Special Purpose Functions for the intraday multiplicative component sGARCH model
################################################################################
.daily2intraday = function(x, v)
{
	Z = as.POSIXct(format(index(x), format="%Y-%m-%d"))
	Y = xts(rep(0, NROW(x)), index(x))
	T = index(v)
	for(i in 1:length(T)){
		idx = which(Z==T[i])
		Y[idx] = v[i]
	}
	return(Y)
}

.intraday2daily = function(v)
{
	idx = unique(format(index(v), "%Y-%m-%d"))
	idx1 = format(index(v),"%Y-%m-%d")
	idx2 = sapply(idx, function(x) min(which(idx1==x)))
	ans = xts(as.numeric(v[idx2]), as.POSIXct(idx))
	return(ans)
}
.unique_intraday = function(x)
{
	Z = format(index(x), format="%H:%M:%S")
	return(sort(unique(Z)))
}
.unique_time = function(x)
{
	Z = format(index(x), format="%H:%M:%S")
	Y = sort(unique(Z))
	idx = vector(mode="list", length=length(Y))
	for(i in 1:length(Y)){
		idx[[i]] = which(Z==Y[i])
	}
	return(idx)
}
.stime = function(x)
{
	Z = format(index(x), format="%H:%M:%S")
	Y = sort(unique(Z))
	U = as.POSIXct(format(index(x), format="%Y-%m-%d"))
	UX = unique(U)
	idx = vector(mode="list", length=length(UX))
	for(i in 1:length(UX)){
		tmp = which(U==UX[i])
		idx[[i]] = tmp[match(Y, Z[tmp])]
	}
	xidx = lapply(idx, function(i) ifelse(!is.na(i), 1, NA))
	return(xidx)
}
# idx2 = .stime(residuals)
# idx1 = .unique_time(residuals)
# v = .daily2intraday(residuals, dailyvar)
# x = residuals
.diurnal_series_aligned = function(x, v, idx1, idx2)
{
	s = sapply(idx1, function(i) median((x[i]^2)/v[i]))
	sx = as.numeric(unlist(sapply(idx2, function(x) na.omit(s*x))))
	return(sx)
}

.diurnal_series = function(x, v, idx1)
{
	s = sapply(idx1, function(i) median(x[i]^2/v[i]))
	return(s)
}