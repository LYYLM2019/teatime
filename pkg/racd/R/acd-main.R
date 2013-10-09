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
.acdfit = function(spec, data, solver = "ucminf", out.sample = 0, solver.control = list(), 
		fit.control = list(stationarity = 0, fixed.se = 0, scale = 0, n.sim = 2000), 
		skew0 = NULL, shape0 = NULL, cluster = NULL, ...)
{
	tic = Sys.time()
	vmodel = spec@model$vmodel$model
	
	if(is.null(solver.control$trace)) trace = 0 else trace = solver.control$trace
	# default for stationarity is off for ACD models
	if(is.null(fit.control$stationarity)) fit.control$stationarity = FALSE
	if(is.null(fit.control$fixed.se)) fit.control$fixed.se = FALSE
	if(is.null(fit.control$scale)) fit.control$scale = FALSE
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
	# create a temporary environment to store values (deleted at end of function)
	garchenv = new.env(hash = TRUE)
	arglist = list()
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
	if(fit.control$scale) dscale = sd(data) else dscale = 1
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
				return(acdfilter(data = data, spec = spec, out.sample = out.sample))
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
	
	fun = switch(vmodel,
			sGARCH   = .sacdLLH,
			csGARCH  = .csacdLLH)
	
	#fun = switch(vmodel,
	#		sGARCH   = racd:::.sacdLLH,
	#		csGARCH  = racd:::.csacdLLH)
	
	fname = switch(vmodel,
			sGARCH   = "sACD",
			csGARCH  = "csACD")
	
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
	fit$skhEst = arglist$skhEst
	ans = new("ACDfit",
			fit = fit,
			model = model)
	rm(garchenv)
	return(ans)
}


.arfimaxfilteracd = function(model, pars, idx, mexdata, h, tskew, tshape, data, N, garchenv)
{
	#if(model[1] == 0) pars[1,1] = 0
	m = as.integer(N[1])
	T = as.integer(N[2])
	if(length(h) <= 1) {
		h = double(length = T)
	} else{
		h = as.double(h)
	}
	if(length(tskew)<=1){
		tskew = double(T)
	} else{
		tskew = as.double(tskew)
	}
	if(length(tshape)<=1){
		tshape = double(T)
	} else{
		tshape = as.double(tshape)
	}
	
	data = as.double(data)
	# flatten exogenous matrix
	if(model[8]>0){
		xmxreg = matrix( pars[idx[8,1]:idx[8,2]], ncol = model[8] )
		imx =  xmxreg %*%t( matrix( mexdata, ncol = model[8] ) )
		mexdata = as.double(as.vector(mexdata))
	} else{
		mexdata = double(T)
		imx = 0
	}
	if(model[5]>0){
		imx = imx + pars[idx[5]]*(h)
	}
	if(model[6]>0){
		imx = imx + pars[idx[6]]*(tskew)
	}
	if(model[7]>0){
		imx = imx + pars[idx[7]]*(tshape)
	}
	res = double(length = T)
	# this routine is used for the mean residuals to initiate the recursion
	# so we ignore arfima before
	zrf = double(length = T)
	constm = double(length = T)
	condm = double(length = T)
	ans = list()
	if(model[2]>0 | model[3]>0){
		ans = try(.C("arfimaxacdfilterC", model = as.integer(model), pars = as.double(pars), 
						idx = as.integer(idx-1), x = data, res = res, mexdata = mexdata, 
						zrf = zrf, constm = constm, condm = condm, h = h, 
						tskew = tskew, tshape = tshape, m = m, T = T,  
						PACKAGE = "racd"), silent = TRUE)
		if(inherits(ans, "try-error")){
			assign(".csol", 1, envir = garchenv)
			assign(".filtermessage", ans, envir = garchenv)
			res = data - pars[idx[1,1]]
			ans$res = res
			if(model[4]>0)
			{
				ans$zrf = rugarch:::.fracdiff(c(1,rep(0,length(data)-1)), darfima = pars[idx[4]])
				ans$res = rugarch:::.fracdiff(ans$res, darfima = pars[idx[4]])
			}
			if(any(is.na(res))) res[which(is.na(res))]=0
			return(ans)
		} else{
			assign(".csol", 0, envir = garchenv)
			if(model[4]>0)
			{
				ans$zrf = rugarch:::.fracdiff(c(1,rep(0,length(data)-1)), darfima = pars[idx[4]])
				ans$res = rugarch:::.fracdiff(ans$res, darfima = pars[idx[4]])
			}
			if(any(is.na(ans$res))) res[which(is.na(ans$res))]=0
			return(ans)
		}
	} else{
		ans = list()
		ans$res = data -  pars[idx[1,1]] - imx
		ans$zrf = zrf
		if(model[4]>0)
		{
			ans$zrf = rugarch:::.fracdiff(c(1,rep(0,length(data)-1)), darfima = pars[idx[4]])
			ans$res = rugarch:::.fracdiff(ans$res, darfima = pars[idx[4]])
		}
		if(any(is.na(ans$res))) res[which(is.na(ans$res))]=0
		return(ans)
	}
}

.acdfilter = function(spec, data, out.sample = 0, n.old = NULL, skew0 = NULL, shape0 = NULL, ...)
{
	# n.old is optional and indicates the length of the original dataseries (in
	# cases when this represents a dataseries augmented by newer data). The reason
	# for using this is so that the old and new datasets agree since the original
	# recursion uses the sum of the residuals to start the recursion and therefore
	# is influenced by new data. For a small augmentation the values converge after
	# x periods, but it is sometimes preferable to have this option so that there is
	# no forward looking information contaminating the study.
	tic = Sys.time()
	vmodel = spec@model$vmodel$model
	xdata = rugarch:::.extractdata(data)
	data = xdata$data
	index = xdata$index
	period = xdata$period
	origdata = data
	origindex = xdata$index
	T = length(origdata)  - out.sample
	
	if(!is.null(n.old) && n.old>T) stop("\nn.old cannot be greater than length data - out.sample!")
	
	data = origdata[1:T]
	index = origindex[1:T]
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
	ans  = switch(vmodel,
			sGARCH   = .sacdLLH(pars, arglist),
			csGARCH  = .csacdLLH(pars, arglist))
	
	filter = list()
	filter$z = ans$z
	filter$sigma = sqrt(ans$h)
	if(vmodel == "csGARCH") filter$q = ans$q
	filter$residuals = ans$res
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

#---------------------------------------------------------------------------------
# SECTION ACD forecast
#---------------------------------------------------------------------------------
.acdforecast1 = function(fit, n.ahead = 10, n.roll = 0, 
		external.forecasts = list(mregfor = NULL, vregfor = NULL, skxregfor = NULL,
				shxregfor = NULL), m.sim = 1000, cluster = NULL, ...)
{
	data = fit@model$modeldata$data
	Nor = length(as.numeric(data))
	index = fit@model$modeldata$index
	period = fit@model$modeldata$period
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
	fspec = acdspec(variance.model = list(model = model$vmodel$model, 
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
	} else{
		tmp =  xts(c(data[1:(N + fcreq)],0), c(index[1:(N + fcreq)], index[(N + fcreq)]+1))	
	}
	flt = acdfilter(spec = fspec, data = tmp, n.old = N, skew0 = fit@fit$tskew[1], shape0 = fit@fit$tshape[1])
	sigmafilter 	= flt@filter$sigma
	resfilter 		= flt@filter$residuals
	zfilter 		= flt@filter$z
	tskewfilter 	= flt@filter$tskew
	tshapefilter 	= flt@filter$tshape
	tempskewfilter 	= flt@filter$tempskew
	tempshapefilter = flt@filter$tempshape
	
	seriesfor = sigmafor = matrix(NA, ncol = n.roll+1, nrow = n.ahead)
	tskewfor = tshapefor = tempshafor = tempskewfor = matrix(NA, ncol = n.roll+1, nrow = n.ahead)
	seriesfor[1,] = fitted(flt)[(N+1):(N+n.roll+1)]
	sigmafor[1,]  =  sigma(flt)[(N+1):(N+n.roll+1)]
	tskewfor[1,]  =  skew(flt)[(N+1):(N+n.roll+1)]
	tshapefor[1,] = shape(flt)[(N+1):(N+n.roll+1)]
	# n.roll x n.ahead (n.ahead=1 generted by model)
	# n.ahead>1 by simulation
	colnames(seriesfor) = colnames(sigmafor) = as.character(index[N:(N+n.roll)])
	colnames(tskewfor) = colnames(tshapefor) = as.character(index[N:(N+n.roll)])
	rownames(seriesfor) = rownames(sigmafor) = paste("T+", 1:n.ahead, sep="")
	rownames(tskewfor) = rownames(tshapefor) = paste("T+", 1:n.ahead, sep="")
	mx = model$maxOrder
	for(i in 1:(n.roll+1)){
		if(n.ahead>1){
			spec = fspec
			presig     = tail(sigmafilter[1:(N+i-1)],  mx)
			preskew    = tail(tskewfilter[1:(N+i-1)],  mx)
			preshape   = tail(tshapefilter[1:(N+i-1)], mx)
			prereturns = tail(data[1:(N+i-1)],         mx)
			for(j in 2:n.ahead){
				# external forecasts passed as fixed:
				if(modelinc[8]>0 && !is.null(external.forecasts$mregfor)){
					mxsim = matrix(mxf[N+(i+j-1),], ncol = modelinc[8], nrow = m.sim, byrow = TRUE)
					mxsim = replicate(m.sim, mxsim, simplify =  FALSE)
					spec@model$modeldata$mexdata = fspec@model$modeldata$mexdata[1:(N+i-1),,drop=FALSE]
				} else{
					mxsim = NULL
				}
				if(modelinc[17]>0 && !is.null(external.forecasts$vregfor)){
					vxsim = matrix(vxf[N+(i+j-1),], ncol = modelinc[17], nrow = m.sim, byrow = TRUE)
					vxsim = replicate(m.sim, vxsim, simplify =  FALSE)
					spec@model$modeldata$vexdata = fspec@model$modeldata$vexdata[1:(N+i-1),,drop=FALSE]
					
				} else{
					vxsim = NULL
				}
				if(modelinc[25]>0 && !is.null(external.forecasts$skxregfor)){
					skxsim = matrix(skxf[N+(i+j-1),], ncol = modelinc[25], nrow = m.sim, byrow = TRUE)
					skxsim = replicate(m.sim, skxsim, simplify =  FALSE)
					spec@model$modeldata$skxdata = fspec@model$modeldata$skxdata[1:(N+i-1),,drop=FALSE]	
				} else{
					skxsim = NULL
				}
				if(modelinc[31]>0 && !is.null(external.forecasts$shxregfor)){
					shxsim = matrix(shxf[N+(i+j-1),], ncol = modelinc[31], nrow = m.sim, byrow = TRUE)
					shxsim = replicate(m.sim, shxsim, simplify =  FALSE)
					spec@model$modeldata$shxdata = fspec@model$modeldata$shxdata[1:(N+i-1),,drop=FALSE]
				} else{
					shxsim = NULL
				}
				sim = acdpath(fspec, n.sim = 1, n.start = 0, m.sim = m.sim, 
						presigma = presig, preskew = preskew, 
						preshape = preshape, prereturns = prereturns, 
						preresiduals = NA, rseed = NA, 
						mexsimdata = mxsim, vexsimdata = vxsim, 
						skxsimdata = skxsim, shxsimdata = shxsim,
						cluster = cluster)
				seriesfor[j, i] = apply(sim@path$seriesSim, 1, "mean")
				sigmafor[j,  i] = sqrt(apply(sim@path$sigmaSim^2, 1, "mean"))
				tskewfor[j,  i]	= apply(sim@path$skewSim, 1, "mean")
				tshapefor[j, i]	= apply(sim@path$shapeSim,  1, "mean")
				# update the previous values:
				prereturns[mx] = seriesfor[j, i]
				presig[mx]   = sigmafor[j, i]
				preskew[mx]  = tskewfor[j, i]
				preshape[mx] = tshapefor[j, i]
			}
		}
	}
	fcst = list()
	fcst$n.ahead = n.ahead
	fcst$n.roll = n.roll
	fcst$N = N+ns
	fcst$n.start = ns
	fcst$seriesFor = seriesfor
	fcst$sigmaFor  = sigmafor
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


.acdforecast2 = function(spec, data = NULL, n.ahead = 10, n.roll = 0, out.sample = 0,
		external.forecasts = list(mregfor = NULL, vregfor = NULL, skxregfor = NULL,
				shxregfor = NULL), m.sim = 1000, cluster = NULL, skew0 = NULL, 
		shape0 = NULL, ...)
{
	vmodel = spec@model$vmodel$model
	xdata = rugarch:::.extractdata(data)
	data = xdata$data
	index = xdata$index
	period = xdata$period
	Nor = length(as.numeric(data))
	ns = out.sample
	N = Nor - ns
	model = spec@model
	ipars = model$pars
	pars = unlist(model$fixed.pars)
	parnames = names(pars)
	modelnames = rugarch:::.checkallfixed(spec)
	if(is.na(all(match(modelnames, parnames), 1:length(modelnames)))) {
		cat("\nacdforecast-->error: parameters names do not match specification\n")
		cat("Expected Parameters are: ")
		cat(paste(modelnames))
		cat("\n")
		stop("Exiting", call. = FALSE)
	}
	# once more into the spec
	# NB Any changes made to the spec are not preserved once we apply set fixed
	setfixed(spec)<-as.list(pars)
	
	
	modelinc = model$modelinc
	idx = model$pidx
	if( n.roll > ns ) stop("\nacdforecast-->error: n.roll must not be greater than out.sample!")

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
	fspec = spec
	fspec@model$modeldata$mexdata = mxf
	fspec@model$modeldata$vexdata = vxf
	fspec@model$modeldata$skxdata = skxf
	fspec@model$modeldata$shxdata = shxf
	# Generate the 1 extra ahead forecast
	if((n.ahead+n.roll)<=out.sample){
		tmp =  xts(data[1:(N + fcreq)], index[1:(N + fcreq)])
	} else{
		tmp =  xts(c(data[1:(N + fcreq)],0), c(index[1:(N + fcreq)], index[(N + fcreq)]+1))	
	}
	flt = acdfilter(spec = fspec, data = tmp, n.old = N, skew0 = skew0, shape0 = shape0)
	sigmafilter 	= flt@filter$sigma
	resfilter 		= flt@filter$residuals
	zfilter 		= flt@filter$z
	tskewfilter 	= flt@filter$tskew
	tshapefilter 	= flt@filter$tshape
	tempskewfilter 	= flt@filter$tempskew
	tempshapefilter = flt@filter$tempshape
	
	seriesfor = sigmafor = matrix(NA, ncol = n.roll+1, nrow = n.ahead)
	tskewfor = tshapefor = tempshafor = tempskewfor = matrix(NA, ncol = n.roll+1, nrow = n.ahead)
	seriesfor[1,] = fitted(flt)[(N+1):(N+n.roll+1)]
	sigmafor[1,]  =  sigma(flt)[(N+1):(N+n.roll+1)]
	tskewfor[1,]  =  skew(flt)[(N+1):(N+n.roll+1)]
	tshapefor[1,] = shape(flt)[(N+1):(N+n.roll+1)]
	# n.roll x n.ahead (n.ahead=1 generted by model)
	# n.ahead>1 by simulation
	colnames(seriesfor) = colnames(sigmafor) = as.character(index[N:(N+n.roll)])
	colnames(tskewfor) = colnames(tshapefor) = as.character(index[N:(N+n.roll)])
	rownames(seriesfor) = rownames(sigmafor) = paste("T+", 1:n.ahead, sep="")
	rownames(tskewfor) = rownames(tshapefor) = paste("T+", 1:n.ahead, sep="")
	mx = model$maxOrder
	for(i in 1:(n.roll+1)){
		if(n.ahead>1){
			spec = fspec
			presig     = tail(sigmafilter[1:(N+i-1)],  mx)
			preskew    = tail(tskewfilter[1:(N+i-1)],  mx)
			preshape   = tail(tshapefilter[1:(N+i-1)], mx)
			prereturns = tail(data[1:(N+i-1)],         mx)
			for(j in 2:n.ahead){
				# external forecasts passed as fixed:
				if(modelinc[8]>0 && !is.null(external.forecasts$mregfor)){
					mxsim = matrix(mxf[N+(i+j-1),], ncol = modelinc[8], nrow = m.sim, byrow = TRUE)
					mxsim = replicate(m.sim, mxsim, simplify =  FALSE)
					spec@model$modeldata$mexdata = fspec@model$modeldata$mexdata[1:(N+i-1),,drop=FALSE]
				} else{
					mxsim = NULL
				}
				if(modelinc[17]>0 && !is.null(external.forecasts$vregfor)){
					vxsim = matrix(vxf[N+(i+j-1),], ncol = modelinc[17], nrow = m.sim, byrow = TRUE)
					vxsim = replicate(m.sim, vxsim, simplify =  FALSE)
					spec@model$modeldata$vexdata = fspec@model$modeldata$vexdata[1:(N+i-1),,drop=FALSE]
					
				} else{
					vxsim = NULL
				}
				if(modelinc[25]>0 && !is.null(external.forecasts$skxregfor)){
					skxsim = matrix(skxf[N+(i+j-1),], ncol = modelinc[25], nrow = m.sim, byrow = TRUE)
					skxsim = replicate(m.sim, skxsim, simplify =  FALSE)
					spec@model$modeldata$skxdata = fspec@model$modeldata$skxdata[1:(N+i-1),,drop=FALSE]	
				} else{
					skxsim = NULL
				}
				if(modelinc[31]>0 && !is.null(external.forecasts$shxregfor)){
					shxsim = matrix(shxf[N+(i+j-1),], ncol = modelinc[31], nrow = m.sim, byrow = TRUE)
					shxsim = replicate(m.sim, shxsim, simplify =  FALSE)
					spec@model$modeldata$shxdata = fspec@model$modeldata$shxdata[1:(N+i-1),,drop=FALSE]
				} else{
					shxsim = NULL
				}
				sim = acdpath(fspec, n.sim = 1, n.start = 0, m.sim = m.sim, 
						presigma = presig, preskew = preskew, 
						preshape = preshape, prereturns = prereturns, 
						preresiduals = NA, rseed = NA, 
						mexsimdata = mxsim, vexsimdata = vxsim, 
						skxsimdata = skxsim, shxsimdata = shxsim,
						cluster = cluster)
				seriesfor[j, i] = apply(sim@path$seriesSim, 1, "mean")
				sigmafor[j,  i] = sqrt(apply(sim@path$sigmaSim^2, 1, "mean"))
				tskewfor[j,  i]	= apply(sim@path$skewSim, 1, "mean")
				tshapefor[j, i]	= apply(sim@path$shapeSim,  1, "mean")
				# update the previous values:
				prereturns[mx] = seriesfor[j, i]
				presig[mx]   = sigmafor[j, i]
				preskew[mx]  = tskewfor[j, i]
				preshape[mx] = tshapefor[j, i]
			}
		}
	}
	
	fcst = list()
	fcst$n.ahead = n.ahead
	fcst$n.roll = n.roll
	fcst$N = N+ns
	fcst$n.start = ns
	fcst$seriesFor = seriesfor
	fcst$sigmaFor  = sigmafor
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
#---------------------------------------------------------------------------------
# SECTION ACD roll
#---------------------------------------------------------------------------------
.acdroll = function(spec, data, n.ahead = 1, forecast.length = 500, n.start = NULL, 
    refit.every = 25, refit.window = c("recursive", "moving"), 
    window.size = NULL, solver = "ucminf", fit.control = list(), 
    solver.control = list(), calculate.VaR = TRUE, VaR.alpha = c(0.01, 
        0.05), cluster = NULL, keep.coef = TRUE, fixARMA = TRUE, fixGARCH = TRUE, 
	fixUBShape = TRUE, UBShapeAdd = 0, fixGHlambda = TRUE, compareGARCH = c("LL", "none"),
	...)
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
	# distribution
	distribution = spec@model$dmodel$model
	if(any(distribution==c("snorm","sstd","sged","jsu","nig","ghyp","ghst"))) skewed = TRUE else skewed = FALSE
	if(any(distribution==c("std","sstd","ged","sged","jsu","nig","ghyp","ghst"))) shaped = TRUE else shaped = FALSE
	if(any(distribution==c("ghyp"))) ghyp = TRUE else ghyp = FALSE
	
	gspec = .spec2GARCH(spec)
	if( !is.null(cluster) ){
		clusterEvalQ(cl = cluster, library(racd))
		clusterExport(cluster, c("data", "index", "s","refit.every", 
						"keep.coef", "shaped", "skewed", "ghyp", "gspec", "fixARMA",
						"fixGARCH", "fixUBShape", "UBShapeAdd", "fixGHlambda","compareGARCH",
						"rollind", "spec", "out.sample", "mex", "vex", "skex", "shex",
						"solver", "solver.control", "fit.control"), envir = environment())
		if(mex)  clusterExport(cluster, c("mexdata"), envir = environment())
		if(vex)  clusterExport(cluster, c("vexdata"), envir = environment())
		if(skex) clusterExport(cluster, c("skdata"), envir = environment())
		if(shex) clusterExport(cluster, c("shdata"), envir = environment())
		tmp = parLapplyLB(cl = cluster, 1:m, fun = function(i){
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
					gfit = ugarchfit(xspec, data[rollind[[i]]], out.sample = out.sample[i], 
							solver = "hybrid")
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
					fit = try(acdfit(zspec, xts(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i], 
									solver = solver, solver.control = solver.control, 
									fit.control = fit.control, shape0 = shape0, skew0 = skew0), silent=TRUE)
					if(inherits(fit, 'try-error') || convergence(fit)!=0 || is.null(fit@fit$cvar)){
						ans = list(y = NA, cf = NA, converge = FALSE, lik = c(NA, glik))
					} else{
						# compare GARCH likelihood with ACD model and reject if lik less than
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
							ans = list(y = y, cf = cf, converge = TRUE, lik = c(likelihood(fit)[1], glik))
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
			gfit = ugarchfit(xspec, data[rollind[[i]]], out.sample = out.sample[i], 
					solver = "hybrid")
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
			fit = try(acdfit(zspec, xts(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i], 
							solver = solver, solver.control = solver.control, 
							fit.control = fit.control, shape0 = shape0, skew0 = skew0), silent=TRUE)
			if(inherits(fit, 'try-error') || convergence(fit)!=0 || is.null(fit@fit$cvar)){
				tmp[[i]] = list(y = NA, cf = NA, converge = FALSE, lik = c(NA, glik))
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
					if(skewed) skw = as.numeric(skew(f)) else skw = rep(0, out.sample[i])
					if(ghyp) shpgig = rep(coef(fit)["ghlambda"], out.sample[i]) else shpgig = rep(0, out.sample[i])
					rlz = tail(data[rollind[[i]]], out.sample[i])
					# use xts for indexing the forecasts
					y = as.data.frame(cbind(ret, sig, skw, shp, shpgig, rlz))
					rownames(y) = tail(as.character(index[rollind[[i]]]), out.sample[i])
					colnames(y) = c("Mu", "Sigma", "Skew", "Shape", "Shape(GIG)", "Realized")
					if(keep.coef) cf = fit@fit$robust.matcoef else cf = NA
					tmp[[i]] = list(y = y, cf = cf, converge = TRUE, lik = c(likelihood(fit)[1], glik))
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


.acdresumeroll = function(object, spec = NULL, solver = "ucminf", fit.control = list(), 
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
		# distribution
		distribution = model$spec@model$dmodel$model
		if(any(distribution==c("snorm","sstd","sged","jsu","nig","ghyp","ghst"))) skewed = TRUE else skewed = FALSE
		if(any(distribution==c("std","sstd","ged","sged","jsu","nig","ghyp","ghst"))) shaped = TRUE else shaped = FALSE
		if(any(distribution==c("ghyp"))) ghyp = TRUE else ghyp = FALSE
		if( !is.null(cluster) ){
			clusterEvalQ(cl = cluster, library(racd))
			clusterExport(cluster, c("data", "index","s","refit.every",
							"keep.coef", "shaped", "skewed", "ghyp", "gspec", "fixARMA",
							"fixGARCH", "fixUBShape", "UBShapeAdd", "fixGHlambda","compareGARCH",
							"rollind", "spec", "out.sample", "mex", "vex", "skex", "shex",
							"noncidx", "solver", "solver.control", "fit.control"),
					envir = environment())
			if(mex)  clusterExport(cluster,  c("mexdata"), envir = environment())
			if(vex)  clusterExport(cluster, c("vexdata"), envir = environment())
			if(skex) clusterExport(cluster, c("skdata"), envir = environment())
			if(shex) clusterExport(cluster, c("shdata"), envir = environment())
			tmp = parLapplyLB(cl = cluster, as.list(noncidx), fun = function(i){
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
						gfit = ugarchfit(xspec, data[rollind[[i]]], out.sample = out.sample[i], 
								solver = "hybrid")
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
						fit = try(acdfit(spec, xts(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i], 
										solver = solver, solver.control = solver.control, 
										fit.control = fit.control, shape0 = shape0, skew0 = skew0), silent=TRUE)
						if(inherits(fit, 'try-error') || convergence(fit)!=0 || is.null(fit@fit$cvar)){
							ans = list(y = NA, cf = NA, converge = FALSE, lik = c(NA, glik))
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
								ans = list(y = y, cf = cf, converge = TRUE, lik = c(likelihood(fit)[1], glik))
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
						gfit = ugarchfit(xspec, data[rollind[[i]]], out.sample = out.sample[i], 
								solver = "hybrid")
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
						fit = try(acdfit(spec, xts(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i], 
										solver = solver, solver.control = solver.control, 
										fit.control = fit.control, shape0 = shape0, skew0 = skew0), silent=TRUE)
						if(inherits(fit, 'try-error') || convergence(fit)!=0 || is.null(fit@fit$cvar)){
							ans = list(y = NA, cf = NA, converge = FALSE, lik = c(NA, glik))
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
								ans = list(y = y, cf = cf, converge = TRUE, lik = c(likelihood(fit)[1], glik))
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

.spec2GARCH = function(spec){
	modelinc = spec@model$modelinc
	gspec = ugarchspec(mean.model=list(include.mean = as.logical(modelinc[1]), 
					armaOrder = modelinc[2:3], external.regressors = spec@model$modeldata$mexdata,
					archm = spec@model$mmodel$archm, arfima = spec@model$mmodel$arfima),
			variance.model=list(model = spec@model$vmodel$model, garchOrder = modelinc[10:11],
					variance.targeting = spec@model$vmodel$variance.targeting,
					external.regressors = spec@model$modeldata$vexdata),
			distribution.model = spec@model$dmodel$model, fixed.pars = spec@model$fixed.pars)
	setbounds(gspec)<-list(skew = spec@model$sbounds[1:2], shape = spec@model$sbounds[3:4])
	return(gspec)
}
