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
.csacdLLH = function(pars, arglist)
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
	if( !is.null(arglist$n.old) ){
		rx = .arfimaxfilteracd(modelinc, ipars[,1], idx, mexdata = arglist$mexdata[1:Nx, , drop=FALSE], 
				h = hm, tskew = 0, tshape = 0, data = data[1:Nx], N = c(m, Nx), arglist$garchenv)
		res2 = rx$res
		res2[is.na(res2) | !is.finite(res2) | is.nan(res2)] = 0		
		mvar = mean(res2*res2)
	} else{
		mvar = mean(res*res)
	}
	# sgarch persistence value
	mexdata = as.double(as.vector(arglist$mexdata))
	vexdata = as.double(as.vector(arglist$vexdata))
	skxdata = as.double(as.vector(arglist$skxdata))
	shxdata = as.double(as.vector(arglist$shxdata))
	persist = (sum(ipars[idx["alpha",1]:idx["alpha",2],1]) + sum(ipars[idx["beta",1]:idx["beta",2],1]))
	# unconditional sigma value
	mvar = mean(res*res)
	# no variance targeting in component GARCH
	ipars[idx["omega",1],1] = max(eps, ipars[idx["omega",1],1]) 
	hEst = mvar
	if(is.na(hEst) | !is.finite(hEst) | is.nan(hEst)) hEst = var(data) * (1 - persist)
	assign("racd_ipars", ipars, envir = arglist$garchenv)
	if(fit.control$stationarity == 1 && modelinc[17] == 0){
		xcond = .csacdcon(pars, arglist)
		if(any(xcond<0)){
			if(arglist$pmode!=1){
				return(llh = get("racd_llh", arglist$garchenv) + 0.1*(abs(get("racd_llh", arglist$garchenv))))
			} else{
				return(llh=1e10)
			}
		}
	}
	#if(modelinc[8]>0) mexdata = as.double(as.vector(mexdata)) else mexdata = double(1)
	#if(modelinc[17]>0) vexdata = as.double(as.vector(vexdata)) else vexdata = double(1)
	#if(modelinc[25]>0) skxdata = as.double(as.vector(skxdata)) else skxdata = double(1)
	#if(modelinc[30]>0) shxdata = as.double(as.vector(shxdata)) else shxdata = double(1)
	
	sbounds = model$sbounds
	skhEst = model$skhEst
	
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
	
	ans = try(.C("csacdfilterC",
					model = as.integer(modelinc), 
					pars = as.double(ipars[,1]), 
					idx = as.integer(idx[,1]-1), 
					hEst = as.double(hEst), 
					x = as.double(data), 
					res = as.double(res), 
					e = double(T), 
					mexdata = as.double(mexdata), 
					vexdata = as.double(vexdata), 
					zrf = as.double(zrf),
					constm = double(T), 
					condm = double(T), 
					m = as.integer(m), 
					T = as.integer(T),
					h = double(T),
					q = double(T),
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
	q = ans$q
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
			all = list(llh = llh, h = h, q = q, res = res, z = z, kappa = kappa, 
					tskew = tskew, tshape = tshape, tempshape = tempshape, 
					tempskew = tempskew, LHT = LHT))
	return( ans )
}

.csacdpath = function(spec, n.sim = 1000, n.start = 0, m.sim = 1, presigma = NA, 
		prereturns = NA, preresiduals = NA, preskew = NA, preshape = NA, rseed = NA,  
		mexsimdata = NULL, vexsimdata = NULL, skxsimdata = NULL, shxsimdata = NULL, 
		cluster = NULL, ...)
{
	if(is.na(rseed[1])){
		sseed = as.integer(runif(m.sim,0,as.integer(Sys.time())))
	} else{
		if(length(rseed) != m.sim) stop("\uacdsim-->error: rseed must be of length m.sim!\n")
		sseed = rseed[1:m.sim]
	}
	n = n.sim + n.start
	model = spec@model
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
	if(modelinc[30]>0) {
		shexdata = matrix(model$modeldata$shxdata, ncol = modelinc[30]) 
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
		stop("\nacdpath-->error: presigma cannot be NA.")
	}
	if(!is.na(prereturns[1])){
		prereturns = as.vector(prereturns)
		if(length(prereturns)<m) stop(paste("\nuacdsim-->error: prereturns must be of length ", m, sep=""))
	} else{
		prereturns = as.numeric(rep(uncmean(spec), m))
	}
	if(!is.na(preresiduals[1])){
		preresiduals = as.vector(preresiduals)
		if(length(preresiduals)<m) stop(paste("\nuacdsim-->error: preresiduals must be of length ", m, sep=""))
		preres = matrix(preresiduals, nrow = m, ncol = m.sim)
	}
	preq = rep(ipars[idx["omega",1],1]/(1 - ipars[idx["eta1",1],1]), m)
	
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
			pretempskew = rep(ipars[idx["skcons",1],1], m)
			pretskew = logtransform(pretempskew, sbounds[1], sbounds[2])
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
	if(model$modelinc[26]>0)
	{
		if( is.na(preshape[1]) ){
			# The tempshape[1] is the transformed shape parameter of the 
			# non-time varying model from which we initiated the original fit.
			pretempshape = rep(ipars[idx["shcons",1],1], m)
			pretshape = exptransform(pretempshape, sbounds[3], sbounds[4], rate = sbounds[5])
		} else{
			pretempshape = exptransform(tail(preshape, m), sbounds[3], sbounds[4], rate = sbounds[5], inverse = TRUE)
			pretshape = tail(preshape, m)
		}
	}
	if(model$modelinc[17]>0){
		tshape = rep(ipars["shape", 1], n+m)
	} else{
		tshape = c(pretshape, rep(0, n))
	}
	# input vectors/matrices
	h = c(presigma^2, rep(0, n))
	q = c(preq, rep(0, n))
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
	if(is.na(preresiduals[1])){
		preres = z[1:m, , drop = FALSE]*matrix(presigma, ncol = m.sim, nrow = m)
	}
	res = rbind( preres, matrix(0, ncol = m.sim, nrow = n) )
	
	# outpus matrices
	sigmaSim  = matrix(0, ncol = m.sim, nrow = n.sim)
	qSim 	  = matrix(0, ncol = m.sim, nrow = n.sim)
	seriesSim = matrix(0, ncol = m.sim, nrow = n.sim)
	residSim  = matrix(0, ncol = m.sim, nrow = n.sim)
	skewSim 		= matrix(0, ncol = m.sim, nrow = n.sim)
	tempskewSim 	= matrix(0, ncol = m.sim, nrow = n.sim)
	shapeSim 		= matrix(0, ncol = m.sim, nrow = n.sim)
	tempshapeSim 	= matrix(0, ncol = m.sim, nrow = n.sim)
	zSim 			= matrix(0, ncol = m.sim, nrow = n.sim)
	
	
	if(!is.null(cluster)){
		parallel::clusterEvalQ(cluster, require(racd))
		parallel::clusterExport(cluster, c("modelinc", "ipars", "idx", "h", "q", "res",
						"tmpskew", "tmpshape", "tskew", "tshape", "sbounds", "sseed",
						"vexsim", "skexsim", "shexsim", "n", "m", "constm", "mexsim"), envir = environment())
		S = parallel::parLapply(cluster, 1:m.sim, function(i){
					set.seed(sseed[i])
					tmp = try(.C("csacdsimC", model = as.integer(modelinc), pars = as.double(ipars[,1]), 
									idx = as.integer(idx[,1]-1), h = as.double(h), q = as.double(q), 
									z = as.double(z[,i]), res = as.double(res[,i]), e = as.double(res[,i]*res[,i]), 
									tempskew = as.double(tmpskew), tempshape = as.double(tmpshape), 
									tskew = as.double(tskew), tshape = as.double(tshape), 
									sbounds = as.double(sbounds), vexdata = as.double(vexsim[[i]]), 
									skxreg = as.double(skexsim[[i]]), shxreg = as.double(shexsim[[i]]),
									T = as.integer(n+m), m = as.integer(m),
									PACKAGE = "racd"), silent = TRUE)
					if(modelinc[8]>0){
						mxreg = matrix( ipars[idx["mxreg",1]:idx["mxreg",2], 1], ncol = modelinc[8] )
						constm[,i] = constm[,i] + mxreg %*%t( matrix( mexsim[[i]], ncol = modelinc[8] ) )
					}
					if(modelinc[5]>0) constm[,i] = constm[,i] + ipars[idx["archm",1]:idx["archm",2], 1]*(sqrt(tmp$h)^modelinc[5])
					if(modelinc[4]>0){
						fres = c(tmp$res[(m+1):(n+m)], if(modelinc[3]>0) rep(0, modelinc[3]) else NULL)
						ans2 = rugarch:::.arfimaxsim(modelinc[1:5], ipars, idx, constm[1:n, i], fres, T = n)
						seriesSim = head(ans2$series, n.sim)
					} else{
						ans2 = rugarch:::.armaxsim(modelinc[1:5], ipars = ipars, idx = idx, constm = constm[,i],  
								x = x, res = tmp$res, T = n + m, m)
						seriesSim = ans2$x[(n.start + m + 1):(n+m)]
					}
					ret = cbind(seriesSim,  tail(sqrt(tmp$h), n.sim), tail(tmp$q, n.sim), 
							tail(tmp$res, n.sim), tail(tmp$tskew, n.sim), 
							tail(tmp$tshape, n.sim), tail(tmp$z, n.sim))
					return(ret)
				})
		seriesSim = sapply(S, function(x) x[,1])
		sigmaSim = sapply(S, function(x) x[,2])
		qSim[,i] = sapply(S, function(x) x[,3])
		residSim = sapply(S, function(x) x[,4])
		skewSim = sapply(S, function(x) x[,5])
		shapeSim = sapply(S, function(x) x[,6])
		zSim = sapply(S, function(x) x[,7])
	} else{
		for(i in 1:m.sim){
			set.seed(sseed[i])
			tmp = try(.C("csacdsimC", model = as.integer(modelinc), pars = as.double(ipars[,1]), 
							idx = as.integer(idx[,1]-1), h = as.double(h), q = as.double(q), 
							z = as.double(z[,i]), res = as.double(res[,i]), 
							e = as.double(res[,i]*res[,i]), 
							tempskew = as.double(tmpskew), tempshape = as.double(tmpshape), 
							tskew = as.double(tskew), tshape = as.double(tshape), 
							sbounds = as.double(sbounds), vexdata = as.double(vexsim[[i]]), 
							skxreg = as.double(skexsim[[i]]), shxreg = as.double(shexsim[[i]]),
							T = as.integer(n+m), m = as.integer(m),
							PACKAGE = "racd"), silent = TRUE)
			if(modelinc[8]>0){
				mxreg = matrix( ipars[idx["mxreg",1]:idx["mxreg",2], 1], ncol = modelinc[8] )
				constm[,i] = constm[,i] + mxreg %*%t( matrix( mexsim[[i]], ncol = modelinc[8] ) )
			}
			if(modelinc[5]>0) constm[,i] = constm[,i] + ipars[idx["archm",1]:idx["archm",2], 1]*(sqrt(tmp$h)^modelinc[5])
			if(modelinc[4]>0){
				fres = c(tmp$res[(m+1):(n+m)], if(modelinc[3]>0) rep(0, modelinc[3]) else NULL)
				ans2 = rugarch:::.arfimaxsim(modelinc[1:5], ipars, idx, constm[1:n, i], fres, T = n)
				seriesSim[,i] = head(ans2$series, n.sim)
			} else{
				ans2 = rugarch:::.armaxsim(modelinc[1:5], ipars = ipars, idx = idx, constm = constm[,i],  
						x = x, res = tmp$res, T = n + m, m)
				seriesSim[,i] = ans2$x[(n.start + m + 1):(n+m)]
			}
			sigmaSim[,i] 	 = tail(sqrt(tmp$h), n.sim)
			qSim[,i]         = tail(tmp$q, n.sim)
			residSim[,i] 	 = tail(tmp$res, n.sim)
			skewSim[,i] 	 = tail(tmp$tskew, n.sim)
			shapeSim[,i] 	 = tail(tmp$tshape, n.sim)
			zSim[,i] 		 = tail(tmp$z, n.sim)
		}
	}
	sim = list(sigmaSim = sigmaSim, qSim = qSim, seriesSim = seriesSim, 
			residSim = residSim, skewSim = skewSim, shapeSim = shapeSim, 
			zSim = zSim)
	sim$n.sim  = n.sim
	sim$m.sim  = m.sim
	
	sol = new("ACDpath",
			path = sim,
			model = model,
			seed = as.integer(sseed))
	return(sol)
}


.csacdsim = function(fit, n.sim = 1000, n.start = 0, m.sim = 1, presigma = NA, 
		prereturns = NA, preresiduals = NA, preskew = NA, preshape = NA, rseed = NA,  
		mexsimdata = NULL, vexsimdata = NULL, skxsimdata = NULL, shxsimdata = NULL, ...)
{
	if(is.na(rseed[1])){
		sseed = as.integer(runif(m.sim,0,as.integer(Sys.time())))
	} else{
		if(length(rseed) != m.sim) stop("\uacdsim-->error: rseed must be of length m.sim!\n")
		sseed = rseed[1:m.sim]
	}
	n = n.sim + n.start
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
	if(modelinc[30]>0) {
		shexdata = matrix(model$modeldata$shxdata, ncol = modelinc[30]) 
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
	}
	preq = tail(fit@fit$q, m)
	
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
	if(model$modelinc[26]>0)
	{
		if( is.na(preshape[1]) ){
			# The tempshape[1] is the transformed shape parameter of the 
			# non-time varying model from which we initiated the original fit.
			pretempshape = tail(fit@fit$tempshape, m)
			pretshape = tail(fit@fit$tshape, m)
		} else{
			pretempshape = exptransform(tail(preshape, m), sbounds[3], sbounds[4], rate = sbounds[5], inverse = TRUE)
			pretshape = tail(preshape, m)
		}
	}
	if(model$modelinc[19]>0){
		tshape = rep(ipars["shape", 1], n+m)
	} else{
		tshape = c(pretshape, rep(0, n))
	}
	# input vectors/matrices
	h = c(presigma^2, rep(0, n))
	q = c(preq, rep(0, n))
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
	if(is.na(preresiduals[1])){
		preres = z[1:m, , drop = FALSE]*matrix(presigma, ncol = m.sim, nrow = m)
	}
	res = rbind( preres, matrix(0, ncol = m.sim, nrow = n) )
	
	# outpus matrices
	sigmaSim  = matrix(0, ncol = m.sim, nrow = n.sim)
	qSim 	  = matrix(0, ncol = m.sim, nrow = n.sim)
	seriesSim = matrix(0, ncol = m.sim, nrow = n.sim)
	residSim  = matrix(0, ncol = m.sim, nrow = n.sim)
	skewSim 		= matrix(0, ncol = m.sim, nrow = n.sim)
	tempskewSim 	= matrix(0, ncol = m.sim, nrow = n.sim)
	shapeSim 		= matrix(0, ncol = m.sim, nrow = n.sim)
	tempshapeSim 	= matrix(0, ncol = m.sim, nrow = n.sim)
	zSim 			= matrix(0, ncol = m.sim, nrow = n.sim)
	
	for(i in 1:m.sim){
		set.seed(sseed[i])
		tmp = try(.C("csacdsimC", model = as.integer(modelinc), pars = as.double(ipars[,1]), 
						idx = as.integer(idx[,1]-1), h = as.double(h), q = as.double(q), 
						z = as.double(z[,i]), res = as.double(res[,i]), 
						e = as.double(res[,i]*res[,i]), 
						tempskew = as.double(tmpskew), tempshape = as.double(tmpshape), 
						tskew = as.double(tskew), tshape = as.double(tshape), 
						sbounds = as.double(sbounds), vexdata = as.double(vexsim[[i]]), 
						skxreg = as.double(skexsim[[i]]), shxreg = as.double(shexsim[[i]]),
						T = as.integer(n+m), m = as.integer(m),
						PACKAGE = "racd"), silent = TRUE)
		if(modelinc[8]>0){
			mxreg = matrix( ipars[idx["mxreg",1]:idx["mxreg",2], 1], ncol = modelinc[8] )
			constm[,i] = constm[,i] + mxreg %*%t( matrix( mexsim[[i]], ncol = modelinc[8] ) )
		}
		if(modelinc[5]>0) constm[,i] = constm[,i] + ipars[idx["archm",1]:idx["archm",2], 1]*(sqrt(tmp$h)^modelinc[5])
		if(modelinc[4]>0){
			fres = c(tmp$res[(m+1):(n+m)], if(modelinc[3]>0) rep(0, modelinc[3]) else NULL)
			ans2 = rugarch:::.arfimaxsim(modelinc[1:5], ipars, idx, constm[1:n, i], fres, T = n)
			seriesSim[,i] = head(ans2$series, n.sim)
		} else{
			ans2 = rugarch:::.armaxsim(modelinc[1:5], ipars = ipars, idx = idx, constm = constm[,i],  
					x = x, res = tmp$res, T = n + m, m)
			seriesSim[,i] = ans2$x[(n.start + m + 1):(n+m)]
		}
		sigmaSim[,i] 	 = tail(sqrt(tmp$h), n.sim)
		qSim[,i]         = tail(tmp$q, n.sim)
		residSim[,i] 	 = tail(tmp$res, n.sim)
		skewSim[,i] 	 = tail(tmp$tskew, n.sim)
		shapeSim[,i] 	 = tail(tmp$tshape, n.sim)
		zSim[,i] 		 = tail(tmp$z, n.sim)
	}
	sim = list(sigmaSim = sigmaSim, qSim = qSim, seriesSim = seriesSim, 
			residSim = residSim, skewSim = skewSim, shapeSim = shapeSim, 
			zSim = zSim)
	sim$n.sim  = n.sim
	sim$m.sim  = m.sim
	
	sol = new("ACDsim",
			simulation = sim,
			model = model,
			seed = as.integer(sseed))
	return(sol)
}

# Constraints from Section 3.1 of Engle and Lee Paper
.csacdcon = function(pars, arglist){
	ipars = arglist$ipars
	estidx = arglist$estidx
	idx = arglist$model$pidx
	ipars[estidx, 1] = pars
	modelinc = arglist$model$modelinc
	eta1 = ipars["eta11",1]
	eta2 = ipars["eta21",1]
	if(modelinc[11]>0) beta = sum(ipars[idx["beta",1]:idx["beta",2],1]) else beta = 0
	if(modelinc[10]>0) alpha = sum(ipars[idx["alpha",1]:idx["alpha",2],1]) else alpha = 0
	# \rho - alpha - beta > 0
	c1 = eta1 - alpha - beta
	# beta - \phi > 0
	c2 = beta - eta2
	return(c(c1, c2))
	#return( 1-  ( (alpha + beta)*(1-eta1) + eta1 ) )
}
