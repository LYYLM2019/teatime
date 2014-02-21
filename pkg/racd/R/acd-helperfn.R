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
.acdmakefitmodel = function(acdmodel, f, T, m, timer, convergence, message, hess, arglist)
{
	# Turn Stationarity Off for numerical derivative calculation
	fit.control = arglist$fit.control
	fit.control$stationarity = arglist$fit.control$stationarity = 0
	data = arglist$data
	ipars = arglist$ipars
	model = arglist$model
	estidx = arglist$estidx
	idx = model$pidx
	arglist$returnType = "llh"
	fit = vector(mode = "list")
	if(is.null(hess)){
		#fit$hessian = rugarch:::.hessian2sidedcpp(f, ipars[estidx, 1], arglist = arglist)
		fit$hessian = numDeriv::hessian(f, x = ipars[estidx, 1], arglist = arglist)
		E = eigen(fit$hessian)$values
		# approx. number of decimal places lost to roundoff/numerical estimation error
		condH = log10(max(E)/min(E))
	} else{
		fit$hessian = hess
		E = eigen(fit$hessian)$values
		condH = log10(max(E)/min(E))
	}
	fit$cvar = try(solve(fit$hessian), silent = TRUE)
	# (might also fail in which case user will see an error about the inversion failure)
	if(inherits(fit$cvar, "try-error")){
		#fit$hessian = numDeriv::hessian(f, ipars[estidx,1], method = "Richardson", arglist = arglist)
		warning("\nracd-->warning: failed to invert hessian\n")
		fit$cvar = NULL
	}
	arglist$returnType = "all"
	temp = f(pars = ipars[estidx, 1],	arglist = arglist)
	if(acdmodel == "apACD"){
		fit$var = temp$h^2
		fit$sigma = temp$h
	} else{
		fit$var = abs(temp$h)
		fit$sigma = sqrt(abs(temp$h))
	}
	if(acdmodel == "csACD") fit$q = temp$q
	fit$condH = condH
	fit$z = temp$z
	fit$LLH = -temp$llh
	fit$log.likelihoods = temp$LHT
	fit$residuals = temp$res
	fit$tskew = temp$tskew
	fit$tshape = temp$tshape
	fit$tempskew = temp$tempskew
	fit$tempshape = temp$tempshape
	
	if(sum(ipars[,2])>0){
		pall = ipars[estidx | as.logical(ipars[,2]==1), 1]
		fixed = match(rownames(ipars[ipars[,2]==1, , drop = FALSE]), names(pall))
		fixedn = length(fixed)
		fNA = rep(NA, fixedn)
		nfixedn = length(pall) - fixedn
		fit$coef = pall
		if(is.null(fit$cvar)){
			fit$se.coef = rep(NA, nfixedn)
			fit$tval = rep(NA, nfixedn)
			# change here
			fit$matcoef = matrix(NA, ncol = 4, nrow = length(pall))
			fit$matcoef[-fixed,] = cbind(ipars[estidx, 1], fit$se.coef,
					fit$tval, rep(NA, nfixedn))
			fit$matcoef[fixed,] = cbind(fit$coef[fixed], fNA, fNA, fNA)
			fit$robust.se.coef = rep(NA, nfixedn)
			fit$robust.tval = rep(NA, nfixedn)
			fit$robust.matcoef = matrix(NA, ncol = 4, nrow = length(fit$coef))
			fit$robust.matcoef[-fixed,] = cbind(fit$coef[-fixed], fit$robust.se.coef, 
					fit$robust.tval, rep(NA, nfixedn))
			fit$robust.matcoef[fixed,] = cbind(fit$coef[fixed], fNA, fNA, fNA)
			fit$hessian.message = "failed to invert hessian"
		} else{
			arglist$returnType = "LHT"
			tmp = rugarch:::robustvcv(fun = f, pars = ipars[estidx, 1], nlag = 0, hess = fit$hessian, n = T, arglist = arglist)
			fit$robust.cvar = tmp$vcv
			fit$scores = numDeriv::jacobian(func = f, x = ipars[estidx, 1], method="Richardson", method.args=list(), arglist = arglist) 
			colnames(fit$scores) = names(ipars[estidx, 1])
			fit$se.coef = sqrt(diag(abs(fit$cvar)))
			fit$tval = fit$coef[-fixed]/fit$se.coef
			# change here
			fit$matcoef = matrix(NA, ncol = 4, nrow = length(fit$coef))
			fit$matcoef[-fixed,] = cbind(fit$coef[-fixed], fit$se.coef,
					fit$tval, 2*(1-pnorm(abs(fit$tval))))
			fit$matcoef[fixed,] = cbind(fit$coef[fixed], fNA,fNA, fNA)
			fit$robust.se.coef = sqrt(diag(fit$robust.cvar))
			fit$robust.tval = fit$coef[-fixed]/fit$robust.se.coef
			# change here
			fit$robust.matcoef = matrix(NA, ncol = 4, nrow = length(fit$coef))
			fit$robust.matcoef[-fixed,] = cbind(fit$coef[-fixed], fit$robust.se.coef,
					fit$robust.tval, 2*(1-pnorm(abs(fit$robust.tval))))
			fit$robust.matcoef[fixed,] = cbind(fit$coef[fixed], fNA, fNA, fNA)
			fit$hessian.message = NULL
		}
		if(model$modelinc[9] == 0){
			vtomega = ipars[idx["omega", 1], 1]
			names(vtomega) = "omega"
			fit$coef = c(fit$coef, vtomega)
			names(fit$coef)[length(fit$coef)] = "omega"
			fit$se.coef = c(fit$se.coef, NA)
			fit$tval = c(fit$tval, NA)
			# change here
			fit$matcoef = rbind(fit$matcoef, c(vtomega, NA, NA, NA))
			
			fit$robust.se.coef = c(fit$robust.se.coef, NA)
			fit$robust.tval = c(fit$robust.tval, NA)
			# change here
			fit$robust.matcoef = rbind(fit$robust.matcoef, c(vtomega, NA, NA, NA))
		}
	} else{
		fit$coef = ipars[estidx, 1]
		if(is.null(fit$cvar)){
			fit$se.coef = rep(NA,length(fit$coef))
			fit$tval = rep(NA,length(fit$coef))
			# change here
			fit$matcoef = cbind(fit$coef, fit$se.coef,
					fit$tval, rep(NA,length(fit$coef)))
			fit$robust.se.coef = rep(NA,length(fit$coef))
			fit$robust.tval = rep(NA,length(fit$coef))
			# change here
			fit$robust.matcoef = cbind(fit$coef, fit$robust.se.coef,fit$robust.tval, 
					rep(NA,length(fit$coef)))
			fit$hessian.message = "failed to invert hessian"
		} else{
			nlag=min(floor(1.2*(T)^(1/3)),(T))
			arglist$returnType = "LHT"
			tmp = rugarch:::robustvcv(fun = f, pars = ipars[estidx,1], nlag = nlag, hess = fit$hessian, n = T, arglist = arglist)
			fit$robust.cvar = tmp$vcv
			fit$scores = numDeriv::jacobian(func = f, x = ipars[estidx, 1], method="Richardson", method.args=list(), arglist = arglist) 
			colnames(fit$scores) = names(ipars[estidx, 1])
			fit$se.coef = sqrt(diag(abs(fit$cvar)))
			fit$tval = fit$coef/fit$se.coef
			# change here
			fit$matcoef = cbind(fit$coef, fit$se.coef,
					fit$tval, 2*(1-pnorm(abs(fit$tval))))
			fit$robust.se.coef = sqrt(diag(fit$robust.cvar))
			fit$robust.tval = fit$coef/fit$robust.se.coef
			# change here
			fit$robust.matcoef = cbind(fit$coef, fit$robust.se.coef,fit$robust.tval, 
					2*(1-pnorm(abs(fit$robust.tval))))
			fit$hessian.message = NULL
		}
		# variance targeting case
		if(model$modelinc[9] == 0){
			vtomega = ipars[idx["omega", 1], 1]
			names(vtomega) = "omega"
			fit$coef = c(fit$coef, vtomega)
			names(fit$coef)[length(fit$coef)] = "omega"
			fit$se.coef = c(fit$se.coef, NA)
			fit$tval = c(fit$tval, NA)
			# change here
			fit$matcoef = rbind(fit$matcoef, c(vtomega, NA, NA, NA))
			
			fit$robust.se.coef = c(fit$robust.se.coef, NA)
			fit$robust.tval = c(fit$robust.tval, NA)
			# change here
			fit$robust.matcoef = rbind(fit$robust.matcoef, c(vtomega, NA, NA, NA))
		}
	}
	# return the correct indicators to the ipars (changed in the case of fixed pars and fixed.se = TRUE)
	ipars[,2:4] = model$pars[,2:4]
	dimnames(fit$matcoef) = list(names(fit$coef), c(" Estimate",
					" Std. Error", " t value", "Pr(>|t|)"))
	dimnames(fit$robust.matcoef) = list(names(fit$coef), c(" Estimate",
					" Std. Error", " t value", "Pr(>|t|)"))
	fit$fitted.values = data-fit$residuals
	fit$convergence = convergence
	fit$message = message
	fit$timer = timer
	fit$ipars = ipars
	return(fit)
}



.acdforcregressors = function(model, mregfor, vregfor, skregfor, shregfor, n.ahead, N, out.sample, n.roll)
{
	# N is the original length
	treq = n.ahead + n.roll
	mxn = model$modelinc[8]
	vxn = model$modelinc[17]
	skxn = model$modelinc[25]
	shxn = model$modelinc[31]
	
	if(mxn>0){
		if(!is.null(mregfor)){
			nmex = NROW(as.matrix(mregfor))
			mmex = NCOL(as.matrix(mregfor))
		} else{
			nmex = 0
			mmex = 0
		}
		if(!is.null(mregfor) && mmex != mxn)
		{
			cat("\nacdforecast-->error: Column dimension of external mean forecast matrix is wrong.")
			cat(paste("\nModel has ", mxn, " external regressors but forecast matrix has ", mmex, sep = ""))
			stop("\n...exiting\n")
		}
		
		if(!is.null(mregfor) && nmex < treq)
		{
			cat(paste("\nacdforecast-->error: You requested ", treq ," actual forecasts (including the rolling periods) but external mean forecasts provided have only ", nmex, " rows",sep=""))
			cat(paste("\nA minimum of ", treq," rows are required", sep=""))
			stop("\n...exiting")
		}
		if(is.null(mregfor)){
			mxf = rbind(as.matrix(model$modeldata$mexdata)[1:(N-out.sample), ,drop = FALSE], matrix(0, ncol = mxn, nrow = treq))
		} else {
			mxf = rbind(as.matrix(model$modeldata$mexdata)[1:(N-out.sample), ,drop = FALSE], as.matrix(mregfor)[1:treq,,drop = FALSE])
		}
	} else{
		mxf = NULL
	}
	if(vxn>0){
		if(!is.null(skregfor)){
			nvex = NROW(as.matrix(vregfor))
			mvex = NCOL(as.matrix(vregfor))
		} else{
			nvex = 0
			mvex = 0
		}
		if(!is.null(vregfor) && mvex != vxn)
		{
			cat("\nacdforecast-->error: Column dimension of external variance forecast matrix is wrong.")
			cat(paste("\nModel has ",vxn," external regressors but forecast matrix has", mvex, sep = ""))
			stop("\n...exiting\n")
		}
		# N is the original length
		if(!is.null(vregfor) && nvex < treq)
		{
			cat(paste("\nacdforecast-->error: You requested ", treq ," actual forecasts (including the rolling periods) but external variance forecasts provided have only ", nvex, " rows",sep=""))
			cat(paste("\nA minimum of ", treq," rows are required", sep=""))
			stop("\n...exiting")
		}
		if(is.null(vregfor)){
			vxf = rbind(as.matrix(model$modeldata$vexdata)[1:(N-out.sample), ,drop = FALSE], matrix(0, ncol = vxn, nrow = treq))
		} else {
			vxf = rbind(as.matrix(model$modeldata$vexdata)[1:(N-out.sample), ,drop = FALSE], as.matrix(vregfor)[1:treq,,drop = FALSE])
		}
	} else{
		vxf = NULL
	}
	
	if(skxn>0){
		if(!is.null(skregfor)){
			nvex = NROW(as.matrix(skregfor))
			mvex = NCOL(as.matrix(skregfor))
		} else{
			nvex = 0
			mvex = 0
		}
		if(!is.null(skregfor) && mvex != skxn)
		{
			cat("\nacdforecast-->error: Column dimension of external skew forecast matrix is wrong.")
			cat(paste("\nModel has ", skxn," external regressors but forecast matrix has", mvex, sep = ""))
			stop("\n...exiting\n")
		}
		# N is the original length
		if(!is.null(skregfor) && nvex < treq)
		{
			cat(paste("\nacdforecast-->error: You requested ", treq ," actual forecasts (including the rolling periods) but external skew forecasts provided have only ", nvex, " rows",sep=""))
			cat(paste("\nA minimum of ", treq," rows are required", sep=""))
			stop("\n...exiting")
		}
		if(is.null(skregfor)){
			skxf = rbind(as.matrix(model$modeldata$skxdata)[1:(N-out.sample), ,drop = FALSE], matrix(0, ncol = skxn, nrow = treq))
		} else {
			skxf = rbind(as.matrix(model$modeldata$skxdata)[1:(N-out.sample), ,drop = FALSE], as.matrix(skregfor)[1:treq,,drop = FALSE])
		}
	} else{
		skxf = NULL
	}
	
	if(shxn>0){
		if(!is.null(shregfor)){
			nvex = NROW(as.matrix(shregfor))
			mvex = NCOL(as.matrix(shregfor))
		} else{
			nvex = 0
			mvex = 0
		}
		if(!is.null(shregfor) && mvex != shxn)
		{
			cat("\nacdforecast-->error: Column dimension of external shape forecast matrix is wrong.")
			cat(paste("\nModel has ", shxn," external regressors but forecast matrix has", mvex, sep = ""))
			stop("\n...exiting\n")
		}
		# N is the original length
		if(!is.null(shregfor) && nvex < treq)
		{
			cat(paste("\nacdforecast-->error: You requested ", treq ," actual forecasts (including the rolling periods) but external shape forecasts provided have only ", nvex, " rows",sep=""))
			cat(paste("\nA minimum of ", treq," rows are required", sep=""))
			stop("\n...exiting")
		}
		if(is.null(shregfor)){
			shxf = rbind(as.matrix(model$modeldata$shxdata)[1:(N-out.sample), ,drop = FALSE], matrix(0, ncol = shxn, nrow = treq))
		} else {
			shxf = rbind(as.matrix(model$modeldata$shxdata)[1:(N-out.sample), ,drop = FALSE], as.matrix(shregfor)[1:treq,,drop = FALSE])
		}
	} else{
		shxf = NULL
	}
	return(list(mxf = mxf, vxf = vxf, skxf = skxf, shxf = shxf))
}

.acdsimregressors = function(model, mexsimdata, vexsimdata, skexsimdata, shexsimdata, N, n, m.sim, m)
{
	mxn = model$modelinc[8]
	vxn = model$modelinc[17]
	skxn = model$modelinc[25]
	shxn = model$modelinc[31]
	
	if(mxn>0){
		if(is.null(mexsimdata)){
			mexsimdata = vector(mode = "list", length = m.sim)
			for(i in 1:m.sim) mexsimdata[[i]] = matrix(0, ncol = mxn, nrow = n)
		}
		if(!is.null(mexsimdata))
		{
			if(!is.list(mexsimdata)) stop("\nacdhsim-->error: mexsimdata should be a list of length m.sim")
			if(length(mexsimdata) != m.sim){
				msd = vector(mode = "list", length = m.sim)
				for(i in 1:m.sim) msd[[i]] = as.matrix(mexsimdata[[1]])
				mexsimdata = msd
				warning("\nacdhsim-->warning: length of mexsimdata list not equal to m.sim...\nreplicating first list element m.sim times.\n")
			}
			for(i in 1:m.sim){
				if(dim(as.matrix(mexsimdata[[i]]))[2] != mxn ) 
					stop(paste("\nacdhsim-->error: mexsimdata ", i," has wrong no. of column", sep=""))
				if(dim(as.matrix(mexsimdata[[i]]))[1] != n )
					stop(paste("\nacdhsim-->error: mexsimdata ", i," has wrong no. of rows", sep=""))
			}		
		}
		mexsimlist = vector(mode = "list", length = m.sim)
		for(i in 1:m.sim){
			premexdata = model$modeldata$mexdata[(N-m+1):N,,drop=FALSE]
			mexsimlist[[i]] = matrix(rbind(premexdata, mexsimdata[[i]]), ncol = mxn)
		}
	} else{
		mexsimlist = vector(mode = "list", length = m.sim)
		for(i in 1:m.sim) mexsimlist[[i]]=0
	}
	if(vxn>0){
		if(is.null(vexsimdata)){
			vexsimdata = vector(mode = "list", length = m.sim)
			for(i in 1:m.sim) vexsimdata[[i]] = matrix(0, ncol = vxn, nrow = n)
		}
		if(!is.null(vexsimdata))
		{
			if(!is.list(vexsimdata)) 
				stop("\nacdsim-->error: vexsimdata should be a list of length m.sim")
			if(length(vexsimdata) != m.sim){
				vsd = vector(mode = "list", length = m.sim)
				for(i in 1:m.sim) vsd[[i]] = as.matrix(vexsimdata[[1]])
				vexsimdata = vsd
				warning("\nacdhsim-->warning: length of vexsimdata list not equal to m.sim...\nreplicating first list element m.sim times.\n")
			}
			for(i in 1:m.sim){
				if(dim(as.matrix(vexsimdata[[i]]))[2] != vxn ) 
					stop(paste("\nacdhsim-->error: vexsimdata ", i," has wrong no. of column", sep=""))
				if(dim(as.matrix(vexsimdata[[i]]))[1] != n )
					stop(paste("\nacdhsim-->error: vexsimdata ", i," has wrong no. of rows", sep=""))
			}		
		}
		vexsimlist = vector(mode = "list", length = m.sim)
		for(i in 1:m.sim){
			prevexdata = model$modeldata$vexdata[(N-m+1):N,,drop=FALSE]
			vexsimlist[[i]] = matrix(rbind(prevexdata, vexsimdata[[i]]), ncol = vxn)
		}
	} else{
		vexsimlist = vector(mode = "list", length = m.sim)
		for(i in 1:m.sim) vexsimlist[[i]]=0
	}
	
	if(skxn>0){
		if(is.null(skexsimdata)){
			skexsimdata = vector(mode = "list", length = m.sim)
			for(i in 1:m.sim) skexsimdata[[i]] = matrix(0, ncol = skxn, nrow = n)
		}
		if(!is.null(skexsimdata))
		{
			if(!is.list(skexsimdata)) 
				stop("\nacdsim-->error: skexsimdata should be a list of length m.sim")
			if(length(skexsimdata) != m.sim){
				sksd = vector(mode = "list", length = m.sim)
				for(i in 1:m.sim) sksd[[i]] = as.matrix(skexsimdata[[1]])
				skexsimdata = sksd
				warning("\nacdhsim-->warning: length of skexsimdata list not equal to m.sim...\nreplicating first list element m.sim times.\n")
			}
			for(i in 1:m.sim){
				if(dim(as.matrix(skexsimdata[[i]]))[2] != skxn ) 
					stop(paste("\nacdhsim-->error: skexsimdata ", i," has wrong no. of column", sep=""))
				if(dim(as.matrix(skexsimdata[[i]]))[1] != n )
					stop(paste("\nacdhsim-->error: skexsimdata ", i," has wrong no. of rows", sep=""))
			}
		}
		skexsimlist = vector(mode = "list", length = m.sim)
		for(i in 1:m.sim){
			preskexdata = model$modeldata$skxdata[(N-m+1):N,,drop=FALSE]
			skexsimlist[[i]] = matrix(rbind(preskexdata, skexsimdata[[i]]), ncol = skxn)
		}
	} else{
		skexsimlist = vector(mode = "list", length = m.sim)
		for(i in 1:m.sim) skexsimlist[[i]]=0
	}
	
	if(shxn>0){
		if(is.null(shexsimdata)){
			shexsimdata = vector(mode = "list", length = m.sim)
			for(i in 1:m.sim) shexsimdata[[i]] = matrix(0, ncol = shxn, nrow = n)
		}
		if(!is.null(shexsimdata))
		{
			if(!is.list(shexsimdata)) 
				stop("\nacdsim-->error: skexsimdata should be a list of length m.sim")
			if(length(shexsimdata) != m.sim){
				shsd = vector(mode = "list", length = m.sim)
				for(i in 1:m.sim) shsd[[i]] = as.matrix(shexsimdata[[1]])
				shexsimdata = shsd
				warning("\nacdhsim-->warning: length of shexsimdata list not equal to m.sim...\nreplicating first list element m.sim times.\n")
			}
			for(i in 1:m.sim){
				if(dim(as.matrix(shexsimdata[[i]]))[2] != shxn ) 
					stop(paste("\nacdhsim-->error: shexsimdata ", i," has wrong no. of column", sep=""))
				if(dim(as.matrix(shexsimdata[[i]]))[1] != n )
					stop(paste("\nacdhsim-->error: shexsimdata ", i," has wrong no. of rows", sep=""))
			}
		}
		shexsimlist = vector(mode = "list", length = m.sim)
		for(i in 1:m.sim){
			preshexdata = model$modeldata$shxdata[(N-m+1):N,,drop=FALSE]
			shexsimlist[[i]] = matrix(rbind(preshexdata, shexsimdata[[i]]), ncol = shxn)
		}
	} else{
		shexsimlist = vector(mode = "list", length = m.sim)
		for(i in 1:m.sim) shexsimlist[[i]]=0
	}
	
	
	return(list(mexsimlist = mexsimlist, vexsimlist = vexsimlist, 
					skexsimlist = skexsimlist, shexsimlist = shexsimlist))
}
