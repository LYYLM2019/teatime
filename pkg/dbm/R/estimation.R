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

# y: xts (multivariate) object (no NAs)
# x.vars: column names of the x variables in the y object
# x.lags: (length of x.vars) lags per  x variable
# --> this allows to include the same variable with different lag and have non-
# consecutive lags
# arp: autoregressive order of probability
# arq: autoregressive order of y Indicator variable (regressand)
# ecm: impose restrictions on value of arq coefficient (1-value(arp coefficient))
# constant: intercept
# link: gaussian (probit) or logistic (logit)
# fixed.pars: optional named vector of fixed parameters (if all are fixed then
# we filter the data without estimation)
# start.pars: optional starting vector of parameters
# solver: either optim or gosolnp
# control: control parameters passed to the solver
# ...: not currently used
# idx [omega, length(x.vars), arp, arq, use_ecm, link, regularize]
dbm = function(y, x.vars = NULL, x.lags = 1, arp = 1, arq = 0, ecm = FALSE, 
		constant = TRUE, link = "gaussian", fixed.pars = NULL, 
		solver = "optim", control=list(), parsearch = TRUE, parsim = 5000, 
		method = "Nelder-Mead", regularization = FALSE, reg.cost = 1, ...)
{
	call <- match.call()
	if(!is(y, "xts")) stop("\ny must be an xts object")
	if(any(is.na(y))){
		stop("\nNA's found in y. Remove and resubmit.")
	}
	#x.vars = colnames(x)
	
	modelnames = c(if(constant) "omega" else NULL, if(!is.null(x.vars)) paste("beta[",1:length(x.vars),"]",sep="") else NULL,
			if(arp>0) paste("alpha[",1:arp,"]",sep="") else NULL, if(arq>0 && !ecm) paste("delta[",1:arq,"]",sep="") else NULL,
			if(link=="glogistic") "skew[k]" else NULL)
	idx = rep(0, 7)
	if(constant){
		idx[1] = 1
		omega  =  0
		pnames = "omega"
	} else{
		pnames = NULL
		omega  = NULL
	}
	model = list()
	model$y = y
	model$x.vars = x.vars
	if(length(x.lags)<length(x.vars)) x.lags = rep(x.lags, length(x.vars)) else x.lags = x.lags[1:length(x.vars)]
	model$x.lags = x.lags
	model$arp = arp
	model$arq = arq
	model$link = link
	model$ecm = ecm
	model$constant = constant
	model$regularization = regularization
	model$reg.cost = reg.cost
	if(!is.null(x.vars)){
		len_x = length(x.vars)
		tmp = match(x.vars, colnames(y))
		if(any(is.na(tmp))) stop("\nx.vars names not found in colnames(y)")
		x = coredata(y[,tmp])
		yname = colnames(y)[-tmp][1]
		y = coredata(y[,-tmp][,1])
		beta = rep(0, len_x)
		pnames = c(pnames, paste("beta[", 1:len_x, "]",sep=""))
		xidx = x.lags
		idx[2] = len_x
	} else{
		x = 0
		beta = NULL
		len_x  = 0
		xidx = 0
		y = coredata(y)[,1]
		yname = colnames(y)[1]
	}
	model$yname = yname
	# initialize when arp is allowed to be > 1
	# if(ecm && arp>1) stop("\necm only supported for ARp=1")
	if(as.integer(arp)>0){
		if(arp>1) stop("\nonly arp=1 currently supported")
		alpha = rep(0, arp)
		pnames = c(pnames, paste("alpha[", 1:arp, "]",sep=""))
		idx[3] = arp
	} else{
		arp = 0
		alpha = NULL
	}
	if(!ecm){
		if(as.integer(arq)>0){
			if(arq>1) stop("\nonly arq=1 currently supported")
			delta = rep(0, arq)
			pnames = c(pnames, paste("delta[", 1:arq, "]",sep=""))
			idx[4] = arq
		} else{
			arq = 0
			delta = NULL
		}
	} else{
		idx[5] = 1
		arq  = 1
		delta = NULL
	}
	if(link=="glogistic"){
		k = 1
		pnames = c(pnames, "skew[k]")
		idx[6] = as.integer(3)
	} else{
		k = NULL
		idx[6] = as.integer(ifelse(link=="gaussian", 1, 2))
	}
	if(regularization) idx[7] = 1
	n = as.integer(NROW(y))
	pmu = double(n)
	
	if( length(pnames)==0 ) stop("\nNo coefficients to estimate!")
	
	pars = as.double(c(omega, alpha, beta, delta, k))
	
	if(!is.null(fixed.pars)){
		fnames = names(fixed.pars)
		chk = match(fnames, pnames)
		if(any(is.na(chk))){
			cx = which(is.na(chk))
			cat("\nunidentified fixed parameter names in fixed.pars\n")
			cat("\nexpected:",pnames)
			cat("\ngot:",fnames)
			cat("\n")
			stop()
		} else{
			fpars = fixed.pars
			fidx = chk
			# if all parameters fixed return filtered object
			if(length(fpars)==length(pars)){
				model$idx = idx
				model$xidx = xidx
				sol = .filter.dbm(pars = fpars, model, call = call)
				return(sol)
			}
			pars = pars[-fidx]
			pnames = pnames[-fidx]
		}
	} else{
		fpars = NULL
		fnames = NULL
		fidx = NULL
	}
	
	arglist = list()
	arglist$idx = idx
	arglist$xidx = xidx
	arglist$y = y
	arglist$x = x
	arglist$pnames = pnames
	arglist$fnames = fnames
	arglist$modelnames = modelnames
	arglist$fpars = fpars
	arglist$fidx = fidx
	arglist$type = "llh"
	arglist$derivarg = 1
	arglist$regularization = regularization
	if(idx[7]==1) arglist$Cost = reg.cost else arglist$Cost = 0
	dbmenv = new.env(hash = TRUE)
	assign("dbm_llh", 1e4, envir = dbmenv)
	arglist$dbmenv <- dbmenv
	if(parsearch){
		arglist$transform = FALSE
		LB = rep(-25, length(pars))
		UB = rep( 25, length(pars))
		if(solver!="optim"){
			if(idx[3]>0) LB[which(substr(pnames, 1,5)=="alpha")]=-0.99999
			if(idx[3]>0) UB[which(substr(pnames, 1,5)=="alpha")]= 0.99999
			if(idx[6]==3) LB[which(substr(pnames, 1,4)=="skew")] = 0.1
			if(idx[6]==3) UB[which(substr(pnames, 1,4)=="skew")] = 100
		} else{
			arglist$transform = TRUE
		}
		dlist = vector(mode="list", length=length(pars))
		spars = startpars(pars  = pars, fun = dbmlik, LB = LB, UB = UB, arglist = arglist,
				n.sim = parsim, bestN = 1, distr.opt = dlist)
		spars = as.numeric(spars)
		spars = spars[-length(spars)]
	} else{
		spars = runif(length(pars))
	}
	if(ecm){
		gr = NULL
	} else{
		gr = switch(as.character(idx[6]),
				"1" = dbmderiv2,
				"2" = dbmderiv1,
				"3" = dbmderiv3)
	}
	if(solver=="optim"){
		arglist$transform = TRUE
		if(is.null(control$maxit)) control$maxit=50000
		if(idx[6]==1){
			sol = optim(par = spars, fn = dbmlik, gr = gr, arglist = arglist,
					control = control, hessian = TRUE, method = method)
		} else if(idx[6]==2){
			sol = optim(par = spars, fn = dbmlik, gr = gr, arglist = arglist,
					control = control, hessian = TRUE, method = method)
		} else{
			sol = optim(par = spars, fn = dbmlik, gr = gr, arglist = arglist,
					control = control, hessian = TRUE, method = method)
		}
		pars = sol$par
	} else if(solver == "gosolnp"){
		arglist$transform = FALSE
		LB = -25*abs(spars)
		UB =  25*abs(spars)
		if(idx[3]>0) LB[which(substr(pnames, 1,5)=="alpha")]=-0.99999
		if(idx[3]>0) UB[which(substr(pnames, 1,5)=="alpha")]= 0.99999
		if(idx[6]==3) LB[which(substr(pnames, 1,4)=="skew")] = 0.1
		if(idx[6]==3) UB[which(substr(pnames, 1,4)=="skew")] = 100
		sol = gosolnp(pars = spars, fun = dbmlik, arglist = arglist, control = control, 
				LB = LB, UB = UB, ...)
		pars = sol$pars
	} else{
		arglist$transform = FALSE
		LB = -25*abs(spars)
		UB = 25*abs(spars)
		if(idx[3]>0) LB[which(substr(pnames, 1,5)=="alpha")]=-0.99999
		if(idx[3]>0) UB[which(substr(pnames, 1,5)=="alpha")]= 0.99999
		if(idx[6]==3) LB[which(substr(pnames, 1,4)=="skew")] = 0.1
		if(idx[6]==3) UB[which(substr(pnames, 1,4)=="skew")] = 100
		if(idx[6]==1){
			sol = nloptr(x0 = spars,  eval_f = dbmlik, 
					eval_grad_f = gr, eval_g_ineq = NULL, 
					eval_jac_g_ineq = NULL, eval_g_eq = NULL, 
    				eval_jac_g_eq = NULL, opts = control,
					arglist = arglist, lb = LB, ub = UB)
		} else if(idx[6]==2){
			sol = nloptr(x0 = spars,  eval_f = dbmlik, 
					eval_grad_f = gr, eval_g_ineq = NULL, 
					eval_jac_g_ineq = NULL, eval_g_eq = NULL, 
    				eval_jac_g_eq = NULL, opts = control,
					arglist = arglist, lb = LB, ub = UB)
		} else{
			sol = nloptr(x0 = spars,  eval_f = dbmlik, 
					eval_grad_f = gr, eval_g_ineq = NULL, 
					eval_jac_g_ineq = NULL, eval_g_eq = NULL, 
    				eval_jac_g_eq = NULL, opts = control,
					arglist = arglist, lb = LB, ub = UB)
		}
		pars = sol$solution
	}
	names(pars)<-pnames
	if(idx[3]>0){
		if(arglist$transform && !any(fnames=="alpha[1]")) pars[paste("alpha[",1:idx[3],"]",sep="")] = logtransform(pars[paste("alpha[",1:idx[3], "]",sep="")],-0.99999,0.99999)
	}
	if(arglist$transform && idx[6]==3 && !any(fnames=="skew[k]")){
		pars["skew[k]"] = logtransform(pars["skew[k]"], 0.01, 100)
	}
	arglist$transform = FALSE
	fit = .postestimate(f = dbmlik, pars, arglist)
	model$parnames = pnames
	model$modelnames = modelnames
	model$fixed.pars = fpars
	model$fidx = fidx
	model$idx = idx
	model$xidx = xidx
	model$transform = arglist$transform
	if(idx[6]==1){
		gfit<-eval(parse(text=paste("glm(",yname,"~1, family=binomial(probit), data=as.data.frame(y))",sep="")))
		null.deviance = deviance(gfit)
		null.epcp = epcp.default(fitted(gfit), y = as.numeric(arglist$y))
	} else if(idx[6]==2){
		gfit<-eval(parse(text=paste("glm(",yname,"~1, family=binomial(logit), data=as.data.frame(y))",sep="")))
		null.deviance = deviance(gfit)
		null.epcp = epcp.default(fitted(gfit), y = as.numeric(arglist$y))
	} else{
		# need to create a function for this
		if(!is.null(x.vars) | arp>0 | arq >0){
			newmod = dbm(model$y[,yname], x.vars = NULL, x.lags = 0, arp = 0, arq = 0, 
					ecm = FALSE, constant = TRUE, link = "glogistic", fixed.pars = NULL, 
					solver = "optim", control=list(maxit=50000,trace=0), 
					parsearch = FALSE, parsim = 5000, method = "Nelder-Mead")
			null.deviance = deviance(newmod)
			null.epcp = epcp.dbm(newmod)
		} else {
			null.deviance = NA
			null.epcp = NA
		}
	}
	fit$null.epcp = null.epcp
	fit$null.deviance = null.deviance
	fit$model.epcp = epcp.default(fit$fitted.values, as.numeric(arglist$y))
	model$estimation = "maxlik"
	model$solver = solver
	if(solver=="optim") model$method = method
	out = list(model = model, fit = fit, call = call)
	class(out) <- "dbm"
	return(out)
}
# idx [omega, length(x.vars), arp, arq, use_ecm, link]

dbmlik = function(pars, arglist)
{
	# extract parameters and prepare vectors
	idx = arglist$idx
	dbmenv<-arglist$dbmenv
	names(pars)<-arglist$pnames
	Cost = 0
	if(!is.null(arglist$fidx)){
		pars = c(pars, arglist$fpars)
	}
	if(idx[1]>0){
		omega = pars["omega"]
		Cost = Cost+omega^2
	} else{
		omega = 0
	}
	if(idx[3]>0){
		alpha = pars[paste("alpha[",1:idx[3], "]",sep="")]
		if(arglist$transform && !any(arglist$fnames=="alpha[1]")) alpha = logtransform(alpha, -0.99999, 0.99999)
		Cost = Cost+sum(alpha^2)
	} else{
		alpha = 0
	}
	if(idx[5]>0){
		delta = 0
	} else{
		if(idx[4]>0) delta = pars[paste("delta[",1:idx[4], "]",sep="")] else delta = 0
		Cost = Cost+sum(delta^2)
	}
	if(idx[2]>0) beta = pars[paste("beta[",1:idx[2], "]",sep="")] else beta = 0
	Cost = Cost+sum(beta^2)
	if(idx[6]==3){
		k = pars["skew[k]"]
		if(arglist$transform && !any(arglist$fnames=="skew[k]")) k = logtransform(k, 0.01, 100) else k = max(0.01, k)
		Cost = Cost+sum(k^2)
	} else{
		k = 1
	}
	x = arglist$x
	y = arglist$y
	n = length(y)
	mpu = double(n)
	xidx = arglist$xidx
	# recursion initialization
	rcs = omega + mean(y)*delta[1]
	if(idx[2]>0) rcs = rcs + sum(colMeans(x)*beta)
	if(idx[3]>0) rcs = rcs/(1-alpha[1])
	# Likelihood Evaluation C
	Cost = (0.5*arglist$Cost * Cost)/n
	tmp = try(.C("c_dbmestimate", y = as.double(y), x = as.double(as.vector(x)), 
					mpu = as.double(mpu), mpuinit = as.double(rcs),
					omega = as.double(omega), alpha = as.double(alpha), 
					delta = as.double(delta), beta = as.double(beta), 
					k = as.double(k), Cost = as.double(Cost),
					lik = double(n), llh = double(1), idx = as.integer(idx), 
					xidx = as.integer(xidx), T = as.integer(n), 
					PACKAGE="dbm"), silent = TRUE)
	# Check and return
	if(is.na(tmp$llh) | is.nan(tmp$llh) | !is.finite(tmp$llh)){
		return(llh = get("dbm_llh", dbmenv) + 0.25*(abs(get("dbm_llh", dbmenv))))
	} else{
		assign("dbm_llh", tmp$llh, envir = dbmenv)
	}
	ans = switch(arglist$type,
			llh = tmp$llh,
			lik = tmp$lik,
			ALL = list(llh = tmp$llh, lik = tmp$lik, mpu = tmp$mpu))
	return(ans)
}

dbmderiv1 = function(pars, arglist)
{
	x = arglist$x
	y = arglist$y
	n = length(y)
	idx = arglist$idx
	dbmenv<-arglist$dbmenv
	names(pars)<-arglist$pnames
	if(!is.null(arglist$fidx)){
		pars = c(pars, arglist$fpars)
	}
	if(idx[1]>0){
		omega = pars["omega"]
    } else{
		omega = 0
		
	}
	dpomega = double(n)
	dvomega = double(n)
	domega = 0
	
	if(idx[2]>0){
		beta = pars[paste("beta[",1:idx[2], "]",sep="")]
		dvbeta = double(n*idx[2])
		dpbeta = double(n*idx[2])
	} else{
		beta = 0
		dvbeta = double(n)
		dpbeta = double(n)
	}
	dbeta = rep(0, length(beta))
	
	if(idx[3]>0){
		alpha = pars[paste("alpha[",1:idx[3], "]",sep="")]
		if(arglist$transform && !any(arglist$fnames=="alpha[1]")) alpha = logtransform(alpha, -0.99999, 0.99999)
		dalpha = rep(0, length(alpha))
		dvalpha = double(n*idx[3])
		dpalpha = double(n*idx[3])
	} else{
		alpha = 0
		dalpha = 0
		dvalpha = double(n)
		dpalpha = double(n)
	}
	if(idx[5]>0){
		delta = 0
		ddelta = 0
		dvdelta = double(n)
		dpdelta = double(n)
	} else{
		if(idx[4]>0){
			delta = pars[paste("delta[",1:idx[4], "]",sep="")]
			ddelta = rep(0, length(delta))
			dvdelta = double(n*idx[4])
			dpdelta = double(n*idx[4])
		} else{
			delta  = 0
			ddelta = 0
			dvdelta = double(n)
			dpdelta = double(n)
		}
	}
	mpu = double(n)
	xidx = arglist$xidx
	# recursion initialization
	if(idx[2]>0){
		meanx = colMeans(x)
	} else{
		meanx = 0
	}
	meany = mean(y)
	Cost = arglist$Cost/n
	
	tmp = try(.C("c_dbmderiv1", y = as.double(y), x = as.double(as.vector(x)), 
					mpu = as.double(mpu), meanx = as.double(meanx), meany = as.double(meany),
					omega = as.double(omega), alpha = as.double(alpha), 
					delta = as.double(delta), beta = as.double(beta), 
					domega = as.double(domega), dalpha = as.double(dalpha),
					ddelta = as.double(ddelta), dbeta = as.double(dbeta),
					dpomega = as.double(dpomega), dpalpha = as.double(dpalpha), 
					dpdelta = as.double(dpdelta), dpbeta = as.double(dpbeta), 
					dvomega = as.double(dvomega), dvalpha = as.double(dvalpha), 
					dvdelta = as.double(dvdelta), dvbeta = as.double(dvbeta), 
					Cost = as.double(Cost), idx = as.integer(idx), xidx = as.integer(xidx), 
					T = as.integer(n), 
					PACKAGE="dbm"), silent = TRUE)
	# Check and return
	if(arglist$derivarg == 1){
		if(inherits(tmp, 'try-error')){
			ans = rep(1e5, length(arglist$pnames))
		} else{
			# NOTE: multiplication by -1 (minimization)
			ans = -1*c(if(idx[1]>0) tmp$domega else NULL, if(idx[2]>0) tmp$dbeta else NULL, 
					if(idx[3]>0) tmp$dalpha else NULL, if(idx[5]==0 && idx[4]>0) tmp$ddelta else NULL)
			names(ans)<-arglist$modelnames
			ans = ans[arglist$pnames]
		}
	} else{
		if(!inherits(tmp, 'try-error')){
			# NOTE: NO multiplication by -1 (this returns the scores of the 
			# maximized likelihood)
			ans = matrix(NA, ncol = 1, nrow = n)
			if(idx[1]>0) ans  = cbind(ans, tmp$dvomega)
			if(idx[2]>0) ans = cbind(ans, matrix(tmp$dvbeta, ncol = idx[2], nrow=n))
			if(idx[3]>0) ans = cbind(ans, matrix(tmp$dvalpha, ncol = idx[3], nrow=n))
			if(idx[5]==0 && idx[4]>0) ans = cbind(ans, matrix(tmp$dvdelta, ncol = idx[4], nrow=n))
			ans = ans[,-1]
			colnames(ans)<-arglist$modelnames
			ans = ans[,arglist$pnames]
		} else{
			ans = tmp
		}
	}
	return(ans)
}

dbmderiv2 = function(pars, arglist)
{
	x = arglist$x
	y = arglist$y
	n = length(y)
	idx = arglist$idx
	dbmenv<-arglist$dbmenv
	names(pars)<-arglist$pnames
	if(!is.null(arglist$fidx)){
		pars = c(pars, arglist$fpars)
	}
	if(idx[1]>0){
		omega = pars["omega"]
    } else{
		omega = 0
		
	}
	dpomega = double(n)
	dvomega = double(n)
	domega = 0
	
	if(idx[2]>0){
		beta = pars[paste("beta[",1:idx[2], "]",sep="")]
		dvbeta = double(n*idx[2])
		dpbeta = double(n*idx[2])
	} else{
		beta = 0
		dvbeta = double(n)
		dpbeta = double(n)
	}
	dbeta = rep(0, length(beta))
	
	if(idx[3]>0){
		alpha = pars[paste("alpha[",1:idx[3], "]",sep="")]
		if(arglist$transform && !any(arglist$fnames=="alpha[1]")) alpha = logtransform(alpha, -0.99999, 0.99999)
		dalpha = rep(0, length(alpha))
		dvalpha = double(n*idx[3])
		dpalpha = double(n*idx[3])
	} else{
		alpha = 0
		dalpha = 0
		dvalpha = double(n)
		dpalpha = double(n)
	}
	if(idx[5]>0){
		delta = 0
		ddelta = 0
		dvdelta = double(n)
		dpdelta = double(n)
	} else{
		if(idx[4]>0){
			delta = pars[paste("delta[",1:idx[4], "]",sep="")]
			ddelta = rep(0, length(delta))
			dvdelta = double(n*idx[4])
			dpdelta = double(n*idx[4])
		} else{
			delta  = 0
			ddelta = 0
			dvdelta = double(n)
			dpdelta = double(n)
		}
	}
	mpu = double(n)
	xidx = arglist$xidx
	# recursion initialization
	if(idx[2]>0){
		meanx = colMeans(x)
	} else{
		meanx = 0
	}
	meany = mean(y)
	Cost = arglist$Cost/n
	
	tmp = try(.C("c_dbmderiv2", y = as.double(y), x = as.double(as.vector(x)), 
					mpu = as.double(mpu), meanx = as.double(meanx), meany = as.double(meany),
					omega = as.double(omega), alpha = as.double(alpha), 
					delta = as.double(delta), beta = as.double(beta), 
					domega = as.double(domega), dalpha = as.double(dalpha),
					ddelta = as.double(ddelta), dbeta = as.double(dbeta),
					dpomega = as.double(dpomega), dpalpha = as.double(dpalpha), 
					dpdelta = as.double(dpdelta), dpbeta = as.double(dpbeta), 
					dvomega = as.double(dvomega), dvalpha = as.double(dvalpha), 
					dvdelta = as.double(dvdelta), dvbeta = as.double(dvbeta), 
					Cost = as.double(Cost), idx = as.integer(idx), xidx = as.integer(xidx), 
					T = as.integer(n), 
					PACKAGE="dbm"), silent = TRUE)
	# Check and return
	if(arglist$derivarg == 1){
		if(inherits(tmp, 'try-error')){
			ans = rep(1e5, length(arglist$pnames))
		} else{
			ans = -1*c(if(idx[1]>0) tmp$domega else NULL, if(idx[2]>0) tmp$dbeta else NULL, 
					if(idx[3]>0) tmp$dalpha else NULL, if(idx[5]==0 && idx[4]>0) tmp$ddelta else NULL)
			names(ans)<-arglist$modelnames
			ans = ans[arglist$pnames]
		}
	} else{
		if(!inherits(tmp, 'try-error')){
			ans = matrix(NA, ncol = 1, nrow = n)
			if(idx[1]>0) ans  = cbind(ans, tmp$dvomega)
			if(idx[2]>0) ans = cbind(ans, matrix(tmp$dvbeta, ncol = idx[2], nrow=n))
			if(idx[3]>0) ans = cbind(ans, matrix(tmp$dvalpha, ncol = idx[3], nrow=n))
			if(idx[5]==0 && idx[4]>0) ans = cbind(ans, matrix(tmp$dvdelta, ncol = idx[4], nrow=n))
			ans = ans[,-1]
			colnames(ans)<-arglist$modelnames
			ans = ans[,arglist$pnames]
		} else{
			ans = tmp
		}
	}
	return(ans)
}

dbmderiv3 = function(pars, arglist)
{
	x = arglist$x
	y = arglist$y
	n = length(y)
	idx = arglist$idx
	dbmenv<-arglist$dbmenv
	names(pars)<-arglist$pnames
	if(!is.null(arglist$fidx)){
		pars = c(pars, arglist$fpars)
	}
	if(idx[1]>0){
		omega = pars["omega"]
    } else{
		omega = 0
		
	}
	dpomega = double(n)
	dvomega = double(n)
	domega = 0
	
	if(idx[2]>0){
		beta = pars[paste("beta[",1:idx[2], "]",sep="")]
		dvbeta = double(n*idx[2])
		dpbeta = double(n*idx[2])
	} else{
		beta = 0
		dvbeta = double(n)
		dpbeta = double(n)
	}
	dbeta = rep(0, length(beta))
	
	if(idx[3]>0){
		alpha = pars[paste("alpha[",1:idx[3], "]",sep="")]
		if(arglist$transform && !any(arglist$fnames=="alpha[1]")) alpha = logtransform(alpha, -0.99999, 0.99999)
		dalpha = rep(0, length(alpha))
		dvalpha = double(n*idx[3])
		dpalpha = double(n*idx[3])
	} else{
		alpha = 0
		dalpha = 0
		dvalpha = double(n)
		dpalpha = double(n)
	}
	if(idx[5]>0){
		delta = 0
		ddelta = 0
		dvdelta = double(n)
		dpdelta = double(n)
	} else{
		if(idx[4]>0){
			delta = pars[paste("delta[",1:idx[4], "]",sep="")]
			ddelta = rep(0, length(delta))
			dvdelta = double(n*idx[4])
			dpdelta = double(n*idx[4])
		} else{
			delta  = 0
			ddelta = 0
			dvdelta = double(n)
			dpdelta = double(n)
		}
	}
	if(arglist$transform && !any(arglist$fnames=="skew[k]")) k = logtransform(pars["skew[k]"], 0.01, 100) else k = max(0.01, pars["skew[k]"])
	dk = double(1)
	dvk = double(n)
	dpk = double(n)
	mpu = double(n)
	xidx = arglist$xidx
	# recursion initialization
	if(idx[2]>0){
		meanx = colMeans(x)
	} else{
		meanx = 0
	}
	meany = mean(y)
	Cost = arglist$Cost/n
	
	tmp = try(.C("c_dbmderiv3", y = as.double(y), x = as.double(as.vector(x)), 
					mpu = as.double(mpu), meanx = as.double(meanx), meany = as.double(meany),
					omega = as.double(omega), alpha = as.double(alpha), 
					delta = as.double(delta), beta = as.double(beta), k = as.double(k),
					domega = as.double(domega), dalpha = as.double(dalpha),
					ddelta = as.double(ddelta), dbeta = as.double(dbeta), dk = as.double(dk),
					dpomega = as.double(dpomega), dpalpha = as.double(dpalpha), 
					dpdelta = as.double(dpdelta), dpbeta = as.double(dpbeta), dpk = as.double(dpk),
					dvomega = as.double(dvomega), dvalpha = as.double(dvalpha), 
					dvdelta = as.double(dvdelta), dvbeta = as.double(dvbeta), dvk = as.double(dvk),
					Cost = as.double(Cost), idx = as.integer(idx), xidx = as.integer(xidx), T = as.integer(n), 
					PACKAGE="dbm"), silent = TRUE)
	# Check and return
	if(arglist$derivarg == 1){
		if(inherits(tmp, 'try-error')){
			ans = rep(1e5, length(arglist$pnames))
		} else{
			ans = -1*c(if(idx[1]>0) tmp$domega else NULL, if(idx[2]>0) tmp$dbeta else NULL, 
					if(idx[3]>0) tmp$dalpha else NULL, if(idx[5]==0 && idx[4]>0) tmp$ddelta else NULL,
					tmp$dk)
			names(ans)<-arglist$modelnames
			ans = ans[arglist$pnames]
		}
	} else{
		if(!inherits(tmp, 'try-error')){
			ans = matrix(NA, ncol = 1, nrow = n)
			if(idx[1]>0) ans  = cbind(ans, tmp$dvomega)
			if(idx[2]>0) ans = cbind(ans, matrix(tmp$dvbeta, ncol = idx[2], nrow=n))
			if(idx[3]>0) ans = cbind(ans, matrix(tmp$dvalpha, ncol = idx[3], nrow=n))
			if(idx[5]==0 && idx[4]>0) ans = cbind(ans, matrix(tmp$dvdelta, ncol = idx[4], nrow=n))
			ans = cbind(ans, tmp$dvk)
			ans = ans[,-1]
			colnames(ans)<-arglist$modelnames
			ans = ans[,arglist$pnames]
		} else{
			ans = tmp
		}
	}
	return(ans)
}


# filter data with a fixed set of parameters (despatched from dbm)
.filter.dbm = function(pars, model, call)
{
	if(model$constant){
		omega  =  0
		pnames = "omega"
	} else{
		pnames = NULL
		omega  = NULL
	}
	if(!is.null(model$x.vars)){
		len_x = length(model$x.vars)
		tmp = match(model$x.vars, colnames(model$y))
		if(any(is.na(tmp))) stop("\nx.vars names not found in colnames(y)")
		x = coredata(model$y[,tmp])
		yname = colnames(model$y)[-tmp]
		y = coredata(model$y[,-tmp])
		beta = rep(0, len_x)
		pnames = c(pnames, paste("beta[", 1:len_x, "]",sep=""))
	} else{
		x = 0
		beta = NULL
		len_x  = 0
		y = coredata(model$y)[,1]
		yname = colnames(model$y)[1]
	}
	if(as.integer(model$arp)>0){
		alpha = rep(0, model$arp)
		pnames = c(pnames, paste("alpha[", 1:model$arp, "]",sep=""))
	} else{
		arp = 0
		alpha = NULL
	}
	if(!model$ecm){
		if(as.integer(model$arq)>0){
			delta = rep(0, model$arq)
			pnames = c(pnames, paste("delta[", 1:model$arq, "]",sep=""))
		} else{
			arq = 0
			delta = NULL
		}
	} else{
		delta = NULL
	}
	if(model$link=="glogistic"){
		k = 1
		pnames = c(pnames, "skew[k]")
	} else{
		k = NULL
	}
	
	n = as.integer(NROW(y))
	pmu = double(n)
	if( length(pnames)==0 ) stop("\nmodel has not coefficients!")
	
	fnames = names(pars)
	chk = match(fnames, pnames)
	if(any(is.na(chk))){
		cx = which(is.na(chk))
		cat("\nunidentified parameter names in pars\n")
		cat("\nexpected:",pnames)
		cat("\ngot:",fnames)
		cat("\n")
		stop()
	}
	idx = model$idx
	xidx = model$xidx
	if(idx[1]>0) omega = pars["omega"] else omega = 0
	if(idx[3]>0){
		alpha = pars[paste("alpha[",1:idx[3], "]",sep="")]
	} else{
		alpha = 0
	}
	if(idx[5]>0){
		delta = 0
	} else{
		if(idx[4]>0) delta = pars[paste("delta[",1:idx[4], "]",sep="")] else delta = 0
	}
	if(idx[2]>0) beta = pars[paste("beta[",1:idx[2], "]",sep="")] else beta = 0
	if(model$link=="glogistic"){
		idx[6] = 1
		k = pars["skew[k]"]
	} else{
		k = 1
	}
	
	n = length(y)
	mpu = double(n)
	# recursion initialization
	rcs = mean(y)*delta[1]
	if(idx[2]>0) rcs = rcs + sum(colMeans(x)*beta)
	if(idx[3]>0) rcs = rcs/(1-alpha[1])
	mpu[1] = as.double(rcs)
	
	temp = try(.C("c_dbmfilter", y = as.double(y), x = as.double(as.vector(x)), 
					mpu = as.double(mpu), 
					omega = as.double(omega), alpha = as.double(alpha), 
					delta = as.double(delta), beta = as.double(beta), 
					k = as.double(k), idx = as.integer(idx), xidx = as.integer(xidx), 
					lik = double(n), T = as.integer(c(0,n)), 
					PACKAGE="dbm"), silent = TRUE)
	fit = list()
	if(idx[6]==1) fit$fitted.values = pnorm(temp$mpu) else fit$fitted.values = plogis(temp$mpu)
	fit$log.likelihoods = temp$lik
	fit$LLH = sum(temp$lik)
	fit$coef = pars
	fit$matcoef = pars
	fit$mpu = temp$mpu
	fit$residuals = (y - fit$mpu)/(fit$mpu*(1-fit$mpu))
	if(idx[6]==1){
		gfit<-eval(parse(text=paste("glm(",yname,"~1, family=binomial(probit), data=as.data.frame(y))",sep="")))
		null.deviance = deviance(gfit)
		null.epcp = epcp.default(fitted(gfit), y = as.numeric(model$y))
	} else if(idx[6]==2){
		gfit<-eval(parse(text=paste("glm(",yname,"~1, family=binomial(logit), data=as.data.frame(y))",sep="")))
		null.deviance = deviance(gfit)
		null.epcp = epcp.default(fitted(gfit), y = as.numeric(model$y))
	} else{
		# need to create a function for this
		if(!is.null(model$x.vars) | arp>0 | arq >0){
			newmod = dbm(model$y[,yname], x.vars = NULL, x.lags = 0, arp = 0, arq = 0, 
					ecm = FALSE, constant = TRUE, link = "glogistic", fixed.pars = NULL, 
					solver = "optim", control=list(maxit=50000,trace=0), 
					parsearch = FALSE, parsim = 5000, method = "Nelder-Mead")
			null.deviance = deviance(newmod)
			null.epcp = epcp.dbm(newmod)
		} else {
			null.deviance = NA
			null.epcp = NA
		}
	}
	fit$null.epcp = null.epcp
	fit$null.deviance = null.deviance
	fit$model.epcp = epcp.default(fit$fitted.values, as.numeric(model$y))
	model$estimation = "filter"
	out = list(model = model, fit = fit, call = call)
	class(out) <- "dbm"
	return(out)
}