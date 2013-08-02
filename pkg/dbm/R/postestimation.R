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
.postestimate = function(f, pars, arglist)
{
	# Turn Stationarity Off for numerical derivative calculation
	fit = list()
	flag = 0
	arglist$type = "llh"
	if(arglist$idx[6]==1){
		gr = dbmderiv2
	} else if(arglist$idx[6]==2){
		gr = dbmderiv1
	} else{
		gr = dbmderiv3
	}
	fit$hessian = optimHess(pars, fn = f, gr = gr, arglist=arglist)
	if(any(is.nan(fit$hessian)) | any(!is.finite(fit$hessian)) ){
		warning("\ndbm: hessian ill conditioned..returning only partial results")
		flag = 1
	} else{
		fit$cvar = try(solve(fit$hessian), silent = TRUE)
		# (might also fail in which case user will see an error about the inversion failure)
		if(inherits(fit$cvar, "try-error")){
			zz = try(solve(.hessian2sided(f, pars, arglist = arglist)), silent=TRUE)
			if(inherits(zz, "try-error")) {
				fit$cvar = NULL
				warning("\ndbm-->warning: failed to invert hessian\n")
				flag = 1
			} else{
				fit$cvar = zz
			}
		}
	}
	arglist$type = "ALL"
	temp = f(pars,	arglist)
	if(arglist$idx[6]==1){
		fit$fitted.values = pnorm(temp$mpu) 
	} else if(arglist$idx[6]==2){
		fit$fitted.values = plogis(temp$mpu)
	} else{
		zpars = c(pars, arglist$fpars)
		fit$fitted.values = plogis(temp$mpu)^zpars["skew[k]"]
	}
	fit$mpu = temp$mpu
	fit$LLH = -temp$llh
	fit$log.likelihoods = temp$lik
	fit$coef = pars
	fit$residuals = (arglist$y - fit$mpu)/(fit$mpu*(1-fit$mpu))
	names(fit$coef) <- arglist$pnames
	if(is.null(fit$cvar) | flag==1){
		m = length(pars)
		fit$se.coef = rep(NA, m)
		fit$tval = rep(NA, m)
		fit$matcoef = cbind(pars, fit$se.coef, fit$tval, rep(NA, m))
		fit$robust.se.coef = rep(NA, m)
		fit$robust.tval = rep(NA, m)
		fit$robust.matcoef = cbind(fit$coef, fit$robust.se.coef, fit$robust.tval, rep(NA, m))
		fit$hessian.message = "failed to invert hessian"
	} else{
		n = length(arglist$y)
		# return returns the maximized likelihood (unlike llh which is already
		# premultiplied by -1 for the minimization)
		arglist$type = "lik"
		tmp = try(robustvcv(fun = f, pars = pars, nlag = floor(4*(n/100)^(2/9)), hess = fit$hessian, n = n, 
				arglist = arglist), silent=TRUE)
		fit$robust.cvar = tmp$vcv
		fit$numeric.scores = jacobian(func = f, x = pars, method="Richardson", method.args=list(), arglist = arglist) 
		colnames(fit$numeric.scores) = arglist$pnames
		fit$se.coef = sqrt(diag(abs(fit$cvar)))
		fit$tval = fit$coef/fit$se.coef
		fit$matcoef = matrix(NA, ncol = 4, nrow = length(fit$coef))
		fit$matcoef = cbind(fit$coef, fit$se.coef, fit$tval, 2*(1-pnorm(abs(fit$tval))))
		if(!is.null(arglist$fidx)) fit$matcoef = rbind(fit$matcoef, 
					cbind(arglist$fpars, rep(NA, length(arglist$fpars)),
    						rep(NA, length(arglist$fpars)), rep(NA, length(arglist$fpars))))
		
		fit$robust.se.coef = sqrt(diag(fit$robust.cvar))
		fit$robust.tval = fit$coef/fit$robust.se.coef
		fit$robust.matcoef = cbind(fit$coef, fit$robust.se.coef,
					fit$robust.tval, 2*(1-pnorm(abs(fit$robust.tval))))
		if(!is.null(arglist$fidx)) fit$robust.matcoef = rbind(fit$robust.matcoef, 
					cbind(arglist$fpars, rep(NA, length(arglist$fpars)),
							rep(NA, length(arglist$fpars)), rep(NA, length(arglist$fpars))))
		fit$hessian.message = NULL
		if(arglist$idx[5]==1){
			fit$analytic.scores = NULL
		} else{
			arglist$derivarg=2
			if(arglist$idx[6]==1){
				fit$analytic.scores = dbmderiv2(pars, arglist)
			} else if(arglist$idx[6]==2){
				fit$analytic.scores = dbmderiv1(pars, arglist)
			} else{
				fit$analytic.scores = dbmderiv3(pars, arglist)
			}
			colnames(fit$analytic.scores) = arglist$pnames
		}
	}
	dimnames(fit$matcoef) = list(c(names(fit$coef), arglist$fnames), c(" Estimate",
					" Std. Error", " z value", "Pr(>|z|)"))
	dimnames(fit$robust.matcoef) = list(c(names(fit$coef), arglist$fnames), c(" Estimate",
					" Std. Error", " z value", "Pr(>|z|)"))

	return(fit)
}
# Fitting procedure required functions:
.hessian2sided = function(f, x, ...)
{
	n = length(x)
	fx = f(x, ...)
	eps = .Machine$double.eps
	# Compute the stepsize (h)
	# h = eps^(1/3)*apply(as.data.frame(x), 1,FUN = function(z) max(abs(z), 1e-4))
	h = apply(as.data.frame(x), 1,FUN = function(z) max(abs(z)*10e-4, 1e-9))
	xh = x+h
	h = xh-x
	if(length(h) == 1) ee = matrix(h, ncol = 1, nrow = 1) else ee = as.matrix(diag(h))
	
	# Compute forward and backward steps
	gp = vector(mode = "numeric", length = n)
	gp = apply(ee, 2, FUN = function(z) f(x+z, ...))
	gm = vector(mode="numeric",length=n)
	gm = apply(ee, 2, FUN = function(z) f(x-z, ...))
	H = h%*%t(h)
	Hm = H
	Hp = H
	# Compute "double" forward and backward steps
	for(i in 1:n){
		for(j in  i:n){
			Hp[i,j] = f(x+ee[,i]+ee[,j], ...)
			Hp[j,i] = Hp[i,j]
			Hm[i,j] = f(x-ee[,i]-ee[,j], ...)
			Hm[j,i] = Hm[i,j]
		}
	}
	#Compute the hessian
	for(i in 1:n){
		for(j in  i:n){
			H[i,j] = ( (Hp[i,j]-gp[i]-gp[j]+fx+fx-gm[i]-gm[j]+Hm[i,j]) /H[i,j] )/2
			H[j,i] = H[i,j]
		}
	}
	return(H)
}


.hessianan = function(f, x, ...)
{
	n = length(x)
	fx = f(x, ...)
	eps = .Machine$double.eps
	# Compute the stepsize (h)
	# h = eps^(1/3)*apply(as.data.frame(x), 1,FUN = function(z) max(abs(z), 1e-4))
	h = apply(as.data.frame(x), 1,FUN = function(z) max(abs(z)*10e-4, 1e-9))
	xh = x+h
	h = xh-x
	if(length(h) == 1) ee = matrix(h, ncol = 1, nrow = 1) else ee = as.matrix(diag(h))
	
	# Compute forward and backward steps
	gp = vector(mode = "numeric", length = n)
	gp = apply(ee, 2, FUN = function(z) f(x+z, ...))
	gm = vector(mode="numeric",length=n)
	gm = apply(ee, 2, FUN = function(z) f(x-z, ...))
	H = h%*%t(h)
	Hm = H
	Hp = H
	# Compute "double" forward and backward steps
	for(i in 1:n){
		for(j in  i:n){
			Hp[i,j] = f(x+ee[,i]+ee[,j], ...)
			Hp[j,i] = Hp[i,j]
			Hm[i,j] = f(x-ee[,i]-ee[,j], ...)
			Hm[j,i] = Hm[i,j]
		}
	}
	#Compute the hessian
	for(i in 1:n){
		for(j in  i:n){
			H[i,j] = ( (Hp[i,j]-gp[i]-gp[j]+fx+fx-gm[i]-gm[j]+Hm[i,j]) /H[i,j] )/2
			H[j,i] = H[i,j]
		}
	}
	return(H)
}

# The following functions are based on on Kevin Sheppard's MFE toolbox
#---------------------------------------------------------------------
neweywestcv = function(data, nlag = NULL, center = TRUE)
{
	# Long-run covariance estimation using Newey-West (Bartlett) weights
	#  if nlag empty=NULL NLAG=min(floor(1.2*T^(1/3)),T)
	N = dim(as.matrix(data))[1]
	if(is.null(nlag)) nlag=min(floor(1.2*N^(1/3)),N)
	if(center) data = apply(data, 2, FUN = function(x) scale(x, center = TRUE, scale = FALSE))
	# weights
	bw = (nlag+1-(seq(0,nlag,by=1)))/(nlag+1)
	cv = 1/N * t(data)%*%data
	for(i in 1:nlag){
		gmi = 1/N * (t(data[(i+1):N,])%*%data[1:(N-i),])
		gpp = gmi + t(gmi)
		cv = cv + bw[i+1]*gpp
	}
	return(cv)
}

robustvcv = function(fun, pars, nlag = 0, hess, n, ...)
{
	arglist = list(...)$arglist
	arglist$derivarg = 2
	scores = switch(as.character(arglist$idx[6]),
			"1" = dbmderiv2(pars, arglist),
			"2" = dbmderiv1(pars, arglist),
			"3" = dbmderiv3(pars, arglist))
	A = hess/n
	hess = A
	Ainv = try( solve(A), silent = TRUE )
	if( inherits(Ainv, "try-error")){
		info = 1
		vcv = NA
	} else{
		if(nlag>0){
			B = neweywestcv(scores, nlag)
			vcv = (Ainv%*%B%*%Ainv)/n
		} else{
			B = cov(scores)
			vcv = (Ainv%*%B%*%Ainv)/n
		}
		info = 0
	}
	return(list(vcv = vcv, scores = scores, info = info))
}




.boxcoxtransform = function(x, lambda)
{
	if(lambda!=0) ret = (x^lambda - 1)/lambda else ret = log(x)
	return(ret)
}

repmat = function(a, n, m)
{
	kronecker(matrix(1, n, m), a)
}
size = function(x, n = NULL)
{
	x = as.matrix(x)
	if(missing(n)) sol = c(n = dim(x)[1], m = dim(x)[2]) else sol = dim(x)[n]
	return(sol)
}

zeros = function(n = 1, m = 1)
{
	if(missing(m)) m = n
	sol = matrix(0, nrow = n, ncol = m)
	return(sol)
}

ones = function(n = 1, m = 1)
{
	if(missing(m)) m = n
	sol = matrix(1, nrow = n, ncol = m)
	return(sol)
}

newlagmatrix = function(x,nlags,xc)
{
	nlags = nlags+1
	xt = size(x, 1);
	newX = rbind(x, zeros(nlags, 1))
	lagmatrix = repmat(newX, nlags, 1)
	lagmatrix = matrix(lagmatrix[1:(size(lagmatrix,1)-nlags)], nrow = (xt+nlags-1), ncol = nlags)
	lagmatrix = lagmatrix[nlags:xt,]
	y = lagmatrix[,1]
	x = lagmatrix[,2:nlags]
	if(xc == 1) x = cbind(ones(size(x,1), 1), x)
	return(data.frame(y = y, x = x))
}