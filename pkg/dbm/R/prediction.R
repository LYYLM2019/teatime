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
# iterated forecast routine
# expand.grid(c(1,0),c(1,0),c(1,0))
# eval(parse(text=paste(c("expand.grid(",rep("c(1,0),",k-1),"c(1,0))"),sep=" ", collapse="")))
# idx [omega, x.lags, length(x.vars), arp, arq, use_ecm, link]

predict.dbm = function(object, newdata = NULL, n.ahead=1, ...)
{
	x.vars = object$model$x.vars
	idx = object$model$idx
	xidx = object$model$xidx
	if(n.ahead>15) stop("\nn.ahead>15 is likely to lead to memory problems...aborting.")
	pars = object$fit$matcoef[,1]
	if(idx[1]>0) omega = pars["omega"] else omega = 0
	if(idx[3]>0){
		alpha = pars[paste("alpha[",1:idx[3], "]",sep="")]
		# if(object$model$transform) alpha = logtransform(alpha, -1, 1)
	} else{
		alpha = 0
	}
	if(idx[5]>0){
		delta = 0
	} else{
		if(idx[4]>0) delta = pars[paste("delta[",1:idx[4], "]",sep="")] else delta = 0
	}
	if(idx[2]>0) beta = pars[paste("beta[",1:idx[2], "]",sep="")] else beta = 0
	if(idx[6]==3) kappa = pars["skew[k]"] else kappa = 1
	y = object$model$y[,object$model$yname]
	if(idx[2]>0) x = coredata(object$model$y[,object$model$x.vars]) else x = 0
	T0 = NROW(y)
	mpu = object$fit$mpu
	y = as.numeric(y)
	pred = rep(NA, n.ahead)
	for(j in 1:n.ahead){
		if(j==1){
			B = NULL
			TN = T0+1
			zeropad = rep(0, 1)
			yf = c(y, zeropad)
			mpuf = c(mpu, zeropad)
			f = try(.C("c_dbmfilter", y = as.double(y), x = as.double(as.vector(x)), 
					mpu = as.double(mpuf), omega = as.double(omega),
					alpha = as.double(alpha), delta = as.double(delta), 
					beta = as.double(beta), k = as.double(kappa), 
					idx = as.integer(idx), xidx = as.integer(xidx),
					lik = double(TN), T = as.integer(c(T0, TN)), PACKAGE="dbm"), silent = TRUE)
			pred[j] = ifelse(idx[6]==1, pnorm(tail(f$mpu, 1)), plogis(tail(f$mpu, 1)))
		} else{
			TN = T0+j
			zeropad = rep(0, j)
			yf = c(y, zeropad)
			mpuf = c(mpu, zeropad)
			# evaluate all possible paths
			if(j==2){
				B = matrix(c(1,0), nrow=2)
			} else{
				B = as.matrix(eval(parse(text=paste(c("expand.grid(",rep("c(1,0),", j-2),"c(1,0))"),sep=" ", collapse=""))))
			}
			#colnames(B) = paste("T+",1:(j-1),sep="")
			m = nrow(B)
			pk = rep(0, m)
			eq1 = rep(1, m)
			eq2 = rep(0, m)
			for(i in 1:m){
				ytmp = as.numeric(c(y, B[i,]))
				eq2[i] = .Call("dbmpk2", p = mpuf, omega = as.numeric(omega), alpha = as.numeric(alpha), 
						delta = as.numeric(delta), beta = as.numeric(beta), kappa = as.double(kappa),
						y = ytmp, x = x, idx = as.integer(idx),
						xidx = as.integer(xidx), V = as.integer(c(j, T0, idx[6])), PACKAGE="dbm")
				for(k in 1:(j-1)){
					tmp = .Call("dbmpk1", p = mpuf, omega = as.numeric(omega), alpha = as.numeric(alpha), 
							delta = as.numeric(delta), beta = as.numeric(beta), kappa = as.double(kappa),
							y = ytmp, x = x, idx = as.integer(idx),
									xidx = as.integer(xidx), V = as.integer(c(k, T0, idx[6])), PACKAGE="dbm")
					eq1[i] = eq1[i] * (tmp^ytmp[T0+k] * (1-tmp)^(1-ytmp[T0+k]))
				}
			}
			pred[j] = sum(eq1 * eq2)
			B = cbind(B, eq1 * eq2)
			colnames(B) = c(paste("T+",1:(j-1),sep=""), paste("P(y[T+",j,"]=1)",sep=""))
		}
	}
	names(pred)<-paste("T+",1:n.ahead,sep="")
	return(list(prediction = pred, path = B))
}