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
logLik.dbm = function(object, ...)
{
	return( object$fit$LLH )
}

score = function(object, pars = NULL, analytic = TRUE, ... ) { UseMethod("score") }

score.dbm = function(object, pars = NULL, analytic = TRUE, ... )
{
	if(is.null(pars)){
		if(analytic) ans = object$fit$analytic.scores else ans = object$fit$numeric.scores
	} else{
		# check parameters
		mnames = object$model$modelnames
		fnames = names(pars)
		chk = match(fnames, mnames)
		if(any(is.na(chk))){
			cx = which(is.na(chk))
			cat("\nunidentified parameter names in pars\n")
			cat("\nexpected:",mnames)
			cat("\ngot:",fnames)
			cat("\n")
			stop()
		}
		arglist = list()
		arglist$idx = object$model$idx
		arglist$xidx = object$model$xidx
		arglist$y = as.numeric(object$model$y[,object$model$yname])
		arglist$x = coredata(object$model$y[,object$model$x.vars])
		arglist$pnames = fnames
		arglist$fnames = NULL
		arglist$regularization = object$model$regularization
		if(arglist$idx[7]==1) arglist$Cost = object$model$reg.cost else arglist$Cost = 0
		arglist$modelnames = fnames
		arglist$fpars = NULL
		arglist$fidx = NULL
		arglist$type = "lik"
		arglist$derivarg = 2
		dbmenv = new.env(hash = TRUE)
		assign("dbm_llh", 0, envir = dbmenv)
		arglist$dbmenv <- dbmenv
		arglist$transform = FALSE
		if(analytic){
			ans = switch(as.character(object$model$idx[6]),
					"1" = dbmderiv2(pars, arglist),
					"2" = dbmderiv1(pars, arglist),
					"3" = dbmderiv3(pars, arglist))
		} else{
			ans = jacobian(func = dbmlik, pars, arglist = arglist)
		}
	}
	return( ans )
}


deviance.dbm = function(object, null = FALSE, ...)
{
	return( ifelse(null, object$fit$null.deviance, -2*logLik(object) ) )
}

plot.dbm = function(x, ...)
{
	old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
	devAskNewPage(ask = TRUE)
	y = x$model$y
	f = fitted(x)
	ep <- axTicksByTime(index(f))
	par(mar = c(2.5, 2.5, 2, 1))
	xx<-plot(as.numeric(f), type="l", xaxt = "n", ylab="",xlab="",main="dbm: Fitted vs Actual", ylim=c(-0.01, 1), yaxs = "i")
	tmp=y[index(f), x$model$yname]
	tmp = as.numeric(tmp)
	start = which(diff(c(0,tmp))==1)
    end = which(diff(c(0,tmp))==(-1))-1
	if(length(start)==0 && length(end)>0) start=1	
    if(length(end)<length(start)) end = c(end, length(tmp))
    if(length(end)>length(start)) start = c(1, start)
    miny = min(f, na.rm=TRUE)
    maxy = max(f, na.rm=TRUE)
    n = length(start)
    for(i in 1:n){
		rect(start[i], -0.01, end[i],  1.01, col= "WhiteSmoke", border = FALSE)
    }
	lines(as.numeric(f))
	axis(1, at = ep, labels = names(ep), tick = TRUE)
	box()
	abline(h=0.5, col = "steelblue", lty=2)
	y = as.numeric(tmp)
	f = as.numeric(f)
	par(mar = c(5, 4, 4, 2) + 0.1)	
	Epi::ROC(test = f, stat = y, plot = "ROC", PV=TRUE, MX=TRUE, AUC = TRUE, main = "ROC Curve")
	invisible()
}

coef.dbm = function(object, ...)
{
	return( object$fit$coef )
}

vcov.dbm = function(object, robust = FALSE, ...)
{
	if(robust) ans = object$fit$robust.cvar else ans = object$fit$cvar
	colnames(ans) = rownames(ans) = object$model$parnames
	return( ans )
}

fitted.dbm = function(object, ...)
{
	ans = object$fit$fitted.values
	ans = xts(ans, index(object$model$y))
	return(ans)
}

residuals.dbm = function(object, type = c("deviance", "pearson", "std.pearson"), ...)
{
	ans = switch(tolower(type[1]),
			deviance = resdeviance(object),
			pearson  = respearson(object),
			std.pearson = respearson(object))
	if(tolower(type)=="std.pearson"){
		L = hatvalues(object)
		ans = ans/sqrt(1-L)
	}
	ans = xts(ans, index(object$model$y))
	return(ans)
}

BIC.dbm = function(object, ...)
{
	nObs = NROW(object$model$y)
	nPars = length(coef(object))
	LLH = logLik(object)
	ans = (-2*LLH)/nObs + nPars * log(nObs)/nObs
	return( ans )
}

AIC.dbm = function(object, ...)
{
	nObs = NROW(object$model$y)
	nPars = length(coef(object))
	LLH = logLik(object)
	ans = (-2*LLH)/nObs + 2 * nPars/nObs
	return( ans )
}


summary.dbm = function(object, ...)
{
	ans = list()
	ans$coefficients = object$fit$matcoef
	ans$logLik = logLik(object)
	ans$model.deviance = deviance(object)
	ans$null.deviance=  deviance(object, TRUE)
	ans$n.obs = length(fitted(object))
	ans$model.epcp = object$fit$model.epcp
	ans$null.epcp = object$fit$null.epcp
	ans$mfRsq = 1 - (ans$logLik/(-ans$null.deviance/2))
	ans$csRsq = 1 - (ans$logLik/(-ans$null.deviance/2))^(2/ans$n.obs)
	ans$tjRsq = tjur.rsq(object)
	class(ans)<-"summary.dbm"
	return(ans)
}

print.dbm = function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
	y = summary(x)
	print.default(format(y$coefficients, digits = digits), 
            print.gap = 2, quote = FALSE)
	cat("\n")
	cat("\nN.observations:\t\t", y$n.obs)
	cat("\nLog-Likelihood:\t\t", y$logLik)
	cat("\nNull Deviance:\t\t", y$null.deviance)
	cat("\nResidual Deviance:\t", y$model.deviance)
	cat("\nMcFadden pseudo  R^2:", round(y$mfRsq,4))
	cat("\nTjur     pseudo  R^2:", round(y$tjRsq,4))
	# cat("\nCox-Snell pseudo R^2:", round(y$csRsq,5))
	cat("\nE[% Correctly predicted (model)]:", round(x$fit$model.epcp,3))
	cat("\nE[% Correctly predicted  (null)]:", round(x$fit$null.epcp,3))
	cat("\n")
	return(invisible(y))
}

tjur.rsq = function(x)
{
	y = x$model$y
	y = as.numeric(y[, if(!is.null(x$model$yname)) x$model$yname else colnames(y)[1]])
	p = as.numeric(fitted(x))
	ans = mean(p[y==1]) - mean(p[y==0])
	return(ans)
}
epcp.dbm = function(x, y = NULL)
{
	# expected percent correctly predicted	
	y = x$model$y
	y = as.numeric(y[, if(!is.null(x$model$yname)) x$model$yname else colnames(y)[1]])
	p = as.numeric(fitted(x))
	ans = (sum(p[y==1]) + sum(1-p[y==0]))/length(y)
	return(ans)
}

epcp.default = function(x, y = NULL)
{
	ans = (sum(x[y==1]) + sum(1-x[y==0]))/length(y)
	return(ans)
}


model.matrix.dbm = function(object, ...)
{
	
	if(object$model$constant){
		y = object$model$y[,-1]
		i = xts(rep(1, nrow(y)), index(y))
		colnames(i)<-"Intercept"
		ans = cbind(i, y)
	} else{
		ans = object$model$y[,-1]
	}
	idx = index(ans)
	ans = coredata(ans)
	
	if(object$model$arp>0){
		for(i in 1:object$model$arp){
			ans = cbind(ans, .lagx(object$fit$mpu, n.lag=i, pad = object$fit$rec.init))
			colnames(ans)[length(colnames(ans))] = paste("pi[t-",i,"]",sep="")
		}
	}
	if(object$model$arq>0 && !object$model$ecm){
		for(i in 1:object$model$arq){
			ans = cbind(ans, .lagx(coredata(object$model$y[,1]), n.lag=i, pad = 0))
			colnames(ans)[length(colnames(ans))] = paste("y[t-",i,"]",sep="")
		}
	}
	ans = xts(ans, idx)
	return(ans)
}

hat.dbm = function(x, ...)
{
	p = as.numeric(fitted(x))
	xp = p * (1-p)
	W.hat =  diag(xp)
	X = coredata(model.matrix(x))
	W.hat.i = .sqrtm(W.hat)
	H = W.hat.i %*% X %*% solve(t(X) %*% W.hat %*% X) %*% t(X) %*% W.hat.i
	return(H)
}

hatvalues.dbm = function(model, ...)
{
	return(diag(hat.dbm(model)))
}

###############################################################################
# Hosmer-Lemeshow Goodness of Fit (GOF) Test
# Code adapted from the ResourceSelection package of Solymos et al.
hoslem.test <- function(y, x, groups, yname) {
	METHOD <- "Hosmer and Lemeshow goodness of fit (GOF) test"
    qq <- unique(quantile(x, probs=seq(0, 1, 1/groups)))
    cutx <- cut(x, breaks = qq, include.lowest = TRUE)
    observed <- xtabs(cbind("y0" = 1 - y, "y1" = y) ~ cutx)
    expected <- xtabs(cbind("yhat0" = 1 - x, "yhat1" = x) ~ cutx)
    chisq <- sum((observed - expected)^2 / expected)
    PVAL = 1 - pchisq(chisq, groups - 2)
    PARAMETER <- groups - 2
    names(chisq) <- "X-squared"
    names(PARAMETER) <- "df"
    return(structure(list(statistic = chisq, parameter = PARAMETER, 
        			p.value = PVAL, method = METHOD, data = yname, 
					observed = observed, 
        			expected = expected), class = "htest"))
}
