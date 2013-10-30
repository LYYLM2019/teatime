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


mfx = function(model, ...)
{
	UseMethod("mfx")
}

#------------------------------------------------------------------------------
# Functions adapted from the erer package of Changyou Sun
# maBina function (renamed and re-worked to mfx: marginal effects)

mfx.dbm <- function(model, x.mean = TRUE, rev.dum = TRUE, ...)
{
    # if (!inherits(w, "glm")) {stop("Please provide an object from 'glm()'.\n")}
    # link <- w$family$link
	link = model$model$link
    #if (link != "probit" & link != "logit") {
    #	stop("Need a binary probit or logit model'.")}
    #if (is.null(dim(w$x))) {
    #   stop("Please specify 'x = TRUE' in glm().\n")}
 	#x <- as.matrix(w$x)
	if(link=="glogistic") stop("\nglogistic not yet supported")
	x = coredata(model.matrix(model))
    x.bar <- as.matrix(colMeans(x))
    b.est <- as.matrix(coef(model))
    K <- nrow(b.est)
    xb <- t(x.bar) %*% b.est
    if(link == "gaussian") f.xb <- dnorm(xb)    
    if(link == "logistic" ) f.xb <- dlogis(xb)
	if(link == "glogistic" ) f.xb <- dlogis(xb)^coef(model)["skew[k]"]
    if(!x.mean){
    	xb2 <- x %*% b.est
		if (link == "gaussian") f.xb <- mean(dnorm(xb))    
    	if (link == "logistic" ) f.xb <- mean(dlogis(xb))
		if (link == "glogistic" ) f.xb <- mean(dlogis(xb)^coef(model)["skew[k]"])
    }
    me <- f.xb * coef(model)
    bx <- b.est %*% t(x.bar)
    if(link == "gaussian") {  
    	dr <- diag(1, K, K) - as.numeric(xb) * bx
    	va <- as.numeric(f.xb)^2 * dr %*% vcov(model) %*% t(dr)
    } else {
    	pg <- as.numeric(plogis(xb))
    	dr <- diag(1, K, K) + (1 - 2 * pg) * bx
    	va <- (pg*(1-pg))^2 * dr %*% vcov(model) %*% t(dr)      
    }
    se <- sqrt(diag(va))
    if (rev.dum) {
    	for (i in 1:ncol(x)){
        	if (identical(unique(x[,i]), c(0, 1))) {
        		x.d1 <- x.bar; x.d1[i, 1] <- 1
        		x.d0 <- x.bar; x.d0[i, 1] <- 0
        		if (link == "gaussian") {
            		me[i] <- pnorm(t(x.d1) %*% b.est) -  
                    		pnorm(t(x.d0) %*% b.est)        
            		dr2 <- dnorm(t(x.d1) %*% b.est) %*%  t(x.d1) -  
                    		dnorm(t(x.d0) %*% b.est) %*%  t(x.d0) 
            		va2 <- dr2 %*% vcov(model) %*% t(dr2)
            		se[i] <- sqrt(as.numeric(va2))                   
        		}
        		if (link == "logistic") {
            		me[i] <- plogis(t(x.d1) %*% b.est) -  
                    		plogis(t(x.d0) %*% b.est)
            		dr2 <- dlogis(t(x.d1) %*% b.est) %*%  t(x.d1) -  
                    		dlogis(t(x.d0) %*% b.est) %*%  t(x.d0) 
            		va2 <- dr2 %*% vcov(model) %*% t(dr2)
            		se[i] <- sqrt(as.numeric(va2))  
        		}
        	}
    	}
    }
    out <- data.frame(effect=me, error=se)
    out$t.value <- out$effect / out$error
    out$p.value <- 2*(1- pt(abs(out[, 3]), model$fit$df.residual))
    result <- list(link=link, f.xb=f.xb, coef=coef(model), 
			x = model.matrix(model), yname = model$model$yname, out=out)
    class(result) <- "mfx"
    return(result)
}

# maTrend function:
mfxTrend =  function(object, ...)
{
	UseMethod("mfxTrend")
}

mfxTrend.mfx <- function(object, n = 300, cregressor, bregressor, ...)
{
    if(!inherits(object, "mfx")) {stop("Need an object from 'mfx()'.\n")}
    if(missing(cregressor)) {stop("Need a regressors name'.\n")}
	x = coredata(object$x)
    if (identical(sort(unique(x[, cregressor])), c(0, 1))) {
    	stop("cregressor must be a continuous variable.") }
    if (!missing(bregressor)) {
    	if (!identical(sort(unique(x[, bregressor])), c(0, 1))) {
        	stop("bregressor must be a binary variable.") } }
    link <- object$link
    b.est <- as.matrix(object$coef)
    result <- list(x = x, cregressor = cregressor)
    mm <- matrix(colMeans(x), ncol=ncol(x), nrow=n, byrow=TRUE)
    colnames(mm) <- colnames(x)
    ran <- range(x[, cregressor])
    tre <- seq(from=ran[1], to=ran[2], length.out=n)
    mm[, cregressor] <- tre
    if(link =="gaussian") { pp <- pnorm( mm %*% b.est)}      
    if(link =="logistic") { pp <- plogis(mm %*% b.est)}
	if(link == "glogistic"){ pp <- plogis(mm %*% b.est)^coef(x$model)["skew[k]"] }
    trend <- data.frame(mm[, cregressor], pp)
    colnames(trend) <- c(cregressor, "pr.all")
    result$mm <- mm
	result$yname = object$yname
    if (!missing(bregressor)) {
    	m1 <- mm; m1[, bregressor] <- 1
    	m0 <- mm; m0[, bregressor] <- 0
    	if(link == "gaussian") {
        	p1 <- pnorm(m1 %*% b.est)
        	p0 <- pnorm(m0 %*% b.est)          
    	} else if(link == "logistic"){
        	p1 <- plogis(m1 %*% b.est)
        	p0 <- plogis(m0 %*% b.est)   
    	} else {
			p1 <- plogis(m1 %*% b.est)^object$coef["skew[k]"]
        	p0 <- plogis(m0 %*% b.est)^object$coef["skew[k]"]  
		}
    	trend <- data.frame(mm[, cregressor], pp, p1, p0)    
    	colnames(trend) <- c(cregressor, "pr.all", 
        		paste("pr", bregressor, "d1", sep="."), 
        		paste("pr", bregressor, "d0", sep="."))
    	result$bregressor <- bregressor
    	result$m1 <- m1 
    	result$m0 <- m0 
    }   
    result$trend <- trend
    class(result) <- "mfxTrend"
    return(result)
}

# plot.maTrend
plot.mfxTrend <- function(x, ...)
{ 
    pr <- x$trend  
    plot(pr[, 2] ~ pr[, 1], type="l", lty=1,
    		ylim=c(min(pr[, -1]), max(pr[, -1])), 
    		xlab=toupper(x$cregressor),
    		ylab=paste("Probability (", x$yname, " = 1)", sep=""), ... )
    abline(v=mean(x$x[,x$cregressor]), lty=4)
    grid()
    if (!is.null(x$bregressor)) {  
    	lines(pr[, 3] ~ pr[, 1], type="l", lty=2)
    	lines(pr[, 4] ~ pr[, 1], type="l", lty=3)
    }
	return(invisible())
}

summary.mfx = function(object, ...)
{
	return(object$out)
}

print.mfx = function(x, ...){
	print.data.frame(x$out, ...)
	return(invisible())
}