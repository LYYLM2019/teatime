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
TinY = 1.0e-8
.meqstartpars = function(pars, arglist)
{
	# Case 1 garchInMean yes:
	data = arglist$data
	N = length(as.numeric(unlist(data)))
	model = arglist$model
	modelinc = model$modelinc
	modeldesc = model$modeldesc
	idx = model$pidx
	mxreg = numeric()
	if(modelinc[8] > 0){
		mexdata = model$modeldata$mexdata[1:N,, drop = FALSE]
		# easier to NULL the names of the data rather than search them by name later
		# (search now is 'mexdata'..1:mxn)
		colnames(mexdata) = NULL
	} else{
		mexdata = NULL
	}
	tmph = 0
	# get the sigma vector should it be needed (garchInMean)
	if(modelinc[5] > 0 && is.null(model$start.pars$archm)){
		tmph = sqrt(ewmav(scale(data, scale = FALSE)))
	}
	# arima without garchInMean
	if( (modelinc[5] == 0 ||  !is.null(model$start.pars$archm)) && (modelinc[2]>0 | modelinc[3]>0)){
		ttemp = arima(data, order = c(modelinc[2], 0, modelinc[3]), include.mean = modelinc[1], xreg = mexdata, method = "CSS")
		fit.mean = ttemp$coef
		#res=ttemp$residuals
		if(modelinc[1]>0){
			pars[idx["mu", 1]:idx["mu", 2], 1] = fit.mean["intercept"]
		}
		if(modelinc[2]>0) pars[idx["ar", 1]:idx["ar", 2], 1] = fit.mean[c(paste("ar",1:modelinc[2],sep=""))]
		if(modelinc[3]>0) pars[idx["ma", 1]:idx["ma", 2], 1] = fit.mean[c(paste("ma",1:modelinc[3],sep=""))]		
		if(modelinc[8]>0){
			i = which(substr(names(fit.mean), 1, 7) == "mexdata")
			pars[idx["mxreg", 1]:idx["mxreg", 2], 1] = fit.mean[i]
		}
	}
	
	# arima with garchInMean	
	if((modelinc[5]>0) && (modelinc[2]+modelinc[3])>0 && is.null(model$start.pars$archm)){
		mexdata = cbind(mexdata, tmph^modelinc[5])
		mxn = modelinc[8]+1
		colnames(mexdata) = paste("xreg", 1:mxn,sep="")
		ttemp = arima0(data, order = c(modelinc[2], 0, modelinc[3]), include.mean = modelinc[1], xreg = mexdata)
		fit.mean = ttemp$coef
		if(modelinc[1]>0){
			pars[idx["mu", 1]:idx["mu", 2], 1] = fit.mean["intercept"]
		}
		if(modelinc[2]>0) pars[idx["ar", 1]:idx["ar", 2], 1] = fit.mean[c(paste("ar",1:modelinc[2],sep=""))]
		if(modelinc[3]>0) pars[idx["ma", 1]:idx["ma", 2], 1] = fit.mean[c(paste("ma",1:modelinc[3],sep=""))]	
		i = which(substr(names(fit.mean), 1, 4) == "xreg")
		z = length(i) # at a minimum it is 2 ex+inmean
		if(modelinc[8]>0){
			pars[idx["mxreg", 1]:idx["mxreg", 2], 1] = fit.mean[i[1:(z-1)]]
		}
		pars[idx["archm", 1]:idx["archm", 2], 1] = fit.mean[i[z]]
	}
	
	# lm for garchInMean without arma
	if((modelinc[5]>0) && (modelinc[2] == 0 && modelinc[3] == 0) && is.null(model$start.pars$archm)){
		mexdata = cbind(mexdata, tmph^modelinc[5])
		mxn = modelinc[8]
		y = data
		if(modelinc[1]>0){
			fit.mean = lm(y~mexdata)
			pars[idx["mu", 1]:idx["mu", 2], 1] = fit.mean$coef["(Intercept)"]				
			i = which(substr(names(fit.mean$coef), 1, 7) == "mexdata")
			z = length(i) # at a minimum it is 2 ex+inmean
			if(modelinc[8]>0){
				pars[idx["mxreg", 1]:idx["mxreg", 2], 1] = fit.mean$coef[i[1:(z-1)]]
			}
			pars[idx["archm", 1]:idx["archm", 2], 1] = fit.mean$coef[i[z]]				
			#res=as.numeric(fit.mean$residuals)
		} else{
			fit.mean = lm(y~mexdata-1)
			if(modelinc[8]>0) pars[idx["mxreg", 1]:idx["mxreg", 2], 1] = fit.mean$coef[c(paste("mexdata",1:(mxn-1),sep=""))]			
			pars[idx["archm", 1]:idx["archm", 2], 1] = fit.mean$coef[c(paste("mexdata", mxn,sep=""))]
			#res=as.numeric(fit.mean$residuals)
		}
	} 
	
	if(modelinc[5]==0 && modelinc[2] == 0 && modelinc[3] == 0){
		y = data
		if(modelinc[8]>0){
			mxn = modelinc[8]
			fit.mean = lm(y~mexdata)
			i = which(substr(names(fit.mean$coef), 1, 7) == "mexdata")
			pars[idx["mu", 1]:idx["mu", 2], 1] = fit.mean$coef["(Intercept)"]
			pars[idx["mxreg", 1]:idx["mxreg", 2], 1] = fit.mean$coef[i]
			#res=as.numeric(fit.mean$residuals)
		} else{
			pars[idx["mu", 1]:idx["mu", 2], 1] = 0
			#res=(data-mu)
		}
	}
	arglist$tmph = tmph
	return(list(pars =  pars, arglist = arglist))
}

# common to all specifications is the mean equation:
.meqstart = function(pars, arglist)
{
	dscale = arglist$dscale
	model = arglist$model
	start.pars = model$start.pars
	start.names = names(start.pars)
	
	fixed.pars = model$fixed.pars
	fixed.names = names(fixed.pars)
	
	idx = model$pidx
	modelinc = model$modelinc
	data = arglist$data
	# this is where we fix the bounds for the fixed parameters
	# fill and then the fixed.names to overwrite it...also for the
	# bounds
	
	if(modelinc[1]>0){
		pars[idx["mu", 1]:idx["mu", 2], 5] = -100*abs(mean(data))
		pars[idx["mu", 1]:idx["mu", 2], 6] =  100*abs(mean(data))
		if(!is.null(start.pars$mu)) pars[idx["mu", 1]:idx["mu", 2], 1] = start.pars$mu[1]/dscale
		if(any(substr(fixed.names, 1, 2)=="mu")){
			pars[idx["mu", 1]:idx["mu", 2], 1] = as.numeric(fixed.pars$mu)
			pars[idx["mu", 1]:idx["mu", 2], 5] = fixed.pars$mu
			pars[idx["mu", 1]:idx["mu", 2], 6] = fixed.pars$mu
		}
	}
	
	# ar (we changed the naming of darfima into arfima which creates some extra problems
	# to be caught when using the substr function
	if(modelinc[2]>0){
		arnames = paste("ar",1:modelinc[2],sep="")
		pars[idx["ar", 1]:idx["ar", 2], 5] = -1+TinY
		pars[idx["ar", 1]:idx["ar", 2], 6] =  1-TinY
		if(any(substr(start.names, 1, 2)=="ar")){
			j = which(substr(start.names, 1, 2)=="ar")
			armatch = charmatch(start.names[j],arnames)
			pars[arnames[armatch], 1]=as.numeric(start.pars[j])
		}
		if(any(substr(fixed.names, 1, 2)=="ar")){
			j = which(substr(fixed.names, 1, 2)=="ar")
			armatch = charmatch(fixed.names[j],arnames)
			pars[arnames[armatch], 1] = as.numeric(fixed.pars[j])
			pars[arnames[armatch], 5] = as.numeric(fixed.pars[j])
			pars[arnames[armatch], 6] = as.numeric(fixed.pars[j])
		}	
	}
	# ma
	if(modelinc[3]>0){
		manames = paste("ma",1:modelinc[3],sep="")
		pars[idx["ma", 1]:idx["ma", 2], 5] = -1+TinY
		pars[idx["ma", 1]:idx["ma", 2], 6] =  1-TinY
		if(any(substr(start.names, 1, 2)=="ma")){
			j = which(substr(start.names, 1, 2)=="ma")
			mamatch = charmatch(start.names[j],manames)
			pars[manames[mamatch], 1]=as.numeric(start.pars[j])
		}
		if(any(substr(fixed.names, 1, 2)=="ma")){
			j = which(substr(fixed.names, 1, 2)=="ma")
			mamatch = charmatch(fixed.names[j],manames)
			pars[manames[mamatch], 1] = as.numeric(fixed.pars[j])
			pars[manames[mamatch], 5] = as.numeric(fixed.pars[j])
			pars[manames[mamatch], 6] = as.numeric(fixed.pars[j])
		}
	}
	
	# garch in mean
	if(modelinc[5]>0){
		pars[idx["archm", 1]:idx["archm", 2], 5] = -10
		pars[idx["archm", 1]:idx["archm", 2], 6] =  10
		if(!is.null(start.pars$archm)) pars[idx["archm", 1]:idx["archm", 2], 1] = start.pars$archm[1]
		if(any(substr(fixed.names, 1, 5)=="archm")){
			pars[idx["archm", 1]:idx["archm", 2], 1] = as.numeric(fixed.pars$archm)
			pars[idx["archm", 1]:idx["archm", 2], 5] = fixed.pars$archm
			pars[idx["archm", 1]:idx["archm", 2], 6] = fixed.pars$archm
		}
	}
	
	if(modelinc[6]>0){
		pars[idx["skewm", 1]:idx["skewm", 2], 1] = 0.1
		pars[idx["skewm", 1]:idx["skewm", 2], 5] = -10
		pars[idx["skewm", 1]:idx["skewm", 2], 6] =  10
		if(!is.null(start.pars$skewm)) pars[idx["skewm", 1]:idx["skewm", 2], 1] = start.pars$skewm[1]
		if(any(!is.na(match(fixed.names, "skewm")))){			
			pars[idx["skewm", 1]:idx["skewm", 2], 1] = as.numeric(fixed.pars$skewm)
			pars[idx["skewm", 1]:idx["skewm", 2], 5] = fixed.pars$skewm
			pars[idx["skewm", 1]:idx["skewm", 2], 6] = fixed.pars$skewm
		}
	}
	if(modelinc[7]>0){
		pars[idx["shapem", 1]:idx["shapem", 2], 1] = 0.1
		pars[idx["shapem", 1]:idx["shapem", 2], 5] = -10
		pars[idx["shapem", 1]:idx["shapem", 2], 6] =  10
		if(!is.null(start.pars$shapem)) pars[idx["shapem", 1]:idx["shapem", 2], 1] = start.pars$shapem[1]
		if(any(!is.na(match(fixed.names, "shapem")))){			
			pars[idx["shapem", 1]:idx["shapem", 2], 1] = as.numeric(fixed.pars$shapem)
			pars[idx["shapem", 1]:idx["shapem", 2], 5] = fixed.pars$shapem
			pars[idx["shapem", 1]:idx["shapem", 2], 6] = fixed.pars$shapem
		}
	}
	
	# arfima
	if(modelinc[4]>0){
		pars[idx["arfima", 1]:idx["arfima", 2], 5] = -1
		pars[idx["arfima", 1]:idx["arfima", 2], 6] =  1
		if(is.null(start.pars$arfima)) pars[idx["arfima", 1]:idx["arfima", 2], 1] = rugarch:::.rsfit(data)-0.5 else pars[idx["arfima", 1]:idx["arfima", 2], 1] = start.pars$arfima[1]
		if(any(substr(fixed.names, 1, 6)=="arfima")){
			pars[idx["arfima", 1]:idx["arfima", 2], 1] = as.numeric(fixed.pars$arfima)
			pars[idx["arfima", 1]:idx["arfima", 2], 5] = as.numeric(fixed.pars$arfima)
			pars[idx["arfima", 1]:idx["arfima", 2], 6] = as.numeric(fixed.pars$arfima)
		}
	}
	
	# exogenous regressors
	if(modelinc[8]>0){
		mxnames = paste("mxreg",1:modelinc[8],sep="")
		pars[idx["mxreg", 1]:idx["mxreg", 2], 5] = as.numeric(abs(pars[idx["mxreg", 1]:idx["mxreg", 2], 1]))*-100
		pars[idx["mxreg", 1]:idx["mxreg", 2], 6] = as.numeric(abs(pars[idx["mxreg", 1]:idx["mxreg", 2], 1]))* 100
		if(any(substr(start.names, 1, 5)=="mxreg")){
			j = which(substr(start.names, 1, 5)=="mxreg")
			mxmatch = charmatch(start.names[j],mxnames)
			pars[mxnames[mxmatch], 1] = as.numeric(start.pars[j])
		}
		if(any(substr(fixed.names, 1, 5)=="mxreg")){
			j = which(substr(fixed.names, 1, 5)=="mxreg")
			mxmatch = charmatch(fixed.names[j],mxnames)
			pars[mxnames[mxmatch], 1] = as.numeric(fixed.pars[j])
			pars[mxnames[mxmatch], 5] = as.numeric(fixed.pars[j])
			pars[mxnames[mxmatch], 6] = as.numeric(fixed.pars[j])
		}
	}		
	return( pars )
}

# starting parameters s.t. specification
.acdstart = function(pars, arglist)
{
	tmp = .meqstartpars(pars, arglist)
	pars = tmp$pars
	arglist = tmp$arglist
	ans = switch(arglist$model$vmodel$model,
			sGARCH = .gjracdstart(pars, arglist),
			csGARCH = .csacdstart(pars, arglist),
			mcsGARCH = .mcsacdstart(pars, arglist))
	return(ans)
}

# GARCH model start parameters
.gjracdstart = function(pars, arglist)
{
	model = arglist$model
	dscale = arglist$dscale
	modelinc = model$modelinc
	start.pars = model$start.pars
	fixed.pars = model$fixed.pars
	idx = model$pidx
	fixed.names = names(fixed.pars)
	start.names = names(start.pars)
	pars = .meqstart(pars, arglist)
	shape0 = arglist$shape0
	skew0 = arglist$skew0
	data  = arglist$data
	if(modelinc[9]>0){
		pars[idx["omega", 1]:idx["omega", 2], 5] = var(data)/100000
		pars[idx["omega", 1]:idx["omega", 2], 6] = var(data)*100000
		if(is.null(start.pars$omega)) pars[idx["omega", 1]:idx["omega", 2], 1] = (var(data, na.rm = TRUE))/1000 else 
			pars[idx["omega", 1]:idx["omega", 2], 1] = start.pars$omega[1]/dscale
		if(any(substr(fixed.names, 1, 5) == "omega")){
			pars[idx["omega", 1]:idx["omega", 2], 1] = as.numeric(fixed.pars$omega)
			pars[idx["omega", 1]:idx["omega", 2], 5] = fixed.pars$omega
			pars[idx["omega", 1]:idx["omega", 2], 6] = fixed.pars$omega
		}
	}
	if(modelinc[10]>0){
		gpnames = paste("alpha",1:modelinc[10],sep="")
		pxd = which(is.na(pars[idx["alpha", 1]:idx["alpha", 2], 5]))
		if(length(pxd)>0) pars[(idx["alpha", 1]:idx["alpha", 2])[pxd], 5] = 0
		pxd = which(is.na(pars[idx["alpha", 1]:idx["alpha", 2], 6]))
		if(length(pxd)>0) pars[(idx["alpha", 1]:idx["alpha", 2])[pxd], 6] =  1-TinY
		pars[idx["alpha", 1]:idx["alpha", 2], 1] = rep(0.05/modelinc[10], modelinc[10])
		sp = na.omit(match(start.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gpnames[sp[i]], 1] = as.numeric(start.pars[gpnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gpnames[sp[i]], 1] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 5] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 6] = as.numeric(fixed.pars[gpnames[sp[i]]])
			}
		}
	}
	if(modelinc[11] > 0){
		gqnames = paste("beta",1:modelinc[11],sep="")
		pxd = which(is.na(pars[idx["beta", 1]:idx["beta", 2], 5]))
		if(length(pxd)>0) pars[(idx["beta", 1]:idx["beta", 2])[pxd], 5] = 0
		pxd = which(is.na(pars[idx["beta", 1]:idx["beta", 2], 6]))
		if(length(pxd)>0) pars[(idx["beta", 1]:idx["beta", 2])[pxd], 6] =  1-TinY
		pars[idx["beta", 1]:idx["beta", 2], 1] = rep(0.7/modelinc[11], modelinc[11])
		sp = na.omit(match(start.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gqnames[sp[i]], 1] = as.numeric(start.pars[gqnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gqnames[sp[i]], 1] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 5] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 6] = as.numeric(fixed.pars[gqnames[sp[i]]])
			}
		}
	}
	
	if(modelinc[12] > 0){
		gqnames = paste("gamma",1:modelinc[12],sep="")
		pxd = which(is.na(pars[idx["gamma", 1]:idx["gamma", 2], 5]))
		if(length(pxd)>0) pars[(idx["gamma", 1]:idx["gamma", 2])[pxd], 5] = -1+TinY
		pxd = which(is.na(pars[idx["gamma", 1]:idx["gamma", 2], 6]))
		if(length(pxd)>0) pars[(idx["gamma", 1]:idx["gamma", 2])[pxd], 6] =  1-TinY
		pars[idx["gamma", 1]:idx["gamma", 2], 1] = rep(0.7/modelinc[12], modelinc[12])
		sp = na.omit(match(start.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gqnames[sp[i]], 1] = as.numeric(start.pars[gqnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gqnames[sp[i]], 1] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 5] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 6] = as.numeric(fixed.pars[gqnames[sp[i]]])
			}
		}
	}
	
	if(modelinc[17]>0){
		vxnames = paste("vxreg",1:modelinc[17],sep="")
		pxd = which(is.na(pars[idx["vxreg", 1]:idx["vxreg", 2], 5]))
		if(length(pxd)>0) pars[(idx["vxreg", 1]:idx["vxreg", 2])[pxd], 5] = 0
		pxd = which(is.na(pars[idx["vxreg", 1]:idx["vxreg", 2], 6]))
		if(length(pxd)>0) pars[(idx["vxreg", 1]:idx["vxreg", 2])[pxd], 6] = 100		
		pars[idx["vxreg", 1]:idx["vxreg", 2], 1] = rep(TinY, modelinc[17])
		sp = na.omit(match(start.names, vxnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[vxnames[sp[i]], 1] = as.numeric(start.pars[vxnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, vxnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[vxnames[sp[i]], 1] = as.numeric(fixed.pars[vxnames[sp[i]]])
				pars[vxnames[sp[i]], 5] = as.numeric(fixed.pars[vxnames[sp[i]]])
				pars[vxnames[sp[i]], 6] = as.numeric(fixed.pars[vxnames[sp[i]]])
			}
		}
	}
	tmp = .distributionstart(pars, model, arglist)
	pars = tmp$pars
	arglist = tmp$arglist
	return( list(pars  = pars, arglist = arglist))
}

.csacdstart = function(pars, arglist)
{
	model = arglist$model
	dscale = arglist$dscale
	modelinc = model$modelinc
	start.pars = model$start.pars
	fixed.pars = model$fixed.pars
	idx = model$pidx
	fixed.names = names(fixed.pars)
	start.names = names(start.pars)
	pars = .meqstart(pars, arglist)
	shape0 = arglist$shape0
	skew0 = arglist$skew0
	data  = arglist$data
	
	if(modelinc[9]>0){
		if(is.na(pars[idx["omega", 1]:idx["omega", 2], 5])) pars[idx["omega", 1]:idx["omega", 2], 5] = 1e-12
		if(is.na(pars[idx["omega", 1]:idx["omega", 2], 6])) pars[idx["omega", 1]:idx["omega", 2], 6] = 1
		if(is.null(start.pars$omega)) pars[idx["omega", 1]:idx["omega", 2], 1] = (var(data, na.rm = TRUE))/1000 else 
			pars[idx["omega", 1]:idx["omega", 2], 1] = start.pars$omega[1]/dscale
		if(any(substr(fixed.names, 1, 5) == "omega")){
			pars[idx["omega", 1]:idx["omega", 2], 1] = as.numeric(fixed.pars$omega)
			pars[idx["omega", 1]:idx["omega", 2], 5] = fixed.pars$omega
			pars[idx["omega", 1]:idx["omega", 2], 6] = fixed.pars$omega
		}
	}
	if(modelinc[10]>0){
		gpnames = paste("alpha",1:modelinc[10],sep="")
		pxd = which(is.na(pars[idx["alpha", 1]:idx["alpha", 2], 5]))
		if(length(pxd)>0) pars[(idx["alpha", 1]:idx["alpha", 2])[pxd], 5] = 0
		pxd = which(is.na(pars[idx["alpha", 1]:idx["alpha", 2], 6]))
		if(length(pxd)>0) pars[(idx["alpha", 1]:idx["alpha", 2])[pxd], 6] =  1-TinY
		pars[idx["alpha", 1]:idx["alpha", 2], 1] = rep(0.05/modelinc[10], modelinc[10])
		sp = na.omit(match(start.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gpnames[sp[i]], 1] = as.numeric(start.pars[gpnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gpnames[sp[i]], 1] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 5] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 6] = as.numeric(fixed.pars[gpnames[sp[i]]])
			}
		}
	}
	if(modelinc[11] > 0){
		gqnames = paste("beta",1:modelinc[11],sep="")
		pxd = which(is.na(pars[idx["beta", 1]:idx["beta", 2], 5]))
		if(length(pxd)>0) pars[(idx["beta", 1]:idx["beta", 2])[pxd], 5] = 0
		pxd = which(is.na(pars[idx["beta", 1]:idx["beta", 2], 6]))
		if(length(pxd)>0) pars[(idx["beta", 1]:idx["beta", 2])[pxd], 6] =  1-TinY
		pars[idx["beta", 1]:idx["beta", 2], 1] = rep(0.7/modelinc[11], modelinc[11])
		sp = na.omit(match(start.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gqnames[sp[i]], 1] = as.numeric(start.pars[gqnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gqnames[sp[i]], 1] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 5] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 6] = as.numeric(fixed.pars[gqnames[sp[i]]])
			}
		}
	}
	# \rho in the paper notation
	if(modelinc[13] > 0){
		gqnames = "eta11"
		if(is.na(pars[idx["eta1", 1]:idx["eta1", 2], 5])) pars[idx["eta1", 1]:idx["eta1", 2], 5] = TinY
		if(is.na(pars[idx["eta1", 1]:idx["eta1", 2], 6])) pars[idx["eta1", 1]:idx["eta1", 2], 6] = 1-TinY
		pars[idx["eta1", 1]:idx["eta1", 2], 1] = 0.98
		if(any(substr(start.names, 1, 4) == "eta1")){
			j = which(substr(start.names, 1, 4) == "eta1")
			gqmatch = charmatch(start.names[j],gqnames)
			pars[gqnames[gqmatch], 1] = as.numeric(start.pars[j])
		}
		if(any(substr(fixed.names, 1, 4) == "eta1")){
			j = which(substr(fixed.names, 1, 4) == "eta1")
			gqmatch = charmatch(fixed.names[j],gqnames)
			pars[gqnames[gqmatch], 1] = as.numeric(fixed.pars[j])
			pars[gqnames[gqmatch], 5] = as.numeric(fixed.pars[j])
			pars[gqnames[gqmatch], 6] = as.numeric(fixed.pars[j])
		}
	}
	# \phi in the paper notation
	if(modelinc[14] > 0){
		gqnames = "eta21"
		if(is.na(pars[idx["eta2", 1]:idx["eta2", 2], 5])) pars[idx["eta2", 1]:idx["eta2", 2], 5] = TinY
		if(is.na(pars[idx["eta2", 1]:idx["eta2", 2], 6])) pars[idx["eta2", 1]:idx["eta2", 2], 6] = 1-TinY
		pars[idx["eta2", 1]:idx["eta2", 2], 1] = 0.05
		if(any(substr(start.names, 1, 4) == "eta2")){
			j = which(substr(start.names, 1, 4) == "eta2")
			gqmatch = charmatch(start.names[j],gqnames)
			pars[gqnames[gqmatch], 1] = as.numeric(start.pars[j])
		}
		if(any(substr(fixed.names, 1, 4) == "eta2")){
			j = which(substr(fixed.names, 1, 4) == "eta2")
			gqmatch = charmatch(fixed.names[j],gqnames)
			pars[gqnames[gqmatch], 1] = as.numeric(fixed.pars[j])
			pars[gqnames[gqmatch], 5] = as.numeric(fixed.pars[j])
			pars[gqnames[gqmatch], 6] = as.numeric(fixed.pars[j])
		}
	}
	if(modelinc[17]>0){
		vxnames = paste("vxreg",1:modelinc[17],sep="")
		pxd = which(is.na(pars[idx["vxreg", 1]:idx["vxreg", 2], 5]))
		if(length(pxd)>0) pars[(idx["vxreg", 1]:idx["vxreg", 2])[pxd], 5] = 0
		pxd = which(is.na(pars[idx["vxreg", 1]:idx["vxreg", 2], 6]))
		if(length(pxd)>0) pars[(idx["vxreg", 1]:idx["vxreg", 2])[pxd], 6] = 100		
		pars[idx["vxreg", 1]:idx["vxreg", 2], 1] = rep(TinY, modelinc[17])
		sp = na.omit(match(start.names, vxnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[vxnames[sp[i]], 1] = as.numeric(start.pars[vxnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, vxnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[vxnames[sp[i]], 1] = as.numeric(fixed.pars[vxnames[sp[i]]])
				pars[vxnames[sp[i]], 5] = as.numeric(fixed.pars[vxnames[sp[i]]])
				pars[vxnames[sp[i]], 6] = as.numeric(fixed.pars[vxnames[sp[i]]])
			}
		}
	}
	tmp = .distributionstart(pars, model, arglist)
	pars = tmp$pars
	arglist = tmp$arglist
	return( list(pars  = pars, arglist = arglist))
}


# multiplicative component apARCH model start parameters
.mcsacdstart = function(pars, arglist)
{
	eps = 1e-12
	data = arglist$data
	model = arglist$model
	dscale = arglist$dscale
	modelinc = model$modelinc
	start.pars = model$start.pars
	fixed.pars = model$fixed.pars
	idx = model$pidx
	fixed.names = names(fixed.pars)
	start.names = names(start.pars)
	pars = .meqstart(pars, arglist)
	shape0 = arglist$shape0
	skew0 = arglist$skew0
	data  = arglist$data
	
	if(modelinc[9]>0){
		if(is.na(pars[idx["omega", 1]:idx["omega", 2], 5])) pars[idx["omega", 1]:idx["omega", 2], 5] = eps
		if(is.na(pars[idx["omega", 1]:idx["omega", 2], 6])) pars[idx["omega", 1]:idx["omega", 2], 6] = 5
		if(is.null(start.pars$omega)) pars[idx["omega", 1]:idx["omega", 2], 1] = 0.05 else 
			pars[idx["omega", 1]:idx["omega", 2], 1] = start.pars$omega[1]/dscale
		if(any(substr(fixed.names, 1, 5) == "omega")){
			pars[idx["omega", 1]:idx["omega", 2], 1] = as.numeric(fixed.pars$omega)
			pars[idx["omega", 1]:idx["omega", 2], 5] = fixed.pars$omega
			pars[idx["omega", 1]:idx["omega", 2], 6] = fixed.pars$omega
		}
	}
	if(modelinc[10]>0){
		gpnames = paste("alpha",1:modelinc[10],sep="")
		pxd = which(is.na(pars[idx["alpha", 1]:idx["alpha", 2], 5]))
		if(length(pxd)>0) pars[(idx["alpha", 1]:idx["alpha", 2])[pxd], 5] = 0
		pxd = which(is.na(pars[idx["alpha", 1]:idx["alpha", 2], 6]))
		if(length(pxd)>0) pars[(idx["alpha", 1]:idx["alpha", 2])[pxd], 6] =  1-TinY
		pars[idx["alpha", 1]:idx["alpha", 2], 1] = rep(0.05/modelinc[10], modelinc[10])
		sp = na.omit(match(start.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gpnames[sp[i]], 1] = as.numeric(start.pars[gpnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gpnames[sp[i]], 1] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 5] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 6] = as.numeric(fixed.pars[gpnames[sp[i]]])
			}
		}
	}
	
	if(modelinc[11] > 0){
		gqnames = paste("beta",1:modelinc[11],sep="")
		pxd = which(is.na(pars[idx["beta", 1]:idx["beta", 2], 5]))
		if(length(pxd)>0) pars[(idx["beta", 1]:idx["beta", 2])[pxd], 5] = 0
		pxd = which(is.na(pars[idx["beta", 1]:idx["beta", 2], 6]))
		if(length(pxd)>0) pars[(idx["beta", 1]:idx["beta", 2])[pxd], 6] =  1-TinY
		pars[idx["beta", 1]:idx["beta", 2], 1] = rep(0.9/modelinc[11], modelinc[11])
		sp = na.omit(match(start.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gqnames[sp[i]], 1] = as.numeric(start.pars[gqnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gqnames[sp[i]], 1] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 5] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 6] = as.numeric(fixed.pars[gqnames[sp[i]]])
			}
		}
	}
	if(modelinc[17]>0){
		vxnames = paste("vxreg",1:modelinc[17],sep="")
		pxd = which(is.na(pars[idx["vxreg", 1]:idx["vxreg", 2], 5]))
		if(length(pxd)>0) pars[(idx["vxreg", 1]:idx["vxreg", 2])[pxd], 5] = 0
		pxd = which(is.na(pars[idx["vxreg", 1]:idx["vxreg", 2], 6]))
		if(length(pxd)>0) pars[(idx["vxreg", 1]:idx["vxreg", 2])[pxd], 6] = 100
		pars[idx["vxreg", 1]:idx["vxreg", 2], 1] = rep(TinY, modelinc[17])
		sp = na.omit(match(start.names, vxnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[vxnames[sp[i]], 1] = as.numeric(start.pars[vxnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, vxnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[vxnames[sp[i]], 1] = as.numeric(fixed.pars[vxnames[sp[i]]])
				pars[vxnames[sp[i]], 5] = as.numeric(fixed.pars[vxnames[sp[i]]])
				pars[vxnames[sp[i]], 6] = as.numeric(fixed.pars[vxnames[sp[i]]])
			}
		}
	}
	tmp = .distributionstart(pars, model, arglist)
	pars = tmp$pars
	arglist = tmp$arglist
	return( list(pars  = pars, arglist = arglist))
}
.distributionstart = function(pars, model, arglist){
	skew0 = arglist$skew0
	shape0 = arglist$shape0 
	model = arglist$model
	modelinc = model$modelinc
	start.pars = model$start.pars
	fixed.pars = model$fixed.pars
	fixed.names = names(fixed.pars)
	start.names = names(start.pars)
	idx = model$pidx
	garchenv = arglist$garchenv
	assign("garchLL", NA, envir = garchenv)
	
	dbounds = arglist$sbounds
	# to be used for starting the recursion in case of problems
	midskew  = ifelse(abs( (dbounds[1] + dbounds[2])/2 ) < 1e-4, (dbounds[1] + dbounds[2]/2)/2, (dbounds[1] + dbounds[2])/2)
	midshape = ifelse(abs( (dbounds[3] + dbounds[4])/2 ) < 1e-4, (dbounds[3] + dbounds[4]/2)/2, (dbounds[3] + dbounds[4])/2)
	
	if(modelinc[18]>0){
		if(is.na(pars[idx["skew", 1]:idx["skew", 2], 5])) pars[idx["skew", 1]:idx["skew", 2], 5] = dbounds[1]
		if(is.na(pars[idx["skew", 1]:idx["skew", 2], 6])) pars[idx["skew", 1]:idx["skew", 2], 6] = dbounds[2]		
		if(is.null(start.pars$skew)) pars[idx["skew", 1]:idx["skew", 2], 1] = midskew else pars[idx["skew", 1]:idx["skew", 2], 1] = start.pars$skew[1]
		if(any(substr(fixed.names, 1, 4) == "skew")){
			pars[idx["skew", 1]:idx["skew", 2], 1] = as.numeric(fixed.pars$skew)
			pars[idx["skew", 1]:idx["skew", 2], 5] = as.numeric(fixed.pars$skew)
			pars[idx["skew", 1]:idx["skew", 2], 6] = as.numeric(fixed.pars$skew)
		}
	}
	if(modelinc[19]>0){
		if(is.na(pars[idx["shape", 1]:idx["shape", 2], 5])) pars[idx["shape", 1]:idx["shape", 2], 5] = dbounds[3]
		if(is.na(pars[idx["shape", 1]:idx["shape", 2], 6])) pars[idx["shape", 1]:idx["shape", 2], 6] = dbounds[4]		
		if(is.null(start.pars$shape)) pars[idx["shape", 1]:idx["shape", 2], 1] = midshape else pars[idx["shape", 1]:idx["shape", 2], 1] = start.pars$shape[1]
		if(any(substr(fixed.names, 1, 5) == "shape")){
			pars[idx["shape", 1]:idx["shape", 2], 1] = as.numeric(fixed.pars$shape)
			pars[idx["shape", 1]:idx["shape", 2], 5] = as.numeric(fixed.pars$shape)
			pars[idx["shape", 1]:idx["shape", 2], 6] = as.numeric(fixed.pars$shape)
		}
	}
	if(modelinc[20]>0){
		if(is.na(pars[idx["ghlambda", 1]:idx["ghlambda", 2], 5])) pars[idx["ghlambda", 1]:idx["ghlambda", 2], 5] = -4
		if(is.na(pars[idx["ghlambda", 1]:idx["ghlambda", 2], 6])) pars[idx["ghlambda", 1]:idx["ghlambda", 2], 6] =  4
		if(is.null(start.pars$ghlambda)) pars[idx["ghlambda", 1]:idx["ghlambda", 2], 1] = -0.5 else pars[idx["ghlambda", 1]:idx["ghlambda", 2], 1] = start.pars$ghlambda[1]
		if(any(substr(fixed.names, 1, 8) == "ghlambda")){
			pars[idx["ghlambda", 1]:idx["ghlambda", 2], 1] = as.numeric(fixed.pars$ghlambda)
			pars[idx["ghlambda", 1]:idx["ghlambda", 2], 5] = as.numeric(fixed.pars$ghlambda)
			pars[idx["ghlambda", 1]:idx["ghlambda", 2], 6] = as.numeric(fixed.pars$ghlambda)
		}
	}
	
	if(modelinc[21]>0 && is.null(skew0)){
		gspex = ugarchspec(mean.model = list(armaOrder = model$modelinc[2:3], include.mean = as.logical(model$modelinc[1])),
				variance.model = list(model = model$vmodel$model, garchOrder = model$modelinc[10:11], 
						variance.targeting = ifelse(model$modelinc[9]==0, TRUE, FALSE)),
				distribution.model = model$dmodel$model, fixed.pars = model$fixed.pars)
		if(model$vmodel$model=="mcsGARCH"){
			gfit = try(ugarchfit(gspex, arglist$data, solver = "hybrid", solver.control = list(trace=0), DailyVar = arglist$DailyVar), silent = TRUE)
		} else{
			gfit = try(ugarchfit(gspex, arglist$data*arglist$dscale, fit.control=list(scale=1), solver = "hybrid", solver.control = list(trace=0)), silent = TRUE)
		}
		if(inherits(gfit, "try-error")){
			garchLL = 1e10
			assign("garchLL", garchLL, envir = garchenv)
			uncskew = midskew
			uncshape = ifelse(is.null(shape0), midshape, shape0)
		} else{
			if(gfit@fit$convergence!=0){
				garchLL = 1e10
				assign("garchLL", garchLL, envir = garchenv)
				uncskew = midskew
				uncshape = ifelse(is.null(shape0), midshape, shape0)
			} else{
				garchLL = likelihood(gfit)
				assign("garchLL", garchLL, envir = garchenv)
				uncskew = coef(gfit)["skew"]
				uncshape = ifelse(is.null(shape0), coef(gfit)["shape"], shape0)
			}
		}
		skew0 = uncskew
		shape0 = uncshape
	}
	if(modelinc[27]>0 && is.null(shape0)){
		gspex = ugarchspec(mean.model = list(armaOrder = model$modelinc[2:3], include.mean = as.logical(model$modelinc[1])),
				variance.model = list(model = model$vmodel$model, garchOrder = model$modelinc[10:11], 
						variance.targeting = ifelse(model$modelinc[9]==0, TRUE, FALSE)),
				distribution.model = model$dmodel$model, fixed.pars = model$fixed.pars)
		if(model$vmodel$model=="mcsGARCH"){
			gfit = try(ugarchfit(gspex, arglist$data, solver = "hybrid", solver.control = list(trace=0), DailyVar = arglist$DailyVar), silent = TRUE)
		} else{
			gfit = try(ugarchfit(gspex, arglist$data*arglist$dscale, fit.control=list(scale=1), solver = "hybrid", solver.control = list(trace=0)), silent = TRUE)
		}
		if(inherits(gfit, "try-error")){
			garchLL = 1e10
			assign("garchLL", garchLL, envir = garchenv)
			uncskew = ifelse(is.null(skew0), midskew, skew0)
			uncshape = midshape
		} else{
			if(gfit@fit$convergence!=0){
				garchLL = 1e10
				assign("garchLL", garchLL, envir = garchenv)
				uncskew = ifelse(is.null(skew0), midskew, skew0)
				uncshape = midshape
			} else{
				garchLL = likelihood(gfit)
				assign("garchLL", garchLL, envir = garchenv)
				uncskew = ifelse(is.null(skew0), coef(gfit)["skew"], skew0)
				uncshape = coef(gfit)["shape"]
			}
		}
		skew0 = uncskew
		shape0 = uncshape
	}
	uncskew = skew0
	uncshape = shape0
	arglist$skhEst = c(0, 0)
	if(modelinc[21]>0){
		xskew = .acdskewbounds(modelinc[22:24], uncskew, model$dmodel$model, dbounds[1:2])
		arglist$skhEst[1] = xskew$sk0
		if(is.na(pars[idx["skcons", 1]:idx["skcons", 2], 5])) pars[idx["skcons", 1]:idx["skcons", 2], 5] = xskew$skewpar.LB[1]
		if(is.na(pars[idx["skcons", 1]:idx["skcons", 2], 6])) pars[idx["skcons", 1]:idx["skcons", 2], 6] = xskew$skewpar.UB[1]
		if(is.null(start.pars$skcons)) pars[idx["skcons", 1]:idx["skcons", 2], 1] = xskew$skewpars[1] else pars[idx["skcons", 1]:idx["skcons", 2], 1] = start.pars$skcons[1]
		if(any(substr(fixed.names, 1, 6) == "skcons")){
			pars[idx["skcons", 1]:idx["skcons", 2], 1] = as.numeric(fixed.pars$skcons)
			pars[idx["skcons", 1]:idx["skcons", 2], 5] = as.numeric(fixed.pars$skcons)
			pars[idx["skcons", 1]:idx["skcons", 2], 6] = as.numeric(fixed.pars$skcons)
		}
		if(modelinc[22]>0){
			sknames = paste("skalpha",1:modelinc[22],sep="")
			pxd = which(is.na(pars[idx["skalpha", 1]:idx["skalpha", 2], 5]))
			if(length(pxd)>0) pars[(idx["skalpha", 1]:idx["skalpha", 2])[pxd], 5] = (xskew$skewpar.LB[2:(2+modelinc[22]-1)])[pxd]
			pxd = which(is.na(pars[idx["skalpha", 1]:idx["skalpha", 2], 6]))
			if(length(pxd)>0) pars[(idx["skalpha", 1]:idx["skalpha", 2])[pxd], 6] = (xskew$skewpar.UB[2:(2+modelinc[22]-1)])[pxd]
			pars[idx["skalpha", 1]:idx["skalpha", 2], 1] = xskew$skewpars[2:(2+modelinc[22]-1)]
			if(any(substr(start.names, 1, 7) == "skalpha")){
				j = which(substr(start.names, 1, 7) == "skalpha")
				skmatch = charmatch(start.names[j],sknames)
				pars[sknames[skmatch], 1] = as.numeric(start.pars[j])
			}
			if(any(substr(fixed.names, 1, 7) == "skalpha")){
				j = which(substr(fixed.names, 1, 7) == "skalpha")
				skmatch = charmatch(fixed.names[j],sknames)
				pars[sknames[skmatch], 1] = as.numeric(fixed.pars[j])
				pars[sknames[skmatch], 5] = as.numeric(fixed.pars[j])
				pars[sknames[skmatch], 6] = as.numeric(fixed.pars[j])
			}
		}
		if(modelinc[23]>0){
			sknames = paste("skgamma",1:modelinc[23],sep="")
			
			pxd = which(is.na(pars[idx["skgamma", 1]:idx["skgamma", 2], 5]))
			if(length(pxd)>0) pars[(idx["skgamma", 1]:idx["skgamma", 2])[pxd], 5] = (xskew$skewpar.LB[(1+modelinc[22]+1):(1+modelinc[22]+modelinc[23])])[pxd]
			pxd = which(is.na(pars[idx["skgamma", 1]:idx["skgamma", 2], 6]))
			if(length(pxd)>0) pars[(idx["skgamma", 1]:idx["skgamma", 2])[pxd], 6] = (xskew$skewpar.UB[(1+modelinc[22]+1):(1+modelinc[22]+modelinc[23])])[pxd]
			pars[idx["skgamma", 1]:idx["skgamma", 2], 1] = xskew$skewpars[(1+modelinc[22]+1):(1+modelinc[22]+modelinc[23])]
			if(any(substr(start.names, 1, 7) == "skgamma")){
				j = which(substr(start.names, 1, 7) == "skgamma")
				skmatch = charmatch(start.names[j],sknames)
				pars[sknames[skmatch], 1] = as.numeric(start.pars[j])
			}
			if(any(substr(fixed.names, 1, 7) == "skgamma")){
				j = which(substr(fixed.names, 1, 7) == "skgamma")
				skmatch = charmatch(fixed.names[j],sknames)
				pars[sknames[skmatch], 1] = as.numeric(fixed.pars[j])
				pars[sknames[skmatch], 5] = as.numeric(fixed.pars[j])
				pars[sknames[skmatch], 6] = as.numeric(fixed.pars[j])
			}
		}
		if(modelinc[24]>0){
			sknames = paste("skbeta",1:modelinc[24],sep="")
			pxd = which(is.na(pars[idx["skbeta", 1]:idx["skbeta", 2], 5]))
			if(length(pxd)>0) pars[(idx["skbeta", 1]:idx["skbeta", 2])[pxd], 5] = (xskew$skewpar.LB[(1+modelinc[22]+modelinc[23]+1):(1+modelinc[22]+modelinc[23]+modelinc[24])])[pxd]
			pxd = which(is.na(pars[idx["skbeta", 1]:idx["skbeta", 2], 6]))
			if(length(pxd)>0) pars[(idx["skbeta", 1]:idx["skbeta", 2])[pxd], 6] = (xskew$skewpar.UB[(1+modelinc[22]+modelinc[23]+1):(1+modelinc[22]+modelinc[23]+modelinc[24])])[pxd]
			pars[idx["skbeta", 1]:idx["skbeta", 2], 1] = xskew$skewpars[(1+modelinc[22]+modelinc[23]+1):(1+modelinc[22]+modelinc[23]+modelinc[24])]
			if(any(substr(start.names, 1, 6) == "skbeta")){
				j = which(substr(start.names, 1, 6) == "skbeta")
				skmatch = charmatch(start.names[j],sknames)
				pars[sknames[skmatch], 1] = as.numeric(start.pars[j])
			}
			if(any(substr(fixed.names, 1, 6) == "skbeta")){
				j = which(substr(fixed.names, 1, 6) == "skbeta")
				skmatch = charmatch(fixed.names[j],sknames)
				pars[sknames[skmatch], 1] = as.numeric(fixed.pars[j])
				pars[sknames[skmatch], 5] = as.numeric(fixed.pars[j])
				pars[sknames[skmatch], 6] = as.numeric(fixed.pars[j])
			}
		}
		if(modelinc[25]>0){
			sknames = paste("skxreg",1:modelinc[25],sep="")
			pxd = which(is.na(pars[idx["skxreg", 1]:idx["skxreg", 2], 5]))
			if(length(pxd)>0) pars[(idx["skxreg", 1]:idx["skxreg", 2])[pxd], 5] = -100
			pxd = which(is.na(pars[idx["skxreg", 1]:idx["skxreg", 2], 6]))
			if(length(pxd)>0) pars[(idx["skxreg", 1]:idx["skxreg", 2])[pxd], 6] = 100
			pars[idx["skxreg", 1]:idx["skxreg", 2], 1] =  0.1
			if(any(substr(start.names, 1, 6) == "skxreg")){
				j = which(substr(start.names, 1, 6) == "skxreg")
				skmatch = charmatch(start.names[j],sknames)
				pars[sknames[skmatch], 1] = as.numeric(start.pars[j])
			}
			if(any(substr(fixed.names, 1, 6) == "skxreg")){
				j = which(substr(fixed.names, 1, 6) == "skxreg")
				skmatch = charmatch(fixed.names[j],sknames)
				pars[sknames[skmatch], 1] = as.numeric(fixed.pars[j])
				pars[sknames[skmatch], 5] = as.numeric(fixed.pars[j])
				pars[sknames[skmatch], 6] = as.numeric(fixed.pars[j])
			}
		}
		if(modelinc[26]>0){
			if(is.na(pars[idx["thskew", 1]:idx["thskew", 2], 5])) pars[idx["thskew", 1]:idx["thskew", 2], 5] = -1.5
			if(is.na(pars[idx["thskew", 1]:idx["thskew", 2], 6])) pars[idx["thskew", 1]:idx["thskew", 2], 6] =  1.5
			if(is.null(start.pars$thskew)) pars[idx["thskew", 1]:idx["thskew", 2], 1] = 0 else pars[idx["thskew", 1]:idx["thskew", 2], 1] = start.pars$thskew[1]
			if(any(substr(fixed.names, 1, 6) == "thskew")){
				pars[idx["thskew", 1]:idx["thskew", 2], 1] = as.numeric(fixed.pars$thskew)
				pars[idx["thskew", 1]:idx["thskew", 2], 5] = as.numeric(fixed.pars$thskew)
				pars[idx["thskew", 1]:idx["thskew", 2], 6] = as.numeric(fixed.pars$thskew)
			}
		}
	}
	
	
	if(modelinc[27]>0){
		xshape = .acdshapebounds(modelinc[28:30], uncshape, model$dmodel$model, dbounds[3:5])
		arglist$skhEst[2] = xshape$sh0
		pxd = which(is.na(pars[idx["shcons", 1]:idx["shcons", 2], 5]))
		if(length(pxd)>0) pars[(idx["shcons", 1]:idx["shcons", 2])[pxd], 5] = xshape$shapepar.LB[1]
		pxd = which(is.na(pars[idx["shcons", 1]:idx["shcons", 2], 6]))
		if(length(pxd)>0) pars[(idx["shcons", 1]:idx["shcons", 2])[pxd], 6] = xshape$shapepar.UB[1]
		if(is.null(start.pars$shcons)) pars[idx["shcons", 1]:idx["shcons", 2], 1] = xshape$shapepars[1] else pars[idx["shcons", 1]:idx["shcons", 2], 1] = start.pars$shcons[1]
		if(any(substr(fixed.names, 1, 6) == "shcons")){
			pars[idx["shcons", 1]:idx["shcons", 2], 1] = as.numeric(fixed.pars$shcons)
			pars[idx["shcons", 1]:idx["shcons", 2], 5] = as.numeric(fixed.pars$shcons)
			pars[idx["shcons", 1]:idx["shcons", 2], 6] = as.numeric(fixed.pars$shcons)
		}
		
		if(modelinc[28]>0){
			shnames = paste("shalpha",1:modelinc[28],sep="")
			pxd = which(is.na(pars[idx["shalpha", 1]:idx["shalpha", 2], 5]))
			if(length(pxd)>0) pars[(idx["shalpha", 1]:idx["shalpha", 2])[pxd], 5] = (xshape$shapepar.LB[2:(1+modelinc[28])])[pxd]
			pxd = which(is.na(pars[idx["shalpha", 1]:idx["shalpha", 2], 6]))
			if(length(pxd)>0) pars[(idx["shalpha", 1]:idx["shalpha", 2])[pxd], 6] = (xshape$shapepar.UB[2:(1+modelinc[28])])[pxd]
			pars[idx["shalpha", 1]:idx["shalpha", 2], 1] = xshape$shapepars[2:(1+modelinc[28])]
			if(any(substr(start.names, 1, 7) == "shalpha")){
				j = which(substr(start.names, 1, 7) == "shalpha")
				shmatch = charmatch(start.names[j],shnames)
				pars[shnames[shmatch], 1] = as.numeric(start.pars[j])
			}
			if(any(substr(fixed.names, 1, 7) == "shalpha")){
				j = which(substr(fixed.names, 1, 7) == "shalpha")
				shmatch = charmatch(fixed.names[j],shnames)
				pars[shnames[shmatch], 1] = as.numeric(fixed.pars[j])
				pars[shnames[shmatch], 5] = as.numeric(fixed.pars[j])
				pars[shnames[shmatch], 6] = as.numeric(fixed.pars[j])
			}
		}
		if(modelinc[29]>0){
			shnames = paste("shgamma",1:modelinc[29],sep="")
			pxd = which(is.na(pars[idx["shgamma", 1]:idx["shgamma", 2], 5]))
			if(length(pxd)>0) pars[(idx["shgamma", 1]:idx["shgamma", 2])[pxd], 5] = (xshape$shapepar.LB[(1+modelinc[28]+1):(1+modelinc[28]+modelinc[29])])[pxd]
			pxd = which(is.na(pars[idx["shgamma", 1]:idx["shgamma", 2], 6]))
			if(length(pxd)>0) pars[(idx["shgamma", 1]:idx["shgamma", 2])[pxd], 6] = (xshape$shapepar.UB[(1+modelinc[28]+1):(1+modelinc[28]+modelinc[29])])[pxd]
			pars[idx["shgamma", 1]:idx["shgamma", 2], 1] = xshape$shapepars[(1+modelinc[28]+1):(1+modelinc[28]+modelinc[29])]
			if(any(substr(start.names, 1, 7) == "shgamma")){
				j = which(substr(start.names, 1, 7) == "shgamma")
				shmatch = charmatch(start.names[j],shnames)
				pars[shnames[shmatch], 1] = as.numeric(start.pars[j])
			}
			if(any(substr(fixed.names, 1, 7) == "shgamma")){
				j = which(substr(fixed.names, 1, 7) == "shgamma")
				shmatch = charmatch(fixed.names[j],shnames)
				pars[shnames[shmatch], 1] = as.numeric(fixed.pars[j])
				pars[shnames[shmatch], 5] = as.numeric(fixed.pars[j])
				pars[shnames[shmatch], 6] = as.numeric(fixed.pars[j])
			}
		}
		if(modelinc[30]>0){
			shnames = paste("shbeta",1:modelinc[30],sep="")
			pxd = which(is.na(pars[idx["shbeta", 1]:idx["shbeta", 2], 5]))
			if(length(pxd)>0) pars[(idx["shbeta", 1]:idx["shbeta", 2])[pxd], 5] = (xshape$shapepar.LB[(1+modelinc[28]+modelinc[29]+1):(1+modelinc[28]+modelinc[29]+modelinc[30])])[pxd]
			pxd = which(is.na(pars[idx["shbeta", 1]:idx["shbeta", 2], 6]))
			if(length(pxd)>0) pars[(idx["shbeta", 1]:idx["shbeta", 2])[pxd], 6] = (xshape$shapepar.UB[(1+modelinc[28]+modelinc[29]+1):(1+modelinc[28]+modelinc[29]+modelinc[30])])[pxd]
			pars[idx["shbeta", 1]:idx["shbeta", 2], 1] = xshape$shapepars[(1+modelinc[28]+modelinc[29]+1):(1+modelinc[28]+modelinc[29]+modelinc[30])]
			if(any(substr(start.names, 1, 6) == "shbeta")){
				j = which(substr(start.names, 1, 6) == "shbeta")
				shmatch = charmatch(start.names[j],shnames)
				pars[shnames[shmatch], 1] = as.numeric(start.pars[j])
			}
			if(any(substr(fixed.names, 1, 6) == "shbeta")){
				j = which(substr(fixed.names, 1, 6) == "shbeta")
				shmatch = charmatch(fixed.names[j],shnames)
				pars[shnames[shmatch], 1] = as.numeric(fixed.pars[j])
				pars[shnames[shmatch], 5] = as.numeric(fixed.pars[j])
				pars[shnames[shmatch], 6] = as.numeric(fixed.pars[j])
			}
		}
		if(modelinc[31]>0){
			shnames = paste("shxreg",1:modelinc[31],sep="")
			pxd = which(is.na(pars[idx["shxreg", 1]:idx["shxreg", 2], 5]))
			if(length(pxd)>0) pars[(idx["shxreg", 1]:idx["shxreg", 2])[pxd], 5] = 0
			pxd = which(is.na(pars[idx["shxreg", 1]:idx["shxreg", 2], 6]))
			if(length(pxd)>0) pars[(idx["shxreg", 1]:idx["shxreg", 2])[pxd], 6] = 100
			pars[idx["shxreg", 1]:idx["shxreg", 2], 1] =  0.1
			if(any(substr(start.names, 1, 6) == "shxreg")){
				j = which(substr(start.names, 1, 6) == "shxreg")
				shmatch = charmatch(start.names[j],shnames)
				pars[shnames[shmatch], 1] = as.numeric(start.pars[j])
			}
			if(any(substr(fixed.names, 1, 6) == "shxreg")){
				j = which(substr(fixed.names, 1, 6) == "shxreg")
				shmatch = charmatch(fixed.names[j],shnames)
				pars[shnames[shmatch], 1] = as.numeric(fixed.pars[j])
				pars[shnames[shmatch], 5] = as.numeric(fixed.pars[j])
				pars[shnames[shmatch], 6] = as.numeric(fixed.pars[j])
			}
		}
		if(modelinc[32]>0){
			if(is.na(pars[idx["thshape", 1]:idx["thshape", 2], 5])) pars[idx["thshape", 1]:idx["thshape", 2], 5] = -1.5
			if(is.na(pars[idx["thshape", 1]:idx["thshape", 2], 6])) pars[idx["thshape", 1]:idx["thshape", 2], 6] =  1.5
			if(is.null(start.pars$thshape)) pars[idx["thshape", 1]:idx["thskew", 2], 1] = 0 else pars[idx["thshape", 1]:idx["thshape", 2], 1] = start.pars$thshape[1]
			if(any(substr(fixed.names, 1, 7) == "thshape")){
				pars[idx["thshape", 1]:idx["thshape", 2], 1] = as.numeric(fixed.pars$thshape)
				pars[idx["thshape", 1]:idx["thshape", 2], 5] = as.numeric(fixed.pars$thshape)
				pars[idx["thshape", 1]:idx["thshape", 2], 6] = as.numeric(fixed.pars$thshape)
			}
		}
	}
	return(list(pars = pars, arglist = arglist))
}
ewmav = function(X, beta = 0.96){
	p = NROW(X)
	h = sqrt(filter((1-beta)*c(0,X[1:(p-1)]^2), filter = beta, method = "recursive", sides = 1, init = c(mean(X^2))))
	return(h)
}