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

#----------------------------------------------------------------------------------
acdspec = function(
		variance.model = list(model = "sGARCH", garchOrder = c(1,1), 
				external.regressors = NULL, variance.targeting = FALSE), 
		mean.model = list(armaOrder = c(1,1), include.mean = TRUE, archm = FALSE, 
				arfima = FALSE, external.regressors = NULL), 
		distribution.model = list(model = "snorm", 
				skewOrder = c(1, 1, 1), skewshock = 1, skewshocktype = 1, 
				skewmodel = "quad", skew.regressors = NULL,
				shapeOrder = c(0, 1, 1), shapeshock = 1, shapeshocktype = 1, 
				shapemodel = "quad", shape.regressors = NULL, exp.rate = 1), 
		start.pars = list(), fixed.pars = list())
{
	UseMethod("acdspec")
}
# shock types: 
# 1 in z^2
# 2, in resids^2
# 3 in abs(z)
# 4 in abs(resid)
.expand.model = function(model){
	modelnames = NULL
	for(i in 1:32){
		if(model[i]>0){
			if(any(c(2,3,8,10,11,12,13,14,17,22,23,24,25,28,29,30,31) == i)){
				modelnames = c(modelnames, paste(names(model)[i], 1:model[i], sep = ""))
			} else{
				modelnames = c(modelnames, names(model)[i])
			}
		}
	}
	return( modelnames )
}

.acdspec = function(
		variance.model = list(model = "sGARCH", garchOrder = c(1,1), 
				external.regressors = NULL, variance.targeting = FALSE), 
		mean.model = list(armaOrder = c(1,1), include.mean = TRUE, archm = FALSE, 
				arfima = FALSE, external.regressors = NULL), 
		distribution.model = list(model = "snorm", 
				skewOrder = c(1, 1, 1), skewshock = 1, skewshocktype = 1, 
				skewmodel = "quad", skew.regressors = NULL,
				shapeOrder = c(0, 1, 1), shapeshock = 1, shapeshocktype = 1, 
				shapemodel = "quad", shape.regressors = NULL, exp.rate=1), 
		start.pars = list(), fixed.pars = list())
{
	# not supported:
	mean.model$skewm = FALSE
	mean.model$shapem = FALSE
	# some checks and preparation to be passed on to specific models by switch
	# at present, external regressors in the skew and shape dynamics are not
	# included, but provision has been made for future expansion (skxreg and shxreg)
	#---------------------------------------------------------------------------
	modelinc = rep(0, 41)
	names(modelinc) = c("mu", "ar", "ma", "arfima", "archm", "skewm", "shapem", "mxreg", 
			"omega", "alpha", "beta", "gamma", "eta1", "eta2", "delta", "lambda", "vxreg", 
			"skew", "shape", "ghlambda", 
			"skcons", "skalpha", "skgamma", "skbeta", "skxreg", "thskew",
			"shcons", "shalpha", "shgamma", "shbeta", "shxreg", "thshape",
			"skewmodel", "shapemodel", "skshock", "shshock", 
			"maxskew", "maxshape", "aux", "aux", "aux")
	
	modeldesc = list()
	modeldata = list()
	
	#---------------------------------------------------------------------------
	mm = match(names(mean.model), c("armaOrder", "include.mean", "archm", "skewm", "shapem", "arfima", "external.regressors"))
	if(any(is.na(mm))){
		idx = which(is.na(mm))
		enx = NULL
		for(i in 1:length(idx)) enx = c(enx, names(mean.model)[idx[i]])
		warning(paste(c("unidentified option(s) in mean.model:\n", enx), sep="", collapse=" "), call. = FALSE, domain = NULL)
	}
	vm = match(names(variance.model), c("model", "garchOrder", "external.regressors", "variance.targeting"))
	if(any(is.na(vm))){
		idx = which(is.na(vm))
		enx = NULL
		for(i in 1:length(idx)) enx = c(enx, names(variance.model)[idx[i]])
		warning(paste(c("unidentified option(s) in variance.model:\n", enx), sep="", collapse=" "), call. = FALSE, domain = NULL)
	}
	dm = match(names(distribution.model), c("model", "skewOrder", "skewshock", "skewshocktype", 
					"skewmodel", "skew.regressors", "shapeOrder", "shapeshock", "shapeshocktype", 
					"shapemodel", "shape.regressors", "exp.rate"))
	if(any(is.na(dm))){
		idx = which(is.na(dm))
		enx = NULL
		for(i in 1:length(idx)) enx = c(enx, names(distribution.model)[idx[i]])
		warning(paste(c("unidentified option(s) in distribution.model:\n", enx), sep="", collapse=" "), call. = FALSE, domain = NULL)
	}
	
	#---------------------------------------------------------------------------
	vmodel = list(model = "sGARCH", garchOrder = c(1,1), external.regressors = NULL, variance.targeting = FALSE)
	idx = na.omit(match(names(variance.model), names(vmodel)))
	if(length(idx)>0) for(i in 1:length(idx)) vmodel[idx[i]] = variance.model[i]
	valid.model = c("sGARCH", "csGARCH", "mcsGARCH")
	if(!any(vmodel$model == valid.model)) 
		stop("\nacdpec-->error: the garch model does not appear to be a valid choice.\n", call. = FALSE)
	modelinc[10] = vmodel$garchOrder[1]
	modelinc[11] = vmodel$garchOrder[2]
	if( vmodel$model == "csGARCH" ){
		modelinc[14] = modelinc[13] = 1
		vmodel$variance.targeting = FALSE
	}
	if(!is.null(vmodel$external.regressors) && !is.matrix(vmodel$external.regressors))
		stop("\nacdspec-->error: external.regressors (variance) must be a matrix.\n", call. = FALSE)
	modeldata$vexdata = vmodel$external.regressors
	if(!is.null(vmodel$external.regressors)) modelinc[17] = NCOL( vmodel$external.regressors )
	if( is.null(vmodel$variance.targeting) ) modelinc[9] = 1 else modelinc[9] = as.integer( 1-vmodel$variance.targeting )
	
	#---------------------------------------------------------------------------
	mmodel = list(armaOrder = c(1,1), include.mean = TRUE, archm = FALSE, skewm = FALSE, 
			shapem = FALSE, arfima = FALSE, external.regressors = NULL)
	idx = na.omit(match(names(mean.model), names(mmodel)))
	if(length(idx)>0) for(i in 1:length(idx)) mmodel[idx[i]] = mean.model[i]
	
	modelinc[2] = mmodel$armaOrder[1]
	modelinc[3] = mmodel$armaOrder[2]
	if(is.null(mmodel$include.mean)) modelinc[1] = 1 else modelinc[1] = as.integer( mmodel$include.mean )
	
	if( as.logical(mmodel$archm) ){
		modelinc[5] = 1
	}
	modelinc[4] = as.integer( as.logical(mmodel$arfima) )
	
	if(!is.null(mmodel$external.regressors) && !is.matrix(mmodel$external.regressors))
		stop("\nacdspec-->error: external.regressors (mean) must be a matrix.\n", call. = FALSE)
	modeldata$mexdata = mmodel$external.regressors
	if( !is.null(mmodel$external.regressors) ) modelinc[8] = NCOL( mmodel$external.regressors )
	#---------------------------------------------------------------------------
	dmodel =  list(model = "snorm", skewOrder = c(1, 1, 1), skewshock = 1, skewshocktype = 1, 
			skewmodel = "quad", skew.regressors = NULL, shapeOrder = c(1, 1, 1), 
			shapeshock = 1, shapeshocktype = 1, shapemodel = "quad", shape.regressors = NULL, 
			exp.rate=1)
	idx = na.omit(match(names(distribution.model), names(dmodel)))
	if(length(idx)>0) for(i in 1:length(idx)) dmodel[idx[i]] = distribution.model[i]
	
	valid.distribution = c("snorm", "std", "sstd", "ged", "sged", "nig", "ghyp", "jsu", "ghst")
	if(!any(dmodel$model==valid.distribution)) 
		stop("\nacdspec-->error: the cond.distribution does not appear to be a valid choice.")
	
	if(dmodel$skewmodel  == "pwl"){
		# skewshocktype not relevant for these type of dynamics
		modelinc[33] = 1
	} else if(dmodel$skewmodel == "xar"){
		if(dmodel$skewshocktype==1) modelinc[33] = 2 else modelinc[33] = 4
		dmodel$skewshock = 2
	} else if(dmodel$skewmodel == "tar"){
		modelinc[33] = 5
		modelinc[26] = 1
	} else{
		# default (quad)
		if(dmodel$skewshocktype==1) modelinc[33] = 0 else modelinc[33] = 3
	}
	
	if(dmodel$shapemodel == "pwl"){
		if(dmodel$shapeshocktype==1) modelinc[34] = 1 else modelinc[34] = 5
	} else if(dmodel$shapemodel == "xar"){
		if(dmodel$shapeshocktype==1) modelinc[34] = 2 else modelinc[34] = 6		
		distribution.model$shapeOrder[1] = 0
		dmodel$shapeshock = 2
	} else if(dmodel$shapemodel == "tar"){
		if(dmodel$shapeshocktype==1) modelinc[34] = 3 else modelinc[34] = 7
		modelinc[32] = 1
	} else{
		if(dmodel$shapeshocktype==1) modelinc[34] = 0 else modelinc[34] = 4
	}
	
	if(dmodel$skewshock  == 1) modelinc[35] = 1
	if(dmodel$shapeshock == 1) modelinc[36] = 1
	
	# because we exlucde the normal, we add 1 to the value (for c code)
	modeldesc$distno = which(dmodel$model == valid.distribution)+1
	di = .DistributionBounds(dmodel$model)
	sbounds = rep(0, 5)
	skmax = shmax = vmax = armax = 0
	# check if time-varying first
	if(is.null(dmodel$skewOrder)){
		modelinc[21:24] = 0
		modelinc[18] = di$include.skew
		if(modelinc[18]>0){
			sbounds[1] = di$skew.LB
			sbounds[2] = di$skew.UB
		}
	} else{
		if(di$include.skew){
			if(sum(dmodel$skewOrder)>0){
				modelinc[21] = 1
				modelinc[22] = dmodel$skewOrder[1]
				modelinc[23] = dmodel$skewOrder[2]
				modelinc[24] = dmodel$skewOrder[3]
				modelinc[18] = 0
				modelinc[37] = skmax = max(dmodel$skewOrder[1:3])
				# if time varying skewness, check for skew-in-mean use
				if(as.logical(mmodel$skewm)) modelinc[6] = as.integer(as.logical(mmodel$skewm))
				
				if(!is.null(dmodel$skew.regressors) && !is.matrix(dmodel$skew.regressors))
					stop("\nacdspec-->error: skew.regressors must be a matrix.\n", call. = FALSE)
				modeldata$skxdata = dmodel$skew.regressors
				if(!is.null(dmodel$skew.regressors)) modelinc[25] = NCOL( dmodel$skew.regressors )
			} else{
				modelinc[21:24] = 0
				modelinc[18] = di$include.skew
				modelinc[26] = 0
			}
			sbounds[1] = di$skew.LB
			sbounds[2] = di$skew.UB
		} else{
			modelinc[21:24] = 0
			modelinc[18] = di$include.skew
		}
	}
	
	if(is.null(dmodel$shapeOrder)){
		modelinc[27:30] = 0
		modelinc[19] = di$include.shape
		if(modelinc[19]>0){
			sbounds[3] = di$shape.LB
			sbounds[4] = di$shape.UB
		}
	} else{
		if(di$include.shape){
			if(sum(dmodel$shapeOrder)>0){
				modelinc[27] = 1
				modelinc[28] = distribution.model$shapeOrder[1]
				modelinc[29] = distribution.model$shapeOrder[2]
				modelinc[30] = distribution.model$shapeOrder[3]
				modelinc[19] = 0
				modelinc[38] = shmax = max(dmodel$shapeOrder[1:3])
				# if time varying shape, check for shape-in-mean use
				if(as.logical(mmodel$shapem)) modelinc[7] = as.integer(as.logical(mmodel$shapem))
				if(!is.null(dmodel$shape.regressors) && !is.matrix(dmodel$shape.regressors))
					stop("\nacdspec-->error: shape.regressors must be a matrix.\n", call. = FALSE)
				modeldata$shxdata = dmodel$shape.regressors
				if(!is.null(dmodel$shape.regressors)) modelinc[31] = NCOL( dmodel$shape.regressors )
			} else{
				modelinc[27:30] = 0
				modelinc[19] = di$include.shape
				modelinc[32] = 0
			}
			sbounds[3] = di$shape.LB
			sbounds[4] = di$shape.UB
		} else{
			modelinc[27:30] = 0
			modelinc[19] = di$include.shape
		}
	}
	modelinc[20] = di$include.ghlambda
	# the last aux value is the distribution number
	modelinc[41] = modeldesc$distno
	
	
	maxOrder = max(c(max(mmodel$armaOrder), max(vmodel$garchOrder), max(dmodel$skewOrder), max(dmodel$shapeOrder)))
	modelnames = .expand.model(modelinc)
	
	pos = 1
	pos.matrix = matrix(0, ncol = 3, nrow = 32)
	colnames(pos.matrix) = c("start", "stop", "include")
	rownames(pos.matrix) = c("mu", "ar", "ma", "arfima", "archm", "skewm", "shapem", "mxreg", 
			"omega", "alpha", "beta", "gamma", "eta1", "eta2", "delta", "lambda", "vxreg", 
			"skew", "shape", "ghlambda", 
			"skcons", "skalpha", "skgamma", "skbeta", "skxreg", "thskew",
			"shcons", "shalpha", "shgamma", "shbeta", "shxreg", "thshape")
	for(i in 1:32){
		if( modelinc[i] > 0 ){
			pos.matrix[i,1:3] = c(pos, pos+modelinc[i]-1, 1)
			pos = max(pos.matrix[1:i,2]+1)
		}
	}
	nn = length(modelnames)
	modelmatrix = matrix(0, ncol = 3, nrow = nn)
	rownames(modelmatrix) = modelnames
	colnames(modelmatrix) = c("opt", "fixed", "start")
	fixed.names = names(fixed.pars)
	fp = charmatch(fixed.names, modelnames)
	
	if(!is.null(fixed.names) && any(!is.na(fp))){
		fixed = fp[!is.na(fp)]
		modelmatrix[fixed,2] = 1
		fz = charmatch(modelnames, fixed.names)
		fz = fz[!is.na(fz)]
		fixed.pars = fixed.pars[fz]
		names(fixed.pars) = fixed.names[fz]
	} else{
		fixed.pars = NULL
	}
	modelmatrix[,1] = 1 - modelmatrix[,2]
	start.names = names(start.pars)
	sp = charmatch(start.names, modelnames)
	if(!is.null(start.names) && any(!is.na(sp))){
		start = sp[!is.na(sp)]
		modelmatrix[start,3] = 1
		sz = charmatch(modelnames, start.names)
		sz = sz[!is.na(sz)]
		start.pars = start.pars[sz]
	} else{
		start.pars = NULL
	}
	
	
	##################################################################
	# Parameter Matrix
	mm = sum(modelinc[c(2,3,8,10,11,12,13,14,17,22,23,24,25,27,28,29,30,31,32)])
	mm = mm - length( which(modelinc[c(2,3,8,10,11,12,13,14,17,22,23,24,25,27,28,29,30,31,32)]>0) )
	pars = matrix(0, ncol = 6, nrow = 32 + mm)
	colnames(pars) = c("Level", "Fixed", "Include", "Estimate", "LB", "UB")
	pidx = matrix(NA, nrow = 32, ncol = 2)
	colnames(pidx) = c("begin", "end")
	rownames(pidx) =  c("mu", "ar", "ma", "arfima", "archm", "skewm", "shapem", "mxreg", 
			"omega", "alpha", "beta", "gamma", "eta1", "eta2", "delta", "lambda", "vxreg", 
			"skew", "shape", "ghlambda", 
			"skcons", "skalpha", "skgamma", "skbeta", "skxreg", "thskew",
			"shcons", "shalpha", "shgamma", "shbeta", "shxreg", "thshape")
	fixed.names = names(fixed.pars)
	pnames = NULL
	nx = 0
	if(pos.matrix[1,3]==1){
		pars[1, 3] = 1
		pars[1, 1] = 0
		if(any(substr(fixed.names, 1, 2)=="mu")) pars[1,2] = 1 else pars[1,4] = 1
	}
	pidx[1,1] = 1
	pidx[1,2] = 1
	pnames = c(pnames, "mu")
	nx = 1
	pn = 1
	pidx[2,1] = 2
	if(pos.matrix[2,3] == 1){
		pn = length( seq(pos.matrix[2,1], pos.matrix[2,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			nnx = paste("ar", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "ar")
	}
	pidx[2,2] = 1+pn
	
	nx = nx + pn
	pn = 1
	pidx[3,1] = nx+1
	if(pos.matrix[3,3] == 1){
		pn = length( seq(pos.matrix[3,1], pos.matrix[3,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			nnx = paste("ma", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "ma")
	}
	pidx[3,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[4,1] = nx+1
	if(pos.matrix[4,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		if(any(!is.na(match(fixed.names, "arfima")))) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	}
	pnames = c(pnames, "arfima")
	pidx[4,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[5,1] = nx+1
	if(pos.matrix[5,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		if(any(!is.na(match(fixed.names, "archm")))) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	}
	pnames = c(pnames, "archm")
	pidx[5,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[6,1] = nx+1
	if(pos.matrix[6,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		if(any(!is.na(match(fixed.names, "skewm")))) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	}
	pnames = c(pnames, "skewm")
	pidx[6,2] = nx+pn
	
	
	nx = nx + pn
	pn = 1
	pidx[7,1] = nx+1
	if(pos.matrix[7,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		if(any(!is.na(match(fixed.names, "shapem")))) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	}
	pnames = c(pnames, "shapem")
	pidx[7,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[8,1] = nx+1
	if(pos.matrix[8,3]==1){
		pn = length( seq(pos.matrix[8,1], pos.matrix[8,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			nnx = paste("mxreg", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "mxreg")
	}
	pidx[8,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[9,1] = nx+1
	if(pos.matrix[9,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		if(any(!is.na(match(fixed.names, "omega")))) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	}
	pnames = c(pnames, "omega")
	pidx[9,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[10,1] = nx+1
	if(pos.matrix[10,3]==1){
		pn = length( seq(pos.matrix[10,1], pos.matrix[10,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			nnx = paste("alpha", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "alpha")
	}
	pidx[10,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[11,1] = nx+1
	if(pos.matrix[11,3]==1){
		pn = length( seq(pos.matrix[11,1], pos.matrix[11,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			nnx = paste("beta", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
		#-------------------------------------------
		# special consideration for the iGARCH model
		#-------------------------------------------
		if(vmodel$model == "iGARCH"){
			# last beta not estimated
			pars[nx+pn, 4] = 0
			nnx = paste("beta", pn, sep="")
			# do not allow the last beta to be fixed
			if(any(substr(fixed.names, 1, nchar(nnx))==nnx)) pars[(nx+pn), 2] = 0
		}
	} else{
		pnames = c(pnames, "beta")
	}
	pidx[11,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[12,1] = nx+1
	
	if(pos.matrix[12,3]==1){
		pn = length( seq(pos.matrix[12,1], pos.matrix[12,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			nnx = paste("gamma", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "gamma")
	}
	pidx[12,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[13,1] = nx+1
	
	if(pos.matrix[13,3]==1){
		pn = length( seq(pos.matrix[13,1], pos.matrix[13,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			nnx = paste("eta1", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "eta1")
	}
	pidx[13,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[14,1] = nx+1
	
	if(pos.matrix[14,3]==1){
		pn = length( seq(pos.matrix[14,1], pos.matrix[14,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			nnx = paste("eta2", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "eta2")
	}
	pidx[14,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[15,1] = nx+1
	
	if(pos.matrix[15,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		if(any(!is.na(match(fixed.names, "delta")))) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	} else{
		#-------------------------------------------
		# special consideration for the fGARCH model
		#-------------------------------------------
		#if(vmodel$model == "fGARCH")
		#{
		#	pars[nx+pn, 3] = 1
		#	pars[nx+pn, 1] = fmodel$fpars$delta
		#	pars[nx+pn, 4] = pars[nx+pn, 2] = 0
		#}
	}
	pidx[15,2] = nx+pn
	
	pnames = c(pnames, "delta")
	
	nx = nx + pn
	pn = 1
	pidx[16,1] = nx+1
	
	if(pos.matrix[16,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		if(any(!is.na(match(fixed.names, "lambda")))) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	} else{
		#-------------------------------------------
		# special consideration for the fGARCH model
		#-------------------------------------------
		#if(vmodel$model == "fGARCH")
		#{
		#	pars[nx+pn, 3] = 1
		#	pars[nx+pn, 1] = fmodel$fpars$lambda
		#	pars[nx+pn, 4] = pars[nx+pn, 2] = 0
		#}
	}
	pidx[16,2] = nx+pn
	
	pnames = c(pnames, "lambda")
	
	nx = nx + pn
	pn = 1
	pidx[17,1] = nx+1
	
	if(pos.matrix[17,3]==1){
		pn = length( seq(pos.matrix[17,1], pos.matrix[17,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			nnx = paste("vxreg", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "vxreg")
		
	}
	pidx[17,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[18,1] = nx+1
	
	if(pos.matrix[18,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		if(any(!is.na(match(fixed.names, "skew")))) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	}
	pidx[18,2] = nx+pn
	
	pnames = c(pnames, "skew")
	
	nx = nx + pn
	pn = 1
	pidx[19,1] = nx+1
	
	if(pos.matrix[19,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		if(any(!is.na(match(fixed.names, "shape")))) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	}
	pnames = c(pnames, "shape")
	pidx[19,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[20,1] = nx+1
	
	if(pos.matrix[20,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		if(any(!is.na(match(fixed.names, "ghlambda")))) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	}
	pidx[20,2] = nx+pn
	# Once more for fgarch model pars
	pnames = c(pnames, "ghlambda")
	
	nx = nx + pn
	pn = 1
	pidx[21,1] = nx+1
	
	if(pos.matrix[21,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		if(any(!is.na(match(fixed.names, "skcons")))) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	}
	pidx[21,2] = nx+pn
	pnames = c(pnames, "skcons")
	
	nx = nx + pn
	pn = 1
	pidx[22,1] = nx+1
	
	if(pos.matrix[22,3]==1){
		pn = length( seq(pos.matrix[22,1], pos.matrix[22,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			nnx = paste("skalpha", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "skalpha")
		
	}
	pidx[22,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[23,1] = nx+1
	
	if(pos.matrix[23,3]==1){
		pn = length( seq(pos.matrix[23,1], pos.matrix[23,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			nnx = paste("skgamma", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "skgamma")
		
	}
	pidx[23,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[24,1] = nx+1
	
	if(pos.matrix[24,3]==1){
		pn = length( seq(pos.matrix[24,1], pos.matrix[24,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			nnx = paste("skbeta", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "skbeta")
		
	}
	pidx[24,2] = nx+pn
	
	
	nx = nx + pn
	pn = 1
	pidx[25,1] = nx+1
	
	if(pos.matrix[25,3]==1){
		pn = length( seq(pos.matrix[25,1], pos.matrix[25,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			nnx = paste("skxreg", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "skxreg")
		
	}
	pidx[25,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[26,1] = nx+1
	
	if(pos.matrix[26,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		if(any(!is.na(match(fixed.names, "thskew")))) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	}
	pidx[26,2] = nx+pn
	pnames = c(pnames, "thskew")
	
	nx = nx + pn
	pn = 1
	pidx[27,1] = nx+1
	
	if(pos.matrix[27,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		if(any(!is.na(match(fixed.names, "shcons")))) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	}
	pidx[27,2] = nx+pn
	# Once more for fgarch model pars
	pnames = c(pnames, "shcons")
	
	nx = nx + pn
	pn = 1
	pidx[28,1] = nx+1
	
	if(pos.matrix[28,3]==1){
		pn = length( seq(pos.matrix[28,1], pos.matrix[28,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			nnx = paste("shalpha", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "shalpha")
		
	}
	pidx[28,2] = nx+pn
	
	
	nx = nx + pn
	pn = 1
	pidx[29,1] = nx+1
	
	if(pos.matrix[29,3]==1){
		pn = length( seq(pos.matrix[29,1], pos.matrix[29,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			nnx = paste("shgamma", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "shgamma")
		
	}
	pidx[29,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[30,1] = nx+1
	
	if(pos.matrix[30,3]==1){
		pn = length( seq(pos.matrix[30,1], pos.matrix[30,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			nnx = paste("shbeta", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "shbeta")
		
	}
	pidx[30,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[31,1] = nx+1
	
	if(pos.matrix[31,3]==1){
		pn = length( seq(pos.matrix[31,1], pos.matrix[31,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			nnx = paste("shxreg", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "shxreg")
		
	}
	pidx[31,2] = nx+pn	
	nx = nx + pn
	pn = 1
	pidx[32,1] = nx+1
	
	if(pos.matrix[32,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		if(any(!is.na(match(fixed.names, "thshape")))) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	}
	pidx[32,2] = nx+pn
	pnames = c(pnames, "thshape")
	
	rownames(pars) = pnames
	
	zf = match(fixed.names, rownames(pars))
	if( length(zf)>0 ) pars[zf, 1] = unlist(fixed.pars)
	pars[,"LB"] = NA
	pars[,"UB"] = NA
	vmodel$external.regressors = NULL
	mmodel$external.regressors = NULL
	sbounds[5] = dmodel$exp.rate
	model = list(modelinc = modelinc, modeldata = modeldata, pars = pars,
			sbounds = sbounds, start.pars = start.pars, fixed.pars = fixed.pars, 
			maxOrder = maxOrder, pos.matrix = pos.matrix, pidx = pidx,
			vmodel = vmodel, mmodel = mmodel, dmodel = dmodel)
	ans = new("ACDspec", model = model)	
	return(ans)
}

setMethod(f = "acdspec", definition = .acdspec)
#----------------------------------------------------------------------------------
.getspecacd = function(object)
{
	vmodel = object@model$vmodel
	mmodel = object@model$mmodel
	dmodel = object@model$dmodel
	vmodel$external.regressors = object@model$modeldata$vexdata
	mmodel$external.regressors = object@model$modeldata$mexdata
	spec = acdspec( variance.model = vmodel, mean.model = mmodel,
			distribution.model = dmodel, start.pars = object@model$start.pars, 
			fixed.pars = object@model$fixed.pars)
	return(spec)
}

setMethod(f = "getspec", signature(object = "ACDfit"), definition = .getspecacd)
#----------------------------------------------------------------------------------
# fixed parameters
.setfixedacd = function(object, value){
	# get parameter values
	oldfixed = unlist(object@model$fixed.pars)
	model = object@model
	ipars = model$pars
	pars = unlist(value)
	names(pars) = parnames = tolower(names(pars))
	# included parameters in model
	modelnames = rownames(ipars[which(ipars[,4]==1 | ipars[,2]==1), ])
	inc = NULL
	for(i in seq_along(parnames)){
		if(is.na(match(parnames[i], modelnames))){
			warning( (paste("Unrecognized Parameter in Fixed Values: ", parnames[i], "...Ignored", sep = "")))
		} else{
			inc = c(inc, i)
		}
	}
	fixed.pars = pars[inc]
	names(fixed.pars) = tolower(names(pars[inc]))
	# need to check for duplicates
	xidx = match(names(fixed.pars), names(oldfixed))
	if(any(!is.na(xidx))){
		oldfixed=oldfixed[-na.omit(xidx)]
	}
	if(!is.null(oldfixed) && length(oldfixed)>0) fixed.pars = c(fixed.pars, oldfixed)
	# set parameter values
	vmodel = object@model$vmodel
	mmodel = object@model$mmodel
	dmodel = object@model$dmodel
	vmodel$external.regressors = object@model$modeldata$vexdata
	mmodel$external.regressors = object@model$modeldata$mexdata
	dmodel$shape.regressors = object@model$modeldata$shxdata
	dmodel$skew.regressors = object@model$modeldata$skxdata
	tmp = acdspec( variance.model = vmodel, mean.model = mmodel,
			distribution.model = dmodel, start.pars  = model$start.pars, 
			fixed.pars = as.list(fixed.pars))
	tmp@model$pars[tmp@model$pars[,2]==0,5:6] = object@model$pars[tmp@model$pars[,2]==0,5:6]
	setbounds(tmp)<-list(shape=object@model$sbounds[3:4], skew = object@model$sbounds[1:2])
	
	return(tmp)
}
setMethod(f="setfixed<-", signature= c(object = "ACDspec", value = "vector"), definition = .setfixedacd)
#----------------------------------------------------------------------------------
# starting parameters
.setstartacd = function(object, value){
	# get parameter values
	model = object@model
	ipars = model$pars
	pars = unlist(value)
	names(pars) = parnames = tolower(names(pars))
	# included parameters in model
	modelnames = rownames(ipars[which(ipars[,4]==1 | ipars[,2]==1), ])
	inc = NULL
	for(i in seq_along(parnames)){
		if(is.na(match(parnames[i], modelnames))){
			warning( (paste("Unrecognized Parameter in Start Values: ", parnames[i], "...Ignored", sep = "")))
		} else{
			inc = c(inc, i)
		}
	}
	start.pars = pars[inc]
	names(start.pars) = tolower(names(pars[inc]))
	# set parameter values
	# set parameter values
	vmodel = object@model$vmodel
	mmodel = object@model$mmodel
	dmodel = object@model$dmodel
	vmodel$external.regressors = object@model$modeldata$vexdata
	mmodel$external.regressors = object@model$modeldata$mexdata
	dmodel$shape.regressors = object@model$modeldata$shxdata
	dmodel$skew.regressors = object@model$modeldata$skxdata
	
	tmp = acdspec( variance.model = vmodel, mean.model = mmodel,
			distribution.model = dmodel, fixed.pars  = model$fixed.pars, 
			start.pars = as.list(start.pars))
	tmp@model$pars[tmp@model$pars[,2]==0,5:6] = object@model$pars[tmp@model$pars[,2]==0,5:6]
	setbounds(tmp)<-list(shape=object@model$sbounds[3:4], skew = object@model$sbounds[1:2])
	return(tmp)
}

setMethod(f="setstart<-", signature= c(object = "ACDspec", value = "vector"), definition = .setstartacd)
#----------------------------------------------------------------------------------
.checkallfixed = function( spec ){
	# check that a given spec with fixed parameters
	model = spec@model
	pars = model$pars
	pnames = rownames(pars)
	estpars = pnames[as.logical(pars[,2] * pars[,3] + pars[,3] * pars[,4])]
	return( estpars )
}
#----------------------------------------------------------------------------------
# set parameters bounds
# Set the lower and upper bounds
# value is a list with names parameters taking 2 values (lower and upper)
# e.g. value  = list(alpha1 = c(0, 0.1), beta1 = c(0.9, 0.99))
# special attention given to skew/shape in case of dynamics
.acdsetbounds = function(object, value){
	model = object@model
	ipars = model$pars
	parnames = tolower(names(value))
	if(model$modelinc[21]>0 && any(parnames=="skew")){
		idx = which(parnames=="skew")
		if(!is.na(value[[idx]][1])) object@model$sbounds[1] = value[[idx]][1]
		if(!is.na(value[[idx]][2])) object@model$sbounds[2] = value[[idx]][2]
	}
	if(model$modelinc[27]>0 && any(parnames=="shape")){
		idx = which(parnames=="shape")
		if(!is.na(value[[idx]][1])) object@model$sbounds[3] = value[[idx]][1]
		if(!is.na(value[[idx]][2])) object@model$sbounds[4] = value[[idx]][2]
	}
	# included parameters in model
	modelnames = rownames(ipars[which(ipars[,4] == 1), ])
	sp = na.omit(match(parnames, modelnames))
	if(length(sp)>0){
		for(i in 1:length(sp)){
			#if(length(value[[modelnames[sp[i]]]])!=2)
			ipars[modelnames[sp[i]], 5] = as.numeric(value[[modelnames[sp[i]]]][1])
			ipars[modelnames[sp[i]], 6] = as.numeric(value[[modelnames[sp[i]]]][2])
		}
	}
	object@model$pars = ipars
	return(object)
}
setReplaceMethod(f="setbounds", signature= c(object = "ACDspec", value = "vector"), definition = .acdsetbounds)
#----------------------------------------------------------------------------------
# estimation
acdfit = function(spec, data, solver = "ucminf", out.sample = 0, solver.control = list(), 
		fit.control = list(stationarity = 0, fixed.se = 0, scale = 0, n.sim = 2000), 
		skew0 = NULL, shape0 = NULL, cluster = NULL, ...)
{
	UseMethod("acdfit")
}
.acdfitswitch = function(spec, data, solver = "ucminf", out.sample = 0, solver.control = list(), 
		fit.control = list(stationarity = 0, fixed.se = 0, scale = 0, n.sim = 2000), 
		skew0 = NULL, shape0 = NULL, cluster = NULL, ...)
{
	switch(spec@model$vmodel$model, 
			sGARCH = .acdfit(spec = spec, data = data, out.sample = out.sample,
					solver = solver, solver.control = solver.control, fit.control = fit.control, 
					skew0 = skew0, shape0 = shape0, cluster = cluster, ...),
			csGARCH = .acdfit(spec = spec, data = data, out.sample = out.sample,
					solver = solver, solver.control = solver.control, fit.control = fit.control, 
					skew0 = skew0, shape0 = shape0, cluster = cluster, ...),
			mcsGARCH = .mcsacdfit(spec = spec, data = data, out.sample = out.sample,
					solver = solver, solver.control = solver.control, fit.control = fit.control, 
					skew0 = skew0, shape0 = shape0, cluster = cluster, ...))
}
setMethod("acdfit", signature(spec = "ACDspec"), .acdfitswitch)
#----------------------------------------------------------------------------------
# filter
acdfilter = function(spec, data, out.sample = 0,  n.old = NULL, skew0 = NULL, shape0 = NULL, ...)
{
	UseMethod("acdfilter")
}

.acdfilterswitch = function(spec, data, out.sample = 0,  n.old = NULL, skew0 = NULL, shape0 = NULL, ...)
{
	switch(spec@model$vmodel$model, 
			sGARCH = .acdfilter(spec = spec, data = data, out.sample = out.sample,
					n.old = n.old, skew0 = skew0, shape0 = shape0, ...),
			csGARCH = .acdfilter(spec = spec, data = data, out.sample = out.sample,
					n.old = n.old, skew0 = skew0, shape0 = shape0, ...),
			mcsGARCH = .mcsacdfilter(spec = spec, data = data, out.sample = out.sample,
					n.old = n.old, skew0 = skew0, shape0 = shape0, ...))
}
setMethod("acdfilter", signature(spec = "ACDspec"), .acdfilterswitch)
#----------------------------------------------------------------------------------
# forecast
acdforecast = function(fitORspec, data = NULL, n.ahead = 10, n.roll = 0, out.sample = 0,
		external.forecasts = list(mregfor = NULL, vregfor = NULL, 
				skregfor = NULL, shregfor = NULL), m.sim = 1000, cluster = NULL, 
		skew0 = NULL, shape0 = NULL, ...)
{
	UseMethod("acdforecast")
}

.acdforecastswitch1 = function(fitORspec, n.ahead = 10, n.roll = 0, 
		external.forecasts = list(mregfor = NULL, vregfor = NULL, 
				skregfor = NULL, shregfor = NULL), m.sim = 1000, cluster = NULL, ...)
{
	fit = fitORspec
	switch(fit@model$vmodel$model, 
			sGARCH = .acdforecast(fit = fit, n.ahead = n.ahead, 
					n.roll = n.roll, 
					external.forecasts = external.forecasts, m.sim = m.sim, 
					cluster = cluster, ...),
			csGARCH = .acdforecast(fit = fit, n.ahead = n.ahead, 
					n.roll = n.roll,  
					external.forecasts = external.forecasts, m.sim = m.sim, 
					cluster = cluster, ...),
			mcsGARCH = .mcsacdforecast(fit = fit, n.ahead = n.ahead, 
					n.roll = n.roll, 
					external.forecasts = external.forecasts, 
					m.sim = m.sim, cluster = cluster, ...))
}

.acdforecastswitch2 = function(fitORspec, data = NULL, n.ahead = 10, n.roll = 0, out.sample = 0,
		external.forecasts = list(mregfor = NULL, vregfor = NULL, 
				skregfor = NULL, shregfor = NULL), m.sim = 1000, 
		cluster = NULL, skew0 = NULL, shape0 = NULL, ...)
{
	spec = fitORspec
	if(spec@model$vmodel$model=="mcsGARCH") stop("\nacdforecast-->error: mcsGARCH model does not support specification dispatch method for forecast.")
	switch(spec@model$vmodel$model, 
			sGARCH = .acdforecast2(spec = spec, data = data, n.ahead = n.ahead, 
					n.roll = n.roll, out.sample = out.sample,
					external.forecasts = external.forecasts, m.sim = m.sim, 
					cluster = cluster, skew0 = skew0, shape0 = shape0, ...),
			csGARCH = .acdforecast2(spec = spec, data = data, n.ahead = n.ahead, 
					n.roll = n.roll,  out.sample = out.sample,
					external.forecasts = external.forecasts, m.sim = m.sim, 
					cluster = cluster, skew0 = skew0, shape0 = shape0, ...))
}
setMethod("acdforecast", signature(fitORspec = "ACDfit"), .acdforecastswitch1)
setMethod("acdforecast", signature(fitORspec = "ACDspec"), .acdforecastswitch2)

#----------------------------------------------------------------------------------
# simulation
acdsim = function(fit, n.sim = 1000, n.start = 0, m.sim = 1, 
		presigma = NA, prereturns = NA, preresiduals = NA, preskew = NA, 
		preshape = NA, rseed = NA, mexsimdata = NULL, vexsimdata = NULL, 
		skxsimdata = NULL, shxsimdata = NULL, ...)
{
	UseMethod("acdsim")
}

.acdsim = function(fit, n.sim = 1000, n.start = 0, m.sim = 1, 
		presigma = NA, prereturns = NA, preresiduals = NA, preskew = NA, 
		preshape = NA, rseed = NA,  mexsimdata = NULL, vexsimdata = NULL, 
		skxsimdata = NULL, shxsimdata = NULL, ...){
	ans = switch(fit@model$vmodel$model,
			sGARCH = .sacdsim(fit = fit, n.sim = n.sim, n.start = n.start, 
					m.sim = m.sim, presigma = presigma, prereturns = prereturns, 
					preresiduals = preresiduals, preskew = preskew, 
					preshape = preshape, rseed = rseed, mexsimdata = mexsimdata, 
					vexsimdata = vexsimdata, skxsimdata = skxsimdata, 
					shxsimdata = shxsimdata, ...),
			csGARCH = .csacdsim(fit = fit, n.sim = n.sim, n.start = n.start, 
					m.sim = m.sim, presigma = presigma, prereturns = prereturns, 
					preresiduals = preresiduals, preskew = preskew, 
					preshape = preshape, rseed = rseed, mexsimdata = mexsimdata, 
					vexsimdata = vexsimdata, skxsimdata = skxsimdata, 
					shxsimdata = shxsimdata, ...),
			mcsGARCH = .mcsacdsim(fit = fit, n.sim = n.sim, n.start = n.start, 
					m.sim = m.sim, presigma = presigma, prereturns = prereturns, 
					preresiduals = preresiduals, preskew = preskew, 
					preshape = preshape, rseed = rseed, mexsimdata = mexsimdata, 
					vexsimdata = vexsimdata, skxsimdata = skxsimdata, 
					shxsimdata = shxsimdata, ...))
	return(ans)
}
setMethod("acdsim", signature(fit = "ACDfit"), .acdsim)
#----------------------------------------------------------------------------------
# path
acdpath = function(spec, n.sim = 1000, n.start = 0, m.sim = 1, presigma = NA, 
		prereturns = NA, preresiduals = NA, preskew = NA, preshape = NA, rseed = NA,  
		mexsimdata = NULL, vexsimdata = NULL, skxsimdata = NULL, shxsimdata = NULL, 
		cluster = NULL, ...)
{
	UseMethod("acdpath")
}

.acdpath = function(spec, n.sim = 1000, n.start = 0, m.sim = 1, presigma = NA, 
		prereturns = NA, preresiduals = NA, preskew = NA, preshape = NA, rseed = NA,  
		mexsimdata = NULL, vexsimdata = NULL, skxsimdata = NULL, shxsimdata = NULL, 
		cluster = NULL, ...){
	ans = switch(spec@model$vmodel$model,
			sGARCH = .sacdpath(spec = spec, n.sim = n.sim, n.start = n.start, 
					m.sim = m.sim, presigma = presigma, prereturns = prereturns, 
					preresiduals = preresiduals, preskew = preskew, 
					preshape = preshape, rseed = rseed, mexsimdata = mexsimdata, 
					vexsimdata = vexsimdata, skxsimdata = skxsimdata, 
					shxsimdata = shxsimdata, cluster = cluster, ...),
			csGARCH = .csacdpath(spec = spec, n.sim = n.sim, n.start = n.start, 
					m.sim = m.sim, presigma = presigma, prereturns = prereturns, 
					preresiduals = preresiduals, preskew = preskew, 
					preshape = preshape, rseed = rseed, mexsimdata = mexsimdata, 
					vexsimdata = vexsimdata, skxsimdata = skxsimdata, 
					shxsimdata = shxsimdata, cluster = cluster, ...),
			mcsGARCH = .sacdpath(spec = spec, n.sim = n.sim, n.start = n.start, 
					m.sim = m.sim, presigma = presigma, prereturns = prereturns, 
					preresiduals = preresiduals, preskew = preskew, 
					preshape = preshape, rseed = rseed, mexsimdata = mexsimdata, 
					vexsimdata = vexsimdata, skxsimdata = skxsimdata, 
					shxsimdata = shxsimdata, cluster = cluster, ...))
	return(ans)
}
setMethod("acdpath", signature(spec = "ACDspec"), .acdpath)
#----------------------------------------------------------------------------------
# rolling estimation/forecast
acdroll = function(spec, data, n.ahead = 1, forecast.length = 500, 
		n.start = NULL, refit.every = 25, refit.window = c("recursive", "moving"), 
		window.size = NULL, solver = "ucminf", fit.control = list(), solver.control = list(),
		calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05), cluster = NULL,
		keep.coef = TRUE, fixARMA = TRUE, fixGARCH = TRUE, 
		fixUBShape = TRUE, UBShapeAdd = 0, fixGHlambda = TRUE, 
		compareGARCH=c("LL", "none"),...)
{
	UseMethod("acdroll")
}

.acdrollswitch = function(spec, data, n.ahead = 1, forecast.length = 500, n.start = NULL, 
		refit.every = 25, refit.window = c("recursive", "moving"), 
		window.size = NULL, solver = "ucminf", fit.control = list(), 
		solver.control = list(), calculate.VaR = TRUE, VaR.alpha = c(0.01, 
				0.05), cluster = NULL, keep.coef = TRUE, fixARMA = TRUE, fixGARCH = TRUE, 
		fixUBShape = TRUE, UBShapeAdd = 0, fixGHlambda = TRUE, compareGARCH = c("LL", "none"),
		...)
{
	ans = switch(spec@model$vmodel$model,
			sGARCH = .acdroll(spec = spec, data = data, n.ahead = n.ahead, 
					forecast.length = forecast.length, n.start = n.start, 
					refit.every = refit.every, refit.window = refit.window[1], 
					window.size = window.size, solver = solver, fit.control = fit.control, 
					solver.control = solver.control, calculate.VaR = calculate.VaR, 
					VaR.alpha = VaR.alpha, cluster = cluster, keep.coef = keep.coef, 
					fixARMA = fixARMA, fixGARCH = fixGARCH, fixUBShape = fixUBShape, 
					UBShapeAdd = UBShapeAdd, fixGHlambda = fixGHlambda, 
					compareGARCH = compareGARCH),
			csGARCH = .acdroll(spec = spec, data = data, n.ahead = n.ahead, 
					forecast.length = forecast.length, n.start = n.start, 
					refit.every = refit.every, refit.window = refit.window[1], 
					window.size = window.size, solver = solver, fit.control = fit.control, 
					solver.control = solver.control, calculate.VaR = calculate.VaR, 
					VaR.alpha = VaR.alpha, cluster = cluster, keep.coef = keep.coef, 
					fixARMA = fixARMA, fixGARCH = fixGARCH, fixUBShape = fixUBShape, 
					UBShapeAdd = UBShapeAdd, fixGHlambda = fixGHlambda, 
					compareGARCH = compareGARCH),
			mcsGARCH = .acdroll.mcs(spec = spec, data = data, n.ahead = n.ahead, 
					forecast.length = forecast.length, n.start = n.start, 
					refit.every = refit.every, refit.window = refit.window[1], 
					window.size = window.size, solver = solver, fit.control = fit.control, 
					solver.control = solver.control, calculate.VaR = calculate.VaR, 
					VaR.alpha = VaR.alpha, cluster = cluster, keep.coef = keep.coef, 
					fixARMA = fixARMA, fixGARCH = fixGARCH, fixUBShape = fixUBShape, 
					UBShapeAdd = UBShapeAdd, fixGHlambda = fixGHlambda, 
					compareGARCH = compareGARCH, ...)
	)
	return(ans)
}
setMethod("acdroll", signature(spec = "ACDspec"),  .acdrollswitch)

.acdresumerollswitch = function(object, spec = NULL, solver = "ucminf", fit.control = list(), 
		solver.control = list(), cluster = NULL, fixARMA = NULL, fixGARCH = NULL, 
		fixUBShape = NULL, UBShapeAdd = NULL, fixGHlambda = NULL, compareGARCH = NULL)
{
	ans = switch(object@model$spec@model$vmodel$model,
			sGARCH = .acdresumeroll(object = object, spec = spec, solver = solver, 
					fit.control = fit.control, solver.control = solver.control, 
					cluster = cluster, fixARMA = fixARMA, fixGARCH = fixGARCH, 
					fixUBShape = fixUBShape, UBShapeAdd = UBShapeAdd, 
					fixGHlambda = fixGHlambda, compareGARCH = compareGARCH),
			csGARCH = .acdresumeroll(object = object, spec = spec, solver = solver, 
					fit.control = fit.control, solver.control = solver.control, 
					cluster = cluster, fixARMA = fixARMA, fixGARCH = fixGARCH, 
					fixUBShape = fixUBShape, UBShapeAdd = UBShapeAdd, 
					fixGHlambda = fixGHlambda, compareGARCH = compareGARCH),
			mcsGARCH = .acdresumeroll.mcs(object = object, spec = spec, solver = solver, 
					fit.control = fit.control, solver.control = solver.control, 
					cluster = cluster, fixARMA = fixARMA, fixGARCH = fixGARCH, 
					fixUBShape = fixUBShape, UBShapeAdd = UBShapeAdd, 
					fixGHlambda = fixGHlambda, compareGARCH = compareGARCH)
	)
	return(ans)
}
setMethod("resume", signature(object = "ACDroll"),  .acdresumerollswitch)
#-------------------------------------------------------------------------------
# convergence method
.acdconvergence = function(object)
{
	return( object@fit$convergence )
}
setMethod("convergence", signature(object = "ACDfit"),  .acdconvergence)

.acdconvergenceroll = function(object){
	nonc = object@model$noncidx
	if(is.null(nonc)){
		ans = 0 
	} else{
		ans = 1
		attr(ans, 'nonconverged')<-nonc
	}
	return(ans)
}

setMethod("convergence", signature(object = "ACDroll"),  definition = .acdconvergenceroll)

.acdmulticonvergence = function(object)
{
	return( sapply(object@fit, function(x) convergence(x) ) )
}
setMethod("convergence", signature(object = "ACDmultifit"),  .acdmulticonvergence)


#----------------------------------------------------------------------------------
# quantile method
.acdfitquantile = function(x, probs = c(0.01, 0.05))
{
	if(class(x)=="ACDroll"){
		d = x@model$spec@model$dmodel$model
		skew = x@forecast$density[,"Skew"]
		shape = x@forecast$density[,"Shape"]
		lambda = x@forecast$density[,"Shape(GIG)"]
		s = x@forecast$density[,"Sigma"]
		m = x@forecast$density[,"Mu"]
		Q = matrix(NA, ncol = length(probs), nrow = length(s))
		for(i in 1:length(probs)) Q[,i] = as.numeric(m) + qdist(d, probs[i], skew = skew, shape = shape, 
					lambda = lambda) * as.numeric(s)
		colnames(Q) = paste("q[", probs,"]", sep="")
		Q = xts(Q, as.POSIXct(rownames(x@forecast$density)))
	} else{
		d = x@model$dmodel$model
		mu = fitted(x)
		sig = sigma(x)
		tskew = skew(x)
		tshape = shape(x)
		if(d=="ghyp") lambda = coef(x)["ghlambda"] else lambda = 0
		Q = matrix(NA, nrow = length(sig), ncol = length(probs))
		for(i in seq_along(probs)) Q[,i] = mu + sig*qdist(d, probs[i], skew = tskew, shape = tshape)
		colnames(Q) = paste("q[", probs,"]", sep="")
		Q = xts(Q, index(sig))
	}
	return(Q)
}

setMethod("quantile", signature(x = "ACDfit"),  .acdfitquantile)
setMethod("quantile", signature(x = "ACDfilter"),  .acdfitquantile)
setMethod("quantile", signature(x = "ACDroll"),  .acdfitquantile)

#-------------------------------------------------------------------------------
# pit method
.acdfitpit = function(object)
{
	if(class(object)=="ACDroll"){
		d = object@model$spec@model$dmodel$model
		skew = skew(object)
		shape = shape(object)
		lambda = object@forecast$density[,"Shape(GIG)"]
		s = sigma(object)
		m = fitted(object)
		r = object@forecast$density[,"Realized"]
		ans =  pdist(d, q = r, mu = as.numeric(m), sigma = as.numeric(s), 
				skew = as.numeric(skew), shape = as.numeric(shape), lambda = as.numeric(ghlambda))
		ans = xts(ans, as.POSIXct(rownames(object@forecast$density)))
		colnames(ans) = "pit"
	} else{
		d = object@model$dmodel$model
		di = .DistributionBounds(d)
		if(di$include.skew)  skew  = skew(object) else skew = 0
		if(di$include.shape) shape = shape(object) else shape = 0
		if(di$include.ghlambda) ghlambda = object@model$pars["ghlambda",1] else ghlambda = 0
		s = sigma(object)
		m = fitted(object)
		r = object@model$modeldata$data[1:object@model$modeldata$T]
		ans =  pdist(d, q = r, mu = as.numeric(m), sigma = as.numeric(s), 
				skew = as.numeric(skew), shape = as.numeric(shape), 
				lambda = as.numeric(ghlambda))
		ans = xts(ans, index(s))
		colnames(ans) = "pit"
	}
	return(ans)
}
setMethod("pit", signature(object = "ACDfit"), .acdfitpit)
setMethod("pit", signature(object = "ACDfilter"), .acdfitpit)
setMethod("pit", signature(object = "ACDroll"), .acdfitpit)
#-------------------------------------------------------------------------------
# coef method
.acdfitcoef = function(object)
{
	if(is(object, "ACDfit")){
		return(object@fit$coef)
	} else{
		return(object@model$pars[object@model$pars[,3]==1, 1])
	}
}

setMethod("coef", signature(object = "ACDfit"), .acdfitcoef)
setMethod("coef", signature(object = "ACDfilter"), .acdfitcoef)

.acdrollcoef = function(object)
{
	if(!is.null(object@model$noncidx)) stop("\nObject contains non-converged estimation windows.")
	return(object@model$coef) 
}

.acdmulticoef = function(object)
{
	if(is(object, "ACDmultifit")){
		ans = sapply(object@fit, function(x) coef(x))
	} else{
		ans = sapply(object@filter, function(x) coef(x))
	}
	return(ans)
}

setMethod("coef", signature(object = "ACDroll"), .acdrollcoef)
setMethod("coef", signature(object = "ACDmultifit"), .acdmulticoef)
setMethod("coef", signature(object = "ACDmultifilter"), .acdmulticoef)

#-------------------------------------------------------------------------------
# conditional mean (fitted method)
.acdfitted = function(object)
{
	if(class(object)[1] == "ACDfit" | class(object)[1] == "ACDfilter"){
		D = object@model$modeldata$index[1:object@model$modeldata$T]
	}
	ans = switch(class(object)[1],
			ACDfit = xts(object@fit$fitted.values, D),
			ACDfilter = xts(object@filter$fitted.values, D),
			ACDforecast = object@forecast$seriesFor,
			ACDsim = {
				ans = object@simulation$seriesSim
				rownames(ans) = paste("T+",1:NROW(object@simulation$seriesSim), sep="")
				return(ans)
			},
			ACDpath ={
				ans = object@path$seriesSim
				rownames(ans) = paste("T+",1:NROW(object@path$seriesSim), sep="")
				return(ans)
			},
			ACDroll = as.xts(object@forecast$density[,"Mu",drop=FALSE]))
	return(ans)
}

.acdmultifitted = function(object)
{
	if(object@desc$type == "equal"){
		if(is(object, "ACDmultifit")){
			ans = sapply(object@fit, FUN = function(x){
						if(convergence(x) == 0) fitted(x) else rep(NA, length = length(x@fit$data))
					}, simplify = TRUE)
		} else if(is(object, "ACDmultifilter")){
			ans = sapply(object@fit, function(x) fitted(x))
		} else{
			# assuming n.ahead=1 else will
			if(object@desc$n.ahead==1){
				ans = sapply(object@forecast, function(x) fitted(x))
			} else{
				ans = lapply(object@forecast, function(x) fitted(x))
			}
		}
	} else{
		if(is(object, "ACDmultifit")){
			ans = lapply(object@fit, FUN = function(x){
						if(convergence(x) == 0) fitted(x) else rep(NA, length = length(x@fit$data))
					}, simplify = TRUE)
		} else if(is(object, "ACDmultifilter")){
			ans = lapply(object@fit, function(x) fitted(x))
		} else{
			ans = lapply(object@forecast, function(x) fitted(x))
		}
	}
	return(ans)
}
setMethod("fitted", signature(object = "ACDfit"), .acdfitted)
setMethod("fitted", signature(object = "ACDfilter"), .acdfitted)
setMethod("fitted", signature(object = "ACDforecast"), .acdfitted)
setMethod("fitted", signature(object = "ACDpath"), .acdfitted)
setMethod("fitted", signature(object = "ACDsim"), .acdfitted)
setMethod("fitted", signature(object = "ACDroll"), .acdfitted)
setMethod("fitted", signature(object = "ACDmultifit"), .acdmultifitted)
setMethod("fitted", signature(object = "ACDmultifilter"), .acdmultifitted)
setMethod("fitted", signature(object = "ACDmultiforecast"), .acdmultifitted)

#-------------------------------------------------------------------------------
# residuals method
.acdfitresids = function(object, standardize = FALSE)
{
	if(is(object, "ACDfit")){
		if(standardize) ans = object@fit$residuals/object@fit$sigma else  ans = object@fit$residuals
	} else{
		if(standardize) ans = object@filter$residuals/object@filter$sigma else  ans = object@filter$residuals
	}
	return(ans)
}

.acdmultifitresids = function(object, standardize = FALSE)
{
	if(object@desc$type == "equal"){
		if(is(object, "ACDmultifit")){
			ans = sapply(object@fit, FUN = function(x){
						if(convergence(x) == 0) residuals(x, standardize) else rep(NA, length = length(x@fit$data))
					}, simplify = TRUE)
		} else {
			ans = sapply(object@fit, function(x) residuals(x, standardize))
		}
	} else{
		if(is(object, "ACDmultifit")){
			ans = lapply(object@fit, FUN = function(x){
						if(convergence(x) == 0) residuals(x, standardize) else rep(NA, length = length(x@fit$data))
					}, simplify = TRUE)
		} else{
			ans = lapply(object@fit, function(x) residuals(x, standardize))
		}
	}
	return(ans)
}

setMethod("residuals", signature(object = "ACDfit"), .acdfitresids)
setMethod("residuals", signature(object = "ACDfilter"), .acdfitresids)
setMethod("residuals", signature(object = "ACDmultifilter"), .acdmultifitresids)
setMethod("residuals", signature(object = "ACDmultifit"), .acdmultifitresids)


#-------------------------------------------------------------------------------
# conditional sigma
.acdsigma = function(object)
{
	if(class(object)[1] == "ACDfit" | class(object)[1] == "ACDfilter"){
		D = object@model$modeldata$index[1:object@model$modeldata$T]
	}
	ans = switch(class(object)[1],
			ACDfit = xts(object@fit$sigma, D),
			ACDfilter = xts(object@filter$sigma, D),
			ACDforecast = object@forecast$sigmaFor,
			ACDsim = {
				ans = object@simulation$sigmaSim
				rownames(ans) = paste("T+",1:NROW(object@simulation$sigmaSim), sep="")
				return(ans)
			},
			ACDpath ={
				ans = object@path$sigmaSim
				rownames(ans) = paste("T+",1:NROW(object@path$sigmaSim), sep="")
				return(ans)
			},
			ACDroll = as.xts(object@forecast$density[,"Sigma",drop=FALSE]))
	return(ans)
}
# The multifit should be array with dates?
.acdmultifitsigma = function(object)
{
	if(object@desc$type == "equal"){
		if(is(object, "ACDmultifit")){
			ans = sapply(object@fit, FUN = function(x){
						if(convergence(x) == 0) sigma(x) else rep(NA, length = length(x@fit$data))
					}, simplify = TRUE)
		} else if(is(object, "ACDmultifilter")){
			ans = sapply(object@filter, function(x) sigma(x))
		} else{
			# assuming n.ahead=1 else will
			if(object@desc$n.ahead==1){
				ans = sapply(object@forecast, function(x) sigma(x))
			} else{
				ans = lapply(object@forecast, function(x) sigma(x))
			}
		}
	} else{
		if(is(object, "ACDmultifit")){
			ans = lapply(object@fit, FUN = function(x){
						if(convergence(x) == 0) sigma(x) else rep(NA, length = length(x@fit$data))
					}, simplify = TRUE)
		} else if(is(object, "ACDmultifilter")){
			ans = lapply(object@filter, function(x) sigma(x))
		} else{
			ans = lapply(object@forecast, function(x) sigma(x))
		}
	}
	return(ans)
}



setMethod("sigma", signature(object = "ACDfit"), .acdsigma)
setMethod("sigma", signature(object = "ACDfilter"), .acdsigma)
setMethod("sigma", signature(object = "ACDforecast"), .acdsigma)
setMethod("sigma", signature(object = "ACDpath"), .acdsigma)
setMethod("sigma", signature(object = "ACDsim"), .acdsigma)
setMethod("sigma", signature(object = "ACDroll"), .acdsigma)
setMethod("sigma", signature(object = "ACDmultifit"), .acdmultifitsigma)
setMethod("sigma", signature(object = "ACDmultifilter"), .acdmultifitsigma)
setMethod("sigma", signature(object = "ACDmultiforecast"), .acdmultifitsigma)

#-------------------------------------------------------------------------------
# conditional skew
skew = function(object, transformed = TRUE, ...)
{
	UseMethod("skew")
}

.acdskew = function(object, transformed = TRUE)
{
	if(class(object)[1] == "ACDfit" | class(object)[1] == "ACDfilter"){
		D = object@model$modeldata$index[1:object@model$modeldata$T]
	}
	if(is(object, "ACDfit")){
		if(transformed) ans = xts(object@fit$tskew, D) else ans = xts(object@fit$tempskew, D)
	} else if(is(object, "ACDfilter")){
		if(transformed) ans = xts(object@filter$tskew, D) else ans = xts(object@filter$tempskew, D)
	} else if(is(object, "ACDforecast")){
		if(transformed){
			ans = object@forecast$tskewFor
		} else{
			ans = apply(object@forecast$tskewFor, 2, function(x) logtransform(x, object@model$sbounds[1], object@model$sbounds[2], inverse = TRUE))
		}
	} else if(is(object, "ACDpath")){
		if(transformed){
			ans = object@path$skewSim
			rownames(ans) = paste("T+",1:NROW(object@path$skewSim), sep="")
		} else{
			ans = apply(object@path$skewSim, 2, function(x) logtransform(x, object@model$sbounds[1], object@model$sbounds[2], inverse = TRUE))
			rownames(ans) = paste("T+",1:NROW(object@path$skewSim), sep="")
			
		}
	} else if(is(object, "ACDsim")){
		if(transformed){
			ans = object@simulation$skewSim
			rownames(ans) = paste("T+",1:NROW(object@simulation$sigmaSim), sep="")
		} else{
			ans = apply(object@simulation$skewSim, 2, function(x) logtransform(x, object@model$sbounds[1], object@model$sbounds[2], inverse = TRUE))
			rownames(ans) = paste("T+",1:NROW(object@simulation$sigmaSim), sep="")
		}
	} else if(is(object, "ACDroll")){
		if(transformed){
			ans = as.xts(object@forecast$density[,"Skew",drop=FALSE])
		} else{
			stop("\nACDroll object does return the (possibly) dynamic bounds needed for the transformation.")
		}
	} else{
		ans = NA
	}
	return(ans)
}

.acdmultiskew = function(object, transformed = TRUE)
{
	if(object@desc$type == "equal"){
		if(is(object, "ACDmultifit")){
			ans = sapply(object@fit, FUN = function(x){
						if(convergence(x) == 0) skew(x, transformed) else rep(NA, length = length(x@fit$data))
					}, simplify = TRUE)
		} else if(is(object, "ACDmultifilter")){
			ans = sapply(object@filter, function(x) skew(x, transformed))
		} else{
			# assuming n.ahead=1 else will
			if(object@desc$n.ahead==1){
				ans = sapply(object@forecast, function(x) skew(x, transformed))
			} else{
				ans = lapply(object@forecast, function(x) skew(x, transformed))
			}
		}
	} else{
		if(is(object, "ACDmultifit")){
			ans = lapply(object@fit, FUN = function(x){
						if(convergence(x) == 0) skew(x, transformed) else rep(NA, length = length(x@fit$data))
					}, simplify = TRUE)
		} else if(is(object, "ACDmultifilter")){
			ans = lapply(object@filter, function(x) skew(x, transformed))
		} else{
			ans = lapply(object@forecast, function(x) skew(x, transformed))
		}
	}
	return(ans)
}


setMethod("skew", signature(object = "ACDfit"), .acdskew)
setMethod("skew", signature(object = "ACDfilter"), .acdskew)
setMethod("skew", signature(object = "ACDforecast"), .acdskew)
setMethod("skew", signature(object = "ACDpath"), .acdskew)
setMethod("skew", signature(object = "ACDsim"), .acdskew)
setMethod("skew", signature(object = "ACDroll"), .acdskew)

setMethod("skew", signature(object = "ACDmultifit"), .acdmultiskew)
setMethod("skew", signature(object = "ACDmultifilter"), .acdmultiskew)
setMethod("skew", signature(object = "ACDmultiforecast"), .acdmultiskew)


#-------------------------------------------------------------------------------
# conditional shape
shape = function(object, transformed = TRUE, ...)
{
	UseMethod("shape")
}

.acdshape = function(object, transformed = TRUE)
{
	if(class(object)[1] == "ACDfit" | class(object)[1] == "ACDfilter"){
		D = object@model$modeldata$index[1:object@model$modeldata$T]
	}
	if(is(object, "ACDfit")){
		if(transformed) ans = xts(object@fit$tshape, D) else ans = xts(object@fit$tempshape, D)
	} else if(is(object, "ACDfilter")){
		if(transformed) ans = xts(object@filter$tshape, D) else ans = xts(object@filter$tempshape, D)
	} else if(is(object, "ACDforecast")){
		if(transformed){
			ans = object@forecast$tshapeFor
		} else{
			ans = apply(object@forecast$tshapeFor, 2, function(x) exptransform(x, object@model$sbounds[3], object@model$sbounds[4], object@model$sbounds[5], inverse = TRUE))
		}
	} else if(is(object, "ACDpath")){
		if(transformed){
			ans = object@path$shapeSim
			rownames(ans) = paste("T+",1:NROW(object@path$shapeSim), sep="")
		} else{
			ans = apply(object@path$shapeSim, 2, function(x) exptransform(x, object@model$sbounds[3], object@model$sbounds[4], object@model$sbounds[5], inverse = TRUE))
			rownames(ans) = paste("T+",1:NROW(object@path$shapeSim), sep="")
		}
	}
	else if(is(object, "ACDsim")){
		if(transformed){
			ans = object@simulation$shapeSim
			rownames(ans) = paste("T+",1:NROW(object@simulation$shapeSim), sep="")
		} else{
			ans = apply(object@simulation$shapeSim, 2, function(x) exptransform(x, object@model$sbounds[3], object@model$sbounds[4], object@model$sbounds[5], inverse = TRUE))
			rownames(ans) = paste("T+",1:NROW(object@simulation$shapeSim), sep="")
		}
	} else if(is(object, "ACDroll")){
		if(transformed){
			ans = as.xts(object@forecast$density[,"Shape",drop=FALSE])
		} else{
			stop("\nACDroll object does return the (possibly) dynamic bounds needed for the transformation.")
		}
	} else{
		ans = NA
	}
	return(ans)
}

.acdmultishape = function(object, transformed = TRUE)
{
	if(object@desc$type == "equal"){
		if(is(object, "ACDmultifit")){
			ans = sapply(object@fit, FUN = function(x){
						if(convergence(x) == 0) shape(x, transformed) else rep(NA, length = length(x@fit$data))
					}, simplify = TRUE)
		} else if(is(object, "ACDmultifilter")){
			ans = sapply(object@filter, function(x) shape(x, transformed))
		} else{
			# assuming n.ahead=1 else will
			if(object@desc$n.ahead==1){
				ans = sapply(object@forecast, function(x) shape(x, transformed))
			} else{
				ans = lapply(object@forecast, function(x) shape(x, transformed))
			}
		}
	} else{
		if(is(object, "ACDmultifit")){
			ans = lapply(object@fit, FUN = function(x){
						if(convergence(x) == 0) shape(x, transformed) else rep(NA, length = length(x@fit$data))
					}, simplify = TRUE)
		} else if(is(object, "ACDmultifilter")){
			ans = lapply(object@filter, function(x) shape(x, transformed))
		} else{
			ans = lapply(object@forecast, function(x) shape(x, transformed))
		}
	}
	return(ans)
}

setMethod("shape", signature(object = "ACDfit"), .acdshape)
setMethod("shape", signature(object = "ACDfilter"), .acdshape)
setMethod("shape", signature(object = "ACDforecast"), .acdshape)
setMethod("shape", signature(object = "ACDpath"), .acdshape)
setMethod("shape", signature(object = "ACDsim"), .acdshape)
setMethod("shape", signature(object = "ACDroll"), .acdshape)


setMethod("shape", signature(object = "ACDmultifit"), .acdmultishape)
setMethod("shape", signature(object = "ACDmultifilter"), .acdmultishape)
setMethod("shape", signature(object = "ACDmultiforecast"), .acdmultishape)

#-------------------------------------------------------------------------------
# conditional skewness
skewness = function(object, ...)
{
	UseMethod("skewness")
}

.acdskewness = function(object, ...)
{
	dist = switch(class(object),
			ACDroll = object@model$spec@model$dmodel$model,
			object@model$dmodel$model)
	if(is(object, "ACDfit") | is(object, "ACDfilter") | is(object, "ACDroll")){
		if(dist=="ghyp"){
			if(is(object, "ACDroll")) lambda = object@forecast$density[,5] else lambda = coef(object)["ghlambda"]
		} else{
			lambda = 0
		}
		S = dskewness(distribution = dist, skew = as.numeric(skew(object)), shape = as.numeric(shape(object)), lambda = lambda)
		if(class(object)[1] == "ACDfit" | class(object)[1] == "ACDfilter"){
			D = object@model$modeldata$index[1:object@model$modeldata$T]
			S = xts(S, D)
		} else{
			S = xts(S, as.POSIXct(rownames(object@forecast$density)))
		}
	} else {
		skew = skew(object)
		shape = shape(object)
		lambda = object@model$pars["ghlambda",1]
		m = NCOL(skew)
		n = NROW(skew)
		S = matrix(NA, ncol = m, nrow = n)
		for(i in 1:m){ S[,i] = dskewness(distribution = dist, skew = skew[,i], shape = shape[,i], lambda = lambda) }
	}
	return(S)
}

setMethod("skewness", signature(object = "ACDfit"), .acdskewness)
setMethod("skewness", signature(object = "ACDfilter"), .acdskewness)
setMethod("skewness", signature(object = "ACDforecast"), .acdskewness)
setMethod("skewness", signature(object = "ACDpath"), .acdskewness)
setMethod("skewness", signature(object = "ACDsim"), .acdskewness)
setMethod("skewness", signature(object = "ACDroll"), .acdskewness)
#-------------------------------------------------------------------------------
# conditional excess kurtosis
kurtosis = function(object, ...)
{
	UseMethod("kurtosis")
}

.acdkurtosis = function(object, ...)
{
	dist = switch(class(object),
			ACDroll = object@model$spec@model$dmodel$model,
			object@model$dmodel$model)
	if(is(object, "ACDfit") | is(object, "ACDfilter") | is(object, "ACDroll")){
		if(dist=="ghyp"){
			if(is(object, "ACDroll")) lambda = object@forecast$density[,5] else lambda = coef(object)["ghlambda"]
		} else{
			lambda = 0
		}
		K = dkurtosis(distribution = dist, skew = as.numeric(skew(object)), shape = as.numeric(shape(object)), lambda=lambda )
		if(class(object)[1] == "ACDfit" | class(object)[1] == "ACDfilter"){
			D = object@model$modeldata$index[1:object@model$modeldata$T]
			K = xts(K, D)
		} else{
			K = xts(K, as.POSIXct(rownames(object@forecast$density)))
		}
	} else {
		skew = skew(object)
		shape = shape(object)
		lambda = object@model$pars["ghlambda",1]
		m = NCOL(skew)
		n = NROW(skew)
		K = matrix(NA, ncol = m, nrow = n)
		for(i in 1:m){ K[,i] = dkurtosis(distribution = dist, skew = skew[,i], shape = shape[,i], lambda = lambda) }
	}
	return(K)
}

setMethod("kurtosis", signature(object = "ACDfit"), .acdkurtosis)
setMethod("kurtosis", signature(object = "ACDfilter"), .acdkurtosis)
setMethod("kurtosis", signature(object = "ACDforecast"), .acdkurtosis)
setMethod("kurtosis", signature(object = "ACDpath"), .acdkurtosis)
setMethod("kurtosis", signature(object = "ACDsim"), .acdkurtosis)
setMethod("kurtosis", signature(object = "ACDroll"), .acdkurtosis)
#-------------------------------------------------------------------------------
# likelihood method
.acdfitLikelihood = function(object)
{
	if(is(object, "ACDfit")){
		return(c("ACD"=object@fit$LLH, "GARCH" = object@model$garchLL))
	} else if(is(object, "ACDfilter")){
		return(c("ACD"=object@filter$LLH, "GARCH" = object@model$garchLL))
	} else{
		if(!is.null(object@model$noncidx)) stop("\nroll object contains non-converged windows...use resume method first.\n")
		tmp = sapply(object@model$LL, function(x) x$log.lik)
		colnames(tmp) = sapply(object@model$LL, function(x) as.character(x$date))
		rownames(tmp) = c("ACD", "GARCH")
		return(tmp)
	}
}

.acdmultifitLikelihood = function(object)
{
	if(is(object, "ACDmultifit")){
		return(sapply(object@fit, function(x) likelihood(x)))
	} else{
		return(sapply(object@filter, function(x) likelihood(x)))
	}
}
setMethod("likelihood", signature(object = "ACDfit"), .acdfitLikelihood)
setMethod("likelihood", signature(object = "ACDfilter"), .acdfitLikelihood)

setMethod("likelihood", signature(object = "ACDmultifit"), .acdmultifitLikelihood)
setMethod("likelihood", signature(object = "ACDmultifilter"), .acdmultifitLikelihood)

#-------------------------------------------------------------------------------
# infocriteria method
.acdinfocriteria = function(object)
{
	if(is(object, "ACDfilter")){
		# np = sum(object@filter$ipars[,2])
		# all parameters fixed
		acdnp = 0
	} else{
		garchnp = sum(object@fit$ipars[1:21,4])
		acdnp = sum(object@fit$ipars[,4])
	}
	acditest = rugarch:::.information.test(likelihood(object)[1], nObs = length(fitted(object)), 
			nPars = acdnp)
	garchitest = rugarch:::.information.test(likelihood(object)[2], nObs = length(fitted(object)), 
			nPars = garchnp)
	itestm = matrix(0, ncol = 2, nrow = 4)
	itestm[1,1] = acditest$AIC
	itestm[2,1] = acditest$BIC
	itestm[3,1] = acditest$SIC
	itestm[4,1] = acditest$HQIC
	
	itestm[1,2] = garchitest$AIC
	itestm[2,2] = garchitest$BIC
	itestm[3,2] = garchitest$SIC
	itestm[4,2] = garchitest$HQIC
	colnames(itestm) = c("ACD", "GARCH")
	rownames(itestm) = c("Akaike", "Bayes", "Shibata", "Hannan-Quinn")
	return(itestm)
}

setMethod("infocriteria", signature(object = "ACDfit"), .acdinfocriteria)
#-------------------------------------------------------------------------------
# show methods

setMethod("show",
		signature(object = "ACDspec"),
		function(object){
			vmodel = object@model$vmodel$model
			model = object@model
			modelinc = object@model$modelinc
			cat(paste("\n*---------------------------------*", sep = ""))
			cat(paste("\n*          ACD Model Spec         *", sep = ""))
			cat(paste("\n*---------------------------------*", sep = ""))
			cat("\n\nConditional Variance Dynamics \t")
			cat(paste("\n-----------------------------------", sep = ""))
			cat(paste("\nGARCH Model\t: ", vmodel,"(", modelinc[10], ",", modelinc[11], ")\n", sep=""))
			if(vmodel == "fGARCH"){
				cat(paste("fGARCH Sub-Model\t: ", model$vmodel$vsubmodel, "\n", sep = ""))
			}
			cat("Mean Model\t: ARFIMA(", modelinc[2],",",ifelse(modelinc[4]>0, "d", 0),",",modelinc[3],")\n", sep = "")
			cat("Distribution\t:", model$dmodel$model,"\n")
			if(sum(object@model$dmodel$skewOrder)>0){
				cat("\nConditional Skew Dynamics \t")
				cat(paste("\n-----------------------------------", sep = ""))
				cat(paste("\nACD Skew Model\t: ", object@model$dmodel$skewmodel,"(", modelinc[22], ",", modelinc[23], ",", modelinc[24],")", sep=""))
				cat(paste("\nShock type:\t", c("Std.Residuals", "Residuals")[object@model$dmodel$skewshock],"\n",sep=""))
			}
			if(sum(object@model$dmodel$shapeOrder)>0){
				cat("\nConditional Shape Dynamics \t")
				cat(paste("\n-----------------------------------", sep = ""))
				cat(paste("\nACD Shape Model\t: ", object@model$dmodel$shapemodel,"(", modelinc[28], ",", modelinc[29], ",", modelinc[30],")", sep=""))
				cat(paste("\nShock type:\t", c("Std.Residuals", "Residuals")[object@model$dmodel$shapeshock],sep=""))
			}
			cat("\n")
			invisible(object)
		})	

setMethod("show",
		signature(object = "ACDfit"),
		function(object){
			vmodel = object@model$vmodel$model
			model = object@model
			modelinc = object@model$modelinc
			cat(paste("\n*---------------------------------*", sep = ""))
			cat(paste("\n*          ACD Model Fit          *", sep = ""))
			cat(paste("\n*---------------------------------*", sep = ""))
			cat("\n\nConditional Variance Dynamics \t")
			cat(paste("\n-----------------------------------", sep = ""))
			cat(paste("\nGARCH Model\t: ", vmodel,"(", modelinc[10], ",", modelinc[11], ")\n", sep=""))
			if(vmodel == "fGARCH"){
				cat(paste("fGARCH Sub-Model\t: ", model$vmodel$vsubmodel, "\n", sep = ""))
			}
			cat("Mean Model\t: ARFIMA(", modelinc[2],",",ifelse(modelinc[4]>0, "d", 0),",",modelinc[3],")\n", sep = "")
			cat("Distribution\t:", model$dmodel$model,"\n")
			if(sum(object@model$dmodel$skewOrder)>0){
				cat("\nConditional Skew Dynamics \t")
				cat(paste("\n-----------------------------------", sep = ""))
				cat(paste("\nACD Skew Model\t: ", object@model$dmodel$skewmodel,"(", modelinc[22], ",", modelinc[23], ",", modelinc[24],")", sep=""))
				cat(paste("\nShock type:\t", c("Std.Residuals", "Residuals")[object@model$dmodel$skewshock],"\n",sep=""))
			}
			if(sum(object@model$dmodel$shapeOrder)>0){
				cat("\nConditional Shape Dynamics \t")
				cat(paste("\n-----------------------------------", sep = ""))
				cat(paste("\nACD Shape Model\t: ", object@model$dmodel$shapemodel,"(", modelinc[28], ",", modelinc[29], ",", modelinc[30],")", sep=""))
				cat(paste("\nShock type:\t", c("Std.Residuals", "Residuals")[object@model$dmodel$shapeshock],sep=""))
			}
			if(object@fit$convergence == 0){
				cat("\n\nOptimal Parameters")
				cat(paste("\n------------------------------------\n",sep=""))
				print(round(object@fit$matcoef,6), digits = 5)
				cat("\nRobust Standard Errors:\n")
				print(round(object@fit$robust.matcoef,6), digits = 5)
				if(!is.null(object@fit$hessian.message)){
					cat(paste("\n", object@fit$hessian.message))
				}
				cat("\nLogLikelihood :", unname(object@fit$LLH[1]), "\n")
				stdresid = object@fit$residuals/object@fit$sigma
				itestm = infocriteria(object)
				cat("\nInformation Criteria")
				cat(paste("\n------------------------------------\n",sep=""))
				print(itestm,digits=5)
				cat("\nQ-Statistics on Standardized Residuals")
				cat(paste("\n------------------------------------\n",sep=""))
				tmp1 = rugarch:::.box.test(stdresid, p = 1, df = sum(modelinc[2:3]))
				print(tmp1, digits = 4)
				cat(paste("d.o.f=", sum(modelinc[2:3]), sep=""))
				cat("\nH0 : No serial correlation\n")
				cat("\nQ-Statistics on Standardized Squared Residuals")
				cat(paste("\n------------------------------------\n",sep=""))
				tmp2 = rugarch:::.box.test(stdresid, p = 2, df = sum(modelinc[10:11]))
				print(tmp2, digits = 4)
				cat(paste("d.o.f=", sum(modelinc[10:11]), sep=""))
				cat("\n\nARCH LM Tests")
				cat(paste("\n------------------------------------\n",sep=""))
				L2 = rugarch:::.archlmtest(stdresid, lags = 2)
				L5 = rugarch:::.archlmtest(stdresid, lags = 5)
				L10 = rugarch:::.archlmtest(stdresid, lags = 10)
				alm = matrix(0,ncol = 3,nrow = 3)
				alm[1,1:3] = c(L2$statistic, L2$parameter, L2$p.value)
				alm[2,1:3] = c(L5$statistic, L5$parameter, L5$p.value)
				alm[3,1:3] = c(L10$statistic, L10$parameter, L10$p.value)
				colnames(alm) = c("Statistic", "DoF", "P-Value")
				rownames(alm) = c("ARCH Lag[2]", "ARCH Lag[5]", "ARCH Lag[10]")
				print(alm,digits = 4)
				nyb = rugarch:::.nyblomTest(object)
				if(is.character(nyb$JointCritical)){
					colnames(nyb$IndividualStat)<-""
					cat("\nNyblom stability test")
					cat(paste("\n------------------------------------\n",sep=""))
					cat("Joint Statistic: ", "no.parameters>20 (not available)")
					cat("\nIndividual Statistics:")
					print(nyb$IndividualStat, digits = 4)
					cat("\nAsymptotic Critical Values (10% 5% 1%)")
					cat("\nIndividual Statistic:\t", round(nyb$IndividualCritical, 2))
					cat("\n\n")
				} else{
					colnames(nyb$IndividualStat)<-""
					cat("\nNyblom stability test")
					cat(paste("\n------------------------------------\n",sep=""))
					cat("Joint Statistic: ", round(nyb$JointStat,4))
					cat("\nIndividual Statistics:")
					print(nyb$IndividualStat, digits = 4)
					cat("\nAsymptotic Critical Values (10% 5% 1%)")
					cat("\nJoint Statistic:     \t", round(nyb$JointCritical, 3))
					cat("\nIndividual Statistic:\t", round(nyb$IndividualCritical, 2))
					cat("\n\n")
				}
				#cat("Sign Bias Test")
				#cat(paste("\n------------------------------------\n",sep=""))
				#sgtest = signbias(object)
				#print(sgtest, digits = 4)
				#cat("\n")
				#cat("\nAdjusted Pearson Goodness-of-Fit Test:")
				#cat(paste("\n------------------------------------\n",sep=""))
				#gofm = gof(object,c(20, 30, 40, 50))
				#print(gofm, digits = 4)
				#cat("\n")
				cat("\nElapsed time :", object@fit$timer,"\n\n")
			} else{
				cat("\nConvergence Problem:")
				cat("\nSolver Message:", object@fit$message,"\n\n")
				
			}
			invisible(object)
		})
setMethod("show",
		signature(object = "ACDfilter"),
		function(object){
			vmodel = object@model$vmodel$model
			model = object@model
			modelinc = object@model$modelinc
			cat(paste("\n*---------------------------------*", sep = ""))
			cat(paste("\n*          ACD Model Filter       *", sep = ""))
			cat(paste("\n*---------------------------------*", sep = ""))
			cat("\n\nConditional Variance Dynamics \t")
			cat(paste("\n-----------------------------------", sep = ""))
			cat(paste("\nGARCH Model\t: ", vmodel,"(", modelinc[10], ",", modelinc[11], ")\n", sep=""))
			if(vmodel == "fGARCH"){
				cat(paste("fGARCH Sub-Model\t: ", model$vmodel$vsubmodel, "\n", sep = ""))
			}
			cat("Mean Model\t: ARFIMA(", modelinc[2],",",ifelse(modelinc[4]>0, "d", 0),",",modelinc[3],")\n", sep = "")
			cat("Distribution\t:", model$dmodel$model,"\n")
			if(sum(object@model$dmodel$skewOrder)>0){
				cat("\nConditional Skew Dynamics \t")
				cat(paste("\n-----------------------------------", sep = ""))
				cat(paste("\nACD Skew Model\t: ", object@model$dmodel$skewmodel,"(", modelinc[22], ",", modelinc[23], ",", modelinc[24],")", sep=""))
				cat(paste("\nShock type:\t", c("Std.Residuals", "Residuals")[object@model$dmodel$skewshock],"\n",sep=""))
			}
			if(sum(object@model$dmodel$shapeOrder)>0){
				cat("\nConditional Shape Dynamics \t")
				cat(paste("\n-----------------------------------", sep = ""))
				cat(paste("\nACD Shape Model\t: ", object@model$dmodel$shapemodel,"(", modelinc[28], ",", modelinc[29], ",", modelinc[30],")", sep=""))
				cat(paste("\nShock type:\t", c("Std.Residuals", "Residuals")[object@model$dmodel$shapeshock],sep=""))
			}
			cat("\n\nFilter Parameters")
			cat(paste("\n------------------------------------\n",sep=""))
			print(matrix(coef(object), ncol=1, dimnames = list(names(coef(object)), "")), digits = 5)
			cat("\nLogLikelihood :", unname(likelihood(object)), "\n")
			stdresid = object@filter$residuals/object@filter$sigma
			cat("\nQ-Statistics on Standardized Residuals")
			cat(paste("\n---------------------------------------\n",sep=""))
			tmp1 = rugarch:::.box.test(stdresid, p = 1, df = sum(modelinc[2:3]))
			print(tmp1, digits = 4)
			cat(paste("d.o.f=", sum(modelinc[2:3]), sep=""))
			cat("\nH0 : No serial correlation\n")
			cat("\nQ-Statistics on Standardized Squared Residuals")
			cat(paste("\n---------------------------------------\n",sep=""))
			tmp2 = rugarch:::.box.test(stdresid, p = 2, df = sum(modelinc[8:9]))
			print(tmp2, digits = 4)
			cat(paste("d.o.f=", sum(modelinc[8:9]), sep=""))
			cat("\n\nARCH LM Tests")
			cat(paste("\n---------------------------------------\n",sep=""))
			L2 = rugarch:::.archlmtest(stdresid, lags = 2)
			L5 = rugarch:::.archlmtest(stdresid, lags = 5)
			L10 = rugarch:::.archlmtest(stdresid, lags = 10)
			alm = matrix(0,ncol = 3,nrow = 3)
			alm[1,1:3] = c(L2$statistic, L2$parameter, L2$p.value)
			alm[2,1:3] = c(L5$statistic, L5$parameter, L5$p.value)
			alm[3,1:3] = c(L10$statistic, L10$parameter, L10$p.value)
			colnames(alm) = c("Statistic", "DoF", "P-Value")
			rownames(alm) = c("ARCH Lag[2]", "ARCH Lag[5]", "ARCH Lag[10]")
			print(alm,digits = 4)
			cat("\nElapsed time :", object@filter$timer,"\n\n")
			invisible(object)
		})

setMethod("show",
		signature(object = "ACDforecast"),
		function(object){
			vmodel = object@model$vmodel$model
			model = object@model
			cat(paste("\n*------------------------------------*", sep = ""))
			cat(paste("\n*       ACD Model Forecast           *", sep = ""))
			cat(paste("\n*------------------------------------*", sep = ""))
			cat(paste("\nGARCH Model  : ", vmodel, sep = ""))
			n.ahead = object@forecast$n.ahead
			cat(paste("\nHorizon      : ", n.ahead, sep = ""))
			cat(paste("\nRoll Steps   : ", object@forecast$n.roll, sep = ""))
			n.start = object@forecast$n.start
			if(n.start>0) infor = ifelse(n.ahead>n.start, n.start, n.ahead) else infor = 0
			cat(paste("\nOut of Sample: ", infor, "\n", sep = ""))
			cat(paste("\n0-roll forecast [T0=", as.character(object@model$modeldata$index[object@model$modeldata$T]), "]:\n", sep=""))
			zz = cbind(object@forecast$seriesFor[,1], object@forecast$sigmaFor[,1],
					object@forecast$tskewFor[,1], object@forecast$tshapeFor[,1])
			colnames(zz) = c("series", "sigma", "skew", "shape")
			print(zz, digits = 4)
			cat("\n\n")
		})

setMethod("show",
		signature(object = "ACDroll"),
		function(object){
			if(!is.null(object@model$noncidx)){
				cat("\nObject contains non-converged estimation windows. Use resume method to re-estimate.\n")
				invisible(object)
			} else{
				cat(paste("\n*-------------------------------------*", sep = ""))
				cat(paste("\n*              ACD Roll               *", sep = ""))
				cat(paste("\n*-------------------------------------*", sep = ""))
				N = object@model$n.refits
				gmodel = object@model$spec@model$vmodel$model
				model = object@model$spec@model
				modelinc = object@model$spec@model$modelinc
				cat("\nNo.Refits\t\t:", N)
				cat("\nRefit Horizon\t:", object@model$refit.every)
				cat("\nNo.Forecasts\t:", NROW(object@forecast$density))
				cat(paste("\nGARCH Model\t\t: ", gmodel, "(",modelinc[10],",",modelinc[11],")\n", sep = ""))
				cat("Mean Model\t\t: ARFIMA(", modelinc[2],",",ifelse(modelinc[4]>0, "d", 0),",",modelinc[3],")\n", sep = "")
				if(sum(model$dmodel$skewOrder)>0){
					cat(paste("\nACD Skew Model\t: ", model$dmodel$skewmodel,"(", modelinc[22], ",", modelinc[23], ",", modelinc[24],")", sep=""))
				}
				if(sum(model$dmodel$shapeOrder)>0){
					cat(paste("\nACD Shape Model\t: ", model$dmodel$shapemodel,"(", modelinc[28], ",", modelinc[29], ",", modelinc[30],")", sep=""))
				}
				cat("\nDistribution\t:", model$dmodel$model,"\n")
				cat("\nForecast Density:\n")
				print(round(head(object@forecast$density),4))
				cat("\n..........................\n")
				print(round(tail(object@forecast$density),4))
				cat("\nElapsed:", format(object@model$elapsed))
				cat("\n")				
				invisible(object)
			}
		})
#----------------------------------------------------------------------------------
# plot methods
setMethod(f = "plot", signature(x = "ACDfit", y = "missing"), .plotacdfit)
setMethod(f = "plot", signature(x = "ACDfilter", y = "missing"), .plotacdfit)