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

distbounds = function(distribution){
	distribution = tolower(distribution[1])
	distribution = match.arg(distribution, c("ged","std","snorm","sged","sstd","nig","ghyp","ghst","jsu"))
	ans = .DistributionBounds(distribution)
	if(distribution == "ghyp"){
		ans = c(ans$skew.LB, ans$skew.UB, ans$shape.LB, ans$shape.UB, ans$ghlambda.LB, ans$ghlambda.UB)
		names(ans) = c("skew(LB)", "skew(UB)", "shape(LB)", "shape(UB)","shape_GIG(LB)","shape_GIG(UB)")
	} else{
		ans = c(ans$skew.LB, ans$skew.UB, ans$shape.LB, ans$shape.UB)
		names(ans) = c("skew(LB)", "skew(UB)", "shape(LB)", "shape(UB)")
	}
	return(ans)
}

logtransform = function(x, lower, upper, inverse = FALSE){
	if(!inverse){
		ans = lower + (upper-lower)/(1+exp(-1*x))
	} else{
		ans = -1*log(-(upper-x)/(-x+lower))
	}
	return(ans)
}


exptransform = function(x, lower, upper, rate = 1, inverse = FALSE){
	if(!inverse){
		ans = lower+exp(-rate*x)*upper
	} else{
		ans =-(1/rate)*log((x - lower)/upper)
	}
	return(ans)
}

.DistributionBounds = function(distribution)
{
	ghlambda = 0
	ghlambda.LB = 0
	ghlambda.UB = 0
	if (distribution == "ged"){
		skew 	= 0
		skew.LB = 0
		skew.UB = 10
		shape 	= 2
		shape.LB = 0.1
		shape.UB = 50
	}
	if (distribution == "std"){
		skew 	= 0
		skew.LB = 0
		skew.UB = 0
		shape 	= 4
		shape.LB = 2.05
		shape.UB = 25
	}
	if (distribution == "snorm"){
		skew 	= 0.9
		skew.LB = 0.1
		skew.UB = 10
		shape 	= 0
		shape.LB = 0
		shape.UB = 0
	}
	if (distribution == "sged"){
		skew 	= 1
		skew.LB	= 0.01
		skew.UB	= 5
		shape 	= 1
		shape.LB = 0.7
		shape.UB = 4
	}
	if (distribution == "sstd"){
		skew 	= 1
		skew.LB = 0.1
		skew.UB = 30
		shape 	= 4
		shape.LB = 4.05
		shape.UB = 25
	}
	if (distribution == "nig"){
		skew 	= 0.2
		skew.LB =-0.99999
		skew.UB = 0.99999
		shape 	= 0.5
		shape.LB = 0.1
		shape.UB = 15
	}
	if(distribution == "ghyp"){
		skew 	= 0.2
		skew.LB = -0.99
		skew.UB	= 0.99
		shape 	= 0.6
		shape.LB = 0.25
		shape.UB = 25
		ghlambda = -0.5
		ghlambda.LB = -4
		ghlambda.UB = 2
	}
	if(distribution == "jsu"){
		skew 	= 0
		skew.LB =-10
		skew.UB = 10
		shape 	= 1
		shape.LB = 0.8
		shape.UB = 5
	}
	# johnson has 2 shape parameters. The second one we model with the "skew"
	# representation in rgarch
	if(distribution == "ghst"){
		skew 	= 0
		skew.LB	= -5
		skew.UB	= 5
		shape 	= 8.2
		shape.LB = 8.1
		shape.UB = 25
	}
	skewed.dists = c("snorm", "sged", "sstd", "nig", "ghyp", "jsu", "ghst")
	shaped.dists = c("ged", "sged", "std", "sstd", "nig", "ghyp", "jsu", "ghst")
	skew0  = 0
	shape0 = 0

	if(any(skewed.dists == distribution)) include.skew = TRUE else include.skew = FALSE
	if(any(shaped.dists == distribution)) include.shape = TRUE else include.shape = FALSE
	if(distribution == "ghyp") include.ghlambda = TRUE else include.ghlambda = FALSE
	
	ans = list(shape = shape, shape.LB = shape.LB, shape.UB = shape.UB, skew = skew,
			skew.LB = skew.LB, skew.UB = skew.UB, include.skew = include.skew, 
			include.shape = include.shape, skew0 = skew0, shape0 = shape0, 
			include.ghlambda = include.ghlambda, ghlambda = ghlambda, ghlambda.LB = ghlambda.LB, 
			ghlambda.UB = ghlambda.UB)
	return(ans)
}

.acdskewbounds = function(acdOrder, unconpar, distribution, dbounds)
{
	.eps = .Machine$double.eps
	skew.LB = dbounds[1]
	skew.UB = dbounds[2]
	par = par.LB = par.UB = numeric()
	# intercept
	par[1]	= .invlogtransform(unconpar, skew.LB, skew.UB)
	par.LB[1] = .invlogtransform(ifelse(skew.LB<0,skew.LB*0.9999,skew.LB*1.00001), skew.LB, skew.UB)
	par.UB[1] = .invlogtransform(skew.UB*0.9999, skew.LB, skew.UB)
	# alpha1 and alpha2 with lower/upper bounds
	if(acdOrder[1] > 0){
		par = c(par, rep(0.1, acdOrder[1]))
		par.LB = c(par.LB, rep(-5 + .eps, acdOrder[1]))
		par.UB = c(par.UB, rep( 5 - .eps, acdOrder[1]))
	}
	if(acdOrder[2] > 0){
		par = c(par, rep(0.1, acdOrder[2]))
		par.LB = c(par.LB, rep(-5 + .eps, acdOrder[2]))
		par.UB = c(par.UB, rep( 5 - .eps, acdOrder[2]))	
	}
	if(acdOrder[3] > 0){
		par = c(par, rep(0.8/acdOrder[3], acdOrder[3]))
		par.LB = c(par.LB, rep( 0 + .eps, acdOrder[3]))
		par.UB = c(par.UB, rep( 1 - .eps, acdOrder[3]))
	}
	return(list(skewpars = par, skewpar.LB = par.LB, skewpar.UB = par.UB, sk0 = par[1]))
}

.acdshapebounds = function(acdOrder, unconpar, distribution, dbounds)
{
	.eps = .Machine$double.eps
	shape.LB = dbounds[1]
	shape.UB = dbounds[2]
	par = par.LB = par.UB = numeric()
	# intercept
	par[1]	= .invexptransform(unconpar, shape.LB, shape.UB, dbounds[3])
	par.LB[1] = 0.0001
	par.UB[1] = shape.UB
	# alpha1 and alpha2 with lower/upper bounds
	if(acdOrder[1] > 0){
		par = c(par, rep(0.1, acdOrder[1]))
		par.LB = c(par.LB, rep( 0 + .eps, acdOrder[1]))
		par.UB = c(par.UB, rep( 1 - .eps, acdOrder[1]))
	}
	if(acdOrder[2] > 0){
		par = c(par, rep(0.1, acdOrder[2]))
		par.LB = c(par.LB, rep( 0 + .eps, acdOrder[2]))
		par.UB = c(par.UB, rep( 1 - .eps, acdOrder[2]))		
	}
	if(acdOrder[3] > 0){
		par = c(par, rep(0.5/acdOrder[3], acdOrder[3]))
		par.LB = c(par.LB, rep( 0 + .eps, acdOrder[3]))
		par.UB = c(par.UB, rep( 1 - .eps, acdOrder[3]))
	}
	return(list(shapepars = par, shapepar.LB = par.LB, shapepar.UB = par.UB, sh0 = par[1]))
}

# logistic transformation and inverse transformation
.invlogtransform = function(y, LB, UB)
{
	x = -1*log(-(UB-y)/(-y+LB))
	return(x)
}

.logtransform = function(x, LB, UB)
{
	y = LB + (UB-LB)/(1+exp(-1*x))
	return(y)
}

.sinetransform = function(x, LB, UB){
	z = (sin(pi*x)*(UB - LB) + UB+LB)/2
	z = pmax(LB, pmin(UB, z))
	return(z)
}

.exptransform = function(x, lower, upper, rate=1){
	lower+exp(-rate*x)*upper
}

.invexptransform = function(x, lower, upper, rate=1)
{
	-(1/rate)*log((x - lower)/upper)
}