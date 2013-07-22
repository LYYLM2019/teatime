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
#-------------------------------------------------------------------------------
# Package Highest Level Virtual Class
#-------------------------------------------------------------------------------
setClass("rACD","VIRTUAL")
#-------------------------------------------------------------------------------
setClass("ACDspec",
		representation(model = "vector"),
		contains = "rACD")
#-------------------------------------------------------------------------------
setClass("ACDfit",
		representation(fit = "vector",
				model = "vector"),
		contains = "rACD")
#-------------------------------------------------------------------------------
setClass("ACDfilter",
		representation(filter = "vector",
				model = "vector"),
		contains = "rACD")
#-------------------------------------------------------------------------------
setClass("ACDforecast",
		representation(forecast = "vector",
				model = "vector"),
		contains = "rACD")
#-------------------------------------------------------------------------------
setClass("ACDsim",
		representation(simulation = "vector",
				model = "vector",
				seed = "integer"),
		contains = "rACD")
#-------------------------------------------------------------------------------
setClass("ACDpath",
		representation(path = "vector",
				model = "vector",
				seed = "integer"),
		contains = "rACD")

setClass("ACDroll",
		representation(
				model = "vector",
				forecast = "vector"),
		contains = "rACD")

#-------------------------------------------------------------------------------
# multiple spec/fit/filter/forecast garch methods (used in 2-stage extension 
# multivariate models)
#-------------------------------------------------------------------------------

# Multiple Spec List Class
setClass("ACDmultispec", 
		representation(spec = "vector",
				type = "character"))

.validspeclist = function(object){
	all(unlist(lapply(object@spec, FUN = function(x) is(x, "ACDspec"))))
}

setValidity("ACDmultispec", .validspeclist)

# Multiple Fit ACD List Class
setClass("ACDmultifit", 
		representation(fit = "vector",
				desc = "vector"))

.validfitlist = function(object){
	all(unlist(lapply(object@fit, FUN = function(x) is(x, "ACDfit"))))
}

setValidity("ACDmultifit", .validfitlist)

# Multiple Fit ACD List Class
setClass("ACDmultifilter", 
		representation(filter = "vector",
				desc = "vector"))

.validfilterlist = function(object){
	all(unlist(lapply(object@filter, FUN = function(x) is(x, "ACDfilter"))))
}

setValidity("ACDmultifilter", .validfilterlist)

# Multiple Forecast ACD List Class
setClass("ACDmultiforecast", 
		representation(forecast = "vector",
				desc = "vector"))

.validforecastlist = function(object){
	all(unlist(lapply(object@forecast, FUN = function(x) is(x, "ACDforecast"))))
}

setValidity("ACDmultiforecast", .validforecastlist)
