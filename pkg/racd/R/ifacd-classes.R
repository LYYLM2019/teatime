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
# IFACD Model
setClass("IFACDspec",
		representation(model = "vector",
				umodel = "vector"))

setClass("IFACDfit",
		representation(mfit = "vector",
				model = "vector"))

setClass("IFACDfilter",
		representation(mfilter = "vector",
				model = "vector"))

setClass("IFACDforecast",
		representation(mforecast = "vector",
				model = "vector"))

setClass("IFACDsim",
		representation(msim = "vector",
				model = "vector"))

setClass("IFACDroll",
		representation(forecast = "vector",
				model = "vector"))

setClass("IFACDfft",
		representation(dist = "vector", model = "vector"))
