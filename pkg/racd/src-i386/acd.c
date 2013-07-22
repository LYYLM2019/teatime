/*################################################################################
##
##   R package rarcd by Alexios Ghalanos Copyright (C) 2012,2013
##   This file is part of the R package rarcd.
##
##   The R package rarcd is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package rarcd is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################*/
# include <R.h>
# include <math.h>
# include "distributions.h"
# include "acd.h"

void acdskewfilter(int *model, double *pars, int *idx, double *z, double *tempskew, double *tskew,
		double *sbounds, double *skewreg, double *h, int i, int T)
{
	int k;
	if( model[20]>0 )
	{
		for(k=0;k<model[36];k++)
		{
			tempskew[i]+=skewdynamics(model, idx, pars, z[i-1-k], fabs(h[i-1-k]), tempskew[i-1-k], k);
		}
		if(model[24]>0)
		{
			for( k=0; k<model[24]; k++ )
			{
				tempskew[i]+=pars[idx[24]+k]*skewreg[i + ( T * k )];
			}
		}
		tskew[i] = logmap1d(sbounds[0], sbounds[1], tempskew[i]);
	}
	else
	{
		tskew[i] = pars[idx[17]];
	}
}

void acdshapefilter(int *model, double *pars, int *idx, double *z, double *tempshape, double *tshape,
		double *sbounds, double *shapereg, double *h, int i, int T)
{
	int k;
	if( model[26]>0 )
	{
		for(k=0;k<model[37];k++)
		{
			tempshape[i]+=shapedynamics(model, idx, pars, z[i-1-k], fabs(h[i-1-k]), tempshape[i-1-k], k);
		}
		if(model[30]>0)
		{
			for( k=0; k<model[30]; k++ )
			{
				tempshape[i]+=pars[idx[30]+k]*shapereg[i + ( T * k )];
			}
		}
		tshape[i] = expmap1d(sbounds[2], sbounds[3], tempshape[i], sbounds[4]);
	}
	else
	{
		tshape[i] = pars[idx[18]];
	}
}

double skewdynamics(const int *model, const int *idx, const double *pars, const double z, const double h, const double lagval, const int lag)
{
	double res=0.0;
	switch(model[32])
	{
		case 0:
		{
			// quadratic model
			res = pars[idx[20]];
			if(model[21]>lag) res+=pars[idx[21]+lag]*z;
			if(model[22]>lag) res+=pars[idx[22]+lag]*(z*z);
			if(model[23]>lag) res+=pars[idx[23]+lag]*lagval;
			break;
		}
		case 1:
		{
			// piece-wise linear threshold model (pwl)
			res = pars[idx[20]];
			double z1 = (z <  0)? z : 0;
			double z2 = (z >= 0)? z : 0;
			if(model[21]>lag) res+=pars[idx[21]+lag]*z1;
			if(model[22]>lag) res+=pars[idx[22]+lag]*z2;
			if(model[23]>lag) res+=pars[idx[23]+lag]*lagval;
			break;
		}
		case 2:
		{
			res = pars[idx[20]];
			double tmp = (z*z) - fabs(h);
			double z1 = (tmp>=0 && z<0 )?z:0;
			double z2 = (tmp>=0 && z>=0)?z:0;
			if(model[21]>lag) res+=pars[idx[21]+lag]*z1;
			if(model[22]>lag) res+=pars[idx[22]+lag]*z2;
			if(model[23]>lag) res+=pars[idx[23]+lag]*lagval;
			break;
		}
		case 3:
		{
			// quadratic model
			res = pars[idx[20]];
			if(model[21]>lag) res+=pars[idx[21]+lag]*z;
			if(model[22]>lag) res+=pars[idx[22]+lag]*fabs(z);
			if(model[23]>lag) res+=pars[idx[23]+lag]*lagval;
			break;
		}
		case 4:
		{
			res = pars[idx[20]];
			double tmp = fabs(z) - sqrt(fabs(h));
			double z1 = (tmp>=0 && z<0 )?z:0;
			double z2 = (tmp>=0 && z>=0)?z:0;
			if(model[21]>lag) res+=pars[idx[21]+lag]*z1;
			if(model[22]>lag) res+=pars[idx[22]+lag]*z2;
			if(model[23]>lag) res+=pars[idx[23]+lag]*lagval;
			break;
		}
		case 5:
		{
			// tar model
			res = pars[idx[20]];
			double z1 = (z <  pars[idx[25]])? z : 0;
			double z2 = (z >= pars[idx[25]])? z : 0;
			if(model[21]>lag) res+=pars[idx[21]+lag]*z1;
			if(model[22]>lag) res+=pars[idx[22]+lag]*z2;
			if(model[23]>lag) res+=pars[idx[23]+lag]*lagval;
			break;
		}
	}
	return res;
}

double shapedynamics(const int *model, const int *idx, const double *pars, const double z, const double h, const double lagval, const int lag)
{
	double res=0.0;
	switch(model[33])
	{
		case 0:
		{
			res = pars[idx[26]];
			if(model[27]>lag) res+=pars[idx[27]+lag]*z;
			if(model[28]>lag) res+=pars[idx[28]+lag]*(z*z);
			if(model[29]>lag) res+=pars[idx[29]+lag]*lagval;
			break;
		}
		case 1:
		{
			res = pars[idx[26]];
			double z1 = (z <  0)? (z*z) : 0;
			double z2 = (z >= 0)? (z*z) : 0;
			if(model[27]>lag) res+=pars[idx[27]+lag]*z1;
			if(model[28]>lag) res+=pars[idx[28]+lag]*z2;
			if(model[29]>lag) res+=pars[idx[29]+lag]*lagval;
			break;
		}
		case 2:
		{
			res = pars[idx[26]];
			double tmp = (z*z) - fabs(h);
			double z1 = (tmp >= 0 ) ? (z*z) : 0;
			if(model[28]>lag) res+=pars[idx[28]+lag]*z1;
			if(model[29]>lag) res+=pars[idx[29]+lag]*lagval;
			break;
		}
		case 3:
		{
			res = pars[idx[26]];
			double z1 = (z <   pars[idx[31]])? (z*z) : 0;
			double z2 = (z >=  pars[idx[31]])? (z*z) : 0;
			if(model[27]>lag) res+=pars[idx[27]+lag]*z1;
			if(model[28]>lag) res+=pars[idx[28]+lag]*z2;
			if(model[29]>lag) res+=pars[idx[29]+lag]*lagval;
			break;
		}
		case 4:
		{
			res = pars[idx[26]];
			if(model[27]>lag) res+=pars[idx[27]+lag]*z;
			if(model[28]>lag) res+=pars[idx[28]+lag]*fabs(z);
			if(model[29]>lag) res+=pars[idx[29]+lag]*lagval;
			break;
		}
		case 5:
		{
			res = pars[idx[26]];
			double z1 = (z < 0)? fabs(z) : 0;
			double z2 = (z >=0)? fabs(z) : 0;
			if(model[27]>lag) res+=pars[idx[27]+lag]*z1;
			if(model[28]>lag) res+=pars[idx[28]+lag]*z2;
			if(model[29]>lag) res+=pars[idx[29]+lag]*lagval;
			break;
		}
		case 6:
		{
			res = pars[idx[26]];
			double tmp = fabs(z) - sqrt(fabs(h));
			double z1 = (tmp >= 0 ) ? fabs(z) : 0;
			if(model[28]>lag) res+=pars[idx[28]+lag]*z1;
			if(model[29]>lag) res+=pars[idx[29]+lag]*lagval;
			break;
		}
		case 7:
		{
			res = pars[idx[26]];
			double z1 = (z <   pars[idx[31]])? fabs(z) : 0;
			double z2 = (z >=  pars[idx[31]])? fabs(z) : 0;
			if(model[27]>lag) res+=pars[idx[27]+lag]*z1;
			if(model[28]>lag) res+=pars[idx[28]+lag]*z2;
			if(model[29]>lag) res+=pars[idx[29]+lag]*lagval;
			break;
		}
	}
	return res;
}

double logmap1d(const double LB, const double UB, const double val)
{
	double res=0.0;
	res=LB+(UB-LB)/(1.0 + exp(-1*val));
	return res;
}

double invlogmap1d(const double LB, const double UB, const double res)
{
	double val=0.0;
	val=-1*log(-(UB-res)/(-res+LB));
	return val;
}

double expmap1d(const double LB, const double UB, const double val, const double rate)
{
	double res=0.0;
	res = LB + exp(-1.0*rate*val)*UB;
	return res;
}

double invexpmap1d(const double LB, const double UB, const double res, const double rate)
{
	double val=0.0;
	val=-(1.0/rate)*log((res - LB)/UB);
	return val;
}
