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
# include "racd.h"
# include "filters.h"
# include "acd.h"
# include "distributions.h"

void arfimaxacdfilterC(int *model, double *pars, int *idx, double *x, double *res, double *mexdata,
		double *zrf, double *constm, double *condm, double *h, double *tskew, double *tshape,
		int *m, int *T)
{
	int i;
	double hm = 0;
	double sk = 0;
	double sh = 0;
	for(i=0; i<*T; i++)
	{
		hm = h[i];
		sk = tskew[i];
		sh = tshape[i];
		arfimaxfilteracd(model, pars, idx, x, res, mexdata, zrf, constm, condm,
				hm, sk, sh, *m, i, *T);
	}
}

void sacdfilterC(int *model, double *pars, int *idx, double *hEst, double *x,
		double *res, double *e, double *mexdata, double *vexdata, double *zrf,
		double *constm, double *condm, int *m, int *T, double *h, double *z,
		double *tempskew, double *tempshape, double *skhEst, double *tskew,
		double *tshape, double *skxreg, double *shxreg, double *sbounds,
		double *llh, double *LHT)
{

	int i;
	double lk=0;
	double hm = 0;
	for(i=0; i<*m; i++)
	{
		if(model[20]>0)
		{
			tempskew[i] = skhEst[0];
			tskew[i] = logmap1d(sbounds[0], sbounds[1], tempskew[i]);
		} else{
			tskew[i] = pars[idx[17]];
		}
		if(model[26]>0)
		{
			tempshape[i] = skhEst[1];
			tshape[i] = expmap1d(sbounds[2], sbounds[3], tempshape[i], sbounds[4]);
		} else{
			tshape[i] = pars[idx[18]];
		}
		h[i] = *hEst;
		arfimaxfilteracd(model, pars, idx, x, res, mexdata, zrf, constm, condm, sqrt(fabs(*hEst)), tskew[i], tshape[i], *m, i, *T);
		e[i] = res[i] * res[i];
		LHT[i] = log(garchdistribution(z[i], sqrt(fabs(h[i])), tskew[i], tshape[i], pars[idx[19]], model[40]));
		lk = lk - LHT[i];
	}
	for (i=*m; i<*T; i++)
	{
		sgarchfilteracd(model, pars, idx, vexdata, e, *T, i, h);
		hm = sqrt(fabs(h[i]));
		if(model[34]==1){
			acdskewfilter(model, pars, idx, z, tempskew, tskew, sbounds, skxreg, h, i, *T);
		} else{
			acdskewfilter(model, pars, idx, res, tempskew, tskew, sbounds, skxreg, h, i, *T);
		}
		if(model[35]==1){
			acdshapefilter(model, pars, idx, z, tempshape, tshape, sbounds, shxreg, h, i, *T);
		} else{
			acdshapefilter(model, pars, idx, res, tempshape, tshape, sbounds, shxreg, h, i, *T);
		}
		arfimaxfilteracd(model, pars, idx, x, res, mexdata, zrf, constm, condm, hm, tskew[i], tshape[i], *m, i, *T);
		e[i] = res[i] * res[i];
		z[i] = res[i]/sqrt(fabs(h[i]));
		LHT[i] = log(garchdistribution(z[i], sqrt(fabs(h[i])), tskew[i], tshape[i], pars[idx[19]], model[40]));
		lk = lk - LHT[i];
	}
	*llh=lk;
}

void csacdfilterC(int *model, double *pars, int *idx, double *hEst, double *x,
		double *res, double *e, double *mexdata, double *vexdata, double *zrf,
		double *constm, double *condm, int *m, int *T, double *h, double *q,
		double *z, double *tempskew, double *tempshape, double *skhEst,
		double *tskew, double *tshape, double *skxreg, double *shxreg,
		double *sbounds, double *llh, double *LHT)
{

	int i;
	double lk=0;
	double hm = 0;
	for(i=0; i<*m; i++)
	{
		if(model[20]>0)
		{
			tempskew[i] = skhEst[0];
			tskew[i] = logmap1d(sbounds[0], sbounds[1], tempskew[i]);
		} else{
			tskew[i] = pars[idx[17]];
		}
		if(model[26]>0)
		{
			tempshape[i] = skhEst[1];
			tshape[i] = expmap1d(sbounds[2], sbounds[3], tempshape[i], sbounds[4]);
		} else{
			tshape[i] = pars[idx[18]];
		}
		h[i] = *hEst;
		q[i] = pars[idx[8]]/(1.0-pars[idx[12]]);
		h[i] = h[i] + q[i];
		arfimaxfilteracd(model, pars, idx, x, res, mexdata, zrf, constm, condm, sqrt(fabs(*hEst)), tskew[i], tshape[i], *m, i, *T);
		e[i] = res[i] * res[i];
		LHT[i] = log(garchdistribution(z[i], sqrt(fabs(h[i])), tskew[i], tshape[i], pars[idx[19]], model[40]));
		lk = lk - LHT[i];
	}
	for (i=*m; i<*T; i++)
	{
		csgarchfilteracd(model, pars, idx, e, vexdata, *T, i, h, q);
		hm = sqrt(fabs(h[i]));
		if(model[34]==1){
			acdskewfilter(model, pars, idx, z, tempskew, tskew, sbounds, skxreg, h, i, *T);
		} else{
			acdskewfilter(model, pars, idx, res, tempskew, tskew, sbounds, skxreg, h, i, *T);
		}
		if(model[35]==1){
			acdshapefilter(model, pars, idx, z, tempshape, tshape, sbounds, shxreg, h, i, *T);
		} else{
			acdshapefilter(model, pars, idx, res, tempshape, tshape, sbounds, shxreg, h, i, *T);
		}
		arfimaxfilteracd(model, pars, idx, x, res, mexdata, zrf, constm, condm, hm, tskew[i], tshape[i], *m, i, *T);
		e[i] = res[i] * res[i];
		z[i] = res[i]/sqrt(fabs(h[i]));
		LHT[i] = log(garchdistribution(z[i], sqrt(fabs(h[i])), tskew[i], tshape[i], pars[idx[19]], model[40]));
		lk = lk - LHT[i];
	}
	*llh=lk;
}

void mcsacdfilterC(int *model, double *pars, int *idx, double *hEst, double *res, double *e,
		double *eres, double *s, double *v, double *vexdata, int *m, int *T, double *h, double *z,
		double *tempskew, double *tempshape, double *skhEst, double *tskew, double *tshape,
		double *skxreg, double *shxreg, double *sbounds, double *llh, double *LHT)
{
	int i;
	double lk=0;
	double hm = 0;
	for(i=0; i<*m; i++)
	{
		if(model[20]>0)
		{
			tempskew[i] = skhEst[0];
			tskew[i] = logmap1d(sbounds[0], sbounds[1], tempskew[i]);
		} else{
			tskew[i] = pars[idx[17]];
		}
		if(model[26]>0)
		{
			tempshape[i] = skhEst[1];
			tshape[i] = expmap1d(sbounds[2], sbounds[3], tempshape[i], sbounds[4]);
		} else{
			tshape[i] = pars[idx[18]];
		}
		h[i] = *hEst;
		hm = sqrt(fabs(h[i])*s[i]*v[i]);
		LHT[i] = log(garchdistribution(z[i], hm, tskew[i], tshape[i], pars[idx[19]], model[40]));
		lk = lk - LHT[i];
	}
	for(i=*m; i<*T; i++)
	{
		sgarchfilteracd(model, pars, idx, vexdata, e, *T, i, h);
		hm = sqrt(fabs(h[i])*s[i]*v[i]);
		if(model[34]==1){
			acdskewfilter(model, pars, idx, z, tempskew, tskew, sbounds, skxreg, h, i, *T);
		} else{
			acdskewfilter(model, pars, idx, res, tempskew, tskew, sbounds, skxreg, h, i, *T);
		}
		if(model[35]==1){
			acdshapefilter(model, pars, idx, z, tempshape, tshape, sbounds, shxreg, h, i, *T);
		} else{
			acdshapefilter(model, pars, idx, res, tempshape, tshape, sbounds, shxreg, h, i, *T);
		}
		// equivalently: res[i]/hm;
		z[i] = eres[i]/sqrt(fabs(h[i]));
		LHT[i] = log(garchdistribution(z[i], hm, tskew[i], tshape[i], pars[idx[19]], model[40]));
		lk = lk - LHT[i];
	}
	*llh=lk;
}

void sacdsimC(int *model, double *pars, int *idx, double *h, double *z,
		double *res, double *e, double *tempskew, double *tempshape,
		double *tskew, double *tshape, double *sbounds, double *vexdata,
		double *skxreg, double *shxreg, int *T, int *m)
{
	int i;
	// set.seed() and put.seed()
	GetRNGstate();
	for(i=*m;i<*T;i++)
	{
		sgarchfilteracd(model, pars, idx, vexdata, e, *T, i, h);
		if(model[34]==1){
			acdskewfilter(model, pars, idx, z, tempskew, tskew, sbounds, skxreg, h, i, *T);
		} else{
			acdskewfilter(model, pars, idx, res, tempskew, tskew, sbounds, skxreg, h, i, *T);
		}
		if(model[35]==1){
			acdshapefilter(model, pars, idx, z, tempshape, tshape, sbounds, shxreg, h, i, *T);
		} else{
			acdshapefilter(model, pars, idx, res, tempshape, tshape, sbounds, shxreg, h, i, *T);
		}
		z[i] = rgarchdist(tshape[i], tskew[i], pars[idx[19]], model[40]);
		res[i] = sqrt(h[i])*z[i];
		e[i] = res[i]*res[i];
	}
	PutRNGstate();
}

void csacdsimC(int *model, double *pars, int *idx, double *h, double *q, double *z,
		double *res, double *e, double *tempskew, double *tempshape,
		double *tskew, double *tshape, double *sbounds, double *vexdata,
		double *skxreg, double *shxreg, int *T, int *m)
{
	int i;
	// set.seed() and put.seed()
	GetRNGstate();
	for ( i=*m; i<*T; i++ )
	{
		csgarchfilteracd(model, pars, idx, e, vexdata, *T, i, h, q);
		if(model[34]==1){
			acdskewfilter(model, pars, idx, z, tempskew, tskew, sbounds, skxreg, h, i, *T);
		} else{
			acdskewfilter(model, pars, idx, res, tempskew, tskew, sbounds, skxreg, h, i, *T);
		}
		if(model[35]==1){
			acdshapefilter(model, pars, idx, z, tempshape, tshape, sbounds, shxreg, h, i, *T);
		} else{
			acdshapefilter(model, pars, idx, res, tempshape, tshape, sbounds, shxreg, h, i, *T);
		}
		z[i] = rgarchdist(tshape[i], tskew[i], pars[idx[19]], model[40]);
		res[i] = sqrt(h[i])*z[i];
		e[i] = res[i]*res[i];
	}
	PutRNGstate();
}

void mcsacdsimC(int *model, double *pars, int *idx, double *h, double *s, double *v,
		double *z, double *res, double *eres, double *e,
		double *tempskew, double *tempshape, double *tskew, double *tshape,
		double *sbounds, double *vexdata, double *skxreg, double *shxreg,
		int *T, int *m)
{
	int i;
	// set.seed() and put.seed()
	GetRNGstate();
	for(i=*m;i<*T;i++)
	{
		sgarchfilteracd(model, pars, idx, vexdata, e, *T, i, h);
		if(model[34]==1){
			acdskewfilter(model, pars, idx, z, tempskew, tskew, sbounds, skxreg, h, i, *T);
		} else{
			acdskewfilter(model, pars, idx, res, tempskew, tskew, sbounds, skxreg, h, i, *T);
		}
		if(model[35]==1){
			acdshapefilter(model, pars, idx, z, tempshape, tshape, sbounds, shxreg, h, i, *T);
		} else{
			acdshapefilter(model, pars, idx, res, tempshape, tshape, sbounds, shxreg, h, i, *T);
		}
		z[i] = rgarchdist(tshape[i], tskew[i], pars[idx[19]], model[40]);
		eres[i] = pow(h[i], 0.5)*z[i];
		res[i] = eres[i]*s[i]*v[i];
		e[i] = eres[i]*eres[i];
	}
	PutRNGstate();
}
