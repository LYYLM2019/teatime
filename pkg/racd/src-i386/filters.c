/*################################################################################
##
##   R package rgarch by Alexios Ghalanos Copyright (C) 2009
##   This file is part of the R package rgarch.
##
##   The R package rgarch is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package rgarch is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################*/
# include <R.h>
# include <math.h>
# include "filters.h"
# include "acd.h"

void sgarchfilteracd(int *model, double *pars, int *idx, double *vexdata, double *e, int T, int i, double *h)
{
	int j;
	h[i] = h[i] + pars[idx[8]];
	if( model[16]>0 )
	{
		int ind=0;
		for( j=0; j<model[16]; j++ )
		{
			ind = i + ( T * j );
			h[i] = h[i] + pars[idx[16]+j]*vexdata[ind];
		}
	}
	for( j=0; j<model[9]; j++ )
	{
		h[i] = h[i] + pars[idx[9]+j]*e[i-(j+1)];
	}
	for( j=0; j<model[10]; j++ )
	{
		h[i] = h[i] + pars[idx[10]+j]*h[i-(j+1)];
	}
}

void csgarchfilteracd(int *model, double *pars, int *idx, double *e, double *vexdata, int T, int i, double *h, double *q)
{
	int j, ind;
	q[i] = pars[idx[8]] + pars[idx[12]]*q[i-1] + pars[idx[13]]*(e[i-1] - h[i-1]);
	// External Regressors are added to the Permanent Component
	if( model[16]>0 )
	{
		for( j=0; j<model[16]; j++ )
		{
			ind = i + ( T * j );
			q[i] = q[i] + pars[idx[16]+j]*vexdata[ind];
		}
	}
	h[i] = h[i] + q[i];
	for( j=0; j<model[9]; j++ )
	{
		h[i] = h[i] + pars[idx[9]+j]*(e[i-(j+1)] - q[i-(j+1)]);
	}
	for( j=0; j<model[10]; j++ )
	{
		h[i] = h[i] + pars[idx[10]+j]*(h[i-(j+1)] - q[i-(j+1)]);
	}
}

void arfimaxfilteracd(int* model, double *pars, int *idx, double *x, double *res,
		double *mexdata, double *zrf, double *constm, double *condm, double h,
		double sk, double sh, int m, int i, int T)
{
/* --------------------------------------------------------------------------------
 * ARFIMA Process :
 * (1-L)^(-darfima).e[t] = phi(1-L)(y[t] - mu[t]) - psi(L).e[t]
 * where mu[t] includes constant, external and arch-in-mean
 * L is the lag operator
 * --------------------------------------------------------------------------------
 * */
	/*0 constm, 1 condm, 2 res*/
	int j, k;
	constm[i] = pars[0];
	// GARCH-In-Mean Initialization (h is always the sigma, not sigma^2 so that h^model[4] is correct)
	if(model[4]>0)
	{
		constm[i]+=pars[idx[4]]*pow(h, model[4]);
	}
	if(model[5]>0)
	{
		constm[i]+=pars[idx[5]]*sk;
	}
	if(model[6]>0)
	{
		constm[i]+=pars[idx[6]]*sh;
	}
	// Exogenous Regressor Initialization
	if(model[7]>0)
	{
		int ind=0;
		for(k=0;k<model[7];k++)
		{
			ind=i+(T*k);
			constm[i]+=pars[idx[7]+k]*mexdata[ind];
		}
	}
	condm[i]+=constm[i];
	//ARMA initialization
	if(model[1]>0 || model[2]>0)
	{
		if(i>=model[1])
		{
			if(model[1]>0)
			{
				for(j=0; j<model[1];j++)
				{
					condm[i]+=pars[idx[1]+j]*(x[i-(j+1)]-constm[i-(j+1)]);
				}
			}
			if(model[2]>0)
			{
				for(j=0; j<model[2];j++)
				{
					if(i-j-1>=0)
					{
						condm[i]+=pars[idx[2]+j]*(x[i-(j+1)]-condm[i-(j+1)]);
					}
				}
			}
		}
	}
	res[i]=x[i]-condm[i];
	//arfima initialization
	if(model[3]>0)
	{
		if(i>0 && i<m)
		{
			double tmp=0;
			for(k=1;k<=i;k++)
			{
				tmp+=(zrf[i-k+1]*res[k-1]) ;
			}
			res[i]=-1.0 * tmp;
		}
		if(i>0 && i>=m)
		{
			double tmp=0;
			for(k=i;k>0;k--)
			{
				// quicker to count down
				tmp+=zrf[k]*(x[i-k] - condm[i-k]) ;
			}
			res[i]+= tmp;
		}
	}
}

void armaxsimacd(int *model, double *pars, int *idx, double *x, double *res, double *constm, int *m, int *T)
{
	int j, i;
	for(i=*m; i<*T; i++)
	{
		x[i] = constm[i];
		for( j=0; j<model[1]; j++ )
		{
			x[i]+= pars[idx[1]+j] * (x[i-(j+1)] - constm[i-(j+1)]);
		}
		for ( j=0; j<model[2]; j++ )
		{
			x[i]+= pars[idx[2]+j] * res[i-(j+1)];
		}
		x[i]+= res[i];
	}
}
