/*################################################################################
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
#################################################################################*/
# include <R.h>
# include <limits.h>
# include <math.h>
# include <Rmath.h>
# include "dbm.h"

double pburr2(const double p, const double k)
{
	double ans = pow(1.0/(1.0+exp(-1.0*p)), k);
	return ans;
}
double dlink(const double p, const double y, const int link, const double k)
{
	double tmp, ans;
	if(link==1){
		tmp = pnorm(p, 0, 1, 1, 0);
	} else if(link==2){
		tmp = plogis(p, 0, 1, 1, 0);
	} else{
		tmp = pburr2(p, k);
	}
	ans = y*log(tmp) + (1-y)*log(1-tmp);
	return ans;
}
void c_dbmestimate(double *y, double *x, double *mpu, double *mpuinit, double *omega,
		double *alpha, double *delta, double *beta, double *k,
		double *lik, double *llh, int *idx, int *xidx, int *T)
{
	int i, j, ind;
	for(i=0;i<*T;i++){
		// include intercept?
		if(idx[0]>0){
			mpu[i] += *omega;
		}
		// include arp? (limited to 1-lag)
		if(idx[2]>0){
			if(i>0){
				mpu[i] += *alpha * mpu[i-1];
			} else{
				mpu[i] += *alpha * mpuinit[0];
			}
		}
		// include arq?
		// ecm model only allowed for lag=1
		if(idx[4]>0){
			if(i>0){
				mpu[i] += (1.0 - *alpha)*y[i-1];
			}
		} else{
			if(idx[3]>0){
				for(j=0;j<idx[3];j++){
					if(i>0){
						if(i>j){
							mpu[i] += delta[j] * y[i-(j+1)];
						}
					}
				}
			}
		}
		// '<=' used since we have zero indexing
		if(idx[1]>0 && i>0){
			for(j=0;j<idx[1];j++){
				if(i>=xidx[j]){
					ind = (i-xidx[j]) + ( *T * j );
					mpu[i] += beta[j]*x[ind];
				}
			}
		}
		lik[i] = dlink(mpu[i], y[i], idx[5], k[0]);
		*llh += lik[i];
	}
	*llh *= -1.0;
}

// logistic function derivative
// also accounts for gradient wrt startup value equation
void c_dbmderiv1(double *y, double *x, double *mpu,
		double *meanx, double *meany, double *omega,
		double *alpha, double *delta, double *beta, double *domega,
		double *dalpha, double *ddelta, double *dbeta,
		double *dpomega, double *dpalpha, double *dpdelta, double *dpbeta,
		double *dvomega, double *dvalpha, double *dvdelta, double *dvbeta,
		int *idx, int *xidx, int *T)
{
	int i, j, ind;
	double xtmp = 0.0;
	double mpuinit = 0.0;
	for(i=0;i<*T;i++){
		// intercept
		if(idx[0]>0){
			mpu[i] += omega[0];
			if(i==0){
				mpuinit += omega[0];
			}
		}
		// arq
		// ecm model only allowed for lag=1
		if(idx[4]>0){
			if(i>0){
				mpu[i] += (1.0 - alpha[0])*y[i-1];
			} else{
				mpuinit += (1.0 - alpha[0]) * meany[0];
			}
		} else{
			if(idx[3]>0){
				for(j=0;j<idx[3];j++){
					if(i>0){
						if(i>j){
							mpu[i] += delta[j]*y[i-(j+1)];
						}
					} else{
						mpuinit += delta[j] * meany[0];
					}
				}
			}
		}
		// explanatory variables
		if(idx[1]>0){
			for(j=0;j<idx[1];j++){
				if(i>0){
					if(i>=xidx[j]){
						ind = (i-xidx[j]) + ( *T * j );
						mpu[i] += beta[j]*x[ind];
					}
				} else{
					mpuinit += beta[j] * meanx[j];
				}
			}
		}
		// arp
		if(idx[2]>0){
			if(i>0){
				mpu[i] += alpha[0] * mpu[i-1];
			} else{
				mpuinit *= 1.0/(1.0 - alpha[0]);
				mpu[i] += alpha[0] * mpuinit;
			}
		}
		// common equation in gradient equation
		xtmp = ((y[i]*exp(-1.0*mpu[i])+y[i]-1.0))/(1.0+exp(-1.0*mpu[i]));

		if(idx[0]>0){
			if(i>0){
				dpomega[i] = 1.0 + alpha[0] * dpomega[i-1];
			} else{
				dpomega[i] = 1.0 + alpha[0]/(1.0-alpha[0]);
			}
			*domega += xtmp * dpomega[i];
			dvomega[i] = xtmp * dpomega[i];
		}
		if(idx[2]>0){
			if(i>0){
				dpalpha[i] = mpu[i-1] + alpha[0] * dpalpha[i-1];
				dvalpha[i] = xtmp * dpalpha[i];
				*dalpha += dvalpha[i];
			} else{
				dpalpha[i] = mpuinit + (alpha[0]*mpuinit)/(1-alpha[0]);
				dvalpha[i] = xtmp * dpalpha[i];
				*dalpha += dvalpha[i];
			}
		}
		if(idx[4]==0){
			if(idx[3]>0){
				for(j=0;j<idx[3];j++){
					if(i>0){
						if(i>j){
							dpdelta[i+(*T * j)] = y[i-(j+1)] + alpha[0] * dpdelta[i+(*T * j)-1];
							dvdelta[i+(*T * j)] = xtmp*dpdelta[i+(*T * j)];
							ddelta[j] += dvdelta[i+(*T * j)];
						}
					} else{
						dpdelta[i+(*T * j)] = (alpha[0] * meany[0])/(1 - alpha[0]);
						dvdelta[i+(*T * j)] = xtmp*dpdelta[i+(*T * j)];
						ddelta[j] += dvdelta[i+(*T * j)];
					}
				}
			}
		}
		if(idx[1]>0){
			for(j=0;j<idx[1];j++){
				if(i>0){
					ind = (i-xidx[j]) + ( *T * j );
					if(i>=xidx[j]){
						dpbeta[i+(*T * j)] = x[ind] + alpha[0] * dpbeta[i+(*T * j)-1];
						dvbeta[i+(*T * j)] =  xtmp * dpbeta[i+(*T * j)];
						dbeta[j] += dvbeta[i+(*T * j)] ;
					}
				} else{
					dpbeta[i+(*T * j)] = (alpha[0] * meanx[j])/(1 - alpha[0]);
					dvbeta[i+(*T * j)] =  xtmp * dpbeta[i+(*T * j)];
					dbeta[j] += dvbeta[i+(*T * j)] ;
				}

			}
		}
	}
}

//Gaussian Link
void c_dbmderiv2(double *y, double *x, double *mpu,
		double *meanx, double *meany, double *omega,
		double *alpha, double *delta, double *beta, double *domega,
		double *dalpha, double *ddelta, double *dbeta,
		double *dpomega, double *dpalpha, double *dpdelta, double *dpbeta,
		double *dvomega, double *dvalpha, double *dvdelta, double *dvbeta,
		int *idx, int *xidx, int *T)
{
	int i, j, ind;
	double xtmp = 0.0;
	double mpuinit = 0.0;
	double numer = 0.0;
	double denom = 0.0;
	double sqrt2 = sqrt(2);
	double sqrtpi = sqrt(PI);
	for(i=0;i<*T;i++){
		// intercept
		if(idx[0]>0){
			mpu[i] += omega[0];
			if(i==0){
				mpuinit += omega[0];
			}
		}
		// arq
		// ecm model only allowed for lag=1
		if(idx[4]>0){
			if(i>0){
				mpu[i] += (1.0 - alpha[0])*y[i-1];
			} else{
				mpuinit += (1.0 - alpha[0]) * meany[0];
			}
		} else{
			if(idx[3]>0){
				for(j=0;j<idx[3];j++){
					if(i>0){
						if(i>j){
							mpu[i] += delta[j]*y[i-(j+1)];
						}
					} else{
						mpuinit += delta[j] * meany[0];
					}
				}
			}
		}
		// explanatory variables
		if(idx[1]>0){
			for(j=0;j<idx[1];j++){
				if(i>0){
					if(i>=xidx[j]){
						ind = (i-xidx[j]) + ( *T * j );
						mpu[i] += beta[j]*x[ind];
					}
				} else{
					mpuinit += beta[j] * meanx[j];
				}
			}
		}
		// arp
		if(idx[2]>0){
			if(i>0){
				mpu[i] += alpha[0] * mpu[i-1];
			} else{
				mpuinit *= 1.0/(1.0 - alpha[0]);
				mpu[i] += alpha[0] * mpuinit;
			}
		}
		// common equation in gradient equation
		numer = exp(-0.5*mpu[i]*mpu[i])*sqrt2*(erf(0.5*sqrt2*mpu[i])-2.0*y[i]+1.0);
		denom = sqrtpi*(pow(erf(-0.5*sqrt2*mpu[i]), 2.0)-1.0);
		xtmp = numer/denom;;

		if(idx[0]>0){
			if(i>0){
				dpomega[i] = 1.0 + alpha[0] * dpomega[i-1];
			} else{
				dpomega[i] = 1.0 + alpha[0]/(1.0-alpha[0]);
			}
			*domega += xtmp * dpomega[i];
			dvomega[i] = xtmp * dpomega[i];
		}
		if(idx[2]>0){
			if(i>0){
				dpalpha[i] = mpu[i-1] + alpha[0] * dpalpha[i-1];
				dvalpha[i] = xtmp * dpalpha[i];
				*dalpha += dvalpha[i];
			} else{
				dpalpha[i] = mpuinit + (alpha[0]*mpuinit)/(1-alpha[0]);
				dvalpha[i] = xtmp * dpalpha[i];
				*dalpha += dvalpha[i];
			}
		}
		if(idx[4]==0){
			if(idx[3]>0){
				for(j=0;j<idx[3];j++){
					if(i>0){
						if(i>j){
							dpdelta[i+(*T * j)] = y[i-(j+1)] + alpha[0] * dpdelta[i+(*T * j)-1];
							dvdelta[i+(*T * j)] = xtmp*dpdelta[i+(*T * j)];
							ddelta[j] += dvdelta[i+(*T * j)];
						}
					} else{
						dpdelta[i+(*T * j)] = (alpha[0] * meany[0])/(1 - alpha[0]);
						dvdelta[i+(*T * j)] = xtmp*dpdelta[i+(*T * j)];
						ddelta[j] += dvdelta[i+(*T * j)];
					}
				}
			}
		}
		if(idx[1]>0){
			for(j=0;j<idx[1];j++){
				if(i>0){
					ind = (i-xidx[j]) + ( *T * j );
					if(i>=xidx[j]){
						dpbeta[i+(*T * j)] = x[ind] + alpha[0] * dpbeta[i+(*T * j)-1];
						dvbeta[i+(*T * j)] =  xtmp * dpbeta[i+(*T * j)];
						dbeta[j] += dvbeta[i+(*T * j)] ;
					}
				} else{
					dpbeta[i+(*T * j)] = (alpha[0] * meanx[j])/(1 - alpha[0]);
					dvbeta[i+(*T * j)] =  xtmp * dpbeta[i+(*T * j)];
					dbeta[j] += dvbeta[i+(*T * j)] ;
				}

			}
		}
	}
}

// Burr Type 2 Link
void c_dbmderiv3(double *y, double *x, double *mpu, double *meanx, double *meany,
		double *omega, double *alpha, double *delta, double *beta, double *k,
		double *domega, double *dalpha, double *ddelta, double *dbeta, double *dk,
		double *dpomega, double *dpalpha, double *dpdelta, double *dpbeta, double *dpk,
		double *dvomega, double *dvalpha, double *dvdelta, double *dvbeta, double *dvk,
		int *idx, int *xidx, int *T)
{
	int i, j, ind;
	double xtmp = 0.0;
	double ktmp = 0.0;
	double mpuinit = 0.0;
	double denom = 0.0;
	double numer = 0.0;
	for(i=0;i<*T;i++){
		// intercept
		if(idx[0]>0){
			mpu[i] += omega[0];
			if(i==0){
				mpuinit += omega[0];
			}
		}
		// arq
		// ecm model only allowed for lag=1
		if(idx[4]>0){
			if(i>0){
				mpu[i] += (1.0 - alpha[0])*y[i-1];
			} else{
				mpuinit += (1.0 - alpha[0]) * meany[0];
			}
		} else{
			if(idx[3]>0){
				for(j=0;j<idx[3];j++){
					if(i>0){
						if(i>j){
							mpu[i] += delta[j]*y[i-(j+1)];
						}
					} else{
						mpuinit += delta[j] * meany[0];
					}
				}
			}
		}
		// explanatory variables
		if(idx[1]>0){
			for(j=0;j<idx[1];j++){
				if(i>0){
					if(i>=xidx[j]){
						ind = (i-xidx[j]) + ( *T * j );
						mpu[i] += beta[j]*x[ind];
					}
				} else{
					mpuinit += beta[j] * meanx[j];
				}
			}
		}
		// arp
		if(idx[2]>0){
			if(i>0){
				mpu[i] += alpha[0] * mpu[i-1];
			} else{
				mpuinit *= 1.0/(1.0 - alpha[0]);
				mpu[i] += alpha[0] * mpuinit;
			}
		}
		// common equation in gradient equation
		numer = k[0]*exp(-1.0*mpu[i])*(exp(k[0]*mpu[i])- 1.0*y[i]*pow(1.0+exp(mpu[i]), k[0]));
		denom = (1.0+exp(-1.0*mpu[i]))*(exp(mpu[i]*k[0]) - pow(1.0+exp(mpu[i]), k[0]));
		ktmp = (log(1.0+exp(-1.0*mpu[i]))*(exp(k[0]*mpu[i])- 1.0*y[i]*pow(1.0+exp(mpu[i]), k[0])))/(exp(mpu[i]*k[0]) - pow(1.0+exp(mpu[i]), k[0]));
		xtmp = numer/denom;
		if(idx[0]>0){
			if(i>0){
				dpomega[i] = 1.0 + alpha[0] * dpomega[i-1];
			} else{
				dpomega[i] = 1.0 + alpha[0]/(1.0-alpha[0]);
			}
			*domega += xtmp * dpomega[i];
			dvomega[i] = xtmp * dpomega[i];
		}
		if(idx[2]>0){
			if(i>0){
				dpalpha[i] = mpu[i-1] + alpha[0] * dpalpha[i-1];
				dvalpha[i] = xtmp * dpalpha[i];
				*dalpha += dvalpha[i];
			} else{
				dpalpha[i] = mpuinit + (alpha[0]*mpuinit)/(1-alpha[0]);
				dvalpha[i] = xtmp * dpalpha[i];
				*dalpha += dvalpha[i];
			}
		}
		if(idx[4]==0){
			if(idx[3]>0){
				for(j=0;j<idx[3];j++){
					if(i>0){
						if(i>j){
							dpdelta[i+(*T * j)] = y[i-(j+1)] + alpha[0] * dpdelta[i+(*T * j)-1];
							dvdelta[i+(*T * j)] = xtmp*dpdelta[i+(*T * j)];
							ddelta[j] += dvdelta[i+(*T * j)];
						}
					} else{
						dpdelta[i+(*T * j)] = (alpha[0] * meany[0])/(1 - alpha[0]);
						dvdelta[i+(*T * j)] = xtmp*dpdelta[i+(*T * j)];
						ddelta[j] += dvdelta[i+(*T * j)];
					}
				}
			}
		}
		if(idx[1]>0){
			for(j=0;j<idx[1];j++){
				if(i>0){
					ind = (i-xidx[j]) + ( *T * j );
					if(i>=xidx[j]){
						dpbeta[i+(*T * j)] = x[ind] + alpha[0] * dpbeta[i+(*T * j)-1];
						dvbeta[i+(*T * j)] =  xtmp * dpbeta[i+(*T * j)];
						dbeta[j] += dvbeta[i+(*T * j)] ;
					}
				} else{
					dpbeta[i+(*T * j)] = (alpha[0] * meanx[j])/(1 - alpha[0]);
					dvbeta[i+(*T * j)] =  xtmp * dpbeta[i+(*T * j)];
					dbeta[j] += dvbeta[i+(*T * j)] ;
				}

			}
		}
		dpk[i] = -1.0;
		dvk[i] = ktmp * dpk[i];
		*dk += dvk[i];
	}
}


void c_dbmfilter(double *y, double *x, double *mpu, double *omega,
		double *alpha, double *delta, double *beta, double *k,
		int *idx, int *xidx, double *lik, int *T)
{
	int i, j, ind;
	for(i=T[0];i<T[1];i++){
		// include intercept?
		if(idx[0]>0){
			mpu[i] += *omega;
		}
		// include arp?
		if(idx[2]>0){
			for(j=0;j<idx[2];j++){
				// make sure i>lag (constrained to lag=1 for now)
				if(j<i){
					mpu[i] += alpha[j]*mpu[i-(j+1)];
				}
			}
		}
		// include arq?
		// ecm model only allowed for lag=1
		if(idx[4]>0){
			if(0<i){
				mpu[i] += (1.0-alpha[0])*y[i-1];
			}
		} else{
			if(idx[3]>0){
				for(j=0;j<idx[3];j++){
					if(j<i){
						mpu[i] += delta[j]*y[i-(j+1)];
					}
				}
			}
		}
		// '<=' used since we have zero indexing
		if(idx[1]>0 && i>0){
			for(j=0;j<idx[1];j++){
				if(i>=xidx[j]){
					ind = (i-xidx[j]) + ( T[1] * j );
					mpu[i] += beta[j]*x[ind];
				}
			}
		}
		lik[i] = dlink(mpu[i], y[i], idx[5], k[0]);
	}
}
