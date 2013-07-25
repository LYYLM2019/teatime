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
#ifndef DBM_H
#define DBM_H
double pburr2(const double , const double );
double dlink(const double , const double , const int, const double );
void c_dbmestimate(double *y, double *x, double *mpu, double *mpuinit, double *omega,
		double *alpha, double *delta, double *beta, double *k,
		double *lik, double *llh, int *idx, int *xidx, int *T);
void c_dbmfilter(double *y, double *x, double *mpu, double *omega,
		double *alpha, double *delta, double *beta, double *k,
		int *idx, int *xidx, double *lik, int *T);
void c_dbmderiv1(double *y, double *x, double *mpu,
		double *meanx, double *meany, double *omega,
		double *alpha, double *delta, double *beta, double *domega,
		double *dalpha, double *ddelta, double *dbeta,
		double *dpomega, double *dpalpha, double *dpdelta, double *dpbeta,
		double *dvomega, double *dvalpha, double *dvdelta, double *dvbeta,
		int *idx, int *xidx, int *T);
void c_dbmderiv2(double *y, double *x, double *mpu,
		double *meanx, double *meany, double *omega,
		double *alpha, double *delta, double *beta, double *domega,
		double *dalpha, double *ddelta, double *dbeta,
		double *dpomega, double *dpalpha, double *dpdelta, double *dpbeta,
		double *dvomega, double *dvalpha, double *dvdelta, double *dvbeta,
		int *idx, int *xidx, int *T);
void c_dbmderiv3(double *y, double *x, double *mpu, double *meanx, double *meany,
		double *omega, double *alpha, double *delta, double *beta, double *k,
		double *domega, double *dalpha, double *ddelta, double *dbeta, double *dk,
		double *dpomega, double *dpalpha, double *dpdelta, double *dpbeta, double *dpk,
		double *dvomega, double *dvalpha, double *dvdelta, double *dvbeta, double *dvk,
		int *idx, int *xidx, int *T);
#endif /* DBM_H */
