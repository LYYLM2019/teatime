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

#include "dbmpredict.h"
using namespace Rcpp;

SEXP dbmpk1(SEXP p, SEXP omega, SEXP alpha, SEXP delta, SEXP beta, SEXP kappa, SEXP y, SEXP x, SEXP idx, SEXP xidx, SEXP V)
{
	try{
		Rcpp::NumericVector xp(p);
		Rcpp::NumericVector xy(y);
		Rcpp::NumericMatrix xx(x);
		double *xomega = REAL(omega);
		Rcpp::NumericVector xalpha(alpha);
		Rcpp::NumericVector xdelta(delta);
		Rcpp::NumericVector xbeta(beta);
		double *xk = REAL(kappa);
		Rcpp::IntegerVector zid(idx);
		Rcpp::IntegerVector zidx(xidx);
		Rcpp::IntegerVector ix(V);
		int k = ix[0];
		//because of zero indexing we subtract 1
		int T = ix[1];
		int link = ix[2];
		int j, i;
		Rcpp::NumericVector zk(k, 0.0);
		Rcpp::NumericVector zy(k, 0.0);
		double tmp2 = 0.0;
		double ans = 0.0;
		if(zid[1]>0){
			for(j=0;j<k;j++){
				tmp2 = 0.0;
				for(i=0;i<zid[1];i++){
					if(zidx[i]>(j+1)){
						tmp2+=xbeta[i]*xx(T-zidx[i]+j,i);
					}
				}
				zk[j] = tmp2;
			}
		}
		if(zid[4]>0){
			for(j=0;j<k;j++){
				zy[j] = (1.0-xalpha[0])*xy[T-1+j];
			}
		} else{
			if(zid[3]>0){
				for(j=0;j<k;j++){
					tmp2 = 0.0;
					for(i=0;i<zid[3];i++){
						tmp2 += xdelta[i]*xy[T+j-i-1];
					}
					zy[j] = tmp2;
				}
			}
		}
		tmp2 = 0.0;
		for(j=0;j<k;j++){
			tmp2+=pow(xalpha[0],j*1.0)*(*xomega+zy[j]+zk[j]);
		}
		tmp2+=pow(xalpha[0],k*1.0)*xp[T-1];
		if(link==1){
			ans = pnorm(NumericVector(1,tmp2), 0.0, 1.0)[0];
		} else if(link==2){
			ans = plogis(NumericVector(1,tmp2), 0.0, 1.0)[0];
		} else{
			ans = pow(1.0/(1+exp(-1.0*tmp2)), *xk);
		}
		return wrap(ans);
		} catch( std::exception &ex ) {
			forward_exception_to_r( ex );
		} catch(...) {
					::Rf_error( "dbm-->dbmpredict c++ exception (unknown reason)" );
				}
		return R_NilValue;
}

SEXP dbmpk2(SEXP p, SEXP omega, SEXP alpha, SEXP delta, SEXP beta, SEXP kappa, SEXP y, SEXP x, SEXP idx, SEXP xidx, SEXP V)
{
	try{
		Rcpp::NumericVector xp(p);
		Rcpp::NumericVector xy(y);
		Rcpp::NumericMatrix xx(x);
		double *xomega = REAL(omega);
		Rcpp::NumericVector xalpha(alpha);
		Rcpp::NumericVector xdelta(delta);
		Rcpp::NumericVector xbeta(beta);
		double *xk = REAL(kappa);
		Rcpp::IntegerVector zid(idx);
		Rcpp::IntegerVector zidx(xidx);
		Rcpp::IntegerVector ix(V);
		int k = ix[0];
		int T = ix[1];
		int link = ix[2];
		int j, i;
		Rcpp::NumericVector zk(k, 0.0);
		Rcpp::NumericVector zy(k, 0.0);
		double tmp2 = 0.0;
		double ans = 0.0;
		if(zid[1]>0){
			for(j=0;j<k;j++){
				tmp2 = 0.0;
				for(i=0;i<zid[1];i++){
					if((k-zidx[i])<(j+1)){
						tmp2+=xbeta[i]*xx(T+k-zidx[i]-j-1,i);
					}
				}
				zk[j] = tmp2;
			}
		}
		if(zid[4]>0){
			for(j=0;j<k;j++){
				zy[j] = (1.0-xalpha[0])*xy[T+k-j-2];
			}
		} else{
			if(zid[3]>0){
				for(j=0;j<k;j++){
					tmp2 = 0.0;
					for(i=0;i<zid[3];i++){
						if(i>=j){
							tmp2 += xdelta[i]*xy[T+k-j-i-2];
						}
					}
					zy[j] = tmp2;
				}
			}
		}
		tmp2 = 0.0;
		for(j=0;j<k;j++){
			tmp2+=pow(xalpha[0],j*1.0)*(*xomega+zy[j]+zk[j]);
		}
		tmp2+=pow(xalpha[0],k*1.0)*xp[T-1];

		if(link==1){
			ans = pnorm(NumericVector(1,tmp2), 0.0, 1.0)[0];
		} else if(link==2){
			ans = plogis(NumericVector(1,tmp2), 0.0, 1.0)[0];
		} else{
			ans = pow(1.0/(1+exp(-1.0*tmp2)), *xk);
		}
		return wrap(ans);
		} catch( std::exception &ex ) {
			forward_exception_to_r( ex );
		} catch(...) {
			::Rf_error( "dbm-->dbmpredict c++ exception (unknown reason)" );
		}
		return R_NilValue;
}

