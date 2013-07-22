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
#ifndef __FILTERS_H__
#define __FILTERS_H__
void sgarchfilteracd(int *model, double *pars, int *idx, double *vexdata, double *e, int T, int i, double *h);
void csgarchfilteracd(int *model, double *pars, int *idx, double *e, double *vexdata, int T, int i, double *h, double *q);
void arfimaxfilteracd(int* model, double *pars, int *idx, double *x, double *res, double *mexdata, double *zrf,
		double *constm, double *condm, double h, double sk, double sh, int m, int i, int T);
void armaxsimacd(int *model, double *pars, int *idx, double *x, double *res, double *constm, int *m, int *T);
#endif
