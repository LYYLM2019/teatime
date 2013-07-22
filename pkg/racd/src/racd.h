/*################################################################################
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
#################################################################################*/
#ifndef RACD_H
#define RACD_H
void arfimaxacdfilterC(int *model, double *pars, int *idx, double *x, double *res, double *mexdata,
		double *zrf, double *constm, double *condm, double *h, double *tskew, double *tshape,
		int *m, int *T);
void sacdfilterC(int *model, double *pars, int *idx, double *hEst, double *x,
		double *res, double *e, double *mexdata, double *vexdata, double *zrf,
		double *constm, double *condm, int *m, int *T, double *h, double *z,
		double *tempskew, double *tempshape, double *skhEst, double *tskew,
		double *tshape, double *skxreg, double *shxreg, double *sbounds,
		double *llh, double *LHT);
void mcsacdfilterC(int *model, double *pars, int *idx, double *hEst, double *res, double *e,
		double *eres, double *s, double *v, double *vexdata, int *m, int *T, double *h, double *z,
		double *tempskew, double *tempshape, double *skhEst, double *tskew, double *tshape,
		double *skxreg, double *shxreg, double *sbounds, double *llh, double *LHT);
void sacdsimC(int *model, double *pars, int *idx, double *h, double *z,
		double *res, double *e, double *tempskew, double *tempshape,
		double *tskew, double *tshape, double *sbounds, double *vexdata,
		double *skxreg, double *shxreg, int *T, int *m);
void csacdfilterC(int *model, double *pars, int *idx, double *hEst, double *x,
		double *res, double *e, double *mexdata, double *vexdata, double *zrf,
		double *constm, double *condm, int *m, int *T, double *h, double *q,
		double *z, double *tempskew, double *tempshape, double *skhEst,
		double *tskew, double *tshape, double *skxreg, double *shxreg,
		double *sbounds, double *llh, double *LHT);
void csacdsimC(int *model, double *pars, int *idx, double *h, double *q, double *z,
		double *res, double *e, double *tempskew, double *tempshape,
		double *tskew, double *tshape, double *sbounds, double *vexdata,
		double *skxreg, double *shxreg, int *T, int *m);
void mcsacdsimC(int *model, double *pars, int *idx, double *h, double *s, double *v,
		double *z, double *res, double *eres, double *e,
		double *tempskew, double *tempshape, double *tskew, double *tshape,
		double *sbounds, double *vexdata, double *skxreg, double *shxreg,
		int *T, int *m);
#endif /* RACD_H */
