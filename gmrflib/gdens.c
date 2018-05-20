
/* gdens.c
 * 
 * Copyright (C) 2001-2006 Havard Rue
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * The author's contact information:
 *
 *        Haavard Rue
 *        CEMSE Division
 *        King Abdullah University of Science and Technology
 *        Thuwal 23955-6900, Saudi Arabia
 *        Email: haavard.rue@kaust.edu.sa
 *        Office: +966 (0)12 808 0640
 *
 */

/*!
  \file gdens.c
  \brief Functions for constructing and working on log-linear and log-quadratic approximations to
  densities. [WARNING: OUTDATED...]
  
  \note This functionality is replace by \file density.c
  
  \note This code could benefite from rewriting, such that each element is switched to a linear
   instead of quadratic it probability for the cell is low. The code is not that efficient nor
   robust as we would like. Linear extrapolation in the tails should be optional.
*/

#include <stddef.h>
#include <float.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

#define TWODIVSQRTPI    1.1283791670955125738961589031215452   /* 2/sqrt(pi) */

#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

/* Pre-hg-Id: $Id: gdens.c,v 1.30 2009/08/26 06:19:51 hrue Exp $ */

/* 
   note: the representation is:  exp(c*x^2+b*x) 
 */

double GMRFLib_erfi(double x)
{
	return erfi_(&x);
}
double GMRFLib_lerf_diff(double x1, double x0)
{
	/*
	 * compute log(erf(x1)-erf(x0)) safely, well, at least when x0 is high and x1-x0 is small, or -x0 is high and x1-x0 is
	 * small.
	 * 
	 * we know that if x0 is high, then x1-x0 are close (hopefully), so if we're out of the limits for the 'erf' function,
	 * just before it breaks down numerically, we switch to an approximation, for high x0 and small (x1-x0).
	 * 
	 */
	double limit = 5.0, xmin, dx, value;

	if (x1 == x0)
		return -DBL_MAX;
	else {
		if (ABS(x1) < limit || ABS(x0) < limit)
			return log((double) (gsl_sf_erf(x1) - gsl_sf_erf(x0)));
		else {
			xmin = DMIN(x0, x1);
			dx = DMAX(x0, x1) - xmin;
			value = -SQR(xmin) + log((1.0 - exp(-2. * dx * xmin)) / (2. * xmin)) + log(TWODIVSQRTPI);
		}
	}
	return value;
}
double GMRFLib_lerfi_diff(double x1, double x0)
{
	/*
	 * same comment apply as for lerf_diff 
	 */
	double limit = 25.0, xmin, dx, value;

	if (x1 == x0) {
		return -DBL_MAX;
	} else {
		if (ABS(x1) < limit || ABS(x0) < limit)
			return log(GMRFLib_erfi(x1) - GMRFLib_erfi(x0));
		else {
			xmin = DMIN(x0, x1);
			dx = DMAX(x0, x1) - xmin;
			value = SQR(xmin) + log((exp(2. * dx * xmin) - 1.0) / (2. * xmin)) + log(TWODIVSQRTPI);
		}
	}
	return value;
}
int GMRFLib_gdens_ElmCompare(const void *e, const void *ee)
{
	const GMRFLib_gdens_Elm *elm, *eelm;

	elm = (const GMRFLib_gdens_Elm *) e;
	eelm = (const GMRFLib_gdens_Elm *) ee;

	if (elm->lp > eelm->lp)
		return -1;
	if (elm->lp < eelm->lp)
		return 1;

	return 0;
}

int GMRFLib_gdens_Free(GMRFLib_gdens_tp * ptr)
{
	if (ptr) {
		Free(ptr->elm);
		Free(ptr);
	}
	return 0;
}

int GMRFLib_gdens_1Update(GMRFLib_gdens_tp * ptr)
{
	/*
	 * update ptr 
	 */
	double lfmax, lpmax, lpsum, psum;
	int i, n;

	n = ptr->n;

	/*
	 * this is just for scaling to prevent overflow/underflow 
	 */
	for (lfmax = ptr->elm[0].lfl, i = 0; i < ptr->n; i++) {
		lfmax = DMAX(lfmax, ptr->elm[i].lfr);
		lfmax = DMAX(lfmax, ptr->elm[i].lfl);
	}

	/*
	 * compute various constants 
	 */
	for (i = 0; i < n; i++)
		ptr->elm[i].b = (ptr->elm[i].lfr - ptr->elm[i].lfl) / ptr->elm[i].dx;

	for (i = 0; i < n; i++) {
		if (ISZERO(ptr->elm[i].b))
			ptr->elm[i].lp = ptr->elm[i].lfl - lfmax + log(ptr->elm[i].dx);
		else {
			double arg = (exp(ptr->elm[i].dx * ptr->elm[i].b) - 1.0) / ptr->elm[i].b;

			if (ISZERO(arg) || arg < 0.)
				ptr->elm[i].lp = ptr->elm[i].lfl - lfmax + log(ptr->elm[i].dx);
			else
				ptr->elm[i].lp = ptr->elm[i].lfl - lfmax + log(arg);
		}
	}

	for (lpmax = ptr->elm[0].lp, i = 1; i < n; i++)
		lpmax = DMAX(lpmax, ptr->elm[i].lp);
	for (i = 0; i < n; i++)
		ptr->elm[i].lp -= lpmax;
	for (i = 0, psum = 0.0; i < n; i++)
		psum += exp(ptr->elm[i].lp);
	for (lpsum = log(psum), i = 0; i < n; i++)
		ptr->elm[i].lp -= lpsum;

	/*
	 * speedup the search for regions and sampling 
	 */
	qsort(ptr->elm, (size_t) ptr->n, sizeof(GMRFLib_gdens_Elm), GMRFLib_gdens_ElmCompare);

	return 0;
}

int GMRFLib_gdens_2Update(GMRFLib_gdens_tp * ptr)
{
	/*
	 * update ptr 
	 */
	double lpmax, lpsum, psum;
	int i, n = ptr->n;

	/*
	 * compute various constants 
	 */
	for (i = 0; i < n; i++) {
		double s, x0 = ptr->elm[i].xl, x1 = ptr->elm[i].xm, x2 = ptr->elm[i].xr, f0 = ptr->elm[i].lfl, f1 = ptr->elm[i].lfm, f2 =
		    ptr->elm[i].lfr;

		ptr->elm[i].lfmax = f0;			       /* scale and shift */
		f1 -= f0;
		f2 -= f0;
		f0 = 0.0;
		x1 -= x0;
		x2 -= x0;
		x0 = 0.0;

		s = 1. / ((-SQR(x2) + x2 * x1) * x1);
		ptr->elm[i].b = -s * (-f2 * SQR(x1) + f1 * SQR(x2));
		ptr->elm[i].c = s * (-f2 * x1 + f1 * x2);
	}

	/*
	 * compute the integral over each cell, unnormalized 
	 */
	for (i = 0; i < n; i++) {
		if (ISZERO(ptr->elm[i].c)) {
			if (ISZERO(ptr->elm[i].b))
				ptr->elm[i].lnormc = log(ptr->elm[i].xr - ptr->elm[i].xl);
			else {
				double tmp = (ptr->elm[i].xr - ptr->elm[i].xl) * ptr->elm[i].b;

				if (fabs(tmp) > FLT_EPSILON) {
					/*
					 * all ok 
					 */
					ptr->elm[i].lnormc = log(fabs(exp((ptr->elm[i].xr - ptr->elm[i].xl) * ptr->elm[i].b) - 1.0))
					    - log(fabs(ptr->elm[i].b));
				} else {
					/*
					 * this is the expansion for small 'tmp' 
					 */
					ptr->elm[i].lnormc = log(fabs(tmp)) - log(fabs(ptr->elm[i].b));
				}
			}
		} else {
			/*
			 * this can (and does) go wrong as the 'erf' and 'erfi' function may fail for high values of the
			 * arguments. perhaps we need to intergrate directly int(exp(x^2),x=low...high) ? 
			 */

			double aa = 0.5 * ptr->elm[i].b / sqrt(fabs(ptr->elm[i].c));

			if (ptr->elm[i].c <= 0.0)
				ptr->elm[i].lnormc = -0.25 * SQR(ptr->elm[i].b) / ptr->elm[i].c + log(1. / (TWODIVSQRTPI * sqrt(-ptr->elm[i].c)))
				    + GMRFLib_lerf_diff(sqrt(-ptr->elm[i].c) * (ptr->elm[i].xr - ptr->elm[i].xl) - aa, -aa);
			else
				ptr->elm[i].lnormc = -0.25 * SQR(ptr->elm[i].b) / ptr->elm[i].c + log(1. / (TWODIVSQRTPI * sqrt(ptr->elm[i].c)))
				    + GMRFLib_lerfi_diff(sqrt(ptr->elm[i].c) * (ptr->elm[i].xr - ptr->elm[i].xl) + aa, aa);
		}
		if (ISINF(ptr->elm[i].lnormc)) {
			ptr->elm[i].lnormc = ptr->elm[i].b = ptr->elm[i].c = 0;
			ptr->elm[i].lp = ptr->elm[i].lfl = ptr->elm[i].lfm = ptr->elm[i].lfr = ptr->elm[i].lfmax = -FLT_MAX;
		} else {
			ptr->elm[i].lp = ptr->elm[i].lnormc + ptr->elm[i].lfmax;
		}
	}

	for (lpmax = ptr->elm[0].lp, i = 1; i < n; i++)
		lpmax = DMAX(lpmax, ptr->elm[i].lp);
	for (i = 0; i < n; i++)
		ptr->elm[i].lp -= lpmax;
	for (i = 0, psum = 0.0; i < n; i++)
		psum += exp(ptr->elm[i].lp);
	for (lpsum = log(psum), i = 0; i < n; i++)
		ptr->elm[i].lp -= lpsum;

	/*
	 * speedup the search for regions and sampling 
	 */
	qsort(ptr->elm, (size_t) ptr->n, sizeof(GMRFLib_gdens_Elm), GMRFLib_gdens_ElmCompare);

	if (0) {
		for (i = 0; i < n; i++) {
			printf("Gdens2 region %d\n", i);
			printf("\txl %f xr %f\n", ptr->elm[i].xl, ptr->elm[i].xr);
			printf("\tlnormc %f lp %f\n", ptr->elm[i].lnormc, ptr->elm[i].lp);
		}
	}

	return 0;
}

GMRFLib_gdens_tp *GMRFLib_gdens_1Init(double xmin, double xmax, int n, int n_min, GMRFLib_gdens_Func_tp * lfunc, void *lfunc_arg)
{
	/*
	 * arguments are:
	 * 
	 * xmin : left limit xmax : right limit n : total number of interior points, which makes n+1 cells n_min: make first
	 * n_min interior points by equal division of (xmax-xmin), then n-n_min divisions by dividing the cell with the largest 
	 * probability in half. lfunc: the function returning log(f(x)), unormalized
	 * 
	 */

	GMRFLib_gdens_tp *ptr;
	int i, nn;
	double xnew, f, df, dx;

	ptr = Calloc(1, GMRFLib_gdens_tp);
	ptr->order = 1;
	ptr->xmin = DMIN(xmin, xmax);
	ptr->xmax = DMAX(xmin, xmax);
	ptr->n = IMAX(2, n_min);
	nn = IMAX(ptr->n, n);
	ptr->elm = Calloc(nn + 1, GMRFLib_gdens_Elm);

	/*
	 * make first equal spaced points 
	 */
	dx = (ptr->xmax - ptr->xmin) / ((double) n_min);
	for (i = 0; i < ptr->n; i++) {
		ptr->elm[i].x = ptr->xmin + i * dx;
		ptr->elm[i].dx = dx;
		ptr->elm[i].lfl = lfunc(ptr->elm[i].x, lfunc_arg);
	}
	for (i = 0; i < ptr->n - 1; i++)
		ptr->elm[i].lfr = ptr->elm[i + 1].lfl;
	ptr->elm[ptr->n - 1].lfr = lfunc(ptr->elm[ptr->n - 1].x + ptr->elm[ptr->n - 1].dx, lfunc_arg);

	GMRFLib_gdens_1Update(ptr);

	/*
	 * decide now the rest of the points by dividing the cell with largest probability in half 
	 */
	for (i = n_min; i < n; i++) {
		/*
		 * split the first elm, and add one at the end 
		 */
		if (ABS(ptr->elm[0].b) < FLT_EPSILON)
			f = 0.5;
		else {
			df = ptr->elm[0].b * ptr->elm[0].dx;
			f = log(0.5 * (exp(df) + 1.0)) / df;
		}

		xnew = ptr->elm[0].x + ptr->elm[0].dx * f;
		ptr->elm[ptr->n].lfr = ptr->elm[0].lfr;
		ptr->elm[ptr->n].x = xnew;
		ptr->elm[ptr->n].dx = ptr->elm[0].dx * (1. - f);
		ptr->elm[0].lfr = ptr->elm[ptr->n].lfl = lfunc(xnew, lfunc_arg);
		ptr->elm[0].dx *= f;

		ptr->n++;
		GMRFLib_gdens_1Update(ptr);
	}

	return ptr;
}

GMRFLib_gdens_tp *GMRFLib_gdens_2Init(double xmin, double xmax, int n, int n_min, GMRFLib_gdens_Func_tp * lfunc, void *lfunc_arg)
{
	/*
	 * arguments are:
	 * 
	 * xmin : left limit xmax : right limit n : total number of interior points, which makes n+1 cells n_min: make first
	 * n_min interior points by equal division of (xmax-xmin), then n-n_min divisions by dividing the cell with the largest 
	 * probability in half. lfunc: the function returning log(f(x,lfunc_arg)), unormalized 
	 */

	GMRFLib_gdens_tp *ptr;
	int i, nn;
	double dx, df, f;

	ptr = Calloc(1, GMRFLib_gdens_tp);
	ptr->order = 2;
	ptr->xmin = DMIN(xmin, xmax);
	ptr->xmax = DMAX(xmin, xmax);
	ptr->n = IMAX(1, n_min);
	nn = IMAX(ptr->n, n);
	ptr->elm = Calloc(nn + 1, GMRFLib_gdens_Elm);

	/*
	 * make first equal spaced points 
	 */
	dx = (ptr->xmax - ptr->xmin) / ((double) n_min);
	for (i = 0; i < ptr->n; i++) {
		ptr->elm[i].xl = ptr->xmin + i * dx;
		ptr->elm[i].xr = ptr->elm[i].xl + dx;
		ptr->elm[i].lfl = (i == 0 ? lfunc(ptr->elm[i].xl, lfunc_arg) : ptr->elm[i - 1].lfr);
		ptr->elm[i].lfr = lfunc(ptr->elm[i].xr, lfunc_arg);

		/*
		 * now, try to decide a nice place to place the midpoint. currently, i try to split the region into
		 * (approximately) equal integrals 
		 */
		df = ptr->elm[i].lfr - ptr->elm[i].lfl;
		f = (ISZERO(df) ? 0.5 : log(0.5 * (exp(df) + 1.0)) / df);

		/*
		 * need also a more safe option than ISZERO 
		 */
		if (ISZERO(f) || ISZERO(1. - f))
			f = 0.5;
		if (f < FLT_EPSILON || 1. - f < FLT_EPSILON)
			f = 0.5;

		if (f == 0.5) {
			/*
			 * this will emulate a log-linear spline 
			 */
			ptr->elm[i].xm = ptr->elm[i].xl + dx * f;
			ptr->elm[i].lfm = 0.5 * (ptr->elm[i].lfl + ptr->elm[i].lfr);
		} else {
			ptr->elm[i].xm = ptr->elm[i].xl + dx * f;
			ptr->elm[i].lfm = lfunc(ptr->elm[i].xm, lfunc_arg);
		}
	}
	for (i = 0; i < ptr->n - 1; i++)
		ptr->elm[i].lfr = ptr->elm[i + 1].lfl;
	ptr->elm[ptr->n - 1].lfr = lfunc(ptr->elm[ptr->n - 1].xr, lfunc_arg);

	GMRFLib_gdens_2Update(ptr);

	/*
	 * decide now the rest of the points by dividing the cell with largest probability in half. try to split the region
	 * into (approximately) equal integrals 
	 */
	for (i = n_min; i < nn; i++) {
		/*
		 * split the first elm, and add one at the end 
		 */
		ptr->elm[ptr->n].lfl = ptr->elm[0].lfm;
		ptr->elm[ptr->n].xl = ptr->elm[0].xm;
		ptr->elm[ptr->n].lfr = ptr->elm[0].lfr;
		ptr->elm[ptr->n].xr = ptr->elm[0].xr;
		df = ptr->elm[ptr->n].lfr - ptr->elm[ptr->n].lfl;
		f = (ISZERO(df) ? 0.5 : log(0.5 * (exp(df) + 1.0)) / df);
		if (ISZERO(f) || ISZERO(1. - f))
			f = 0.5;
		if (f < FLT_EPSILON || 1. - f < FLT_EPSILON)
			f = 0.5;
		if (f == 0.5) {
			/*
			 * this will emulate a log-linear spline 
			 */
			ptr->elm[ptr->n].xm = ptr->elm[ptr->n].xl + f * (ptr->elm[ptr->n].xr - ptr->elm[ptr->n].xl);
			ptr->elm[ptr->n].lfm = 0.5 * (ptr->elm[ptr->n].lfl + ptr->elm[ptr->n].lfr);
		} else {
			ptr->elm[ptr->n].xm = ptr->elm[ptr->n].xl + f * (ptr->elm[ptr->n].xr - ptr->elm[ptr->n].xl);
			ptr->elm[ptr->n].lfm = lfunc(ptr->elm[ptr->n].xm, lfunc_arg);
		}

		ptr->elm[0].lfr = ptr->elm[0].lfm;
		ptr->elm[0].xr = ptr->elm[0].xm;
		df = ptr->elm[0].lfr - ptr->elm[0].lfl;
		f = (ISZERO(df) ? 0.5 : log(0.5 * (exp(df) + 1.0)) / df);
		if (ISZERO(f) || ISZERO(1. - f))
			f = 0.5;
		if (f < FLT_EPSILON || 1. - f < FLT_EPSILON)
			f = 0.5;
		if (f == 0.5) {
			/*
			 * this will emulate a log-linear spline 
			 */
			ptr->elm[0].xm = ptr->elm[0].xl + f * (ptr->elm[0].xr - ptr->elm[0].xl);
			ptr->elm[0].lfm = 0.5 * (ptr->elm[0].lfl + ptr->elm[0].lfr);
		} else {
			ptr->elm[0].xm = ptr->elm[0].xl + f * (ptr->elm[0].xr - ptr->elm[0].xl);
			ptr->elm[0].lfm = lfunc(ptr->elm[0].xm, lfunc_arg);
		}

		ptr->n++;
		GMRFLib_gdens_2Update(ptr);
	}

	return ptr;
}
double *GMRFLib_cut_points(double fac, int n)
{
	/*
	 * return the n+1 cutpoints between -fac...fac, that makes n regions of equal probability 
	 */

	static double fac_save = 0.0;

#pragma omp threadprivate(fac_save)
	static int n_save = 0;

#pragma omp threadprivate (n_save)
	static double *cp = NULL;

#pragma omp threadprivate(cp)

	double p_left, p, q, mean = 0.0, sd = 1.0;
	int i;

	if (fac == fac_save && n == n_save)
		return cp;

	fac_save = fac;
	n_save = n;

	Free(cp);
	cp = Calloc(n + 1, double);

	p = gsl_cdf_ugaussian_P((fac - mean) / sd);
	q = 1.0 - p;
	p_left = q;

	cp[0] = -fac;
	cp[n] = fac;
	for (i = 1; i < n; i++) {
		p = p_left + i / ((double) n) * (1.0 - 2.0 * p_left);
		cp[i] = mean + sd * gsl_cdf_ugaussian_Pinv(p);
	}

	if (0) {
		for (i = 0; i < n + 1; i++)
			printf("cutpoint %d %f\n", i, cp[i]);
	}
	return cp;
}

GMRFLib_gdens_tp *GMRFLib_gdens_2InitNew(double fac, double mean, double stdev, int n, GMRFLib_gdens_Func_tp * lfunc, void *lfunc_arg)
{
	/*
	 * arguments are:
	 * 
	 * limits are mean +- fac*stdev, and n-regions.
	 * 
	 * lfunc: the function returning log(f(x,lfunc_arg)), unormalized 
	 */

	GMRFLib_gdens_tp *ptr;
	int i;
	double dx, df, f, *cpoints;

	fac = ABS(fac);

	ptr = Calloc(1, GMRFLib_gdens_tp);
	ptr->xmin = mean - fac * stdev;
	ptr->xmax = mean + fac * stdev;
	ptr->n = n;
	ptr->elm = Calloc(n + 1, GMRFLib_gdens_Elm);

	cpoints = GMRFLib_cut_points(fac, n);

	for (i = 0; i < ptr->n; i++) {
		dx = stdev * (cpoints[i + 1] - cpoints[i]);
		ptr->elm[i].xl = mean + cpoints[i] * stdev;
		ptr->elm[i].xr = ptr->elm[i].xl + dx;
		ptr->elm[i].lfl = (i == 0 ? lfunc(ptr->elm[i].xl, lfunc_arg) : ptr->elm[i - 1].lfr);
		ptr->elm[i].lfr = lfunc(ptr->elm[i].xr, lfunc_arg);

		/*
		 * now, try to decide a nice place to place the midpoint. currently, i try to split the region into
		 * (approximately) equal integrals 
		 */
		df = ptr->elm[i].lfr - ptr->elm[i].lfl;
		f = (ISZERO(df) ? 0.5 : log(0.5 * (exp(df) + 1.0)) / df);

		/*
		 * need also a more safe option than ISZERO 
		 */
		if (ISZERO(f) || ISZERO(1. - f))
			f = 0.5;
		if (f < FLT_EPSILON || 1. - f < FLT_EPSILON)
			f = 0.5;

		if (f == 0.5) {
			/*
			 * this will emulate a log-linear spline 
			 */
			ptr->elm[i].xm = ptr->elm[i].xl + dx * f;
			ptr->elm[i].lfm = 0.5 * (ptr->elm[i].lfl + ptr->elm[i].lfr);
		} else {
			ptr->elm[i].xm = ptr->elm[i].xl + dx * f;
			ptr->elm[i].lfm = lfunc(ptr->elm[i].xm, lfunc_arg);
		}
	}
	for (i = 0; i < ptr->n - 1; i++)
		ptr->elm[i].lfr = ptr->elm[i + 1].lfl;
	ptr->elm[ptr->n - 1].lfr = lfunc(ptr->elm[ptr->n - 1].xr, lfunc_arg);

	GMRFLib_gdens_2Update(ptr);

	return ptr;
}

GMRFLib_gdens_tp *GMRFLib_gdens_InitNew(double fac, double mean, double stdev, int n, GMRFLib_gdens_Func_tp * lfunc, void *lfunc_arg)
{
	return GMRFLib_gdens_2InitNew(fac, mean, stdev, n, lfunc, lfunc_arg);
}

GMRFLib_gdens_tp *GMRFLib_gdens_Init(double xmin, double xmax, int n, int n_min, int order, GMRFLib_gdens_Func_tp * lfunc, void *lfunc_arg)
{
	if (!(order == 1 || order == 2)) {
		GMRFLib_ERROR_NO_RETURN(GMRFLib_EPARAMETER);
		return NULL;
	}

	if (order == 1)
		return GMRFLib_gdens_1Init(xmin, xmax, n, n_min, lfunc, lfunc_arg);
	else
		return GMRFLib_gdens_2Init(xmin, xmax, n, n_min, lfunc, lfunc_arg);
}
int GMRFLib_gdens_1Eval(GMRFLib_gdens_tp * ptr, double *x, double *ldens, int flag)
{
	int i, r;
	double psum, p, xx;

	if (flag) {					       /* sample, draw first segment then within the segment */
		p = GMRFLib_uniform();
		r = ptr->n - 1;
		for (i = 0, psum = 0.0; i < ptr->n - 1; i++)
			if ((psum += exp(ptr->elm[i].lp)) > p) {
				r = i;
				break;
			}
		p = GMRFLib_uniform();
		if (ISZERO(ptr->elm[r].b))
			*x = ptr->elm[r].x + p * ptr->elm[r].dx;
		else
			*x = ptr->elm[r].x + log(1.0 - p + p * exp(ptr->elm[r].b * ptr->elm[r].dx)) / ptr->elm[r].b;
	} else {					       /* evaluate the density, find the corresponding segment */

		r = ptr->n - 1;
		for (i = 0; i < ptr->n - 1; i++)
			if (*x >= ptr->elm[i].x && *x < ptr->elm[i].x + ptr->elm[i].dx) {
				r = i;
				break;
			}
	}

	if (*x < ptr->xmin || *x > ptr->xmax) {
		fprintf(stdout, "\nGMRFLib_gdens_1: x=%f is not in [%f,%f], truncate!\n", *x, ptr->xmin, ptr->xmax);
		xx = DMAX(ptr->xmin, DMIN(ptr->xmax, *x));
	} else
		xx = *x;

	if (ISZERO(ptr->elm[r].b))
		*ldens = ptr->elm[r].lp - log(ptr->elm[r].dx);
	else {
		double arg = ptr->elm[r].b * exp(ptr->elm[r].b * (xx - ptr->elm[r].x))
		    / (exp(ptr->elm[r].b * ptr->elm[r].dx) - 1.0);

		if (ISZERO(arg) || arg <= 0.0)
			*ldens = ptr->elm[r].lp - log(ptr->elm[r].dx);
		else
			*ldens = ptr->elm[r].lp + log(arg);
	}

	return 0;
}
int GMRFLib_gdens_2Eval(GMRFLib_gdens_tp * ptr, double *x, double *ldens, int flag)
{
	int i, r;
	double psum, p, xx, xlocalmax, lfmax, b, dx;

	if (flag) {					       /* sample, draw first segment then within the segment */
		p = GMRFLib_uniform();
		r = ptr->n - 1;
		for (i = 0, psum = 0.0; i < ptr->n - 1; i++)
			if ((psum += exp(ptr->elm[i].lp)) > p) {
				r = i;
				break;
			}
		/*
		 * here i use a rejection sampling algorithm with two proposals: a uniform and a linear in log-scale. pick your 
		 * choice ;-) 
		 */
		if (1) {
			/*
			 * uniform 
			 */
			dx = ptr->elm[r].xr - ptr->elm[r].xl;
			if (ISZERO(ptr->elm[r].c))
				lfmax = DMAX(ptr->elm[r].lfr, ptr->elm[r].lfl) - ptr->elm[r].lfmax;
			else {
				xlocalmax = -ptr->elm[r].b / (2.0 * ptr->elm[r].c);
				xlocalmax = DMAX(0.0, xlocalmax);
				xlocalmax = DMIN(dx, xlocalmax);
				lfmax = ptr->elm[r].c * SQR(xlocalmax) + ptr->elm[r].b * xlocalmax;
			}
			do
				xx = GMRFLib_uniform() * dx;
			while (GMRFLib_uniform() > exp(ptr->elm[r].c * SQR(xx) + ptr->elm[r].b * xx - lfmax));
			*x = xx + ptr->elm[r].xl;
		} else {
			/*
			 * linear 
			 */
			dx = ptr->elm[r].xr - ptr->elm[r].xl;
			b = (ptr->elm[r].lfr - ptr->elm[r].lfl) / dx;

			if (ISZERO(ptr->elm[r].c)) {
				p = GMRFLib_uniform();
				xx = log(1.0 - p + p * exp(b * dx)) / b;	/* are on a straight line */
			} else {
				xlocalmax = -(ptr->elm[r].b - b) / (2.0 * ptr->elm[r].c);
				if (ptr->elm[r].c < 0.0) {
					xlocalmax = DMAX(0.0, xlocalmax);
					xlocalmax = DMIN(dx, xlocalmax);
					lfmax = ptr->elm[r].c * SQR(xlocalmax) + (ptr->elm[r].b - b) * xlocalmax;
				} else {
					double lf_left, lf_right;

					lf_left = 0.0;
					lf_right = ptr->elm[r].c * SQR(dx) + (ptr->elm[r].b - b) * dx;
					lfmax = (lf_right > lf_left ? lf_right : lf_left);
				}

				do {
					p = GMRFLib_uniform();
					if (ISZERO(b))
						xx = dx * p;
					else
						xx = log(1.0 - p + p * exp(b * dx)) / b;
				}
				while (GMRFLib_uniform() > exp(ptr->elm[r].c * SQR(xx) + (ptr->elm[r].b - b) * xx - lfmax));
			}
			*x = xx + ptr->elm[r].xl;
		}
	} else {					       /* evaluate the density, find the corresponding segment */

		r = ptr->n - 1;
		for (i = 0; i < ptr->n - 1; i++)
			if (*x >= ptr->elm[i].xl && *x < ptr->elm[i].xr) {
				r = i;
				break;
			}
	}

	if (*x < ptr->xmin || *x > ptr->xmax) {
		if (0)
			fprintf(stdout, "\nGMRFLib_gdens_2: x=%f is not in [%f,%f], truncate!\n", *x, ptr->xmin, ptr->xmax);
		xx = DMAX(ptr->xmin, DMIN(ptr->xmax, *x));     /* only happens in eval-mode */
	} else {
		xx = *x;
	}

	xx -= ptr->elm[r].xl;				       /* relative value within the cell */

	*ldens = ptr->elm[r].lp + ptr->elm[r].c * SQR(xx) + ptr->elm[r].b * xx - ptr->elm[r].lnormc;

	return 0;
}
int GMRFLib_gdens_Eval(GMRFLib_gdens_tp * ptr, double *x, double *ldens, int flag)
{
	if (ptr->order == 1)
		return GMRFLib_gdens_1Eval(ptr, x, ldens, flag);
	else
		return GMRFLib_gdens_2Eval(ptr, x, ldens, flag);
}
int GMRFLib_gdens_1Sample(GMRFLib_gdens_tp * ptr, double *x, double *ldens)
{
	return GMRFLib_gdens_1Eval(ptr, x, ldens, 1);
}
int GMRFLib_gdens_2Sample(GMRFLib_gdens_tp * ptr, double *x, double *ldens)
{
	return GMRFLib_gdens_2Eval(ptr, x, ldens, 1);
}
int GMRFLib_gdens_Sample(GMRFLib_gdens_tp * ptr, double *x, double *ldens)
{
	if (ptr->order == 1)
		return GMRFLib_gdens_1Eval(ptr, x, ldens, 1);
	else
		return GMRFLib_gdens_2Eval(ptr, x, ldens, 1);
}
int GMRFLib_gdens_1LDens(GMRFLib_gdens_tp * ptr, double *x, double *ldens)
{
	return GMRFLib_gdens_1Eval(ptr, x, ldens, 0);
}
int GMRFLib_gdens_2LDens(GMRFLib_gdens_tp * ptr, double *x, double *ldens)
{
	return GMRFLib_gdens_2Eval(ptr, x, ldens, 0);
}
int GMRFLib_gdens_LDens(GMRFLib_gdens_tp * ptr, double *x, double *ldens)
{
	if (ptr->order == 1)
		return GMRFLib_gdens_1Eval(ptr, x, ldens, 0);
	else
		return GMRFLib_gdens_2Eval(ptr, x, ldens, 0);
}
