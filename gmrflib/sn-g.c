
/* sn-g.c
 * 
 * Copyright (C) 2024-2024 Havard Rue
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

#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"

#define ORDER 5
static GMRFLib_spline_tp **cxs = NULL;
static GMRFLib_spline_tp **icxs = NULL;
static double skew_lim = 0.0;				       /* will be replaced by max(skew) */
static const double upper_lim = 2.807033768;		       /* set to 0.0 to disable truncation */

int GMRFLib_sn_g_get_order(void)
{
	return (ORDER);
}

double *GMRFLib_sn_g_get_coof(double skew, double *cx)
{
	// return a ALLOC'ed vector of polynomial coefficients for SKEW if CX=NULL, otherwise, return and overwrite CX

	if (!cxs) {
#pragma omp critical (Name_6f9a5166f27bbd78fbd805e2a3ed48da6b5e2cfa)
		if (!cxs) {
#define TABLE_C
#include "GMRFLib/sn-g-tables.h"
#undef TABLE_C
			assert(ORDER == 5);
			GMRFLib_spline_tp **ss = Calloc(ORDER + 1, GMRFLib_spline_tp *);
			int n = sizeof(table_skew) / sizeof(double);
			ss[0] = GMRFLib_spline_create(table_skew, table_c0, n);
			ss[1] = GMRFLib_spline_create(table_skew, table_c1, n);
			ss[2] = GMRFLib_spline_create(table_skew, table_c2, n);
			ss[3] = GMRFLib_spline_create(table_skew, table_c3, n);
			ss[4] = GMRFLib_spline_create(table_skew, table_c4, n);
			ss[5] = GMRFLib_spline_create(table_skew, table_c5, n);
			if (ISZERO(skew_lim)) {
				skew_lim = table_skew[n - 1];
			}
			cxs = ss;
		}
	}

	double *cxx = (cx ? cx : Calloc(ORDER + 1, double));
	double askew = DMIN(ABS(skew), skew_lim);
	for (int i = 0; i < ORDER + 1; i++) {
		cxx[i] = GMRFLib_spline_eval(askew, cxs[i]);
	}

	if (skew < 0.0) {
		for (int i = 0; i < ORDER + 1; i += 2) {
			cxx[i] = -cxx[i];
		}
	}

	return (cxx);
}

double *GMRFLib_sn_ginv_get_coof(double skew, double *cx)
{
	// return a ALLOC'ed vector of polynomial coefficients for SKEW if CX=NULL, otherwise, return and overwrite CX

	if (!icxs) {
#pragma omp critical (Name_cae1655d239527e4f5b52231b8bb3b6e7c511c5e)
		if (!icxs) {
#define TABLE_IC
#include "GMRFLib/sn-g-tables.h"
#undef TABLE_IC
			assert(ORDER == 5);
			GMRFLib_spline_tp **ss = Calloc(ORDER + 1, GMRFLib_spline_tp *);
			int n = sizeof(table_skew) / sizeof(double);
			ss[0] = GMRFLib_spline_create(table_skew, table_ic0, n);
			ss[1] = GMRFLib_spline_create(table_skew, table_ic1, n);
			ss[2] = GMRFLib_spline_create(table_skew, table_ic2, n);
			ss[3] = GMRFLib_spline_create(table_skew, table_ic3, n);
			ss[4] = GMRFLib_spline_create(table_skew, table_ic4, n);
			ss[5] = GMRFLib_spline_create(table_skew, table_ic5, n);
			if (ISZERO(skew_lim)) {
				skew_lim = table_skew[n - 1];
			}
			icxs = ss;
		}
	}

	double *cxx = (cx ? cx : Calloc(ORDER + 1, double));
	double askew = DMIN(ABS(skew), skew_lim);
	for (int i = 0; i < ORDER + 1; i++) {
		cxx[i] = GMRFLib_spline_eval(askew, icxs[i]);
	}

	if (skew < 0.0) {
		for (int i = 0; i < ORDER + 1; i += 2) {
			cxx[i] = -cxx[i];
		}
	}

	return (cxx);
}

double GMRFLib_sn_g_eval(double x, double *cx)
{
	// return g(x) for given polynomial coefficients cx

	double res = 0.0;
	if (upper_lim && ABS(x) > upper_lim) {
		double x0 = DSIGN(x) * upper_lim;
		double val = GMRFLib_sn_g_eval(x0, cx);
		double deriv = GMRFLib_sn_g_eval_deriv(x0, cx);
		res = val + deriv * (x - x0);
	} else {
		double fact = 1.0, pow = 1.0;
		res = cx[0];
		for (int i = 1; i < ORDER + 1; i++) {
			fact *= i;
			pow *= x;
			res += cx[i] * pow / fact;
		}
	}
	return (res);
}

double GMRFLib_sn_g_eval_deriv(double x, double *cx)
{
	// return first derivative in x, g'(x), for given polynomial coefficients cx

	double res = 0.0;
	if (upper_lim && ABS(x) > upper_lim) {
		double x0 = DSIGN(x) * upper_lim;
		res = GMRFLib_sn_g_eval_deriv(x0, cx);
	} else {
		double fact = 1.0, pow = 1.0;
		res = cx[1];
		for (int i = 2; i < ORDER + 1; i++) {
			fact *= (i - 1);
			pow *= x;
			res += cx[i] * pow / fact;
		}
	}
	return (res);
}
