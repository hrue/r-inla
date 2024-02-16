
/* blockupdate.c
 * 
 * Copyright (C) 2001-2024 Havard Rue
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

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

forceinline double PROD_DIFF(double a, double b, double c, double d)
{
	// return a*b-c*d
	// reference: https://pharr.org/matt/blog/2019/11/03/difference-of-floats 
	double cd = c * d;
	double err = fma(-c, d, cd);
	double dop = fma(a, b, -cd);
	return dop + err;
}

int GMRFLib_default_blockupdate_param(GMRFLib_blockupdate_param_tp **blockupdate_par)
{
	GMRFLib_ASSERT(blockupdate_par, GMRFLib_EINVARG);

	*blockupdate_par = Calloc(1, GMRFLib_blockupdate_param_tp);
	(*blockupdate_par)->modeoption = GMRFLib_MODEOPTION_MODE;
	(*blockupdate_par)->fp = NULL;
	(*blockupdate_par)->step_len = GSL_ROOT4_DBL_EPSILON;
	(*blockupdate_par)->stencil = 5;

	return GMRFLib_SUCCESS;
}


int GMRFLib_2order_taylor(int thread_id, double *a, double *b, double *c, double *dd, double d, double x0, int indx,
			  double *x_vec, GMRFLib_logl_tp *loglFunc, void *loglFunc_arg, double *step_len, int *stencil)
{
	/*
	 * compute a,b,c in the taylor expansion around x0 of d*loglFunc(x0,...)
	 * 
	 * a + b*(x-x0) + 0.5*c*(x-x0)^2 + 1/6*dd*(x-x0)^3
	 * 
	 */
	double f0 = 0.0, df = 0.0, ddf = 0.0, dddf = 0.0;

	if (ISZERO(d)) {
		f0 = df = ddf = 0.0;
	} else {
		GMRFLib_2order_approx_core(thread_id, &f0, &df, &ddf, (dd ? &dddf : NULL), x0, indx, x_vec, loglFunc, loglFunc_arg, step_len,
					   stencil);
	}

	if (a)
		*a = f0;
	if (b)
		*b = df;
	if (c)
		*c = ddf;
	if (dd)
		*dd = dddf;

	if (d != 1.0) {
		if (a)
			*a *= d;
		if (b)
			*b *= d;
		if (c)
			*c *= d;
		if (dd)
			*dd *= d;
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_2order_approx(int thread_id, double *a, double *b, double *c, double *dd, double d, double x0, int indx,
			  double *x_vec, GMRFLib_logl_tp *loglFunc, void *loglFunc_arg, double *step_len, int *stencil, double *cmin)
{
	/*
	 * compute a,b,c in the taylor expansion around x0 of d*loglFunc(x0,...)
	 * 
	 * a + b*x - 0.5*c*x^2 + 1/6*dd*x^3
	 *
	 * where cmin is the minimum value of c.
	 */

	/*
	 * > A:=collect(expand(a + b * (x-x0) + 1/2 \                                      
	 * > * c * (x-x0)^2 + 1/6 * d * (x-x0)^3), [x,x^2, x^3]);
	 *          3                    /            2    \             2              3
	 *       d x   /      d x0\  2   |        d x0     |         c x0           d x0
	 * A := ---- + |c/2 - ----| x  + |-c x0 + ----- + b| x + a + ----- - b x0 - -----
	 *             \       2  /      \          2      /           2              6
	 *
	 * > coeff(A,x);                                                                   
	 *             2
	 *         d x0
	 * -c x0 + ----- + b
	 *           2
	 * 
	 * > coeff(A,x^2);
	 *       d x0
	 * c/2 - ----
	 *        2
	 * 
	 * > coeff(A,x^3);
	 * d/6
	 * 
	 */

#define INVALID(x_) (ISNAN(x_) || ISINF(x_))

	double f0 = 0.0, df = 0.0, ddf = 0.0, dddf = 0.0;
	int rescue = 0;

	GMRFLib_2order_approx_core(thread_id, &f0, &df, &ddf, (dd ? &dddf : NULL), x0, indx, x_vec, loglFunc, loglFunc_arg,
				   step_len, stencil);
	if (INVALID(ddf)) {
		//if (INVALID(x0) || INVALID(f0) || INVALID(df) || INVALID(ddf)) {
		fprintf(stderr, "GMRFLib_2order_approx: rescue NAN/INF values in logl for idx=%1d\n", indx);
		f0 = df = 0.0;
		ddf = -1.0;			       /* we try with this */
		if (dd) {
			dddf = 0.0;
		}
		rescue = 1;
	} else {
		if (cmin) {
			ddf = DMIN(-(*cmin), ddf);
		}
	}

	if (rescue) {
		if (a)
			*a = 0.0;
		if (b)
			*b = 0.0;
		if (c)
			*c = -d * ddf;
		if (dd)
			*dd = 0.0;
	} else {
		if (a)
			*a = f0 + x0 * (-df + 0.5 * x0 * (ddf + 0.3333333333333333333 * dddf * x0));
		if (b)
			*b = df + x0 * (-ddf + 0.5 * dddf * x0);
		if (c)
			*c = -ddf + dddf * x0;
		if (dd)
			*dd = dddf;

		if (d != 1.0) {
			if (a)
				*a *= d;
			if (b)
				*b *= d;
			if (c)
				*c *= d;
			if (dd)
				*dd *= d;
		}
	}

#undef INVALID
	return GMRFLib_SUCCESS;
}

forceinline int GMRFLib_2order_approx_core(int thread_id, double *a, double *b, double *c, double *dd, double x0, int indx,
					   double *x_vec, GMRFLib_logl_tp *loglFunc, void *loglFunc_arg, double *step_len, int *stencil)
{
	double step, df = 0.0, ddf = 0.0, dddf = 0.0, xx[9], f[9], f0 = 0.0, x00;

	if (step_len && *step_len < 0.0) {
		/*
		 * for internal use only! 
		 */
		step = -(*step_len);

		xx[0] = x0 - 2 * step;
		xx[1] = x0 - step;
		xx[2] = x0;
		xx[3] = x0 + step;
		xx[4] = x0 + 2 * step;

		loglFunc(thread_id, f, xx, 5, indx, x_vec, NULL, loglFunc_arg, NULL);

		f0 = f[2];
		df = (1.0 / 12.0 * f[4] - 2.0 / 3.0 * f[3] + 0.0 * f[2] + 2.0 / 3.0 * f[1] - 1.0 / 12.0 * f[0]) / step;
		ddf = (-1.0 / 12.0 * f[4] - 4.0 / 3.0 * f[3] - 5.0 / 2.0 * f[2] + 4.0 / 3.0 * f[1] - 1.0 / 12.0 * f[0]) / SQR(step);
		dddf = (-1.0 / 2.0 * f[4] + 1.0 * f[3] + 0.0 * f[2] - 1.0 * f[1] + 1.0 / 2.0 * f[0]) / POW3(step);
	} else if (stencil && *stencil == 3) {
		// special implementation for *stencil==3, which is ONLY used for initial values
		step = 1.0e-4;
		int n = 3;
		xx[0] = x0 - step;
		xx[1] = x0;
		xx[2] = x0 + step;

		loglFunc(thread_id, f, xx, n, indx, x_vec, NULL, loglFunc_arg, NULL);

		f0 = f[1];
		df = 0.5 * (-f[0] + f[2]);
		ddf = f[0] - 2.0 * f[1] + f[2];

		df /= step;
		ddf /= SQR(step);
		dddf = 0.0;
	} else {
		// DOES NOT WORK FOR n==3!

		int num_points = (stencil ? *stencil : 5);

		if (!step_len || ISZERO(*step_len)) {
			static double ref = GSL_DBL_EPSILON / 2.220446049e-16;
			step = ref * (*stencil == 5 || *stencil == 3 ? 1.0e-4 : (*stencil == 7 ? 5.0e-4 : 1.0e-3));
		} else {
			step = *step_len;
		}

		// see https://en.wikipedia.org/wiki/Finite_difference_coefficients
		int n, nn, wlength;

		/*
		 * static double wf3[] = { // -0.5, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, // 1.0, -2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, // 0.0, 0.0, 
		 * 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }; 
		 */
		static double wf5[] = {
			// 
			0.08333333333333333, -0.6666666666666666, 0.0, 0.6666666666666666, -0.08333333333333333, 0.0, 0.0, 0.0,
			// 
			-0.08333333333333333, 1.333333333333333, -2.5, 1.333333333333333, -0.08333333333333333, 0.0, 0.0, 0.0,
			// 
			-0.5, 1.0, 0.0, -1.0, 0.5, 0.0, 0.0, 0.0
		};

		static double wf7[] = {
			// 
			-0.01666666666666667, 0.15, -0.75, 0.0, 0.75, -0.15, 0.01666666666666667, 0.0,
			// 
			0.01111111111111111, -0.15, 1.5, -2.722222222222222, 1.5, -0.15, 0.01111111111111111, 0.0,
			// 
			0.125, -1.0, 1.625, 0.0, -1.625, 1.0, -0.125, 0.0
		};

		static double wf9[] = {
			// 
			0.003571428571428571, -0.0380952380952381, 0.2, -0.8, 0.0, 0.8, -0.2, 0.0380952380952381, -0.003571428571428571, 0.0, 0.0,
			0.0, 0.0, 0.0, 0.0, 0.0,
			// 
			-0.001785714285714286, 0.0253968253968254, -0.2, 1.6, -2.847222222222222, 1.6, -0.2, 0.0253968253968254,
			-0.001785714285714286, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			// 
			-0.02916666666666667, 0.3, -1.408333333333333, 2.033333333333333, 0.0, -2.033333333333333, 1.408333333333333, -0.3,
			0.02916666666666667, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
		};

		double *wf = NULL;

		switch (num_points) {
			/*
			 * case 3: n = 3; nn = 1; wlength = 8; wf = wf3; break; 
			 */
		case 5:
			n = 5;
			nn = 2;
			wlength = 8;
			wf = wf5;
			break;
		case 7:
			n = 7;
			nn = 3;
			wlength = 8;
			wf = wf7;
			break;
		case 9:
			n = 9;
			nn = 4;
			wlength = 16;
			wf = wf9;
			break;
		default:
			GMRFLib_ASSERT(num_points == 5 || num_points == 7 || num_points == 9, GMRFLib_EINVARG);
			abort();
		}

		double *wff = wf + wlength;
		double *wfff = wf + 2 * wlength;

		x00 = x0 - nn * step;
#pragma omp simd
		for (int i = 0; i < n; i++) {
			xx[i] = x00 + i * step;
		}

		loglFunc(thread_id, f, xx, n, indx, x_vec, NULL, loglFunc_arg, NULL);
		f0 = f[nn];

		if (0) {
			// this is the plain code
			df = GMRFLib_ddot(n, wf, f);
			ddf = GMRFLib_ddot(n, wff, f);
			if (dd) {
				dddf = GMRFLib_ddot(n, wfff, f);
			}
		} else {
			// better code to protect for rounding errors, as we know how the signs change.
			// see www-reference in PROD_DIFF function

			int iref = n / 2L;
			double *f_ref = f + iref;
			double *wf_ref = wf + iref;
			double *wff_ref = wff + iref;

			// we know the sign is alternating with the same abs(coof)
			df = 0.0;
#pragma omp simd reduction(+: df)
			for (int i = 1; i < n - iref; i++) {
				double _a = f_ref[i], _b = f_ref[-i];
				df += PROD_DIFF(wf_ref[i], _a, wf_ref[i], _b);
			}

			// abs(coof) is the same but with oposite sign
			ddf = wff_ref[0] * f_ref[0];
			switch (n) {
			case 5L:
				ddf += PROD_DIFF(wff_ref[1], f_ref[-1] + f_ref[1], -wff_ref[2], f_ref[-2] + f_ref[2]);
				break;
			case 7L:
				ddf += PROD_DIFF(wff_ref[1], f_ref[-1] + f_ref[1], -wff_ref[2], f_ref[-2] + f_ref[2]) +
				    PROD_DIFF(wff_ref[3], f_ref[3], -wff_ref[3], f_ref[-3]);
				break;
			case 9L:
				ddf += PROD_DIFF(wff_ref[1], f_ref[-1] + f_ref[1], -wff_ref[2], f_ref[-2] + f_ref[2]) +
				    PROD_DIFF(wff_ref[3], f_ref[-3] + f_ref[3], -wff_ref[4], f_ref[-4] + f_ref[4]);
				break;
			default:
				assert(0 == 1);
			}

			if (dd) {
				double *wfff_ref = wfff + iref;
				// we know the sign is alternating with the same abs(coof)
				dddf = 0.0;
#pragma omp simd reduction(+: dddf)
				for (int i = 1; i < n - iref; i++) {
					double _a = f_ref[i], _b = f_ref[-i];
					dddf += PROD_DIFF(wfff_ref[i], _a, wfff_ref[i], _b);
				}
			}
		}

		df /= step;
		ddf /= SQR(step);
		if (dd) {
			dddf /= POW3(step);
		}
	}

	*a = f0;
	*b = df;
	*c = ddf;
	if (dd) {
		*dd = dddf;
	}

	return GMRFLib_SUCCESS;
}
