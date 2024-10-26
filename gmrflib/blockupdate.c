
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

#pragma omp declare simd
static double GMRFLib_prod_diff(double a, double b, double c, double d)
{
	// return a*b-c*d , see https://pharr.org/matt/blog/2019/11/03/difference-of-floats 
	double cd = c * d;
	return fma(a, b, -cd) + fma(-c, d, cd);
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

	*a = f0;
	*b = df;
	*c = ddf;
	if (dd)
		*dd = dddf;

	if (d != 1.0) {
		*a *= d;
		*b *= d;
		*c *= d;
		if (dd) {
			*dd *= d;
		}
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

	GMRFLib_2order_approx_core(thread_id, &f0, &df, &ddf, (dd ? &dddf : NULL), x0, indx, x_vec, loglFunc, loglFunc_arg, step_len, stencil);
	if (INVALID(ddf)) {
		// if (INVALID(x0) || INVALID(f0) || INVALID(df) || INVALID(ddf)) {
		fprintf(stderr, " *** WARNING *** GMRFLib_2order_approx: rescue NAN/INF values in logl for idx=%1d\n", indx);
		f0 = df = 0.0;
		ddf = -1.0;				       /* we try with this */
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
		*a = 0.0;
		*b = 0.0;
		*c = -d * ddf;
		if (dd) {
			*dd = 0.0;
		}
	} else {
		*a = f0 + x0 * (-df + 0.5 * x0 * (ddf + 0.3333333333333333333 * dddf * x0));
		*b = df + x0 * (-ddf + 0.5 * dddf * x0);
		*c = -ddf + dddf * x0;
		if (dd) {
			*dd = dddf;
		}

		if (d != 1.0) {
			*a *= d;
			*b *= d;
			*c *= d;
			if (dd) {
				*dd *= d;
			}
		}
	}

#undef INVALID
	return GMRFLib_SUCCESS;
}

forceinline int GMRFLib_2order_approx_core(int thread_id, double *a, double *b, double *c, double *dd, double x0, int indx,
					   double *x_vec, GMRFLib_logl_tp *loglFunc, void *loglFunc_arg, double *step_len, int *stencil)
{
	// default step-size is determined using test=151. stencil=9 does not bring much...

	double step, df = 0.0, ddf = 0.0, dddf = 0.0, xx[9], f[9], f0 = 0.0, x00;
	int stenc = (stencil ? *stencil : 5);

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
	} else {

		// this is the plain code
		// df=GMRFLib_ddot(n, wf, f);
		// ddf=GMRFLib_ddot(n, wff, f);

		switch (stenc) {
		case 3:
		{
			// special implementation: ONLY used for initial values
			step = 1.0e-4;
			int n = 3;
			xx[0] = x0 - step;
			xx[1] = x0;
			xx[2] = x0 + step;

			loglFunc(thread_id, f, xx, n, indx, x_vec, NULL, loglFunc_arg, NULL);

			f0 = f[1];
			df = 0.5 * (-f[0] + f[2]);
			ddf = f[0] - 2.0 * f[1] + f[2];
		}
			break;

		case 5:
		{
			if (!step_len || ISZERO(*step_len)) {
				static double ref = GSL_DBL_EPSILON / 2.220446049e-16;
				step = ref * 5.0e-4;
			} else {
				step = *step_len;
			}

			// see https://en.wikipedia.org/wiki/Finite_difference_coefficients
			static double wf5[] = {
				// 
				0.0833333333333333333333333, -0.666666666666666666666667, 0.,
				0.666666666666666666666667, -0.0833333333333333333333333, 0.0, 0.0, 0.0,
				// 
				-0.0833333333333333333333333, 1.33333333333333333333333, -2.50000000000000000000000, 1.33333333333333333333333,
				-0.0833333333333333333333333, 0.0, 0.0, 0.0,
				// 
				-0.5, 1.0, 0.0, -1.0, 0.5, 0.0, 0.0, 0.0
			};

			int n = 5;
			int nn = 2;
			int wlength = 8;
			double *wf = wf5;
			double *wff = wf + wlength;
			double *wfff = wf + 2 * wlength;

			x00 = x0 - nn * step;
#pragma omp simd
			for (int i = 0; i < n; i++) {
				xx[i] = x00 + i * step;
			}

			loglFunc(thread_id, f, xx, n, indx, x_vec, NULL, loglFunc_arg, NULL);
			f0 = f[nn];

			int iref = n / 2L;
			double *f_ref = f + iref;
			double *wf_ref = wf + iref;
			double *wff_ref = wff + iref;

			df = GMRFLib_prod_diff(wf_ref[1], f_ref[1] - f_ref[-1], -wf_ref[2], f_ref[2] - f_ref[-2]);
			ddf = GMRFLib_prod_diff(wff_ref[1], f_ref[-1] + f_ref[1], -wff_ref[2], f_ref[-2] + f_ref[2]);
			ddf = fma(wff_ref[0], f_ref[0], ddf);

			if (dd) {
				double *wfff_ref = wfff + iref;
				dddf = GMRFLib_prod_diff(wfff_ref[1], f_ref[1] - f_ref[-1], -wfff_ref[2], f_ref[2] - f_ref[-2]);
			}
		}
			break;

		case 7:
		{
			if (!step_len || ISZERO(*step_len)) {
				static double ref = GSL_DBL_EPSILON / 2.220446049e-16;
				step = ref * 100.0e-4;
			} else {
				step = *step_len;
			}

			static double wf7[] = {
				// 
				-0.0166666666666666666666667, 0.150000000000000000000000, -0.750000000000000000000000, 0.0,
				0.750000000000000000000000, -0.150000000000000000000000, 0.016666666666666666666666, 0.0,
				// 
				0.0111111111111111111111111, -0.150000000000000000000000, 1.50000000000000000000000, -2.72222222222222222222222,
				1.50000000000000000000000, -0.150000000000000000000000, 0.0111111111111111111111111, 0.0,
				// 
				0.125, -1.0, 1.625, 0.0, -1.625, 1.0, -0.125, 0.0
			};

			int n = 7;
			int nn = 3;
			int wlength = 8;
			double *wf = wf7;
			double *wff = wf + wlength;
			double *wfff = wf + 2 * wlength;

			x00 = x0 - nn * step;
#pragma omp simd
			for (int i = 0; i < n; i++) {
				xx[i] = x00 + i * step;
			}

			loglFunc(thread_id, f, xx, n, indx, x_vec, NULL, loglFunc_arg, NULL);
			f0 = f[nn];

			int iref = n / 2L;
			double *f_ref = f + iref;
			double *wf_ref = wf + iref;
			double *wff_ref = wff + iref;

			df = GMRFLib_prod_diff(wf_ref[1], f_ref[1] - f_ref[-1], -wf_ref[2], f_ref[2] - f_ref[-2]);
			df = fma(wf_ref[3], f_ref[3] - f_ref[-3], df);

			ddf = GMRFLib_prod_diff(wff_ref[0], f_ref[0], -wff_ref[1], f_ref[1] + f_ref[-1]) +
			    GMRFLib_prod_diff(wff_ref[2], f_ref[2] + f_ref[-2], -wff_ref[3], f_ref[3] + f_ref[-3]);

			if (dd) {
				double *wfff_ref = wfff + iref;
				dddf = GMRFLib_prod_diff(wfff_ref[1], f_ref[1] - f_ref[-1], -wfff_ref[2], f_ref[2] - f_ref[-2]);
				dddf = fma(wfff_ref[3], f_ref[3] - f_ref[-3], dddf);
			}
		}
			break;

		case 9:
		{
			if (!step_len || ISZERO(*step_len)) {
				static double ref = GSL_DBL_EPSILON / 2.220446049e-16;
				step = ref * 250.0e-4;
			} else {
				step = *step_len;
			}

			static double wf9[] = {
				// 
				0.00357142857142857142857143, -0.0380952380952380952380952, 0.200000000000000000000000, -0.800000000000000000000000,
				0.,
				0.800000000000000000000000, -0.200000000000000000000000, 0.0380952380952380952380952, -0.00357142857142857142857143,
				0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,

				// 
				-0.00178571428571428571428571, 0.0253968253968253968253968, -0.200000000000000000000000, 1.60000000000000000000000,
				-2.84722222222222222222222, 1.60000000000000000000000, -0.200000000000000000000000, 0.0253968253968253968253968,
				-0.00178571428571428571428571, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,

				// 
				-0.02916666666666667, 0.3, -1.408333333333333, 2.033333333333333, 0.0, -2.033333333333333, 1.408333333333333, -0.3,
				0.02916666666666667, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
				    - 0.0291666666666666666666667, 0.300000000000000000000000, -1.40833333333333333333333,
				2.03333333333333333333333, 0.0,
				-2.03333333333333333333333, 1.40833333333333333333333, -0.300000000000000000000000, 0.0291666666666666666666667,
				0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
			};

			int n = 9;
			int nn = 4;
			int wlength = 16;
			double *wf = wf9;
			double *wff = wf + wlength;
			double *wfff = wf + 2 * wlength;

			x00 = x0 - nn * step;
#pragma omp simd
			for (int i = 0; i < n; i++) {
				xx[i] = x00 + i * step;
			}

			loglFunc(thread_id, f, xx, n, indx, x_vec, NULL, loglFunc_arg, NULL);
			f0 = f[nn];

			int iref = n / 2L;
			double *f_ref = f + iref;
			double *wf_ref = wf + iref;
			double *wff_ref = wff + iref;

			df = GMRFLib_prod_diff(wf_ref[1], f_ref[1] - f_ref[-1], -wf_ref[2], f_ref[2] - f_ref[-2]) +
			    GMRFLib_prod_diff(wf_ref[3], f_ref[3] - f_ref[-3], -wf_ref[4], f_ref[4] - f_ref[-4]);

			ddf = GMRFLib_prod_diff(wff_ref[1], f_ref[-1] + f_ref[1], -wff_ref[2], f_ref[-2] + f_ref[2]) +
			    GMRFLib_prod_diff(wff_ref[3], f_ref[-3] + f_ref[3], -wff_ref[4], f_ref[-4] + f_ref[4]);
			ddf = fma(wff_ref[0], f_ref[0], ddf);

			if (dd) {
				double *wfff_ref = wfff + iref;
				dddf = GMRFLib_prod_diff(wfff_ref[1], f_ref[1] - f_ref[-1], -wfff_ref[2], f_ref[2] - f_ref[-2]) +
				    GMRFLib_prod_diff(wfff_ref[3], f_ref[3] - f_ref[-3], -wfff_ref[4], f_ref[4] - f_ref[-4]);
			}
		}
			break;

		default:
			assert(0 == 1);
		}
	}

	double istep = 1.0 / step;
	df *= istep;
	ddf *= SQR(istep);
	if (dd) {
		dddf *= POW3(istep);
	}

	*a = f0;
	*b = df;
	*c = ddf;
	if (dd) {
		*dd = dddf;
	}

	return GMRFLib_SUCCESS;
}
