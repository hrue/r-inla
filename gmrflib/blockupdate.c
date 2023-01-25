
/* blockupdate.c
 * 
 * Copyright (C) 2001-2023 Havard Rue
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

int GMRFLib_default_blockupdate_param(GMRFLib_blockupdate_param_tp ** blockupdate_par)
{
	GMRFLib_ASSERT(blockupdate_par, GMRFLib_EINVARG);

	*blockupdate_par = Calloc(1, GMRFLib_blockupdate_param_tp);
	(*blockupdate_par)->modeoption = GMRFLib_MODEOPTION_MODE;
	(*blockupdate_par)->fp = NULL;
	(*blockupdate_par)->step_len = GMRFLib_eps(0.25);
	(*blockupdate_par)->stencil = 5;

	return GMRFLib_SUCCESS;
}


int GMRFLib_2order_taylor(int thread_id, double *a, double *b, double *c, double *dd, double d, double x0, int indx,
			  double *x_vec, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, double *step_len, int *stencil)
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

	if (a) {
		*a = d * f0;
	}
	if (b) {
		*b = d * df;
	}
	if (c) {
		*c = d * ddf;
	}
	if (dd) {
		*dd = d * dddf;
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_2order_approx(int thread_id, double *a, double *b, double *c, double *dd, double d, double x0, int indx,
			  double *x_vec, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, double *step_len, int *stencil, double *cmin)
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

	if (ISZERO(d)) {
		f0 = df = ddf = dddf = 0.0;
	} else {
		if (!INVALID(x0)) {
			GMRFLib_2order_approx_core(thread_id, &f0, &df, &ddf, (dd ? &dddf : NULL), x0, indx, x_vec, loglFunc, loglFunc_arg,
						   step_len, stencil);
		}

		if (INVALID(x0) || INVALID(ddf) || INVALID(df) || INVALID(f0)) {
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
	}

	if (rescue) {
		if (a) {
			*a = 0.0;
		}
		if (b) {
			*b = 0.0;
		}
		if (c) {
			*c = -d * ddf;
		}
		if (dd) {
			*dd = 0.0;
		}
	} else {
		if (a) {
			*a = d * (f0 - df * x0 + 0.5 * ddf * SQR(x0) + 1.0 / 6.0 * dddf * gsl_pow_3(x0));
		}
		if (b) {
			*b = d * (df - ddf * x0 + dddf / 2.0 * SQR(x0));
		}
		if (c) {
			*c = -d * (ddf - dddf * x0);
		}
		if (dd) {
			*dd = d * dddf;
		}
	}

#undef INVALID
	return GMRFLib_SUCCESS;
}

int GMRFLib_2order_approx_core(int thread_id, double *a, double *b, double *c, double *dd, double x0, int indx,
			       double *x_vec, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, double *step_len, int *stencil)
{

#define ERR if (dd) {							\
		fprintf(stderr, "2order_approx_core: 3rd derivative requested but there are to few points in the stencil\n"); \
		assert(dd == NULL);					\
		exit(1);						\
	}

	double step, df = 0.0, ddf = 0.0, dddf = 0.0, xx[9], f[9], f0 = 0.0;

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
		dddf = (-1.0 / 2.0 * f[4] + 1.0 * f[3] + 0.0 * f[2] - 1.0 * f[1] + 1.0 / 2.0 * f[0]) / gsl_pow_3(step);
	} else {
		int num_points = (stencil ? *stencil : 5);
		step = (step_len && *step_len > 0.0 ? *step_len : GMRFLib_eps(1.0 / 3.9134));
		switch (num_points) {
			/*
			 * see https://en.wikipedia.org/wiki/Finite_difference_coefficients
			 */
		case 3:
		{
			xx[0] = x0 - step;
			xx[1] = x0;
			xx[2] = x0 + step;

			loglFunc(thread_id, f, xx, 3, indx, x_vec, NULL, loglFunc_arg, NULL);
			f0 = f[1];
			df = 0.5 * (f[2] - f[0]) / step;
			ddf = (f[2] - 2.0 * f[1] + f[0]) / SQR(step);
			ERR;
		}
			break;

		case 5:
		{
			double wf[] = { 1.0 / 12.0, -2.0 / 3.0, 0.0, 2.0 / 3.0, -1.0 / 12.0 };
			double wff[] = { -1.0 / 12.0, 4.0 / 3.0, -5.0 / 2.0, 4.0 / 3.0, -1.0 / 12.0 };
			double wfff[] = { -1.0 / 2.0, 1.0, 0.0, -1.0, 1.0 / 2.0 };

			xx[0] = x0 - 2.0 * step;
			xx[1] = x0 - step;
			xx[2] = x0;
			xx[3] = x0 + step;
			xx[4] = x0 + 2.0 * step;

			loglFunc(thread_id, f, xx, 5, indx, x_vec, NULL, loglFunc_arg, NULL);
			f0 = f[2];
			df = (wf[0] * f[0] + wf[1] * f[1] + wf[2] * f[2] + wf[3] * f[3] + wf[4] * f[4]) / step;
			ddf = (wff[0] * f[0] + wff[1] * f[1] + wff[2] * f[2] + wff[3] * f[3] + wff[4] * f[4]) / SQR(step);
			if (dd) {
				dddf = (wfff[0] * f[0] + wfff[1] * f[1] + wfff[2] * f[2] + wfff[3] * f[3] + wfff[4] * f[4]) / gsl_pow_3(step);
			}
		}
			break;

		case 7:
		{
			double wf[] = { -1.0 / 60.0, 3.0 / 20.0, -3.0 / 4.0, 0.0, 3.0 / 4.0, -3.0 / 20.0, 1.0 / 60.0 };
			double wff[] = { 1.0 / 90.0, -3.0 / 20.0, 3.0 / 2.0, -49.0 / 18.0, 3.0 / 2.0, -3.0 / 20.0, 1.0 / 90.0 };
			double wfff[] = { 1.0 / 8.0, -1.0, 13.0 / 8.0, 0.0, -13.0 / 8.0, 1.0, -1.0 / 8.0 };

			xx[0] = x0 - 3.0 * step;
			xx[1] = x0 - 2.0 * step;
			xx[2] = x0 - step;
			xx[3] = x0;
			xx[4] = x0 + step;
			xx[5] = x0 + 2.0 * step;
			xx[6] = x0 + 3.0 * step;

			loglFunc(thread_id, f, xx, 7, indx, x_vec, NULL, loglFunc_arg, NULL);
			f0 = f[3];
			df = (wf[0] * f[0] + wf[1] * f[1] + wf[2] * f[2] + wf[3] * f[3] + wf[4] * f[4] + wf[5] * f[5] + wf[6] * f[6]) / step;
			ddf = (wff[0] * f[0] + wff[1] * f[1] + wff[2] * f[2] + wff[3] * f[3] + wff[4] * f[4] + wff[5] * f[5] +
			       wff[6] * f[6]) / SQR(step);
			if (dd) {
				dddf = (wfff[0] * f[0] + wfff[1] * f[1] + wfff[2] * f[2] + wfff[3] * f[3] + wfff[4] * f[4] + wfff[5] * f[5] +
					wfff[6] * f[6]) / gsl_pow_3(step);
			}
		}
			break;

		case 9:
		{
			double wf[] = { 1.0 / 280.0, -4.0 / 105.0, 1.0 / 5.0, -4.0 / 5.0, 0.0, 4.0 / 5.0, -1.0 / 5.0, 4.0 / 105.0, -1.0 / 280.0 };
			double wff[] =
			    { -1.0 / 560.0, 8.0 / 315.0, -1.0 / 5.0, 8.0 / 5.0, -205.0 / 72.0, 8.0 / 5.0, -1.0 / 5.0, 8.0 / 315.0, -1.0 / 560.0 };

			double wfff[] =
			    { -7.0 / 240.0, 3.0 / 10.0, -169.0 / 120.0, 61.0 / 30.0, 0.0, -61.0 / 30.0, 169.0 / 120.0, -3.0 / 10.0, 7.0 / 240.0 };

			xx[0] = x0 - 4.0 * step;
			xx[1] = x0 - 3.0 * step;
			xx[2] = x0 - 2.0 * step;
			xx[3] = x0 - step;
			xx[4] = x0;
			xx[5] = x0 + step;
			xx[6] = x0 + 2.0 * step;
			xx[7] = x0 + 3.0 * step;
			xx[8] = x0 + 4.0 * step;

			loglFunc(thread_id, f, xx, 9, indx, x_vec, NULL, loglFunc_arg, NULL);
			f0 = f[4];
			df = (wf[0] * f[0] + wf[1] * f[1] + wf[2] * f[2] + wf[3] * f[3] + wf[4] * f[4] + wf[5] * f[5] + wf[6] * f[6] +
			      wf[7] * f[7] + wf[8] * f[8]) / step;
			ddf = (wff[0] * f[0] + wff[1] * f[1] + wff[2] * f[2] + wff[3] * f[3] + wff[4] * f[4] + wff[5] * f[5] + wff[6] * f[6] +
			       wff[7] * f[7] + wff[8] * f[8]) / SQR(step);
			if (dd) {
				dddf = (wfff[0] * f[0] + wfff[1] * f[1] + wfff[2] * f[2] + wfff[3] * f[3] + wfff[4] * f[4] + wfff[5] * f[5] +
					wfff[6] * f[6] + wfff[7] * f[7] + wfff[8] * f[8]) / gsl_pow_3(step);
			}
		}
			break;

		default:
			GMRFLib_ASSERT(num_points == 3 || num_points == 5 || num_points == 7 || num_points == 9, GMRFLib_EINVARG);
			abort();
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
