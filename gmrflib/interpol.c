
/* interpol.c
 * 
 * Copyright (C) 2011-2022 Havard Rue
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
#ifndef GITCOMMIT
#define GITCOMMIT
#endif

#include <assert.h>
#include <stddef.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>
#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

static const char GitID[] = "file: " __FILE__ "  " GITCOMMIT;

GMRFLib_spline_tp *GMRFLib_spline_create(double *x, double *y, int n)
{
	return GMRFLib_spline_create_x(x, y, n, GMRFLib_INTPOL_TRANS_NONE);
}

GMRFLib_spline_tp *GMRFLib_spline_create_x(double *x, double *y, int n,
					 GMRFLib_intpol_transform_tp trans)
{
	/*
	 * Return a spline interpolant for {(x,y)} 
	 */
	int nn = n;
	double *xx = NULL, *yy = NULL;
	double eps = GMRFLib_eps(0.5);
	GMRFLib_spline_tp *s = Calloc(1, GMRFLib_spline_tp);

	assert(n > 0);
	Calloc_init(2 * n);
	xx = Calloc_get(n);
	yy = Calloc_get(n);
	Memcpy(xx, x, n * sizeof(double));
	Memcpy(yy, y, n * sizeof(double));

	if (trans == GMRFLib_INTPOL_TRANS_P) {
		for(int i = 0; i < n; i++) {
			yy[i] = TRUNCATE(yy[i], eps, 1.0 - eps);
			yy[i] = log(yy[i]/(1.0 - yy[i]));
		}
	} else if (trans == GMRFLib_INTPOL_TRANS_Pinv) {
		for(int i = 0; i < n; i++) {
			xx[i] = TRUNCATE(xx[i], eps, 1.0 - eps);
			xx[i] = log(xx[i]/(1.0 - xx[i]));
		}
	}
	
	GMRFLib_qsorts(xx, (size_t) n, sizeof(double), yy, sizeof(double), NULL, 0, GMRFLib_dcmp);
	GMRFLib_unique_relative2(&nn, xx, yy, eps);

	s->trans = trans;
	s->xmin = xx[0];
	s->xmax = xx[nn - 1];
	s->accel = gsl_interp_accel_alloc();
	s->spline = gsl_spline_alloc(GMRFLib_density_interp_type(nn), (unsigned int) nn);
	gsl_spline_init(s->spline, xx, yy, (unsigned int) nn);

	Calloc_free();
	return s;
}

GMRFLib_spline_tp *GMRFLib_spline_create_from_matrix(GMRFLib_matrix_tp * M)
{
	/*
	 * Return a spline interpolant between {(x,y)} where x is the first column in M and y is the second column of M. 
	 */

	GMRFLib_ASSERT_RETVAL(M->ncol == 2, GMRFLib_EINVARG, NULL);
	GMRFLib_ASSERT_RETVAL(M->nrow > 0, GMRFLib_EINVARG, NULL);
	GMRFLib_ASSERT_RETVAL(M->A, GMRFLib_EINVARG, NULL);

	double *x = M->A;
	double *y = M->A + M->nrow;			       /* column-based storage */

	return GMRFLib_spline_create(x, y, M->nrow);
}

double GMRFLib_spline_eval(double x, GMRFLib_spline_tp * s)
{
	/*
	 * Evaluate a spline 's' in point 'x' 
	 */

	int extrapolate = 1;
	double xx, val;
	double eps = GMRFLib_eps(0.5);

	if (s->trans == GMRFLib_INTPOL_TRANS_Pinv) {
		xx = TRUNCATE(x, eps, 1.0-eps);
		xx = log(xx/(1.0 - xx));
	} else {
		xx = x;
	}
	
	if (xx < s->xmin || xx > s->xmax) {
		if (extrapolate) {
			// maybe I should put this into the GMRFLib_spline_tp as a parameter...
			double deriv;
			if (xx > s->xmax) {
				deriv = GMRFLib_spline_eval_deriv_x(s->xmax, s);
				val = GMRFLib_spline_eval(s->xmax, s) + deriv * (xx - s->xmax);
			} else if (xx < s->xmin) {
				deriv = GMRFLib_spline_eval_deriv_x(s->xmin, s);
				val = GMRFLib_spline_eval(s->xmin, s) + deriv * (xx - s->xmin);
			} else {
				assert(0 == 1);
			}
		} else {
			val = NAN;
		}
	} else {
		val = gsl_spline_eval(s->spline, xx, s->accel);
	}

	if (s->trans == GMRFLib_INTPOL_TRANS_P){
		val = 1.0/(1.0 + exp(-val));
	} 

	return val;
}

double GMRFLib_spline_eval_deriv(double x, GMRFLib_spline_tp * s)
{
	/*
	 * Evaluate the derivative of the spline 's' in point 'x' 
	 */

	assert(s->trans == GMRFLib_INTPOL_TRANS_NONE);
	double val;
	if (x < s->xmin || x > s->xmax) {
		val = NAN;
	} else {
		val = gsl_spline_eval_deriv(s->spline, x, s->accel);
	}

	return val;
}

double GMRFLib_spline_eval_deriv2(double x, GMRFLib_spline_tp * s)
{
	/*
	 * Evaluate the 2.derivative of the spline 's' in point 'x' 
	 */

	assert(s->trans == GMRFLib_INTPOL_TRANS_NONE);
	double val;
	if (x < s->xmin || x > s->xmax) {
		val = NAN;
	} else {
		val = gsl_spline_eval_deriv2(s->spline, x, s->accel);
	}

	return val;
}

double GMRFLib_spline_eval_deriv_x(double x, GMRFLib_spline_tp * s)
{
	double val;
	if (x < s->xmin || x > s->xmax) {
		val = NAN;
	} else {
		val = gsl_spline_eval_deriv(s->spline, x, s->accel);
	}
	return val;
}

double GMRFLib_spline_eval_deriv2_x(double x, GMRFLib_spline_tp * s)
{
	double val;
	if (x < s->xmin || x > s->xmax) {
		val = NAN;
	} else {
		val = gsl_spline_eval_deriv2(s->spline, x, s->accel);
	}
	return val;
}

int GMRFLib_spline_free(GMRFLib_spline_tp * s)
{
	/*
	 * Free spline in 's' including 's' iteself. 
	 */

	if (s) {
		gsl_spline_free(s->spline);
		gsl_interp_accel_free(s->accel);
		Free(s);
	}

	return GMRFLib_SUCCESS;
}
