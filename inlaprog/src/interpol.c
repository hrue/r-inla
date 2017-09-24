
/* interpol.c
 * 
 * Copyright (C) 2011-2014 Havard Rue
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
#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = HGVERSION;

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

#include "interpol.h"

GMRFLib_spline_tp *inla_spline_create(double *x, double *y, int n)
{
	/*
	 * Return a spline interpolant for {(x,y)} 
	 */
	int nn = n;
	double *work, *xx, *yy;
	GMRFLib_spline_tp *s = Calloc(1, GMRFLib_spline_tp);

	work = Calloc(2 * n, double);
	xx = work;
	yy = work + n;
	memcpy(xx, x, n * sizeof(double));
	memcpy(yy, y, n * sizeof(double));

	GMRFLib_qsorts(xx, (size_t) n, sizeof(double), yy, sizeof(double), NULL, 0, GMRFLib_dcmp);
	GMRFLib_unique_relative2(&nn, xx, yy, GMRFLib_eps(1. / 3.0));

	s->xmin = xx[0];
	s->xmax = xx[nn - 1];
	s->accel = gsl_interp_accel_alloc();
	s->spline = gsl_spline_alloc(GMRFLib_density_interp_type(nn), (unsigned int) nn);
	gsl_spline_init(s->spline, xx, yy, (unsigned int) nn);

	Free(work);

	return s;
}

GMRFLib_spline_tp *inla_spline_create_from_matrix(GMRFLib_matrix_tp * M)
{
	/*
	 * Return a spline interpolant between {(x,y)} where x is the first column in M and y is the second column of M. 
	 */

	GMRFLib_ASSERT_RETVAL(M->ncol == 2, GMRFLib_EINVARG, NULL);
	GMRFLib_ASSERT_RETVAL(M->nrow > 0, GMRFLib_EINVARG, NULL);
	GMRFLib_ASSERT_RETVAL(M->A, GMRFLib_EINVARG, NULL);

	double *x = M->A;
	double *y = M->A + M->nrow;			       /* column-based storage */

	return inla_spline_create(x, y, M->nrow);
}

double inla_spline_eval(double x, GMRFLib_spline_tp * s)
{
	/*
	 * Evaluate a spline 's' in point 'x' 
	 */

	int extrapolate = 1;
	double val;

	if (x < s->xmin || x > s->xmax) {
		if (extrapolate) {
			// maybe I should put this into the GMRFLib_spline_tp as a parameter...
			double deriv;
			if (x > s->xmax) {
				deriv = inla_spline_eval_deriv(s->xmax, s);
				val = inla_spline_eval(s->xmax, s) + deriv * (x - s->xmax);
			} else if (x < s->xmin) {
				deriv = inla_spline_eval_deriv(s->xmin, s);
				val = inla_spline_eval(s->xmin, s) + deriv * (x - s->xmin);
			} else {
				assert(0 == 1);
			}
		} else {
			val = NAN;
		}
	} else {
		val = gsl_spline_eval(s->spline, x, s->accel);
	}

	return val;
}

double inla_spline_eval_deriv(double x, GMRFLib_spline_tp * s)
{
	/*
	 * Evaluate the derivative of the spline 's' in point 'x' 
	 */

	double val;

	if (x < s->xmin || x > s->xmax) {
		val = NAN;
	} else {
		val = gsl_spline_eval_deriv(s->spline, x, s->accel);
	}

	return val;
}

double inla_spline_eval_deriv2(double x, GMRFLib_spline_tp * s)
{
	/*
	 * Evaluate the 2.derivative of the spline 's' in point 'x' 
	 */

	double val;

	if (x < s->xmin || x > s->xmax) {
		val = NAN;
	} else {
		val = gsl_spline_eval_deriv2(s->spline, x, s->accel);
	}

	return val;
}

int inla_spline_free(GMRFLib_spline_tp * s)
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
