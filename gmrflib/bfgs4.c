
/* bfgs4.c
 *  This is a modified version of the VECTOR_BFGS2 optimiser in GSL 
 */

/* multimin/vector_bfgs2.c
 * 
 * Copyright (C) 2007 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA.
 */

/* vector_bfgs2.c -- Fletcher's implementation of the BFGS method,
   from R.Fletcher, "Practical Method's of Optimization", Second
   Edition, ISBN 0471915475.  Algorithms 2.6.2 and 2.6.4. */

/* Thanks to Alan Irwin irwin@beluga.phys.uvic.ca. for suggesting this
   algorithm and providing sample fortran benchmarks */

#include <errno.h>
#include <stddef.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_randist.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"
#include "GMRFLib/bfgs4.h"

static const int debug = 0;

/* Find a minimum in x=[0,1] of the interpolating quadratic through
 * (0,f0) (1,f1) with derivative fp0 at x=0.  The interpolating
 * polynomial is q(x) = f0 + fp0 * z + (f1-f0-fp0) * z^2
 */

static double interp_quad(double f0, double fp0, double f1, double zl, double zh)
{
	double fl = f0 + zl * (fp0 + zl * (f1 - f0 - fp0));
	double fh = f0 + zh * (fp0 + zh * (f1 - f0 - fp0));
	double c = 2 * (f1 - f0 - fp0);			       /* curvature */

	double zmin = zl, fminn = fl;

	if (fh < fminn) {
		zmin = zh;
		fminn = fh;
	}

	if (c > 0) {					       /* positive curvature required for a minimum */
		double z = -fp0 / c;			       /* location of minimum */
		if (z > zl && z < zh) {
			double f = f0 + z * (fp0 + z * (f1 - f0 - fp0));
			if (f < fminn) {
				zmin = z;
				fminn = f;
			};
		}
	}

	return zmin;
}

/* Find a minimum in x=[0,1] of the interpolating cubic through
 * (0,f0) (1,f1) with derivatives fp0 at x=0 and fp1 at x=1.
 *
 * The interpolating polynomial is:
 *
 * c(x) = f0 + fp0 * z + eta * z^2 + xi * z^3
 *
 * where eta=3*(f1-f0)-2*fp0-fp1, xi=fp0+fp1-2*(f1-f0). 
 */

static double cubic(double c0, double c1, double c2, double c3, double z)
{
	return c0 + z * (c1 + z * (c2 + z * c3));
}

static void check_extremum(double c0, double c1, double c2, double c3, double z, double *zmin, double *fminn)
{
	/*
	 * could make an early return by testing curvature >0 for minimum 
	 */

	double y = cubic(c0, c1, c2, c3, z);

	if (y < *fminn) {
		*zmin = z;				       /* accepted new point */
		*fminn = y;
	}
}

static double interp_cubic(double f0, double fp0, double f1, double fp1, double zl, double zh)
{
	double eta = 3 * (f1 - f0) - 2 * fp0 - fp1;
	double xi = fp0 + fp1 - 2 * (f1 - f0);
	double c0 = f0, c1 = fp0, c2 = eta, c3 = xi;
	double zmin, fminn;
	double z0, z1;

	zmin = zl;
	fminn = cubic(c0, c1, c2, c3, zl);
	check_extremum(c0, c1, c2, c3, zh, &zmin, &fminn);

	{
		int n = gsl_poly_solve_quadratic(3 * c3, 2 * c2, c1, &z0, &z1);

		if (n == 2) {				       /* found 2 roots */
			if (z0 > zl && z0 < zh)
				check_extremum(c0, c1, c2, c3, z0, &zmin, &fminn);
			if (z1 > zl && z1 < zh)
				check_extremum(c0, c1, c2, c3, z1, &zmin, &fminn);
		} else if (n == 1) {			       /* found 1 root */
			if (z0 > zl && z0 < zh)
				check_extremum(c0, c1, c2, c3, z0, &zmin, &fminn);
		}
	}

	return zmin;
}

static double interpolate(double a, double fa, double fpa, double b, double fb, double fpb, double xmin, double xmax, int order)
{
	/*
	 * Map [a,b] to [0,1] 
	 */
	double z, alpha, zmin, zmax;

	zmin = (xmin - a) / (b - a);
	zmax = (xmax - a) / (b - a);

	if (zmin > zmax) {
		double tmp = zmin;
		zmin = zmax;
		zmax = tmp;
	};

	if (order > 2 && GSL_IS_REAL(fpb)) {
		z = interp_cubic(fa, fpa * (b - a), fb, fpb * (b - a), zmin, zmax);
	} else {
		z = interp_quad(fa, fpa * (b - a), fb, zmin, zmax);
	}

	alpha = a + z * (b - a);

	return alpha;
}

static void moveto(double alpha, bfgs4_wrapper_t *w)
{
	if (alpha == w->x_cache_key) {			       /* using previously cached position */
		return;
	}

	/*
	 * set x_alpha = x + alpha * p 
	 */

	gsl_vector_memcpy(w->x_alpha, w->x);
	gsl_blas_daxpy(alpha, w->p, w->x_alpha);

	w->x_cache_key = alpha;
}

static double slope(bfgs4_wrapper_t *w)
{							       /* compute gradient . direction */
	double df;
	gsl_blas_ddot(w->g_alpha, w->p, &df);
	return df;
}

static double wrap_f(double alpha, void *params)
{
	bfgs4_wrapper_t *w = (bfgs4_wrapper_t *) params;
	if (alpha == w->f_cache_key) {			       /* using previously cached f(alpha) */
		return w->f_alpha;
	}

	moveto(alpha, w);

	w->f_alpha = GSL_MULTIMIN_FN_EVAL_F(w->fdf, w->x_alpha);
	w->f_cache_key = alpha;

	return w->f_alpha;
}

static double wrap_df(double alpha, void *params)
{
	bfgs4_wrapper_t *w = (bfgs4_wrapper_t *) params;
	if (alpha == w->df_cache_key) {			       /* using previously cached df(alpha) */
		return w->df_alpha;
	}

	moveto(alpha, w);

	if (alpha != w->g_cache_key) {
		GSL_MULTIMIN_FN_EVAL_DF(w->fdf, w->x_alpha, w->g_alpha);
		w->g_cache_key = alpha;
	}

	w->df_alpha = slope(w);
	w->df_cache_key = alpha;

	return w->df_alpha;
}

static void wrap_fdf(double alpha, void *params, double *f, double *df)
{
	bfgs4_wrapper_t *w = (bfgs4_wrapper_t *) params;

	/*
	 * Check for previously cached values 
	 */

	if (alpha == w->f_cache_key && alpha == w->df_cache_key) {
		*f = w->f_alpha;
		*df = w->df_alpha;
		return;
	}

	if (alpha == w->f_cache_key || alpha == w->df_cache_key) {
		*f = wrap_f(alpha, params);
		*df = wrap_df(alpha, params);
		return;
	}

	moveto(alpha, w);
	GSL_MULTIMIN_FN_EVAL_F_DF(w->fdf, w->x_alpha, &w->f_alpha, w->g_alpha);
	w->f_cache_key = alpha;
	w->g_cache_key = alpha;

	w->df_alpha = slope(w);
	w->df_cache_key = alpha;

	*f = w->f_alpha;
	*df = w->df_alpha;
}

static void
prepare_wrapper(bfgs4_wrapper_t *w, gsl_multimin_function_fdf *fdf,
		const gsl_vector *x, double f, const gsl_vector *g, const gsl_vector *p, gsl_vector *x_alpha, gsl_vector *g_alpha)
{
	w->fdf_linear.f = &wrap_f;
	w->fdf_linear.df = &wrap_df;
	w->fdf_linear.fdf = &wrap_fdf;
	w->fdf_linear.params = (void *) w;		       /* pointer to "self" */

	w->fdf = fdf;

	w->x = x;
	w->g = g;
	w->p = p;

	w->x_alpha = x_alpha;
	w->g_alpha = g_alpha;

	gsl_vector_memcpy(w->x_alpha, w->x);
	w->x_cache_key = 0.0;

	w->f_alpha = f;
	w->f_cache_key = 0.0;

	gsl_vector_memcpy(w->g_alpha, w->g);
	w->g_cache_key = 0.0;

	w->df_alpha = slope(w);
	w->df_cache_key = 0.0;
}

static void update_position(bfgs4_wrapper_t *w, double alpha, gsl_vector *x, double *f, gsl_vector *g)
{
	/*
	 * ensure that everything is fully cached 
	 */
	{
		double f_alpha, df_alpha;
		wrap_fdf(alpha, w, &f_alpha, &df_alpha);
	};

	*f = w->f_alpha;
	gsl_vector_memcpy(x, w->x_alpha);
	gsl_vector_memcpy(g, w->g_alpha);
}

static void change_direction(bfgs4_wrapper_t *w)
{
	/*
	 * Convert the cache values from the end of the current minimisation to those needed for the start of the next minimisation, alpha=0 
	 */

	/*
	 * The new x_alpha for alpha=0 is the current position 
	 */
	gsl_vector_memcpy(w->x_alpha, w->x);
	w->x_cache_key = 0.0;

	/*
	 * The function value does not change 
	 */
	w->f_cache_key = 0.0;

	/*
	 * The new g_alpha for alpha=0 is the current gradient at the endpoint 
	 */
	gsl_vector_memcpy(w->g_alpha, w->g);
	w->g_cache_key = 0.0;

	/*
	 * Calculate the slope along the new direction vector, p 
	 */
	w->df_alpha = slope(w);
	w->df_cache_key = 0.0;
}

static int vector_bfgs4_alloc(void *vstate, size_t n)
{
	vector_bfgs4_state_t *state = (vector_bfgs4_state_t *) vstate;

	state->p = gsl_vector_calloc(n);

	if (state->p == 0) {
		GSL_ERROR("failed to allocate space for p", GSL_ENOMEM);
	}

	state->x0 = gsl_vector_calloc(n);

	if (state->x0 == 0) {
		gsl_vector_free(state->p);
		GSL_ERROR("failed to allocate space for g0", GSL_ENOMEM);
	}

	state->g0 = gsl_vector_calloc(n);

	if (state->g0 == 0) {
		gsl_vector_free(state->x0);
		gsl_vector_free(state->p);
		GSL_ERROR("failed to allocate space for g0", GSL_ENOMEM);
	}

	state->dx0 = gsl_vector_calloc(n);

	if (state->dx0 == 0) {
		gsl_vector_free(state->g0);
		gsl_vector_free(state->x0);
		gsl_vector_free(state->p);
		GSL_ERROR("failed to allocate space for g0", GSL_ENOMEM);
	}

	state->dg0 = gsl_vector_calloc(n);

	if (state->dg0 == 0) {
		gsl_vector_free(state->dx0);
		gsl_vector_free(state->g0);
		gsl_vector_free(state->x0);
		gsl_vector_free(state->p);
		GSL_ERROR("failed to allocate space for g0", GSL_ENOMEM);
	}

	state->x_alpha = gsl_vector_calloc(n);

	if (state->x_alpha == 0) {
		gsl_vector_free(state->dg0);
		gsl_vector_free(state->dx0);
		gsl_vector_free(state->g0);
		gsl_vector_free(state->x0);
		gsl_vector_free(state->p);
		GSL_ERROR("failed to allocate space for g0", GSL_ENOMEM);
	}

	state->g_alpha = gsl_vector_calloc(n);

	if (state->g_alpha == 0) {
		gsl_vector_free(state->x_alpha);
		gsl_vector_free(state->dg0);
		gsl_vector_free(state->dx0);
		gsl_vector_free(state->g0);
		gsl_vector_free(state->x0);
		gsl_vector_free(state->p);
		GSL_ERROR("failed to allocate space for g0", GSL_ENOMEM);
	}

	return GSL_SUCCESS;
}

static int vector_bfgs4_set(void *vstate, gsl_multimin_function_fdf *fdf, const gsl_vector *x, double *f, gsl_vector *gradient, double step_size,
			    double tol)
{
	vector_bfgs4_state_t *state = (vector_bfgs4_state_t *) vstate;

	state->iter = 0;
	state->step = step_size;
	state->delta_f = 0;

	GSL_MULTIMIN_FN_EVAL_F_DF(fdf, x, f, gradient);

	/*
	 * Use the gradient as the initial direction 
	 */

	gsl_vector_memcpy(state->x0, x);
	gsl_vector_memcpy(state->g0, gradient);
	state->g0norm = gsl_blas_dnrm2(state->g0);

	gsl_vector_memcpy(state->p, gradient);
	gsl_blas_dscal(-1 / state->g0norm, state->p);
	state->pnorm = gsl_blas_dnrm2(state->p);	       /* should be 1 */
	state->fp0 = -state->g0norm;

	/*
	 * Prepare the wrapper 
	 */

	prepare_wrapper(&state->wrap, fdf, state->x0, *f, state->g0, state->p, state->x_alpha, state->g_alpha);

	/*
	 * Prepare 1d minimisation parameters 
	 */
	state->rho = 0.01;
	state->sigma = tol;
	state->tau1 = 9;
	state->tau2 = 0.05;
	state->tau3 = 0.5;
	state->order = 3;				       /* use cubic interpolation where possible */

	return GSL_SUCCESS;
}

static void vector_bfgs4_free(void *vstate)
{
	vector_bfgs4_state_t *state = (vector_bfgs4_state_t *) vstate;

	gsl_vector_free(state->x_alpha);
	gsl_vector_free(state->g_alpha);
	gsl_vector_free(state->dg0);
	gsl_vector_free(state->dx0);
	gsl_vector_free(state->g0);
	gsl_vector_free(state->x0);
	gsl_vector_free(state->p);
}

static int vector_bfgs4_restart(void *vstate)
{
	vector_bfgs4_state_t *state = (vector_bfgs4_state_t *) vstate;

	state->iter = 0;
	return GSL_SUCCESS;
}

static int vector_bfgs4_iterate(void *vstate, gsl_multimin_function_fdf *UNUSED(fdf),
				gsl_vector *x, double *f, gsl_vector *gradient, gsl_vector *dx)
{
	vector_bfgs4_state_t *state = (vector_bfgs4_state_t *) vstate;
	double alpha = 0.0, alpha1;
	gsl_vector *x0 = state->x0;
	gsl_vector *g0 = state->g0;
	gsl_vector *p = state->p;

	double g0norm = state->g0norm;
	double pnorm = state->pnorm;
	double delta_f = state->delta_f;
	double pg, dir;
	int status;

	double f0 = *f;

	if (pnorm == 0.0 || g0norm == 0.0 || state->fp0 == 0) {
		gsl_vector_set_zero(dx);
		return GSL_ENOPROG;
	}

	if (delta_f < 0) {
		double del = GSL_MAX_DBL(-delta_f, 10 * GSL_DBL_EPSILON * fabs(f0));
		alpha1 = GSL_MIN_DBL(1.0, 2.0 * del / (-state->fp0));
	} else {
		alpha1 = fabs(state->step);
	}

	/*
	 * line minimisation, with cubic interpolation (order = 3) 
	 */
	if (debug)
		printf("...call minimize()\n");
	status = minimize(&state->wrap.fdf_linear, state, state->rho, state->sigma, state->tau1, alpha1, &alpha);
	if (debug)
		printf("...end minimize()\n");

	if (status != GSL_SUCCESS) {
		update_position(&(state->wrap), alpha, x, f, gradient);	/* YES! hrue */
		return status;
	}

	update_position(&(state->wrap), alpha, x, f, gradient);

	state->delta_f = *f - f0;

	/*
	 * Choose a new direction for the next step 
	 */

	{
		/*
		 * This is the BFGS update: 
		 */
		/*
		 * p' = g1 - A dx - B dg 
		 */
		/*
		 * A = - (1+ dg.dg/dx.dg) B + dg.g/dx.dg 
		 */
		/*
		 * B = dx.g/dx.dg 
		 */

		gsl_vector *dx0 = state->dx0;
		gsl_vector *dg0 = state->dg0;

		double dxg, dgg, dxdg, dgnorm, A, B;

		/*
		 * dx0 = x - x0 
		 */
		gsl_vector_memcpy(dx0, x);
		gsl_blas_daxpy(-1.0, x0, dx0);

		gsl_vector_memcpy(dx, dx0);		       /* keep a copy */

		/*
		 * dg0 = g - g0 
		 */
		gsl_vector_memcpy(dg0, gradient);
		gsl_blas_daxpy(-1.0, g0, dg0);

		gsl_blas_ddot(dx0, gradient, &dxg);
		gsl_blas_ddot(dg0, gradient, &dgg);
		gsl_blas_ddot(dx0, dg0, &dxdg);

		dgnorm = gsl_blas_dnrm2(dg0);

		if (dxdg != 0) {
			B = dxg / dxdg;
			A = -(1.0 + dgnorm * dgnorm / dxdg) * B + dgg / dxdg;
		} else {
			B = 0;
			A = 0;
		}

		gsl_vector_memcpy(p, gradient);
		gsl_blas_daxpy(-A, dx0, p);
		gsl_blas_daxpy(-B, dg0, p);
	}

	gsl_vector_memcpy(g0, gradient);
	gsl_vector_memcpy(x0, x);
	state->g0norm = gsl_blas_dnrm2(g0);
	state->pnorm = gsl_blas_dnrm2(p);

	/*
	 * update direction and fp0 
	 */

	gsl_blas_ddot(p, gradient, &pg);
	dir = (pg >= 0.0) ? -1.0 : +1.0;
	gsl_blas_dscal(dir / state->pnorm, p);
	state->pnorm = gsl_blas_dnrm2(p);
	gsl_blas_ddot(p, g0, &state->fp0);

	change_direction(&state->wrap);

	return GSL_SUCCESS;
}

static const gsl_multimin_fdfminimizer_type vector_bfgs4_type = {
	"vector_bfgs4",					       /* name */
	sizeof(vector_bfgs4_state_t),
	&vector_bfgs4_alloc,
	&vector_bfgs4_set,
	&vector_bfgs4_iterate,
	&vector_bfgs4_restart,
	&vector_bfgs4_free
};

const gsl_multimin_fdfminimizer_type *gsl_multimin_fdfminimizer_vector_bfgs4 = &vector_bfgs4_type;

int bfgs4_dofit(const gsl_multifit_robust_type *T, const gsl_matrix *X, const gsl_vector *y, gsl_vector *c, gsl_matrix *cov)
{
	int s;
	gsl_multifit_robust_workspace *work = gsl_multifit_robust_alloc(T, X->size1, X->size2);
	s = gsl_multifit_robust(X, y, c, cov, work);
	gsl_multifit_robust_free(work);
	return s;
}

int gsl_bfgs4_test1(size_t n)
{
	// test-example from the GSL documentation

	size_t i;
	const size_t p = 2;				       /* linear fit */
	gsl_matrix *X, *cov;
	gsl_vector *x, *y, *c;
	const double a = 1.45;				       /* slope */
	const double b = 3.88;				       /* intercept */
	gsl_rng *r;
	X = gsl_matrix_alloc(n, p);
	x = gsl_vector_alloc(n);
	y = gsl_vector_alloc(n);
	c = gsl_vector_alloc(p);
	cov = gsl_matrix_alloc(p, p);
	r = gsl_rng_alloc(gsl_rng_default);

	for (i = 0; i < n - 3; i++) {
		double dx = 10.0 / (n - 1.0);
		double ei = gsl_rng_uniform(r);
		double xi = -5.0 + i * dx;
		double yi = a * xi + b;
		gsl_vector_set(x, i, xi);
		gsl_vector_set(y, i, yi + ei);
	}

	gsl_vector_set(x, n - 3, 4.7);
	gsl_vector_set(y, n - 3, -8.3);
	gsl_vector_set(x, n - 2, 3.5);
	gsl_vector_set(y, n - 2, -6.7);
	gsl_vector_set(x, n - 1, 4.1);
	gsl_vector_set(y, n - 1, -6.0);

	for (i = 0; i < n; ++i) {
		double xi = gsl_vector_get(x, i);
		gsl_matrix_set(X, i, 0, 1.0);
		gsl_matrix_set(X, i, 1, xi);
	}

	bfgs4_dofit(gsl_multifit_robust_bisquare, X, y, c, cov);

	for (i = 0; i < n; ++i) {
		double xi = gsl_vector_get(x, i);
		double yi = gsl_vector_get(y, i);
		gsl_vector_view v = gsl_matrix_row(X, i);
		double y_rob, y_err;
		gsl_multifit_robust_est(&v.vector, c, cov, &y_rob, &y_err);
		printf("%g %g %g\n", xi, yi, y_rob);
	}

#define C(i) (gsl_vector_get(c,(i)))
#define COV(i,j) (gsl_matrix_get(cov,(i),(j)))
	{
		printf(" best fit: Y = %g + %g X\n", C(0), C(1));
		printf(" covariance matrix:\n");
		printf(" [ %.5e, %.5e\n", COV(0, 0), COV(0, 1));
		printf("   %.5e, %.5e\n", COV(1, 0), COV(1, 1));
	}
	gsl_matrix_free(X);
	gsl_vector_free(x);
	gsl_vector_free(y);
	gsl_vector_free(c);
	gsl_matrix_free(cov);
	gsl_rng_free(r);
	return 0;
}

int bfgs4_robust_minimize(double *xmin, double *ymin, int nn, double *x, double *y, int mm, double *xd, double *yd, int order)
{
	// input n pairs of (x_i, y_i), fit a robust regression model of given order
	// and return the x* and optional y* that minimize the fitted model.

	size_t n = (size_t) nn, m = (size_t) mm, i, j, p = order + 1;
	size_t nm = n + m, idx = 0;
	gsl_matrix *X, *cov;
	gsl_vector *yy, *c;
	double xtmp, xtmp2, ytmp, x_min = 0.0, y_min = 0.0;

	X = gsl_matrix_alloc(nm, p);
	yy = gsl_vector_alloc(nm);
	c = gsl_vector_alloc(p);
	cov = gsl_matrix_alloc(p, p);

	for (i = 0; i < n; ++i) {
		gsl_vector_set(yy, idx + i, y[i]);
		double xi = x[i], xxi = xi;
		gsl_matrix_set(X, idx + i, 0, 1.0);
		for (j = 1; j < p; j++) {
			gsl_matrix_set(X, idx + i, j, xxi);
			xxi *= xi;
		}
	}

	idx = n;
	for (i = 0; i < m; ++i) {
		gsl_vector_set(yy, idx + i, yd[i]);
		double xi = xd[i], xxi = 1.0;
		gsl_matrix_set(X, idx + i, 0, 0.0);
		for (j = 1; j < p; j++) {
			gsl_matrix_set(X, idx + i, j, j * xxi);
			xxi *= xi;
		}
	}

	int err;
	gsl_set_error_handler_off();

	err = bfgs4_dofit(gsl_multifit_robust_bisquare, X, yy, c, cov);
	// err = bfgs4_dofit(gsl_multifit_robust_fair, X, yy, c, cov);
	// err = bfgs4_dofit(gsl_multifit_robust_huber, X, yy, c, cov);
	// err = bfgs4_dofit(gsl_multifit_robust_welsch, X, yy, c, cov);
	// err = bfgs4_dofit(gsl_multifit_robust_cauchy, X, yy, c, cov);

	gsl_set_error_handler(NULL);
	if (err == GSL_EMAXITER) {
		int iidx;
		GMRFLib_min_value(y, nn, &iidx);
		*xmin = x[iidx];
		if (ymin) {
			*ymin = y[iidx];
		}
		return GMRFLib_SUCCESS;
	}

	size_t mp = 50;
	double dx = 0.0, xx_max = 0.0, xx_min = 0.0;

	if (xd) {
		xx_max = DMAX(GMRFLib_max_value(x, nn, NULL), GMRFLib_max_value(xd, mm, NULL));
		xx_min = DMIN(GMRFLib_min_value(x, nn, NULL), GMRFLib_min_value(xd, mm, NULL));
	} else {
		xx_max = GMRFLib_max_value(x, nn, NULL);
		xx_min = GMRFLib_min_value(x, nn, NULL);
	}
	dx = (xx_max - xx_min) / (mp - 1.0);

	for (i = 0; i < mp; ++i) {
		xtmp = xtmp2 = xx_min + dx * i;
		ytmp = gsl_vector_get(c, 0);
		for (j = 1; j < p; j++) {
			ytmp += gsl_vector_get(c, j) * xtmp2;
			xtmp2 *= xtmp;
		}
		if (i == 0 || ytmp < y_min) {
			x_min = xtmp;
			y_min = ytmp;
		}
	}

	int max_iter = 10;
	for (int iter = 0; iter < max_iter; iter++) {
		double val, grad;

		xtmp = x_min;
		val = gsl_vector_get(c, 1);
		for (j = 2; j < p; j++) {
			val += j * gsl_vector_get(c, j) * xtmp;
			xtmp *= x_min;
		}

		xtmp = x_min;
		grad = 2.0 * gsl_vector_get(c, 2);
		for (j = 3; j < p; j++) {
			grad += j * (j - 1.0) * gsl_vector_get(c, j) * xtmp;
			xtmp *= x_min;
		}

		double ddx = -val / grad;
		x_min += ddx;
		if (iter == max_iter - 1 || ABS(ddx) < GSL_ROOT3_DBL_EPSILON) {
			break;
		}
	}

	y_min = gsl_vector_get(c, 0);
	xtmp = x_min;
	for (j = 1; j < p; j++) {
		y_min += gsl_vector_get(c, j) * xtmp;
		xtmp *= x_min;
	}

	gsl_matrix_free(X);
	gsl_vector_free(yy);
	gsl_vector_free(c);
	gsl_matrix_free(cov);

	*xmin = x_min;
	if (ymin) {
		*ymin = y_min;
	}

	return GMRFLib_SUCCESS;
}

static int minimize(gsl_function_fdf *fn, vector_bfgs4_state_t *state, double rho, double sigma, double tau1, double alpha1, double *alpha_new)
{
	double f0, fp0, falpha, falpha_prev, fpalpha, fpalpha_prev, delta, alpha_next;
	double alpha = alpha1, alpha_prev = 0.0;
	double a, b;
	const size_t bracket_iters = 25;
	size_t i = 0, j;

	double *dfalphas = Calloc(bracket_iters + 1, double);
	double *dfunval = Calloc(bracket_iters + 1, double);
	size_t ndfunval = 0;

	if (debug)
		printf("...enter minimize() sigma = %.12g\n", sigma);

	GSL_FN_FDF_EVAL_F_DF(fn, 0.0, &f0, &fp0);
	if (debug)
		printf("..eval F_DF %g %g \n", f0, fp0);

	dfalphas[ndfunval] = 0.0;
	dfunval[ndfunval++] = fp0;

	falpha_prev = f0;
	fpalpha_prev = fp0;

	/*
	 * Avoid uninitialized variables morning 
	 */
	a = 0.0;
	b = alpha;

	/*
	 * Begin bracketing 
	 */

	while (i++ < bracket_iters) {
		if (debug)
			printf("...begin bracketing\n");

		falpha = GSL_FN_FDF_EVAL_F(fn, alpha);
		if (debug)
			printf("...begin bracketing: eval f %.12g\n", falpha);

		/*
		 * Fletcher's rho test 
		 */

		if (falpha > f0 + alpha * rho * fp0 || falpha >= falpha_prev) {
			a = alpha_prev;
			b = alpha;
			break;				       /* goto sectioning */
		}

		fpalpha = GSL_FN_FDF_EVAL_DF(fn, alpha);
		if (debug)
			printf("...begin bracketing: eval df %.12g\n", falpha);

		dfalphas[ndfunval] = alpha;
		dfunval[ndfunval++] = fpalpha;

		/*
		 * Fletcher's sigma test 
		 */
		if (fabs(fpalpha) <= -sigma * fp0) {
			*alpha_new = alpha;
			return GSL_SUCCESS;
		}

		if (fpalpha >= 0) {
			a = alpha;
			b = alpha_prev;
			break;				       /* goto sectioning */
		}

		delta = alpha - alpha_prev;

		{
			double lower = alpha + delta;
			double upper = alpha + tau1 * delta;
			alpha_next = interpolate(alpha_prev, falpha_prev, fpalpha_prev, alpha, falpha, fpalpha, lower, upper, state->order);
		}

		alpha_prev = alpha;
		falpha_prev = falpha;
		fpalpha_prev = fpalpha;
		alpha = alpha_next;
	}

	/*
	 * Sectioning of bracket [a,b] 
	 */

	size_t dim = state->p->size;
	size_t num_threads = GMRFLib_openmp->max_threads_outer;
	size_t na = IMIN(16, IMAX(num_threads, 8));	       /* since 'zero' is already computed */

	double *aa = Calloc(na, double);
	double *fun = Calloc(na, double);
	double **thetas = Calloc(na, double *);
	double *pos = Calloc(na, double);

	// layout points on [0,1], adding one left point outside for stability. make points close to each other close to 0 compared to 1. since
	// the gradient information is available now at alpha=0, we can reduce the two points to the left into one point to the left
	int idx_zero = 1;
	double dzero = (double) idx_zero;

	pos[0] = -0.15;
	pos[idx_zero] = 0.0;
	for (i = idx_zero + 1; i < na; i++) {
		pos[i] = SQR((i - dzero) / ((na - 1.0) - dzero));
	}
	GMRFLib_scale_vector(pos, na);

	for (i = 0; i < na; i++) {
		aa[i] = a + (b - a) * pos[i];
		thetas[i] = Calloc(dim, double);
		for (j = 0; j < dim; j++) {
			thetas[i][j] = gsl_vector_get(state->x0, j) + aa[i] * gsl_vector_get(state->p, j);
		}
	}

	int ierr = 0;
	GMRFLib_opt_f_omp(thetas, na, fun, &ierr);

	// remove possible Inf values (do nan's in the same round...)
	for (i = j = 0; i < na; i++) {
		if (!(ISINF(fun[i]) || ISNAN(fun[i]))) {
			if (j < i) {
				fun[j] = fun[i];
				Memcpy(thetas[j], thetas[i], dim * sizeof(double));
			}
			j++;
		}
	}
	na = j;
	assert(na > 2);

	// include derivative information within the range of alphas
	double amax = GMRFLib_max_value(aa, na, NULL);
	double amin = GMRFLib_min_value(aa, na, NULL);
	int nd = 0;
	for (i = j = 0; i < ndfunval; i++) {
		if ((amin <= dfalphas[i]) && (dfalphas[i] <= amax)) {
			dfalphas[j] = dfalphas[i];
			dfunval[j] = dfunval[i];
			j++;
			nd++;
		}
	}

	if (1 || debug) {
		printf("\tparallel linesearch:\n");
		for (i = 0; i < na; i++) {
			printf("\t\tfun  i=%1zu aa=%6.2f ", i, aa[i]);
			for (j = 0; j < dim; j++) {
				printf(" %6.2f", thetas[i][j]);
			}
			printf(" %10.6f\n", fun[i]);
		}
		for (i = 0; i < (size_t) nd; i++) {
			printf("\t\tdfun i=%1zu aa=%6.2f ", i, dfalphas[i]);
			for (j = 0; j < dim; j++) {
				printf(" %6s", "");
			}
			printf("  %10.6f\n", dfunval[i]);
		}
	}

	double fmin, aa_min;
	int robust_regression = 1, order = 2;

	bfgs4_robust_minimize(&aa_min, &fmin, na, aa, fun, nd, dfalphas, dfunval, order);
	if (aa_min > amax || aa_min < amin) {
		int idx_min;
		GMRFLib_min_value(fun, na, &idx_min);
		aa_min = aa[idx_min];
		robust_regression = 0;
	}

	if (1 || debug) {
		printf("\tamin %f fmin %f (%s)\n", aa_min, fmin,
		       (robust_regression ? "robust regression, internal minimum" : "enable emergency mode"));
	}
	*alpha_new = aa_min;

	for (i = 0; i < na; i++) {
		Free(thetas[i]);
	}
	Free(thetas);
	Free(aa);
	Free(fun);
	Free(pos);

	Free(dfalphas);
	Free(dfunval);

	return GSL_SUCCESS;
}
