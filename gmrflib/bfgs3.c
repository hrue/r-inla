
/* bfgs3.c
 *
 * EXPERIMENTAL ONLY!!!
 *
 *
 *  This is a modified version of the VECTOR_BFGS2 optimiser in GSL 
 */

#ifndef HGVERSION
#define HGVERSION
#endif
//static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

/* Pre-hg-Id: $Id: bfgs3.c,v 1.8 2009/12/15 12:26:03 hrue Exp $ */


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
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_poly.h>
#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"
#include "GMRFLib/bfgs3.h"
static int debug = 0;


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

/* recommended values from Fletcher are 
   rho = 0.01, sigma = 0.1, tau1 = 9, tau2 = 0.05, tau3 = 0.5 */

static int minimize(gsl_function_fdf * fn, double rho, double sigma, double tau1, double tau2, double tau3, int order, double alpha1,
		    double *alpha_new)
{
	double f0, fp0, falpha, falpha_prev, fpalpha, fpalpha_prev, delta, alpha_next;
	double alpha = alpha1, alpha_prev = 0.0;
	double a, b, fa, fb, fpa, fpb;
	const size_t bracket_iters = 100, section_iters = 100;
	size_t i = 0;

	if (debug)
		printf("...enter minimize() sigma = %.12g\n", sigma);

	GSL_FN_FDF_EVAL_F_DF(fn, 0.0, &f0, &fp0);
	if (debug)
		printf("..eval F_DF %g %g \n", f0, fp0);

	falpha_prev = f0;
	fpalpha_prev = fp0;

	/*
	 * Avoid uninitialized variables morning 
	 */
	a = 0.0;
	b = alpha;
	fa = f0;
	fb = 0.0;
	fpa = fp0;
	fpb = 0.0;

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

		if (falpha > f0 + alpha * rho * fp0 || falpha >= falpha_prev || GMRFLib_request_optimiser_to_stop) {
			a = alpha_prev;
			fa = falpha_prev;
			fpa = fpalpha_prev;
			b = alpha;
			fb = falpha;
			fpb = GSL_NAN;
			break;				       /* goto sectioning */
		}

		fpalpha = GSL_FN_FDF_EVAL_DF(fn, alpha);
		if (debug)
			printf("...begin bracketing: eval df %.12g\n", falpha);

		/*
		 * Fletcher's sigma test 
		 */

		if (fabs(fpalpha) <= -sigma * fp0 || GMRFLib_request_optimiser_to_stop) {
			*alpha_new = alpha;
			return GSL_SUCCESS;
		}

		if (fpalpha >= 0 || GMRFLib_request_optimiser_to_stop) {
			a = alpha;
			fa = falpha;
			fpa = fpalpha;
			b = alpha_prev;
			fb = falpha_prev;
			fpb = fpalpha_prev;
			break;				       /* goto sectioning */
		}

		delta = alpha - alpha_prev;

		{
			double lower = alpha + delta;
			double upper = alpha + tau1 * delta;

			alpha_next = interpolate(alpha_prev, falpha_prev, fpalpha_prev, alpha, falpha, fpalpha, lower, upper, order);
		}

		alpha_prev = alpha;
		falpha_prev = falpha;
		fpalpha_prev = fpalpha;
		alpha = alpha_next;

		if (GMRFLib_request_optimiser_to_stop)
			break;
	}

	if (GMRFLib_request_optimiser_to_stop) {
		*alpha_new = alpha;
		return GSL_SUCCESS;			       /* terminate */
	}

	/*
	 * Sectioning of bracket [a,b] 
	 */

	while (i++ < section_iters) {

		if (debug)
			printf("...sectioning of bracket b %.12g a %.12g alpha %.12g\n", b, a, alpha);

		delta = b - a;

		{
			double lower = a + tau2 * delta;
			double upper = b - tau3 * delta;

			alpha = interpolate(a, fa, fpa, b, fb, fpb, lower, upper, order);
		}

		falpha = GSL_FN_FDF_EVAL_F(fn, alpha);
		if (debug)
			printf("...eval F %.12g\n", falpha);

		if (debug)
			printf("... roundoff check %.12g\n", (a - alpha) * fpa);

		if ((a - alpha) * fpa <= GMRFLib_eps(2.0 / 3.0)) {	/* hrue */
			/*
			 * roundoff prevents progress 
			 */
			if (debug)
				printf("BFGS3: minimizer: abort search as we do not seem to get any longer...\n");

			*alpha_new = alpha;		       /* hrue */
			return GSL_ENOPROG;
		}

		if (debug)
			printf("...TEST %.12g > %.12g || %.12g >= %.12g\n", falpha, f0 + rho * alpha * fp0, falpha, fa);

		if (GMRFLib_request_optimiser_to_stop) {
			*alpha_new = alpha;
			return GSL_SUCCESS;		       /* terminate */
		}

		if (falpha > f0 + rho * alpha * fp0 || falpha >= fa) {
			/*
			 * a_next = a; 
			 */
			b = alpha;
			fb = falpha;
			fpb = GSL_NAN;
		} else {
			fpalpha = GSL_FN_FDF_EVAL_DF(fn, alpha);
			if (debug)
				printf("...eval DF %.12g\n", fpalpha);

			if (debug)
				printf("... TEST %.12g <= %.12g\n", fabs(fpalpha), -sigma * fp0);
			if (fabs(fpalpha) <= -sigma * fp0 || GMRFLib_request_optimiser_to_stop) {
				*alpha_new = alpha;
				return GSL_SUCCESS;	       /* terminate */
			}

			if (((b - a) >= 0 && fpalpha >= 0) || ((b - a) <= 0 && fpalpha <= 0)) {
				b = a;
				fb = fa;
				fpb = fpa;
				a = alpha;
				fa = falpha;
				fpa = fpalpha;
			} else {
				a = alpha;
				fa = falpha;
				fpa = fpalpha;
			}
		}
	}

	return GSL_SUCCESS;
}

static void moveto(double alpha, wrapper_t * w)
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

static double slope(wrapper_t * w)
{							       /* compute gradient . direction */
	double df;
	gsl_blas_ddot(w->g_alpha, w->p, &df);
	return df;
}

static double wrap_f(double alpha, void *params)
{
	wrapper_t *w = (wrapper_t *) params;
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
	wrapper_t *w = (wrapper_t *) params;
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
	wrapper_t *w = (wrapper_t *) params;

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
prepare_wrapper(wrapper_t * w, gsl_multimin_function_fdf * fdf,
		const gsl_vector * x, double f, const gsl_vector * g, const gsl_vector * p, gsl_vector * x_alpha, gsl_vector * g_alpha)
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

static void update_position(wrapper_t * w, double alpha, gsl_vector * x, double *f, gsl_vector * g)
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

static void change_direction(wrapper_t * w)
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

static int vector_bfgs3_alloc(void *vstate, size_t n)
{
	vector_bfgs3_state_t *state = (vector_bfgs3_state_t *) vstate;

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

static int vector_bfgs3_set(void *vstate, gsl_multimin_function_fdf * fdf, const gsl_vector * x, double *f, gsl_vector * gradient, double step_size,
			    double tol)
{
	vector_bfgs3_state_t *state = (vector_bfgs3_state_t *) vstate;

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

static void vector_bfgs3_free(void *vstate)
{
	vector_bfgs3_state_t *state = (vector_bfgs3_state_t *) vstate;

	gsl_vector_free(state->x_alpha);
	gsl_vector_free(state->g_alpha);
	gsl_vector_free(state->dg0);
	gsl_vector_free(state->dx0);
	gsl_vector_free(state->g0);
	gsl_vector_free(state->x0);
	gsl_vector_free(state->p);
}

static int vector_bfgs3_restart(void *vstate)
{
	vector_bfgs3_state_t *state = (vector_bfgs3_state_t *) vstate;

	state->iter = 0;
	return GSL_SUCCESS;
}

static int vector_bfgs3_iterate(void *vstate, gsl_multimin_function_fdf * fdf, gsl_vector * x, double *f, gsl_vector * gradient, gsl_vector * dx)
{
	vector_bfgs3_state_t *state = (vector_bfgs3_state_t *) vstate;
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
	status = minimize(&state->wrap.fdf_linear, state->rho, state->sigma, state->tau1, state->tau2, state->tau3, state->order, alpha1, &alpha);
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

static const gsl_multimin_fdfminimizer_type vector_bfgs3_type = {
	"vector_bfgs3",					       /* name */
	sizeof(vector_bfgs3_state_t),
	&vector_bfgs3_alloc,
	&vector_bfgs3_set,
	&vector_bfgs3_iterate,
	&vector_bfgs3_restart,
	&vector_bfgs3_free
};

const gsl_multimin_fdfminimizer_type *gsl_multimin_fdfminimizer_vector_bfgs3 = &vector_bfgs3_type;
