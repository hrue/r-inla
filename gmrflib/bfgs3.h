
/* bfgs3.h
 *
 *
 * EXPERIMENTAL ONLY!!!
 *
 * 
 * Copyright (C) 2009 Havard Rue
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
 *       H{\aa}vard Rue
 *       Department of Mathematical Sciences
 *       The Norwegian University of Science and Technology
 *       N-7491 Trondheim, Norway
 *       Voice: +47-7359-3533    URL  : http://www.math.ntnu.no/~hrue  
 *       Fax  : +47-7359-3524    Email: havard.rue@math.ntnu.no
 *
 * RCSId: $Id: bfgs3.h,v 1.3 2009/05/05 17:49:06 hrue Exp $
 *
 */

#ifndef __GMRFLib_BFGS3_H__
#define __GMRFLib_BFGS3_H__

#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS					       /* empty */
#define __END_DECLS					       /* empty */
#endif

__BEGIN_DECLS

/* 
 */
    typedef struct {
	gsl_function_fdf fdf_linear;
	gsl_multimin_function_fdf *fdf;
	/*
	 * fixed values 
	 */
	const gsl_vector *x;
	const gsl_vector *g;
	const gsl_vector *p;

	/*
	 * cached values, for x(alpha) = x + alpha * p 
	 */
	double f_alpha;
	double df_alpha;
	gsl_vector *x_alpha;
	gsl_vector *g_alpha;

	/*
	 * cache "keys" 
	 */
	double f_cache_key;
	double df_cache_key;
	double x_cache_key;
	double g_cache_key;
} wrapper_t;
typedef struct {
	int iter;
	double step;
	double g0norm;
	double pnorm;
	double delta_f;
	double fp0;					       /* f'(0) for f(x-alpha*p) */
	gsl_vector *x0;
	gsl_vector *g0;
	gsl_vector *p;
	/*
	 * work space 
	 */
	gsl_vector *dx0;
	gsl_vector *dg0;
	gsl_vector *x_alpha;
	gsl_vector *g_alpha;
	/*
	 * wrapper function 
	 */
	wrapper_t wrap;
	/*
	 * minimization parameters 
	 */
	double rho;
	double sigma;
	double tau1;
	double tau2;
	double tau3;
	int order;
} vector_bfgs3_state_t;

static double interp_quad(double f0, double fp0, double f1, double zl, double zh);
static double cubic(double c0, double c1, double c2, double c3, double z);
static void check_extremum(double c0, double c1, double c2, double c3, double z, double *zmin, double *fminn);
static double interp_cubic(double f0, double fp0, double f1, double fp1, double zl, double zh);
static double interpolate(double a, double fa, double fpa, double b, double fb, double fpb, double xmin, double xmax, int order);
static int minimize(gsl_function_fdf * fn, double rho, double sigma, double tau1, double tau2, double tau3, int order, double alpha1, double *alpha_new);
static void moveto(double alpha, wrapper_t * w);
static double slope(wrapper_t * w);
static double wrap_f(double alpha, void *params);
static double wrap_df(double alpha, void *params);
static void wrap_fdf(double alpha, void *params, double *f, double *df);
static void prepare_wrapper(wrapper_t * w, gsl_multimin_function_fdf * fdf,
			    const gsl_vector * x, double f, const gsl_vector * g, const gsl_vector * p, gsl_vector * x_alpha, gsl_vector * g_alpha);
static void update_position(wrapper_t * w, double alpha, gsl_vector * x, double *f, gsl_vector * g);
static void change_direction(wrapper_t * w);
static int vector_bfgs3_alloc(void *vstate, size_t n);
static int vector_bfgs3_set(void *vstate, gsl_multimin_function_fdf * fdf, const gsl_vector * x, double *f, gsl_vector * gradient, double step_size, double tol);
static void vector_bfgs3_free(void *vstate);
static int vector_bfgs3_restart(void *vstate);
static int vector_bfgs3_iterate(void *vstate, gsl_multimin_function_fdf * fdf, gsl_vector * x, double *f, gsl_vector * gradient, gsl_vector * dx);

__END_DECLS
#endif
