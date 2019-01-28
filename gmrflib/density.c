
/* density.c
 * 
 * Copyright (C) 2006-2013 Havard Rue
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
  \file density.c
  \brief Functions to represent densities and compute properties of them

  These functions represent densities which are reasonable close a Gaussian, by fitting a Gaussian, Skew-Normal distribution, or
  a more general construct; a Gaussian times exp of a smooth spline.  The use in GMRFLib is to represent posterior marginals for
  for approximate inference (see \ref approx-inference.c).

  These routines assumes that the variable itself is standardised, with mean \c GMRFLib_density_tp::std_mean and standard deviation
  \c GMRFLib_density_tp::std_stdev, and all properties computed or extracted from a \c GMRFLib_density_tp -object, is with respect
  to this standardisation.

  Example: Let \f$\mu\f$ and \f$\sigma\f$ be the \c GMRFLib_density_tp::std_mean and \c GMRFLib_density_tp::std_stdev used for
  standardisation. Then the \c GMRFLib_density_tp -object represent the density for \f$\tilde{X}\f$, where
  \f[ \tilde{X} = \frac{X - \mu}{\sigma}. \f]

  Functions that extract or compute properties from a \c GMRFLib_density_tp -object, are the following.
  - \c GMRFLib_evaluate_logdensity()  Compute the log of the density at a given position.
  - \c GMRFLib_evaluate_density() Compute the density at a given position.
  - \c GMRFLib_density_Pinv() Compute quantiles 
  - \c GMRFLib_density_P() The cummulative density function
  
  additionally, we can free an object, using \c GMRFLib_free_density() and print a summary of the density using \c
  GMRFLib_density_printf().

  There are two functions which compare two densities, by computing the Kullback-Leibler divergence between them, it is
  - \c GMRFLib_kld() Compute the Kullback-Leibler divergence between two densities
  - \c GMRFLib_kld_sym() Compute the symmetric Kullback-Leibler divergence between two densities
  - \c GMRFLib_mkld() Compute the Kullback-Leibler divergence between two densities using only the first two moments
  - \c GMRFLib_mkld_sym() Compute the symmetric Kullback-Leibler divergence between two densities using only the first two moments
  
  Properties of the density is also available in the \c GMRFLib_density_tp -object, as
  - \c GMRFLib_density_tp::mean The mean in the standardised scale
  - \c GMRFLib_density_tp::stdev The standard deviation in the standardised scale
  - \c GMRFLib_density_tp::kld The Kullback-Leibler divergence from the standard normal
  - \c GMRFLib_density_tp::user_mean The mean in the users scale, ie \f$\mbox{E}(X)\f$ in the example above
  - \c GMRFLib_density_tp::user_stdev The standard deviation in the users scale, ie \f$\mbox{Stdev}(X)\f$ in the example above
  - \c GMRFLib_density_tp::std_mean The location parameter used for the standardisation
  - \c GMRFLib_density_tp::std_stdev The scale parameter used for the standardisation

  These functions convert from the the user scale to the internal standardised scale, and back,
  - \c GMRFLib_density_user2std()
  - \c GMRFLib_density_std2user()
  where the standarisation is given by the \c GMRFLib_density_tp object.
  
*/

#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

/* Pre-hg-Id: $Id: density.c,v 1.147 2009/11/04 18:24:30 hrue Exp $ */

#include <stddef.h>
#include <float.h>
#include <math.h>
#include <strings.h>
#include <string.h>
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

static double WEIGHT_PROB = 0.98;

#define CONST_1 0.6266570686577500604			       // sqrt(M_PI/8.0);
#define CONST_2 (-0.69314718055994528623)		       // log(0.5);

int GMRFLib_density_prune_weights(int *n_idx, int *idx, double *weights, int n)
{
	// make a list of the largest scaled weights so that the cummulative sum is at least WEIGHT_PROB

	int i, debug = 0;
	double w_sum = 0.0;
	double *ww = Calloc(n, double);

	memcpy(ww, weights, n * sizeof(double));
	for (i = 0, w_sum = 0.0; i < n; i++) {
		w_sum += ww[i];
	}
	w_sum = 1.0 / w_sum;
	for (i = 0; i < n; i++) {
		ww[i] *= w_sum;
		idx[i] = i;
	}
	GMRFLib_qsorts((void *) ww, (size_t) n, sizeof(double), (void *) idx, sizeof(int), NULL, NULL, GMRFLib_dcmp_r);
	for (i = 0, *n_idx = 0, w_sum = 0.0; i < n; i++) {
		w_sum += ww[i];
		(*n_idx)++;
		if (w_sum > WEIGHT_PROB)
			break;
	}

	if (debug) {
		w_sum = 0.0;
		for (i = 0; i < *n_idx; i++) {
			w_sum += ww[i];
			printf("i %1d idx %1d n_idx %1d n %1d ww %g w_sum %g\n", i, idx[i], *n_idx, n, ww[i], w_sum);
		}
	}
	Free(ww);

	return GMRFLib_SUCCESS;
}
int GMRFLib_sn_density(double *dens, double x, void *param)
{
	/*
	 * compute the sn-density 
	 */
	double ldens;

	GMRFLib_EWRAP0(GMRFLib_sn_logdensity(&ldens, x, param));
	*dens = exp(ldens);

	return GMRFLib_SUCCESS;
}
double GMRFLib_log_gsl_cdf_ugaussian_P(double z)
{
	/*
	 * compute log(gsl_cdf_ugaussian_P(z)) for large |z| as well 
	 */
	if (ABS(z) < 8.0) {
		if (0) {
			// faster option. see also the code with some doc in inla.c. I havn't yet merged these functions
			double val;
			if (z > 0.0) {
				// val = 0.5 + 0.5 * sqrt(1.0 - exp(-sqrt(M_PI / 8.0) * SQR(z)));
				val = 0.5 + 0.5 * sqrt(1.0 - exp(-CONST_1 * SQR(z)));
			} else {
				val = 1.0 - (0.5 + 0.5 * sqrt(1.0 - exp(-CONST_1 * SQR(z))));
			}
			return (log(val));
		} else {
			return log(gsl_cdf_ugaussian_P(z));
		}
	} else {
		if (z > 0) {
			if (z > 37.0) {
				return 0.0;
			} else {
				double t1, t3, t8, t13;

				t1 = sqrt(M_PI);
				t3 = sqrt(0.2e1);
				t8 = exp(z * z / 2.0);
				t13 = -0.1e1 / t1 * t3 / z / t8 / 0.2e1;

				return t13;
			}
		} else {
			return -0.5 * SQR(z);
		}
	}
	abort();

	return 0.0;
}
int GMRFLib_sn_logdensity(double *ldens, double x, void *param)
{
	/*
	 * compute the log sn-density 
	 */
	static const double log_norm_const_gaussian = -0.918938533204672741780329736407;	/* log(1.0/sqrt(2.0*M_PI)) */
	GMRFLib_sn_param_tp *p = (GMRFLib_sn_param_tp *) param;
	double z = (x - p->xi) / p->omega;

	*ldens = M_LN2 + log_norm_const_gaussian - 0.5 * SQR(z) + GMRFLib_log_gsl_cdf_ugaussian_P(p->alpha * z) - log(p->omega);

	return GMRFLib_SUCCESS;
}
double GMRFLib_sn_logdensity_diff_xi(double x, void *param)
{
	GMRFLib_sn_param_tp *p = (GMRFLib_sn_param_tp *) param;
	double xi = p->xi;
	double omega = p->omega;
	double alpha = p->alpha;

	double t1, t2, t4, t9, t16, t17, t19, t20, t23, t25, t37;

	t1 = sqrt(0.2e1);
	t2 = M_PI * t1;
	t4 = x - xi;
	t9 = gsl_sf_erf(0.1e1 / omega * t4 * alpha * t1 / 0.2e1);
	t16 = alpha * alpha;
	t17 = t4 * t4;
	t19 = omega * omega;
	t20 = 0.1e1 / t19;
	t23 = exp(-t20 * t17 * t16 / 0.2e1);
	t25 = sqrt(M_PI);
	t37 =
	    0.1e1 / (t9 + 1.0) / 0.3141592654e1 * t20 * t1 * (t9 * x * t2 + x * t2 - t9 * xi * t2 - xi * t2 -
							      0.2e1 * omega * t25 * alpha * t23) / 0.2e1;

	return t37;
}
double GMRFLib_sn_logdensity_diff_omega(double x, void *param)
{
	GMRFLib_sn_param_tp *p = (GMRFLib_sn_param_tp *) param;
	double xi = p->xi;
	double omega = p->omega;
	double alpha = p->alpha;

	double t1, t2, t3, t5, t10, t14, t20, t24, t25, t27, t31, t32, t33, t34, t41, t55;

	t1 = sqrt(0.2e1);
	t2 = M_PI * t1;
	t3 = x * x;
	t5 = x - xi;
	t10 = gsl_sf_erf(0.1e1 / omega * t5 * alpha * t1 / 0.2e1);
	t14 = xi * x;
	t20 = xi * xi;
	t24 = alpha * alpha;
	t25 = t5 * t5;
	t27 = omega * omega;
	t31 = exp(-0.1e1 / t27 * t25 * t24 / 0.2e1);
	t32 = alpha * t31;
	t33 = sqrt(M_PI);
	t34 = omega * t33;
	t41 = t27 * t1;
	t55 = 0.1e1 / (t10 + 1.0) / M_PI / t27 / omega * t1
	    * (t10 * t3 * t2 + t3 * t2 - 0.2e1 * t10 * t14 * t2 - 0.2e1 * t14 * t2 + t10 * t20 * t2
	       + t20 * t2 - 0.2e1 * x * t34 * t32 + 0.2e1 * xi * t34 * t32 - t10 * M_PI * t41 - M_PI * t41) / 0.2e1;

	return t55;
}
double GMRFLib_sn_logdensity_diff_alpha(double x, void *param)
{
	GMRFLib_sn_param_tp *p = (GMRFLib_sn_param_tp *) param;
	double xi = p->xi;
	double omega = p->omega;
	double alpha = p->alpha;

	double t1, t2, t3, t5, t9, t11, t13, t15, t21, t25;

	t1 = alpha * alpha;
	t2 = x - xi;
	t3 = t2 * t2;
	t5 = omega * omega;
	t9 = exp(-0.1e1 / t5 * t3 * t1 / 0.2e1);
	t11 = sqrt(0.2e1);
	t13 = sqrt(M_PI);
	t15 = 0.1e1 / omega;
	t21 = gsl_sf_erf(t15 * t2 * alpha * t11 / 0.2e1);
	t25 = 0.1e1 / (t21 + 1.0) * t15 / t13 * t11 * t2 * t9;

	return t25;
}
int GMRFLib_sn_moments(double *mean, double *stdev, double *skewness, GMRFLib_sn_param_tp * p)
{
	/*
	 * compute two first moments of the sn 
	 */
	if (p) {
		double delta = p->alpha / sqrt(1.0 + SQR(p->alpha));

		if (mean) {
			*mean = p->xi + p->omega * sqrt(2.0 / M_PI) * delta;
		}
		if (stdev) {
			*stdev = p->omega * sqrt(1.0 - 2.0 * SQR(delta) / M_PI);
		}
		if (skewness) {
			/*
			 * compute the skewness of the sn (https://en.wikipedia.org/wiki/Skew_normal_distribution)
			 */
			*skewness = (4.0 - M_PI) / 2.0 * pow(delta * sqrt(2.0 / M_PI), 3.0) / pow(1.0 - 2 * SQR(delta) / M_PI, 3.0 / 2.0);
		}
	}
	return GMRFLib_SUCCESS;
}
int GMRFLib_sn_fit(GMRFLib_sn_param_tp * param, double *fval, double *x, double *log_density, int n)
{
	/*
	 * fit a skew-normal, return the parameters in 'param' and 'error' in fval 
	 */

	return GMRFLib_sn_fit__intern((void *) param, fval, x, log_density, (size_t) n, (size_t) 4);
}
int GMRFLib_normal_fit(double *mean, double *variance, double *fval, double *x, double *log_density, int n)
{
	/*
	 * fit a normal, return the mean, variance and the 'error' fval, [if !NULL]
	 * 
	 * we use the same code as the skew-normal, but with alpha fixed to be zero. 
	 */

	int retval;
	GMRFLib_sn_param_tp param;

	retval = GMRFLib_sn_fit__intern((void *) &param, fval, x, log_density, (size_t) n, (size_t) 3);
	if (retval == GMRFLib_SUCCESS) {
		if (mean) {
			*mean = param.xi;
		}
		if (variance) {
			*variance = SQR(param.omega);
		}
	}
	return retval;
}
int GMRFLib_sn_fit__intern(void *param, double *fval, double *x, double *log_density, size_t n, size_t m)
{
	/*
	 * fit the skew-normal to a density using the `n' observations
	 * 
	 * (x[0], log_density[0]), ...., (x[n-1], log_density[n-1])
	 * 
	 * returning the solution in 'param'.
	 * 
	 * the density is assumed to not normalized.
	 * 
	 * make 'alpha' the last parameters so we can use the same code to fit the normal later 
	 */

#define MAXIT (500)					       /* maximum iterations */
#define print_state(iter, s) if (debug) { \
		printf ("iter: %3u x = %.16f %.16f %.16f %.16f " "|f(x)| = %g\n", \
			iter, gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1), \
			gsl_vector_get (s->x, 2), (m == 4 ? gsl_vector_get (s->x, 3) : 0.0), \
			gsl_blas_dnrm2 (s->f));				\
		printf ("         dx = %.16f %.16f %.16f %.16f " "|f(x)| = %g\n", \
			gsl_vector_get (s->dx, 0), gsl_vector_get (s->dx, 1), \
			gsl_vector_get (s->dx, 2), (m == 4 ? gsl_vector_get (s->dx, 3) : 0.0), \
			gsl_blas_dnrm2 (s->f));}

	gsl_multifit_fdfsolver *s = NULL;
	gsl_multifit_function_fdf f;
	gsl_vector *xx = NULL;

	GMRFLib_sn_fit_data_tp data;
	GMRFLib_sn_param_tp *p = (GMRFLib_sn_param_tp *) param;

	int iter = 0, status, debug = 0, i, imax;
	double eps = GMRFLib_eps(1. / 3.), *log_density_scaled;

	GMRFLib_EWRAP0_GSL_PTR(xx = gsl_vector_alloc(m));

	/*
	 * scale the log_density, so the maximum is 0. 
	 */
	imax = 0;					       /* imax = index with max log_density */
	for (i = 1; i < (int) n; i++) {
		if (log_density[i] > log_density[imax]) {
			imax = i;
		}
	}
	log_density_scaled = Calloc(n, double);

	for (i = 0; i < (int) n; i++) {
		log_density_scaled[i] = log_density[i] - log_density[imax];
	}

	data.n = n;
	data.m = m;					       /* need this information in f, df etc */
	data.x = x;
	data.log_density = log_density_scaled;

	/*
	 * set initial values 
	 */
	gsl_vector_set(xx, 0, x[imax]);			       /* xi */
	gsl_vector_set(xx, 1, 1.0);			       /* omega */
	gsl_vector_set(xx, 2, 0.9189385332046725);	       /* equal to 0.5*log(2.0 * M_PI) */
	if (m == 4) {
		gsl_vector_set(xx, 3, 0.0);		       /* alpha */
	}

	f.f = GMRFLib_sn_fit_f;
	f.df = GMRFLib_sn_fit_df;
	f.fdf = GMRFLib_sn_fit_fdf;
	f.n = n;
	f.p = m;
	f.params = (void *) &data;

	GMRFLib_EWRAP0_GSL_PTR(s = gsl_multifit_fdfsolver_alloc(gsl_multifit_fdfsolver_lmsder, n, m));
	GMRFLib_EWRAP0_GSL(gsl_multifit_fdfsolver_set(s, &f, xx));

	print_state(iter, s);
	do {
		iter++;
		status = gsl_multifit_fdfsolver_iterate(s);
		print_state(iter, s);
		if (status) {
			break;
		}
		status = gsl_multifit_test_delta(s->dx, s->x, eps, eps);
	}
	while (status == GSL_CONTINUE && iter < MAXIT);

	if (iter >= MAXIT) {
		GMRFLib_ERROR(GMRFLib_ESN);
	}

	/*
	 * write the solution back 
	 */
	p->xi = gsl_vector_get(s->x, 0);
	p->omega = gsl_vector_get(s->x, 1);
	p->alpha = (m == 4 ? gsl_vector_get(s->x, 3) : 0.0);

	if (fval) {
		*fval = gsl_blas_dnrm2(s->f) / sqrt((double) n);
	}

	gsl_multifit_fdfsolver_free(s);
	gsl_vector_free(xx);
	Free(log_density_scaled);

	return GMRFLib_SUCCESS;
#undef print_state
#undef MAXIT
}
int GMRFLib_sn_fit_f(const gsl_vector * param, void *data, gsl_vector * f)
{
	GMRFLib_sn_fit_data_tp *d = (GMRFLib_sn_fit_data_tp *) data;
	GMRFLib_sn_param_tp p;
	double constant, ldens;
	size_t i;

	p.xi = gsl_vector_get(param, 0);
	p.omega = gsl_vector_get(param, 1);
	constant = gsl_vector_get(param, 2);
	p.alpha = (d->m == 4 ? gsl_vector_get(param, 3) : 0.0);

	if (p.omega <= 0.0) {
		return GSL_EFAILED;
	}

	for (i = 0; i < d->n; i++) {
		GMRFLib_EWRAP0(GMRFLib_sn_logdensity(&ldens, d->x[i], (void *) &p));
		gsl_vector_set(f, i, (d->log_density[i] - (ldens + constant)) * exp(d->log_density[i] * GMRFLib_SN_WEIGHT_SCALING));
	}

	return GSL_SUCCESS;
}
int GMRFLib_sn_fit_df(const gsl_vector * param, void *data, gsl_matrix * J)
{
	GMRFLib_sn_fit_data_tp *d = (GMRFLib_sn_fit_data_tp *) data;
	GMRFLib_sn_param_tp p;
	size_t i;

	p.xi = gsl_vector_get(param, 0);
	p.omega = gsl_vector_get(param, 1);
	p.alpha = (d->m == 4 ? gsl_vector_get(param, 3) : 0.0);

	if (p.omega <= 0.0) {
		return GSL_EFAILED;
	}

	for (i = 0; i < d->n; i++) {
		/*
		 * compute d f_i /d param_j
		 * 
		 * where f_i = [ y_i - (sn_dens() + constant) ] * exp(y_i)
		 * 
		 * recall that y_i is scaled so that maximum is 0. 
		 */
		double factor = exp(d->log_density[i] * GMRFLib_SN_WEIGHT_SCALING);

		gsl_matrix_set(J, i, 0, -GMRFLib_sn_logdensity_diff_xi(d->x[i], (void *) &p) * factor);
		gsl_matrix_set(J, i, 1, -GMRFLib_sn_logdensity_diff_omega(d->x[i], (void *) &p) * factor);
		gsl_matrix_set(J, i, 2, -factor);

		if (d->m == 4) {
			gsl_matrix_set(J, i, 3, -GMRFLib_sn_logdensity_diff_alpha(d->x[i], (void *) &p) * factor);
		}
	}

	return GSL_SUCCESS;
}
int GMRFLib_sn_fit_fdf(const gsl_vector * param, void *data, gsl_vector * f, gsl_matrix * J)
{
	int retval;

	retval = GMRFLib_sn_fit_f(param, data, f);
	if (retval != GSL_SUCCESS) {
		return retval;
	}

	retval = GMRFLib_sn_fit_df(param, data, J);
	if (retval != GSL_SUCCESS) {
		return retval;
	}

	return GSL_SUCCESS;
}
int GMRFLib_init_density(GMRFLib_density_tp * density, int lookup_tables)
{
	/*
	 * initialize 'density': compute the mean, stdev and the norm_const (for the log spline fit) 
	 */
	int i, k, np = 2 * GMRFLib_faster_integration_np, npm = 2 * np, debug = 0;
	double result, error, eps = GMRFLib_eps(1. / 2.), tmp, low = 0.0, high = 0.0, xval, ldens_max = -FLT_MAX, *xpm =
	    NULL, *ldm = NULL, *xp = NULL, integral, w[2] = {
	4.0, 2.0}, dx = 0.0, m1, m2, m3, x0, x1, d0, d1;

	if (!density) {
		return GMRFLib_SUCCESS;
	}

	if (density->type == GMRFLib_DENSITY_TYPE_GAUSSIAN) {
		density->mean = density->mean_gaussian;
		density->stdev = density->stdev_gaussian;
		density->skewness = 0.0;
		density->x_min = -GMRFLib_DENSITY_INTEGRATION_LIMIT * density->stdev + density->mean;
		density->x_max = GMRFLib_DENSITY_INTEGRATION_LIMIT * density->stdev + density->mean;
	} else {
		if (density->type == GMRFLib_DENSITY_TYPE_SKEWNORMAL) {
			/*
			 * for the skew-normal we know the moments 
			 */
			GMRFLib_sn_moments(&(density->mean), &(density->stdev), &(density->skewness), density->sn_param);
			density->x_min = -GMRFLib_DENSITY_INTEGRATION_LIMIT * density->stdev + density->mean;
			density->x_max = GMRFLib_DENSITY_INTEGRATION_LIMIT * density->stdev + density->mean;
		} else {
			GMRFLib_ASSERT(density->type == GMRFLib_DENSITY_TYPE_SCGAUSSIAN, GMRFLib_ESNH);

			low = density->x_min;
			high = density->x_max;
			dx = (high - low) / (npm - 1.0);

			if (debug) {
				P(low);
				P(high);
				P(dx);
			}

			if (GMRFLib_faster_integration) {
				double ldmax, log_integral;

				xpm = Calloc(2 * npm, double);
				ldm = xpm + npm;
				for (xval = low, i = 0; i < npm; xval += dx, i++) {
					xpm[i] = xval;
				}
				density->log_norm_const = 0.0;
				GMRFLib_evaluate_nlogdensity(ldm, xpm, npm, density);

				if (debug) {
					for (i = 0; i < npm; i++) {
						printf("INIT: i %d xpm %g ldm %g\n", i, xpm[i], ldm[i]);
					}
				}

				ldmax = GMRFLib_max_value(ldm, npm, NULL);
				GMRFLib_adjust_vector(ldm, npm);	/* so its well-behaved... */

				integral = exp(ldm[0]) + exp(ldm[npm - 1]);
				for (i = 1, k = 0; i < npm - 1; i++, k = (k + 1) % 2) {
					integral += exp(ldm[i]) * w[k];
				}
				integral *= dx / 3.0;

				log_integral = log(integral);
				for (i = 0; i < npm; i++) {
					ldm[i] -= log_integral;
				}
				density->log_norm_const = log_integral + ldmax;

				d0 = exp(ldm[0]);
				d1 = exp(ldm[npm - 1]);
				x0 = xpm[0];
				x1 = xpm[npm - 1];
				m1 = x0 * d0 + x1 * d1;
				m2 = SQR(x0) * d0 + SQR(x1) * d1;
				m3 = gsl_pow_3(x0) * d0 + gsl_pow_3(x1) * d1;
				for (i = 1, k = 0; i < npm - 1; i++, k = (k + 1) % 2) {
					double d, x, x2, x3;

					d = exp(ldm[i]) * w[k];
					x = xpm[i];
					x2 = x * x;
					x3 = x * x2;

					m1 += x * d;
					m2 += x2 * d;
					m3 += x3 * d;
				}
				m1 *= dx / 3.0;
				m2 *= dx / 3.0;
				m3 *= dx / 3.0;

				if (debug) {
					P(m1);
					P(m2);
					P(m3);
				}

				density->mean = m1;
				density->stdev = sqrt(DMAX(0.0, m2 - SQR(m1)));
				density->skewness = (m3 - 3.0 * m1 * SQR(density->stdev) - gsl_pow_3(m1)) / gsl_pow_3(density->stdev);
			} else {
				double ldens;
				GMRFLib_density_properties_tp prop;
				gsl_function F;

				low = density->x_min;
				high = density->x_max;
				dx = (high - low) / (npm - 1.0);

				prop.density = density;
				F.function = GMRFLib_evaluate_density__intern;
				F.params = (void *) &prop;

				/*
				 * the __intern function *use* the norm_const, so first we need to set it temporary to 1.0 
				 */
				density->log_norm_const = 0.0;
				for (xval = low; xval < high; xval += (high - low) / 20.0) {
					GMRFLib_evaluate_logdensity(&ldens, xval, density);
					if (xval == low || ldens > ldens_max) {
						ldens_max = ldens;
					}
				}
				density->log_norm_const = ldens_max;
				GMRFLib_gsl_integration_wrapper(&F, low, high, eps, eps, &result, &error);
				density->log_norm_const = log(result) + ldens_max;

				/*
				 * ...and from now on, the density is normalised 
				 */
				F.function = GMRFLib_evaluate_density_power__intern;	/* use this function for f(x)*x^power */

				prop.power = 1;
				GMRFLib_gsl_integration_wrapper(&F, low, high, eps, eps, &result, &error);
				density->mean = result;

				prop.power = 2;
				GMRFLib_gsl_integration_wrapper(&F, low, high, eps, eps, &result, &error);
				tmp = result - SQR(density->mean);
				density->stdev = sqrt(DMAX(tmp, 0.0));


				/*
				 *  I do not care to add the skewness computations right now, as this code part is obsolete!
				 */
				fprintf(stderr, "\n\n\nTODO: added computation of skewness for a density, HERE!!!\n\n\n");
				abort();

			}
		}
	}

	/*
	 * for convenience, here the mean and the stdev in the users scale 
	 */
	density->user_mean = density->std_stdev * density->mean + density->std_mean;
	density->user_stdev = density->std_stdev * density->stdev;
	density->user_mode = NAN;			       /* yes, this is the value if its not computed */

	/*
	 * new style speedup 
	 */
	if (lookup_tables) {
		/*
		 * build fast lookup tables for P() and Pinv() calculations using linear interpolation 
		 */
		double *p = NULL, *work = NULL, *val = NULL, *dens;
		int some_nans = 0;

		if (!xpm) {
			/*
			 * in case we haven't done this before; same as above 
			 */
			low = density->x_min;
			high = density->x_max;
			dx = (high - low) / (npm - 1.0);

			xpm = Calloc(2 * npm, double);
			ldm = xpm + npm;
			for (xval = low, i = 0; i < npm; xval += dx, i++) {
				xpm[i] = xval;
			}
			GMRFLib_evaluate_nlogdensity(ldm, xpm, npm, density);
			GMRFLib_adjust_vector(ldm, npm);       /* so its well-behaved... */
		}

		/*
		 * find the mode fitting a quadratic around the best point, like a one-step Newton-Raphson, since we have so excellent initial value
		 */
		int imax;
		GMRFLib_max_value(ldm, npm - 1, &imax);
		density->user_mode = density->std_mean + density->std_stdev *
		    (xpm[imax] - (((ldm[imax + 1] - ldm[imax - 1]) / (2.0 * dx)) / ((ldm[imax + 1] - 2.0 * ldm[imax] + ldm[imax - 1]) / SQR(dx))));

		work = Calloc(4 * np, double);
		dens = work;
		val = work + np;
		p = work + 2 * np;
		xp = work + 3 * np;

		for (i = 0; i < np; i++) {
			xp[i] = xpm[2 * i];
			dens[i] = exp(ldm[2 * i]);
		}
		val[0] = 0.0;
		integral = val[0];
		for (i = 1; i < np; i++) {
			val[i] = (dens[i - 1] + dens[i] + 4.0 * exp(ldm[2 * i - 1])) * dx / 6.0;
			val[i] = DMAX(0.0, val[i]);
			integral += val[i];
		}

		p[0] = 0.0;
		for (i = 1; i < np - 1; i++) {
			p[i] = p[i - 1] + val[i] / integral;
			some_nans = some_nans || isnan(p[i]);
		}
		p[np - 1] = 1.0;

		if (!some_nans) {
			density->P = Calloc(1, GMRFLib_spline_tp);
			density->Pinv = Calloc(1, GMRFLib_spline_tp);
			GMRFLib_EWRAP0_GSL_PTR(density->P->accel = gsl_interp_accel_alloc());
			GMRFLib_EWRAP0_GSL_PTR(density->P->spline = gsl_spline_alloc(gsl_interp_linear, (unsigned int) np));
			GMRFLib_EWRAP0_GSL(gsl_spline_init(density->P->spline, xp, p, (unsigned int) np));

			GMRFLib_unique_additive2(&np, p, xp, GMRFLib_eps(1. / 2.));
			GMRFLib_EWRAP0_GSL_PTR(density->Pinv->accel = gsl_interp_accel_alloc());
			GMRFLib_EWRAP0_GSL_PTR(density->Pinv->spline = gsl_spline_alloc(gsl_interp_linear, (unsigned int) np));
			GMRFLib_EWRAP0_GSL(gsl_spline_init(density->Pinv->spline, p, xp, (unsigned int) np));
		} else {
			density->P = NULL;
		}
		Free(work);
	} else {
		density->P = NULL;
		density->Pinv = NULL;
	}

	Free(xpm);

	return GMRFLib_SUCCESS;
}

/*!
  \brief Compute the log of the density at location x

  This function compute the log of the density

  \param[out] logdens At output \c *logdens contains the log of the density at \f$x\f$
  \param[in] x The argument for which the density is evaluated. Note that \f$x\f$ is in standardised scale.
  \param[in] density A pointer to the density of type \c GMRFLib_density_tp
  
  \sa GMRFLib_evaluate_density()
*/
int GMRFLib_evaluate_logdensity(double *logdens, double x, GMRFLib_density_tp * density)
{
	return GMRFLib_evaluate_nlogdensity(logdens, &x, 1, density);
}
int GMRFLib_evaluate_nlogdensity(double *logdens, double *x, int n, GMRFLib_density_tp * density)
{
	/*
	 * evaluate the log-density-object. Note that x is in *standardised scale* . 
	 */
	static double log_norm_const_gaussian = -0.918938533204672741780329736407;	/* log(1.0/sqrt(2.0*M_PI)) */
	int i;

	switch (density->type) {
	case GMRFLib_DENSITY_TYPE_GAUSSIAN:

		for (i = 0; i < n; i++) {
			logdens[i] = log_norm_const_gaussian - log(density->stdev) - 0.5 * SQR(x[i] - density->mean) / SQR(density->stdev);
		}
		break;

	case GMRFLib_DENSITY_TYPE_SKEWNORMAL:
	{
		/*
		 * just inline GMRFLib_sn_logdensity() 
		 */

		GMRFLib_sn_param_tp *p = density->sn_param;
		double local_const_1 = M_LN2 + log_norm_const_gaussian - log(p->omega);
		double a = 1.0 / p->omega;
		double b = -p->xi / p->omega;
		double z, zz, val;

		for (i = 0; i < n; i++) {
			// inline the most important case of log(Phi(...)) and use a very good approximation
			z = a * x[i] + b;
			zz = p->alpha * z;
			if (ABS(zz) < 8.0) {
				if (zz > 0.0) {
					// val = log(0.5 + 0.5 * sqrt(1.0 - exp(- CONST_1 * SQR(zz))));
					val = CONST_2 + log(1.0 + sqrt(1.0 - exp(-CONST_1 * SQR(zz))));
				} else {
					// val = log(0.5 - 0.5 * sqrt(1.0 - exp(- CONST_1 * SQR(zz))));
					val = CONST_2 + log(1.0 - sqrt(1.0 - exp(-CONST_1 * SQR(zz))));
				}
			} else {
				// use the more complitated asympt expression, which we do here (for which the code in
				// the prev {} is a copy of
				val = GMRFLib_log_gsl_cdf_ugaussian_P(zz);
			}
			logdens[i] = local_const_1 - 0.5 * SQR(z) + val;
		}
		if (0) {
			// OLD code
			for (i = 0; i < n; i++) {
				z = (x[i] - p->xi) / p->omega;
				logdens[i] =
				    M_LN2 + log_norm_const_gaussian - 0.5 * SQR(z) + GMRFLib_log_gsl_cdf_ugaussian_P(p->alpha * z) - log(p->omega);
			}
		}
		break;
	}

	case GMRFLib_DENSITY_TYPE_SCGAUSSIAN:
	{
		for (i = 0; i < n; i++) {
			double xmax = density->log_correction->spline->interp->xmax, xmin = density->log_correction->spline->interp->xmin;

			if (x[i] >= xmin && x[i] <= xmax) {
				logdens[i] = gsl_spline_eval(density->log_correction->spline, x[i], density->log_correction->accel)
				    - 0.5 * SQR(x[i]) - density->log_norm_const;
			} else {
				double diff, cor, f0;
				if (x[i] >= xmax - DBL_EPSILON) {
					f0 = gsl_spline_eval(density->log_correction->spline, xmax, density->log_correction->accel);
					diff = gsl_spline_eval_deriv(density->log_correction->spline, xmax, density->log_correction->accel);
					cor = f0 + DMIN(0.0, diff) * (x[i] - xmax);
				} else {
					f0 = gsl_spline_eval(density->log_correction->spline, xmin, density->log_correction->accel);
					diff = gsl_spline_eval_deriv(density->log_correction->spline, xmin, density->log_correction->accel);
					cor = f0 + DMAX(0.0, diff) * (x[i] - xmin);
				}
				logdens[i] = cor - 0.5 * SQR(x[i]) - density->log_norm_const;
			}
		}
		break;
	}

	default:

		GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
		break;
	}

	return GMRFLib_SUCCESS;
}

/*!
  \brief Compute the density at location x

  This function compute the density

  \param[out] dens At output \c *dens contains the value of the density at \f$x\f$
  \param[in] x The argument for which the density is evaluated. Note that \f$x\f$ is in standardised scale.
  \param[in] density A pointer to the density of type \c GMRFLib_density_tp
  
  \sa GMRFLib_evaluate_logdensity()
*/
int GMRFLib_evaluate_density(double *dens, double x, GMRFLib_density_tp * density)
{
	double ldens;

	GMRFLib_evaluate_logdensity(&ldens, x, density);
	*dens = exp(ldens);

	return GMRFLib_SUCCESS;
}
int GMRFLib_evaluate_ndensity(double *dens, double *x, int n, GMRFLib_density_tp * density)
{
	int i;

	GMRFLib_evaluate_nlogdensity(dens, x, n, density);
	for (i = 0; i < n; i++) {
		dens[i] = exp(dens[i]);
	}
	return GMRFLib_SUCCESS;
}
double GMRFLib_evaluate_density__intern(double x, void *param)
{
	double dens = 0.0;
	GMRFLib_density_properties_tp *prop = (GMRFLib_density_properties_tp *) param;

	GMRFLib_evaluate_density(&dens, x, prop->density);

	return dens;
}
double GMRFLib_evaluate_density_power__intern(double x, void *param)
{
	double dens = 0.0;
	GMRFLib_density_properties_tp *prop = (GMRFLib_density_properties_tp *) param;

	GMRFLib_evaluate_density(&dens, x, prop->density);

	return dens * gsl_pow_int(x, prop->power);
}

/*!
  \brief Free a density-object

  This function free a density object

  \param[in] density A pointer to the density object to be free'd.
*/
int GMRFLib_free_density(GMRFLib_density_tp * density)
{
#define FreeSpline(_s) \
	if (_s) {						\
		if (_s->spline)					\
			gsl_spline_free(_s->spline);		\
		if (_s->accel)					\
			gsl_interp_accel_free(_s->accel);	\
		Free(_s);					\
	}


	if (density) {
		switch (density->type) {
		case GMRFLib_DENSITY_TYPE_GAUSSIAN:
			break;

		case GMRFLib_DENSITY_TYPE_SKEWNORMAL:
			break;

		case GMRFLib_DENSITY_TYPE_SCGAUSSIAN:
			FreeSpline(density->log_correction);
			break;

		default:
			GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
			break;
		}
		FreeSpline(density->P);
		FreeSpline(density->Pinv);
		Free(density);
	}
#undef FreeSpline
	return GMRFLib_SUCCESS;
}
double GMRFLib_evaluate_density_kld__intern(double x, void *param)
{
	static double norm_const_gaussian = 0.398942280401432677939946059934;	/* 1.0/sqrt(2.0*M_PI) */
	static double log_norm_const_gaussian = -0.918938533204672741780329736407;	/* log(1.0/sqrt(2.0*M_PI)) */
	GMRFLib_density_properties_tp *prop = (GMRFLib_density_properties_tp *) param;
	double ldens;

	GMRFLib_evaluate_logdensity(&ldens, x, prop->density);

	return norm_const_gaussian * exp(-0.5 * SQR(x)) * (log_norm_const_gaussian - 0.5 * SQR(x) - ldens);
}
double GMRFLib_density_Pinv_f(double x, void *param)
{
	GMRFLib_density_properties_tp *prop = (GMRFLib_density_properties_tp *) param;
	double px = NAN;

	GMRFLib_density_P(&px, x, prop->density);

	return px - prop->prob;
}
double GMRFLib_density_Pinv_df(double x, void *param)
{
	GMRFLib_density_properties_tp *prop = (GMRFLib_density_properties_tp *) param;
	double dens = 0.0;

	GMRFLib_evaluate_density(&dens, x, prop->density);

	return dens;
}
void GMRFLib_density_Pinv_fdf(double x, void *param, double *f, double *df)
{
	*f = GMRFLib_density_Pinv_f(x, param);
	*df = GMRFLib_density_Pinv_df(x, param);

	return;
}

/*!
  \brief Compute the quantiles of a density

  This function compute the quantiles for a density by solving the equation \f[ \mbox{Prob}( \tilde{X} \le x_p) = p\f] for
  \f$x_p\f$.  Note that \f$\tilde{X}\f$ is in standardised scale.

  \param[out] xp  At output, then \f$*xp\f$ contains the quantile
  \param[in] p The probability
  \param[in] density The density itself

  \note If the density is not Gaussian, then if \f$p=0\f$ or \f$p=1\f$ then the lower or higher limit for the numerical
  integration itself is returned. In this way, the quantiles are continous with respect to \f$p\f$.
  
  \sa GMRFLib_density_P()
*/
int GMRFLib_density_Pinv(double *xp, double p, GMRFLib_density_tp * density)
{
	/*
	 * solve x such that Prob (X <= xp) = p, for the `density'
	 * 
	 * NOTE that 'xp' is in standarized scale. 
	 */
	GMRFLib_ASSERT(p >= 0 && p <= 1, GMRFLib_EPARAMETER);

	/*
	 * if the density is Gaussian, then its easy 
	 */
	if (density->type == GMRFLib_DENSITY_TYPE_GAUSSIAN) {
		*xp = density->mean + density->stdev * gsl_cdf_ugaussian_Pinv(p);
	} else {
		GMRFLib_ASSERT(density->Pinv, GMRFLib_ESNH);
		if (density->Pinv->spline) {
			*xp = gsl_spline_eval(density->Pinv->spline, p, density->Pinv->accel);
		} else {
			GMRFLib_ASSERT(density->Pinv->spline != NULL, GMRFLib_ESNH);
		}
	}
	return GMRFLib_SUCCESS;
}

/*!
  \brief The cummulative distribution function

  This function computes \f[ \mbox{Prob}( \tilde{X} \le x) \f] for
  given \f$x\f$.  Note that \f$\tilde{X}\f$ is in standardised scale.

  \param[out] px  At output, the probability for \f$\tilde{X} \le x\f$.
  \param[in] x The argument
  \param[in] density The density itself

  \note If the density is not Gaussian: When \f$x\f$ is lower or higher than the limit for the numerical integration, then 0 or
  1 is returned. In this way, \c GMRFLib_density_P() is consistent with \c GMRFLib_density_Pinv()
  
  \sa GMRFLib_density_Pinv()
*/
int GMRFLib_density_P(double *px, double x, GMRFLib_density_tp * density)
{
	/*
	 * compute px = Prob (X <= x) for the `density'
	 * 
	 * NOTE that 'x' is in standarized scale. 
	 */
	double result;

	GMRFLib_ENTER_ROUTINE;

	/*
	 * if the density is Gaussian, then its easy 
	 */
	if (density->type == GMRFLib_DENSITY_TYPE_GAUSSIAN) {
		result = gsl_cdf_ugaussian_P((x - density->mean) / density->stdev);
	} else {
		if (density->P->spline) {
			if (x <= density->P->spline->interp->xmin) {
				result = 0.0;
			} else if (x >= density->P->spline->interp->xmax) {
				result = 1.0;
			} else {
				result = gsl_spline_eval(density->P->spline, x, density->P->accel);
			}
		} else {
			GMRFLib_ASSERT(density->P->spline != NULL, GMRFLib_ESNH);
		}
	}

	*px = result;
	GMRFLib_LEAVE_ROUTINE;

	return GMRFLib_SUCCESS;
}
int GMRFLib_evaluate_densities(double *dens, double x_user, int n, GMRFLib_density_tp ** densities, double *weights)
{
	/*
	 * evaluate the density in ***USER SCALE*** at `x_user', where the density is given as a mixture
	 * 
	 * \sum_{i=0}^{n-1} weights[i]*densities[i]
	 * 
	 * the weights need not to be scaled. 
	 */
	int i, j, *idx = NULL, n_idx;
	double w_sum = 0.0, d_tmp = 0.0, d = 0.0, x_std;

	idx = Calloc(n, int);
	GMRFLib_density_prune_weights(&n_idx, idx, weights, n);

	for (j = 0; j < n_idx; j++) {
		i = idx[j];
		x_std = GMRFLib_density_user2std(x_user, densities[i]);
		GMRFLib_evaluate_density(&d_tmp, x_std, densities[i]);
		d += weights[i] * d_tmp / densities[i]->std_stdev;
		w_sum += weights[i];
	}
	*dens = d / w_sum;
	Free(idx);

	return GMRFLib_SUCCESS;
}
int GMRFLib_evaluate_ndensities(double *dens, int nd, double *x_user, int nx, GMRFLib_density_tp ** densities, double *weights)
{
	/*
	 * evaluate the density in ***USER SCALE*** at `x_user[j]', j=0...nx, where the density is given as a mixture
	 * 
	 * \sum_{i=0}^{n-1} weights[i]*densities[i]
	 * 
	 * the weights need not to be scaled. 
	 */
	int i, j, k, n_idx, *idx = NULL, n_alloc = IMAX(nd, nx);
	double w_sum = 0.0, *d_tmp, *d = NULL, *x_std;

	d = Calloc(4 * n_alloc, double);
	d_tmp = &d[n_alloc];
	x_std = &d[2 * n_alloc];
	idx = (int *) &d[3 * n_alloc];

	GMRFLib_density_prune_weights(&n_idx, idx, weights, nd);

	for (k = 0; k < n_idx; k++) {
		i = idx[k];
		w_sum += weights[i];

		GMRFLib_density_user2std_n(x_std, x_user, densities[i], nx);
		if (0) {
			// Old code
			for (j = 0; j < nx; j++)
				x_std[j] = GMRFLib_density_user2std(x_user[j], densities[i]);
		}
		GMRFLib_evaluate_ndensity(d_tmp, x_std, nx, densities[i]);
		for (j = 0; j < nx; j++) {
			d[j] += weights[i] * d_tmp[j] / densities[i]->std_stdev;
		}
	}

	w_sum /= w_sum;
	for (j = 0; j < nx; j++) {
		dens[j] = d[j] * w_sum;
	}

	Free(d);
	return GMRFLib_SUCCESS;
}
int GMRFLib_evaluate_gdensities(double *dens, double x_user, int n, GMRFLib_density_tp ** densities, double *weights)
{
	/*
	 * evaluate the gaussian densities in ***USER SCALE*** at `x_user', where the density is given as a gaussian
	 * mixture. the mean and variances for each term in the mixture are given by the std_mean and std_stdev.
	 * 
	 * density = \sum_{i=0}^{n-1} weights[i]*gdensities[i] / C, where C = sum weights[i]
	 * 
	 * the weights need not to be scaled. 
	 */
	int i, j, *idx = NULL, n_idx;
	double w_sum = 0.0, d_tmp = 0.0, d = 0.0, x_std;

	idx = Calloc(n, int);
	GMRFLib_density_prune_weights(&n_idx, idx, weights, n);

	for (j = 0; j < n_idx; j++) {
		i = idx[j];
		x_std = GMRFLib_density_user2std(x_user, densities[i]);
		d_tmp = exp(-0.5 * SQR(x_std));
		d += weights[i] * d_tmp / densities[i]->std_stdev;
		w_sum += weights[i];
	}
	*dens = d / w_sum;
	Free(idx);

	return GMRFLib_SUCCESS;
}
const gsl_interp_type *GMRFLib_density_interp_type(int n)
{
	/*
	 * return the interpolation type depending on the number of points 
	 */

	if (n >= (int) gsl_interp_cspline->min_size) {	       /* n >= 3 */
		return gsl_interp_cspline;
	} else if (n >= (int) gsl_interp_linear->min_size) {   /* n >= 2 */
		return gsl_interp_linear;
	} else {
		GMRFLib_ASSERT_RETVAL(n >= (int) gsl_interp_linear->min_size, GMRFLib_EPARAMETER, NULL);
	}

	return NULL;
}
int GMRFLib_density_duplicate(GMRFLib_density_tp ** density_to, GMRFLib_density_tp * density_from)
{
	int n = 1;
	double weights = 1.0;

	GMRFLib_density_combine(density_to, NULL, n, &density_from, &weights);
	(*density_to)->flags = density_from->flags;

	return GMRFLib_SUCCESS;
}
int GMRFLib_density_combine(GMRFLib_density_tp ** density, GMRFLib_density_tp ** gdensity, int n, GMRFLib_density_tp ** densities, double *weights)
{
	/*
	 * make a new spline-corrected-gaussian density out of a weighted sum of densities and return this in DENSITY.  make a
	 * new spline-corrected-gaussian density out of a weighted sum of gaussian densities and return this in GDENSITY.
	 *
	 * both DENSITY and GDENSITY is optional, either of them can be NULL.
	 * 
	 * \sum_{i=0}^{n-1} weights[i]*densities[i]
	 * 
	 * the weights need not to be scaled. 
	 */

	int i, j, n_points = 30, np, np_g, np_max, nf, minp = 3;
	double mean, stdev, mean_g, stdev_g, *x_points = NULL, *x_points_g = NULL,
	    *log_dens = NULL, *log_dens_g = NULL, dens, x_real, m1, m2, sum_w, *ptr = NULL, m, sd, xx,
		f[] = { 0, 0.1, -0.1, 0.25, -0.25, 0.5, -0.5, 0.75, -0.75, 1.0, -1.0, 1.5, -1.5, 2.0, -2.0, 3.0, -3.0 };

	GMRFLib_ENTER_ROUTINE;
	nf = sizeof(f) / sizeof(double);
	if (n == 0) {
		if (density) {
			*density = NULL;
		}
		if (gdensity) {
			*gdensity = NULL;
		}
		GMRFLib_LEAVE_ROUTINE;
		return GMRFLib_SUCCESS;
	}

	/*
	 * compute the mean and variance in the user-scale 
	 */
	m1 = m2 = sum_w = 0.0;
	for (i = 0; i < n; i++) {
		m1 += weights[i] * densities[i]->user_mean;
		m2 += weights[i] * (SQR(densities[i]->user_stdev) + SQR(densities[i]->user_mean));
		sum_w += weights[i];
	}
	mean = m1 / sum_w;
	stdev = sqrt(DMAX(0.0, m2 / sum_w - SQR(mean)));

	m1 = m2 = sum_w = 0.0;
	for (i = 0; i < n; i++) {
		m1 += weights[i] * densities[i]->std_mean;
		m2 += weights[i] * (SQR(densities[i]->std_stdev) + SQR(densities[i]->std_mean));
		sum_w += weights[i];
	}
	mean_g = m1 / sum_w;
	stdev_g = sqrt(DMAX(0.0, m2 / sum_w - SQR(mean_g)));

	if (1) {
		// new code. only use the mean/stdev to layout points
		// FIXME1("COMBINE USING NEW CODE");
		np_max = n_points + nf;
		x_points = Calloc(np_max, double);
		x_points_g = (gdensity ? Calloc(np_max, double) : NULL);

		GMRFLib_ghq_abscissas(&ptr, n_points);
		memcpy(x_points, ptr, n_points * sizeof(double));
		memcpy(x_points + n_points, f, nf * sizeof(double));
		np = np_max;
		if (gdensity) {
			memcpy(x_points_g, ptr, n_points * sizeof(double));
			memcpy(x_points_g + n_points, f, nf * sizeof(double));
			np_g = np_max;
		} else {
			np_g = 0;
		}

		/*
		 * sort and remove ties or points to close. the _additive option is more 'pratical', whereas the _relative option is
		 * very conservative. We need to ensure that we are not ending up with fewer than minp points which is minimum for the
		 * spline-interpolant (or 5 for the akima-interpolant).
		 */
		qsort(x_points, (size_t) np, sizeof(double), GMRFLib_dcmp);
		if (gdensity) {
			qsort(x_points_g, (size_t) np_g, sizeof(double), GMRFLib_dcmp);
		}

		double *x_points_tmp = NULL;
		int np_tmp;

		x_points_tmp = Calloc(np, double);
		np_tmp = np;
		memcpy(x_points_tmp, x_points, np * sizeof(double));
		GMRFLib_unique_additive(&np_tmp, x_points_tmp, GMRFLib_eps(1. / 4.0));
		if (np_tmp >= minp) {			       /* then its ok */
			np = np_tmp;
			memcpy(x_points, x_points_tmp, np * sizeof(double));
		} else {
			GMRFLib_unique_relative(&np_tmp, x_points_tmp, GMRFLib_eps(1. / 3.0));
			if (np_tmp >= minp) {		       /* then its ok */
				np = np_tmp;
				memcpy(x_points, x_points_tmp, np * sizeof(double));
			}
		}
		Free(x_points_tmp);

		if (gdensity) {
			x_points_tmp = Calloc(np_g, double);
			np_tmp = np_g;
			memcpy(x_points_tmp, x_points_g, np_g * sizeof(double));
			GMRFLib_unique_additive(&np_tmp, x_points_tmp, GMRFLib_eps(1. / 4.0));
			if (np_tmp >= minp) {		       /* then its ok */
				np_g = np_tmp;
				memcpy(x_points_g, x_points_tmp, np_g * sizeof(double));
			} else {
				GMRFLib_unique_relative(&np_tmp, x_points_tmp, GMRFLib_eps(1. / 3.0));
				if (np_tmp >= minp) {	       /* then its ok */
					np_g = np_tmp;
					memcpy(x_points_g, x_points_tmp, np_g * sizeof(double));
				}
			}
			Free(x_points_tmp);
		}
	} else {
		// old code

		/*
		 * chose abscissas as for the GHQ + the mean +/- 0, 0.25, 0.5 and 1.0*stdev, for each density in the mixture 
		 */
		nf = sizeof(f) / sizeof(double);
		np_max = n_points + nf * n;		       /* maximum value of np */
		x_points = Calloc(np_max, double);
		x_points_g = (gdensity ? Calloc(np_max, double) : NULL);

		GMRFLib_ghq_abscissas(&ptr, n_points);
		memcpy(x_points, ptr, n_points * sizeof(double));
		if (gdensity) {
			memcpy(x_points_g, ptr, n_points * sizeof(double));
		}

		/*
		 * compute the new points. it is required that the new points are within the INTEGRATION_LIMIT, otherwise its density
		 * is essentially zero. this loop also determine `np' 
		 */
		np = n_points;
		for (i = 0; i < n; i++) {
			m = densities[i]->user_mean;
			sd = densities[i]->user_stdev;
			for (j = 0; j < nf; j++) {
				xx = (m + f[j] * sd - mean) / stdev;
				if (ABS(xx) < GMRFLib_DENSITY_INTEGRATION_LIMIT) {
					x_points[np++] = xx;
				}
			}
		}
		if (gdensity) {
			np_g = n_points;
			for (i = 0; i < n; i++) {
				m = densities[i]->std_mean;
				sd = densities[i]->std_stdev;
				for (j = 0; j < nf; j++) {
					xx = (m + f[j] * sd - mean_g) / stdev_g;
					if (ABS(xx) < GMRFLib_DENSITY_INTEGRATION_LIMIT) {
						x_points_g[np_g++] = xx;
					}
				}
			}
		} else {
			np_g = 0;
		}

		/*
		 * sort and remove ties or points to close. the _additive option is more 'pratical', whereas the _relative option is
		 * very conservative. We need to ensure that we are not ending up with fewer than minp points which is minimum for the
		 * spline-interpolant (or 5 for the akima-interpolant).
		 */
		qsort(x_points, (size_t) np, sizeof(double), GMRFLib_dcmp);
		if (gdensity) {
			qsort(x_points_g, (size_t) np_g, sizeof(double), GMRFLib_dcmp);
		}

		double *x_points_tmp = NULL;
		int np_tmp;

		x_points_tmp = Calloc(np, double);
		np_tmp = np;
		memcpy(x_points_tmp, x_points, np * sizeof(double));
		GMRFLib_unique_additive(&np_tmp, x_points_tmp, GMRFLib_eps(1. / 4.0));
		if (np_tmp >= minp) {			       /* then its ok */
			np = np_tmp;
			memcpy(x_points, x_points_tmp, np * sizeof(double));
		} else {
			GMRFLib_unique_relative(&np_tmp, x_points_tmp, GMRFLib_eps(1. / 3.0));
			if (np_tmp >= minp) {		       /* then its ok */
				np = np_tmp;
				memcpy(x_points, x_points_tmp, np * sizeof(double));
			}
		}
		Free(x_points_tmp);

		if (gdensity) {
			x_points_tmp = Calloc(np_g, double);
			np_tmp = np_g;
			memcpy(x_points_tmp, x_points_g, np_g * sizeof(double));
			GMRFLib_unique_additive(&np_tmp, x_points_tmp, GMRFLib_eps(1. / 4.0));
			if (np_tmp >= minp) {		       /* then its ok */
				np_g = np_tmp;
				memcpy(x_points_g, x_points_tmp, np_g * sizeof(double));
			} else {
				GMRFLib_unique_relative(&np_tmp, x_points_tmp, GMRFLib_eps(1. / 3.0));
				if (np_tmp >= minp) {	       /* then its ok */
					np_g = np_tmp;
					memcpy(x_points_g, x_points_tmp, np_g * sizeof(double));
				}
			}
			Free(x_points_tmp);
		}
	}

	log_dens = Calloc(np, double);
	log_dens_g = (gdensity ? Calloc(np_g, double) : NULL);

	/*
	 * compute the weighted density. note that we have to go through the user/real-scale to get this right 
	 */
	if (1) {
		/*
		 * new improved code; use inline 
		 */
		double *xx_real = NULL, *ddens = NULL;

		xx_real = Calloc(2 * np, double);
		ddens = &xx_real[np];
		for (i = 0; i < np; i++) {
			xx_real[i] = x_points[i] * stdev + mean;
		}
		GMRFLib_evaluate_ndensities(ddens, n, xx_real, np, densities, weights);
		for (i = 0; i < np; i++) {
			log_dens[i] = (ddens[i] > 0.0 ? log(ddens[i]) : -FLT_MAX);
		}
		Free(xx_real);
	} else {
		/*
		 * old version 
		 */
		for (i = 0; i < np; i++) {
			x_real = x_points[i] * stdev + mean;
			GMRFLib_evaluate_densities(&dens, x_real, n, densities, weights);
			log_dens[i] = (dens > 0.0 ? log(dens) : -FLT_MAX);
		}
	}
	GMRFLib_adjust_vector(log_dens, np);

	if (gdensity) {
		for (i = 0; i < np_g; i++) {
			x_real = x_points_g[i] * stdev_g + mean_g;
			GMRFLib_evaluate_gdensities(&dens, x_real, n, densities, weights);
			log_dens_g[i] = (dens > 0.0 ? log(dens) : -FLT_MAX);
		}
		GMRFLib_adjust_vector(log_dens_g, np_g);
	}

	if (density) {
		GMRFLib_density_create(density, GMRFLib_DENSITY_TYPE_SCGAUSSIAN, np, x_points, log_dens, mean, stdev, GMRFLib_TRUE);
	}
	if (gdensity) {
		GMRFLib_density_create(gdensity, GMRFLib_DENSITY_TYPE_SCGAUSSIAN, np_g, x_points_g, log_dens_g, mean_g, stdev_g, GMRFLib_TRUE);
	}

	Free(x_points);
	Free(x_points_g);
	Free(log_dens);
	Free(log_dens_g);

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}
int GMRFLib_density_create_normal(GMRFLib_density_tp ** density, double mean, double stdev, double std_mean, double std_stdev)
{
	/*
	 * create a normal density,
	 * 
	 * (X - std_mean)/std_stdev ~ N (mean, stdev^2)
	 * 
	 */

	(*density) = Calloc(1, GMRFLib_density_tp);
	(*density)->type = GMRFLib_DENSITY_TYPE_GAUSSIAN;
	(*density)->std_mean = std_mean;
	(*density)->std_stdev = std_stdev;
	(*density)->mean_gaussian = mean;
	(*density)->stdev_gaussian = stdev;

	GMRFLib_EWRAP0(GMRFLib_init_density(*density, GMRFLib_TRUE));

	return GMRFLib_SUCCESS;
}
int GMRFLib_density_create_sn(GMRFLib_density_tp ** density, GMRFLib_sn_param_tp sn_param, double std_mean, double std_stdev, int lookup_tables)
{
	/*
	 * create a skew-normal density,
	 * 
	 * (X - std_mean)/std_stdev ~ SN(param) 
	 */

	(*density) = Calloc(1, GMRFLib_density_tp);
	(*density)->type = GMRFLib_DENSITY_TYPE_SKEWNORMAL;
	(*density)->std_mean = std_mean;
	(*density)->std_stdev = std_stdev;
	(*density)->sn_param = Calloc(1, GMRFLib_sn_param_tp);
	memcpy((void *) (*density)->sn_param, (const void *) &sn_param, sizeof(GMRFLib_sn_param_tp));

	GMRFLib_EWRAP0(GMRFLib_init_density(*density, lookup_tables));

	return GMRFLib_SUCCESS;
}
int GMRFLib_density_adjust_vector(double *ldens, int n)
{
	int i;
	double maxdev = 2.0 * log(DBL_EPSILON);

	GMRFLib_adjust_vector(ldens, n);
	for (i = 0; i < n; i++) {
		ldens[i] = DMAX(maxdev, ldens[i]);
	}
	return GMRFLib_SUCCESS;
}
int GMRFLib_density_create(GMRFLib_density_tp ** density, int type, int n, double *x, double *logdens, double std_mean,
			   double std_stdev, int lookup_tables)
{
	/*
	 * create a density of type 'type', so that
	 * 
	 * (X - std_mean)/std_stdev
	 * 
	 * is either Normal, Skew-normal or a spline-corrected standard Gaussian.
	 * 
	 * given n points (x_i, logdens_i), for i=0...n-1, of the density of (X - std_mean)/std_stdev
	 *
	 * make lookup_tables if LOOKUP_TABLES is TRUE
	 */
	int i, debug = 0;
	double *xx = NULL, *ldens = NULL, g_mean, g_var;
	GMRFLib_sn_param_tp sn_param;

	xx = Calloc(n, double);
	ldens = Calloc(n, double);
	memcpy(xx, x, (size_t) n * sizeof(double));
	memcpy(ldens, logdens, (size_t) n * sizeof(double));

	/*
	 * sort xx and remove ties
	 */
	GMRFLib_qsorts(xx, (size_t) n, sizeof(double), ldens, sizeof(double), NULL, 0, GMRFLib_dcmp);
	GMRFLib_unique_relative2(&n, xx, ldens, GMRFLib_eps(1. / 3.0));
	GMRFLib_adjust_vector(ldens, n);

	if (debug) {
		int ii;

		printf("%s: Create density n %d type %d\n", __GMRFLib_FuncName, n, type);
		for (ii = 0; ii < n; ii++) {
			printf("\tx %f ldens %f\n", xx[ii], ldens[ii]);
		}
	}

	/*
	 * special option 
	 */
	if (n == 1) {
		GMRFLib_EWRAP0(GMRFLib_density_create_normal(density, 0.0, 1.0, std_mean, std_stdev));
	} else {
		switch (type) {
		case GMRFLib_DENSITY_TYPE_GAUSSIAN:
			/*
			 * fit a gaussian 
			 */
			GMRFLib_EWRAP0(GMRFLib_normal_fit(&g_mean, &g_var, NULL, xx, ldens, n));
			GMRFLib_EWRAP0(GMRFLib_density_create_normal(density, g_mean, sqrt(g_var), std_mean, std_stdev));
			break;

		case GMRFLib_DENSITY_TYPE_SKEWNORMAL:
			/*
			 * fit skew-normal 
			 */
			GMRFLib_EWRAP0(GMRFLib_sn_fit(&sn_param, NULL, xx, ldens, n));
			GMRFLib_EWRAP0(GMRFLib_density_create_sn(density, sn_param, std_mean, std_stdev, lookup_tables));
			break;

		case GMRFLib_DENSITY_TYPE_SCGAUSSIAN:
			/*
			 * fit spline-corrected gaussian 
			 */
			(*density) = Calloc(1, GMRFLib_density_tp);
			(*density)->type = GMRFLib_DENSITY_TYPE_SCGAUSSIAN;
			(*density)->std_mean = std_mean;
			(*density)->std_stdev = std_stdev;
			(*density)->x_min = GMRFLib_min_value(xx, n, NULL);
			(*density)->x_max = GMRFLib_max_value(xx, n, NULL);

			GMRFLib_density_adjust_vector(ldens, n);	/* prevent extreme cases for the spline */
			for (i = 0; i < n; i++) {
				ldens[i] += 0.5 * SQR(xx[i]);  /* ldens is now the correction */
			}

			(*density)->log_correction = Calloc(1, GMRFLib_spline_tp);
			GMRFLib_EWRAP0_GSL_PTR((*density)->log_correction->accel = gsl_interp_accel_alloc());
			GMRFLib_EWRAP0_GSL_PTR((*density)->log_correction->spline =
					       gsl_spline_alloc(GMRFLib_density_interp_type(n), (unsigned int) n));
			GMRFLib_EWRAP0_GSL(gsl_spline_init((*density)->log_correction->spline, xx, ldens, (unsigned int) n));
			GMRFLib_EWRAP0(GMRFLib_init_density(*density, lookup_tables));
			/*
			 * to be sure, we reset them here
			 */
			(*density)->x_min = (*density)->log_correction->spline->interp->xmin;
			(*density)->x_max = (*density)->log_correction->spline->interp->xmax;

			break;

		default:
			GMRFLib_ASSERT(1 == 0, GMRFLib_ESNH);
		}
	}

	Free(xx);
	Free(ldens);

	return GMRFLib_SUCCESS;
}
int GMRFLib_density_new_mean(GMRFLib_density_tp ** new_density, GMRFLib_density_tp * density, double new_mean)
{
	/*
	 * return a new density, which the density with the given new mean 
	 */

#define N (30)
#define M (4)
	int i, n = N + 2 * M;
	double *x, *ld, alpha, eps[M] = { 1e-6, 1e-5, 1e-4, 1e-3 };

	x = Calloc(n, double);
	for (i = 0; i < M; i++) {
		GMRFLib_density_Pinv(&x[i], eps[i], density);
		GMRFLib_density_Pinv(&x[M + i], 1.0 - eps[i], density);
	}
	for (i = 0; i < N; i++) {
		alpha = 0.01 + 0.98 * (1.0 / (double) N) * i;
		GMRFLib_density_Pinv(&x[2 * M + i], alpha, density);
	}
	ld = Calloc(n, double);

	GMRFLib_evaluate_nlogdensity(ld, x, n, density);
	GMRFLib_density_create(new_density, GMRFLib_DENSITY_TYPE_SCGAUSSIAN, n, x, ld, new_mean, density->std_stdev, GMRFLib_TRUE);

	Free(x);
	Free(ld);

	return GMRFLib_SUCCESS;
}

/*!
  \brief Print a summary of a density

  \param[in] fp The file-pointer to write to. If \c NULL, then \c stdout is used.
  \param[in] density The density itself
*/
int GMRFLib_density_printf(FILE * fp, GMRFLib_density_tp * density)
{
	/*
	 * print summary statistics for a density 
	 */

	if (density) {
		if (!fp) {
			fp = stdout;
		}
		fprintf(fp, "%-35s %16.10f\n", "Mean", density->mean);
		fprintf(fp, "%-35s %16.10f\n", "Stdev", density->stdev);
		fprintf(fp, "%-35s %16.10f\n", "User mean", density->user_mean);
		fprintf(fp, "%-35s %16.10f\n", "User stdev", density->user_stdev);
		fprintf(fp, "%-35s %16.10f\n", "User mode", density->user_mode);
		fprintf(fp, "%-35s %16.10f\n", "Standarisation: mean", density->std_mean);
		fprintf(fp, "%-35s %16.10f\n", "Standarisation: stdev", density->std_stdev);
		switch (density->type) {
		case GMRFLib_DENSITY_TYPE_GAUSSIAN:
			fprintf(fp, "%-35s %-20s\n", "Density type", "Gaussian");
			fprintf(fp, "     %-30s %16.10f\n", "mean", density->mean_gaussian);
			fprintf(fp, "     %-30s %16.10f\n", "stdev", density->stdev_gaussian);
			break;
		case GMRFLib_DENSITY_TYPE_SKEWNORMAL:
			fprintf(fp, "%-35s %-20s\n", "Density type", "Skew Normal");
			fprintf(fp, "\t%-35s %16.10f\n", "xi", density->sn_param->xi);
			fprintf(fp, "     %-30s %16.10f\n", "omega", density->sn_param->omega);
			fprintf(fp, "     %-30s %16.10f\n", "alpha", density->sn_param->alpha);
			break;
		case GMRFLib_DENSITY_TYPE_SCGAUSSIAN:
			fprintf(fp, "%-26s %-30s\n", "Density type", "Spline corrected Gaussian");
			fprintf(fp, "     %-30s %16.10f\n", "Log normalisation constant", density->log_norm_const);
			fprintf(fp, "     %-30s %16.10f\n", "x_min", density->x_min);
			fprintf(fp, "     %-30s %16.10f\n", "x_max", density->x_max);
			break;
		}

		int i;
		for (i = 0; i < 8; i++) {
			fprintf(fp, "Flag %1d = %1d\n", i, GMRFLib_getbit(density->flags, i));
		}

		fflush(fp);
	}
	return GMRFLib_SUCCESS;
}
double GMRFLib_evaluate_density_kld2__intern(double x, void *param)
{
	/*
	 * x is in user-scale. compute the log(f(x)/g(x))f(x) in user-scale 
	 */
	double x_std0, x_std1, ld0, ld1, ld0u, ld1u;
	GMRFLib_density_tp **d = (GMRFLib_density_tp **) param;

	x_std0 = GMRFLib_density_user2std(x, d[0]);
	x_std1 = GMRFLib_density_user2std(x, d[1]);

	GMRFLib_evaluate_logdensity(&ld0, x_std0, d[0]);
	GMRFLib_evaluate_logdensity(&ld1, x_std1, d[1]);
	ld0u = ld0 - log(d[0]->std_stdev);
	ld1u = ld1 - log(d[1]->std_stdev);

	return (ld0u - ld1u) * exp(ld0u);
}

/*!
  \brief Compute the Kullback-Leibler divergence between two densities

  This function computes the Kullback-Leibler divergence from \c density (\f$f(x)\f$) to \c ddensity (\f$g(x)\f$), as \f[ \int
  \log\left(\frac{f(x)}{g(x)} \right) f(x) dx \f]

  \param[out] kld At output, then \c *kld contains the Kullback-Leibler divergence
  \param[in] density The reference density, \f$f(x)\f$
  \param[in] ddensity The other density density, \f$g(x)\f$

  \sa GMRFLib_kld_sym(), GMRFLib_mkld()
*/
int GMRFLib_kld(double *kld, GMRFLib_density_tp * density, GMRFLib_density_tp * ddensity)
{
	/*
	 * compute the KLD from density to ddensity, \int log(density/ddensity)*density dx 
	 */
	double result, error, eps = GMRFLib_eps(1. / 2.), low0, low1, high0, high1, low, high;
	gsl_function F;
	GMRFLib_density_tp *d[2];

	GMRFLib_ENTER_ROUTINE;

	low0 = -GMRFLib_DENSITY_INTEGRATION_LIMIT * density->user_stdev + density->user_mean;
	high0 = GMRFLib_DENSITY_INTEGRATION_LIMIT * density->user_stdev + density->user_mean;
	low1 = -GMRFLib_DENSITY_INTEGRATION_LIMIT * ddensity->user_stdev + ddensity->user_mean;
	high1 = GMRFLib_DENSITY_INTEGRATION_LIMIT * ddensity->user_stdev + ddensity->user_mean;

	low = DMIN(low0, low1);
	high = DMAX(high0, high1);

	d[0] = density;
	d[1] = ddensity;

	F.function = GMRFLib_evaluate_density_kld2__intern;
	F.params = (void *) d;

	GMRFLib_gsl_integration_wrapper(&F, low, high, eps, eps, &result, &error);
	*kld = DMAX(0.0, result);			       /* known to be positive */

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

/*!
  \brief Compute the Kullback-Leibler divergence between two densities using only the first two moments

  This function computes the Kullback-Leibler divergence from \c density (\f$f(x)\f$) to \c ddensity (\f$g(x)\f$), as \f[ \int
  \log\left(\frac{f(x)}{g(x)} \right) f(x) dx \f], approximating both densities with a Gaussian fitting the two first moments.

  \param[out] mkld At output, then \c *kld contains the Kullback-Leibler divergence
  \param[in] density The reference density, \f$f(x)\f$
  \param[in] ddensity The other density density, \f$g(x)\f$

  \sa GMRFLib_kld(), GMRFLib_mkld_sym()
*/
int GMRFLib_mkld(double *mkld, GMRFLib_density_tp * density, GMRFLib_density_tp * ddensity)
{
	/*
	 * Maple: simplify(int(phi((x-m1)/s1)*(-(1/2)*(x-m1)^2/s1^2-log(s1)+(1/2)*log(s2)+(1/2)*(x-m2)^2/s2^2)/s1, x = -infinity .. 
	 * infinity)) 
	 */
	if (density && ddensity) {
		double m1 = density->user_mean, m2 = ddensity->user_mean, cg1 = density->user_stdev, cg2 = ddensity->user_stdev;

		*mkld = -(0.2e1 * log(cg1) * pow(cg2, 0.2e1) - 0.2e1 * log(cg2) * pow(cg2, 0.2e1) - (double) (m2 * m2)
			  + (double) (2 * m1 * m2) - (double) (m1 * m1)) * pow(cg2, -0.2e1) / 0.2e1
		    + (-pow(cg2, 0.2e1) + pow(cg1, 0.2e1)) * pow(cg2, -0.2e1) / 0.2e1;
	} else {
		*mkld = 0.0;
	}

	return GMRFLib_SUCCESS;
}

/*!
  \brief Compute the symmetric Kullback-Leibler divergence between two densities

  This function computes the symmetric Kullback-Leibler divergence between \c density (\f$f(x)\f$) and \c ddensity (\f$g(x)\f$),
  as \f[\frac{1}{2}\int \log\left(\frac{f(x)}{g(x)} \right) f(x) dx +\frac{1}{2}\int \log\left(\frac{g(x)}{f(x)} \right) g(x) dx
  \f]

  \param[out] kld_sym At output, then \c *kld_sym contains the symmetric Kullback-Leibler divergence
  \param[in] density First density, \f$f(x)\f$
  \param[in] ddensity Second density, \f$g(x)\f$

  \sa GMRFLib_kld()
*/
int GMRFLib_kld_sym(double *kld_sym, GMRFLib_density_tp * density, GMRFLib_density_tp * ddensity)
{
	/*
	 * compute the symmetric KLD distance between density and ddensity 
	 */
	double kld0, kld1;

	GMRFLib_kld(&kld0, density, ddensity);
	GMRFLib_kld(&kld1, ddensity, density);
	*kld_sym = DMAX(0.0, (kld0 + kld1) / 2.0);

	return GMRFLib_SUCCESS;
}

/*!
  \brief Compute the symmetric Kullback-Leibler divergence between two densities using only the first two moments

  This function computes the symmetric Kullback-Leibler divergence between \c density (\f$f(x)\f$) and \c ddensity (\f$g(x)\f$),
  as \f[\frac{1}{2}\int \log\left(\frac{f(x)}{g(x)} \right) f(x) dx +\frac{1}{2}\int \log\left(\frac{g(x)}{f(x)} \right) g(x) dx
  \f] using only the first two moments

  \param[out] mkld_sym At output, then \c *kld_sym contains the symmetric Kullback-Leibler divergence
  \param[in] density First density, \f$f(x)\f$
  \param[in] ddensity Second density, \f$g(x)\f$

  \sa GMRFLib_mkld(), GMRFLib_kld_sym()
*/
int GMRFLib_mkld_sym(double *mkld_sym, GMRFLib_density_tp * density, GMRFLib_density_tp * ddensity)
{
	double mkld0, mkld1;

	GMRFLib_mkld(&mkld0, density, ddensity);
	GMRFLib_mkld(&mkld1, ddensity, density);
	*mkld_sym = DMAX(0.0, (mkld0 + mkld1) / 2.0);

	return GMRFLib_SUCCESS;
}

/*!
  \brief Convert a value \c x in the standarised scale to the user scale

  \param[in] x The value in the standarised scale as defined in \c density
  \param[in] density The density which defines the standarisation

  The value in the user scale is returned through the function itself.

  \sa GMRFLib_density_user2std()
*/
double GMRFLib_density_std2user(double x, GMRFLib_density_tp * density)
{
	return density->std_mean + x * density->std_stdev;
}
double GMRFLib_density_std2user_n(double *x_user, double *x, int n, GMRFLib_density_tp * density)
{
	int i;

	for (i = 0; i < n; i++) {
		x_user[i] = density->std_mean + x[i] * density->std_stdev;
	}
	return GMRFLib_SUCCESS;
}

/*!
  \brief Convert a value \c x in the user scale to the standarised scale as defined in \c density

  \param[in] x The value in the user scale
  \param[in] density The density which defines the standarisation

  The value in the standarised scale is returned through the function itself.s

  \sa GMRFLib_density_user2std()
*/
double GMRFLib_density_user2std(double x, GMRFLib_density_tp * density)
{
	return (x - density->std_mean) / density->std_stdev;
}
int GMRFLib_density_user2std_n(double *x_std, double *x, GMRFLib_density_tp * density, int n)
{
	// the vectorised version
	int i;
	double a = 1.0 / density->std_stdev;
	double b = -density->std_mean / density->std_stdev;
	for (i = 0; i < n; i++) {
		x_std[i] = a * x[i] + b;
	}
	return GMRFLib_SUCCESS;
}
int GMRFLib_gsl_integration_fix_limits(double *new_lower, double *new_upper, gsl_function * F, double lower, double upper)
{
	/*
	 *  compute new and improver lower and upper limits for the integration. skip parts with |feval| < eps. 
	 */
	double eps = 1.0e-10, min_step = GMRFLib_eps(1. / 3.), step, x, val;
	int nstep = 20, debug = 0, i;

	if (lower == upper) {
		*new_lower = lower;
		*new_upper = upper;

		return GMRFLib_SUCCESS;
	}
	if (lower > upper) {
		return GMRFLib_gsl_integration_fix_limits(new_upper, new_lower, F, upper, lower);
	}

	step = (upper - lower) / nstep;
	if (debug) {
		printf("%s: lower %g upper %g, seach with step_size = %g\n", __GMRFLib_FuncName, lower, upper, step);
		fflush(stdout);
	}

	*new_lower = lower;
	if (step > min_step) {
		for (i = 0; i < nstep; i++) {
			x = lower + i * step;
			val = GSL_FN_EVAL(F, x);
			if (debug) {
				printf("\tsearch lower x %.12g F %.12g\n", x, val);
				fflush(stdout);
			}
			if (ABS(val) < eps) {
				*new_lower = x;
			} else {
				break;
			}
		}
	}

	*new_upper = upper;
	if (step > min_step) {
		for (i = 0; i < nstep; i++) {
			x = upper - i * step;
			val = GSL_FN_EVAL(F, x);
			if (debug) {
				printf("\tsearch upper x %.12g F %.12g\n", x, val);
				fflush(stdout);
			}
			if (ABS(val) < eps) {
				*new_upper = x;
			} else {
				break;
			}
		}
	}

	if (*new_upper < *new_lower) {
		/*
		 * then the function is all below < eps. In this case, do nothing.
		 */
		*new_lower = lower;
		*new_upper = upper;
	}
	if (debug) {
		printf("\tfound new integration limits: %.12g %.12g\n", *new_lower, *new_upper);
		fflush(stdout);
	}

	return GMRFLib_SUCCESS;
}
int GMRFLib_density_layout_x(double **x_vec, int *len_x, GMRFLib_density_tp * density)
{
	/*
	 * return points for printing the marginals. this is on a standarised scale, so the SD is one. 
	 */
	double p_few[] = {
		// 0.0000001,
		// 0.000001,
		// 0.00001,
		// 0.0001,
		// 0.0005,
		// 0.001,
		// 0.005,
		0.01,
		// 0.025, 
		0.05,
		// 0.075, 
		0.10,
		// 0.125, 
		0.15,
		// 0.175, 
		0.2,
		// 0.225,
		0.25,
		// 0.275,
		0.30,
		// 0.325,
		0.35,
		// 0.375,
		0.40,
		// 0.425,
		0.45,
		// 0.46,
		// 0.47,
		0.475,
		// 0.48,
		// 0.49,
		0.50,
		// 0.51,
		// 0.52,
		0.525,
		// 0.53,
		// 0.54,
		0.55,
		// 0.575,
		0.60,
		// 0.625,
		0.65,
		// 0.675,
		0.70,
		// 0.725,
		0.75,
		// 0.775,
		0.80,
		// 0.825,
		0.85,
		// 0.875,
		0.9,
		// 0.925, 
		0.95,
		// 0.975,
		0.99,
		// 0.995, 
		// 0.999,
		// 0.9995,
		// 0.9999,
		// 0.99999,
		// 0.999999
		// 0.9999999
	};

	double p_many[] = {
		0.0000001,
		0.000001,
		0.00001,
		0.0001,
		0.0005,
		0.001,
		0.005,
		0.01,
		0.025,
		0.05,
		0.075,
		0.10,
		0.125,
		0.15,
		0.175,
		0.2,
		0.225,
		0.25,
		0.275,
		0.30,
		0.325,
		0.35,
		0.375,
		0.40,
		0.425,
		0.45,
		0.46,
		0.47,
		0.475,
		0.48,
		0.49,
		0.50,
		0.51,
		0.52,
		0.525,
		0.53,
		0.54,
		0.55,
		0.575,
		0.60,
		0.625,
		0.65,
		0.675,
		0.70,
		0.725,
		0.75,
		0.775,
		0.80,
		0.825,
		0.85,
		0.875,
		0.9,
		0.925,
		0.95,
		0.975,
		0.99,
		0.995,
		0.999,
		0.9995,
		0.9999,
		0.99999,
		0.999999,
		0.9999999
	};

	double x_add[] = {
		-10.0, -8.0, -6.0, -5.0, -4.0, -3.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0
	};

	int i, m, use_many;

	if ((GMRFLib_density_storage_strategy == GMRFLib_DENSITY_STORAGE_STRATEGY_DEFAULT ||
	     GMRFLib_density_storage_strategy == GMRFLib_DENSITY_STORAGE_STRATEGY_HIGH)) {
		use_many = GMRFLib_TRUE;
	} else {
		use_many = GMRFLib_FALSE;
	}

	double *p_ptr = (use_many ? p_many : p_few);

	*len_x = (use_many ? sizeof(p_many) + sizeof(x_add) : sizeof(p_few) + sizeof(x_add)) / sizeof(double);
	m = (use_many ? sizeof(p_many) : sizeof(p_few)) / sizeof(double);
	*x_vec = Calloc(*len_x, double);
	for (i = 0; i < m; i++) {
		GMRFLib_density_Pinv(&((*x_vec)[i]), p_ptr[i], density);
	}
	for (i = m; i < *len_x; i++) {
		(*x_vec)[i] = x_add[i - m];
	}
	qsort(*x_vec, (size_t) (*len_x), sizeof(double), GMRFLib_dcmp);
	GMRFLib_unique_additive(len_x, *x_vec, GMRFLib_eps(0.5));

	return GMRFLib_SUCCESS;
}
