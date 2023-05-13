
/* density.c
 * 
 * Copyright (C) 2006-2023 Havard Rue
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

#include <stddef.h>
#include <float.h>
#include <math.h>
#include <strings.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

#define CONST_1 0.6266570686577500604			       // sqrt(M_PI/8.0);
#define CONST_2 (-0.69314718055994528623)		       // log(0.5);
#define CONST_3 (-0.918938533204672741780329736407)	       // log(1.0/sqrt(2.0*M_PI))
#define CONST_4 (0.398942280401432677939946059934)	       // 1.0/sqrt(2.0*M_PI)

int GMRFLib_sn_par2moments(double *mean, double *stdev, double *skewness, GMRFLib_sn_param_tp *p)
{
	/*
	 * compute two first moments of the sn 
	 */
	double delta, c1, c2, c3;

	delta = p->alpha / sqrt(1.0 + SQR(p->alpha));
	c1 = 1.0 - 2.0 * SQR(delta) / M_PI;
	c2 = 0.79788456080286535588;			       // sqrt(2.0/M_PI)
	c3 = 0.4292036732051033808;			       // (4 - M_PI)/2.0

	// *stdev = p->omega * sqrt(1.0 - 2.0 * SQR(delta) / M_PI);
	// *skewness = (4.0 - M_PI) / 2.0 * pow(delta * sqrt(2.0 / M_PI), 3.0) / pow(1.0 - 2 * SQR(delta) / M_PI, 3.0 / 2.0);

	*mean = p->xi + p->omega * c2 * delta;
	*stdev = p->omega * sqrt(c1);
	*skewness = c3 * pow(delta * c2, 3.0) / pow(c1, 1.5);

	return GMRFLib_SUCCESS;
}

int GMRFLib_sn_moments2par(GMRFLib_sn_param_tp *p, double *mean, double *stdev, double *skew)
{
	double delta, c1, c2, c3, ss = pow(ABS(*skew), 2.0 / 3.0);

	c1 = 0.63661977236758134306;			       // 2.0 / M_PI
	c2 = 0.79788456080286535588;			       // sqrt(2.0/M_PI)
	c3 = 0.4292036732051033808;			       // (4 - M_PI)/2.0

	delta = (*skew < 0 ? -1.0 : 1.0) * sqrt(M_PI_2 * ss / (ss + pow(c3, 2.0 / 3.0)));
	p->alpha = delta / sqrt(1.0 - SQR(delta));
	p->omega = *stdev / sqrt(1.0 - c1 * SQR(delta));
	p->xi = *mean - p->omega * delta * c2;

	return GMRFLib_SUCCESS;
}


GMRFLib_idxval_tp *GMRFLib_density_prune_weights(double *weights, int n, double prob)
{
	// return an idxval with some of the weights pruned off, so that the sum is at least prob

	size_t one = 1;
	int nn;
	double ww_sum = 0.0;
	double *ww = Calloc(n, double);

	Memcpy(ww, weights, n * sizeof(double));
	GMRFLib_normalize(n, ww);

	size_t *perm = Calloc(n, size_t);
	gsl_sort_index(perm, ww, one, (size_t) n);

	prob = TRUNCATE(prob, 0.0, 1.0);
	ww_sum = 0;
	nn = 0;
	for (int i = n - 1; i >= 0; i--) {		       /* since 'perm' is increasing */
		int j = perm[i];
		if (ww_sum < prob) {
			ww_sum += ww[j];
			nn++;
		} else {
			ww[j] = 0.0;
		}
	}
	GMRFLib_normalize(n, ww);

	GMRFLib_idxval_tp *idxval = NULL;
	GMRFLib_idxval_create_x(&idxval, nn);

	for (int i = 0; i < n; i++) {
		if (!ISZERO(ww[i])) {
			GMRFLib_idxval_add(&idxval, i, ww[i]);
		}
	}

	Free(ww);
	Free(perm);

	return idxval;
}

int GMRFLib_sn_density(double *dens, double x, void *param)
{
	/*
	 * compute the sn-density 
	 */
	double ldens = 0.0;

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
	GMRFLib_sn_param_tp *p = (GMRFLib_sn_param_tp *) param;
	double z = (x - p->xi) / p->omega;

	*ldens = M_LN2 + CONST_3 - 0.5 * SQR(z) + GMRFLib_log_gsl_cdf_ugaussian_P(p->alpha * z) - log(p->omega);

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

int GMRFLib_sn_fit(GMRFLib_sn_param_tp *param, double *fval, double *x, double *log_density, int n)
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
	GMRFLib_sn_param_tp param = { 0.0, 0.0, 0.0 };

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
		printf ("iter: %3d x = %.16f %.16f %.16f %.16f " "|f(x)| = %g\n", \
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

	const int debug = 0;
	int iter = 0, status, i, imax;
	double eps = GSL_ROOT3_DBL_EPSILON, *log_density_scaled;

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

	double ld_max = log_density[imax];
#pragma GCC ivdep
	for (i = 0; i < (int) n; i++) {
		log_density_scaled[i] = log_density[i] - ld_max;
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

int GMRFLib_sn_fit_f(const gsl_vector *param, void *data, gsl_vector *f)
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

int GMRFLib_sn_fit_df(const gsl_vector *param, void *data, gsl_matrix *J)
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

int GMRFLib_sn_fit_fdf(const gsl_vector *param, void *data, gsl_vector *f, gsl_matrix *J)
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

int GMRFLib_init_density(GMRFLib_density_tp *density, int lookup_tables)
{
	/*
	 * initialize 'density': compute the mean, stdev and the norm_const (for the log spline fit) 
	 */

	int np = GMRFLib_INT_NUM_POINTS;
	int npm = GMRFLib_INT_NUM_INTERPOL * np - (GMRFLib_INT_NUM_INTERPOL - 1);
	int i;
	double low = 0.0, high = 0.0, xval, *xpm = NULL, *ld = NULL, *ldm = NULL, *pm = NULL, *xp = NULL, dx = 0.0, dxm = 0.0, d0, d1;

	if (!density) {
		return GMRFLib_SUCCESS;
	}

	// GMRFLib_ENTER_ROUTINE;
	Calloc_init(4 * npm + 2 * np, 6);

	if (density->type == GMRFLib_DENSITY_TYPE_GAUSSIAN) {
		// density->mean = density->mean_gaussian;
		// density->stdev = density->stdev_gaussian;
		density->skewness = 0.0;
		density->x_min = -GMRFLib_DENSITY_INTEGRATION_LIMIT * density->stdev + density->mean;
		density->x_max = GMRFLib_DENSITY_INTEGRATION_LIMIT * density->stdev + density->mean;
	} else if (density->type == GMRFLib_DENSITY_TYPE_SKEWNORMAL) {
		/*
		 * for the skew-normal we know the moments 
		 */
		double mom[3] = { 0, 0, 0 };
		GMRFLib_sn_par2moments(&mom[0], &mom[1], &mom[2], density->sn_param);
		density->mean = mom[0];
		density->stdev = mom[1];
		density->skewness = mom[2];
		density->x_min = -GMRFLib_DENSITY_INTEGRATION_LIMIT * density->stdev + density->mean;
		density->x_max = GMRFLib_DENSITY_INTEGRATION_LIMIT * density->stdev + density->mean;
	}

	/*
	 * the mean and the stdev in the users scale 
	 */
	density->user_mean = density->std_stdev * density->mean + density->std_mean;
	density->user_stdev = density->std_stdev * density->stdev;

	if (density->type == GMRFLib_DENSITY_TYPE_GAUSSIAN) {
		density->user_mode = density->user_mean;
	} else {
		density->user_mode = NAN;		       /* yes, this is the value if its not computed */
	}

	if (density->type == GMRFLib_DENSITY_TYPE_GAUSSIAN) {
		Calloc_free();
		density->P = density->Pinv = NULL;
		// GMRFLib_LEAVE_ROUTINE;
		return GMRFLib_SUCCESS;
	}

	if (!lookup_tables && density->type != GMRFLib_DENSITY_TYPE_SCGAUSSIAN) {
		Calloc_free();
		density->P = density->Pinv = NULL;
		// GMRFLib_LEAVE_ROUTINE;
		return GMRFLib_SUCCESS;
	}

	low = density->x_min;
	high = density->x_max;
	dx = (high - low) / (np - 1.0);
	dxm = (high - low) / (npm - 1.0);

	double x_max;
	double ldmax;
	double w[2] = { 4.0, 2.0 };

	xp = Calloc_get(np);
	ld = Calloc_get(np);

	for (xval = low, i = 0; i < np; xval += dx, i++) {
		xp[i] = xval;
	}
	density->log_norm_const = 0.0;
	GMRFLib_evaluate_nlogdensity(ld, xp, np, density);
	ldmax = GMRFLib_max_value(ld, np, NULL);
	GMRFLib_adjust_vector(ld, np);

	// interpolate
	xpm = Calloc_get(npm);
	ldm = Calloc_get(npm);
	pm = Calloc_get(npm);

	if (GMRFLib_INT_NUM_INTERPOL == 3) {
		const double div3 = 1.0 / 3.0;
#pragma GCC ivdep
		for (i = 0; i < np - 1; i++) {
			xpm[3 * i + 0] = xp[i];
			xpm[3 * i + 1] = (2.0 * xp[i] + xp[i + 1]) * div3;
			xpm[3 * i + 2] = (xp[i] + 2.0 * xp[i + 1]) * div3;
		}
#pragma GCC ivdep
		for (i = 0; i < np - 1; i++) {
			ldm[3 * i + 0] = ld[i];
			ldm[3 * i + 1] = (2.0 * ld[i] + ld[i + 1]) * div3;
			ldm[3 * i + 2] = (ld[i] + 2.0 * ld[i + 1]) * div3;
		}
		xpm[3 * (np - 2) + 3] = xp[np - 1];
		ldm[3 * (np - 2) + 3] = ld[np - 1];
		assert(3 * (np - 2) + 3 == npm - 1);
	} else if (GMRFLib_INT_NUM_INTERPOL == 2) {
		const double div2 = 0.5;
#pragma GCC ivdep
		for (i = 0; i < np - 1; i++) {
			xpm[2 * i + 0] = xp[i];
			xpm[2 * i + 1] = (xp[i] + xp[i + 1]) * div2;
		}
#pragma GCC ivdep
		for (i = 0; i < np - 1; i++) {
			ldm[2 * i + 0] = ld[i];
			ldm[2 * i + 1] = (ld[i] + ld[i + 1]) * div2;
		}
		xpm[2 * (np - 2) + 2] = xp[np - 1];
		ldm[2 * (np - 2) + 2] = ld[np - 1];
		assert(2 * (np - 2) + 2 == npm - 1);
	} else {
		assert(GMRFLib_INT_NUM_INTERPOL == 2 || GMRFLib_INT_NUM_INTERPOL == 3);
	}

	// convert scale
#if defined(INLA_LINK_WITH_MKL)
	vdExp(npm, ldm, ldm);
#else
#pragma GCC ivdep
	for (i = 0; i < npm; i++) {
		ldm[i] = exp(ldm[i]);
	}
#endif

	int idx_max = 0;
	GMRFLib_max_value(ldm, npm, &idx_max);
	if (idx_max == 0) {
		x_max = xpm[0];
	} else if (idx_max == npm - 1) {
		x_max = xpm[npm - 1];
	} else {
		double *xx = xpm + idx_max - 1;
		double *tld = ldm + idx_max - 1;
		// see inla.c and 'inla_integrate_func'
		x_max = (tld[0] * xx[1] * xx[1] - tld[0] * xx[2] * xx[2] - tld[1] * xx[0] * xx[0] +
			 tld[1] * xx[2] * xx[2] + tld[2] * xx[0] * xx[0] - tld[2] * xx[1] * xx[1]) /
		    (tld[0] * xx[1] - tld[0] * xx[2] - tld[1] * xx[0] + tld[1] * xx[2] + tld[2] * xx[0] - xx[1] * tld[2]) / 0.2e1;
	}

	// compute moments
	double mm[4] = { 0.0, 0.0, 0.0, 0.0 };
	double xx[4] = { 0.0, 0.0, 0.0, 0.0 };

	d0 = ldm[0];
	d1 = ldm[npm - 1];
	xx[0] = xpm[0];
	xx[1] = xpm[npm - 1];
	mm[0] = d0 + d1;
	mm[1] = xx[0] * d0 + xx[1] * d1;
	mm[2] = SQR(xx[0]) * d0 + SQR(xx[1]) * d1;
	mm[3] = gsl_pow_3(xx[0]) * d0 + gsl_pow_3(xx[1]) * d1;

	for (i = 1; i < npm - 1; i++) {
		double d = ldm[i] * w[(i - 1) % 2];
		xx[1] = xpm[i];
		xx[2] = SQR(xx[1]);
		xx[3] = xx[2] * xx[1];

		mm[0] += d;
		mm[1] += xx[1] * d;
		mm[2] += xx[2] * d;
		mm[3] += xx[3] * d;
	}
	mm[1] /= mm[0];
	mm[2] /= mm[0];
	mm[3] /= mm[0];

	if (density->type == GMRFLib_DENSITY_TYPE_SCGAUSSIAN) {
		density->log_norm_const = log(mm[0] * dxm / 3.0) + ldmax;
		density->mean = mm[1];
		density->stdev = sqrt(DMAX(0.0, mm[2] - SQR(mm[1])));
		density->skewness = (mm[3] - 3.0 * mm[1] * SQR(density->stdev) - gsl_pow_3(mm[1])) / gsl_pow_3(density->stdev);
		density->x_min = -GMRFLib_DENSITY_INTEGRATION_LIMIT * density->stdev + density->mean;
		density->x_max = GMRFLib_DENSITY_INTEGRATION_LIMIT * density->stdev + density->mean;
		density->user_mean = density->std_stdev * density->mean + density->std_mean;
		density->user_stdev = density->std_stdev * density->stdev;
		density->user_mode = x_max * density->std_stdev + density->std_mean;
	}

	if (lookup_tables) {
		pm = Calloc_get(npm);
		pm[0] = 0.5 * ldm[0];
		for (i = 1; i < npm; i++) {
			pm[i] = pm[i - 1] + 0.5 * (ldm[i] + ldm[i - 1]);
		}

		double cc = 1.0 / (pm[npm - 1] + 0.5 * ldm[npm - 1]);

		// for (i = 0; i < npm; i++) pm[i] *= cc;
		GMRFLib_dscale(npm, cc, pm);

		// shrink before creating the spline
		int k = 0;
		for (i = 0; i < npm; i += GMRFLib_INT_NUM_INTERPOL) {
			xpm[k] = xpm[i];
			pm[k] = pm[i];
			k++;
		}
		npm = k;
		density->Pinv = GMRFLib_spline_create_x(pm, xpm, npm, GMRFLib_INTPOL_TRANS_Pinv);
		//density->P = GMRFLib_spline_create_x(xpm, pm, npm, GMRFLib_INTPOL_TRANS_P);
		density->P = NULL;
	}

	Calloc_free();
	// GMRFLib_LEAVE_ROUTINE;

	return GMRFLib_SUCCESS;
}

int GMRFLib_evaluate_logdensity(double *logdens, double x, GMRFLib_density_tp *density)
{
	return GMRFLib_evaluate_nlogdensity(logdens, &x, 1, density);
}

int GMRFLib_evaluate_nlogdensity(double *logdens, double *x, int n, GMRFLib_density_tp *density)
{
	/*
	 * evaluate the log-density-object. Note that x is in *standardised scale* . 
	 */
	int i;

	assert(density);

	switch (density->type) {
	case GMRFLib_DENSITY_TYPE_GAUSSIAN:
	{
		double c1 = CONST_3 - log(density->stdev);
		double c2 = -0.5 / SQR(density->stdev);
		double m = density->mean;
#pragma GCC ivdep
		for (i = 0; i < n; i++) {
			// logdens[i] = log_norm_const_gaussian - log(density->stdev) - 0.5 * SQR(x[i] - density->mean) / SQR(density->stdev);
			logdens[i] = c1 + c2 * SQR(x[i] - m);
		}
	}
		break;

	case GMRFLib_DENSITY_TYPE_SKEWNORMAL:
	{
		/*
		 * just inline GMRFLib_sn_logdensity() 
		 */

		GMRFLib_sn_param_tp *p = density->sn_param;
		double local_const_1 = M_LN2 + CONST_3 - log(p->omega);
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
					val = CONST_2 + log1p(0.0 + sqrt(1.0 - exp(-CONST_1 * SQR(zz))));
				} else {
					// val = log(0.5 - 0.5 * sqrt(1.0 - exp(- CONST_1 * SQR(zz))));
					val = CONST_2 + log1p(0.0 - sqrt(1.0 - exp(-CONST_1 * SQR(zz))));
				}
			} else {
				// use the more complitated asympt expression, which we do here (for which the code in
				// the prev {} is a copy of
				val = GMRFLib_log_gsl_cdf_ugaussian_P(zz);
			}
			logdens[i] = local_const_1 - 0.5 * SQR(z) + val;
		}
	}
		break;

	case GMRFLib_DENSITY_TYPE_SCGAUSSIAN:
	{
		for (i = 0; i < n; i++) {
			logdens[i] = GMRFLib_spline_eval(x[i], density->log_correction) - 0.5 * SQR(x[i]) - density->log_norm_const;
		}
	}
		break;

	default:

		FIXME("Unknown type");
		P(density->type);
		GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
		break;
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_evaluate_density(double *dens, double x, GMRFLib_density_tp *density)
{
	double ldens = 0.0;

	GMRFLib_evaluate_logdensity(&ldens, x, density);
	*dens = exp(ldens);

	return GMRFLib_SUCCESS;
}

int GMRFLib_evaluate_ndensity(double *dens, double *x, int n, GMRFLib_density_tp *density)
{
	assert(dens);

	GMRFLib_evaluate_nlogdensity(dens, x, n, density);

#if defined(INLA_LINK_WITH_MKL)
	vdExp(n, dens, dens);
#else
#pragma GCC ivdep
	for (int i = 0; i < n; i++) {
		dens[i] = exp(dens[i]);
	}
#endif
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

int GMRFLib_free_density(GMRFLib_density_tp *density)
{
	if (density) {
		switch (density->type) {
		case GMRFLib_DENSITY_TYPE_GAUSSIAN:
			break;

		case GMRFLib_DENSITY_TYPE_SKEWNORMAL:
		{
			Free(density->sn_param);
		}
			break;

		case GMRFLib_DENSITY_TYPE_SCGAUSSIAN:
		{
			GMRFLib_spline_free(density->log_correction);
		}
			break;

		default:
			GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
			break;
		}
		if (density->P) {
			GMRFLib_spline_free(density->P);
		}
		if (density->Pinv) {
			GMRFLib_spline_free(density->Pinv);
		}
		Free(density);
	}

	return GMRFLib_SUCCESS;
}

double GMRFLib_evaluate_density_kld__intern(double x, void *param)
{
	GMRFLib_density_properties_tp *prop = (GMRFLib_density_properties_tp *) param;
	double ldens = 0.0;
	GMRFLib_evaluate_logdensity(&ldens, x, prop->density);
	return CONST_4 * exp(-0.5 * SQR(x)) * (CONST_3 - 0.5 * SQR(x) - ldens);
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

int GMRFLib_density_Pinv(double *xp, double p, GMRFLib_density_tp *density)
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
			*xp = GMRFLib_spline_eval(p, density->Pinv);
		} else {
			GMRFLib_ASSERT(density->Pinv->spline != NULL, GMRFLib_ESNH);
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_density_P(double *px, double x, GMRFLib_density_tp *density)
{
	/*
	 * compute px = Prob (X <= x) for the `density'
	 * 
	 * NOTE that 'x' is in standarized scale. 
	 */
	double result = 0.0;

	GMRFLib_ENTER_ROUTINE;

	/*
	 * if the density is Gaussian, then its easy 
	 */
	if (density->type == GMRFLib_DENSITY_TYPE_GAUSSIAN) {
		result = gsl_cdf_ugaussian_P((x - density->mean) / density->stdev);
	} else {
		result = GMRFLib_spline_eval(x, density->P);
	}

	*px = result;
	GMRFLib_LEAVE_ROUTINE;

	return GMRFLib_SUCCESS;
}

int GMRFLib_evaluate_densities(double *dens, double x_user, int n, GMRFLib_density_tp **densities, double *weights)
{
	/*
	 * evaluate the density in ***USER SCALE*** at `x_user', where the density is given as a mixture
	 * 
	 * \sum_{i=0}^{n-1} weights[i]*densities[i]
	 * 
	 * the weights need not to be scaled. 
	 */
	int i, j;
	double d_tmp = 0.0, d = 0.0, x_std, p;
	GMRFLib_idxval_tp *probs = GMRFLib_density_prune_weights(weights, n, GMRFLib_weight_prob);

	for (j = 0; j < probs->n; j++) {
		i = probs->idx[j];
		p = probs->val[j];
		x_std = GMRFLib_density_user2std(x_user, densities[i]);
		GMRFLib_evaluate_density(&d_tmp, x_std, densities[i]);
		d += p * d_tmp / densities[i]->std_stdev;
	}
	*dens = d;
	GMRFLib_idxval_free(probs);

	return GMRFLib_SUCCESS;
}

int GMRFLib_evaluate_ndensities(double *dens, double *x_user, int nx, GMRFLib_density_tp **densities, GMRFLib_idxval_tp *probs)
{
	/*
	 * evaluate the density in ***USER SCALE*** at `x_user[j]', j=0...nx, where the density is given as a mixture
	 * 
	 * \sum_{i=0}^{n-1} weights[i]*densities[i]
	 * 
	 * the weights need not to be scaled. 
	 */
	int nd = (probs ? probs->n : 1);
	int i, k, n_alloc = IMAX(nd, nx);
	double *d_tmp = NULL, *x_std = NULL, p;

	Calloc_init(2 * n_alloc, 2);
	d_tmp = Calloc_get(n_alloc);
	x_std = Calloc_get(n_alloc);

	Memset(dens, 0, nx * sizeof(double));
	for (k = 0; k < probs->n; k++) {
		i = probs->idx[k];
		p = probs->val[k];
		GMRFLib_density_user2std_n(x_std, x_user, densities[i], nx);
		GMRFLib_evaluate_ndensity(d_tmp, x_std, nx, densities[i]);

		double a = p / densities[i]->std_stdev;
		int inc = 1;
		daxpy_(&nx, &a, d_tmp, &inc, dens, &inc);

		// for (j = 0; j < nx; j++) 
		// dens[j] += p * d_tmp[j] / densities[i]->std_stdev;
	}

	Calloc_free();
	return GMRFLib_SUCCESS;
}

int GMRFLib_density_duplicate(GMRFLib_density_tp **density_to, GMRFLib_density_tp *density_from)
{
	GMRFLib_density_combine(density_to, &density_from, NULL);
	(*density_to)->flags = density_from->flags;

	return GMRFLib_SUCCESS;
}

int GMRFLib_density_combine(GMRFLib_density_tp **density, GMRFLib_density_tp **densities, GMRFLib_idxval_tp *probs)
{
	/*
	 * make a new spline-corrected-gaussian density out of a weighted sum of densities and return this in DENSITY.  make a
	 * new spline-corrected-gaussian density out of a weighted sum of gaussian densities and return this in GDENSITY.
	 *
	 * \sum_{i=0}^{n-1} weights[i]*densities[i]
	 * 
	 * the weights need not to be scaled. 
	 */

	// probs == NULL means n=1
	int n = (probs ? probs->n : 1);

	if (n == 0) {
		if (density) {
			*density = NULL;
		}
		return GMRFLib_SUCCESS;
	}

	// GMRFLib_ENTER_ROUTINE;

	// this actually happens like for 'eb' and is also how 'duplicate' is implemented
	if (n == 1) {
		if ((*densities)->type == GMRFLib_DENSITY_TYPE_GAUSSIAN) {
			GMRFLib_density_create_normal(density, (*densities)->mean, (*densities)->stdev,
						      (*densities)->std_mean, (*densities)->std_stdev,
						      ((*densities)->P || (*densities)->Pinv ? 1 : 0));
		} else if ((*densities)->type == GMRFLib_DENSITY_TYPE_SKEWNORMAL) {
			GMRFLib_density_create_sn(density, *((*densities)->sn_param),
						  (*densities)->std_mean, (*densities)->std_stdev, ((*densities)->P || (*densities)->Pinv ? 1 : 0));
		} else if ((*densities)->type == GMRFLib_DENSITY_TYPE_SCGAUSSIAN) {
			int m = 0;
			double *x = NULL, *ld = NULL;

			GMRFLib_density_layout_x(NULL, &m, NULL);
			Calloc_init(2 * m, 2);
			x = Calloc_get(m);
			ld = Calloc_get(m);
			GMRFLib_density_layout_x(x, &m, *densities);
			GMRFLib_evaluate_nlogdensity(ld, x, m, *densities);
			GMRFLib_density_create(density, GMRFLib_DENSITY_TYPE_SCGAUSSIAN, m, x, ld,
					       (*densities)->std_mean, (*densities)->std_stdev, ((*densities)->P || (*densities)->Pinv ? 1 : 0));
			Calloc_free();
		} else {
			FIXME("Unknown type");
			P((*densities)->type);
			assert(0 == 1);
		}
		// GMRFLib_LEAVE_ROUTINE;
		return GMRFLib_SUCCESS;
	}

	double mean, stdev, *ddens = NULL, *log_dens = NULL, *xx_real = NULL, m1, m2, sum_w, p;
	double xx[] = { -5.0, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.25, -1.0, -0.75, -0.5, -0.25, -0.125, 0.0,
		0.125, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0
	};
	int nx = sizeof(xx) / sizeof(double);

	/*
	 * compute the mean and variance in the user-scale 
	 */
	m1 = m2 = sum_w = 0.0;
	for (int ii = 0; ii < probs->n; ii++) {
		int i = probs->idx[ii];
		p = probs->val[ii];
		m1 += p * densities[i]->user_mean;
		m2 += p * (SQR(densities[i]->user_stdev) + SQR(densities[i]->user_mean));
		sum_w += p;
	}
	mean = m1 / sum_w;
	stdev = sqrt(DMAX(0.0, m2 / sum_w - SQR(mean)));

	Calloc_init(2 * nx, 2);
	/*
	 * compute the weighted density. note that we have to go through the user/real-scale to get this right 
	 */
	xx_real = Calloc_get(nx);
	ddens = Calloc_get(nx);
	log_dens = ddens;				       /* same storage */
#pragma GCC ivdep
	for (int i = 0; i < nx; i++) {
		xx_real[i] = xx[i] * stdev + mean;
	}

	GMRFLib_evaluate_ndensities(ddens, xx_real, nx, densities, probs);
#pragma GCC ivdep
	for (int i = 0; i < nx; i++) {
		log_dens[i] = (ddens[i] > 0.0 ? log(ddens[i]) : -FLT_MAX);
	}
	GMRFLib_density_create(density, GMRFLib_DENSITY_TYPE_SCGAUSSIAN, nx, xx, log_dens, mean, stdev, GMRFLib_TRUE);

	Calloc_free();

	// GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}


int GMRFLib_density_create_normal(GMRFLib_density_tp **density, double mean, double stdev, double std_mean, double std_stdev, int lookup_tables)
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
	(*density)->mean = mean;
	(*density)->stdev = stdev;

	GMRFLib_init_density(*density, lookup_tables);

	return GMRFLib_SUCCESS;
}

int GMRFLib_density_create_sn(GMRFLib_density_tp **density, GMRFLib_sn_param_tp sn_param, double std_mean, double std_stdev, int lookup_tables)
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
	Memcpy((void *) (*density)->sn_param, (const void *) &sn_param, sizeof(GMRFLib_sn_param_tp));

	GMRFLib_init_density(*density, lookup_tables);

	return GMRFLib_SUCCESS;
}

int GMRFLib_density_adjust_vector(double *ldens, int n)
{
	double maxdev = 2.0 * log(DBL_EPSILON);

	GMRFLib_adjust_vector(ldens, n);
	for (int i = 0; i < n; i++) {
		ldens[i] = DMAX(maxdev, ldens[i]);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_density_create(GMRFLib_density_tp **density, int type, int n, double *x, double *logdens, double std_mean,
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
	int i, j;
	const int debug = 0;
	double *xx = NULL, *ldens = NULL, g_mean = 0.0, g_var = 1.0;
	GMRFLib_sn_param_tp sn_param = { 0, 0, 0 };

	Calloc_init(2 * n, 2);
	xx = Calloc_get(n);
	ldens = Calloc_get(n);

	Memcpy(xx, x, (size_t) n * sizeof(double));
	Memcpy(ldens, logdens, (size_t) n * sizeof(double));

	/*
	 * sort xx and remove ties. that that we need to sort first. In most cases we do not
	 */
	int is_sorted = 1;
	for (i = 1; i < n && is_sorted; i++) {
		is_sorted = (xx[i] > xx[i - 1]);
	}
	if (!is_sorted) {
		// gsl_sort2(xx, (size_t) 1, ldens, (size_t) 1, (size_t) n);
		my_sort2_dd(xx, ldens, n);
	}
	GMRFLib_unique_relative2(&n, xx, ldens, GSL_SQRT_DBL_EPSILON);
	GMRFLib_adjust_vector(ldens, n);

	if (debug) {
		printf("%s: Create density n %d type %d\n", __GMRFLib_FuncName, n, type);
		for (int ii = 0; ii < n; ii++) {
			printf("\tx %f ldens %f\n", xx[ii], ldens[ii]);
		}
	}

	if (n == 1) {
		GMRFLib_density_create_normal(density, 0.0, 1.0, std_mean, std_stdev, lookup_tables);
	} else {
		switch (type) {
		case GMRFLib_DENSITY_TYPE_GAUSSIAN:
		{
			/*
			 * fit a gaussian 
			 */
			GMRFLib_normal_fit(&g_mean, &g_var, NULL, xx, ldens, n);
			GMRFLib_density_create_normal(density, g_mean, sqrt(g_var), std_mean, std_stdev, lookup_tables);
		}
			break;

		case GMRFLib_DENSITY_TYPE_SKEWNORMAL:
		{
			/*
			 * fit skew-normal 
			 */
			GMRFLib_sn_fit(&sn_param, NULL, xx, ldens, n);
			GMRFLib_density_create_sn(density, sn_param, std_mean, std_stdev, lookup_tables);
		}
			break;

		case GMRFLib_DENSITY_TYPE_SCGAUSSIAN:
		{
			/*
			 * fit spline-corrected gaussian 
			 */
			for (i = j = 0; i < n; i++) {	       /* this could happen */
				if (!ISINF(ldens[i])) {
					ldens[j] = ldens[i];
					xx[j] = xx[i];
					j++;
				}
			}
			n = j;

			(*density) = Calloc(1, GMRFLib_density_tp);
			(*density)->type = GMRFLib_DENSITY_TYPE_SCGAUSSIAN;
			(*density)->std_mean = std_mean;
			(*density)->std_stdev = std_stdev;
			(*density)->x_min = GMRFLib_min_value(xx, n, NULL);
			(*density)->x_max = GMRFLib_max_value(xx, n, NULL);


#pragma GCC ivdep
			for (i = 0; i < n; i++) {
				ldens[i] += 0.5 * SQR(xx[i]);  /* ldens is now the correction */
			}
			(*density)->log_correction = GMRFLib_spline_create(xx, ldens, n);
			GMRFLib_init_density(*density, lookup_tables);
		}
			break;

		default:
			GMRFLib_ASSERT(1 == 0, GMRFLib_ESNH);
		}
	}

	Calloc_free();
	return GMRFLib_SUCCESS;
}

int GMRFLib_density_new_mean(GMRFLib_density_tp **new_density, GMRFLib_density_tp *density, double new_mean)
{
	/*
	 * return a new density, which the density with the given new mean 
	 */
#define N (60)
#define M (4)
	int i, n = N + 2 * M;
	double *x, *ld, alpha, eps[M] = { 1e-6, 1e-5, 1e-4, 1e-3 };

	Calloc_init(2 * n, 2);
	x = Calloc_get(n);
	ld = Calloc_get(n);

	for (i = 0; i < M; i++) {
		GMRFLib_density_Pinv(&x[i], eps[i], density);
		GMRFLib_density_Pinv(&x[M + i], 1.0 - eps[i], density);
	}
	for (i = 0; i < N; i++) {
		double delta = eps[M - 1] * 2.0;
		alpha = delta + (1.0 - 2.0 * delta) * i * (1.0 / (double) (N - 1.0));
		GMRFLib_density_Pinv(&x[2 * M + i], alpha, density);
	}
	GMRFLib_evaluate_nlogdensity(ld, x, n, density);
	GMRFLib_density_create(new_density, GMRFLib_DENSITY_TYPE_SCGAUSSIAN, n, x, ld, new_mean, density->std_stdev, GMRFLib_TRUE);

	Calloc_free();
#undef N
#undef M

	return GMRFLib_SUCCESS;
}

int GMRFLib_density_new_meansd(GMRFLib_density_tp **new_density, GMRFLib_density_tp *density, double new_mean, double new_stdev)
{
	/*
	 * return a new density, which the density with the given new mean and stdev
	 */

#define N (30)
#define M (4)
	int i, n = N + 2 * M;
	double *x, *ld, alpha, eps[M] = { 1e-6, 1e-5, 1e-4, 1e-3 };

	Calloc_init(2 * n, 2);
	x = Calloc_get(n);
	ld = Calloc_get(n);

	for (i = 0; i < M; i++) {
		GMRFLib_density_Pinv(&x[i], eps[i], density);
		GMRFLib_density_Pinv(&x[M + i], 1.0 - eps[i], density);
	}
	for (i = 0; i < N; i++) {
		alpha = 0.01 + 0.98 * (1.0 / (double) N) * i;
		GMRFLib_density_Pinv(&x[2 * M + i], alpha, density);
	}
	GMRFLib_evaluate_nlogdensity(ld, x, n, density);
	GMRFLib_density_create(new_density, GMRFLib_DENSITY_TYPE_SCGAUSSIAN, n, x, ld, new_mean, new_stdev, GMRFLib_TRUE);

	Calloc_free();
#undef N
#undef M
	return GMRFLib_SUCCESS;
}

int GMRFLib_density_new_user_mean(GMRFLib_density_tp *density, double new_user_mean)
{
	if (density) {
		double diff = new_user_mean - density->user_mean;
		density->std_mean += diff;
		density->user_mean = new_user_mean;
		if (!ISNAN(density->user_mode)) {
			density->user_mode += diff;
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_density_new_user_stdev(GMRFLib_density_tp *density, double new_user_stdev)
{
	assert(density->type == GMRFLib_DENSITY_TYPE_GAUSSIAN);
	if (density) {
		double diff = new_user_stdev / density->user_stdev;
		density->stdev *= diff;
		density->user_stdev *= diff;
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_density_printf(FILE *fp, GMRFLib_density_tp *density)
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
		fprintf(fp, "%-35s %16.10f\n", "Skewness", density->skewness);
		fprintf(fp, "%-35s %16.10f\n", "User mean", density->user_mean);
		fprintf(fp, "%-35s %16.10f\n", "User stdev", density->user_stdev);
		fprintf(fp, "%-35s %16.10f\n", "User mode", density->user_mode);
		fprintf(fp, "%-35s %16.10f\n", "Standarisation: mean", density->std_mean);
		fprintf(fp, "%-35s %16.10f\n", "Standarisation: stdev", density->std_stdev);
		switch (density->type) {
		case GMRFLib_DENSITY_TYPE_GAUSSIAN:
		{
			fprintf(fp, "%-35s %-20s\n", "Density type", "Gaussian");
			fprintf(fp, "     %-30s %16.10f\n", "mean", density->mean);
			fprintf(fp, "     %-30s %16.10f\n", "stdev", density->stdev);
		}
			break;

		case GMRFLib_DENSITY_TYPE_SKEWNORMAL:
		{
			fprintf(fp, "%-35s %-20s\n", "Density type", "Skew Normal");
			fprintf(fp, "\t%-35s %16.10f\n", "xi", density->sn_param->xi);
			fprintf(fp, "     %-30s %16.10f\n", "omega", density->sn_param->omega);
			fprintf(fp, "     %-30s %16.10f\n", "alpha", density->sn_param->alpha);
		}
			break;

		case GMRFLib_DENSITY_TYPE_SCGAUSSIAN:
		{
			fprintf(fp, "%-26s %-30s\n", "Density type", "Spline corrected Gaussian");
			fprintf(fp, "     %-30s %16.10f\n", "Log normalisation constant", density->log_norm_const);
			fprintf(fp, "     %-30s %16.10f\n", "x_min", density->x_min);
			fprintf(fp, "     %-30s %16.10f\n", "x_max", density->x_max);
		}
			break;

		default:
			break;
		}

		for (int i = 0; i < 8; i++) {
			fprintf(fp, "Flag %1d = %1d\n", i, GMRFLib_getbit(density->flags, i));
		}
	}
	return GMRFLib_SUCCESS;
}

double GMRFLib_evaluate_density_kld2__intern(double x, void *param)
{
	/*
	 * x is in user-scale. compute the log(f(x)/g(x))f(x) in user-scale 
	 */
	double x_std0, x_std1, ld0 = 0.0, ld1 = 0.0, ld0u, ld1u;
	GMRFLib_density_tp **d = (GMRFLib_density_tp **) param;

	x_std0 = GMRFLib_density_user2std(x, d[0]);
	x_std1 = GMRFLib_density_user2std(x, d[1]);

	GMRFLib_evaluate_logdensity(&ld0, x_std0, d[0]);
	GMRFLib_evaluate_logdensity(&ld1, x_std1, d[1]);
	ld0u = ld0 - log(d[0]->std_stdev);
	ld1u = ld1 - log(d[1]->std_stdev);

	return (ld0u - ld1u) * exp(ld0u);
}

int GMRFLib_kld(double *kld, GMRFLib_density_tp *density, GMRFLib_density_tp *ddensity)
{
	/*
	 * compute the KLD from density to ddensity, \int log(density/ddensity)*density dx 
	 */
	double result = 0.0, error, eps = GSL_SQRT_DBL_EPSILON, low0, low1, high0, high1, low, high;
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

int GMRFLib_mkld(double *mkld, GMRFLib_density_tp *density, GMRFLib_density_tp *ddensity)
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

int GMRFLib_kld_sym(double *kld_sym, GMRFLib_density_tp *density, GMRFLib_density_tp *ddensity)
{
	/*
	 * compute the symmetric KLD distance between density and ddensity 
	 */
	double kld0 = NAN, kld1 = NAN;
	GMRFLib_kld(&kld0, density, ddensity);
	GMRFLib_kld(&kld1, ddensity, density);
	*kld_sym = DMAX(0.0, (kld0 + kld1) / 2.0);

	return GMRFLib_SUCCESS;
}

int GMRFLib_mkld_sym(double *mkld_sym, GMRFLib_density_tp *density, GMRFLib_density_tp *ddensity)
{
	double mkld0, mkld1;
	GMRFLib_mkld(&mkld0, density, ddensity);
	GMRFLib_mkld(&mkld1, ddensity, density);
	*mkld_sym = DMAX(0.0, (mkld0 + mkld1) / 2.0);

	return GMRFLib_SUCCESS;
}

double GMRFLib_density_std2user(double x, GMRFLib_density_tp *density)
{
	return density->std_mean + x * density->std_stdev;
}

double GMRFLib_density_std2user_n(double *__restrict x_user, double *__restrict x, int n, GMRFLib_density_tp *__restrict density)
{
	double m = density->std_mean;
	double s = density->std_stdev;

	if (1) {
		GMRFLib_daxpb(n, s, x, m, x_user);
	} else {
#pragma GCC ivdep
		for (int i = 0; i < n; i++) {
			x_user[i] = m + x[i] * s;
		}
	}
	return GMRFLib_SUCCESS;
}

double GMRFLib_density_user2std(double x, GMRFLib_density_tp *density)
{
	return (x - density->std_mean) / density->std_stdev;
}

int GMRFLib_density_user2std_n(double *__restrict x_std, double *__restrict x, GMRFLib_density_tp *__restrict density, int n)
{
	double a = 1.0 / density->std_stdev;
	double b = -density->std_mean / density->std_stdev;

	if (1) {
		GMRFLib_daxpb(n, a, x, b, x_std);
	} else {
#pragma GCC ivdep
		for (int i = 0; i < n; i++) {
			x_std[i] = a * x[i] + b;
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_gsl_integration_fix_limits(double *new_lower, double *new_upper, gsl_function *F, double lower, double upper)
{
	/*
	 *  compute new and improver lower and upper limits for the integration. skip parts with |feval| < eps. 
	 */
	double eps = 1.0e-10, min_step = GSL_ROOT3_DBL_EPSILON, step, x, val;
	int nstep = 20, i;
	const int debug = 0;

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

int GMRFLib_density_layout_x(double *x_vec, int *len_x, GMRFLib_density_tp *density)
{
	GMRFLib_ENTER_ROUTINE;

	// this one must be increasing
	double p[] = { 0.00001, 0.0001, 0.001, 0.01, 0.025,
		0.05, 0.10, 0.15, 0.2, 0.25, 0.30,
		0.35, 0.375, 0.40, 0.425, 0.45, 0.46, 0.47, 0.475, 0.48, 0.49, 0.50, 0.51,
		0.52, 0.525, 0.53, 0.54, 0.55, 0.575, 0.60, 0.625, 0.65, 0.70,
		0.75, 0.80, 0.85, 0.9, 0.95, 0.975, 0.99, 0.999, 0.9999, 0.99999
	};

	*len_x = sizeof(p) / sizeof(double);
	if (x_vec) {
		for (int i = 0; i < *len_x; i++) {
			GMRFLib_density_Pinv(&(x_vec[i]), p[i], density);
		}
	}

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

/*

  Both functions below modified by Haavard Rue, April 2021
  
  Repo is https://github.com/thomasluu/sncdfinv/

  Copyright 2015-2016 Thomas Luu

  Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the
  License. You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

  Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS"
  BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language
  governing permissions and limitations under the License.

  plog: Computation of the Lambert W-function by Halley's Method.

  Initial guesses based on:

  D.A. Barry, J.-Y. Parlange, L. Li, H. Prommer, C.J. Cunningham, and F. Stagnitti. Analytical approximations for real values of
  the Lambert W-function. Mathematics and Computers in Simulation, 53(1):95-103, 2000.

  D.A. Barry, J.-Y. Parlange, L. Li, H. Prommer, C.J. Cunningham, and F. Stagnitti. Erratum to analytical approximations for
  real values of the Lambert W-function. Mathematics and computers in simulation, 59(6):543-543, 2002.

  GMRFLib_sn_Pinv:
  
  Based on: Luu, T; (2016) Fast and accurate parallel computation of quantile functions for random number generation. Doctoral
  thesis, UCL (University College London). http://discovery.ucl.ac.uk/1482128/
*/

double plog(double x)
{
	if (x == 0.0) {
		return 0.0;
	}

	double w0, w1;
	if (x > 0.0) {
		w0 = log(1.2 * x / log(2.4 * x / log1p(2.4 * x)));
	} else {
		double v = 1.4142135623730950488 * sqrt(1 + 2.7182818284590452354 * x);
		double N2 = 10.242640687119285146 + 1.9797586132081854940 * v;
		double N1 = 0.29289321881345247560 * (1.4142135623730950488 + N2);
		w0 = -1.0 + v * (N2 + v) / (N2 + v + N1 * v);
	}

	while (1) {
		double e = exp(w0);
		double f = w0 * e - x;
		w1 = w0 + ((f + f) * (1.0 + w0)) / (f * (2.0 + w0) - (e + e) * (1.0 + w0) * (1.0 + w0));
		if (ABS(w0 / w1 - 1.0) < 1.4901161193847656e-8) {
			break;
		}
		w0 = w1;
	}

	return w1;
}

double GMRFLib_sn_Pinv(double u, double a)
{
	double tol = 0.01;

	if (u == 0.0) {
		return -INFINITY;
	}

	if (u == 1.0) {
		return INFINITY;
	}

	/*
	 * Change of variable + special cases
	 */
	if (a == 1.0)
		u = sqrt(u);
	if (a == -1.0)
		u = sqrt(1.0 - u);
	double z = GMRFLib_Phi_inv(u);

	if (a == 0)
		return z;
	else if (a == 1)
		return z;
	else if (a == -1)
		return -z;
	if (a < 0)
		z = -z;

	double A = ABS(a);
	double right_limit = GMRFLib_erf(GMRFLib_erfc_inv(2 * tol) / A);

	/*
	 * Tails
	 */
	if (a > 0 && u > right_limit) {
		return 1.4142135623730950488 * GMRFLib_erf_inv(u);
	} else if (a < 0 && (1 - u) > right_limit) {
		return -1.4142135623730950488 * GMRFLib_erfc_inv(u);
	}

	double x = GMRFLib_Phi_inv(0.5 - 0.31830988618379067154 * atan(A));
	double expon = exp(0.5 * (-x * x));
	double errfn = 1.0;
	double efder = expon * 0.79788456080286535588 * A / errfn;

	double c0 = 0;
	double c1 = expon / errfn;
	double c2 = -expon * (efder + errfn * x) / (2 * errfn * errfn);
	double c3 = 0.16666666666666666667 *
	    expon * (3 * efder * efder + errfn * errfn * (-1 + x * x) + expon * expon + efder * (3 * errfn * x)) / (errfn * errfn * errfn);
	double c4 = -0.041666666666666666667 *
	    expon * (15 * efder * efder * efder + errfn * errfn * errfn * x * (-3 + x * x) +
		     6 * errfn * expon * expon * x + efder * efder * (18 * errfn * x) +
		     efder * (errfn * errfn * (-4 + 7 * x * x) + expon * expon * (7 - A * A))) / (errfn * errfn * errfn * errfn);
	double c5 = 0.0083333333333333333333 *
	    expon * (105 * efder * efder * efder * efder + errfn * errfn * errfn * errfn * (3 - 6 * x * x + x * x * x * x) +
		     5 * errfn * errfn * expon * expon * (-2 + 5 * x * x) + expon * expon * expon * expon * (7) +
		     15 * efder * efder * efder * (10 * errfn * x) +
		     efder * (5 * errfn * errfn * errfn * x * (-5 + 3 * x * x) + 10 * errfn * expon * expon * x * (7 - A * A)) +
		     5 * efder * efder * (3 * errfn * errfn * (-2 + 5 * x * x) + expon * expon * (-3 * (-4 + A * A)))) /
	    (errfn * errfn * errfn * errfn * errfn);

	double h = 0.75 * pow(ABS(tol / c5), 0.2);
	double left_limit = x - h;
	if (z < left_limit) {
		if (a > 0) {
			return -sqrt(2 * plog(1 / (6.2831853071795864769 * u * a)) / (1 + a * a));
		} else {
			return sqrt(2 * plog(1 / (6.2831853071795864769 * (1 - u) * ABS(a))) / (1 + a * a));
		}
	}
	h = z - x;
	double res = c0 + h * (c1 + h * (c2 + h * (c3 + h * (c4 + h * c5))));
	return (a < 0 ? -res : res);
}

double GMRFLib_sn_mode(double skew)
{
	// return the mode for a skew-normal with moments=c(0,1,skew)

	static GMRFLib_spline_tp **spline = NULL;
	int idx = 0;
	GMRFLib_CACHE_SET_ID(idx);

	if (!spline) {
#pragma omp critical (Name_f7f083055f5255ebc7e4aae6b1b8f1baa3991d18)
		{
			if (!spline) {
				spline = Calloc(GMRFLib_CACHE_LEN, GMRFLib_spline_tp *);
			}
		}
	}

	if (!spline[idx]) {
#pragma omp critical (Name_948ad03ea9be172f0a7c5a6c9a6445f21830dd19)
		{
			if (!spline[idx]) {

				// find.sn.mode <- function(skew, mode.initial = NULL) {
				// dp <- unlist(INLA:::inla.sn.reparam(moments = c(0, 1, skew)))
				// if (is.null(mode.initial)) {
				// mode.initial <- qsn(0.5, dp = dp)
				// }
				// res <- optim(mode.initial, dsn, log = TRUE, dp = dp,
				// method = "BFGS",
				// control = list(fnscale = -1, ndeps = 1e-5))
				// return(res$par)
				// }
				// ss <- seq(-0.988, 0.988, by = 0.005)
				// mm <- numeric(length(ss))

				double skews[] = {
					-0.988, -0.983, -0.978, -0.973, -0.968, -0.963, -0.958, -0.953, -0.948, -0.943, -0.938, -0.933, -0.928,
					-0.923, -0.918, -0.913, -0.908, -0.903, -0.898, -0.893, -0.888, -0.883, -0.878, -0.873, -0.868, -0.863,
					-0.858, -0.853, -0.848, -0.843, -0.838, -0.833, -0.828, -0.823, -0.818, -0.813, -0.808, -0.803, -0.798,
					-0.793, -0.788, -0.783, -0.778, -0.773, -0.768, -0.763, -0.758, -0.753, -0.748, -0.743, -0.738, -0.733,
					-0.728, -0.723, -0.718, -0.713, -0.708, -0.703, -0.698, -0.693, -0.688, -0.683, -0.678, -0.673, -0.668,
					-0.663, -0.658, -0.653, -0.648, -0.643, -0.638, -0.633, -0.628, -0.623, -0.618, -0.613, -0.608, -0.603,
					-0.598, -0.593, -0.588, -0.583, -0.578, -0.573, -0.568, -0.563, -0.558, -0.553, -0.548, -0.543, -0.538,
					-0.533, -0.528, -0.523, -0.518, -0.513, -0.508, -0.503, -0.498, -0.493, -0.488, -0.483, -0.478, -0.473,
					-0.468, -0.463, -0.458, -0.453, -0.448, -0.443, -0.438, -0.433, -0.428, -0.423, -0.418, -0.413, -0.408,
					-0.403, -0.398, -0.393, -0.388, -0.383, -0.378, -0.373, -0.368, -0.363, -0.358, -0.353, -0.348, -0.343,
					-0.338, -0.333, -0.328, -0.323, -0.318, -0.313, -0.308, -0.303, -0.298, -0.293, -0.288, -0.283, -0.278,
					-0.273, -0.268, -0.263, -0.258, -0.253, -0.248, -0.243, -0.238, -0.233, -0.228, -0.223, -0.218, -0.213,
					-0.208, -0.203, -0.198, -0.193, -0.188, -0.183, -0.178, -0.173, -0.168, -0.163, -0.158, -0.153, -0.148,
					-0.143, -0.138, -0.133, -0.128, -0.123, -0.118, -0.113, -0.108, -0.103, -0.098, -0.093, -0.088, -0.083,
					-0.078, -0.073, -0.068, -0.063, -0.058, -0.053, -0.048, -0.043, -0.038, -0.033, -0.028, -0.023, -0.018,
					-0.013, -0.008, -0.003, 0.002, 0.007, 0.012, 0.017, 0.022, 0.027, 0.032, 0.037, 0.042, 0.047, 0.052,
					0.057, 0.062, 0.067, 0.072, 0.077, 0.082, 0.087, 0.092, 0.097, 0.102, 0.107, 0.112, 0.117, 0.122, 0.127,
					0.132, 0.137, 0.142, 0.147, 0.152, 0.157, 0.162, 0.167, 0.172, 0.177, 0.182, 0.187, 0.192, 0.197, 0.202,
					0.207, 0.212, 0.217, 0.222, 0.227, 0.232, 0.237, 0.242, 0.247, 0.252, 0.257, 0.262, 0.267, 0.272, 0.277,
					0.282, 0.287, 0.292, 0.297, 0.302, 0.307, 0.312, 0.317, 0.322, 0.327, 0.332, 0.337, 0.342, 0.347, 0.352,
					0.357, 0.362, 0.367, 0.372, 0.377, 0.382, 0.387, 0.392, 0.397, 0.402, 0.407, 0.412, 0.417, 0.422, 0.427,
					0.432, 0.437, 0.442, 0.447, 0.452, 0.457, 0.462, 0.467, 0.472, 0.477, 0.482, 0.487, 0.492, 0.497, 0.502,
					0.507, 0.512, 0.517, 0.522, 0.527, 0.532, 0.537, 0.542, 0.547, 0.552, 0.557, 0.562, 0.567, 0.572, 0.577,
					0.582, 0.587, 0.592, 0.597, 0.602, 0.607, 0.612, 0.617, 0.622, 0.627, 0.632, 0.637, 0.642, 0.647, 0.652,
					0.657, 0.662, 0.667, 0.672, 0.677, 0.682, 0.687, 0.692, 0.697, 0.702, 0.707, 0.712, 0.717, 0.722, 0.727,
					0.732, 0.737, 0.742, 0.747, 0.752, 0.757, 0.762, 0.767, 0.772, 0.777, 0.782, 0.787, 0.792, 0.797, 0.802,
					0.807, 0.812, 0.817, 0.822, 0.827, 0.832, 0.837, 0.842, 0.847, 0.852, 0.857, 0.862, 0.867, 0.872, 0.877,
					0.882, 0.887, 0.892, 0.897, 0.902, 0.907, 0.912, 0.917, 0.922, 0.927, 0.932, 0.937, 0.942, 0.947, 0.952,
					0.957, 0.962, 0.967, 0.972, 0.977, 0.982, 0.987, 0.988
				};
				double modes[] = {
					1.114585842, 1.065449627, 1.027881735, 0.9967178902, 0.9698704599, 0.9460064999, 0.924449043, 0.9047137376,
					0.8864606977, 0.8694413827, 0.8534679246, 0.838394446, 0.8241051319, 0.8105063076, 0.7975209973,
					0.7850850806, 0.77314451, 0.7616532528, 0.7505669292, 0.7398569764, 0.7294983525, 0.7194584508,
					0.70971416, 0.7002448982, 0.691032261, 0.6820597188, 0.6733123661, 0.6647791261, 0.6564424868,
					0.6482942012, 0.6403240413, 0.6325226312, 0.6248813534, 0.6173922672, 0.6100480376, 0.6028418732,
					0.595767473, 0.5888189779, 0.5819909307, 0.5752782376, 0.568676137, 0.5621801698, 0.5557861541,
					0.5494901601, 0.5432884944, 0.5371776746, 0.5311544183, 0.5252156249, 0.519358363, 0.5135798576,
					0.5078774801, 0.5022487363, 0.4966912553, 0.491202792, 0.4857811891, 0.4804244217, 0.4751305445,
					0.4698976905, 0.4647241047, 0.4596080836, 0.4545480101, 0.4495423455, 0.4445895978, 0.439688347,
					0.4348372386, 0.4300349576, 0.4252802518, 0.4205719153, 0.4159087892, 0.4112897584, 0.4067137503,
					0.4021797296, 0.3976867142, 0.3932337235, 0.3888198353, 0.3844441705, 0.3801058661, 0.3758040772,
					0.3715380058, 0.3673068744, 0.3631099308, 0.3589464475, 0.3548157203, 0.3507170671, 0.3466498275,
					0.3426133613, 0.3386070477, 0.3346302812, 0.3306824894, 0.3267630941, 0.3228715489, 0.3190073162,
					0.3151698844, 0.3113587467, 0.307573406, 0.3038134057, 0.3000782529, 0.2963675003, 0.2926807374,
					0.2890175074, 0.2853774035, 0.2817600117, 0.2781649603, 0.27459184, 0.2710402762, 0.2675099159,
					0.2640004004, 0.2605113685, 0.2570424857, 0.2535934195, 0.2501638451, 0.2467534452, 0.2433619099,
					0.2399889362, 0.2366342279, 0.2332974954, 0.2299784556, 0.2266768315, 0.2233923518, 0.2201247514,
					0.2168737704, 0.2136391546, 0.2104206522, 0.2072180245, 0.2040310321, 0.2008594358, 0.1977030062,
					0.1945615218, 0.1914347725, 0.1883225074, 0.1852245535, 0.182140665, 0.1790706684, 0.1760143481,
					0.1729715153, 0.1699419538, 0.1669255036, 0.1639219537, 0.1609311221, 0.1579528439, 0.1549869184,
					0.1520331932, 0.1490914996, 0.1461616047, 0.143243418, 0.1403367373, 0.1374414108, 0.1345572957,
					0.1316841929, 0.1288219537, 0.1259704665, 0.1231295378, 0.1202990176, 0.1174787906, 0.1146687001,
					0.1118685727, 0.1090783003, 0.1062977307, 0.1035267231, 0.1007651394, 0.09801284307, 0.09526969859,
					0.09253557171, 0.08981032931, 0.08709383929, 0.08438594199, 0.0816865599, 0.07899554989, 0.07631276676,
					0.07363811092, 0.07097140669, 0.06831253016, 0.06566141674, 0.06301785793, 0.06038178142, 0.05775305322,
					0.05513154798, 0.05251709206, 0.04990961362, 0.04730895685, 0.04471499446, 0.04212760487, 0.0395466322,
					0.0369719521, 0.0344034263, 0.03184090653, 0.02928426998, 0.02673334651, 0.02418798502, 0.02164803605,
					0.01911335275, 0.01658371003, 0.01405894938, 0.01153887215, 0.009023262093, 0.006511879629,
					0.004004449754, 0.001500681424, -0.001000307827, -0.003503410089, -0.006010082189, -0.008520640322,
					-0.01103537855, -0.01355453761, -0.01607835198, -0.01860699845, -0.02114067591, -0.02367956564,
					-0.02622382112, -0.02877361694, -0.03132911105, -0.03389043752, -0.03645773759, -0.03903118297,
					-0.0416108849, -0.0441969904, -0.04678962705, -0.04938893724, -0.05199504357, -0.0546080752,
					-0.05722816835, -0.05985544782, -0.06249004493, -0.06513207628, -0.0677816799, -0.07043897937,
					-0.07310410191, -0.0757771798, -0.07845833153, -0.08114768987, -0.0838453847, -0.08655154333,
					-0.08926629573, -0.09198977692, -0.09472211286, -0.09746344124, -0.1002138963, -0.1029736078,
					-0.1057427193, -0.1085213634, -0.111309681, -0.1141078727, -0.1169159751, -0.1197341571, -0.1225625954,
					-0.1254014237, -0.1282508141, -0.1311108642, -0.1339817971, -0.1368636975, -0.1397567726, -0.1426611666,
					-0.1455770415, -0.1485045576, -0.1514438942, -0.1543952083, -0.1573586763, -0.160334473, -0.1633227763,
					-0.1663237667, -0.1693376219, -0.1723645409, -0.1754047079, -0.1784583238, -0.1815255712, -0.1846066475,
					-0.1877017759, -0.1908111674, -0.1939350109, -0.197073523, -0.2002269548, -0.2033954854, -0.2065793824,
					-0.2097788798, -0.2129941731, -0.2162255601, -0.2194732486, -0.2227374896, -0.2260185723, -0.2293167426,
					-0.2326322762, -0.2359654522, -0.2393165423, -0.2426858395, -0.2460736384, -0.2494802361, -0.2529059551,
					-0.2563510972, -0.2598159903, -0.2633009654, -0.2668063553, -0.270332526, -0.2738798062, -0.2774485925,
					-0.2810392259, -0.284652121, -0.2882876404, -0.2919462095, -0.2956282505, -0.2993341523, -0.3030643892,
					-0.3068194044, -0.310599629, -0.3144055821, -0.3182377011, -0.322096533, -0.3259825867, -0.3298963543,
					-0.3338384406, -0.3378093445, -0.3418097033, -0.3458400898, -0.3499011239, -0.3539934418, -0.3581177001,
					-0.3622745762, -0.3664647698, -0.3706890039, -0.3749480255, -0.3792426069, -0.3835735397, -0.3879416727,
					-0.3923478397, -0.3967929349, -0.4012778767, -0.405803608, -0.410371146, -0.414981487, -0.4196357055,
					-0.4243349078, -0.4290802439, -0.4338729099, -0.4387141507, -0.4436052629, -0.4485475972, -0.4535425625,
					-0.4585916284, -0.4636963296, -0.4688582697, -0.4740791254, -0.4793606512, -0.4847046847, -0.4901131518,
					-0.4955880726, -0.5011315682, -0.5067458676, -0.5124333151, -0.5181963792, -0.5240376614, -0.5299599061,
					-0.5359660124, -0.5420590452, -0.5482422497, -0.5545190652, -0.5608931421, -0.5673683598, -0.5739488472,
					-0.5806390049, -0.5874435304, -0.5943674466, -0.6014161335, -0.6085958898, -0.6159119207, -0.6233713921,
					-0.6309815213, -0.6387501172, -0.6466856486, -0.6547973232, -0.6630951796, -0.6715901933, -0.6802954699,
					-0.6892225117, -0.6983867413, -0.7078044105, -0.7174936184, -0.7274746142, -0.7377701665, -0.7484060197,
					-0.7593983179, -0.7708029105, -0.7826528429, -0.7949843753, -0.8078533955, -0.8213216825, -0.8354631758,
					-0.8503675374, -0.8661451515, -0.8829343369, -0.9009120832, -0.9203106208, -0.941496803, -0.9648408635,
					-0.9910483779, -1.021200258, -1.05723934, -1.103375097, -1.114585842
				};

				spline[idx] = GMRFLib_spline_create(skews, modes, (int) sizeof(skews) / sizeof(double));
			}
		}
	}

	skew = TRUNCATE(skew, -GMRFLib_SN_SKEWMAX, GMRFLib_SN_SKEWMAX);
	return (GMRFLib_spline_eval(skew, spline[idx]));
}

double GMRFLib_sn_d3_to_skew(double d3)
{
	// find the skewness for a given third order derivative at the model, with mean=0 and var=1.

	static GMRFLib_spline_tp **spline = NULL;
	int idx = 0;
	GMRFLib_CACHE_SET_ID(idx);

	if (!spline) {
#pragma omp critical (Name_06501c73f0089b8702336f89a7e7c85e10465cf1)
		{
			if (!spline) {
				spline = Calloc(GMRFLib_CACHE_LEN, GMRFLib_spline_tp *);
			}
		}
	}

	if (!spline[idx]) {
#pragma omp critical (Name_0b84ff96beeead8255ba5c20d89eda3f5360f93c)
		{
			if (!spline[idx]) {
				// both skew and d3 have been POWER13 transformed, as this gives a much better function to interpolate
				double skew3s[] = {
					-0.9959838925, -0.9943009155, -0.9926122218, -0.9909177627, -0.9892174886, -0.9875113495, -0.9857992945,
					-0.9840812721, -0.9823572299, -0.9806271149, -0.9788908735, -0.977148451, -0.9753997922, -0.973644841,
					-0.9718835404, -0.9701158327, -0.9683416593, -0.9665609608, -0.9647736769, -0.9629797462, -0.9611791067,
					-0.9593716954, -0.957557448, -0.9557362998, -0.9539081846, -0.9520730354, -0.9502307842, -0.9483813619,
					-0.9465246982, -0.944660722, -0.9427893606, -0.9409105407, -0.9390241873, -0.9371302245, -0.9352285752,
					-0.9333191608, -0.9314019016, -0.9294767164, -0.927543523, -0.9256022375, -0.9236527746, -0.9216950477,
					-0.9197289687, -0.9177544479, -0.915771394, -0.9137797144, -0.9117793146, -0.9097700985, -0.9077519683,
					-0.9057248245, -0.9036885658, -0.901643089, -0.8995882891, -0.897524059, -0.8954502899, -0.8933668708,
					-0.8912736887, -0.8891706283, -0.8870575722, -0.884934401, -0.8828009925, -0.8806572225, -0.8785029644,
					-0.8763380887, -0.8741624639, -0.8719759553, -0.869778426, -0.8675697359, -0.8653497422, -0.8631182992,
					-0.8608752582, -0.8586204672, -0.8563537711, -0.8540750116, -0.8517840269, -0.8494806516, -0.8471647169,
					-0.84483605, -0.8424944747, -0.8401398104, -0.8377718728, -0.8353904732, -0.8329954186, -0.8305865115,
					-0.8281635499, -0.825726327, -0.8232746311, -0.8208082453, -0.8183269477, -0.8158305107, -0.8133187014,
					-0.8107912808, -0.8082480041, -0.8056886203, -0.8031128718, -0.8005204946, -0.7979112176, -0.7952847628,
					-0.7926408444, -0.7899791695, -0.7872994366, -0.7846013365, -0.7818845511, -0.7791487536, -0.7763936077,
					-0.7736187677, -0.7708238778, -0.768008572, -0.7651724731, -0.7623151931, -0.7594363318, -0.7565354772,
					-0.7536122043, -0.750666075, -0.747696637, -0.7447034238, -0.7416859539, -0.7386437295, -0.7355762368,
					-0.7324829445, -0.729363303, -0.726216744, -0.7230426793, -0.7198404996, -0.7166095742, -0.713349249,
					-0.7100588459, -0.7067376615, -0.7033849656, -0.7, -0.6965819768, -0.6931300768, -0.6896434481,
					-0.6861212036, -0.6825624197, -0.6789661336, -0.6753313417, -0.6716569962, -0.6679420032, -0.6641852195,
					-0.6603854498, -0.6565414427, -0.6526518879, -0.6487154117, -0.6447305727, -0.6406958577, -0.636609676,
					-0.6324703543, -0.6282761305, -0.6240251469, -0.6197154435, -0.6153449494, -0.6109114744, -0.6064126995,
					-0.6018461655, -0.597209262, -0.5924992137, -0.5877130659, -0.5828476683, -0.5778996565, -0.5728654316,
					-0.5677411371, -0.5625226328, -0.5572054656, -0.5517848353, -0.5462555571, -0.5406120176, -0.5348481241,
					-0.5289572473, -0.5229321532, -0.5167649252, -0.5104468722, -0.50396842, -0.4973189833, -0.4904868132,
					-0.4834588127, -0.4762203156, -0.4687548148, -0.4610436292, -0.4530654896, -0.4447960181, -0.4362070671,
					-0.4272658682, -0.4179339196, -0.4081655102, -0.3979057208, -0.3870876641, -0.3756285754, -0.3634241186,
					-0.350339806, -0.3361975407, -0.320753433, -0.3036588972, -0.284386698, -0.2620741394, -0.2351334688,
					-0.2, -0.144224957, 0.125992105, 0.1912931183, 0.2289428485, 0.2571281591, 0.2802039331, 0.3,
					0.3174802104, 0.3332221852, 0.3476026645, 0.360882608, 0.3732511157, 0.3848501131, 0.395789161,
					0.40615481, 0.4160167646, 0.4254320865, 0.4344481486, 0.4431047622, 0.4514357435, 0.4594700892,
					0.4672328728, 0.4747459399, 0.4820284528, 0.4890973247, 0.4959675664, 0.5026525695, 0.509164337,
					0.5155136735, 0.5217103446, 0.5277632088, 0.5336803297, 0.5394690712, 0.5451361778, 0.5506878446,
					0.5561297767, 0.5614672408, 0.5667051108, 0.5718479065, 0.5768998281, 0.5818647867, 0.5867464308,
					0.59154817, 0.5962731958, 0.6009245007, 0.6055048947, 0.61001702, 0.6144633651, 0.6188462762,
					0.6231679684, 0.6274305357, 0.6316359598, 0.635786118, 0.639882791, 0.6439276696, 0.6479223603,
					0.6518683915, 0.6557672186, 0.6596202284, 0.6634287437, 0.6671940272, 0.6709172852, 0.6745996712,
					0.6782422886, 0.6818461941, 0.6854124002, 0.6889418775, 0.6924355573, 0.6958943337, 0.6993190657,
					0.7027105788, 0.7060696671, 0.7093970945, 0.7126935967, 0.7159598825, 0.7191966348, 0.7224045124,
					0.7255841507, 0.7287361631, 0.731861142, 0.7349596597, 0.7380322692, 0.7410795055, 0.7441018861,
					0.7470999115, 0.7500740668, 0.7530248212, 0.7559526299, 0.7588579338, 0.7617411603, 0.7646027242,
					0.7674430279, 0.7702624618, 0.7730614053, 0.7758402264, 0.7785992832, 0.7813389232, 0.7840594846,
					0.786761296, 0.7894446773, 0.7921099395, 0.7947573855, 0.7973873099, 0.8, 0.8025957353, 0.8051747881,
					0.8077374241, 0.8102839019, 0.8128144739, 0.8153293862, 0.8178288788, 0.8203131859, 0.8227825361,
					0.8252371525, 0.8276772529, 0.8301030501, 0.8325147517, 0.8349125609, 0.837296676, 0.8396672908,
					0.8420245948, 0.8443687734, 0.8467000076, 0.8490184749, 0.8513243484, 0.853617798, 0.8558989894,
					0.8581680854, 0.8604252449, 0.8626706237, 0.8649043743, 0.867126646, 0.8693375853, 0.8715373356,
					0.8737260372, 0.875903828, 0.8780708428, 0.8802272141, 0.8823730714, 0.8845085422, 0.8866337511,
					0.8887488205, 0.8908538706, 0.8929490191, 0.8950343817, 0.8971100718, 0.8991762009, 0.9012328783,
					0.9032802112, 0.9053183053, 0.9073472639, 0.9093671888, 0.9113781798, 0.9133803351, 0.9153737512,
					0.9173585227, 0.9193347428, 0.9213025029, 0.9232618931, 0.9252130018, 0.927155916, 0.9290907211,
					0.9310175012, 0.9329363391, 0.934847316, 0.9367505121, 0.938646006, 0.9405338751, 0.9424141957,
					0.9442870428, 0.9461524903, 0.9480106107, 0.9498614756, 0.9517051555, 0.9535417196, 0.9553712362,
					0.9571937726, 0.9590093948, 0.9608181683, 0.962620157, 0.9644154244, 0.9662040328, 0.9679860436,
					0.9697615172, 0.9715305133, 0.9732930906, 0.9750493072, 0.9767992199, 0.9785428852, 0.9802803585,
					0.9820116944, 0.9837369469, 0.9854561691, 0.9871694135, 0.9888767317, 0.9905781747, 0.9922737928,
					0.9939636356, 0.9956477521
				};

				double d33s[] = {
					-4.90997749, -4.235163358, -3.831583864, -3.547746883, -3.332973363, -3.160324901, -3.017221905,
					-2.895548811,
					-2.790066112, -2.697210047, -2.614453695, -2.539944085, -2.472284424, -2.410397163, -2.353434582,
					-2.300718341,
					-2.251697909, -2.205920556, -2.16299577, -2.122625515, -2.0845534, -2.048532342, -2.014365658, -1.981882874,
					-1.950935145, -1.921391799, -1.893137844, -1.866075308, -1.840104336, -1.815149269, -1.791137009,
					-1.76800217,
					-1.7456857, -1.724133961, -1.70329812, -1.683134081, -1.663601211, -1.644662346, -1.62628348, -1.608433297,
					-1.591082775, -1.574205395, -1.557776378, -1.541772604, -1.526173054, -1.510957782, -1.49610839,
					-1.481607671,
					-1.467439178, -1.453588179, -1.440040309, -1.426782431, -1.413802075, -1.401087196, -1.388627006,
					-1.376411093,
					-1.364429234, -1.352672575, -1.341131967, -1.329799057, -1.318666164, -1.307725428, -1.296970153,
					-1.286392916,
					-1.275987752, -1.265748278, -1.255668558, -1.24574318, -1.23596678, -1.226334184, -1.216840275,
					-1.207480843,
					-1.19825095, -1.189146434, -1.180163165, -1.171297341, -1.162545175, -1.153902491, -1.14536661,
					-1.136933328,
					-1.128600033, -1.120363508, -1.112220068, -1.104167381, -1.096202572, -1.088322712, -1.080525264,
					-1.072807717,
					-1.065167579, -1.057602307, -1.050110115, -1.042687974, -1.035334021, -1.028046025, -1.020822233,
					-1.013660521,
					-1.006558526, -0.9995149145, -0.9925276209, -0.9855944034, -0.9787143046, -0.9718845834, -0.9651047846,
					-0.9583721229, -0.9516852592, -0.94504301, -0.9384434168, -0.9318849065, -0.9253665945, -0.9188856964,
					-0.9124419599, -0.9060334706, -0.8996586778, -0.8933171883, -0.8870056476, -0.8807243502, -0.8744710544,
					-0.8682453127, -0.8620444646, -0.8558679377, -0.8497146348, -0.8435826853, -0.837470988, -0.8313782538,
					-0.8253036296, -0.8192445368, -0.813200221, -0.8071700339, -0.8011518856, -0.7951442923, -0.789146817,
					-0.7831572906, -0.7771742019, -0.7711968435, -0.7652225055, -0.7592515183, -0.753280373, -0.7473097709,
					-0.7413366722, -0.7353597784, -0.729376837, -0.723387945, -0.7173892731, -0.7113808689, -0.7053596713,
					-0.6993248627, -0.6932734026, -0.6872039892, -0.6811137176, -0.6750019388, -0.6688656108, -0.6627027124,
					-0.6565092451, -0.6502840527, -0.6440244769, -0.6377277393, -0.6313896271, -0.6250081411, -0.6185798621,
					-0.6121007705, -0.6055680735, -0.5989765651, -0.59232087, -0.5856003558, -0.5788063504, -0.571934422,
					-0.5649804903,
					-0.5579380628, -0.5507979586, -0.5435554761, -0.536201069, -0.5287278908, -0.5211252141, -0.5133820422,
					-0.5054876945, -0.497428898, -0.4891903532, -0.4807558969, -0.4721099083, -0.4632276881, -0.4540847412,
					-0.4446564596, -0.4349092489, -0.4248027143, -0.4142941186, -0.4033294175, -0.3918411722, -0.3797456696,
					-0.3669439254, -0.3533001657, -0.3386323889, -0.3227044594, -0.3051719475, -0.2855042223, -0.2628471099,
					-0.2356039192, -0.2002255401, -0.1442893202, 0.1260279293, 0.1914846493, 0.229357556, 0.2578358508,
					0.2812476452,
					0.3014329649, 0.3193424917, 0.3355609642, 0.350455108, 0.3642880366, 0.3772456037, 0.3894740454,
					0.4010753747,
					0.4121392262, 0.422735397, 0.4329188381, 0.4427324699, 0.4522233133, 0.4614203074, 0.470353076,
					0.4790447723,
					0.4875204088, 0.4957964526, 0.5038892848, 0.5118153296, 0.5195881717, 0.5272191117, 0.5347160472,
					0.5420933562,
					0.549357132, 0.5565168535, 0.5635795583, 0.570550925, 0.5774381407, 0.5842468045, 0.5909828774,
					0.5976501001,
					0.6042533692, 0.6107985143, 0.6172891075, 0.6237270514, 0.6301177642, 0.636464308, 0.6427686185,
					0.6490356785,
					0.6552668416, 0.6614659993, 0.6676351308, 0.6737766914, 0.6798938983, 0.6859874074, 0.6920608253,
					0.6981159443,
					0.704154038, 0.7101774175, 0.7161887287, 0.7221887265, 0.7281793524, 0.7341630542, 0.7401413104,
					0.7461150308,
					0.7520862394, 0.7580573248, 0.7640283143, 0.7700013424, 0.7759780151, 0.7819600214, 0.7879483138,
					0.7939441399,
					0.7999495046, 0.805965797, 0.811993654, 0.8180346895, 0.8240903013, 0.8301621473, 0.8362510257,
					0.8423585585,
					0.8484860827, 0.8546356553, 0.8608073444, 0.867002689, 0.8732239223, 0.8794716125, 0.8857473446,
					0.8920523675,
					0.8983878707, 0.9047559478, 0.9111575416, 0.9175941067, 0.9240673444, 0.930578244, 0.9371283888,
					0.9437198638,
					0.9503535033, 0.9570310466, 0.9637542057, 0.9705248702, 0.9773441424, 0.9842142727, 0.9911366703,
					0.9981130466,
					1.005145479, 1.012235327, 1.019384755, 1.026596309, 1.033871066, 1.041211662, 1.048619711, 1.056098059,
					1.063648624,
					1.071273658, 1.078975473, 1.086756775, 1.094619681, 1.102567448, 1.110602178, 1.118727158, 1.126945175,
					1.135258946,
					1.14367172, 1.152187143, 1.160807775, 1.169537802, 1.17838077, 1.187340262, 1.196420021, 1.205624547,
					1.21495778,
					1.224424229, 1.234029014, 1.243776125, 1.25367142, 1.263719654, 1.273926791, 1.284298355, 1.294840451,
					1.305559848,
					1.316462921, 1.327556779, 1.338849174, 1.350347444, 1.362060275, 1.373996231, 1.386164561, 1.398575122,
					1.411238192,
					1.42416458, 1.437365931, 1.450854732, 1.464643906, 1.478747605, 1.493180797, 1.507959177, 1.523099882,
					1.538621003,
					1.55454213, 1.570884322, 1.587670199, 1.604923932, 1.622671928, 1.64094257, 1.659766579, 1.679177759,
					1.699212761,
					1.719910761, 1.741315544, 1.763474802, 1.786441191, 1.810272718, 1.835033312, 1.86079436, 1.887635175,
					1.915646392,
					1.944924709, 1.97558397, 2.007751546, 2.041572744, 2.077214182, 2.114867631, 2.154755773, 2.197098356,
					2.24226701,
					2.29061247, 2.342542229, 2.398596352, 2.459423563, 2.525831304, 2.598841699, 2.679772629, 2.770360879,
					2.872955193,
					2.990831769, 3.129309187, 3.294994665, 3.500289878, 3.766951589, 4.140012105, 4.73787115
				};

				int n = (int) sizeof(skew3s) / sizeof(double);
				spline[idx] = GMRFLib_spline_create(d33s, skew3s, n);
			}
		}
	}

#define POWER13(x_) ((x_) <  0.0 ? - pow(ABS(x_), 1.0/3.0) : pow(ABS(x_), 1.0/3.0))
#define iPOWER13(x_) gsl_pow_3(x_)

	double skew;

	d3 = POWER13(d3);
	skew = GMRFLib_spline_eval(d3, spline[idx]);
	skew = iPOWER13(skew);
	skew = TRUNCATE(skew, -GMRFLib_SN_SKEWMAX, GMRFLib_SN_SKEWMAX);

#undef POWER13
#undef iPOWER13

	return (skew);
}
