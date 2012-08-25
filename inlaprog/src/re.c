/* re.c
 * 
 * Copyright (C) 2012 Havard Rue
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
 */
#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = HGVERSION;

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

#include "re.h"

#define Pq(_arg) (exp(0.25)/sqrt(8.0*M_PI) * (gsl_sf_bessel_Knu(((_arg)+1.0)/2.0, 0.25) + gsl_sf_bessel_Knu(ABS(((_arg)-1.0)/2.0), 0.25)))
#define M1(_epsilon, _delta) (sinh((_epsilon)/(_delta))*Pq(1.0/(_delta)))
#define M2(_epsilon, _delta) (0.5*(cosh(2.0*(_epsilon)/(_delta))*Pq(2.0/(_delta)) -1.0))
#define M3(_epsilon, _delta) (0.25*(sinh(3.0*(_epsilon)/(_delta))*Pq(3.0/(_delta)) - \
				    3.0*sinh((_epsilon)/(_delta))*Pq(1.0/(_delta))))
#define M4(_epsilon, _delta) (0.125*(cosh(4.0*(_epsilon)/(_delta))*Pq(4.0/(_delta)) - \
				     4.0*cosh(2.0*(_epsilon)/(_delta))*Pq(2.0/(_delta)) + 3.0))

int re_valid_skew_kurt(double *dist, double skew, double kurt)
{
#define KURT_LIMIT(s) (2.15 + SQR((s)/0.8))
	int retval;

	retval = (kurt > KURT_LIMIT(skew) ? GMRFLib_TRUE : GMRFLib_FALSE);
	if (dist) {
		if (retval == GMRFLib_FALSE) {
			*dist = ABS(kurt - KURT_LIMIT(skew));
		} else {
			*dist = 0.0;
		}
	}

#undef KURT_LIMIT
	return retval;
}
int re_shash_skew_kurt(double *skew, double *kurt, double epsilon, double delta)
{
	double m1 = M1(epsilon, delta);
	double m2 = M2(epsilon, delta);
	double m3 = M3(epsilon, delta);
	double m4 = M4(epsilon, delta);

	if (skew) {
		*skew = (m3 - 3.0 * m1 * m2 + 3.0 * gsl_pow_3(m1) - gsl_pow_3(m1)) / gsl_pow_3(sqrt(m2 - SQR(m1)));
	}
	if (kurt) {
		*kurt = (m4 - 4.0 * m1 * m3 + 6.0 * SQR(m1) * m2 - 4.0 * gsl_pow_4(m1) + gsl_pow_4(m1)) / gsl_pow_4(sqrt(m2 - SQR(m1)));
	}

	return GMRFLib_SUCCESS;
}


int re_shash_fit_parameters(re_shash_param_tp * param, double *mean, double *prec, double *skew, double *kurt)
{
	/*
	 * for given mean, prec, skew and kurt, fit the parameters in the shash-model. Either none or both of 'skew' and 'kurt' can be NULL. 
	 */

	int use_lookup = 1;

	static double **pin = NULL;
	static double **pout = NULL;
	static int *perr = NULL;

	assert((!skew && !kurt) || (skew && kurt));	       /* either non or both */

	if (use_lookup) {
		if (!pin && !pout && !perr) {
#pragma omp critical
			{
				if (!pin && !pout && !perr) {
					pin = Calloc(GMRFLib_MAX_THREADS, double *);
					pout = Calloc(GMRFLib_MAX_THREADS, double *);
					perr = Calloc(GMRFLib_MAX_THREADS, int);

					int i;
					for(i=0; i<GMRFLib_MAX_THREADS; i++){
						pin[i] = Calloc(2, double);
						pout[i] = Calloc(2, double);
					}
				}
			}
		}
	}

	double npin[2], pout_tmp[2];
	int thread_id, err = 0, debug = 0;

	thread_id = omp_get_thread_num();
	npin[0] = (skew ? *skew : 0.0);
	npin[1] = (kurt ? *kurt : 3.0);

	if (use_lookup && memcmp(npin, pin[thread_id], sizeof(npin)) == 0) {
		/*
		 * we have already computed these values
		 */
		memcpy(pout_tmp, pout[thread_id], sizeof(pout_tmp));
	} else if (!skew && !kurt) {
		pout_tmp[0] = 0.0;				       /* epsilon */
		pout_tmp[1] = 1.0;				       /* delta */
	} else {
		double epsilon = 0.0;
		double delta = 1.0;

		gsl_vector *x = gsl_vector_alloc(2);
		gsl_vector_set(x, 0, 0.0);		       /* initial values */
		gsl_vector_set(x, 1, 1.0);

		double target[2];
		target[0] = *skew;
		target[1] = *kurt;

		const gsl_multifit_fdfsolver_type *T;
		gsl_multifit_fdfsolver *s;
		int n = 2;
		int p = 2;
		int status;

		gsl_multifit_function_fdf f;
		f.f = re_shash_f;
		f.df = re_shash_df;
		f.fdf = re_shash_fdf;
		f.n = n;
		f.p = p;
		f.params = (void *) target;


		T = gsl_multifit_fdfsolver_lmsder;
		s = gsl_multifit_fdfsolver_alloc(T, n, p);
		gsl_multifit_fdfsolver_set(s, &f, (const gsl_vector *) x);

		int iter = 0, iter_max = 1000;
		double eps = 1e-4;
		do {
			iter++;
			status = gsl_multifit_fdfsolver_iterate(s);
			if (status) {
				break;
			}
			if (debug) {
				printf("status = %s\n", gsl_strerror(status));
				printf("iter %1d: x %g %g |f| = %g\n", iter, gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1), gsl_blas_dnrm2(s->f));
			}
			status = gsl_multifit_test_delta(s->dx, s->x, eps, eps);
		}
		while (status == GSL_CONTINUE && iter < iter_max);

		err = (gsl_blas_dnrm2(s->f) > 1e-6);
		if (err) {
			if (re_valid_skew_kurt(NULL, target[0], target[1]) == GMRFLib_TRUE){
				fprintf(stderr, "SHASH fail to fit target skew=%g kurt=%g err=%g, but VALID FAIL!\n",
					target[0], target[1], gsl_blas_dnrm2(s->f));
			}
		}

		pout_tmp[0] = gsl_vector_get(s->x, 0);
		pout_tmp[1] = gsl_vector_get(s->x, 1);

		gsl_vector_free(x);
		gsl_multifit_fdfsolver_free(s);
	}

	double m1, m2, mu, stdev, ss;
	
	m1 = M1(pout_tmp[0], pout_tmp[1]);
	m2 = M2(pout_tmp[0], pout_tmp[1]);
	ss = sqrt(m2 - SQR(m1));
	mu = *mean - m1 / sqrt(*prec) / ss;
	stdev = 1.0 / sqrt(*prec) / ss;

	param->mu = mu;
	param->stdev = stdev;
	param->epsilon = pout_tmp[0];
	param->delta = pout_tmp[1];

	if (use_lookup){
		perr[thread_id] = err;
		memcpy(pin[thread_id], npin, sizeof(npin));
		memcpy(pout[thread_id], pout_tmp, sizeof(pout_tmp));
	}
		
	return (err ? !GMRFLib_SUCCESS : GMRFLib_SUCCESS);
}
int re_shash_f(const gsl_vector * x, void *data, gsl_vector * f)
{
	double epsilon, delta, skew, sskew, kurt, kkurt, *ddata = (double *) data;

	epsilon = gsl_vector_get(x, 0);
	delta = gsl_vector_get(x, 1);
	sskew = ddata[0];
	kkurt = ddata[1];

	re_shash_skew_kurt(&skew, &kurt, epsilon, delta);
	gsl_vector_set(f, 0, skew - sskew);
	gsl_vector_set(f, 1, kurt - kkurt);

	return GSL_SUCCESS;
}
int re_shash_df(const gsl_vector * x, void *data, gsl_matrix * J)
{
	double h = GMRFLib_eps(1. / 3.), epsilon, delta, skew_h, skew_mh, kurt_h, kurt_mh; 

	epsilon = gsl_vector_get(x, 0);
	delta = gsl_vector_get(x, 1);

	re_shash_skew_kurt(&skew_h, &kurt_h, epsilon + h, delta);
	re_shash_skew_kurt(&skew_mh, &kurt_mh, epsilon - h, delta);
	gsl_matrix_set(J, 0, 0, (skew_h - skew_mh) / 2.0 / h);
	gsl_matrix_set(J, 1, 0, (kurt_h - kurt_mh) / 2.0 / h);

	re_shash_skew_kurt(&skew_h, &kurt_h, epsilon, delta + h);
	re_shash_skew_kurt(&skew_mh, &kurt_mh, epsilon, delta - h);
	gsl_matrix_set(J, 0, 1, (skew_h - skew_mh) / 2.0 / h);
	gsl_matrix_set(J, 1, 1, (kurt_h - kurt_mh) / 2.0 / h);

	return GSL_SUCCESS;
}
int re_shash_fdf(const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J)
{
	re_shash_f(x, data, f);
	re_shash_df(x, data, J);

	return GSL_SUCCESS;
}
