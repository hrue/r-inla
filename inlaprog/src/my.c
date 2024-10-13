
/* my.c
 * 
 * Copyright (C) 2016-2024 Havard Rue
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

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "inla.h"
#include "my.h"
#include "my-fix.h"
#include "GMRFLib/GMRFLibP.h"
#include "GMRFLib/GMRFLib.h"

int my_file_exists(const char *filename)
{
	struct stat sb;
	return ((stat(filename, &sb) == 0 && S_ISREG(sb.st_mode)) ? INLA_OK : !INLA_OK);
}

int my_dir_exists(const char *dirname)
{
	struct stat sb;
	return ((stat(dirname, &sb) == 0 && S_ISDIR(sb.st_mode)) ? INLA_OK : !INLA_OK);
}

int my_setenv(char *str, int prefix)
{
	/*
	 * set a variable in the enviroment; if PREFIX prepend with inla_, so that a=b yields inla_a=b. 
	 */
	char *p = NULL, *var = NULL;
	const int debug = 0;

	if (debug)
		printf("enter my_setenv with [%s]\n", str);

	p = strchr(str, '=');
	if (!p) {
		fprintf(stderr, "*** Error: Argument is void [%s]\n", str);
		exit(EXIT_FAILURE);
	}
	*p = '\0';
#if defined(WINDOWS)
	if (prefix) {
		GMRFLib_sprintf(&var, "inla_%s=%s", str, p + 1);
	} else {
		GMRFLib_sprintf(&var, "%s=%s", str, p + 1);
	}
	putenv(var);
	if (debug)
		printf("putenv \t%s\n", var);
#else
	if (prefix) {
		GMRFLib_sprintf(&var, "inla_%s", str);
	} else {
		GMRFLib_sprintf(&var, "%s", str);
	}
	setenv(var, p + 1, 1);
	if (debug)
		printf("\tsetenv %s=%s\n", var, p + 1);
#endif
	Free(var);
	return INLA_OK;
}

double my_gsl_sf_lnfact(int x)
{
	static int first = 1;
	static int nmax = 1048576L;
	static double *lng = NULL;

	if (first) {
#pragma omp critical (Name_764ffe066cfbd16ba7b1096b9b762ad9b1f8e669)
		if (first) {
			lng = Calloc(nmax, double);
			lng[0] = 0.0;
			for (int i = 1; i < nmax; i++) {
				lng[i] = lng[i - 1] + log((double) i);
			}
			first = 0;
		}
	}

	if (x >= nmax) {
		return gsl_sf_lnfact((unsigned int) x);
	} else {
		return lng[x];
	}

	assert(0 == 1);
	return NAN;
}

double my_gsl_sf_lngamma(double x)
{
	if (round(x) != x) {
		return gsl_sf_lngamma(x);
	} else {
		// x is an int, then use the cached values

		static int first = 1;
		static int nmax = 1048576;
		static double *lng = NULL;

		if (first) {
#pragma omp critical (Name_72a7f789baa1bbf55989513ddf777ec4ee6c91df)
			if (first) {
				lng = Calloc(nmax, double);
				lng[0] = NAN;
				lng[1] = 0.0;
				for (int i = 2; i < nmax; i++) {
					lng[i] = lng[i - 1] + log((double) (i - 1));
				}
				first = 0;
			}
		}
		if (x >= nmax) {
			return gsl_sf_lngamma(x);
		} else {
			return lng[(int) round(x)];
		}
	}

	assert(0 == 1);
	return NAN;
}

int my_gsl_sf_lnfact_e(const unsigned int n, gsl_sf_result *result)
{
	// copy of gsl_sf_lnfact_e

	result->val = my_gsl_sf_lnfact((int) n);
	result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
	return GSL_SUCCESS;
}

int my_gsl_sf_lnchoose_e(unsigned int n, unsigned int m, gsl_sf_result *result)
{
	// copy of gsl_sf_lnchoose_e

	assert(m <= n);
	if (m == n || m == 0) {
		result->val = 0.0;
		result->err = 0.0;
	} else {
		gsl_sf_result nf;
		gsl_sf_result mf;
		gsl_sf_result nmmf;
		if (m * 2 > n)
			m = n - m;
		my_gsl_sf_lnfact_e(n, &nf);
		my_gsl_sf_lnfact_e(m, &mf);
		my_gsl_sf_lnfact_e(n - m, &nmmf);
		result->val = nf.val - mf.val - nmmf.val;
		result->err = nf.err + mf.err + nmmf.err;
		result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
	}
	return GMRFLib_SUCCESS;
}

double my_gsl_sf_lnchoose(unsigned int n, unsigned int m)
{
	gsl_sf_result r = { 0, 0 };
	my_gsl_sf_lnchoose_e(n, m, &r);
	return (r.val);
}

double my_gsl_sf_lnbeta(double a, double b)
{
	double ab_min, ab_max;

	if (b < a) {
		ab_min = b;
		ab_max = a;
	} else {
		ab_min = a;
		ab_max = b;
	}

	ab_min = DMAX(DBL_EPSILON, ab_min);
	double n = ab_max / ab_min;
	const double threshold = 100.0 / DBL_EPSILON;

	if (n > threshold) {
		// log(Beta(a,a/n)) = as n->infinity + symmetry

/*			
 *			asympt(log(Beta(a,a/n)), n,2);
 *                                                         (Psi(a~) + gamma) a~      1
 *                                        ln(n) - ln(a~) - -------------------- + O(----)
 *                                                                  n                 2
 *                                                                                   n
 */

		return (log(n) - log(ab_max) - (gsl_sf_psi(ab_max) + 0.5772156649015328606065120) * ab_max / n);
	} else {
		return (gsl_sf_lnbeta(a, b));
	}
}

double my_betabinomial_helper4(int n, double a, double *work)
{
	const int roll = 4L;
	double s0 = 0.0;
	div_t d = div(n, roll);
	int m = d.quot * roll;
	int nn = d.quot;

#pragma omp simd
	for (int i = 0; i < nn; i++) {
		double aa = i * roll + a;
		work[i] = aa * (aa + 1) * (aa + 2) * (aa + 3);
	}

	GMRFLib_log(nn, work, work);
	s0 = GMRFLib_dsum(nn, work);

	if (d.rem) {
		double aa = m + a;
		double s = aa;
#pragma omp simd reduction(*: s)
		for (int i = 1; i < d.rem; i++) {
			s *= (aa + i);
		}
		s0 += log(s);
	}

	return (s0);
}

double my_betabinomial_helper8(int n, double a, double *work)
{
	const int roll = 8L;
	double s0 = 0.0;
	div_t d = div(n, roll);
	int m = d.quot * roll;
	int nn = d.quot;

#pragma omp simd
	for (int i = 0; i < nn; i++) {
		double aa = i * roll + a;
		work[i] = aa * (aa + 1) * (aa + 2) * (aa + 3) * (aa + 4) * (aa + 5) * (aa + 6) * (aa + 7);
	}

	GMRFLib_log(nn, work, work);
	s0 = GMRFLib_dsum(nn, work);

	if (d.rem) {
		double aa = m + a;
		double s = aa;
#pragma omp simd reduction(*: s)
		for (int i = 1; i < d.rem; i++) {
			s *= (aa + i);
		}
		s0 += log(s);
	}

	return (s0);
}

double my_betabinomial_helper16(int n, double a, double *work)
{
	const int roll = 16L;
	double s0 = 0.0;
	div_t d = div(n, roll);
	int m = d.quot * roll;
	int nn = d.quot;

#pragma omp simd
	for (int i = 0; i < nn; i++) {
		double aa = i * roll + a;
		work[i] = aa * (aa + 1) * (aa + 2) * (aa + 3) * (aa + 4) * (aa + 5) * (aa + 6) * (aa + 7)
		    * (aa + 8) * (aa + 9) * (aa + 10) * (aa + 11) * (aa + 12) * (aa + 13) * (aa + 14) * (aa + 15);
	}

	GMRFLib_log(nn, work, work);
	s0 = GMRFLib_dsum(nn, work);

	if (d.rem) {
		double aa = m + a;
		double s = aa;
#pragma omp simd reduction(*: s)
		for (int i = 1; i < d.rem; i++) {
			s *= (aa + i);
		}
		s0 += log(s);
	}

	return (s0);
}

double my_betabinomial_helper_core(int n, double a, double *work, int roll)
{
	div_t d = div(n, roll);
	int m = d.quot * roll;
	int nn = d.quot;

	for (int i = 0; i < nn; i++) {
		double aa = i * roll + a;
		double s = 1.0;
#pragma omp simd reduction(*: s)
		for (int j = 0; j < roll; j++) {
			s *= (aa + j);
		}
		work[i] = s;
	}

	GMRFLib_log(nn, work, work);
	double s0 = GMRFLib_dsum(nn, work);

	if (d.rem) {
		double aa = m + a;
		double s = aa;
#pragma omp simd reduction(*: s)
		for (int j = 1; j < d.rem; j++) {
			s *= (aa + j);
		}
		s0 += log(s);
	}

	return (s0);
}

double my_betabinomial_helper(int n, double a, double *work)
{
	return (n <= 104L ? my_betabinomial_helper8(n, a, work) : my_betabinomial_helper16(n, a, work));
}

double my_betabinomial(int y, int n, double a, double b, double *work)
{
	// WORK needs to be >= n
	double s1 = my_betabinomial_helper(y, a, work);
	double s2 = my_betabinomial_helper(n - y, b, work);
	double s3 = my_betabinomial_helper(n, a + b, work);

	return (s1 + s2 - s3);
}

double my_betabinomial2(int y, int n, double a, double b)
{
	double work[n];

	// using Gamma(1+z)=z*Gamma(z), we can get this
	double ladd = 0.0;
	while (a > 1.0) {
		a--;
		ladd += log((y + a) * (a + b) / (n + a + b) / a);
	}
	while (b > 1.0) {
		b--;
		ladd += log((n - y + b) * (a + b) / (n + a + b) / b);
	}

	// here we have 0<a<1, 0<b<1, but NOT a+b<1.
	// this could be helpful creating approximations
	double s1 = my_betabinomial_helper(y, a, work);
	double s2 = my_betabinomial_helper(n - y, b, work);
	double s3 = my_betabinomial_helper(n, a + b, work);
	return (s1 + s2 - s3 + ladd);
}

double my_lambert_W0(double y)
{
	double val = 0.0;
	my_lambert_W0s(1, &y, &val);

	return val;
}

void my_lambert_W0s(int m, double *y, double *res)
{
	// solve for x, so that x*exp(x)=y. This is lambert_W0(y). cache and interpolate the function in log-log scale

	static double logy_lim[2] = { -10.0, 20.0 };
	static GMRFLib_spline_tp *spline_lambert_W0 = NULL;

	if (!spline_lambert_W0) {
#pragma omp critical (Name_8d6da91249d12a4b187c29c3d315ea60b99212fb)
		{
			if (!spline_lambert_W0) {
				int n0 = 140;
				int n = n0 + 1;
				double dy = (logy_lim[1] - logy_lim[0]) / n0;
				double *work = Calloc(2 * n, double);
				double *xx = work;
				double *yy = work + n;

				for (int i = 0; i < n; i++) {
					yy[i] = logy_lim[0] + i * dy;
					xx[i] = log(gsl_sf_lambert_W0(exp(yy[i])));
				}
				spline_lambert_W0 = GMRFLib_spline_create(yy, xx, n);
				Free(work);
			}
		}
	}

	for (int k = 0; k < m; k++) {
		if (y[k] > 0.0) {
			double log_y = log(y[k]);
			if (log_y < logy_lim[1]) {
				// this version adds an extra Newton-R correction step. then we can do the caching less accurate
				double theta = GMRFLib_spline_eval(log_y, spline_lambert_W0);
				double exp_theta = exp(theta);
				double err = theta + exp_theta - log_y;
				double t1 = 1.0 + exp_theta;
				theta -= err / (t1 + err * exp_theta / t1);
				res[k] = exp(theta);
			} else {
				res[k] = gsl_sf_lambert_W0(y[k]);
			}
		} else {
			if (ISZERO(y[k]) || ISNAN(y[k])) {
				res[k] = y[k];
			} else {
				assert(0 == 1);
			}
		}
	}
}

double my_lbell(int y)
{
	// return log(Bell(y)/y!)

	static double *lbell = NULL;
	static int ymax = 256L;

	if (!lbell || y > ymax) {
#pragma omp critical (Name_47ba9e44a455e20c4bf966b4712c6ece203cb147)
		{
			if (!lbell || y > ymax) {
				int n = IMAX(2 * ymax, y);
				double *lbell_tmp = my_compute_lbell(n);

				ymax = n;
				lbell = lbell_tmp;
			}
		}
	}

	return ((y < 0) ? NAN : lbell[y]);
}

double *my_compute_lbell(int nmax)
{
	// return a vector with log(Bell(y)/y!) for i=0...nmax, computed using the recursive result: B_{n+1} = \sum_{k=0}^n Choose(n,k) B_k.
	// Since we cache these numbers, the time used is less of an issue.

	nmax = IMAX(7, nmax);

	double *log_bell = Calloc(nmax + 1, double);
	double *terms = Calloc(nmax, double);

	log_bell[0] = log_bell[1] = 0.0;
	for (int n = 2; n <= nmax; n++) {
		int n1 = n - 1;

#pragma omp parallel for if(n > 1024L)
		for (int k = 0; k < n; k++) {
			terms[k] = my_gsl_sf_lnchoose((unsigned int) n1, (unsigned int) k) + log_bell[k];
		}
		QSORT_FUN((void *) terms, (size_t) n, sizeof(double), GMRFLib_dcmp);

		// need to compute log(exp(terms[0]) + ... + exp(terms[n1])), do this the obvious way: summing the smallest terms first using
		// the largest element (the last one) as scaling.
		double sum = 0.0;
		for (int k = 0; k < n1; k++) {
			sum += exp(terms[k] - terms[n1]);
		}
		log_bell[n] = terms[n1] + log1p(sum);
	}

	for (int n = 0; n <= nmax; n++) {
		// adjust with n! as this is they way it is used in the Bell-distribution
		log_bell[n] -= my_gsl_sf_lnfact(n);
	}

	Free(terms);
	return (log_bell);
}
