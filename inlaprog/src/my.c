
/* my.c
 * 
 * Copyright (C) 2016-2022 Havard Rue
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

#include <string.h>
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "inla.h"
#include "my.h"
#include "my-fix.h"
#include "GMRFLib/GMRFLibP.h"
#include "GMRFLib/GMRFLib.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-const-variable"
static const char GitID[] = "file: " __FILE__ "  " GITCOMMIT;
#pragma GCC diagnostic pop

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
	static int nmax = 1048576;
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

int my_gsl_sf_lnfact_e(const unsigned int n, gsl_sf_result * result)
{
	// copy of gsl_sf_lnfact_e

	result->val = my_gsl_sf_lnfact((int) n);
	result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
	return GSL_SUCCESS;
}

int my_gsl_sf_lnchoose_e(unsigned int n, unsigned int m, gsl_sf_result * result)
{
	// copy of gsl_sf_lnchoose_e

	if (m > n) {
		assert(0 == 1);
	} else if (m == n || m == 0) {
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

double my_betabinomial_helper(int n, double a)
{
	const int roll = 4L;				       /* dont change this */
	double s0 = 0.0;
	div_t d = div(n, roll);
	int m = d.quot * roll;

#pragma GCC ivdep
	for (int i = 0; i < m; i += roll) {
		double aa = i + a;
		s0 += log(aa * (aa + 1.0) * (aa + 2.0) * (aa + 3.0));
	}

	if (d.rem) {
		double aa = m + a;
		switch (d.rem) {
		case 1:
			s0 += log(aa);
			break;
		case 2:
			s0 += log(aa * (aa + 1.0));
			break;
		case 3:
			s0 += log(aa * (aa + 1.0) * (aa + 2.0));
			break;
		}
	}

	return (s0);
}
double my_betabinomial(int y, int n, double a, double b)
{
	double s1 = my_betabinomial_helper(y, a);
	double s2 = my_betabinomial_helper(n - y, b);
	double s3 = my_betabinomial_helper(n, a + b);
	return (s1 + s2 - s3);
}

double my_betabinomial2(int y, int n, double a, double b)
{
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
	double s1 = my_betabinomial_helper(y, a);
	double s2 = my_betabinomial_helper(n - y, b);
	double s3 = my_betabinomial_helper(n, a + b);
	return (s1 + s2 - s3 + ladd);
}
