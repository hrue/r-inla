
/* tweedie.c
 * 
 * Copyright (C) 2021-2021 Havard Rue
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
static const char GitID[] = GITCOMMIT;

#include <math.h>
#include <strings.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>
#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

#include "tweedie.h"

// This code is modified from tweedie.c in the (GPL'ed) cplm_0.7-9.tar.gz
// package of Wayne Zhang <actuary_zhang@hotmail.com>

#define TWEEDIE_DROP 37.0
#define TWEEDIE_INCRE 5
#define TWEEDIE_NTERM 20000

/**
 * Compute the log density for the tweedie compound Poisson distribution.
 * This is based on the dtweedie.series function in R, but bounds
 * are not determined using all observations because this could result in 
 * dramatically slower computation in certain circumstances. 
 *
 * @param n the number of observations
 * @param y the observation
 * @param mu the vector of means
 * @param phi scalar: the dispersion parameter
 * @param p scalar: the index parameter
 * @param ldens the vector that stores the computed log density
 */

void dtweedie(int n, double y, double *mu, double phi, double p, double *ldens)
{
	static struct {
		int n_max;
		int nn_max;
		double *work;
		double *wwork;
	} store = { -1, -1, NULL, NULL };
#pragma omp threadprivate(store)

	if (store.n_max < 0) {
		store.n_max = 8;
		store.work = Calloc(store.n_max * 4, double);
		store.nn_max = 128;
		store.wwork = Calloc(store.nn_max, double);
	}

	double p1 = p - 1.0, p2 = 2.0 - p;
	double a = -p2 / p1, a1 = 1.0 / p1;
	double cc, j, w, sum_ww = 0.0, ww_max = 0.0;

	if (y == 0.0) {
		for (int i = 0; i < n; i++) {
			ldens[i] = -pow(mu[i], p2) / (phi * p2);
		}
		return;
	}
	// only need the lower bound and the # terms to be stored 
	int jh, *jl, *jd, jd_max = 0, jl_low = TWEEDIE_NTERM;
	double *jmax, *logz;

	if (n > store.n_max) {
		store.n_max = n * 2;
		Free(store.work);
		store.work = Calloc(store.n_max * 4, double);
	}

	jl = (int *) (store.work + 0);
	jd = (int *) (store.work + n);
	jmax = store.work + 2 * n;
	logz = store.work + 3 * n;

	cc = a * log(p1) - log(p2);
	for (int i = 0; i < n; i++) {
		jmax[i] = DMAX(1.0, pow(y, p2) / (phi * p2));
		logz[i] = -a * log(y) - a1 * log(phi) + cc;
		P(logz[i]);
	}

	// find bounds in the summation 
	for (int i = 0; i < n; i++) {
		// locate upper bound 
		cc = logz[i] + a1 + a * log(-a);
		j = jmax[i];
		w = a1 * j;
		while (1) {
			j += TWEEDIE_INCRE;
			if (j * (cc - a1 * log(j)) < (w - TWEEDIE_DROP))
				break;
		}
		jh = ceil(j);
		// locate lower bound 
		j = jmax[i];
		while (1) {
			j -= TWEEDIE_INCRE;
			if (j < 1 || j * (cc - a1 * log(j)) < w - TWEEDIE_DROP)
				break;
		}
		jl[i] = IMAX(1, floor(j));
		jd[i] = jh - jl[i] + 1;
		jl_low = IMIN(jl_low, jd[i]);
		jd_max = IMAX(jd_max, jd[i]);
	}

	// set limit for # terms in the sum 
	int nterms = IMIN(jd_max, TWEEDIE_NTERM);

	if (nterms > store.nn_max) {
		store.nn_max = nterms * 2;
		Free(store.wwork);
		store.wwork = Calloc(store.nn_max, double);
	}

	// this term does not depend on mu[i], so we only need to do this once, as only 'mu[i]' varies with 'i'.
	for (int k = 0; k < nterms; k++) {
		j = k + jl_low;
		store.wwork[k] = j * logz[0] - my_gsl_sf_lngamma(1 + j) - my_gsl_sf_lngamma(-a * j);
		ww_max = (k == 0 ? store.wwork[k] : DMAX(ww_max, store.wwork[k]));
	}
	sum_ww = 0.0;
	for (int k = 0; k < nterms; k++) {
		sum_ww += exp(store.wwork[k] - ww_max);
	}

	// evaluate series using the finite sum
	for (int i = 0; i < n; i++) {
		ldens[i] = -pow(mu[i], p2) / (phi * p2);
		ldens[i] += -y / (phi * p1 * pow(mu[i], p1)) - log(y) + log(sum_ww) + ww_max;
	}

	return;
}

double ptweedie(double y, double mu, double phi, double p)
{
	// compute Prob(Y <= y)

	double lambda = pow(mu, 2.0 - p) / (phi * (2.0 - p));
	double alpha = (2.0 - p) / (p - 1.0);
	double gamma = phi * (p - 1.0) * pow(mu, p - 1.0);
	double plim = 0.999;
	double retval, prob, pacc;

	retval = pacc = gsl_ran_poisson_pdf((unsigned int) 0, lambda);
	for (int n = 1; pacc < plim; n++, pacc += prob) {
		prob = gsl_ran_poisson_pdf((unsigned int) n, lambda);
		retval += prob * gsl_cdf_gamma_P(y, n * alpha, gamma);
	}

	return (retval);
}
