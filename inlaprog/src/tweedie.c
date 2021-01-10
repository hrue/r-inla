
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

#include "rmath.h"
#undef ISNAN

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"
#include "tweedie.h"

// the dtweedie-code is modified from tweedie.c in the (GPL'ed) cplm_0.7-9.tar.gz package of Wayne Zhang
// <actuary_zhang@hotmail.com>, to tailor it to the INLA use.

#define TWEEDIE_DROP 40.0
#define TWEEDIE_INCRE 4
#define TWEEDIE_NTERMS_ADD 8
#define TWEEDIE_MAX_IDX 16384

/**
 * n the number of observations
 * y the observation
 * mu the vector of means
 * phi scalar: the dispersion parameter
 * p scalar: the index parameter
 * ldens the vector that stores the computed log density
 */

void dtweedie(int n, double y, double *mu, double phi, double p, double *ldens)
{
	static struct {
		int nterms;
		double save_p;
		double *work;
		double *wwork;
	} cache = { -1, -9999.9999, NULL, NULL };
#pragma omp threadprivate(cache)

	static size_t cache_count[] = { 0, 0, 0 };
	static double sum_nterms = 0.0;

	double p1 = p - 1.0, p2 = 2.0 - p;
	double a = -p2 / p1, a1 = 1.0 / p1;
	double cc, j, w, sum_ww = 0.0, ww_max = 0.0;

	int use_interpolation = 1, nterms, k, i, one = 1, k_low = -1, reuse = 0;
	int verbose = 1, show_stat = 0;

	if (cache.nterms < 0) {
		if (verbose) {
			printf("\tdtweedie: initialize cache. use_interpolation= %1d\n", use_interpolation);
		}
		cache.nterms = 0;
		cache.work = Calloc(2 * TWEEDIE_MAX_IDX, double);
		cache.wwork = cache.work + TWEEDIE_MAX_IDX;
	}

	if (ISZERO(y)) {
		for (int i = 0; i < n; i++) {
			ldens[i] = -pow(mu[i], p2) / (phi * p2);
		}
		return;
	}


	int jmax;
	double logz, logz_stripped;

	cc = a * log(p1) - log(p2);
	jmax = DMAX(1.0, pow(y, p2) / (phi * p2));
	logz = -a * log(y) - a1 * log(phi) + cc;
	logz_stripped = cc;

	cc = logz + a1 + a * log(-a);
	j = jmax;
	w = a1 * j;
	while (1) {
		j += TWEEDIE_INCRE;
		if (j * (cc - a1 * log(j)) < (w - TWEEDIE_DROP))
			break;
	}
	nterms = IMIN(TWEEDIE_MAX_IDX, ceil(j));
	sum_nterms += nterms;

	if (!(p == cache.save_p)) {
		reuse = 0;
		k_low = 0;
		cache_count[0]++;
		show_stat = 1;
	} else {
		if (nterms > cache.nterms) {
			// if only the 'nterms' have increased, we can add the terms.
			nterms = IMIN(TWEEDIE_MAX_IDX, nterms + TWEEDIE_NTERMS_ADD);
			reuse = 0;
			k_low = cache.nterms;
			cache_count[1]++;
			if (0 && verbose) {
				printf("\tdtweedie: increase cache from %d to %d\n", cache.nterms, nterms);
			}
		} else {
			// this is the case where we are fine and can reuse what we have
			reuse = 1;
			cache_count[2]++;
		}
	}

	if (!reuse) {
		if (!use_interpolation) {
			for (k = k_low; k < nterms; k++) {
				j = k + one;
				cache.wwork[k] = j * logz_stripped - inla_lgamma_fast(1.0 + j) - inla_lgamma_fast(-a * j);
			}
		} else {
			for (k = k_low; k < nterms + 1; k += 2) {
				j = k + one;
				cache.wwork[k] = j * logz_stripped - inla_lgamma_fast(1.0 + j) - inla_lgamma_fast(-a * j);

				if (k > 0) {
					double estimate = 0.5 * (cache.wwork[k] + cache.wwork[k - 2]);
					double correction = (1.0 - a) * (1.0 / j + 1.0 / SQR(j));
					double limit = 0.05;
					if (ABS(correction) > limit) {
						double jj = j - 1.0;
						cache.wwork[k - 1] = jj * logz_stripped - inla_lgamma_fast(1.0 + jj) - inla_lgamma_fast(-a * jj);
					} else {
						cache.wwork[k - 1] = estimate + correction;
					}
				}
			}
		}
		cache.save_p = p;
		cache.nterms = nterms;
	}
	double term_removed = -a * log(y) -a1 * log(phi);      // the terms we have removed from 'logz'
	double lim = -20.72326584;			       // log(1.0e-9)
	int idx_max = -1;

	sum_ww = 0.0;
	ww_max = cache.wwork[0] + one * term_removed;
	for (k = 0; k < nterms; k++) {
		j = k + one;
		cache.work[k] = cache.wwork[k] + term_removed * j;
		if (cache.work[k] > ww_max) {
			ww_max = cache.work[k];
			idx_max = k;
		}
	}

	for (k = idx_max; k < nterms; k++) {		       /* assume there is one mode */
		double tmp = cache.work[k] - ww_max;
		if (tmp > lim) {
			sum_ww += exp(tmp);
		} else {
			break;
		}
	}
	for (k = idx_max - 1; k >= 0; k--) {
		double tmp = cache.work[k] - ww_max;
		if (tmp > lim) {
			sum_ww += exp(tmp);
		} else {
			break;
		}
	}

	for (i = 0; i < n; i++) {
		ldens[i] = -pow(mu[i], p2) / (phi * p2) - y / (phi * p1 * pow(mu[i], p1)) - log(y) + log(sum_ww) + ww_max;
	}

	if (verbose) {
		static size_t count = 0;
		size_t ntot = cache_count[0] + cache_count[1] + cache_count[2];

		count++;
		if (show_stat) {
			printf("\tdtweedie: ntimes=%zu rebuild=%.3f%%  adjust=%.3f%% reuse=%.3f%% mean.nterms=%1d\n",
			       ntot,
			       100.0 * (double) cache_count[0] / (double) ntot,
			       100.0 * (double) cache_count[1] / (double) ntot,
			       100.0 * (double) cache_count[2] / (double) ntot, (int) (sum_nterms / (double) count));
		}
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

	// result from a test: Rmath takes 30% more time than the GSL call to 'pgamma'
	// retval += prob * MATHLIB_FUN(pgamma)(y, n * alpha, gamma, 1, 0);

	// when sd/mean < low, we basically are in the Gaussian regime. this is for the case Gamma(n*alpha, gamma), so mean=
	// n*alpha*gamma, and sd= sqrt(n*alpha)*gamma, so gamma cancel.
	double low = 0.25;
	int n, n_gauss = (int) (1.0 / (alpha * SQR(low)) + 0.5);

	// stride should scale with stdev sqrt(lambda). the below choice gives (about)
	// stride=1 for lam < 5
	// stride=2 for 5 <= lam < 13
	// stride=3 for 13 <= lam < 25
	// stride=4 for 25 <= lam < 41, etc...
	int stride = IMAX(1, (int) (sqrt(lambda / 2.0) + 0.5));

	retval = pacc = prob = gsl_ran_poisson_pdf((unsigned int) 0, lambda);

	if (stride == 1) {
		for (n = 1; pacc < plim; n++, pacc += prob) {
			// prob = gsl_ran_poisson_pdf((unsigned int) n, lambda);
			prob *= lambda / n;

			if (n < n_gauss) {
				retval += prob * gsl_cdf_gamma_P(y, n * alpha, gamma);
			} else {
				retval += prob * inla_Phi_fast((y - (n * alpha * gamma)) / (sqrt(n * alpha) * gamma));
			}
		}
	} else {
		// use log-linear interpolation for the cdf's between each stride
		double cdf, lcdf, lcdf_prev;

		// we have to start manually
		n = 1;
		prob *= lambda / n;
		pacc += prob;
		cdf = gsl_cdf_gamma_P(y, n * alpha, gamma);
		lcdf_prev = log(cdf);
		retval += prob * cdf;

		for (n = stride + 1; pacc < plim; n += stride) {
			if (n < n_gauss) {
				cdf = gsl_cdf_gamma_P(y, n * alpha, gamma);
			} else {
				cdf = inla_Phi_fast((y - (n * alpha * gamma)) / (sqrt(n * alpha) * gamma));
			}
			lcdf = lcdf_prev = log(cdf);

			for (int k = stride - 1; k >= 1; k--) {
				prob *= lambda / (n - k);
				pacc += prob;
				retval += prob * exp((k * lcdf_prev + (stride - k) * lcdf) / (double) stride);
			}

			prob *= lambda / n;
			pacc += prob;
			retval += prob * cdf;
		}
	}

	return (retval);
}
