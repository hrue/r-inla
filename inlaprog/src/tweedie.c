
/* tweedie.c
 * 
 * Copyright (C) 2021-2022 Havard Rue
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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-const-variable"
static const char GitID[] = "file: " __FILE__ "  " GITCOMMIT;
#pragma GCC diagnostic pop

// the 'dtweedie'-code is inspired from tweedie.c in the (GPL'ed) cplm_0.7-9.tar.gz package of Wayne Zhang
// <actuary_zhang@hotmail.com>, but got largly rewritten to incorporate caching and interpolation and other optimization
// features. This to optimize it to how its used within the INLA context.

#define TWEEDIE_DROP 40.0
#define TWEEDIE_INCRE 1.2
#define TWEEDIE_MAX_IDX 16384

// this is the G.Nemes (2007) approximation from https://en.wikipedia.org/wiki/Stirling's_approximation
#define LGAMMA_FAST(_x) ((_x) < 1.0 ? gsl_sf_lngamma(_x) : 0.5 * (1.837877066409345 - log(_x)) + (_x) * (log((_x) + 1.0/( 12.0*(_x) - 0.1/(_x))) - 1.0))

/**
 * n scalar length of mu
 * y scalar the observation
 * mu vector means
 * phi scalar dispersion parameter
 * p scalar index parameter
 * ldens vector the computed log densities
 */

typedef struct {
	int nterms;
	int interpolation_ok;
	double save_p;
	double *work;
	double *wwork;
} dtweedie_cache_tp;

static dtweedie_cache_tp **cache = NULL;
static double *lgammas = NULL;

// must be called before use to initialize cache
void dtweedie_init_cache(void)
{
	if (!cache) {
		cache = Calloc(GMRFLib_CACHE_LEN, dtweedie_cache_tp *);
		for (int i = 0; i < GMRFLib_CACHE_LEN; i++) {
			cache[i] = Calloc(GMRFLib_CACHE_LEN, dtweedie_cache_tp);
			cache[i]->nterms = -1;
			cache[i]->interpolation_ok = 0;
			cache[i]->save_p = -9999.9999;
			cache[i]->work = NULL;
			cache[i]->wwork = NULL;
		}
		lgammas = Calloc(TWEEDIE_MAX_IDX, double);
		for (int i = 0; i < TWEEDIE_MAX_IDX; i++) {
			lgammas[i] = my_gsl_sf_lngamma(1.0 + i);	/* factorials for integers */
		}
	}
}

void dtweedie(int n, double y, double *mu, double phi, double p, double *ldens)
{
	// this function cache wrt 'p' only

	static size_t cache_count[] = { 0, 0, 0 };
	static double sum_nterms = 0.0;

	double p1 = p - 1.0, p2 = 2.0 - p;
	double a = -p2 / p1, a1 = 1.0 / p1;
	double cc, w, sum_ww = 0.0, ww_max = 0.0, lsum_ww, ly;
	double jmax, logz, logz_stripped;
	int id = -1, use_interpolation = 1, nterms, k, i, j, one = 1, k_low = -1, reuse = 0, verbose = 0, show_stat = 0;

	GMRFLib_CACHE_SET_ID(id);
	dtweedie_cache_tp *cache_ptr = cache[id];

	if (cache_ptr->nterms < 0) {
		if (verbose) {
			printf("\tdtweedie: initialize cache[%1d] use_interpolation= %1d\n", id, use_interpolation);
		}
		cache_ptr->nterms = 0;
		cache_ptr->work = Calloc(2 * TWEEDIE_MAX_IDX, double);
		cache_ptr->wwork = cache_ptr->work + TWEEDIE_MAX_IDX;
	}

	if (ISZERO(y)) {
		for (i = 0; i < n; i++) {
			ldens[i] = -pow(mu[i], p2) / (phi * p2);
		}
		return;
	}

	ly = log(y);
	cc = a * log(p1) - log(p2);
	jmax = DMAX(1.0, pow(y, p2) / (phi * p2));
	jmax = DMAX(jmax, 10.0);
	logz = -a * ly - a1 * log(phi) + cc;
	logz_stripped = cc;

	cc = logz + a1 + a * log(-a);
	w = a1 * jmax;
	double ljmax = log(jmax);
	double ljmax_add = log((double) TWEEDIE_INCRE);
	while (1) {
		jmax *= TWEEDIE_INCRE;
		ljmax += ljmax_add;
		if (jmax * (cc - a1 * ljmax) < (w - TWEEDIE_DROP))
			break;
	}

	nterms = IMIN(TWEEDIE_MAX_IDX, (int) ceil(jmax));
	if (nterms == TWEEDIE_MAX_IDX) {
		fprintf(stderr, "\ndtweedie: Reached upper limit. Increase TWEEDIE_MAX_IDX or scale your data!\n");
	}
	if (verbose) {
#pragma omp atomic
		sum_nterms += nterms;
	}

	if (!(p == cache_ptr->save_p)) {
		reuse = 0;
		cache_ptr->interpolation_ok = 0;	       /* yes, we need to start again */
		k_low = 0;
		if (verbose) {
#pragma omp atomic
			cache_count[0]++;
			show_stat = 1;
		}

	} else {
		if (nterms > cache_ptr->nterms) {
			// if only the 'nterms' have increased, we can add the terms.
			nterms = IMIN(TWEEDIE_MAX_IDX, nterms);
			reuse = 0;
			k_low = cache_ptr->nterms;
			if (verbose) {
#pragma omp atomic
				cache_count[1]++;
				if (verbose) {
					printf("\tdtweedie: increase cache from %d to %d\n", cache_ptr->nterms, nterms);
				}
			}
		} else {
			// this is the case where we are fine and can reuse what we have
			reuse = 1;
			if (verbose) {
#pragma omp atomic
				cache_count[2]++;
			}
		}
	}

	if (!reuse) {
		if (!use_interpolation) {
			for (k = k_low; k < nterms; k++) {
				double xx = -a * j;
				j = k + one;
				cache_ptr->wwork[k] = j * logz_stripped - lgammas[j] - LGAMMA_FAST(xx);
			}
		} else {
			// correction term has expansion c[0]/j + c[1]/j^2 + c[2]/j^3 + c[3]/j^4
			double coofs[] = {
				0.5 - 0.5 * a,
				0.5 - 0.5 * a,
				-(1.0 + (-8.0 + 7.0 * a) * a) / (12.0 * a),
				-(1.0 + (-4.0 + 3.0 * a) * a) / (4.0 * a)
			};

			for (k = k_low; k < nterms + 1; k += 2) {
				j = k + one;
				cache_ptr->wwork[k] = j * logz_stripped - lgammas[j] - LGAMMA_FAST(-a * j);

				if (k > k_low) {
					// no need to check before around here
					if (k >= 18) {
						// for ever increasing 'j', we check if the error in the interpolation is small, and if
						// so, we use interpolation for all j' > j (for this cache). any time the cache is
						// rebuilt, then we start over again
						double estimate = 0.5 * (cache_ptr->wwork[k] + cache_ptr->wwork[k - 2]);
						double inv_j = 1.0 / (double) j;
						double correction = inv_j * (coofs[0] + inv_j * (coofs[1] + inv_j * (coofs[2] + inv_j * coofs[3])));
						double limit = 1.0e-5;

						if (cache_ptr->interpolation_ok) {
							cache_ptr->wwork[k - 1] = estimate + correction;
						} else {
							int jj = j - 1.0;
							cache_ptr->wwork[k - 1] = jj * logz_stripped - lgammas[jj] - LGAMMA_FAST(-a * jj);
							if (ABS(cache_ptr->wwork[k - 1] - (estimate + correction)) < limit) {
								// from here on, we can safely use interpolation
								cache_ptr->interpolation_ok = 1;
								if (verbose) {
									printf("\tdtweedie: set interpolation_ok=1 at k= %d\n", k);
								}
							}
						}
					} else {
						int jj = j - 1.0;
						cache_ptr->wwork[k - 1] = jj * logz_stripped - lgammas[jj] - LGAMMA_FAST(-a * jj);
					}
				}
			}
		}
		cache_ptr->save_p = p;
		cache_ptr->nterms = nterms;
	}
	double term_removed = -a * ly - a1 * log(phi);	       // the terms we have removed from 'logz'
	double lim = -20.72326584;			       // log(1.0e-9)

	int idx_max = 0;
	sum_ww = 0.0;
	ww_max = cache_ptr->wwork[0] + one * term_removed;
	for (k = 0; k < nterms; k++) {
		j = k + one;
		cache_ptr->work[k] = cache_ptr->wwork[k] + term_removed * j;
		if (cache_ptr->work[k] > ww_max) {
			ww_max = cache_ptr->work[k];
			idx_max = k;
		}
	}

	for (k = idx_max; k < nterms; k++) {		       /* assume there is one mode */
		double tmp = cache_ptr->work[k] - ww_max;
		sum_ww += exp(tmp);
		if (tmp < lim) {
			break;
		}
	}
	for (k = idx_max - 1; k >= 0; k--) {
		double tmp = cache_ptr->work[k] - ww_max;
		sum_ww += exp(tmp);
		if (tmp < lim) {
			break;
		}
	}
	lsum_ww = log(sum_ww);

	for (i = 0; i < n; i++) {
		ldens[i] = -pow(mu[i], p2) / (phi * p2) - y / (phi * p1 * pow(mu[i], p1)) - ly + lsum_ww + ww_max;
	}

	if (verbose) {
		static size_t count = 0;
		size_t ntot = cache_count[0] + cache_count[1] + cache_count[2];

#pragma omp atomic
		count++;

		if (show_stat) {
#pragma omp critical (Name_a5c97756f6b00566ba28b5dce73bc4447f991240)
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
