#include <assert.h>
#include <math.h>
#include <strings.h>
#include <stdio.h>

#include "rmath.h"
#if defined(ISNAN)
#       undef ISNAN
#endif

#include "GMRFLib/GMRFLib.h"
#include "fast-math/special-functions.h"
#include "tweedie2.h"

// this code is very insprired from tweedie.c in cplm_0.7-9.tar.gz package of Wayne Zhang <actuary_zhang@hotmail.com>, but rewritten
// to incorporate caching to optimize it to how its used within the INLA context

#define TWEEDIE_DROP 30.0
#define TWEEDIE_INC 5

/**
 * n scalar length of mu
 * y scalar the observation
 * mu vector means
 * phi scalar dispersion parameter
 * p scalar index parameter
 * ldens vector the computed log densities
 */

typedef struct {
	int len;
	int lower;
	int upper;
	double save_p;
	double *w;
	double *arg;
	double *res;
	double *lgam_terms;
	double *lnfact;
} dtweedie_cache_tp;

static dtweedie_cache_tp **cache = NULL;
static int verbose = 0;

// must be called before use to initialize cache
// env INLA_DEBUG=dtweedie2_init_cache inla.mkl.work -v -t1 Model.ini
void dtweedie2_init_cache_idx(int idx)
{
#define LEN 64
	if (cache && cache[idx]) {
		return;
	}

	cache[idx] = Calloc(1, dtweedie_cache_tp);
	cache[idx]->save_p = -9999.9999;
	cache[idx]->len = LEN;
	cache[idx]->w = Calloc(LEN, double);
	cache[idx]->arg = Calloc(LEN, double);
	cache[idx]->res = Calloc(LEN, double);
	cache[idx]->lgam_terms = Calloc(LEN, double);
	cache[idx]->lnfact = Calloc(LEN, double);
	for (int j = 1; j < LEN; j++) {
		cache[idx]->lnfact[j] = cache[idx]->lnfact[j - 1] + log(j);
	}
#undef LEN
}
void dtweedie2_init_cache(void)
{
	if (!cache) {
#pragma omp critical (Name_92509c30f7c8ce2ff56520888da767c88a1ae7d4)
		if (!cache) {
			dtweedie_cache_tp **ccache = Calloc(GMRFLib_CACHE_LEN(), dtweedie_cache_tp *);
			verbose = 0; //GMRFLib_DEBUG_IF_TRUE();
			cache = ccache;
		}
	}
}

static void dtweedie2_adjust_cache(int idx, int nlen)
{
#define MINLEN 64
	if (nlen > cache[idx]->len) {
		int olen = cache[idx]->len;
		cache[idx]->len = IMIN(MINLEN, cache[idx]->len);
		while (cache[idx]->len < nlen) {
			cache[idx]->len *= 2;
		}
		nlen = cache[idx]->len;

		// we cache do new alloc here, no need for realloc
		Free(cache[idx]->w);
		Free(cache[idx]->arg);
		Free(cache[idx]->res);
		cache[idx]->w = Calloc(nlen, double);
		cache[idx]->arg = Calloc(nlen, double);
		cache[idx]->res = Calloc(nlen, double);

		// these needs realloc
		cache[idx]->lgam_terms = Realloc(cache[idx]->lgam_terms, nlen, double);
		cache[idx]->lnfact = Realloc(cache[idx]->lnfact, nlen, double);
		for (int j = olen; j < nlen; j++) {
			cache[idx]->lnfact[j] = cache[idx]->lnfact[j - 1] + log((double) j);
		}

		if (verbose) {
			// we ignore small numbers here...
			static double total_cache_size = 0.0;
			double change = (nlen - olen) * 5 * sizeof(double);
#pragma omp atomic
			total_cache_size += change;
			printf("\ttweedie2: extend cache[%1d] from len=%1d to %1d [total.size=%.2fMb]\n", idx, olen, nlen,
			       total_cache_size / SQR(1024.));
		}
	}
#undef MINLEN
}

void dtweedie2(int n, double y, double *mu, double phi, double p, double *ldens)
{
	double p1 = p - 1.0, p2 = 2.0 - p;
	double alpha = -p2 / p1, a1 = 1.0 / p1;

	if (ISZERO(y)) {
		for (int i = 0; i < n; i++) {
			ldens[i] = -pow(mu[i], p2) / (phi * p2);
		}
		return;
	}
#undef MEASURE_TIME
#if defined(MEASURE_TIME)
	static double tref[10] = { 0 };
	static double trefc = 0.0;
	tref[0] -= GMRFLib_timer();
#endif
	double ly = log(y);

	int id;
	GMRFLib_CACHE_SET_IDX(id);
	dtweedie2_init_cache_idx(id);
	dtweedie_cache_tp *c = cache[id];

	double cc = alpha * log(p1) - log(p2);
	double jmax = DMAX(10.0, pow(y, p2) / (phi * p2));
	double logz = -alpha * ly - a1 * log(phi) + cc;

	// lower and upper bound in the sum
	cc = logz + a1 + alpha * log(-alpha);
	double w = a1 * jmax;
	int jj = jmax;
	int inc = TWEEDIE_INC;
	while (1) {
		jj += inc;
		if (jj * (cc - a1 * log(jj)) < (w - TWEEDIE_DROP))
			break;
		inc += 2;				       /* speed it up a little */
	}
	int upper = jj;
	dtweedie2_adjust_cache(id, upper + 1);

	jj = jmax;
	inc = TWEEDIE_INC;
	while (1) {
		jj -= inc;
		if (jj < 1 || jj * (cc - a1 * log(jj)) < w - TWEEDIE_DROP)
			break;
		inc += 2;				       /* speed it up a little */
	}
	int lower = IMAX(1, floor(jj));

#if defined(MEASURE_TIME)
	tref[0] += GMRFLib_timer();
	tref[1] -= GMRFLib_timer();
#endif

	if (c->save_p != p) {
		for (int j = lower; j <= upper; j++) {
			c->arg[j] = -alpha * j;
		}
		LGAMMAfn_m(upper - lower + 1, c->arg + lower, c->res + lower);
#pragma omp simd
		for (int j = lower; j <= upper; j++) {
			c->lgam_terms[j] = c->lnfact[j] + c->res[j];
		}
		c->save_p = p;
		c->lower = lower;
		c->upper = upper;
	} else {
		// save_p == p
		if (lower < c->lower) {
			// only add lower term
#pragma omp simd
			for (int j = lower; j < c->lower; j++) {
				c->arg[j] = -alpha * j;
			}
			LGAMMAfn_m(c->lower - lower, c->arg + lower, c->res + lower);
#pragma omp simd
			for (int j = lower; j < c->lower; j++) {
				c->lgam_terms[j] = c->lnfact[j] + c->res[j];
			}
			c->lower = lower;
		}
		if (upper > c->upper) {
			// include ->upper as its cleaner code
#pragma omp simd
			for (int j = c->upper; j <= upper; j++) {
				c->arg[j] = -alpha * j;
			}
			LGAMMAfn_m(upper - c->upper + 1, c->arg + c->upper, c->res + c->upper);
#pragma omp simd
			for (int j = c->upper; j <= upper; j++) {
				c->lgam_terms[j] = c->lnfact[j] + c->res[j];
			}
			c->upper = upper;
		}
	}

#if defined(MEASURE_TIME)
	tref[1] += GMRFLib_timer();
	tref[2] -= GMRFLib_timer();
#endif

	c->w[lower] = lower * logz - c->lgam_terms[lower];
	double w_max = c->w[lower];
#pragma omp simd
	for (int j = lower + 1; j <= upper; j++) {
		// w[j] = j * logz - lgamma(1 + j) - lgamma(-alpha * j);
		c->w[j] = j * logz - c->lgam_terms[j];
		w_max = DMAX(w_max, c->w[j]);
	}

#pragma omp simd
	for (int j = lower; j <= upper; j++) {
		c->w[j] -= w_max;
	}
	GMRFLib_exp(upper - lower + 1, c->w + lower, c->res + lower);
	double sum_w = GMRFLib_dsum(upper - lower + 1, c->res + lower);

	for (int i = 0; i < n; i++) {
		ldens[i] = -pow(mu[i], p2) / (phi * p2);       // y == 0
		ldens[i] += -y / (phi * p1 * pow(mu[i], p1)) - ly + log(sum_w) + w_max;
	}

#if defined(MEASURE_TIME)
	tref[2] += GMRFLib_timer();
	trefc++;
	if ((int) trefc % 1000 == 0) {
		double s = 1.0 / GMRFLib_dsum(3, tref);
		for (int i = 0; i < 3; i++) {
			printf("chunk %1d: %.3f ", i, tref[i] * s);
		}
		printf("\n");
	}
#endif

#if 0
	// verify against the old version?
	static int first = 1;
	void dtweedie_init_cache(void);
	if (first)
		dtweedie_init_cache();
	first = 0;
	void dtweedie(int n, double y, double *mu, double phi, double p, double *ldens);
	for (int i = 0; i < n; i++) {
		double ld = 0;
		dtweedie(1, y, &(mu[i]), phi, p, &ld);
		if (ABS(ldens[i] - ld) > 0.001) {
			printf("i %d %.8f %.8f %.12f\n", i, ldens[i], ld, ldens[i] - ld);
			abort();
		}
	}
#endif

#if defined(MEASURE_TIME)
#       undef MEASURE_TIME
#endif
	return;
}

// this function is more or less a copy from the old version in tweedie.c
double ptweedie2(double y, double mu, double phi, double p)
{
	// compute Prob(Y <= y)

#define LOGNFACTORIAL(y_) (y_ < c->len ? c->lnfact[y_] : my_gsl_sf_lnfact(y_))
#define LOG_PDF_POISSON(y_) ((y_)*log_lambda - lambda - LOGNFACTORIAL((int) (y_)))
#define CDF(n_) ((n_) < n_gauss ?					\
		 gsl_cdf_gamma_P(y, (n_) * alpha, gamma) :		\
		 inla_cdf_normal_fast((y - ((n_) * c1)) / (sqrt((n_)) * c2)))
#define LOG_CDF(n_) ((n_) < n_gauss ?					\
		     MATHLIB_FUN(pgamma) (y, (n_) * alpha, gamma, 1, 1) : \
		     inla_logcdf_normal_fast((y - ((n_) * c1)) / (sqrt((n_)) * c2)))

	int id;
	GMRFLib_CACHE_SET_IDX(id);
	dtweedie_cache_tp *c = cache[id];

	double lambda = pow(mu, 2.0 - p) / (phi * (2.0 - p));
	double log_lambda = log(lambda);
	double alpha = (2.0 - p) / (p - 1.0);
	double gamma = phi * (p - 1.0) * pow(mu, p - 1.0);
	double plim = 0.999;
	double c1 = alpha * gamma;
	double c2 = sqrt(alpha) * gamma;
	double retval, prob, lprob, lprob_max, pacc, diff, lower_diff;

	// when sd/mean < low, we basically are in the Gaussian regime. this is for the case Gamma(n*alpha, gamma), so mean=
	// n*alpha*gamma, and sd= sqrt(n*alpha)*gamma, so gamma cancel.
	double low = 0.25;
	int n, n_gauss = (int) round(1.0 + 1.0 / (alpha * SQR(low)));

	// stride should scale with stdev which is sqrt(lambda)
	int stride = IMAX(1, (int) (0.707107 * sqrt(lambda)));
	int nfirst = 0;

	lower_diff = log(1.0E-6);
	lprob = LOG_PDF_POISSON(nfirst);
	lprob_max = LOG_PDF_POISSON(round(lambda));
	diff = lprob - lprob_max;
	retval = pacc = prob = exp(lprob);

	if (diff < lower_diff) {
		// find first pdf such that diff > lower_diff using a binary search
		int llow = 0;
		int hhigh = (int) (lambda - 4.0 * sqrt(lambda));
		hhigh = IMAX(nfirst + 1, hhigh);
		while (1) {
			int mmid = (llow + hhigh) / 2;
			diff = LOG_PDF_POISSON(mmid) - lprob_max;
			if (diff > lower_diff) {
				hhigh = mmid;
			} else {
				llow = mmid;
			}
			if (hhigh - llow <= 1)
				break;
		}
		nfirst = hhigh;
		lprob = LOG_PDF_POISSON(nfirst);
		pacc = prob = exp(lprob);
		retval = prob * CDF(nfirst);
	}
	// as we already have accounted for nfirst, we need to 'nfirst' to be the next one in the sum 
	nfirst++;

	if (stride == 1) {
		for (n = nfirst; pacc < plim; n++) {
			lprob += log_lambda - log(n);
			prob = exp(lprob);
			pacc += prob;
			retval += prob * CDF(n);
		}
	} else {
		// use interpolation of gamma_P() between each stride. use Rmath implementation
		// as it compute log(gamma_P()) directly
		double lcdf_left = LOG_CDF(nfirst);
		double lcdf_right;
		for (n = nfirst; pacc < plim; n += stride) {
			int nn = n + stride;
			lcdf_right = LOG_CDF(nn);
			for (int k = (n == nfirst ? 0 : 1); k <= stride; k++) {
				double w = k / (double) stride;
				double est = (1.0 - w) * lcdf_left + w * lcdf_right;
				lprob += log_lambda - log(n + k);
				pacc += exp(lprob);
				retval += exp(lprob + est);
			}
			lcdf_left = lcdf_right;
		}
	}

#undef LOG_PDF_POISSON
#undef CDF
#undef LOG_CDF
	return (retval);
}
