#include <assert.h>
#include <float.h>
#include <math.h>
#include <omp.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/hashP.h"

int GMRFLib_ghq__intern(double *x, double *w, int n)
{
	int i, its, j, m;
	double p1, p2, p3, pp, z = 0, z1;
	double eps = 3.0e-14;				       /* shown to work fine for n < 199 */

	/*
	 * this function breaks down for n >= 199 
	 */
	GMRFLib_ASSERT(n < 199, GMRFLib_EINVARG);

	x--;						       /* indices from 1...n */
	if (w) {
		w--;					       /* in this routine */
	}

	m = (n + 1) / 2;
	for (i = 1; i <= m; i++) {
		switch (i) {
		case 1:
		{
			z = sqrt(2.0 * n + 1.0) - 1.85575 * pow(2.0 * n + 1.0, -0.16667);
		}
			break;
		case 2:
		{
			z -= 1.14 * pow((double) n, 0.426) / z;
		}
			break;
		case 3:
		{
			z = 1.86 * z - 0.86 * x[1];
		}
			break;
		case 4:
		{
			z = 1.91 * z - 0.91 * x[2];
		}
			break;
		default:
			z = 2.0 * z - x[i - 2];
		}
		its = 0;
		while (1) {
			its++;

			p1 = 0.7511255444649425;
			p2 = 0.0;
			for (j = 1; j <= n; j++) {
				p3 = p2;
				p2 = p1;
				p1 = z * sqrt(2.0 / j) * p2 - sqrt((j - 1.0) / j) * p3;
			}
			pp = sqrt(2.0 * n) * p2;
			z1 = z;
			z = z1 - p1 / pp;
			if (ABS(z - z1) <= eps)
				break;			       /* converged */
			GMRFLib_ASSERT(its <= 500, GMRFLib_ESNH);
		}

		x[i] = z;
		x[n + 1 - i] = -z;

		if (w) {
			w[i] = 2.0 / (pp * pp);
			w[n + 1 - i] = w[i];
		}
	}
	return 0;
}

int GMRFLib_ghq_abscissas(double **xp, int n)
{
	return GMRFLib_ghq(xp, NULL, n);
}

int GMRFLib_ghq_weights(double **wp, int n)
{
	return GMRFLib_ghq(NULL, wp, n);
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((target_clones(INLA_CLONE_TARGETS "default")))
int GMRFLib_ghq_ms(double **xp, double **wp, int n, double mean, double stdev)
{
	// the same for a given mean and stdev. Allocated new memory for xp and wp
	int i;
	double *xxp = NULL, *wwp = NULL;
	GMRFLib_ghq(&xxp, &wwp, n);
	assert(xxp);
	assert(wwp);

	*xp = aMalloc(n, double);
	*wp = aMalloc(n, double);
	assert(*xp);
	assert(*wp);
	for (i = 0; i < n; i++) {
		(*xp)[i] = xxp[i] * stdev + mean;
		(*wp)[i] = wwp[i];
	}
	return GMRFLib_SUCCESS;
}
#pragma GCC diagnostic pop

int GMRFLib_ghq(double **xp, double **wp, int n)
{
	/*
	 * return a pointer to the abscissas and the weights for the Gauss-Hermite quadrature for given `n' with kernel
	 * exp(-x^2/2)/sqrt(2*pi)
	 * 
	 * both abscissas and the weights are optional
	 * 
	 * an internal storage is used to store previously computed abscissas and weights 
	 */

	static map_ivp **abscissas = NULL;		       /* keep previous computed elements here */
	static map_ivp **weights = NULL;		       /* keep previous computed elements here */

	if (!abscissas) {
#pragma omp critical (Name_57dd787c76d8e98b908fbe47a9af2b183bc7a84a)
		if (!abscissas) {
			weights = Calloc(GMRFLib_CACHE_LEN(), map_ivp *);
			map_ivp **tmp = Calloc(GMRFLib_CACHE_LEN(), map_ivp *);
			abscissas = tmp;
		}
	}
	int idx = 0;
	GMRFLib_CACHE_SET_IDX(idx);

	if (!abscissas[idx]) {
#pragma omp critical (Name_144aa75e163ba7c9b9b2548a2c758b4aa2b19808)
		if (!abscissas[idx]) {
			map_ivp *tmp = Calloc(1, map_ivp);
			weights[idx] = Calloc(1, map_ivp);
			map_ivp_init(tmp);
			map_ivp_init(weights[idx]);
			abscissas[idx] = tmp;
		}
	}

	double *x = NULL, *w = NULL;
	void *ptr = NULL, *pptr = NULL;

	GMRFLib_ASSERT(n > 0, GMRFLib_EINVARG);
	if ((ptr = map_ivp_ptr(abscissas[idx], n))) {
		/*
		 * use previously computed. note that map_ivp_ptr returns a ptr to the stored ptr. 
		 */
		if (xp) {
			*xp = *((double **) ptr);
		}

		/*
		 * the weights should now be stored as well 
		 */
		if (wp) {
			pptr = map_ivp_ptr(weights[idx], n);
			GMRFLib_ASSERT(pptr, GMRFLib_ESNH);
			*wp = *((double **) pptr);
		}
	} else {
		/*
		 * compute new ones 
		 */

		// this storage is never free'd
		x = aMalloc(n, double);
		w = aMalloc(n, double);
		GMRFLib_ghq__intern(x, w, n);

		/*
		 * the Gauss-Hermite is with kernel exp(-x^2), transform to kernel exp(-x^2/2)/sqrt(2*pi)
		 */
		double s = 1.0 / sqrt(2.0 * M_PI);
		GMRFLib_dscale(n, M_SQRT2, x);
		GMRFLib_dscale(n, M_SQRT2 * s, w);

		/*
		 * reverse the order so its small to large, and sort the weights along 
		 */
		// GMRFLib_qsort2((void *) x, (size_t) n, sizeof(double), w, sizeof(double), NULL, 0, GMRFLib_dcmp);
		my_sort2_dd(x, w, n);
		map_ivp_set(abscissas[idx], n, (void *) x);
		map_ivp_set(weights[idx], n, (void *) w);

		if (xp) {
			*xp = x;			       /* return a ptr only */
		}
		if (wp) {
			*wp = w;			       /* return a ptr only */
		}
	}
	return GMRFLib_SUCCESS;
}
