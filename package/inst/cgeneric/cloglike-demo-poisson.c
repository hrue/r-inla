#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <strings.h>

#include "cgeneric.h"

#define Calloc(n_, type_)  (type_ *)calloc((n_), sizeof(type_))
#define Malloc(n_, type_)  (type_ *)malloc((n_) * sizeof(type_))
#define SQR(x_) ((x_)*(x_))

double poisson_cdf(int y, double lambda)
{
	// simple standalone code. of course it is better to compute it using the incomplete Gamma-distribution/external library.
	double res = 0.0;
	if (y >= 0) {
		double p = exp(-lambda);
		res = p;
		for (int yy = 1; yy <= y; yy++) {
			p *= lambda / yy;
			res += p;
		}
	}
	return res;
}

double *inla_cloglike_poisson(inla_cloglike_cmd_tp cmd, double *theta,
			      inla_cgeneric_data_tp *data, int ny, double *y, int nx, double *x, double *result)
{
	double *ret = NULL;

	switch (cmd) {
	case INLA_CLOGLIKE_INITIAL:
	{
		ret = Malloc(1, double);
		ret[0] = 0;
	}
		break;

	case INLA_CLOGLIKE_LOG_PRIOR:
	{
		// not in use
		assert(0 == 1);
	}
		break;

	case INLA_CLOGLIKE_LOGLIKE:
	{
		// if ny=2, then the 2nd column is log(y[0]!), otherwise, we ignore the constant normalizing constant (which only
		// involves the marginal likelihood calculations)
		if (ny == 1) {
#pragma omp simd
			for (int i = 0; i < nx; i++) {
				result[i] = x[i] * y[0] - exp(x[i]);
			}
		} else if (ny == 2) {
#pragma omp simd
			for (int i = 0; i < nx; i++) {
				result[i] = x[i] * y[0] - exp(x[i]) - y[1];
			}
		} else {
			assert(0 == 1);
		}
	}
		break;

	case INLA_CLOGLIKE_CDF:
	{
		for (int i = 0; i < nx; i++) {
			result[i] = poisson_cdf(y[0], exp(x[i]));
		}
	}
		break;

	case INLA_CLOGLIKE_QUIT:
		break;
	}

	return (ret);
}

double *inla_cloglike_poisson_cache(inla_cloglike_cmd_tp cmd, double *theta,
				    inla_cgeneric_data_tp *data, int ny, double *y, int nx, double *x, double *result)
{
	// this code shows how to 'cache' expensive calculations, in this case log(y[0]!)

	typedef struct {
		int ymax;
		double *lfactorial;
	} Cache_tp;

	if (!(data->cache)) {
		// use a random name
#pragma omp critical (Name_e2814d0ff0cb393dee01d0eb049e6e976f56cce8)
		if (!(data->cache)) {
			Cache_tp *c = Malloc(1, Cache_tp);
			c->ymax = 1024;			       /* or something */
			c->lfactorial = Malloc(c->ymax + 1, double);
			c->lfactorial[0] = 0.0;
			for (int k = 1; k <= c->ymax; k++) {
				c->lfactorial[k] = c->lfactorial[k - 1] + log(k);
			}
			// IMPORTANT: need to assign data->cache as the last expression within this 'if(!(data->cache))',
			// otherwise, we can run into race-conditions
			*((Cache_tp **) (&data->cache)) = c;
		}
	}
	Cache_tp *cache = *((Cache_tp **) (&data->cache));

	double *ret = NULL;
	switch (cmd) {
	case INLA_CLOGLIKE_INITIAL:
	{
		// no hyperparameters
		ret = Malloc(1, double);
		ret[0] = 0;
	}
		break;

	case INLA_CLOGLIKE_LOG_PRIOR:
	{
		// not in use
		assert(0 == 1);
	}
		break;

	case INLA_CLOGLIKE_LOGLIKE:
	{
		// if y[0] is to large, we have to rebuild the cache...
		int iy = (int) y[0];
		assert(iy <= cache->ymax);
		double lfac = cache->lfactorial[iy];

		if (iy == 0) {
#pragma omp simd
			for (int i = 0; i < nx; i++) {
				result[i] = -exp(x[i]) - lfac;
			}
		} else {
#pragma omp simd
			for (int i = 0; i < nx; i++) {
				result[i] = x[i] * iy - exp(x[i]) - lfac;
			}
		}
	}
		break;

	case INLA_CLOGLIKE_CDF:
	{
		for (int i = 0; i < nx; i++) {
			result[i] = poisson_cdf(y[0], exp(x[i]));
		}
	}
		break;

	case INLA_CLOGLIKE_QUIT:
		break;
	}

	return (ret);
}
