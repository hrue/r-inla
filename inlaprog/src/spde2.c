
/* spde2.c
 * 
 * Copyright (C) 2011-2023  Havard Rue
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


#include <stdlib.h>
#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"
#include "inla.h"
#include "spde2.h"

extern G_tp G;						       /* import some global parametes from inla */

#define Sqr(x_) ((x_)*(x_))

double inla_spde2_Qfunction_orig(int thread_id, int ii, int jj, double *UNUSED(values), void *arg)
{
	if (jj < 0) {
		return NAN;
	}

	int i, j;
	if (ii <= jj) {
		i = ii;
		j = jj;
	} else {
		i = jj;
		j = ii;
	}

	inla_spde2_tp *model = (inla_spde2_tp *) arg;
	double value = 0.0;
	double d_i[3], d_j[3];

	int nc = model->B[0]->ncol;
	double *vals = GMRFLib_vmatrix_get(model->vmatrix, i, j);

	double *vals_i0 = vals;
	double *vals_i1 = vals + nc;
	double *vals_i2 = vals + 2 * nc;

	d_i[0] = vals_i0[0];
	d_i[1] = vals_i1[0];
	d_i[2] = vals_i2[0];

	if (i == j) {
#pragma GCC ivdep
		for (int k = 1; k < nc; k++) {
			double theta = model->theta[k - 1][thread_id][0];
			d_i[0] += vals_i0[k] * theta;
			d_i[1] += vals_i1[k] * theta;
			d_i[2] += vals_i2[k] * theta;
		}
#pragma GCC ivdep
		for (int k = 0; k < 2; k++) {
			d_i[k] = exp(d_i[k]);
		}

		if (model->transform != SPDE2_TRANSFORM_IDENTITY) {
			switch (model->transform) {
			case SPDE2_TRANSFORM_LOG:
			{
				d_i[2] = 2 * exp(d_i[2]) - 1.0;
			}
				break;

			case SPDE2_TRANSFORM_LOGIT:
			{
				d_i[2] = cos(M_PI * map_probability(d_i[2], MAP_FORWARD, NULL));
			}
				break;

			case SPDE2_TRANSFORM_IDENTITY:
				break;

			default:
				assert(0 == 1);
			}
		}

		double *v = vals + 3 * nc;
		value = SQR(d_i[0]) * (SQR(d_i[1]) * v[0] + d_i[2] * d_i[1] * (v[1] + v[2]) + v[3]);
	} else {
		double *vals_j0 = vals + 3 * nc;
		double *vals_j1 = vals + 4 * nc;
		double *vals_j2 = vals + 5 * nc;

		d_j[0] = vals_j0[0];
		d_j[1] = vals_j1[0];
		d_j[2] = vals_j2[0];

#pragma GCC ivdep
		for (int k = 1; k < nc; k++) {
			double theta = model->theta[k - 1][thread_id][0];
			d_i[0] += vals_i0[k] * theta;
			d_i[1] += vals_i1[k] * theta;
			d_i[2] += vals_i2[k] * theta;

			d_j[0] += vals_j0[k] * theta;
			d_j[1] += vals_j1[k] * theta;
			d_j[2] += vals_j2[k] * theta;
		}

#pragma GCC ivdep
		for (int k = 0; k < 2; k++) {
			d_i[k] = exp(d_i[k]);
			d_j[k] = exp(d_j[k]);
		}

		if (model->transform != SPDE2_TRANSFORM_IDENTITY) {
			switch (model->transform) {
			case SPDE2_TRANSFORM_LOG:
			{
				d_i[2] = 2 * exp(d_i[2]) - 1.0;
				d_j[2] = 2 * exp(d_j[2]) - 1.0;
			}
				break;

			case SPDE2_TRANSFORM_LOGIT:
			{
				d_i[2] = cos(M_PI * map_probability(d_i[2], MAP_FORWARD, NULL));
				d_j[2] = cos(M_PI * map_probability(d_j[2], MAP_FORWARD, NULL));
			}
				break;

			case SPDE2_TRANSFORM_IDENTITY:
				break;

			default:
				assert(0 == 1);
			}
		}

		double *v = vals + 6 * nc;
		value = d_i[0] * d_j[0] * (d_i[1] * d_j[1] * v[0] + d_i[2] * d_i[1] * v[1] + d_j[1] * d_j[2] * v[2] + v[3]);
	}

	return value;
}

double inla_spde2_Qfunction(int thread_id, int ii, int jj, double *values, void *arg)
{
	// use simple caching for this function. only cache calculations for one 'i', that can be used for all (i,j) with the same i. recall to
	// enable init in 'build_model below' before enable this function.

	if (jj < 0) {
		// it is benefitial do so a such approach, as we then ensure that 'i' is cached and we can bypass the computation of the
		// cache-index

		int idx = -1;
		GMRFLib_CACHE_SET_ID(idx);

		double fake_values[2];
		fake_values[0] = idx;			       /* transfer the cache-index */
		fake_values[1] = 1;			       /* tell that we know that 'i' is in the cache */

		inla_spde2_tp *model = (inla_spde2_tp *) arg;
		values[0] = inla_spde2_Qfunction(thread_id, ii, ii, fake_values, arg);
		for(int k = 0; k < model->graph->lnnbs[ii]; k++){
			values[1 + k] = inla_spde2_Qfunction(thread_id, ii, model->graph->lnbs[ii][k], fake_values, arg);
		}

		return 0.0;
	}

	const int debug = 0;
	const int debug_details = 0;

	int i, j;
	if (ii <= jj) {
		i = ii;
		j = jj;
	} else {
		i = jj;
		j = ii;
	}

	int idx = -1;
	int i_in_cache = 0; 

	if (values) {
		// this is a hack. these values are 'fake' and is part of the trick described above
		idx = (int) values[0];
		i_in_cache = (int) values[1];
	} else {
		GMRFLib_CACHE_SET_ID(idx);
	}

	inla_spde2_tp *model = (inla_spde2_tp *) arg;
	int nc = model->B[0]->ncol;
	double *vals = GMRFLib_vmatrix_get(model->vmatrix, i, j);

	if (!(model->cache[idx])) {
#pragma omp critical (Name_096287ed3ed383c234e780a2ee2897e5fede0116)
		{
			if (!(model->cache[idx])) {
				if (debug) {
					printf("spde2: init cache for idx = %1d\n", idx);
				}
				model->cache[idx] = Calloc(1, spde2_cache_tp);
				model->cache[idx]->i = -1;
				model->cache[idx]->need_transform = (model->transform != SPDE2_TRANSFORM_IDENTITY);

				double *work = Calloc(3 + nc, double);
				model->cache[idx]->theta = work;
				model->cache[idx]->vals = work + nc;
				// add theta[0] = 1.0 here, and append other 'thetas' after, so we can make cleaner loops
				model->cache[idx]->theta[0] = 1.0;
			}
		}
	}

	spde2_cache_tp *cache = model->cache[idx];
	int need_transform = cache->need_transform;

	double value;
	double d_i0 = 0.0;
	double d_i1 = 0.0;
	double d_i2 = 0.0;
	double d_j0 = 0.0;
	double d_j1 = 0.0;
	double d_j2 = 0.0;
	
	if (i == j) {
		double *vals_i0 = vals;
		double *vals_i1 = vals + nc;
		double *vals_i2 = vals + 2 * nc;

#pragma omp simd
		for (int k = 0; k < nc - 1; k++) {
			cache->theta[k + 1] = model->theta[k][thread_id][0];
		}
		double *theta = cache->theta;

#pragma omp simd reduction(+: d_i0, d_i1, d_i2)
		for (int k = 0; k < nc; k++) {
			double th = theta[k];
			d_i0 += vals_i0[k] * th;
			d_i1 += vals_i1[k] * th;
			d_i2 += vals_i2[k] * th;
		}
		d_i0 = exp(d_i0);
		d_i1 = exp(d_i1);

		if (need_transform) {
			switch (model->transform) {
			case SPDE2_TRANSFORM_IDENTITY:
				break;

			case SPDE2_TRANSFORM_LOG:
				d_i2 = 2 * exp(d_i2) - 1.0;
				break;

			case SPDE2_TRANSFORM_LOGIT:
				d_i2 = cos(M_PI * map_probability(d_i2, MAP_FORWARD, NULL));
				break;

			default:
				assert(0 == 1);
			}
		}

		double *v = vals + 3 * nc;
		value = Sqr(d_i0) * (Sqr(d_i1) * v[0] + d_i2 * d_i1 * (v[1] + v[2]) + v[3]);

		// store in cache. 'theta' is done already
		cache->vals[0] = d_i0;
		cache->vals[1] = d_i1;
		cache->vals[2] = d_i2;

		if (debug) {
#pragma omp critical (Name_81be0810e04979398ed477ac445e942f3056221b)
			{
				printf("spde2: store cache for idx=%1d i=%1d\n", idx, i);
				if (debug_details) {
					for (int k = 1; k < nc; k++) {
						printf("\ttheta[%1d] = %.12f\n", k, cache->theta[k]);
					}
					printf("\td_i0 = %.12f\n", d_i0);
					printf("\td_i1 = %.12f\n", d_i1);
					printf("\td_i2 = %.12f\n", d_i2);
				}
			}
		}
	} else {
		// check if we have the 'i' value in cache
		int in_cache = i_in_cache;

		// with this hack, we know its there if i_in_cache is TRUE, otherwise, we need to check
		if (!in_cache) {
			in_cache = (i == cache->i);
			if (in_cache) {
				for (int k = 0; k < nc - 1; k++) {
					if (cache->theta[1 + k] != model->theta[k][thread_id][0]) {
						in_cache = 0;
						break;
					}
				}
			}
		}
		
		if (debug) {
#pragma omp critical (Name_170ba04061977b3f6655c4b7b4bb6c086b3a7c68)
			{
				if (in_cache) {
					printf("spde2: use cache for idx=%1d i=%1d j=%1d\n", idx, i, j);
					if (debug_details) {
						for (int k = 1; k < nc; k++) {
							printf("\ttheta[%1d] = %.12f\n", k, cache->theta[k]);
						}
						printf("\td_i0 = %.12f\n", cache->vals[0]);
						printf("\td_i1 = %.12f\n", cache->vals[1]);
						printf("\td_i2 = %.12f\n", cache->vals[2]);
					}
				} else {
					printf("spde2: not in cache for idx = %1d, i = %1d, j = %1d\n", idx, i, j);
				}
			}
		}
		// check hit/miss rates... might be useful in the future again, so I keep it here
		if (debug_details) {
			static double cache_hit = 0.0;
			static double cache_miss = 0.0;
			static double cache_count = 0.0;
#pragma omp critical (Name_8a2bd7f139bfd60665130abf233f21903ff83ad9)
			{
				cache_count++;
				cache_hit += in_cache;
				cache_miss += (1 - in_cache);

				if (!((int) cache_count % 100000))
					printf("hit %.3f miss %.3f\n", cache_hit / (cache_hit + cache_miss), cache_miss / (cache_hit + cache_miss));
			}
		}

		double *vals_j0 = vals + 3 * nc;
		double *vals_j1 = vals + 4 * nc;
		double *vals_j2 = vals + 5 * nc;

		if (in_cache) {
			d_i0 = cache->vals[0];
			d_i1 = cache->vals[1];
			d_i2 = cache->vals[2];

			double *theta = cache->theta;

#pragma omp simd reduction(+: d_j0, d_j1, d_j2)
			for (int k = 0; k < nc; k++) {
				double th = theta[k];
				d_j0 += vals_j0[k] * th;
				d_j1 += vals_j1[k] * th;
				d_j2 += vals_j2[k] * th;
			}

			d_j0 = exp(d_j0);
			d_j1 = exp(d_j1);

			if (need_transform) {
				switch (model->transform) {
				case SPDE2_TRANSFORM_IDENTITY:
					break;

				case SPDE2_TRANSFORM_LOG:
					d_j2 = 2.0 * exp(d_j2) - 1.0;
					break;

				case SPDE2_TRANSFORM_LOGIT:
					d_j2 = cos(M_PI * map_probability(d_j2, MAP_FORWARD, NULL));
					break;

				default:
					assert(0 == 1);
				}
			}

		} else {
			// not in_cache
			double *vals_i0 = vals;
			double *vals_i1 = vals + nc;
			double *vals_i2 = vals + 2 * nc;

			for (int k = 0; k < nc - 1; k++) {
				cache->theta[1 + k] = model->theta[k][thread_id][0];
			}
			double *theta = cache->theta;

#pragma omp simd reduction(+: d_i0, d_i1, d_i2)
			for (int k = 0; k < nc; k++) {
				double th = theta[k];
				d_i0 += vals_i0[k] * th;
				d_i1 += vals_i1[k] * th;
				d_i2 += vals_i2[k] * th;
			}
#pragma omp simd reduction(+: d_j0, d_j1, d_j2)
			for (int k = 0; k < nc; k++) {
				double th = theta[k];
				d_j0 += vals_j0[k] * th;
				d_j1 += vals_j1[k] * th;
				d_j2 += vals_j2[k] * th;
			}
			d_i0 = exp(d_i0);
			d_i1 = exp(d_i1);
			d_j0 = exp(d_j0);
			d_j1 = exp(d_j1);

			if (need_transform) {
				switch (model->transform) {
				case SPDE2_TRANSFORM_IDENTITY:
					break;

				case SPDE2_TRANSFORM_LOG:
					d_i2 = 2.0 * exp(d_i2) - 1.0;
					d_j2 = 2.0 * exp(d_j2) - 1.0;
					break;

				case SPDE2_TRANSFORM_LOGIT:
					d_i2 = cos(M_PI * map_probability(d_i2, MAP_FORWARD, NULL));
					d_j2 = cos(M_PI * map_probability(d_j2, MAP_FORWARD, NULL));
					break;

				default:
					assert(0 == 1);
				}
			}
		}

		double *v = vals + 6 * nc;
		value = d_i0 * d_j0 * (d_i1 * d_j1 * v[0] + d_i2 * d_i1 * v[1] + d_j1 * d_j2 * v[2] + v[3]);

		if (!in_cache) {
			// cache this value
			cache->i = i;
			cache->vals[0] = d_i0;
			cache->vals[1] = d_i1;
			cache->vals[2] = d_i2;

			if (debug) {
#pragma omp critical (Name_ad1bbe7257f1c3db0c64e992dc63fdd484fda64d)
				{
					printf("spde2: store cache for idx=%1d i=%1d\n", idx, i);
					if (debug_details) {
						for (int k = 1; k < nc; k++) {
							printf("\ttheta[%1d] = %.12f\n", k, cache->theta[k]);
						}
						printf("\td_i0 = %.12f\n", cache->vals[0]);
						printf("\td_i1 = %.12f\n", cache->vals[1]);
						printf("\td_i2 = %.12f\n", cache->vals[2]);
					}
				}
			}
		}
	}

	return value;
}

double inla_spde2_Qfunction_old(int thread_id, int ii, int jj, double *UNUSED(values), void *arg)
{
	if (jj < 0) {
		return NAN;
	}

	int i, j;
	if (ii <= jj) {
		i = ii;
		j = jj;
	} else {
		i = jj;
		j = ii;
	}

	inla_spde2_tp *model = (inla_spde2_tp *) arg;
	int nc = model->B[0]->ncol;

	// manual inline
	// double *vals = GMRFLib_vmatrix_get(model->vmatrix, i, j);
	double *vals = (double *) *map_ivp_ptr(&(model->vmatrix->vmat[i]), j);


	double value;
	double d_i[6];
	double *d_j = d_i + 3;

	if (i == j) {

#pragma omp simd
		for (int k = 0; k < 3; k++) {
			d_i[k] = vals[k * nc];
		}
		for (int k = 1; k < nc; k++) {
			double th = model->theta[k - 1][thread_id][0];
			double *v = vals + k;
#pragma omp simd
			for (int kk = 0; kk < 3; kk++) {
				d_i[kk] += v[kk * nc] * th;
			}
		}

#pragma omp simd
		for (int k= 0; k < 2; k++) {
			d_i[k] = exp(d_i[k]);
		}

		if (model->transform != SPDE2_TRANSFORM_IDENTITY) {
			switch (model->transform) {
			case SPDE2_TRANSFORM_IDENTITY:
				break;

			case SPDE2_TRANSFORM_LOG:
			{
				d_i[2] = 2 * exp(d_i[2]) - 1.0;
			}
				break;

			case SPDE2_TRANSFORM_LOGIT:
			{
				d_i[2] = cos(M_PI * map_probability(d_i[2], MAP_FORWARD, NULL));
			}
				break;

			default:
				assert(0 == 1);
			}
		}

		double *v = vals + 3 * nc;
		value = Sqr(d_i[0]) * (Sqr(d_i[1]) * v[0] + d_i[2] * d_i[1] * (v[1] + v[2]) + v[3]);

	} else {
#pragma omp simd
		for (int k = 0; k < 6; k++) {
			d_i[k] = vals[k * nc];
		}

		for (int k = 1; k < nc; k++) {
			double th = model->theta[k - 1][thread_id][0];
			double *v = vals + k;
#pragma omp simd
			for (int kk = 0; kk < 6; kk++) {
				d_i[kk] += v[kk * nc] * th;
			}
		}

#pragma omp simd
		for (int k = 0; k < 2; k++) {
			d_i[k] = exp(d_i[k]);
			d_j[k] = exp(d_j[k]);
		}

		if (model->transform != SPDE2_TRANSFORM_IDENTITY) {
			switch (model->transform) {
			case SPDE2_TRANSFORM_IDENTITY:
				break;

			case SPDE2_TRANSFORM_LOG:
			{
				d_i[2] = 2.0 * exp(d_i[2]) - 1.0;
				d_j[2] = 2.0 * exp(d_j[2]) - 1.0;
			}
				break;

			case SPDE2_TRANSFORM_LOGIT:
			{
				d_i[2] = cos(M_PI * map_probability(d_i[2], MAP_FORWARD, NULL));
				d_j[2] = cos(M_PI * map_probability(d_j[2], MAP_FORWARD, NULL));
			}
				break;

			default:
				assert(0 == 1);
			}
		}

		double *v = vals + 6 * nc;
		value = d_i[0] * d_j[0] * (d_i[1] * d_j[1] * v[0] + d_i[2] * d_i[1] * v[1] + d_j[1] * d_j[2] * v[2] + v[3]);
	}

	return value;
}

int inla_spde2_build_model(int UNUSED(thread_id), inla_spde2_tp **smodel, const char *prefix, const char *transform)
{
	int i;
	const int debug = 0;
	inla_spde2_tp *model = NULL;
	char *fnm = NULL;

	model = Calloc(1, inla_spde2_tp);

	if (strcasecmp(transform, "logit") == 0) {
		model->transform = SPDE2_TRANSFORM_LOGIT;
	} else if (strcasecmp(transform, "log") == 0) {
		model->transform = SPDE2_TRANSFORM_LOG;
	} else if (strcasecmp(transform, "identity") == 0) {
		model->transform = SPDE2_TRANSFORM_IDENTITY;
	} else {
		assert(0 == 1);
	}

	model->B = Calloc(3, GMRFLib_matrix_tp *);
	model->M = Calloc(3, GMRFLib_matrix_tp *);
	for (i = 0; i < 3; i++) {
		GMRFLib_sprintf(&fnm, "%s%s%1d", prefix, "B", i);
		model->B[i] = GMRFLib_read_fmesher_file((const char *) fnm, 0, -1);

		GMRFLib_sprintf(&fnm, "%s%s%1d", prefix, "M", i);
		model->M[i] = GMRFLib_read_fmesher_file((const char *) fnm, 0, -1);
	}

	for (i = 1; i < 3; i++) {
		/*
		 * all need the same dimensions n x (p+1) 
		 */
		assert(model->B[0]->nrow == model->B[i]->nrow);
		assert(model->B[0]->ncol == model->B[i]->ncol);

		/*
		 * all are square with the same dimension n x n 
		 */
		assert(model->M[i]->nrow == model->M[i]->ncol);
		assert(model->M[0]->nrow == model->M[i]->nrow);

		/*
		 * and the number of rows must be the same 
		 */
		assert(model->B[i]->nrow == model->M[i]->nrow);
	}

	model->n = model->M[0]->nrow;
	model->ntheta = model->B[0]->ncol - 1;

	GMRFLib_sprintf(&fnm, "%s%s", prefix, "BLC");
	if (GMRFLib_is_fmesher_file((const char *) fnm, 0L, -1) == GMRFLib_SUCCESS) {
		model->BLC = GMRFLib_read_fmesher_file((const char *) fnm, 0, -1);
	} else {
		model->BLC = NULL;
	}

	if (debug) {
		P(model->n);
		P(model->ntheta);
		P(model->ntheta_used);
	}

	/*
	 * I need to build the graph. Need to add both M_ij and M_ji as M1 can be non-symmetric. 
	 */
	GMRFLib_ged_tp *ged = NULL;
	GMRFLib_ged_init(&ged, NULL);

#define ADD_GRAPH(_G)							\
	{								\
		int i_, j_, jj_;					\
		for(i_ = 0; i_ < _G->n; i_++) {				\
			GMRFLib_ged_add(ged, i_, i_);			\
			for(jj_ = 0; jj_ < _G->nnbs[i_]; jj_++){	\
				j_ = _G->nbs[i_][jj_];			\
				GMRFLib_ged_add(ged, i_, j_);		\
				GMRFLib_ged_add(ged, j_, i_);		\
			}						\
		}							\
	}

	for (i = 0; i < 3; i++) {
		ADD_GRAPH(model->M[i]->graph);
	}

#undef ADD_GRAPH

	GMRFLib_ged_build(&(model->graph), ged);
	assert(model->n == model->graph->n);
	GMRFLib_ged_free(ged);

	model->Qfunc = inla_spde2_Qfunction;
	model->Qfunc_arg = (void *) model;

	HYPER_NEW2(model->theta, 0.0, model->ntheta);
	*smodel = model;

	// add better storage
	int nc = model->B[0]->ncol;
	GMRFLib_vmatrix_init(&(model->vmatrix), model->n, model->graph);
	for (i = 0; i < model->n; i++) {
		int j = i;

		double *v = Calloc(3 * nc + 4, double);
		assert(v);

		GMRFLib_matrix_get_row(v + 0 * nc, i, model->B[0]);
		GMRFLib_matrix_get_row(v + 1 * nc, i, model->B[1]);
		GMRFLib_matrix_get_row(v + 2 * nc, i, model->B[2]);

		double *vv = v + 3 * nc;
		vv[0] = GMRFLib_matrix_get(i, j, model->M[0]);
		vv[1] = GMRFLib_matrix_get(i, j, model->M[1]);
		vv[2] = GMRFLib_matrix_get(j, i, model->M[1]);
		vv[3] = GMRFLib_matrix_get(i, j, model->M[2]);

		GMRFLib_vmatrix_set(model->vmatrix, i, j, v);

		for (int jj = 0; jj < model->graph->lnnbs[i]; jj++) {
			j = model->graph->lnbs[i][jj];
			v = Calloc(6 * nc + 4, double);
			assert(v);

			GMRFLib_matrix_get_row(v + 0 * nc, i, model->B[0]);
			GMRFLib_matrix_get_row(v + 1 * nc, i, model->B[1]);
			GMRFLib_matrix_get_row(v + 2 * nc, i, model->B[2]);

			GMRFLib_matrix_get_row(v + 3 * nc, j, model->B[0]);
			GMRFLib_matrix_get_row(v + 4 * nc, j, model->B[1]);
			GMRFLib_matrix_get_row(v + 5 * nc, j, model->B[2]);

			vv = v + 6 * nc;
			vv[0] = GMRFLib_matrix_get(i, j, model->M[0]);
			vv[1] = GMRFLib_matrix_get(i, j, model->M[1]);
			vv[2] = GMRFLib_matrix_get(j, i, model->M[1]);
			vv[3] = GMRFLib_matrix_get(i, j, model->M[2]);

			GMRFLib_vmatrix_set(model->vmatrix, i, j, v);
		}
	}

	// recall to enable this if using cache
	model->cache = Calloc(GMRFLib_MAX_THREADS(), spde2_cache_tp *);

	return INLA_OK;
}

double *inla_spde2_userfunc2(int number, double *theta, int nhyper, double *covmat, void *arg)
{
	/*
	 * compute the marginals for BLC*theta
	 */

	/*
	 * Cov_spde2: Covariance between the local theta's in this spde2 model. Cov: covariance between all theta's in the full model. Theta_spde2: the theta's in
	 * this spde2 model. Theta: the theta's in the full model. 
	 */
#define Cov_spde2(_i, _j) covmat[(idx_offset + (_i)) + nhyper * (idx_offset + (_j))]
#define Cov(_i, _j)       covmat[(_i) + nhyper * (_j)]
#define Theta_spde2(_i)   theta[idx_offset + (_i)]
#define Theta(_i)         theta[(_i)]

	// these are the useful new ones
#define CovNew(_i, _j) covmat_new[(_i) + nhyper_new * (_j)]
#define ThetaNew(_i)   theta_new[(_i)]

	if (!covmat || nhyper == 0) {
		return NULL;
	}

	int use_new_version = 1, i, ii, jj, k, kk, *idx_map = NULL, nhyper_new = 0;
	const int debug = 0;
	double *theta_new = NULL, *covmat_new = NULL, *stdev_corr_pos_new = NULL, *stdev_corr_neg_new = NULL;
	gsl_vector *sqrt_eigen_values_new = NULL;
	gsl_matrix *eigen_vectors_new = NULL;
	GMRFLib_userfunc2_arg_tp *a = (GMRFLib_userfunc2_arg_tp *) arg;
	inla_spde2_tp *model = (inla_spde2_tp *) GMRFLib_ai_INLA_userfunc2_args[number];

	if (model->BLC == NULL) {
		return NULL;
	}

	int idx_offset = model->theta_first_idx;
	int nrow = model->BLC->nrow;
	int ncol = model->BLC->ncol;
	double *row_spde2 = Calloc(ncol, double);	       /* yes, one row has length ncol. */
	double *row = NULL;				       /* set later */
	if (debug) {
		P(nhyper);
		P(model->ntheta);
		P(model->ntheta_used);
	}

	assert(ncol - 1 == model->ntheta);

	/*
	 * build new theta vector and covariance matrix by inserting fixed values
	 */
	nhyper_new = nhyper + model->ntheta - model->ntheta_used;
	idx_map = Calloc(nhyper_new, int);
	for (k = kk = 0; k < nhyper_new; k++) {
		if ((k >= idx_offset) && (k < idx_offset + model->ntheta) && model->fixed[k - idx_offset]) {
			idx_map[k] = -1;
		} else {
			idx_map[k] = kk;
			kk++;
		}
	}

	theta_new = Calloc(nhyper_new, double);
	stdev_corr_pos_new = Calloc(nhyper_new, double);
	stdev_corr_neg_new = Calloc(nhyper_new, double);
	sqrt_eigen_values_new = gsl_vector_alloc(nhyper_new);
	for (k = kk = 0; k < nhyper_new; k++) {
		if (idx_map[k] >= 0) {
			theta_new[k] = theta[idx_map[k]];
			stdev_corr_pos_new[k] = a->stdev_corr_pos[idx_map[k]];
			stdev_corr_neg_new[k] = a->stdev_corr_neg[idx_map[k]];
			gsl_vector_set(sqrt_eigen_values_new, k, gsl_vector_get(a->sqrt_eigen_values, idx_map[k]));
		} else {
			theta_new[k] = model->fixed_values[kk];
			if (ISNAN(theta_new[k])) {
				// otherwise other things will break
				theta_new[k] = 0.0;
			}
			stdev_corr_pos_new[k] = 1.0;
			stdev_corr_neg_new[k] = 1.0;
			gsl_vector_set(sqrt_eigen_values_new, k, 0.0);
			kk++;
		}
	}

	covmat_new = Calloc(ISQR(nhyper_new), double);
	eigen_vectors_new = gsl_matrix_calloc(nhyper_new, nhyper_new);
	for (k = 0; k < nhyper_new; k++) {
		for (kk = 0; kk < nhyper_new; kk++) {
			if (idx_map[k] >= 0 && idx_map[kk] >= 0) {
				CovNew(k, kk) = Cov(idx_map[k], idx_map[kk]);
				gsl_matrix_set(eigen_vectors_new, k, kk, gsl_matrix_get(a->eigen_vectors, idx_map[k], idx_map[kk]));
			} else {
				CovNew(k, kk) = 0.0;
				gsl_matrix_set(eigen_vectors_new, k, kk, 0.0);
			}
		}
	}

	if (debug) {
		printf("inla_spde2_userfunc2:\n");
		for (k = 0; k < nhyper_new; k++) {
			printf("\tidx_map[%1d] = %d\n", k, idx_map[k]);
		}
		for (k = 0; k < nhyper; k++) {
			printf("\ttheta[%1d] = %g\n", k, theta[k]);
		}
		for (k = 0; k < nhyper_new; k++) {
			printf("\ttheta_new[%1d] = %g\n", k, theta_new[k]);
		}
		printf("\tcovmat\n");
		for (k = 0; k < nhyper; k++) {
			printf("\t\t");
			for (kk = 0; kk < nhyper; kk++) {
				printf(" %8.4f", Cov(k, kk));
			}
			printf("\n");
		}
		printf("\n");
		printf("\tcovmat_new\n");
		for (k = 0; k < nhyper_new; k++) {
			printf("\t\t");
			for (kk = 0; kk < nhyper_new; kk++) {
				printf(" %8.4f", CovNew(k, kk));
			}
			printf("\n");
		}
		printf("\n");
	}

	if (!GMRFLib_ai_INLA_userfunc2_len) {
		GMRFLib_ai_INLA_userfunc2_len = Calloc(GMRFLib_ai_INLA_userfunc2_n, int);
	}
	GMRFLib_ai_INLA_userfunc2_len[number] = nrow;
	if (!GMRFLib_ai_INLA_userfunc2_density) {
		GMRFLib_ai_INLA_userfunc2_density = Calloc(GMRFLib_ai_INLA_userfunc2_n, GMRFLib_density_tp **);
	}
	GMRFLib_ai_INLA_userfunc2_density[number] = Calloc(GMRFLib_ai_INLA_userfunc2_len[number], GMRFLib_density_tp *);

	row = Calloc(nhyper_new + 1, double);
	if (use_new_version) {
		for (i = 0; i < nrow; i++) {
			/*
			 * insert it into a larger vector to get the the full 'row' 
			 */
			GMRFLib_matrix_get_row(row_spde2, i, model->BLC);
			Memset(row, 0, (nhyper_new + 1) * sizeof(double));
			row[0] = row_spde2[0];
			Memcpy(row + idx_offset + 1, row_spde2 + 1, model->ntheta * sizeof(double));

			/*
			 * Sigma * a, a = row
			 */
			double *Sigma_a = Calloc(nhyper_new, double);
			for (ii = 0; ii < nhyper_new; ii++) {
				for (jj = 0; jj < nhyper_new; jj++) {
					Sigma_a[ii] += CovNew(ii, jj) * row[1 + jj];
				}
			}

			/*
			 * Get the mean and the variance using the covariance and the mode 
			 */
			double mean = row[0];
			double var = 0.0;
			for (ii = 0; ii < nhyper_new; ii++) {
				mean += ThetaNew(ii) * row[1 + ii];
				var += Sigma_a[ii] * row[1 + ii];
			}

			GMRFLib_ai_integrator_arg_tp *iarg = Calloc(1, GMRFLib_ai_integrator_arg_tp);

			iarg->nhyper = nhyper_new;
			iarg->idx = -1;
			iarg->return_log = GMRFLib_TRUE;       /* return log(density), yes. */
			iarg->hyper_count = -1;
			iarg->hyper_z = NULL;
			iarg->hyper_ldens = NULL;
			iarg->theta_mode = &ThetaNew(0);
			iarg->sqrt_eigen_values = sqrt_eigen_values_new;
			iarg->eigen_vectors = eigen_vectors_new;
			iarg->z = Calloc(nhyper_new, double);
			iarg->theta = Calloc(nhyper_new, double);
			iarg->stdev_corr_pos = stdev_corr_pos_new;
			iarg->stdev_corr_neg = stdev_corr_neg_new;
			iarg->dz = -1;
			iarg->interpolator = GMRFLib_AI_INTERPOLATOR_CCD;

			int npoints = 51;
			double *x = Calloc(nhyper_new, double), *xx = NULL, *xxx = Calloc(npoints, double), *ldens_values = Calloc(npoints, double);

			GMRFLib_ghq_abscissas(&xx, npoints);
			Memcpy(xxx, xx, npoints * sizeof(double));
			xxx[0] = DMIN(xxx[0], -GMRFLib_DENSITY_INTEGRATION_LIMIT);
			xxx[npoints - 1] = DMAX(xxx[npoints - 1], GMRFLib_DENSITY_INTEGRATION_LIMIT);

			if (var > 0.0) {
				for (ii = 0; ii < npoints; ii++) {
					for (jj = 0; jj < nhyper_new; jj++) {
						x[jj] = ThetaNew(jj) + Sigma_a[jj] * xxx[ii] / sqrt(var);
					}
					ldens_values[ii] = GMRFLib_ai_integrator_func(nhyper_new, x, iarg);
				}
				GMRFLib_density_create(&(GMRFLib_ai_INLA_userfunc2_density[number][i]),
						       GMRFLib_DENSITY_TYPE_SCGAUSSIAN, npoints, xxx, ldens_values, mean, sqrt(var), GMRFLib_TRUE);
			} else {
				// return a point-mass
				GMRFLib_density_create_normal(&(GMRFLib_ai_INLA_userfunc2_density[number][i]), 0.0, 1.0, mean, DBL_EPSILON,
							      GMRFLib_TRUE);
			}
			Free(Sigma_a);
			Free(x);
			Free(xxx);
			Free(ldens_values);
			Free(iarg->z);
			Free(iarg->theta);
			Free(iarg);
		}
	} else {
		// Not yet modified to deal with the case there ntheta_used < ntheta. hrue 26/7/2016.
		assert(model->ntheta == model->ntheta_used);
		for (i = 0; i < nrow; i++) {
			double mean = 0.0, var = 0.0;

			GMRFLib_matrix_get_row(row_spde2, i, model->BLC);
			mean = row_spde2[0];
			for (ii = 0; ii < model->ntheta; ii++) {
				mean += Theta_spde2(ii) * row_spde2[1 + ii];
				for (jj = 0; jj < model->ntheta; jj++) {
					var += row_spde2[1 + ii] * row_spde2[1 + jj] * Cov_spde2(ii, jj);	/* yes the first column is a
														 * constant offset */
				}
			}
			GMRFLib_density_create_normal(&(GMRFLib_ai_INLA_userfunc2_density[number][i]), 0.0, 1.0, mean,
						      (var > 0 ? sqrt(var) : DBL_EPSILON), GMRFLib_TRUE);
		}
	}

	Free(row_spde2);
	Free(row);
	Free(idx_map);
	Free(covmat_new);
	Free(theta_new);
	Free(stdev_corr_neg_new);
	Free(stdev_corr_pos_new);
	gsl_vector_free(sqrt_eigen_values_new);
	gsl_matrix_free(eigen_vectors_new);

#undef Cov_spde2
#undef Cov
#undef Theta_spde2
#undef Theta
#undef CovNew
#undef ThetaNew
	return NULL;
}

#undef Sqr
