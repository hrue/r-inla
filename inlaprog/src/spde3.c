
/* spde3.c
 * 
 * Copyright (C) 2014  Havard Rue
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
static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;


#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"
#include "inla.h"
#include "spde3.h"


extern G_tp G;						       /* import some global parametes from inla */

double inla_spde3_Qfunction(int i, int j, void *arg)
{
	inla_spde3_tp *model = (inla_spde3_tp *) arg;
	double value, phi_i[3], phi_j[3], d_i[3], d_j[3];
	int k, kk;

	/*
	 * to hold the i'th and j'th row of the B-matrices. use one storage only 
	 */
	double *row_i = Calloc(2 * model->B[0]->ncol, double);
	double *row_j = &row_i[model->B[0]->ncol];

	for (k = 0; k < 3; k++) {
		if (i == j) {
			/*
			 * some savings for i == j 
			 */
			GMRFLib_matrix_get_row(row_i, i, model->B[k]);
			phi_i[k] = row_i[0];
			for (kk = 1; kk < model->B[k]->ncol; kk++) {
				/*
				 * '-1' is the correction for the first intercept column in B 
				 */
				phi_i[k] += row_i[kk] * model->theta[kk - 1][GMRFLib_thread_id][0];
			}
			phi_j[k] = phi_i[k];		       /* they are equal in this case */
		} else {
			/*
			 * i != j 
			 */
			GMRFLib_matrix_get_row(row_i, i, model->B[k]);
			GMRFLib_matrix_get_row(row_j, j, model->B[k]);
			phi_i[k] = row_i[0];
			phi_j[k] = row_j[0];
			for (kk = 1; kk < model->B[k]->ncol; kk++) {
				/*
				 * '-1' is the correction for the first intercept column in B 
				 */
				phi_i[k] += row_i[kk] * model->theta[kk - 1][GMRFLib_thread_id][0];
				phi_j[k] += row_j[kk] * model->theta[kk - 1][GMRFLib_thread_id][0];
			}
		}
	}
	Free(row_i);

	for (k = 0; k < 2; k++) {
		d_i[k] = exp(phi_i[k]);
		d_j[k] = exp(phi_j[k]);
	}

	/*
	 * change this later on, need an option here for various 'link' functions. some savings possible for i==j.
	 */
	if (i == j) {
		switch (model->transform) {
		case SPDE3_TRANSFORM_LOGIT:
			d_i[2] = cos(M_PI * map_probability(phi_i[2], MAP_FORWARD, NULL));
			break;
		case SPDE3_TRANSFORM_LOG:
			d_i[2] = 2 * exp(phi_i[2]) - 1.0;
			break;
		case SPDE3_TRANSFORM_IDENTITY:
			d_i[2] = phi_i[2];
			break;
		default:
			assert(0 == 1);
		}
		d_j[2] = d_i[2];
	} else {
		switch (model->transform) {
		case SPDE3_TRANSFORM_LOGIT:
			d_i[2] = cos(M_PI * map_probability(phi_i[2], MAP_FORWARD, NULL));
			d_j[2] = cos(M_PI * map_probability(phi_j[2], MAP_FORWARD, NULL));
			break;
		case SPDE3_TRANSFORM_LOG:
			d_i[2] = 2 * exp(phi_i[2]) - 1.0;
			d_j[2] = 2 * exp(phi_j[2]) - 1.0;
			break;
		case SPDE3_TRANSFORM_IDENTITY:
			d_i[2] = phi_i[2];
			d_j[2] = phi_j[2];
			break;
		default:
			assert(0 == 1);
		}
	}

	value = d_i[0] * d_j[0] * (d_i[1] * d_j[1] * GMRFLib_matrix_get(i, j, model->M[0]) +
				   d_i[2] * d_i[1] * GMRFLib_matrix_get(i, j, model->M[1]) +
				   d_j[1] * d_j[2] * GMRFLib_matrix_get(j, i, model->M[1]) + GMRFLib_matrix_get(i, j, model->M[2]));

	return value;
}

int inla_spde3_build_model(inla_spde3_tp ** smodel, const char *prefix, const char *transform)
{
	int i, debug = 0;
	inla_spde3_tp *model = NULL;
	char *fnm = NULL;

	model = Calloc(1, inla_spde3_tp);

	if (strcasecmp(transform, "logit") == 0) {
		model->transform = SPDE3_TRANSFORM_LOGIT;
	} else if (strcasecmp(transform, "log") == 0) {
		model->transform = SPDE3_TRANSFORM_LOG;
	} else if (strcasecmp(transform, "identity") == 0) {
		model->transform = SPDE3_TRANSFORM_IDENTITY;
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

	model->Qfunc = inla_spde3_Qfunction;
	model->Qfunc_arg = (void *) model;

	HYPER_NEW2(model->theta, 0.0, model->ntheta);
	*smodel = model;

	return INLA_OK;
}

int inla_spde3_userfunc3(int number, double *theta, int nhyper, double *covmat, void *arg)
{
	/*
	 * compute the marginals for BLC*theta
	 */

	GMRFLib_userfunc3_arg_tp *a = (GMRFLib_userfunc3_arg_tp *) arg;

	/*
	 * Cov_spde3: Covariance between the local theta's in this spde3 model. Cov: covariance between all theta's in the full model. Theta_spde3: the theta's in
	 * this spde3 model. Theta: the theta's in the full model. 
	 */
#define Cov_spde3(_i, _j) covmat[(idx_offset + (_i)) + nhyper * (idx_offset + (_j))]
#define Cov(_i, _j)       covmat[(_i) + nhyper * (_j)]
#define Theta_spde3(_i)   theta[idx_offset + (_i)]
#define Theta(_i)         theta[(_i)]

	if (!covmat || nhyper == 0) {
		return INLA_OK;
	}

	int use_new_version = 1;
	int i, ii, jj;
	inla_spde3_tp *model = (inla_spde3_tp *) GMRFLib_ai_INLA_userfunc3_args[number];

	if (model->BLC == NULL) {
		return INLA_OK;
	}

	int idx_offset = model->theta_first_idx;
	int nrow = model->BLC->nrow;
	int ncol = model->BLC->ncol;
	double *row_spde3 = Calloc(ncol, double);	       /* yes, one row has length ncol. */
	double *row = Calloc(nhyper + 1, double);

	assert(ncol - 1 == model->ntheta);

	if (!GMRFLib_ai_INLA_userfunc3_len) {
		GMRFLib_ai_INLA_userfunc3_len = Calloc(GMRFLib_ai_INLA_userfunc3_n, int);
	}
	GMRFLib_ai_INLA_userfunc3_len[number] = nrow;

	if (!GMRFLib_ai_INLA_userfunc3_density) {
		GMRFLib_ai_INLA_userfunc3_density = Calloc(GMRFLib_ai_INLA_userfunc3_n, GMRFLib_density_tp **);
	}
	GMRFLib_ai_INLA_userfunc3_density[number] = Calloc(nrow, GMRFLib_density_tp *);

	if (use_new_version) {
		for (i = 0; i < nrow; i++) {

			GMRFLib_matrix_get_row(row_spde3, i, model->BLC);

			/*
			 * insert it into a larger vector to get the the full 'row' 
			 */
			memset(row, 0, (nhyper + 1) * sizeof(double));
			row[0] = row_spde3[0];
			memcpy(row + idx_offset + 1, row_spde3 + 1, model->ntheta * sizeof(double));

			/*
			 * Sigma * a, a = row
			 */
			double *Sigma_a = Calloc(nhyper, double);
			for (ii = 0; ii < nhyper; ii++) {
				for (jj = 0; jj < nhyper; jj++) {
					Sigma_a[ii] += Cov(ii, jj) * row[1 + jj];
				}
			}

			/*
			 * Get the mean and the variance using the covariance and the mode 
			 */
			double mean = row[0];
			double var = 0.0;
			for (ii = 0; ii < nhyper; ii++) {
				mean += Theta(ii) * row[1 + ii];
				var += Sigma_a[ii] * row[1 + ii];
			}

			GMRFLib_ai_integrator_arg_tp *iarg = Calloc(1, GMRFLib_ai_integrator_arg_tp);

			iarg->nhyper = nhyper;
			iarg->idx = -1;
			iarg->return_log = GMRFLib_TRUE;       /* return log(density), yes. */
			iarg->hyper_count = -1;
			iarg->hyper_z = NULL;
			iarg->hyper_ldens = NULL;
			iarg->theta_mode = &Theta(0);
			iarg->sqrt_eigen_values = a->sqrt_eigen_values;
			iarg->eigen_vectors = a->eigen_vectors;
			iarg->z = Calloc(nhyper, double);
			iarg->theta = Calloc(nhyper, double);
			iarg->stdev_corr_pos = a->stdev_corr_pos;
			iarg->stdev_corr_neg = a->stdev_corr_neg;
			iarg->dz = -1;
			iarg->interpolator = GMRFLib_AI_INTERPOLATOR_CCD;

			int npoints = 51;
			double *x = Calloc(nhyper, double), *xx = NULL, *xxx = Calloc(npoints, double), *ldens_values = Calloc(npoints, double);

			GMRFLib_ghq_abscissas(&xx, npoints);
			memcpy(xxx, xx, npoints * sizeof(double));
			xxx[0] = DMIN(xxx[0], -GMRFLib_DENSITY_INTEGRATION_LIMIT);
			xxx[npoints - 1] = DMAX(xxx[npoints - 1], GMRFLib_DENSITY_INTEGRATION_LIMIT);

			for (ii = 0; ii < npoints; ii++) {
				for (jj = 0; jj < nhyper; jj++) {
					x[jj] = Theta(jj) + Sigma_a[jj] * xxx[ii] / sqrt(var);
				}
				ldens_values[ii] = GMRFLib_ai_integrator_func(nhyper, x, iarg);
			}
			GMRFLib_density_create(&(GMRFLib_ai_INLA_userfunc3_density[number][i]),
					       GMRFLib_DENSITY_TYPE_SCGAUSSIAN, npoints, xxx, ldens_values, mean, sqrt(var), GMRFLib_TRUE);

			Free(Sigma_a);
			Free(x);
			Free(xxx);
			Free(ldens_values);
			Free(iarg->z);
			Free(iarg->theta);
			Free(iarg);
		}
	} else {
		for (i = 0; i < nrow; i++) {
			double mean = 0.0, var = 0.0;

			GMRFLib_matrix_get_row(row_spde3, i, model->BLC);
			mean = row_spde3[0];
			for (ii = 0; ii < model->ntheta; ii++) {
				mean += Theta_spde3(ii) * row_spde3[1 + ii];
				for (jj = 0; jj < model->ntheta; jj++) {
					var += row_spde3[1 + ii] * row_spde3[1 + jj] * Cov_spde3(ii, jj);	/* yes the first column is a constant offset */
				}
			}
			GMRFLib_density_create_normal(&(GMRFLib_ai_INLA_userfunc3_density[number][i]), 0.0, 1.0, mean, (var > 0 ? sqrt(var) : DBL_EPSILON));
		}
	}

	Free(row_spde3);
	Free(row);

#undef Cov_spde3
#undef Cov
#undef Theta_spde3
#undef Theta
	return INLA_OK;
}
