
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
 *        Haavard Rue
 *        CEMSE Division
 *        King Abdullah University of Science and Technology
 *        Thuwal 23955-6900, Saudi Arabia
 *        Email: haavard.rue@kaust.edu.sa
 *        Office: +966 (0)12 808 0640
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

int inla_spde3_build_model(inla_spde3_tp ** smodel, const char *prefix, const char *transform)
{
	int i, debug = 0;
	inla_spde3_tp *model = NULL;
	char *fnm = NULL;

	model = Calloc(1, inla_spde3_tp);

	if (!strcasecmp(transform, "identity")) {
		model->transform = SPDE3_TRANSFORM_IDENTITY;
	} else if (!strcasecmp(transform, "log")) {
		model->transform = SPDE3_TRANSFORM_LOG;
	} else if (!strcasecmp(transform, "shiftedlog")) {
		model->transform = SPDE3_TRANSFORM_SHIFTEDLOG;
	} else if (!strcasecmp(transform, "logit")) {
		model->transform = SPDE3_TRANSFORM_LOGIT;
	} else if (!strcasecmp(transform, "oldlogit")) {
		model->transform = SPDE3_TRANSFORM_OLDLOGIT;
	} else {
		assert(0 == 1);
	}

	model->B = Calloc(4, GMRFLib_matrix_tp *);
	model->M = Calloc(4, GMRFLib_matrix_tp *);
	for (i = 0; i < 4; i++) {
		GMRFLib_sprintf(&fnm, "%s%s%1d", prefix, "B", i);
		model->B[i] = GMRFLib_read_fmesher_file((const char *) fnm, 0, -1);

		GMRFLib_sprintf(&fnm, "%s%s%1d", prefix, "M", i);
		model->M[i] = GMRFLib_read_fmesher_file((const char *) fnm, 0, -1);
	}

	for (i = 1; i < 4; i++) {
		/*
		 * all need the same dimensions n x (p+1) , except for B[3]
		 */
		if (i != 3) {
			assert(model->B[0]->nrow == model->B[i]->nrow);
		}
		assert(model->B[0]->ncol == model->B[i]->ncol);

		/*
		 * all are square with the same dimension n x n, except M[3]
		 */
		if (i != 3) {
			assert(model->M[i]->nrow == model->M[i]->ncol);
		}
		assert(model->M[0]->ncol == model->M[i]->ncol);

		/*
		 * and the number of rows of B must be the same as the number of cols of M
		 */
		assert(model->B[i]->nrow == model->M[i]->nrow);
	}

	model->n = model->M[0]->nrow;
	model->n3 = model->M[3]->nrow;
	model->ntheta = model->B[0]->ncol - 1;

	/*
	 * big savings possible if M[3] is a matrix with only zeros 
	 */
	if (model->M[3]->elems == 0) {
		GMRFLib_matrix_free(model->M[3]);
		model->M[3] = NULL;
	}
	// since this matrix is non-symmetric, we need the graph for for the transpose as well. the
	// graph gives the entries for for each row.
	if (model->M[3]) {
		model->M3transpose = GMRFLib_matrix_transpose(model->M[3]);
	} else {
		model->M3transpose = NULL;
	}

	GMRFLib_sprintf(&fnm, "%s%s", prefix, "BLC");
	if (GMRFLib_is_fmesher_file((const char *) fnm, 0L, -1) == GMRFLib_SUCCESS) {
		model->BLC = GMRFLib_read_fmesher_file((const char *) fnm, 0, -1);
	} else {
		model->BLC = NULL;
	}

	if (debug) {
		P(model->n);
		P(model->n3);
		P(model->ntheta);
	}

	/*
	 * I need to build the graph. Need to add both M_ij and M_ji as M3 can be non-symmetric. 
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

	if (model->M[3]) {
		/*
		 * add the graph-contributions for the \sum_k d_k M_ki M_kj = \sum_k d_k M^T_ik M_kj - term for M^(3). M^(3) is non-symmetric.  the graph is
		 * row-wise, so g->nnbs[i] is the number of neighours to row i. so M_ki is M^T_ik.
		 */
		int j, jj, k, kk;
		GMRFLib_graph_tp *g = model->M3transpose->graph;
		GMRFLib_graph_tp *gg = model->M[3]->graph;

		for (i = 0; i < model->n3; i++) {
			// case k = i. then we get M_ii*M_ij, so we need to add all M3-neighbours to i.
			for (kk = 0; kk < gg->nnbs[i]; kk++) {
				j = gg->nbs[i][kk];
				GMRFLib_ged_add(ged, i, j);
			}
			// case i ~ k
			for (kk = 0; kk < g->nnbs[i]; kk++) {
				k = g->nbs[i][kk];
				// we have now M_ik, and k != i. then we need to add j for which j = k ...
				j = k;
				GMRFLib_ged_add(ged, i, j);
				// ... and then all M3-neighbours to k
				for (jj = 0; jj < gg->nnbs[k]; jj++) {
					j = gg->nbs[k][jj];
					GMRFLib_ged_add(ged, i, j);
				}
			}
		}
	}

	GMRFLib_ged_build(&(model->graph), ged);
	assert(model->n == model->graph->n);
	GMRFLib_ged_free(ged);

	model->Qfunc = inla_spde3_Qfunction;
	model->Qfunc_arg = (void *) model;

	model->store = Calloc(ISQR(GMRFLib_MAX_THREADS), inla_spde3_d3store_tp *);
	for (i = 0; i < ISQR(GMRFLib_MAX_THREADS); i++) {
		int j;

		model->store[i] = Calloc(1, inla_spde3_d3store_tp);
		model->store[i]->theta = Calloc(model->ntheta, double);
		for (j = 0; j < model->ntheta; j++) {
			model->store[i]->theta[j] = GMRFLib_uniform();
		}
		model->store[i]->d3 = NULL;		       /* alloc this one when it is used */
	}

	HYPER_NEW2(model->theta, 0.0, model->ntheta);
	*smodel = model;

	return INLA_OK;
}
double inla_spde3_Qfunction(int i, int j, void *arg)
{
	inla_spde3_tp *model = (inla_spde3_tp *) arg;
	double value, phi_i[3], phi_j[3], d_i[3], d_j[3];
	int k, kk, use_store = 1, debug = 0;

	/*
	 * to hold the i'th and j'th and k'th row of the B-matrices. use one storage only
	 */
	double *row_i = Calloc(3 * model->B[0]->ncol, double), *row_j, *row_k;
	row_j = &row_i[model->B[0]->ncol];
	row_k = &row_i[2 * model->B[0]->ncol];

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

	for (k = 0; k < 2; k++) {
		d_i[k] = exp(phi_i[k]);
		d_j[k] = exp(phi_j[k]);
	}

	/*
	 * change this later on, need an option here for various 'link' functions. some savings possible for i==j.
	 */
	if (i == j) {
		switch (model->transform) {
		case SPDE3_TRANSFORM_IDENTITY:
			d_i[2] = phi_i[2];
			break;
		case SPDE3_TRANSFORM_LOG:
			d_i[2] = exp(phi_i[2]);
			break;
		case SPDE3_TRANSFORM_SHIFTEDLOG:
			d_i[2] = 2.0 * exp(phi_i[2]) - 1.0;
			break;
		case SPDE3_TRANSFORM_LOGIT:
			d_i[2] = (1.0 - exp(phi_i[2])) / (1.0 + exp(phi_i[2]));
			break;
		case SPDE3_TRANSFORM_OLDLOGIT:
			d_i[2] = cos(M_PI / (1.0 + exp(-phi_i[2])));
			break;
		default:
			assert(0 == 1);
		}
		d_j[2] = d_i[2];
	} else {
		switch (model->transform) {
		case SPDE3_TRANSFORM_IDENTITY:
			d_i[2] = phi_i[2];
			d_j[2] = phi_j[2];
			break;
		case SPDE3_TRANSFORM_LOG:
			d_i[2] = exp(phi_i[2]);
			d_j[2] = exp(phi_j[2]);
			break;
		case SPDE3_TRANSFORM_SHIFTEDLOG:
			d_i[2] = 2.0 * exp(phi_i[2]) - 1.0;
			d_j[2] = 2.0 * exp(phi_j[2]) - 1.0;
			break;
		case SPDE3_TRANSFORM_LOGIT:
			d_i[2] = (1.0 - exp(phi_i[2])) / (1.0 + exp(phi_i[2]));
			d_j[2] = (1.0 - exp(phi_j[2])) / (1.0 + exp(phi_j[2]));
			break;
		case SPDE3_TRANSFORM_OLDLOGIT:
			d_i[2] = cos(M_PI / (1.0 + exp(-phi_i[2])));
			d_j[2] = cos(M_PI / (1.0 + exp(-phi_j[2])));
			break;
		default:
			assert(0 == 1);
		}
	}

	value = d_i[0] * d_j[0] * (d_i[1] * d_j[1] * GMRFLib_matrix_get(i, j, model->M[0]) +
				   d_i[2] * d_i[1] * GMRFLib_matrix_get(i, j, model->M[1]) +
				   d_j[1] * d_j[2] * GMRFLib_matrix_get(j, i, model->M[1]) + GMRFLib_matrix_get(i, j, model->M[2]));

	// Add the new M^(3) term
#define COMPUTE_D3(_d3, _k) if (1)					\
	{								\
		double _tmp = NAN;					\
		int _kk;						\
		GMRFLib_matrix_get_row(row_k, _k, model->B[3]);		\
		_tmp = row_k[0];					\
		for (_kk = 1; _kk < model->B[3]->ncol; _kk++) {		\
			_tmp += row_k[_kk] * model->theta[_kk - 1][GMRFLib_thread_id][0]; \
		}							\
		_d3 = exp(_tmp);					\
	}

	if (model->M[3]) {
		int id = GMRFLib_thread_id + omp_get_thread_num() * GMRFLib_MAX_THREADS;

		if (use_store) {
			// check if we need to recompute storage
			int recompute = 0;
			for (k = 0; k < model->ntheta; k++) {
				if (model->theta[k][GMRFLib_thread_id][0] != model->store[id]->theta[k]) {
					recompute = 1;
					break;
				}
			}
			if (recompute) {
				if (debug) {
					fprintf(stderr, "recompute id=%1d\n", id);
				}
				for (k = 0; k < model->ntheta; k++) {
					model->store[id]->theta[k] = model->theta[k][GMRFLib_thread_id][0];
				}
				if (!(model->store[id]->d3)) {
					model->store[id]->d3 = Calloc(model->n3, double);
				}
				for (k = 0; k < model->n3; k++) {
					COMPUTE_D3(model->store[id]->d3[k], k);
				}
			}
		}
		// For the M3-term, it is easier to work with M3transpose. We need to be careful since M3 is generally not
		// symmetric.
		double m3_value = 0.0, d3 = NAN, tmp;
		GMRFLib_graph_tp *g = model->M3transpose->graph;
		double *m3t_row_i = Calloc(2 * g->n, double);
		double *m3t_row_j = m3t_row_i + g->n;
		int i_min;

		GMRFLib_matrix_get_row(m3t_row_i, i, model->M3transpose);
		GMRFLib_matrix_get_row(m3t_row_j, j, model->M3transpose);
		i_min = (g->nnbs[i] <= g->nnbs[j] ? i : j);    /* loop over the one with fewest neigbours */

		if (use_store) {
			for (kk = 0; kk <= g->nnbs[i_min]; kk++) {	/* yes! we want '<=' */
				k = (kk < g->nnbs[i_min] ? g->nbs[i_min][kk] : i_min);	/* since we include the case k=i here */
				if (k < model->n3) {	       /* it might be that M3_ii is not there */
					m3_value += m3t_row_i[k] * m3t_row_j[k] * model->store[id]->d3[k];
				}
			}
		} else {
			for (kk = 0; kk <= g->nnbs[i_min]; kk++) {	/* yes! we want '<=' */
				k = (kk < g->nnbs[i_min] ? g->nbs[i_min][kk] : i_min);	/* since we include the case k=i here */
				if (k < model->n3) {	       /* it might be that M3_ii is not there */
					tmp = m3t_row_i[k] * m3t_row_j[k];
					if (tmp) {	       /* this term is often zero, and computing D3 is expensive */
						COMPUTE_D3(d3, k);
						m3_value += tmp * d3;
					}
				}
			}
		}
		Free(m3t_row_i);
		value += d_i[0] * d_j[0] * m3_value;
	}
#undef D3
	Free(row_i);

	return value;
}
double *inla_spde3_userfunc3(int number, double *theta, int nhyper, double *covmat, void *arg)
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
		return NULL;
	}

	int use_new_version = 1;
	int i, ii, jj;
	inla_spde3_tp *model = (inla_spde3_tp *) GMRFLib_ai_INLA_userfunc3_args[number];

	if (model->BLC == NULL) {
		return NULL;
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
					var += row_spde3[1 + ii] * row_spde3[1 + jj] * Cov_spde3(ii, jj);	/* yes the first column is a
														 * constant offset */
				}
			}
			GMRFLib_density_create_normal(&(GMRFLib_ai_INLA_userfunc3_density[number][i]), 0.0, 1.0, mean,
						      (var > 0 ? sqrt(var) : DBL_EPSILON));
		}
	}

	Free(row_spde3);
	Free(row);

#undef Cov_spde3
#undef Cov
#undef Theta_spde3
#undef Theta
	return NULL;
}
