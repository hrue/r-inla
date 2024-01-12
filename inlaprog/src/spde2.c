
/* spde2.c
 * 
 * Copyright (C) 2011-2024  Havard Rue
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

double inla_spde2_Qfunction(int thread_id, int ii, int jj, double *values, void *arg)
{
	if (jj < 0) {
		inla_spde2_tp *model = (inla_spde2_tp *) arg;
		int nc = model->B[0]->ncol;

		int nb = model->graph->lnnbs[ii];
		double *V = model->row_V[ii];
		double *v = model->row_v[ii];
		double theta[nc];

		theta[0] = 1.0;
		for (int k = 1; k < nc; k++) {
			theta[k] = model->theta[k - 1][thread_id][0];
		}

		// dij = V %*% theta
		int m = nc;
		int lda = nc;
		int n = (1 + nb) * 3;
		int inc = 1;
		double alpha = 1.0;
		double beta = 0.0;
		double dij[(1 + nb) * 3];
		dgemv_("T", &m, &n, &alpha, V, &lda, theta, &inc, &beta, dij, &inc, F_ONE);

		GMRFLib_exp_inc(1 + nb, dij, 3, dij);
		GMRFLib_exp_inc(1 + nb, dij + 1, 3, dij + 1);

		if (model->transform != SPDE2_TRANSFORM_IDENTITY) {
			switch (model->transform) {
			case SPDE2_TRANSFORM_LOG:
			{
				for (int kk = 0, off = 2; kk < 1 + nb; kk++, off += 3) {
					dij[off] = 2.0 * exp(dij[off]) - 1.0;
				}
			}
				break;

			case SPDE2_TRANSFORM_LOGIT:
			{
				for (int kk = 0, off = 2; kk < 1 + nb; kk++, off += 3) {
					dij[off] = cos(M_PI / (1.0 + exp(-dij[off])));
				}
			}
				break;

			case SPDE2_TRANSFORM_IDENTITY:
				break;

			default:
				assert(0 == 1);
			}
		}

		double *d_i = dij;
		double *vv = v;
		double *d_j = dij;
		for (int kk = 0; kk < 1 + nb; kk++, vv += 4, d_j += 3) {
			values[kk] = d_i[0] * d_j[0] * (d_i[1] * d_j[1] * vv[0] + d_i[2] * d_i[1] * vv[1] + d_j[1] * d_j[2] * vv[2] + vv[3]);
		}
		return 0.0;
	} else {
		return inla_spde2_Qfunction_ij(thread_id, IMIN(ii, jj), IMAX(ii, jj), values, arg);
	}
}

double inla_spde2_Qfunction_ij(int thread_id, int ii, int jj, double *UNUSED(values), void *arg)
{
	// do not use directly. need ``if (jj < 0)'' code

	int use_ddot_lim = 16;
	inla_spde2_tp *model = (inla_spde2_tp *) arg;
	int nc = model->B[0]->ncol, nc2 = 2 * nc;
	double d_i[6] = { 0, 0, 0, 0, 0, 0 };
	double *d_j = d_i + 3;
	double *vals_i = model->row_V[ii];

	double theta[nc];
	theta[0] = 1.0;
	for (int k = 1; k < nc; k++) {
		theta[k] = model->theta[k - 1][thread_id][0];
	}

	if (nc < use_ddot_lim) {
		// better to inline the code, as the vectors are pretty short
		if (0) {
			d_i[0] = exp(GMRFLib_ddot(nc, vals_i, theta));
			d_i[1] = exp(GMRFLib_ddot(nc, vals_i + nc, theta));
			d_i[2] = GMRFLib_ddot(nc, vals_i + 2 * nc, theta);
		} else {
			double d0 = 0.0, d1 = 0.0, d2 = 0.0;
#pragma omp simd reduction(+: d0, d1, d2)
			for (int k = 0; k < nc; k++) {
				d0 += vals_i[k] * theta[k];
				d1 += vals_i[k + nc] * theta[k];
				d2 += vals_i[k + nc2] * theta[k];
			}
			d_i[0] = exp(d0);
			d_i[1] = exp(d1);
			d_i[2] = d2;
		}
	} else {
		int m = nc;
		int lda = nc;
		int n = 3;
		int inc = 1;
		double alpha = 1.0;
		double beta = 0.0;
		dgemv_("T", &m, &n, &alpha, vals_i, &lda, theta, &inc, &beta, d_i, &inc, F_ONE);
		d_i[0] = exp(d_i[0]);
		d_i[1] = exp(d_i[1]);
	}

	if (model->transform != SPDE2_TRANSFORM_IDENTITY) {
		switch (model->transform) {
		case SPDE2_TRANSFORM_LOG:
			d_i[2] = 2 * exp(d_i[2]) - 1.0;
			break;

		case SPDE2_TRANSFORM_LOGIT:
			d_i[2] = cos(M_PI / (1.0 + exp(-d_i[2])));
			break;

		case SPDE2_TRANSFORM_IDENTITY:
			break;

		default:
			assert(0 == 1);
		}
	}

	if (ii == jj) {
		double *v = model->row_v[ii];
		double value = SQR(d_i[0]) * (SQR(d_i[1]) * v[0] + d_i[2] * d_i[1] * (v[1] + v[2]) + v[3]);
		return value;
	}

	spde2_vV_tp *vals_j_p = (spde2_vV_tp *) * map_ivp_ptr(&(model->Vmatrix->vmat[ii]), jj);
	double *vals_j = vals_j_p->V;

	if (nc < use_ddot_lim) {
		// better to inline the code, as the vectors are pretty short
		if (0) {
			d_j[0] = exp(GMRFLib_ddot(nc, vals_j, theta));
			d_j[1] = exp(GMRFLib_ddot(nc, vals_j + nc, theta));
			d_j[2] = GMRFLib_ddot(nc, vals_j + 2 * nc, theta);
		} else {
			double d0 = 0.0, d1 = 0.0, d2 = 0.0;
#pragma omp simd reduction(+: d0, d1, d2)
			for (int k = 0; k < nc; k++) {
				d0 += vals_j[k] * theta[k];
				d1 += vals_j[k + nc] * theta[k];
				d2 += vals_j[k + nc2] * theta[k];
			}
			d_j[0] = exp(d0);
			d_j[1] = exp(d1);
			d_j[2] = d2;
		}
	} else {
		int m = nc;
		int lda = nc;
		int n = 3;
		int inc = 1;
		double alpha = 1.0;
		double beta = 0.0;
		dgemv_("T", &m, &n, &alpha, vals_j, &lda, theta, &inc, &beta, d_j, &inc, F_ONE);
		d_j[0] = exp(d_j[0]);
		d_j[1] = exp(d_j[1]);
	}

	if (model->transform != SPDE2_TRANSFORM_IDENTITY) {
		switch (model->transform) {
		case SPDE2_TRANSFORM_LOG:
			d_j[2] = 2.0 * exp(d_j[2]) - 1.0;
			break;

		case SPDE2_TRANSFORM_LOGIT:
			d_j[2] = cos(M_PI / (1.0 + exp(-d_j[2])));
			break;

		case SPDE2_TRANSFORM_IDENTITY:
			break;

		default:
			assert(0 == 1);
		}
	}

	double *v = vals_j_p->v;
	double value = d_i[0] * d_j[0] * (d_i[1] * d_j[1] * v[0] + d_i[2] * d_i[1] * v[1] + d_j[1] * d_j[2] * v[2] + v[3]);

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

	model->row_V = Calloc(model->n, double *);
	model->row_v = Calloc(model->n, double *);
	GMRFLib_vmatrix_init(&(model->Vmatrix), model->n, model->graph);
	int nc = model->B[0]->ncol;

	for (i = 0; i < model->n; i++) {
		int j = i;

		int nb = model->graph->lnnbs[i];
		int msize = 3 * nc;
		int vsize = 4;
		double *V = Calloc((1 + nb) * msize, double);
		double *v = Calloc((1 + nb) * vsize, double);

		model->row_V[i] = V;
		model->row_v[i] = v;

		GMRFLib_matrix_get_row(V + 0 * nc, i, model->B[0]);
		GMRFLib_matrix_get_row(V + 1 * nc, i, model->B[1]);
		GMRFLib_matrix_get_row(V + 2 * nc, i, model->B[2]);

		v[0] = GMRFLib_matrix_get(i, j, model->M[0]);
		v[1] = GMRFLib_matrix_get(i, j, model->M[1]);
		v[2] = GMRFLib_matrix_get(j, i, model->M[1]);
		v[3] = GMRFLib_matrix_get(i, j, model->M[2]);

		spde2_vV_tp *vV = Calloc(1, spde2_vV_tp);
		vV->v = v;
		vV->V = V;
		map_ivp_set(&(model->Vmatrix->vmat[i]), i, (void *) vV);

		for (int jj = 0; jj < model->graph->lnnbs[i]; jj++) {
			j = model->graph->lnbs[i][jj];

			V += msize;
			GMRFLib_matrix_get_row(V + 0 * nc, j, model->B[0]);
			GMRFLib_matrix_get_row(V + 1 * nc, j, model->B[1]);
			GMRFLib_matrix_get_row(V + 2 * nc, j, model->B[2]);

			v += vsize;
			v[0] = GMRFLib_matrix_get(i, j, model->M[0]);
			v[1] = GMRFLib_matrix_get(i, j, model->M[1]);
			v[2] = GMRFLib_matrix_get(j, i, model->M[1]);
			v[3] = GMRFLib_matrix_get(i, j, model->M[2]);

			vV = Calloc(1, spde2_vV_tp);
			vV->v = v;
			vV->V = V;
			map_ivp_set(&(model->Vmatrix->vmat[i]), j, (void *) vV);
		}
	}

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
