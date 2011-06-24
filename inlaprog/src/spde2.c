
/* spde2.c
 * 
 * Copyright (C) 2011  Havard Rue
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
#include "spde2.h"


extern G_tp G;						       /* import some global parametes from inla */

double inla_spde2_Qfunction(int i, int j, void *arg)
{
	inla_spde2_tp *model = (inla_spde2_tp *) arg;
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
			for (kk = 1; kk < B->ncol; kk++) {
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
			for (kk = 1; kk < B->ncol; kk++) {
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
		case SPDE2_TRANSFORM_LOGIT:
			d_i[2] = cos(M_PI * map_probability(phi_i[2], MAP_FORWARD, NULL));
			break;
		case SPDE2_TRANSFORM_LOG:
			d_i[2] = 2 * exp(phi_i[2]) - 1.0;
			break;
		case SPDE2_TRANSFORM_IDENTITY:
			d_i[2] = phi_i[2];
			break;
		default:
			assert(0 == 1);
		}
		d_j[2] = d_i[2];
	} else {
		switch (model->transform) {
		case SPDE2_TRANSFORM_LOGIT:
			d_i[2] = cos(M_PI * map_probability(phi_i[2], MAP_FORWARD, NULL));
			d_j[2] = cos(M_PI * map_probability(phi_j[2], MAP_FORWARD, NULL));
			break;
		case SPDE2_TRANSFORM_LOG:
			d_i[2] = 2 * exp(phi_i[2]) - 1.0;
			d_j[2] = 2 * exp(phi_j[2]) - 1.0;
			break;
		case SPDE2_TRANSFORM_IDENTITY:
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

int inla_spde2_build_model(inla_spde2_tp ** smodel, const char *prefix, const char *transform)
{
	int i, debug = 0;
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
	// model->theta_extra = Calloc(GMRFLib_MAX_THREADS, double *);

	*smodel = model;

	return INLA_OK;
}
