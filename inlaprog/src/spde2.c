
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

	for (k = 0; k < 3; k++) {
		GMRFLib_matrix_tp *B;

		switch (k) {
		case 0:
			B = model->B0;
			break;
		case 1:
			B = model->B1;
			break;
		case 2:
			B = model->B2;
			break;
		default:
			assert(0 == 1);
		}

		if (i == j) {
			/*
			 * some savings for i == j 
			 */
			phi_i[k] = GMRFLib_matrix_get(i, 0, B);
			for (kk = 1; kk < B->ncol; kk++) {
				/*
				 * '-1' is the correction for the first intercept column in B 
				 */
				phi_i[k] += GMRFLib_matrix_get(i, kk, B) * model->theta[kk - 1][GMRFLib_thread_id][0];
			}
			phi_j[k] = phi_i[k];		       /* they are equal in this case */
		} else {
			/*
			 * i != j 
			 */
			phi_i[k] = GMRFLib_matrix_get(i, 0, B);
			phi_j[k] = GMRFLib_matrix_get(j, 0, B);
			for (kk = 1; kk < B->ncol; kk++) {
				/*
				 * '-1' is the correction for the first intercept column in B 
				 */
				phi_i[k] += GMRFLib_matrix_get(i, kk, B) * model->theta[kk - 1][GMRFLib_thread_id][0];
				phi_j[k] += GMRFLib_matrix_get(j, kk, B) * model->theta[kk - 1][GMRFLib_thread_id][0];
			}
		}
	}

	for (k = 0; k < 2; k++) {
		d_i[k] = exp(phi_i[k]);
		d_j[k] = exp(phi_j[k]);
	}

	/*
	 * change this later on, need an option here for various 'link' functions. 
	 */
	d_i[2] = cos(M_PI * map_probability(phi_i[2], MAP_FORWARD, NULL));
	d_j[2] = cos(M_PI * map_probability(phi_j[2], MAP_FORWARD, NULL));

	value = d_i[0] * d_j[0] * (d_i[1] * d_j[1] * GMRFLib_matrix_get(i, j, model->M0) +
				   d_i[2] * d_j[1] * GMRFLib_matrix_get(i, j, model->M1) +
				   d_i[1] * d_j[2] * GMRFLib_matrix_get(j, i, model->M1) + GMRFLib_matrix_get(i, j, model->M2));

	return value;
}

int inla_spde2_build_model(inla_spde2_tp ** smodel, const char *prefix)
{
	int i, j, debug = 1;
	inla_spde2_tp *model = NULL;
	char *fnm = NULL;

	model = Calloc(1, inla_spde2_tp);

	GMRFLib_sprintf(&fnm, "%s%s", prefix, "B0");
	model->B0 = GMRFLib_read_fmesher_file((const char *) fnm, 0, -1);

	GMRFLib_sprintf(&fnm, "%s%s", prefix, "B1");
	model->B1 = GMRFLib_read_fmesher_file((const char *) fnm, 0, -1);

	GMRFLib_sprintf(&fnm, "%s%s", prefix, "B2");
	model->B2 = GMRFLib_read_fmesher_file((const char *) fnm, 0, -1);

	GMRFLib_sprintf(&fnm, "%s%s", prefix, "M0");
	model->M0 = GMRFLib_read_fmesher_file((const char *) fnm, 0, -1);

	GMRFLib_sprintf(&fnm, "%s%s", prefix, "M1");
	model->M1 = GMRFLib_read_fmesher_file((const char *) fnm, 0, -1);

	GMRFLib_sprintf(&fnm, "%s%s", prefix, "M2");
	model->M2 = GMRFLib_read_fmesher_file((const char *) fnm, 0, -1);

	assert(model->B0->nrow == model->B1->nrow);
	assert(model->B0->nrow == model->B2->nrow);
	assert(model->B0->ncol == model->B1->ncol);
	assert(model->B0->ncol == model->B2->ncol);

	assert(model->M0->nrow == model->M0->ncol);
	assert(model->M1->nrow == model->M1->ncol);
	assert(model->M2->nrow == model->M2->ncol);
	assert(model->M0->nrow == model->M1->nrow);
	assert(model->M0->nrow == model->M2->nrow);

	model->n = model->M0->nrow;
	model->ntheta = model->B0->ncol - 1;

	GMRFLib_sprintf(&fnm, "%s%s", prefix, "BLC");
	model->BLC = GMRFLib_read_fmesher_file((const char *) fnm, 0, -1);

	if (debug) {
		P(model->n);
		P(model->ntheta);
	}

	/*
	 * I need to build the graph. Need to add both M_ij and M_ji as M1 can be non-symmetric. 
	 */
	GMRFLib_ged_tp *ged = NULL;
	GMRFLib_ged_init(&ged, NULL);

#define ADD_GRAPH(_G)					\
	for(i = 0; i< _G->n; i++) {			\
		GMRFLib_ged_add(ged, i, i);		\
		int jj;					\
		for(jj = 0; jj < _G->nnbs[i]; jj++){	\
			j = _G->nbs[i][jj];		\
			GMRFLib_ged_add(ged, i, j);	\
			GMRFLib_ged_add(ged, j, i);	\
		}					\
	}

	ADD_GRAPH(model->M0->graph);
	ADD_GRAPH(model->M1->graph);
	ADD_GRAPH(model->M2->graph);

#undef ADD_GRAPH

	GMRFLib_ged_build(&(model->graph), ged);
	assert(model->n == model->graph->n);
	GMRFLib_ged_free(ged);

	model->Qfunc = inla_spde2_Qfunction;
	model->Qfunc_arg = (void *) model;

	HYPER_NEW2(model->theta, 0.0, model->ntheta);
	// model->theta_extra = Calloc(GMRFLib_MAX_THREADS, double *);

	if (0) {
		FILE *fp;

		FIXME("write graph to file spde2-graph.dat");
		fp = fopen("spde2-graph.dat", "w");
		GMRFLib_print_graph(fp, model->graph);
		fclose(fp);
	}

	if (0) {
		int k;
		for (k = 1; k < 8; k++) {
			GMRFLib_problem_tp *problem;
			GMRFLib_reorder = k;
			// GMRFLib_optimize_reorder(model->graph, NULL);
			GMRFLib_init_problem(&problem, NULL, NULL, NULL, NULL, model->graph, model->Qfunc, model->Qfunc_arg, NULL, NULL, GMRFLib_NEW_PROBLEM);
			char *nm;
			GMRFLib_sprintf(&nm, "Qspde2-%1d", k);
			GMRFLib_bitmap_problem(nm, problem);
			GMRFLib_free_problem(problem);
		}
	}

	*smodel = model;

	return INLA_OK;
}
