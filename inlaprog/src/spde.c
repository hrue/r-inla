
/* spde.c
 * 
 * Copyright (C) 2007-2010 Havard Rue
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

/* Pre-hg-Id: $Id: spde.c,v 1.41 2010/03/01 17:43:07 hrue Exp $ */

#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"
#include "inla.h"
#include "spde.h"


extern G_tp G;						       /* import some global parametes from inla */

/* 
   This is to get the E() and Var() of the deformations itself.  This is a bit dirty, but for now on... I need to find a better way to do this.
 */
static inla_spde_tp *func_smodel = NULL;

inla_spde_points_tp *inla_spde_set_points(GMRFLib_matrix_tp * M)
{
	assert(M->nrow > 0);
	assert(M->ncol > 0);

	int i, j;
	double *hold = NULL;
	inla_spde_points_tp *p = NULL;

	hold = Calloc(M->nrow * M->ncol, double);
	p = Calloc(1, inla_spde_points_tp);
	p->n = M->nrow;
	p->dim = M->ncol;
	p->s = Calloc(M->nrow, double *);
	for (i = 0; i < M->nrow; i++) {
		p->s[i] = &hold[i * M->ncol];
	}

	for (i = 0; i < M->nrow; i++) {
		for (j = 0; j < M->ncol; j++) {
			p->s[i][j] = GMRFLib_matrix_get(i, j, M);
		}
	}

	return p;
}
int inla_spde_free_points(inla_spde_points_tp * p)
{
	if (p) {
		Free(p->s);
		Free(p);
	}
	return INLA_OK;
}
double inla_spde_Qfunction(int node, int nnode, void *arg)
{
	inla_spde_tp *model = (inla_spde_tp *) arg;
	double value;

	/*
	 * Store computed values of the ocillating coeff... 
	 */
	typedef struct {
		double theta;
		double oc;
	} OC_tp;
	static OC_tp *OC = NULL;
	OC_tp *a;

	if (!OC) {
#pragma omp critical
		{
			if (!OC) {
				OC = Calloc(GMRFLib_MAX_THREADS, OC_tp);
			}
		}
	}

	a = &(OC[GMRFLib_thread_id]);
	// argument is cos(theta*PI)
	if (a->theta != model->oc[GMRFLib_thread_id][0]) {
		a->theta = model->oc[GMRFLib_thread_id][0];
		a->oc = cos(M_PI * map_probability(a->theta, MAP_FORWARD, NULL));
	}

	if (node == nnode) {
		double G1ii, G2ii, K, T;

		G1ii = model->G1->Qfunc(node, node, (void *) model->G1->Qfunc_arg);
		G2ii = model->G2->Qfunc(node, node, (void *) model->G2->Qfunc_arg);

		if (model->K) {
			K = model->K[node];
		} else {
			K = inla_spde_KT_model_eval(model->Kmodel, node);
		}

		if (model->T) {
			T = model->T[node];
		} else {
			T = inla_spde_KT_model_eval(model->Tmodel, node);
		}

		value = SQR(T) * (SQR(K) * model->C[node] + a->oc * 2.0 * K * G1ii + G2ii);
	} else {
		double G1ij, G2ij, K, KK, T, TT;

		if (GMRFLib_is_neighb(node, nnode, model->G1_graph)) {
			G1ij = model->G1->Qfunc(node, nnode, (void *) model->G1->Qfunc_arg);
		} else {
			G1ij = 0.0;
		}
		G2ij = model->G2->Qfunc(node, nnode, (void *) model->G2->Qfunc_arg);

		if (model->K) {
			K = model->K[node];
			KK = model->K[nnode];
		} else {
			inla_spde_KT_model_eval2(&K, &KK, model->Kmodel, node, nnode);
		}

		if (model->T) {
			T = model->T[node];
			TT = model->T[nnode];
		} else {
			inla_spde_KT_model_eval2(&T, &TT, model->Tmodel, node, nnode);
		}

		value = T * (a->oc * (K * G1ij + KK * G1ij) + G2ij) * TT;
	}
	return value;
}
int inla_spde_KT_model_init(inla_spde_theta_tp * theta_model, GMRFLib_matrix_tp * basis)
{
	double ***theta, *hold;
	int i, j;

	theta_model->ntheta = basis->ncol;
	assert(theta_model->ntheta > 0);

	HYPER_NEW2(theta, 0.0, theta_model->ntheta);
	theta_model->theta = theta;
	theta_model->theta_extra = Calloc(GMRFLib_MAX_THREADS, double *);

	theta_model->n = basis->nrow;
	theta_model->basis = Calloc(basis->nrow, double *);
	hold = Calloc(basis->elems, double);

	for (i = 0; i < basis->nrow; i++) {
		theta_model->basis[i] = &(hold[basis->ncol * i]);
		for (j = 0; j < basis->ncol; j++) {
			theta_model->basis[i][j] = GMRFLib_matrix_get(i, j, basis);
		}
	}
	return INLA_OK;
}
double inla_spde_KT_model_eval(inla_spde_theta_tp * theta_model, int idx)
{
	int i;
	double value;

	if (theta_model->theta_extra && theta_model->theta_extra[GMRFLib_thread_id]) {
		/*
		 * use the one from extra 
		 */
		// printf("From extra %g\n", theta_model->theta_extra[GMRFLib_thread_id][0]);
		for (i = 0, value = 0.0; i < theta_model->ntheta; i++) {
			value += theta_model->theta_extra[GMRFLib_thread_id][i] * theta_model->basis[idx][i];
		}
	} else {
		/*
		 * use the common one 
		 */
		// printf("From general %g\n", theta_model->theta[0][GMRFLib_thread_id][0]);
		for (i = 0, value = 0.0; i < theta_model->ntheta; i++) {
			value += theta_model->theta[i][GMRFLib_thread_id][0] * theta_model->basis[idx][i];
		}
	}
	return exp(value);
}
int inla_spde_KT_model_eval2(double *value0, double *value1, inla_spde_theta_tp * theta_model, int idx, int iidx)
{
	/*
	 * for speedup; do this for two locations at the same time 
	 */

	int i;
	double t, *b0, *b1;

	*value0 = 0.0;
	*value1 = 0.0;
	b0 = theta_model->basis[idx];
	b1 = theta_model->basis[iidx];

	if (theta_model->theta_extra && theta_model->theta_extra[GMRFLib_thread_id]) {
		/*
		 * use the theta defined from extra() 
		 */
		for (i = 0; i < theta_model->ntheta; i++) {
			t = theta_model->theta_extra[GMRFLib_thread_id][i];
			*value0 += t * b0[i];
			*value1 += t * b1[i];
		}
	} else {
		/*
		 * use the common theta 
		 */
		for (i = 0; i < theta_model->ntheta; i++) {
			t = theta_model->theta[i][GMRFLib_thread_id][0];
			*value0 += t * b0[i];
			*value1 += t * b1[i];
		}
	}
	*value0 = exp(*value0);
	*value1 = exp(*value1);

	return INLA_OK;
}
int inla_spde_build_model(inla_spde_tp ** smodel, const char *prefix)
{
	int n, i, j;
	inla_spde_tp *model;
	char *fnm;
	GMRFLib_matrix_tp *M;

	model = Calloc(1, inla_spde_tp);

	/*
	 * ocillating coeff 
	 */
	model->oc = Calloc(GMRFLib_MAX_THREADS, double *);
	for (i = 0; i < GMRFLib_MAX_THREADS; i++)
		model->oc[i] = Calloc(1, double);

	GMRFLib_sprintf(&fnm, "%s%s", prefix, "c0");
	M = GMRFLib_read_fmesher_file((const char *) fnm, 0, -1);
	n = M->nrow;
	model->C = GMRFLib_matrix_get_diagonal(M);
	GMRFLib_matrix_free(M);
	Free(fnm);

	if (0) {
		for (i = 0; i < n; i++)
			printf("%d %g\n", i, model->C[i]);
		exit(0);
	}

	GMRFLib_sprintf(&fnm, "%s%s", prefix, "s");
	M = GMRFLib_read_fmesher_file((const char *) fnm, 0, -1);
	model->s = inla_spde_set_points(M);
	assert(model->s->n == n);
	GMRFLib_matrix_free(M);
	Free(fnm);
	if (0) {
		for (i = 0; i < n; i++) {
			printf("%d ", i);
			for (j = 0; j < model->s->dim; j++)
				printf(" %g", model->s->s[i][j]);
			printf("\n");
		}
	}

	GMRFLib_sprintf(&fnm, "%s%s", prefix, "g1");
	M = GMRFLib_read_fmesher_file((const char *) fnm, 0, -1);
	assert(M->nrow == n);
	if (0) {
		int imax = 0;
		int jmax = 0;
		for (i = 0; i < M->elems; i++) {
			printf("%d %d %d %g\n", i, M->i[i], M->j[i], M->values[i]);
			imax = IMAX(imax, M->i[i]);
			jmax = IMAX(jmax, M->j[i]);
		}
		P(M->elems);
		P(imax);
		P(jmax);
	}

	GMRFLib_tabulate_Qfunc_from_list(&(model->G1), &(model->G1_graph), M->elems, M->i, M->j, M->values, n, NULL, NULL, NULL);

	GMRFLib_matrix_free(M);
	Free(fnm);

	GMRFLib_sprintf(&fnm, "%s%s", prefix, "g2");
	M = GMRFLib_read_fmesher_file((const char *) fnm, 0, -1);
	assert(M->nrow == n);
	GMRFLib_tabulate_Qfunc_from_list(&(model->G2), &(model->G2_graph), M->elems, M->i, M->j, M->values, n, NULL, NULL, NULL);
	GMRFLib_matrix_free(M);
	Free(fnm);

	GMRFLib_sprintf(&fnm, "%s%s", prefix, "basisT");
	if (GMRFLib_file_exists(fnm, "rb") == GMRFLib_SUCCESS) {
		M = GMRFLib_read_fmesher_file((const char *) fnm, 0, -1);
	} else {
		M = GMRFLib_matrix_1(n);
	}
	model->Tmodel = Calloc(1, inla_spde_theta_tp);
	inla_spde_KT_model_init(model->Tmodel, M);
	GMRFLib_matrix_free(M);
	Free(fnm);

	GMRFLib_sprintf(&fnm, "%s%s", prefix, "basisK");
	if (GMRFLib_file_exists(fnm, "rb") == GMRFLib_SUCCESS) {
		M = GMRFLib_read_fmesher_file((const char *) fnm, 0, -1);
	} else {
		M = GMRFLib_matrix_1(n);
	}
	model->Kmodel = Calloc(1, inla_spde_theta_tp);
	inla_spde_KT_model_init(model->Kmodel, M);
	GMRFLib_matrix_free(M);
	Free(fnm);

	model->Qfunc = inla_spde_Qfunction;
	model->Qfunc_arg = (void *) model;
	model->graph = model->G2_graph;
	model->n = model->graph->n;
	*smodel = model;
	func_smodel = model;				       /* store it also here */

	if (0) {
		FILE *fp;

		FIXME("write graph to file spde-graph.dat");
		fp = fopen("spde-graph.dat", "w");
		GMRFLib_print_graph(fp, model->graph);
		fclose(fp);
	}

	if (0) {
		int k;
		for (k = 1; k < 8; k++) {
			GMRFLib_problem_tp *problem;
			GMRFLib_reorder = k;
			// GMRFLib_optimize_reorder(model->graph, NULL);
			GMRFLib_init_problem(&problem, NULL, NULL, NULL, NULL, model->graph, model->Qfunc, model->Qfunc_arg, NULL,
					     NULL, GMRFLib_NEW_PROBLEM);
			char *nm;
			GMRFLib_sprintf(&nm, "Qspde%1d", k);
			GMRFLib_bitmap_problem(nm, problem);
			GMRFLib_free_problem(problem);
		}
	}

	return INLA_OK;
}
double *inla_spde_userfunc0(GMRFLib_problem_tp * problem, double *theta, int nhyper)
{
	/*
	 * return the log(deformations). First the T's so the K's
	 */
	assert(func_smodel);
	assert(func_smodel->Tmodel->theta_extra[GMRFLib_thread_id] == NULL);
	assert(func_smodel->Kmodel->theta_extra[GMRFLib_thread_id] == NULL);

	int i, n = func_smodel->n, debug = 0;
	double *deformations;

	GMRFLib_ai_INLA_userfunc0_dim = 2 * n;		       /* set this variable, yes */
	deformations = Calloc(GMRFLib_ai_INLA_userfunc0_dim, double);

	for (i = 0; i < n; i++) {
		deformations[i] = log(inla_spde_KT_model_eval(func_smodel->Tmodel, i));	/* First the T's */
		deformations[i + n] = log(inla_spde_KT_model_eval(func_smodel->Kmodel, i));	/* Then the K's */

		if (debug) {
			printf("deformations[%1d]: T= %g K= %g\n", i, deformations[i], deformations[i + n]);
		}
	}

	return deformations;
}
double *inla_spde_userfunc1(double *theta, int nhyper, double *covmat)
{
	/*
	 * compute the marginals for the log(deformations). First the T's so the K's. using the mode and the covariance at the mode
	 */

#define DO_COMPUTE							\
	for (i = 0, mean = var = 0.0; i < t->ntheta; i++) {		\
		mean += t->theta[i][GMRFLib_thread_id][0] * t->basis[idx][i]; \
		for (j = 0; j < t->ntheta; j++) {			\
			var += t->basis[idx][i] * t->basis[idx][j] *	\
				covmat[ (i + where_to_start + offset_theta) + (j + where_to_start + offset_theta) * nhyper ]; \
		}							\
	}								\
	GMRFLib_density_create_normal(&(GMRFLib_ai_INLA_userfunc1_density[idx+offset]), 0.0, 1.0, mean, (var > 0 ? sqrt(var) : DBL_EPSILON))


	if (!covmat) {
		return NULL;
	}

	assert(func_smodel);
	assert(func_smodel->Tmodel->theta_extra[GMRFLib_thread_id] == NULL);
	assert(func_smodel->Kmodel->theta_extra[GMRFLib_thread_id] == NULL);

	int idx, i, j, n = func_smodel->n, where_to_start;
	double mean, var;

	where_to_start = GMRFLib_ai_INLA_userfunc1_dim;

	GMRFLib_ai_INLA_userfunc1_dim = 2 * n;		       /* set this variable, yes */
	GMRFLib_ai_INLA_userfunc1_density = Calloc(GMRFLib_ai_INLA_userfunc1_dim, GMRFLib_density_tp *);

	for (idx = 0; idx < n; idx++) {
		inla_spde_theta_tp *t;
		int offset, offset_theta;

		t = func_smodel->Tmodel;
		offset = 0;
		offset_theta = 0;
		DO_COMPUTE;

		t = func_smodel->Kmodel;
		offset = n;
		offset_theta = func_smodel->Tmodel->ntheta;
		DO_COMPUTE;
	}

#undef DO_COMPUTE
	return NULL;
}
