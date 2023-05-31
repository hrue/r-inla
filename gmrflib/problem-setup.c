
/* problem-setup.c
 * 
 * Copyright (C) 2001-2023 Havard Rue
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
#include <math.h>
#include <string.h>
#include <stdio.h>

#include "GMRFLib/sha.h"
#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

int error_check_validate_constr1 = 0;

#if defined(INLA_WINDOWS32)
static int constr_store_use = 0;			       /* do not use it as the sha is not prepared for it */
#else
static int constr_store_use = 1;
#endif
static map_strvp constr_store;
static int constr_store_must_init = 1;
static int constr_store_debug = 0;

int GMRFLib_init_constr_store(void)
{
	GMRFLib_ENTER_ROUTINE;
	constr_store_debug = GMRFLib_DEBUG_IF_TRUE();
	if (constr_store_use) {
		if (constr_store_must_init) {
			map_strvp_init_hint(&constr_store, 128);
			constr_store.alwaysdefault = 1;
			constr_store_must_init = 0;
			if (constr_store_debug) {
				printf("constr_store: init storage\n");
			}
		}
	}
	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}


#if defined(INLA_WINDOWS32)
static int constr_store_logdet_use = 0;			       /* do not use it as the sha is not prepared for it */
#else
static int constr_store_logdet_use = 1;
#endif
static map_strd constr_store_logdet;
static int constr_store_logdet_must_init = 1;
static int constr_store_logdet_debug = 0;

int GMRFLib_init_constr_store_logdet(void)
{
	GMRFLib_ENTER_ROUTINE;
	constr_store_logdet_debug = GMRFLib_DEBUG_IF_TRUE();
	if (constr_store_logdet_use) {
		if (constr_store_logdet_must_init) {
			map_strd_init_hint(&constr_store_logdet, 128);
			constr_store_logdet.alwaysdefault = 1;
			constr_store_logdet_must_init = 0;
			if (constr_store_logdet_debug) {
				printf("constr_store_logdet: init storage\n");
			}
		}
	}
	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

double GMRFLib_Qfunc_wrapper(int thread_id, int sub_node, int sub_nnode, double *values, void *arguments)
{
	int node, nnode;
	double val;
	GMRFLib_Qfunc_arg_tp *args = NULL;
	args = (GMRFLib_Qfunc_arg_tp *) arguments;

	// we know that the mapping is identity in this case and that the sub_graph is the same as the graph
	// this is also validated in the problem-setup

	node = sub_node;
	nnode = sub_nnode;

	if (nnode >= 0) {
		// the normal case, nothing spesific to do
		if (node == nnode) {
			val = (*(args->user_Qfunc)) (thread_id, node, nnode, NULL, args->user_Qfunc_args) + args->diagonal_adds[node];
		} else {
			val = (*(args->user_Qfunc)) (thread_id, node, nnode, NULL, args->user_Qfunc_args);
		}
	} else {
		// this is the multi-case
		val = (*(args->user_Qfunc)) (thread_id, node, -1, values, args->user_Qfunc_args);
		if (ISNAN(val)) {
			// the Qfunction does not support it, move on doing it manually
			int j, jj, k = 0;
			values[k++] = (*(args->user_Qfunc)) (thread_id, node, node, NULL, args->user_Qfunc_args) + args->diagonal_adds[node];
			for (jj = 0; jj < args->graph->lnnbs[node]; jj++) {
				j = args->graph->lnbs[node][jj];
				values[k++] = (*(args->user_Qfunc)) (thread_id, node, j, NULL, args->user_Qfunc_args);
			}
		} else {
			values[0] += args->diagonal_adds[node];
		}
		val = 0.0;
	}
	return val;
}

int validate_constr1(GMRFLib_constr_tp *constr, int n)
{
	GMRFLib_constr_tp *nnew = NULL;
	GMRFLib_graph_tp *g = Calloc(1, GMRFLib_graph_tp);
	g->n = n;

	GMRFLib_duplicate_constr(&nnew, constr, g);
	assert(nnew);
	for (int j = 0; j < constr->nc; j++) {
		if ((nnew->jfirst[j] != constr->jfirst[j]) || (nnew->jlen[j] != constr->jlen[j])) {
			printf("CONSTR jfirst/jlen ERROR: i= %d nnew->jfirst= %d old->jfirst= %d nnew->jlen= %d old->jlen= %d\n",
			       j, nnew->jfirst[j], constr->jfirst[j], nnew->jlen[j], constr->jlen[j]);
			abort();
		}
	}
	GMRFLib_free_constr(nnew);
	Free(g);

	return GMRFLib_SUCCESS;
}

int dgemm_special(int m, int n, double *C, double *UNUSED(A), double *B, GMRFLib_constr_tp *constr)
{
	if (error_check_validate_constr1)
		validate_constr1(constr, n);

	// compute C=A*B, where A is the constr matrix, and we know that C is symmetric.
	// see below where this is used.

	typedef struct {
		int m;
		int K;
		int *ii;
		int *jj;
	} storage_t;

	static storage_t **storage = NULL;

	if (!storage) {
#pragma omp critical (Name_e4b888063524c281c8bef772bc2579873731fa49)
		{
			if (!storage) {
				storage = Calloc(GMRFLib_CACHE_LEN(), storage_t *);
			}
		}
	}

	int id = 0;
	GMRFLib_CACHE_SET_ID(id);

	if (!storage[id]) {
		storage[id] = Calloc(1, storage_t);
	}

	if (storage[id]->m != m) {
		Free(storage[id]->ii);

		int K = m * (m + 1) / 2;
		int *ii = Calloc(2 * K, int);
		int *jj = ii + K;
		for (int k = 0, i = 0; i < m; i++) {
			for (int j = i; j < m; j++) {
				ii[k] = i;
				jj[k] = j;
				k++;
			}
		}
		storage[id]->m = m;
		storage[id]->K = K;
		storage[id]->ii = ii;
		storage[id]->jj = jj;
	}
#define CODE_BLOCK							\
	for (int k = 0; k < storage[id]->K; k++) {			\
		int i = storage[id]->ii[k], j = storage[id]->jj[k];	\
		double value;						\
		value = GMRFLib_dot_product(constr->idxval[i], B + j * n); \
		C[i + j * m] = C[j + i * m] = value;			\
	}

	int nt = 0;
	if (omp_get_level() > 2) {
		nt = 1;
	} else {
		nt = (m > GMRFLib_MAX_THREADS()? GMRFLib_MAX_THREADS() : 1);
	}

	RUN_CODE_BLOCK(nt, 0, 0);
#undef CODE_BLOCK

	return GMRFLib_SUCCESS;
}

int dgemm_special2(int m, double *C, double *A, GMRFLib_constr_tp *constr)
{
	// compute C=A*A', where A is the constr matrix. C is symmetric. see below where this is used.

	typedef struct {
		int m;
		int K;
		int *ii;
		int *jj;
	} storage_t;

	static storage_t **storage = NULL;

	if (!storage) {
#pragma omp critical (Name_756fb3e47fb1205cd5e595775867506d805d78f6)
		{
			if (!storage) {
				storage = Calloc(GMRFLib_CACHE_LEN(), storage_t *);
			}
		}
	}

	int id = 0;
	GMRFLib_CACHE_SET_ID(id);

	if (!storage[id]) {
		storage[id] = Calloc(1, storage_t);
	}

	if (storage[id]->m != m) {
		Free(storage[id]->ii);

		int K = m * (m + 1) / 2;
		int *ii = Calloc(2 * K, int);
		int *jj = ii + K;
		for (int k = 0, i = 0; i < m; i++) {
			for (int j = i; j < m; j++) {
				ii[k] = i;
				jj[k] = j;
				k++;
			}
		}
		storage[id]->m = m;
		storage[id]->K = K;
		storage[id]->ii = ii;
		storage[id]->jj = jj;
	}
#define CODE_BLOCK							\
	for (int k = 0; k < storage[id]->K; k++) {			\
		int i = storage[id]->ii[k], j = storage[id]->jj[k], incx = m, jf, je, jlen; \
		double value;						\
		jf = IMAX(constr->jfirst[i], constr->jfirst[j]);	\
		je = IMIN(constr->jfirst[i] + constr->jlen[i], constr->jfirst[j] + constr->jlen[j]); \
		jlen = je - jf;						\
		if (jlen > 0) {						\
			value = ddot_(&jlen, &(A[i + m * jf]), &incx, &(A[j + m * jf]), &incx);	\
		} else {						\
			value = 0.0;					\
		}							\
		C[i + j * m] = C[j + i * m] = value;			\
	}

	int nt = 0;
	if (omp_get_level() > 2) {
		nt = 1;
	} else {
		nt = (m > GMRFLib_MAX_THREADS()? GMRFLib_MAX_THREADS() : 1);
	}

	RUN_CODE_BLOCK(nt, 0, 0);
#undef CODE_BLOCK

	return GMRFLib_SUCCESS;
}

int dgemv_special(double *res, double *x, GMRFLib_constr_tp *constr)
{
	// compute 'res = A %*% x'

	int nc = constr->nc;
	int nt = 0;

	// can be used on third level
	if (omp_get_level() > 2) {
		nt = 1;
	} else {
		int f = 64;
		nt = IMAX(1, IMIN(nc / f, GMRFLib_MAX_THREADS()));
	}

#define CODE_BLOCK							\
	for (int i = 0; i < nc; i++) {					\
		res[i] = GMRFLib_dot_product(constr->idxval[i], x);	\
	}

	RUN_CODE_BLOCK(nt, 0, 0);
#undef CODE_BLOCK

	return GMRFLib_SUCCESS;
}

int GMRFLib_Qsolve(double *x, double *b, GMRFLib_problem_tp *problem, int idx)
{
	// solve Q x = b, but correct for constraints, like eq 2.30 in the GMRF-book, x := x - Q^-1A^T(AQ^-1A^T)^-1 (Ax-e).

	// if IDX >=0 then assume only b[idx] != 0, if IDX < 0, then assume a general B

	GMRFLib_ENTER_ROUTINE;

	static double **wwork = NULL;
	static int *wwork_len = NULL;
	if (!wwork) {
#pragma omp critical (Name_415e5dc8137b4c7fc43dd511591be71a0710d398)
		{
			if (!wwork) {
				wwork_len = Calloc(GMRFLib_CACHE_LEN(), int);
				wwork = Calloc(GMRFLib_CACHE_LEN(), double *);
			}
		}
	}

	int n = problem->sub_graph->n;
	int nc = (problem->sub_constr && problem->sub_constr->nc > 0 ? problem->sub_constr->nc : 0);
	double *xx = NULL;
	int cache_idx = 0;

	GMRFLib_CACHE_SET_ID(cache_idx);
	if (n + nc > wwork_len[cache_idx]) {
		Free(wwork[cache_idx]);
		wwork_len[cache_idx] = n + nc;
		wwork[cache_idx] = Calloc(wwork_len[cache_idx], double);
	}
	xx = wwork[cache_idx];

	Memcpy(xx, b, n * sizeof(double));
	if (idx >= 0) {
		GMRFLib_solve_llt_sparse_matrix_special(xx, &(problem->sub_sm_fact), problem->sub_graph, idx);
	} else {
		GMRFLib_solve_llt_sparse_matrix(xx, 1, &(problem->sub_sm_fact), problem->sub_graph);
	}

	if ((problem->sub_constr && problem->sub_constr->nc > 0)) {
		int nnc = problem->sub_constr->nc, inc = 1;
		double alpha = -1.0, beta = 1.0, *t_vector = wwork[cache_idx] + n;
		GMRFLib_eval_constr0(t_vector, NULL, xx, problem->sub_constr, problem->sub_graph);
		dgemv_("N", &n, &nnc, &alpha, problem->constr_m, &n, t_vector, &inc, &beta, xx, &inc, F_ONE);
	}

	Memcpy(x, xx, n * sizeof(double));

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_init_problem(int thread_id, GMRFLib_problem_tp **problem,
			 double *x,
			 double *b,
			 double *c, double *mean, GMRFLib_graph_tp *graph, GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_args, GMRFLib_constr_tp *constr)
{
	GMRFLib_ENTER_ROUTINE;
	GMRFLib_EWRAP1(GMRFLib_init_problem_store(thread_id, problem, x, b, c, mean, graph, Qfunc, Qfunc_args, constr, NULL));
	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_init_problem_store(int thread_id,
			       GMRFLib_problem_tp **problem,
			       double *x,
			       double *b,
			       double *c,
			       double *mean,
			       GMRFLib_graph_tp *graph,
			       GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_args, GMRFLib_constr_tp *constr, GMRFLib_store_tp *store)
{
	double *bb = NULL;
	int sub_n, free_x = 0, retval;
	GMRFLib_smtp_tp smtp;

	int store_store_sub_graph = 0, store_use_sub_graph = 0;
	int store_store_remap = 0, store_use_remap = 0;
	int store_store_symb_fact = 0, store_use_symb_fact = 0;

	GMRFLib_Qfunc_tp *sub_Qfunc;
	GMRFLib_Qfunc_arg_tp *sub_Qfunc_arg;

	GMRFLib_ENTER_ROUTINE;

	GMRFLib_ASSERT(problem, GMRFLib_EINVARG);
	GMRFLib_ASSERT(graph, GMRFLib_EINVARG);
	GMRFLib_ASSERT(Qfunc, GMRFLib_EINVARG);

	if (store) {
		store_store_sub_graph = (store->sub_graph ? 0 : 1);
		store_use_sub_graph = !store_store_sub_graph;
		store_store_remap = (store->remap ? 0 : 1);
		store_use_remap = !store_store_remap;
	}

	/*
	 * if the smtp is defined in the store, we use that one. otherwise, we use the default one. 
	 */
	if (store && GMRFLib_valid_smtp((int) store->smtp) == GMRFLib_TRUE) {
		smtp = store->smtp;
	} else {
		smtp = GMRFLib_smtp;
	}
	if (store) {
		if (smtp == GMRFLib_SMTP_TAUCS) {
			store_store_symb_fact = (store && store->TAUCS_symb_fact ? 0 : 1);
			store_use_symb_fact = !store_store_symb_fact;
		} else if (smtp == GMRFLib_SMTP_PARDISO) {
			store_store_symb_fact = 0;
			store_use_symb_fact = !store_store_symb_fact;
		} else {
			store_store_symb_fact = 0;
			store_use_symb_fact = 0;
		}
	}

	/*
	 * if x = NULL, make it zeros 
	 */
	if (!x) {
		x = Calloc(graph->n, double);
		free_x = 1;
	}

	*problem = Calloc(1, GMRFLib_problem_tp);

	/*
	 * this is the sparse-matrix method 
	 */
	(*problem)->sub_sm_fact.smtp = smtp;

	/*
	 * compute the internal stuff if needed 
	 */
	if (constr) {
		GMRFLib_EWRAP1(GMRFLib_prepare_constr(constr, graph, 0));
	}

	if (store_use_sub_graph) {
		/*
		 * copy from store 
		 */
		GMRFLib_graph_duplicate(&((*problem)->sub_graph), store->sub_graph);
	} else {
		/*
		 * compute it 
		 */
		GMRFLib_graph_comp_subgraph(&((*problem)->sub_graph), graph, NULL, NULL);

		/*
		 * store a copy, if requested 
		 */
		if (store_store_sub_graph) {
			GMRFLib_graph_duplicate(&(store->sub_graph), (*problem)->sub_graph);
		}
	}

	sub_n = (*problem)->sub_graph->n;
	if (sub_n == 0) {				       /* fast return if there is nothing todo */
		GMRFLib_graph_free((*problem)->sub_graph);
		Free(*problem);
		if (free_x) {
			Free(x);
		}
		GMRFLib_LEAVE_ROUTINE;
		return GMRFLib_SUCCESS;
	}

	/*
	 * compute the reordering which is stored in sub_sm_fact 
	 */
	if (store_use_remap) {
		/*
		 * use the reordering in store 
		 */
		(*problem)->sub_sm_fact.remap = Calloc(sub_n, int);
		Memcpy((*problem)->sub_sm_fact.remap, store->remap, sub_n * sizeof(int));
		if (smtp == GMRFLib_SMTP_BAND) {
			(*problem)->sub_sm_fact.bandwidth = store->bandwidth;
		}
	} else {
		/*
		 * compute it 
		 */
		GMRFLib_EWRAP1(GMRFLib_compute_reordering(&((*problem)->sub_sm_fact), (*problem)->sub_graph, NULL));

		/*
		 * store a copy, if requested 
		 */
		if (store_store_remap) {
			if ((*problem)->sub_sm_fact.remap != NULL) {
				store->remap = Calloc(sub_n, int);
				Memcpy(store->remap, (*problem)->sub_sm_fact.remap, sub_n * sizeof(int));
				if (smtp == GMRFLib_SMTP_BAND) {
					store->bandwidth = (*problem)->sub_sm_fact.bandwidth;
				}
			} else {
				store->remap = NULL;
			}
		}
	}

	/*
	 * setup space 
	 */
	(*problem)->n = graph->n;
	(*problem)->sample = Calloc(graph->n, double);
	(*problem)->mean = Calloc(graph->n, double);
	(*problem)->mean_constr = Calloc(graph->n, double);
	(*problem)->sub_sample = Calloc(sub_n, double);
	(*problem)->sub_mean = Calloc(sub_n, double);
	(*problem)->sub_mean_constr = Calloc(sub_n, double);

	Memcpy((*problem)->sample, x, graph->n * sizeof(double));
	Memcpy((*problem)->mean, x, graph->n * sizeof(double));
	Memcpy((*problem)->mean_constr, x, graph->n * sizeof(double));

	/*
	 * tabulate the Qfunc on the sub_graph. first, make the Qfunc, then tabulate it. 
	 */

	sub_Qfunc = GMRFLib_Qfunc_wrapper;
	sub_Qfunc_arg = Calloc(1, GMRFLib_Qfunc_arg_tp);
	sub_Qfunc_arg->diagonal_adds = Calloc(sub_n, double);
	sub_Qfunc_arg->graph = (*problem)->sub_graph;

	if (c) {
		Memcpy(sub_Qfunc_arg->diagonal_adds, c, sub_n * sizeof(double));
	}
	sub_Qfunc_arg->user_Qfunc = Qfunc;
	sub_Qfunc_arg->user_Qfunc_args = Qfunc_args;
	GMRFLib_EWRAP1(GMRFLib_tabulate_Qfunc(thread_id, &((*problem)->tab), (*problem)->sub_graph, sub_Qfunc, (void *) sub_Qfunc_arg, NULL));

	Free(sub_Qfunc_arg->diagonal_adds);
	Free(sub_Qfunc_arg);

	/*
	 * now compute the new 'effective' b, and then the mean 
	 */

	/*
	 * i use bb as name, pointing to the same storage as sub_mean. this to avoid yet another arraw. but i assume
	 * that sub_mean[] is zero, which is the case for !keep, but NOT the case otherwise. hence, i have to set it to 
	 * zero unless it already is so. [this was a nasty bug...] 
	 */
	bb = (*problem)->sub_mean;

	if (b) {
		Memcpy(bb, b, sub_n * sizeof(double));
	}

	if (mean) {
		double *tmp = Calloc(sub_n, double);
		GMRFLib_Qx(thread_id, tmp, mean, (*problem)->sub_graph, (*problem)->tab->Qfunc, (*problem)->tab->Qfunc_arg);
		GMRFLib_daddto(sub_n, tmp, bb);
		Free(tmp);
	}

	/*
	 * now compute the Cholesky-factorisation
	 * 
	 * solve to obtain the mean (recall that bb=(*problem)->sub_mean) later! 
	 */
	/*
	 * 13/5/2005. the next is a small hack, to be fixed properly later, well, it is not decided (yet,) if this is
	 * really needed.
	 * 
	 * in the case where graph is fixed and we're using TAUCS, we can reuse the symbolic factorisation. the
	 * symbolic, is roughly, about 10% of the costs of the numeric factorisation for large problems. however, the
	 * 'free_fact_sparse_...' does free both, and there is currently no options to free only the numerical one. so
	 * therefore i do this here. 
	 */

	if (store_use_symb_fact && (smtp == GMRFLib_SMTP_TAUCS)) {
		(*problem)->sub_sm_fact.TAUCS_symb_fact = GMRFLib_sm_fact_duplicate_TAUCS(store->TAUCS_symb_fact);
	}

	if (store_use_symb_fact && (smtp == GMRFLib_SMTP_PARDISO)) {
		// FIXME1("ADDED NEW EXPERIMENTAL CODE");
		GMRFLib_pardiso_store_tp *s = Calloc(1, GMRFLib_pardiso_store_tp);
		s->graph = (*problem)->sub_graph;
		// use the internal cached storage
		GMRFLib_duplicate_pardiso_store(&((*problem)->sub_sm_fact.PARDISO_fact), s, GMRFLib_FALSE, GMRFLib_FALSE);
		Free(s);
	}

	int ret;
	ret = GMRFLib_build_sparse_matrix(thread_id, &((*problem)->sub_sm_fact), (*problem)->tab->Qfunc,
					  (char *) ((*problem)->tab->Qfunc_arg), (*problem)->sub_graph);
	if (ret != GMRFLib_SUCCESS) {
		return ret;
	}

	ret = GMRFLib_factorise_sparse_matrix(&((*problem)->sub_sm_fact), (*problem)->sub_graph);
	if (ret != GMRFLib_SUCCESS) {
		return ret;
	}

	if (store_store_symb_fact && (smtp == GMRFLib_SMTP_TAUCS)) {
		store->TAUCS_symb_fact = GMRFLib_sm_fact_duplicate_TAUCS((*problem)->sub_sm_fact.TAUCS_symb_fact);
	}

	/*
	 * the next step is to initialize the constraints and take that part into account. 
	 */
	if (constr && constr->nc > 0) {
		int nc;
		double *aat_m, alpha, beta, *b_add = NULL, *aqat_m;

		/*
		 * first make new constraints on the subgraph. ok to call this function if fixed_values is NULL, then
		 * we just get a copy of constr back. 
		 */
		b_add = Calloc((*problem)->sub_graph->n, double);

		if (sub_n == graph->n) {
			GMRFLib_duplicate_constr(&((*problem)->sub_constr), constr, graph);
		} else {
			GMRFLib_recomp_constr(&((*problem)->sub_constr), constr, x, b_add, NULL, graph, (*problem)->sub_graph);
		}

		/*
		 * if we should keep the mean, then do not add the correction-terms 
		 */
		GMRFLib_daddto((*problem)->sub_graph->n, b_add, (*problem)->sub_mean);
		Free(b_add);

		if ((*problem)->sub_constr && (*problem)->sub_constr->nc > 0) {
			/*
			 * go further only if the constraint is still there: it might go away!!! 
			 */
			nc = (*problem)->sub_constr->nc;       /* shortname */

			// we assume this is ok for INLA so we turn this off.

			double *qi_at_m_store = NULL;	       /* possible reuse old results */

			(*problem)->qi_at_m = Calloc(nc * sub_n, double);

			if (qi_at_m_store == NULL) {
				/*
				 * compute it as usual 
				 */
				for (int k = 0; k < nc; k++) {
					int kk = k * sub_n;
					double *yy = (*problem)->qi_at_m + kk;
					double *xx = (*problem)->sub_constr->a_matrix + k; 
#pragma GCC ivdep
					for (int i = 0, j = 0; i < sub_n; i++, j += nc) {
						yy[i] = xx[j];
					}
				}
				GMRFLib_EWRAP1(GMRFLib_solve_llt_sparse_matrix
					       ((*problem)->qi_at_m, nc, &((*problem)->sub_sm_fact), (*problem)->sub_graph));
			} else {
				/*
				 * reuse 
				 */
				Memcpy((*problem)->qi_at_m, qi_at_m_store, (nc - 1) * sub_n * sizeof(double));
				for (int k = nc - 1; k < nc; k++) {
					int kk = k * sub_n;
					double *yy = (*problem)->qi_at_m + kk;
					double *xx = (*problem)->sub_constr->a_matrix + k; 
#pragma GCC ivdep
					for (int i = 0, j = 0; i < sub_n; i++, j += nc) {
						yy[i] = xx[j];
					}
				}
				GMRFLib_EWRAP1(GMRFLib_solve_llt_sparse_matrix(&((*problem)->qi_at_m[(nc - 1) * sub_n]), 1,
									       &((*problem)->sub_sm_fact), (*problem)->sub_graph));
			}
			Free(qi_at_m_store);

			/*
			 * compute l_aqat_m = chol(AQ^{-1}A^T)^{-1}) = chol(A qi_at_m)^{-1}, size = nc x nc 
			 */
			aqat_m = Calloc(nc * nc, double);
			alpha = 1.0;
			beta = 0.0;
			if (GMRFLib_faster_constr) {
				dgemm_special(nc, sub_n, aqat_m, (*problem)->sub_constr->a_matrix, (*problem)->qi_at_m, (*problem)->sub_constr);
			} else {
				dgemm_("N", "N", &nc, &nc, &sub_n, &alpha, (*problem)->sub_constr->a_matrix, &nc,
				       (*problem)->qi_at_m, &sub_n, &beta, aqat_m, &nc, F_ONE, F_ONE);
			}

			if (GMRFLib_aqat_m_diag_add > 0.0) {
				for (int i = 0; i < nc; i++) {
					aqat_m[i + i * nc] += GMRFLib_aqat_m_diag_add;
				}
			}

			/*
			 * compute chol(aqat_m), recall that GMRFLib_comp_chol_general returns a new malloced L 
			 */
			retval = GMRFLib_comp_chol_general(&((*problem)->l_aqat_m), aqat_m, nc, &((*problem)->logdet_aqat), GMRFLib_ESINGCONSTR);
			if (retval != GMRFLib_SUCCESS) {
				GMRFLib_WARNING("Matrix AQA^t is numerical singular, remove singularity and move on");
				GMRFLib_ensure_spd(aqat_m, nc, GSL_SQRT_DBL_EPSILON, NULL);
				GMRFLib_EWRAP1(GMRFLib_comp_chol_general
					       (&((*problem)->l_aqat_m), aqat_m, nc, &((*problem)->logdet_aqat), GMRFLib_ESINGCONSTR));
			}
			Free(aqat_m);

			/*
			 * ...and the constr-matrix Q^-1A^T inv(AQ^{-1}A^T + Sigma) 
			 */
			(*problem)->constr_m = Calloc(sub_n * nc, double);
			double *tmp_vector = Calloc(sub_n * nc, double);

			for (int j = 0; j < nc; j++) {
				double *yy = tmp_vector + j;
				double *xx = (*problem)->qi_at_m + j * sub_n;
#pragma GCC ivdep
				for (int i = 0, k = 0; i < sub_n; i++, k += nc) {
					yy[k] = xx[i];
				}
			}

			GMRFLib_solveAxb_posdef(tmp_vector, (*problem)->l_aqat_m, tmp_vector, nc, sub_n);

			for (int j = 0; j < nc; j++) {
				double *yy = (*problem)->constr_m + j * sub_n;
				double *xx = tmp_vector + j;
#pragma GCC ivdep
				for (int i = 0, k = 0; i < sub_n; i++, k += nc) {
					yy[i] = xx[k];
				}
			}
				
			Free(tmp_vector);

			GMRFLib_constr_tp *con = (*problem)->sub_constr;
			double *p = NULL;

			if (con->sha && constr_store_logdet_use) {
				p = map_strd_ptr(&constr_store_logdet, (char *) con->sha);
				if (p) {
					(*problem)->logdet_aat = *p;
				}
				if (constr_store_logdet_debug) {
					if (p) {
						printf("constr_store_logdet: constr found in store= %f\n", *p);
					} else {
						printf("constr_store_logdet: constr not found in store\n");
					}
				}
			}
			if (!p) {
				/*
				 * compute |A*A'| 
				 */
				alpha = 1.0;
				beta = 0.0;
				aat_m = Calloc(nc * nc, double);

				if (GMRFLib_faster_constr) {
					dgemm_special2(nc, aat_m, (*problem)->sub_constr->a_matrix, (*problem)->sub_constr);
				} else {
					dgemm_("N", "T", &nc, &nc, &sub_n, &alpha, (*problem)->sub_constr->a_matrix,
					       &nc, (*problem)->sub_constr->a_matrix, &nc, &beta, aat_m, &nc, F_ONE, F_ONE);
				}
				tmp_vector = NULL;
				retval = GMRFLib_comp_chol_general(&tmp_vector, aat_m, nc, &((*problem)->logdet_aat), GMRFLib_ESINGCONSTR2);
				if (retval != GMRFLib_SUCCESS) {
					GMRFLib_WARNING("Matrix AA^t is numerical singular, remove singularity and move on");
					GMRFLib_ensure_spd(aat_m, nc, GSL_SQRT_DBL_EPSILON, NULL);
					GMRFLib_EWRAP1(GMRFLib_comp_chol_general
						       (&tmp_vector, aat_m, nc, &((*problem)->logdet_aat), GMRFLib_ESINGCONSTR2));
				}
				Free(aat_m);
				Free(tmp_vector);

				if (constr_store_logdet_use) {
					if (!(con->sha)) {
						if (constr_store_logdet_debug) {
							printf("constr_store_logdet: value computed %f, but not set\n", (*problem)->logdet_aat);
						}
					} else if (con->sha) {
						if (constr_store_logdet_debug) {
							printf("constr_store_logdet: store value %f\n", (*problem)->logdet_aat);
						}
#pragma omp critical (Name_8c313c5cb0ba5eb20ede5a81e455580200ca1348)
						{
							map_strd_set(&constr_store_logdet, GMRFLib_strdup((char *) con->sha),
								     (*problem)->logdet_aat);
						}
					}
				}
			}
		}
	}

	GMRFLib_EWRAP1(GMRFLib_solve_llt_sparse_matrix((*problem)->sub_mean, 1, &((*problem)->sub_sm_fact), (*problem)->sub_graph));
	if (!((*problem)->sub_mean_constr)) {
		(*problem)->sub_mean_constr = Calloc(sub_n, double);
	}
	Memcpy((*problem)->sub_mean_constr, (*problem)->sub_mean, sub_n * sizeof(double));

	if (((*problem)->sub_constr && (*problem)->sub_constr->nc > 0)) {
		/*
		 * compute the mean after correcting for the constraint. this is the same as if the sample itself is
		 * sub_mean! i have copied parts of code from GMRFLib_sample into here...
		 */
		int nc = constr->nc, inc = 1;
		double alpha, beta, *t_vector;

		Free((*problem)->sub_constr_value);
		(*problem)->sub_constr_value = t_vector = Calloc(nc, double);

		GMRFLib_eval_constr(t_vector, NULL, (*problem)->sub_mean, (*problem)->sub_constr, (*problem)->sub_graph);

		/*
		 * sub_mean_constr is pr.default equal to sub_mean 
		 */
		alpha = -1.0;
		beta = 1.0;				       /* mean_constr = mean - cond_m*t_vector */
		dgemv_("N", &sub_n, &nc, &alpha, (*problem)->constr_m, &sub_n, t_vector, &inc, &beta, (*problem)->sub_mean_constr, &inc, F_ONE);
	}

	/*
	 * make the ``mean variables'', which are easier accessible for the user. 
	 */
	Memcpy((*problem)->mean, (*problem)->sub_mean, sub_n * sizeof(double));
	Memcpy((*problem)->mean_constr, (*problem)->sub_mean_constr, sub_n * sizeof(double));

	/*
	 * evaluate the log-likelihood and compute the constants therein 
	 */
	GMRFLib_EWRAP1(GMRFLib_evaluate__intern(*problem, 1));

	if (free_x) {
		Free(x);
	}

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_sample(GMRFLib_problem_tp *problem)
{
	int i, n;
	double sqrterm;

	if (!problem)
		return GMRFLib_SUCCESS;
	GMRFLib_ENTER_ROUTINE;

	n = problem->sub_graph->n;

	/*
	 * just solve L^Tx=z, then add mean and store in correct place 
	 */
	for (i = 0, sqrterm = 0.0; i < n; i++) {
		double z = GMRFLib_stdnormal();
		sqrterm += SQR(z);
		problem->sub_sample[i] = z;
	}

	GMRFLib_EWRAP1(GMRFLib_solve_lt_sparse_matrix(problem->sub_sample, 1, &(problem->sub_sm_fact), problem->sub_graph));
	GMRFLib_daddto(n, problem->sub_mean, problem->sub_sample);
	Memcpy(problem->sample, problem->sub_sample, n * sizeof(double));

	/*
	 * if there is no constraints, then we can evaluate the log-density and return 
	 */
	if (!problem->sub_constr || (problem->sub_constr && problem->sub_constr == 0)) {
		problem->sub_logdens = -0.5 * n * log(2.0 * M_PI) + problem->log_normc - 0.5 * sqrterm;
		GMRFLib_LEAVE_ROUTINE;
		return GMRFLib_SUCCESS;
	} else {
		/*
		 * otherwise, we must account for the constraint, BUT, we need to call GMRFLib_evaluate() to compute the
		 * log-density.
		 * 
		 * FIXME: There is a slight overhead, since the quadratic form can be faster evaluated knowing the random
		 * bits... 
		 */

		double alpha, beta, *t_vector;
		int inc = 1, nc = problem->sub_constr->nc;

		Free(problem->sub_constr_value);
		problem->sub_constr_value = Calloc(nc, double);
		t_vector = Calloc(nc, double);		       /* t_vector = Ax-e */

		GMRFLib_EWRAP1(GMRFLib_eval_constr(problem->sub_constr_value, NULL, problem->sub_sample, problem->sub_constr, problem->sub_graph));
		Memcpy(t_vector, problem->sub_constr_value, nc * sizeof(double));

		alpha = -1.0;
		beta = 1.0;				       /* sample := sample - cond_m*t_vector */
		dgemv_("N", &n, &nc, &alpha, problem->constr_m, &n, t_vector, &inc, &beta, problem->sub_sample, &inc, F_ONE);
		Free(t_vector);
		Memcpy(problem->sample, problem->sub_sample, n * sizeof(double));
		GMRFLib_EWRAP1(GMRFLib_evaluate(problem));     /* to compute the log-density */
	}

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_evaluate(GMRFLib_problem_tp *problem)
{
	GMRFLib_ENTER_ROUTINE;
	GMRFLib_EWRAP1(GMRFLib_evaluate__intern(problem, 0));
	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_evaluate__intern(GMRFLib_problem_tp *problem, int compute_const)
{
	/*
	 * evaluate the log-density in point 'sample' in the problem definition 
	 */

	int i, n;
	int thread_id = -1;
	double sqrterm, *xx = NULL, *yy = NULL, *work = NULL;

	if (!problem) {
		return GMRFLib_SUCCESS;
	}

	n = problem->sub_graph->n;
	work = Calloc(2 * n, double);
	xx = work;
	yy = work + n;

	/*
	 * user has altered the 'sample', put the correct subset into sub_sample and compute (x-\mu)^TQ(x-\mu)
	 */

	Memcpy(problem->sub_sample, problem->sample, n * sizeof(double));
#pragma omp simd private(i)
	for (i = 0; i < n; i++) {
		xx[i] = problem->sub_sample[i] - problem->sub_mean[i];
	}
	GMRFLib_Qx(thread_id, yy, xx, problem->sub_graph, problem->tab->Qfunc, (void *) problem->tab->Qfunc_arg);
	sqrterm = GMRFLib_ddot(n, yy, xx);

	Free(work);

	/*
	 * evaluate the normalization constant and add up 
	 */
	if (compute_const) {
		GMRFLib_EWRAP0(GMRFLib_log_determinant(&(problem->log_normc), &(problem->sub_sm_fact), problem->sub_graph));
		problem->log_normc /= 2.0;		       /* |Q|^1/2 */
	}

	problem->sub_logdens = -0.5 * n * log(2.0 * M_PI) + problem->log_normc - 0.5 * sqrterm;

	/*
	 * now correct for the constraint, if any 
	 */
	if (!problem->sub_constr || problem->sub_constr->nc == 0) {
		return GMRFLib_SUCCESS;
	} else {
		/*
		 * deterministic constraints 
		 */
		int nc = problem->sub_constr->nc;

		if (compute_const) {
			/*
			 * t_vector = A mu-b tt_vector = i_aqat_m*t_vector 
			 */
			double *t_vector = NULL, *tt_vector = NULL;

			Free(problem->sub_constr_value);
			problem->sub_constr_value = t_vector = Calloc(nc, double);
			tt_vector = Calloc(nc, double);

			GMRFLib_EWRAP0(GMRFLib_eval_constr(t_vector, NULL, problem->sub_mean, problem->sub_constr, problem->sub_graph));
			GMRFLib_EWRAP0(GMRFLib_solveAxb_posdef(tt_vector, problem->l_aqat_m, t_vector, nc, 1));
			problem->exp_corr = GMRFLib_ddot(nc, t_vector, tt_vector);
			Free(tt_vector);
		}

		/*
		 * [x|Ax] = [x] [Ax|x] / [Ax] 
		 */
		problem->sub_logdens += (-0.5 * problem->logdet_aat)	/* [Ax|x] */
		    -(-0.5 * nc * log(2.0 * M_PI) - 0.5 * problem->logdet_aqat - 0.5 * problem->exp_corr);	/* [Ax] */
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_free_problem(GMRFLib_problem_tp *problem)
{
	/*
	 * free all malloced stuff in 'problem' 
	 */
	if (!problem) {
		return GMRFLib_SUCCESS;
	}
	GMRFLib_free_Qinv(problem);			       /* do this here, so that `n' is still there... */

	Free(problem->sample);
	Free(problem->mean);
	Free(problem->mean_constr);
	Free(problem->sub_sample);
	Free(problem->sub_mean);
	Free(problem->sub_mean_constr);
	GMRFLib_free_fact_sparse_matrix(&(problem->sub_sm_fact));
	GMRFLib_free_reordering(&(problem->sub_sm_fact));
	Free(problem->sub_constr_value);
	Free(problem->constr_m);
	Free(problem->l_aqat_m);
	Free(problem->inv_aqat_m);
	Free(problem->qi_at_m);
	GMRFLib_graph_free(problem->sub_graph);
	GMRFLib_free_tabulate_Qfunc(problem->tab);

	GMRFLib_free_constr(problem->sub_constr);
	problem->sub_constr = NULL;

	Free(problem);
	return GMRFLib_SUCCESS;
}

int GMRFLib_free_Qinv(GMRFLib_problem_tp *problem)
{
	if (problem && problem->sub_inverse) {
		int i, n = problem->sub_graph->n;

		for (i = 0; i < n; i++) {
			map_id_free(problem->sub_inverse->Qinv[i]);
			Free(problem->sub_inverse->Qinv[i]);
		}
		Free(problem->sub_inverse->Qinv);

		Free(problem->sub_inverse->mapping);
		Free(problem->sub_inverse);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_free_store(GMRFLib_store_tp *store)
{
	/*
	 * free contents in store 
	 */
	if (!store) {
		return GMRFLib_SUCCESS;
	}

	Free(store->remap);
	if (store->copy_ptr) {
		/*
		 * do nothing 
		 */
	} else {
		GMRFLib_graph_free(store->sub_graph);
		if (store->TAUCS_symb_fact) {
			taucs_supernodal_factor_free(store->TAUCS_symb_fact);
		}
	}

	if (store->copy_pardiso_ptr) {
		/*
		 * do nothing 
		 */
	} else {
		if (store->PARDISO_fact) {
			GMRFLib_pardiso_free(&(store->PARDISO_fact));
			store->PARDISO_fact = NULL;
		}
	}

	store->sub_graph = NULL;
	store->TAUCS_symb_fact = NULL;

	/*
	 * free the diag and sub-store. its of the same type, therefore we can do this recursively. 
	 */
	GMRFLib_free_store(store->sub_store);
	GMRFLib_free_store(store->diag_store);

	if (store->store_problems) {
		GMRFLib_free_problem(store->problem_old2new);
		GMRFLib_free_problem(store->problem_new2old);
		Free(store->old_logdens);
		Free(store->new_logdens);
	}

	/*
	 * free the ptr itself 
	 */
	Free(store);

	return GMRFLib_SUCCESS;
}

int GMRFLib_Qinv(GMRFLib_problem_tp *problem)
{
	if (problem) {
		GMRFLib_EWRAP1(GMRFLib_compute_Qinv((void *) problem));
	}
	return GMRFLib_SUCCESS;
}

double *GMRFLib_Qinv_get(GMRFLib_problem_tp *problem, int i, int j)
{
	int ii = problem->sub_inverse->mapping[i];
	int jj = problem->sub_inverse->mapping[j];
	return map_id_ptr(problem->sub_inverse->Qinv[IMIN(ii, jj)], IMAX(ii, jj));
}

double GMRFLib_Qinv_get0(GMRFLib_problem_tp *problem, int i, int j)
{
	int ii = problem->sub_inverse->mapping[i];
	int jj = problem->sub_inverse->mapping[j];
	double *d = map_id_ptr(problem->sub_inverse->Qinv[IMIN(ii, jj)], IMAX(ii, jj));
	if (d == NULL)
		printf("i j NULL %d %d\n", i, j);
	return (d ? *d : 0.0);
}

int GMRFLib_make_empty_constr(GMRFLib_constr_tp **constr)
{
	if (constr) {
		*constr = Calloc(1, GMRFLib_constr_tp);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_free_constr(GMRFLib_constr_tp *constr)
{
	if (!constr) {
		return GMRFLib_SUCCESS;
	}

	if (constr_store_use && constr->sha) {
		void *p;
		p = map_strvp_ptr(&constr_store, (char *) constr->sha);
		if (constr_store_debug) {
			if (p) {
				printf("\t[%1d] constr_store: constr is found in store: do not free\n", omp_get_thread_num());
			} else {
				printf("\t[%1d] constr_store: constr is not found in store: free\n", omp_get_thread_num());
			}
		}
		if (p) {
			return GMRFLib_SUCCESS;
		}
	}

	Free(constr->a_matrix);
	Free(constr->e_vector);
	Free(constr->jfirst);

	if (constr->idxval) {
		for (int i = 0; i < constr->nc; i++) {
			GMRFLib_idxval_free(constr->idxval[i]);
		}
		Free(constr->idxval);
	}

	Free(constr->sha);
	Free(constr);

	return GMRFLib_SUCCESS;
}

int GMRFLib_printf_constr(FILE *fp, GMRFLib_constr_tp *constr, GMRFLib_graph_tp *graph)
{
	int i, j;
	FILE *fpp = NULL;

	if (!constr) {
		return GMRFLib_SUCCESS;
	}
	fpp = (fp ? fp : stdout);

	fprintf(fpp, "n_constr %d  is_scaled %d\n", constr->nc, constr->is_scaled);
	for (i = 0; i < constr->nc; i++) {
		fprintf(fpp, "constraint %d, e= %f\na[%1d, ] = ", i, constr->e_vector[i], i);
		for (j = 0; j < graph->n; j++) {
			fprintf(fpp, " %9.6f", constr->a_matrix[i + j * constr->nc]);
		}
		fprintf(fpp, "\n");
	}
	for (i = 0; i < constr->nc; i++) {
		fprintf(fpp, "constraint %d, jfirst= %d jlen= %d\n", i, constr->jfirst[i], constr->jlen[i]);
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_prepare_constr(GMRFLib_constr_tp *constr, GMRFLib_graph_tp *graph, int scale_constr)
{
	/*
	 * prepare a constraint
	 */
	int nc, n;

	if (!constr || constr->is_prepared) {
		return GMRFLib_SUCCESS;
	}

	nc = constr->nc;
	n = graph->n;

	if (scale_constr && !(constr->is_scaled)) {
		/*
		 * scale the constraints so that max(|A[i,]|)=1 
		 */
		for (int k = 0; k < nc; k++) {
			int idx = idamax_(&n, constr->a_matrix + k, &nc) - 1;
			double a = 1.0 / ABS(constr->a_matrix[idx * nc + k]);
			if (!ISEQUAL(a, 1.0)) {
				dscal_(&n, &a, constr->a_matrix + k, &nc);
				constr->e_vector[k] *= a;
			}
		}
		constr->is_scaled = 1;
	}

	constr->jfirst = Calloc(2 * nc, int);
	constr->jlen = constr->jfirst + nc;

	for (int i = 0; i < nc; i++) {
		double *a = constr->a_matrix + i;
		for (int j = 0; j < n; j++) {
			if (!ISZERO(a[j * nc])) {
				constr->jfirst[i] = j;
				break;
			}
		}
		for (int j = n - 1; j >= 0; j--) {
			if (!ISZERO(a[j * nc])) {
				constr->jlen[i] = j - constr->jfirst[i] + 1;
				break;
			}
		}
	}

	constr->idxval = Calloc(nc, GMRFLib_idxval_tp *);
	for (int i = 0; i < nc; i++) {
		double *a = constr->a_matrix;
		for (int j = 0; j < n; j++) {
			int k = i + j * nc;
			if (!ISZERO(a[k])) {
				GMRFLib_idxval_add(&(constr->idxval[i]), j, a[k]);
			}
		}
	}
	GMRFLib_idxval_prepare(constr->idxval, nc, 1);

	GMRFLib_constr_add_sha(constr, graph);
	constr->is_prepared = 1;

	return GMRFLib_SUCCESS;
}

int GMRFLib_constr_add_sha(GMRFLib_constr_tp *constr, GMRFLib_graph_tp *graph)
{
	GMRFLib_SHA_TP c;
	unsigned char *md = Calloc(GMRFLib_SHA_DIGEST_LEN + 1, unsigned char);

	Memset(md, 0, GMRFLib_SHA_DIGEST_LEN + 1);
	GMRFLib_SHA_Init(&c);

	GMRFLib_SHA_DUPDATE(constr->a_matrix, graph->n * constr->nc);
	GMRFLib_SHA_DUPDATE(constr->e_vector, constr->nc);
	GMRFLib_SHA_Final(md, &c);
	md[GMRFLib_SHA_DIGEST_LEN] = '\0';
	constr->sha = md;

	return GMRFLib_SUCCESS;
}

int GMRFLib_eval_constr(double *value, double *sqr_value, double *x, GMRFLib_constr_tp *constr, GMRFLib_graph_tp *graph)
{
	/*
	 * eval a constraint `constr' at x-value `x':
	 * 
	 * if value, *value = Ax-e if sqr_value, *sqr_value = (Ax-e)'Q(Ax-e) 
	 */
	int nc = constr->nc;
	Calloc_init(2 * nc, 2);

	double *t_vector = Calloc_get(nc);
	double *res = Calloc_get(nc);
	Memcpy(t_vector, constr->e_vector, nc * sizeof(double));

	if (GMRFLib_faster_constr) {
		dgemv_special(res, x, constr);
		if (nc >= GMRFLib_SIMD_LIM) {
#pragma omp simd
			for (int i = 0; i < nc; i++) {
				t_vector[i] = res[i] - t_vector[i];
			}
		} else {
			for (int i = 0; i < nc; i++) {
				t_vector[i] = res[i] - t_vector[i];
			}
		}
	} else {
		int inc = 1;
		double alpha = 1.0;
		double beta = -1.0;
		dgemv_("N", &nc, &(graph->n), &alpha, constr->a_matrix, &nc, x, &inc, &beta, t_vector, &inc, F_ONE);
	}

	if (value) {
		Memcpy(value, t_vector, nc * sizeof(double));
	}
	if (sqr_value) {
		*sqr_value = 0.0;
	}
	Calloc_free();

	return GMRFLib_SUCCESS;
}

int GMRFLib_eval_constr0(double *value, double *sqr_value, double *x, GMRFLib_constr_tp *constr, GMRFLib_graph_tp *graph)
{
	/*
	 * eval a constraint `constr' at x-value `x', setting e=0
	 */
	int nc = constr->nc;
	Calloc_init(2 * nc, 2);

	double *t_vector = Calloc_get(nc);
	double *res = Calloc_get(nc);

	if (GMRFLib_faster_constr) {
		dgemv_special(res, x, constr);
		if (nc >= GMRFLib_SIMD_LIM) {
#pragma omp simd
			for (int i = 0; i < nc; i++) {
				t_vector[i] = res[i] - t_vector[i];
			}
		} else {
			for (int i = 0; i < nc; i++) {
				t_vector[i] = res[i] - t_vector[i];
			}
		}
	} else {
		int inc = 1;
		double alpha = 1.0;
		double beta = -1.0;
		dgemv_("N", &nc, &(graph->n), &alpha, constr->a_matrix, &nc, x, &inc, &beta, t_vector, &inc, F_ONE);
	}

	if (value) {
		Memcpy(value, t_vector, nc * sizeof(double));
	}
	if (sqr_value) {
		*sqr_value = 0.0;
	}
	Calloc_free();

	return GMRFLib_SUCCESS;
}

int GMRFLib_duplicate_constr(GMRFLib_constr_tp **new_constr, GMRFLib_constr_tp *constr, GMRFLib_graph_tp *graph)
{
	if (!constr) {
		*new_constr = NULL;
		return GMRFLib_SUCCESS;
	}

	if (constr_store_use && constr->sha) {
		void **p;
		p = map_strvp_ptr(&constr_store, (char *) constr->sha);
		if (constr_store_debug) {
			if (p) {
				printf("\t[%1d] constr_store: constr is found in store: do not duplicate.\n", omp_get_thread_num());
			} else {
				printf("\t[%1d] constr_store: constr is not found in store: duplicate.\n", omp_get_thread_num());
			}
		}
		if (p) {
			*new_constr = (GMRFLib_constr_tp *) * p;
			return GMRFLib_SUCCESS;
		}
	}

	GMRFLib_recomp_constr(new_constr, constr, NULL, NULL, NULL, graph, NULL);

	if (constr_store_use && constr->sha) {
		if (constr_store_debug) {
			printf("\t[%1d] constr_store: store constr 0x%p\n", omp_get_thread_num(), (void *) *new_constr);
		}
#pragma omp critical (Name_94faae67756d65c5760e5596c1b377f8844e3f00)
		{
			map_strvp_set(&constr_store, (char *) (*new_constr)->sha, (void *) *new_constr);
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_recomp_constr(GMRFLib_constr_tp **new_constr, GMRFLib_constr_tp *constr, double *x,
			  double *b_add, char *mask, GMRFLib_graph_tp *graph, GMRFLib_graph_tp *sub_graph)
{
	/*
	 * remap the constaints Ax=e, on graph to Ax=e on sub-graph
	 * 
	 * if stochastic constraints, then compute also ``b_add'', which is the additional terms to the 'b'-term: b^Tx.
	 * 
	 */

	GMRFLib_ENTER_ROUTINE;

	int i, k, kk, n, ns, *in_use = NULL, *cmap = NULL, nc = 0;

	if (!constr) {
		if (new_constr) {
			*new_constr = NULL;
		}
		GMRFLib_LEAVE_ROUTINE;
		return GMRFLib_SUCCESS;
	}

	if (!mask) {
		/*
		 * new_constr is equal to the old one, make a copy 
		 */
		GMRFLib_make_empty_constr(new_constr);
		(*new_constr)->nc = constr->nc;
		(*new_constr)->is_scaled = constr->is_scaled;

		(*new_constr)->a_matrix = Calloc(graph->n * constr->nc, double);
		Memcpy((*new_constr)->a_matrix, constr->a_matrix, graph->n * constr->nc * sizeof(double));

		(*new_constr)->e_vector = Calloc(constr->nc, double);
		Memcpy((*new_constr)->e_vector, constr->e_vector, constr->nc * sizeof(double));

		// this will add jfirst and jlen
		GMRFLib_prepare_constr(*new_constr, graph, 0);

		GMRFLib_LEAVE_ROUTINE;
		return GMRFLib_SUCCESS;
	}

	n = graph->n;
	ns = sub_graph->n;
	in_use = Calloc(constr->nc, int);
	cmap = Calloc(constr->nc, int);

	GMRFLib_make_empty_constr(new_constr);
	if (b_add) {
		Memset(b_add, 0, sub_graph->n * sizeof(double));
	}

	/*
	 * find those constrs that are in use, part I 
	 */
	for (k = 0; k < constr->nc; k++) {
		for (i = 0; i < ns && in_use[k] == 0; i++) {
			if (constr->a_matrix[k + constr->nc * i]) {
				in_use[k] = 1;
			}
		}
	}

	/*
	 * find those constrs that are in use, part III 
	 */
	for (k = nc = 0; k < constr->nc; k++) {
		if (in_use[k]) {
			cmap[nc++] = k;			       /* the mapping */
		}
	}
	(*new_constr)->nc = nc;

	if ((*new_constr)->nc == 0) {			       /* we do not have any constraints in sub_graph */
		Free(in_use);
		Free(cmap);
		Free(*new_constr);
		return GMRFLib_SUCCESS;
	}

	(*new_constr)->a_matrix = Calloc(nc * ns, double);
	(*new_constr)->e_vector = Calloc(nc, double);

	for (k = 0; k < nc; k++) {
		kk = cmap[k];
		(*new_constr)->e_vector[k] = constr->e_vector[kk];
		for (i = 0; i < n; i++) {
			if (mask[i]) {
				(*new_constr)->e_vector[k] -= constr->a_matrix[kk + constr->nc * i] * x[i];
			}
		}
		for (i = 0; i < ns; i++) {
			(*new_constr)->a_matrix[k + nc * i] = constr->a_matrix[kk + constr->nc * i];
		}
	}

	GMRFLib_prepare_constr(*new_constr, sub_graph, 1);

	Free(in_use);
	Free(cmap);
	GMRFLib_LEAVE_ROUTINE;

	return GMRFLib_SUCCESS;
}

int GMRFLib_info_problem(FILE *fp, GMRFLib_problem_tp *problem)
{
	if (problem) {
		return GMRFLib_fact_info_report(fp, &(problem->sub_sm_fact));
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_fact_info_report(FILE *fp, GMRFLib_sm_fact_tp *sm_fact)
{
	GMRFLib_fact_info_tp f;
	double possible_fillins;
	FILE *ffp = NULL;
	const char *sep = "-------------------------------------------------------------------------";

	if (!sm_fact)
		return GMRFLib_SUCCESS;

	ffp = (fp ? fp : stdout);
	f = sm_fact->finfo;

	fprintf(ffp, "\n\nGMRFLib report on factorisation\n%s\n", sep);

	fprintf(ffp, "Size ................................: %8d\n", f.n);
	if (f.n > 0.0) {
		fprintf(ffp, "Number of non-zeros in Q ............: %8d  (In percentage %.6f )\n", f.nnzero,
			(100.0 * f.nnzero / (SQR((double) f.n))));
		fprintf(ffp, "Number of non-zeros in L ............: %8d  (In percentage %.6f )\n",
			(f.nnzero - f.n) / 2 + f.n + f.nfillin,
			(100.0 * ((f.nnzero - f.n) / 2. + f.n + f.nfillin)) / (SQR((double) f.n) / 2. + (double) f.n / 2.));

		possible_fillins = SQR((double) f.n) / 2. + (double) f.n / 2. - f.nnzero;	/* terms in L - those in Q */

		if (possible_fillins > 0.0) {
			fprintf(ffp, "Number of fillins ...................: %8d  (In percentage %.6f )\n", f.nfillin,
				(100.0 * f.nfillin) / possible_fillins);
		}
		fprintf(ffp, "Terms in L / Minimum terms in L .....: %8.3f\n",
			((f.nnzero - f.n) / 2. + f.n + f.nfillin) / ((double) (f.nnzero - f.n) / 2. + f.n));
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_print_problem(FILE *fp, GMRFLib_problem_tp *problem)
{
	/*
	 * for verification purposes; just print the contents of 'problem' to 'fp' 
	 */

	int n, ns, nc;
	FILE *fpp = NULL;

	if (!problem) {
		return GMRFLib_SUCCESS;
	}

	n = problem->n;					       /* full graph */
	ns = problem->sub_graph->n;			       /* sub_graph */
	nc = (problem->sub_constr ? problem->sub_constr->nc : 0);
	fpp = (fp ? fp : stdout);

	fprintf(fpp, "\nContents of problem:\n");
	GMRFLib_print_darray(fpp, problem->sample, n, "sample");
	GMRFLib_print_darray(fpp, problem->mean, n, "mean");
	GMRFLib_print_darray(fpp, problem->mean_constr, n, "mean_constr");
	GMRFLib_print_darray(fpp, &(problem->sub_logdens), 1, "sub_logdens");
	GMRFLib_print_darray(fpp, problem->sub_sample, ns, "sub_sample");
	GMRFLib_print_darray(fpp, problem->sub_mean, ns, "sub_mean");
	GMRFLib_print_darray(fpp, problem->sub_mean_constr, ns, "sub_mean_constr");
	GMRFLib_print_darray(fpp, problem->constr_m, ns * nc, "constr_m");
	GMRFLib_print_darray(fpp, problem->l_aqat_m, nc * nc, "l_aqat_m");
	GMRFLib_print_darray(fpp, problem->qi_at_m, ns * nc, "qi_at_m");
	GMRFLib_print_darray(fpp, &(problem->logdet_aqat), 1, "logdet_aqat");
	GMRFLib_print_darray(fpp, &(problem->log_normc), 1, "log_normc");
	GMRFLib_print_darray(fpp, &(problem->exp_corr), 1, "exp_corr");

	GMRFLib_printf_constr(fpp, problem->sub_constr, problem->sub_graph);

	fflush(fpp);
	return GMRFLib_SUCCESS;
}

GMRFLib_problem_tp *GMRFLib_duplicate_problem(GMRFLib_problem_tp *problem, int skeleton, int copy_ptr, int copy_pardiso_ptr)
{
	/*
	 * duplicate a problem 
	 */

#define DUPLICATE(name, len, tp, skeleton_)				\
	if (1) {							\
		if (problem->name && ((len)>0) && !skeleton_){		\
			np->name = Calloc((len), tp);			\
			Memcpy(np->name, problem->name, (len)*sizeof(tp)); \
		} else {						\
			np->name = NULL;				\
		}							\
	}
#define COPY(name) np->name = problem->name

	if (!problem) {
		return NULL;
	}

	GMRFLib_problem_tp *np = Calloc(1, GMRFLib_problem_tp);

	int i;
	int n = problem->n;				       /* full graph */
	int ns = problem->sub_graph->n;			       /* sub_graph */
	int nc = (problem->sub_constr ? problem->sub_constr->nc : 0);

	DUPLICATE(sample, n, double, skeleton);
	DUPLICATE(mean, n, double, skeleton);
	DUPLICATE(mean_constr, n, double, skeleton);

	COPY(n);
	COPY(sub_logdens);
	DUPLICATE(sub_sample, ns, double, skeleton);
	DUPLICATE(sub_mean, ns, double, skeleton);
	DUPLICATE(sub_mean_constr, ns, double, skeleton);

	/*
	 * duplicate the sparse-matrix factorisation 
	 */
	DUPLICATE(sub_sm_fact.remap, ns, int, 0);
	DUPLICATE(sub_sm_fact.bchol, ns * (problem->sub_sm_fact.bandwidth + 1), double, 0);

	COPY(sub_sm_fact.bandwidth);
	COPY(sub_sm_fact.smtp);

	if (problem->sub_sm_fact.TAUCS_L && !skeleton) {
		np->sub_sm_fact.TAUCS_L = GMRFLib_L_duplicate_TAUCS(problem->sub_sm_fact.TAUCS_L, problem->sub_sm_fact.TAUCS_L->flags);
	} else {
		np->sub_sm_fact.TAUCS_L = NULL;
	}

	if (problem->sub_sm_fact.TAUCS_L_inv_diag && !skeleton) {
		DUPLICATE(sub_sm_fact.TAUCS_L_inv_diag, ns, double, skeleton);
	} else {
		np->sub_sm_fact.TAUCS_L_inv_diag = NULL;
	}
	np->sub_sm_fact.TAUCS_symb_fact = GMRFLib_sm_fact_duplicate_TAUCS(problem->sub_sm_fact.TAUCS_symb_fact);
	COPY(sub_sm_fact.finfo);

	if (problem->sub_sm_fact.PARDISO_fact) {
		GMRFLib_duplicate_pardiso_store(&(np->sub_sm_fact.PARDISO_fact), problem->sub_sm_fact.PARDISO_fact, copy_ptr, copy_pardiso_ptr);
	}

	/*
	 * then the constraint 
	 */
	if (skeleton) {
		np->sub_constr = NULL;
	} else {
		if (problem->sub_constr) {
			// this will make use of the cache
			GMRFLib_duplicate_constr(&(np->sub_constr), problem->sub_constr, problem->sub_graph);
		} else {
			COPY(sub_constr);
		}
	}

	DUPLICATE(sub_constr_value, nc, double, skeleton);
	DUPLICATE(constr_m, ns * nc, double, skeleton);
	DUPLICATE(l_aqat_m, nc * nc, double, skeleton);
	DUPLICATE(inv_aqat_m, nc * nc, double, skeleton);
	DUPLICATE(qi_at_m, ns * nc, double, skeleton);

	COPY(logdet_aat);
	COPY(logdet_aqat);
	COPY(log_normc);
	COPY(exp_corr);
	GMRFLib_graph_duplicate(&(np->sub_graph), problem->sub_graph);

	/*
	 * copy the tab 
	 */
	if (problem->tab && !skeleton) {
		GMRFLib_tabulate_Qfunc_tp *tab = Calloc(1, GMRFLib_tabulate_Qfunc_tp);
		GMRFLib_tabulate_Qfunc_arg_tp *Qfunc_arg = Calloc(1, GMRFLib_tabulate_Qfunc_arg_tp);
		GMRFLib_tabulate_Qfunc_arg_tp *tmp = (GMRFLib_tabulate_Qfunc_arg_tp *) (problem->tab->Qfunc_arg);

		tab->Qfunc = problem->tab->Qfunc;
		Qfunc_arg->n = tmp->n;
		if (tmp->log_prec_omp) {
			int tmax = GMRFLib_MAX_THREADS();
			Qfunc_arg->log_prec_omp = Calloc(tmax, double *);
			for (i = 0; i < tmax; i++) {
				Qfunc_arg->log_prec_omp[i] = tmp->log_prec_omp[i];
			}
		}
		if (tmp->values) {
			Qfunc_arg->values = Calloc(ns, map_id *);
			for (i = 0; i < ns; i++) {
				Qfunc_arg->values[i] = GMRFLib_duplicate_map_id(tmp->values[i]);
			}
		} else {
			Qfunc_arg->values = NULL;
		}
		if (tmp->Q) {
			GMRFLib_csr_duplicate(&(Qfunc_arg->Q), tmp->Q, 1);
		}

		tab->Qfunc_arg = (void *) Qfunc_arg;
		np->tab = tab;
	} else {
		np->tab = NULL;
	}

	/*
	 * copy the sub_inverse 
	 */
	if (problem->sub_inverse && !skeleton) {
		np->sub_inverse = Calloc(1, GMRFLib_Qinv_tp);
		map_id **Qinv = Calloc(n, map_id *);

		for (i = 0; i < n; i++) {
			Qinv[i] = GMRFLib_duplicate_map_id(problem->sub_inverse->Qinv[i]);
		}
		np->sub_inverse->Qinv = Qinv;
		np->sub_inverse->mapping = Calloc(n, int);
		Memcpy(np->sub_inverse->mapping, problem->sub_inverse->mapping, n * sizeof(int));
	} else {
		np->sub_inverse = NULL;
	}

#undef DUPLICATE
#undef COPY
	return np;
}

GMRFLib_store_tp *GMRFLib_duplicate_store(GMRFLib_store_tp *store, int skeleton, int copy_ptr, int copy_pardiso_ptr)
{
	/*
	 * duplicate STORE 
	 */
#define DUPLICATE(name, len, tp, skeleton_)				\
	if (1) {							\
		if (store->name && ((len)>0) && !skeleton_){		\
			new_store->name = Calloc((len), tp);		\
			Memcpy(new_store->name, store->name, (len)*sizeof(tp)); \
		} else {						\
			new_store->name = NULL;				\
		}							\
	}
#define COPY(name) new_store->name = store->name

	if (!store) {
		return NULL;
	}

	int ns = store->sub_graph->n;
	GMRFLib_store_tp *new_store = Calloc(1, GMRFLib_store_tp);

	COPY(bandwidth);
	DUPLICATE(remap, ns, int, 0);

	if (copy_ptr == GMRFLib_TRUE) {
		/*
		 * just copy ptr's; read only 
		 */
		new_store->sub_graph = store->sub_graph;
		new_store->TAUCS_symb_fact = store->TAUCS_symb_fact;
	} else {
		GMRFLib_graph_duplicate(&(new_store->sub_graph), store->sub_graph);
		new_store->TAUCS_symb_fact = GMRFLib_sm_fact_duplicate_TAUCS(store->TAUCS_symb_fact);
	}
	new_store->copy_ptr = copy_ptr;
	new_store->copy_pardiso_ptr = copy_pardiso_ptr;
	if (store->PARDISO_fact) {
		GMRFLib_duplicate_pardiso_store(&(new_store->PARDISO_fact), store->PARDISO_fact, copy_ptr, copy_pardiso_ptr);
	}

	char *tmp = Calloc(1, char);
	Free(tmp);

	COPY(store_problems);
	COPY(fixed_hyperparameters);
	COPY(decision);
	DUPLICATE(old_logdens, 1, double, skeleton);
	DUPLICATE(new_logdens, 1, double, skeleton);

	if (!skeleton) {
		new_store->problem_old2new = GMRFLib_duplicate_problem(store->problem_old2new, skeleton, copy_ptr, copy_pardiso_ptr);
		new_store->problem_new2old = GMRFLib_duplicate_problem(store->problem_new2old, skeleton, copy_ptr, copy_pardiso_ptr);
	} else {
		new_store->problem_new2old = NULL;
		new_store->problem_old2new = NULL;
	}

	if (store->diag_store) {
		new_store->diag_store = GMRFLib_duplicate_store(store->diag_store, skeleton, copy_ptr, copy_pardiso_ptr);
	}
	if (store->sub_store) {
		new_store->sub_store = GMRFLib_duplicate_store(store->sub_store, skeleton, copy_ptr, copy_pardiso_ptr);
	}
#undef DUPLICATE
#undef COPY
	return new_store;
}

double GMRFLib_Qfunc_generic(int UNUSED(thread_id), int i, int j, double *UNUSED(values), void *arg)
{
	if (j < 0) {
		return NAN;
	}

	if (i != j) {
		return -1.0;
	} else {
		GMRFLib_graph_tp *g = (GMRFLib_graph_tp *) arg;
		return g->n;
	}
}

int GMRFLib_optimize_reorder(GMRFLib_graph_tp *graph, size_t *nnz_opt, int *use_global, GMRFLib_global_node_tp *gn)
{
	if (!graph) {
		if (nnz_opt)
			*nnz_opt = 0;
		return GMRFLib_SUCCESS;
	}

	if (GMRFLib_smtp == GMRFLib_SMTP_BAND) {
		GMRFLib_reorder = GMRFLib_REORDER_DEFAULT;
		*nnz_opt = 0;
	} else if (GMRFLib_smtp == GMRFLib_SMTP_PARDISO) {
		GMRFLib_reorder = GMRFLib_REORDER_PARDISO;
		*nnz_opt = 0;
	} else {
		static int debug = 0;
		size_t *nnzs = NULL, nnz_best;
		int k, n = -1, nk, r, i, ne = 0, use_global_nodes;
		GMRFLib_reorder_tp rs[] = { GMRFLib_REORDER_METIS }; // , GMRFLib_REORDER_AMDBARC };
		taucs_ccs_matrix *Q = NULL;
		char *fixed = NULL;
		double *cputime = NULL;

		n = graph->n;

		/*
		 * build the Q-matrix; just symbolically, so I set Q_ij = 1. this matrix is common for all reordering
		 * schemes. 
		 */
		int ic, kk, j, nnz;

		//GMRFLib_printf_graph(stdout, graph);
		
		nnz = n + GMRFLib_isum(n, graph->nnbs);
		Q = taucs_ccs_create(n, n, nnz, TAUCS_DOUBLE);
		GMRFLib_ASSERT(Q, GMRFLib_EMEMORY);
		Q->flags = (TAUCS_DOUBLE | TAUCS_SYMMETRIC | TAUCS_TRIANGULAR | TAUCS_LOWER | TAUCS_PATTERN);
		Q->colptr[0] = 0;

		for (i = 0, ic = 0; i < n; i++) {
			Q->rowind[ic] = i;
			Q->values.d[ic] = 1.0;
			ic++;
			ne = 1;

			for (kk = 0; kk < graph->nnbs[i]; kk++) {
				j = graph->nbs[i][kk];
				if (j > i) {
					break;
				}
				Q->rowind[ic] = j;
				Q->values.d[ic] = 1.0;
				ic++;
				ne++;
			}
			Q->colptr[i + 1] = Q->colptr[i] + ne;
		}

		/*
		 * do this for all reorderings, first with global node detection and then without global nodes detection 
		 */
		nk = 2 * (int) (sizeof(rs) / sizeof(int));     /* yes, twice... */
		nnzs = Calloc(nk, size_t);
		cputime = Calloc(nk, double);

		for (k = 0; k < nk; k++) {
			int *iperm = NULL, *perm = NULL, ii, kkk;
			supernodal_factor_matrix *TAUCS_symb_fact = NULL;
			taucs_ccs_matrix *L = NULL;

			GMRFLib_global_node_tp lgn;
			if (gn) {
				Memcpy((void *) &lgn, (void *) gn, sizeof(GMRFLib_global_node_tp));
			} else {
				Memcpy((void *) &lgn, (void *) &GMRFLib_global_node, sizeof(GMRFLib_global_node_tp));
			}

			/*
			 * first half is with, second half is without, the detection of global nodes 
			 */
			kkk = k % (nk / 2);
			use_global_nodes = (k < nk / 2);

			if (!use_global_nodes) {
				/*
				 * currently disable 
				 */
				lgn.factor = 2.0;
				lgn.degree = INT_MAX;
			}

			/*
			 * if we run with global nodes, it can be that the settings are so that they are not in effect. we check for this here, and if so, we do not
			 * need to try it, as it will be tried again in the second half where the global_nodes are disabled.
			 */
			if (!use_global_nodes || !(use_global_nodes && (lgn.factor > 1.0) && (lgn.degree > graph->n - 1))) {

				cputime[k] = GMRFLib_cpu();
				GMRFLib_compute_reordering_TAUCS(&iperm, graph, rs[kkk], &lgn);

				perm = Calloc(n, int);
				for (ii = 0; ii < n; ii++) {
					perm[iperm[ii]] = ii;
				}

				L = taucs_ccs_permute_symmetrically(Q, perm, iperm);	/* permute the matrix */
				TAUCS_symb_fact = (supernodal_factor_matrix *) taucs_ccs_factor_llt_symbolic(L);
				nnzs[k] = GMRFLib_sm_fact_nnz_TAUCS(TAUCS_symb_fact);
				Free(perm);
				Free(iperm);
				taucs_ccs_free(L);
				taucs_supernodal_factor_free(TAUCS_symb_fact);

				cputime[k] = GMRFLib_cpu() - cputime[k];

				if (debug) {
#pragma omp critical (Name_4dc800d9b856792e63baa1e9a01d82865c857322)
					{
						printf("%s: reorder=[%s] \tnnz=%zu \tUseGlobalNodes=%1d cpu=%.4f\n",
						       __GMRFLib_FuncName, GMRFLib_reorder_name(rs[kkk]), nnzs[k], use_global_nodes, cputime[k]);
					}
				}

			} else {
				nnzs[k] = UINT_MAX;
			}
		}

		/*
		 * find the best one 
		 */
		r = 0;
		nnz_best = nnzs[r];
		for (k = 1; k < nk; k++) {
			if (nnzs[k] < nnz_best) {
				r = k;
				nnz_best = nnzs[r];
			}
		}

		/*
		 * find out which one this is 
		 */
		use_global_nodes = (r < (nk / 2));
		r = r % (nk / 2);

		if (!use_global_nodes) {
			GMRFLib_global_node_tp g = { 2.0, INT_MAX };
			GMRFLib_global_node = g;
			if (gn) {
				/*
				 * if so, copy the result over there as well
				 */
				Memcpy((void *) gn, (void *) &g, sizeof(GMRFLib_global_node));
			}
		} else {
			if (gn) {
				Memcpy((void *) &GMRFLib_global_node, (void *) gn, sizeof(GMRFLib_global_node));
			}
		}
		GMRFLib_reorder = rs[r];

		if (debug) {
			printf("%s: best reordering=[%s] UseGlobalNodes=%1d\n", __GMRFLib_FuncName,
			       GMRFLib_reorder_name(GMRFLib_reorder), use_global_nodes);
		}
		taucs_ccs_free(Q);
		if (nnz_opt)
			*nnz_opt = nnz_best;
		if (use_global)
			*use_global = use_global_nodes;

		Free(nnzs);
		Free(fixed);
	}
	return GMRFLib_SUCCESS;
}
