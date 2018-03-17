/* smtp-pardiso.c
 * 
 * Copyright (C) 2018 Havard Rue
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

/*!
  \file smtp-pardiso.c
  \brief The interface to the PARDISO library
*/

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#if !defined(__FreeBSD__)
#include <malloc.h>
#endif

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

#define WARNING(_msg) fprintf(stderr, "\n\n%s:%1d: %s\n\n", __FILE__, __LINE__, _msg)

static struct 
{
	int verbose;
	int s_verbose;
	int csr_check;
	int mnum;
	int *busy;
	GMRFLib_pardiso_store_tp **static_pstores;
}
	S = { 0, 0, 0, 0, NULL, NULL};
	
//#define PSTORES_NUM() ISQR(GMRFLib_MAX_THREADS)
#define PSTORES_NUM() 1024

int GMRFLib_free_csr(GMRFLib_csr_tp ** csr)
{
	if (*csr) {
		Free((*csr)->ia);
		Free((*csr)->ja);
		Free((*csr)->a);
		Free(*csr);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_duplicate_csr(GMRFLib_csr_tp ** csr_to, GMRFLib_csr_tp * csr_from)
{
	if (csr_from == NULL) {
		*csr_to = NULL;
		return GMRFLib_SUCCESS;
	}

	*csr_to = Calloc(1, GMRFLib_csr_tp);
	(*csr_to)->base = csr_from->base;
	(*csr_to)->n = csr_from->n;
	(*csr_to)->na = csr_from->na;

	(*csr_to)->ia = Calloc(csr_from->n + 1, int);
	memcpy((void *) ((*csr_to)->ia), (void *) (csr_from->ia), (size_t) (csr_from->n + 1) * sizeof(int));

	(*csr_to)->ja = Calloc(csr_from->na, int);
	memcpy((void *) ((*csr_to)->ja), (void *) (csr_from->ja), (size_t) (csr_from->na) * sizeof(int));

	(*csr_to)->a = Calloc(csr_from->na, double);
	memcpy((void *) ((*csr_to)->a), (void *) (csr_from->a), (size_t) (csr_from->na) * sizeof(double));

	return GMRFLib_SUCCESS;
}

int GMRFLib_csr_base(int base, GMRFLib_csr_tp * M)
{
	if (M->base != base) {
		GMRFLib_csr_convert(M);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_csr_convert(GMRFLib_csr_tp * M)
{
	int i;
	if (M->base == 1) {
		M->base = 0;
		for (i = 0; i < M->n + 1; i++) {
			M->ia[i] -= 1;
		}
		for (i = 0; i < M->na; i++) {
			M->ja[i] -= 1;
		}
	} else if (M->base == 0) {
		M->base = 1;
		for (i = 0; i < M->n + 1; i++) {
			M->ia[i] += 1;
		}
		for (i = 0; i < M->na; i++) {
			M->ja[i] += 1;
		}
	} else {
		assert(0 == 1);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_Q2csr_check(GMRFLib_csr_tp * M)
{
	int mtype = GMRFLib_PARDISO_MTYPE, error = 0;
	GMRFLib_csr_tp *Q;

	GMRFLib_duplicate_csr(&Q, M);
	GMRFLib_csr_base(1, Q);
	pardiso_chkmatrix(&mtype, &(Q->n), Q->a, Q->ia, Q->ja, &error);
	GMRFLib_free_csr(&Q);
	if (error != 0) {
		GMRFLib_ERROR(GMRFLib_ESNH);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_Q2csr(GMRFLib_csr_tp ** csr, GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg)
{
	// create a upper triangular csr matrix from Q
#define M (*csr)
	int i, j, jj, k, n, na, nnz;

	M = Calloc(1, GMRFLib_csr_tp);
	M->base = 0;
	n = graph->n;
	GMRFLib_nQelm(&nnz, graph);			       // symmetric
	na = (nnz - n) / 2 + n;				       // only upper triangular. yes, integer division
	M->na = na;
	M->n = n;
	M->a = Calloc(na, double);
	M->ja = Calloc(na, int);
	M->ia = Calloc(n + 1, int);

	M->ia[0] = 0;
	for (i = 0; i < n; i++) {
		k = 1;
		for (jj = 0; jj < graph->nnbs[i]; jj++) {
			if (graph->nbs[i][jj] > i)
				k++;
		}
		M->ia[i + 1] = M->ia[i] + k;
	}
	assert(M->ia[n] == na);

	for (i = k = 0; i < n; i++) {
		M->ja[k] = i;
		M->a[k] = Qfunc(i, i, Qfunc_arg);
		k++;
		for (jj = 0; jj < graph->nnbs[i]; jj++) {
			j = graph->nbs[i][jj];
			if (j > i) {
				M->ja[k] = j;
				M->a[k] = Qfunc(i, j, Qfunc_arg);
				k++;
			}
		}
	}
#undef M

	return GMRFLib_SUCCESS;
}

int GMRFLib_print_csr(FILE * fp, GMRFLib_csr_tp * csr)
{
	int i, j, k, jj, nnb;
	double value;

	fprintf(fp, "Q->base = %1d, Q->n = %1d, Q->na = %1d\n", csr->base, csr->n, csr->na);
	for (i = k = 0; i < csr->n; i++) {
		nnb = csr->ia[i + 1] - csr->ia[i];
		for (jj = 0; jj < nnb; jj++) {
			j = csr->ja[k];
			value = csr->a[k];
			fprintf(fp, "%sQ[ %1d , %1d ] = %.12f\n", (jj > 0 ? "\t" : ""), i + csr->base, j, value);
			k++;
		}
	}
	fflush(fp);

	return GMRFLib_SUCCESS;
}

int GMRFLib_csr2Q(GMRFLib_tabulate_Qfunc_tp ** Qtab, GMRFLib_graph_tp ** graph, GMRFLib_csr_tp * csr)
{
	int i, j, k, jj, nnb;

	int *iarr = Calloc(csr->na, int);
	int *jarr = Calloc(csr->na, int);
	double *arr = Calloc(csr->na, double);

	for (i = k = 0; i < csr->n; i++) {
		nnb = csr->ia[i + 1] - csr->ia[i];
		for (jj = 0; jj < nnb; jj++) {
			j = csr->ja[k];
			iarr[k] = i - csr->base;
			jarr[k] = j - csr->base;
			arr[k] = csr->a[k];
			k++;
		}
	}
	assert(k == csr->na);
	GMRFLib_tabulate_Qfunc_from_list(Qtab, graph, csr->na, iarr, jarr, arr, csr->n, NULL, NULL, NULL);

	Free(iarr);
	Free(jarr);
	Free(arr);

	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_init(GMRFLib_pardiso_store_tp ** store)
{
	GMRFLib_ENTER_ROUTINE;

	int error = 0;
	GMRFLib_pardiso_store_tp *s = Calloc(1, GMRFLib_pardiso_store_tp);
	
	if (S.s_verbose)
		PP("_pardiso_init()", s);

	s->maxfct = ISQR(GMRFLib_MAX_THREADS);
	s->pstore = Calloc(s->maxfct, GMRFLib_pardiso_store_pr_thread_tp *); 
	for(int i = 0; i < s->maxfct; i++){
		s->pstore[i] = Calloc(1, GMRFLib_pardiso_store_pr_thread_tp);
	}
	s->mtype = GMRFLib_PARDISO_MTYPE;
	s->msglvl = 0;
	s->solver = 0;
	s->iparm_default = Calloc(GMRFLib_PARDISO_PLEN, int);
	s->dparm_default = Calloc(GMRFLib_PARDISO_PLEN, double);
	s->iparm_default[0] = 0;			       /* use default values */
	s->iparm_default[2] = GMRFLib_PARDISO_NUM_PROC_DEFAULT;
	P1(s->iparm_default[2]);
	s->iparm_default[4] = 0;			       /* user internal reordering */
	
	pardisoinit(s->pt, &(s->mtype), &(s->solver), s->iparm_default, s->dparm_default, &error);
	if (error != 0) {
		GMRFLib_ERROR(GMRFLib_EPARDISO_INTERNAL_ERROR);
	}

	s->done_with_init = GMRFLib_TRUE;
	*store = s;

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_setparam(GMRFLib_pardiso_flag_tp flag, GMRFLib_pardiso_store_tp * store)
{
	static int check_install = GMRFLib_TRUE;
	if (check_install == GMRFLib_TRUE) {
#pragma omp critical
		if (check_install == GMRFLib_TRUE) {
			GMRFLib_pardiso_check_install(1, 0);
			check_install = GMRFLib_FALSE;
		}
	}

	assert(store->done_with_init == GMRFLib_TRUE);
	memcpy((void *) (store->pstore[S.mnum]->iparm), (void *) (store->iparm_default), GMRFLib_PARDISO_PLEN * sizeof(int));
	memcpy((void *) (store->pstore[S.mnum]->dparm), (void *) (store->dparm_default), GMRFLib_PARDISO_PLEN * sizeof(double));

	store->pstore[S.mnum]->nrhs = 0;
	store->pstore[S.mnum]->err_code = 0;

	switch (flag) {
	case GMRFLib_PARDISO_FLAG_REORDER:
		store->pstore[S.mnum]->phase = 11;			       // analysis
		store->pstore[S.mnum]->iparm[4] = 0;			       /* 0 = compute the permutation */
		break;

	case GMRFLib_PARDISO_FLAG_SYMFACT:
		store->pstore[S.mnum]->phase = 11;			       // analysis
		break;

	case GMRFLib_PARDISO_FLAG_CHOL:
		store->pstore[S.mnum]->phase = 12;			       // Analysis, numerical factorization
		store->pstore[S.mnum]->iparm[32] = 1;			       /* determinant */
		break;

	case GMRFLib_PARDISO_FLAG_QINV:
		store->pstore[S.mnum]->phase = -22;
		store->pstore[S.mnum]->iparm[35] = 1;			       /* do not overwrite internal factor L with selected inversion */
		store->pstore[S.mnum]->iparm[36] = 0;			       /* return upper triangular Qinv */
		break;

	case GMRFLib_PARDISO_FLAG_SOLVE_L:
		store->pstore[S.mnum]->phase = 33;			       // solve
		store->pstore[S.mnum]->nrhs = 1;
		store->pstore[S.mnum]->iparm[7] = 0;			       /* Max numbers of iterative refinement steps. */
		store->pstore[S.mnum]->iparm[25] = 1;
		break;

	case GMRFLib_PARDISO_FLAG_SOLVE_LT:
		store->pstore[S.mnum]->phase = 33;			       // solve
		store->pstore[S.mnum]->nrhs = 1;
		store->pstore[S.mnum]->iparm[7] = 0;			       /* Max numbers of iterative refinement steps. */
		store->pstore[S.mnum]->iparm[25] = 2;
		break;

	case GMRFLib_PARDISO_FLAG_SOLVE_LLT:
		store->pstore[S.mnum]->phase = 33;			       // solve
		store->pstore[S.mnum]->nrhs = 1;
		store->pstore[S.mnum]->iparm[7] = 0;			       /* Max numbers of iterative refinement steps. */
		store->pstore[S.mnum]->iparm[25] = 0;
		break;

	default:
		GMRFLib_ERROR(GMRFLib_ESNH);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_check_install(int quiet, int no_err)
{
	int *iparm = Calloc(GMRFLib_PARDISO_PLEN, int);
	int mtype = GMRFLib_PARDISO_MTYPE, err_code = 0, solver = 0;
	double *dparm = Calloc(GMRFLib_PARDISO_PLEN, double);
	void **pt = Calloc(GMRFLib_PARDISO_PLEN, void *);

	STDOUT_TO_DEV_NULL_START(quiet);
	pardisoinit(pt, &mtype, &solver, iparm, dparm, &err_code);
	STDOUT_TO_DEV_NULL_END;

	Free(pt);
	Free(dparm);
	Free(iparm);

	if (!no_err) {
		if (err_code != 0) {
			if (err_code == -10) {
				GMRFLib_ERROR(GMRFLib_EPARDISO_LICENSE_NOTFOUND);
			} else if (err_code == -11) {
				GMRFLib_ERROR(GMRFLib_EPARDISO_LICENSE_EXPIRED);
			} else if (err_code == -12) {
				GMRFLib_ERROR(GMRFLib_EPARDISO_LICENSE_ERR_USERNAME);
			} else {
				GMRFLib_ERROR(GMRFLib_ESNH);
			}
			assert(0 == 1);
		}
	}
	return (err_code == 0 ? GMRFLib_SUCCESS : !GMRFLib_SUCCESS);
}

double GMRFLib_pardiso_Qfunc_default(int i, int j, void *arg)
{
	GMRFLib_graph_tp *g = (GMRFLib_graph_tp *) arg;
	return (i == j ? g->n + 2.0 * g->nnbs[i] : -1.0);
}

int GMRFLib_pardiso_reorder(GMRFLib_pardiso_store_tp * store, GMRFLib_graph_tp * graph)
{
	assert(store != NULL);
	assert(store->done_with_init == GMRFLib_TRUE);

	if (store->done_with_reorder == GMRFLib_TRUE) {
		return GMRFLib_SUCCESS;
	}

	GMRFLib_ENTER_ROUTINE;

	int i, mnum1 = 1;
	GMRFLib_csr_tp *Q = NULL;

	GMRFLib_copy_graph(&(store->graph), graph);
	GMRFLib_pardiso_setparam(GMRFLib_PARDISO_FLAG_REORDER, store);
	GMRFLib_Q2csr(&Q, store->graph, GMRFLib_pardiso_Qfunc_default, (void *) store->graph);
	GMRFLib_csr_base(1, Q);
	if (S.csr_check) {
		GMRFLib_Q2csr_check(Q);
	}
	if (S.verbose) {
		GMRFLib_print_csr(stdout, Q);
	}

	pardiso(store->pt, &(store->maxfct), &mnum1, &(store->mtype),
		&(store->pstore[S.mnum]->phase),
		&(Q->n), Q->a, Q->ia, Q->ja, NULL, 
		&(store->pstore[S.mnum]->nrhs), store->pstore[S.mnum]->iparm,
		&(store->msglvl), &(store->pstore[S.mnum]->dummy), &(store->pstore[S.mnum]->dummy),
		&(store->pstore[S.mnum]->err_code), store->pstore[S.mnum]->dparm);

	if (store->pstore[S.mnum]->err_code) {
		GMRFLib_ERROR(GMRFLib_EPARDISO_INTERNAL_ERROR);
	}

	store->done_with_reorder = GMRFLib_TRUE;
	store->pstore[S.mnum]->L_nnz = store->pstore[S.mnum]->iparm[17] - 1;
	GMRFLib_free_csr(&Q);

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_symfact(GMRFLib_pardiso_store_tp * store)
{
	assert(store->done_with_init == GMRFLib_TRUE);
	assert(store->done_with_reorder == GMRFLib_TRUE);

	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_build(GMRFLib_pardiso_store_tp * store, 
			  GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg)
{
	GMRFLib_ENTER_ROUTINE;

	assert(store->done_with_init == GMRFLib_TRUE);
	assert(store->done_with_reorder == GMRFLib_TRUE);

	if (store->pstore[S.mnum]->done_with_build == GMRFLib_TRUE && store->pstore[S.mnum]->Q) {
		GMRFLib_free_csr(&(store->pstore[S.mnum]->Q));
	}
	GMRFLib_Q2csr(&(store->pstore[S.mnum]->Q), graph, Qfunc, (void *) Qfunc_arg);
	GMRFLib_csr_base(1, store->pstore[S.mnum]->Q);
	if (S.csr_check) {
		GMRFLib_Q2csr_check(store->pstore[S.mnum]->Q);
	}
	if (S.verbose)
		GMRFLib_print_csr(stdout, store->pstore[S.mnum]->Q);

	store->pstore[S.mnum]->done_with_build = GMRFLib_TRUE;
	GMRFLib_LEAVE_ROUTINE;

	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_chol(GMRFLib_pardiso_store_tp * store)
{
	GMRFLib_ENTER_ROUTINE;

	assert(store->done_with_init == GMRFLib_TRUE);
	assert(store->done_with_reorder == GMRFLib_TRUE);
	assert(store->pstore[S.mnum]->done_with_build == GMRFLib_TRUE);
	assert(store->pstore[S.mnum]->Q != NULL);

	GMRFLib_pardiso_setparam(GMRFLib_PARDISO_FLAG_CHOL, store);
	
	int mnum1 = S.mnum+1;
	pardiso(store->pt, &(store->maxfct), &mnum1, &(store->mtype), &(store->pstore[S.mnum]->phase),
		&(store->pstore[S.mnum]->Q->n),
		store->pstore[S.mnum]->Q->a, store->pstore[S.mnum]->Q->ia, store->pstore[S.mnum]->Q->ja,
		&(store->pstore[S.mnum]->idummy), &(store->pstore[S.mnum]->nrhs),
		store->pstore[S.mnum]->iparm, &(store->msglvl), NULL, NULL, &(store->pstore[S.mnum]->err_code),
		store->pstore[S.mnum]->dparm);

	if (store->pstore[S.mnum]->err_code != 0) {
		GMRFLib_ERROR(GMRFLib_EPOSDEF);
	}

	if (store->pstore[S.mnum]->iparm[22] > 0) {
		printf("\nERROR not pos def matrix, #neg.eigen %d", store->pstore[S.mnum]->iparm[22]);
		GMRFLib_ERROR(GMRFLib_EPOSDEF);
	}

	store->pstore[S.mnum]->log_det_Q = store->pstore[S.mnum]->dparm[32];
	store->pstore[S.mnum]->done_with_chol = GMRFLib_TRUE;
	GMRFLib_LEAVE_ROUTINE;

	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_solve_core(GMRFLib_pardiso_store_tp * store, 
			       GMRFLib_pardiso_flag_tp flag, double *x, double *b)
{
	assert(store->done_with_init == GMRFLib_TRUE);
	assert(store->done_with_reorder == GMRFLib_TRUE);
	assert(store->pstore[S.mnum]->done_with_build == GMRFLib_TRUE);
	assert(store->pstore[S.mnum]->done_with_chol == GMRFLib_TRUE);

	// this is so that the RHS can be overwritten
	double *xx = Calloc(store->pstore[S.mnum]->Q->n, double);

	GMRFLib_pardiso_setparam(flag, store);
	int mnum1 = S.mnum+1;
	pardiso(store->pt, &(store->maxfct), &mnum1, &(store->mtype), &(store->pstore[S.mnum]->phase),
		&(store->pstore[S.mnum]->Q->n),
		store->pstore[S.mnum]->Q->a, store->pstore[S.mnum]->Q->ia, store->pstore[S.mnum]->Q->ja,
		NULL,
		&(store->pstore[S.mnum]->nrhs), store->pstore[S.mnum]->iparm, &(store->msglvl), b, xx,
		&(store->pstore[S.mnum]->err_code), store->pstore[S.mnum]->dparm);

	if (store->pstore[S.mnum]->err_code != 0) {
		GMRFLib_ERROR(GMRFLib_EPARDISO_INTERNAL_ERROR);
	}

	memcpy((void *) x, (void *) xx, store->pstore[S.mnum]->Q->n * sizeof(double));
	Free(xx);

	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_solve_L(GMRFLib_pardiso_store_tp * store, double *x, double *b)
{
	GMRFLib_ENTER_ROUTINE;
	int res = GMRFLib_pardiso_solve_core(store, GMRFLib_PARDISO_FLAG_SOLVE_L, x, b);
	GMRFLib_LEAVE_ROUTINE;

	return res;
}

int GMRFLib_pardiso_solve_LT(GMRFLib_pardiso_store_tp * store, double *x, double *b)
{
	GMRFLib_ENTER_ROUTINE;
	int res = GMRFLib_pardiso_solve_core(store, GMRFLib_PARDISO_FLAG_SOLVE_LT, x, b);
	GMRFLib_LEAVE_ROUTINE;

	return res;
}

int GMRFLib_pardiso_solve_LLT(GMRFLib_pardiso_store_tp * store, double *x, double *b)
{
	GMRFLib_ENTER_ROUTINE;
	int res = GMRFLib_pardiso_solve_core(store, GMRFLib_PARDISO_FLAG_SOLVE_LLT, x, b);
	GMRFLib_LEAVE_ROUTINE;

	return res;
}

double GMRFLib_pardiso_logdet(GMRFLib_pardiso_store_tp * store)
{
	return (store->pstore[S.mnum]->log_det_Q);
}

int GMRFLib_pardiso_bitmap(void)
{
	/*
	 * Not available for PARDISO, USE TAUCS
	 */
	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_Qinv_INLA(GMRFLib_problem_tp * problem) 
{
	if (problem == NULL) {
		return GMRFLib_SUCCESS;
	}

	GMRFLib_ENTER_ROUTINE;
	GMRFLib_pardiso_Qinv(problem->sub_sm_fact.PARDISO_fact);

	GMRFLib_csr_tp *Qi = problem->sub_sm_fact.PARDISO_fact->pstore[S.mnum]->Qinv;
	int n = Qi->n, i, j, jj, k, kk;
	double value;
	map_id **Qinv = Calloc(n, map_id *);

	assert(Qi->base == 0);

	for (i = k = 0; i < n; i++) {
		int nnb;

		nnb = Qi->ia[i + 1] - Qi->ia[i];
		Qinv[i] = Calloc(1, map_id);
		map_id_init_hint(Qinv[i], nnb);
		for (jj = 0; jj < nnb; jj++) {
			j = Qi->ja[k];
			map_id_set(Qinv[i], j, Qi->a[k]);
			k++;
		}
	}

	if (problem->sub_constr && problem->sub_constr->nc > 0) {
#pragma omp parallel for private(i, k, j, kk, value)
		for (i = 0; i < n; i++) {
			for (k = -1; (k = (int) map_id_next(Qinv[i], k)) != -1;) {
				j = Qinv[i]->contents[k].key;
				map_id_get(Qinv[i], j, &value);
				for (kk = 0; kk < problem->sub_constr->nc; kk++) {
					value -= problem->constr_m[i + kk * n] * problem->qi_at_m[j + kk * n];
				}
				map_id_set(Qinv[i], j, value);
			}
		}
	}

	GMRFLib_Qinv_tp *subQinv = Calloc(1, GMRFLib_Qinv_tp);

	subQinv->Qinv = Qinv;
	subQinv->mapping = Calloc(1, map_ii);
	map_ii_init_hint(subQinv->mapping, n);
	for (i = 0; i < n; i++) {
		//map_ii_set(subQinv->mapping, i, i);
		//printf("Mapping %d %d\n", i, problem->sub_graph->mothergraph_idx[i]);
		map_ii_set(subQinv->mapping, problem->sub_graph->mothergraph_idx[i], problem->sub_graph->mothergraph_idx[i]);
	}

	problem->sub_inverse = subQinv;
	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_Qinv(GMRFLib_pardiso_store_tp * store)
{
	GMRFLib_ENTER_ROUTINE;

	assert(store->done_with_init == GMRFLib_TRUE);
	assert(store->done_with_reorder == GMRFLib_TRUE);
	assert(store->pstore[S.mnum]->done_with_build == GMRFLib_TRUE);
	assert(store->pstore[S.mnum]->done_with_chol == GMRFLib_TRUE);

	if (store->pstore[S.mnum]->Qinv) {
		GMRFLib_free_csr(&(store->pstore[S.mnum]->Qinv));
	}
	
	GMRFLib_duplicate_csr(&(store->pstore[S.mnum]->Qinv), store->pstore[S.mnum]->Q);
	GMRFLib_csr_base(1, store->pstore[S.mnum]->Qinv);

	GMRFLib_pardiso_setparam(GMRFLib_PARDISO_FLAG_QINV, store);
	int mnum1 = S.mnum+1;
	pardiso(store->pt, &(store->maxfct), &mnum1, &(store->mtype), &(store->pstore[S.mnum]->phase),
		&(store->pstore[S.mnum]->Qinv->n),
		store->pstore[S.mnum]->Qinv->a, store->pstore[S.mnum]->Qinv->ia, store->pstore[S.mnum]->Qinv->ja,
		NULL, &(store->pstore[S.mnum]->nrhs), store->pstore[S.mnum]->iparm, &(store->msglvl), NULL, NULL,
		&(store->pstore[S.mnum]->err_code), store->pstore[S.mnum]->dparm);

	if (store->pstore[S.mnum]->err_code != 0) {
		GMRFLib_ERROR(GMRFLib_EPARDISO_INTERNAL_ERROR);
	}

	GMRFLib_csr_base(0, store->pstore[S.mnum]->Qinv);
	if (S.verbose)
		GMRFLib_print_csr(stdout, store->pstore[S.mnum]->Qinv);

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_free(GMRFLib_pardiso_store_tp ** store)
{
	if (store == NULL || *store == NULL) {
		return GMRFLib_SUCCESS;
	}
	GMRFLib_ENTER_ROUTINE;

	if (S.s_verbose)
		PP("free: old=", *store);

#pragma omp critical
	{
		int found = 0, i;
		if (S.static_pstores != NULL) {
			for(i = 0; i < PSTORES_NUM(); i++){
				if (S.busy[i] && S.s_verbose){
					printf("in store: i=%1d s=%lx\n", i, (unsigned long int) ((void *)S.static_pstores[i]));
				}
			}
			for(i = 0; i < PSTORES_NUM() && !found; i++){
				if (S.busy[i] && (S.static_pstores[i] == *store)) {

					if (S.s_verbose) {
						P(S.busy[i]);
						PP("S.static_pstores[i]", S.static_pstores[i]);
						PP("*store", *store);
					}

					found = 1;
					S.busy[i] = 0;
					if (S.s_verbose)
						printf("==> free store[%1d]\n", i);
				}
			}
		}

		if (!found) {
			if (S.s_verbose)
				printf("==> free manually as not found\n");

			int mnums1 = 0;
			if ((*store)->pstore) {
				(*store)->pstore[S.mnum]->phase = -1;
				int mnum1 = S.mnum+1;
				pardiso((*store)->pt, &((*store)->maxfct), &mnum1, &((*store)->mtype),
					&((*store)->pstore[S.mnum]->phase),
					&((*store)->pstore[S.mnum]->idummy),
					&((*store)->pstore[S.mnum]->dummy), &((*store)->pstore[S.mnum]->idummy),
					&((*store)->pstore[S.mnum]->idummy),
					&((*store)->pstore[S.mnum]->idummy),
					&((*store)->pstore[S.mnum]->nrhs), (*store)->pstore[S.mnum]->iparm, &((*store)->msglvl),
					NULL, NULL, &((*store)->pstore[S.mnum]->err_code), (*store)->pstore[S.mnum]->dparm);

				for(i = 0; i < (*store)->maxfct; i++){
					GMRFLib_free_csr(&((*store)->pstore[i]->Q));
					GMRFLib_free_csr(&((*store)->pstore[i]->Qinv));
				}
				Free((*store)->pstore);
			}
	
			GMRFLib_free_graph((*store)->graph);
			Free((*store)->iparm_default);
			Free((*store)->dparm_default);
			Free(*store);
		}
	}
	
	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_duplicate_pardiso_store(GMRFLib_pardiso_store_tp ** new, GMRFLib_pardiso_store_tp * old)
{
	if (old == NULL) {
		*new = NULL;
		return GMRFLib_SUCCESS;
	}
	
	if (0) {
		FIXME("-->duplicate by creating a new one each time");
		GMRFLib_pardiso_init(new);
		GMRFLib_pardiso_reorder(*new, old->graph);
		return GMRFLib_SUCCESS;
	}
	GMRFLib_ENTER_ROUTINE;

#pragma omp critical
	{
	
		int i;
		if (S.static_pstores == NULL) {
			S.static_pstores = Calloc(PSTORES_NUM(), GMRFLib_pardiso_store_tp *);
			S.busy = Calloc(PSTORES_NUM(), int);
			if (S.s_verbose)
				printf("==> init static_pstores\n");
		}

		int found;
		for(i = 0, found = 0; i < PSTORES_NUM() && !found; i++){
			if (!S.busy[i]) {
				if (S.static_pstores[i]) {
					*new = S.static_pstores[i];
					if (S.s_verbose)
						printf("==> reuse store[%1d]\n", i);
				} else {
					GMRFLib_pardiso_init(&(S.static_pstores[i]));
					GMRFLib_pardiso_reorder(S.static_pstores[i], old->graph);
					*new = S.static_pstores[i];
					if (S.s_verbose)
						printf("==> new store[%1d]\n", i);
				}
				found = 1;
				S.busy[i] = 1;
			}
		}
		assert(found == 1);
		if (S.s_verbose) {
			printf("duplicate: new=%x old=%x i=%1d\n", *((void **)new), ((void *)old), i);
		}
	}


	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

//
int my_pardiso_test(void)
{
	int err = 0;
	
	if (0) {
		err = GMRFLib_pardiso_check_install(1, 0);
		if (err == GMRFLib_SUCCESS) {
			printf("PARDISO OK\n");
		} else {
			printf("PARDISO FAIL\n");
		}
	}

	GMRFLib_tabulate_Qfunc_tp *Qtab;
	GMRFLib_graph_tp *g;

	GMRFLib_tabulate_Qfunc_from_file(&Qtab, &g, "Qdense.txt", -1, NULL, NULL, NULL);

	GMRFLib_csr_tp *csr, *csr2;

	GMRFLib_Q2csr(&csr, g, Qtab->Qfunc, Qtab->Qfunc_arg);
	//GMRFLib_print_csr(stdout, csr);
	P(csr->n);
	
	GMRFLib_duplicate_csr(&csr2, csr);
	// GMRFLib_print_csr(stdout, csr2);

	GMRFLib_csr2Q(&Qtab, &g, csr2);
	// GMRFLib_print_Qfunc(stdout, g, Qtab->Qfunc, Qtab->Qfunc_arg);

	int *perm = NULL;
	int i, j;

	//GMRFLib_print_graph(stdout, g);
	GMRFLib_pardiso_store_tp *store = NULL;
	GMRFLib_pardiso_init(&store);
	GMRFLib_pardiso_reorder(store, g);
	GMRFLib_pardiso_build(store, g, Qtab->Qfunc, Qtab->Qfunc_arg);
	GMRFLib_pardiso_chol(store);
	GMRFLib_pardiso_Qinv(store);

	P(GMRFLib_openmp->max_threads_outer);
	P(GMRFLib_openmp->max_threads_inner);
	
	//S.s_verbose=1;
#pragma omp parallel for private(j) num_threads(2)
	for (j = 0; j < 100; j++) {
		printf("this is j= %d from thread %d\n", j, omp_get_thread_num());

		GMRFLib_pardiso_store_tp *store2 = NULL;
		if (0) {
			GMRFLib_pardiso_init(&store2);
			GMRFLib_pardiso_reorder(store2, g);
		} else {
			GMRFLib_duplicate_pardiso_store(&store2, store);
		}
		GMRFLib_pardiso_build(store2, g, Qtab->Qfunc, Qtab->Qfunc_arg);
		GMRFLib_pardiso_chol(store2);

		int view = 0;
		double *x = Calloc(g->n, double);
		double *b = Calloc(g->n, double);

		for (i = 0; i < g->n; i++) {
			b[i] = ISQR(i + 1);
		}
		GMRFLib_pardiso_solve_LLT(store2, x, b);
		if (view) {
			printf("solve LLT\n");
			for (i = 0; i < g->n; i++) {
				printf("i %d  x= %f  b=%f\n", i, x[i], b[i]);
			}
		}

		for (i = 0; i < g->n; i++) {
			b[i] = ISQR(i + 1);
		}
		GMRFLib_pardiso_solve_LT(store2, x, b);
		if (view) {
			printf("solve LT\n");
			for (i = 0; i < g->n; i++) {
				printf("i %d  x= %f  b=%f\n", i, x[i], b[i]);
			}
		}

		for (i = 0; i < g->n; i++) {
			b[i] = ISQR(i + 1);
		}
		GMRFLib_pardiso_solve_L(store2, x, b);
		if (view) {
			printf("solve L\n");
			for (i = 0; i < g->n; i++) {
				printf("i %d  x= %f  b=%f\n", i, x[i], b[i]);
			}
		}

		GMRFLib_pardiso_Qinv(store2);
		GMRFLib_pardiso_free(&store2);
		Free(x);
		Free(b);
	}
	
	printf("call Qinv\n");
	GMRFLib_pardiso_Qinv(store);

	printf("call free\n");
	GMRFLib_free_tabulate_Qfunc(Qtab);
	GMRFLib_pardiso_free(&store);
	GMRFLib_free_csr(&csr);
	GMRFLib_free_csr(&csr2);
	GMRFLib_free_graph(g);
	Free(perm);

	exit(0);
}

#undef WARNING
