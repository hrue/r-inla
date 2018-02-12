/* smtp-pardiso.c
 * 
 * Copyright (C) 2018 Havard Rue & Alexander Litvinenko
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
  \brief The implementation of the interface towards PARDISO library.
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

#define WARN(_msg) fprintf(stderr, "\n\n%s:%1d: %s\n\n", __FILE__, __LINE__, _msg)

int GMRFLib_free_csr(GMRFLib_csr_tp **csr)
{
	if (*csr) {
		Free((*csr)->ia);
		Free((*csr)->ja);
		Free((*csr)->a);
		Free(*csr);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_duplicate_csr(GMRFLib_csr_tp ** csr_to, GMRFLib_csr_tp *csr_from)
{
	if (csr_from == NULL)
	{
		*csr_to = NULL;
		return GMRFLib_SUCCESS;
	}

	*csr_to = Calloc(1, GMRFLib_csr_tp);
	(*csr_to)->base = csr_from->base;
	(*csr_to)->n = csr_from->n;
	(*csr_to)->na = csr_from->na;

	(*csr_to)->ia = Calloc(csr_from->n + 1, int);
	memcpy((void *)((*csr_to)->ia), (void *)(csr_from->ia), (size_t) (csr_from->n+1) * sizeof(int));

	(*csr_to)->ja = Calloc(csr_from->na, int);
	memcpy((void *)((*csr_to)->ja), (void *)(csr_from->ja), (size_t) (csr_from->na) * sizeof(int));

	(*csr_to)->a = Calloc(csr_from->na, double);
	memcpy((void *)((*csr_to)->a), (void *)(csr_from->a), (size_t) (csr_from->na) * sizeof(double));
	
	return GMRFLib_SUCCESS;
}

int GMRFLib_csr_base(int base, GMRFLib_csr_tp *M)
{
	if (M->base != base){
		GMRFLib_csr_convert(M);
	}
	return GMRFLib_SUCCESS;
}
int GMRFLib_csr_convert(GMRFLib_csr_tp *M)
{
	int i;
	if (M->base == 1) {
		M->base = 0;
		for (i = 0; i < M->n+1; i++) {
			M->ia[i] -= 1;
		}
		for (i = 0; i < M->na; i++) {
			M->ja[i] -= 1;
		}
	} else if (M->base == 0) {
		M->base = 1;
		for (i = 0; i < M->n+1; i++) {
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

int GMRFLib_Q2csr_check(GMRFLib_csr_tp *M)
{
	int mtype = -2, error=0;
	GMRFLib_csr_tp *Q;

	GMRFLib_duplicate_csr(&Q, M);
	GMRFLib_csr_base(1, Q);
	pardiso_chkmatrix(&mtype , &(Q->n), Q->a, Q->ia, Q->ja, &error);
	GMRFLib_free_csr(&Q);
	if (error != 0) {
		fprintf(stderr, "\n\nMatrix failed check\n");
		exit(1);
	}
	return (error == 0 ? GMRFLib_SUCCESS : !GMRFLib_SUCCESS);
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
	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_setparam(GMRFLib_pardiso_flag_tp flag, GMRFLib_pardiso_store_tp *store)
{
	static int check_install = 1;
	if (check_install) {
		int silent = 0;
		GMRFLib_pardiso_check_install(silent);
		check_install = 0;
	}

	assert(store->init_done == 1);

	store->mtype = -2;
	store->maxfct = 1;				       /* Maximum number of numerical factorizations.  */
	store->mnum   = 1;				       /* Which factorization to use. */
	store->msglvl = 1;				       /* Print statistical information  */
	store->nrhs = 0;
	store->err_code = 0;				       /* Initialize err_code flag */

	switch(flag) {
	case GMRFLib_PARDISO_FLAG_REORDER:
		if (store->Q) WARN("->Q exists");
		store->Q = Calloc(1, GMRFLib_csr_tp);
		store->phase = 11;			       // analysis
		store->iparm[4] = 0;			       /* compute the permutation */
		break;
	case GMRFLib_PARDISO_FLAG_SYMFACT: 
		store->phase = 11;			       // analysis
		break;
	case GMRFLib_PARDISO_FLAG_CHOL: 
		store->phase = 12;			       // Analysis, numerical factorization
		store->iparm[32] = 1;			       /* determinant */
		break;
	case GMRFLib_PARDISO_FLAG_INV: 
		if (store->Qinv) WARN("->Qinv exists");
		store->Qinv = Calloc(1, GMRFLib_csr_tp);
		store->phase = -22;			  /* do not overwrite internal factor L with selected inversion*/ 
		store->iparm[35]  = 1;			  /* do not overwrite internal factor L with selected inversion*/ 
		break;
	case GMRFLib_PARDISO_FLAG_SOLVE: 
		store->iparm[7] = 1;			       /* Max numbers of iterative refinement steps. */
		store->phase = 33;			       // solve
		store->nrhs = 1;
		break;
	default:
		GMRFLib_ERROR(GMRFLib_ESNH);
	}
}

int GMRFLib_pardiso_check_install(int quiet)
{
	void **pt = Calloc(64, void *);
	double *dparm = Calloc(64, double);
	int *iparm = Calloc(64, int);
	int mtype = -2, err_code = 0, solver = 0;

	STDOUT_TO_DEV_NULL_START(1);
	pardisoinit(pt, &mtype, &solver, iparm, dparm, &err_code);
	STDOUT_TO_DEV_NULL_END;

	Free(pt);
	Free(dparm);
	Free(iparm);
	
	if (quiet) {
		return (err_code == 0 ? GMRFLib_SUCCESS : !GMRFLib_SUCCESS);
	} else {
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
			assert(0==1);
		} else {
			return (GMRFLib_SUCCESS);
		}
	}
}

double GMRFLib_pardiso_Qfunc_default(int i, int j, void *arg) 
{
	GMRFLib_graph_tp *g = (GMRFLib_graph_tp *) arg;
	return (i == j ? g->n + 2.0*g->nnbs[i] : -1.0);
}

int GMRFLib_pardiso_init(GMRFLib_pardiso_store_tp *store)
{
	assert(store->init_done == 0);
	
	int error = 0, solver = 0;
	store->mtype = -2;
	store->init_done = 1;

	store->iparm[2] = 4;				       /* num_proc: FIXME LATER */
	pardisoinit(store->pt, &(store->mtype), &solver, store->iparm, store->dparm, &error);
	if (error != 0) {
		P(error);
		exit(1);
	}
	
	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_reorder(GMRFLib_pardiso_store_tp *store, GMRFLib_graph_tp *graph, int *reordering)
{
	assert(store->init_done == 1);

	int i;
	GMRFLib_pardiso_setparam(GMRFLib_PARDISO_FLAG_REORDER, store);
	GMRFLib_Q2csr(&(store->Q), graph, GMRFLib_pardiso_Qfunc_default, (void *) graph);
	GMRFLib_csr_base(1, store->Q);
	GMRFLib_print_csr(stdout, store->Q);
	GMRFLib_Q2csr_check(store->Q);

	if (reordering) {
		store->my_perm = Calloc(store->Q->n, int);
		memcpy((void *) store->my_perm, (void *) reordering, (store->Q->n)*sizeof(int));

		for(i = 0; i< store->Q->n; i++){
			store->my_perm[i]++;		       /* base=1 */
			printf("my_perm[%1d] = %1d\n", i, store->my_perm[i]);
		}
		store->iparm[4] = 1;			       /* use this one */
	}

	pardiso(store->pt, &(store->maxfct), &(store->mnum), &(store->mtype), &(store->phase),
		&(store->Q->n), store->Q->a, store->Q->ia, store->Q->ja, store->my_perm, 
		&(store->nrhs), store->iparm, &(store->msglvl), &(store->dummy), &(store->dummy),
		&(store->err_code), store->dparm);

	P(store->err_code);
	if (store->err_code) {
		P(store->err_code);
		exit(1);
	}
	
	store->L_nnz = store->iparm[17] -1;
	P(store->L_nnz);
	
	GMRFLib_free_csr(&(store->Q));

	
	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_symfact(GMRFLib_pardiso_store_tp *store)
{
	assert(store->init_done == 1);
	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_chol(GMRFLib_pardiso_store_tp *store, GMRFLib_graph_tp *graph,
			 GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_arg)
{
	GMRFLib_pardiso_setparam(GMRFLib_PARDISO_FLAG_CHOL, store);
	GMRFLib_Q2csr(&(store->Q), graph, Qfunc, (void *) Qfunc_arg);
	GMRFLib_csr_base(1, store->Q);
	GMRFLib_Q2csr_check(store->Q);
	GMRFLib_print_csr(stdout, store->Q);

	store->iparm[4] = 1;				       /* use the stored permutation */
	pardiso(store->pt, &(store->maxfct), &(store->mnum), &(store->mtype), &(store->phase),
		&(store->Q->n), store->Q->a, store->Q->ia, store->Q->ja, store->my_perm, &(store->nrhs),
		store->iparm, &(store->msglvl), NULL, NULL, &(store->err_code), store->dparm);
	
	if (store->err_code != 0) {
		printf("\nERROR during numerical factorization: %d", store->err_code);
		exit(2);
	}

	if (store->iparm[22] > 0 || 1) {
		printf("\nERROR not pos def matrix, #neg.eigen %d", store->iparm[22]);
	}
		
	store->log_det_Q = store->dparm[32];
	printf("\nLogDeterminant = %.12g ...\n ", store->log_det_Q);

	return GMRFLib_SUCCESS;
}





//
int pardiso_test()
{
	int err=0;

	err = GMRFLib_pardiso_check_install(1);
	if (err == GMRFLib_SUCCESS) {
		printf("PARDISO OK\n");
	} else {
		printf("PARDISO FAIL\n");
	}

	GMRFLib_tabulate_Qfunc_tp *Qtab;
	GMRFLib_graph_tp *g;

	GMRFLib_tabulate_Qfunc_from_file(&Qtab, &g, "Q1.txt", -1, NULL, NULL, NULL);

	GMRFLib_csr_tp *csr, *csr2;

	GMRFLib_Q2csr(&csr, g, Qtab->Qfunc, Qtab->Qfunc_arg);
	GMRFLib_print_csr(stdout, csr);

	GMRFLib_duplicate_csr(&csr2, csr);
	//GMRFLib_print_csr(stdout, csr2);
	
	GMRFLib_csr2Q(&Qtab, &g, csr2);
	//GMRFLib_print_Qfunc(stdout, g, Qtab->Qfunc, Qtab->Qfunc_arg);

	GMRFLib_pardiso_store_tp *store = Calloc(1, GMRFLib_pardiso_store_tp);

	int *perm = Calloc(g->n, int);
	int i;
	
	GMRFLib_print_graph(stdout, g);
	for(i=0; i<g->n; i++) perm[i] = i;
	GMRFLib_pardiso_init(store);
	GMRFLib_pardiso_reorder(store, g, perm);

	GMRFLib_pardiso_chol(store, g, Qtab->Qfunc,  Qtab->Qfunc_arg);
	
	exit(0);
}

#undef WARN
