
/* smtp-pardiso.c
 * 
 * Copyright (C) 2018-2020 Havard Rue
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

#include <time.h>
#include <math.h>
#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#if !defined(__FreeBSD__)
#include <malloc.h>
#endif

#include "GMRFLib/hashP.h"
#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

#ifndef GITCOMMIT
#define GITCOMMIT
#endif
static const char GitID[] = "file: " __FILE__ "  " GITCOMMIT;

// do not change: also inlaprog/src/libpardiso.c uses this code
#define NOLIB_ECODE (270465)

typedef struct {
	int verbose;
	int s_verbose;
	int debug;
	int csr_check;
	int mtype;
	int msglvl;
	int nrhs_max;
	int parallel_reordering;
	int *busy;
	GMRFLib_pardiso_store_tp **static_pstores;
} GMRFLib_static_pardiso_tp;

GMRFLib_static_pardiso_tp S = {
	0,						       // verbose
	0,						       // s_verbose
	0,						       // debug
	0,						       // csr_check
	-2,						       // mtype (-2 = sym, 2 = sym pos def)
	0,						       // msg-level (0: no, 1: yes)
	1,						       // maximum number of rhs
	0,						       // parallel reordering? no
	NULL,						       // busy
	NULL
};


#define PSTORES_NUM (16384)

int GMRFLib_pardiso_set_parallel_reordering(int value)
{
	if (value > 0) {
		S.parallel_reordering = (value ? 1 : 0);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_set_nrhs(int nrhs)
{
	if (nrhs > 0) {
		S.nrhs_max = nrhs;
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_set_verbose(int verbose)
{
	if (verbose > 0) {
		S.verbose = (verbose ? 1 : 0);
		S.msglvl = (verbose ? 1 : 0);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_set_debug(int debug)
{
	if (debug) {
		S.s_verbose = (debug ? 1 : 0);
		S.debug = (debug ? 1 : 0);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_csr_free(GMRFLib_csr_tp ** csr)
{
	if (*csr) {
		Free((*csr)->ia);
		Free((*csr)->ja);
		if ((*csr)->copy_only) {
			// do not free
		} else {
			Free((*csr)->a);
		}
		Free(*csr);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_csr_duplicate(GMRFLib_csr_tp ** csr_to, GMRFLib_csr_tp * csr_from)
{
	if (csr_from == NULL) {
		*csr_to = NULL;
		return GMRFLib_SUCCESS;
	}

	*csr_to = Calloc(1, GMRFLib_csr_tp);
	(*csr_to)->base = csr_from->base;
	(*csr_to)->n = csr_from->n;
	(*csr_to)->na = csr_from->na;
	(*csr_to)->copy_only = csr_from->copy_only;

	(*csr_to)->ia = Calloc(csr_from->n + 1, int);
	memcpy((void *) ((*csr_to)->ia), (void *) (csr_from->ia), (size_t) (csr_from->n + 1) * sizeof(int));

	(*csr_to)->ja = Calloc(csr_from->na, int);
	memcpy((void *) ((*csr_to)->ja), (void *) (csr_from->ja), (size_t) (csr_from->na) * sizeof(int));

	if (csr_from->copy_only) {
		(*csr_to)->a = csr_from->a;
	} else {
		(*csr_to)->a = Calloc(csr_from->na, double);
		memcpy((void *) ((*csr_to)->a), (void *) (csr_from->a), (size_t) (csr_from->na) * sizeof(double));
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_csr_base(int base, GMRFLib_csr_tp * M)
{
	int imin = GMRFLib_imin_value(M->ia, M->n + 1, NULL);
	int jmin = GMRFLib_imin_value(M->ja, M->na, NULL);
	int ijmin = IMIN(imin, jmin);

	if (ijmin == base) {
		// all ok. just make sure the base _was_ correct
		assert(base == M->base);
	} else {
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

int GMRFLib_csr_check(GMRFLib_csr_tp * M)
{
	int mtype = S.mtype, error = 0;
	GMRFLib_csr_tp *Q;

	GMRFLib_csr_duplicate(&Q, M);
	GMRFLib_csr_base(1, Q);
	pardiso_chkmatrix(&mtype, &(Q->n), Q->a, Q->ia, Q->ja, &error);
	GMRFLib_csr_free(&Q);
	if (error != 0) {
		GMRFLib_ERROR(GMRFLib_ESNH);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_Q2csr(GMRFLib_csr_tp ** csr, GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg)
{
	// create a upper triangular csr matrix from Q
#define M (*csr)
	int i, j, k, n, na, nnz, nan_error = 0;

	M = Calloc(1, GMRFLib_csr_tp);
	M->base = 0;
	n = graph->n;
	GMRFLib_graph_nnodes(&nnz, graph);		       // symmetric
	na = (nnz - n) / 2 + n;				       // only upper triangular. yes, integer division
	M->na = na;
	M->n = n;
	M->ja = Calloc(na, int);
	M->ia = Calloc(n + 1, int);

	M->ia[0] = 0;
	for (i = k = 0; i < n; i++) {
		M->ja[k++] = i;
		if (graph->lnnbs[i]) {
			memcpy(&(M->ja[k]), graph->lnbs[i], graph->lnnbs[i] * sizeof(int));
			k += graph->lnnbs[i];
		}
		M->ia[i + 1] = M->ia[i] + (1 + graph->lnnbs[i]);
	}
	assert(M->ia[n] == na);

	// when this is true, we can just copy the pointer to the matrix.
	int used_fast_tab = 0;
	if (Qfunc == GMRFLib_tabulate_Qfunction) {
		GMRFLib_tabulate_Qfunc_arg_tp *arg = (GMRFLib_tabulate_Qfunc_arg_tp *) Qfunc_arg;
		if (arg->Q) {
			M->a = arg->Q->a;
			// mark this a copy only, not to be free'd.
			M->copy_only = 1;
			used_fast_tab = 1;
		}
	}
	if (!used_fast_tab) {
		M->a = Calloc(na, double);
	}

	int id_save = GMRFLib_thread_id;
	if (!used_fast_tab) {
		// a bit more manual work
		double val = Qfunc(0, -1, &(M->a[0]), Qfunc_arg);
		if (ISNAN(val)) {
#pragma omp parallel for private(i, k, j) num_threads(GMRFLib_openmp->max_threads_inner)
			for (i = 0; i < n; i++) {
				GMRFLib_thread_id = id_save;
				for (k = M->ia[i]; k < M->ia[i + 1]; k++) {
					j = M->ja[k];
					M->a[k] = Qfunc(i, j, NULL, Qfunc_arg);
				}
			}
		} else {
#pragma omp parallel for private(i, k) num_threads(GMRFLib_openmp->max_threads_inner)
			for (i = 0; i < n; i++) {
				GMRFLib_thread_id = id_save;
				double v;
				k = M->ia[i];
				v = Qfunc(i, -1, &(M->a[k]), Qfunc_arg);
				assert(!(ISNAN(v)));
			}
		}
	}
	GMRFLib_thread_id = id_save;

	for (i = 0; i < na; i++) {
		GMRFLib_STOP_IF_NAN_OR_INF(M->a[i], i, -1);
		if (nan_error)
			break;
	}

	if (nan_error) {
		return !GMRFLib_SUCCESS;
	}
#undef M
	return GMRFLib_SUCCESS;
}

int GMRFLib_csr_write(char *filename, GMRFLib_csr_tp * csr)
{
	// write to file
	GMRFLib_io_tp *io = NULL;
	GMRFLib_io_open(&io, filename, "wb");
	GMRFLib_io_write(io, (const void *) &(csr->n), sizeof(int));
	GMRFLib_io_write(io, (const void *) &(csr->na), sizeof(int));
	GMRFLib_io_write(io, (const void *) &(csr->base), sizeof(int));
	GMRFLib_io_write(io, (const void *) (csr->ia), sizeof(int) * (csr->n + 1));
	GMRFLib_io_write(io, (const void *) (csr->ja), sizeof(int) * csr->na);
	GMRFLib_io_write(io, (const void *) (csr->a), sizeof(double) * csr->na);
	GMRFLib_io_close(io);

	return GMRFLib_SUCCESS;
}

int GMRFLib_csr_read(char *filename, GMRFLib_csr_tp ** csr)
{
	// write to file
	GMRFLib_io_tp *io = NULL;
#define M (*csr)
	M = Calloc(1, GMRFLib_csr_tp);

	GMRFLib_io_open(&io, filename, "rb");
	GMRFLib_io_read(io, (void *) &(M->n), sizeof(int));
	GMRFLib_io_read(io, (void *) &(M->na), sizeof(int));
	GMRFLib_io_read(io, (void *) &(M->base), sizeof(int));

	M->ia = Calloc(M->n + 1, int);
	GMRFLib_io_read(io, (void *) (M->ia), sizeof(int) * (M->n + 1));

	M->ja = Calloc(M->na, int);
	GMRFLib_io_read(io, (void *) (M->ja), sizeof(int) * M->na);

	M->a = Calloc(M->na, double);
	GMRFLib_io_read(io, (void *) (M->a), sizeof(double) * M->na);

	GMRFLib_io_close(io);
#undef M
	return GMRFLib_SUCCESS;
}


int GMRFLib_csr_print(FILE * fp, GMRFLib_csr_tp * csr)
{
	int i, j, k, jj, nnb;
	double value;

	fprintf(fp, "Q->base = %1d, Q->n = %1d, Q->na = %1d copy_only= %1d \n", csr->base, csr->n, csr->na, csr->copy_only);
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
			iarr[k] = i;
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
	int inla_ncpu(void);				       /* external function */
	GMRFLib_pardiso_store_tp *s = Calloc(1, GMRFLib_pardiso_store_tp);

	if (S.s_verbose) {
		PP("_pardiso_init()", s);
	}

	s->maxfct = 1;
	s->pstore = Calloc(1, GMRFLib_pardiso_store_pr_thread_tp);
	assert(S.mtype == -2 || S.mtype == 2);
	s->mtype = S.mtype;
	s->msglvl = S.msglvl;
	s->solver = 0;
	s->iparm_default = Calloc(GMRFLib_PARDISO_PLEN, int);
	s->dparm_default = Calloc(GMRFLib_PARDISO_PLEN, double);
	s->iparm_default[0] = 0;			       /* use default values */
	s->iparm_default[2] = GMRFLib_PARDISO_MAX_NUM_THREADS;

	if (S.s_verbose) {
		PPg("_pardiso_init(): num_threads", (double) (s->iparm_default[2]));
	}

	pardisoinit(s->pt, &(s->mtype), &(s->solver), s->iparm_default, s->dparm_default, &error);

	s->iparm_default[1] = 2;			       /* use this so we can have identical solutions */
	s->iparm_default[4] = 0;			       /* use internal reordering */
	s->iparm_default[7] = 4;			       /* maximum number of refinement steps */
	s->iparm_default[10] = 0;			       /* These are the default, but... */
	s->iparm_default[12] = 0;			       /* I need these for the divided LDL^Tx=b solver to work */
	s->iparm_default[20] = 0;			       /* Diagonal pivoting, and... */
	s->iparm_default[23] = 1;			       /* two level scheduling, and... */
	s->iparm_default[24] = 0;			       /* use parallel solve, as... */
	s->iparm_default[27] = S.parallel_reordering;	       /* parallel reordering? */
	s->iparm_default[33] = 1;			       /* I want identical solutions (require iparm_default[1]=2) */

	if (error != 0) {
		if (error == NOLIB_ECODE) {
			GMRFLib_ERROR(GMRFLib_EPARDISO_NO_LIBRARY);
		} else if (error == -10) {
			GMRFLib_ERROR(GMRFLib_EPARDISO_LICENSE_NOTFOUND);
		} else if (error == -11) {
			GMRFLib_ERROR(GMRFLib_EPARDISO_LICENSE_EXPIRED);
		} else if (error == -12) {
			GMRFLib_ERROR(GMRFLib_EPARDISO_LICENSE_ERR_USERNAME);
		} else {
			GMRFLib_ERROR(GMRFLib_EPARDISO_INTERNAL_ERROR);
		}
	}
	s->done_with_init = GMRFLib_TRUE;
	*store = s;

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_setparam(GMRFLib_pardiso_flag_tp flag, GMRFLib_pardiso_store_tp * store)
{
	assert(store->done_with_init == GMRFLib_TRUE);
	memcpy((void *) (store->pstore->iparm), (void *) (store->iparm_default), GMRFLib_PARDISO_PLEN * sizeof(int));
	memcpy((void *) (store->pstore->dparm), (void *) (store->dparm_default), GMRFLib_PARDISO_PLEN * sizeof(double));

	store->pstore->nrhs = 0;
	store->pstore->err_code = 0;
	store->pstore->iparm[2] = GMRFLib_PARDISO_MAX_NUM_THREADS;

	switch (flag) {
	case GMRFLib_PARDISO_FLAG_REORDER:
	case GMRFLib_PARDISO_FLAG_SYMFACT:
		store->pstore->phase = 11;		       // analysis
		store->pstore->iparm[4] = 0;		       /* 0 = compute the permutation */
		store->pstore->iparm[39] = 1;		       /* 1 = return the permutation */
		store->pstore->nrhs = S.nrhs_max;	       /* this is how it is, apparently */
		break;

	case GMRFLib_PARDISO_FLAG_CHOL:
		store->pstore->phase = 22;		       // numerical factorization
		store->pstore->iparm[32] = 1;		       /* determinant */
		store->pstore->iparm[39] = 1;		       /* 1 = return the permutation (does not do that anymore) */
		break;

	case GMRFLib_PARDISO_FLAG_QINV:
		store->pstore->phase = -22;
		store->pstore->iparm[35] = 1;		       /* do not overwrite internal factor L with selected inversion */
		store->pstore->iparm[36] = 0;		       /* return upper triangular Qinv */
		break;

	case GMRFLib_PARDISO_FLAG_SOLVE_L:
		store->pstore->phase = 33;		       // solve
		store->pstore->iparm[25] = (S.mtype == 2 ? 1 : -12);
		break;

	case GMRFLib_PARDISO_FLAG_SOLVE_LT:
		store->pstore->phase = 33;		       // solve
		store->pstore->iparm[25] = (S.mtype == 2 ? 2 : -23);
		break;

	case GMRFLib_PARDISO_FLAG_SOLVE_LLT:
		store->pstore->phase = 33;		       // solve
		store->pstore->iparm[25] = 0;
		break;

	default:
		GMRFLib_ERROR(GMRFLib_ESNH);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_check_install(int quiet, int no_err)
{
	int *iparm = Calloc(GMRFLib_PARDISO_PLEN, int);
	int mtype = S.mtype, err_code = 0, solver = 0;
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
			if (err_code == NOLIB_ECODE) {
				GMRFLib_ERROR(GMRFLib_EPARDISO_NO_LIBRARY);
			} else if (err_code == -10) {
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

double GMRFLib_pardiso_Qfunc_default(int i, int j, double *UNUSED(values), void *arg)
{
	if (i >= 0 && j < 0) {
		return NAN;
	}

	GMRFLib_graph_tp *g = (GMRFLib_graph_tp *) arg;
	return (i == j ? g->n + 2.0 * g->nnbs[i] : -1.0);
}

int GMRFLib_pardiso_reorder(GMRFLib_pardiso_store_tp * store, GMRFLib_graph_tp * graph)
{
	int debug = S.debug;

	assert(store != NULL);
	assert(store->done_with_init == GMRFLib_TRUE);

	if (store->done_with_reorder == GMRFLib_TRUE) {
		return GMRFLib_SUCCESS;
	}
	GMRFLib_ENTER_ROUTINE;

	int i, n, mnum1 = 1;
	GMRFLib_csr_tp *Q = NULL;

	GMRFLib_graph_duplicate(&(store->graph), graph);
	GMRFLib_pardiso_setparam(GMRFLib_PARDISO_FLAG_REORDER, store);
	GMRFLib_Q2csr(&Q, store->graph, GMRFLib_pardiso_Qfunc_default, (void *) store->graph);
	GMRFLib_csr_base(1, Q);

	if (S.csr_check) {
		GMRFLib_csr_check(Q);
	}
	if (0 && S.debug) {
		GMRFLib_csr_print(stdout, Q);
	}

	n = Q->n;
	store->pstore->perm = Calloc(n, int);
	store->pstore->iperm = Calloc(n, int);

	pardiso(store->pt, &(store->maxfct), &mnum1, &(store->mtype),
		&(store->pstore->phase),
		&(Q->n), Q->a, Q->ia, Q->ja, store->pstore->perm,
		&(store->pstore->nrhs), store->pstore->iparm,
		&(store->msglvl), &(store->pstore->dummy), &(store->pstore->dummy), &(store->pstore->err_code), store->pstore->dparm);

	// Just fill it with a dummy (identity) reordering
	for (i = 0; i < n; i++) {
		store->pstore->perm[i] = i;
		store->pstore->iperm[store->pstore->perm[i]] = i;
	}

	if (0 && debug) {
		for (i = 0; i < n; i++) {
			printf("perm[%1d] = %1d | iperm[%1d] = %1d\n", i, store->pstore->perm[i], i, store->pstore->iperm[i]);
		}
	}

	if (store->pstore->err_code) {
		GMRFLib_ERROR(GMRFLib_EPARDISO_INTERNAL_ERROR);
	}

	store->done_with_reorder = GMRFLib_TRUE;
	store->pstore->L_nnz = store->pstore->iparm[17] - 1;
	GMRFLib_csr_free(&Q);

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_perm(double *x, int m, GMRFLib_pardiso_store_tp * store)
{
	return GMRFLib_pardiso_perm_core(x, m, store, 1);
}
int GMRFLib_pardiso_iperm(double *x, int m, GMRFLib_pardiso_store_tp * store)
{
	return GMRFLib_pardiso_perm_core(x, m, store, 0);
}
int GMRFLib_pardiso_perm_core(double *x, int m, GMRFLib_pardiso_store_tp * store, int direction)
{
	int i, j, k, n, *permutation;
	double *xx;

	n = store->graph->n;
	xx = Calloc(n * m, double);
	memcpy(xx, x, n * m * sizeof(double));
	permutation = (direction ? store->pstore->perm : store->pstore->iperm);
	assert(permutation);
	assert(m > 0);

	for (j = 0; j < m; j++) {
		k = j * n;
		for (i = 0; i < n; i++) {
			x[k + i] = xx[k + permutation[i]];
		}
	}
	Free(xx);

	return GMRFLib_SUCCESS;
}
int GMRFLib_pardiso_symfact(GMRFLib_pardiso_store_tp * store)
{
	assert(store->done_with_init == GMRFLib_TRUE);
	assert(store->done_with_reorder == GMRFLib_TRUE);

	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_build(GMRFLib_pardiso_store_tp * store, GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg)
{
	GMRFLib_ENTER_ROUTINE;

	assert(store != NULL);
	assert(store->done_with_init == GMRFLib_TRUE);
	assert(store->done_with_reorder == GMRFLib_TRUE);

	if (store->pstore->done_with_build == GMRFLib_TRUE && store->pstore->Q) {
		GMRFLib_csr_free(&(store->pstore->Q));
	}
	GMRFLib_Q2csr(&(store->pstore->Q), graph, Qfunc, (void *) Qfunc_arg);
	GMRFLib_csr_base(1, store->pstore->Q);

	if (S.csr_check) {
		GMRFLib_csr_check(store->pstore->Q);
	}
	if (0 && S.debug) {
		GMRFLib_csr_print(stdout, store->pstore->Q);
	}

	store->pstore->done_with_build = GMRFLib_TRUE;
	GMRFLib_LEAVE_ROUTINE;

	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_chol(GMRFLib_pardiso_store_tp * store)
{
	int debug = S.debug;

	GMRFLib_ENTER_ROUTINE;
	assert(store->done_with_init == GMRFLib_TRUE);
	assert(store->done_with_reorder == GMRFLib_TRUE);
	assert(store->pstore->done_with_build == GMRFLib_TRUE);
	assert(store->pstore->Q != NULL);

	int mnum1 = 1, n = store->pstore->Q->n, i;
	GMRFLib_pardiso_setparam(GMRFLib_PARDISO_FLAG_CHOL, store);

	if (debug) {
		printf("CHOL: level %d NUM_THREADS %d iparm[2] %d\n", omp_get_level(), omp_get_num_threads(), store->pstore->iparm[2]);
	}

	if (GMRFLib_openmp->adaptive && omp_get_level() == 0) {
		// this is the exception of the rule, as we want to run this in parallel if we are in adaptive model and
		// level=0.
		omp_set_num_threads(store->pstore->iparm[2]);
	} else {
		omp_set_num_threads(GMRFLib_openmp->max_threads_inner);
		assert(GMRFLib_openmp->max_threads_inner <= store->pstore->iparm[2]);
	}

	pardiso(store->pt, &(store->maxfct), &mnum1, &(store->mtype), &(store->pstore->phase),
		&n, store->pstore->Q->a, store->pstore->Q->ia, store->pstore->Q->ja,
		store->pstore->perm, &(store->pstore->nrhs),
		store->pstore->iparm, &(store->msglvl), NULL, NULL, &(store->pstore->err_code), store->pstore->dparm);

	if (debug) {
		printf("Average number of non-zeros in L per row %.2f\n", store->pstore->iparm[17] / (double) n);
	}
	// Have to check if we need to revert back to C indexing
	int perm_min = GMRFLib_imin_value(store->pstore->perm, n, NULL);
	assert(perm_min == 0 || perm_min == 1);		       /* must either be C or F... */
	if (perm_min == 1) {
		for (i = 0; i < n; i++) {
			store->pstore->perm[i]--;	       /* back to C indexing */
		}
	}
	for (i = 0; i < n; i++) {
		store->pstore->iperm[store->pstore->perm[i]] = i;
	}

	if (0 && debug) {
		for (i = 0; i < n; i++) {
			printf("perm[%1d] = %1d | iperm[%1d] = %1d\n", i, store->pstore->perm[i], i, store->pstore->iperm[i]);
		}
	}

	if (store->pstore->err_code != 0 || store->pstore->iparm[22] > 0) {
		if (store->pstore->iparm[22] > 0) {
			printf("\n*** PARDISO ERROR(%1d): not pos.def matrix: %1d eigenvalues are negative.\n",
			       store->pstore->err_code, store->pstore->iparm[22]);
		} else {
			printf("\n*** PARDISO ERROR(%1d): check manual.\n", store->pstore->err_code);
		}
		printf("*** PARDISO ERROR: I will try to work around the problem...\n\n");
		return GMRFLib_EPOSDEF;
	}

	store->pstore->log_det_Q = store->pstore->dparm[32];
	store->pstore->done_with_chol = GMRFLib_TRUE;

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_solve_core(GMRFLib_pardiso_store_tp * store, GMRFLib_pardiso_flag_tp flag, double *x, double *b, int nrhs)
{
	assert(nrhs > 0);
	assert(store->done_with_init == GMRFLib_TRUE);
	assert(store->done_with_reorder == GMRFLib_TRUE);
	assert(store->pstore->done_with_build == GMRFLib_TRUE);
	assert(store->pstore->done_with_chol == GMRFLib_TRUE);

	// this is so that the RHS can be overwritten
	int n = store->graph->n, mnum1 = 1, i, offset, nblock, reminder, local_nrhs, max_nrhs, idum = 0;
	div_t d;
	double *bb = NULL;

	max_nrhs = S.nrhs_max;
	d = div(nrhs, max_nrhs);
	nblock = d.quot;
	reminder = (d.rem != 0);
	bb = Calloc(nrhs * n, double);
	memcpy((void *) bb, (void *) b, n * nrhs * sizeof(double));

	// we might want to tweak the number of threads here, even do this in parallel for many rhs when the version is
	// thread-safe
	GMRFLib_pardiso_setparam(flag, store);
	for (i = 0; i < nblock + reminder; i++) {
		offset = i * n * max_nrhs;
		local_nrhs = (i < nblock ? max_nrhs : (int) d.rem);
		pardiso(store->pt, &(store->maxfct), &mnum1, &(store->mtype), &(store->pstore->phase),
			&n, store->pstore->Q->a, store->pstore->Q->ia, store->pstore->Q->ja,
			&idum, &local_nrhs, store->pstore->iparm, &(store->msglvl),
			bb + offset, x + offset, &(store->pstore->err_code), store->pstore->dparm);
		if (store->pstore->err_code != 0) {
			GMRFLib_ERROR(GMRFLib_EPARDISO_INTERNAL_ERROR);
		}
	}
	Free(bb);

	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_solve_L(GMRFLib_pardiso_store_tp * store, double *x, double *b, int nrhs)
{
	GMRFLib_ENTER_ROUTINE;
	int res = GMRFLib_pardiso_solve_core(store, GMRFLib_PARDISO_FLAG_SOLVE_L, x, b, nrhs);
	GMRFLib_pardiso_iperm(x, nrhs, store);
	GMRFLib_LEAVE_ROUTINE;

	return res;
}

int GMRFLib_pardiso_solve_LT(GMRFLib_pardiso_store_tp * store, double *x, double *b, int nrhs)
{
	GMRFLib_ENTER_ROUTINE;
	GMRFLib_pardiso_perm(b, nrhs, store);
	int res = GMRFLib_pardiso_solve_core(store, GMRFLib_PARDISO_FLAG_SOLVE_LT, x, b, nrhs);
	if (x != b) {
		GMRFLib_pardiso_iperm(b, nrhs, store);
	}
	GMRFLib_LEAVE_ROUTINE;

	return res;
}

int GMRFLib_pardiso_solve_LLT(GMRFLib_pardiso_store_tp * store, double *x, double *b, int nrhs)
{
	GMRFLib_ENTER_ROUTINE;
	int res = GMRFLib_pardiso_solve_core(store, GMRFLib_PARDISO_FLAG_SOLVE_LLT, x, b, nrhs);
	GMRFLib_LEAVE_ROUTINE;

	return res;
}

double GMRFLib_pardiso_logdet(GMRFLib_pardiso_store_tp * store)
{
	return (store->pstore->log_det_Q);
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

	GMRFLib_csr_tp *Qi = problem->sub_sm_fact.PARDISO_fact->pstore->Qinv;
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
		// map_ii_set(subQinv->mapping, i, i);
		// printf("Mapping %d %d\n", i, problem->sub_graph->mothergraph_idx[i]);
		map_ii_set(subQinv->mapping, problem->sub_graph->mothergraph_idx[i], problem->sub_graph->mothergraph_idx[i]);
	}

	problem->sub_inverse = subQinv;
	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_Qinv(GMRFLib_pardiso_store_tp * store)
{
	GMRFLib_ENTER_ROUTINE;

	assert(store->done_with_reorder == GMRFLib_TRUE);
	assert(store->pstore->done_with_build == GMRFLib_TRUE);
	assert(store->pstore->done_with_chol == GMRFLib_TRUE);

	if (store->pstore->Qinv) {
		GMRFLib_csr_free(&(store->pstore->Qinv));
	}

	GMRFLib_csr_duplicate(&(store->pstore->Qinv), store->pstore->Q);
	GMRFLib_csr_base(1, store->pstore->Qinv);

	GMRFLib_pardiso_setparam(GMRFLib_PARDISO_FLAG_QINV, store);
	int mnum1 = 1;

	if (GMRFLib_openmp->adaptive && omp_get_level() == 0) {
		// this is the exception of the rule, as we want to run this in parallel if we are in adaptive model and
		// level=0.
		omp_set_num_threads(store->pstore->iparm[2]);
	} else {
		omp_set_num_threads(GMRFLib_openmp->max_threads_inner);
		assert(GMRFLib_openmp->max_threads_inner <= store->pstore->iparm[2]);
	}

	pardiso(store->pt, &(store->maxfct), &mnum1, &(store->mtype), &(store->pstore->phase),
		&(store->pstore->Qinv->n),
		store->pstore->Qinv->a, store->pstore->Qinv->ia, store->pstore->Qinv->ja,
		NULL, &(store->pstore->nrhs), store->pstore->iparm, &(store->msglvl), NULL, NULL, &(store->pstore->err_code), store->pstore->dparm);

	if (store->pstore->err_code != 0) {
		GMRFLib_ERROR(GMRFLib_EPARDISO_INTERNAL_ERROR);
	}

	GMRFLib_csr_base(0, store->pstore->Qinv);
	if (0 && S.debug) {
		GMRFLib_csr_print(stdout, store->pstore->Qinv);
	}

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_free(GMRFLib_pardiso_store_tp ** store)
{
	if (store == NULL || *store == NULL) {
		return GMRFLib_SUCCESS;
	}
	GMRFLib_ENTER_ROUTINE;

	if ((*store)->copy_pardiso_ptr) {
		// this is special
		if (S.s_verbose) {
			FIXME("Free pardiso store with copy_pardiso_ptr = 1");
		}
		Free((*store)->pstore->perm);
		Free((*store)->pstore->iperm);
		Free((*store)->pstore);
		Free((*store));
		*store = NULL;

		GMRFLib_LEAVE_ROUTINE;
		return GMRFLib_SUCCESS;
	}

	if (S.s_verbose) {
		PP("free: old=", *store);
	}

	int found = 0, i;
	if (S.static_pstores != NULL) {
		if (S.s_verbose) {
			for (i = 0; i < PSTORES_NUM; i++) {
				if (S.busy[i]) {
					printf("in store: i=%1d s=%p\n", i, (void *) S.static_pstores[i]);
				}
			}
		}
		for (i = 0; i < PSTORES_NUM && !found; i++) {

			if (S.static_pstores[i] == *store) {
				found = 1;
				if (S.busy[i]) {
					S.busy[i] = 0;
					if (S.s_verbose) {
						printf("==> S.busy[%1d] = 1\n", i);
						PP("S.static_pstores[i]", S.static_pstores[i]);
						PP("*store", *store);
						printf("==> free store[%1d]\n", i);
					}
				} else {
					if (S.s_verbose) {
						printf("==> this one is already free [%1d]. ignore\n", i);
					}
				}
			}
		}
	}
	if (!found) {
		if (S.s_verbose) {
			printf("==> free manually as not found\n");
		}

		if ((*store)->pstore) {
			(*store)->pstore->phase = -1;
			int mnum1 = 1;
			pardiso((*store)->pt, &((*store)->maxfct), &mnum1, &((*store)->mtype),
				&((*store)->pstore->phase),
				&((*store)->pstore->idummy),
				&((*store)->pstore->dummy), &((*store)->pstore->idummy),
				&((*store)->pstore->idummy),
				&((*store)->pstore->idummy),
				&((*store)->pstore->nrhs), (*store)->pstore->iparm, &((*store)->msglvl),
				NULL, NULL, &((*store)->pstore->err_code), (*store)->pstore->dparm);

			GMRFLib_csr_free(&((*store)->pstore->Q));
			GMRFLib_csr_free(&((*store)->pstore->Qinv));
			Free((*store)->pstore);
		}

		GMRFLib_graph_free((*store)->graph);
		Free((*store)->iparm_default);
		Free((*store)->dparm_default);
		Free(*store);
	}

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_duplicate_pardiso_store(GMRFLib_pardiso_store_tp ** new, GMRFLib_pardiso_store_tp * old, int UNUSED(copy_ptr), int copy_pardiso_ptr)
{
	// if copy_pardiso_ptr, then copy the ptr to read-only objects. 'copy_ptr' is NOT USED

	int debug = S.debug, failsafe_mode = 0;
	if (old == NULL) {
		*new = NULL;
		return GMRFLib_SUCCESS;
	}

	if (failsafe_mode) {
		FIXME("-->duplicate by creating a new one each time");
		GMRFLib_pardiso_init(new);
		GMRFLib_pardiso_reorder(*new, old->graph);
		return GMRFLib_SUCCESS;
	}
	GMRFLib_ENTER_ROUTINE;

	if (copy_pardiso_ptr) {
		int i;

#define CP(_what) dup->_what = old->_what
#define CP2(_what) dup->pstore->_what = old->pstore->_what
#define CPv(_what, type, len)						\
		if (old->pstore->_what) {				\
			dup->pstore->_what = Calloc(len, type);		\
			memcpy((void *) (dup->pstore->_what), (void *) (old->pstore->_what), (len) * sizeof(type)); \
		} else {						\
			dup->pstore->_what = NULL;			\
		}							\

		GMRFLib_pardiso_store_tp *dup = Calloc(1, GMRFLib_pardiso_store_tp);
		dup->copy_pardiso_ptr = 1;		       /* YES! */
		for (i = 0; i < GMRFLib_PARDISO_PLEN; i++) {
			CP(pt[i]);
		}
		CP(iparm_default);
		CP(dparm_default);
		CP(maxfct);
		CP(done_with_init);
		CP(done_with_reorder);
		CP(msglvl);
		CP(mtype);
		CP(solver);
		CP(graph);

		dup->pstore = Calloc(1, GMRFLib_pardiso_store_pr_thread_tp);
		for (i = 0; i < GMRFLib_PARDISO_PLEN; i++) {
			CP2(iparm[i]);
			CP2(dparm[i]);
		}
		CP2(done_with_build);
		CP2(done_with_chol);
		CP2(dummy);
		CP2(err_code);
		CP2(idummy);
		CP2(nrhs);
		CP2(phase);
		CP2(L_nnz);

		CPv(perm, int, old->graph->n);		       // CP2(perm);
		CPv(iperm, int, old->graph->n);		       // CP2(iperm);

		CP2(log_det_Q);
		CP2(Q);
		CP2(Qinv);
#undef CP
#undef CP2
#undef CPv
		*new = dup;
		return GMRFLib_SUCCESS;
	}

	if (S.static_pstores == NULL) {
#pragma omp critical
		{
			if (S.static_pstores == NULL) {
				S.static_pstores = Calloc(PSTORES_NUM, GMRFLib_pardiso_store_tp *);
				S.busy = Calloc(PSTORES_NUM, int);
				if (S.s_verbose) {
					printf("==> init static_pstores\n");
				}
			}
		}
	}

	int found = 0, idx = -1, ok = 0;
#pragma omp critical
	{
		for (int i = 0; i < PSTORES_NUM && !found; i++) {
			if (!S.busy[i]) {
				S.busy[i] = 1;
				idx = i;
				found = 1;
			}
		}
	}

	assert(found == 1);
	if (S.static_pstores[idx]) {
		if (debug) {
			printf("%s:%1d: static_pstores...iparm[2] = %1d\n", __FILE__, __LINE__, S.static_pstores[idx]->pstore->iparm[2]);
			printf("%s:%1d: level %d max_threads_nested = %1d\n", __FILE__, __LINE__, omp_get_level(),
			       GMRFLib_openmp->max_threads_nested[1]);
		}
		ok = (S.static_pstores[idx]->pstore->iparm[2] >= GMRFLib_openmp->max_threads_nested[1]);
	} else {
		ok = 1;
	}
	if (!ok) {
		P(S.static_pstores[idx]->pstore->iparm[2]);
		P(GMRFLib_openmp->max_threads_nested[1]);
		FIXME("THIS IS NOT TRUE: iparm[2] >= threads_nested[1]");
	}

	if (S.static_pstores[idx] && ok) {
		*new = S.static_pstores[idx];
		if (S.s_verbose) {
			printf("==> reuse store[%1d]\n", idx);
		}
	} else {
		GMRFLib_pardiso_init(&(S.static_pstores[idx]));
		GMRFLib_pardiso_reorder(S.static_pstores[idx], old->graph);
		*new = S.static_pstores[idx];
		if (S.s_verbose) {
			printf("==> new store[%1d]\n", idx);
		}
	}

	if (S.s_verbose) {
		printf("duplicate: new=%p old=%p i=%1d\n", *((void **) new), ((void *) old), idx);
	}

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}



// **********************************************************************
// **********************************************************************
// **********************************************************************
// *** 
// *** test-programs
// *** 
// **********************************************************************
// **********************************************************************
// **********************************************************************


double my_pardiso_test_Q(int i, int j, double *UNUSED(values), void *arg)
{
	GMRFLib_graph_tp *graph = (GMRFLib_graph_tp *) arg;
	return (i == j ? graph->n + i : -1.0);
}

int my_pardiso_test1(void)
{
	int err = 0;

	if (1) {
		err = GMRFLib_pardiso_check_install(1, 0);
		if (err == GMRFLib_SUCCESS) {
			printf("PARDISO OK\n");
		} else {
			printf("PARDISO FAIL\n");
		}
	}

	GMRFLib_tabulate_Qfunc_tp *Qtab;
	GMRFLib_graph_tp *g;

	GMRFLib_openmp->strategy = GMRFLib_OPENMP_STRATEGY_PARDISO;
	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_OPTIMIZE, NULL, NULL);

	// GMRFLib_tabulate_Qfunc_from_file(&Qtab, &g, "Qdense.txt", -1, NULL, NULL, NULL);
	// GMRFLib_tabulate_Qfunc_from_file(&Qtab, &g, "I5.txt", -1, NULL, NULL, NULL);
	GMRFLib_tabulate_Qfunc_from_file(&Qtab, &g, "Q.txt", -1, NULL, NULL, NULL);

	GMRFLib_csr_tp *csr, *csr2;
	GMRFLib_Q2csr(&csr, g, Qtab->Qfunc, Qtab->Qfunc_arg);
	GMRFLib_csr_print(stdout, csr);

	GMRFLib_csr_duplicate(&csr2, csr);
	// GMRFLib_csr_print(stdout, csr2);

	GMRFLib_csr2Q(&Qtab, &g, csr2);
	// GMRFLib_Qfunc_print(stdout, g, Qtab->Qfunc, Qtab->Qfunc_arg);

	int *perm = NULL;
	int i, k, nrhs;

	// GMRFLib_graph_printf(stdout, g);
	GMRFLib_pardiso_store_tp *store = NULL;
	GMRFLib_pardiso_init(&store);
	GMRFLib_pardiso_reorder(store, g);
	GMRFLib_pardiso_build(store, g, Qtab->Qfunc, Qtab->Qfunc_arg);
	GMRFLib_pardiso_chol(store);
	GMRFLib_pardiso_Qinv(store);

	P(GMRFLib_openmp->max_threads_outer);
	P(GMRFLib_openmp->max_threads_inner);

	nrhs = 1;
	// S.s_verbose=1;
#pragma omp parallel for private(k) num_threads(GMRFLib_openmp->max_threads_outer)
	for (k = 0; k < 100; k++) {
		printf("this is k= %d from thread %d\n", k, omp_get_thread_num());

		GMRFLib_pardiso_store_tp *store2 = NULL;
		if (0) {
			GMRFLib_pardiso_init(&store2);
			GMRFLib_pardiso_reorder(store2, g);
			GMRFLib_pardiso_build(store2, g, Qtab->Qfunc, Qtab->Qfunc_arg);
			GMRFLib_pardiso_chol(store2);
		} else {
			GMRFLib_duplicate_pardiso_store(&store2, store, -1, 1);
		}

		int view = 1;
		double *x = Calloc(g->n * nrhs, double);
		double *b = Calloc(g->n * nrhs, double);

		for (i = 0; i < g->n * nrhs; i++) {
			b[i] = ISQR(i + 1);
		}
		GMRFLib_pardiso_solve_LLT(store2, x, b, nrhs);
		if (view) {
			printf("solve LLT\n");
			for (int ii = 0; ii < nrhs; ii++) {
				for (i = 0; i < g->n; i++) {
					printf("ii %d i %d  x= %f  b=%f\n", ii, i, x[i + ii * g->n], b[i + ii * g->n]);
				}
			}
		}

		for (i = 0; i < g->n * nrhs; i++) {
			b[i] = ISQR(i + 1);
		}
		GMRFLib_pardiso_solve_LT(store2, x, b, 1);
		if (view) {
			printf("solve LT\n");
			for (int ii = 0; ii < nrhs; ii++) {
				for (i = 0; i < g->n; i++) {
					printf("ii %d i %d  x= %f  b=%f\n", ii, i, x[i + ii * g->n], b[i + ii * g->n]);
				}
			}
		}

		for (i = 0; i < g->n * nrhs; i++) {
			b[i] = ISQR(i + 1);
		}
		GMRFLib_pardiso_solve_L(store2, x, b, 1);
		if (view) {
			printf("solve L\n");
			for (int ii = 0; ii < nrhs; ii++) {
				for (i = 0; i < g->n; i++) {
					printf("ii %d i %d  x= %f  b=%f\n", ii, i, x[i + ii * g->n], b[i + ii * g->n]);
				}
			}
		}

		GMRFLib_pardiso_free(&store2);
		Free(x);
		Free(b);
	}

	printf("call free\n");
	GMRFLib_free_tabulate_Qfunc(Qtab);
	GMRFLib_pardiso_free(&store);
	GMRFLib_csr_free(&csr);
	GMRFLib_csr_free(&csr2);
	GMRFLib_graph_free(g);
	Free(perm);

	exit(0);
}

int my_pardiso_test2(void)
{
	int n = 5, m = 1, nc = 1, i;
	double *var;

	GMRFLib_graph_tp *graph = NULL;
	GMRFLib_graph_mk_linear(&graph, n, m, 0);
	GMRFLib_problem_tp *problem = NULL;
	GMRFLib_constr_tp *constr = NULL;
	GMRFLib_make_empty_constr(&constr);

	constr->nc = nc;
	constr->a_matrix = Calloc(n * nc, double);
	for (i = 0; i < n * nc; i++)
		constr->a_matrix[i] = GMRFLib_uniform();
	constr->e_vector = Calloc(nc, double);
	for (i = 0; i < nc; i++)
		constr->e_vector[i] = GMRFLib_uniform();
	GMRFLib_prepare_constr(constr, graph, 1);

	// GMRFLib_openmp->strategy = GMRFLib_OPENMP_STRATEGY_PARDISO;
	// GMRFLib_smtp = GMRFLib_SMTP_PARDISO;
	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_BUILD_MODEL, NULL, NULL);

	double *x = Calloc(n, double);
	double *b = Calloc(n, double);
	double *c = Calloc(n, double);
	double *mean = Calloc(n, double);
	for (i = 0; i < n; i++) {
		x[i] = GMRFLib_uniform();
		b[i] = GMRFLib_uniform();
		c[i] = exp(GMRFLib_uniform());
		mean[i] = GMRFLib_uniform();
	}

	GMRFLib_init_problem(&problem, x, b, c, mean, graph, my_pardiso_test_Q, (void *) graph, NULL, constr, GMRFLib_NEW_PROBLEM);
	GMRFLib_evaluate(problem);
	GMRFLib_Qinv(problem, GMRFLib_QINV_ALL);

	for (i = 0; i < n; i++) {
		var = GMRFLib_Qinv_get(problem, i, i);
		printf("Qinv[%1d,%1d] = %g\n", i, i, *var);
	}
	GMRFLib_free_problem(problem);

	return 0;
}

int my_pardiso_test3(void)
{
	int err = 0;

	FIXME("this is test3");
	if (1) {
		err = GMRFLib_pardiso_check_install(1, 0);
		if (err == GMRFLib_SUCCESS) {
			printf("PARDISO OK\n");
		} else {
			printf("PARDISO FAIL\n");
		}
	}

	GMRFLib_tabulate_Qfunc_tp *Qtab;
	GMRFLib_graph_tp *g;

	GMRFLib_openmp->strategy = GMRFLib_OPENMP_STRATEGY_PARDISO;
	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_OPTIMIZE, NULL, NULL);

	// GMRFLib_tabulate_Qfunc_from_file(&Qtab, &g, "Qdense.txt", -1, NULL, NULL, NULL);
	// GMRFLib_tabulate_Qfunc_from_file(&Qtab, &g, "I5.txt", -1, NULL, NULL, NULL);
	GMRFLib_tabulate_Qfunc_from_file(&Qtab, &g, "Q.txt", -1, NULL, NULL, NULL);

	GMRFLib_csr_tp *csr;
	GMRFLib_Q2csr(&csr, g, Qtab->Qfunc, Qtab->Qfunc_arg);
	GMRFLib_csr_print(stdout, csr);
	P(csr->n);

	// GMRFLib_graph_printf(stdout, g);
	GMRFLib_pardiso_store_tp *store = NULL;
	GMRFLib_pardiso_init(&store);
	GMRFLib_pardiso_reorder(store, g);
	GMRFLib_pardiso_build(store, g, Qtab->Qfunc, Qtab->Qfunc_arg);
	GMRFLib_pardiso_chol(store);

	GMRFLib_pardiso_Qinv(store);

	P(GMRFLib_openmp->max_threads_outer);
	P(GMRFLib_openmp->max_threads_inner);

	int nrhs = S.nrhs_max, k;
#pragma omp parallel for private(k) num_threads(GMRFLib_openmp->max_threads_outer)
	for (k = 0; k < 1000; k++) {

		int i;
		GMRFLib_pardiso_store_tp *local_store = NULL;

		printf("this is k= %d from thread %d\n", k, omp_get_thread_num());
		GMRFLib_duplicate_pardiso_store(&local_store, store, -1, 1);

		int view = 1;
		double *x = Calloc(g->n * nrhs, double);
		double *b = Calloc(g->n * nrhs, double);

		for (i = 0; i < g->n * nrhs; i++) {
			b[i] = ISQR(i + 1);
		}
		GMRFLib_pardiso_solve_LLT(local_store, x, b, nrhs);
		if (view) {
			printf("solve LLT\n");
			for (int ii = 0; ii < nrhs; ii++) {
				for (i = 0; i < g->n; i++) {
					printf("ii %d i %d  x= %f  b=%f\n", ii, i, x[i + ii * g->n], b[i + ii * g->n]);
				}
			}
		}

		if (1) {
			for (i = 0; i < g->n; i++) {
				b[i] = ISQR(i + 1);
				x[i] = 0.0;
			}

			double *x2 = Calloc(g->n, double);
			GMRFLib_pardiso_solve_L(local_store, x2, b, 1);
			GMRFLib_pardiso_solve_LT(local_store, x, x2, 1);
			if (view) {
				printf("solve L D LT\n");
				for (int ii = 0; ii < nrhs; ii++) {
					for (i = 0; i < g->n; i++) {
						printf("ii %d i %d  x= %f  b=%f\n", ii, i, x[i + ii * g->n], b[i + ii * g->n]);
					}
				}
			}
		}


		for (i = 0; i < g->n * nrhs; i++) {
			b[i] = ISQR(i + 1);
		}
		GMRFLib_pardiso_solve_LT(local_store, x, b, 1);
		if (view) {
			printf("solve LT\n");
			for (int ii = 0; ii < nrhs; ii++) {
				for (i = 0; i < g->n; i++) {
					printf("ii %d i %d  x= %f  b=%f\n", ii, i, x[i + ii * g->n], b[i + ii * g->n]);
				}
			}
		}

		for (i = 0; i < g->n * nrhs; i++) {
			b[i] = ISQR(i + 1);
		}
		GMRFLib_pardiso_solve_L(local_store, x, b, 1);
		if (view) {
			printf("solve L\n");
			for (int ii = 0; ii < nrhs; ii++) {
				for (i = 0; i < g->n; i++) {
					printf("ii %d i %d  x= %f  b=%f\n", ii, i, x[i + ii * g->n], b[i + ii * g->n]);
				}
			}
		}

		GMRFLib_pardiso_free(&local_store);
		Free(x);
		Free(b);
	}

	GMRFLib_free_tabulate_Qfunc(Qtab);
	GMRFLib_pardiso_free(&store);
	GMRFLib_csr_free(&csr);
	GMRFLib_graph_free(g);

	return GMRFLib_SUCCESS;
}

int my_pardiso_test4(void)
{
	int k;
	GMRFLib_tabulate_Qfunc_tp *Qtab;
	GMRFLib_graph_tp *g;
	GMRFLib_csr_tp *csr;

	GMRFLib_smtp = GMRFLib_SMTP_PARDISO;
	GMRFLib_openmp->strategy = GMRFLib_OPENMP_STRATEGY_PARDISO;
	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_OPTIMIZE, NULL, NULL);
	P(GMRFLib_openmp->max_threads_inner);

	GMRFLib_tabulate_Qfunc_from_file(&Qtab, &g, "Q-problem-ijformat.txt", -1, NULL, NULL, NULL);
	GMRFLib_Q2csr(&csr, g, Qtab->Qfunc, Qtab->Qfunc_arg);

	// GMRFLib_csr_write("Q-problem-csr.dat", csr);
	// GMRFLib_csr_read("Q-problem-csr.dat", &csr);

	GMRFLib_pardiso_store_tp *store = NULL;
	GMRFLib_pardiso_init(&store);
	GMRFLib_pardiso_reorder(store, g);
	GMRFLib_pardiso_build(store, g, Qtab->Qfunc, Qtab->Qfunc_arg);

	for (k = 0; k < 1000; k++) {
		GMRFLib_pardiso_chol(store);
		printf("k %d logdet %.12f\n", k, GMRFLib_pardiso_logdet(store));
	}

	exit(0);
}

int my_pardiso_test5(void)
{
	S.msglvl = 1;
	S.csr_check = 1;
	GMRFLib_openmp->strategy = GMRFLib_OPENMP_STRATEGY_PARDISO;
	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_DEFAULT, NULL, NULL);

	P(GMRFLib_openmp->max_threads_outer);
	P(GMRFLib_openmp->max_threads_inner);
	P(omp_get_nested());

	int k;
#pragma omp parallel for private(k) num_threads(GMRFLib_openmp->max_threads_outer)
	for (k = 0; k < 1000; k++) {

		// I do not free anything here...

		P(k);
		GMRFLib_tabulate_Qfunc_tp *Qtab = NULL;
		GMRFLib_graph_tp *g = NULL;

		GMRFLib_tabulate_Qfunc_from_file(&Qtab, &g, "Qsparse2.dat", -1, NULL, NULL, NULL);
		GMRFLib_csr_tp *csr = NULL;
		GMRFLib_Q2csr(&csr, g, Qtab->Qfunc, Qtab->Qfunc_arg);

		GMRFLib_pardiso_store_tp *store = NULL;
		GMRFLib_pardiso_init(&store);
		GMRFLib_pardiso_reorder(store, g);
		GMRFLib_pardiso_build(store, g, Qtab->Qfunc, Qtab->Qfunc_arg);
		GMRFLib_pardiso_chol(store);
		GMRFLib_pardiso_Qinv(store);
	}
	exit(0);
}

int my_pardiso_test6(GMRFLib_ai_store_tp * ai_store, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg, double *c)
{
	int n = ai_store->problem->sub_graph->n;
	int i;

	assert(omp_get_num_threads() == 1);

	GMRFLib_thread_id = 0;
	GMRFLib_tabulate_Qfunc_tp *tab = NULL;
	GMRFLib_problem_tp *problem = ai_store->problem;
	GMRFLib_pardiso_store_tp *pardiso_store = problem->sub_sm_fact.PARDISO_fact;

	GMRFLib_tabulate_Qfunc(&tab, ai_store->problem->sub_graph, Qfunc, Qfunc_arg, NULL, NULL, NULL);
	P(GMRFLib_openmp->max_threads_outer);

#pragma omp parallel for private(i) num_threads(GMRFLib_openmp->max_threads_outer)
	for (i = 0; i < n; i++) {
		int *iparm = Calloc(GMRFLib_PARDISO_PLEN, int);
		double *dparm = Calloc(GMRFLib_PARDISO_PLEN, double);

		memcpy(iparm, pardiso_store->iparm_default, GMRFLib_PARDISO_PLEN * sizeof(int));
		memcpy(dparm, pardiso_store->dparm_default, GMRFLib_PARDISO_PLEN * sizeof(double));

		iparm[7] = 0;
		iparm[25] = 0;

		int j;
		int error = 0;
		int maxfct = 1;
		int mnum = 1;
		int msglvl = 0;
		int mtype = -2;
		int phase = 33;
		int idum = 0;
		int one = 1;

		double *work = Calloc(3 * n, double);
		double *b = work;
		double *x = work + n;
		double *res = work + 2 * n;
		double err;
		double fake_a = 0.0;
		int fake_ia = 0;
		int fake_ja = 0;

		// if I set b[10]=1, then it works in parallel, but not if they are different
		b[i] = 1.0;

		assert(GMRFLib_openmp->max_threads_inner == iparm[2]);
		omp_set_num_threads(iparm[2]);

		pardiso(pardiso_store->pt, &maxfct, &mnum, &mtype, &phase, &n,
			// not in use
			&fake_a,			       // pardiso_store->pstore->Q->a,
			&fake_ia,			       // pardiso_store->pstore->Q->ia,
			&fake_ja,			       // pardiso_store->pstore->Q->ja,
			// 
			&idum, &one, iparm, &msglvl, b, x, &error, dparm);

		// res = Q x
		GMRFLib_Qx2(res, x, problem->sub_graph, tab->Qfunc, tab->Qfunc_arg, c);

		for (err = 0.0, j = 0; j < n; j++) {
			err += SQR(res[j] - b[j]);
		}
		err /= n;
		if (err > 1E-4)
			printf("i %d err %g\n", i, err);

		Free(work);
	}
	exit(0);
}

int my_pardiso_test7(void)
{
	S.msglvl = 0;
	S.csr_check = 1;
	GMRFLib_openmp->strategy = GMRFLib_OPENMP_STRATEGY_PARDISO;
	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_DEFAULT, NULL, NULL);

	P(GMRFLib_openmp->max_threads_outer);
	P(GMRFLib_openmp->max_threads_inner);
	P(omp_get_nested());

	int k;
#pragma omp parallel for private(k) num_threads(GMRFLib_openmp->max_threads_outer)
	for (k = 0; k < 1000; k++) {

		// I do not free anything here...

		omp_set_num_threads(GMRFLib_openmp->max_threads_inner);
		P(k);
		GMRFLib_tabulate_Qfunc_tp *Qtab = NULL;
		GMRFLib_graph_tp *g = NULL;

		if (k < 500)
			GMRFLib_tabulate_Qfunc_from_file(&Qtab, &g, "Q.txt", -1, NULL, NULL, NULL);
		else
			GMRFLib_tabulate_Qfunc_from_file(&Qtab, &g, "I5.txt", -1, NULL, NULL, NULL);
		GMRFLib_csr_tp *csr = NULL;
		GMRFLib_Q2csr(&csr, g, Qtab->Qfunc, Qtab->Qfunc_arg);

		GMRFLib_pardiso_store_tp *store = NULL;
		GMRFLib_pardiso_init(&store);
		GMRFLib_pardiso_reorder(store, g);
		GMRFLib_pardiso_build(store, g, Qtab->Qfunc, Qtab->Qfunc_arg);
		GMRFLib_pardiso_chol(store);
		GMRFLib_pardiso_Qinv(store);
	}
	exit(0);
}
