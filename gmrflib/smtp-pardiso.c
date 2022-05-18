
/* smtp-pardiso.c
 * 
 * Copyright (C) 2018-2022 Havard Rue
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
	-1,						       // maximum number of rhs
	1,						       // parallel reordering? yes
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
	if (nrhs != 0) {
		S.nrhs_max = nrhs;
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_get_nrhs(void)
{
	return S.nrhs_max;
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
		Free((*csr)->iwork);
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

	int n = csr_from->n;
	int na = csr_from->na;
	int len = n + 1 + na;

	*csr_to = Calloc(1, GMRFLib_csr_tp);
	(*csr_to)->n = n;
	(*csr_to)->na = na;
	(*csr_to)->copy_only = csr_from->copy_only;

	(*csr_to)->iwork = Calloc(2 * len, int);
	(*csr_to)->ia = (*csr_to)->iwork;
	(*csr_to)->ja = (*csr_to)->iwork + n + 1;
	(*csr_to)->ia1 = (*csr_to)->iwork + len;
	(*csr_to)->ja1 = (*csr_to)->iwork + len + n + 1;
	Memcpy((void *) ((*csr_to)->iwork), (void *) (csr_from->iwork), (size_t) (2 * len) * sizeof(int));

	if (csr_from->copy_only) {
		(*csr_to)->a = csr_from->a;
	} else {
		(*csr_to)->a = Calloc(csr_from->na, double);
		Memcpy((void *) ((*csr_to)->a), (void *) (csr_from->a), (size_t) (csr_from->na) * sizeof(double));
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_csr_check(GMRFLib_csr_tp * M)
{
	int mtype = S.mtype, error = 0;

	assert(M);
	pardiso_chkmatrix(&mtype, &(M->n), M->a, M->ia1, M->ja, &error);
	if (error != 0) {
		GMRFLib_ERROR(GMRFLib_ESNH);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_Q2csr(int thread_id, GMRFLib_csr_tp ** csr, GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg)
{
	// create a upper triangular csr matrix from Q
#define M (*csr)
	int i, k, n, na, nan_error = 0, len;

	M = Calloc(1, GMRFLib_csr_tp);
	n = graph->n;
	na = graph->nnz / 2 + n;			       // only upper triangular. yes, integer division
	len = n + 1 + na;
	M->na = na;
	M->n = n;
	M->iwork = Calloc(2 * len, int);
	M->ia = M->iwork;
	M->ja = M->iwork + n + 1;
	M->ia1 = M->iwork + len;
	M->ja1 = M->iwork + len + n + 1;

	// new code. by doing it in two steps we can do the second one in parallel, and this is the one that take time.
	int *k_arr = Calloc(n, int);
	M->ia[0] = 0;
	for (i = k = 0; i < n; i++) {
		M->ja[k++] = i;
		k_arr[i] = k;
		k += graph->lnnbs[i];
		M->ia[i + 1] = M->ia[i] + (1 + graph->lnnbs[i]);
	}
	assert(M->ia[n] == na);

#define CODE_BLOCK							\
	for (int i = 0; i < n; i++) {					\
		if (graph->lnnbs[i]) {					\
			int k = k_arr[i];				\
			Memcpy(&(M->ja[k]), graph->lnbs[i], graph->lnnbs[i] * sizeof(int)); \
		}							\
	}

	RUN_CODE_BLOCK(1, 0, 0);			       /* TO QUICK TO DO IN PARALLEL */
#undef CODE_BLOCK
	Free(k_arr);

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

	if (!used_fast_tab) {
		// a bit more manual work
		double val = Qfunc(thread_id, 0, -1, &(M->a[0]), Qfunc_arg);
		if (ISNAN(val)) {
#define CODE_BLOCK							\
			for (int i = 0; i < n; i++) {			\
				for (int k = M->ia[i]; k < M->ia[i + 1]; k++) { \
					int j = M->ja[k];	\
					M->a[k] = Qfunc(thread_id, i, j, NULL, Qfunc_arg);	\
				}					\
			}

			RUN_CODE_BLOCK((GMRFLib_Qx_strategy ? GMRFLib_MAX_THREADS() : 1), 0, 0);
#undef CODE_BLOCK
		} else {
#define CODE_BLOCK							\
			for (int i = 0; i < n; i++) {			\
				int k = M->ia[i];			\
				Qfunc(thread_id, i, -1, &(M->a[k]), Qfunc_arg);	\
			}

			RUN_CODE_BLOCK((GMRFLib_Qx_strategy ? GMRFLib_MAX_THREADS() : 1), 0, 0);
#undef CODE_BLOCK
		}
	}

	for (i = 0; i < n + 1; i++) {
		M->ia1[i] = M->ia[i] + 1;
	}
	for (i = 0; i < na; i++) {
		M->ja1[i] = M->ja[i] + 1;
	}

	for (i = 0; i < na; i++) {
		GMRFLib_STOP_IF_NAN_OR_INF(M->a[i], i, -1);
		if (nan_error) {
			break;
		}
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
	GMRFLib_io_write(io, (const void *) (csr->ia), sizeof(int) * (csr->n + 1));
	GMRFLib_io_write(io, (const void *) (csr->ja), sizeof(int) * csr->na);
	GMRFLib_io_write(io, (const void *) (csr->ia1), sizeof(int) * (csr->n + 1));
	GMRFLib_io_write(io, (const void *) (csr->ja1), sizeof(int) * csr->na);
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

	int len = M->n + 1 + M->na;
	M->iwork = Calloc(2 * len, int);
	M->ia = M->iwork;
	GMRFLib_io_read(io, (void *) (M->ia), sizeof(int) * (M->n + 1));
	M->ja = M->iwork + M->n + 1;
	GMRFLib_io_read(io, (void *) (M->ja), sizeof(int) * M->na);

	M->ia1 = M->iwork + len;
	GMRFLib_io_read(io, (void *) (M->ia1), sizeof(int) * (M->n + 1));
	M->ja1 = M->iwork + len + M->n + 1;
	GMRFLib_io_read(io, (void *) (M->ja1), sizeof(int) * M->na);

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

	fprintf(fp, "Q->n = %1d, Q->na = %1d copy_only= %1d \n", csr->n, csr->na, csr->copy_only);
	for (i = k = 0; i < csr->n; i++) {
		nnb = csr->ia[i + 1] - csr->ia[i];
		for (jj = 0; jj < nnb; jj++) {
			j = csr->ja[k];
			value = csr->a[k];
			fprintf(fp, "%sQ[ %1d , %1d ] = %.12f\n", (jj > 0 ? "\t" : ""), i, j, value);
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
			jarr[k] = j;
			arr[k] = csr->a[k];
			k++;
		}
	}
	assert(k == csr->na);
	GMRFLib_tabulate_Qfunc_from_list(Qtab, graph, csr->na, iarr, jarr, arr, csr->n, NULL);

	Free(iarr);
	Free(jarr);
	Free(arr);

	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_init(GMRFLib_pardiso_store_tp ** store)
{
	int error = 0;
	int inla_ncpu(void);				       /* external function */
	GMRFLib_pardiso_store_tp *s = Calloc(1, GMRFLib_pardiso_store_tp);

	if (S.s_verbose) {
		PP("_pardiso_init()", s);
	}

	s->maxfct = 1;
	s->pstore = Calloc(GMRFLib_MAX_THREADS(), GMRFLib_pardiso_store_pr_thread_tp *);
	for (int i = 0; i < GMRFLib_MAX_THREADS(); i++) {
		s->pstore[i] = Calloc(1, GMRFLib_pardiso_store_pr_thread_tp);
	}
	assert(S.mtype == -2 || S.mtype == 2);
	s->mtype = S.mtype;
	s->msglvl = S.msglvl;
	s->solver = 0;
	s->iparm_default = Calloc(GMRFLib_PARDISO_PLEN, int);
	s->dparm_default = Calloc(GMRFLib_PARDISO_PLEN, double);
	s->iparm_default[0] = 0;			       /* use default values */
	s->iparm_default[2] = GMRFLib_PARDISO_MAX_NUM_THREADS();

	if (S.s_verbose) {
		PPg("_pardiso_init(): num_threads", (double) (s->iparm_default[2]));
	}

	pardisoinit(s->pt, &(s->mtype), &(s->solver), s->iparm_default, s->dparm_default, &error);

	s->iparm_default[1] = 3;			       /* use METIS5 */
	s->iparm_default[4] = 0;			       /* use internal reordering */
	s->iparm_default[7] = 0;			       /* maximum number of refinement steps */
	s->iparm_default[10] = 0;			       /* These are the default, but... */
	s->iparm_default[12] = 0;			       /* I need these for the divided LDL^Tx=b solver to work */
	s->iparm_default[20] = 0;			       /* Diagonal pivoting, and... */
	s->iparm_default[23] = 1;			       /* two level scheduling, and... */
	s->iparm_default[24] = 0;			       /* use parallel solve? */
	s->iparm_default[27] = S.parallel_reordering;	       /* parallel reordering? */
	s->iparm_default[33] = 1;			       /* want identical solutions */

	// options for METIS5; see manual. Pays off to do a good reordering
	s->iparm_default[57] = 2;			       /* default 1 */
	s->iparm_default[60] = 200;			       /* default */
	s->iparm_default[61] = 5;			       /* default 1 */

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

	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_setparam(GMRFLib_pardiso_flag_tp flag, GMRFLib_pardiso_store_tp * store, int *thread_num)
{
	int tnum = (thread_num ? *thread_num : omp_get_thread_num());

	assert(store->done_with_init == GMRFLib_TRUE);
	Memcpy((void *) (store->pstore[tnum]->iparm), (void *) (store->iparm_default), GMRFLib_PARDISO_PLEN * sizeof(int));
	Memcpy((void *) (store->pstore[tnum]->dparm), (void *) (store->dparm_default), GMRFLib_PARDISO_PLEN * sizeof(double));

	store->pstore[tnum]->nrhs = 0;
	store->pstore[tnum]->err_code = 0;
	store->pstore[tnum]->iparm[2] = store->pstore[GMRFLib_PSTORE_TNUM_REF]->iparm[2] = GMRFLib_PARDISO_MAX_NUM_THREADS();

	switch (flag) {
	case GMRFLib_PARDISO_FLAG_REORDER:
	case GMRFLib_PARDISO_FLAG_SYMFACT:
		store->pstore[tnum]->phase = 11;	       // analysis
		store->pstore[tnum]->iparm[4] = 0;	       /* 0 = compute the permutation */
		store->pstore[tnum]->iparm[39] = 1;	       /* 1 = return the permutation */
		store->pstore[tnum]->nrhs = S.nrhs_max;	       /* this is how it is, apparently */
		break;

	case GMRFLib_PARDISO_FLAG_CHOL:
		store->pstore[tnum]->phase = 22;	       // numerical factorization
		store->pstore[tnum]->iparm[32] = 1;	       /* determinant */
		store->pstore[tnum]->iparm[39] = 1;	       /* 1 = return the permutation (does not do that anymore) */
		break;

	case GMRFLib_PARDISO_FLAG_QINV:
		store->pstore[tnum]->phase = -22;
		store->pstore[tnum]->iparm[35] = 1;	       /* do not overwrite internal factor L with selected inversion */
		store->pstore[tnum]->iparm[36] = 0;	       /* return upper triangular Qinv */
		break;

	case GMRFLib_PARDISO_FLAG_SOLVE_L:
		store->pstore[tnum]->phase = 33;	       // solve
		store->pstore[tnum]->iparm[25] = (S.mtype == 2 ? 1 : -12);
		break;

	case GMRFLib_PARDISO_FLAG_SOLVE_LT:
		store->pstore[tnum]->phase = 33;	       // solve
		store->pstore[tnum]->iparm[25] = (S.mtype == 2 ? 2 : -23);
		break;

	case GMRFLib_PARDISO_FLAG_SOLVE_LLT:
		store->pstore[tnum]->phase = 33;	       // solve
		store->pstore[tnum]->iparm[25] = 0;
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

double GMRFLib_pardiso_Qfunc_default(int UNUSED(thread_id), int i, int j, double *UNUSED(values), void *arg)
{
	if (j < 0) {
		return NAN;
	}

	GMRFLib_graph_tp *g = (GMRFLib_graph_tp *) arg;
	return (i == j ? g->n + 2.0 * g->nnbs[i] : -1.0);
}

int GMRFLib_pardiso_reorder(GMRFLib_pardiso_store_tp * store, GMRFLib_graph_tp * graph)
{
	int tnum = omp_get_thread_num();
	int debug = S.debug;

	assert(store != NULL);
	assert(store->done_with_init == GMRFLib_TRUE);

	if (store->done_with_reorder == GMRFLib_TRUE) {
		return GMRFLib_SUCCESS;
	}

	int i, n, mnum1 = 1;
	GMRFLib_csr_tp *Q = NULL;

	GMRFLib_graph_duplicate(&(store->graph), graph);
	GMRFLib_pardiso_setparam(GMRFLib_PARDISO_FLAG_REORDER, store, NULL);
	GMRFLib_Q2csr(0, &Q, store->graph, GMRFLib_pardiso_Qfunc_default, (void *) store->graph);
	assert(Q);

	if (S.csr_check) {
		GMRFLib_csr_check(Q);
	}
	if (0 && S.debug) {
		GMRFLib_csr_print(stdout, Q);
	}


	n = Q->n;
	store->pstore[GMRFLib_PSTORE_TNUM_REF]->perm = Calloc(n, int);
	store->pstore[GMRFLib_PSTORE_TNUM_REF]->iperm = Calloc(n, int);

	if (S.parallel_reordering) {
		if (GMRFLib_openmp->adaptive && omp_get_level() == 0) {
			// this is the exception of the rule, as we want to run this in parallel if we are in adaptive model and
			// level=0.
			omp_set_num_threads(store->pstore[GMRFLib_PSTORE_TNUM_REF]->iparm[2]);
		} else {
			omp_set_num_threads(GMRFLib_openmp->max_threads_inner);
			assert(GMRFLib_openmp->max_threads_inner <= store->pstore[GMRFLib_PSTORE_TNUM_REF]->iparm[2]);
		}
	}

	pardiso(store->pt, &(store->maxfct), &mnum1, &(store->mtype),
		&(store->pstore[tnum]->phase),
		&(Q->n), Q->a, Q->ia1, Q->ja1, store->pstore[GMRFLib_PSTORE_TNUM_REF]->perm,
		&(store->pstore[tnum]->nrhs), store->pstore[tnum]->iparm,
		&(store->msglvl), &(store->pstore[tnum]->dummy), &(store->pstore[tnum]->dummy), &(store->pstore[tnum]->err_code),
		store->pstore[tnum]->dparm);

	// Just fill it with a dummy (identity) reordering
	for (i = 0; i < n; i++) {
		store->pstore[GMRFLib_PSTORE_TNUM_REF]->perm[i] = i;
		store->pstore[GMRFLib_PSTORE_TNUM_REF]->iperm[store->pstore[GMRFLib_PSTORE_TNUM_REF]->perm[i]] = i;
	}

	if (0 && debug) {
		for (i = 0; i < n; i++) {
			printf("perm[%1d] = %1d | iperm[%1d] = %1d\n", i, store->pstore[GMRFLib_PSTORE_TNUM_REF]->perm[i], i,
			       store->pstore[GMRFLib_PSTORE_TNUM_REF]->iperm[i]);
		}
	}

	if (store->pstore[tnum]->err_code) {
		GMRFLib_ERROR(GMRFLib_EPARDISO_INTERNAL_ERROR);
	}

	store->done_with_reorder = GMRFLib_TRUE;
	store->pstore[tnum]->L_nnz = store->pstore[tnum]->iparm[17] - 1;
	GMRFLib_csr_free(&Q);

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
	int n, *permutation = NULL;
	double *xx = NULL;

	assert(store);
	assert(store->graph);
	n = store->graph->n;
	xx = Calloc(n * m, double);
	Memcpy(xx, x, n * m * sizeof(double));
	permutation = (direction ? store->pstore[GMRFLib_PSTORE_TNUM_REF]->perm : store->pstore[GMRFLib_PSTORE_TNUM_REF]->iperm);
	assert(permutation);
	assert(m > 0);


#define CODE_BLOCK						\
	for (int j = 0; j < m; j++) {				\
		int k = j * n;					\
		for (int i = 0; i < n; i++) {			\
			x[k + i] = xx[k + permutation[i]];	\
		}						\
	}

	RUN_CODE_BLOCK((m > 8 ? m : 1), 0, 0);
#undef CODE_BLOCK

	Free(xx);
	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_symfact(GMRFLib_pardiso_store_tp * store)
{
	assert(store->done_with_init == GMRFLib_TRUE);
	assert(store->done_with_reorder == GMRFLib_TRUE);

	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_build(int thread_id, GMRFLib_pardiso_store_tp * store, GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg)
{
	assert(store != NULL);
	assert(store->done_with_init == GMRFLib_TRUE);
	assert(store->done_with_reorder == GMRFLib_TRUE);

	if (store->pstore[GMRFLib_PSTORE_TNUM_REF]->done_with_build == GMRFLib_TRUE && store->pstore[GMRFLib_PSTORE_TNUM_REF]->Q) {
		GMRFLib_csr_free(&(store->pstore[GMRFLib_PSTORE_TNUM_REF]->Q));
	}
	GMRFLib_Q2csr(thread_id, &(store->pstore[GMRFLib_PSTORE_TNUM_REF]->Q), graph, Qfunc, (void *) Qfunc_arg);

	if (S.csr_check) {
		GMRFLib_csr_check(store->pstore[GMRFLib_PSTORE_TNUM_REF]->Q);
	}
	if (0 && S.debug) {
		GMRFLib_csr_print(stdout, store->pstore[GMRFLib_PSTORE_TNUM_REF]->Q);
	}

	store->pstore[GMRFLib_PSTORE_TNUM_REF]->done_with_build = GMRFLib_TRUE;

	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_chol(GMRFLib_pardiso_store_tp * store)
{
	int tnum = omp_get_thread_num();
	int debug = S.debug;

	assert(store->done_with_init == GMRFLib_TRUE);
	assert(store->done_with_reorder == GMRFLib_TRUE);
	assert(store->pstore[GMRFLib_PSTORE_TNUM_REF]->done_with_build == GMRFLib_TRUE);
	assert(store->pstore[GMRFLib_PSTORE_TNUM_REF]->Q != NULL);

	int mnum1 = 1, n = store->pstore[GMRFLib_PSTORE_TNUM_REF]->Q->n, i;
	GMRFLib_pardiso_setparam(GMRFLib_PARDISO_FLAG_CHOL, store, NULL);

	if (debug) {
		printf("CHOL: level %d NUM_THREADS %d iparm[2] %d\n", omp_get_level(), omp_get_num_threads(),
		       store->pstore[GMRFLib_PSTORE_TNUM_REF]->iparm[2]);
	}

	if (GMRFLib_openmp->adaptive && omp_get_level() == 0) {
		// this is the exception of the rule, as we want to run this in parallel if we are in adaptive model and
		// level=0.
		omp_set_num_threads(store->pstore[GMRFLib_PSTORE_TNUM_REF]->iparm[2]);
	} else {
		omp_set_num_threads(GMRFLib_openmp->max_threads_inner);
		assert(GMRFLib_openmp->max_threads_inner <= store->pstore[GMRFLib_PSTORE_TNUM_REF]->iparm[2]);
	}

	GMRFLib_csr_tp *Q = store->pstore[GMRFLib_PSTORE_TNUM_REF]->Q;
	pardiso(store->pt, &(store->maxfct), &mnum1, &(store->mtype), &(store->pstore[tnum]->phase),
		&n, Q->a, Q->ia1, Q->ja1, store->pstore[GMRFLib_PSTORE_TNUM_REF]->perm, &(store->pstore[tnum]->nrhs),
		store->pstore[tnum]->iparm, &(store->msglvl), NULL, NULL, &(store->pstore[tnum]->err_code), store->pstore[tnum]->dparm);

	if (debug) {
		printf("Average number of non-zeros in L per row %.2f\n", store->pstore[tnum]->iparm[17] / (double) n);
	}
	// Revert back to C indexing ?
	if (GMRFLib_imin_value(store->pstore[GMRFLib_PSTORE_TNUM_REF]->perm, n, NULL) == 1) {
		for (i = 0; i < n; i++) {
			store->pstore[GMRFLib_PSTORE_TNUM_REF]->perm[i]--;
		}
	}
	for (i = 0; i < n; i++) {
		store->pstore[GMRFLib_PSTORE_TNUM_REF]->iperm[store->pstore[GMRFLib_PSTORE_TNUM_REF]->perm[i]] = i;
	}

	if (0 && debug) {
		for (i = 0; i < n; i++) {
			printf("perm[%1d] = %1d | iperm[%1d] = %1d\n", i, store->pstore[GMRFLib_PSTORE_TNUM_REF]->perm[i], i,
			       store->pstore[GMRFLib_PSTORE_TNUM_REF]->iperm[i]);
		}
	}

	if (store->pstore[tnum]->err_code != 0 || store->pstore[tnum]->iparm[22] > 0) {
		if (store->pstore[tnum]->iparm[22] > 0) {
			printf("\n*** PARDISO ERROR(%1d): not pos.def matrix: %1d eigenvalues are negative.\n",
			       store->pstore[tnum]->err_code, store->pstore[tnum]->iparm[22]);
		} else {
			printf("\n*** PARDISO ERROR(%1d): check manual.\n", store->pstore[tnum]->err_code);
		}
		printf("*** PARDISO ERROR: I will try to work around the problem...\n\n");
		return GMRFLib_EPOSDEF;
	}

	store->pstore[tnum]->log_det_Q = store->pstore[tnum]->dparm[32];
	store->pstore[GMRFLib_PSTORE_TNUM_REF]->done_with_chol = GMRFLib_TRUE;

	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_solve_core(GMRFLib_pardiso_store_tp * store, GMRFLib_pardiso_flag_tp flag, double *x, double *b, int nrhs)
{
	// note that 'x' = 'b' !!!!!!!!!!!!!!!

	assert(store);
	assert(nrhs > 0);
	assert(store->done_with_init == GMRFLib_TRUE);
	assert(store->done_with_reorder == GMRFLib_TRUE);
	assert(store->pstore[GMRFLib_PSTORE_TNUM_REF]->done_with_build == GMRFLib_TRUE);
	assert(store->pstore[GMRFLib_PSTORE_TNUM_REF]->done_with_chol == GMRFLib_TRUE);

	// this is so that the RHS can be overwritten
	int n = store->graph->n, mnum1 = 1, nblock, block_nrhs, err_code = 0, debug = 0;

	GMRFLib_pardiso_setparam(flag, store, NULL);
	int need_workaround = (GMRFLib_openmp->max_threads_inner > 1) && (store->pstore[omp_get_thread_num()]->iparm[7] > 0);
	assert(need_workaround == 0);			       /* disable this for the moment */

	int max_nt = GMRFLib_MAX_THREADS();
	double **ddparm = NULL;
	int **iiparm = NULL;
	int *factor = NULL;
	int *init_done = NULL;
	void **ppt = NULL;

	if (need_workaround) {
		size_t siz = sizeof(void *);		       /* to prevent compiler warning */
		ppt = Calloc(max_nt, void *);
		init_done = Calloc(max_nt, int);
		iiparm = Calloc(max_nt, int *);
		ddparm = Calloc(max_nt, double *);
		factor = Calloc(max_nt, int);
		ppt[0] = Calloc(max_nt * GMRFLib_PARDISO_PLEN, void *);
		iiparm[0] = Calloc(max_nt * GMRFLib_PARDISO_PLEN, int);
		ddparm[0] = Calloc(max_nt * GMRFLib_PARDISO_PLEN, double);
		for (int i = 1; i < max_nt; i++) {
			ppt[i] = ppt[i - 1] + GMRFLib_PARDISO_PLEN * siz;
			iiparm[i] = iiparm[i - 1] + GMRFLib_PARDISO_PLEN;
			ddparm[i] = ddparm[i - 1] + GMRFLib_PARDISO_PLEN;
			factor[i] = 0;
		}
	}

	int nt = 1;
	int nsolve;
	div_t d;

	if (nrhs > 1) {
		if (GMRFLib_openmp->adaptive && omp_get_level() == 0) {
			// this is the exception of the rule, as we want to run this in parallel if we are in adaptive model and
			// level=0.
			nt = store->pstore[GMRFLib_PSTORE_TNUM_REF]->iparm[2];
		} else {
			nt = GMRFLib_openmp->max_threads_inner;
			assert(GMRFLib_openmp->max_threads_inner <= store->pstore[GMRFLib_PSTORE_TNUM_REF]->iparm[2]);
		}
	}
	// printf("possible nt %d S.nrhs_max %d\n", nt, S.nrhs_max);

	if (nrhs == 1) {
		// in this case we always set block_nrhs=1 and nt=1
		block_nrhs = 1;
		nt = 1;
	} else if (S.nrhs_max < 0) {
		// this is the adaptive choice
		d = div(nrhs, nt);
		if (nrhs >= nt) {
			// if this is true, we need to divide the work between the nt cores
			block_nrhs = d.quot + (d.rem != 0);
		} else {
			// else we use one core pr rhs
			block_nrhs = 1;
			nt = IMIN(nt, nrhs);
		}
	} else if (nrhs <= S.nrhs_max) {
		// then we do all of them in one block
		block_nrhs = nrhs;
		nt = 1;
	} else {
		// S.nrhs_max define the max nrhs, so then we divide the work
		block_nrhs = S.nrhs_max;
		d = div(nrhs, block_nrhs);
		nt = IMIN(nt, d.quot + (d.rem != 0));
	}

	if (nt > 1) {
		omp_set_num_threads(nt);
	}

	d = div(nrhs, block_nrhs);
	nblock = d.quot;
	nsolve = nblock + (d.rem != 0);

	// printf("nrhs %d max_nrhs %d nt %d nsolve %d\n", nrhs, max_nrhs, nt, nsolve);

	GMRFLib_csr_tp *Q = store->pstore[GMRFLib_PSTORE_TNUM_REF]->Q;
#define CODE_BLOCK							\
	for (int i = 0; i < nsolve; i++) {				\
		CODE_BLOCK_ALL_WORK_ZERO();				\
		int idum = 0;						\
		int tnum = omp_get_thread_num();			\
		int offset = i * n * block_nrhs;			\
		int local_nrhs = (i < nblock ? block_nrhs : (int) d.rem); \
		double *bb = CODE_BLOCK_WORK_PTR(0);			\
									\
		GMRFLib_pardiso_setparam(flag, store, &tnum);		\
		Memcpy((void *) bb, (void *) (b + offset), n * local_nrhs * sizeof(double));\
		if (need_workaround) {					\
			if (!init_done[tnum]) {				\
				init_done[tnum] = 1;			\
				pardiso_copy_symbolic_factor_single(store->pt, ppt[tnum], \
								    store->pstore[tnum]->iparm, iiparm[tnum], \
								    store->pstore[tnum]->dparm, ddparm[tnum], \
								    &n, &(store->maxfct), &(factor[tnum])); \
			}						\
			pardiso(ppt[tnum], &(store->maxfct), &mnum1, &(store->mtype), &(store->pstore[tnum]->phase), \
				&n, Q->a, Q->ia1, Q->ja1, &idum, &local_nrhs, iiparm[tnum], &(store->msglvl), \
				bb, x + offset, &(store->pstore[tnum]->err_code), ddparm[tnum]); \
		} else {						\
			pardiso(store->pt, &(store->maxfct), &mnum1, &(store->mtype), &(store->pstore[tnum]->phase), \
				&n, Q->a, Q->ia1, Q->ja1, &idum, &local_nrhs, store->pstore[tnum]->iparm, &(store->msglvl), \
				bb, x + offset, &(store->pstore[tnum]->err_code), store->pstore[tnum]->dparm); \
		}							\
									\
		if (store->pstore[tnum]->err_code != 0) {		\
			err_code = GMRFLib_EPARDISO_INTERNAL_ERROR;	\
		}							\
									\
		if (debug) {						\
			double *yy = Calloc(local_nrhs * n,  double);	\
			for (int j = 0; j < local_nrhs; j++) {		\
				double normb, normr;			\
				pardiso_residual(&(store->mtype), &n,	\
						 Q->a, Q->ia1, Q->ja1, bb + j * n, x + offset + j * n, yy + j * n, &normb, &normr); \
				printf("ni j %d %d The norm of the residual is %e \n ", i, j, normr / normb); \
			}						\
			Free(yy);					\
		}							\
	}

	RUN_CODE_BLOCK(nt, 1, block_nrhs * n);
#undef CODE_BLOCK

	if (need_workaround) {
		for (int i = 0; i < max_nt; i++) {
			if (ppt[i]) {
				pardiso_delete_symbolic_factor_single(ppt[i], iiparm[i], &(factor[i]));
			}
		}
		Free(iiparm[0]);
		Free(init_done);
		Free(iiparm);
		Free(ddparm[0]);
		Free(ddparm);
		Free(ppt);
		Free(factor);
	}

	if (err_code) {
		GMRFLib_ERROR(err_code);
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_solve_L(GMRFLib_pardiso_store_tp * store, double *x, double *b, int nrhs)
{
	int res = GMRFLib_pardiso_solve_core(store, GMRFLib_PARDISO_FLAG_SOLVE_L, x, b, nrhs);
	GMRFLib_pardiso_iperm(x, nrhs, store);

	return res;
}

int GMRFLib_pardiso_solve_LT(GMRFLib_pardiso_store_tp * store, double *x, double *b, int nrhs)
{
	GMRFLib_pardiso_perm(b, nrhs, store);
	int res = GMRFLib_pardiso_solve_core(store, GMRFLib_PARDISO_FLAG_SOLVE_LT, x, b, nrhs);
	if (x != b) {
		GMRFLib_pardiso_iperm(b, nrhs, store);
	}

	return res;
}

int GMRFLib_pardiso_solve_LLT(GMRFLib_pardiso_store_tp * store, double *x, double *b, int nrhs)
{
	int res = GMRFLib_pardiso_solve_core(store, GMRFLib_PARDISO_FLAG_SOLVE_LLT, x, b, nrhs);
	return res;
}

double GMRFLib_pardiso_logdet(GMRFLib_pardiso_store_tp * store)
{
	int tnum = omp_get_thread_num();
	return (store->pstore[tnum]->log_det_Q);
}

int GMRFLib_pardiso_bitmap(void)
{
	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_Qinv_INLA(GMRFLib_problem_tp * problem)
{
	if (problem == NULL) {
		return GMRFLib_SUCCESS;
	}

	GMRFLib_pardiso_Qinv(problem->sub_sm_fact.PARDISO_fact);

	GMRFLib_csr_tp *Qi = problem->sub_sm_fact.PARDISO_fact->pstore[GMRFLib_PSTORE_TNUM_REF]->Qinv;
	int n = Qi->n, i, j, jj, k;
	map_id **Qinv = Calloc(n, map_id *);

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
#define CODE_BLOCK							\
		for (int i = 0; i < n; i++) {				\
			for (int k = -1; (k = (int) map_id_next(Qinv[i], k)) != -1;) { \
				double value = 0.0;			\
				int j = Qinv[i]->contents[k].key;	\
				map_id_get(Qinv[i], j, &value);		\
				for (int kk = 0; kk < problem->sub_constr->nc; kk++) { \
					value -= problem->constr_m[i + kk * n] * problem->qi_at_m[j + kk * n]; \
				}					\
				map_id_set(Qinv[i], j, value);		\
			}						\
		}

		RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 0, 0);
#undef CODE_BLOCK
	}

	GMRFLib_Qinv_tp *subQinv = Calloc(1, GMRFLib_Qinv_tp);

	subQinv->Qinv = Qinv;
	subQinv->mapping = Calloc(n, int);
	for (i = 0; i < n; i++) {
		subQinv->mapping[i] = i;
	}
	problem->sub_inverse = subQinv;

	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_Qinv(GMRFLib_pardiso_store_tp * store)
{
	int tnum = omp_get_thread_num();
	assert(store->done_with_reorder == GMRFLib_TRUE);
	assert(store->pstore[GMRFLib_PSTORE_TNUM_REF]->done_with_build == GMRFLib_TRUE);
	assert(store->pstore[GMRFLib_PSTORE_TNUM_REF]->done_with_chol == GMRFLib_TRUE);

	if (store->pstore[GMRFLib_PSTORE_TNUM_REF]->Qinv) {
		GMRFLib_csr_free(&(store->pstore[GMRFLib_PSTORE_TNUM_REF]->Qinv));
	}

	GMRFLib_csr_duplicate(&(store->pstore[GMRFLib_PSTORE_TNUM_REF]->Qinv), store->pstore[GMRFLib_PSTORE_TNUM_REF]->Q);
	GMRFLib_pardiso_setparam(GMRFLib_PARDISO_FLAG_QINV, store, NULL);
	int mnum1 = 1;

	if (GMRFLib_openmp->adaptive && omp_get_level() == 0) {
		// this is the exception of the rule, as we want to run this in parallel if we are in adaptive model and
		// level=0.
		omp_set_num_threads(store->pstore[GMRFLib_PSTORE_TNUM_REF]->iparm[2]);
	} else {
		omp_set_num_threads(GMRFLib_openmp->max_threads_inner);
		assert(GMRFLib_openmp->max_threads_inner <= store->pstore[GMRFLib_PSTORE_TNUM_REF]->iparm[2]);
	}
	GMRFLib_csr_tp *Qinv = store->pstore[GMRFLib_PSTORE_TNUM_REF]->Qinv;
	pardiso(store->pt, &(store->maxfct), &mnum1, &(store->mtype), &(store->pstore[tnum]->phase),
		&(Qinv->n), Qinv->a, Qinv->ia1, Qinv->ja1, NULL,
		&(store->pstore[tnum]->nrhs), store->pstore[tnum]->iparm, &(store->msglvl),
		NULL, NULL, &(store->pstore[tnum]->err_code), store->pstore[tnum]->dparm);

	if (store->pstore[tnum]->err_code != 0) {
		GMRFLib_ERROR(GMRFLib_EPARDISO_INTERNAL_ERROR);
	}

	if (0 && S.debug) {
		GMRFLib_csr_print(stdout, Qinv);
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_exit(void)
{
	if (S.static_pstores != NULL) {
		for (int i = 0; i < PSTORES_NUM; i++) {
			if (S.static_pstores[i]) {
				GMRFLib_pardiso_free(&(S.static_pstores[i]));
			}
		}
	}
	Free(S.busy);
	Free(S.static_pstores);
	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_free(GMRFLib_pardiso_store_tp ** store)
{
	int tnum = omp_get_thread_num();
	if (store == NULL || *store == NULL) {
		return GMRFLib_SUCCESS;
	}

	if ((*store)->copy_pardiso_ptr) {
		// this is special
		if (S.s_verbose) {
			FIXME("Free pardiso store with copy_pardiso_ptr = 1");
		}
		if ((*store)->pstore) {
			for (int i = 0; i < GMRFLib_MAX_THREADS(); i++) {
				if ((*store)->pstore[i]) {
					Free((*store)->pstore[i]->perm);
					Free((*store)->pstore[i]->iperm);
					Free((*store)->pstore[i]);
				}
			}
		}
		Free((*store)->pstore);
		Free((*store));
		*store = NULL;

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

		if ((*store)->pstore[tnum]) {
			(*store)->pstore[tnum]->phase = -1;
			int mnum1 = 1;
			pardiso((*store)->pt, &((*store)->maxfct), &mnum1, &((*store)->mtype),
				&((*store)->pstore[tnum]->phase),
				&((*store)->pstore[tnum]->idummy),
				&((*store)->pstore[tnum]->dummy), &((*store)->pstore[tnum]->idummy),
				&((*store)->pstore[tnum]->idummy),
				&((*store)->pstore[tnum]->idummy),
				&((*store)->pstore[tnum]->nrhs), (*store)->pstore[tnum]->iparm, &((*store)->msglvl),
				NULL, NULL, &((*store)->pstore[tnum]->err_code), (*store)->pstore[tnum]->dparm);

			GMRFLib_csr_free(&((*store)->pstore[GMRFLib_PSTORE_TNUM_REF]->Q));
			GMRFLib_csr_free(&((*store)->pstore[GMRFLib_PSTORE_TNUM_REF]->Qinv));
			Free((*store)->pstore[tnum]);
		}

		GMRFLib_graph_free((*store)->graph);
		Free((*store)->iparm_default);
		Free((*store)->dparm_default);
		Free(*store);
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_duplicate_pardiso_store(GMRFLib_pardiso_store_tp ** new, GMRFLib_pardiso_store_tp * old, int UNUSED(copy_ptr), int copy_pardiso_ptr)
{
	int tnum = omp_get_thread_num();
	// if copy_pardiso_ptr, then copy the ptr to read-only objects. 'copy_ptr' is NOT USED
	int debug = S.debug, failsafe_mode = 0;
	if (old == NULL) {
		*new = NULL;
		return GMRFLib_SUCCESS;
	}

	if (failsafe_mode) {
		// FIXME("-->duplicate by creating a new one each time");
		GMRFLib_pardiso_init(new);
		GMRFLib_pardiso_reorder(*new, old->graph);
		return GMRFLib_SUCCESS;
	}
	GMRFLib_ENTER_ROUTINE;

	if (copy_pardiso_ptr) {
		int i;

#define CP(_what) dup->_what = old->_what
#define CP2(_what) dup->pstore[tnum]->_what = old->pstore[tnum]->_what
#define CP2_ref(_what) dup->pstore[GMRFLib_PSTORE_TNUM_REF]->_what = old->pstore[GMRFLib_PSTORE_TNUM_REF]->_what
#define CPv_ref(_what, type, len)					\
		if (old->pstore[GMRFLib_PSTORE_TNUM_REF]->_what) {			\
			dup->pstore[GMRFLib_PSTORE_TNUM_REF]->_what = Calloc(len, type); \
			Memcpy((void *) (dup->pstore[GMRFLib_PSTORE_TNUM_REF]->_what), (void *) (old->pstore[GMRFLib_PSTORE_TNUM_REF]->_what), (len) * sizeof(type)); \
		} else {						\
			dup->pstore[GMRFLib_PSTORE_TNUM_REF]->_what = NULL;		\
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

		if (!dup->pstore)
			dup->pstore = Calloc(GMRFLib_MAX_THREADS(), GMRFLib_pardiso_store_pr_thread_tp *);
		if (!dup->pstore[tnum])
			dup->pstore[tnum] = Calloc(1, GMRFLib_pardiso_store_pr_thread_tp);
		if (!dup->pstore[GMRFLib_PSTORE_TNUM_REF])
			dup->pstore[GMRFLib_PSTORE_TNUM_REF] = Calloc(1, GMRFLib_pardiso_store_pr_thread_tp);
		for (i = 0; i < GMRFLib_PARDISO_PLEN; i++) {
			CP2_ref(iparm[i]);
			CP2_ref(dparm[i]);
		}
		CP2_ref(done_with_build);
		CP2_ref(done_with_chol);
		CP2(dummy);
		CP2(err_code);
		CP2(idummy);
		CP2(nrhs);
		CP2(phase);
		CP2(L_nnz);

		CPv_ref(perm, int, old->graph->n);
		CPv_ref(iperm, int, old->graph->n);

		CP2(log_det_Q);
		CP2_ref(Q);
		CP2_ref(Qinv);
#undef CP
#undef CP2
#undef CP2_ref
#undef CPv
#undef CPv_ref
		*new = dup;
		return GMRFLib_SUCCESS;
	}

	if (S.static_pstores == NULL) {
#pragma omp critical
		{
			if (S.static_pstores == NULL) {
				if (S.s_verbose) {
					printf("==> init static_pstores\n");
				}
				S.busy = Calloc(PSTORES_NUM, int);
				S.static_pstores = Calloc(PSTORES_NUM, GMRFLib_pardiso_store_tp *);
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
			printf("%s:%1d: static_pstores...iparm[2] = %1d\n", __FILE__, __LINE__,
			       S.static_pstores[idx]->pstore[GMRFLib_PSTORE_TNUM_REF]->iparm[2]);
			printf("%s:%1d: level %d max_threads_nested = %1d\n", __FILE__, __LINE__, omp_get_level(),
			       GMRFLib_openmp->max_threads_nested[1]);
		}
		ok = (S.static_pstores[idx]->pstore[GMRFLib_PSTORE_TNUM_REF]->iparm[2] >= GMRFLib_openmp->max_threads_nested[1]);
	} else {
		ok = 1;
	}
	if (!ok) {
		P(S.static_pstores[idx]->pstore[GMRFLib_PSTORE_TNUM_REF]->iparm[2]);
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
