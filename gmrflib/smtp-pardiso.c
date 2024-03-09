
/* smtp-pardiso.c
 * 
 * Copyright (C) 2018-2024 Havard Rue
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

#if defined(WIN32) || defined(WINDOWS)
#include <io.h>
#endif

#include "GMRFLib/hashP.h"
#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

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

static int csr_store_use = 1;
static map_strvp csr_store;
static int csr_store_must_init = 1;
static int csr_store_debug = 0;

int GMRFLib_csr_init_store(void)
{
	GMRFLib_ENTER_ROUTINE;
	csr_store_debug = GMRFLib_DEBUG_IF_TRUE();

	if (csr_store_use) {
		if (csr_store_must_init) {
			map_strvp_init_hint(&csr_store, 128);
			csr_store_must_init = 0;
			if (csr_store_debug) {
				printf("\tcsr_store: init storage\n");
			}
		}
	}
	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

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

int GMRFLib_csr_free(GMRFLib_csr_tp **csr)
{
	if (*csr) {
		// Free((*csr)->s->iwork);
		if ((*csr)->copy_only) {
			// do not free
		} else {
			Free((*csr)->a);
		}
		Free(*csr);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_csr_duplicate(GMRFLib_csr_tp **csr_to, GMRFLib_csr_tp *csr_from, int skeleton)
{
	// skeleton: if TRUE, point to the copy of 'a' if its a copy, otherwise alloc. if FALSE, always alloc.
	if (csr_from == NULL) {
		*csr_to = NULL;
		return GMRFLib_SUCCESS;
	}

	*csr_to = Calloc(1, GMRFLib_csr_tp);
	void **p = NULL;

	if (csr_store_use && csr_from->s->sha) {
		p = map_strvp_ptr(&csr_store, (char *) csr_from->s->sha);
		if (csr_store_debug) {
			if (p) {
				printf("\t[%1d] csr_duplicate: crs is found in store: do not duplicate\n", omp_get_thread_num());
			} else {
				printf("\t[%1d] csr_duplicate: crs is not found in store: duplicate\n", omp_get_thread_num());
			}
		}
	} else {
		if (csr_store_debug) {
			printf("\t[%1d] csr_duplicate: no sha\n", omp_get_thread_num());
		}
	}

	if (p) {
		(*csr_to)->s = (GMRFLib_csr_skeleton_tp *) * p;
	} else {
		int n = csr_from->s->n;
		int na = csr_from->s->na;
		int n1 = GMRFLib_align(n + 1, sizeof(int));
		int len = n1 + na;
		int llen = GMRFLib_align(len, sizeof(int));

		(*csr_to)->copy_only = csr_from->copy_only;
		(*csr_to)->s = Calloc(1, GMRFLib_csr_skeleton_tp);
		(*csr_to)->s->sha = NULL;
		(*csr_to)->s = Calloc(1, GMRFLib_csr_skeleton_tp);
		(*csr_to)->s->n = n;
		(*csr_to)->s->na = na;
		(*csr_to)->s->iwork = Calloc(2 * llen, int);
		(*csr_to)->s->ia = (*csr_to)->s->iwork;
		(*csr_to)->s->ja = (*csr_to)->s->iwork + n1;
		(*csr_to)->s->ia1 = (*csr_to)->s->iwork + llen;
		(*csr_to)->s->ja1 = (*csr_to)->s->iwork + llen + n1;
		Memcpy((void *) ((*csr_to)->s->iwork), (void *) (csr_from->s->iwork), (size_t) (2 * llen) * sizeof(int));
	}

	if (csr_from->a) {
		if (skeleton) {
			if (csr_from->copy_only) {
				(*csr_to)->a = csr_from->a;
			} else {
				(*csr_to)->a = Calloc(csr_from->s->na, double);
				Memcpy((void *) ((*csr_to)->a), (void *) (csr_from->a), (size_t) (csr_from->s->na) * sizeof(double));
			}
		} else {
			(*csr_to)->copy_only = 0;
			(*csr_to)->a = Calloc(csr_from->s->na, double);
			Memcpy((void *) ((*csr_to)->a), (void *) (csr_from->a), (size_t) (csr_from->s->na) * sizeof(double));
		}
	} else {
		assert(0 == 1);
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_csr_check(GMRFLib_csr_tp *M)
{
	int mtype = S.mtype, error = 0;

	assert(M);
	pardiso_chkmatrix(&mtype, &(M->s->n), M->a, M->s->ia1, M->s->ja1, &error);
	if (error != 0) {
		GMRFLib_ERROR(GMRFLib_ESNH);
	}
	return GMRFLib_SUCCESS;
}

GMRFLib_csr_skeleton_tp *GMRFLib_csr_skeleton(GMRFLib_graph_tp *graph)
{
	GMRFLib_csr_skeleton_tp *Ms = NULL;

	int n, na, len, n1, llen;

	if (csr_store_use && graph->sha) {
		void **p = NULL;
		p = map_strvp_ptr(&csr_store, (char *) graph->sha);
		if (csr_store_debug) {
			if (p) {
				printf("\t[%1d] csr_skeleton: crs is found in store\n", omp_get_thread_num());
			} else {
				printf("\t[%1d] csr_skeleton: crs is not found in store\n", omp_get_thread_num());
			}
		}
		if (p) {
			Ms = (GMRFLib_csr_skeleton_tp *) * p;
			return (Ms);
		}
	} else {
		if (csr_store_debug) {
			printf("\t[%1d] csr_skeleton: no sha\n", omp_get_thread_num());
		}
	}

	Ms = Calloc(1, GMRFLib_csr_skeleton_tp);
	if (graph->sha) {
		Ms->sha = Calloc(GMRFLib_SHA_DIGEST_LEN + 1, unsigned char);
		Memcpy(Ms->sha, graph->sha, GMRFLib_SHA_DIGEST_LEN + 1);
	}
	n = graph->n;
	na = graph->nnz / 2 + n;			       // only upper triangular. yes, integer division
	n1 = GMRFLib_align(n + 1, sizeof(int));
	len = n1 + na;
	llen = GMRFLib_align(len, sizeof(int));
	Ms->na = na;
	Ms->n = n;
	Ms->iwork = Calloc(2 * llen, int);
	Ms->ia = Ms->iwork;
	Ms->ja = Ms->iwork + n1;
	Ms->ia1 = Ms->iwork + llen;
	Ms->ja1 = Ms->iwork + llen + n1;

	// new code. by doing it in two steps we can do the second one in parallel, and this is the one that take time.
	int *k_arr = Calloc(n, int);
	Ms->ia[0] = 0;
	for (int i = 0, k = 0; i < n; i++) {
		Ms->ja[k++] = i;
		k_arr[i] = k;
		k += graph->lnnbs[i];
		Ms->ia[i + 1] = Ms->ia[i] + (1 + graph->lnnbs[i]);
	}
	assert(Ms->ia[n] == na);

#define CODE_BLOCK							\
	for (int i = 0; i < n; i++) {					\
		if (graph->lnnbs[i]) {					\
			int k = k_arr[i];				\
			Memcpy(&(Ms->ja[k]), graph->lnbs[i], graph->lnnbs[i] * sizeof(int)); \
		}							\
	}

	RUN_CODE_BLOCK(1, 0, 0);			       /* TO QUICK TO DO IN PARALLEL */
#undef CODE_BLOCK
	Free(k_arr);

#pragma GCC ivdep
	for (int i = 0; i < n + 1; i++) {
		Ms->ia1[i] = Ms->ia[i] + 1;
	}
#pragma GCC ivdep
	for (int i = 0; i < na; i++) {
		Ms->ja1[i] = Ms->ja[i] + 1;
	}

	if (csr_store_use && graph->sha) {
		if (csr_store_debug) {
			printf("\t[%1d] csr_store: store crs 0x%p\n", omp_get_thread_num(), (void *) Ms);
		}
#pragma omp critical (Name_488cde57983063f09dc9ddf7c473d77c79ea2929)
		{
			map_strvp_set(&csr_store, (char *) graph->sha, (void *) Ms);
		}
	}

	return (Ms);
}

int GMRFLib_Q2csr(int thread_id, GMRFLib_csr_tp **csr, GMRFLib_graph_tp *graph, GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_arg)
{
	GMRFLib_ENTER_ROUTINE;

#define M (*csr)

	M = Calloc(1, GMRFLib_csr_tp);
	M->s = GMRFLib_csr_skeleton(graph);

	// when this is true, we can just copy the pointer to the matrix.
	int used_fast_tab = 0;
	if (Qfunc == GMRFLib_tabulate_Qfunction_std) {
		GMRFLib_tabulate_Qfunc_arg_tp *arg = (GMRFLib_tabulate_Qfunc_arg_tp *) Qfunc_arg;
		if (arg->Q) {
			M->a = arg->Q->a;
			// mark this a copy only, not to be free'd.
			M->copy_only = 1;
			used_fast_tab = 1;
		}
	}

	if (!used_fast_tab) {
		M->a = Calloc(M->s->na, double);
		// a bit more manual work
		double val = Qfunc(thread_id, 0, -1, &(M->a[0]), Qfunc_arg);
		if (ISNAN(val)) {
#define CODE_BLOCK							\
			for (int i = 0; i < M->s->n; i++) {		\
				for (int k = M->s->ia[i]; k < M->s->ia[i + 1]; k++) { \
					int j = M->s->ja[k];		\
					M->a[k] = Qfunc(thread_id, i, j, NULL, Qfunc_arg); \
				}					\
			}

			RUN_CODE_BLOCK((GMRFLib_Qx_strategy ? GMRFLib_MAX_THREADS() : 1), 0, 0);
#undef CODE_BLOCK
		} else {
#define CODE_BLOCK							\
			for (int i = 0; i < M->s->n; i++) {		\
				int k = M->s->ia[i];			\
				Qfunc(thread_id, i, -1, &(M->a[k]), Qfunc_arg);	\
			}

			RUN_CODE_BLOCK((GMRFLib_Qx_strategy ? GMRFLib_MAX_THREADS() : 1), 0, 0);
#undef CODE_BLOCK
		}
	}

	int nan_error = 0;
	for (int i = 0; i < M->s->na; i++) {
		GMRFLib_STOP_IF_NAN_OR_INF(M->a[i], i, -1);
		if (nan_error) {
			GMRFLib_LEAVE_ROUTINE;
			return !GMRFLib_SUCCESS;
		}
	}

#undef M
	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_csr_write(char *filename, GMRFLib_csr_tp *csr)
{
	// write to file
	GMRFLib_io_tp *io = NULL;
	GMRFLib_io_open(&io, filename, "wb");
	GMRFLib_io_write(io, (const void *) &(csr->s->n), sizeof(int));
	GMRFLib_io_write(io, (const void *) &(csr->s->na), sizeof(int));
	GMRFLib_io_write(io, (const void *) (csr->s->ia), sizeof(int) * (csr->s->n + 1));
	GMRFLib_io_write(io, (const void *) (csr->s->ja), sizeof(int) * csr->s->na);
	GMRFLib_io_write(io, (const void *) (csr->s->ia1), sizeof(int) * (csr->s->n + 1));
	GMRFLib_io_write(io, (const void *) (csr->s->ja1), sizeof(int) * csr->s->na);
	GMRFLib_io_write(io, (const void *) (csr->a), sizeof(double) * csr->s->na);
	GMRFLib_io_close(io);

	return GMRFLib_SUCCESS;
}

int GMRFLib_csr_read(char *filename, GMRFLib_csr_tp **csr)
{
	// write to file
	GMRFLib_io_tp *io = NULL;
#define M (*csr)
	M = Calloc(1, GMRFLib_csr_tp);
	M->s = Calloc(1, GMRFLib_csr_skeleton_tp);

	GMRFLib_io_open(&io, filename, "rb");
	GMRFLib_io_read(io, (void *) &(M->s->n), sizeof(int));
	GMRFLib_io_read(io, (void *) &(M->s->na), sizeof(int));

	int len = M->s->n + 1 + M->s->na;
	M->s->iwork = Calloc(2 * len, int);
	M->s->ia = M->s->iwork;
	GMRFLib_io_read(io, (void *) (M->s->ia), sizeof(int) * (M->s->n + 1));
	M->s->ja = M->s->iwork + M->s->n + 1;
	GMRFLib_io_read(io, (void *) (M->s->ja), sizeof(int) * M->s->na);

	M->s->ia1 = M->s->iwork + len;
	GMRFLib_io_read(io, (void *) (M->s->ia1), sizeof(int) * (M->s->n + 1));
	M->s->ja1 = M->s->iwork + len + M->s->n + 1;
	GMRFLib_io_read(io, (void *) (M->s->ja1), sizeof(int) * M->s->na);

	M->a = Calloc(M->s->na, double);
	GMRFLib_io_read(io, (void *) (M->a), sizeof(double) * M->s->na);

	GMRFLib_io_close(io);
#undef M
	return GMRFLib_SUCCESS;
}

int GMRFLib_csr_print(FILE *fp, GMRFLib_csr_tp *csr)
{
	int i, j, k, jj, nnb;
	double value;

	fprintf(fp, "Q->s->n = %1d, Q->s->na = %1d copy_only= %1d \n", csr->s->n, csr->s->na, csr->copy_only);
	for (i = k = 0; i < csr->s->n; i++) {
		nnb = csr->s->ia[i + 1] - csr->s->ia[i];
		for (jj = 0; jj < nnb; jj++) {
			j = csr->s->ja[k];
			value = csr->a[k];
			fprintf(fp, "%sQ[ %1d , %1d ] = %.12f\n", (jj > 0 ? "\t" : ""), i, j, value);
			k++;
		}
	}
	fflush(fp);

	return GMRFLib_SUCCESS;
}

int GMRFLib_csr2Q(GMRFLib_tabulate_Qfunc_tp **Qtab, GMRFLib_graph_tp **graph, GMRFLib_csr_tp *csr)
{
	int i, j, k, jj, nnb;

	int *iarr = Calloc(csr->s->na, int);
	int *jarr = Calloc(csr->s->na, int);
	double *arr = Calloc(csr->s->na, double);

	for (i = k = 0; i < csr->s->n; i++) {
		nnb = csr->s->ia[i + 1] - csr->s->ia[i];
		for (jj = 0; jj < nnb; jj++) {
			j = csr->s->ja[k];
			iarr[k] = i;
			jarr[k] = j;
			arr[k] = csr->a[k];
			k++;
		}
	}
	assert(k == csr->s->na);
	GMRFLib_tabulate_Qfunc_from_list(Qtab, graph, csr->s->na, iarr, jarr, arr, csr->s->n, NULL);

	Free(iarr);
	Free(jarr);
	Free(arr);

	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_init(GMRFLib_pardiso_store_tp **store)
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
	s->iparm_default[53] = 0;			       /* number of expected dense columns (>0 gives memory failure) */

	// options for METIS5; see manual. Pays off to do a good reordering
	s->iparm_default[54] = 1;			       /* use options */
	s->iparm_default[55] = 3;			       /* default */
	s->iparm_default[56] = 1;			       /* default */
	s->iparm_default[57] = 1;			       /* !default (gives memory failure) */
	s->iparm_default[57] = 2;			       /* default */
	s->iparm_default[58] = 0;			       /* default */
	s->iparm_default[59] = 0;			       /* default */
	s->iparm_default[59] = 3;			       /* !default */
	s->iparm_default[60] = 200;			       /* default */
	s->iparm_default[60] = 150;			       /* !default */
	s->iparm_default[61] = 1;			       /* default */
	s->iparm_default[61] = 5;			       /* !default */

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

int GMRFLib_pardiso_setparam(GMRFLib_pardiso_flag_tp flag, GMRFLib_pardiso_store_tp *store, int *thread_num)
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
	{
		store->pstore[tnum]->phase = 11;	       // analysis
		store->pstore[tnum]->iparm[4] = 0;	       /* 0 = compute the permutation */
		store->pstore[tnum]->iparm[39] = 1;	       /* 1 = return the permutation */
		store->pstore[tnum]->nrhs = S.nrhs_max;	       /* this is how it is, apparently */
	}
		break;

	case GMRFLib_PARDISO_FLAG_CHOL:
	{
		store->pstore[tnum]->phase = 22;	       // numerical factorization
		store->pstore[tnum]->iparm[32] = 1;	       /* determinant */
		store->pstore[tnum]->iparm[39] = 1;	       /* 1 = return the permutation (does not do that anymore) */
	}
		break;

	case GMRFLib_PARDISO_FLAG_QINV:
	{
		store->pstore[tnum]->phase = -22;
		store->pstore[tnum]->iparm[35] = 1;	       /* do not overwrite internal factor L with selected inversion */
		store->pstore[tnum]->iparm[36] = 0;	       /* return upper triangular Qinv */
	}
		break;

	case GMRFLib_PARDISO_FLAG_SOLVE_L:
	{
		store->pstore[tnum]->phase = 33;	       // solve
		store->pstore[tnum]->iparm[25] = (S.mtype == 2 ? 1 : -12);
	}
		break;

	case GMRFLib_PARDISO_FLAG_SOLVE_LT:
	{
		store->pstore[tnum]->phase = 33;	       // solve
		store->pstore[tnum]->iparm[25] = (S.mtype == 2 ? 2 : -23);
	}
		break;

	case GMRFLib_PARDISO_FLAG_SOLVE_LLT:
	{
		store->pstore[tnum]->phase = 33;	       // solve
		store->pstore[tnum]->iparm[25] = 0;
	}
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

int GMRFLib_pardiso_reorder(GMRFLib_pardiso_store_tp *store, GMRFLib_graph_tp *graph)
{
	int tnum = omp_get_thread_num();
	int debug = S.debug;

	assert(store != NULL);
	assert(store->done_with_init == GMRFLib_TRUE);

	if (store->done_with_reorder == GMRFLib_TRUE) {
		return GMRFLib_SUCCESS;
	}

	int n, mnum1 = 1;
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

	assert(Q->s->n == store->graph->n);
	n = store->graph->n;
	store->pstore[GMRFLib_PSTORE_TNUM_REF]->perm = Calloc(n, int);
	store->pstore[GMRFLib_PSTORE_TNUM_REF]->iperm = Calloc(n, int);

	if (S.s_verbose) {
		printf("Init store perm and iperm with n = %1d\n", n);
	}

	if (S.parallel_reordering) {
		if (GMRFLib_openmp->adaptive && omp_get_level() == 0) {
			// this is the exception of the rule, as we want to run this in parallel if we are in adaptive model and
			// level=0.
			omp_set_num_threads(store->pstore[GMRFLib_PSTORE_TNUM_REF]->iparm[2]);
		} else {
			omp_set_num_threads(IMIN(GMRFLib_openmp->max_threads_inner, store->pstore[GMRFLib_PSTORE_TNUM_REF]->iparm[2]));
			// assert(GMRFLib_openmp->max_threads_inner <= store->pstore[GMRFLib_PSTORE_TNUM_REF]->iparm[2]);
		}
	}

	pardiso(store->pt, &(store->maxfct), &mnum1, &(store->mtype),
		&(store->pstore[tnum]->phase),
		&(Q->s->n), Q->a, Q->s->ia1, Q->s->ja1, store->pstore[GMRFLib_PSTORE_TNUM_REF]->perm,
		&(store->pstore[tnum]->nrhs), store->pstore[tnum]->iparm,
		&(store->msglvl), &(store->pstore[tnum]->dummy), &(store->pstore[tnum]->dummy), &(store->pstore[tnum]->err_code),
		store->pstore[tnum]->dparm);

	// Just fill it with a dummy (identity) reordering
	for (int i = 0; i < n; i++) {
		store->pstore[GMRFLib_PSTORE_TNUM_REF]->perm[i] = i;
		store->pstore[GMRFLib_PSTORE_TNUM_REF]->iperm[store->pstore[GMRFLib_PSTORE_TNUM_REF]->perm[i]] = i;
	}
	store->pstore[GMRFLib_PSTORE_TNUM_REF]->perm_identity = 1;

	if (0 && debug) {
		for (int i = 0; i < n; i++) {
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

int GMRFLib_pardiso_perm(double *x, int m, GMRFLib_pardiso_store_tp *store)
{
	return GMRFLib_pardiso_perm_core(x, m, store, 1);
}

int GMRFLib_pardiso_iperm(double *x, int m, GMRFLib_pardiso_store_tp *store)
{
	return GMRFLib_pardiso_perm_core(x, m, store, 0);
}

int GMRFLib_pardiso_perm_core(double *x, int m, GMRFLib_pardiso_store_tp *store, int direction)
{
	int n, *permutation = NULL;
	double *xx = NULL;

	assert(store);
	assert(store->graph);

	if (store->pstore[GMRFLib_PSTORE_TNUM_REF]->perm_identity) {
		return GMRFLib_SUCCESS;
	}

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

int GMRFLib_pardiso_symfact(GMRFLib_pardiso_store_tp *store)
{
	assert(store->done_with_init == GMRFLib_TRUE);
	assert(store->done_with_reorder == GMRFLib_TRUE);

	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_build(int thread_id, GMRFLib_pardiso_store_tp *store, GMRFLib_graph_tp *graph, GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_arg)
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

int GMRFLib_pardiso_chol(GMRFLib_pardiso_store_tp *store)
{
	int tnum = omp_get_thread_num();
	int debug = S.debug;

	assert(store->done_with_init == GMRFLib_TRUE);
	assert(store->done_with_reorder == GMRFLib_TRUE);
	assert(store->pstore[GMRFLib_PSTORE_TNUM_REF]->done_with_build == GMRFLib_TRUE);
	assert(store->pstore[GMRFLib_PSTORE_TNUM_REF]->Q != NULL);

	int mnum1 = 1, n = store->pstore[GMRFLib_PSTORE_TNUM_REF]->Q->s->n;
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
		omp_set_num_threads(IMIN(store->pstore[GMRFLib_PSTORE_TNUM_REF]->iparm[2], GMRFLib_openmp->max_threads_inner));
		// assert(GMRFLib_openmp->max_threads_inner <= store->pstore[GMRFLib_PSTORE_TNUM_REF]->iparm[2]);
	}

	GMRFLib_csr_tp *Q = store->pstore[GMRFLib_PSTORE_TNUM_REF]->Q;
	pardiso(store->pt, &(store->maxfct), &mnum1, &(store->mtype), &(store->pstore[tnum]->phase),
		&n, Q->a, Q->s->ia1, Q->s->ja1, store->pstore[GMRFLib_PSTORE_TNUM_REF]->perm, &(store->pstore[tnum]->nrhs),
		store->pstore[tnum]->iparm, &(store->msglvl), NULL, NULL, &(store->pstore[tnum]->err_code), store->pstore[tnum]->dparm);

	if (debug) {
		printf("Average number of non-zeros in L per row %.2f\n", store->pstore[tnum]->iparm[17] / (double) n);
	}
	// Revert back to C indexing ?
	if (GMRFLib_imin_value(store->pstore[GMRFLib_PSTORE_TNUM_REF]->perm, n, NULL) == 1) {
#pragma omp simd
		for (int i = 0; i < n; i++) {
			store->pstore[GMRFLib_PSTORE_TNUM_REF]->perm[i]--;
		}
	}
	for (int i = 0; i < n; i++) {
		store->pstore[GMRFLib_PSTORE_TNUM_REF]->iperm[store->pstore[GMRFLib_PSTORE_TNUM_REF]->perm[i]] = i;
	}

	int identity = 1;
	for (int i = 0; i < n && identity; i++) {
		identity = (store->pstore[GMRFLib_PSTORE_TNUM_REF]->iperm[i] == i);
	}
	store->pstore[GMRFLib_PSTORE_TNUM_REF]->perm_identity = identity;

	if (0 && debug) {
		for (int i = 0; i < n; i++) {
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

int GMRFLib_pardiso_solve_core(GMRFLib_pardiso_store_tp *store, GMRFLib_pardiso_flag_tp flag, double *x, double *b, int nrhs)
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
	int nt = 1;
	int nsolve;
	div_t d;

	if (nrhs > 1) {
		if (GMRFLib_openmp->adaptive && omp_get_level() == 0) {
			// this is the exception of the rule, as we want to run this in parallel if we are in adaptive model and
			// level=0.
			nt = store->pstore[GMRFLib_PSTORE_TNUM_REF]->iparm[2];
		} else {
			// this can be used in third level
			if (omp_get_level() > 2) {
				nt = 1;
			} else {
				nt = GMRFLib_openmp->max_threads_inner;
			}
			// assert(GMRFLib_openmp->max_threads_inner <= store->pstore[GMRFLib_PSTORE_TNUM_REF]->iparm[2]);
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
		GMRFLib_pardiso_setparam(flag, store, &tnum);		\
		Memcpy((void *) bb, (void *) (b + offset), n * local_nrhs * sizeof(double));\
		/* if diagonal matrix, do this manually as there is an _error_ or _feature_ in the library */ \
		if (n == Q->s->na) {					\
			assert(store->pstore[GMRFLib_PSTORE_TNUM_REF]->perm_identity); \
			double *xx = x + offset, *qq = Q->a;		\
			if (store->pstore[tnum]->iparm[25] != 0) {	\
				/* in this case, Lx=b, L^Tx=b, are both the same */ \
				for(int k = 0; k < n; k++) {		\
					xx[k] = bb[k] / sqrt(qq[k]); \
				}					\
			} else {					\
				for(int k = 0; k < n; k++) {		\
					xx[k] = bb[k] / qq[k];		\
				}					\
			}						\
		} else {						\
			pardiso(store->pt, &(store->maxfct), &mnum1, &(store->mtype), &(store->pstore[tnum]->phase), \
				&n, Q->a, Q->s->ia1, Q->s->ja1, &idum, &local_nrhs, store->pstore[tnum]->iparm, &(store->msglvl), \
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
						 Q->a, Q->s->ia1, Q->s->ja1, bb + j * n, x + offset + j * n, yy + j * n, &normb, &normr); \
				printf("ni j %d %d The norm of the residual is %e \n ", i, j, normr / normb); \
			}						\
			Free(yy);					\
		}							\
	}

	RUN_CODE_BLOCK(nt, 1, block_nrhs * n);
#undef CODE_BLOCK

	if (err_code) {
		GMRFLib_ERROR(err_code);
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_solve_L(GMRFLib_pardiso_store_tp *store, double *x, double *b, int nrhs)
{
	int res = GMRFLib_pardiso_solve_core(store, GMRFLib_PARDISO_FLAG_SOLVE_L, x, b, nrhs);
	GMRFLib_pardiso_iperm(x, nrhs, store);

	return res;
}

int GMRFLib_pardiso_solve_LT(GMRFLib_pardiso_store_tp *store, double *x, double *b, int nrhs)
{
	GMRFLib_pardiso_perm(b, nrhs, store);
	int res = GMRFLib_pardiso_solve_core(store, GMRFLib_PARDISO_FLAG_SOLVE_LT, x, b, nrhs);
	if (x != b) {
		GMRFLib_pardiso_iperm(b, nrhs, store);
	}

	return res;
}

int GMRFLib_pardiso_solve_LLT(GMRFLib_pardiso_store_tp *store, double *x, double *b, int nrhs)
{
	int res = GMRFLib_pardiso_solve_core(store, GMRFLib_PARDISO_FLAG_SOLVE_LLT, x, b, nrhs);
	return res;
}

double GMRFLib_pardiso_logdet(GMRFLib_pardiso_store_tp *store)
{
	int tnum = omp_get_thread_num();
	return (store->pstore[tnum]->log_det_Q);
}

int GMRFLib_pardiso_bitmap(void)
{
	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_Qinv_INLA(GMRFLib_problem_tp *problem)
{
	if (problem == NULL) {
		return GMRFLib_SUCCESS;
	}

	GMRFLib_pardiso_Qinv(problem->sub_sm_fact.PARDISO_fact);

	GMRFLib_csr_tp *Qi = problem->sub_sm_fact.PARDISO_fact->pstore[GMRFLib_PSTORE_TNUM_REF]->Qinv;
	int n = Qi->s->n;
	map_id **Qinv = Calloc(n, map_id *);

	for (int i = 0, k = 0; i < n; i++) {
		int nnb = Qi->s->ia[i + 1] - Qi->s->ia[i];
		Qinv[i] = Calloc(1, map_id);
		map_id_init_hint(Qinv[i], nnb);
		for (int jj = 0; jj < nnb; jj++) {
			int j = Qi->s->ja[k];
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
#pragma GCC ivdep
	for (int i = 0; i < n; i++) {
		subQinv->mapping[i] = i;
	}
	problem->sub_inverse = subQinv;

	return GMRFLib_SUCCESS;
}

int GMRFLib_pardiso_Qinv(GMRFLib_pardiso_store_tp *store)
{
	int tnum = omp_get_thread_num();
	assert(store->done_with_reorder == GMRFLib_TRUE);
	assert(store->pstore[GMRFLib_PSTORE_TNUM_REF]->done_with_build == GMRFLib_TRUE);
	assert(store->pstore[GMRFLib_PSTORE_TNUM_REF]->done_with_chol == GMRFLib_TRUE);

	if (store->pstore[GMRFLib_PSTORE_TNUM_REF]->Qinv) {
		GMRFLib_csr_free(&(store->pstore[GMRFLib_PSTORE_TNUM_REF]->Qinv));
	}

	GMRFLib_csr_duplicate(&(store->pstore[GMRFLib_PSTORE_TNUM_REF]->Qinv), store->pstore[GMRFLib_PSTORE_TNUM_REF]->Q, 0);
	GMRFLib_pardiso_setparam(GMRFLib_PARDISO_FLAG_QINV, store, NULL);
	int mnum1 = 1;

	if (GMRFLib_openmp->adaptive && omp_get_level() == 0) {
		// this is the exception of the rule, as we want to run this in parallel if we are in adaptive model and
		// level=0.
		omp_set_num_threads(store->pstore[GMRFLib_PSTORE_TNUM_REF]->iparm[2]);
	} else {
		omp_set_num_threads(IMIN(GMRFLib_openmp->max_threads_inner, store->pstore[GMRFLib_PSTORE_TNUM_REF]->iparm[2]));
		// assert(GMRFLib_openmp->max_threads_inner <= store->pstore[GMRFLib_PSTORE_TNUM_REF]->iparm[2]);
	}
	GMRFLib_csr_tp *Qinv = store->pstore[GMRFLib_PSTORE_TNUM_REF]->Qinv;

	pardiso(store->pt, &(store->maxfct), &mnum1, &(store->mtype), &(store->pstore[tnum]->phase),
		&(Qinv->s->n), Qinv->a, Qinv->s->ia1, Qinv->s->ja1, NULL,
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

int GMRFLib_pardiso_pstores_do(GMRFLib_pardiso_store_tp **store)
{
	// if STORE and *STORE is !NULL, then free STORE and its storage.
	// is STORE is NULL, return IDX for a new storage

	int retval = -1;
	int tnum = 0;
	GMRFLib_CACHE_SET_ID(tnum);

	// need critical, as we operate on the global variable S.static_...
#pragma omp critical
	{
		if (store && *store) {
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
				retval = -1;
			} else {
				if (S.s_verbose) {
					PP("free: old=", *store);
				}

				int found = 0;
				if (S.static_pstores != NULL) {
					if (S.s_verbose) {
						for (int i = 0; i < PSTORES_NUM; i++) {
							if (S.busy[i]) {
								printf("in store: i=%1d s=%p\n", i, (void *) S.static_pstores[i]);
							}
						}
					}
					for (int i = 0; i < PSTORES_NUM && !found; i++) {
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
			}
		} else {
			// find a new IDX

			int idx = -1;
			for (int i = 0; i < PSTORES_NUM; i++) {
				if (!S.busy[i]) {
					S.busy[i] = 1;
					idx = i;
					break;
				}
			}
			retval = idx;
		}
	}

	return retval;
}

int GMRFLib_pardiso_free(GMRFLib_pardiso_store_tp **store)
{
	return GMRFLib_pardiso_pstores_do(store);
}

int GMRFLib_duplicate_pardiso_store(GMRFLib_pardiso_store_tp **nnew, GMRFLib_pardiso_store_tp *old, int UNUSED(copy_ptr),
				    int copy_pardiso_ptr, GMRFLib_graph_tp *graph)
{
	int tnum = omp_get_thread_num();
	// if copy_pardiso_ptr, then copy the ptr to read-only objects. 'copy_ptr' is NOT USED
	int debug = S.debug, failsafe_mode = 0;
	if (old == NULL) {
		*nnew = NULL;
		return GMRFLib_SUCCESS;
	}

	if (failsafe_mode) {
		FIXME1("-->duplicate by creating a new one each time");
		GMRFLib_pardiso_init(nnew);
		GMRFLib_pardiso_reorder(*nnew, old->graph);
		return GMRFLib_SUCCESS;
	}
	GMRFLib_ENTER_ROUTINE;

	if (copy_pardiso_ptr) {
#define CP(_what) dup->_what = old->_what
#define CP2(_what) dup->pstore[tnum]->_what = old->pstore[tnum]->_what
#define CP2_ref(_what) dup->pstore[GMRFLib_PSTORE_TNUM_REF]->_what = old->pstore[GMRFLib_PSTORE_TNUM_REF]->_what
#define CPv_ref(_what, type, len)					\
		if (old->pstore[GMRFLib_PSTORE_TNUM_REF]->_what) {	\
			dup->pstore[GMRFLib_PSTORE_TNUM_REF]->_what = Calloc(len, type); \
			Memcpy((void *) (dup->pstore[GMRFLib_PSTORE_TNUM_REF]->_what), \
			       (void *) (old->pstore[GMRFLib_PSTORE_TNUM_REF]->_what), (len) * sizeof(type)); \
		} else {						\
			dup->pstore[GMRFLib_PSTORE_TNUM_REF]->_what = NULL; \
		}							\

		GMRFLib_pardiso_store_tp *dup = Calloc(1, GMRFLib_pardiso_store_tp);
		dup->copy_pardiso_ptr = 1;		       /* YES! */
		for (int i = 0; i < GMRFLib_PARDISO_PLEN; i++) {
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
		for (int i = 0; i < GMRFLib_PARDISO_PLEN; i++) {
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
		CP2(perm_identity);

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
		*nnew = dup;
		return GMRFLib_SUCCESS;
	}

	if (S.static_pstores == NULL) {
#pragma omp critical (Name_046c40f5fd2e202479d5c486dfdf986558c6e681)
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

	int idx = GMRFLib_pardiso_pstores_do(NULL);
	int ok = 1;

	if (S.s_verbose) {
		printf("pstores_do return idx=%d\n", idx);
	}

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
		int redo_reordering = 0;
		*nnew = S.static_pstores[idx];
		if (graph && (graph->sha == NULL || S.static_pstores[idx]->graph->sha == NULL ||
			      strcmp((const char *) S.static_pstores[idx]->graph->sha, (const char *) graph->sha) != 0)) {
			S.static_pstores[idx]->done_with_reorder = GMRFLib_FALSE;
			GMRFLib_pardiso_reorder(S.static_pstores[idx], graph);
			redo_reordering = 1;
		}
		if (S.s_verbose) {
			printf("==> reuse store[%1d]  redo_reordering[%1d]\n", idx, redo_reordering);
		}
	} else {
		if (S.static_pstores[idx]) {
			GMRFLib_pardiso_free(&(S.static_pstores[idx]));
		}
		GMRFLib_pardiso_init(&(S.static_pstores[idx]));
		GMRFLib_pardiso_reorder(S.static_pstores[idx], (graph && graph->n > 0 ? graph : old->graph));
		*nnew = S.static_pstores[idx];
		if (S.s_verbose) {
			printf("==> new store[%1d]\n", idx);
		}
	}

	if (S.s_verbose) {
		printf("duplicate: new=%p old=%p i=%1d\n", *((void **) nnew), ((void *) old), idx);
	}

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}
