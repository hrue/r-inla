#include <assert.h>
#include <math.h>
#include <omp.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#if defined(_WIN32)
#       include <io.h>
#endif

#include "GMRFLib/hashP.h"
#include "GMRFLib/GMRFLib.h"

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
	-2,						       // mtype (-2 = sym, 2 = sym pos def) DO NOT CHANGE!
	0,						       // msg-level (0: no, 1: yes)
	-1,						       // maximum number of rhs
	1,						       // parallel reordering? yes
	NULL,						       // busy
	NULL
};

#define PSTORES_NUM (16384)


int GMRFLib_pardiso_set_debug(int debug)
{
	if (debug) {
		S.s_verbose = (debug ? 1 : 0);
		S.debug = (debug ? 1 : 0);
	}
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
	s->iparm_default[2] = IMAX(GMRFLib_openmp->max_threads_nested[1], GMRFLib_ADAPTIVE_NUM_THREADS());

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
	store->pstore[tnum]->iparm[2] = store->pstore[GMRFLib_PSTORE_TNUM_REF]->iparm[2] = GMRFLib_openmp->max_threads_nested[2];

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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((target_clones(INLA_CLONE_TARGETS "default")))
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
			omp_set_num_threads(GMRFLib_openmp->max_threads_inner);
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
#pragma GCC diagnostic pop

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
		CODE_BLOCK_INIT();				\
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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((target_clones(INLA_CLONE_TARGETS "default")))
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
#pragma GCC diagnostic pop

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((target_clones(INLA_CLONE_TARGETS "default")))
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
		CODE_BLOCK_INIT();					\
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
#pragma GCC diagnostic pop

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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((target_clones(INLA_CLONE_TARGETS "default")))
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
			CODE_BLOCK_INIT();				\
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
#pragma GCC diagnostic pop

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
		omp_set_num_threads(GMRFLib_openmp->max_threads_inner);
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
	GMRFLib_CACHE_SET_IDX(tnum);

	// need critical, as we operate on the global variable S.static_...
#pragma omp critical (Name_62af846ab39772752ebc33a713cb5579f0d0904e)
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
	GMRFLib_ENTER_FUNCTION;

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
		if (S.static_pstores == NULL) {
			if (S.s_verbose) {
				printf("==> init static_pstores\n");
			}
			S.busy = Calloc(PSTORES_NUM, int);
			S.static_pstores = Calloc(PSTORES_NUM, GMRFLib_pardiso_store_tp *);
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
		P(GMRFLib_openmp->max_threads_nested[0]);
		P(GMRFLib_openmp->max_threads_nested[1]);
		P(GMRFLib_openmp->max_threads_nested[2]);
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

	GMRFLib_LEAVE_FUNCTION;
	return GMRFLib_SUCCESS;
}
