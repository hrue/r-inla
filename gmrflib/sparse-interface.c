#include <assert.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>

#include "GMRFLib/GMRFLib.h"

static int csr_store_use = 1;
static map_strvp csr_store;
static int csr_store_must_init = 1;
static int csr_store_debug = 0;

int GMRFLib_csr_init_store(void)
{
	GMRFLib_ENTER_FUNCTION;
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
	GMRFLib_LEAVE_FUNCTION;
	return GMRFLib_SUCCESS;
}

int GMRFLib_csr_free(GMRFLib_csr_tp **csr)
{
	if (*csr) {
		// Free((*csr)->s->iwork);
		// Free((*csr)->s->iwork1);
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
		int n1 = GMRFLib_align_len(n + 1, sizeof(int));
		int len = n1 + na;
		int llen = GMRFLib_align_len(len, sizeof(int));

		(*csr_to)->copy_only = csr_from->copy_only;
		(*csr_to)->s = Calloc(1, GMRFLib_csr_skeleton_tp);
		(*csr_to)->s->sha = NULL;
		(*csr_to)->s = Calloc(1, GMRFLib_csr_skeleton_tp);
		(*csr_to)->s->n = n;
		(*csr_to)->s->na = na;
		(*csr_to)->s->iwork = Malloc(llen, int);
		(*csr_to)->s->ia = (*csr_to)->s->iwork;
		(*csr_to)->s->ja = (*csr_to)->s->iwork + n1;
		Memcpy((void *) ((*csr_to)->s->iwork), (void *) (csr_from->s->iwork), (size_t) llen * sizeof(int));
#if defined(INLA_WITH_PARDISO)
		(*csr_to)->s->iwork1 = Malloc(llen, int);
		(*csr_to)->s->ia1 = (*csr_to)->s->iwork1;
		(*csr_to)->s->ja1 = (*csr_to)->s->iwork1 + n1;
		Memcpy((void *) ((*csr_to)->s->iwork1), (void *) (csr_from->s->iwork1), (size_t) llen * sizeof(int));
#else
		(*csr_to)->s->iwork1 = (*csr_to)->s->ia1 = (*csr_to)->s->ja1 = NULL;
#endif
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
#if defined(INLA_WITH_PARDISO)
	assert(M);
	int mtype = -2;
	int error = 0;
	pardiso_chkmatrix(&mtype, &(M->s->n), M->a, M->s->ia1, M->s->ja1, &error);
	if (error != 0) {
		GMRFLib_ERROR(GMRFLib_ESNH);
	}
#else
	assert(M->s->ia1 == NULL);
	assert(M->s->ja1 == NULL);
#endif
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
		Ms->sha = Malloc(GMRFLib_SHA_DIGEST_LEN + 1, unsigned char);
		Memcpy(Ms->sha, graph->sha, GMRFLib_SHA_DIGEST_LEN + 1);
	}
	n = graph->n;
	na = graph->nnz / 2 + n;			       // only upper triangular. yes, integer division
	n1 = GMRFLib_align_len(n + 1, sizeof(int));
	len = n1 + na;
	llen = GMRFLib_align_len(len, sizeof(int));
	Ms->na = na;
	Ms->n = n;
	Ms->iwork = Malloc(llen, int);
	Ms->ia = Ms->iwork;
	Ms->ja = Ms->iwork + n1;
#if defined(INLA_WITH_PARDISO)
	Ms->iwork1 = Malloc(llen, int);
	Ms->ia1 = Ms->iwork1;
	Ms->ja1 = Ms->iwork1 + n1;
#else
	Ms->iwork1 = Ms->ia1 = Ms->ja1 = NULL;
#endif
	// new code. by doing it in two steps we can do the second one in parallel, and this is the one that take time.
	int *k_arr = Malloc(n, int);
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
		CODE_BLOCK_INIT();					\
		if (graph->lnnbs[i]) {					\
			int k = k_arr[i];				\
			Memcpy(&(Ms->ja[k]), graph->lnbs[i], graph->lnnbs[i] * sizeof(int)); \
		}							\
	}

	RUN_CODE_BLOCK(1, 0, 0);			       /* TO QUICK TO DO IN PARALLEL */
#undef CODE_BLOCK
	Free(k_arr);

#if defined(INLA_WITH_PARDISO)
#       pragma omp simd
	for (int i = 0; i < n + 1; i++) {
		Ms->ia1[i] = Ms->ia[i] + 1;
	}
#       pragma omp simd
	for (int i = 0; i < na; i++) {
		Ms->ja1[i] = Ms->ja[i] + 1;
	}
#endif

	if (csr_store_use && graph->sha) {
		if (csr_store_debug) {
			printf("\t[%1d] csr_store: store crs 0x%p\n", omp_get_thread_num(), (void *) Ms);
		}
#pragma omp critical (Name_488cde57983063f09dc9ddf7c473d77c79ea2929)
		map_strvp_set(&csr_store, (char *) graph->sha, (void *) Ms);
	}

	return (Ms);
}

int GMRFLib_Q2csr(int thread_id, GMRFLib_csr_tp **csr, GMRFLib_graph_tp *graph, GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_arg)
{
	GMRFLib_ENTER_FUNCTION;

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
		M->a = Malloc(M->s->na, double);
		// a bit more manual work
		double val = Qfunc(thread_id, 0, -1, &(M->a[0]), Qfunc_arg);
		if (ISNAN(val)) {
			static char *tag = NULL;
			if (!tag) {
#pragma omp critical (Name_7600f798b7727e8eb5fbed77a2db305e4db69365)
				if (!tag) {
					GMRFLib_sprintf(&tag, "%s:%1d", __FILE__, __LINE__);
				}
			}
			int level = omp_get_level();
			int tnum = omp_get_thread_num();
			int num_threads = (level == 0 ? GMRFLib_ADAPTIVE_NUM_THREADS() : GMRFLib_openmp->max_threads_inner);
			int nt_loc = GMRFLib_adapt_nt_get(tag, tnum, level, num_threads);
			double tref = -GMRFLib_timer();

#define CODE_BLOCK							\
			for (int i = 0; i < M->s->n; i++) {		\
				CODE_BLOCK_INIT();			\
				for (int k = M->s->ia[i]; k < M->s->ia[i + 1]; k++) { \
					int j = M->s->ja[k];		\
					M->a[k] = Qfunc(thread_id, i, j, NULL, Qfunc_arg); \
				}					\
			}

			RUN_CODE_BLOCK(nt_loc, 0, 0);
#undef CODE_BLOCK
			tref += GMRFLib_timer();
			GMRFLib_adapt_nt_update(tag, tnum, level, tref);

		} else {

			static char *tag = NULL;
			if (!tag) {
#pragma omp critical (Name_5d44f84bdfc2d2b324a71dbddd46ed4c72f4fba7)
				if (!tag) {
					GMRFLib_sprintf(&tag, "%s:%1d", __FILE__, __LINE__);
				}
			}
			int level = omp_get_level();
			int tnum = omp_get_thread_num();
			int num_threads = (level == 0 ? GMRFLib_ADAPTIVE_NUM_THREADS() : GMRFLib_openmp->max_threads_inner);
			int nt_loc = GMRFLib_adapt_nt_get(tag, tnum, level, num_threads);
			double tref = -GMRFLib_timer();

#define CODE_BLOCK							\
			for (int i = 0; i < M->s->n; i++) {		\
				CODE_BLOCK_INIT();			\
				int k = M->s->ia[i];			\
				Qfunc(thread_id, i, -1, &(M->a[k]), Qfunc_arg);	\
			}

			RUN_CODE_BLOCK(nt_loc, 0, 0);
#undef CODE_BLOCK
			tref += GMRFLib_timer();
			GMRFLib_adapt_nt_update(tag, tnum, level, tref);
		}
	}

	int nan_error = 0;
	for (int i = 0; i < M->s->na; i++) {
		GMRFLib_STOP_IF_NAN_OR_INF(M->a[i], i, -1);
		if (nan_error) {
			GMRFLib_LEAVE_FUNCTION;
			return !GMRFLib_SUCCESS;
		}
	}

#undef M
	GMRFLib_LEAVE_FUNCTION;
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
#if defined(INLA_WITH_PARDISO)
	GMRFLib_io_write(io, (const void *) (csr->s->ia1), sizeof(int) * (csr->s->n + 1));
	GMRFLib_io_write(io, (const void *) (csr->s->ja1), sizeof(int) * csr->s->na);
#else
	csr->s->ia1 = csr->s->ja1 = NULL;
#endif
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
	M->s->iwork = Malloc(len, int);
	M->s->ia = M->s->iwork;
	GMRFLib_io_read(io, (void *) (M->s->ia), sizeof(int) * (M->s->n + 1));
	M->s->ja = M->s->iwork + M->s->n + 1;
	GMRFLib_io_read(io, (void *) (M->s->ja), sizeof(int) * M->s->na);
#if defined(INLA_WITH_PARDISO)
	M->s->iwork1 = Malloc(len, int);
	M->s->ia1 = M->s->iwork1;
	GMRFLib_io_read(io, (void *) (M->s->ia1), sizeof(int) * (M->s->n + 1));
	M->s->ja1 = M->s->iwork1 + M->s->n + 1;
	GMRFLib_io_read(io, (void *) (M->s->ja1), sizeof(int) * M->s->na);
#else
	M->s->iwork1 = M->s->ia1 = M->s->ja1 = NULL;
#endif

	M->a = Malloc(M->s->na, double);
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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((target_clones(INLA_CLONE_TARGETS "default")))
int GMRFLib_csr2Q(GMRFLib_tabulate_Qfunc_tp **Qtab, GMRFLib_graph_tp **graph, GMRFLib_csr_tp *csr)
{
	int i, j, k, jj, nnb;

	int *iarr = Malloc(csr->s->na, int);
	int *jarr = Malloc(csr->s->na, int);
	double *arr = Malloc(csr->s->na, double);

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
#pragma GCC diagnostic pop

int GMRFLib_compute_reordering(GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph, GMRFLib_global_node_tp *gn)
{
	GMRFLib_ENTER_FUNCTION;
	GMRFLib_global_node_tp lgn, *gn_ptr = NULL;

	if (gn) {
		/*
		 * then this defines the global node definition 
		 */
		gn_ptr = (GMRFLib_global_node_tp *) gn;
	} else {
		lgn = GMRFLib_global_node;
		gn_ptr = &lgn;
	}

	GMRFLib_reorder_tp r = GMRFLib_reorder;
	if (sm_fact->smtp == GMRFLib_SMTP_PARDISO || sm_fact->smtp == GMRFLib_SMTP_STILES) {
		r = GMRFLib_REORDER_DEFAULT;
	} else if ((sm_fact->smtp == GMRFLib_SMTP_TAUCS || sm_fact->smtp == GMRFLib_SMTP_BAND) &&
		   (r == GMRFLib_REORDER_STILES || r == GMRFLib_REORDER_PARDISO)) {
		r = GMRFLib_REORDER_DEFAULT;
	}

	switch (r) {
	case GMRFLib_REORDER_DEFAULT:
	{
		/*
		 * this choice depends on the sparse-solver 
		 */
		switch (sm_fact->smtp) {
		case GMRFLib_SMTP_BAND:
		{
			GMRFLib_EWRAP1(GMRFLib_compute_reordering_BAND(&(sm_fact->remap), graph));
		}
			break;

		case GMRFLib_SMTP_TAUCS:
		{
			GMRFLib_EWRAP1(GMRFLib_compute_reordering_TAUCS(&(sm_fact->remap), graph, r, gn_ptr));
		}
			break;

		case GMRFLib_SMTP_PARDISO:
		{
			if (sm_fact->PARDISO_fact == NULL) {
				GMRFLib_pardiso_init(&(sm_fact->PARDISO_fact));
			}
			GMRFLib_pardiso_reorder(sm_fact->PARDISO_fact, graph);
			sm_fact->remap = Malloc(graph->n, int);
			Memcpy((void *) sm_fact->remap, (void *) sm_fact->PARDISO_fact->pstore[GMRFLib_PSTORE_TNUM_REF]->perm,
			       graph->n * sizeof(int));
		}
			break;

		case GMRFLib_SMTP_STILES:
		{
			int k = -1;
			GMRFLib_stiles_store_tp *p = (GMRFLib_stiles_store_tp *) GMRFLib_stiles_get_store_ptr();
			for (int i = 0; i < p->graphs->n; i++) {
				GMRFLib_graph_tp *g = (GMRFLib_graph_tp *) (p->graphs->ptr[i]);
				if (strcmp((const char *) graph->sha, (const char *) g->sha) == 0) {
					k = i;
					break;
				}
			}
			assert(k >= 0);

			GMRFLib_stiles_idx_tp stiles_idx = { k, 0, 0 };
			sm_fact->remap = Malloc(graph->n, int);
			Memcpy((void *) sm_fact->remap, (void *) GMRFLib_stiles_get_perm(&stiles_idx), graph->n * sizeof(int));
		}
			break;

		default:
			GMRFLib_ASSERT(1 == 0, GMRFLib_ESNH);
			break;
		}
	}
		break;

	case GMRFLib_REORDER_PARDISO:
	{
		switch (sm_fact->smtp) {
		case GMRFLib_SMTP_PARDISO:		       /* same code as above */
		{
			if (sm_fact->PARDISO_fact == NULL) {
				GMRFLib_pardiso_init(&(sm_fact->PARDISO_fact));
			}
			GMRFLib_pardiso_reorder(sm_fact->PARDISO_fact, graph);
			sm_fact->remap = Malloc(graph->n, int);
			Memcpy((void *) sm_fact->remap, (void *) sm_fact->PARDISO_fact->pstore[GMRFLib_PSTORE_TNUM_REF]->perm,
			       graph->n * sizeof(int));
		}
			break;

		default:
			GMRFLib_ASSERT(1 == 0, GMRFLib_ESNH);
			break;
		}
	}
		break;

	case GMRFLib_REORDER_BAND:
		GMRFLib_EWRAP1(GMRFLib_compute_reordering_BAND(&(sm_fact->remap), graph));
		break;

		/*
		 * all the remaining ones are treated by the _TAUCS routine 
		 */
	case GMRFLib_REORDER_IDENTITY:
	case GMRFLib_REORDER_REVERSE_IDENTITY:
	case GMRFLib_REORDER_METIS:
	case GMRFLib_REORDER_GENMMD:
	case GMRFLib_REORDER_AMD:
	case GMRFLib_REORDER_AMDC:
	case GMRFLib_REORDER_AMDBAR:
	case GMRFLib_REORDER_AMDBARC:
	case GMRFLib_REORDER_MD:
	case GMRFLib_REORDER_MMD:
	{
		GMRFLib_EWRAP1(GMRFLib_compute_reordering_TAUCS(&(sm_fact->remap), graph, r, gn_ptr));
	}
		break;

	case GMRFLib_REORDER_STILES:
	{
		GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
	}
		break;

	default:
	{
		GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
	}
		break;
	}

	if (sm_fact->remap) {				       /* need this still for the wa-routines. FIXME */
		GMRFLib_EWRAP1(GMRFLib_graph_comp_bw(&(sm_fact->bandwidth), graph, sm_fact->remap));
	} else {
		sm_fact->bandwidth = -1;
	}

	GMRFLib_LEAVE_FUNCTION;

	return GMRFLib_SUCCESS;
}

int GMRFLib_free_reordering(GMRFLib_sm_fact_tp *sm_fact)
{
	if (sm_fact) {
		Free(sm_fact->remap);
		sm_fact->bandwidth = 0;
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_build_sparse_matrix(int thread_id, GMRFLib_sm_fact_tp *sm_fact, GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_arg,
				GMRFLib_graph_tp *graph, GMRFLib_problem_tp *problem)
{
	GMRFLib_ENTER_FUNCTION;
	int ret = 0;

	switch (sm_fact->smtp) {
	case GMRFLib_SMTP_BAND:
	{
		ret = GMRFLib_build_sparse_matrix_BAND(thread_id, &(sm_fact->bchol), Qfunc, Qfunc_arg, graph, sm_fact->remap, sm_fact->bandwidth);
		if (ret != GMRFLib_SUCCESS) {
			GMRFLib_LEAVE_FUNCTION;
			return ret;
		}
	}
		break;

	case GMRFLib_SMTP_TAUCS:
	{
		ret = GMRFLib_build_sparse_matrix_TAUCS(thread_id, &(sm_fact->TAUCS_L), Qfunc, Qfunc_arg, graph, sm_fact->remap);
		if (ret != GMRFLib_SUCCESS) {
			GMRFLib_LEAVE_FUNCTION;
			return ret;
		}
	}
		break;

	case GMRFLib_SMTP_PARDISO:
	{
		if (sm_fact->PARDISO_fact == NULL) {
			GMRFLib_pardiso_init(&(sm_fact->PARDISO_fact));
			GMRFLib_pardiso_reorder(sm_fact->PARDISO_fact, graph);
		}
		ret = GMRFLib_pardiso_build(thread_id, sm_fact->PARDISO_fact, graph, Qfunc, Qfunc_arg);
		if (ret != GMRFLib_SUCCESS) {
			GMRFLib_LEAVE_FUNCTION;
			return ret;
		}
	}
		break;

	case GMRFLib_SMTP_STILES:
	{
		ret = GMRFLib_stiles_build(problem->stiles_idx, thread_id, Qfunc, Qfunc_arg);
		if (ret != GMRFLib_SUCCESS) {
			GMRFLib_LEAVE_FUNCTION;
			return ret;
		}
	}
		break;

	default:
		GMRFLib_ASSERT(1 == 0, GMRFLib_ESNH);
		break;
	}

	GMRFLib_LEAVE_FUNCTION;
	return GMRFLib_SUCCESS;
}

int GMRFLib_factorise_sparse_matrix(GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph, GMRFLib_problem_tp *problem)
{
	int ret;
	GMRFLib_ENTER_FUNCTION;

	switch (sm_fact->smtp) {
	case GMRFLib_SMTP_BAND:
	{
		ret = GMRFLib_factorise_sparse_matrix_BAND(sm_fact->bchol, &(sm_fact->finfo), graph, sm_fact->bandwidth);
		if (ret != GMRFLib_SUCCESS) {
			GMRFLib_LEAVE_FUNCTION;
			return ret;
		}
	}
		break;

	case GMRFLib_SMTP_TAUCS:
	{
		ret =
		    GMRFLib_factorise_sparse_matrix_TAUCS(&(sm_fact->TAUCS_L), &(sm_fact->TAUCS_symb_fact), &(sm_fact->TAUCS_cache),
							  &(sm_fact->finfo));
		if (ret != GMRFLib_SUCCESS) {
			GMRFLib_LEAVE_FUNCTION;
			return ret;
		}
	}
		break;

	case GMRFLib_SMTP_PARDISO:
	{
		ret = GMRFLib_pardiso_chol(sm_fact->PARDISO_fact);
		if (ret != GMRFLib_SUCCESS) {
			GMRFLib_LEAVE_FUNCTION;
			return ret;
		}
	}
		break;

	case GMRFLib_SMTP_STILES:
	{
		ret = GMRFLib_stiles_chol(problem->stiles_idx);
		if (ret != GMRFLib_SUCCESS) {
			GMRFLib_LEAVE_FUNCTION;
			return ret;
		}
	}
		break;

	default:
		GMRFLib_ASSERT(1 == 0, GMRFLib_ESNH);
		break;
	}

	GMRFLib_LEAVE_FUNCTION;
	return GMRFLib_SUCCESS;
}

int GMRFLib_free_fact_sparse_matrix(GMRFLib_sm_fact_tp *sm_fact)
{
	if (sm_fact) {
		switch (sm_fact->smtp) {
		case GMRFLib_SMTP_BAND:
		{
			GMRFLib_EWRAP1(GMRFLib_free_fact_sparse_matrix_BAND(sm_fact->bchol));
			sm_fact->bchol = NULL;
		}
			break;

		case GMRFLib_SMTP_TAUCS:
		{
			GMRFLib_free_fact_sparse_matrix_TAUCS(sm_fact->TAUCS_L, sm_fact->TAUCS_LL, sm_fact->TAUCS_symb_fact);
			GMRFLib_taucs_cache_free(sm_fact->TAUCS_cache);
			sm_fact->TAUCS_L = NULL;
			sm_fact->TAUCS_symb_fact = NULL;
			sm_fact->TAUCS_symb_fact = NULL;
			sm_fact->TAUCS_cache = NULL;
		}
			break;

		case GMRFLib_SMTP_PARDISO:
		{
			if (sm_fact->PARDISO_fact) {
				GMRFLib_pardiso_free(&(sm_fact->PARDISO_fact));
				sm_fact->PARDISO_fact = NULL;
			}
		}
			break;

		case GMRFLib_SMTP_STILES:
			break;

		default:
			GMRFLib_ASSERT(1 == 0, GMRFLib_ESNH);
			break;
		}
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_solve_l_sparse_matrix(double *rhs, int nrhs, GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph, GMRFLib_problem_tp *problem)
{
	/*
	 * rhs in real world. solve L x=rhs, rhs is overwritten by the solution 
	 */
	GMRFLib_ENTER_FUNCTION;

	switch (sm_fact->smtp) {
	case GMRFLib_SMTP_BAND:
	{
		omp_set_num_threads(GMRFLib_openmp->max_threads_inner);
#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_inner) schedule(static)
		for (int i = 0; i < nrhs; i++) {
			GMRFLib_solve_l_sparse_matrix_BAND(&rhs[i * graph->n], sm_fact->bchol, graph, sm_fact->remap, sm_fact->bandwidth);
		}
	}
		break;

	case GMRFLib_SMTP_TAUCS:
	{
		omp_set_num_threads(GMRFLib_openmp->max_threads_inner);
#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_inner) schedule(static)
		for (int i = 0; i < nrhs; i++) {
			GMRFLib_solve_l_sparse_matrix_TAUCS(&rhs[i * graph->n], sm_fact->TAUCS_L, graph, sm_fact->remap);
		}
	}
		break;

	case GMRFLib_SMTP_PARDISO:
	{
		GMRFLib_pardiso_solve_L(sm_fact->PARDISO_fact, rhs, rhs, nrhs);
	}
		break;

	case GMRFLib_SMTP_STILES:
	{
		GMRFLib_stiles_idx_tp stiles_idx = { problem->stiles_idx->in_group, 0, 0 };
		int err = GMRFLib_stiles_set_idx(&stiles_idx, nrhs);
		if (err == GMRFLib_SUCCESS) {
			GMRFLib_stiles_solve_L(&stiles_idx, rhs);
		} else {
			GMRFLib_stiles_set_idx(&stiles_idx, 1);
			for (int k = 0; k < nrhs; k++) {
				GMRFLib_stiles_solve_L(&stiles_idx, rhs + k * graph->n);
			}
		}
	}
		break;

	default:
		GMRFLib_ERROR(GMRFLib_ESNH);
		break;
	}

	GMRFLib_LEAVE_FUNCTION;

	return GMRFLib_SUCCESS;
}

int GMRFLib_solve_lt_sparse_matrix(double *rhs, int nrhs, GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph, GMRFLib_problem_tp *problem)
{
	/*
	 * rhs in real world. solve L^Tx=rhs, rhs is overwritten by the solution 
	 */
	GMRFLib_ENTER_FUNCTION;

	switch (sm_fact->smtp) {
	case GMRFLib_SMTP_BAND:
	{
		omp_set_num_threads(GMRFLib_openmp->max_threads_inner);
#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_inner) schedule(static)
		for (int i = 0; i < nrhs; i++) {
			GMRFLib_solve_lt_sparse_matrix_BAND(&rhs[i * graph->n], sm_fact->bchol, graph, sm_fact->remap, sm_fact->bandwidth);
		}
	}
		break;

	case GMRFLib_SMTP_TAUCS:
	{
		omp_set_num_threads(GMRFLib_openmp->max_threads_inner);
#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_inner) schedule(static)
		for (int i = 0; i < nrhs; i++) {
			GMRFLib_solve_lt_sparse_matrix_TAUCS(&rhs[i * graph->n], sm_fact->TAUCS_L, graph, sm_fact->remap);
		}
	}
		break;

	case GMRFLib_SMTP_PARDISO:
	{
		GMRFLib_pardiso_solve_LT(sm_fact->PARDISO_fact, rhs, rhs, nrhs);
	}
		break;

	case GMRFLib_SMTP_STILES:
	{
		GMRFLib_stiles_idx_tp stiles_idx = { problem->stiles_idx->in_group, 0, 0 };
		int err = GMRFLib_stiles_set_idx(&stiles_idx, nrhs);
		if (err == GMRFLib_SUCCESS) {
			GMRFLib_stiles_solve_LT(&stiles_idx, rhs);
		} else {
			GMRFLib_stiles_set_idx(&stiles_idx, 1);
			for (int k = 0; k < nrhs; k++) {
				GMRFLib_stiles_solve_LT(&stiles_idx, rhs + k * graph->n);
			}
		}
	}
		break;

	default:
		GMRFLib_ERROR(GMRFLib_ESNH);
		break;
	}

	GMRFLib_LEAVE_FUNCTION;

	return GMRFLib_SUCCESS;
}

int GMRFLib_solve_llt_sparse_matrix(double *rhs, int nrhs, GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph,
				    GMRFLib_problem_tp *problem, GMRFLib_stiles_idx_tp *stiles_idx)
{
	// if (nrhs > 1) P(nrhs);

	/*
	 * rhs in real world. solve Q x=rhs, where Q=L L^T 
	 */
	GMRFLib_ENTER_FUNCTION;

	static double **wwork = NULL;
	static int *wwork_len = NULL;
	if (!wwork) {
#pragma omp critical (Name_c02cfe7c85f984ba167d3d158f5219787998c27f)
		if (!wwork) {
			wwork_len = Calloc(GMRFLib_CACHE_LEN(), int);
			double **tmp = Calloc(GMRFLib_CACHE_LEN(), double *);
			wwork = tmp;
		}
	}

	int cache_idx = 0;
	GMRFLib_CACHE_SET_IDX(cache_idx);

	int nw = graph->n * nrhs;
	if (nw > wwork_len[cache_idx]) {
		int numa_node = GMRFLib_numa_get_node();
		if (wwork[cache_idx] && wwork_len[cache_idx] > 0) {
			GMRFLib_numa_free(wwork[cache_idx], wwork_len[cache_idx]);
		}
		wwork_len[cache_idx] = nw;
		wwork[cache_idx] = (double *) GMRFLib_numa_alloc_onnode(wwork_len[cache_idx] * sizeof(double), numa_node);
	}
	double *work = wwork[cache_idx];

	if (sm_fact->smtp == GMRFLib_SMTP_BAND) {
		omp_set_num_threads(GMRFLib_openmp->max_threads_inner);
#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_inner) schedule(static)
		for (int i = 0; i < nrhs; i++) {
			int offset = i * graph->n;
			GMRFLib_solve_llt_sparse_matrix_BAND(rhs + offset, sm_fact->bchol, graph, sm_fact->remap,
							     sm_fact->bandwidth, work + offset);
		}
	} else if (sm_fact->smtp == GMRFLib_SMTP_TAUCS) {
		int min_block_size = GMRFLib_taucs_get_min_block_size();
		//int block_size = GMRFLib_taucs_get_block_size();
		int ntt = (omp_get_level() == 0 ? GMRFLib_openmp->max_threads_outer :
			   (omp_get_level() == 1 ? GMRFLib_openmp->max_threads_inner : 1));
		if (nrhs == 1) {
			// simple entry, part 1
			GMRFLib_solve_llt_sparse_matrix_TAUCS(rhs, sm_fact->TAUCS_L, sm_fact->TAUCS_LL, graph, sm_fact->remap, work);
		} else if (ntt ==  1) {
			// simple entry, part 2
			GMRFLib_solve_llt_sparse_matrix2_TAUCS(rhs, sm_fact->TAUCS_L, graph, sm_fact->remap, nrhs, work);
		} else {
			int len = GMRFLib_align_len(ntt, sizeof(int));
			int *iwork = Calloc(2 * len, int);
			int *csize = iwork;
			int *offset = iwork + len;

			GMRFLib_ifill(ntt, 0, csize);
			int done = 0;
			while (!done) {
				for (int i = 0; i < ntt && !done; i++) {
					int off = IMIN(nrhs - GMRFLib_isum(ntt, csize), min_block_size);
					csize[i] += off;
					if (GMRFLib_isum(ntt, csize) == nrhs) {
						done = 1;
					}
				}
			}

			// cleanup if a thread as few nrhs, then move those to the earlier one
			for (int i = ntt-1; i > 0; i--) {
				if (csize[i] < min_block_size / 2) {
					csize[i-1] += csize[i];
					csize[i] = 0;
				}
			}

			// int ntt_orig = ntt;
			int new_ntt = 0;
			for (int i = 0; i < ntt; i++) {
				new_ntt +=  (csize[i] > 0);
			}
			ntt = new_ntt;

			// printf("CALL with nrhs %1d and ntt_orig %1d ntt %1d\n", nrhs, ntt_orig, ntt);
			// for(int i = 0; i < ntt_orig; i++) printf("\tthread %1d nrhs %1d\n", i, csize[i]);
			
			assert(GMRFLib_isum(ntt, csize) == nrhs);
			offset[0] = 0;
			for (int i = 1; i < ntt; i++) {
				offset[i] = offset[i - 1] + csize[i - 1] * graph->n;
			}

#pragma omp parallel for num_threads(ntt) schedule(static)
			for (int i = 0; i < ntt; i++) {
				GMRFLib_solve_llt_sparse_matrix2_TAUCS(rhs + offset[i], sm_fact->TAUCS_L, graph,
								       sm_fact->remap, csize[i], work + offset[i]);
			}
			Free(iwork);
			
		}
	} else if (sm_fact->smtp == GMRFLib_SMTP_PARDISO) {
		GMRFLib_EWRAP1(GMRFLib_pardiso_solve_LLT(sm_fact->PARDISO_fact, rhs, rhs, nrhs));
	} else if (sm_fact->smtp == GMRFLib_SMTP_STILES) {
		GMRFLib_stiles_idx_tp s_idx = { 0, 0, 0 };
		if (stiles_idx) {
			s_idx.in_group = stiles_idx->in_group;
		} else {
			s_idx.in_group = problem->stiles_idx->in_group;
		}

		GMRFLib_stiles_set_idx(&s_idx, nrhs);
		GMRFLib_stiles_solve_LLT(&s_idx, rhs);
	} else {
		GMRFLib_ERROR(GMRFLib_ESNH);
	}
	GMRFLib_LEAVE_FUNCTION;

	return GMRFLib_SUCCESS;
}

int GMRFLib_solve_llt_sparse_matrix_special(double *rhs, GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph, int idx, GMRFLib_problem_tp *problem)
{
	/*
	 * rhs in real world. solve Q x=rhs, where Q=L L^T. BUT, here we know that rhs is 0 execpt for a 1 at index idx.
	 */
	GMRFLib_ENTER_FUNCTION;

	switch (sm_fact->smtp) {
	case GMRFLib_SMTP_BAND:
	{
		GMRFLib_solve_llt_sparse_matrix_special_BAND(rhs, sm_fact->bchol, graph, sm_fact->remap, sm_fact->bandwidth, idx);
	}
		break;

	case GMRFLib_SMTP_TAUCS:
	{
		GMRFLib_solve_llt_sparse_matrix_special_TAUCS(rhs, sm_fact->TAUCS_L, graph, sm_fact->remap, idx);
	}
		break;

	case GMRFLib_SMTP_PARDISO:
	{
		GMRFLib_pardiso_solve_LLT(sm_fact->PARDISO_fact, rhs, rhs, 1);
	}
		break;

	case GMRFLib_SMTP_STILES:
	{
		GMRFLib_stiles_idx_tp stiles_idx = { problem->stiles_idx->in_group, 0, 0 };
		GMRFLib_stiles_set_idx(&stiles_idx, 1);
		GMRFLib_stiles_solve_LLT(&stiles_idx, rhs);
	}
		break;

	default:
		GMRFLib_ERROR(GMRFLib_ESNH);
		break;
	}

	GMRFLib_LEAVE_FUNCTION;

	return GMRFLib_SUCCESS;
}

int GMRFLib_solve_lt_sparse_matrix_special(double *rhs, GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph, int findx, int toindx,
					   int remapped, GMRFLib_problem_tp *problem)
{
	/*
	 * rhs in real world, bchol in mapped world. solve L^Tx=b backward only from rhs[findx] up to rhs[toindx]. note that
	 * findx and toindx is in mapped world. if remapped, do not remap/remap-back the rhs before solving.
	 */
	GMRFLib_ENTER_FUNCTION;
	switch (sm_fact->smtp) {
	case GMRFLib_SMTP_BAND:
	{
		GMRFLib_EWRAP0(GMRFLib_solve_lt_sparse_matrix_special_BAND
			       (rhs, sm_fact->bchol, graph, sm_fact->remap, sm_fact->bandwidth, findx, toindx, remapped));
	}
		break;

	case GMRFLib_SMTP_TAUCS:
	{
		GMRFLib_EWRAP0(GMRFLib_solve_lt_sparse_matrix_special_TAUCS(rhs, sm_fact->TAUCS_L, graph, sm_fact->remap, findx, toindx, remapped));
	}
		break;

	case GMRFLib_SMTP_PARDISO:
	{
		GMRFLib_pardiso_solve_LT(sm_fact->PARDISO_fact, rhs, rhs, 1);
	}
		break;

	case GMRFLib_SMTP_STILES:
	{
		GMRFLib_stiles_idx_tp stiles_idx = { problem->stiles_idx->in_group, 0, 0 };
		GMRFLib_stiles_set_idx(&stiles_idx, 1);
		GMRFLib_stiles_solve_LT(&stiles_idx, rhs);
	}
		break;

	default:
	{
		GMRFLib_ERROR(GMRFLib_ESNH);
	}
		break;
	}
	GMRFLib_LEAVE_FUNCTION;

	return GMRFLib_SUCCESS;
}

int GMRFLib_solve_l_sparse_matrix_special(double *rhs, GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph, int findx, int toindx,
					  int remapped, GMRFLib_problem_tp *problem)
{
	/*
	 * rhs in real world, bchol in mapped world. solve Lx=b backward only from rhs[findx] up to rhs[toindx]. note that
	 * findx and toindx is in mapped world. if remapped, do not remap/remap-back the rhs before solving.
	 */

	GMRFLib_ENTER_FUNCTION;
	switch (sm_fact->smtp) {
	case GMRFLib_SMTP_BAND:
	{
		GMRFLib_EWRAP0(GMRFLib_solve_l_sparse_matrix_special_BAND
			       (rhs, sm_fact->bchol, graph, sm_fact->remap, sm_fact->bandwidth, findx, toindx, remapped));
	}
		break;

	case GMRFLib_SMTP_TAUCS:
	{
		GMRFLib_EWRAP0(GMRFLib_solve_l_sparse_matrix_special_TAUCS(rhs, sm_fact->TAUCS_L, graph, sm_fact->remap, findx, toindx, remapped));
	}
		break;

	case GMRFLib_SMTP_PARDISO:
	{
		if (remapped) {
			GMRFLib_pardiso_perm(rhs, 1, sm_fact->PARDISO_fact);
		}
		GMRFLib_pardiso_solve_L(sm_fact->PARDISO_fact, rhs, rhs, 1);
	}
		break;

	case GMRFLib_SMTP_STILES:
	{
		if (remapped) {
			FIXME1("UNSURE ABOUT THIS THING");
			int *perm = GMRFLib_stiles_get_perm(problem->stiles_idx);
			double *y = Malloc(graph->n, double);
			Memcpy(y, rhs, graph->n * sizeof(double));
			GMRFLib_pack(graph->n, y, perm, rhs);
			Free(y);
		}
		GMRFLib_stiles_idx_tp stiles_idx = { problem->stiles_idx->in_group, 0, 0 };
		GMRFLib_stiles_set_idx(&stiles_idx, 1);
		GMRFLib_stiles_solve_L(&stiles_idx, rhs);
	}
		break;

	default:
		GMRFLib_ERROR(GMRFLib_ESNH);
		break;
	}

	GMRFLib_LEAVE_FUNCTION;
	return GMRFLib_SUCCESS;
}

int GMRFLib_log_determinant(double *logdet, GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph, GMRFLib_problem_tp *problem)
{
	switch (sm_fact->smtp) {
	case GMRFLib_SMTP_BAND:
	{
		GMRFLib_EWRAP0(GMRFLib_log_determinant_BAND(logdet, sm_fact->bchol, graph, sm_fact->bandwidth));
	}
		break;

	case GMRFLib_SMTP_TAUCS:
	{
		GMRFLib_EWRAP0(GMRFLib_log_determinant_TAUCS(logdet, sm_fact->TAUCS_L));
	}
		break;

	case GMRFLib_SMTP_PARDISO:
	{
		*logdet = GMRFLib_pardiso_logdet(sm_fact->PARDISO_fact);
	}
		break;

	case GMRFLib_SMTP_STILES:
	{
		*logdet = GMRFLib_stiles_logdet(problem->stiles_idx);
	}
		break;

	default:
		GMRFLib_ERROR(GMRFLib_ESNH);
		break;
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_comp_cond_meansd(double *cmean, double *csd, int indx, double *x, int remapped, GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph)
{
	GMRFLib_ENTER_FUNCTION;
	switch (sm_fact->smtp) {
	case GMRFLib_SMTP_BAND:
	{
		GMRFLib_EWRAP1(GMRFLib_comp_cond_meansd_BAND
			       (cmean, csd, indx, x, remapped, sm_fact->bchol, graph, sm_fact->remap, sm_fact->bandwidth));
	}
		break;

	case GMRFLib_SMTP_TAUCS:
	{
		GMRFLib_EWRAP1(GMRFLib_comp_cond_meansd_TAUCS(cmean, csd, indx, x, remapped, sm_fact->TAUCS_L, graph, sm_fact->remap));
	}
		break;

	case GMRFLib_SMTP_PARDISO:
	case GMRFLib_SMTP_STILES:
	{
		assert(0 == 1);
	}
		break;

	default:
		GMRFLib_ERROR(GMRFLib_ESNH);
		break;
	}
	GMRFLib_LEAVE_FUNCTION;

	return GMRFLib_SUCCESS;
}

int GMRFLib_bitmap_factorisation(const char *filename_body, GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph)
{
	switch (sm_fact->smtp) {
	case GMRFLib_SMTP_BAND:
	{
		GMRFLib_EWRAP1(GMRFLib_bitmap_factorisation_BAND(filename_body, sm_fact->bchol, graph, sm_fact->remap, sm_fact->bandwidth));
	}
		break;

	case GMRFLib_SMTP_TAUCS:
	{
		GMRFLib_EWRAP1(GMRFLib_bitmap_factorisation_TAUCS(filename_body, sm_fact->TAUCS_L));
	}
		break;

	case GMRFLib_SMTP_PARDISO:
	{
		GMRFLib_EWRAP1(GMRFLib_pardiso_bitmap());
	}
		break;

	case GMRFLib_SMTP_STILES:
		break;

	default:
		GMRFLib_ERROR(GMRFLib_ESNH);
		break;
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_compute_Qinv(void *problem)
{
	GMRFLib_ENTER_FUNCTION;
	GMRFLib_problem_tp *p = (GMRFLib_problem_tp *) problem;

	switch (p->sub_sm_fact.smtp) {
	case GMRFLib_SMTP_BAND:
	{
		GMRFLib_EWRAP0(GMRFLib_compute_Qinv_BAND(p));
	}
		break;

	case GMRFLib_SMTP_TAUCS:
	{
		GMRFLib_EWRAP0(GMRFLib_compute_Qinv_TAUCS(p));
	}
		break;

	case GMRFLib_SMTP_PARDISO:
	{
		GMRFLib_pardiso_Qinv_INLA(p);
	}
		break;

	case GMRFLib_SMTP_STILES:
	{
		GMRFLib_stiles_Qinv_INLA(p);
	}
		break;

	default:
		GMRFLib_ERROR(GMRFLib_ESNH);
		break;
	}
	GMRFLib_LEAVE_FUNCTION;
	return GMRFLib_SUCCESS;
}

int GMRFLib_valid_smtp(int smtp)
{
	return ((smtp == GMRFLib_SMTP_BAND || smtp == GMRFLib_SMTP_TAUCS || smtp == GMRFLib_SMTP_PARDISO
		 || smtp == GMRFLib_SMTP_STILES) ? GMRFLib_TRUE : GMRFLib_FALSE);
}

const char *GMRFLib_reorder_name(GMRFLib_reorder_tp r)
{
	switch (r) {
	case GMRFLib_REORDER_DEFAULT:
		return "default";
	case GMRFLib_REORDER_IDENTITY:
		return "identity";
	case GMRFLib_REORDER_REVERSE_IDENTITY:
		return "reverseidentity";
	case GMRFLib_REORDER_BAND:
		return "band";
	case GMRFLib_REORDER_METIS:
		return "metis";
	case GMRFLib_REORDER_GENMMD:
		return "genmmd";
	case GMRFLib_REORDER_AMD:
		return "amd";
	case GMRFLib_REORDER_MD:
		return "md";
	case GMRFLib_REORDER_MMD:
		return "mmd";
	case GMRFLib_REORDER_AMDBAR:
		return "amdbar";
	case GMRFLib_REORDER_AMDC:
		return "amdc";
	case GMRFLib_REORDER_AMDBARC:
		return "amdbarc";
	case GMRFLib_REORDER_PARDISO:
		return "pardiso";
	case GMRFLib_REORDER_STILES:
		return "stiles";
	default:
		fprintf(stderr, "\n\t*** ERROR *** Reordering [%d] not defined.\n", r);
		GMRFLib_ASSERT_RETVAL(0 == 1, GMRFLib_EPARAMETER, "(unknown reording)");
	}

	return "(unknown reording)";
}

int GMRFLib_reorder_id(const char *name)
{
	if (!strcasecmp(name, "default") || !strcasecmp(name, "auto"))
		return GMRFLib_REORDER_DEFAULT;
	else if (!strcasecmp(name, "identity"))
		return GMRFLib_REORDER_IDENTITY;
	else if (!strcasecmp(name, "reverseidentity"))
		return GMRFLib_REORDER_REVERSE_IDENTITY;
	else if (!strcasecmp(name, "band"))
		return GMRFLib_REORDER_BAND;
	else if (!strcasecmp(name, "metis"))
		return GMRFLib_REORDER_METIS;
	else if (!strcasecmp(name, "genmmd"))
		return GMRFLib_REORDER_GENMMD;
	else if (!strcasecmp(name, "amd"))
		return GMRFLib_REORDER_AMD;
	else if (!strcasecmp(name, "amdbar"))
		return GMRFLib_REORDER_AMDBAR;
	else if (!strcasecmp(name, "md"))
		return GMRFLib_REORDER_MD;
	else if (!strcasecmp(name, "mmd"))
		return GMRFLib_REORDER_MMD;
	else if (!strcasecmp(name, "amdc"))
		return GMRFLib_REORDER_AMDC;
	else if (!strcasecmp(name, "amdbarc"))
		return GMRFLib_REORDER_AMDBARC;
	else if (!strcasecmp(name, "pardiso"))
		return GMRFLib_REORDER_PARDISO;
	else if (!strcasecmp(name, "stiles"))
		return GMRFLib_REORDER_STILES;
	else {
		fprintf(stderr, "\n\t*** ERROR *** Reordering [%s] not defined.\n", name);
		GMRFLib_ASSERT_RETVAL(0 == 1, GMRFLib_EPARAMETER, -1);
	}

	return -1;
}
