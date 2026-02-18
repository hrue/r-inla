#include <time.h>
#include <assert.h>
#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "GMRFLib/GMRFLib.h"

static GMRFLib_stiles_ctl_tp *ctl = NULL;
static GMRFLib_stiles_store_tp *store = NULL;
static GMRFLib_ptr_tp *free_ptrs = NULL;

#define RESCALE_GROUP() (store->ng2)

int GMRFLib_stiles_setup(GMRFLib_stiles_setup_tp *setup)
{
	GMRFLib_STOP_IF_NOT_SERIAL();
	GMRFLib_ENTER_FUNCTION;

	if (store) {
		GMRFLib_stiles_quit();
	}

	double tref = GMRFLib_timer();
	GMRFLib_ptr_tp *graphs = setup->graphs;
	GMRFLib_idx_tp *nrhss = setup->nrhss;
	int nt_outer = GMRFLib_openmp->max_threads_nested[0];
	int nt_inner = GMRFLib_openmp->max_threads_nested[1];
	int nt_special = GMRFLib_ADAPTIVE_NUM_THREADS();
	int nt_max_threads = GMRFLib_MAX_THREADS();
	int ng = graphs->n;
	int ng2 = 2 * ng;
	int ngt = ng2 + 1;

	assert(nt_outer > 0);
	assert(nt_inner > 0);
	assert(ng > 0);

	store = Calloc(1, GMRFLib_stiles_store_tp);
	store->n = Calloc(ngt, int);
	store->nnz = Calloc(ngt, int);

	// copy graphs
	store->graphs = NULL;
	for (int i = 0; i < ng2; i++) {
		GMRFLib_graph_tp *g = NULL;
		GMRFLib_graph_duplicate(&g, (GMRFLib_graph_tp *) (graphs->ptr[i % ng]));
		GMRFLib_ptr_add(&(store->graphs), g);
		store->n[i] = g->n;
		store->nnz[i] = g->nnz;
	}
	for (int i = 0; i < 1; i++) {
		GMRFLib_graph_tp *g = NULL;
		GMRFLib_graph_duplicate(&g, (GMRFLib_graph_tp *) (graphs->ptr[i]));
		GMRFLib_ptr_add(&(store->graphs), g);
		store->n[ng2] = g->n;
		store->nnz[ng2] = g->nnz;
	}

	int *calls_g = Malloc(ngt, int);
	int *cores_g = Malloc(ngt, int);
	int *zeros = Calloc(ngt, int);
	GMRFLib_ifill(ng, nt_outer, calls_g);
	GMRFLib_ifill(ng, nt_inner * nt_outer, cores_g);       /* yes, the total number of threads for one group of matrices */
	GMRFLib_ifill(ng, 1, calls_g + ng);
	GMRFLib_ifill(ng, nt_special * 1, cores_g + ng);
	GMRFLib_ifill(1, nt_max_threads, calls_g + ng2);
	GMRFLib_ifill(1, 1 * nt_max_threads, cores_g + ng2);

	bool *inv = Malloc(ngt, bool);
	GMRFLib_bfill(ng2, true, inv);
	GMRFLib_bfill(1, false, inv + ng2);

	int nn = 3 + (nrhss ? nrhss->n : 0);
	assert(nn > 0);
	int *nrhs = Calloc(nn, int);
	nrhs[0] = nn;					       /* the first element is the number of different size of rhs's */
	nrhs[1] = 1;					       /* always */
	nrhs[2] = IMAX(1, GMRFLib_stiles_get_block_size());    /* always */

	if (nrhss && nrhss->n) {
		Memcpy(nrhs + 3, nrhss->idx, nrhss->n * sizeof(int));
		GMRFLib_sort_i(nrhs + 1, nn - 1);
		int nnew = 0, *nrhsnew = NULL;
		GMRFLib_iuniques(&nnew, &nrhsnew, nrhs + 1, nn - 1);
		GMRFLib_sort_i(nrhsnew, nnew);
		nrhs[0] = nnew;
		Memcpy(nrhs + 1, nrhsnew, nnew * sizeof(int));
	}
	store->nrhss = nrhs[0];
	store->rhss = nrhs + 1;

	int *rescale = Calloc(ngt, int);
	rescale[ng2] = 1;
	sTiles_set_rescale_cores(rescale, ngt); 

	//FIXME("FIX LATER");
	//sTiles_set_tile_type_mode(0);

	sTiles_create(&(store->obj), ngt, calls_g, cores_g, zeros, inv);
	store->ng = ng;
	store->ng2 = ng2;
	store->ngt = ngt;
	store->n_in_group = ngt;
	store->n_within_group = calls_g;
	store->n_cores_group = cores_g;
	store->nt_outer = nt_outer;
	store->nt_inner = nt_inner;
	store->nt_special = nt_special;
	store->rescale_on = 0;
	store->nt_max_threads = nt_max_threads;
	store->Qinv_done = Malloc(ngt, bool *);
	store->bind_done = Malloc(ngt, bool *);
	store->chol_done = Malloc(ngt, bool *);
	for (int i = 0; i < ngt; i++) {
		store->Qinv_done[i] = Malloc(store->n_within_group[i], bool);
		store->bind_done[i] = Malloc(store->n_within_group[i], bool);
		store->chol_done[i] = Malloc(store->n_within_group[i], bool);
		GMRFLib_bfill(store->n_within_group[i], false, store->Qinv_done[i]);
		GMRFLib_bfill(store->n_within_group[i], false, store->bind_done[i]);
		GMRFLib_bfill(store->n_within_group[i], false, store->chol_done[i]);
	}

	// prepare for parallel, but it seems not worth it at the moment.
	// anyway, we build first, then store
	int **sidx_i = Calloc(ng, int *);
	int **sidx_j = Calloc(ng, int *);

	for (int ig = 0; ig < ng; ig++) {
		GMRFLib_graph_tp *g = (GMRFLib_graph_tp *) (graphs->ptr[ig]);
		int nz = g->n + g->nnz / 2;
		int *idx_i = Malloc(nz, int);
		int *idx_j = Malloc(nz, int);
		int k = 0;
		for (int i = 0; i < g->n; i++) {
			idx_i[k] = i;
			idx_j[k++] = i;
			int m = g->lnnbs[i];
			if (m) {
				GMRFLib_ifill(m, i, idx_i + k);
				Memcpy(idx_j + k, g->lnbs[i], m * sizeof(int));
				k += m;
			}
		}
		assert(k == nz);
		sidx_i[ig] = idx_i;
		sidx_j[ig] = idx_j;
	}

	// ptr's are stored, contents not copied, so we need to free them later.
	// since sTiles wants the lower triangular matrix, we just swap 'idx_i' and 'idx_j'
	for (int ig = 0; ig < ng; ig++) {
		GMRFLib_graph_tp *g = (GMRFLib_graph_tp *) (graphs->ptr[ig]);
		int nz = g->n + g->nnz / 2;
		sTiles_assign_graph(ig, &(store->obj), g->n, nz, sidx_j[ig], sidx_i[ig]);	/* oops, yes we swap */
		sTiles_assign_graph(ig + ng, &(store->obj), g->n, nz, sidx_j[ig], sidx_i[ig]);	/* oops, yes we swap */
		GMRFLib_ptr_add(&free_ptrs, sidx_i[ig]);
		GMRFLib_ptr_add(&free_ptrs, sidx_j[ig]);

		if (ig == 0) {
			sTiles_assign_graph(ng2, &(store->obj), g->n, nz, sidx_j[ig], sidx_i[ig]);
		}
	}
	Free(sidx_i);
	Free(sidx_j);

	sTiles_init(&(store->obj));

	store->perm = Calloc(ngt, int *);
	store->iperm = Calloc(ngt, int *);
	for (int i = 0; i < ng2; i++) {
		int *p = sTiles_return_perm_vec(i, &(store->obj));
		int *pi = sTiles_return_iperm_vec(i, &(store->obj));
		store->perm[i] = Malloc(store->n[i], int);
		store->iperm[i] = Malloc(store->n[i], int);
		Memcpy(store->perm[i], p, store->n[i] * sizeof(int));
		Memcpy(store->iperm[i], pi, store->n[i] * sizeof(int));

		if (i == 0) {
			store->perm[ng2] = Malloc(store->n[i], int);
			store->iperm[ng2] = Malloc(store->n[i], int);
			Memcpy(store->perm[ng2], p, store->n[i] * sizeof(int));
			Memcpy(store->iperm[ng2], pi, store->n[i] * sizeof(int));
		}
	}

	Free(zeros);
	Free(inv);
	store->wtime = GMRFLib_timer() - tref;

	if (ctl->verbose) {
		GMRFLib_stiles_print(stdout);
	}

	GMRFLib_LEAVE_FUNCTION;
	return GMRFLib_SUCCESS;
}

void GMRFLib_stiles_rescale_start(void) 
{
	if (store) {
		sTiles_turn_on_rescale(0, &(store->obj));
		store->rescale_on = 1;
	}
}

void GMRFLib_stiles_rescale_end(void) 
{
	if (store) {
		sTiles_turn_off_rescale(0, &(store->obj));
		store->rescale_on = 0;
	}
}

int GMRFLib_stiles_is_rescale(void) 
{
	return (store && store->rescale_on ? 1 : 0);
}

int GMRFLib_stiles_rescale_group(void) 
{
	if (store) {
		assert(GMRFLib_stiles_is_rescale());
		return RESCALE_GROUP();
	} else {
		return 0;
	}
}

void GMRFLib_stiles_quit(void)
{
	// not sure if all is free'ed here...
	
	if (!store) {
		return;
	}

	GMRFLib_STOP_IF_NOT_SERIAL();

	GMRFLib_stiles_unbind_all();
	sTiles_quit();

	if (free_ptrs) {
		for (int i = 0; i < free_ptrs->n; i++) {
			Free(free_ptrs->ptr[i]);
		}
		GMRFLib_ptr_free(free_ptrs);
		free_ptrs = NULL;
	}
	GMRFLib_ptr_free(store->graphs);

	if (store->perm) {
		for (int i = 0; i < store->n_in_group; i++) {
			Free(store->perm[i]);
		}
	}
	if (store->iperm) {
		for (int i = 0; i < store->n_in_group; i++) {
			Free(store->iperm[i]);
		}
	}

	if (store->Qinv_done) {
		for (int i = 0; i < store->n_in_group; i++) {
			Free(store->Qinv_done[i]);
		}
		Free(store->Qinv_done);
	}

	if (store->bind_done) {
		for (int i = 0; i < store->n_in_group; i++) {
			Free(store->bind_done[i]);
		}
		Free(store->bind_done);
	}

	if (store->chol_done) {
		for (int i = 0; i < store->n_in_group; i++) {
			Free(store->chol_done[i]);
		}
		Free(store->chol_done);
	}

	Free(store->perm);
	Free(store->iperm);
	Free(store->n_within_group);
	Free(store->n_cores_group);
	Free(store->n);
	Free(store->nnz);
	free(store->rhss - 1);				       /* yes, since we shift when creating */
	store->rhss = NULL;
	Free(store);
}

void GMRFLib_stiles_print_idx(GMRFLib_stiles_idx_tp *stiles_idx, FILE *fp)
{
	fprintf(fp, "%s:%1d stiles_idx: in_group=%1d within_group=%1d nrhs=%1d\n",
		__FILE__, __LINE__, stiles_idx->in_group, stiles_idx->within_group, stiles_idx->nrhs);
}

int GMRFLib_stiles_set_idx(GMRFLib_stiles_idx_tp *stiles_idx, int nrhs)
{
	if (!store) {
		return GMRFLib_SUCCESS;
	}
	
	if (GMRFLib_smtp == GMRFLib_SMTP_STILES) {
		if (GMRFLib_OPENMP_IN_SERIAL()) {
			if (stiles_idx->in_group < store->ng) {
				stiles_idx->in_group += store->ng;
			}
			stiles_idx->within_group = 0;
		} else {
			stiles_idx->within_group = omp_get_thread_num();
		}
	}
	stiles_idx->nrhs = nrhs;
	return GMRFLib_SUCCESS;
}

void *GMRFLib_stiles_get_store_ptr(void)
{
	return (void *) store;
}

void GMRFLib_stiles_print(FILE *fp)
{
	if (!store)
		return;

#pragma omp critical (Name_4c8dac87b14702b8de3511c972d6b27af33cc04c)
	{
		fprintf(fp, "\n\ncontent of 'store' (computed in %.3fs):\n", store->wtime);
		fprintf(fp, "\tngroup[%1d] verbose[%1d] ng[%1d] ng2[%1d] ngt[%1d]\n", store->n_in_group, ctl->verbose,
			store->ng, store->ng2, store->ngt);
		fprintf(fp, "\tnt_outer[%1d] nt_inner[%1d] nt_special[%1d] block.size[%1d]\n", store->nt_outer, store->nt_inner, store->nt_special,
			GMRFLib_stiles_get_block_size());
		// GMRFLib_stiles_print_ctl_param(fp, "\t\t");
		for (int i = 0; i < store->n_in_group; i++) {
			fprintf(fp, "\tgroup[%1d]: n[%1d] nnz[%1d] n_within_group[%1d] n_cores_group[%1d] rescaled_group[%s]\n",
				i, store->n[i], store->nnz[i], store->n_within_group[i], store->n_cores_group[i], (i == RESCALE_GROUP() ? "yes" : "no"));

			fprintf(fp, "\t\tnrhs = [ ");
			for (int j = 0; j < store->nrhss; j++) {
				fprintf(fp, "%1d ", store->rhss[j]);
			}
			fprintf(fp, "]\n");
			int perm_identity = 1, preview = 8;
			for (int j = 0; j < store->n[i]; j++) {
				if (store->perm[i]) {
					if (store->perm[i][j] != j || store->iperm[i][j] != j) {
						fprintf(fp, "\t\tperm[%1d][%1d] = %5d iperm[%1d][%1d] = %5d\n", i, j,
							store->perm[i][j], i, j, store->iperm[i][j]);
						preview--;
						perm_identity = 0;
					}
					if (preview <= 0)
						break;
				}
			}

			if (store->perm[i]) {
				if (perm_identity) {
					fprintf(fp, "\t\tperm[%1d] = identity\n", i);
				}
			}

			printf("\t\tQinv_done: ");
			for (int j = 0; j < store->n_within_group[i]; j++) {
				printf("%1d:%1d ", j, store->Qinv_done[i][j]);
			}
			printf("\n");

			printf("\t\tbind_done: ");
			for (int j = 0; j < store->n_within_group[i]; j++) {
				printf("%1d:%1d ", j, store->bind_done[i][j]);
			}
			printf("\n");

			printf("\t\tchol_done: ");
			for (int j = 0; j < store->n_within_group[i]; j++) {
				printf("%1d:%1d ", j, store->chol_done[i][j]);
			}
			printf("\n");
		}
		fprintf(fp, "\n");
	}
}

int *GMRFLib_stiles_get_perm(GMRFLib_stiles_idx_tp *stiles_idx)
{
	return (store ? store->perm[stiles_idx->in_group] : NULL);
}

int *GMRFLib_stiles_get_iperm(GMRFLib_stiles_idx_tp *stiles_idx)
{
	return (store ? store->iperm[stiles_idx->in_group] : NULL);
}

int GMRFLib_stiles_set_ctl(int verbose, int block_size, int len, int *param)
{
	GMRFLib_STOP_IF_NOT_SERIAL();

	if (ctl) {
		Free(ctl->param);
		Free(ctl);
	}
	ctl = Calloc(1, GMRFLib_stiles_ctl_tp);
	ctl->verbose = (verbose >= 0 ? verbose : 0);
	ctl->block_size = (block_size > 0 ? block_size : 32);
	if (len > 0) {
		ctl->param_len = len;
		ctl->param = Malloc(len, int);
		Memcpy(ctl->param, param, len * sizeof(int));
	} else {
		ctl->param_len = 32;
		ctl->param = Malloc(ctl->param_len, int);
		GMRFLib_ifill(ctl->param_len, -1, ctl->param);
	}

	sTiles_expert_user();
	for(int i = 0; i < ctl->param_len; i++) {
		sTiles_set_control_param(i, ctl->param[i]);
	}
	for(int i = 0; i < ctl->param_len; i++) {
		ctl->param[i] = sTiles_get_control_param(i);
	}

	return GMRFLib_SUCCESS;
}

GMRFLib_stiles_setup_tp *GMRFLib_stiles_get_setup(void *mb)
{
	GMRFLib_STOP_IF_NOT_SERIAL();

	// need to call it in the inla_ parts as need need access to inla_tp. it is a mess to make 'inla_tp' available here
	GMRFLib_stiles_setup_tp *inla_stiles_get_setup(void *);
	return inla_stiles_get_setup(mb);
}

void GMRFLib_stiles_free_setup(GMRFLib_stiles_setup_tp *setup)
{
	if (!setup)
		return;

	GMRFLib_ptr_free(setup->graphs);
	GMRFLib_idx_free(setup->nrhss);
	Free(setup);
}

int GMRFLib_stiles_chol(GMRFLib_stiles_idx_tp *stiles_idx)
{
#if 0
	// FIXME("CHOL ENTER");
	double tref = -GMRFLib_timer();
#endif
	assert(!GMRFLib_stiles_is_rescale());

	GMRFLib_stiles_idx_tp lidx = *stiles_idx;
	lidx.within_group = omp_get_thread_num();

	GMRFLib_stiles_bind(&lidx);
	int status = sTiles_chol(lidx.in_group, lidx.within_group, &(store->obj));
	GMRFLib_stiles_unbind(&lidx);

	if (status) {
		fprintf(stderr, "\n\n*** ERROR *** sTiles_chol %d \n\n", status);
		fflush(stderr);
		store->chol_done[lidx.in_group][lidx.within_group] = false;
	} else {
		store->chol_done[lidx.in_group][lidx.within_group] = true;
	}
#if 0
//#pragma omp critical (Name_a59d65352b63a2cd6aac7d155e2f7f307080c4d0)
	{
		tref += GMRFLib_timer();
		printf("CHOL LEAVE thread %d num_threads %d time %f\n", omp_get_thread_num(), omp_get_num_threads(), tref);
	}
#endif

	return (status ? GMRFLib_EPOSDEF : GMRFLib_SUCCESS);
}

double GMRFLib_stiles_logdet(GMRFLib_stiles_idx_tp *stiles_idx)
{
	GMRFLib_stiles_idx_tp lidx = *stiles_idx;
	lidx.within_group = omp_get_thread_num();

	assert(!GMRFLib_stiles_is_rescale());
	assert(store->chol_done[lidx.in_group][lidx.within_group] == true);

	GMRFLib_stiles_bind(&lidx);
	double ldet = sTiles_get_logdet(lidx.in_group, lidx.within_group, &(store->obj));
	GMRFLib_stiles_unbind(&lidx);

	return ldet;
}

void GMRFLib_stiles_Qinv(GMRFLib_stiles_idx_tp *stiles_idx)
{
	GMRFLib_stiles_idx_tp lidx = *stiles_idx;
	lidx.within_group = omp_get_thread_num();
	assert(!GMRFLib_stiles_is_rescale());
	assert(store->chol_done[lidx.in_group][lidx.within_group] == true);

	GMRFLib_stiles_bind(&lidx);
	if (store->Qinv_done[lidx.in_group][lidx.within_group]) {
		sTiles_clear_selinv(lidx.in_group, lidx.within_group, &(store->obj));
	}

	sTiles_selinv(lidx.in_group, lidx.within_group, &(store->obj));
	GMRFLib_stiles_unbind(&lidx);
	store->Qinv_done[lidx.in_group][lidx.within_group] = true;
}

int GMRFLib_stiles_build(GMRFLib_stiles_idx_tp *stiles_idx, int thread_id, GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_arg)
{
	GMRFLib_stiles_idx_tp lidx = *stiles_idx;
	lidx.within_group = omp_get_thread_num();
	// return a malloc'ed vector of Qi. values
	assert(LEGAL(lidx.in_group, store->n_in_group));
	assert(LEGAL(lidx.within_group, store->n_within_group[lidx.in_group]));

	GMRFLib_graph_tp *graph = (GMRFLib_graph_tp *) (store->graphs->ptr[lidx.in_group]);
	int n = graph->n;

	int N = graph->n + graph->nnz / 2;
	double *x = Malloc(N, double);

	GMRFLib_tabulate_Qfunc_arg_tp *arg = (GMRFLib_tabulate_Qfunc_arg_tp *) Qfunc_arg;
	int fast_copy = (Qfunc == GMRFLib_tabulate_Qfunction_std && arg->Q);

	if (fast_copy) {
		Memcpy(x, arg->Q->a, N * sizeof(double));
	} else {
		double *values = Calloc(1 + graph->max_nnbs, double);
		if (ISNAN(Qfunc(thread_id, 0, -1, values, Qfunc_arg))) {
			for (int i = 0, k = 0; i < n; i++) {
				x[k++] = Qfunc(thread_id, i, i, NULL, Qfunc_arg);
				for (int jj = 0; jj < graph->lnnbs[i]; jj++) {
					int j = graph->lnbs[i][jj];
					x[k++] = Qfunc(thread_id, i, j, NULL, Qfunc_arg);
				}
			}
		} else {
			for (int i = 0, k = 0; i < n; i++) {
				Qfunc(thread_id, i, -1, values, Qfunc_arg);
				int len = 1 + graph->lnnbs[i];
				Memcpy(x + k, values, len * sizeof(double));
				k += len;
			}
		}
		Free(values);
	}

	sTiles_assign_values(lidx.in_group, lidx.within_group, &(store->obj), x);
	Free(x);

	store->Qinv_done[lidx.in_group][lidx.within_group] = false;
	store->chol_done[lidx.in_group][lidx.within_group] = false;

	return GMRFLib_SUCCESS;
}

#define CHOL_DONE_CHECK(msg_)						\
	if (!store->chol_done[llidx.in_group][llidx.within_group]) {	\
		printf(" *** ERROR *** solve %s in_group=%1d within_group=%1d BUT chol_done=%1d\n", \
		       #msg_,  llidx.in_group, llidx.within_group,	\
		       store->chol_done[llidx.in_group][llidx.within_group]); \
		fflush(stdout);						\
		assert(store->chol_done[llidx.in_group][llidx.within_group]); \
	}

int GMRFLib_stiles_solve_LLT(GMRFLib_stiles_idx_tp *stiles_idx, double *rhs)
{
	assert(stiles_idx->nrhs);
	GMRFLib_stiles_idx_tp lidx = *stiles_idx;
	lidx.within_group = omp_get_thread_num();
	GMRFLib_stiles_idx_tp llidx = lidx;
	if (GMRFLib_stiles_is_rescale() && (lidx.in_group == RESCALE_GROUP())) {
		llidx.in_group =  0;
		llidx.within_group =  0;
	}

	CHOL_DONE_CHECK(LLT);
	
#if 0
	static double tref = 0;
#pragma omp threadprivate(tref)
	static double trefc = 0;
#pragma omp threadprivate(trefc)
	if (stiles_idx->nrhs > 1) {
		tref += -GMRFLib_timer();
	}
#endif
	
	GMRFLib_stiles_bind(&lidx);
	if (GMRFLib_stiles_is_rescale()) {
		sTiles_solve_LLT_rescale(llidx.in_group, llidx.within_group, &(store->obj), rhs, llidx.nrhs, lidx.in_group, lidx.within_group);
	} else {
		sTiles_solve_LLT(llidx.in_group, llidx.within_group, &(store->obj), rhs, llidx.nrhs);
	}
	GMRFLib_stiles_unbind(&lidx);

#if 0
	if (stiles_idx->nrhs > 1) {
		tref += GMRFLib_timer();
		trefc += stiles_idx->nrhs;
		printf("solve %1d rhs using %.6f * E-6 each\n", (int) trefc, 1.0E6 * tref / trefc);
	}
#endif
	return GMRFLib_SUCCESS;
}

int GMRFLib_stiles_solve_L(GMRFLib_stiles_idx_tp *stiles_idx, double *rhs)
{
	GMRFLib_stiles_idx_tp lidx = *stiles_idx;
	lidx.within_group = omp_get_thread_num();
	GMRFLib_stiles_idx_tp llidx = lidx;
	if (GMRFLib_stiles_is_rescale() && (lidx.in_group == RESCALE_GROUP())) {
		llidx.in_group =  0;
		llidx.within_group =  0;
	}
	
	CHOL_DONE_CHECK(L);
	
	GMRFLib_stiles_bind(&lidx);
	sTiles_solve_L(llidx.in_group, llidx.within_group, &(store->obj), rhs, llidx.nrhs);
	GMRFLib_stiles_unbind(&lidx);

	return GMRFLib_SUCCESS;
}

int GMRFLib_stiles_solve_LT(GMRFLib_stiles_idx_tp *stiles_idx, double *rhs)
{
	GMRFLib_stiles_idx_tp lidx = *stiles_idx;
	lidx.within_group = omp_get_thread_num();
	GMRFLib_stiles_idx_tp llidx = lidx;
	if (GMRFLib_stiles_is_rescale() && (lidx.in_group == RESCALE_GROUP())) {
		llidx.in_group =  0;
		llidx.within_group =  0;
	}
	
	CHOL_DONE_CHECK(LT);

	GMRFLib_stiles_bind(&lidx);
	sTiles_solve_LT(llidx.in_group, llidx.within_group, &(store->obj), rhs, llidx.nrhs);
	GMRFLib_stiles_unbind(&lidx);

	return GMRFLib_SUCCESS;
}
#undef CHOL_DONE_CHECK

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((target_clones(INLA_CLONE_TARGETS "default")))
int GMRFLib_stiles_Qinv_INLA(GMRFLib_problem_tp *problem)
{
#define QINV_GET(i_, j_) sTiles_get_selinv_elm(problem->stiles_idx->in_group, \
					       problem->stiles_idx->within_group, \
					       i_, j_, &(store->obj))
	assert(!GMRFLib_stiles_is_rescale());

	if (problem == NULL) {
		return GMRFLib_SUCCESS;
	}

	int in_group = problem->stiles_idx->in_group;
	int n = store->n[in_group];

	GMRFLib_stiles_Qinv(problem->stiles_idx);
	GMRFLib_graph_tp *graph = (GMRFLib_graph_tp *) (store->graphs->ptr[in_group]);
	map_id **Qinv = Calloc(n, map_id *);

	for (int i = 0; i < n; i++) {
		int nnb = graph->lnnbs[i];
		Qinv[i] = Calloc(1, map_id);
		map_id_init_hint(Qinv[i], 1 + nnb);
		map_id_set(Qinv[i], i, QINV_GET(i, i));
		for (int jj = 0; jj < nnb; jj++) {
			int j = graph->lnbs[i][jj];
			map_id_set(Qinv[i], j, QINV_GET(i, j));
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
#pragma omp simd
	for (int i = 0; i < n; i++) {
		subQinv->mapping[i] = i;
	}
	problem->sub_inverse = subQinv;

#undef QINV_GET
	return GMRFLib_SUCCESS;
}
#pragma GCC diagnostic pop

void GMRFLib_stiles_bind(GMRFLib_stiles_idx_tp *stiles_idx)
{
	if (!store)
		return;
	
	GMRFLib_ENTER_FUNCTION;
	bool *p = &(store->bind_done[stiles_idx->in_group][stiles_idx->within_group]);
	if (*p == false) {
		sTiles_bind(stiles_idx->in_group, stiles_idx->within_group, &(store->obj));
		*p = true;
	}
	GMRFLib_LEAVE_FUNCTION;
}

void GMRFLib_stiles_unbind(GMRFLib_stiles_idx_tp *stiles_idx)
{
	if (!store)
		return;

	GMRFLib_ENTER_FUNCTION;
	bool *p = &(store->bind_done[stiles_idx->in_group][stiles_idx->within_group]);
	if (*p == true) {
		sTiles_unbind(stiles_idx->in_group, stiles_idx->within_group, &(store->obj));
		*p = false;
	}
	GMRFLib_LEAVE_FUNCTION;
}

void GMRFLib_stiles_unbind_group(int in_group)
{
	if (!store)
		return;

	// need to do this one in parallel
#pragma omp parallel for num_threads(store->n_within_group[in_group])
	for (int j = 0; j < store->n_within_group[in_group]; j++) {
		GMRFLib_stiles_idx_tp stiles_idx = { in_group, j, 0 };
		GMRFLib_stiles_unbind(&stiles_idx);
	}
}

void GMRFLib_stiles_unbind_all(void)
{
	if (!store)
		return;

	for (int g = 0; g < store->n_in_group; g++) {
		GMRFLib_stiles_unbind_group(g);
	}
}

int GMRFLib_stiles_get_verbose()
{
	return (ctl ? ctl->verbose : 0);
}

GMRFLib_stiles_ctl_tp * GMRFLib_stiles_get_ctl(void) 
{
	return ctl;
}

int GMRFLib_stiles_get_tile_size(void)
{
	return (ctl && ctl->param[1] > 0 ? ctl->param[1] : sTiles_get_auto_tile_size(0));
}

int GMRFLib_stiles_get_block_size(void)
{
	return (ctl && ctl->block_size > 0 ? ctl->block_size : 32);
}

void GMRFLib_stiles_print_ctl_param(FILE *UNUSED(fp), char *UNUSED(suf)) 
{
	sTiles_print_params();
}

