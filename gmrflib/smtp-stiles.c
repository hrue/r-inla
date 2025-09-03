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
#include "GMRFLib/GMRFLibP.h"

static GMRFLib_stiles_ctl_tp *ctl = NULL;
static GMRFLib_stiles_store_tp *store = NULL;
static GMRFLib_ptr_tp *free_ptrs = NULL;

int GMRFLib_stiles_setup(GMRFLib_stiles_setup_tp *setup)
{
	double tref = GMRFLib_timer();
	GMRFLib_STOP_IF_NOT_SERIAL();

	GMRFLib_ptr_tp *graphs = setup->graphs;
	GMRFLib_idx_tp *nrhss = setup->nrhss;
	if (store) {
		GMRFLib_stiles_quit();
	}

	if (!ctl) {
		GMRFLib_stiles_set_ctl(0, 0);
	}

	int nt_outer = GMRFLib_openmp->max_threads_nested[0];
	int nt_inner = GMRFLib_openmp->max_threads_nested[1];
	int nt_special = GMRFLib_ADAPTIVE_NUM_THREADS();
	int ng = graphs->n;
	int ng2 = 2 * ng;

	assert(nt_outer > 0);
	assert(nt_inner > 0);
	assert(ng > 0);

	store = Calloc(1, GMRFLib_stiles_store_tp);
	store->n = Calloc(ng2, int);
	store->nnz = Calloc(ng2, int);

	// copy graphs
	store->graphs = NULL;
	for (int i = 0; i < ng2; i++) {
		GMRFLib_graph_tp *g = NULL;
		GMRFLib_graph_duplicate(&g, (GMRFLib_graph_tp *) (graphs->ptr[i % ng]));
		GMRFLib_ptr_add(&(store->graphs), g);
		store->n[i] = g->n;
		store->nnz[i] = g->nnz;
	}

	int *calls_g = Malloc(ng2, int);
	int *cores_g = Malloc(ng2, int);
	int *zeros = Calloc(ng2, int);
	GMRFLib_ifill(ng, nt_outer, calls_g);
	GMRFLib_ifill(ng, nt_inner * nt_outer, cores_g);       /* yes, the total number of threads for one group of matrices */
	GMRFLib_ifill(ng, 1, calls_g + ng);
	GMRFLib_ifill(ng, nt_special * 1, cores_g + ng);

	bool *inv = Malloc(ng2, bool);
	GMRFLib_bfill(ng2, true, inv);

	int nn = 3 + (nrhss ? nrhss->n : 0);
	assert(nn > 0);
	int *nrhs = Calloc(nn, int);
	nrhs[0] = nn;					       /* the first element is the number of different size of rhs's */
	nrhs[1] = 1;					       /* always */
	nrhs[2] = IMAX(1, GMRFLib_stiles_get_tile_size());     /* always */

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

	sTiles_create(&(store->obj), ng2, calls_g, cores_g, zeros, inv, nrhs);
	store->ng = ng;
	store->ng2 = ng2;
	store->n_in_group = ng2;
	store->n_within_group = calls_g;
	store->n_cores_group = cores_g;
	store->nt_outer = nt_outer;
	store->nt_inner = nt_inner;
	store->nt_special = nt_special;
	store->Qinv_done = Malloc(ng2, bool *);
	store->bind_done = Malloc(ng2, bool *);
	for (int i = 0; i < ng2; i++) {
		store->Qinv_done[i] = Malloc(store->n_within_group[i], bool);
		store->bind_done[i] = Malloc(store->n_within_group[i], bool);
		GMRFLib_bfill(store->n_within_group[i], false, store->Qinv_done[i]);
		GMRFLib_bfill(store->n_within_group[i], false, store->bind_done[i]);
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
	}
	Free(sidx_i);
	Free(sidx_j);

	sTiles_init(&(store->obj));

	store->perm = Calloc(ng2, int *);
	store->iperm = Calloc(ng2, int *);
	for (int i = 0; i < ng2; i++) {
		int *p = sTiles_return_perm_vec(i, &(store->obj));
		int *pi = sTiles_return_iperm_vec(i, &(store->obj));
		store->perm[i] = Malloc(store->n[i], int);
		store->iperm[i] = Malloc(store->n[i], int);
		Memcpy(store->perm[i], p, store->n[i] * sizeof(int));
		Memcpy(store->iperm[i], pi, store->n[i] * sizeof(int));
	}

	Free(zeros);
	Free(inv);
	store->wtime = GMRFLib_timer() - tref;

	if (ctl->verbose) {
		GMRFLib_stiles_print(stdout);
	}

	return GMRFLib_SUCCESS;
}

void GMRFLib_stiles_quit(void)
{
	if (!store)
		return;

	GMRFLib_STOP_IF_NOT_SERIAL();

	// int nt_outer = store->nt_outer;
	// int nt_inner = store->nt_inner;

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

int GMRFLib_stiles_set_idx_copy(GMRFLib_stiles_idx_tp *stiles_idx, int nrhs)
{
	return GMRFLib_stiles_set_idx(stiles_idx, nrhs);
}

int GMRFLib_stiles_set_idx_special(GMRFLib_stiles_idx_tp *stiles_idx, int nrhs)
{
	GMRFLib_stiles_idx_tp lidx;
	Memcpy(&lidx, stiles_idx, sizeof(GMRFLib_stiles_idx_tp));
	lidx.in_group += store->ng;
	return GMRFLib_stiles_set_idx(&lidx, nrhs);
}

int GMRFLib_stiles_set_idx(GMRFLib_stiles_idx_tp *stiles_idx, int nrhs)
{

	// rewrite ->within_group using omp_get_thread_num(), keep in_group fixed
	if (!store)
		return GMRFLib_SUCCESS;

	int nt = omp_get_thread_num();
	if (GMRFLib_smtp == GMRFLib_SMTP_STILES) {
		if (GMRFLib_OPENMP_IN_SERIAL()) {
			if (stiles_idx->in_group < store->ng) {
				stiles_idx->in_group += store->ng;
			}
			assert(nt == 0);
			stiles_idx->within_group = 0;
		} else {
			stiles_idx->within_group = nt;
		}
	}

	stiles_idx->within_group = (nt % store->n_within_group[stiles_idx->in_group]);
	stiles_idx->nrhs = nrhs;

	return GMRFLib_SUCCESS;
}

void *GMRFLib_stiles_get_store_ptr(void)
{
	return (void *) store;
}

void GMRFLib_stiles_print(FILE *fp)
{
#pragma omp critical (Name_4c8dac87b14702b8de3511c972d6b27af33cc04c)
	{
		fprintf(fp, "\n\ncontent of 'store' (computed in %.3fs):\n", store->wtime);
		fprintf(fp, "\tngroup[%1d] verbose[%1d] tile.size[%1d] ng[%1d] ng2[%1d]\n", store->n_in_group, ctl->verbose,
			GMRFLib_stiles_get_tile_size(), store->ng, store->ng2);
		fprintf(fp, "\tnt_outer[%1d] nt_inner[%1d] nt_special[%1d]\n", store->nt_outer, store->nt_inner, store->nt_special);

		for (int i = 0; i < store->n_in_group; i++) {
			fprintf(fp, "\tgroup[%1d]: n[%1d] nnz[%1d] n_within_group[%1d] n_cores_group[%1d]\n",
				i, store->n[i], store->nnz[i], store->n_within_group[i], store->n_cores_group[i]);

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
		}
		fprintf(fp, "\n");
	}
}

int *GMRFLib_stiles_get_perm(GMRFLib_stiles_idx_tp *stiles_idx)
{
	return store->perm[stiles_idx->in_group];
}

int *GMRFLib_stiles_get_iperm(GMRFLib_stiles_idx_tp *stiles_idx)
{
	return store->iperm[stiles_idx->in_group];
}

int GMRFLib_stiles_set_ctl(int verbose, int tile_size)
{
	GMRFLib_STOP_IF_NOT_SERIAL();
	Free(ctl);

	ctl = Calloc(1, GMRFLib_stiles_ctl_tp);
	ctl->verbose = (verbose >= 0 ? verbose : 0);
	ctl->tile_size = IMAX(0, tile_size);
#if defined(INLA_WITH_STILES)
	if (ctl->tile_size == 0) {
		ctl->tile_size = get_auto_tile_size();
	}
#endif
	if (ctl->tile_size > 0) {
		sTiles_set_tile_size(ctl->tile_size);
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
	if (setup) {
		GMRFLib_ptr_free(setup->graphs);
		GMRFLib_idx_free(setup->nrhss);
		Free(setup);
	}
}

int GMRFLib_stiles_chol(GMRFLib_stiles_idx_tp *stiles_idx)
{
#if 0
	// FIXME("CHOL ENTER");
	double tref = -GMRFLib_timer();
#endif

	GMRFLib_stiles_idx_tp lidx;
	Memcpy(&lidx, stiles_idx, sizeof(GMRFLib_stiles_idx_tp));
	lidx.within_group = omp_get_thread_num();

	GMRFLib_stiles_bind(&lidx);
	int status = sTiles_chol(lidx.in_group, lidx.within_group, &(store->obj));
	GMRFLib_stiles_unbind(&lidx);

	if (status) {
		fprintf(stderr, "\n\n*** ERROR *** sTiles_chol %d \n\n", status);
		fflush(stderr);
	}
#if 0
//#pragma omp critical (Name_a59d65352b63a2cd6aac7d155e2f7f307080c4d0)
	{
		tref += GMRFLib_timer();
		printf("CHOL LEAVE thread %d num_threads %d time %f\n", omp_get_thread_num(), omp_get_num_threads(), tref);
	}
#endif

	return (status ? !GMRFLib_SUCCESS : GMRFLib_SUCCESS);
}

double GMRFLib_stiles_logdet(GMRFLib_stiles_idx_tp *stiles_idx)
{
	GMRFLib_stiles_idx_tp lidx;
	Memcpy(&lidx, stiles_idx, sizeof(GMRFLib_stiles_idx_tp));
	lidx.within_group = omp_get_thread_num();

	GMRFLib_stiles_bind(&lidx);
	double ldet = sTiles_get_logdet(lidx.in_group, lidx.within_group, &(store->obj));
	GMRFLib_stiles_unbind(&lidx);

	return ldet;
}

void GMRFLib_stiles_Qinv(GMRFLib_stiles_idx_tp *stiles_idx)
{
	GMRFLib_stiles_idx_tp lidx;
	Memcpy(&lidx, stiles_idx, sizeof(GMRFLib_stiles_idx_tp));
	lidx.within_group = omp_get_thread_num();

	GMRFLib_stiles_bind(&lidx);
	if (store->Qinv_done[lidx.in_group][lidx.within_group]) {
		sTiles_clear_selinv(lidx.in_group, lidx.within_group, &(store->obj));
	} else {
		store->Qinv_done[lidx.in_group][lidx.within_group] = true;
	}
	sTiles_selinv(lidx.in_group, lidx.within_group, &(store->obj));
	GMRFLib_stiles_unbind(&lidx);
}

double GMRFLib_stiles_Qinv_get(int i, int j, GMRFLib_stiles_idx_tp *stiles_idx)
{
	int in_group = stiles_idx->in_group;
	int within_group = omp_get_thread_num();
	return sTiles_get_selinv_elm(in_group, within_group, i, j, &(store->obj));
}

int GMRFLib_stiles_build(GMRFLib_stiles_idx_tp *stiles_idx, int thread_id, GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_arg)
{
	int in_group = stiles_idx->in_group;
	int within_group = omp_get_thread_num();

	// return a malloc'ed vector of Qi. values

	assert(LEGAL(in_group, store->n_in_group));
	assert(LEGAL(within_group, store->n_within_group[in_group]));

	GMRFLib_graph_tp *graph = (GMRFLib_graph_tp *) (store->graphs->ptr[in_group]);
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

	sTiles_assign_values(in_group, within_group, &(store->obj), x);
	Free(x);

	return GMRFLib_SUCCESS;
}

int GMRFLib_stiles_solve_LLT(GMRFLib_stiles_idx_tp *stiles_idx, double *rhs)
{
	GMRFLib_stiles_idx_tp lidx;
	Memcpy(&lidx, stiles_idx, sizeof(GMRFLib_stiles_idx_tp));
	lidx.within_group = omp_get_thread_num();
	GMRFLib_stiles_bind(&lidx);
	sTiles_solve_LLT(stiles_idx->in_group, stiles_idx->within_group, &(store->obj), rhs, stiles_idx->nrhs);
	GMRFLib_stiles_unbind(&lidx);

	return GMRFLib_SUCCESS;
}

int GMRFLib_stiles_solve_L(GMRFLib_stiles_idx_tp *stiles_idx, double *rhs)
{
	GMRFLib_stiles_idx_tp lidx;
	Memcpy(&lidx, stiles_idx, sizeof(GMRFLib_stiles_idx_tp));
	lidx.within_group = omp_get_thread_num();
	GMRFLib_stiles_bind(&lidx);
	sTiles_solve_L(stiles_idx->in_group, stiles_idx->within_group, &(store->obj), rhs, stiles_idx->nrhs);
	GMRFLib_stiles_unbind(&lidx);

	return GMRFLib_SUCCESS;
}

int GMRFLib_stiles_solve_LT(GMRFLib_stiles_idx_tp *stiles_idx, double *rhs)
{
	GMRFLib_stiles_idx_tp lidx;
	Memcpy(&lidx, stiles_idx, sizeof(GMRFLib_stiles_idx_tp));
	lidx.within_group = omp_get_thread_num();
	GMRFLib_stiles_bind(&lidx);
	sTiles_solve_LT(stiles_idx->in_group, stiles_idx->within_group, &(store->obj), rhs, stiles_idx->nrhs);
	GMRFLib_stiles_unbind(&lidx);

	return GMRFLib_SUCCESS;
}

int GMRFLib_stiles_Qinv_INLA(GMRFLib_problem_tp *problem)
{
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
		map_id_set(Qinv[i], i, GMRFLib_stiles_Qinv_get(i, i, problem->stiles_idx));
		for (int jj = 0; jj < nnb; jj++) {
			int j = graph->lnbs[i][jj];
			map_id_set(Qinv[i], j, GMRFLib_stiles_Qinv_get(i, j, problem->stiles_idx));
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

	return GMRFLib_SUCCESS;
}

void GMRFLib_stiles_bind(GMRFLib_stiles_idx_tp *stiles_idx)
{
#if 0
	if (!tref) {
#pragma omp critical (Name_32d219aff2336da13ace9454f552486ecea90e5c)
		if (!tref) {
			tref = Calloc(GMRFLib_MAX_THREADS(), double);
		}
	}
	double tt = -GMRFLib_timer();
	assert(store);
#endif

	bool *p = &(store->bind_done[stiles_idx->in_group][stiles_idx->within_group]);
	if (*p == false) {
		sTiles_bind(stiles_idx->in_group, stiles_idx->within_group, &(store->obj));
		*p = true;
	}
#if 0
	int tnum = omp_get_thread_num();
	tt += GMRFLib_timer();
	tref[tnum] += tt;
	if (tnum == 0) {
#pragma omp critical
		{
			printf("BIND ");
			for (int i = 0; i < GMRFLib_MAX_THREADS(); i++) {
				printf(" %.3f", tref[i]);
			}
			printf("\n");
		}
	}
#endif
}

void GMRFLib_stiles_unbind(GMRFLib_stiles_idx_tp *stiles_idx)
{
	bool *p = &(store->bind_done[stiles_idx->in_group][stiles_idx->within_group]);
	if (*p == true) {
		sTiles_unbind(stiles_idx->in_group, stiles_idx->within_group, &(store->obj));
		*p = false;
	}
}

void GMRFLib_stiles_unbind_group(int in_group)
{
	// need to do this one in parallel
#pragma omp parallel for num_threads(store->n_within_group[in_group])
	for (int j = 0; j < store->n_within_group[in_group]; j++) {
		GMRFLib_stiles_idx_tp stiles_idx = { in_group, j, 0 };
		GMRFLib_stiles_unbind(&stiles_idx);
	}
}

void GMRFLib_stiles_unbind_all(void)
{
	for (int g = 0; g < store->n_in_group; g++) {
		GMRFLib_stiles_unbind_group(g);
	}
}

int GMRFLib_stiles_get_verbose()
{
	return (ctl ? ctl->verbose : 0);
}

int GMRFLib_stiles_get_tile_size(void)
{
#if defined(INLA_WITH_STILES)
	int get_auto_tile_size(void);
	if (!ctl) {
		return get_auto_tile_size();
	} else {
		return (ctl->tile_size > 0 ? ctl->tile_size : get_auto_tile_size());
	}
#else
	return (ctl ? ctl->tile_size : 0);
#endif
}
