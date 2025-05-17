#include <stddef.h>
#include <assert.h>
#include <strings.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

static GMRFLib_stiles_ctl_tp *ctl = NULL;
static GMRFLib_stiles_store_tp *store = NULL;
static GMRFLib_idxptr_tp *free_ptrs = NULL;

int GMRFLib_stiles_setup(GMRFLib_stiles_setup_tp *setup)
{
	GMRFLib_idxptr_tp *graphs = setup->graphs;
	GMRFLib_idx_tp *nrhss = setup->nrhss;

	GMRFLib_STOP_IF_NOT_SERIAL();

	if (!ctl) {
		GMRFLib_stiles_set_ctl(0, 0);
	}

	int nt_outer = GMRFLib_openmp->max_threads_nested[0];
	int nt_inner = GMRFLib_openmp->max_threads_nested[1];
	int ng = graphs->n;
	int ng2 = ng + nt_outer;			       /* the copies */

	assert(nt_outer > 0);
	assert(nt_inner > 0);
	assert(ng > 0);

	if (store) {
		GMRFLib_stiles_quit();
	}
	store = Calloc(1, GMRFLib_stiles_store_tp);
	store->n = Calloc(ng2, int);
	store->nnz = Calloc(ng2, int);

	// copy graphs, not the copy-one
	store->graphs = NULL;
	for (int i = 0; i < ng; i++) {
		GMRFLib_graph_tp *g = NULL;
		GMRFLib_graph_duplicate(&g, (GMRFLib_graph_tp *) (graphs->ptr[i]));
		GMRFLib_idxptr_add(&(store->graphs), g);
		store->n[i] = g->n;
		store->nnz[i] = g->nnz;
	}

	int *calls_g = Malloc(ng2, int);
	int *cores_g = Malloc(ng2, int);
	int *zeros = Calloc(ng2, int);
	GMRFLib_ifill(ng2, nt_outer, calls_g);
	GMRFLib_ifill(ng2, nt_inner * nt_outer, cores_g);      /* yes, the total number of threads for one group of matrices */

	bool *inv = Malloc(ng2, bool);
	GMRFLib_bfill(ng2, true, inv);

	int nn = 1 + (nrhss ? nrhss->n : 0);
	assert(nn > 0);
	int *nrhs = Malloc(nn + 1, int);
	nrhs[0] = nn;					       /* the first element is the number of different size of rhs's */
	nrhs[1] = 1;					       /* always */
	if (nrhss && nrhss->n) {
		Memcpy(nrhs + 2, nrhss->idx, nrhss->n * sizeof(int));
		GMRFLib_sort_i(nrhs + 1, nn);

		int nnew = 0, *nrhsnew = NULL;
		GMRFLib_iuniques(&nnew, &nrhsnew, nrhs + 1, nn);
		GMRFLib_sort_i(nrhsnew, nnew);
		nrhs[0] = nnew;
		Memcpy(nrhs + 1, nrhsnew, nnew * sizeof(int));
	}
	store->nrhss = nrhs[0];
	store->rhss = nrhs + 1;

	sTiles_create(&(store->obj), ng2, calls_g, cores_g, zeros, inv, nrhs);
	store->n_in_group = ng2;
	store->offset_copy = ng;
	store->n_within_group = calls_g;
	store->n_cores_group = cores_g;
	store->nt_outer = nt_outer;
	store->nt_inner = nt_inner;
	store->Qinv_done = Malloc(ng2, bool *);
	store->bind_done = Malloc(ng2, bool *);
	for (int i = 0; i < ng2; i++) {
		store->Qinv_done[i] = Malloc(store->n_within_group[i], bool);
		store->bind_done[i] = Malloc(store->n_within_group[i], bool);
		GMRFLib_bfill(store->n_within_group[i], false, store->Qinv_done[i]);
		GMRFLib_bfill(store->n_within_group[i], false, store->bind_done[i]);
	}

	// skip the copy ones
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

		// ptr's are stored, contents not copied, so we need to free them later.
		// since sTiles wants the lower triangular matrix, we just swap 'idx_i' and 'idx_j'

		sTiles_assign_graph(ig, store->obj, g->n, nz, idx_j, idx_i);	/* oops, yes we swap */
		GMRFLib_idxptr_add(&free_ptrs, idx_i);
		GMRFLib_idxptr_add(&free_ptrs, idx_j);

		if (0) {
			printf("\nUpper triangular format\n");
			for (int i = 0; i < nz; i++) {
				printf("idx[%1d] = (%1d, %1d)\n", i, idx_i[i], idx_j[i]);
			}
		}
	}

	// we initialize all but the copy-group's
	for (int i = 0; i < ng; i++) {
		sTiles_init_group(i, &(store->obj));
	}

	// this is the copy ones
	for (int k = 0; k < calls_g[0]; k++) {
		int kk = store->offset_copy + k;
		for (int i = 0; i < calls_g[kk]; i++) {
			sTiles_map_group_call_to_group_call(&(store->obj), kk, i, 0, k);
		}
	}

	// no need for the copy-one
	store->perm = Calloc(ng2, int *);
	store->iperm = Calloc(ng2, int *);
	for (int i = 0; i < ng; i++) {
		int *p;

		p = sTiles_return_perm_vec(i, &(store->obj));
		store->perm[i] = Malloc(store->n[i], int);
		Memcpy(store->perm[i], p, store->n[i] * sizeof(int));

		p = sTiles_return_iperm_vec(i, &(store->obj));
		store->iperm[i] = Malloc(store->n[i], int);
		Memcpy(store->iperm[i], p, store->n[i] * sizeof(int));
	}

	Free(zeros);
	Free(inv);

	if (ctl->verbose) {
		GMRFLib_stiles_print(stdout);
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_stiles_get_offset_copy(void)
{
	return (store->offset_copy);
}

void GMRFLib_stiles_quit(void)
{
	if (!store)
		return;

	GMRFLib_STOP_IF_NOT_SERIAL();

	// int nt_outer = store->nt_outer;
	// int nt_inner = store->nt_inner;

	for (int g = 0; g < store->n_in_group; g++) {
		for (int k = 0; k < store->n_within_group[g]; k++) {
			GMRFLib_stiles_idx_tp stiles_idx = { g, k, -1 };
			GMRFLib_stiles_unbind(&stiles_idx);
		}
	}

	sTiles_quit();
	if (free_ptrs) {
		for (int i = 0; i < free_ptrs->n; i++) {
			Free(free_ptrs->ptr[i]);
		}
		GMRFLib_idxptr_free(free_ptrs);
		free_ptrs = NULL;
	}
	GMRFLib_idxptr_free(store->graphs);

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
	fprintf(fp, "%s:%1d stiles_idx: in_group=%1d within_group=%1d sidx=%1d\n",
		__FILE__, __LINE__, stiles_idx->in_group, stiles_idx->within_group, stiles_idx->sidx);
}

int GMRFLib_stiles_set_idx(GMRFLib_stiles_idx_tp *stiles_idx, int nrhs)
{
	// rewrite ->within_group using thread_num(), and set sidx
	stiles_idx->within_group = (omp_get_thread_num() % store->n_within_group[stiles_idx->in_group]);
	if (nrhs >= 0) {
		if ((stiles_idx->sidx = GMRFLib_find_ivalue(store->rhss, store->nrhss, 1, nrhs)) < 0) {
			fprintf(stderr, "\n%s: %1d: nrhs = %1d not found. continue with nrhs=1\n\n", __FILE__, __LINE__, nrhs);
			stiles_idx->sidx = GMRFLib_find_ivalue(store->rhss, store->nrhss, 1, nrhs);
		}
	} else {
		stiles_idx->sidx = -1;
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_stiles_set_idx_copy(GMRFLib_stiles_idx_tp *stiles_idx, int nrhs)
{
	// rewrite ->in_group and ->within_group into the corresponding copy, and set sidx
	stiles_idx->in_group = store->offset_copy + stiles_idx->within_group;
	stiles_idx->within_group = (omp_get_thread_num() % store->n_within_group[stiles_idx->in_group]);
	if (nrhs >= 0) {
		if ((stiles_idx->sidx = GMRFLib_find_ivalue(store->rhss, store->nrhss, 1, nrhs)) < 0) {
			fprintf(stderr, "\n%s: %1d: nrhs = %1d not found\n\n", __FILE__, __LINE__, nrhs);
			stiles_idx->sidx = GMRFLib_find_ivalue(store->rhss, store->nrhss, 1, nrhs);
		}
	} else {
		stiles_idx->sidx = -1;
	}
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
		fprintf(fp, "\n\ncontent of 'store':\n");
		fprintf(fp, "\t\tngroup[%1d] verbose[%1d] tile.size[%1d]\n", store->n_in_group, ctl->verbose, GMRFLib_stiles_get_tile_size());
		fprintf(fp, "\t\tnt_outer[%1d] nt_inner[%1d] offset_copy[%1d]\n", store->nt_outer, store->nt_inner, store->offset_copy);

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
						fprintf(fp, "\t\tperm[%1d][%5d] = %5d  iperm[%1d][%5d] = %5d\n", i, j,
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
			} else {
				fprintf(fp, "\t\tthis group is 'copy' of in_group[0] within_group[%1d]\n", i - store->offset_copy);
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

int GMRFLib_stiles_get_tile_size(void)
{
	return sTiles_return_tile_size();
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
	if (ctl->tile_size) {
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
		GMRFLib_idxptr_free(setup->graphs);
		GMRFLib_idx_free(setup->nrhss);
		Free(setup);
	}
}

int GMRFLib_stiles_chol(GMRFLib_stiles_idx_tp *stiles_idx)
{
	int in_group = stiles_idx->in_group;
	int within_group = stiles_idx->within_group;

	int status = sTiles_chol(in_group, within_group, &(store->obj));
	if (status) {
		fprintf(stderr, "\n\n*** ERROR *** sTiles_chol %d \n\n", status);
		fflush(stderr);
	}
	return (status ? !GMRFLib_SUCCESS : GMRFLib_SUCCESS);
}

double GMRFLib_stiles_logdet(GMRFLib_stiles_idx_tp *stiles_idx)
{
	int in_group = stiles_idx->in_group;
	int within_group = stiles_idx->within_group;

	return sTiles_get_logdet(in_group, within_group, &(store->obj));
}

void GMRFLib_stiles_Qinv(GMRFLib_stiles_idx_tp *stiles_idx)
{
	int in_group = stiles_idx->in_group;
	int within_group = stiles_idx->within_group;

	if (store->Qinv_done[in_group][within_group]) {
		sTiles_clear_selinv(in_group, within_group, &(store->obj));
	} else {
		store->Qinv_done[in_group][within_group] = true;
	}
	sTiles_selinv(in_group, within_group, &(store->obj));
}

double GMRFLib_stiles_Qinv_get(int i, int j, GMRFLib_stiles_idx_tp *stiles_idx)
{
	int in_group = stiles_idx->in_group;
	int within_group = stiles_idx->within_group;

	return sTiles_get_selinv_elm(in_group, within_group, i, j, &(store->obj));
}

int GMRFLib_stiles_build(GMRFLib_stiles_idx_tp *stiles_idx, int thread_id, GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_arg)
{
	int in_group = stiles_idx->in_group;
	int within_group = stiles_idx->within_group;

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
		FIXME1("FAST COPY");
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
	assert(stiles_idx->sidx >= 0);
	sTiles_solve_LLT(stiles_idx->in_group, stiles_idx->within_group, &(store->obj), rhs, stiles_idx->sidx);
	return GMRFLib_SUCCESS;
}

int GMRFLib_stiles_solve_L(GMRFLib_stiles_idx_tp *stiles_idx, double *rhs)
{
	assert(stiles_idx->sidx >= 0);
	sTiles_solve_L(stiles_idx->in_group, stiles_idx->within_group, &(store->obj), rhs, stiles_idx->sidx);
	return GMRFLib_SUCCESS;
}

int GMRFLib_stiles_solve_LT(GMRFLib_stiles_idx_tp *stiles_idx, double *rhs)
{
	assert(stiles_idx->sidx >= 0);
	sTiles_solve_LT(stiles_idx->in_group, stiles_idx->within_group, &(store->obj), rhs, stiles_idx->sidx);
	return GMRFLib_SUCCESS;
}

int GMRFLib_stiles_Qinv_INLA(GMRFLib_problem_tp *problem)
{
	if (problem == NULL) {
		return GMRFLib_SUCCESS;
	}

	int in_group = problem->stiles_idx->in_group;
	// int within_group = problem->stiles_idx->within_group;
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
	assert(store);
	bool *p = &(store->bind_done[stiles_idx->in_group][stiles_idx->within_group]);
	if (*p == false) {
		sTiles_bind(stiles_idx->in_group, stiles_idx->within_group, &(store->obj));
		*p = true;
	}
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
#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_nested[0])
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

int GMRFLib_stiles_verbose()
{
	return (ctl ? ctl->verbose : 0);
}

int GMRFLib_stiles_set_tile_size(void)
{
	return (ctl ? ctl->tile_size : 0);
}


//
//
// TEST FUNCTIONS GOES HERE
//
//

double GMRFLib_stiles_test_Qfunc(int UNUSED(thread_id), int i, int j, double *UNUSED(values), void *args)
{
	if (j < 0)
		return NAN;

	GMRFLib_graph_tp *g = (GMRFLib_graph_tp *) args;
	return (double) ((i == j) ? ISQR(g->n) + 2 * i : -(1 + IMIN(i, j)));
}

double GMRFLib_stiles_test_Qfunc2(int thread_id, int i, int j, double *values, void *args)
{
	GMRFLib_graph_tp *g = (GMRFLib_graph_tp *) args;
	if (j < 0) {
		int k = 0;
		values[k++] = GMRFLib_stiles_test_Qfunc(thread_id, i, i, NULL, args);
		for (int jj = 0; jj < g->lnnbs[i]; jj++) {
			j = g->lnbs[i][jj];
			values[k++] = GMRFLib_stiles_test_Qfunc(thread_id, i, j, NULL, args);
		}
	} else {
		return GMRFLib_stiles_test_Qfunc(thread_id, i, j, NULL, args);
	}
	return 0.0;
}

double GMRFLib_stiles_test_Qfunc3(int thread_id, int i, int j, double *values, void *args)
{
	GMRFLib_graph_tp *g = (GMRFLib_graph_tp *) args;
	int n = g->n;
	double n2 = n * 2.0;

	if (j < 0) {
		int k = 0;
		values[k++] = (2.0 * n + 1.0 + thread_id) / n2;
		for (int jj = 0; jj < g->lnnbs[i]; jj++) {
			j = g->lnbs[i][jj];
			values[k++] = -1.0 / n2;
		}
	} else {
		if (i == j) {
			return (2.0 * n + 1.0 + thread_id) / n2;
		} else {
			return -1.0 / n2;
		}
	}
	return 0.0;
}

int GMRFLib_stiles_test(void)
{
	GMRFLib_STOP_IF_NOT_SERIAL();

	GMRFLib_stiles_idx_tp *stiles_idx = Calloc(1, GMRFLib_stiles_idx_tp);

	GMRFLib_stiles_set_ctl(0, 0);

	int n[4] = { 7, 7, 7, 7 };
	GMRFLib_graph_tp *g[4] = { NULL, NULL, NULL, NULL };

	GMRFLib_graph_mk_linear(&(g[0]), n[0], 1, 0);
	GMRFLib_graph_mk_linear(&(g[1]), n[1], 1, 0);
	GMRFLib_graph_mk_linear(&(g[2]), n[2], 1, 0);
	GMRFLib_graph_mk_linear(&(g[3]), n[3], 1, 0);

	GMRFLib_idxptr_tp *graphs = NULL;
	GMRFLib_idxptr_add(&graphs, g[0]);
	GMRFLib_idxptr_add(&graphs, g[1]);
	GMRFLib_idxptr_add(&graphs, g[2]);
	GMRFLib_idxptr_add(&graphs, g[3]);

	GMRFLib_stiles_setup_tp setup = { graphs, NULL };
	GMRFLib_stiles_setup(&setup);

	if (0) {
		GMRFLib_printf_graph(stdout, g[0]);
		GMRFLib_printf_graph(stdout, g[1]);
		GMRFLib_printf_graph(stdout, g[2]);
		GMRFLib_printf_graph(stdout, g[3]);
	}

	GMRFLib_tabulate_Qfunc_tp *tab[2] = { NULL, NULL };
	GMRFLib_tabulate_Qfunc(0, &(tab[0]), g[0], GMRFLib_stiles_test_Qfunc, (void *) g[0], NULL);
	GMRFLib_tabulate_Qfunc(0, &(tab[1]), g[1], GMRFLib_stiles_test_Qfunc2, (void *) g[1], NULL);

	stiles_idx->in_group = 0;
	GMRFLib_stiles_build(stiles_idx, 0, GMRFLib_stiles_test_Qfunc, (void *) g[0]);
	stiles_idx->in_group = 1;
	GMRFLib_stiles_build(stiles_idx, 0, GMRFLib_stiles_test_Qfunc2, (void *) g[1]);
	stiles_idx->in_group = 2;
	GMRFLib_stiles_build(stiles_idx, 0, tab[0]->Qfunc, tab[0]->Qfunc_arg);
	stiles_idx->in_group = 3;
	GMRFLib_stiles_build(stiles_idx, 0, tab[1]->Qfunc, tab[1]->Qfunc_arg);

	double *rhs = Malloc(n[0], double);

	stiles_idx->in_group = 0;
	sTiles_bind(stiles_idx->in_group, stiles_idx->within_group, &(store->obj));
	GMRFLib_stiles_chol(stiles_idx);

	for (int i = 0; i < n[0]; i++) {
		rhs[i] = i;
	}

	stiles_idx->in_group = 0;
	GMRFLib_stiles_solve_LLT(stiles_idx, rhs);
	for (int i = 0; i < n[0]; i++) {
		printf("solve_LLT[%1d] = %f\n", i, rhs[i]);
	}

	for (int i = 0; i < n[0]; i++) {
		rhs[i] = i;
		printf("rhs[%1d] = %f\n", i, rhs[i]);
	}
	GMRFLib_stiles_solve_L(stiles_idx, rhs);
	for (int i = 0; i < n[0]; i++) {
		printf("solve_L[%1d] = %f\n", i, rhs[i]);
	}

	for (int i = 0; i < n[0]; i++) {
		rhs[i] = i;
		printf("rhs[%1d] = %f\n", i, rhs[i]);
	}
	GMRFLib_stiles_solve_LT(stiles_idx, rhs);
	for (int i = 0; i < n[0]; i++) {
		printf("solve_LT[%1d] = %f\n", i, rhs[i]);
	}

	stiles_idx->in_group = 1;
	sTiles_bind(stiles_idx->in_group, stiles_idx->within_group, &(store->obj));
	GMRFLib_stiles_chol(stiles_idx);

	stiles_idx->in_group = 2;
	sTiles_bind(stiles_idx->in_group, stiles_idx->within_group, &(store->obj));
	GMRFLib_stiles_chol(stiles_idx);

	stiles_idx->in_group = 3;
	sTiles_bind(stiles_idx->in_group, stiles_idx->within_group, &(store->obj));
	GMRFLib_stiles_chol(stiles_idx);

	for (int i = 0; i < n[0]; i++) {
		printf("perm[%1d] = %1d iperm = %1d\n", i, store->perm[0][i], store->iperm[0][i]);
	}

	stiles_idx->in_group = 0;
	printf("Q[0]->logdet %f\n", GMRFLib_stiles_logdet(stiles_idx));
	stiles_idx->in_group = 1;
	printf("Q[1]->logdet %f\n", GMRFLib_stiles_logdet(stiles_idx));
	stiles_idx->in_group = 2;
	printf("Q[2]->logdet %f\n", GMRFLib_stiles_logdet(stiles_idx));
	stiles_idx->in_group = 3;
	printf("Q[3]->logdet %f\n", GMRFLib_stiles_logdet(stiles_idx));

	stiles_idx->in_group = 0;
	GMRFLib_stiles_Qinv(stiles_idx);
	for (int i = 0; i < g[0]->n; i++) {
		printf("Qinv[%1d, %1d] = %.8f\n", i, i, GMRFLib_stiles_Qinv_get(i, i, stiles_idx));
		for (int jj = 0; jj < g[0]->nnbs[i]; jj++) {
			int j = g[0]->nbs[i][jj];
			printf("\tQinv[%1d, %1d] = %.8f\n", i, j, GMRFLib_stiles_Qinv_get(i, j, stiles_idx));
		}
	}

	FILE *fp = fopen("Q1.txt", "w");
	GMRFLib_printf_Qfunc(0, fp, g[0], GMRFLib_stiles_test_Qfunc, (void *) g[0]);
	fclose(fp);

	GMRFLib_stiles_print(stdout);

	GMRFLib_stiles_quit();
	return GMRFLib_SUCCESS;
}

int GMRFLib_stiles_test2(void)
{
	GMRFLib_STOP_IF_NOT_SERIAL();

	GMRFLib_graph_tp *graph = NULL;
	GMRFLib_graph_read(&graph, "germany.graph");
	int n = graph->n;

	GMRFLib_constr_tp *constr = NULL;
	GMRFLib_make_empty_constr(&constr);
	int nc = 5;
	constr->a_matrix = Calloc(nc * n, double);
	constr->e_vector = Calloc(nc, double);
	for (int i = 0; i < n * nc; i++) {
		constr->a_matrix[i] = GMRFLib_uniform();
	}
	for (int i = 0; i < nc; i++) {
		constr->e_vector[i] = GMRFLib_uniform();
	}
	GMRFLib_prepare_constr(constr, graph, GMRFLib_TRUE);

	int thread_id = 0;

	GMRFLib_idxptr_tp *graphs = NULL;
	GMRFLib_idxptr_add(&graphs, graph);
	GMRFLib_stiles_setup_tp setup = { graphs, NULL };
	GMRFLib_stiles_setup(&setup);

	GMRFLib_stiles_idx_tp *stiles_idx;
	stiles_idx = Calloc(1, GMRFLib_stiles_idx_tp);
	stiles_idx->in_group = 0;
	stiles_idx->within_group = 0;

	double *b = Malloc(n, double);
	for (int i = 0; i < n; i++) {
		b[i] = n * (i - n / 2.0);		       // / (double) n;
	}

	GMRFLib_smtp_tp smtp = GMRFLib_SMTP_TAUCS;
	GMRFLib_problem_tp *problem = NULL;
	GMRFLib_init_problem(thread_id, &problem, NULL, b, NULL, NULL, graph, GMRFLib_stiles_test_Qfunc, (void *) graph, constr, stiles_idx, &smtp);
	GMRFLib_compute_Qinv_TAUCS(problem);
	double *v = Malloc(n, double);
	for (int i = 0; i < n; i++) {
		v[i] = *GMRFLib_Qinv_get(problem, i, i);
	}

	smtp = GMRFLib_SMTP_STILES;
	GMRFLib_problem_tp *problem2 = NULL;

	GMRFLib_stiles_bind(stiles_idx);
	GMRFLib_init_problem(thread_id, &problem2, NULL, b, NULL, NULL, graph, GMRFLib_stiles_test_Qfunc, (void *) graph, constr,
			     stiles_idx, &smtp);

	GMRFLib_compute_Qinv(problem2);
	P(sTiles_get_selinv_elm(stiles_idx->in_group, stiles_idx->within_group, 0, 0, &(store->obj)));
	P(*GMRFLib_Qinv_get(problem2, 0, 0));

	GMRFLib_compute_Qinv(problem2);
	P(sTiles_get_selinv_elm(stiles_idx->in_group, stiles_idx->within_group, 0, 0, &(store->obj)));
	P(*GMRFLib_Qinv_get(problem2, 0, 0));


	double *v2 = Malloc(n, double);
	for (int i = 0; i < n; i++) {
		v2[i] = *GMRFLib_Qinv_get(problem2, i, i);
	}

	double err_diff_mean = 0.0;
	double err_diff_var = 0.0;
	for (int i = 0; i < n; i++) {
		err_diff_mean += SQR(problem->mean_constr[i] - problem2->mean_constr[i]);
		err_diff_var += SQR(1.0 - v2[i] / v[i]);
	}

	P(sqrt(err_diff_mean / n));
	P(sqrt(err_diff_var / n));

	for (int k = 0; k < 100; k++) {
		b = Malloc(n, double);

		GMRFLib_stiles_bind(stiles_idx);
		GMRFLib_init_problem(thread_id, &problem2, NULL, b, NULL, NULL, graph, GMRFLib_stiles_test_Qfunc, (void *) graph,
				     constr, stiles_idx, &smtp);
		GMRFLib_compute_Qinv(problem2);
		P(*GMRFLib_Qinv_get(problem2, 0, 0));

		for (int i = 0; i < n; i++)
			b[i] = GMRFLib_uniform();
		GMRFLib_solve_llt_sparse_matrix(b, 1, &(problem2->sub_sm_fact), graph, problem2, NULL);
		GMRFLib_solve_lt_sparse_matrix(b, 1, &(problem2->sub_sm_fact), graph, problem2);
		GMRFLib_solve_l_sparse_matrix(b, 1, &(problem2->sub_sm_fact), graph, problem2);
		for (int i = 0; i < n; i++)
			assert(!ISNAN(b[i]));
		Free(b);

		GMRFLib_free_problem(problem2);
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_stiles_test3(void)
{
	int nh = 0;
	assert(GMRFLib_smtp == GMRFLib_SMTP_STILES);
	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_OPTIMIZE, (void *) &nh, &GMRFLib_smtp);
	// omp_set_max_active_levels(2);

	int nt = GMRFLib_openmp->max_threads_nested[0];
	int n = 5000;
	int mm = 40;
	int m = 1 * mm * nt;

	P(n);
	P(mm);
	P(m);

	GMRFLib_graph_tp *graph = NULL;
	GMRFLib_graph_mk_linear(&graph, n, n - 1, 1);

	GMRFLib_idxptr_tp *graphs = NULL;
	GMRFLib_idx_tp *rh = NULL;

	GMRFLib_idxptr_add(&graphs, graph);
	GMRFLib_idx_add(&rh, mm);
	GMRFLib_idx_add(&rh, m);
	GMRFLib_stiles_setup_tp setup = { graphs, rh };
	GMRFLib_stiles_setup(&setup);
	GMRFLib_stiles_print(stdout);

	double *b = Malloc(n, double);
	for (int i = 0; i < n; i++) {
		b[i] = (i - n / 2.0) / (double) n;
	}

	P(GMRFLib_openmp->max_threads_nested[0]);
	P(GMRFLib_openmp->max_threads_nested[1]);

#pragma omp parallel for num_threads(nt)
	for (int k = 0; k < nt; k++) {
		GMRFLib_stiles_idx_tp stiles_idx = { 0, 0, 0 };
		GMRFLib_stiles_set_idx(&stiles_idx, 1);
		GMRFLib_stiles_bind(&stiles_idx);

		int thread_id = k;

		printf("running with k %1d tnum %d within_group %1d\n", k, omp_get_thread_num(), stiles_idx.within_group);
		fflush(stdout);

		GMRFLib_problem_tp *problem = NULL;
		GMRFLib_init_problem(thread_id, &problem, NULL, b, NULL, NULL, graph, GMRFLib_stiles_test_Qfunc3, (void *) graph,
				     NULL, &stiles_idx, NULL);
		// GMRFLib_compute_Qinv(problem);
		// printf("k %1d Qinv[,] = %.12f %.12f\n", k, *GMRFLib_Qinv_get(problem, 0, 0), *GMRFLib_Qinv_get(problem, n - 1, n 
		// 
		// 
		// - 1));
		GMRFLib_free_problem(problem);
	}
	GMRFLib_stiles_unbind_group(0);

	double *B = Malloc(n * m, double);
	double *BB = Malloc(n * m, double);
	for (int i = 0; i < m * n; i++) {
		B[i] = BB[i] = GMRFLib_uniform() - 0.5;
	}

	FIXME("solve m start");
	double tref = -GMRFLib_timer();
#pragma omp parallel for num_threads(nt)
	for (int i = 0; i < 1; i++) {
		GMRFLib_stiles_idx_tp stiles_idx = { 0, -1, -1 };
		GMRFLib_stiles_set_idx(&stiles_idx, m);
		GMRFLib_stiles_bind(&stiles_idx);
		GMRFLib_stiles_solve_LLT(&stiles_idx, B);
	}
	GMRFLib_stiles_unbind_group(0);
	tref += GMRFLib_timer();
	FIXME("solve m end");
	P(tref);

	assert((m / (mm * nt)) * (mm * nt) == m);

	FIXME("solve mm start");
	tref = -GMRFLib_timer();
	for (int j = 0; j < m / (mm * nt); j++) {
#pragma omp parallel for num_threads(nt)
		for (int i = 0; i < nt; i++) {
			GMRFLib_stiles_idx_tp stiles_idx = { 0, i, -1 };
			GMRFLib_stiles_set_idx_copy(&stiles_idx, mm);
			GMRFLib_stiles_bind(&stiles_idx);
			GMRFLib_stiles_solve_LLT(&stiles_idx, BB + (j * mm * nt + i * mm) * n);
		}
	}
	GMRFLib_stiles_unbind_group(store->offset_copy);
	FIXME("solve mm end");

	tref += GMRFLib_timer();
	P(tref);

	double err = 0.0;
	for (int i = 0; i < n * m; i++) {
		err += SQR(B[i] - BB[i]);
		if (0 && err > FLT_EPSILON) {
			printf("i %d err %.8f\n", i, sqrt(err / (i + 1)));
			abort();
		}
	}
	P(sqrt(err / n / m));
	GMRFLib_stiles_unbind_all();
	GMRFLib_stiles_quit();
	return GMRFLib_SUCCESS;
}
