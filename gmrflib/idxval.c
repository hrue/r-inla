#include <assert.h>
#include <float.h>
#include <math.h>
#include <omp.h>
#include <signal.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <time.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/hashP.h"

#if defined(INLA_WITH_ARMPL)
#       include "armpl_sparse.h"
#endif

#define IDX_ALLOC_INITIAL 8
#define IDX_ALLOC_ADD     64
#define IDX_ALLOC_NDIV    4
#define IDX_ALLOC_INCREASE (IMAX(IDX_ALLOC_ADD, (*hold)->n / IDX_ALLOC_NDIV))

int GMRFLib_idx_create(GMRFLib_idx_tp **hold)
{
	return GMRFLib_idx_create_x(hold, IDX_ALLOC_INITIAL);
}

int GMRFLib_idx_create_x(GMRFLib_idx_tp **hold, int len)
{
	len = IMAX(1, len);
	*hold = Calloc(1, GMRFLib_idx_tp);
	(*hold)->idx = Calloc(len, int);
	(*hold)->n_alloc = len;
	(*hold)->n = 0;

	return GMRFLib_SUCCESS;
}

int GMRFLib_str_create(GMRFLib_str_tp **hold)
{
	return GMRFLib_str_create_x(hold, IDX_ALLOC_INITIAL);
}

int GMRFLib_str_create_x(GMRFLib_str_tp **hold, int len)
{
	len = IMAX(1, len);
	*hold = Calloc(1, GMRFLib_str_tp);
	(*hold)->str = Calloc(len, char *);
	(*hold)->n_alloc = len;
	(*hold)->n = 0;

	return GMRFLib_SUCCESS;
}

int GMRFLib_idx2_create(GMRFLib_idx2_tp **hold)
{
	return GMRFLib_idx2_create_x(hold, IDX_ALLOC_INITIAL);
}

int GMRFLib_idx2_create_x(GMRFLib_idx2_tp **hold, int len)
{
	len = IMAX(1, len);
	*hold = Calloc(1, GMRFLib_idx2_tp);
	(*hold)->idx = Calloc(2, int *);
	(*hold)->idx[0] = Calloc(len, int);
	(*hold)->idx[1] = Calloc(len, int);
	(*hold)->n_alloc = len;
	(*hold)->n = 0;

	return GMRFLib_SUCCESS;
}

int GMRFLib_val_create(GMRFLib_val_tp **hold)
{
	*hold = Calloc(1, GMRFLib_val_tp);
	(*hold)->val = Calloc(IDX_ALLOC_INITIAL, double);
	(*hold)->n_alloc = IDX_ALLOC_INITIAL;
	(*hold)->n = 0;

	return GMRFLib_SUCCESS;
}

int GMRFLib_ptr_create(GMRFLib_ptr_tp **hold)
{
	return GMRFLib_ptr_create_x(hold, IDX_ALLOC_INITIAL);
}

int GMRFLib_ptr_create_x(GMRFLib_ptr_tp **hold, int len)
{
	*hold = Calloc(1, GMRFLib_ptr_tp);
	(*hold)->ptr = Malloc(len, void *);
	(*hold)->n_alloc = len;
	(*hold)->n = 0;

	return GMRFLib_SUCCESS;
}

int GMRFLib_idxval_create(GMRFLib_idxval_tp **hold)
{
	return GMRFLib_idxval_create_x(hold, IDX_ALLOC_INITIAL);
}

int GMRFLib_idxval_create_x(GMRFLib_idxval_tp **hold, int len)
{
	len = IMAX(1, len);
	*hold = Calloc(1, GMRFLib_idxval_tp);
	(*hold)->idx = Calloc(len, int);
	(*hold)->val = Calloc(len, double);
	(*hold)->n_alloc = len;
	(*hold)->n = 0;
	(*hold)->iaddto = 0;

	return GMRFLib_SUCCESS;
}

GMRFLib_idx_tp **GMRFLib_idx_ncreate(int n)
{
	if (n > 0) {
		GMRFLib_idx_tp **a = Calloc(n, GMRFLib_idx_tp *);
		for (int i = 0; i < n; i++) {
			GMRFLib_idx_create(&(a[i]));
		}
		return a;
	} else {
		return NULL;
	}
}

GMRFLib_idx_tp **GMRFLib_idx_ncreate_x(int n, int len)
{
	if (n > 0) {
		GMRFLib_idx_tp **a = Calloc(n, GMRFLib_idx_tp *);
		for (int i = 0; i < n; i++) {
			GMRFLib_idx_create_x(&(a[i]), len);
		}
		return a;
	} else {
		return NULL;
	}
}

GMRFLib_str_tp **GMRFLib_str_ncreate(int n)
{
	if (n > 0) {
		GMRFLib_str_tp **a = Calloc(n, GMRFLib_str_tp *);
		for (int i = 0; i < n; i++) {
			GMRFLib_str_create(&(a[i]));
		}
		return a;
	} else {
		return NULL;
	}
}

GMRFLib_str_tp **GMRFLib_str_ncreate_x(int n, int len)
{
	if (n > 0) {
		GMRFLib_str_tp **a = Calloc(n, GMRFLib_str_tp *);
		for (int i = 0; i < n; i++) {
			GMRFLib_str_create_x(&(a[i]), len);
		}
		return a;
	} else {
		return NULL;
	}
}

GMRFLib_idx2_tp **GMRFLib_idx2_ncreate(int n)
{
	if (n > 0) {
		GMRFLib_idx2_tp **a = Calloc(n, GMRFLib_idx2_tp *);
		for (int i = 0; i < n; i++) {
			GMRFLib_idx2_create(&(a[i]));
		}
		return a;
	} else {
		return NULL;
	}
}

GMRFLib_idx2_tp **GMRFLib_idx2_ncreate_x(int n, int len)
{
	if (n > 0) {
		GMRFLib_idx2_tp **a = Calloc(n, GMRFLib_idx2_tp *);
		for (int i = 0; i < n; i++) {
			GMRFLib_idx2_create_x(&(a[i]), len);
		}
		return a;
	} else {
		return NULL;
	}
}

GMRFLib_val_tp **GMRFLib_val_ncreate(int n)
{
	if (n > 0) {
		GMRFLib_val_tp **a = Calloc(n, GMRFLib_val_tp *);
		for (int i = 0; i < n; i++) {
			GMRFLib_val_create(&(a[i]));
		}
		return a;
	} else {
		return NULL;
	}
}

GMRFLib_ptr_tp **GMRFLib_ptr_ncreate(int n)
{
	if (n > 0) {
		GMRFLib_ptr_tp **a = Calloc(n, GMRFLib_ptr_tp *);
		for (int i = 0; i < n; i++) {
			GMRFLib_ptr_create(&(a[i]));
		}
		return a;
	} else {
		return NULL;
	}
}

GMRFLib_idxval_tp **GMRFLib_idxval_ncreate(int n)
{
	if (n > 0) {
		GMRFLib_idxval_tp **a = Calloc(n, GMRFLib_idxval_tp *);
		for (int i = 0; i < n; i++) {
			GMRFLib_idxval_create(&(a[i]));
		}
		return a;
	} else {
		return NULL;
	}
}

GMRFLib_idxval_tp **GMRFLib_idxval_ncreate_x(int n, int len, int num_threads)
{
	if (n > 0) {
		GMRFLib_idxval_tp **a = Calloc(n, GMRFLib_idxval_tp *);
		if (num_threads > 0) {
#pragma omp parallel for num_threads(num_threads)
			for (int i = 0; i < n; i++) {
				GMRFLib_idxval_create_x(&(a[i]), len);
			}
		} else {
			for (int i = 0; i < n; i++) {
				GMRFLib_idxval_create_x(&(a[i]), len);
			}
		}
		return a;
	} else {
		return NULL;
	}
}

int GMRFLib_idx_printf(FILE *fp, GMRFLib_idx_tp *hold, const char *msg)
{
	if (hold) {
		fprintf(fp, "[%s] n = %1d  nalloc = %1d\n", msg, hold->n, hold->n_alloc);
		for (int i = 0; i < hold->n; i++) {
			fprintf(fp, "\tidx[%1d] = %1d\n", i, hold->idx[i]);
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_str_printf(FILE *fp, GMRFLib_str_tp *hold, const char *msg)
{
	if (hold) {
		fprintf(fp, "[%s] n = %1d  nalloc = %1d\n", msg, hold->n, hold->n_alloc);
		for (int i = 0; i < hold->n; i++) {
			fprintf(fp, "\tidx[%1d] = %s\n", i, hold->str[i]);
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idx2_printf(FILE *fp, GMRFLib_idx2_tp *hold, const char *msg)
{
	if (hold) {
		fprintf(fp, "[%s] n = %1d  nalloc = %1d\n", msg, hold->n, hold->n_alloc);
		for (int i = 0; i < hold->n; i++) {
			fprintf(fp, "\tidx[][%1d] = %1d %1d\n", i, hold->idx[0][i], hold->idx[1][i]);
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_val_printf(FILE *fp, GMRFLib_val_tp *hold, const char *msg)
{
	if (hold) {
		fprintf(fp, "[%s] n = %1d  nalloc = %1d\n", msg, hold->n, hold->n_alloc);
		for (int i = 0; i < hold->n; i++) {
			fprintf(fp, "\tval[%1d] = %g\n", i, hold->val[i]);
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_ptr_printf(FILE *fp, GMRFLib_ptr_tp *hold, const char *msg)
{
	if (hold) {
		fprintf(fp, "[%s] n = %1d  nalloc = %1d\n", msg, hold->n, hold->n_alloc);
		for (int i = 0; i < hold->n; i++) {
			fprintf(fp, "\tptr[%1d] = %p\n", i, hold->ptr[i]);
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idxval_printf(FILE *fp, GMRFLib_idxval_tp *hold, const char *msg)
{
	if (hold) {
		int show_details = 0;
		fprintf(fp, "[%s] n = %1d  nalloc = %1d iaddto = %1d\n", msg, hold->n, hold->n_alloc, hold->iaddto);
		if (show_details) {
			for (int i = 0; i < hold->n; i++) {
				fprintf(fp, "\t(idx, val)[%1d] = (%d, %g)\n", i, hold->idx[i], hold->val[i]);
			}
		}
		if (hold->g_n) {
			fprintf(fp, "\tg_n = %1d\n", hold->g_n);
			for (int g = 0; g < hold->g_n; g++) {
				fprintf(fp, "\tgroup %d has length %d (one=%s) (aligned=%s:%s)\n", g, hold->g_len[g],
					(hold->g_1 && hold->g_1[g] ? "TRUE" : "FALSE"),
					(GMRFLib_is_aligned(hold->g_idx[g][0]) ? "TRUE" : "FALSE"), (GMRFLib_is_aligned(hold->g_val[g][0]) ? "TRUE" : "FALSE"));
				if (show_details) {
					fprintf(fp, "\t\t");
					for (int k = 0; k < IABS(hold->g_len[g]); k++) {
						fprintf(fp, " %1d(%g)", hold->g_idx[g][k], hold->g_val[g][k]);

					}
					fprintf(fp, "\n");
				}
			}
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idxval_info_printf(FILE *fp, GMRFLib_idxval_tp *hold, const char *msg)
{
	if (hold) {
		fprintf(fp, "[%s] n=%1d  nalloc=%1d iaddto=%1d\n", (msg ? msg : ""), hold->n, hold->n_alloc, hold->iaddto);
		if (hold->g_n) {
			for (int g = 0; g < hold->g_n; g++) {
				fprintf(fp, "\tgroup %1d/%1d: len=%1d one=%s\n", g, hold->g_n, hold->g_len[g],
					(hold->g_1 && hold->g_1[g] ? "TRUE" : "FALSE"));
			}
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idx_nprune(GMRFLib_idx_tp **a, int n)
{
	if (a) {
		for (int i = 0; i < n; i++) {
			GMRFLib_idx_prune(a[i]);
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_str_nprune(GMRFLib_str_tp **a, int n)
{
	if (a) {
		for (int i = 0; i < n; i++) {
			GMRFLib_str_prune(a[i]);
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idx_sort(GMRFLib_idx_tp *hold)
{
	if (hold) {
		GMRFLib_sort_i(hold->idx, hold->n);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idx_nsort(GMRFLib_idx_tp **a, int n, int nt)
{
#define CODE_BLOCK							\
	for(int i = 0; i < n; i++) {					\
		CODE_BLOCK_INIT();					\
		if (a[i] && a[i]->n > 1) {				\
			GMRFLib_sort_i(a[i]->idx, a[i]->n);		\
		}							\
	}

	RUN_CODE_BLOCK(nt, 0, 0);
#undef CODE_BLOCK

	return GMRFLib_SUCCESS;
}

int GMRFLib_idx_uniq(GMRFLib_idx_tp *hold)
{
	if (hold && hold->n > 1) {
		int i, j;

		GMRFLib_idx_sort(hold);
		for (j = 0, i = 0; i < hold->n; i++) {
			if (hold->idx[j] != hold->idx[i]) {
				hold->idx[++j] = hold->idx[i];
			}
		}
		hold->n = j + 1;
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idx_nuniq(GMRFLib_idx_tp **a, int n, int nt)
{
#define CODE_BLOCK				\
	for (int i = 0; i < n; i++) {		\
		CODE_BLOCK_INIT();		\
		GMRFLib_idx_uniq(a[i]);		\
	}

	RUN_CODE_BLOCK(nt, 0, 0);
#undef CODE_BLOCK

	return GMRFLib_SUCCESS;
}

int GMRFLib_idx2_nprune(GMRFLib_idx2_tp **a, int n)
{
	if (a) {
		for (int i = 0; i < n; i++) {
			GMRFLib_idx2_prune(a[i]);
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_val_nprune(GMRFLib_val_tp **a, int n)
{
	if (a) {
		for (int i = 0; i < n; i++) {
			GMRFLib_val_prune(a[i]);
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_ptr_nprune(GMRFLib_ptr_tp **a, int n)
{
	if (a) {
		for (int i = 0; i < n; i++) {
			GMRFLib_ptr_prune(a[i]);
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idxval_nprune(GMRFLib_idxval_tp **a, int n, int nt)
{
#define CODE_BLOCK					\
	for (int i = 0; i < n; i++) {			\
		CODE_BLOCK_INIT();			\
		GMRFLib_idxval_prune(a[i]);		\
	}

	RUN_CODE_BLOCK(nt, 0, 0);
#undef CODE_BLOCK

	return GMRFLib_SUCCESS;
}

int GMRFLib_idx_prune(GMRFLib_idx_tp *hold)
{
	if (hold) {
		if (hold->n_alloc - hold->n > IDX_ALLOC_INITIAL) {
			hold->idx = Realloc(hold->idx, IMAX(1, hold->n), int);
			hold->n_alloc = IMAX(1, hold->n);
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_str_prune(GMRFLib_str_tp *hold)
{
	if (hold) {
		if (hold->n_alloc - hold->n > IDX_ALLOC_INITIAL) {
			hold->str = Realloc(hold->str, IMAX(1, hold->n), char *);
			hold->n_alloc = IMAX(1, hold->n);
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idx2_prune(GMRFLib_idx2_tp *hold)
{
	if (hold) {
		if (hold->n_alloc - hold->n > IDX_ALLOC_INITIAL) {
			hold->idx[0] = Realloc(hold->idx[0], IMAX(1, hold->n), int);
			hold->idx[1] = Realloc(hold->idx[1], IMAX(1, hold->n), int);
			hold->n_alloc = IMAX(1, hold->n);
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_val_prune(GMRFLib_val_tp *hold)
{
	if (hold) {
		if (hold->n_alloc - hold->n > IDX_ALLOC_INITIAL) {
			hold->val = Realloc(hold->val, IMAX(1, hold->n), double);
			hold->n_alloc = IMAX(1, hold->n);
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_ptr_prune(GMRFLib_ptr_tp *hold)
{
	if (hold) {
		if (hold->n_alloc - hold->n > IDX_ALLOC_INITIAL) {
			hold->ptr = Realloc(hold->ptr, IMAX(1, hold->n), void *);
			hold->n_alloc = IMAX(1, hold->n);
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idxval_prune(GMRFLib_idxval_tp *hold)
{
	if (hold) {
		if (hold->n_alloc - hold->n > IDX_ALLOC_INITIAL) {
			int n = IMAX(1, hold->n);
			hold->idx = Realloc(hold->idx, n, int);
			hold->val = Realloc(hold->val, n, double);
			hold->n_alloc = n;
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idxval_sort(GMRFLib_idxval_tp *hold)
{
	return GMRFLib_idxval_nsort(&hold, 1, 1);
}

int GMRFLib_idxval_nsort(GMRFLib_idxval_tp **hold, int n, int nt)
{
	return GMRFLib_idxval_nsort_x(hold, n, nt, 0, 1);
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((target_clones(INLA_CLONE_TARGETS "default")))
int GMRFLib_idxval_nsort_x_core(GMRFLib_idxval_tp *h, double *x, int prepare, int accumulate)
{
	// x is a test vector
	const int limit_merge = 8L, limit_sequential = 4L;
#if 0
	static int limit_merge = 0, limit_sequential = 0;
	if (!limit_merge)
		limit_merge = atoi(getenv("LIMIT_MERGE"));
	if (!limit_sequential)
		limit_sequential = atoi(getenv("LIMIT_SEQUENTIAL"));
#endif

	int debug = 0;
	if (GMRFLib_testit_mode && GMRFLib_testit_debug) {
		debug = 1;
	}

	if (!h || h->n <= 0) {
		return GMRFLib_SUCCESS;
	}
	// sort
	if (h->n > 1) {
		if (!GMRFLib_is_sorted_iinc(h->n, h->idx)) {
			my_sort2_id(h->idx, h->val, h->n);
		}
	}
	// unique
	if (h->n > 1) {
		int all_unique = 1;
		for (int j = 1; all_unique && j < h->n; j++) {
			all_unique = (h->idx[j] > h->idx[j - 1]);
		}

		if (!all_unique) {
			int k = 0;
			for (int j = 1; j < h->n; j++) {
				if (h->idx[j] != h->idx[k]) {
					k++;
					h->idx[k] = h->idx[j];
					h->val[k] = h->val[j];
				} else {
					if (accumulate) {
						h->val[k] += h->val[j];
					}
				}
			}
			if (debug && (h->n > k + 1)) {
				printf("Make unique: accumulate=[%1d], reduce length from %d to %dn", accumulate, h->n, k + 1);
			}
			h->n = k + 1;
		}
	}

	if (h->n <= GMRFLib_DOT_GROUP_NLIM || !prepare || !GMRFLib_internal_opt) {
		// check special cases: check if ddot or dsum can be used directly

		int is_sequential = 1;
		for (int i = 1; i < h->n && is_sequential; i++) {
			is_sequential = (h->idx[i] == h->idx[i - 1] + 1);
		}

		if (is_sequential) {
			int all_one = (h->val[0] == 1.0);
			for (int i = 1; i < h->n && all_one; i++) {
				all_one = (h->val[i] == 1.0);
			}

			if (all_one) {
				// special case, all_one and sequential, call _dsum directly
				h->dot_product_func = (GMRFLib_dot_product_tp *) GMRFLib_sparse_ddot_sum_;
			} else {
				// special case, sequential, call _ddot directly
				h->dot_product_func = (GMRFLib_dot_product_tp *) GMRFLib_sparse_ddot_ddot_;
			}
			h->preference = IDXVAL_SERIAL;
			return GMRFLib_SUCCESS;
		} else {
#if defined(INLA_WITH_ARMPL)
			armpl_status_t info = armpl_spvec_create_d(&(h->spvec), 0, h->idx[h->n - 1], h->n, h->idx, h->val, 0);
			assert(info == ARMPL_STATUS_SUCCESS);
			h->preference = IDXVAL_SERIAL_ARMPL;
#else
			h->preference = IDXVAL_SERIAL;
#endif
			h->dot_product_func = (GMRFLib_dot_product_tp *) NULL;	// GMRFLib_sparse_ddot_;
			return GMRFLib_SUCCESS;
		}
	}
	// an upper bound for the number of groups for memory allocation
	int ng = 1;
	int i = 1;
	while (i < h->n) {
		while (i < h->n && h->idx[i] == h->idx[i - 1] + 1)
			i++;
		while (i < h->n && h->idx[i] > h->idx[i - 1] + 1)
			i++;
		ng += 2;
	}

	int *g_istart = Calloc(ng, int);
	int *g_len = Calloc(ng, int);
	int *g_1 = Calloc(ng, int);
	int **g_idx = Calloc(ng, int *);
	double **g_val = Calloc(ng, double *);

	// collect sequential groups, starting from group=1, as group=0 is reserved for the irregular one
	ng = 1;
	g_istart[0] = 0;
	g_len[0] = 0;
	i = 1;
	while (i < h->n) {
		if (h->idx[i] == h->idx[i - 1] + 1) {
			while (i < h->n && h->idx[i] == h->idx[i - 1] + 1)
				i++;
			g_len[ng] = i - g_istart[ng];
			ng++;
		} else {
			while (i < h->n && h->idx[i] > h->idx[i - 1] + 1)
				i++;
			g_istart[ng] = i - 1;
		}
	}

	// keep only sequential groups with length >= limit_sequential
	int ggg = 1;
	for (int g = 1; g < ng; g++) {
		if (g_len[g] >= limit_sequential) {
			g_istart[ggg] = g_istart[g];
			g_len[ggg] = g_len[g];
			ggg++;
		}
	}
	ng = ggg;

	g_istart[ng] = h->n;
	g_len[ng] = 0;

	if (debug) {
		if (0) {
			for (i = 0; i < h->n; i++) {
				printf("idx[%1d] =  %1d\n", i, h->idx[i]);
			}
		}
		printf("ng = %1d\n", ng);
		for (int g = 1; g < ng; g++) {
			printf("group %1d start %1d len %1d\n", g, g_istart[g], g_len[g]);
			if (0) {
				for (i = 0; i < g_len[g]; i++) {
					printf("\t\t\tidx %1d\n", h->idx[g_istart[g] + i]);
				}
			}
		}
	}

	int irr_len = 1;
	int seq_len = 1;
	for (int g = 1; g < ng + 1; g++) {
		int istart = g_istart[g - 1] + g_len[g - 1];
		irr_len += g_istart[g] - istart + 1;
		seq_len += g_len[g];
	}

	int align_i = GMRFLib_MEM_ALIGN / sizeof(int);
	int align_d = GMRFLib_MEM_ALIGN / sizeof(double);
	int len_i = irr_len + seq_len + (ng + 1) * limit_merge + 2 * (ng + 1) * align_i;
	int len_d = irr_len + seq_len + (ng + 1) * limit_merge + 2 * (ng + 1) * align_d;
	int *new_idx = Calloc(len_i, int);
	double *new_val = Calloc(len_d, double);

	// build the irregular group
	int k = 0;
	for (int g = 1; g < ng + 1; g++) {
		int istart = g_istart[g - 1] + g_len[g - 1];
		int len = g_istart[g] - istart;
		if (len) {
			assert(k + len <= len_i);
			assert(k + len <= len_d);
			Memcpy(new_idx + k, h->idx + istart, len * sizeof(int));
			Memcpy(new_val + k, h->val + istart, len * sizeof(double));
			k += len;
		}
	}
	g_len[0] = k;
	g_idx[0] = new_idx;
	g_val[0] = new_val;

	if (debug) {
		if (0) {
			for (i = 0; i < g_len[0]; i++) {
				printf("nidx[%1d] %1d  val %g\n", i, new_idx[i], new_val[i]);
			}
		}
	}
	// copy each sequential group and pad for possible grouping
	int *seq_idx = new_idx + GMRFLib_align_len(irr_len, sizeof(int));
	double *seq_val = new_val + GMRFLib_align_len(irr_len, sizeof(double));

	k = 0;
	for (i = 1; i < ng; i++) {
		int istart = g_istart[i];
		int len = g_len[i];
		if (len) {
			Memcpy(seq_idx + k, h->idx + istart, len * sizeof(int));
			Memcpy(seq_val + k, h->val + istart, len * sizeof(double));
		}
		g_istart[i] = k;

		if (debug) {
			printf("group %d istart %d len %d limit_merge %d\n", i, istart, len,
			       (i < ng - 1 ? IMIN(limit_merge, h->idx[g_istart[i + 1]] - h->idx[istart + len - 1] - 1) : 0));
		}

		int pad = (i < ng - 1 ? IMIN(limit_merge, h->idx[g_istart[i + 1]] - h->idx[istart + len - 1] - 1) : 0);
		k += len;
		int offset = seq_idx[k - 1] + 1;
		for (int j = 0; j < pad; j++) {
			seq_idx[k + j] = offset + j;
		}
		k += pad;
		k = GMRFLib_align_len(k, sizeof(int));
	}

	// setup pointers to each sequential group
	for (int g = 1; g < ng; g++) {
		int istart = g_istart[g];
		g_idx[g] = seq_idx + istart;
		g_val[g] = seq_val + istart;
		g_len[g] *= -1;
	}

	if (debug) {
		if (0) {
			for (i = 0; i < k; i++) {
				printf("i %d new_idx %1d new_val %f\n", i, new_idx[i], new_val[i]);
			}
		}
	}
	// set g_1's
	g_1[0] = 0;
	for (int g = 0; g < ng; g++) {
		int all_one = 1;
		double *val = g_val[g];
		for (i = 0; all_one && i < IABS(g_len[g]); i++) {
			all_one = (val[i] == 1.0);
		}
		g_1[g] = all_one;
	}

	if (debug) {
		printf("NEW\nng = %1d\n", ng);
		for (int g = 0; g < ng; g++) {
			printf("group %1d start %1d len %1d g_1 %1d\n", g, g_istart[g], g_len[g], g_1[g]);
			if (0) {
				for (i = 0; i < IABS(g_len[g]); i++) {
					printf("\t\t\tidx %1d val %g\n", g_idx[g][i], g_val[g][i]);
				}
			}
		}
	}
	// remove groups with zero length
	int g = 0;
	while (1) {
		if (IABS(g_len[g]) == 0) {
			ng--;
			for (int gg = g; gg < ng; gg++) {
				g_len[gg] = g_len[gg + 1];
				g_istart[gg] = g_istart[gg + 1];
				g_1[gg] = g_1[gg + 1];
				g_idx[gg] = g_idx[gg + 1];
				g_val[gg] = g_val[gg + 1];
			}
		} else {
			g++;
		}

		if (g >= ng - 1) {
			break;
		}
	}

	if (debug) {
		printf("BEFORE MERGE\nng = %1d\n", ng);
		for (g = 0; g < ng; g++) {
			printf("group %1d start %1d len %1d g_1 %1d\n", g, g_istart[g], g_len[g], g_1[g]);
		}
	}
	// merge
	if (ng > 2) {
		g = 1;
		while (1) {
			int g1_end = seq_idx[g_istart[g] + IABS(g_len[g])];
			int g2_start = seq_idx[g_istart[g + 1]];
			if (g2_start - g1_end <= limit_merge && (g_1[g] == 0 && g_1[g + 1] == 0)) {
				g_len[g] = g_istart[g + 1] + IABS(g_len[g + 1]) - g_istart[g];
				for (int gg = g + 2; gg < ng; gg++) {
					g_istart[gg - 1] = g_istart[gg];
					g_len[gg - 1] = g_len[gg];
					g_1[gg - 1] = g_1[gg];
					g_idx[gg - 1] = g_idx[gg];
					g_val[gg - 1] = g_val[gg];
				}
				ng--;
			} else {
				g++;
			}
			if (g >= ng - 1) {
				break;
			}
		}

		// check for sequential. if so, set negative len
		for (g = 0; g < ng; g++) {
			if (g_len[g] > 0) {
				if (g_len[g] == g_idx[g][g_len[g] - 1] - g_idx[g][0] + 1) {
					g_len[g] *= -1;
				}
			}
		}
	}

	if (debug) {
		printf("AFTER MERGE\nng = %1d\n", ng);
		for (g = 0; g < ng; g++) {
			printf("group %1d start %1d len %1d g_1 %1d\n", g, g_istart[g], g_len[g], g_1[g]);
		}
	}

	h->g_n = ng;
	h->g_len = g_len;
	h->g_1 = g_1;
	h->g_idx = g_idx;
	h->g_val = g_val;
	h->g_n_mem = 2;
	h->g_mem = Calloc(h->g_n_mem, void *);
	h->g_mem[0] = (void *) new_idx;
	h->g_mem[1] = (void *) new_val;
	h->dot_product_func = NULL;
#if defined(INLA_WITH_ARMPL)
	h->spvec = NULL;
	h->spvec_g = NULL;
#endif
	Free(g_istart);

	int ntimes = 2;
	double treff[2] = { 0.0, 0.0 };
	double value[2] = { 0.0, 0.0 };

	int simple = (ng == 1 && g_len[0] < 0 ? 1 : 0);	       /* this is the simple one, where its either a sum or dot (no sparse dot) */
	double (*ddot_group)(GMRFLib_idxval_tp *, double *);
	ddot_group = (simple ? GMRFLib_sparse_ddot_group_simple_ : GMRFLib_sparse_ddot_group_);

#if defined(INLA_WITH_ARMPL)
	armpl_status_t info = armpl_spvec_create_d(&(h->spvec), 0, h->idx[h->n - 1], h->n, h->idx, h->val, 0);
	assert(info == ARMPL_STATUS_SUCCESS);
	if (!simple && (g_len[0] > 0 && g_1[0] == 0)) {
		info = armpl_spvec_create_d(&(h->spvec_g), 0, h->g_idx[0][h->g_len[0] - 1], h->g_len[0], h->g_idx[0], h->g_val[0], 0);
		assert(info == ARMPL_STATUS_SUCCESS);
	}
#endif
	for (int time = -ntimes; time < ntimes; time++) {
		if (time < 0) {
			GMRFLib_sparse_ddot_(h, x);
			ddot_group(h, x);
		} else {
			treff[0] -= GMRFLib_timer();
			value[0] = GMRFLib_sparse_ddot_(h, x);
			treff[0] += GMRFLib_timer();

			treff[1] -= GMRFLib_timer();
			value[1] = ddot_group(h, x);
			treff[1] += GMRFLib_timer();
		}
	}

	for (k = 0; k < 2; k++) {
		treff[k] /= (double) ntimes;
	}

	for (k = 1; k < 2; k++) {
		if (ABS(value[k] - value[0]) > 1000.0 * FLT_EPSILON * sqrt(h->n)) {
			P(ABS(value[k] - value[0]));
			P(k);
			P(value[0]);
			P(value[1]);

			printf("n %dn", h->n);
			for (i = 0; i < h->n; i++) {
				printf("\tidx[%1d] =  %1d  val = %g\n", i, h->idx[i], h->val[i]);
			}
			printf("ng %dn", h->g_n);
			if (0) {
				for (g = 0; g < h->g_n; g++) {
					printf("\tg = %d g_1 = %d\n", g, h->g_1[g]);
					for (i = 0; i < IABS(h->g_len[g]); i++) {
						printf("\t\tidx[%1d] =  %1d  val = %g\n", i, h->g_idx[g][i], h->g_val[g][i]);
					}
				}
			}

			P(ABS(value[k] - value[0]));
			P(k);
			P(value[0]);
			P(value[1]);

			assert(0 == 1);
		}
	}

	int kmin = -1;
	double tmin = GMRFLib_min_value(treff, 2, &kmin);

	if (debug) {
		double s = 1.0 / (DBL_EPSILON + treff[0] + treff[1]);
		printf("n= %1d chose kmin=%1d [serial= %.3f group= %.3f]\n", h->n, kmin, treff[0] * s, treff[1] * s);
	}

	switch (kmin) {
	case 0:
		h->preference = IDXVAL_SERIAL;
		h->dot_product_func = (GMRFLib_dot_product_tp *) NULL;	// GMRFLib_sparse_ddot_;
		break;
	case 1:
#if defined(INLA_WITH_ARMPL)
		h->preference = IDXVAL_GROUP_ARMPL;
#else
		h->preference = IDXVAL_GROUP;
#endif
		h->dot_product_func = (GMRFLib_dot_product_tp *) ddot_group;
		break;
	default:
		assert(0 == 1);
	}

	if (!GMRFLib_testit_mode) {
		if (h->preference == IDXVAL_SERIAL) {
			h->g_n = 0;
			Free(h->g_idx);
			Free(h->g_val);
			Free(h->g_len);
			Free(h->g_1);
			for (k = 0; k < h->g_n_mem; k++) {
				Free(h->g_mem[k]);
			}
			Free(h->g_mem);
			h->g_n_mem = 0;
#if defined(INLA_WITH_ARMPL)
			if (h->spvec_g) {
				armpl_status_t POSSIBLY_UNUSED(info) = armpl_spvec_destroy(h->spvec_g);
				h->spvec_g = NULL;
			}
#endif
		}
	}

	if (GMRFLib_dot_product_optim_report || GMRFLib_testit_mode) {
		int idx = 0;
		GMRFLib_CACHE_SET_IDX(idx);
		for (k = 0; k < 2; k++) {
			GMRFLib_dot_product_optim_report[idx][k] += treff[k];
		}
		GMRFLib_dot_product_optim_report[idx][2] += tmin;
		GMRFLib_dot_product_optim_report[idx][3 + kmin]++;	/* count... */
	}

	return GMRFLib_SUCCESS;
}
#pragma GCC diagnostic pop

int GMRFLib_idxval_prepare(GMRFLib_idxval_tp **hold, int n, int nt)
{
	return GMRFLib_idxval_nsort_x(hold, n, nt, 1, 1);
}

int GMRFLib_idxval_nsort_x(GMRFLib_idxval_tp **hold, int n, int nt, int prepare, int accumulate)
{
	int nmax = 1;
	for (int i = 0; i < n; i++) {
		GMRFLib_idxval_tp *h = hold[i];
		if (h->n) {
			nmax = IMAX(nmax, h->idx[h->n - 1] + 1);
		}
	}

	static double *x_ran = NULL;
	static int len_x_ran = 0;

	if (nmax > len_x_ran) {
#pragma omp critical (Name_5e61892bb7ffd3f396829c4740024f79b59c9da9)
		if (nmax > len_x_ran) {
			// do not free 'x_ran' as it might be in use. by purpose, we leak memory here
			if (0) {
				if (len_x_ran) {
					printf("[%s:%d] Leak %.2g Mb memory\n", __FILE__, __LINE__, len_x_ran * sizeof(double) / SQR(1024.0));
				}
			}

			int len = IMAX(2 * nmax, ISQR(256));
			double *xx = Calloc(len, double);
			for (int i = 0; i < len; i++) {
				xx[i] = GMRFLib_uniform();
			}
			x_ran = xx;
			len_x_ran = len;
		}
	}
#define CODE_BLOCK							\
	for(int k = 0; k < n; k++) {					\
		CODE_BLOCK_INIT();					\
		GMRFLib_idxval_nsort_x_core(hold[k], x_ran, prepare, accumulate); \
	}

	RUN_CODE_BLOCK(nt, 0, 0);
#undef CODE_BLOCK

	return GMRFLib_SUCCESS;
}

int GMRFLib_idx_free(GMRFLib_idx_tp *hold)
{
	if (hold) {
		Free(hold->idx);
		Free(hold);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_str_free(GMRFLib_str_tp *hold)
{
	if (hold) {
		for (int i = 0; i < hold->n; i++) {
			if (hold->str[i]) {
				Free(hold->str[i]);
			}
		}
		Free(hold->str);
		Free(hold);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idx2_free(GMRFLib_idx2_tp *hold)
{
	if (hold) {
		Free(hold->idx[0]);
		Free(hold->idx[1]);
		Free(hold->idx);
		Free(hold);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_val_free(GMRFLib_val_tp *hold)
{
	if (hold) {
		Free(hold->val);
		Free(hold);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_ptr_free(GMRFLib_ptr_tp *hold)
{
	if (hold) {
		Free(hold->ptr);
		Free(hold);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idxval_free(GMRFLib_idxval_tp *hold)
{
	if (hold) {
		Free(hold->idx);
		Free(hold->val);
		Free(hold->g_idx);
		Free(hold->g_val);
		Free(hold->g_len);
		Free(hold->g_1);
		for (int i = 0; i < hold->g_n_mem; i++) {
			Free(hold->g_mem[i]);
		}
		if (hold->g_mem) {
			Free(hold->g_mem);
		}
#if defined(INLA_WITH_ARMPL)
		if (hold->spvec) {
			armpl_status_t info = armpl_spvec_destroy(hold->spvec);
			assert(info == ARMPL_STATUS_SUCCESS);
		}
		if (hold->spvec_g) {
			armpl_status_t info = armpl_spvec_destroy(hold->spvec_g);
			assert(info == ARMPL_STATUS_SUCCESS);
		}
#endif
		Free(hold);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idx_add(GMRFLib_idx_tp **hold, int idx)
{
	if (*hold == NULL) {
		GMRFLib_idx_create(hold);
	}
	if ((*hold)->n == (*hold)->n_alloc) {
		(*hold)->n_alloc += IDX_ALLOC_INCREASE;
		(*hold)->idx = Realloc((*hold)->idx, (*hold)->n_alloc, int);
	}
	(*hold)->idx[(*hold)->n] = idx;
	(*hold)->n++;

	return GMRFLib_SUCCESS;
}

int GMRFLib_idx_nadd(GMRFLib_idx_tp **hold, int n, int *idx)
{
	if (*hold == NULL) {
		GMRFLib_idx_create_x(hold, n);
	}
	if ((*hold)->n + n > (*hold)->n_alloc) {
		(*hold)->n_alloc += IMAX(n, IDX_ALLOC_INCREASE);
		(*hold)->idx = Realloc((*hold)->idx, (*hold)->n_alloc, int);
		assert((*hold)->idx);
	}
	Memcpy((*hold)->idx + (*hold)->n, idx, n * sizeof(int));
	(*hold)->n += n;

	return GMRFLib_SUCCESS;
}

int GMRFLib_idx_find(int id, GMRFLib_idx_tp *h)
{
	// return index for which h->idx[i]=id, -1 if none
	for (int i = 0; i < h->n; i++) {
		if (h->idx[i] == id) {
			return i;
		}
	}
	return -1;
}

int GMRFLib_idxval_find(int *id, double *val, GMRFLib_idxval_tp *h)
{
	// return index for which h->idx[i]=id and/or h->val[i]=val, -1 if none

	if (id && !val) {
		for (int i = 0; i < h->n; i++) {
			if (h->idx[i] == *id) {
				return i;
			}
		}
	} else if (!id && val) {
		for (int i = 0; i < h->n; i++) {
			if (ISEQUAL(h->val[i], *val)) {
				return i;
			}
		}
	} else if (id && val) {
		for (int i = 0; i < h->n; i++) {
			if (h->idx[i] == *id && ISEQUAL(h->val[i], *val)) {
				return i;
			}
		}
	} else {
		return -1;
	}

	return -1;
}

GMRFLib_idx_tp *GMRFLib_idx_duplicate(GMRFLib_idx_tp *h)
{
	if (!h)
		return NULL;
	GMRFLib_idx_tp *nnew = NULL;
	GMRFLib_idx_create_x(&nnew, IMAX(1, h->n));
	if (h->n > 0) {
		GMRFLib_idx_nadd(&nnew, h->n, h->idx);
	}

	return nnew;
}

GMRFLib_idx2_tp *GMRFLib_idx2_duplicate(GMRFLib_idx2_tp *h)
{
	if (!h)
		return NULL;
	GMRFLib_idx2_tp *nnew = NULL;
	GMRFLib_idx2_create_x(&nnew, IMAX(1, h->n));
	if (h->n > 0) {
		Memcpy(nnew->idx[0], h->idx[0], h->n * sizeof(int));
		Memcpy(nnew->idx[1], h->idx[1], h->n * sizeof(int));
		nnew->n = h->n;
	}

	return nnew;
}

int GMRFLib_idx_overlap(GMRFLib_idx_tp *idx1, GMRFLib_idx_tp *idx2)
{
	// check if any indices are the same, if so return 1, else return 0.
	// ASSUME IDX1 and IDX2 are SORTED!

	if (idx1 == NULL || idx2 == NULL || idx1->n == 0 || idx2->n == 0) {
		return 0;
	}

	const int debug = 1;
	if (debug) {
		FIXME("idx_overlap: check if sorted as debug = 1");
		assert(GMRFLib_is_sorted_iinc(idx1->n, idx1->idx));
		assert(GMRFLib_is_sorted_iinc(idx2->n, idx2->idx));
	}
	// this is fast check MAX(idx1) < MIN(idx2) or opposite
	if ((idx1->idx[idx1->n - 1] < idx2->idx[0]) || (idx2->idx[idx2->n - 1] < idx1->idx[0])) {
		return 0;
	}

	int i1 = 0, i2 = 0;
	while (1) {
		if (idx1->idx[i1] < idx2->idx[i2]) {
			i1++;
		} else if (idx1->idx[i1] > idx2->idx[i2]) {
			i2++;
		} else {
			return 1;			       // then they must be equal
		}

		if (i1 == idx1->n || i2 == idx2->n) {
			return 0;
		}
	}
	return 0;
}

int GMRFLib_idxval_overlap(GMRFLib_idxval_tp *idx1, GMRFLib_idxval_tp *idx2)
{
	if (idx1 == NULL || idx2 == NULL) {
		return 0;
	}

	GMRFLib_idx_tp v, vv;
	v.n = idx1->n;
	v.idx = idx1->idx;
	vv.n = idx2->n;
	vv.idx = idx2->idx;

	return GMRFLib_idx_overlap(&v, &vv);
}

int GMRFLib_str_add(GMRFLib_str_tp **hold, char *s)
{
	if (*hold == NULL) {
		GMRFLib_str_create(hold);
	}
	if ((*hold)->n == (*hold)->n_alloc) {
		(*hold)->n_alloc += IDX_ALLOC_INCREASE;
		(*hold)->str = Realloc((*hold)->str, (*hold)->n_alloc, char *);
	}
	(*hold)->str[(*hold)->n] = Strdup(s);
	(*hold)->n++;

	return GMRFLib_SUCCESS;
}



int GMRFLib_str_is_member(GMRFLib_str_tp *hold, char *s, int case_sensitive, int *idx_match)
{
	if (hold == NULL) {
		return 0;
	}

	int (*cmp)(const char *, const char *) =(case_sensitive ? strcmp : strcasecmp);
	for (int i = 0; i < hold->n; i++) {
		if (cmp(s, hold->str[i]) == 0) {
			if (idx_match) {
				*idx_match = i;
			}
			return 1;
		}
	}
	return 0;
}

int GMRFLib_idx2_add(GMRFLib_idx2_tp **hold, int idx0, int idx1)
{
	if (*hold == NULL) {
		GMRFLib_idx2_create(hold);
	}
	if ((*hold)->n == (*hold)->n_alloc) {
		(*hold)->n_alloc += IDX_ALLOC_INCREASE;
		(*hold)->idx[0] = Realloc((*hold)->idx[0], (*hold)->n_alloc, int);
		(*hold)->idx[1] = Realloc((*hold)->idx[1], (*hold)->n_alloc, int);
	}
	(*hold)->idx[0][(*hold)->n] = idx0;
	(*hold)->idx[1][(*hold)->n] = idx1;
	(*hold)->n++;

	return GMRFLib_SUCCESS;
}

int GMRFLib_val_add(GMRFLib_val_tp **hold, double val)
{
	if (*hold == NULL) {
		GMRFLib_val_create(hold);
	}
	if ((*hold)->n == (*hold)->n_alloc) {
		(*hold)->n_alloc += IDX_ALLOC_INCREASE;
		(*hold)->val = Realloc((*hold)->val, (*hold)->n_alloc, double);
	}
	(*hold)->val[(*hold)->n] = val;
	(*hold)->n++;

	return GMRFLib_SUCCESS;
}

int GMRFLib_ptr_add(GMRFLib_ptr_tp **hold, void *ptr)
{
	if (*hold == NULL) {
		GMRFLib_ptr_create(hold);
	}
	if ((*hold)->n == (*hold)->n_alloc) {
		(*hold)->n_alloc += IDX_ALLOC_INCREASE;
		(*hold)->ptr = Realloc((*hold)->ptr, (*hold)->n_alloc, void *);
	}
	(*hold)->ptr[(*hold)->n] = ptr;
	(*hold)->n++;

	return GMRFLib_SUCCESS;
}

int GMRFLib_idxval_add(GMRFLib_idxval_tp **hold, int idx, double val)
{
	if (*hold == NULL) {
		GMRFLib_idxval_create(hold);
	}
	if ((*hold)->n == (*hold)->n_alloc) {
		(*hold)->n_alloc += IDX_ALLOC_INCREASE;
		(*hold)->idx = Realloc((*hold)->idx, (*hold)->n_alloc, int);
		(*hold)->val = Realloc((*hold)->val, (*hold)->n_alloc, double);
	}
	(*hold)->idx[(*hold)->n] = idx;
	(*hold)->val[(*hold)->n] = val;
	(*hold)->n++;

	return GMRFLib_SUCCESS;
}

int GMRFLib_idxval_addto(GMRFLib_idxval_tp **hold, int idx, double val)
{
	// if idx exists before, add val to value , otherwise just 'add'.
	// if there are two entries of 'idx', then only the first is used.

	if (*hold == NULL) {
		GMRFLib_idxval_create(hold);
	}

	int i;

	// hopefully this is it.
	i = (*hold)->iaddto;
	if ((*hold)->idx[i] == idx) {
		(*hold)->val[i] += val;
		return GMRFLib_SUCCESS;
	}
	// FIXME: this should be improved in general, but I think for the usage its ok. Since we are likely to add with same or increasing idx,
	// then I added this 'iaddto' which recall the last index, and try to be a little smarter.
	for (i = (*hold)->iaddto + 1; i < (*hold)->n; i++) {
		if ((*hold)->idx[i] == idx) {
			(*hold)->val[i] += val;
			(*hold)->iaddto = i;
			return GMRFLib_SUCCESS;
		}
	}
	for (i = 0; i < (*hold)->iaddto; i++) {
		if ((*hold)->idx[i] == idx) {
			(*hold)->val[i] += val;
			(*hold)->iaddto = i;
			return GMRFLib_SUCCESS;
		}
	}

	if ((*hold)->n == (*hold)->n_alloc) {
		(*hold)->n_alloc += IDX_ALLOC_INCREASE;
		(*hold)->idx = Realloc((*hold)->idx, (*hold)->n_alloc, int);
		(*hold)->val = Realloc((*hold)->val, (*hold)->n_alloc, double);
	}
	(*hold)->idx[(*hold)->n] = idx;
	(*hold)->val[(*hold)->n] = val;
	(*hold)->n++;

	return GMRFLib_SUCCESS;
}

int GMRFLib_idx_split_free(GMRFLib_ptr_tp *ptr)
{
	if (!ptr) {
		return GMRFLib_SUCCESS;
	}
	for (int i = 0; i < ptr->n; i++) {
		GMRFLib_idx_free((GMRFLib_idx_tp *) (ptr->ptr[i]));
	}
	GMRFLib_ptr_free(ptr);
	return GMRFLib_SUCCESS;
}

GMRFLib_ptr_tp *GMRFLib_idx_split(GMRFLib_idx_tp *sel, int size)
{
	// split IDX in SEL, into groups of SIZE.
	// return a _ptr_tp of GMRFLib_idx_tp's, which can be free'd with _idx_split_free

	if (size <= 0 || !sel || sel->n <= 0) {
		return NULL;
	}

	div_t d = div(sel->n, size);
	int NN = d.quot + (d.rem > 0 ? 1 : 0);
	GMRFLib_ptr_tp *ptr = NULL;
	GMRFLib_ptr_create_x(&ptr, NN);

	for (int i = 0; i < d.quot; i++) {
		GMRFLib_idx_tp *idx = NULL;
		GMRFLib_idx_create_x(&idx, size);
		GMRFLib_idx_nadd(&idx, size, sel->idx + i * size);
		GMRFLib_ptr_add(&ptr, (void *) idx);
	}
	if (d.rem > 0) {
		GMRFLib_idx_tp *idx = NULL;
		GMRFLib_idx_create_x(&idx, d.rem);
		GMRFLib_idx_nadd(&idx, d.rem, sel->idx + d.quot * size);
		GMRFLib_ptr_add(&ptr, (void *) idx);
	}
	assert(ptr->n == NN);

	return ptr;
}
