
/* idxval.c
 * 
 * Copyright (C) 2022-2022 Havard Rue
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

#ifndef GITCOMMIT
#define GITCOMMIT
#endif

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-const-variable"
static const char GitID[] = "file: " __FILE__ "  " GITCOMMIT;
#pragma GCC diagnostic pop

#if !defined(__FreeBSD__)
#include <assert.h>
#include <float.h>
#include <signal.h>
#include <stdarg.h>
#include <math.h>
#include <strings.h>
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stddef.h>
#include <string.h>
#include <malloc.h>
#endif
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"
#include "GMRFLib/hashP.h"

#define IDX_ALLOC_INITIAL 32
#define IDX_ALLOC_ADD     512

// set to 0 to keep the groups anyway
#define IDXVAL_FREE_GROUPS_IF_NOT_BEST() 1

int GMRFLib_idx_create(GMRFLib_idx_tp ** hold)
{
	return GMRFLib_idx_create_x(hold, IDX_ALLOC_INITIAL);
}

int GMRFLib_idx_create_x(GMRFLib_idx_tp ** hold, int len)
{
	len = IMAX(1, len);
	*hold = Calloc(1, GMRFLib_idx_tp);
	(*hold)->idx = Calloc(len, int);
	(*hold)->n_alloc = len;
	(*hold)->n = 0;

	return GMRFLib_SUCCESS;
}

int GMRFLib_str_create(GMRFLib_str_tp ** hold)
{
	return GMRFLib_str_create_x(hold, IDX_ALLOC_INITIAL);
}

int GMRFLib_str_create_x(GMRFLib_str_tp ** hold, int len)
{
	len = IMAX(1, len);
	*hold = Calloc(1, GMRFLib_str_tp);
	(*hold)->str = Calloc(len, char *);
	(*hold)->n_alloc = len;
	(*hold)->n = 0;

	return GMRFLib_SUCCESS;
}

int GMRFLib_idx2_create(GMRFLib_idx2_tp ** hold)
{
	return GMRFLib_idx2_create_x(hold, IDX_ALLOC_INITIAL);
}

int GMRFLib_idx2_create_x(GMRFLib_idx2_tp ** hold, int len)
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

int GMRFLib_val_create(GMRFLib_val_tp ** hold)
{
	*hold = Calloc(1, GMRFLib_val_tp);
	(*hold)->val = Calloc(IDX_ALLOC_INITIAL, double);
	(*hold)->n_alloc = IDX_ALLOC_INITIAL;
	(*hold)->n = 0;

	return GMRFLib_SUCCESS;
}

int GMRFLib_idxval_create(GMRFLib_idxval_tp ** hold)
{
	return GMRFLib_idxval_create_x(hold, IDX_ALLOC_INITIAL);
}

int GMRFLib_idxval_create_x(GMRFLib_idxval_tp ** hold, int len)
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

GMRFLib_idxval_tp **GMRFLib_idxval_ncreate_x(int n, int len)
{
	if (n > 0) {
		GMRFLib_idxval_tp **a = Calloc(n, GMRFLib_idxval_tp *);
		for (int i = 0; i < n; i++) {
			GMRFLib_idxval_create_x(&(a[i]), len);
		}
		return a;
	} else {
		return NULL;
	}
}

int GMRFLib_idx_printf(FILE * fp, GMRFLib_idx_tp * hold, const char *msg)
{
	if (hold) {
		fprintf(fp, "[%s] n = %1d  nalloc = %1d\n", msg, hold->n, hold->n_alloc);
		for (int i = 0; i < hold->n; i++) {
			fprintf(fp, "\tidx[%1d] = %1d\n", i, hold->idx[i]);
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_str_printf(FILE * fp, GMRFLib_str_tp * hold, const char *msg)
{
	if (hold) {
		fprintf(fp, "[%s] n = %1d  nalloc = %1d\n", msg, hold->n, hold->n_alloc);
		for (int i = 0; i < hold->n; i++) {
			fprintf(fp, "\tidx[%1d] = %s\n", i, hold->str[i]);
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idx2_printf(FILE * fp, GMRFLib_idx2_tp * hold, const char *msg)
{
	if (hold) {
		fprintf(fp, "[%s] n = %1d  nalloc = %1d\n", msg, hold->n, hold->n_alloc);
		for (int i = 0; i < hold->n; i++) {
			fprintf(fp, "\tidx[][%1d] = %1d %1d\n", i, hold->idx[0][i], hold->idx[1][i]);
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_val_printf(FILE * fp, GMRFLib_val_tp * hold, const char *msg)
{
	if (hold) {
		fprintf(fp, "[%s] n = %1d  nalloc = %1d\n", msg, hold->n, hold->n_alloc);
		for (int i = 0; i < hold->n; i++) {
			fprintf(fp, "\tval[%1d] = %g\n", i, hold->val[i]);
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idxval_printf(FILE * fp, GMRFLib_idxval_tp * hold, const char *msg)
{
	if (hold) {
		fprintf(fp, "[%s] n = %1d  nalloc = %1d iaddto = %1d\n", msg, hold->n, hold->n_alloc, hold->iaddto);
		for (int i = 0; i < hold->n; i++) {
			fprintf(fp, "\t(idx, val)[%1d] = (%d, %g)\n", i, hold->idx[i], hold->val[i]);
		}
		if (hold->g_i) {
			fprintf(fp, "\tg_n = %1d\n", hold->g_n);
			for (int g = 0; g < hold->g_n; g++) {
				fprintf(fp, "\tgroup %d has length %d and start at index %d (one=%s)\n", g, hold->g_len[g], hold->g_i[g],
					(hold->g_1 && hold->g_1[g] ? "TRUE" : "FALSE"));
				fprintf(fp, "\t\t");
				for (int k = 0; k < IABS(hold->g_len[g]); k++) {
					fprintf(fp, " %1d", hold->idx[hold->g_i[g] + k]);
				}
				fprintf(fp, "\n");
			}
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idxval_info_printf(FILE * fp, GMRFLib_idxval_tp * hold, const char *msg)
{
	if (hold) {
		fprintf(fp, "[%s] n=%1d  nalloc=%1d iaddto=%1d\n", (msg ? msg : ""), hold->n, hold->n_alloc, hold->iaddto);
		if (hold->g_i) {
			for (int g = 0; g < hold->g_n; g++) {
				fprintf(fp, "\tgroup %1d/%1d: len=%1d one=%s\n", g, hold->g_n, hold->g_len[g],
					(hold->g_1 && hold->g_1[g] ? "TRUE" : "FALSE"));
			}
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idx_nprune(GMRFLib_idx_tp ** a, int n)
{
	if (a) {
		for (int i = 0; i < n; i++) {
			GMRFLib_idx_prune(a[i]);
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_str_nprune(GMRFLib_str_tp ** a, int n)
{
	if (a) {
		for (int i = 0; i < n; i++) {
			GMRFLib_str_prune(a[i]);
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idx_sort(GMRFLib_idx_tp * hold)
{
	if (hold) {
		qsort((void *) hold->idx, (size_t) hold->n, sizeof(int), GMRFLib_icmp);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idx_nsort(GMRFLib_idx_tp ** a, int n, int nt)
{
#define CODE_BLOCK							\
	for(int i = 0; i < n; i++) {					\
		if (a[i] && a[i]->n > 1) {				\
			qsort((void *) a[i]->idx, (size_t) a[i]->n,  sizeof(int), GMRFLib_icmp); \
		}							\
	}

	RUN_CODE_BLOCK(nt, 0, 0);
#undef CODE_BLOCK

	return GMRFLib_SUCCESS;
}

int GMRFLib_idx_uniq(GMRFLib_idx_tp * hold)
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

int GMRFLib_idx_nuniq(GMRFLib_idx_tp ** a, int n, int nt)
{
#define CODE_BLOCK				\
	for (int i = 0; i < n; i++) {		\
		GMRFLib_idx_uniq(a[i]);		\
	}

	RUN_CODE_BLOCK(nt, 0, 0);
#undef CODE_BLOCK

	return GMRFLib_SUCCESS;
}

int GMRFLib_idx2_nprune(GMRFLib_idx2_tp ** a, int n)
{
	if (a) {
		for (int i = 0; i < n; i++) {
			GMRFLib_idx2_prune(a[i]);
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_val_nprune(GMRFLib_val_tp ** a, int n)
{
	if (a) {
		for (int i = 0; i < n; i++) {
			GMRFLib_val_prune(a[i]);
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idxval_nuniq(GMRFLib_idxval_tp ** a, int n, int nt)
{
#define CODE_BLOCK					\
	for (int i = 0; i < n; i++) {			\
		GMRFLib_idxval_uniq(a[i]);		\
	}

	RUN_CODE_BLOCK(nt, 0, 0);
#undef CODE_BLOCK
	return GMRFLib_SUCCESS;
}

int GMRFLib_idxval_uniq(GMRFLib_idxval_tp * hold)
{
	// sort idx, and accumulate values and then prune
	if (hold && hold->n > 1) {
		int i, j;

		GMRFLib_idxval_sort(hold);
		for (j = 0, i = 0; i < hold->n; i++) {
			if (hold->idx[j] == hold->idx[i]) {
				if (i > j) {
					hold->val[j] += hold->val[i];
				}
			} else {
				j++;
				hold->idx[j] = hold->idx[i];
				hold->val[j] = hold->val[i];
			}
		}
		hold->n = j + 1;
		GMRFLib_idxval_prune(hold);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idxval_nprune(GMRFLib_idxval_tp ** a, int n, int nt)
{
#define CODE_BLOCK					\
	for (int i = 0; i < n; i++) {			\
		GMRFLib_idxval_prune(a[i]);		\
	}

	RUN_CODE_BLOCK(nt, 0, 0);
#undef CODE_BLOCK

	return GMRFLib_SUCCESS;
}

int GMRFLib_idx_prune(GMRFLib_idx_tp * hold)
{
	if (hold) {
		hold->idx = Realloc(hold->idx, IMAX(1, hold->n), int);
		hold->n_alloc = IMAX(1, hold->n);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_str_prune(GMRFLib_str_tp * hold)
{
	if (hold) {
		hold->str = Realloc(hold->str, IMAX(1, hold->n), char *);
		hold->n_alloc = IMAX(1, hold->n);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idx2_prune(GMRFLib_idx2_tp * hold)
{
	if (hold) {
		hold->idx[0] = Realloc(hold->idx[0], IMAX(1, hold->n), int);
		hold->idx[1] = Realloc(hold->idx[1], IMAX(1, hold->n), int);
		hold->n_alloc = IMAX(1, hold->n);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_val_prune(GMRFLib_val_tp * hold)
{
	if (hold) {
		hold->val = Realloc(hold->val, IMAX(1, hold->n), double);
		hold->n_alloc = IMAX(1, hold->n);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idxval_prune(GMRFLib_idxval_tp * hold)
{
	if (hold) {
		int n = IMAX(1, hold->n);
		hold->idx = Realloc(hold->idx, n, int);
		hold->val = Realloc(hold->val, n, double);
		hold->n_alloc = n;
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idxval_sort(GMRFLib_idxval_tp * hold)
{
	return GMRFLib_idxval_nsort(&hold, 1, 1);
}

int GMRFLib_idxval_nsort(GMRFLib_idxval_tp ** hold, int n, int nt)
{
	return GMRFLib_idxval_nsort_x(hold, n, nt, 0);
}

int GMRFLib_idxval_nsort_x(GMRFLib_idxval_tp ** hold, int n, int nt, int prune_zeros)
{
	GMRFLib_ENTER_ROUTINE;

	// prune_zeros, 0=do nothing, >0 remove duplicate with same index within each group, <0 remove all zeros
	// prune_zeros and 'group'-indexing is limited to vectors with length > 'limit'
	const int limit = 16;
	int debug = GMRFLib_DEBUG_IF_TRUE();

	int nmax = 1;
	for (int i = 0; i < n; i++) {
		GMRFLib_idxval_tp *h = hold[i];
		nmax = IMAX(nmax, h->n);
	}


#define CODE_BLOCK							\
	for (int i = 0; i < n; i++) {					\
		GMRFLib_idxval_tp *h = hold[i];				\
		if (!h) continue;					\
		double tref[5] =  {0, 0, 0, 0, 0};			\
		if (debug) tref[0] -= GMRFLib_cpu();			\
		if (h->n > 1) {						\
			int is_sorted = 1;				\
			for(int j = 1; j < h->n && is_sorted; j++) {	\
				is_sorted = is_sorted && (h->idx[j] >= h->idx[j-1]); \
			}						\
			if (!is_sorted) {				\
				if (0) {				\
					GMRFLib_qsorts((void *) h->idx, (size_t) h->n, sizeof(int), \
						       (void *) h->val, sizeof(double), NULL, (size_t) 0, GMRFLib_icmp); \
				} else {				\
					double *idx_tmp = CODE_BLOCK_WORK_PTR(0); \
					for(int j = 0; j < h->n; j++) {	\
						idx_tmp[j] = (double) h->idx[j]; \
					}				\
					gsl_sort2(idx_tmp, (size_t) 1, h->val, (size_t) 1, (size_t) h->n); \
					for(int j = 0; j < h->n; j++) {	\
						h->idx[j] = (int) idx_tmp[j]; \
					}				\
				}					\
			}						\
		}							\
		if (debug) tref[0] += GMRFLib_cpu();			\
									\
		if (h->n <= limit) {					\
			h->preference = IDXVAL_SERIAL;			\
			continue;					\
		}							\
									\
		/*							\
		 * build basic groups with one group for each sequence and then one for each individual \
		 */							\
		if (debug) tref[1] -= GMRFLib_cpu();			\
		int ng = 1;						\
		for (int j = 1; j < h->n; j++) {			\
			if (h->idx[j] != h->idx[j - 1] + 1) {		\
				ng++;					\
			}						\
		}							\
		int *g_i = Calloc(ng + 1, int);				\
		int *g_len = Calloc(ng + 1, int);			\
									\
		int k = 0;						\
		g_i[0] = 0;						\
		for (int j = 1; j < h->n; j++) {			\
			if (h->idx[j] != h->idx[j - 1] + 1) {		\
				g_len[k] = j - g_i[k];			\
				k++;					\
				g_i[k] = j;				\
			}						\
		}							\
		g_len[ng - 1] = h->n - g_i[ng - 1];			\
		if (debug) tref[1] += GMRFLib_cpu();			\
									\
		if (1 && debug) {					\
			GMRFLib_idxval_printf(stdout, h, "OLD");	\
		}							\
									\
		k = h->n;				       /* so the code below works */ \
		int kg = ng;				       /* so the code below works */ \
									\
		if (debug) tref[2] -= GMRFLib_cpu();			\
		if (prune_zeros) {					\
			/*						\
			 * remove zeros in 'val', either all or all but the first one \
			 */						\
			k = 0;						\
			kg = 0;						\
			int nskip = 0;					\
			for (int g = 0; g < ng; g++) {			\
									\
				int gl = IABS(g_len[g]);		\
				int num_zero = 0;			\
				int idx_zz = -1;			\
				int gn = 0;				\
				int kk = k;				\
									\
				for (int j = 0; j < gl; j++) {		\
					int idx_z = h->idx[g_i[g] + j];	\
									\
					if (prune_zeros < 0) {		\
						if (ISZERO(h->val[g_i[g] + j])) { \
							nskip++;	\
						} else {		\
							h->idx[k] = idx_z; \
							h->val[k] = h->val[g_i[g] + j];	\
							k++;		\
							gn++;		\
						}			\
					} else {			\
						if (idx_z == idx_zz) {	\
							if (ISZERO(h->val[g_i[g] + j])) { \
								num_zero++; \
								if (prune_zeros > 0 && num_zero == 1) {	\
									h->idx[k] = idx_z; \
									h->val[k] = h->val[g_i[g] + j];	\
									k++; \
									gn++; \
								} else { \
									nskip++; \
								}	\
							} else {	\
								h->idx[k] = idx_z; \
								h->val[k] = h->val[g_i[g] + j];	\
								k++;	\
								gn++;	\
							}		\
						} else {		\
							h->idx[k] = idx_z; \
							h->val[k] = h->val[g_i[g] + j];	\
							k++;		\
							gn++;		\
							num_zero = 0;	\
							if (ISZERO(h->val[g_i[g] + j])) { \
								num_zero++; \
							}		\
						}			\
						idx_zz = idx_z;		\
					}				\
				}					\
				if (debug && (g_len[g] != gn)) {	\
					printf("modify group %d from %d to %d, nskip=%1d\n", g, gl, gn, nskip); \
				}					\
				if (gn) {				\
					g_i[kg] = kk;			\
					g_len[kg] = (g_len[g] < 0 ? -gn : gn); \
					kg++;				\
				}					\
			}						\
		}							\
		h->n = k;						\
		h->g_n = kg;						\
		h->g_i = g_i;						\
		h->g_len = g_len;					\
									\
		if (debug) tref[2] += GMRFLib_cpu();			\
									\
		if (1 && debug) {					\
			GMRFLib_idxval_printf(stdout, h, "NEW");	\
		}							\
									\
		/*							\
		 * grow the groups together by merging individual groups, and merging sequential groups \
		 */							\
		if (debug) tref[3] -= GMRFLib_cpu();			\
		int gg = 0;						\
		int irregular = 0;					\
		for (int g = 0; g < h->g_n; g++) {			\
			if (h->g_len[g] == 1) {				\
				if (irregular) {			\
					h->g_len[gg]++;			\
				} else {				\
					h->g_len[gg] = 1;		\
					h->g_i[gg] = h->g_i[g];		\
				}					\
				irregular = 1;				\
			} else if (IABS(h->g_len[g]) > 0) {		\
				if (irregular) {			\
					gg++;				\
				}					\
				irregular = 0;				\
				h->g_len[gg] = -IABS(h->g_len[g]);	\
				h->g_i[gg] = h->g_i[g];			\
				gg++;					\
			} else {					\
				continue;				\
			}						\
		}							\
		h->g_n = gg + irregular;				\
									\
		int fill = 8L;						\
		if (h && h->n > 1) {					\
			int len_extra = 0;				\
			for (int g = 0; g < h->g_n - 1; g++) {		\
				if (h->g_len[g] < 0 && h->g_len[g + 1] < 0) { \
					int last_this = h->idx[h->g_i[g] + IABS(h->g_len[g]) - 1]; \
					int first_next = h->idx[h->g_i[g + 1]];	\
					if (first_next - last_this < fill) { \
						len_extra += first_next - last_this - 1; \
						if (debug) {		\
							printf("\tmerge group %1d and %1d with dist %1d and lengths %1d and %1d\n", \
							       g, g + 1, first_next - last_this, h->g_len[g], h->g_len[g + 1]);	\
						}			\
					}				\
				}					\
			}						\
									\
			if (len_extra > 0) {				\
									\
				int n_new = h->n + len_extra;		\
				int *idx_new = Calloc(n_new, int);	\
				double *val_new = Calloc(n_new, double); \
									\
				int kk = 0;				\
				int gg = 0;				\
									\
				int this_len = h->g_len[gg];		\
				int this_alen = IABS(this_len);		\
									\
				Memcpy(idx_new + kk, h->idx, this_alen * sizeof(int)); \
				Memcpy(val_new + kk, h->val, this_alen * sizeof(double)); \
									\
				kk += this_alen;			\
				h->g_i[gg] = h->g_i[0];			\
				h->g_len[gg] = h->g_len[0];		\
									\
				int gg_len = h->g_len[gg];		\
				int gg_alen = IABS(gg_len);		\
				int gg_first_i = h->g_i[gg];		\
				int gg_last_i = gg_first_i + gg_alen - 1; \
				int gg_last_idx = idx_new[gg_last_i];	\
				int pending = 0;			\
									\
				if (debug) {				\
					printf("n=%1d len_extra=%1d n_new=%1d\n", h->n, len_extra, n_new); \
				}					\
									\
				for (int g = 1; g < h->g_n; g++) {	\
									\
					if (debug) {			\
						printf("\tprocess group=%1d with len=%d kk=%1d\n", g, h->g_len[g], kk);	\
					}				\
									\
					/*				\
					 * merge this group into the current new one or just add it ? \
					 */				\
					int g_len = h->g_len[g];	\
					int g_alen = IABS(g_len);	\
					int g_first_i = h->g_i[g];	\
					int g_first_idx = h->idx[g_first_i]; \
									\
					int idx_diff = g_first_idx - gg_last_idx; \
					if (debug) {			\
						printf("g_len=%1d gg_len=%1d idx_diff=%1d\n", g_len, gg_len, idx_diff);	\
					}				\
									\
					if (g_len < 0 && gg_len < 0 && idx_diff < fill) { \
						if (debug) {		\
							printf("\tmerge group g=%1d into gg=%1d\n", g, gg); \
						}			\
									\
						/*			\
						 * yes, merge it	\
						 */			\
						int len_fill = idx_diff - 1; \
						for (int j = 0; j < len_fill; j++) { \
							idx_new[kk + j] = idx_new[kk - 1] + 1 + j; \
						}			\
						if (debug) {		\
							printf("\t\tlen_fill=%1d\n", len_fill);	\
						}			\
									\
						h->g_len[gg] -= len_fill; \
						kk += len_fill;		\
						Memcpy(idx_new + kk, h->idx + g_first_i, g_alen * sizeof(int));	\
						Memcpy(val_new + kk, h->val + g_first_i, g_alen * sizeof(double)); \
						h->g_len[gg] -= g_alen;	\
						kk += g_alen;		\
									\
						gg_len = h->g_len[gg];	\
						gg_alen = IABS(gg_len);	\
						gg_first_i = h->g_i[gg]; \
						gg_last_i = gg_first_i + gg_alen - 1; \
						gg_last_idx = idx_new[gg_last_i]; \
									\
						if (debug) {		\
							printf("\t\tappend group g=%1d with len %1d\n", g, g_len); \
							printf("\t\tgg_last_idx=%1d kk=%1d\n", gg_last_idx, kk); \
						}			\
						pending = 1;		\
					} else {			\
						/*			\
						 * no, just add a new group \
						 */			\
						gg++;			\
						Memcpy(idx_new + kk, h->idx + g_first_i, g_alen * sizeof(int));	\
						Memcpy(val_new + kk, h->val + g_first_i, g_alen * sizeof(double)); \
						h->g_len[gg] = g_len;	\
						h->g_i[gg] = kk;	\
						kk += g_alen;		\
									\
						if (debug) {		\
							printf("\tmake a new group, gg=%1d with len=%1d\n", gg, g_len);	\
						}			\
									\
						gg_len = h->g_len[gg];	\
						gg_alen = IABS(gg_len);	\
						gg_first_i = h->g_i[gg]; \
						gg_last_i = gg_first_i + gg_alen - 1; \
						gg_last_idx = idx_new[gg_last_i]; \
						pending = (g < h->g_n - 1 ? 0 : 1); \
					}				\
					if (debug) {			\
						printf("gg=%1d kk=%1d\n", gg, kk); \
					}				\
				}					\
									\
				h->n_n = n_new;				\
				h->g_n = gg + pending;			\
				h->g_idx = idx_new;			\
				h->g_val = val_new;			\
				h->free_g_mem = 1;			\
			} else {					\
				h->n_n = h->n;				\
				h->g_idx = h->idx;			\
				h->g_val = h->val;			\
				h->free_g_mem = 0;			\
			}						\
		}							\
		if (debug) tref[3] += GMRFLib_cpu();			\
									\
		/*							\
		 * add a boolean for all(val[]==1), which makes dot-products into sums \
		 */							\
		if (debug) tref[4] -= GMRFLib_cpu();			\
		int *g_1 = Calloc(h->g_n + 1, int);			\
		for (int g = 0; g < h->g_n; g++) {			\
			int all_1 = 1;					\
			for (int j = 0; j < IABS(h->g_len[g]) && all_1; j++) { \
				all_1 = all_1 && ISEQUAL(h->g_val[h->g_i[g] + j], 1.0); \
			}						\
			g_1[g] = (all_1 ? 1 : 0);			\
		}							\
		h->g_1 = g_1;						\
		if (debug) tref[4] += GMRFLib_cpu();			\
									\
		if (debug) {						\
			GMRFLib_idxval_info_printf(stdout, h, "\t");	\
		}							\
		if (debug) {						\
			double sum = tref[0] + tref[1] + tref[2] + tref[3] + tref[4]; \
			printf("%.3f %.3f %.3f %.3f %.3f\n", tref[0]/sum, tref[1]/sum, tref[2]/sum, tref[3]/sum, tref[4]/sum); \
		}							\
	}

	RUN_CODE_BLOCK(nt, 1, nmax);
#undef CODE_BLOCK

	/*
	 * Add a tag about which is faster, the group or the serial algorithm, for each 'idxval'. I'm a little reluctant about doing this within
	 * the parallel loop. Usually, this is a very quick procedure, so it does not really matter...
	 */

	nmax = 1;
	for (int i = 0; i < n; i++) {
		GMRFLib_idxval_tp *h = hold[i];
		if (h->n) {
			nmax = IMAX(nmax, h->idx[h->n - 1]);
		}
	}
	nmax++;

	static double **wwork = NULL;
	static int *wwork_len = NULL;
	if (!wwork) {
#pragma omp critical (Name_5904d3e1eebd4db435c9b26e0854ee01328b78b2)
		{
			if (!wwork) {
				wwork_len = Calloc(GMRFLib_CACHE_LEN, int);
				wwork = Calloc(GMRFLib_CACHE_LEN, double *);
			}
		}
	}

	int cache_idx = 0;
	GMRFLib_CACHE_SET_ID(cache_idx);

	if (nmax > wwork_len[cache_idx]) {
		Free(wwork[cache_idx]);
		wwork_len[cache_idx] = nmax + 1024L;
		wwork[cache_idx] = Calloc(wwork_len[cache_idx], double);
		double *w = wwork[cache_idx];
		for (int j = 0; j < wwork_len[cache_idx]; j++) {
			w[j] = GMRFLib_uniform();
		}
	}
	double *x = wwork[cache_idx];

	double time_min = 0.0;
	double time_max = 0.0;
	int ntimes = 2;

	for (int i = 0; i < n; i++) {

		if (hold[i]->preference != IDXVAL_UNKNOWN) {
			continue;
		}

		double tref[4] = { 0.0, 0.0, 0.0, 0.0 };
		double value[4] = { 0.0, 0.0, 0.0, 0.0 };

		if (debug) {
			printf("start testing for hold[%1d]...\n", i);
		}
		for (int time = -1; time < ntimes; time++) {
			int measure = (time >= 0);
			if (measure) {
				tref[0] -= GMRFLib_cpu();
			}
			value[0] = GMRFLib_dot_product_serial(hold[i], x);
			if (measure) {
				tref[0] += GMRFLib_cpu();
			}
#if defined(INLA_LINK_WITH_MKL)
			if (measure) {
				tref[1] -= GMRFLib_cpu();
			}
			value[1] = GMRFLib_dot_product_serial_mkl(hold[i], x);
			if (measure) {
				tref[1] += GMRFLib_cpu();
			}
#else
			value[1] = value[0];
			tref[1] = tref[0];
#endif

			if (measure) {
				tref[2] -= GMRFLib_cpu();
			}
			value[2] = GMRFLib_dot_product_group(hold[i], x);
			if (measure) {
				tref[2] += GMRFLib_cpu();
			}
#if defined(INLA_LINK_WITH_MKL)
			if (measure) {
				tref[3] -= GMRFLib_cpu();
			}
			value[3] = GMRFLib_dot_product_group_mkl(hold[i], x);
			if (measure) {
				tref[3] += GMRFLib_cpu();
			}
#else
			value[3] = value[2];
			tref[3] = tref[2];
#endif
		}

		for (int k = 1; k < 4; k++) {
			// without cost, we can validate that the results is the same...
			if ((ABS(value[k] - value[0]) / (1.0 + (ABS(value[0]) + ABS(value[k])) / 2.0)) > 1.0E-6) {
				P(ABS(value[k] - value[0]) / (1.0 + (ABS(value[0]) + ABS(value[k])) / 2.0));
				P(k);
				P(value[0]);
				P(value[k]);
				assert(0 == 1);
			}
		}

		int k = -1;
		double tmin = GMRFLib_min_value(tref, 4, &k);
		double tmax = GMRFLib_max_value(tref, 4, NULL);

		if (debug) {
			double s = 1.0 / (tref[0] + tref[1] + tref[2] + tref[3]) / ntimes;
			printf("for h[%1d] with n= %1d chose k=%1d [serial= %.3f serial.mkl= %.3f group= %.3f group.mkl= %.3f]\n",
			       i, hold[i]->n, k, tref[0] * s, tref[1] * s, tref[2] * s, tref[3] * s);
		}

		switch (k) {
		case 0:
			hold[i]->preference = IDXVAL_SERIAL;
			break;
		case 1:
			hold[i]->preference = IDXVAL_SERIAL_MKL;
			break;
		case 2:
			hold[i]->preference = IDXVAL_GROUP;
			break;
		case 3:
			hold[i]->preference = IDXVAL_GROUP_MKL;
			break;
		default:
			assert(0 == 1);
		}

		if (IDXVAL_FREE_GROUPS_IF_NOT_BEST()) {
			if (hold[i]->preference == IDXVAL_SERIAL || hold[i]->preference == IDXVAL_SERIAL_MKL) {
				// no need to keep the group info in the struct
				hold[i]->n_n = 0;
				hold[i]->g_n = 0;
				hold[i]->g_idx = NULL;
				hold[i]->g_val = NULL;
				Free(hold[i]->g_len);
				Free(hold[i]->g_i);
				Free(hold[i]->g_1);
				if (hold[i]->free_g_mem) {
					Free(hold[i]->g_idx);
					Free(hold[i]->g_val);
				}
			}
		}

		if (GMRFLib_dot_product_optim_report) {
			int idx;
			GMRFLib_CACHE_SET_ID(idx);
			for (int k = 0; k < 4; k++) {
				GMRFLib_dot_product_optim_report[idx][k] += tref[k];
			}
			GMRFLib_dot_product_optim_report[idx][4] += tmin;
			GMRFLib_dot_product_optim_report[idx][5 + k]++;	/* count... */
		}

		time_min += tmin / ntimes;
		time_max += tmax / ntimes;
	}

	if (debug) {
		if (time_min > 0.0) {
			printf("idxval opt: saving %.6f seconds/M.evals, %.2f%% improvement\n",
			       (time_max - time_min) * ISQR(1024), 100.0 * (1.0 - time_min / time_max));
		}
	}

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_idx_free(GMRFLib_idx_tp * hold)
{
	if (hold) {
		Free(hold->idx);
		Free(hold);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_str_free(GMRFLib_str_tp * hold)
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

int GMRFLib_idx2_free(GMRFLib_idx2_tp * hold)
{
	if (hold) {
		Free(hold->idx[0]);
		Free(hold->idx[1]);
		Free(hold->idx);
		Free(hold);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_val_free(GMRFLib_val_tp * hold)
{
	if (hold) {
		Free(hold->val);
		Free(hold);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idxval_free(GMRFLib_idxval_tp * hold)
{
	if (hold) {
		Free(hold->idx);
		Free(hold->val);

		Free(hold->g_i);
		Free(hold->g_len);
		Free(hold->g_1);
		if (hold->free_g_mem) {
			Free(hold->g_idx);
			Free(hold->g_val);
		}
		Free(hold);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idx_add(GMRFLib_idx_tp ** hold, int idx)
{
	if (*hold == NULL) {
		GMRFLib_idx_create(hold);
	}
	if ((*hold)->n == (*hold)->n_alloc) {
		(*hold)->n_alloc += IMAX(IDX_ALLOC_ADD, (*hold)->n / 4);
		(*hold)->idx = Realloc((*hold)->idx, (*hold)->n_alloc, int);
	}
	(*hold)->idx[(*hold)->n] = idx;
	(*hold)->n++;

	return GMRFLib_SUCCESS;
}

int GMRFLib_str_add(GMRFLib_str_tp ** hold, char *s)
{
	if (*hold == NULL) {
		GMRFLib_str_create(hold);
	}
	if ((*hold)->n == (*hold)->n_alloc) {
		(*hold)->n_alloc += IMAX(IDX_ALLOC_ADD, (*hold)->n / 4);
		(*hold)->str = Realloc((*hold)->str, (*hold)->n_alloc, char *);
	}
	(*hold)->str[(*hold)->n] = GMRFLib_strdup(s);
	(*hold)->n++;

	return GMRFLib_SUCCESS;
}

int GMRFLib_idx2_add(GMRFLib_idx2_tp ** hold, int idx0, int idx1)
{
	if (*hold == NULL) {
		GMRFLib_idx2_create(hold);
	}
	if ((*hold)->n == (*hold)->n_alloc) {
		(*hold)->n_alloc += IMAX(IDX_ALLOC_ADD, (*hold)->n / 4);
		(*hold)->idx[0] = Realloc((*hold)->idx[0], (*hold)->n_alloc, int);
		(*hold)->idx[1] = Realloc((*hold)->idx[1], (*hold)->n_alloc, int);
	}
	(*hold)->idx[0][(*hold)->n] = idx0;
	(*hold)->idx[1][(*hold)->n] = idx1;
	(*hold)->n++;

	return GMRFLib_SUCCESS;
}

int GMRFLib_val_add(GMRFLib_val_tp ** hold, double val)
{
	if (*hold == NULL) {
		GMRFLib_val_create(hold);
	}
	if ((*hold)->n == (*hold)->n_alloc) {
		(*hold)->n_alloc += IMAX(IDX_ALLOC_ADD, (*hold)->n / 4);
		(*hold)->val = Realloc((*hold)->val, (*hold)->n_alloc, double);
	}
	(*hold)->val[(*hold)->n] = val;
	(*hold)->n++;

	return GMRFLib_SUCCESS;
}

int GMRFLib_idxval_add(GMRFLib_idxval_tp ** hold, int idx, double val)
{
	if (*hold == NULL) {
		GMRFLib_idxval_create(hold);
	}
	if ((*hold)->n == (*hold)->n_alloc) {
		(*hold)->n_alloc += IMAX(IDX_ALLOC_ADD, (*hold)->n / 4);
		(*hold)->idx = Realloc((*hold)->idx, (*hold)->n_alloc, int);
		(*hold)->val = Realloc((*hold)->val, (*hold)->n_alloc, double);
	}
	(*hold)->idx[(*hold)->n] = idx;
	(*hold)->val[(*hold)->n] = val;
	(*hold)->n++;

	return GMRFLib_SUCCESS;
}

int GMRFLib_idxval_addto(GMRFLib_idxval_tp ** hold, int idx, double val)
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
		(*hold)->n_alloc += IMAX(IDX_ALLOC_ADD, (*hold)->n / 4);
		(*hold)->idx = Realloc((*hold)->idx, (*hold)->n_alloc, int);
		(*hold)->val = Realloc((*hold)->val, (*hold)->n_alloc, double);
	}
	(*hold)->idx[(*hold)->n] = idx;
	(*hold)->val[(*hold)->n] = val;
	(*hold)->n++;

	return GMRFLib_SUCCESS;
}


int GMRFLib_str_is_member(GMRFLib_str_tp * hold, char *s, int case_sensitive, int *idx_match)
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
