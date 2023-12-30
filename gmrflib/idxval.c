
/* idxval.c
 * 
 * Copyright (C) 2022-2023 Havard Rue
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

#include <assert.h>
#include <float.h>
#include <signal.h>
#include <stdarg.h>
#include <math.h>
#include <strings.h>
#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"
#include "GMRFLib/hashP.h"

#define IDX_ALLOC_INITIAL 32
#define IDX_ALLOC_ADD     512
#define IDX_ALLOC_SUBMATRIX_INITIAL 8
#define IDX_ALLOC_SUBMATRIX_ADD 32

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

int GMRFLib_idxval_printf(FILE *fp, GMRFLib_idxval_tp *hold, const char *msg)
{
	if (hold) {
		fprintf(fp, "[%s] n = %1d  nalloc = %1d iaddto = %1d\n", msg, hold->n, hold->n_alloc, hold->iaddto);
		for (int i = 0; i < hold->n; i++) {
			fprintf(fp, "\t(idx, val)[%1d] = (%d, %g)\n", i, hold->idx[i], hold->val[i]);
		}
		if (hold->g_n) {
			fprintf(fp, "\tg_n = %1d\n", hold->g_n);
			for (int g = 0; g < hold->g_n; g++) {
				fprintf(fp, "\tgroup %d has length %d (one=%s)\n", g, hold->g_len[g],
					(hold->g_1 && hold->g_1[g] ? "TRUE" : "FALSE"));
				fprintf(fp, "\t\t");
				for (int k = 0; k < IABS(hold->g_len[g]); k++) {
					fprintf(fp, " %1d(%g)", hold->g_idx[g][k], hold->g_val[g][k]);

				}
				fprintf(fp, "\n");
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
		qsort((void *) hold->idx, (size_t) hold->n, sizeof(int), GMRFLib_icmp);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idx_nsort(GMRFLib_idx_tp **a, int n, int nt)
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

int GMRFLib_idxval_nprune(GMRFLib_idxval_tp **a, int n, int nt)
{
#define CODE_BLOCK					\
	for (int i = 0; i < n; i++) {			\
		GMRFLib_idxval_prune(a[i]);		\
	}

	RUN_CODE_BLOCK(nt, 0, 0);
#undef CODE_BLOCK

	return GMRFLib_SUCCESS;
}

int GMRFLib_idx_prune(GMRFLib_idx_tp *hold)
{
	if (hold) {
		hold->idx = Realloc(hold->idx, IMAX(1, hold->n), int);
		hold->n_alloc = IMAX(1, hold->n);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_str_prune(GMRFLib_str_tp *hold)
{
	if (hold) {
		hold->str = Realloc(hold->str, IMAX(1, hold->n), char *);
		hold->n_alloc = IMAX(1, hold->n);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idx2_prune(GMRFLib_idx2_tp *hold)
{
	if (hold) {
		hold->idx[0] = Realloc(hold->idx[0], IMAX(1, hold->n), int);
		hold->idx[1] = Realloc(hold->idx[1], IMAX(1, hold->n), int);
		hold->n_alloc = IMAX(1, hold->n);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_val_prune(GMRFLib_val_tp *hold)
{
	if (hold) {
		hold->val = Realloc(hold->val, IMAX(1, hold->n), double);
		hold->n_alloc = IMAX(1, hold->n);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idxval_prune(GMRFLib_idxval_tp *hold)
{
	if (hold) {
		int n = IMAX(1, hold->n);
		hold->idx = Realloc(hold->idx, n, int);
		hold->val = Realloc(hold->val, n, double);
		hold->n_alloc = n;
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

int GMRFLib_idxval_nsort_x_core(GMRFLib_idxval_tp *h, double *x, int prepare, int accumulate)
{
	const int limit = 16L;
	const int debug = 0;

	if (!h || h->n == 0) {
		return GMRFLib_SUCCESS;
	}
	// sort
	if (h->n > 1) {
		int is_sorted = 1;
		for (int j = 1; is_sorted && j < h->n; j++) {
			is_sorted = (h->idx[j] >= h->idx[j - 1]);
		}
		if (!is_sorted) {
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

	if (h->n <= limit || !prepare || !GMRFLib_internal_opt) {
		h->preference = IDXVAL_GROUP_MKL;
		return GMRFLib_SUCCESS;
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
	g_istart[ng] = h->n;
	g_len[ng] = 0;

	if (debug) {
		for (i = 0; i < h->n; i++) {
			printf("idx[%1d] =  %1d\n", i, h->idx[i]);
		}
		printf("ng = %1d\n", ng);
		for (int g = 1; g < ng; g++) {
			printf("group %1d start %1d len %1d\n", g, g_istart[g], g_len[g]);
			for (i = 0; i < g_len[g]; i++) {
				printf("\t\t\tidx %1d\n", h->idx[g_istart[g] + i]);
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

	int *new_idx = Calloc(irr_len + seq_len + ng * limit, int);
	double *new_val = Calloc(irr_len + seq_len + ng * limit, double);

	// build the irregular group
	int k = 0;
	for (int g = 1; g < ng + 1; g++) {
		int istart = g_istart[g - 1] + g_len[g - 1];
		int len = g_istart[g] - istart;
		if (len) {
			Memcpy(new_idx + k, h->idx + istart, len * sizeof(int));
			Memcpy(new_val + k, h->val + istart, len * sizeof(double));
		}
		k += len;
	}
	int new_len = k;
	g_len[0] = new_len;
	g_idx[0] = new_idx;
	g_val[0] = new_val;

	if (debug) {
		for (i = 0; i < g_len[0]; i++) {
			printf("nidx[%1d] %1d  val %g\n", i, new_idx[i], new_val[i]);
		}
	}
	// copy each sequential group and pad for possible grouping
	int *seq_idx = new_idx + irr_len;
	double *seq_val = new_val + irr_len;

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
			printf("group %d istart %d len %d limit %d\n", i, istart, len,
			       (i < ng - 1 ? IMIN(limit, h->idx[g_istart[i + 1]] - h->idx[istart + len - 1] - 1) : 0));
		}

		int pad = (i < ng - 1 ? IMIN(limit, h->idx[g_istart[i + 1]] - h->idx[istart + len - 1] - 1) : 0);
		k += len;
		int offset = seq_idx[k - 1] + 1;
		for (int j = 0; j < pad; j++) {
			seq_idx[k + j] = offset + j;
		}
		k += pad;
	}

	// setup pointers to each sequential group
	for (int g = 1; g < ng; g++) {
		int istart = g_istart[g];
		g_idx[g] = seq_idx + istart;
		g_val[g] = seq_val + istart;
		g_len[g] *= -1;
	}

	if (debug) {
		for (i = 0; i < k; i++) {
			printf("i %d new_idx %1d new_val %f\n", i, new_idx[i], new_val[i]);
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
			for (i = 0; i < IABS(g_len[g]); i++) {
				printf("\t\t\tidx %1d val %g\n", g_idx[g][i], g_val[g][i]);
			}
		}
	}
	// remove groups with len=0
	int g = 0;
	while (1) {
		if (IABS(g_len[g] == 0)) {
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
			if (g2_start - g1_end <= limit && (g_1[g] == 0 && g_1[g + 1] == 0)) {
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
	Free(g_istart);

	int ntimes = 2;
#if defined(INLA_LINK_WITH_MKL)
	int with_mkl = 1;
#else
	int with_mkl = 0;
#endif

	if (h->preference != IDXVAL_UNKNOWN) {
		return GMRFLib_SUCCESS;
	}

	double treff[4] = { 0.0, 0.0, 0.0, 0.0 };
	double value[4] = { 0.0, 0.0, 0.0, 0.0 };

	for (int time = -1; time < ntimes; time++) {
		if (time < 0) {
			GMRFLib_dot_product_serial(h, x);
			GMRFLib_dot_product_group(h, x);
			if (with_mkl) {
				GMRFLib_dot_product_serial_mkl(h, x);
				GMRFLib_dot_product_group_mkl(h, x);
			}
		} else {
			treff[0] -= GMRFLib_cpu();
			value[0] = GMRFLib_dot_product_serial(h, x);
			treff[0] += GMRFLib_cpu();

			if (with_mkl) {
				treff[1] -= GMRFLib_cpu();
				value[1] = GMRFLib_dot_product_serial_mkl(h, x);
				treff[1] += GMRFLib_cpu();
			} else {
				value[1] = value[0];
				treff[1] = treff[0];
			}

			treff[2] -= GMRFLib_cpu();
			value[2] = GMRFLib_dot_product_group(h, x);
			treff[2] += GMRFLib_cpu();

			if (with_mkl) {
				treff[3] -= GMRFLib_cpu();
				value[3] = GMRFLib_dot_product_group_mkl(h, x);
				treff[3] += GMRFLib_cpu();
			} else {
				value[3] = value[2];
				treff[3] = treff[2];
			}

			if (0) {
				printf("idxval optimisation: length = %1d\n", h->n);
				printf("\tserial     value   = %.16g\n", value[0]);
				printf("\tserial_mkl abs.err = %.16g\n", ABS(value[1] - value[0]));
				printf("\tgroup      abs.err = %.16g\n", ABS(value[2] - value[0]));
				printf("\tgroup_mkl  abs.err = %.16g\n", ABS(value[3] - value[0]));
			}
		}
	}

	for (k = 0; k < 4; k++) {
		treff[k] /= (double) ntimes;
	}

	for (k = 1; k < 4; k++) {
		if (ABS(value[k] - value[0]) > 1000.0 * FLT_EPSILON * sqrt(h->n)) {
			P(ABS(value[k] - value[0]));
			P(k);
			P(value[0]);
			P(value[1]);
			P(value[2]);
			P(value[3]);

			printf("n %d\n", h->n);
			for (i = 0; i < h->n; i++) {
				printf("\tidx[%1d] =  %1d  val = %g\n", i, h->idx[i], h->val[i]);
			}
			printf("ng %d\n", h->g_n);
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
			P(value[2]);
			P(value[3]);

			assert(0 == 1);
		}
	}

	int kmin = -1;
	double tmin = GMRFLib_min_value(treff, 4, &kmin);

	if (debug) {
		double s = 1.0 / (DBL_EPSILON + treff[0] + treff[1] + treff[2] + treff[3]);
		printf("for h with n= %1d chose kmin=%1d [serial= %.3f serial.mkl= %.3f group= %.3f group.mkl= %.3f]\n",
		       h->n, kmin, treff[0] * s, treff[1] * s, treff[2] * s, treff[3] * s);
	}

	switch (kmin) {
	case 0:
		h->preference = IDXVAL_SERIAL;
		break;
	case 1:
		h->preference = IDXVAL_SERIAL_MKL;
		break;
	case 2:
		h->preference = IDXVAL_GROUP;
		break;
	case 3:
		h->preference = IDXVAL_GROUP_MKL;
		break;
	default:
		assert(0 == 1);
	}
	h->cpu_gain = treff[1] - treff[kmin];
	assert(h->cpu_gain >= 0);

	if (h->preference == IDXVAL_SERIAL || h->preference == IDXVAL_SERIAL_MKL) {
		/*
		 * no need to keep the group info in the struct 
		 */
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
	}

	if (GMRFLib_dot_product_optim_report) {
		int idx = 0;
		GMRFLib_CACHE_SET_ID(idx);
		for (k = 0; k < 4; k++) {
			GMRFLib_dot_product_optim_report[idx][k] += treff[k];
		}
		GMRFLib_dot_product_optim_report[idx][4] += tmin;
		GMRFLib_dot_product_optim_report[idx][5 + kmin]++;	/* count... */
	}

	return GMRFLib_SUCCESS;
}

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
		GMRFLib_idxval_nsort_x_core(hold[k], x_ran, prepare, accumulate); \
	}

	RUN_CODE_BLOCK_DYNAMIC(nt, 0, 0);
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

int GMRFLib_idxval_free(GMRFLib_idxval_tp *hold)
{
	if (hold) {
		Free(hold->idx);
		Free(hold->val);

		Free(hold->g_len);
		Free(hold->g_1);
		for (int i = 0; i < hold->g_n_mem; i++) {
			Free(hold->g_mem[i]);
		}
		if (hold->g_mem) {
			Free(hold->g_mem);
		}
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
		(*hold)->n_alloc += IMAX(IDX_ALLOC_ADD, (*hold)->n / 4);
		(*hold)->idx = Realloc((*hold)->idx, (*hold)->n_alloc, int);
	}
	(*hold)->idx[(*hold)->n] = idx;
	(*hold)->n++;

	return GMRFLib_SUCCESS;
}

int GMRFLib_idx_nadd(GMRFLib_idx_tp **hold, int n, int *idx)
{
	if (*hold == NULL) {
		GMRFLib_idx_create(hold);
	}
	if ((*hold)->n + n > (*hold)->n_alloc) {
		(*hold)->n_alloc += IMAX(IDX_ALLOC_ADD, n + (*hold)->n / 4);
		(*hold)->idx = Realloc((*hold)->idx, (*hold)->n_alloc, int);
	}

	Memcpy((*hold)->idx + (*hold)->n, idx, n * sizeof(int));
	(*hold)->n += n;

	return GMRFLib_SUCCESS;
}

GMRFLib_idx_tp *GMRFLib_idx_duplicate(GMRFLib_idx_tp *h)
{
	GMRFLib_idx_tp *new = NULL;
	GMRFLib_idx_create_x(&new, (h ? IMAX(1, h->n) : 1));
	if (h && h->n > 0) {
		GMRFLib_idx_nadd(&new, h->n, h->idx);
	}

	return new;
}

int GMRFLib_str_add(GMRFLib_str_tp **hold, char *s)
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



int GMRFLib_str_is_member(GMRFLib_str_tp *hold, char *s, int case_sensitive, int *idx_match)
{
	if (hold == NULL) {
		return 0;
	}

	int (*cmp)(const char *, const char *) = (case_sensitive ? strcmp : strcasecmp);
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
		(*hold)->n_alloc += IMAX(IDX_ALLOC_ADD, (*hold)->n / 4);
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
		(*hold)->n_alloc += IMAX(IDX_ALLOC_ADD, (*hold)->n / 4);
		(*hold)->val = Realloc((*hold)->val, (*hold)->n_alloc, double);
	}
	(*hold)->val[(*hold)->n] = val;
	(*hold)->n++;

	return GMRFLib_SUCCESS;
}

int GMRFLib_idxval_add(GMRFLib_idxval_tp **hold, int idx, double val)
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
		(*hold)->n_alloc += IMAX(IDX_ALLOC_ADD, (*hold)->n / 4);
		(*hold)->idx = Realloc((*hold)->idx, (*hold)->n_alloc, int);
		(*hold)->val = Realloc((*hold)->val, (*hold)->n_alloc, double);
	}
	(*hold)->idx[(*hold)->n] = idx;
	(*hold)->val[(*hold)->n] = val;
	(*hold)->n++;

	return GMRFLib_SUCCESS;
}
