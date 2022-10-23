
/* dot-product.c
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

#include <stdint.h>
#include <assert.h>
#include <stddef.h>
#include <time.h>
#include <strings.h>
#include <string.h>
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <omp.h>

#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_math.h>
#include "sparse-vec.h"

#define ISZERO(x) (gsl_fcmp(1.0 + (x), 1.0, DBL_EPSILON) == 0)
#define ISEQUAL(x, y) (gsl_fcmp(x, y, DBL_EPSILON) == 0)

#define GMRFLib_ALLOC_SAFE_SIZE(n_, type_) ((size_t)(n_) * sizeof(type_) < PTRDIFF_MAX ? (size_t)(n_) : (size_t)1)
#define Calloc(n, type)         (type *)calloc(GMRFLib_ALLOC_SAFE_SIZE(n, type), sizeof(type))
#define Malloc(n, type)         (type *)malloc(GMRFLib_ALLOC_SAFE_SIZE((n) * sizeof(type), char))
#define Realloc(ptr, n, type)   ((ptr) ? \
                                 (type *)realloc((void *)ptr, GMRFLib_ALLOC_SAFE_SIZE((n) * sizeof(type), char)) : \
                                 (type *)calloc(GMRFLib_ALLOC_SAFE_SIZE(n, type), sizeof(type)))
#define Free(ptr)               if (ptr) {free((void *)(ptr)); ptr=NULL;}
#define Memcpy(dest, src, n)    memcpy((void *) (dest), (void *) (src), GMRFLib_ALLOC_SAFE_SIZE(n, char))
#define Memset(dest, value, n)  memset((void *) (dest), (int) (value), (size_t) (n))

#define ABS(x)   fabs(x)
#define SQR(x)   ((x)*(x))
#define IABS(x)   abs(x)
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define P(x)     if (1) { printf("[%s:%1d] " #x " = [ %.12f ]\n",__FILE__, __LINE__,(double)(x)); }
#define UNIFORM() (rand() / ((double) RAND_MAX + 1.0))
#define SEED(seed_) srand((unsigned int) (seed_))

#define SIMPLE_LOOP_LIMIT 8L
#define IDX_ALLOC_INITIAL 32
#define IDX_ALLOC_ADD     512

#if !defined(SHOW_TEST_OUTPUT)
#define SHOW_TEST_OUTPUT 0
#endif

void mkl_dcsrmv(const char *transa, const int *m, const int *k, const double *alpha,
		const char *matdescra, const double *val, const int *indx,
		const int *pntrb, const int *pntre, const double *x, const double *beta, double *y);
double ddot_(int *len, double *x, int *incx, double *y, int *incy);

#define GMRFLib_cpu() omp_get_wtime()
#define GMRFLib_uniform() ((double) rand() / (1.0 + RAND_MAX))

double GMRFLib_dot_product_group(GMRFLib_idxval_tp * __restrict ELM_, double *__restrict ARR_)
{
	// this uses g_idx and g_val

	double value_ = 0.0;
	for (int g_ = 0; g_ < ELM_->g_n; g_++) {
		int len_ = ELM_->g_len[g_];
		if (len_ == 0)
			continue;

		int *__restrict ii_ = ELM_->g_idx[g_];
		double *__restrict vv_ = ELM_->g_val[g_];

		if (len_ > 0) {
			double *__restrict aa_ = &(ARR_[0]);
			if (ELM_->g_1[g_]) {
				if (len_ < SIMPLE_LOOP_LIMIT) {
#pragma GCC ivdep
					for (int i_ = 0; i_ < len_; i_++) {
						value_ += aa_[ii_[i_]];
					}
				} else {
					value_ += GMRFLib_dsum_idx(len_, aa_, ii_);
				}
			} else {
				if (len_ < SIMPLE_LOOP_LIMIT) {
#pragma GCC ivdep
					for (int i_ = 0; i_ < len_; i_++) {
						value_ += vv_[i_] * aa_[ii_[i_]];
					}
				} else {
					value_ += GMRFLib_ddot_idx(len_, vv_, aa_, ii_);
				}
			}
		} else if (len_ < 0) {
			int llen_ = -len_;
			double *__restrict aa_ = &(ARR_[ii_[0]]);
			if (ELM_->g_1[g_]) {
				if (llen_ < SIMPLE_LOOP_LIMIT) {
#pragma GCC ivdep
					for (int i_ = 0; i_ < llen_; i_++) {
						value_ += aa_[i_];
					}
				} else {
					value_ += GMRFLib_dsum(llen_, aa_);
				}
			} else {
				if (llen_ < SIMPLE_LOOP_LIMIT) {
#pragma GCC ivdep
					for (int i_ = 0; i_ < llen_; i_++) {
						value_ += vv_[i_] * aa_[i_];
					}
				} else {
					value_ += GMRFLib_ddot(llen_, vv_, aa_);
				}
			}
		}
	}
	return (value_);
}

double GMRFLib_dot_product_group_mkl(GMRFLib_idxval_tp * __restrict ELM_, double *__restrict ARR_)
{
	// this uses n_idx and n_val

	double value_ = 0.0;
	for (int g_ = 0; g_ < ELM_->g_n; g_++) {
		int len_ = ELM_->g_len[g_];
		if (len_ == 0)
			continue;

		int *__restrict ii_ = ELM_->g_idx[g_];
		double *__restrict vv_ = ELM_->g_val[g_];

		if (len_ > 0) {
			double *__restrict aa_ = &(ARR_[0]);
			if (ELM_->g_1[g_]) {
				if (len_ < SIMPLE_LOOP_LIMIT) {
#pragma GCC ivdep
					for (int i_ = 0; i_ < len_; i_++) {
						value_ += aa_[ii_[i_]];
					}
				} else {
					value_ += GMRFLib_dsum_idx(len_, aa_, ii_);
				}
			} else {
				if (len_ < SIMPLE_LOOP_LIMIT) {
#pragma GCC ivdep
					for (int i_ = 0; i_ < len_; i_++) {
						value_ += vv_[i_] * aa_[ii_[i_]];
					}
				} else {
					value_ += GMRFLib_ddot_idx_mkl(len_, vv_, aa_, ii_);
				}
			}
		} else if (len_ < 0) {
			int llen_ = -len_;
			double *__restrict aa_ = &(ARR_[ii_[0]]);
			if (ELM_->g_1[g_]) {
				if (llen_ < SIMPLE_LOOP_LIMIT) {
#pragma GCC ivdep
					for (int i_ = 0; i_ < llen_; i_++) {
						value_ += aa_[i_];
					}
				} else {
					value_ += GMRFLib_dsum(llen_, aa_);
				}
			} else {
				if (llen_ < SIMPLE_LOOP_LIMIT) {
#pragma GCC ivdep
					for (int i_ = 0; i_ < llen_; i_++) {
						value_ += vv_[i_] * aa_[i_];
					}
				} else {
					value_ += GMRFLib_ddot(llen_, vv_, aa_);
				}
			}
		}
	}
	return (value_);
}

double GMRFLib_dot_product_serial(GMRFLib_idxval_tp * __restrict ELM_, double *__restrict ARR_)
{
	double value_ = 0.0;
	double *__restrict vv_ = ELM_->val;
	double *__restrict aa_ = ARR_;
	int *__restrict idx_ = ELM_->idx;

	if (ELM_->n < SIMPLE_LOOP_LIMIT) {
#pragma GCC ivdep
		for (int i_ = 0; i_ < ELM_->n; i_++) {
			value_ += vv_[i_] * aa_[idx_[i_]];
		}
	} else {
		value_ += GMRFLib_ddot_idx(ELM_->n, vv_, aa_, idx_);
	}

	return (value_);
}

double GMRFLib_dot_product_serial_mkl(GMRFLib_idxval_tp * __restrict ELM_, double *__restrict ARR_)
{
	double value_ = 0.0;
	double *__restrict vv_ = ELM_->val;
	double *__restrict aa_ = ARR_;
	int *__restrict idx_ = ELM_->idx;

	if (ELM_->n < SIMPLE_LOOP_LIMIT) {
#pragma GCC ivdep
		for (int i_ = 0; i_ < ELM_->n; i_++) {
			value_ += vv_[i_] * aa_[idx_[i_]];
		}
	} else {
		value_ += GMRFLib_ddot_idx_mkl(ELM_->n, vv_, aa_, idx_);
	}

	return (value_);
}

double GMRFLib_dot_product(GMRFLib_idxval_tp * __restrict ELM_, double *__restrict ARR_)
{
	switch (ELM_->preference) {
	case IDXVAL_SERIAL:
		return (GMRFLib_dot_product_serial(ELM_, ARR_));
		break;
	case IDXVAL_SERIAL_MKL:
		return (GMRFLib_dot_product_serial_mkl(ELM_, ARR_));
		break;
	case IDXVAL_GROUP:
		return (GMRFLib_dot_product_group(ELM_, ARR_));
		break;
	case IDXVAL_GROUP_MKL:
		return (GMRFLib_dot_product_group_mkl(ELM_, ARR_));
		break;
	case IDXVAL_UNKNOWN:
		return (GMRFLib_dot_product_serial_mkl(ELM_, ARR_));
		break;
	default:
		assert(0 == 1);
		break;
	}
	assert(0 == 1);

	return NAN;
}


double GMRFLib_dsum(int n, double *x)
{
	double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
	const int roll = 8L;
	div_t d = div(n, roll);
	int m = d.quot * roll;

#pragma GCC ivdep
	for (int i = 0; i < m; i += roll) {
		double *xx = x + i;

		s0 += xx[0];
		s1 += xx[1];
		s2 += xx[2];
		s3 += xx[3];

		s0 += xx[4];
		s1 += xx[5];
		s2 += xx[6];
		s3 += xx[7];
	}

#pragma GCC ivdep
	for (int i = d.quot * roll; i < n; i++) {
		s0 += x[i];
	}

	return (s0 + s1 + s2 + s3);
}

double GMRFLib_ddot(int n, double *__restrict x, double *__restrict y)
{
	int one = 1;
	return ddot_(&n, x, &one, y, &one);
}

double GMRFLib_dsum_idx(int n, double *__restrict a, int *__restrict idx)
{
	const int roll = 8L;
	double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
	div_t d = div(n, roll);

#pragma GCC ivdep
	for (int i = 0; i < d.quot * roll; i += roll) {
		int *iidx = idx + i;

		s0 += a[iidx[0]];
		s1 += a[iidx[1]];
		s2 += a[iidx[2]];
		s3 += a[iidx[3]];

		s0 += a[iidx[4]];
		s1 += a[iidx[5]];
		s2 += a[iidx[6]];
		s3 += a[iidx[7]];
	}

#pragma GCC ivdep
	for (int i = d.quot * roll; i < n; i++) {
		s0 += a[idx[i]];
	}

	return (s0 + s1 + s2 + s3);
}

double GMRFLib_ddot_idx(int n, double *__restrict v, double *__restrict a, int *__restrict idx)
{
	const int roll = 8L;
	double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
	div_t d = div(n, roll);
	int m = d.quot * roll;

#pragma GCC ivdep
	for (int i = 0; i < m; i += roll) {
		double *vv = v + i;
		int *iidx = idx + i;

		s0 += vv[0] * a[iidx[0]];
		s1 += vv[1] * a[iidx[1]];
		s2 += vv[2] * a[iidx[2]];
		s3 += vv[3] * a[iidx[3]];

		s0 += vv[4] * a[iidx[4]];
		s1 += vv[5] * a[iidx[5]];
		s2 += vv[6] * a[iidx[6]];
		s3 += vv[7] * a[iidx[7]];
	}

#pragma GCC ivdep
	for (int i = d.quot * roll; i < n; i++) {
		s0 += v[i] * a[idx[i]];
	}

	return (s0 + s1 + s2 + s3);
}

#if defined(WITH_MKL)

double GMRFLib_ddot_idx_mkl(int n, double *__restrict v, double *__restrict a, int *__restrict idx)
{
	// this is the MKL version, which is done using a sparse '1 x n' matrix.
	// we could include <mkl.h> but we can just do this, as we only need one non-standard function

	int iarr[4] = { 1, 0, n, idx[n - 1] + 1 };
	double darr[3] = { 1.0, 0.0, 0.0 };
	// we need to define this with length 6. the fifth argument is not used, so we use it for the
	// argument 'trans', the first argument in the call, trans='N'
	const char matdescra[6] = { 'G', '.', '.', 'C', 'N', '.' };

	mkl_dcsrmv(matdescra + 4, iarr, iarr + 3, darr, matdescra, v, idx, iarr + 1, iarr + 2, a, darr + 1, darr + 2);
	return (darr[2]);
}

#else							       /* defined(INLA_LINK_WITH_MKL) */
double GMRFLib_ddot_idx_mkl(int n, double *__restrict v, double *__restrict a, int *__restrict idx)
{
	return GMRFLib_ddot_idx(n, v, a, idx);
}

#endif							       /* if defined(INLA_LINK_WITH_MKL) */

int GMRFLib_idxval_free(GMRFLib_idxval_tp * hold)
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
	return 0;
}

int GMRFLib_idxval_create(GMRFLib_idxval_tp ** hold)
{
	return GMRFLib_idxval_create_x(hold, IDX_ALLOC_INITIAL);
}

int GMRFLib_idxval_create_x(GMRFLib_idxval_tp ** hold, int len)
{
	len = MAX(1, len);
	*hold = Calloc(1, GMRFLib_idxval_tp);
	(*hold)->idx = Calloc(len, int);
	(*hold)->val = Calloc(len, double);
	(*hold)->n_alloc = len;
	(*hold)->n = 0;
	(*hold)->iaddto = 0;

	return 0;
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

int GMRFLib_idxval_printf(FILE * fp, GMRFLib_idxval_tp * hold, const char *msg)
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
	return 0;
}

int GMRFLib_idxval_info_printf(FILE * fp, GMRFLib_idxval_tp * hold, const char *msg)
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
	return 0;
}

int GMRFLib_idxval_nuniq(GMRFLib_idxval_tp ** a, int n, int nt)
{
	for (int i = 0; i < n; i++) {
		GMRFLib_idxval_uniq(a[i]);
	}
	return 0;
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
	return 0;
}

int GMRFLib_idxval_nprune(GMRFLib_idxval_tp ** a, int n, int nt)
{
	for (int i = 0; i < n; i++) {
		GMRFLib_idxval_prune(a[i]);
	}

	return 0;
}

int GMRFLib_idxval_prune(GMRFLib_idxval_tp * hold)
{
	if (hold) {
		int n = MAX(1, hold->n);
		hold->idx = Realloc(hold->idx, n, int);
		hold->val = Realloc(hold->val, n, double);
		hold->n_alloc = n;
	}
	return 0;
}

int GMRFLib_idxval_sort(GMRFLib_idxval_tp * hold)
{
	return GMRFLib_idxval_nsort(&hold, 1, 1);
}

int GMRFLib_idxval_nsort(GMRFLib_idxval_tp ** hold, int n, int nt)
{
	return GMRFLib_idxval_nsort_x(hold, n, nt);
}

int GMRFLib_idxval_prepare(GMRFLib_idxval_tp * hold)
{
	return (GMRFLib_idxval_nsort_x(&hold, 1, 1));
}

int GMRFLib_idxval_nsort_x_core(GMRFLib_idxval_tp * h, double *x)
{
	const int limit = 8L;
	const int debug = 0;

	if (!h || h->n == 0) {
		return 0;
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
					h->val[k] += h->val[j];
				}
			}
			if (debug && (h->n > k + 1)) {
				printf("Make unique: reduce length from %d to %dn", h->n, k + 1);
			}
			h->n = k + 1;
		}
	}

	if (h->n <= limit) {
		h->preference = IDXVAL_SERIAL_MKL;
		return 0;
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
		for (int i = 0; i < h->n; i++) {
			printf("idx[%1d] =  %1d\n", i, h->idx[i]);
		}
		printf("ng = %1d\n", ng);
		for (int g = 1; g < ng; g++) {
			printf("group %1d start %1d len %1d\n", g, g_istart[g], g_len[g]);
			for (int i = 0; i < g_len[g]; i++) {
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
		for (int i = 0; i < g_len[0]; i++) {
			printf("nidx[%1d] %1d  val %g\n", i, new_idx[i], new_val[i]);
		}
	}
	// copy each sequential group and pad for possible grouping
	int *seq_idx = new_idx + irr_len;
	double *seq_val = new_val + irr_len;

	k = 0;
	for (int i = 1; i < ng; i++) {
		int istart = g_istart[i];
		int len = g_len[i];
		if (len) {
			Memcpy(seq_idx + k, h->idx + istart, len * sizeof(int));
			Memcpy(seq_val + k, h->val + istart, len * sizeof(double));
		}
		g_istart[i] = k;

		if (debug) {
			printf("group %d istart %d len %d limit %d\n", i, istart, len,
			       (i < ng - 1 ? MIN(limit, h->idx[g_istart[i + 1]] - h->idx[istart + len - 1] - 1) : 0));
		}

		int pad = (i < ng - 1 ? MIN(limit, h->idx[g_istart[i + 1]] - h->idx[istart + len - 1] - 1) : 0);
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
		for (int i = 0; i < k; i++) {
			printf("i %d new_idx %1d new_val %f\n", i, new_idx[i], new_val[i]);
		}
	}
	// set g_1's
	g_1[0] = 0;
	for (int g = 0; g < ng; g++) {
		int all_one = 1;
		double *val = g_val[g];
		for (int i = 0; all_one && i < IABS(g_len[g]); i++) {
			all_one = (val[i] == 1.0);
		}
		g_1[g] = all_one;
	}

	if (debug) {
		printf("NEW\nng = %1d\n", ng);
		for (int g = 0; g < ng; g++) {
			printf("group %1d start %1d len %1d g_1 %1d\n", g, g_istart[g], g_len[g], g_1[g]);
			for (int i = 0; i < IABS(g_len[g]); i++) {
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
		for (int g = 0; g < ng; g++) {
			printf("group %1d start %1d len %1d g_1 %1d\n", g, g_istart[g], g_len[g], g_1[g]);
			if (0)
				for (int i = 0; i < IABS(g_len[g]); i++) {
					printf("\t\t\tidx %1d val %g\n", new_idx[g_istart[g] + i], new_val[g_istart[g] + i]);
				}
		}
	}
	// merge
	if (ng > 2) {
		int g = 1;
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
		for (int g = 0; g < ng; g++) {
			printf("group %1d start %1d len %1d g_1 %1d\n", g, g_istart[g], g_len[g], g_1[g]);
			if (0)
				for (int i = 0; i < IABS(g_len[g]); i++) {
					printf("\t\t\tidx %1d val %g\n", new_idx[g_istart[g] + i], new_val[g_istart[g] + i]);
				}
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

	double time_min = 0.0;
	double time_max = 0.0;
	int ntimes = 2;
#if defined(INLA_LINK_WITH_MKL)
	int with_mkl = 1;
#else
	int with_mkl = 0;
#endif

	if (h->preference != IDXVAL_UNKNOWN) {
		return 0;
	}

	double treff[4] = { 0.0, 0.0, 0.0, 0.0 };
	double value[4] = { 0.0, 0.0, 0.0, 0.0 };

	for (int time = 0; time < ntimes; time++) {
		int measure = (time >= 0);
		if (measure) {
			treff[0] -= GMRFLib_cpu();
		}

		value[0] = GMRFLib_dot_product_serial(h, x);
		if (measure) {
			treff[0] += GMRFLib_cpu();
		}
		if (with_mkl) {
			if (measure) {
				treff[1] -= GMRFLib_cpu();
			}
			value[1] = GMRFLib_dot_product_serial_mkl(h, x);
			if (measure) {
				treff[1] += GMRFLib_cpu();
			}
		} else {
			value[1] = value[0];
			treff[1] = treff[0];
		}

		if (measure) {
			treff[2] -= GMRFLib_cpu();
		}
		value[2] = GMRFLib_dot_product_group(h, x);
		if (measure) {
			treff[2] += GMRFLib_cpu();
		}
		if (with_mkl) {
			if (measure) {
				treff[3] -= GMRFLib_cpu();
			}
			value[3] = GMRFLib_dot_product_group_mkl(h, x);
			if (measure) {
				treff[3] += GMRFLib_cpu();
			}
		} else {
			value[3] = value[2];
			treff[3] = treff[2];
		}
	}

	for (int k = 1; k < 4; k++) {
		if (ABS(value[k] - value[0]) > FLT_EPSILON * sqrt(h->n)) {
			P(ABS(value[k] - value[0]));
			P(k);
			P(value[0]);
			P(value[1]);
			P(value[2]);
			P(value[3]);

			printf("n %d\n", h->n);
			for (int i = 0; i < h->n; i++) {
				printf("\tidx[%1d] =  %1d  val = %g\n", i, h->idx[i], h->val[i]);
			}
			printf("ng %d\n", h->g_n);
			for (int g = 0; g < h->g_n; g++) {
				printf("\tg = %d g_1 = %d\n", g, h->g_1[g]);
				for (int i = 0; i < IABS(h->g_len[g]); i++) {
					printf("\t\tidx[%1d] =  %1d  val = %g\n", i, h->g_idx[g][i], h->g_val[g][i]);
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

	k = -1;
	double tmin = GMRFLib_min_value(treff, 4, &k);
	double tmax = GMRFLib_max_value(treff, 4, NULL);

	if (debug) {
		double s = 1.0 / (treff[0] + treff[1] + treff[2] + treff[3]) / ntimes;
		printf("for h with n= %1d chose k=%1d [serial= %.3f serial.mkl= %.3f group= %.3f group.mkl= %.3f]\n",
		       h->n, k, treff[0] * s, treff[1] * s, treff[2] * s, treff[3] * s);
	}

	switch (k) {
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

	if (h->preference == IDXVAL_SERIAL || h->preference == IDXVAL_SERIAL_MKL) {
		/*
		 * no need to keep the group info in the struct 
		 */
		h->g_n = 0;
		Free(h->g_idx);
		Free(h->g_val);
		Free(h->g_len);
		Free(h->g_1);
		for (int k = 0; k < h->g_n_mem; k++) {
			Free(h->g_mem[k]);
		}
		Free(h->g_mem);
		h->g_n_mem = 0;
	}

	time_min += tmin / ntimes;
	time_max += tmax / ntimes;

	return 0;
}

int GMRFLib_idxval_nsort_x(GMRFLib_idxval_tp ** hold, int n, int nt)
{
	int nmax = 1;
	for (int i = 0; i < n; i++) {
		GMRFLib_idxval_tp *h = hold[i];
		if (h->n) {
			nmax = MAX(nmax, h->idx[h->n - 1] + 1);
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
			int len = MAX(2 * nmax, SQR(256));
			double *xx = Calloc(len, double);
			for (int i = 0; i < len; i++) {
				xx[i] = GMRFLib_uniform();
			}
			x_ran = xx;
			len_x_ran = len;
		}
	}

	if (nt > 0) {
#pragma omp parallel for threads(nt)
		for(int k = 0; k < n; k++) {
			GMRFLib_idxval_nsort_x_core(hold[k], x_ran);
		}
	} else {
		for(int k = 0; k < n; k++) {
			GMRFLib_idxval_nsort_x_core(hold[k], x_ran);
		}
	}

	return 0;
}

int GMRFLib_idxval_add(GMRFLib_idxval_tp ** hold, int idx, double val)
{
	if (*hold == NULL) {
		GMRFLib_idxval_create(hold);
	}
	if ((*hold)->n == (*hold)->n_alloc) {
		(*hold)->n_alloc += MAX(IDX_ALLOC_ADD, (*hold)->n / 4);
		(*hold)->idx = Realloc((*hold)->idx, (*hold)->n_alloc, int);
		(*hold)->val = Realloc((*hold)->val, (*hold)->n_alloc, double);
	}
	(*hold)->idx[(*hold)->n] = idx;
	(*hold)->val[(*hold)->n] = val;
	(*hold)->n++;

	return 0;
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
		return 0;
	}
	// FIXME: this should be improved in general, but I think for the usage its ok. Since we are likely to add with same or increasing idx,
	// then I added this 'iaddto' which recall the last index, and try to be a little smarter.
	for (i = (*hold)->iaddto + 1; i < (*hold)->n; i++) {
		if ((*hold)->idx[i] == idx) {
			(*hold)->val[i] += val;
			(*hold)->iaddto = i;
			return 0;
		}
	}
	for (i = 0; i < (*hold)->iaddto; i++) {
		if ((*hold)->idx[i] == idx) {
			(*hold)->val[i] += val;
			(*hold)->iaddto = i;
			return 0;
		}
	}

	if ((*hold)->n == (*hold)->n_alloc) {
		(*hold)->n_alloc += MAX(IDX_ALLOC_ADD, (*hold)->n / 4);
		(*hold)->idx = Realloc((*hold)->idx, (*hold)->n_alloc, int);
		(*hold)->val = Realloc((*hold)->val, (*hold)->n_alloc, double);
	}
	(*hold)->idx[(*hold)->n] = idx;
	(*hold)->val[(*hold)->n] = val;
	(*hold)->n++;

	return 0;
}

double GMRFLib_min_value(double *x, int n, int *idx)
{
	/*
	 * return the MIN(x[]) 
	 */
	int i, imin;
	double min_val;

	min_val = x[0];
	imin = 0;
	for (i = 1; i < n; i++) {
		if (x[i] < min_val) {
			min_val = x[i];
			imin = i;
		}
	}

	if (idx) {
		*idx = imin;
	}

	return min_val;
}
double GMRFLib_max_value(double *x, int n, int *idx)
{
	/*
	 * return the MAX(x[]), optional idx
	 */
	int i, imax;
	double max_val;

	max_val = x[0];
	imax = 0;
	for (i = 1; i < n; i++) {
		if (x[i] > max_val) {
			max_val = x[i];
			imax = i;
		}
	}

	if (idx) {
		*idx = imax;
	}
	return max_val;
}

void my_downheap2_id(int *__restrict data1, double *__restrict data2, const int N, int k)
{
	int v1 = data1[k];
	double v2 = data2[k];

	while (k <= N / 2) {
		int j = 2 * k;
		if (j < N && data1[j] < data1[j + 1]) {
			j++;
		}

		if (!(v1 < data1[j])) {
			break;
		}

		data1[k] = data1[j];
		data2[k] = data2[j];
		k = j;
	}
	data1[k] = v1;
	data2[k] = v2;
}

void gsl_sort2_id(int *__restrict data1, double *__restrict data2, const int n)
{
	int N, k;

	if (n == 0) {
		return;					       /* No data to sort */
	}

	/*
	 * We have n_data elements, last element is at 'n_data-1', first at '0' Set N to the last element number. 
	 */

	N = n - 1;
	k = N / 2;
	k++;						       /* Compensate the first use of 'k--' */
	do {
		k--;
		my_downheap2_id(data1, data2, N, k);
	} while (k > 0);

	while (N > 0) {
		int tmp1 = data1[0];
		data1[0] = data1[N];
		data1[N] = tmp1;

		double tmp2 = data2[0];
		data2[0] = data2[N];
		data2[N] = tmp2;

		/*
		 * then process the heap 
		 */
		N--;
		my_downheap2_id(data1, data2, N, 0);
	}
}


void my_insertionSort_id(int *__restrict iarr, double *__restrict darr, int n)
{
	if (darr) {
		for (int i = 1; i < n; i++) {
			int key = iarr[i];
			double dkey = darr[i];
			int j = i - 1;
			while (j >= 0 && iarr[j] > key) {
				iarr[j + 1] = iarr[j];
				darr[j + 1] = darr[j];
				j--;
			}
			iarr[j + 1] = key;
			darr[j + 1] = dkey;
		}
	} else {
		for (int i = 1; i < n; i++) {
			int key = iarr[i];
			int j = i - 1;
			while (j >= 0 && iarr[j] > key) {
				iarr[j + 1] = iarr[j];
				j--;
			}
			iarr[j + 1] = key;
		}
	}
}

void my_sort2_id(int *__restrict ix, double *__restrict x, int n)
{
	if (n < 128) {
		my_insertionSort_id(ix, x, n);
	} else {
		gsl_sort2_id(ix, x, n);
	}
}

