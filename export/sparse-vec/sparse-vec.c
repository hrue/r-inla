
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

		int istart_ = ELM_->g_i[g_];
		int *__restrict ii_ = &(ELM_->g_idx[istart_]);
		double *__restrict vv_ = &(ELM_->g_val[istart_]);

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

		int istart_ = ELM_->g_i[g_];
		int *__restrict ii_ = &(ELM_->g_idx[istart_]);
		double *__restrict vv_ = &(ELM_->g_val[istart_]);

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
	return 0;
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
	return GMRFLib_idxval_nsort_x(hold, n, nt, 0);
}

int GMRFLib_idxval_prepare(GMRFLib_idxval_tp * hold)
{
	return (GMRFLib_idxval_nsort_x(&hold, 1, 1, -1));
}

int GMRFLib_idxval_nsort_x(GMRFLib_idxval_tp ** hold, int n, int nt, int prune_zeros)
{
	// prune_zeros, 0=do nothing, >0 remove duplicate with same index within each group, <0 remove all zeros
	// prune_zeros and 'group'-indexing is limited to vectors with length > 'limit'
	const int limit = 16;
	const int debug = 0;

	int nmax = 1;
	for (int i = 0; i < n; i++) {
		GMRFLib_idxval_tp *h = hold[i];
		nmax = MAX(nmax, h->n);
	}

	double *tmp_work = Calloc(nmax, double);

	for (int i = 0; i < n; i++) {
		GMRFLib_idxval_tp *h = hold[i];
		if (!h)
			continue;
		double tref[5] = { 0, 0, 0, 0, 0 };
		if (debug)
			tref[0] -= omp_get_wtime();
		if (h->n > 1) {
			int is_sorted = 1;
			for (int j = 1; j < h->n && is_sorted; j++) {
				is_sorted = is_sorted && (h->idx[j] >= h->idx[j - 1]);
			}
			if (!is_sorted) {
				double *idx_tmp = tmp_work;
				for (int j = 0; j < h->n; j++) {
					idx_tmp[j] = (double) h->idx[j];
				}
				gsl_sort2(idx_tmp, (size_t) 1, h->val, (size_t) 1, (size_t) h->n);
				for (int j = 0; j < h->n; j++) {
					h->idx[j] = (int) idx_tmp[j];
				}
			}
		}
		if (debug)
			tref[0] += omp_get_wtime();

		if (h->n <= limit) {
			h->preference = IDXVAL_SERIAL;
			continue;
		}

		/*
		 * build basic groups with one group for each sequence and then one for each individual 
		 */
		if (debug)
			tref[1] -= omp_get_wtime();
		int ng = 1;
		for (int j = 1; j < h->n; j++) {
			if (h->idx[j] != h->idx[j - 1] + 1) {
				ng++;
			}
		}
		int *g_i = Calloc(ng + 1, int);
		int *g_len = Calloc(ng + 1, int);

		int k = 0;
		g_i[0] = 0;
		for (int j = 1; j < h->n; j++) {
			if (h->idx[j] != h->idx[j - 1] + 1) {
				g_len[k] = j - g_i[k];
				k++;
				g_i[k] = j;
			}
		}
		g_len[ng - 1] = h->n - g_i[ng - 1];
		if (debug)
			tref[1] += omp_get_wtime();

		if (1 && debug) {
			GMRFLib_idxval_printf(stdout, h, "OLD");
		}

		k = h->n;				       /* so the code below works */
		int kg = ng;				       /* so the code below works */

		if (debug)
			tref[2] -= omp_get_wtime();
		if (prune_zeros) {
			/*
			 * remove zeros in 'val', either all or all but the first one 
			 */
			k = 0;
			kg = 0;
			int nskip = 0;
			for (int g = 0; g < ng; g++) {

				int gl = IABS(g_len[g]);
				int num_zero = 0;
				int idx_zz = -1;
				int gn = 0;
				int kk = k;

				for (int j = 0; j < gl; j++) {
					int idx_z = h->idx[g_i[g] + j];

					if (prune_zeros < 0) {
						if (ISZERO(h->val[g_i[g] + j])) {
							nskip++;
						} else {
							h->idx[k] = idx_z;
							h->val[k] = h->val[g_i[g] + j];
							k++;
							gn++;
						}
					} else {
						if (idx_z == idx_zz) {
							if (ISZERO(h->val[g_i[g] + j])) {
								num_zero++;
								if (prune_zeros > 0 && num_zero == 1) {
									h->idx[k] = idx_z;
									h->val[k] = h->val[g_i[g] + j];
									k++;
									gn++;
								} else {
									nskip++;
								}
							} else {
								h->idx[k] = idx_z;
								h->val[k] = h->val[g_i[g] + j];
								k++;
								gn++;
							}
						} else {
							h->idx[k] = idx_z;
							h->val[k] = h->val[g_i[g] + j];
							k++;
							gn++;
							num_zero = 0;
							if (ISZERO(h->val[g_i[g] + j])) {
								num_zero++;
							}
						}
						idx_zz = idx_z;
					}
				}
				if (debug && (g_len[g] != gn)) {
					printf("modify group %d from %d to %d, nskip=%1dn", g, gl, gn, nskip);
				}
				if (gn) {
					g_i[kg] = kk;
					g_len[kg] = (g_len[g] < 0 ? -gn : gn);
					kg++;
				}
			}
		}
		h->n = k;
		h->g_n = kg;
		h->g_i = g_i;
		h->g_len = g_len;

		if (debug)
			tref[2] += omp_get_wtime();

		if (1 && debug) {
			GMRFLib_idxval_printf(stdout, h, "NEW");
		}

		/*
		 * grow the groups together by merging individual groups, and merging sequential groups 
		 */
		if (debug)
			tref[3] -= omp_get_wtime();
		int gg = 0;
		int irregular = 0;
		for (int g = 0; g < h->g_n; g++) {
			if (h->g_len[g] == 1) {
				if (irregular) {
					h->g_len[gg]++;
				} else {
					h->g_len[gg] = 1;
					h->g_i[gg] = h->g_i[g];
				}
				irregular = 1;
			} else if (IABS(h->g_len[g]) > 0) {
				if (irregular) {
					gg++;
				}
				irregular = 0;
				h->g_len[gg] = -IABS(h->g_len[g]);
				h->g_i[gg] = h->g_i[g];
				gg++;
			} else {
				continue;
			}
		}
		h->g_n = gg + irregular;

		int fill = 8L;
		if (h && h->n > 1) {
			int len_extra = 0;
			for (int g = 0; g < h->g_n - 1; g++) {
				if (h->g_len[g] < 0 && h->g_len[g + 1] < 0) {
					int last_this = h->idx[h->g_i[g] + IABS(h->g_len[g]) - 1];
					int first_next = h->idx[h->g_i[g + 1]];
					if (first_next - last_this < fill) {
						len_extra += first_next - last_this - 1;
						if (debug) {
							printf("tmerge group %1d and %1d with dist %1d and lengths %1d and %1dn",
							       g, g + 1, first_next - last_this, h->g_len[g], h->g_len[g + 1]);
						}
					}
				}
			}

			if (len_extra > 0) {

				int n_new = h->n + len_extra;
				int *idx_new = Calloc(n_new, int);
				double *val_new = Calloc(n_new, double);

				int kk = 0;
				int gg = 0;

				int this_len = h->g_len[gg];
				int this_alen = IABS(this_len);

				Memcpy(idx_new + kk, h->idx, this_alen * sizeof(int));
				Memcpy(val_new + kk, h->val, this_alen * sizeof(double));

				kk += this_alen;
				h->g_i[gg] = h->g_i[0];
				h->g_len[gg] = h->g_len[0];

				int gg_len = h->g_len[gg];
				int gg_alen = IABS(gg_len);
				int gg_first_i = h->g_i[gg];
				int gg_last_i = gg_first_i + gg_alen - 1;
				int gg_last_idx = idx_new[gg_last_i];
				int pending = 0;

				if (debug) {
					printf("n=%1d len_extra=%1d n_new=%1dn", h->n, len_extra, n_new);
				}

				for (int g = 1; g < h->g_n; g++) {

					if (debug) {
						printf("tprocess group=%1d with len=%d kk=%1dn", g, h->g_len[g], kk);
					}

					/*
					 * merge this group into the current new one or just add it ? 
					 */
					int g_len = h->g_len[g];
					int g_alen = IABS(g_len);
					int g_first_i = h->g_i[g];
					int g_first_idx = h->idx[g_first_i];

					int idx_diff = g_first_idx - gg_last_idx;
					if (debug) {
						printf("g_len=%1d gg_len=%1d idx_diff=%1dn", g_len, gg_len, idx_diff);
					}

					if (g_len < 0 && gg_len < 0 && idx_diff < fill) {
						if (debug) {
							printf("tmerge group g=%1d into gg=%1dn", g, gg);
						}

						/*
						 * yes, merge it       
						 */
						int len_fill = idx_diff - 1;
						for (int j = 0; j < len_fill; j++) {
							idx_new[kk + j] = idx_new[kk - 1] + 1 + j;
						}
						if (debug) {
							printf("ttlen_fill=%1dn", len_fill);
						}

						h->g_len[gg] -= len_fill;
						kk += len_fill;
						Memcpy(idx_new + kk, h->idx + g_first_i, g_alen * sizeof(int));
						Memcpy(val_new + kk, h->val + g_first_i, g_alen * sizeof(double));
						h->g_len[gg] -= g_alen;
						kk += g_alen;

						gg_len = h->g_len[gg];
						gg_alen = IABS(gg_len);
						gg_first_i = h->g_i[gg];
						gg_last_i = gg_first_i + gg_alen - 1;
						gg_last_idx = idx_new[gg_last_i];

						if (debug) {
							printf("ttappend group g=%1d with len %1dn", g, g_len);
							printf("ttgg_last_idx=%1d kk=%1dn", gg_last_idx, kk);
						}
						pending = 1;
					} else {
						/*
						 * no, just add a new group 
						 */
						gg++;
						Memcpy(idx_new + kk, h->idx + g_first_i, g_alen * sizeof(int));
						Memcpy(val_new + kk, h->val + g_first_i, g_alen * sizeof(double));
						h->g_len[gg] = g_len;
						h->g_i[gg] = kk;
						kk += g_alen;

						if (debug) {
							printf("tmake a new group, gg=%1d with len=%1dn", gg, g_len);
						}

						gg_len = h->g_len[gg];
						gg_alen = IABS(gg_len);
						gg_first_i = h->g_i[gg];
						gg_last_i = gg_first_i + gg_alen - 1;
						gg_last_idx = idx_new[gg_last_i];
						pending = (g < h->g_n - 1 ? 0 : 1);
					}
					if (debug) {
						printf("gg=%1d kk=%1dn", gg, kk);
					}
				}

				h->n_n = n_new;
				h->g_n = gg + pending;
				h->g_idx = idx_new;
				h->g_val = val_new;
				h->free_g_mem = 1;
			} else {
				h->n_n = h->n;
				h->g_idx = h->idx;
				h->g_val = h->val;
				h->free_g_mem = 0;
			}
		}
		if (debug)
			tref[3] += omp_get_wtime();

		/*
		 * add a boolean for all(val[]==1), which makes dot-products into sums 
		 */
		if (debug)
			tref[4] -= omp_get_wtime();
		int *g_1 = Calloc(h->g_n + 1, int);
		for (int g = 0; g < h->g_n; g++) {
			int all_1 = 1;
			for (int j = 0; j < IABS(h->g_len[g]) && all_1; j++) {
				all_1 = all_1 && ISEQUAL(h->g_val[h->g_i[g] + j], 1.0);
			}
			g_1[g] = (all_1 ? 1 : 0);
		}
		h->g_1 = g_1;
		if (debug)
			tref[4] += omp_get_wtime();

		if (debug) {
			GMRFLib_idxval_info_printf(stdout, h, "t");
		}
		if (debug) {
			double sum = tref[0] + tref[1] + tref[2] + tref[3] + tref[4];
			printf("%.3f %.3f %.3f %.3f %.3fn", tref[0] / sum, tref[1] / sum, tref[2] / sum, tref[3] / sum, tref[4] / sum);
		}
	}

	Free(tmp_work);

	/*
	 * Add a tag about which option that is faster, the group or the serial algorithm, for each 'idxval'. I'm a little reluctant about doing
	 * this within the parallel loop. Usually, this is a very quick procedure, so it does not really matter...
	 */

	nmax = 8;
	for (int i = 0; i < n; i++) {
		GMRFLib_idxval_tp *h = hold[i];
		if (h->n) {
			nmax = MAX(nmax, h->idx[h->n - 1] + 1);
		}
	}

	double *x = Calloc(nmax, double);
	assert(x);
	for (int i = 0; i < nmax; i++) {
		x[i] = GMRFLib_uniform();
	}

	double time_min = 0.0;
	double time_max = 0.0;
	int ntimes = 1;

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
			if (ABS(value[k] - value[0]) > FLT_EPSILON) {
				P(ABS(value[k] - value[0]));
				P(k);
				P(value[0]);
				P(value[1]);
				P(value[2]);
				P(value[3]);
				assert(0 == 1);
			}
		}

		int k = -1;
		double tmin = GMRFLib_min_value(tref, 4, &k);
		double tmax = GMRFLib_max_value(tref, 4, NULL);

		if (0) {
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

		if (1) {
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

		time_min += tmin / ntimes;
		time_max += tmax / ntimes;
	}

	if (0) {
		if (time_min > 0.0) {
			printf("idxval opt: saving %.6f seconds/M.evals, %.2f%% improvement\n",
			       (time_max - time_min) * SQR(1024.0), 100.0 * (1.0 - time_min / time_max));
		}
	}

	return 0;
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
