#include <assert.h>
#include <omp.h>
#include <stdlib.h>

#if defined(__linux__) && defined(__AVX2__)
#include <immintrin.h>					       // For AVX/SSE intrinsics
#endif

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"
#include "GMRFLib/dot.h"

double GMRFLib_dot_product_opt(GMRFLib_idxval_tp *__restrict ELM_, double *__restrict ARR_)
{
	// Use __builtin_expect for better branch prediction on hot path
	if (__builtin_expect(ELM_->dot_product_func != NULL, 1)) {
#if defined(INLA_WITH_DEVEL) && !defined(INLA_WITH_ARMPL)
		if (__builtin_expect(GMRFLib_dot_product_gain >= 0.0, 0)) {
#pragma omp atomic
			GMRFLib_dot_product_gain += ELM_->cpu_gain;
		}
#endif
		return (ELM_->dot_product_func(ELM_, ARR_));
	} else {
		return GMRFLib_dot_product_sparse_mkl(ELM_, ARR_);
	}
}

double GMRFLib_dot_product_group_prefetch(GMRFLib_idxval_tp *__restrict ELM_, double *__restrict ARR_)
{
	double value_ = 0.0;
	const int g_n = ELM_->g_n;

	// Prefetch first group data
	if (__builtin_expect(g_n > 0, 1)) {
		__builtin_prefetch(ELM_->g_len, 0, 3);
		__builtin_prefetch(ELM_->g_idx[0], 0, 3);
		__builtin_prefetch(ELM_->g_val[0], 0, 3);
	}

	for (int g_ = 0; g_ < g_n; g_++) {
		// Prefetch next iteration data
		if (__builtin_expect(g_ + 1 < g_n, 1)) {
			__builtin_prefetch(ELM_->g_idx[g_ + 1], 0, 3);
			__builtin_prefetch(ELM_->g_val[g_ + 1], 0, 3);
		}

		const int len_ = ELM_->g_len[g_];
		int *__restrict const ii_ = ELM_->g_idx[g_];
		double *__restrict const vv_ = ELM_->g_val[g_];

		if (__builtin_expect(len_ > 0, 1)) {
			double *__restrict const aa_ = &(ARR_[0]);
			if (__builtin_expect(ELM_->g_1[g_], 1)) {
				value_ += GMRFLib_dsum_idx_opt(len_, aa_, ii_);
			} else {
				value_ += GMRFLib_ddot_idx_opt(len_, vv_, aa_, ii_);
			}
		} else if (__builtin_expect(len_ < 0, 0)) {
			const int llen_ = -len_;
			double *__restrict const aa_ = &(ARR_[ii_[0]]);
			if (__builtin_expect(ELM_->g_1[g_], 1)) {
				value_ += GMRFLib_dsum(llen_, aa_);
			} else {
				value_ += GMRFLib_ddot_opt(llen_, vv_, aa_);
			}
		}
	}
	return value_;
}

double GMRFLib_ddot_opt(int n, double *__restrict x, double *__restrict y)
{
	if (__builtin_expect(n <= 8, 0)) {
		// Small arrays - manual unroll
		aligned_double(r) = 0.0;
		switch (n) {
		case 8:
			r += x[7] * y[7];
			__attribute__((fallthrough));
		case 7:
			r += x[6] * y[6];
			__attribute__((fallthrough));
		case 6:
			r += x[5] * y[5];
			__attribute__((fallthrough));
		case 5:
			r += x[4] * y[4];
			__attribute__((fallthrough));
		case 4:
			r += x[3] * y[3];
			__attribute__((fallthrough));
		case 3:
			r += x[2] * y[2];
			__attribute__((fallthrough));
		case 2:
			r += x[1] * y[1];
			__attribute__((fallthrough));
		case 1:
			r += x[0] * y[0];
		}
		return r;
	} else if (__builtin_expect(n <= 16, 0)) {
		// Medium arrays - OpenMP SIMD
		aligned_double(r) = 0.0;
#pragma omp simd reduction(+: r)
		for (int i = 0; i < n; i++) {
			r += x[i] * y[i];
		}
		return r;
	} else {
		// Large arrays - use BLAS
		int one = 1;
		return ddot_(&n, x, &one, y, &one);
	}
}

double GMRFLib_ddot_idx_opt(int n, double *__restrict v, double *__restrict a, int *__restrict idx)
{
#define ROLL8 8
#define ROLL16 16
	double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
	int i = 0;
	if (n < ROLL16) {
		for (; i + ROLL8 <= n; i += ROLL8) {
			// Prefetch next cache lines
			__builtin_prefetch(v + i + ROLL8, 0, 3);
			__builtin_prefetch(idx + i + ROLL8, 0, 3);

			const double *vv = v + i;
			const int *iidx = idx + i;

			s0 += vv[0] * a[iidx[0]];
			s1 += vv[1] * a[iidx[1]];
			s2 += vv[2] * a[iidx[2]];
			s3 += vv[3] * a[iidx[3]];

			s0 += vv[4] * a[iidx[4]];
			s1 += vv[5] * a[iidx[5]];
			s2 += vv[6] * a[iidx[6]];
			s3 += vv[7] * a[iidx[7]];
		}
		for (; i < n; i++) {
			s0 += v[i] * a[idx[i]];
		}
		return s0 + s1 + s2 + s3;
	} else {
		double s4 = 0.0, s5 = 0.0, s6 = 0.0, s7 = 0.0;
		for (; i + ROLL16 <= n; i += ROLL16) {
			// Prefetch next cache lines
			__builtin_prefetch(v + i + ROLL16, 0, 3);
			__builtin_prefetch(idx + i + ROLL16, 0, 3);

			const double *vv = v + i;
			const int *iidx = idx + i;

			s0 += vv[0] * a[iidx[0]];
			s1 += vv[1] * a[iidx[1]];
			s2 += vv[2] * a[iidx[2]];
			s3 += vv[3] * a[iidx[3]];
			s4 += vv[4] * a[iidx[4]];
			s5 += vv[5] * a[iidx[5]];
			s6 += vv[6] * a[iidx[6]];
			s7 += vv[7] * a[iidx[7]];

			s0 += vv[8] * a[iidx[8]];
			s1 += vv[9] * a[iidx[9]];
			s2 += vv[10] * a[iidx[10]];
			s3 += vv[11] * a[iidx[11]];
			s4 += vv[12] * a[iidx[12]];
			s5 += vv[13] * a[iidx[13]];
			s6 += vv[14] * a[iidx[14]];
			s7 += vv[15] * a[iidx[15]];
		}
		for (; i < n; i++) {
			s0 += v[i] * a[idx[i]];
		}
		return s0 + s1 + s2 + s3 + s4 + s5 + s6 + s7;
	}
#undef ROLL8
#undef ROLL16
}


// Template-like macro for generating optimized versions
#define DEFINE_OPTIMIZED_GROUP_FUNC(SUFFIX, DOT_FUNC) \
	double GMRFLib_dot_product_group_##SUFFIX(GMRFLib_idxval_tp *__restrict ELM_, double *__restrict ARR_) \
	{								\
		double value_ = 0.0;					\
		const int g_n = ELM_->g_n;				\
									\
		_Pragma("omp simd reduction(+: value_)")		\
			for (int g_ = 0; g_ < g_n; g_++) {		\
				const int len_ = ELM_->g_len[g_];	\
				int * __restrict const ii_ = ELM_->g_idx[g_]; \
				double * __restrict const vv_ = ELM_->g_val[g_]; \
				if (__builtin_expect(len_ > 0, 1)) {	\
					double * __restrict const aa_ = &(ARR_[0]); \
					value_ += (ELM_->g_1[g_] ? GMRFLib_dsum_idx_opt(len_, aa_, ii_) : DOT_FUNC(len_, vv_, aa_, ii_)); \
				} else if (__builtin_expect(len_ < 0, 0)) { \
					const int llen_ = -len_;	\
					double * __restrict const aa_ = &(ARR_[ii_[0]]); \
					value_ += (ELM_->g_1[g_] ? GMRFLib_dsum(llen_, aa_) : GMRFLib_ddot_opt(llen_, vv_, aa_)); \
				}					\
				}					\
		return value_;						\
		}

// Generate optimized versions
DEFINE_OPTIMIZED_GROUP_FUNC(mkl_opt, GMRFLib_ddot_idx_mkl);
DEFINE_OPTIMIZED_GROUP_FUNC(sparse_opt, GMRFLib_ddot_idx_opt);

double GMRFLib_dot_product_sparse_opt(GMRFLib_idxval_tp *__restrict ELM_, double *__restrict ARR_)
{
	double *__restrict vv_ = ELM_->val;
	double *__restrict aa_ = ARR_;
	int *__restrict idx_ = ELM_->idx;

	return GMRFLib_ddot_idx_opt(ELM_->n, vv_, aa_, idx_);
}

double GMRFLib_dot_product_sparse_mkl(GMRFLib_idxval_tp *__restrict ELM_, double *__restrict ARR_)
{
	double *__restrict vv_ = ELM_->val;
	double *__restrict aa_ = ARR_;
	int *__restrict idx_ = ELM_->idx;

	return GMRFLib_ddot_idx_mkl(ELM_->n, vv_, aa_, idx_);
}

double GMRFLib_ddot(int n, double *x, double *y)
{
	return GMRFLib_ddot_opt(n, x, y);
}

double GMRFLib_ddot_idx(int n, double *__restrict v, double *__restrict a, int *__restrict idx)
{
#if defined(INLA_WITH_MKL)
	return cblas_ddoti(n, v, idx, a);
#else
	return GMRFLib_ddot_idx_opt(n, v, a, idx);
#endif
}

double GMRFLib_ddot_idx_mkl(int n, double *__restrict v, double *__restrict a, int *__restrict idx)
{
#if defined(INLA_WITH_MKL)
	return cblas_ddoti(n, v, idx, a);
#else
	return GMRFLib_ddot_idx_opt(n, v, a, idx);
#endif
}

double GMRFLib_dot_product_sparse_armpl(GMRFLib_idxval_tp *__restrict ELM_, double *__restrict ARR_)
{
#if defined(INLA_WITH_ARMPL)
	double res = 0.0;
	armpl_status_t info = armpl_spdot_exec_d(ELM_->spvec, ARR_, &res);
	assert(info == ARMPL_STATUS_SUCCESS);
	return (res);
#else
	return (GMRFLib_dot_product_sparse_mkl(ELM_, ARR_));
#endif
}
