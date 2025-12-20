#include <assert.h>
#include <omp.h>
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"

#pragma GCC push_options
#pragma GCC optimize("O3")
double GMRFLib_sparse_ddot(int n, double *__restrict v, double *__restrict a, int *__restrict idx)
{
	if (n == 0) return 0.0;
	
	// sum_i v[i] * a[idx[i]]
#if defined(INLA_WITH_MKL)
	return cblas_ddoti(n, v, idx, a);
#elif defined(__linux__) && defined(__SSE2__) && defined(INLA_WITH_INTRINSICS)
	__m128d sum0 = _mm_setzero_pd();
	__m128d sum1 = _mm_setzero_pd();
	int i = 0;
	int limit = n & ~3;
	for (; i < limit; i += 4) {
		__m128d sparse1 = _mm_load_sd(&a[idx[i]]);
		__m128d sparse2 = _mm_load_sd(&a[idx[i+1]]);
		__m128d sparse3 = _mm_load_sd(&a[idx[i+2]]);
		__m128d sparse4 = _mm_load_sd(&a[idx[i+3]]);
		__m128d sparse_vals0 = _mm_unpacklo_pd(sparse1, sparse2);
		__m128d sparse_vals1 = _mm_unpacklo_pd(sparse3, sparse4);
		__m128d dense_vals0 = _mm_load_pd(&v[i]);
		__m128d dense_vals1 = _mm_load_pd(&v[i + 2]);
		__m128d prod0 = _mm_mul_pd(sparse_vals0, dense_vals0);
		__m128d prod1 = _mm_mul_pd(sparse_vals1, dense_vals1);
		sum0 = _mm_add_pd(sum0, prod0);
		sum1 = _mm_add_pd(sum1, prod1);
	}
	sum0 = _mm_add_pd(sum0, sum1);
	__m128d sum_swapped = _mm_shuffle_pd(sum0, sum0, 1);
	__m128d sum_total = _mm_add_pd(sum0, sum_swapped);
	double result;
	_mm_store_sd(&result, sum_total);
	for (; i < n; i++) {
		result += a[idx[i]] * v[i];
	}
	return result;
#else
	double s0 = 0.0;
#pragma omp simd reduction(+: s0)
	for (int i = 0; i < n; i++) {
		s0 += v[i] * a[idx[i]];
	}
	return s0;
#endif
}
#pragma GCC pop_options

double GMRFLib_sparse_ddot_ddot_(GMRFLib_idxval_tp *__restrict ELM_, double *__restrict ARR_)
{
	// special case: ->idx == sequential
	return (GMRFLib_ddot(ELM_->n, ELM_->val, ARR_ + ELM_->idx[0]));
}

double GMRFLib_sparse_ddot_sum_(GMRFLib_idxval_tp *__restrict ELM_, double *__restrict ARR_)
{
	// special case: ->idx == sequential and all(->val == 1.0)
	return (GMRFLib_dsum(ELM_->n, ARR_ + ELM_->idx[0]));
}

double GMRFLib_sparse_ddot_(GMRFLib_idxval_tp *__restrict ELM_, double *__restrict ARR_)
{
#if defined(INLA_WITH_ARMPL)
	if (ELM_->spvec) {
		double res = 0.0;
		armpl_status_t info = armpl_spdot_exec_d(ELM_->spvec, ARR_, &res);
		assert(info == ARMPL_STATUS_SUCCESS);
		return (res);
	}
#endif
	double *__restrict vv_ = ELM_->val;
	double *__restrict aa_ = ARR_;
	int *__restrict idx_ = ELM_->idx;
	if (ELM_->dot_product_func) {
		return (ELM_->dot_product_func(ELM_, ARR_));
	} else {
		return (GMRFLib_sparse_ddot(ELM_->n, vv_, aa_, idx_));
	}
}

double GMRFLib_sparse_ddot_group_(GMRFLib_idxval_tp *__restrict ELM_, double *__restrict ARR_)
{
	double value = 0.0;
	const int g_n = ELM_->g_n;
#define USE_PREFETCH 1
#if USE_PREFETCH
	if (__builtin_expect(g_n > 0, 1)) {
		__builtin_prefetch(ELM_->g_len, 0, 3);
		__builtin_prefetch(ELM_->g_idx[0], 0, 3);
		__builtin_prefetch(ELM_->g_val[0], 0, 3);
	}
#endif

	for (int g_ = 0; g_ < g_n; g_++) {
#if USE_PREFETCH
		if (__builtin_expect(g_ + 1 < g_n, 1)) {
			__builtin_prefetch(ELM_->g_idx[g_ + 1], 0, 3);
			__builtin_prefetch(ELM_->g_val[g_ + 1], 0, 3);
		}
#endif
		const int len_ = ELM_->g_len[g_];
		int *__restrict const ii_ = ELM_->g_idx[g_];
		double *__restrict const vv_ = ELM_->g_val[g_];

		if (__builtin_expect(len_ > 0, 1)) {
			double *__restrict const aa_ = &(ARR_[0]);
			if (__builtin_expect(ELM_->g_1[g_], 0)) {
				value += GMRFLib_sparse_dsum(len_, aa_, ii_);
			} else {
#if defined(INLA_WITH_ARMPL)
				// the way we store it, then this happens either with g_=0 or does not happens at all
				if ((g_ == 0) && ELM_->spvec_g) {
					double res = 0.0;
					armpl_status_t info = armpl_spdot_exec_d(ELM_->spvec_g, ARR_, &res);
					assert(info == ARMPL_STATUS_SUCCESS);
					value += res;
				} else {
					value += GMRFLib_sparse_ddot(len_, vv_, aa_, ii_);
				}
#else
				value += GMRFLib_sparse_ddot(len_, vv_, aa_, ii_);
#endif
			}
		} else if (__builtin_expect(len_ < 0, 1)) {
			const int llen_ = -len_;
			double *__restrict const aa_ = &(ARR_[ii_[0]]);
			if (__builtin_expect(ELM_->g_1[g_], 1)) {
				value += GMRFLib_dsum(llen_, aa_);
			} else {
				value += GMRFLib_ddot(llen_, vv_, aa_);
			}
		}
	}
#undef USE_PREFETCH
	return value;
}

double GMRFLib_sparse_ddot_group_simple_(GMRFLib_idxval_tp *__restrict ELM_, double *__restrict ARR_)
{
	// in this case, there is only one group giving dense vector calls to 'dot' or 'sum'
	int llen_ = IABS(ELM_->g_len[0]);
	int *__restrict const ii_ = ELM_->g_idx[0];
	double *__restrict const vv_ = ELM_->g_val[0];
	double *__restrict const aa_ = &(ARR_[ii_[0]]);
	double value = 0.0;
	if (__builtin_expect(ELM_->g_1[0], 1)) {
		value = GMRFLib_dsum(llen_, aa_);
	} else {
		value = GMRFLib_ddot(llen_, vv_, aa_);
	}
	return value;
}
