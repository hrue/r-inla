#include <assert.h>
#include <omp.h>
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"

#define SPARSE_DOT()							\
	double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;			\
	int unroll = 8;							\
	int m = n & ~(unroll - 1);					\
	for (int i = 0; i < m; i += unroll) {				\
		s0 += v[i + 0] * a[idx[i + 0]];				\
		s1 += v[i + 1] * a[idx[i + 1]];				\
		s2 += v[i + 2] * a[idx[i + 2]];				\
		s3 += v[i + 3] * a[idx[i + 3]];				\
		s0 += v[i + 4] * a[idx[i + 4]];				\
		s1 += v[i + 5] * a[idx[i + 5]];				\
		s2 += v[i + 6] * a[idx[i + 6]];				\
		s3 += v[i + 7] * a[idx[i + 7]];				\
	}								\
	for (int i = m; i < n; i++) {					\
		s0 += v[i] * a[idx[i]];					\
	}								\
	return (s0 + s1) + (s2 + s3)

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((optimize("O3")))
    __attribute__((target_clones(INLA_CLONE_TARGETS "default")))
double GMRFLib_sparse_ddot(int n, double *__restrict v, double *__restrict a, int *__restrict idx)
{
	// sum_i v[i] * a[idx[i]]
#if defined(INLA_WITH_MKL)
	return cblas_ddoti(n, v, idx, a);
#else
	SPARSE_DOT();
#endif
}
#pragma GCC diagnostic pop
#undef SPARSE_DOT

double GMRFLib_sparse_ddot_ddot_(GMRFLib_idxval_tp *__restrict ELM_, double *__restrict ARR_)
{
	// special case: ->idx == sequential
	return (GMRFLib_ddot(ELM_->n, ELM_->val, ARR_ + ELM_->idx[0]));
}

double GMRFLib_sparse_ddot_sum_(GMRFLib_idxval_tp *__restrict ELM_, double *__restrict ARR_)
{
	// special case: ->idx == sequential and all(->val == 1.0)
	// return (GMRFLib_dsum(ELM_->n, ARR_ + ELM_->idx[0]));
#if defined(INLA_WITH_SIMDE)
	int n = ELM_->n;
	if (n < 16) {
		double *x = ARR_ + ELM_->idx[0];
		// it returns inside here:
#       include "intrinsics/simde/dsum-sse2-small.h"
	} 
#endif
	return (GMRFLib_dsum(ELM_->n, ARR_ + ELM_->idx[0]));
}

// special cases
#pragma GCC diagnostic push
__attribute__((optimize("fast-math")))
double GMRFLib_sparse_ddot_sum1_(GMRFLib_idxval_tp *__restrict ELM_, double *__restrict ARR_)
{
	double * x = ARR_ + ELM_->idx[0];
	return x[0];
}
double GMRFLib_sparse_ddot_sum2_(GMRFLib_idxval_tp *__restrict ELM_, double *__restrict ARR_)
{
	double * x = ARR_ + ELM_->idx[0];
	return x[0] + x[1];
}
double GMRFLib_sparse_ddot_sum3_(GMRFLib_idxval_tp *__restrict ELM_, double *__restrict ARR_)
{
	double * x = ARR_ + ELM_->idx[0];
	return x[0] + x[1] + x[2];
}
double GMRFLib_sparse_ddot_sum4_(GMRFLib_idxval_tp *__restrict ELM_, double *__restrict ARR_)
{
	double * x = ARR_ + ELM_->idx[0];
	return (x[0] + x[1]) + (x[2] + x[3]);
}
double GMRFLib_sparse_ddot_sum5_(GMRFLib_idxval_tp *__restrict ELM_, double *__restrict ARR_)
{
	double * x = ARR_ + ELM_->idx[0];
	return ((x[0] + x[1]) + (x[2] + x[3])) + x[4];
}
double GMRFLib_sparse_ddot_sum6_(GMRFLib_idxval_tp *__restrict ELM_, double *__restrict ARR_)
{
	double * x = ARR_ + ELM_->idx[0];
	return (x[0] + x[1] + x[2]) + (x[3] + x[4] + x[5]);
}
double GMRFLib_sparse_ddot_sum7_(GMRFLib_idxval_tp *__restrict ELM_, double *__restrict ARR_)
{
	double * x = ARR_ + ELM_->idx[0];
	return ((x[0] + x[1] + x[2]) + (x[3] + x[4] + x[5])) + x[6];
}
#pragma GCC diagnostic pop

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
	for (int g_ = 0; g_ < g_n; g_++) {
		const int len_ = ELM_->g_len[g_];
		int *__restrict const ii_ = ELM_->g_idx[g_];
		double *__restrict const vv_ = ELM_->g_val[g_];

		if (len_ > 0) {
			double *__restrict const aa_ = &(ARR_[0]);
			if (ELM_->g_1[g_]) {
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
		} else if (len_ < 0) {
			const int llen_ = -len_;
			double *__restrict const aa_ = &(ARR_[ii_[0]]);
			if (ELM_->g_1[g_]) {
				value += GMRFLib_dsum(llen_, aa_);
			} else {
				value += GMRFLib_ddot(llen_, vv_, aa_);
			}
		}
	}
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
	if (ELM_->g_1[0]) {
		value = GMRFLib_dsum(llen_, aa_);
	} else {
		value = GMRFLib_ddot(llen_, vv_, aa_);
	}
	return value;
}
