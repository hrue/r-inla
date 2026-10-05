
#ifndef __GMRFLib_DOT_H__
#       define __GMRFLib_DOT_H__

#       include <stdlib.h>
#       include <stddef.h>
#       include <math.h>

#       undef __BEGIN_DECLS
#       undef __END_DECLS
#       ifdef __cplusplus
#              define __BEGIN_DECLS extern "C" {
#              define __END_DECLS }
#       else
#              define __BEGIN_DECLS			       /* empty */
#              define __END_DECLS			       /* empty */
#       endif

__BEGIN_DECLS
#       include "GMRFLib/GMRFLibP.h"
#       if defined(INLA_WITH_ARMPL)
#              include "armpl_sparse.h"
#       endif


double GMRFLib_dsum(int n, double *x);
double GMRFLib_dsum_ext(int n, double *x);
double GMRFLib_sparse_ddot(int n, double *RESTRICT v, double *RESTRICT a, int *RESTRICT idx);
double GMRFLib_sparse_ddot_(GMRFLib_idxval_tp * RESTRICT ELM_, double *RESTRICT ARR_);
double GMRFLib_sparse_ddot_ddot_(GMRFLib_idxval_tp * RESTRICT ELM_, double *RESTRICT ARR_);
double GMRFLib_sparse_ddot_group_(GMRFLib_idxval_tp * RESTRICT ELM_, double *RESTRICT ARR_);
double GMRFLib_sparse_ddot_group_simple_(GMRFLib_idxval_tp * RESTRICT ELM_, double *RESTRICT ARR_);
double GMRFLib_sparse_ddot_sum1_(GMRFLib_idxval_tp * RESTRICT ELM_, double *RESTRICT ARR_);
double GMRFLib_sparse_ddot_sum2_(GMRFLib_idxval_tp * RESTRICT ELM_, double *RESTRICT ARR_);
double GMRFLib_sparse_ddot_sum3_(GMRFLib_idxval_tp * RESTRICT ELM_, double *RESTRICT ARR_);
double GMRFLib_sparse_ddot_sum4_(GMRFLib_idxval_tp * RESTRICT ELM_, double *RESTRICT ARR_);
double GMRFLib_sparse_ddot_sum5_(GMRFLib_idxval_tp * RESTRICT ELM_, double *RESTRICT ARR_);
double GMRFLib_sparse_ddot_sum6_(GMRFLib_idxval_tp * RESTRICT ELM_, double *RESTRICT ARR_);
double GMRFLib_sparse_ddot_sum7_(GMRFLib_idxval_tp * RESTRICT ELM_, double *RESTRICT ARR_);
double GMRFLib_sparse_ddot_sum_(GMRFLib_idxval_tp * RESTRICT ELM_, double *RESTRICT ARR_);
double GMRFLib_sparse_dsum(int n, double *RESTRICT a, int *RESTRICT idx);
int GMRFLib_isum(int n, int *ix);

//
static FORCEINLINE double GMRFLib_sparse_dsum_INLINE(int n, double *RESTRICT a, int *RESTRICT idx)
{
	double res = 0.0;
#pragma omp simd reduction(+: res)
	for (int i = 0; i < n; i++) {
		res += a[idx[i]];
	}
	return res;
}

#define SPARSE_DOT()					\
	double res = 0.0;				\
	_Pragma("omp simd reduction(+:res)")		\
	for (int i = 0; i < n; i++) {			\
		res += v[i] * a[idx[i]];		\
	}						\
	return res

static FORCEINLINE double GMRFLib_sparse_ddot_INLINE(int n, double *RESTRICT v, double *RESTRICT a, int *RESTRICT idx)
{
	// sum_i v[i] * a[idx[i]]
#if defined(INLA_WITH_MKL)
	if (n > 256) {
		double cblas_ddoti(const int nz, const double *x, const int *indx, const double *y);
		return cblas_ddoti(n, v, idx, a);
	}
#endif
	SPARSE_DOT();
}

__END_DECLS
#endif
