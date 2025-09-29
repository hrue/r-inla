
#ifndef __GMRFLib_DOT_H__
#define __GMRFLib_DOT_H__

#include <stdlib.h>
#include <stddef.h>
#include <math.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS					       /* empty */
#define __END_DECLS					       /* empty */
#endif

__BEGIN_DECLS
#include "GMRFLib/GMRFLibP.h"
#define GMRFLib_DOT_GROUP_NLIM 256
double GMRFLib_ddot(int n, double *x, double *y);
double GMRFLib_ddot_idx(int n, double *v, double *a, int *idx);
double GMRFLib_ddot_idx_mkl(int n, double *v, double *a, int *idx);
double GMRFLib_dot_product_group_mkl_opt(GMRFLib_idxval_tp * ELM_, double *ARR_);
double GMRFLib_dot_product_group_prefetch(GMRFLib_idxval_tp * __restrict ELM_, double *__restrict ARR_);
double GMRFLib_dot_product_group_sparse_opt(GMRFLib_idxval_tp * ELM_, double *ARR_);
double GMRFLib_dot_product_opt(GMRFLib_idxval_tp * __restrict ELM_, double *__restrict ARR_);
double GMRFLib_dot_product_sparse_mkl(GMRFLib_idxval_tp * ELM_, double *ARR_);
double GMRFLib_dot_product_sparse_opt(GMRFLib_idxval_tp * __restrict ELM_, double *__restrict ARR_);
void GMRFLib_chose_threshold_ddot(void);

double GMRFLib_ddot_idx_avx2(int n, double *__restrict v, double *__restrict a, int *__restrict idx);
double GMRFLib_ddot_idx_opt(int n, double *__restrict v, double *__restrict a, int *__restrict idx);
double GMRFLib_ddot_opt(int n, double *__restrict x, double *__restrict y);
double GMRFLib_dsum_idx_opt(int n, double *__restrict a, int *__restrict idx);

#define GMRFLib_dot_product_INLINE(ans_, v_, a_) ans_ = GMRFLib_dot_product_opt(v_, a_)
#define GMRFLib_dot_product_INLINE_ADDTO(ans_, v_, a_)	ans_ += GMRFLib_dot_product_opt(v_, a_)

#if defined(INLA_WITH_ARMPL)
#include "armpl_sparse.h"
double GMRFLib_dot_product_sparse_armpl(GMRFLib_idxval_tp * __restrict ELM_, double *__restrict ARR_);
#endif

__END_DECLS
#endif
