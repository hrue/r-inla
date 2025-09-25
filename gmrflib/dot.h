
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
double GMRFLib_dot_product_group_serial_opt(GMRFLib_idxval_tp * ELM_, double *ARR_);
double GMRFLib_dot_product_optimized(GMRFLib_idxval_tp * __restrict ELM_, double *__restrict ARR_);
double GMRFLib_dot_product_serial_mkl(GMRFLib_idxval_tp * ELM_, double *ARR_);
double GMRFLib_dot_product_serial_optimized(GMRFLib_idxval_tp * __restrict ELM_, double *__restrict ARR_);
void GMRFLib_chose_threshold_ddot(void);

double GMRFLib_ddot_idx_avx2(int n, double *__restrict v, double *__restrict a, int *__restrict idx);
double GMRFLib_ddot_idx_optimized(int n, double *__restrict v, double *__restrict a, int *__restrict idx);
double GMRFLib_ddot_optimized(int n, double *__restrict x, double *__restrict y);
double GMRFLib_dsum_idx_optimized(int n, double *__restrict a, int *__restrict idx);
double GMRFLib_dsum_optimized(int n, double *__restrict a);

#define NOT_IN_USE____GMRFLib_dot_product_INLINE(ans_, v_, a_)		\
	if (v_->n >= 8L) {						\
		ans_ = GMRFLib_dot_product_optimized(v_, a_);		\
	} else {							\
		double *_v = v_->val;					\
		int *_idx = v_->idx;					\
									\
		switch(v_->n) {						\
		case 0:	ans_ = 0.0; break;				\
		case 1:	ans_ = _v[0] * a_[_idx[0]]; break;		\
		case 2:	ans_ = _v[0] * a_[_idx[0]] + _v[1] * a_[_idx[1]]; break; \
		case 3:	ans_ = _v[0] * a_[_idx[0]] + _v[1] * a_[_idx[1]] + _v[2] * a_[_idx[2]]; break; \
		case 4:	ans_ = _v[0] * a_[_idx[0]] + _v[1] * a_[_idx[1]] + _v[2] * a_[_idx[2]] + _v[3] * a_[_idx[3]]; break; \
		case 5:	ans_ = _v[0] * a_[_idx[0]] + _v[1] * a_[_idx[1]] + _v[2] * a_[_idx[2]] + _v[3] * a_[_idx[3]] + _v[4] * a_[_idx[4]]; break; \
		case 6:	ans_ = _v[0] * a_[_idx[0]] + _v[1] * a_[_idx[1]] + _v[2] * a_[_idx[2]] + _v[3] * a_[_idx[3]] + _v[4] * a_[_idx[4]] + _v[5] * a_[_idx[5]]; break; \
		case 7:	ans_ = _v[0] * a_[_idx[0]] + _v[1] * a_[_idx[1]] + _v[2] * a_[_idx[2]] + _v[3] * a_[_idx[3]] + _v[4] * a_[_idx[4]] + _v[5] * a_[_idx[5]] + _v[6] * a_[_idx[6]];	\
		}							\
	}

#define NOT_IN_USE____GMRFLib_dot_product_INLINE_ADDTO(ans_, v_, a_)	\
	if (v_->n >= 8L) {						\
		ans_ += GMRFLib_dot_product_optimized(v_, a_);		\
	} else {							\
		double *_v = v_->val;					\
		int *_idx = v_->idx;					\
									\
		switch(v_->n) {						\
		case 0:	break;						\
		case 1:	ans_ += _v[0] * a_[_idx[0]]; break;		\
		case 2:	ans_ += _v[0] * a_[_idx[0]] + _v[1] * a_[_idx[1]]; break; \
		case 3:	ans_ += _v[0] * a_[_idx[0]] + _v[1] * a_[_idx[1]] + _v[2] * a_[_idx[2]]; break; \
		case 4:	ans_ += _v[0] * a_[_idx[0]] + _v[1] * a_[_idx[1]] + _v[2] * a_[_idx[2]] + _v[3] * a_[_idx[3]]; break; \
		case 5:	ans_ += _v[0] * a_[_idx[0]] + _v[1] * a_[_idx[1]] + _v[2] * a_[_idx[2]] + _v[3] * a_[_idx[3]] + _v[4] * a_[_idx[4]]; break; \
		case 6:	ans_ += _v[0] * a_[_idx[0]] + _v[1] * a_[_idx[1]] + _v[2] * a_[_idx[2]] + _v[3] * a_[_idx[3]] + _v[4] * a_[_idx[4]] + _v[5] * a_[_idx[5]]; break; \
		case 7:	ans_ += _v[0] * a_[_idx[0]] + _v[1] * a_[_idx[1]] + _v[2] * a_[_idx[2]] + _v[3] * a_[_idx[3]] + _v[4] * a_[_idx[4]] + _v[5] * a_[_idx[5]] + _v[6] * a_[_idx[6]]; \
		}							\
	}


#define GMRFLib_dot_product_INLINE(ans_, v_, a_) ans_ = GMRFLib_dot_product_optimized(v_, a_)
#define GMRFLib_dot_product_INLINE_ADDTO(ans_, v_, a_)	ans_ += GMRFLib_dot_product_optimized(v_, a_)

#if defined(INLA_WITH_ARMPL)
#include "armpl_sparse.h"
double GMRFLib_dot_product_serial_armpl(GMRFLib_idxval_tp * __restrict ELM_, double *__restrict ARR_);
#endif

__END_DECLS
#endif
