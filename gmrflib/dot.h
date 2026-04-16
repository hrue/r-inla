
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
double GMRFLib_sparse_ddot(int n, double *__restrict v, double *__restrict a, int *__restrict idx);
double GMRFLib_sparse_ddot_(GMRFLib_idxval_tp * __restrict ELM_, double *__restrict ARR_);
double GMRFLib_sparse_ddot_group_(GMRFLib_idxval_tp * __restrict ELM_, double *__restrict ARR_);
double GMRFLib_sparse_ddot_group_simple_(GMRFLib_idxval_tp * __restrict ELM_, double *__restrict ARR_);

double GMRFLib_sparse_ddot_ddot_(GMRFLib_idxval_tp * __restrict ELM_, double *__restrict ARR_);
double GMRFLib_sparse_ddot_sum_(GMRFLib_idxval_tp * __restrict ELM_, double *__restrict ARR_);
double GMRFLib_sparse_ddot_sum1_(GMRFLib_idxval_tp * __restrict ELM_, double *__restrict ARR_);

__END_DECLS
#endif
