#ifndef __GMRFLib_FM_H__
#       define __GMRFLib_FM_H__

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
#       include "GMRFLib/GMRFLib.h"
void GMRFLib_abs(int, double *restrict, double *restrict);
void GMRFLib_exp(int, double *restrict, double *restrict);
void GMRFLib_exp_inc(int n, double *restrict x, int inc, double *restrict y);
void GMRFLib_log(int, double *restrict, double *restrict);
void GMRFLib_log1p(int, double *restrict, double *restrict);
void GMRFLib_sqr(int n, double *restrict x, double *restrict y);
void GMRFLib_sqrt(int n, double *restrict x, double *restrict y);
void GMRFLib_add(int n, double *restrict x, double *restrict y, double *restrict z);
void GMRFLib_mul(int n, double *restrict x, double *restrict y, double *restrict z);
void GMRFLib_daddto(int n, double *restrict x, double *restrict y);
void GMRFLib_cdaddto(int n, double *restrict x, double cx, double *restrict y);

__END_DECLS
#endif
