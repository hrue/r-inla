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
void GMRFLib_abs(int, double *RESTRICT, double *RESTRICT);
void GMRFLib_add(int n, double *RESTRICT x, double *RESTRICT y, double *RESTRICT z);
void GMRFLib_cdaddto(int n, double *RESTRICT x, double cx, double *RESTRICT y);
void GMRFLib_daddto(int n, double *RESTRICT x, double *RESTRICT y);
void GMRFLib_exp(int, double *RESTRICT, double *RESTRICT);
void GMRFLib_exp_inc(int n, double *RESTRICT x, int inc, double *RESTRICT y);
void GMRFLib_fm_po_1(int n, double *RESTRICT mask, double *RESTRICT ll, double *RESTRICT wp, double *RESTRICT integral2, double *RESTRICT integral3, double *RESTRICT integral4);
void GMRFLib_log(int, double *RESTRICT, double *RESTRICT);
void GMRFLib_log1p(int, double *RESTRICT, double *RESTRICT);
void GMRFLib_mul(int n, double *RESTRICT x, double *RESTRICT y, double *RESTRICT z);
void GMRFLib_sqr(int n, double *RESTRICT x, double *RESTRICT y);
void GMRFLib_sqrt(int n, double *RESTRICT x, double *RESTRICT y);

__END_DECLS
#endif
