#ifndef __GMRFLib_SIMD_H__
#define __GMRFLib_SIMD_H__

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
#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"
#if defined(INLA_WITH_SIMD)
//#define SLEEF_ENABLE_OMP_SIMD
#include <sleef.h>
#endif
void GMRFLib_exp(int, double *, double *);
void GMRFLib_exp_inc(int n, double *x, int inc, double *y);
void GMRFLib_log(int, double *, double *);
void GMRFLib_log1p(int, double *, double *);
void GMRFLib_sqr(int n, double *x, double *y);
void GMRFLib_sqrt(int n, double *x, double *y);
void GMRFLib_add(int n, double *x, double *y, double *z);
void GMRFLib_mul(int n, double *x, double *y, double *z);
void GMRFLib_daddto(int n, double *x, double *y);
void GMRFLib_cdaddto(int n, double *x, double cx, double *y);

__END_DECLS
#endif
