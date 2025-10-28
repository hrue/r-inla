#ifndef __GMRFLib_FIT_SN_H__
#define __GMRFLib_FIT_SN_H__
#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS					       /* empty */
#define __END_DECLS					       /* empty */
#endif

#include <math.h>

__BEGIN_DECLS void fitsn_ld(int n, double *x, double *param, double *ld);
void fitsn_gradhess(double x, double *param, double *grad, double *hess);
void fitsn_test(void);
void fitsn_test_grad(void);
void fitsn_test_hess(void);



__END_DECLS
#endif
