/*!
  \file utils.h
  \brief Typedefs for \ref utils.c
*/

#ifndef __GMRFLib_SN_G_H__
#define __GMRFLib_SN_G_H__

#include <stdlib.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS					       /* empty */
#define __END_DECLS					       /* empty */
#endif

__BEGIN_DECLS int GMRFLib_sn_g_get_order(void);

//
double *GMRFLib_sn_g_get_coof(double skew, double *cx);
double *GMRFLib_sn_ginv_get_coof(double skew, double *cx);
double GMRFLib_sn_g_eval(double x, double *cx);
double GMRFLib_sn_g_eval_deriv(double x, double *cx);
//

__END_DECLS
#endif
