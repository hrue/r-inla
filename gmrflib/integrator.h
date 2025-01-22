#ifndef __GMRFLib_INTEGRATOR_H__
#define __GMRFLib_INTEGRATOR_H__

#include <math.h>
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

__BEGIN_DECLS

/* 
   dummy comment
 */
typedef double (*integrand)(unsigned ndim, const double *x, void *);

/* Integrate the function f from xmin[dim] to xmax[dim], with at most maxEval function evaluations (0 for no limit), until the
   given absolute or relative error is achieved.  val returns the integral, and err returns the estimate for the absolute error
   in val.  The return value of the function is 0 on success and non-zero if there was an error. */

int adapt_integrate(integrand f, void *fdata, unsigned dim, const double *xmin, const double *xmax,
		    unsigned maxEval, double reqAbsError, double reqRelError, double *val, double *err);

__END_DECLS
#endif
