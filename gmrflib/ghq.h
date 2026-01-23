#ifndef __GMRFLib_GHQ_H__
#       define __GMRFLib_GHQ_H__

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

__BEGIN_DECLS typedef struct {
	int n;
	double skew3;					       /* skewness^(1/3) */
	double *nodes;
	double *w;
	double *w_grad;
	double *w_hess;
} GMRFLib_snq_tp;

/*
 */

#       define GMRFLib_skew_to_skew3(skew_) (DSIGN(skew_) * pow(ABS(skew_), 1.0/3.0))
#       define GMRFLib_skew3_to_skew(skew3_) POW3(skew3_)

GMRFLib_snq_tp *GMRFLib_snq(int n, double skew3);
int GMRFLib_ghq(double **xp, double **wp, int n);
int GMRFLib_ghq__intern(double *x, double *w, int n);
int GMRFLib_ghq_abscissas(double **xp, int n);
int GMRFLib_ghq_ms(double **xp, double **wp, int n, double mean, double stdev);
int GMRFLib_ghq_weights(double **wp, int n);
int GMRFLib_snq_free(GMRFLib_snq_tp * q);

double inla_logcdf_normal(double x);			       /* external function */
__END_DECLS
#endif
