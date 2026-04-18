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

__BEGIN_DECLS

/*
 */
int GMRFLib_ghq(double **xp, double **wp, int n);
int GMRFLib_ghq__intern(double *x, double *w, int n);
int GMRFLib_ghq_abscissas(double **xp, int n);
int GMRFLib_ghq_ms(double **xp, double **wp, int n, double mean, double stdev);
int GMRFLib_ghq_weights(double **wp, int n);

__END_DECLS
#endif
