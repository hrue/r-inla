#ifndef __INLA_EVAL_H__
#       define __INLA_EVAL_H__
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
// 
double inla_eval(char *expression, double *x, double *theta, int ntheta);
double inla_eval_expression(char *expression, double *x, double *theta, int ntheta);
double inla_eval_table(char *expression, double *x, double *theta, int ntheta);

__END_DECLS
#endif
