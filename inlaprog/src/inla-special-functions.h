#ifndef __INLA_SPECIAL_FUNCTIONS_H__
#       define __INLA_SPECIAL_FUNCTIONS_H__
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
 *
 */

double inla_gamma_fast1(double x);
double inla_gamma_fast2(double x);
double inla_lgamma_fast1(double x);
double inla_lgamma_fast2(double x);

#define LGAMMAfn(x_) inla_lgamma_fast2(x_)
#define GAMMAfn(x_) exp(LGAMMAfn(x_))

__END_DECLS
#endif
