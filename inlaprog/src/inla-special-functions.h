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
double inla_lbeta(double a, double b);
double inla_beta(double a, double b);
void inla_lgamma_fast2_m(size_t m, double *restrict x, double *restrict res);
void inla_lbeta_m(size_t m, double *restrict a, double *restrict b, double *restrict llbeta);

#       define LGAMMAfn(x_) inla_lgamma_fast2(x_)
#       define GAMMAfn(x_) inla_gamma_fast2(x_)
#       define LBETAfn(a_, b_) inla_lbeta(a_, b_)


__END_DECLS
#endif
