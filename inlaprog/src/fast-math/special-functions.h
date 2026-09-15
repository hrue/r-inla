#include "GMRFLib/GMRFLib.h"
#ifndef __SPECIAL_FUNCTIONS_H__
#       define __SPECIAL_FUNCTIONS_H__
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

double inla_logcdf_normal(double x);
double inla_cdf_normal(double x);
double inla_cdf_normal_fast(double x);
double inla_logitcdf_normal(double x);
double inla_logcdf_normal_fast(double x);
double inla_lgamma_fast(double x);
double inla_lbeta(double a, double b);
void inla_lbeta_m(size_t m, double *RESTRICT a, double *RESTRICT b, double *RESTRICT llbeta);
void inla_lgamma_fast_m(size_t m, double *RESTRICT x, double *RESTRICT res);
void inla_lgamma_m(size_t m, double *RESTRICT x, double *RESTRICT res);

// we chose here if to use the faster approximation to 'lgamma()', that is slightly less accurate
#       if 1
#              define LGAMMAfn(x_) inla_lgamma_fast(x_)
#              define LGAMMAfn_m(m_, x_, r_) inla_lgamma_fast_m((size_t) (m_), x_, r_)
#              define GAMMAfn(x_) inla_gamma_fast(x_)
#       else
#              define LGAMMAfn(x_) inla_lgamma(x_)
#              define LGAMMAfn_m(m_, x_, r_) inla_lgamma_m((size_t) (m_), x_, r_)
#              define GAMMAfn(x_) inla_gamma(x_)
#       endif

#       define LBETAfn(a_, b_) inla_lbeta(a_, b_)

__END_DECLS
#endif
