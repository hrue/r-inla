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
double inla_lgamma_fast(double x);
double inla_lbeta(double a, double b);
void inla_lbeta_m(size_t m, double *RESTRICT a, double *RESTRICT b, double *RESTRICT llbeta);
void inla_lgamma_fast_m(size_t m, double *RESTRICT x, double *RESTRICT res);
void inla_lgamma_m(size_t m, double *RESTRICT x, double *RESTRICT res);

void inla_llike_nbinomial_1(const int m, const double d1, const double d2, const double d3,
			    const double d4, const double d5, double *RESTRICT x, double *RESTRICT logll);
void inla_llike_nbinomial_2(const int m, const double d1, const double d2, const double d3, const double d4, 
			    double *RESTRICT x, double *RESTRICT logll);
void inla_llike_log1p_exp(const int m, double *RESTRICT x, double *RESTRICT y);
void inla_llike_log1p_exp_1(const int m, const double c1, const double c2, double *RESTRICT x, double *RESTRICT y);

void inla_llike_tweedie2_1(const int n, const double p1, const double p2, const double phi, const double y, const double ly, const double sum_w,
			   const double w_max, double *RESTRICT mu, double *RESTRICT ldens);

// define 'INLA_WITHOUT_FAST_LGAMMA' to use the libm-versions: lgamma() instead of the faster approximation that is slightly less
// accurate
#       if defined(INLA_WITHOUT_FAST_LGAMMA)
#              define LGAMMAfn(x_) inla_lgamma(x_)
#              define LGAMMAfn_m(m_, x_, r_) inla_lgamma_m((size_t) (m_), x_, r_)
#              define GAMMAfn(x_) inla_gamma(x_)
#       else
#              define LGAMMAfn(x_) inla_lgamma_fast(x_)
#              define LGAMMAfn_m(m_, x_, r_) inla_lgamma_fast_m((size_t) (m_), x_, r_)
#              define GAMMAfn(x_) inla_gamma_fast(x_)
#       endif

#       define LBETAfn(a_, b_) inla_lbeta(a_, b_)

__END_DECLS
#endif
