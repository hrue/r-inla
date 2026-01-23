#ifndef __INLA_QUANTILE_REGRESSION_H__
#       define __INLA_QUANTILE_REGRESSION_H__
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
#       include "GMRFLib/density.h"
    struct inla_qgamma_cache_tp {
	double quantile;
	GMRFLib_spline_tp *s;
};

GMRFLib_spline_tp **inla_qcontpois_func(double alpha, int num);
double inla_pcontpois(double y, double lambda);
double inla_pcontpois_deriv(double y, double lambda);
double inla_qcontpois(double quantile, double alpha, double *initial_guess);
double inla_qcontpois_eta(double quantile, double alpha, double *initial_guess);
double inla_qgamma_cache(double shape, double quantile);

__END_DECLS
#endif
