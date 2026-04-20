#ifndef __INLA_PC_POWERLINK_H__
#       define __INLA_PC_POWERLINK_H__

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
    typedef struct {
	GMRFLib_spline_tp *cdf, *icdf;
	double power, xmin, xmax, pmin, pmax, mean, sd;
} inla_powerlink_table_tp;

double map_inv_powerlink_core(double arg, map_arg_tp typ, void *param, double *intercept);

__END_DECLS
#endif
