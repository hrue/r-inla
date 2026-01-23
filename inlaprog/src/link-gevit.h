#ifndef __INLA_LINK_GEVIT_H__
#       define __INLA_LINK_GEVIT_H__

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
double inla_log_pgev(double y, double xi, double *l_xi);
double inla_log_pcgev(double y, double xi, double *l_xi);
double inla_pgev(double y, double xi, double *l_xi);
double inla_pcgev(double y, double xi, double *l_xi);
double inla_inv_pgev(double p, double xi, double *l_xi);
double inla_inv_pcgev(double p, double xi, double *l_xi);
double link_gev(int thread_id, double arg, map_arg_tp typ, void *param, double *cov);
double link_cgev(int thread_id, double arg, map_arg_tp typ, void *param, double *cov);
double link_gev_core(int thread_id, double arg, map_arg_tp typ, void *param, int type);
double link_gev_bound(double xi, double *l_xi);
void link_gev_test(double xi, double intercept);


__END_DECLS
#endif
