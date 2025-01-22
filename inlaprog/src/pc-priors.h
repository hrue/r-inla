#ifndef __INLA_PC_PRIORS_H__
#define __INLA_PC_PRIORS_H__

#define inla_pc_sn_alpha2skew(_skew) inla_pc_sn_core(-1, _skew)
#define inla_pc_sn_skew2alpha(_alpha) inla_pc_sn_core(1, _alpha)

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS					       /* empty */
#define __END_DECLS					       /* empty */
#endif
__BEGIN_DECLS

/* 
 *
 */
    GMRFLib_spline_tp * inla_pc_sn_create_spline(void);
GMRFLib_spline_tp *inla_pcp_dof_create_spline(void);
double inla_pc_h_default(double x, int inverse, int derivative);
double inla_pc_simplex_core_d(double *x, int p, double lambda);
double inla_pc_simplex_d(double *x, double *b, int p, double lambda);
double inla_pc_sn_core(int code, double arg);
double inla_pc_sn_d(double skew, double *deriv);
double inla_pcp_dof_d(double dof);
double inla_pcp_dof_dof(double d);
double inla_pcp_dof_kld_approx(double dof);

__END_DECLS
#endif
