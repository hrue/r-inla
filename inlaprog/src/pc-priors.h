
/* pc-priors.h
 * 
 * Copyright (C) 2014 Havard Rue
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * The author's contact information:
 *
 *        Haavard Rue
 *        CEMSE Division
 *        King Abdullah University of Science and Technology
 *        Thuwal 23955-6900, Saudi Arabia
 *        Email: haavard.rue@kaust.edu.sa
 *        Office: +966 (0)12 808 0640
 *
 */
#ifndef __INLA_PC_PRIORS_H__
#define __INLA_PC_PRIORS_H__
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
double inla_pc_sn_d(double alpha, double *deriv);
double inla_pcp_dof_d(double dof);
double inla_pcp_dof_dof(double d);
double inla_pcp_dof_kld_approx(double dof);


__END_DECLS
#endif
