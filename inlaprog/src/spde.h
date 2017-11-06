
/* spde.h
 * 
 * Copyright (C) 2007-2009 Havard Rue
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
 *
 */
#ifndef __INLA_SPDE_H__
#define __INLA_SPDE_H__
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
#include "GMRFLib/GMRFLib.h"

/* 
   
 */
    typedef struct {
	int n;
	int dim;
	double **s;					       /* n x dim */
} inla_spde_points_tp;

typedef struct {
	int ntheta;
	double ***theta;				       /* theta[i][thread_id][0] */
	double **theta_extra;				       /* theta_extra[thread_id][i] */

	int n;
	double **basis;
} inla_spde_theta_tp;


typedef struct {
	int n;
	inla_spde_points_tp *s;

	double *T;
	double *K;
	double *C;

	GMRFLib_tabulate_Qfunc_tp *G1;
	GMRFLib_graph_tp *G1_graph;

	GMRFLib_tabulate_Qfunc_tp *G2;
	GMRFLib_graph_tp *G2_graph;

	inla_spde_theta_tp *Tmodel;
	inla_spde_theta_tp *Kmodel;

	double **oc;					       /* the hyperparameter: ocillating coeff */

	GMRFLib_Qfunc_tp *Qfunc;
	void *Qfunc_arg;
	GMRFLib_graph_tp *graph;
} inla_spde_tp;


double inla_spde_KT_model_eval(inla_spde_theta_tp * theta_model, int idx);
double inla_spde_Qfunction(int node, int nnode, void *arg);
double *inla_spde_userfunc0(GMRFLib_problem_tp * problem, double *theta, int *nhyper);
int inla_spde_userfunc1(double *theta, int nhyper, double *covmat);
int inla_spde_KT_model_eval2(double *value0, double *value1, inla_spde_theta_tp * theta_model, int idx, int iidx);
int inla_spde_KT_model_init(inla_spde_theta_tp * theta_model, GMRFLib_matrix_tp * basis);
int inla_spde_free_points(inla_spde_points_tp * p);
int inla_spde_basis_eval(int kmax, double *s, double *res_array);
int inla_spde_basis_n(int kmax);
int inla_spde_build_model(inla_spde_tp ** smodel, const char *prefix);
inla_spde_points_tp *inla_spde_set_points(GMRFLib_matrix_tp * M);

__END_DECLS
#endif
