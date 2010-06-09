
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
 *       H{\aa}vard Rue
 *       Department of Mathematical Sciences
 *       The Norwegian University of Science and Technology
 *       N-7491 Trondheim, Norway
 *       Voice: +47-7359-3533    URL  : http://www.math.ntnu.no/~hrue  
 *       Fax  : +47-7359-3524    Email: havard.rue@math.ntnu.no
 *
 *
 */
#ifndef __INLA_SPDE_H__
#define __INLA_SPDE_H__
#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS					       /* empty */
# define __END_DECLS					       /* empty */
#endif
__BEGIN_DECLS

#include "fmesher-io.h"

/* 
   
 */
    typedef enum {
	SPDE_BASIS_GENERAL = 1,
	SPDE_BASIS_ROTSYM = 2
} spde_basis_enum_tp;

typedef struct {
	spde_basis_enum_tp type;
	int order;
	int nparam;
} spde_basis_model_tp;

typedef int spde_basis_n_func_tp(int order);
typedef int spde_basis_eval_func_tp(spde_basis_model_tp * model, double *s, double *res_array);


typedef struct {
	int n;
	int dim;
	double **s;					       /* n x dim */
} inla_spde_points_tp;

typedef struct {
	int order;
	int ntheta;
	double ***theta;				       /* theta[i][thread_id][0] */
	double **theta_extra;				       /* theta_extra[thread_id][i] */

	int n;
	double **basis;

	spde_basis_n_func_tp *basis_func_n;
	spde_basis_eval_func_tp *basis_func_eval;
} inla_spde_theta_tp;


typedef struct {
	int n;
	inla_spde_points_tp *s;

	double *T;
	double *K;
	double *C;

	GMRFLib_tabulate_Qfunc_tp *G;
	GMRFLib_graph_tp *G_graph;

	GMRFLib_tabulate_Qfunc_tp *G2;
	GMRFLib_graph_tp *G2_graph;

	inla_spde_theta_tp *Kmodel;
	inla_spde_theta_tp *Tmodel;

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
int inla_spde_KT_model_init(inla_spde_theta_tp * theta_model, inla_matrix_tp *basis);
int inla_spde_free_points(inla_spde_points_tp * p);
int inla_spde_basis_eval(int kmax, double *s, double *res_array);
int inla_spde_basis_n(int kmax);
int inla_spde_build_model(inla_spde_tp ** smodel, const char *prefix, spde_basis_model_tp * basisT, spde_basis_model_tp * basisK);
int spde_basis_eval_general(spde_basis_model_tp * model, double *s, double *res_array);
int spde_basis_eval_rotsym(spde_basis_model_tp * model, double *s, double *res_array);
int spde_basis_n_general(int kmax);
int spde_basis_n_rotsym(int kmax);
int spde_basis_n(spde_basis_model_tp * model);

inla_spde_points_tp *inla_spde_set_points(inla_matrix_tp *M);

//inla_spde_points_tp *inla_spde_alloc_points(int n);
//int inla_spde_read_diagonal_matrix(const char *filename, int *n, double **x);
//int inla_spde_read_points(const char *filename, inla_spde_points_tp ** points);

__END_DECLS
#endif
