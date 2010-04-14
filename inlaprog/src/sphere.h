
/* sphere.h
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
 * RCSId: $Id: sphere.h,v 1.12 2010/02/09 13:16:37 hrue Exp $
 *
 */
#ifndef __INLA_SPHERE_H__
#define __INLA_SPHERE_H__
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

/* 
   
 */
    typedef enum {
	SPH_BASIS_GENERAL = 1,
	SPH_BASIS_ROTSYM = 2
} sph_basis_enum_tp;

typedef struct {
	sph_basis_enum_tp type;
	int order;
	int nparam;
} sph_basis_model_tp;

typedef int sph_basis_n_func_tp(int order);
typedef int sph_basis_eval_func_tp(sph_basis_model_tp * model, double *s, double *res_array);


typedef struct {
	int n;
	double **s;
} inla_points_tp;

typedef struct {
	int order;
	int ntheta;
	double ***theta;				       /* theta[i][thread_id][0] */
	double **theta_extra;				       /* theta_extra[thread_id][i] */

	int n;
	double **basis;

	sph_basis_n_func_tp *basis_func_n;
	sph_basis_eval_func_tp *basis_func_eval;
} inla_theta_tp;


typedef struct {
	int n;
	inla_points_tp *s;

	double *T;
	double *K;
	double *C;

	GMRFLib_tabulate_Qfunc_tp *G;
	GMRFLib_graph_tp *G_graph;

	GMRFLib_tabulate_Qfunc_tp *G2;
	GMRFLib_graph_tp *G2_graph;

	inla_theta_tp *Kmodel;
	inla_theta_tp *Tmodel;

	double **oc;					       /* the hyperparameter: ocillating coeff */

	GMRFLib_Qfunc_tp *Qfunc;
	void *Qfunc_arg;
	GMRFLib_graph_tp *graph;
} inla_sphere_tp;


double inla_KT_model_eval(inla_theta_tp * theta_model, int idx);
double inla_sphere_Qfunction(int node, int nnode, void *arg);
double *inla_sphere_userfunc0(GMRFLib_problem_tp * problem, double *theta, int *nhyper);
int inla_sphere_userfunc1(double *theta, int nhyper, double *covmat);
inla_points_tp *inla_alloc_points(int n);
int inla_KT_model_eval2(double *value0, double *value1, inla_theta_tp * theta_model, int idx, int iidx);
int inla_KT_model_init(inla_theta_tp * theta_model, sph_basis_model_tp * bmodel, inla_points_tp * points);
int inla_free_points(inla_points_tp * p);
int inla_read_diagonal_matrix(const char *filename, int *n, double **x);
int inla_read_points(const char *filename, inla_points_tp ** points);
int inla_sph_basis_eval(int kmax, double *s, double *res_array);
int inla_sph_basis_n(int kmax);
int inla_sphere_build_model(inla_sphere_tp ** smodel, const char *dir, sph_basis_model_tp * basisT, sph_basis_model_tp * basisK);
int sph_basis_eval_general(sph_basis_model_tp * model, double *s, double *res_array);
int sph_basis_eval_rotsym(sph_basis_model_tp * model, double *s, double *res_array);
int sph_basis_n_general(int kmax);
int sph_basis_n_rotsym(int kmax);
int sph_basis_n(sph_basis_model_tp * model);


__END_DECLS
#endif
