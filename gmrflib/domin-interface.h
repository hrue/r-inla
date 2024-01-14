
/* domin-interface.h
 * 
 * Copyright (C) 2006-2024 Havard Rue
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

#ifndef __GMRFLib_DOMIN_INTERFACE_H__
#define __GMRFLib_DOMIN_INTERFACE_H__

#include <stdlib.h>
#include <stddef.h>
#include <math.h>

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
   
*/
    typedef struct {
	gsl_matrix *A;
	gsl_matrix *tAinv;
	int reset_directions;
	int thread_id;
} opt_dir_params_tp;

typedef struct {
	int *f_count;
	int nhyper;
	int use_directions;
	void **Qfunc_arg;
	void *log_extra_arg;
	void *loglFunc_arg;
	char *compute;
	double ***hyperparam;
	double *b;
	double *c;
	double *d;
	double *mean;
	double *solution;
	double *x;
	double fvalue;
	GMRFLib_Qfunc_tp **Qfunc;
	GMRFLib_ai_log_extra_tp *log_extra;
	GMRFLib_ai_param_tp *ai_par;
	GMRFLib_ai_store_tp *ai_store;
	GMRFLib_bfunc_tp **bfunc;
	GMRFLib_constr_tp *constr;
	GMRFLib_graph_tp *graph;
	GMRFLib_logl_tp *loglFunc;
	gsl_matrix *directions;
	GMRFLib_preopt_tp *preopt;
	int parallel_linesearch;
	GMRFLib_idx_tp *d_idx;
} GMRFLib_opt_arg_tp;

int GMRFLib_opt_setup(double ***hyperparam, int nhyper,
		      GMRFLib_ai_log_extra_tp * log_extra, void *log_extra_arg,
		      char *compute,
		      double *x, double *b, double *c, double *mean,
		      GMRFLib_bfunc_tp ** bfunc,
		      double *d, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg,
		      GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg,
		      GMRFLib_constr_tp * constr, GMRFLib_ai_param_tp * ai_par, GMRFLib_ai_store_tp * ai_store,
		      GMRFLib_preopt_tp * preopt, GMRFLib_idx_tp * d_idx);
int GMRFLib_opt_exit(void);
int GMRFLib_opt_f_intern(int thread_id, double *x, double *fx, int *ierr, GMRFLib_ai_store_tp * ais, GMRFLib_tabulate_Qfunc_tp ** tabQfunc,
			 double **bnew);
int GMRFLib_opt_f(int thread_id, double *x, double *fx, int *ierr, GMRFLib_tabulate_Qfunc_tp ** tabQfunc, double **bnew);
int GMRFLib_opt_f_omp(double **x, int nx, double *f, int *ierr);
int GMRFLib_opt_gradf(double *x, double *gradx, int *ierr);
int GMRFLib_opt_estimate_hessian(double *hessian, double *x, double *log_dens_mode, int count);
int GMRFLib_opt_get_f_count(void);
int GMRFLib_opt_gradf_intern(double *x, double *gradx, double *f0, int *ierr);
int GMRFLib_opt_get_hyper(double *x);
int GMRFLib_opt_get_latent(double *latent);
int GMRFLib_opt_set_hyper(double *x);
int GMRFLib_opt_set_latent(double *latent);

GMRFLib_matrix_tp *GMRFLib_opt_get_directions(void);
double GMRFLib_opt_get_f(void);
double GMRFLib_gsl_f(const gsl_vector * v, void *params);
int GMRFLib_gsl_get_results(double *theta_mode, double *log_dens_mode);
int GMRFLib_opt_get_latent(double *latent);
int GMRFLib_opt_reset_directions(void);
int GMRFLib_gsl_optimize(GMRFLib_ai_param_tp * ai_par);
int GMRFLib_opt_dir_step(double *x, int idx, double h);
int GMRFLib_opt_dir_transform_gradient(double *grad);
int GMRFLib_opt_dir_transform_hessian(double *hessian);
int GMRFLib_opt_turn_off_parallel_linesearch();
void GMRFLib_gsl_df(const gsl_vector * v, void *params, gsl_vector * df);
void GMRFLib_gsl_fdf(const gsl_vector * x, void *params, double *f, gsl_vector * df);

GSL_VAR const gsl_multimin_fdfminimizer_type *gsl_multimin_fdfminimizer_vector_bfgs3;	/* my version of vector_bfgs2() */
GSL_VAR const gsl_multimin_fdfminimizer_type *gsl_multimin_fdfminimizer_vector_bfgs4;	/* my 2nd version of vector_bfgs2() */

void GMRFLib_opt_trace_append(GMRFLib_opt_trace_tp ** otrace, double f, double *theta, int nfunc);
void GMRFLib_opt_trace_free(GMRFLib_opt_trace_tp * otrace);
GMRFLib_opt_trace_tp *GMRFLib_opt_trace_get(void);
void inla_write_state_to_file(double fval, int nfun, int ntheta, double *theta, int nx, double *x);
__END_DECLS
#endif
