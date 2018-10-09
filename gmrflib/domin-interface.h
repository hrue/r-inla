
/* domin-interface.h
 * 
 * Copyright (C) 2006-2008 Havard Rue
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
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif

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
	double ***hyperparam;
	int nhyper;
	GMRFLib_ai_log_extra_tp *log_extra;
	void *log_extra_arg;

	int *f_count;
	char *compute;

	double *solution;
	double fvalue;

	double *x;
	double *b;
	double *c;
	double *mean;
	GMRFLib_bfunc_tp **bfunc;
	double *d;
	GMRFLib_logl_tp *loglFunc;
	void *loglFunc_arg;
	char *fixed_value;
	GMRFLib_graph_tp *graph;
	GMRFLib_Qfunc_tp **Qfunc;
	void **Qfunc_arg;
	GMRFLib_constr_tp *constr;
	GMRFLib_ai_param_tp *ai_par;
	GMRFLib_ai_store_tp *ai_store;
} GMRFLib_domin_arg_tp;

int GMRFLib_domin_setup(double ***hyperparam, int nhyper,
			GMRFLib_ai_log_extra_tp * log_extra, void *log_extra_arg,
			char *compute,
			double *x, double *b, double *c, double *mean,
			GMRFLib_bfunc_tp **bfunc, 
			double *d, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, char *fixed_value,
			GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg,
			GMRFLib_constr_tp * constr, GMRFLib_ai_param_tp * ai_par, GMRFLib_ai_store_tp * ai_store);
int GMRFLib_domin_exit(void);
int GMRFLib_domin_f_intern(double *x, double *fx, int *ierr, GMRFLib_ai_store_tp * ais, GMRFLib_tabulate_Qfunc_tp **tabQfunc, double ** bnew);
int GMRFLib_domin_f(double *x, double *fx, int *ierr, GMRFLib_tabulate_Qfunc_tp **tabQfunc, double ** bnew);
int GMRFLib_domin_f_omp(double **x, int nx, double *f, int *ierr);
int GMRFLib_domin_gradf(double *x, double *gradx, int *ierr);
int GMRFLib_domin_gradf_OLD(double *x, double *gradx, int *ierr);
int GMRFLib_domin_estimate_hessian(double *hessian, double *x, double *log_dens_mode, int count);
int GMRFLib_domin_estimate_hessian_OLD(double *hessian, double *x);
int GMRFLib_domin_get_f_count(void);
int GMRFLib_domin_gradf_intern(double *x, double *gradx, double *f0, int *ierr);

double GMRFLib_gsl_f(const gsl_vector * v, void *params);
void GMRFLib_gsl_df(const gsl_vector * v, void *params, gsl_vector * df);
void GMRFLib_gsl_fdf(const gsl_vector * x, void *params, double *f, gsl_vector * df);
int GMRFLib_gsl_optimize(GMRFLib_ai_param_tp * ai_par);
int GMRFLib_gsl_get_results(double *theta_mode, double *log_dens_mode);

int GMRFLib_test_something____omp(void);

GSL_VAR const gsl_multimin_fdfminimizer_type *gsl_multimin_fdfminimizer_vector_bfgs3;	/* my version of vector_bfgs2() */

__END_DECLS
#endif
