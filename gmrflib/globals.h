
/* globals.h
 * 
 * Copyright (C) 2001-2024 Havard Rue
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

/*!
  \file globals.h
  \brief Typedefs and defines for \ref globals.c

  This file contains definitions of global variables and other typedefs.

  The random number generator (RNG) used in GMRFLib, need to support the following requirements.

  - A function of type \c GMRFLib_uniform_tp() used to return a Uniform(0,1) variable.
  - A function of type \c GMRFLib_uniform_init_tp() to initialize the RNG
  - A function of type \c GMRFLib_uniform_getstate_tp() which returns in a malloc'ed array the state of the RNG
  - A function of type \c GMRFLib_uniform_setstate_tp() which set the state in the RNG.

  GMRFLib's implementation of this RNG is available in the file \c random.c

*/

#ifndef __GMRFLib_GLOBALS_H__
#define __GMRFLib_GLOBALS_H__

#include <stdlib.h>

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
#if !defined(GMRFLib_FALSE)
#define GMRFLib_FALSE (0)
#endif
#if !defined(GMRFLib_TRUE)
#define GMRFLib_TRUE  (1)
#endif
typedef double GMRFLib_uniform_tp(void);
typedef int GMRFLib_uniform_init_tp(unsigned long int seed);
typedef void *GMRFLib_uniform_getstate_tp(size_t *siz);
typedef int GMRFLib_uniform_setstate_tp(void *state);

typedef double *GMRFLib_ai_INLA_userfunc0_tp(int thread_id, GMRFLib_problem_tp * problem, double *theta, int nhyper);
typedef double *GMRFLib_ai_INLA_userfunc1_tp(int thread_id, double *theta, int nhyper, double *covmat);
typedef double *GMRFLib_ai_INLA_userfunc2_tp(int number, double *theta, int nhyper, double *covmat, void *arg);
typedef double *GMRFLib_ai_INLA_userfunc3_tp(int number, double *theta, int nhyper, double *covmat, void *arg);

#ifndef __GMRFLib_DONT_DEFINE_GLOBALS

extern char GMRFLib_path[];

extern GMRFLib_smtp_tp GMRFLib_smtp;
extern GMRFLib_reorder_tp GMRFLib_reorder;

extern gsl_rng *GMRFLib_rng_ptr;
#pragma omp threadprivate(GMRFLib_rng_ptr)

extern GMRFLib_uniform_tp *GMRFLib_uniform;
extern GMRFLib_uniform_init_tp *GMRFLib_uniform_init;
extern GMRFLib_uniform_getstate_tp *GMRFLib_uniform_getstate;
extern GMRFLib_uniform_setstate_tp *GMRFLib_uniform_setstate;

extern GMRFLib_ai_INLA_userfunc0_tp *GMRFLib_ai_INLA_userfunc0;	/* points to the function */
extern int GMRFLib_ai_INLA_userfunc0_dim;		       /* dimension of func() */
extern GMRFLib_density_tp **GMRFLib_ai_INLA_userfunc0_density; /* return the marginal densities here */

extern GMRFLib_ai_INLA_userfunc1_tp *GMRFLib_ai_INLA_userfunc1;	/* points to the function */
extern int GMRFLib_ai_INLA_userfunc1_dim;		       /* dimension of func() */
extern GMRFLib_density_tp **GMRFLib_ai_INLA_userfunc1_density; /* return the marginal densities here */

extern GMRFLib_ai_INLA_userfunc2_tp **GMRFLib_ai_INLA_userfunc2;
extern void **GMRFLib_ai_INLA_userfunc2_args;
extern GMRFLib_density_tp ***GMRFLib_ai_INLA_userfunc2_density;
extern int GMRFLib_ai_INLA_userfunc2_n;
extern int *GMRFLib_ai_INLA_userfunc2_len;
extern char **GMRFLib_ai_INLA_userfunc2_tag;

extern GMRFLib_ai_INLA_userfunc3_tp **GMRFLib_ai_INLA_userfunc3;
extern void **GMRFLib_ai_INLA_userfunc3_args;
extern GMRFLib_density_tp ***GMRFLib_ai_INLA_userfunc3_density;
extern int GMRFLib_ai_INLA_userfunc3_n;
extern int *GMRFLib_ai_INLA_userfunc3_len;
extern char **GMRFLib_ai_INLA_userfunc3_tag;

extern int GMRFLib_bitmap_max_dimension;
extern int GMRFLib_bitmap_swap;
extern GMRFLib_openmp_tp *GMRFLib_openmp;
extern GMRFLib_global_node_tp GMRFLib_global_node;
extern GMRFLib_density_storage_strategy_tp GMRFLib_density_storage_strategy;
extern int GMRFLib_pardiso_ok;

extern int GMRFLib_faster_constr;
extern double GMRFLib_aqat_m_diag_add;
extern int GMRFLib_inla_mode;
extern int GMRFLib_Qx_strategy;				       // 0 = serial, 1 = parallel
extern int GMRFLib_preopt_predictor_strategy;		       // 0 = !data_rich, 1 = data_rich
extern double GMRFLib_weight_prob;
extern double GMRFLib_weight_prob_one;
extern double **GMRFLib_dot_product_optim_report;
extern double GMRFLib_dot_product_gain;			       // set to < 0 to disable it
extern int GMRFLib_internal_opt;
extern int GMRFLib_save_memory;

extern int GMRFLib_sort2_id_cut_off;
extern int GMRFLib_sort2_dd_cut_off;

extern int GMRFLib_write_state;
extern int GMRFLib_gaussian_data;

extern int GMRFLib_testit_mode;
extern int GMRFLib_testit_debug;
extern int GMRFLib_taucs_sort_L;
extern int GMRFLib_opt_solve;

extern int GMRFLib_intern_flag;
extern int GMRFLib_cachelinesize;

extern char *GMRFLib_tmpdir;
#endif
__END_DECLS
#endif
