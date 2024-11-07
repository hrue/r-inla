
/* globals.c
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
 */

/*!
  \file globals.c
  \brief Set values of global variables.
*/

#define __GMRFLib_DONT_DEFINE_GLOBALS
#include <limits.h>
#include <time.h>
#include <stdlib.h>
#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"
#undef __GMRFLib_DONT_DEFINE_GLOBALS

/*!
  \brief Define the path to various GMRFLib data-files.

  This variable define a colon-separated list of paths to various GMRFLib-files, and can be
  preprendend by an environment variable \c GMRFLib_PATH. The initial value is set to the \c PREFIX
  variable in the installation.
*/
const char GMRFLib_path[] = GMRFLib_PATH;

/*!
  \brief Define the sparse solver for the sparse-matrix computations.

  This variable defines the which method is used to factorise the sparse matrices. The current
  implementation includes
  - #GMRFLib_SMTP_BAND, using the band-matrix routines in \c LAPACK
  - #GMRFLib_SMTP_TAUCS, using the multifrontal supernodal factorisation in the \c TAUCS library.

  and its values are define in GMRFLib_smtp_tp.  Default value is #GMRFLib_SMTP_TAUCS.\n\n
*/
GMRFLib_smtp_tp GMRFLib_smtp = GMRFLib_SMTP_TAUCS;

/*!
  \brief Define the reordering routine for sparse-matrix computations.

  This variable defines the which method is used to reorder the sparse matrices.
  
  The current implementation includes
  - #GMRFLib_REORDER_DEFAULT, use the default choice depending on \c GMRFLib_smtp
  - #GMRFLib_REORDER_IDENTITY, do not use any reordering (the identity reordering)
  - #GMRFLib_REORDER_REVERSE_IDENTITY, the reverse identity reordering
  - #GMRFLib_REORDER_BAND, minimise the band-width
  - #GMRFLib_REORDER_METIS, nested dissection
  - #GMRFLib_REORDER_GENMMD, multiple minimum degree
  - #GMRFLib_REORDER_AMD, approximate minimum degree
  - #GMRFLib_REORDER_AMDC, approximate minimum degree (C-version)
  - #GMRFLib_REORDER_AMDBAR, approximate minimum degree, without agressive absorption
  - #GMRFLib_REORDER_AMDBARC, approximate minimum degree, without agressive absorption (C-version)
  - #GMRFLib_REORDER_MD, true minimum degree
  - #GMRFLib_REORDER_MMD, multiple minimum degree

  and its values are define in GMRFLib_reorder_tp.  The default choice, is #GMRFLib_REORDER_BAND if \c GMRFLib_smtp =
  #GMRFLib_SMTP_BAND and #GMRFLib_REORDER_METIS if  \c GMRFLib_smtp = #GMRFLib_SMTP_TAUCS.\n\n

  \sa The routine GMRFLib_optimise_reorder() chose the best reordering for a given graph.
*/
GMRFLib_reorder_tp GMRFLib_reorder = GMRFLib_REORDER_DEFAULT;

int GMRFLib_verify_graph_read_from_disc = GMRFLib_FALSE;

gsl_rng *GMRFLib_rng_ptr = NULL;			       /* this holds the RNG and its state and is avail globally */
#pragma omp threadprivate(GMRFLib_rng_ptr)

GMRFLib_uniform_tp *GMRFLib_uniform = GMRFLib_rng_uniform;
GMRFLib_uniform_init_tp *GMRFLib_uniform_init = GMRFLib_rng_init;
GMRFLib_uniform_getstate_tp *GMRFLib_uniform_getstate = GMRFLib_rng_getstate;
GMRFLib_uniform_setstate_tp *GMRFLib_uniform_setstate = GMRFLib_rng_setstate;

/* 
   function to compute E_{theta|y}(func(...)) in INLA
 */
GMRFLib_ai_INLA_userfunc0_tp *GMRFLib_ai_INLA_userfunc0 = NULL;	/* function-ptr */
int GMRFLib_ai_INLA_userfunc0_dim = 0;			       /* dimension of func() */
GMRFLib_density_tp **GMRFLib_ai_INLA_userfunc0_density;	       /* returning the marginal densities here */

GMRFLib_ai_INLA_userfunc1_tp *GMRFLib_ai_INLA_userfunc1;       /* points to the function */
GMRFLib_density_tp **GMRFLib_ai_INLA_userfunc1_density;	       /* returning the marginal densities here */
int GMRFLib_ai_INLA_userfunc1_dim = 0;			       /* dimension of func() */

GMRFLib_ai_INLA_userfunc2_tp **GMRFLib_ai_INLA_userfunc2 = NULL;
void **GMRFLib_ai_INLA_userfunc2_args = NULL;
GMRFLib_density_tp ***GMRFLib_ai_INLA_userfunc2_density = NULL;
int GMRFLib_ai_INLA_userfunc2_n = 0;
int *GMRFLib_ai_INLA_userfunc2_len = NULL;
char **GMRFLib_ai_INLA_userfunc2_tag = NULL;

GMRFLib_ai_INLA_userfunc3_tp **GMRFLib_ai_INLA_userfunc3 = NULL;
void **GMRFLib_ai_INLA_userfunc3_args = NULL;
GMRFLib_density_tp ***GMRFLib_ai_INLA_userfunc3_density = NULL;
int GMRFLib_ai_INLA_userfunc3_n = 0;
int *GMRFLib_ai_INLA_userfunc3_len = NULL;
char **GMRFLib_ai_INLA_userfunc3_tag = NULL;


int GMRFLib_bitmap_max_dimension = -1;
int GMRFLib_bitmap_swap = 0;

/* 
   Holds the thread strategy
 */
GMRFLib_openmp_tp *GMRFLib_openmp = NULL;

/* 
   define global nodes = {factor, degree}. factor: a node is defined to be global if nneig(i) >= (n-1) *factor degree :node is define to be global if nneig(i) >=
   degree.
 */
GMRFLib_global_node_tp GMRFLib_global_node = { 2.0, INT_MAX };

/* 
   storage strategy for density
 */
GMRFLib_density_storage_strategy_tp GMRFLib_density_storage_strategy = GMRFLib_DENSITY_STORAGE_STRATEGY_DEFAULT;

/* 
   tell if we have a working pardiso library, -1, is for 'not checked yet'
 */
int GMRFLib_pardiso_ok = -1;

int GMRFLib_faster_constr = 1;

int GMRFLib_inla_mode = 0;

// add stability to AQ^-1A^T
double GMRFLib_aqat_m_diag_add = 0.0;

int GMRFLib_Qx_strategy = 0;				       // 0 = serial, 1 = parallel
int GMRFLib_preopt_predictor_strategy = 0;		       // 0 = !data_rich, 1 = data_rich

double GMRFLib_weight_prob = 0.975;			       // for pruning weights for densities
double GMRFLib_weight_prob_one = 0.999;			       // for pruning weights otherwise
double **GMRFLib_dot_product_optim_report = NULL;
double GMRFLib_dot_product_gain = -1.0;
int GMRFLib_sort2_id_cut_off = 50;
int GMRFLib_sort2_dd_cut_off = 70;
int GMRFLib_internal_opt = 1;
int GMRFLib_save_memory = 0;

int GMRFLib_write_state = 0;
int GMRFLib_gaussian_data = 0;
int GMRFLib_testit_mode = 0;
int GMRFLib_testit_debug = 0;
int GMRFLib_taucs_sort_L = 0;
int GMRFLib_opt_solve = 0;
int GMRFLib_intern_flag = 0;
int GMRFLib_cachelinesize = 64;
char *GMRFLib_tmpdir = NULL;
