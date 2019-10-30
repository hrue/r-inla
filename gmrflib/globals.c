
/* globals.c
 * 
 * Copyright (C) 2001-2006 Havard Rue
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

#ifndef HGVERSION
#define HGVERSION
#endif
//static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

/* Pre-hg-Id: $Id: globals.c,v 1.43 2010/02/15 08:26:37 hrue Exp $ */

#define __GMRFLib_DONT_DEFINE_GLOBALS
#include <limits.h>
#include <time.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
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

/*! 
  \brief Set the blas level in the Lapack routines using the band solver.

  When \c GMRFLib_smtp is set to #GMRFLib_SMTP_BAND, then there is an additional option of chosing
  the blas level for the routines used.  The two options are
  - level 2 (#BLAS_LEVEL2),  and
  - level 3 (#BLAS_LEVEL3)
  
  Use level 2 for smaller problems and level 3 for larger problems. Default is #BLAS_LEVEL3.

  \remark For larger problems you may want to use the TAUCS library instead.
 */
int GMRFLib_blas_level = BLAS_LEVEL3;

/*!
  \brief Use internal lookup-tables or not.

  GMRFLib use in the file wa.c internal lookup-tables to speed up the calculations (default).  This
  feature can be turned off by setting \c GMRFLib_use_wa_table_lookup() to \c #GMRFLib_FALSE.\n\n
*/
int GMRFLib_use_wa_table_lookup = GMRFLib_TRUE;

/*! 
  \brief Verify graph after reading.

  GMRFLib can automatically verify a graph after reading it from file. By default this feature is
  off, but can be turned on setting \c GMRFLib_verify_graph_read_from_disc() to \c #GMRFLib_TRUE \n\n
*/
int GMRFLib_verify_graph_read_from_disc = GMRFLib_FALSE;

/*!
  \brief Collect timing statistics

  GMRFLib collects by default various internal statistics of the CPU demanding functions.  This
  statistics can be extracted using the functions \c GMRFLib_timer_full_report() or
  \c GMRFLib_timer_report(). This feature can be turned off by setting 
  \c GMRFLib_collect_timer_statistics()
  to \c #GMRFLib_FALSE. \n\n
*/
int GMRFLib_collect_timer_statistics = GMRFLib_TRUE;

/*!
  \brief The CPU timing function

  Writing a function that (correctly) returns the used CPU time from a fixed reference, is often OS dependent. You may override
  the GMRFLib's implementation by letting \c GMRFLib_cpu() point to your function. \n\n
*/
GMRFLib_cpu_tp *GMRFLib_cpu = GMRFLib_cpu_default;

/* 
   default uniform(0,1) generator and support utilities. these routines are defined in random.c. make the rng_ptr threadprivate,
   so each thread has its own copy.
 */
gsl_rng *GMRFLib_rng_ptr = NULL;			       /* this holds the RNG and its state and is avail globally */
#pragma omp threadprivate(GMRFLib_rng_ptr)

unsigned long int GMRFLib_rng_seed;			       /* this holds a copy of the last seed */
#pragma omp threadprivate(GMRFLib_rng_seed)

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

/* 
   use faster integration than GSL?
 */
int GMRFLib_faster_integration = GMRFLib_TRUE;

/* 
   number of subdivisions of the faster integration
 */
int GMRFLib_faster_integration_np = 80;

/* 
   OpenMP spesifics
 */
int GMRFLib_thread_id = 0;
#pragma omp threadprivate(GMRFLib_thread_id)


/*
  Signal USR2: Stop optimiser and present results
*/
int GMRFLib_request_optimiser_to_stop = GMRFLib_FALSE;

/* 
   Maximum dimension of bitmap-files, or unlimited if <= 0
 */
int GMRFLib_bitmap_max_dimension = -1;

/* 
   Swap bitmap-file ?
 */
int GMRFLib_bitmap_swap = 0;

/* 
   Holds the thread strategy
 */
GMRFLib_openmp_tp *GMRFLib_openmp = NULL;


/* 
   Holds the MemInfo flag
 */
int GMRFLib_meminfo_thread_id = 0;
#pragma omp threadprivate(GMRFLib_meminfo_thread_id)

/* 
   INLA catch error...
 */
int GMRFLib_catch_error_for_inla = GMRFLib_FALSE;


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
   internal use only; for debugging
 */
int GMRFLib_debug_code = 0;

/* 
   tell the pardiso-interface that we're in a thread-safe area
 */
int GMRFLib_pardiso_thread_safe = 1;

/* 
   tell if we have a working pardiso library, -1, is for 'not checked yet'
 */
int GMRFLib_pardiso_ok = -1;
