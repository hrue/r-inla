
/* openmp.h
 * 
 * Copyright (C) 2007-2024 Havard Rue
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
  \file openmp.h
  \brief Typedefs used to do approximative inference
*/

#ifndef __GMRFLib_OPENMP_H__
#define __GMRFLib_OPENMP_H__

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
#include <omp.h>
    typedef enum {
	GMRFLib_OPENMP_STRATEGY_SMALL = 1,
	GMRFLib_OPENMP_STRATEGY_MEDIUM,
	GMRFLib_OPENMP_STRATEGY_LARGE,
	GMRFLib_OPENMP_STRATEGY_HUGE,
	GMRFLib_OPENMP_STRATEGY_PARDISO,
	GMRFLib_OPENMP_STRATEGY_DEFAULT,
	GMRFLib_OPENMP_STRATEGY_NONE
} GMRFLib_openmp_strategy_tp;

#define GMRFLib_OPENMP_STRATEGY_NAME(num)				\
	((num) == GMRFLib_OPENMP_STRATEGY_SMALL ? "small" :		\
	 ((num) == GMRFLib_OPENMP_STRATEGY_MEDIUM ? "medium" :		\
	  ((num) == GMRFLib_OPENMP_STRATEGY_LARGE ? "large" :		\
	   ((num) == GMRFLib_OPENMP_STRATEGY_HUGE ? "huge" :		\
	    ((num) == GMRFLib_OPENMP_STRATEGY_PARDISO ? "pardiso" :	\
	     ((num) == GMRFLib_OPENMP_STRATEGY_DEFAULT ? "default" :	\
	      ((num) == GMRFLib_OPENMP_STRATEGY_NONE ? "none" : "THIS SHOULD NOT HAPPEN")))))))

typedef enum {
	GMRFLib_OPENMP_PLACES_PARSE_MODEL = 1,
	GMRFLib_OPENMP_PLACES_BUILD_MODEL,
	GMRFLib_OPENMP_PLACES_OPTIMIZE,
	GMRFLib_OPENMP_PLACES_HESSIAN,
	GMRFLib_OPENMP_PLACES_HESSIAN_SCALE,
	GMRFLib_OPENMP_PLACES_INTEGRATE_HYPERPAR,
	GMRFLib_OPENMP_PLACES_COMBINE,
	GMRFLib_OPENMP_PLACES_EXTERNAL,
	GMRFLib_OPENMP_PLACES_TIMING,
	GMRFLib_OPENMP_PLACES_GCPO_BUILD,
	GMRFLib_OPENMP_PLACES_SPECIAL,
	GMRFLib_OPENMP_PLACES_DEFAULT,
	GMRFLib_OPENMP_PLACES_NONE
} GMRFLib_openmp_place_tp;

#define GMRFLib_OPENMP_PLACE_NAME(num)					\
	((num) == GMRFLib_OPENMP_PLACES_PARSE_MODEL ? "parse.model" :	\
	 ((num) == GMRFLib_OPENMP_PLACES_BUILD_MODEL ? "build.model" :	\
	  ((num) == GMRFLib_OPENMP_PLACES_OPTIMIZE ? "optimize" :	\
	   ((num) == GMRFLib_OPENMP_PLACES_HESSIAN ? "hessian" :	\
	    ((num) == GMRFLib_OPENMP_PLACES_HESSIAN_SCALE ? "hessian.scale" : \
	     ((num) == GMRFLib_OPENMP_PLACES_INTEGRATE_HYPERPAR ? "integrate.hyperpar" : \
	      ((num) == GMRFLib_OPENMP_PLACES_COMBINE ? "combine" :	\
	       ((num) == GMRFLib_OPENMP_PLACES_EXTERNAL ? "external" :	\
		((num) == GMRFLib_OPENMP_PLACES_GCPO_BUILD ? "gcpo_build" : \
		 ((num) == GMRFLib_OPENMP_PLACES_SPECIAL ? "special" : \
		  ((num) == GMRFLib_OPENMP_PLACES_DEFAULT ? "default" :	\
		   ((num) == GMRFLib_OPENMP_PLACES_NONE ? "none" : "THIS SHOULD NOT HAPPEN"))))))))))))

typedef struct {
	GMRFLib_openmp_place_tp place;
	GMRFLib_openmp_strategy_tp strategy;
	int max_threads;
	int *max_threads_nested;
	int blas_num_threads;
	// for PARDISO, like _outer is the number of threads in the outer loop, while _inner is the number of threads for
	// pardiso. the _inner is only relevant if nested=1.
	int max_threads_outer;
	int max_threads_inner;
	// when this is TRUE, then do PARDISO is parallel if the function call is serial
	int adaptive;
} GMRFLib_openmp_tp;

#define GMRFLib_MAX_THREADS() (GMRFLib_openmp->max_threads)

// Might replace `4' in the generic pardiso control statement later (if that happens)
#define GMRFLib_PARDISO_MAX_NUM_THREADS() (GMRFLib_openmp->adaptive ?	\
					   IMIN(GMRFLib_MAX_THREADS(), GMRFLib_openmp->max_threads_nested[1] * 8) : \
					   GMRFLib_openmp->max_threads_nested[1])

#define GMRFLib_OPENMP_IN_SERIAL()                  ((omp_get_num_threads() == 1) && (omp_get_level() == 0))
#define GMRFLib_OPENMP_IN_PARALLEL()                (!GMRFLib_OPENMP_IN_SERIAL())
#define GMRFLib_OPENMP_IN_PARALLEL_ONE_THREAD()     ((omp_get_num_threads() == 1) && (omp_get_level() == 1))
#define GMRFLib_OPENMP_IN_PARALLEL_ONEPLUS_THREAD() (omp_in_parallel() == 1)

int GMRFLib_set_blas_num_threads(int threads);
int GMRFLib_openmp_nested_fix(void);
int GMRFLib_openmp_implement_strategy(GMRFLib_openmp_place_tp place, void *arg, GMRFLib_smtp_tp * smtp);
int GMRFLib_openmp_implement_strategy_special(int outer, int inner);

#if defined(INLA_LINK_WITH_MKL)
void MKL_Set_Num_Threads(int);
#endif

__END_DECLS
#endif
