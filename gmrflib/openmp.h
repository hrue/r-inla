
/* openmp.h
 * 
 * Copyright (C) 2007 Havard Rue
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
#ifdef _OPENMP
#include <omp.h>
#else
extern void omp_set_num_threads(int);
extern int omp_get_num_threads(void);
extern int omp_get_max_threads(void);
extern int omp_get_thread_num(void);
extern int omp_get_thread_num_(void);
extern int omp_get_num_procs(void);
extern int omp_in_parallel(void);
extern void omp_set_dynamic(int);
extern int omp_get_dynamic(void);
extern void omp_set_nested(int);
extern int omp_get_nested(void);
extern double omp_get_wtime(void);
extern double omp_get_wtick(void);
#endif


typedef enum {
	GMRFLib_OPENMP_STRATEGY_SMALL = 1,
	GMRFLib_OPENMP_STRATEGY_MEDIUM,
	GMRFLib_OPENMP_STRATEGY_LARGE,
	GMRFLib_OPENMP_STRATEGY_HUGE,
	GMRFLib_OPENMP_STRATEGY_PARDISO_SERIAL,
	GMRFLib_OPENMP_STRATEGY_PARDISO_PARALLEL,
	GMRFLib_OPENMP_STRATEGY_DEFAULT, 
	GMRFLib_OPENMP_STRATEGY_NONE
} GMRFLib_openmp_strategy_tp;

#define GMRFLib_OPENMP_STRATEGY_NAME(num) \
	((num) == GMRFLib_OPENMP_STRATEGY_SMALL ? "small" :		\
	 ((num) == GMRFLib_OPENMP_STRATEGY_MEDIUM ? "medium" :	\
	  ((num) == GMRFLib_OPENMP_STRATEGY_LARGE ? "large" :		\
	   ((num) == GMRFLib_OPENMP_STRATEGY_HUGE ? "huge" :		\
	    ((num) == GMRFLib_OPENMP_STRATEGY_PARDISO_SERIAL ? "pardiso.serial" : \
	     ((num) == GMRFLib_OPENMP_STRATEGY_PARDISO_PARALLEL ? "pardiso.parallel" : \
	      ((num) == GMRFLib_OPENMP_STRATEGY_DEFAULT ? "default" :	\
	       ((num) == GMRFLib_OPENMP_STRATEGY_NONE ? "none" : "THIS SHOULD NOT HAPPEN"))))))))

typedef enum {
	GMRFLib_OPENMP_PLACES_BUILD_MODEL = 1,
	GMRFLib_OPENMP_PLACES_OPTIMIZE,
	GMRFLib_OPENMP_PLACES_HESSIAN,
	GMRFLib_OPENMP_PLACES_HESSIAN_SCALE,
	GMRFLib_OPENMP_PLACES_INTEGRATE,
	GMRFLib_OPENMP_PLACES_INTEGRATE_HYPERPAR,
	GMRFLib_OPENMP_PLACES_COMBINE,
	GMRFLib_OPENMP_PLACES_EXTERNAL,
	GMRFLib_OPENMP_PLACES_DEFAULT, 
	GMRFLib_OPENMP_PLACES_NONE
} GMRFLib_openmp_place_tp;

#define GMRFLib_OPENMP_PLACE_NAME(num) \
	((num) == GMRFLib_OPENMP_PLACES_BUILD_MODEL ? "build.model" :	\
	 ((num) == GMRFLib_OPENMP_PLACES_OPTIMIZE ? "optimize" :	\
	  ((num) == GMRFLib_OPENMP_PLACES_HESSIAN ? "hessian" :		\
	   ((num) == GMRFLib_OPENMP_PLACES_HESSIAN_SCALE ? "hessian.scale" : \
	    ((num) == GMRFLib_OPENMP_PLACES_INTEGRATE ? "integrate" :	\
	     ((num) == GMRFLib_OPENMP_PLACES_INTEGRATE_HYPERPAR ? "integrate.hyperpar" : \
	      ((num) == GMRFLib_OPENMP_PLACES_COMBINE ? "combine" : \
	       ((num) == GMRFLib_OPENMP_PLACES_EXTERNAL ? "external" : \
		((num) == GMRFLib_OPENMP_PLACES_DEFAULT ? "default" : \
		 ((num) == GMRFLib_OPENMP_PLACES_NONE ? "none" : "THIS SHOULD NOT HAPPEN"))))))))))

typedef struct {
	GMRFLib_openmp_place_tp place;
	int max_threads;
	// for PARDISO, like _outer is the number of threads in the outer loop, while _inner is the number of threads for
	// pardiso. the _inner is only relevant if nested=1.
	int max_threads_outer;
	int max_threads_inner;
	GMRFLib_openmp_strategy_tp strategy;
} GMRFLib_openmp_tp;

#define GMRFLib_MAX_THREADS (GMRFLib_openmp ? GMRFLib_openmp->max_threads : IMIN(omp_get_max_threads(), omp_get_num_procs()))

int GMRFLib_openmp_implement_strategy(GMRFLib_openmp_place_tp place, void *arg, GMRFLib_smtp_tp *smtp);

__END_DECLS
#endif
