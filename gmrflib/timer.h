
/* timer.h
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
  \file timer.h
  \brief Typedefs and defines for the timer-utility in GMRFLib.
*/

#ifndef __GMRFLib_TIMER_H__
#define __GMRFLib_TIMER_H__

#include <time.h>
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

/* 
   dummy comment
 */
// defined in high-prec-timer.cpp
double GMRFLib_timer_chrono(void);
//#define GMRFLib_timer() omp_get_wtime()
#define GMRFLib_timer() GMRFLib_timer_chrono()

typedef struct {
	char *name;					       /* function name */
	double ntimes;					       /* number of times called */
	double ctime_acc;				       /* accumulated ctime */
	double ctime_ref;				       /* last reference time */

	double ctime_min;				       /* minimum time used */
	double ctime_max;				       /* maximum times used */
	double ctime_acc2;				       /* accumulated ctime^2 */
} GMRFLib_timer_hashval_tp;

#define GMRFLib_ENTER_ROUTINE \
	GMRFLib_DEBUG_INIT();						\
	GMRFLib_TRACE_INIT();						\
	static double trace_cpu_acc_ = 0.0;				\
	_Pragma("omp threadprivate(trace_cpu_acc_)")			\
	static double trace_cpu_ = 0.0;					\
	_Pragma("omp threadprivate(trace_cpu_)")			\
	trace_cpu_ = GMRFLib_timer();					\
	GMRFLib_TRACE_i("Enter, total", trace_count_);

#define GMRFLib_LEAVE_ROUTINE if (1)					\
	{								\
		trace_cpu_acc_ += (GMRFLib_timer() - trace_cpu_);		\
		GMRFLib_TRACE_idd("Leave, count cpu/count*1E6 total", trace_count_, 1.0E6 * trace_cpu_acc_ / (double) trace_count_, trace_cpu_acc_); \
	}

double GMRFLib_timer_windows(void);

__END_DECLS
#endif
