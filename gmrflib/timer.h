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

#define GMRFLib_ENTER_ROUTINE						\
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
