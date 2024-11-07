
/* GMRFLibP.h
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
  \file GMRFLibP.h
  \brief Internal include-file for the GMRFLib source.
*/

#ifndef __GMRFLibP_H__
#define __GMRFLibP_H__

#include <assert.h>
#include <stddef.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>

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
#if __GNUC__ > 7
typedef size_t fortran_charlen_t;
#else
typedef int fortran_charlen_t;
#endif
#define F_ONE ((fortran_charlen_t)1)

// see https://stackoverflow.com/questions/3599160/how-to-suppress-unused-parameter-warnings-in-c
#ifdef __GNUC__
#define UNUSED(x) UNUSED_ ## x __attribute__((__unused__))
#else
#define UNUSED(x) UNUSED_ ## x
#endif

#if defined(NDEBUG)
#error The code assume that NDEBUG is *NOT* defined
#endif

#ifdef __GNUC__
#define POSSIBLY_UNUSED_FUNCTION(x) __attribute__((__unused__)) x
#else
#define POSSIBLY_UNUSED_FUNCTION(x) x
#endif

#pragma omp declare simd
static double POSSIBLY_UNUSED_FUNCTION(SQR) (double x) {
	return (x * x);
}

#pragma omp declare simd
static int POSSIBLY_UNUSED_FUNCTION(ISQR) (int ix) {
	return (ix * ix);
}

#pragma omp declare simd
static double POSSIBLY_UNUSED_FUNCTION(POW3) (double x) {
	return (x * x * x);
}

#pragma omp declare simd
static int POSSIBLY_UNUSED_FUNCTION(IPOW3) (int ix) {
	return (ix * ix * ix);
}

#pragma omp declare simd
static double POSSIBLY_UNUSED_FUNCTION(POW4) (double x) {
	double xx = x * x;
	return (xx * xx);
}

#pragma omp declare simd
static int POSSIBLY_UNUSED_FUNCTION(IPOW4) (int ix) {
	int ixx = ix * ix;
	return (ixx * ixx);
}

//#define DMAX(a,b) GSL_MAX_DBL(a, b)
//#define DMIN(a,b) GSL_MIN_DBL(a, b)
//#define IMAX(a,b) GSL_MAX_INT(a, b)
//#define IMIN(a,b) GSL_MIN_INT(a, b)

#if 1
#define DMAX(a_, b_) fmax(a_, b_)
#define DMIN(a_, b_) fmin(a_, b_)
#else							       // if 1
#pragma omp declare simd
static double POSSIBLY_UNUSED_FUNCTION(DMAX) (double a, double b) {
	return ((a) > (b) ? (a) : (b));
}

#pragma omp declare simd
static double POSSIBLY_UNUSED_FUNCTION(DMIN) (double a, double b) {
	return ((a) < (b) ? (a) : (b));
}
#endif							       // if 1

#pragma omp declare simd
static int POSSIBLY_UNUSED_FUNCTION(IMAX) (int a, int b) {
	return ((a) > (b) ? (a) : (b));
}

#pragma omp declare simd
static int POSSIBLY_UNUSED_FUNCTION(IMIN) (int a, int b) {
	return ((a) < (b) ? (a) : (b));
}

#pragma omp declare simd
static int POSSIBLY_UNUSED_FUNCTION(ITRUNCATE) (int x, int low, int high) {
	// #define ITRUNCATE(x, low, high) IMIN(IMAX(x, low), high)
	return IMIN(IMAX(x, low), high);
}

#pragma omp declare simd
static double POSSIBLY_UNUSED_FUNCTION(TRUNCATE) (double x, double low, double high) {
	// #define TRUNCATE(x, low, high) DMIN( DMAX(x, low), high) 
	return DMIN(DMAX(x, low), high);
}

#define GMRFLib_CACHELINESIZE (GMRFLib_cachelinesize)
#define GMRFLib_CACHELINESIZE_ND (GMRFLib_CACHELINESIZE / sizeof(double))
#define GMRFLib_CACHELINESIZE_NI (GMRFLib_CACHELINESIZE / sizeof(int))
#define GMRFLib_MEM_ALIGN (32L)

typedef enum {
	GMRFLib_MODE_CLASSIC = 1,
	GMRFLib_MODE_COMPACT
} GRMFLib_preopt_mode_tp;

#define GMRFLib_MODE_NAME() (GMRFLib_inla_mode == GMRFLib_MODE_CLASSIC ? "Classic" : \
			     (GMRFLib_inla_mode == GMRFLib_MODE_COMPACT ? "Compact" : "(UNKNOWN MODE)"))

// utility functions for this are mostly in smtp-pardiso.c
typedef struct {
	int n;
	int na;
	int *ia;
	int *ia1;
	int *ja;
	int *ja1;
	int *iwork;
	unsigned char *sha;
} GMRFLib_csr_skeleton_tp;

typedef struct {
	GMRFLib_csr_skeleton_tp *s;
	double *a;
	int copy_only;
} GMRFLib_csr_tp;

typedef struct {
	double *x;
	int free;
} GMRFLib_vec_tp;

typedef enum {
	INLA_B_STRATEGY_SKIP = 0,
	INLA_B_STRATEGY_KEEP = 1
} inla_b_strategy_tp;

/* 
   here are the wrappers for calling functions which return the error-code if it fails
*/
#define GMRFLib_EWRAP__intern(func_call, leave)		\
	if (1) {					\
		int rretval;				\
		rretval = func_call;			\
		if (rretval != GMRFLib_SUCCESS){	\
			GMRFLib_ERROR(rretval);			\
			return rretval;				\
		}						\
	}
#define GMRFLib_EWRAP_MAPKIT__intern(func_call, leave)		\
	if (1){							\
		int rrretval;						\
		rrretval = func_call;					\
		if (rrretval != MAPKIT_OK){				\
			char *msg;					\
			GMRFLib_EWRAP__intern(GMRFLib_sprintf(&msg, "Mapkit-library returned error-code [%1d]", rrretval), leave); \
			GMRFLib_ERROR_MSG(GMRFLib_EMAPKIT, msg);	\
			Free(msg);					\
			return GMRFLib_EMAPKIT;				\
		}							\
	}
#define GMRFLib_EWRAP_GSL__intern(func_call, leave)			\
	if (1){								\
		int rrretval;						\
		gsl_error_handler_t *ehandler = gsl_set_error_handler_off(); \
		rrretval = func_call;					\
		gsl_set_error_handler(ehandler);			\
		if (rrretval != GSL_SUCCESS){				\
			char *msg;					\
			GMRFLib_EWRAP__intern(GMRFLib_sprintf(&msg, "GSL-library returned error-code [%1d]", rrretval), leave); \
			GMRFLib_ERROR_MSG(GMRFLib_EGSL, msg);		\
			Free(msg);					\
			return GMRFLib_EGSL;				\
		}							\
	}
#define GMRFLib_EWRAP_GSL_PTR__intern(func_call, leave)			\
	if (1){								\
		gsl_error_handler_t *ehandler = gsl_set_error_handler_off(); \
		void *retval_ptr = (void *)(func_call);			\
		gsl_set_error_handler(ehandler);			\
		if (retval_ptr == NULL){				\
			char *msg;					\
			GMRFLib_EWRAP__intern(GMRFLib_sprintf(&msg, "GSL-library call returned NULL-pointer"), leave); \
			GMRFLib_ERROR_MSG(GMRFLib_EMEMORY, msg);	\
			Free(msg);					\
			return GMRFLib_EMEMORY;				\
		}							\
	}
#define GMRFLib_EWRAP0(func_call) GMRFLib_EWRAP__intern(func_call, 0)
#define GMRFLib_EWRAP1(func_call) GMRFLib_EWRAP__intern(func_call, 1)
#define GMRFLib_EWRAP0_MAPKIT(func_call) GMRFLib_EWRAP_MAPKIT__intern(func_call, 0)
#define GMRFLib_EWRAP1_MAPKIT(func_call) GMRFLib_EWRAP_MAPKIT__intern(func_call, 1)
#define GMRFLib_EWRAP0_GSL(func_call) GMRFLib_EWRAP_GSL__intern(func_call, 0)
#define GMRFLib_EWRAP1_GSL(func_call) GMRFLib_EWRAP_GSL__intern(func_call, 1)
#define GMRFLib_EWRAP0_GSL_PTR(func_call) GMRFLib_EWRAP_GSL_PTR__intern(func_call, 0)
#define GMRFLib_EWRAP1_GSL_PTR(func_call) GMRFLib_EWRAP_GSL_PTR__intern(func_call, 1)

/* 
   this simply measure the cpu spend in `expressions' and print statistics every occation. there is also a BEGIN/END variant for
   larger blocks.
*/
#define GMRFLib_MEASURE_CPU(msg, expression)		\
	if (1) {					\
		static double _tacc = 0.0;		\
		static int _ntimes = 0;			\
		double _tref;				\
		_tref = GMRFLib_timer();			\
		_ntimes++;				\
		expression;						\
		_tacc += GMRFLib_timer() - _tref;				\
		printf("%s:%s:%d: cpu accumulative [%s] %.6f mean %.8f n %d\n", \
		       __FILE__, __GMRFLib_FuncName, __LINE__, msg, _tacc, _tacc/_ntimes, _ntimes); \
	}

#define GMRFLib_MEASURE_CPU_BEGIN()			\
	if (1) {					\
		static double _tacc = 0.0;		\
		static int _ntimes = 0;			\
		double _tref;				\
		_tref = GMRFLib_timer();			\
		_ntimes++;

#define GMRFLib_MEASURE_CPU_END(msg)					\
	_tacc += GMRFLib_timer() - _tref;					\
	printf("%s:%s:%d: cpu accumulative [%s] %.6f mean %.8f n %d\n",	\
	       __FILE__, __GMRFLib_FuncName, __LINE__, msg, _tacc, _tacc/_ntimes, _ntimes); \
	}

#define GMRFLib_DEBUG_INIT() static int debug_ = -1;			\
	static int debug_count_ = 0;					\
	_Pragma("omp threadprivate(debug_count_)")			\
	debug_count_++;							\
	if (debug_ < 0)	{						\
		debug_ = GMRFLib_debug_functions(__GMRFLib_FuncName);	\
	}

#define GMRFLib_TRACE_INIT() static int trace_ = -1;			\
	static int trace_count_ = 0;					\
	_Pragma("omp threadprivate(trace_count_)")			\
	trace_count_++;							\
	if (trace_ < 0)	{						\
		trace_ = GMRFLib_trace_functions(__GMRFLib_FuncName);	\
	}

#define GMRFLib_DEBUG_IF_TRUE() (debug_)
#define GMRFLib_DEBUG_IF()      (debug_ > 0 && !((debug_count_ - 1) % debug_))

#define GMRFLib_DEBUG(msg_)						\
	if (debug_ && !((debug_count_ - 1) % debug_)) {			\
		printf("\t[%1d] %s (%s:%1d): %s\n", omp_get_thread_num(), GMRFLib_function_name_strip(__GMRFLib_FuncName), __FILE__, __LINE__, msg_); \
	}

#define GMRFLib_DEBUG_i(msg_, i_)					\
	if (debug_ && !((debug_count_ - 1) % debug_)) {			\
		printf("\t[%1d] %s (%s:%1d): %s %d\n", omp_get_thread_num(), GMRFLib_function_name_strip(__GMRFLib_FuncName), __FILE__, __LINE__, msg_, i_); \
	}

#define GMRFLib_DEBUG_ii(msg_, i_, ii_)					\
	if (debug_ && !((debug_count_ - 1) % debug_)) {			\
		printf("\t[%1d] %s (%s:%1d): %s %d %d\n", omp_get_thread_num(), GMRFLib_function_name_strip(__GMRFLib_FuncName), __FILE__, __LINE__, msg_, i_, ii_); \
	}

#define GMRFLib_DEBUG_iii(msg_, i_, ii_, iii_)				\
	if (debug_ && !((debug_count_ - 1) % debug_)) {			\
		printf("\t[%1d] %s (%s:%1d): %s %d %d %d\n", omp_get_thread_num(), GMRFLib_function_name_strip(__GMRFLib_FuncName), __FILE__, __LINE__, msg_, i_, ii_, iii_); \
	}

#define GMRFLib_DEBUG_i_iv(msg_, i_, ii_, iii_, iv_)			\
	if (debug_ && !((debug_count_ - 1) % debug_)) {			\
		printf("\t[%1d] %s (%s:%1d): %s %d %d %d %d\n", omp_get_thread_num(), GMRFLib_function_name_strip(__GMRFLib_FuncName), __FILE__, __LINE__, msg_, i_, ii_, iii_, iv_); \
	}

#define GMRFLib_DEBUG_i_v(msg_, i_, ii_, iii_, iv_, v_)			\
	if (debug_ && !((debug_count_ - 1) % debug_)) {			\
		printf("\t[%1d] %s (%s:%1d): %s %d %d %d %d %d\n", omp_get_thread_num(), GMRFLib_function_name_strip(__GMRFLib_FuncName), __FILE__, __LINE__, msg_, i_, ii_, iii_, iv_, v_); \
	}

#define GMRFLib_DEBUG_d(msg_, d_)					\
	if (debug_ && !((debug_count_ - 1) % debug_)) {			\
		printf("\t[%1d] %s (%s:%1d): %s %.4f\n", omp_get_thread_num(), GMRFLib_function_name_strip(__GMRFLib_FuncName), __FILE__, __LINE__, msg_, d_); \
	}

#define GMRFLib_DEBUG_dd(msg_, d_, dd_)					\
	if (debug_ && !((debug_count_ - 1) % debug_)) {			\
		printf("\t[%1d] %s (%s:%1d): %s %.4f %.4f\n", omp_get_thread_num(), GMRFLib_function_name_strip(__GMRFLib_FuncName), __FILE__, __LINE__, msg_, d_, dd_); \
	}

#define GMRFLib_DEBUG_ddd(msg_, d_, dd_, ddd_)				\
	if (debug_ && !((debug_count_ - 1) % debug_)) {			\
		printf("\t[%1d] %s (%s:%1d): %s %.4f %.4f %.4f\n", omp_get_thread_num(), GMRFLib_function_name_strip(__GMRFLib_FuncName), __FILE__, __LINE__, msg_, d_, dd_, ddd_); \
	}

#define GMRFLib_DEBUG_id(msg_, i_, d_)					\
	if (debug_ && !((debug_count_ - 1) % debug_)) {			\
		printf("\t[%1d] %s (%s:%1d): %s %d %.4f\n", omp_get_thread_num(), GMRFLib_function_name_strip(__GMRFLib_FuncName), __FILE__, __LINE__, msg_, i_, d_); \
	}

#define GMRFLib_DEBUG_idd(msg_, i_, d_, dd_)				\
	if (debug_ && !((debug_count_ - 1) % debug_)) {			\
		printf("\t[%1d] %s (%s:%1d): %s %d %.4f %.4f\n", omp_get_thread_num(), GMRFLib_function_name_strip(__GMRFLib_FuncName), __FILE__, __LINE__, msg_, i_, d_, dd_); \
	}

#define GMRFLib_DEBUG_iddd(msg_, i_, d_, dd_, ddd_)			\
	if (debug_ && !((debug_count_ - 1) % debug_)) {			\
		printf("\t[%1d] %s (%s:%1d): %s %d %.4f %.4f %.4f\n", omp_get_thread_num(), GMRFLib_function_name_strip(__GMRFLib_FuncName), __FILE__, __LINE__, msg_, i_, d_, dd_, ddd_); \
	}

#define GMRFLib_DEBUG_idddd(msg_, i_, d_, dd_, ddd_, dddd_)			\
	if (debug_ && !((debug_count_ - 1) % debug_)) {			\
		printf("\t[%1d] %s (%s:%1d): %s %d %.4f %.4f %.4f %.4f\n", omp_get_thread_num(), GMRFLib_function_name_strip(__GMRFLib_FuncName), __FILE__, __LINE__, msg_, i_, d_, dd_, ddd_, dddd_); \
	}

#define GMRFLib_TRACE_i(msg_, i_)					\
	if (trace_ && !((trace_count_ - 1) % trace_)) {			\
		printf("\t[%1d] %s (%s:%1d): %s %d\n", omp_get_thread_num(), GMRFLib_function_name_strip(__GMRFLib_FuncName), __FILE__, __LINE__, msg_, i_); \
	}

#define GMRFLib_TRACE_idd(msg_, i_, d_, dd_)				\
	if (trace_ && !((trace_count_ - 1) % trace_)) {			\
		printf("\t[%1d] %s (%s:%1d): %s %d %.4f %.4f\n", omp_get_thread_num(), GMRFLib_function_name_strip(__GMRFLib_FuncName), __FILE__, __LINE__, msg_, i_, d_, dd_); \
	}

#define Calloc_init(n_, m_)						\
	size_t calloc_m_ = (m_);					\
	size_t calloc_l1_cacheline_ = GMRFLib_CACHELINESIZE_ND;		\
	size_t calloc_mem_align_ = GMRFLib_MEM_ALIGN / sizeof(double);	\
	size_t calloc_len_ = (size_t)((n_) + calloc_m_ * calloc_mem_align_ * calloc_l1_cacheline_); \
	size_t calloc_offset_ = 0;					\
	size_t calloc_m_count_ = 0;					\
	double *calloc_work_ = Calloc(IMAX(1, calloc_len_), double);	\
	assert(calloc_work_)

#define iCalloc_init(n_, m_)						\
	size_t icalloc_m_ = (m_);					\
	size_t icalloc_l1_cacheline_ = GMRFLib_CACHELINESIZE_NI;		\
	size_t icalloc_mem_align_ = GMRFLib_MEM_ALIGN / sizeof(int);	\
	size_t icalloc_len_ = (size_t)((n_) + icalloc_m_ * icalloc_mem_align_ * icalloc_l1_cacheline_); \
	size_t icalloc_offset_ = 0;					\
	size_t icalloc_m_count_ = 0;					\
	int *icalloc_work_ = Calloc(IMAX(1, icalloc_len_), int); \
	assert(icalloc_work_)

#define Calloc_get(_n)							\
	calloc_work_ + calloc_offset_;					\
	calloc_offset_ += GMRFLib_align((size_t)(_n), sizeof(double));	\
	calloc_m_count_++;						\
	Calloc_check()

#define iCalloc_get(_n)							\
	icalloc_work_ + icalloc_offset_;				\
	icalloc_offset_ += GMRFLib_align((size_t)(_n), sizeof(int));	\
	icalloc_m_count_++;						\
	iCalloc_check()

#define Calloc_check()							\
	if (!(calloc_offset_ <= calloc_len_)) { P(calloc_offset_); P(calloc_len_); }; assert(calloc_offset_ <= calloc_len_); \
	if (!(calloc_m_count_ <= calloc_m_)) { P(calloc_m_); P(calloc_m_count_); }; assert(calloc_m_count_ <= calloc_m_)

#define iCalloc_check()							\
	if (!(icalloc_offset_ <= icalloc_len_)) { P(icalloc_offset_); P(icalloc_len_); }; assert(icalloc_offset_ <= icalloc_len_); \
	if (!(icalloc_m_count_ <= icalloc_m_)) { P(icalloc_m_); P(icalloc_m_count_); }; assert(icalloc_m_count_ <= icalloc_m_)

#define Calloc_free()   if (1) { Calloc_check(); Free(calloc_work_);}
#define iCalloc_free()  if (1) { iCalloc_check(); Free(icalloc_work_); }

#define GMRFLib_ALLOC_SAFE_SIZE(n_, type_) ((size_t)(n_) < PTRDIFF_MAX ? (size_t)(n_) : (size_t)1)
#if 0
#define Calloc(n, type)         (type *)GMRFLib_calloc(GMRFLib_ALLOC_SAFE_SIZE(n, type), sizeof(type), __FILE__, __GMRFLib_FuncName, __LINE__)
#define Malloc(n, type)         (type *)GMRFLib_malloc(GMRFLib_ALLOC_SAFE_SIZE((n) * sizeof(type), char), __FILE__, __GMRFLib_FuncName, __LINE__)
#define Realloc(ptr, n, type)   (type *)GMRFLib_realloc((void *)ptr, GMRFLib_ALLOC_SAFE_SIZE((n)*sizeof(type), char), __FILE__, __GMRFLib_FuncName, __LINE__)
#define Free(ptr)               if (ptr) {GMRFLib_free((void *)(ptr), __FILE__, __GMRFLib_FuncName, __LINE__); ptr=NULL;}
#define Memcpy(dest, src, n)    GMRFLib_memcpy(dest, src, n)
#else
#undef  GMRFLib_TRACE_MEMORY
#define Calloc(n, type)         (type *)calloc(GMRFLib_ALLOC_SAFE_SIZE(n, type), sizeof(type))
#define Malloc(n, type)         (type *)malloc(GMRFLib_ALLOC_SAFE_SIZE((n) * sizeof(type), char))
#define Realloc(ptr, n, type)   ((ptr) ? \
				 (type *)realloc((void *)ptr, GMRFLib_ALLOC_SAFE_SIZE((n) * sizeof(type), char)) : \
				 (type *)calloc(GMRFLib_ALLOC_SAFE_SIZE(n, type), sizeof(type)))
#define Free(ptr)               if (ptr) {free((void *)(ptr)); ptr=NULL;}
#define Memcpy(dest, src, n)    memcpy((void *) (dest), (void *) (src), GMRFLib_ALLOC_SAFE_SIZE(n, char))
#endif
#define Memset(dest, value, n)  memset((void *) (dest), (int) (value), (size_t) (n))

#define likely(x)   __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)

#define ABS(x) fabs(x)
#define FIXME( msg) if (1) { printf("\n{%1d}[%s:%1d] %s: FIXME [%s]\n",  omp_get_thread_num(), __FILE__, __LINE__, __GMRFLib_FuncName,(msg?msg:""));	}
#define FIXME1(msg) if (1) { static int first=1; if (first) { first=0; FIXME(msg); }}
#define FIXME1stderr(msg) if (1) { static int first=1; if (first) { first=0; FIXMEstderr(msg); }}
#define FIXMEstderr( msg) if (1) { fprintf(stderr, "\n{%1d}[%s:%1d] %s: FIXME [%s]\n",  omp_get_thread_num(), __FILE__, __LINE__, __GMRFLib_FuncName,(msg?msg:""));	}
#define IABS(x)   abs(x)
#define ISEQUAL(x, y) (gsl_fcmp(x, y, DBL_EPSILON) == 0)
#define ISEQUAL_x(x, y, eps) (gsl_fcmp(x, y, eps) == 0)
#define ISINF(x) isinf(x)
// same definition in rmath.h
#define ISNAN(x) (isnan(x) != 0)
#define ISSMALL(x) (gsl_fcmp(1.0 + (x), 1.0, DBL_EPSILON) == 0)
#define ISSMALL_x(x, eps) (gsl_fcmp(1.0 + (x), 1.0, eps) == 0)
#define ISZERO(x) (((__typeof (x)) (x)) == 0)
#define LEGAL(i, n) ((i) >= 0 && (i) < (n))
#define MOD(i,n)  (((i)+(n))%(n))
#define OVERLAP(p_, pp_, n_) (!(((pp_) + (n_) - 1 <  (p_)) || ((p_) + (n_) - 1 <  (pp_))))
#define P(x)        if (1) { printf("[%s:%1d] " #x " = [ %.16f ]\n",__FILE__, __LINE__,(double)(x)); }
#define P1(x)       if (1) { static int first=1;  if (first) { printf("[%s:%1d] " #x " = [ %.16f ]\n", __FILE__, __LINE__, (double)(x)); first=0; }}
#define P1stderr(x) if (1) { static int first=1;  if (first) { fprintf(stderr, "[%s:%1d] " #x " = [ %.16f ]\n", __FILE__, __LINE__, (double)(x)); first=0; }}
#define PP(msg,pt)  if (1) { fprintf(stdout, "[%s:%1d] %s ptr " #pt " = %p\n", __FILE__, __LINE__, msg, pt); }
#define PPg(msg,pt) if (1) { fprintf(stdout, "[%s:%1d] %s value " #pt " = %g\n", __FILE__, __LINE__, msg, pt); }
#define PPstderr(msg,pt)  if (1) { fprintf(stderr, "[%s:%1d] %s ptr " #pt " = %p\n", __FILE__, __LINE__, msg, pt); }
#define PPstderrg(msg,pt) if (1) { fprintf(stderr, "[%s:%1d] %s value " #pt " = %g\n", __FILE__, __LINE__, msg, *((double *))pt); }
#define Pstderr(x)  if (1) { fprintf(stderr, "[%s:%1d] " #x " = [ %.16f ]\n",__FILE__, __LINE__,(double)(x)); }
#define ISIGN(x) ((x) >= 0 ? 1 : -1)
#define DSIGN(x) ((x) >= 0.0 ? 1.0 : -1.0)
#define SWAP(x_, y_) if (1) { typeof(x_) tmp___ = x_; x_ = y_; y_ = tmp___; }
#define MAKE_ODD(n_) if (GSL_IS_EVEN(n_)) (n_)++
#define PUSH_AWAY(x_) (DMAX(GSL_DBL_EPSILON, ABS(x_)) * DSIGN(x_))
#define GMRFLib_GLOBAL_NODE(n, gptr) ((int) IMIN((n-1)*(gptr ? (gptr)->factor :  GMRFLib_global_node.factor), \
						 (gptr ? (gptr)->degree : GMRFLib_global_node.degree)))

// https://philippegroarke.com/blog/2017/02/19/quicktip-understanding-16-byte-memory-alignment-detection/
#define SIMD_ALIGNED(ptr_) (((intptr_t)(ptr_) & 0xF) == 0)


#define Orig_GMRFLib_STOP_IF_NAN_OR_INF(value, idx, jdx)		\
	if (ISNAN(value) || ISINF(value)) {				\
		if (!nan_error)						\
			fprintf(stdout,					\
				"\n\tFunction: %s(), Line: %1d, Thread: %1d\n\tVariable evaluates to NAN or INF. idx=(%1d,%1d). I will try to fix it...", \
				__GMRFLib_FuncName, __LINE__, omp_get_thread_num(), idx, jdx); \
		nan_error = 1;						\
	}

#define GMRFLib_STOP_IF_NAN_OR_INF(value, idx, jdx)			\
	if (!gsl_finite(value)) {					\
		if (!nan_error)						\
			fprintf(stdout,					\
				"\n\tFunction: %s(), Line: %1d, Thread: %1d\n\tVariable evaluates to NAN or INF. idx=(%1d,%1d). I will try to fix it...", \
				__GMRFLib_FuncName, __LINE__, omp_get_thread_num(), idx, jdx); \
		nan_error = 1;						\
	}

#define GMRFLib_SET_PREC(arg_) (arg_->log_prec_omp ? exp(*(arg_->log_prec_omp[thread_id])) : 1.0)
#define GMRFLib_SET_RANGE(arg_) (arg_->log_range_omp ? exp(*(arg_->log_range_omp[thread_id])) : 1.0)

#define GMRFLib_CACHE_DELAY() GMRFLib_delay_random(25, 50)
// assume _level() <= 2
#define GMRFLib_CACHE_LEN() (GMRFLib_MAX_THREADS() * (GMRFLib_MAX_THREADS() + 1))
#define GMRFLib_CACHE_SET_ID(__id)					\
	{								\
		int level_ = omp_get_level();				\
		int tnum_ = omp_get_thread_num();			\
		if (level_ <= 1) {					\
			__id =  tnum_;					\
		} else if (level_ == 2) {				\
			int level2_ = omp_get_ancestor_thread_num(level_ -1); \
			__id = IMAX(1, 1 + level2_) * GMRFLib_MAX_THREADS() + tnum_; \
		} else {						\
			assert(0 == 1);					\
		}							\
	}

// this use level1 only. set __id to -1 if we're on level2
#define GMRFLib_CACHE_LEN_LEVEL1_ONLY() (GMRFLib_MAX_THREADS())
#define GMRFLib_CACHE_SET_ID_LEVEL1_ONLY(__id)				\
	__id = (omp_get_level() <= 1 ? omp_get_thread_num() : -1)

// len_work_ * n_work_ >0 will create n_work_ workspaces for all threads, each of (len_work_ * n_work_) doubles. _PTR(i_) will return the ptr to
// the thread spesific workspace index i_ and _ZERO will zero-set it, i_=0,,,n_work_-1. CODE_BLOCK_THREAD_ID must be used to set

#define CODE_BLOCK_WORK_PTR(i_work_) (work__[(nt__ == 1 ? 0 : omp_get_thread_num())] + (size_t) (i_work_) * len_work__)
#define CODE_BLOCK_WORK_ZERO(i_work_) Memset(CODE_BLOCK_WORK_PTR(i_work_), 0, (size_t) len_work__ * sizeof(double))
#define CODE_BLOCK_ALL_WORK_ZERO() if (work__) Memset(CODE_BLOCK_WORK_PTR(0), 0, (size_t) (len_work__ * n_work__ * sizeof(double)))

// this avoids a potential second call to omp_get_thread_num() as if this is known, it can be passed as the 'thread_num_' argument
#define CODE_BLOCK_WORK_PTR_x(i_work_, thread_num_) (work__[thread_num_] + (size_t) (i_work_) * len_work__)
#define CODE_BLOCK_WORK_ZERO_x(i_work_, thread_num_) Memset(CODE_BLOCK_WORK_PTR_x(i_work_, thread_num_), 0, (size_t) len_work__ * sizeof(double))
#define CODE_BLOCK_ALL_WORK_ZERO_x(thread_num_) Memset(CODE_BLOCK_WORK_PTR_x(0, thread_num_), 0, (size_t) (len_work__ * n_work__ * sizeof(double)))


#define RUN_CODE_BLOCK(thread_max_, n_work_, len_work_)			\
	if (1) {							\
		int nt__ = ((GMRFLib_OPENMP_IN_PARALLEL_ONE_THREAD() || GMRFLib_OPENMP_IN_SERIAL()) ? \
			    IMAX(GMRFLib_openmp->max_threads_inner, GMRFLib_openmp->max_threads_outer) : GMRFLib_openmp->max_threads_inner); \
		int tmax__ = thread_max_;				\
		int len_work__ = IMAX(1, len_work_);			\
		int n_work__ = IMAX(1, n_work_);			\
		nt__ = IMAX(1, (tmax__ < 0 ? -tmax__ : IMAX(1, IMIN(nt__, tmax__)))); \
									\
		double ** work__ = Calloc(nt__, double *);		\
		for (int i_ = 0; i_ < nt__; i_++) {			\
			work__[i_] = Calloc(len_work__ * n_work__, double); \
			assert(work__[i_]);				\
		}							\
		assert(work__);						\
									\
		if (nt__ > 1) {						\
                        _Pragma("omp parallel for num_threads(nt__) schedule(guided)") \
				CODE_BLOCK;				\
		} else {						\
			CODE_BLOCK;					\
		}							\
		for (int i_ = 0; i_ < nt__; i_++) {			\
			Free(work__[i_]);				\
		}							\
		Free(work__);						\
        }

#define RUN_CODE_BLOCK_GUIDED(thread_max_, n_work_, len_work_)		\
	if (1) {							\
		int nt__ = ((GMRFLib_OPENMP_IN_PARALLEL_ONE_THREAD() || GMRFLib_OPENMP_IN_SERIAL()) ? \
			    IMAX(GMRFLib_openmp->max_threads_inner, GMRFLib_openmp->max_threads_outer) : GMRFLib_openmp->max_threads_inner); \
		int tmax__ = thread_max_;				\
		int len_work__ = IMAX(1, len_work_);			\
		int n_work__ = IMAX(1, n_work_);			\
		nt__ = IMAX(1, (tmax__ < 0 ? -tmax__ : IMAX(1, IMIN(nt__, tmax__)))); \
									\
		double ** work__ = Calloc(nt__, double *);		\
		for (int i_ = 0; i_ < nt__; i_++) {			\
			work__[i_] = Calloc(len_work__ * n_work__, double); \
			assert(work__[i_]);				\
		}							\
		assert(work__);						\
									\
		if (nt__ > 1) {						\
                        _Pragma("omp parallel for num_threads(nt__) schedule(guided)") \
				CODE_BLOCK;				\
		} else {						\
			CODE_BLOCK;					\
		}							\
		for (int i_ = 0; i_ < nt__; i_++) {			\
			Free(work__[i_]);				\
		}							\
		Free(work__);						\
        }

#define RUN_CODE_BLOCK_DYNAMIC(thread_max_, n_work_, len_work_)		\
	if (1) {							\
		int nt__ = ((GMRFLib_OPENMP_IN_PARALLEL_ONE_THREAD() || GMRFLib_OPENMP_IN_SERIAL()) ? \
			    IMAX(GMRFLib_openmp->max_threads_inner, GMRFLib_openmp->max_threads_outer) : GMRFLib_openmp->max_threads_inner); \
		int tmax__ = thread_max_;				\
		int len_work__ = IMAX(1, len_work_);			\
		int n_work__ = IMAX(1, n_work_);			\
		nt__ = IMAX(1, (tmax__ < 0 ? -tmax__ : IMAX(1, IMIN(nt__, tmax__)))); \
									\
		double ** work__ = Calloc(nt__, double *);		\
		for (int i_ = 0; i_ < nt__; i_++) { \
			work__[i_] = Calloc(len_work__ * n_work__, double); \
			assert(work__[i_]);				\
		}							\
		assert(work__);						\
									\
		if (nt__ > 1) {						\
			_Pragma("omp parallel for num_threads(nt__) schedule(dynamic)") \
				CODE_BLOCK;				\
		} else {						\
			CODE_BLOCK;					\
		}							\
		for (int i_ = 0; i_ < nt__; i_++) {			\
			Free(work__[i_]);				\
		}							\
		Free(work__);						\
        }

#define RUN_CODE_BLOCK_STATIC(thread_max_, n_work_, len_work_)		\
	if (1) {							\
		int nt__ = ((GMRFLib_OPENMP_IN_PARALLEL_ONE_THREAD() || GMRFLib_OPENMP_IN_SERIAL()) ? \
			    IMAX(GMRFLib_openmp->max_threads_inner, GMRFLib_openmp->max_threads_outer) : GMRFLib_openmp->max_threads_inner); \
		int tmax__ = thread_max_;				\
		int len_work__ = IMAX(1, len_work_);			\
		int n_work__ = IMAX(1, n_work_);			\
		nt__ = IMAX(1, (tmax__ < 0 ? -tmax__ : IMAX(1, IMIN(nt__, tmax__)))); \
									\
		double ** work__ = Calloc(nt__, double *);		\
		for (int i_ = 0; i_ < nt__; i_++) { \
			work__[i_] = Calloc(len_work__ * n_work__, double); \
			assert(work__[i_]);				\
		}							\
		assert(work__);						\
									\
		if (nt__ > 1) {						\
			_Pragma("omp parallel for num_threads(nt__) schedule(static)") \
				CODE_BLOCK;				\
		} else {						\
			CODE_BLOCK;					\
		}							\
		for (int i_ = 0; i_ < nt__; i_++) {			\
			Free(work__[i_]);				\
		}							\
		Free(work__);						\
        }

#define CODE_BLOCK_WORK_TP_PTR() work_t__[(nt__ == 1 ? 0 : omp_get_thread_num())]
// CODE_BLOCK_WORK_TP_FREE(ptr_) needs to be defined

#define RUN_CODE_BLOCK_X(thread_max_, n_work_, len_work_, work_tp_)	\
	if (1) {							\
		int nt__ = ((GMRFLib_OPENMP_IN_PARALLEL_ONE_THREAD() || GMRFLib_OPENMP_IN_SERIAL()) ? \
			    IMAX(GMRFLib_openmp->max_threads_inner, GMRFLib_openmp->max_threads_outer) : GMRFLib_openmp->max_threads_inner); \
		int tmax__ = thread_max_;				\
		int len_work__ = IMAX(1, len_work_);			\
		int n_work__ = IMAX(1, n_work_);			\
		nt__ = IMAX(1, (tmax__ < 0 ? -tmax__ : IMAX(1, IMIN(nt__, tmax__)))); \
									\
		work_tp_ ** work_t__ = Calloc(nt__, work_tp_ *);	\
		for (int i_ = 0; i_ < nt__; i_++) {			\
			work_t__[i_] = Calloc(1, work_tp_);		\
		}							\
									\
		double ** work__ = Calloc(nt__, double *);		\
		for (int i_ = 0; i_ < nt__; i_++) {			\
			work__[i_] = Calloc(len_work__ * n_work__, double); \
			assert(work__[i_]);				\
		}							\
		assert(work__);						\
									\
		if (nt__ > 1) {						\
			_Pragma("omp parallel for num_threads(nt__) schedule(static)") \
				CODE_BLOCK;				\
		} else {						\
			CODE_BLOCK;					\
		}							\
									\
		for (int i_ = 0; i_ < nt__; i_++) {			\
			CODE_BLOCK_WORK_TP_FREE(work_t__[i_]);		\
		}							\
		Free(work_t__);						\
		for (int i_ = 0; i_ < nt__; i_++) {			\
			Free(work__[i_]);				\
		}							\
		Free(work__);						\
        }

#define RUN_CODE_BLOCK_PLAIN(thread_max_, n_work_, len_work_)		\
	if (1) {							\
		int nt__ = ((GMRFLib_OPENMP_IN_PARALLEL_ONE_THREAD() || GMRFLib_OPENMP_IN_SERIAL()) ? \
			    IMAX(GMRFLib_openmp->max_threads_inner, GMRFLib_openmp->max_threads_outer) : GMRFLib_openmp->max_threads_inner); \
		int tmax__ = thread_max_;				\
		int len_work__ = IMAX(1, len_work_);			\
		int n_work__ = IMAX(1, n_work_);			\
		nt__ = IMAX(1, (tmax__ < 0 ? -tmax__ : IMAX(1, IMIN(nt__, tmax__)))); \
									\
		double ** work__ = Calloc(nt__, double *);		\
		for (int i_ = 0; i_ < nt__; i_++) { \
			work__[i_] = Calloc(len_work__ * n_work__, double); \
			assert(work__[i_]);				\
		}							\
		assert(work__);						\
									\
		CODE_BLOCK;						\
									\
		for (int i_ = 0; i_ < nt__; i_++) {			\
			Free(work__[i_]);				\
		}							\
		Free(work__);						\
        }


#define GMRFLib_INT_NUM_POINTS   (45)			       /* number of points for integration,... */
#define GMRFLib_INT_NUM_INTERPOL  (3)			       /* ...which are then interpolated: use 2 or 3 */
#define GMRFLib_INT_GHQ_POINTS   (15)			       /* MUST BE ODD!!!! for the quadrature */
#define GMRFLib_INT_GHQ_POINTS_PAD (1)			       /* So the _ALLOC_LEN is aligned well */
#define GMRFLib_INT_GHQ_ALLOC_LEN (GMRFLib_INT_GHQ_POINTS + GMRFLib_INT_GHQ_POINTS_PAD)

#define TIMER_INIT(use_, n_)				\
	static double timer_[2 + (n_)];			\
	static int timer_first_ = 1;			\
	const int timer_use_ = use_;			\
	int timer_n_ = 2 + n_;				\
	int timer_idx_ = -1;				\
	if (timer_use_ && timer_first_)	{		\
		timer_first_ = 0;			\
		GMRFLib_fill(timer_n_, 0.0, timer_);	\
	}

#define TIMER_CHECK							\
	if (timer_use_) {						\
		assert(timer_idx_ < timer_n_);				\
		double tim = GMRFLib_timer();				\
		if (timer_idx_ >= 0)					\
			timer_[timer_idx_] += tim;			\
		timer_[++timer_idx_] -= tim;				\
	}

#define TIMER_SUMMARY							\
	TIMER_CHECK;							\
	if (timer_use_ && timer_idx_ > 0) {				\
		double sum = GMRFLib_dsum(timer_idx_, timer_);		\
		printf("\n%s:%d: (%s) relative ", __FILE__, __LINE__, __GMRFLib_FuncName); \
		for (int i_ = 0; i_ < timer_idx_; i_++) {		\
			printf(" [%1d] %.4f", i_, timer_[i_] / sum);	\
		}							\
		printf("\n\n");						\
	}								\


/* from /usr/include/assert.h. use __GMRFLib_FuncName to define name of current function.

   Version 2.4 and later of GCC define a magical variable `__PRETTY_FUNCTION__' which contains the
   name of the function currently being defined.  This is broken in G++ before version 2.6.  C9x has
   a similar variable called __func__, but prefer the GCC one since it demangles C++ function names.
*/
#ifndef __GNUC_PRERQ
#if defined __GNUC__ && defined __GNUC_MINOR__
#define __GNUC_PREREQ(maj, min) ((__GNUC__ << 16) + __GNUC_MINOR__ >= ((maj) << 16) + (min))
#else
#define __GNUC_PREREQ(maj, min) 0
#endif
#endif
#if defined __GNUC__
#if defined __cplusplus ? __GNUC_PREREQ (2, 6) : __GNUC_PREREQ (2, 4)
#define  __GMRFLib_FuncName   __PRETTY_FUNCTION__
#else
#if defined __STDC_VERSION__ && __STDC_VERSION__ >= 199901L
#define __GMRFLib_FuncName  __func__
#else
#define __GMRFLib_FuncName  ((const char *) "(function-name unavailable)")
#endif
#endif
#else
#if defined(__sun)					       /* it works here */
#define __GMRFLib_FuncName __func__
#else
#define __GMRFLib_FuncName ((const char *) "(function-name unavailable)")
#endif
#endif

// from https://en.wikipedia.org/wiki/Inline_function
#ifdef _MSC_VER
#define forceinline __forceinline
#elif defined(__GNUC__)
#define forceinline inline __attribute__((__always_inline__))
#elif defined(__CLANG__)
#if __has_attribute(__always_inline__)
#define forceinline inline __attribute__((__always_inline__))
#else
#define forceinline inline
#endif
#else
#define forceinline inline
#endif

/* 
   parts taken from /usr/include/tcl.h
 */
#ifdef __STRING
#define __GMRFLib_STRINGIFY(x) __STRING(x)
#else
#ifdef _MSC_VER
#define __GMRFLib_STRINGIFY(x) #x
#else
#ifdef RESOURCE_INCLUDED
#define __GMRFLib_STRINGIFY(x) #x
#else
#ifdef __STDC__
#define __GMRFLib_STRINGIFY(x) #x
#else
#define __GMRFLib_STRINGIFY(x) "x"
#endif
#endif
#endif
#endif


__END_DECLS
#endif
