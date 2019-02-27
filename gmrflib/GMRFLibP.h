
/* GMRFLibP.h
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
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
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
typedef short int GMRFLib_short_int;
typedef long unsigned int GMRFLib_sizeof_tp;


/* 
   here are the wrappers for calling functions which return the error-code if it fails
*/
#define GMRFLib_EWRAP__intern(func_call, leave)		\
	if (1) {					\
		int rretval;				\
		rretval = func_call;			\
		if (rretval != GMRFLib_SUCCESS){	\
			if (leave) GMRFLib_LEAVE_ROUTINE;	\
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
			if (leave) GMRFLib_LEAVE_ROUTINE;		\
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
			if (leave) GMRFLib_LEAVE_ROUTINE;		\
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
			if (leave) GMRFLib_LEAVE_ROUTINE;		\
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
		_tref = GMRFLib_cpu();			\
		_ntimes++;				\
		expression;						\
		_tacc += GMRFLib_cpu() - _tref;				\
		printf("%s:%s:%d: cpu accumulative [%s] %.6f mean %.8f n %d\n", \
		       __FILE__, __GMRFLib_FuncName, __LINE__, msg, _tacc, _tacc/_ntimes, _ntimes); \
	}
#define GMRFLib_MEASURE_CPU_BEGIN			\
	if (1) {					\
		static double _tacc = 0.0;		\
		static int _ntimes = 0;			\
		double _tref;				\
		_tref = GMRFLib_cpu();			\
		_ntimes++;
#define GMRFLib_MEASURE_CPU_END(msg)					\
	_tacc += GMRFLib_cpu() - _tref;					\
	printf("%s:%s:%d: cpu accumulative [%s] %.6f mean %.8f n %d\n",	\
	       __FILE__, __GMRFLib_FuncName, __LINE__, msg, _tacc, _tacc/_ntimes, _ntimes); \
	}

/* 
   some useful macros
*/
#if 1
//#define GMRFLib_TRACE_MEMORY    1000000   // trace memory larger than this ammount. undefine it to disable this feature.
#define Calloc(n, type)         (type *)GMRFLib_calloc((size_t)(n),sizeof(type), __FILE__, __GMRFLib_FuncName, __LINE__, RCSId)
#define Malloc(n, type)         (type *)GMRFLib_malloc((size_t)(n)*sizeof(type), __FILE__, __GMRFLib_FuncName, __LINE__, RCSId)
#define Realloc(ptr, n, type)   (type *)GMRFLib_realloc((void *)ptr, (size_t)(n)*sizeof(type), __FILE__, __GMRFLib_FuncName, __LINE__, RCSId)
#define Free(ptr)               {GMRFLib_free((void *)(ptr), __FILE__, __GMRFLib_FuncName, __LINE__, RCSId); ptr=NULL;}
#else
#undef  GMRFLib_TRACE_MEMORY
#define Calloc(n, type)         (type *)calloc((size_t)(n),sizeof(type))
#define Malloc(n, type)         (type *)malloc((size_t)(n)*sizeof(type))
#define Realloc(ptr, n, type)   (type *)realloc((void *)ptr, (size_t)(n)*sizeof(type))
#define Free(ptr)               {free((void *)(ptr)); ptr=NULL;}
#endif

/* 
   ABS is for double, IABS is for int, and so on.
*/
#define ABS(x)   fabs(x)
#define DMAX(a,b) GSL_MAX_DBL(a, b)
#define DMIN(a,b) GSL_MIN_DBL(a, b)
#define TRUNCATE(x, low, high)  DMIN( DMAX(x, low), high)      /* ensure that x is in the inteval [low,high] */
#define SQR(x) gsl_pow_2(x)
#define IABS(x)   abs(x)
#define IMAX(a,b) GSL_MAX_INT(a, b)
#define IMIN(a,b) GSL_MIN_INT(a, b)
#define ITRUNCATE(x, low, high) IMIN(IMAX(x, low), high)
#define ISQR(x) ((x)*(x))
#define MOD(i,n)  (((i)+(n))%(n))
#define FIXME( msg) if (1) { printf("\n[%1d]:%s:%1d:%s: FIXME [%s]\n",  omp_get_thread_num(), __FILE__, __LINE__, __GMRFLib_FuncName,(msg?msg:""));	}
#define FIXME1(msg) if (1) { static int first=1; if (first) { first=0; FIXME(msg); }}
#define FIXMEstderr( msg) if (1) { fprintf(stderr, "\n[%1d]:%s:%1d:%s: FIXME [%s]\n",  omp_get_thread_num(), __FILE__, __LINE__, __GMRFLib_FuncName,(msg?msg:""));	}
#define FIXME1stderr(msg) if (1) { static int first=1; if (first) { first=0; FIXMEstderr(msg); }}
#define P(x)        if (1) { printf("line[%1d] " #x " = [ %.12f ]\n",__LINE__,(double)(x)); }
#define Pstderr(x)  if (1) { fprintf(stderr, "line[%1d] " #x " = [ %.12f ]\n",__LINE__,(double)(x)); }
#define P1(x)       if (1) { static int first=1;  if (first) { printf("line[%1d] " #x " = [ %.12f ]\n", __LINE__, (double)(x)); first=0; }}
#define P1stderr(x) if (1) { static int first=1;  if (first) { fprintf(stderr, "line[%1d] " #x " = [ %.12f ]\n", __LINE__, (double)(x)); first=0; }}
#define PP(msg,pt) if (1) { fprintf(stdout, "%d: %s ptr " #pt " = 0x%p\n", __LINE__, msg, pt); }
#define PPstderr(msg,pt)  if (1) { fprintf(stderr, "%d: %s ptr " #pt " = 0x%p\n", __LINE__, msg, pt); }
#define PPg(msg,pt) if (1) { fprintf(stdout, "%d: %s value " #pt " = %g\n", __LINE__, msg, pt); }
#define PPstderrg(msg,pt) if (1) { fprintf(stderr, "%d: %s value " #pt " = %g\n", __LINE__, msg, pt); }
#define ISINF(x) gsl_isinf(x)
#define ISNAN(x) gsl_isnan(x)
#define ISZERO(x) (gsl_fcmp(x, 0.0, DBL_EPSILON) == 0)
#define ISEQUAL(x, y) (gsl_fcmp(x, y, DBL_EPSILON) == 0)
#define LEGAL(i, n) ((i) >= 0 && (i) < (n))

#define GMRFLib_GLOBAL_NODE(n, gptr) ((int) IMIN((n-1)*(gptr ? (gptr)->factor :  GMRFLib_global_node.factor), \
						 (gptr ? (gptr)->degree : GMRFLib_global_node.degree)))

#define GMRFLib_STOP_IF_NAN_OR_INF(value, idx, jdx)			\
	if (ISNAN(value) || ISINF(value)) {				\
		if (!nan_error)						\
			fprintf(stdout,					\
				"\n\t%s\n\tFunction: %s(), Line: %1d, Thread: %1d\n\tVariable evaluates to NAN or INF. idx=(%1d,%1d). I will try to fix it...", \
				RCSId, __GMRFLib_FuncName, __LINE__, omp_get_thread_num(), idx, jdx); \
		if (GMRFLib_catch_error_for_inla) {			\
			nan_error = 1;					\
		} else {						\
			if (1)abort();					\
			if (1)exit(1);					\
		}							\
	}

#define GMRFLib_SET_PREC(arg_)						\
	(arg_->prec ? *(arg_->prec)					\
	 : (arg_->log_prec ? exp(*(arg_->log_prec))			\
	    : (arg_->log_prec_omp ? exp(*(arg_->log_prec_omp[GMRFLib_thread_id])) : 1.0)))

#define GMRFLib_SET_RANGE(arg_)						\
	(arg_->range ? *(arg_->range)					\
	 : (arg_->log_range ? exp(*(arg_->log_range))			\
	    : (arg_->log_range_omp ? exp(*(arg_->log_range_omp[GMRFLib_thread_id])) : 1.0)))

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
