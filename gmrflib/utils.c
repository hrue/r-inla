
/* utils.c
 * 
 * Copyright (C) 2006-2010 Havard Rue
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
  \file utils.c
  \brief Various simple utilities.
*/

#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

/* Pre-hg-Id: $Id: utils.c,v 1.97 2010/02/15 08:26:40 hrue Exp $ */

#if !defined(__FreeBSD__)
#include <stddef.h>
#include <string.h>
#include <malloc.h>
#endif
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"
#include "GMRFLib/hashP.h"

/* 
 *  compile with -DGMRFLib_MEMCHECK to enable the internal memcheck utility. Do not work with OPENMP. compile with -DGMRFLib_MEMINFO to enable the meminfo utility.
 *  compute with -DGMRFLib_TRACE_MEMORY to view memory allocation
 */

static map_vpvp memcheck_hash_table;
static int memcheck_verbose = 0;
static int memcheck_first = 1;

#if defined(GMRFLib_MEMINFO)
static GMRFLib_meminfo_tp *MemInfo = NULL;
#define MEMINFO(_size) \
	if (GMRFLib_meminfo_thread_id != 0) {				\
		assert(IABS(GMRFLib_meminfo_thread_id) < 1+omp_get_max_threads()); \
		if (!MemInfo){						\
			_Pragma("omp critial")				\
			{						\
				if (!MemInfo)				\
					MemInfo = (GMRFLib_meminfo_tp *) calloc(1+omp_get_max_threads(), sizeof(GMRFLib_meminfo_tp)); \
			}						\
		}							\
		if (GMRFLib_meminfo_thread_id > 0) {			\
			MemInfo[GMRFLib_meminfo_thread_id].n++;	\
			MemInfo[GMRFLib_meminfo_thread_id].bytes += _size; \
		} else if (GMRFLib_meminfo_thread_id < 0 && MemInfo[-GMRFLib_meminfo_thread_id].n > 0) { \
			_Pragma("omp critial")				\
			{						\
				if (GMRFLib_meminfo_thread_id < 0 && MemInfo[-GMRFLib_meminfo_thread_id].n > 0) { \
					printf("\nMemInfo thread_id %d\n", -GMRFLib_meminfo_thread_id); \
					printf("\tn %g\n", (double) MemInfo[-GMRFLib_meminfo_thread_id].n); \
					printf("\tb %g Mb\n\n", ((double) MemInfo[-GMRFLib_meminfo_thread_id].bytes)/1024.0/1024.0); \
					MemInfo[-GMRFLib_meminfo_thread_id].n = MemInfo[-GMRFLib_meminfo_thread_id].bytes = 0; \
				}					\
			}						\
		}							\
	}
#else
#define MEMINFO(_size) if (0) {}
#endif

char *GMRFLib_rindex(const char *p, int ch)
{
	/*
	 * as Windows does not have it... 
	 */
	char *save, *pp = (char *) p;
	for (save = NULL;; ++pp) {
		if (*pp == ch) {
			save = pp;
		}
		if (!*pp) {
			return (save);
		}
	}
	abort();
	return NULL;
}

int GMRFLib_which(double val, double *array, int len)
{
	/*
	 * return the first index in array such that array[idx] == val, and -1 if not there 
	 */
	int i;

	for (i = 0; i < len; i++) {
		if (ISEQUAL(val, array[i])) {
			return i;
		}
	}
	return -1;
}
int GMRFLib_find_nonzero(double *array, int len, int direction)
{
	/*
	 * return the first/last index in array such that array[idx] != 0, and -1 if not there. direction > 0 : look for first. direction < 0 : look for last
	 */
	int i;

	if (direction >= 0) {
		for (i = 0; i < len; i++) {
			if (array[i] != 0.0)
				return i;
		}
		return -1;
	} else {
		for (i = len - 1; i >= 0; i--) {
			if (array[i] != 0.0)
				return i;
		}
		return -1;
	}

	return -1;
}
double GMRFLib_eps(double power)
{
	/*
	 * Return eps^power, where eps is the smalles number such that 1+eps != eps.  
	 
	 * i had to rewrite this due to compiler behaviour using `-ffast-math' and gcc version 4.1. it seems ok now. 
	 */

	static double eps = -1.0;

	if (eps < 0.0) {
		double val;

		eps = 1.0;
		val = 2.0;
		while (val > 1.0) {
			eps /= 2.0;
			val = 1.0 + eps;
		}
		eps *= 2.0;
	}
	return pow(eps, power);
}
int GMRFLib_print_darray(FILE * fp, double *x, int n, const char *desc)
{
	int i;

	fp = (fp ? fp : stdout);
	fprintf(fp, "Double array with length %1d [%1s]\n", n, (desc ? desc : "no description given"));
	for (i = 0; i < n; i++) {
		fprintf(fp, "\telement %1d = %.10f\n", i, x[i]);
	}
	return GMRFLib_SUCCESS;
}
int GMRFLib_print_iarray(FILE * fp, int *x, int n, const char *desc)
{
	int i;

	fp = (fp ? fp : stdout);
	fprintf(fp, "Integer array with length %1d [%1s]\n", n, (desc ? desc : "no description given"));
	for (i = 0; i < n; i++) {
		fprintf(fp, "\telement %1d = %1d\n", i, x[i]);
	}
	return GMRFLib_SUCCESS;
}
int GMRFLib_adjust_vector(double *x, int n)
{
	/*
	 * x := x - max(x[]) 
	 */
	int i;
	double max_value;

	if (n <= 0 || !x) {
		return GMRFLib_SUCCESS;
	}

	max_value = GMRFLib_max_value(x, n, NULL);
	for (i = 0; i < n; i++) {
		x[i] -= max_value;
	}

	return GMRFLib_SUCCESS;
}
int GMRFLib_scale_vector(double *x, int n)
{
	/*
	 * x := x/max(x)
	 */
	int i;
	double scale;

	if (n <= 0 || !x) {
		return GMRFLib_SUCCESS;
	}

	scale = GMRFLib_max_value(x, n, NULL);
	if (!ISZERO(scale)) {
		scale = 1.0 / scale;
		for (i = 0; i < n; i++) {
			x[i] *= scale;
		}
	}

	return GMRFLib_SUCCESS;
}
double GMRFLib_min_value(double *x, int n, int *idx)
{
	/*
	 * return the MIN(x[]) 
	 */
	int i;
	double min_val;

	min_val = x[0];
	if (idx) {
		*idx = 0;
	}
	for (i = 1; i < n; i++) {
		if (x[i] < min_val && !ISNAN(x[i])) {
			min_val = x[i];
			if (idx) {
				*idx = i;
			}
		}
	}

	return min_val;
}
int GMRFLib_imin_value(int *x, int n, int *idx)
{
	/*
	 * return the IMIN(x[]) 
	 */
	int i;
	int min_val;

	min_val = x[0];
	if (idx) {
		*idx = 0;
	}
	for (i = 1; i < n; i++) {
		if (x[i] < min_val && !ISNAN(x[i])) {
			min_val = x[i];
			if (idx) {
				*idx = i;
			}
		}
	}

	return min_val;
}
double GMRFLib_max_value(double *x, int n, int *idx)
{
	/*
	 * return the MAX(x[]), optional idx
	 */
	int i;
	double max_val;

	max_val = x[0];
	if (idx) {
		*idx = 0;
	}
	for (i = 1; i < n; i++) {
		if (x[i] > max_val && !ISNAN(x[i])) {
			max_val = x[i];
			if (idx) {
				*idx = i;
			}
		}
	}

	return max_val;
}
int GMRFLib_imax_value(int *x, int n, int *idx)
{
	/*
	 * return IMAX(x[]) 
	 */

	int i;
	int max_val;

	max_val = x[0];
	if (idx) {
		*idx = 0;
	}
	for (i = 1; i < n; i++) {
		if (x[i] > max_val && !ISNAN(x[i])) {
			max_val = x[i];
			if (idx) {
				*idx = i;
			}
		}
	}

	return max_val;
}
int GMRFLib_icmp(const void *a, const void *b)
{
	const int *ia = NULL, *ib = NULL;

	ia = (const int *) a;
	ib = (const int *) b;

	if (*ia > *ib) {
		return 1;
	}
	if (*ia < *ib) {
		return -1;
	}

	return 0;
}
int GMRFLib_dcmp(const void *a, const void *b)
{
	const double *da = NULL, *db = NULL;

	da = (const double *) a;
	db = (const double *) b;

	if (*da > *db) {
		return 1;
	}
	if (*da < *db) {
		return -1;
	}

	return 0;
}
int GMRFLib_dcmp_r(const void *a, const void *b)
{
	const double *da = NULL, *db = NULL;

	da = (const double *) a;
	db = (const double *) b;

	if (*da > *db) {
		return -1;
	}
	if (*da < *db) {
		return 1;
	}

	return 0;
}
int GMRFLib_dcmp_abs(const void *a, const void *b)
{
	const double *da = NULL, *db = NULL;

	da = (const double *) a;
	db = (const double *) b;

	/*
	 * sort on ABS() 
	 */
	if (ABS(*da) > ABS(*db)) {
		return 1;
	}
	if (ABS(*da) < ABS(*db)) {
		return -1;
	}

	/*
	 * if they're equal, sort on sign 
	 */
	if ((*da) > (*db)) {
		return 1;
	}
	if ((*da) < (*db)) {
		return -1;
	}

	/*
	 * identical 
	 */
	return 0;
}
int GMRFLib_qsorts(void *x, size_t nmemb, size_t size_x, void *y, size_t size_y, void *z, size_t size_z, int (*compar) (const void *, const void *))
{
	/*
	 * sort x and optionally sort y and z along
	 * 
	 * if z is non-NULL, then y must be so as well. 
	 */

	char *xyz = NULL, *xx = NULL, *yy = NULL, *zz = NULL;
	size_t siz, i, offset;

	if (nmemb == 0) {
		return GMRFLib_SUCCESS;
	}
	if (z) {
		GMRFLib_ASSERT(y, GMRFLib_EINVARG);
	}
	xx = (char *) x;
	yy = (char *) y;
	zz = (char *) z;

	siz = size_x;
	if (y) {
		siz += size_y;
		if (z) {
			siz += size_z;
		}
	}
	xyz = Calloc(nmemb * siz, char);

	for (i = 0; i < nmemb; i++) {
		memcpy((void *) &xyz[i * siz], (void *) &xx[i * size_x], size_x);
	}
	if (y) {
		offset = size_x;
		for (i = 0; i < nmemb; i++) {
			memcpy((void *) &xyz[i * siz + offset], (void *) &yy[i * size_y], size_y);
		}
		if (z) {
			offset = size_x + size_y;
			for (i = 0; i < nmemb; i++) {
				memcpy((void *) &xyz[i * siz + offset], (void *) &zz[i * size_z], size_z);
			}
		}
	}

	qsort((void *) xyz, nmemb, siz, compar);

	for (i = 0; i < nmemb; i++) {
		memcpy((void *) &xx[i * size_x], (void *) &xyz[i * siz], size_x);
	}
	if (y) {
		offset = size_x;
		for (i = 0; i < nmemb; i++) {
			memcpy((void *) &yy[i * size_y], (void *) &xyz[i * siz + offset], size_y);
		}
		if (z) {
			offset = size_x + size_y;
			for (i = 0; i < nmemb; i++) {
				memcpy((void *) &zz[i * size_z], (void *) &xyz[i * siz + offset], size_z);
			}
		}
	}

	Free(xyz);

	return GMRFLib_SUCCESS;
}
double GMRFLib_log_apbex(double a, double b)
{
	/*
	 * try to evaluate log(a + exp(b)) safely 
	 */

	if (a == 0.0)
		return b;

	double B = exp(b);

	if (B > a) {
		return b + log(1.0 + a / B);
	} else {
		return log(a) + log(1.0 + B / a);
	}
}

void *GMRFLib_calloc(size_t nmemb, size_t size, const char *file, const char *funcname, int lineno, const char *id)
{
	void *ptr = NULL;
	char *msg = NULL;

	MEMINFO(nmemb * size);

#if defined(GMRFLib_MEMCHECK) && !defined(_OPENMP)
	ptr = GMRFLib_calloc__(nmemb, size, file, funcname, lineno, id);
#else
	ptr = calloc(nmemb, size);
#endif
#if defined(GMRFLib_TRACE_MEMORY)
	if (nmemb * size > GMRFLib_TRACE_MEMORY)
		printf("%s:%s:%u: calloc %zu x %zu bytes, total %zu\n", file, funcname, lineno, nmemb, size, nmemb * size);
#endif
	if (ptr) {
		return ptr;
	}

	/*
	 * alloc failed 
	 */
	GMRFLib_sprintf(&msg, "Fail to calloc nmemb=%1zu elements of size=%1zu bytes", nmemb, size);
	GMRFLib_handle_error(file, funcname, lineno, id, GMRFLib_EMEMORY, msg);
	abort();

	return NULL;
}
void *GMRFLib_malloc(size_t size, const char *file, const char *funcname, int lineno, const char *id)
{
	void *ptr = NULL;
	char *msg = NULL;

	MEMINFO(size);

#if defined(GMRFLib_MEMCHECK) && !defined(_OPENMP)
	ptr = GMRFLib_malloc__(size, file, funcname, lineno, id);
#else
	ptr = malloc(size);
#endif
#if defined(GMRFLib_TRACE_MEMORY)
	if (size > GMRFLib_TRACE_MEMORY)
		printf("%s:%s:%u: malloc %zu bytes, total %zu\n", file, funcname, lineno, size, size);
#endif
	if (ptr) {
		return ptr;
	}

	/*
	 * alloc failed 
	 */
	GMRFLib_sprintf(&msg, "Fail to malloc size=%1zu bytes", size);
	GMRFLib_handle_error(file, funcname, lineno, id, GMRFLib_EMEMORY, msg);
	abort();

	return NULL;
}
void *GMRFLib_realloc(void *old_ptr, size_t size, const char *file, const char *funcname, int lineno, const char *id)
{
	void *ptr = NULL;
	char *msg = NULL;

	MEMINFO(size);

#if defined(GMRFLib_MEMCHECK) && !defined(_OPENMP)
	ptr = GMRFLib_realloc__(old_ptr, size, file, funcname, lineno, id);
#else
	ptr = realloc(old_ptr, size);
#endif
#if defined(GMRFLib_TRACE_MEMORY)
	if (size > GMRFLib_TRACE_MEMORY)
		printf("%s:%s:%u: realloc %zu bytes, total %zu\n", file, funcname, lineno, size, size);
#endif
	if (ptr) {
		return ptr;
	}

	/*
	 * realloc failed 
	 */
	GMRFLib_sprintf(&msg, "Fail to realloc size=%1zu bytes", size);
	GMRFLib_handle_error(file, funcname, lineno, id, GMRFLib_EMEMORY, msg);
	abort();

	return NULL;
}
void GMRFLib_free(void *ptr, const char *file, const char *funcname, int lineno, const char *id)
{
	if (ptr) {
#if defined(GMRFLib_MEMCHECK) && !defined(_OPENMP)
		GMRFLib_free__(ptr, file, funcname, lineno, id);
#else
		free(ptr);
#endif
	}
}
char *GMRFLib_strdup(const char *ptr)
{
	/*
	 * strdup is not part of ANSI-C 
	 */
#if defined(__STRICT_ANSI__) || defined(GMRFLib_MEMCHECK)
	if (ptr) {
		size_t len = strlen(ptr);
		char *str = Malloc(len + 1, char);

		return strcpy(str, ptr);
	} else {
		return (char *) NULL;
	}
#else
	if (ptr) {
		char *p = strdup(ptr);
		GMRFLib_ASSERT_RETVAL(p, GMRFLib_EMEMORY, (char *) NULL);

		return p;
	} else {
		return NULL;
	}
#endif
}
void *GMRFLib_malloc__(size_t size, const char *file, const char *funcname, int lineno, const char *id)
{
	void *p = NULL;

	p = malloc(size);
	GMRFLib_memcheck_register(p, size, file, funcname, lineno, id);
	return p;
}
void *GMRFLib_calloc__(size_t nmemb, size_t size, const char *file, const char *funcname, int lineno, const char *id)
{
	void *p = NULL;

	p = calloc(nmemb, size);
	GMRFLib_memcheck_register(p, size, file, funcname, lineno, id);
	return p;
}
void *GMRFLib_realloc__(void *old_ptr, size_t size, const char *file, const char *funcname, int lineno, const char *id)
{
	void *p = NULL;

	p = realloc(old_ptr, size);
	GMRFLib_memcheck_remove(old_ptr, file, funcname, lineno, id);
	GMRFLib_memcheck_register(p, size, file, funcname, lineno, id);
	return p;
}
void GMRFLib_free__(void *ptr, const char *file, const char *funcname, int lineno, const char *id)
{
	if (!ptr) {
		return;
	}
	if (memcheck_first) {
		map_vpvp_init_hint(&memcheck_hash_table, 10000);
		memcheck_first = 0;
	}
	if (ptr && !map_vpvp_ptr(&memcheck_hash_table, ptr)) {
		GMRFLib_memcheck_error("try to free a ptr not alloced", ptr, file, funcname, lineno, id);
	}
	GMRFLib_memcheck_remove(ptr, file, funcname, lineno, id);
	free(ptr);

	return;
}
char *GMRFLib_memcheck_make_tag(size_t size, const char *file, const char *funcname, int lineno, const char *id)
{
	char *tag = NULL;

	tag = (char *) calloc(1024 + 1, sizeof(char));
	snprintf(tag, 1024, "size: %zu bytes \tfile: %s \tfuncname: %s \tlineno: %d \tid: %s", size, file, funcname, lineno, id);

	return tag;
}
int GMRFLib_memcheck_error(const char *msg, void *p, const char *file, const char *funcname, int lineno, const char *id)
{
	printf("GMRFLib_memcheck: %s [0x%" PRIxPTR "]\n", msg, (uintptr_t) p);
	printf("called from file %s funcname %s lineno %d id %s\n", file, funcname, lineno, id);
	abort();
	return GMRFLib_SUCCESS;
}
int GMRFLib_memcheck_register(void *p, size_t size, const char *file, const char *funcname, int lineno, const char *id)
{
	void *tag = NULL;

	if (!p) {
		return GMRFLib_SUCCESS;
	}
	if (memcheck_first) {
		map_vpvp_init_hint(&memcheck_hash_table, 10000);
		memcheck_first = 0;
	}
	if (map_vpvp_ptr(&memcheck_hash_table, p)) {
		GMRFLib_memcheck_error("register a ptr not free'd", p, file, funcname, lineno, id);
	} else {
		tag = (void *) GMRFLib_memcheck_make_tag(size, file, funcname, lineno, id);
		map_vpvp_set(&memcheck_hash_table, p, tag);
	}

	if (memcheck_verbose) {
		printf("%s: 0x%" PRIxPTR " %sn\n", __GMRFLib_FuncName, (uintptr_t) p, (char *) tag);
	}

	return GMRFLib_SUCCESS;
}
int GMRFLib_memcheck_remove(void *p, const char *file, const char *funcname, int lineno, const char *id)
{
	void *ptr = NULL;

	if (!p) {
		return GMRFLib_SUCCESS;
	}
	if (memcheck_first) {
		map_vpvp_init_hint(&memcheck_hash_table, 10000);
		memcheck_first = 0;
	}
	ptr = map_vpvp_ptr(&memcheck_hash_table, p);
	if (!ptr) {
		GMRFLib_memcheck_error("remove ptr not registered", p, file, funcname, lineno, id);
	}
	if (memcheck_verbose)
		printf("%s: 0x%" PRIxPTR " %sn\n", __GMRFLib_FuncName, (uintptr_t) p, *((char **) ptr));

	free(*((void **) ptr));
	map_vpvp_remove(&memcheck_hash_table, p);

	return GMRFLib_SUCCESS;
}

int GMRFLib_memcheck_printf(FILE * fp)
{
#if defined(GMRFLib_MEMCHECK) && !defined(_OPENMP)
	int i;

	fp = (fp ? fp : stdout);
	if (!memcheck_first) {
		for (i = -1; (i = map_vpvp_next(&memcheck_hash_table, i)) != -1;) {
			fprintf(fp, "0x%x %s\n", (size_t) memcheck_hash_table.contents[i].key, (char *) memcheck_hash_table.contents[i].value);
		}
	}
#endif
	return GMRFLib_SUCCESS;
}
int GMRFLib_unique_relative(int *n, double *x, double eps)
{
	/*
	 * assume x is sorted, remove ties and change *n accordingly.
	 * 
	 * ties are defined if relative error between x_i and x_j <= eps, roughly, by using the routine gsl_fcmp()
	 * 
	 */
	int i = 0, j = 0, jj = 0;

	while (jj < *n) {
		while (jj + 1 < *n && gsl_fcmp(x[jj + 1], x[j], eps) == 0) {
			jj++;
		}
		x[i] = x[jj];
		i++;
		jj++;
		j = jj;
	}
	*n = i;

	return GMRFLib_SUCCESS;
}
int GMRFLib_unique_relative2(int *n, double *x, double *y, double eps)
{
	/*
	 * assume x is sorted, remove ties and change *n accordingly. make the same changes to y.
	 * 
	 * ties are defined if relative error between x_i and x_j <= eps, roughly, by using the routine gsl_fcmp()
	 * 
	 */

	if (!y) {
		return GMRFLib_unique_relative(n, x, eps);
	}

	int i = 0, j = 0, jj = 0;

	while (jj < *n) {
		while (jj + 1 < *n && gsl_fcmp(x[jj + 1], x[j], eps) == 0) {
			jj++;
		}
		x[i] = x[jj];
		y[i] = y[jj];
		i++;
		jj++;
		j = jj;
	}
	*n = i;

	return GMRFLib_SUCCESS;
}
int GMRFLib_unique_additive(int *n, double *x, double eps)
{
	/*
	 * assume x is sorted, remove ties and change *n accordingly. use the median in each bin
	 * 
	 * ties are defined if |x_i and x_j|  <= eps
	 * 
	 */
	int i = 0, j = 0, jj = 0;

	while (jj < *n) {
		while (jj + 1 < *n && fabs(x[jj + 1] - x[j]) <= eps) {
			jj++;
		}
		x[i] = x[jj];
		i++;
		jj++;
		j = jj;
	}
	*n = i;

	return GMRFLib_SUCCESS;
}
int GMRFLib_unique_additive2(int *n, double *x, double *y, double eps)
{
	/*
	 * assume x is sorted, remove ties and change *n accordingly. use the median in each bin. make the same changes to y.
	 * 
	 * ties are defined if |x_i and x_j|  <= eps
	 * 
	 */

	if (!y) {
		return GMRFLib_unique_additive(n, x, eps);
	}

	int i = 0, j = 0, jj = 0;

	while (jj < *n) {
		while (jj + 1 < *n && fabs(x[jj + 1] - x[j]) <= eps) {
			jj++;
		}
		x[i] = x[jj];
		y[i] = y[jj];
		i++;
		jj++;
		j = jj;
	}
	*n = i;

	return GMRFLib_SUCCESS;
}

int GMRFLib_matrix_fprintf(FILE * fp, double *A, int m, int n)
{
	// A is m x n matrix
#pragma omp critial
	{
		fprintf(fp, "\n\n");
		for (int i = 0; i < m; i++) {
			fprintf(fp, "\t");
			for (int j = 0; j < n; j++)
				fprintf(fp, " %10.6f", A[i + j * m]);
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
		fflush(fp);
	}
	return 0;
}


int GMRFLib_gsl_matrix_fprintf(FILE * fp, gsl_matrix * matrix, const char *format)
{
	size_t i, j;

	for (i = 0; i < matrix->size1; i++) {
		for (j = 0; j < matrix->size2; j++) {
			fprintf(fp, (format ? format : " %g"), gsl_matrix_get(matrix, i, j));
		}
		fprintf(fp, "\n");
	}
	return GMRFLib_SUCCESS;
}
double GMRFLib_signed_pow(double x, double power)
{
	if (ISZERO(x)) {
		return 0.0;
	} else {
		return (x > 0.0 ? 1.0 : -1.0) * pow(fabs(x), power);
	}
}
int GMRFLib_2order_poleq(double *sol1, double *sol2, double a, double b, double c)
{
	/*
	 * solve the equation a*x^2 + b*x + c,
	 * 
	 * return the two solutions in sol1 and sol2. return GMRFLib_EMISC if there is complex solutions. 
	 */

	double tmp = SQR(b) - 4.0 * a * c;

	if (tmp < 0.0) {
		return GMRFLib_EMISC;
	}

	if (ISZERO(a)) {
		if (sol1) {
			*sol1 = -c / b;
		}
		if (sol2) {
			*sol2 = -c / b;
		}
	} else {
		if (sol1) {
			*sol1 = (-b + sqrt(tmp)) / (2.0 * a);
		}
		if (sol2) {
			*sol2 = (-b - sqrt(tmp)) / (2.0 * a);
		}
	}
	return GMRFLib_SUCCESS;
}

mapkit_size_t GMRFLib_nelm_map_ii(map_ii * hash)
{
	/*
	 * return the number of elements in HASH 
	 */

	mapkit_size_t i, nelm = 0;

	for (i = -1; (i = map_ii_next(hash, i)) != -1;) {
		nelm++;
	}
	return nelm;
}

mapkit_size_t GMRFLib_nelm_map_id(map_id * hash)
{
	/*
	 * return the number of elements in HASH 
	 */
	mapkit_size_t i, nelm = 0;

	for (i = -1; (i = map_id_next(hash, i)) != -1;) {
		nelm++;
	}
	return nelm;
}

map_ii *GMRFLib_duplicate_map_ii(map_ii * hash)
{
	/*
	 * return a copy of HASH 
	 */
	if (!hash) {
		return NULL;
	}

	mapkit_size_t i, nelm;
	map_ii *newhash = NULL;

	nelm = GMRFLib_nelm_map_ii(hash);
	newhash = Calloc(1, map_ii);
	map_ii_init_hint(newhash, nelm);

	for (i = -1; (i = map_ii_next(hash, i)) != -1;) {
		map_ii_set(newhash, hash->contents[i].key, hash->contents[i].value);
	}
	newhash->alwaysdefault = hash->alwaysdefault;
	newhash->defaultvalue = hash->defaultvalue;

	return newhash;
}
GMRFLib_sizeof_tp GMRFLib_sizeof_map_ii(map_ii * hash)
{
	if (!hash) {
		return 0;
	}
	GMRFLib_sizeof_tp siz = 0;
	mapkit_size_t nelm = GMRFLib_nelm_map_ii(hash);
	siz += sizeof(map_ii) + nelm * sizeof(int);

	return siz;
}

map_id *GMRFLib_duplicate_map_id(map_id * hash)
{
	/*
	 * return a copy of HASH 
	 */
	if (!hash) {
		return NULL;
	}

	mapkit_size_t i, nelm;
	map_id *newhash = NULL;

	nelm = GMRFLib_nelm_map_id(hash);
	newhash = Calloc(1, map_id);
	map_id_init_hint(newhash, nelm);

	for (i = -1; (i = map_id_next(hash, i)) != -1;) {
		map_id_set(newhash, hash->contents[i].key, hash->contents[i].value);
	}
	newhash->alwaysdefault = hash->alwaysdefault;
	newhash->defaultvalue = hash->defaultvalue;

	return newhash;
}
GMRFLib_sizeof_tp GMRFLib_sizeof_map_id(map_id * hash)
{
	if (!hash) {
		return 0;
	}
	GMRFLib_sizeof_tp siz = 0;
	mapkit_size_t nelm = GMRFLib_nelm_map_id(hash);
	siz += sizeof(map_id) + nelm * sizeof(double);

	return siz;
}
int GMRFLib_is_int(char *str, int *value)
{
	/*
	 * return 1 if an int can be read from STR and 0 if not. return value in *VALUE is non-NULL.
	 */
	if (!str) {
		return 0;
	}

	int x, err;
	err = (sscanf(str, "%d", &x) == 1 ? 1 : 0);
	if (value) {
		*value = x;
	}
	return err;
}
char *GMRFLib_strtok_r(char *s1, const char *s2, char **lasts)
{
	char *ret;

	if (*lasts == NULL && s1 == NULL) {		       /* added this: hrue */
		return NULL;
	}

	if (s1 == NULL) {
		s1 = *lasts;
	}
	while (*s1 && strchr(s2, *s1)) {
		++s1;
	}
	if (*s1 == '\0') {
		return NULL;
	}
	ret = s1;
	while (*s1 && !strchr(s2, *s1)) {
		++s1;
	}
	if (*s1) {
		*s1++ = '\0';
	}
	*lasts = s1;

	return ret;
}

/* 
   define the floating point exception routines
 */
#if defined(__sun)||defined(__FreeBSD__)
#include <ieeefp.h>
#include <floatingpoint.h>
static void fpe_handler(void)
{
	fprintf(stderr, "\n\nFloating-point exception occured.\n");
	abort();
}
static int fpe(void)
{
	// signal(SIGTRAP, (void (*)()) fpe_handler); 
	signal(SIGFPE, (void (*)()) fpe_handler);
	// signal(SIGILL, (void (*)()) fpe_handler); 
	fpsetmask(FP_X_INV | FP_X_DZ | FP_X_OFL);
	return 0;
}
#endif
#ifdef __linux
#ifndef __USE_GNU
#define __USE_GNU 1
#endif
#include <fenv.h>
static int fpe(void)
{
	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
	return 0;
}
#endif
#if !defined(__sun) && !defined(__linux) && !defined(__FreeBSD__)
static int fpe(void)
{
	return 0;
}
#endif
int GMRFLib_fpe(void)
{
	return fpe();
}
int GMRFLib_iuniques(int *nuniques, int **uniques, int *ix, int nx)
{
	/*
	 * return in number of unique entries in ix != 0 and list them in `uniques' 
	 */

	int nu, *un = NULL, i, j, *ixx;

	if (nx <= 0 || !ix) {
		*nuniques = 0;
		if (uniques) {
			*uniques = NULL;
		}
		return GMRFLib_SUCCESS;
	}

	ixx = Calloc(nx, int);
	memcpy(ixx, ix, nx * sizeof(int));
	qsort((void *) ixx, (size_t) nx, sizeof(int), GMRFLib_icmp);

	for (j = nu = i = 0; i < nx; i++) {
		if (ixx[i] && (!i || ixx[i] != ixx[j])) {
			nu++;
			j = i;
		}
	}
	// printf("nu %d\n", nu);
	un = Calloc(nu, int);

	for (j = nu = i = 0; i < nx; i++) {
		if (ixx[i] && (!i || ixx[i] != ixx[j])) {
			// printf("\t un[%1d] = %d\n", nu, ixx[i]);
			un[nu++] = ixx[i];
			j = i;
		}
	}

	*nuniques = nu;
	if (uniques) {
		*uniques = un;
	} else {
		Free(un);
	}
	Free(ixx);

	return GMRFLib_SUCCESS;
}

#undef MEMINFO
