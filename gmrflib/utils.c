
/* utils.c
 * 
 * Copyright (C) 2006-2021 Havard Rue
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

#ifndef GITCOMMIT
#define GITCOMMIT
#endif
static const char GitID[] = "file: " __FILE__ "  " GITCOMMIT;

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

#define IDX_ALLOC_INITIAL 32
#define IDX_ALLOC_ADD     128

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

int GMRFLib_iwhich_sorted(int val, int *ix, int len)
{
	// return the index of iarray for which ix[idx]=val and
	// we KNOW that ix is sorted, and return -1 if not found

	int low, high, mid, n, n_lim = 8, i;
	if (len == 0) {
		return -1;
	}

	low = 0;
	high = len - 1;

	while (1) {
		n = high - low - 1;			       /* n is how many alternatives left */
		if (n <= n_lim) {
			for (i = low; i <= high; i++) {
				if (ix[i] == val) {
					return (i);
				}
			}
			return (-1);
		} else {
			mid = low + (high - low) / 2;	       /* integer division */
			if (ix[mid] > val) {
				high = mid;
			} else {
				low = mid;
			}
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

double GMRFLib_eps2(void)
{
	static double eps2 = -1.0;
	if (eps2 < 0.0) {
		eps2 = GMRFLib_eps(0.5);
	}
	return (eps2);
}

double GMRFLib_eps1(void)
{
	static double eps1 = -1.0;
	if (eps1 < 0.0) {
		eps1 = GMRFLib_eps(1.0);
	}
	return (eps1);
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
		val = 1.0 + eps;
		while (val > 1.0) {
			eps /= 2.0;
			val = 1.0 + eps;
		}
		eps *= 2.0;
	}
	return (pow(eps, power));
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

int GMRFLib_icmp(const void *a, const void *b)
{
	const int *ia, *ib;

	ia = (const int *) a;
	ib = (const int *) b;

	return (*ia - *ib);
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

int GMRFLib_qsorts(void *x, size_t nmemb, size_t size_x, void *y, size_t size_y, void *z, size_t size_z, int (*compar)(const void *, const void *))
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
		Memcpy((void *) &xyz[i * siz], (void *) &xx[i * size_x], size_x);
	}
	if (y) {
		offset = size_x;
		for (i = 0; i < nmemb; i++) {
			Memcpy((void *) &xyz[i * siz + offset], (void *) &yy[i * size_y], size_y);
		}
		if (z) {
			offset = size_x + size_y;
			for (i = 0; i < nmemb; i++) {
				Memcpy((void *) &xyz[i * siz + offset], (void *) &zz[i * size_z], size_z);
			}
		}
	}

	qsort((void *) xyz, nmemb, siz, compar);

	for (i = 0; i < nmemb; i++) {
		Memcpy((void *) &xx[i * size_x], (void *) &xyz[i * siz], size_x);
	}
	if (y) {
		offset = size_x;
		for (i = 0; i < nmemb; i++) {
			Memcpy((void *) &yy[i * size_y], (void *) &xyz[i * siz + offset], size_y);
		}
		if (z) {
			offset = size_x + size_y;
			for (i = 0; i < nmemb; i++) {
				Memcpy((void *) &zz[i * size_z], (void *) &xyz[i * siz + offset], size_z);
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
		return b + log1p(0.0 + a / B);
	} else {
		return log(a) + log1p(0.0 + B / a);
	}
}

void *GMRFLib_memcpy(void *dest, const void *src, size_t n)
{
	assert(n < PTRDIFF_MAX);
	memcpy(dest, src, n);
	return NULL;
}

void *GMRFLib_calloc(size_t nmemb, size_t size, const char *file, const char *funcname, int lineno, const char *id)
{
	void *ptr = NULL;
	char *msg = NULL;

	assert(nmemb * size < PTRDIFF_MAX);
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
	GMRFLib_sprintf(&msg, "Fail to calloc nmemb=%1lu elements of size=%1lu bytes", nmemb, size);
	GMRFLib_handle_error(file, funcname, lineno, id, GMRFLib_EMEMORY, msg);
	abort();

	return NULL;
}

void *GMRFLib_malloc(size_t size, const char *file, const char *funcname, int lineno, const char *id)
{
	void *ptr = NULL;
	char *msg = NULL;

	assert(size < PTRDIFF_MAX);
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
	GMRFLib_sprintf(&msg, "Fail to malloc size=%1lu bytes", size);
	GMRFLib_handle_error(file, funcname, lineno, id, GMRFLib_EMEMORY, msg);
	abort();

	return NULL;
}

void *GMRFLib_realloc(void *old_ptr, size_t size, const char *file, const char *funcname, int lineno, const char *id)
{
	void *ptr = NULL;
	char *msg = NULL;

	assert(size < PTRDIFF_MAX);
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
	GMRFLib_sprintf(&msg, "Fail to realloc size=%1lu bytes", size);
	GMRFLib_handle_error(file, funcname, lineno, id, GMRFLib_EMEMORY, msg);
	abort();

	return NULL;
}

void GMRFLib_free(void *ptr, const char *file, const char *funcname, int lineno, const char *id)
{
	if (ptr) {
#if defined(GMRFLib_MEMCHECK) && !defined(_OPENMP)
		int choice = 1;
#else
		int choice = 0;
#endif
		if (choice) {
			GMRFLib_free__(ptr, file, funcname, lineno, id);
		} else {
			free(ptr);
		}
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
	int choice = 1;
#else
	int choice = 0;
#endif

	if (choice == 1) {
		int i;
		fp = (fp ? fp : stdout);
		if (!memcheck_first) {
			for (i = -1; (i = map_vpvp_next(&memcheck_hash_table, i)) != -1;) {
				fprintf(fp, "0x%zu %s\n", (size_t) memcheck_hash_table.contents[i].key,
					(char *) memcheck_hash_table.contents[i].value);
			}
		}
	}

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

int GMRFLib_printf_matrix(FILE * fp, double *A, int m, int n)
{
	// A is m x n matrix
#pragma omp critical
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

int GMRFLib_printf_gsl_matrix(FILE * fp, gsl_matrix * matrix, const char *format)
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

int GMRFLib_printf_gsl_vector(FILE * fp, gsl_vector * vector, const char *format)
{
	size_t i;

	for (i = 0; i < vector->size; i++) {
		fprintf(fp, (format ? format : " %g"), gsl_vector_get(vector, i));
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

size_t GMRFLib_sizeof_map_ii(map_ii * hash)
{
	if (!hash) {
		return 0;
	}
	size_t siz = 0;
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

size_t GMRFLib_sizeof_map_id(map_id * hash)
{
	if (!hash) {
		return 0;
	}
	size_t siz = 0;
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
	Memcpy(ixx, ix, nx * sizeof(int));
	qsort((void *) ixx, (size_t) nx, sizeof(int), GMRFLib_icmp);

	for (j = nu = i = 0; i < nx; i++) {
		if (ixx[i] && (!i || ixx[i] != ixx[j])) {
			nu++;
			j = i;
		}
	}
	un = Calloc(nu, int);

	for (j = nu = i = 0; i < nx; i++) {
		if (ixx[i] && (!i || ixx[i] != ixx[j])) {
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

int GMRFLib_gsl_vec2plain(double **out, gsl_vector * vec)
{
	if (!vec || vec->size == 0) {
		*out = NULL;
	} else {
		*out = Calloc(vec->size, double);
		for (size_t i = 0; i < vec->size; i++) {
			(*out)[i] = gsl_vector_get(vec, i);
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_gsl_mat2plain(double **out, gsl_matrix * mat)
{
	if (!mat || mat->size1 == 0 || mat->size2 == 0) {
		*out = NULL;
	} else {
		*out = Calloc(mat->size1 * mat->size2, double);
		for (size_t j = 0; j < mat->size2; j++) {
			size_t off = j * mat->size1;
			for (size_t i = 0; i < mat->size1; i++) {
				(*out)[i + off] = gsl_matrix_get(mat, i, j);
			}
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idx_create(GMRFLib_idx_tp ** hold)
{
	*hold = Calloc(1, GMRFLib_idx_tp);
	(*hold)->idx = Calloc(IDX_ALLOC_INITIAL, int);
	(*hold)->n_alloc = IDX_ALLOC_INITIAL;
	(*hold)->n = 0;

	return GMRFLib_SUCCESS;
}

int GMRFLib_idx2_create(GMRFLib_idx2_tp ** hold)
{
	*hold = Calloc(1, GMRFLib_idx2_tp);
	(*hold)->idx = Calloc(2, int *);
	(*hold)->idx[0] = Calloc(IDX_ALLOC_INITIAL, int);
	(*hold)->idx[1] = Calloc(IDX_ALLOC_INITIAL, int);
	(*hold)->n_alloc = IDX_ALLOC_INITIAL;
	(*hold)->n = 0;

	return GMRFLib_SUCCESS;
}

int GMRFLib_val_create(GMRFLib_val_tp ** hold)
{
	*hold = Calloc(1, GMRFLib_val_tp);
	(*hold)->val = Calloc(IDX_ALLOC_INITIAL, double);
	(*hold)->n_alloc = IDX_ALLOC_INITIAL;
	(*hold)->n = 0;

	return GMRFLib_SUCCESS;
}

int GMRFLib_idxval_create(GMRFLib_idxval_tp ** hold)
{
	*hold = Calloc(1, GMRFLib_idxval_tp);
	(*hold)->store = Calloc(IDX_ALLOC_INITIAL, GMRFLib_idxval_elm_tp);
	(*hold)->n_alloc = IDX_ALLOC_INITIAL;
	(*hold)->n = 0;
	(*hold)->iaddto = 0;

	return GMRFLib_SUCCESS;
}

GMRFLib_idx_tp **GMRFLib_idx_ncreate(int n)
{
	GMRFLib_idx_tp **a = Calloc(n, GMRFLib_idx_tp *);
	for (int i = 0; i < n; i++) {
		GMRFLib_idx_create(&(a[i]));
	}
	return a;
}

GMRFLib_idx2_tp **GMRFLib_idx2_ncreate(int n)
{
	GMRFLib_idx2_tp **a = Calloc(n, GMRFLib_idx2_tp *);
	for (int i = 0; i < n; i++) {
		GMRFLib_idx2_create(&(a[i]));
	}
	return a;
}

GMRFLib_val_tp **GMRFLib_val_ncreate(int n)
{
	GMRFLib_val_tp **a = Calloc(n, GMRFLib_val_tp *);
	for (int i = 0; i < n; i++) {
		GMRFLib_val_create(&(a[i]));
	}
	return a;
}

GMRFLib_idxval_tp **GMRFLib_idxval_ncreate(int n)
{
	GMRFLib_idxval_tp **a = Calloc(n, GMRFLib_idxval_tp *);
	for (int i = 0; i < n; i++) {
		GMRFLib_idxval_create(&(a[i]));
	}
	return a;
}

int GMRFLib_idx_printf(FILE * fp, GMRFLib_idx_tp * hold, char *msg)
{
	if (hold) {
		fprintf(fp, "[%s] n = %1d  nalloc = %1d\n", msg, hold->n, hold->n_alloc);
		for (int i = 0; i < hold->n; i++) {
			fprintf(fp, "\tidx[%1d] = %1d\n", i, hold->idx[i]);
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idx2_printf(FILE * fp, GMRFLib_idx2_tp * hold, char *msg)
{
	if (hold) {
		fprintf(fp, "[%s] n = %1d  nalloc = %1d\n", msg, hold->n, hold->n_alloc);
		for (int i = 0; i < hold->n; i++) {
			fprintf(fp, "\tidx[][%1d] = %1d %1d\n", i, hold->idx[0][i], hold->idx[1][i]);
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_val_printf(FILE * fp, GMRFLib_val_tp * hold, char *msg)
{
	if (hold) {
		fprintf(fp, "[%s] n = %1d  nalloc = %1d\n", msg, hold->n, hold->n_alloc);
		for (int i = 0; i < hold->n; i++) {
			fprintf(fp, "\tval[%1d] = %g\n", i, hold->val[i]);
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idxval_printf(FILE * fp, GMRFLib_idxval_tp * hold, char *msg)
{
	if (hold) {
		fprintf(fp, "[%s] n = %1d  nalloc = %1d iaddto = %1d\n", msg, hold->n, hold->n_alloc, hold->iaddto);
		for (int i = 0; i < hold->n; i++) {
			fprintf(fp, "\tstore[%1d] = (%d, %g)\n", i, hold->store[i].idx, hold->store[i].val);
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idx_nprune(GMRFLib_idx_tp ** a, int n)
{
	if (a) {
		for (int i = 0; i < n; i++) {
			GMRFLib_idx_prune(a[i]);
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idx_sort(GMRFLib_idx_tp * hold)
{
	if (hold) {
		qsort((void *) hold->idx, (size_t) hold->n, sizeof(int), GMRFLib_icmp);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idx_nsort(GMRFLib_idx_tp ** a, int n, int nt)
{
#define CODE_BLOCK							\
	for(int i = 0; i < n; i++) {					\
		if (a[i] && a[i]->n > 1) {				\
			qsort((void *) a[i]->idx, (size_t) a[i]->n,  sizeof(int), GMRFLib_icmp); \
		}							\
	}
	RUN_CODE_BLOCK(nt);
#undef CODE_BLOCK

	return GMRFLib_SUCCESS;
}

int GMRFLib_idx_uniq(GMRFLib_idx_tp * hold)
{
	if (hold && hold->n > 1) {
		int i, j;

		GMRFLib_idx_sort(hold);
		for (j = 0, i = 0; i < hold->n; i++) {
			if (hold->idx[j] != hold->idx[i]) {
				hold->idx[++j] = hold->idx[i];
			}
		}
		hold->n = j + 1;
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idx_nuniq(GMRFLib_idx_tp ** a, int n, int nt)
{
#define CODE_BLOCK				\
	for (int i = 0; i < n; i++) {		\
		GMRFLib_idx_uniq(a[i]);		\
	}

	RUN_CODE_BLOCK(nt);
#undef CODE_BLOCK

	return GMRFLib_SUCCESS;
}

int GMRFLib_idx2_nprune(GMRFLib_idx2_tp ** a, int n)
{
	if (a) {
		for (int i = 0; i < n; i++) {
			GMRFLib_idx2_prune(a[i]);
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_val_nprune(GMRFLib_val_tp ** a, int n)
{
	if (a) {
		for (int i = 0; i < n; i++) {
			GMRFLib_val_prune(a[i]);
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idxval_nuniq(GMRFLib_idxval_tp ** a, int n, int nt)
{
#define CODE_BLOCK					\
	for (int i = 0; i < n; i++) {			\
		GMRFLib_idxval_uniq(a[i]);		\
	}
	RUN_CODE_BLOCK(nt);
#undef CODE_BLOCK
	return GMRFLib_SUCCESS;
}

int GMRFLib_idxval_uniq(GMRFLib_idxval_tp * hold)
{
	// sort idx, and accumulate values and then prune
	if (hold && hold->n > 1) {
		int i, j;

		GMRFLib_idxval_sort(hold);
		for (j = 0, i = 0; i < hold->n; i++) {
			if (hold->store[j].idx == hold->store[i].idx) {
				if (i > j) {
					hold->store[j].val += hold->store[i].val;
				}
			} else {
				j++;
				hold->store[j] = hold->store[i];
			}
		}
		hold->n = j + 1;
		GMRFLib_idxval_prune(hold);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idxval_nprune(GMRFLib_idxval_tp ** a, int n, int nt)
{
#define CODE_BLOCK					\
	for (int i = 0; i < n; i++) {			\
		GMRFLib_idxval_prune(a[i]);		\
	}

	RUN_CODE_BLOCK(nt);
#undef CODE_BLOCK

	return GMRFLib_SUCCESS;
}

int GMRFLib_idx_prune(GMRFLib_idx_tp * hold)
{
	if (hold) {
		hold->idx = Realloc(hold->idx, IMAX(1, hold->n), int);
		hold->n_alloc = IMAX(1, hold->n);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idx2_prune(GMRFLib_idx2_tp * hold)
{
	if (hold) {
		hold->idx[0] = Realloc(hold->idx[0], IMAX(1, hold->n), int);
		hold->idx[1] = Realloc(hold->idx[1], IMAX(1, hold->n), int);
		hold->n_alloc = IMAX(1, hold->n);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_val_prune(GMRFLib_val_tp * hold)
{
	if (hold) {
		hold->val = Realloc(hold->val, IMAX(1, hold->n), double);
		hold->n_alloc = IMAX(1, hold->n);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idxval_prune(GMRFLib_idxval_tp * hold)
{
	if (hold) {
		hold->store = Realloc(hold->store, IMAX(1, hold->n), GMRFLib_idxval_elm_tp);
		hold->n_alloc = IMAX(1, hold->n);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idxval_cmp(const void *a, const void *b)
{
	GMRFLib_idxval_elm_tp *aa = (GMRFLib_idxval_elm_tp *) a;
	GMRFLib_idxval_elm_tp *bb = (GMRFLib_idxval_elm_tp *) b;

	return (aa->idx - bb->idx);
}

int GMRFLib_idxval_sort(GMRFLib_idxval_tp * hold)
{
	if (hold) {
		qsort((void *) hold->store, (size_t) hold->n, sizeof(GMRFLib_idxval_elm_tp), GMRFLib_idxval_cmp);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idxval_nsort(GMRFLib_idxval_tp ** hold, int n, int nt)
{
#define CODE_BLOCK							\
	for(int i = 0; i < n; i++) {					\
		if (hold[i] && hold[i]->n > 1) {			\
			qsort((void *) hold[i]->store, (size_t) hold[i]->n,  sizeof(GMRFLib_idxval_elm_tp), GMRFLib_idxval_cmp); \
		}							\
	}

	RUN_CODE_BLOCK(nt);
#undef CODE_BLOCK

	return GMRFLib_SUCCESS;
}

int GMRFLib_idx_free(GMRFLib_idx_tp * hold)
{
	if (hold) {
		Free(hold->idx);
		Free(hold);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idx2_free(GMRFLib_idx2_tp * hold)
{
	if (hold) {
		Free(hold->idx[0]);
		Free(hold->idx[1]);
		Free(hold->idx);
		Free(hold);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_val_free(GMRFLib_val_tp * hold)
{
	if (hold) {
		Free(hold->val);
		Free(hold);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idxval_free(GMRFLib_idxval_tp * hold)
{
	if (hold) {
		Free(hold->store);
		Free(hold);
	}
	return GMRFLib_SUCCESS;
}

forceinline int GMRFLib_idx_add(GMRFLib_idx_tp ** hold, int idx)
{
	if (*hold == NULL) {
		GMRFLib_idx_create(hold);
	}
	if ((*hold)->n == (*hold)->n_alloc) {
		(*hold)->n_alloc += IMAX(IDX_ALLOC_ADD, (*hold)->n / 4);
		(*hold)->idx = Realloc((*hold)->idx, (*hold)->n_alloc, int);
	}
	(*hold)->idx[(*hold)->n] = idx;
	(*hold)->n++;

	return GMRFLib_SUCCESS;
}

forceinline int GMRFLib_idx2_add(GMRFLib_idx2_tp ** hold, int idx0, int idx1)
{
	if (*hold == NULL) {
		GMRFLib_idx2_create(hold);
	}
	if ((*hold)->n == (*hold)->n_alloc) {
		(*hold)->n_alloc += IMAX(IDX_ALLOC_ADD, (*hold)->n / 4);
		(*hold)->idx[0] = Realloc((*hold)->idx[0], (*hold)->n_alloc, int);
		(*hold)->idx[1] = Realloc((*hold)->idx[1], (*hold)->n_alloc, int);
	}
	(*hold)->idx[0][(*hold)->n] = idx0;
	(*hold)->idx[1][(*hold)->n] = idx1;
	(*hold)->n++;

	return GMRFLib_SUCCESS;
}

forceinline int GMRFLib_val_add(GMRFLib_val_tp ** hold, double val)
{
	if (*hold == NULL) {
		GMRFLib_val_create(hold);
	}
	if ((*hold)->n == (*hold)->n_alloc) {
		(*hold)->n_alloc += IMAX(IDX_ALLOC_ADD, (*hold)->n / 4);
		(*hold)->val = Realloc((*hold)->val, (*hold)->n_alloc, double);
	}
	(*hold)->val[(*hold)->n] = val;
	(*hold)->n++;

	return GMRFLib_SUCCESS;
}

forceinline int GMRFLib_idxval_add(GMRFLib_idxval_tp ** hold, int idx, double val)
{
	if (*hold == NULL) {
		GMRFLib_idxval_create(hold);
	}
	if ((*hold)->n == (*hold)->n_alloc) {
		(*hold)->n_alloc += IMAX(IDX_ALLOC_ADD, (*hold)->n / 4);
		(*hold)->store = Realloc((*hold)->store, (*hold)->n_alloc, GMRFLib_idxval_elm_tp);
	}
	(*hold)->store[(*hold)->n].idx = idx;
	(*hold)->store[(*hold)->n].val = val;
	(*hold)->n++;

	return GMRFLib_SUCCESS;
}

forceinline int GMRFLib_idxval_addto(GMRFLib_idxval_tp ** hold, int idx, double val)
{
	// if idx exists before, add val to value , otherwise just 'add'.
	// if there are two entries of 'idx', then only the first is used.

	if (*hold == NULL) {
		GMRFLib_idxval_create(hold);
	}

	int i;

	// hopefully this is it.
	i = (*hold)->iaddto;
	if ((*hold)->store[i].idx == idx) {
		(*hold)->store[i].val += val;
		return GMRFLib_SUCCESS;
	}
	// FIXME: this should be improved in general, but I think for the usage its ok. Since we are likely to add with same or increasing idx,
	// then I added this 'iaddto' which recall the last index, and try to be a little smarter.
	for (i = (*hold)->iaddto + 1; i < (*hold)->n; i++) {
		if ((*hold)->store[i].idx == idx) {
			(*hold)->store[i].val += val;
			(*hold)->iaddto = i;
			return GMRFLib_SUCCESS;
		}
	}
	for (i = 0; i < (*hold)->iaddto; i++) {
		if ((*hold)->store[i].idx == idx) {
			(*hold)->store[i].val += val;
			(*hold)->iaddto = i;
			return GMRFLib_SUCCESS;
		}
	}

	if ((*hold)->n == (*hold)->n_alloc) {
		(*hold)->n_alloc += IMAX(IDX_ALLOC_ADD, (*hold)->n / 4);
		(*hold)->store = Realloc((*hold)->store, (*hold)->n_alloc, GMRFLib_idxval_elm_tp);
	}
	(*hold)->store[(*hold)->n].idx = idx;
	(*hold)->store[(*hold)->n].val = val;
	(*hold)->n++;

	return GMRFLib_SUCCESS;
}

#undef MEMINFO
#undef IDX_ALLOC_ADD
#undef IDX_ALLOC_INITIAL

/////////////////////////////////////////////////////////////////////////

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
	int i, imin;
	double min_val;

	min_val = x[0];
	imin = 0;
	for (i = 1; i < n; i++) {
		if (x[i] < min_val) {
			min_val = x[i];
			imin = i;
		}
	}

	if (idx) {
		*idx = imin;
	}

	return min_val;
}

int GMRFLib_imin_value(int *x, int n, int *idx)
{
	/*
	 * return the IMIN(x[]) 
	 */
	int i, imin, min_val;

	min_val = x[0];
	imin = 0;
	for (i = 1; i < n; i++) {
		if (x[i] < min_val) {
			min_val = x[i];
			imin = i;
		}
	}

	if (idx) {
		*idx = imin;
	}
	return min_val;
}

double GMRFLib_max_value(double *x, int n, int *idx)
{
	/*
	 * return the MAX(x[]), optional idx
	 */
	int i, imax;
	double max_val;

	max_val = x[0];
	imax = 0;
	for (i = 1; i < n; i++) {
		if (x[i] > max_val) {
			max_val = x[i];
			imax = i;
		}
	}

	if (idx) {
		*idx = imax;
	}
	return max_val;
}

int GMRFLib_imax_value(int *x, int n, int *idx)
{
	/*
	 * return IMAX(x[]) 
	 */
	int i, imax, max_val;

	max_val = x[0];
	imax = 0;
	for (i = 1; i < n; i++) {
		if (x[i] > max_val) {
			max_val = x[i];
			imax = i;
		}
	}

	if (idx) {
		*idx = imax;
	}
	return max_val;
}

const char *GMRFLib_debug_functions_strip(const char *name)
{
	char *s = (char *) name;
	if (!strncmp("GMRFLib_", s, 8)) {
		s += 8;
	}
	if (!strncmp("inla_", s, 5)) {
		s += 5;
	}
	return s;
}

int GMRFLib_debug_functions(const char *name)
{
	static int not_defined = 0;
	if (not_defined) {
		return 0;
	}

	static int first = 1;
#pragma omp threadprivate(first)

	static map_stri *defs = NULL;
#pragma omp threadprivate(defs)

	if (first == 1) {
		// format FUN[:N],...
		// prefix's GMRFLib_ and inla_ are removed automatically
		char *def = getenv("INLA_TRACE");
		int verbose = 0;

		if (def) {
			def = GMRFLib_strdup(def);
		}
		if (verbose) {
			printf("\t\tREAD %s\n", def);
		}

		if (!def) {
			not_defined = 1;
			first = 0;
			return 0;
		} else {
			char sep1[] = ",";

			defs = Calloc(1, map_stri);
			map_stri_init_hint(defs, 128);
			char *str = def;
			char *s;

			first = -1;
			while ((s = strtok(str, sep1))) {
				str = NULL;

				int val = 0;
				char *s2 = strchr(s, ':');
				char *ss, *sss;
				if (!s2) {
					ss = s;
					val = 1;
				} else {
					int len = s2 - s + 1;
					ss = Calloc(len + 1, char);
					ss[len] = '\0';
					strncpy(ss, s, len - 1);
					val = atoi(s2 + 1);
					val = IMAX(val, 1);
				}
				// strip leading whitespace
				while (!strncmp(ss, " ", 1))
					ss++;
				// special option that override all others
				if (!strcmp(ss, "*")) {
					first = 2;
				}

				sss = (char *) GMRFLib_debug_functions_strip((const char *)ss);
				char *nm = NULL;
				if (strlen(ss)) {
					GMRFLib_sprintf(&nm, "%s", sss);
					map_stri_set(defs, nm, val);
				}
				if (first != 2) {
					first = 0;
				}

				if (verbose) {
					printf("\t\t[%1d] debug init: ADD [%s]=%1d\n", omp_get_thread_num(), sss, val);
				}
			}
		}
	}

	int *p = map_stri_ptr(defs, (char *) (first == 2 ? "*" : GMRFLib_debug_functions_strip(name)));

	return (p ? *p : 0);
}
