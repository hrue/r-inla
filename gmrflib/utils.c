
/* utils.c
 * 
 * Copyright (C) 2006-2022 Havard Rue
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
#include <assert.h>
#include <float.h>
#include <signal.h>
#include <stdarg.h>
#include <math.h>
#include <strings.h>
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stddef.h>
#include <string.h>
#include <malloc.h>
#endif
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"
#include "GMRFLib/hashP.h"

#define IDX_ALLOC_INITIAL 64
#define IDX_ALLOC_ADD     512

void GMRFLib_delay(int msec)
{
	long pause;
	clock_t now, then;

	pause = msec * (CLOCKS_PER_SEC / 1000);
	now = then = clock();
	while ((now - then) < pause) {
		now = clock();
	}
}

void GMRFLib_delay_random(int msec_low, int msec_high)
{
	GMRFLib_delay(msec_low + (int) ((msec_high - msec_low) * GMRFLib_uniform()));
}

int GMRFLib_sprintf(char **ptr, const char *fmt, ...)
{
	/*
	 * parts of this code is copied from the manual page of snprintf. 
	 */

	int n, size = 128 + 1;
	char *p;
	va_list ap;

	GMRFLib_ASSERT(ptr, GMRFLib_EINVARG);
	GMRFLib_ASSERT(fmt, GMRFLib_EINVARG);

	p = Calloc(size, char);

	while (1) {
		/*
		 * Try to print in the allocated space. 
		 */
		va_start(ap, fmt);
		n = vsnprintf(p, (unsigned int) size, fmt, ap);
		va_end(ap);

		/*
		 * if that worked, return the string, 
		 */
		if (n > -1 && n < size) {
			*ptr = p;
			return GMRFLib_SUCCESS;
		}

		/*
		 * ...else try again with more space 
		 */
		if (n > -1) {
			size = n + 1;
		} else {
			size *= 2;
		}
		p = Realloc(p, size, char);
	}

	return GMRFLib_SUCCESS;
}

void *GMRFLib_memcpy(void *dest, const void *src, size_t n)
{
	// assert(n < PTRDIFF_MAX);
	memcpy(dest, src, n);
	return NULL;
}

void *GMRFLib_calloc(size_t nmemb, size_t size, const char *file, const char *funcname, int lineno, const char *id)
{
	void *ptr = NULL;
	char *msg = NULL;

	assert(nmemb * size < PTRDIFF_MAX);
	ptr = calloc(nmemb, size);

	if (ptr) {
		return ptr;
	}
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
	ptr = malloc(size);
	if (ptr) {
		return ptr;
	}
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
	ptr = realloc(old_ptr, size);
	if (ptr) {
		return ptr;
	}
	GMRFLib_sprintf(&msg, "Fail to realloc size=%1lu bytes", size);
	GMRFLib_handle_error(file, funcname, lineno, id, GMRFLib_EMEMORY, msg);
	abort();

	return NULL;
}

void GMRFLib_free(void *ptr, const char *file, const char *funcname, int lineno, const char *id)
{
	if (ptr) {
		free(ptr);
	} else {
		fprintf(stderr, "%s:%s:%d (%s): Try to free a NULL-ptr\n", file, funcname, lineno, id);
	}
}

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

int GMRFLib_iwhich_sorted(int val, int *ix, int len, int *guess)
{
	// return the index of iarray for which ix[idx]=val and we KNOW that ix is sorted, and return -1 if not found. 'guess' (NULL is not
	// allowed) is an initial guess for [low,high] and automatically updated. initialize with guess[1]=0. 'guess' must be thread-safe

	if (len == 0) {
		return -1;
	}

	int low, high, mid, n, lguess[2] = { 0, 0 };

	if (!guess) {
		guess = lguess;
	}
	// use the guess of [low,high] ? MUST BE INITIALIZED to [0,0]!
	if (guess[1] == 0 || guess[1] >= len || guess[1] <= guess[0]) {
		// invalid values for 'guess', no need to check
		low = 0;
		high = len - 1;
	} else {
		low = (val >= ix[guess[0]] ? guess[0] : 0);
		high = (val <= ix[guess[1]] ? guess[1] : len - 1);
	}

	while (1) {
		n = high - low + 1;			       // 'n' is how many alternatives left 
		if (n <= 8) {
			guess[1] = high + 1;
			for (int i = low; i <= high; i++) {
				if (ix[i] == val) {
					guess[0] = i + 1;
					return i;
				}
			}
			guess[0] = low;
			guess[1] = high;
			return -1;
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

double GMRFLib_eps(double power)
{
	/*
	 * Return eps^power, where eps is the smalles number such that 1+eps != eps.  
	 */

	return (pow(DBL_EPSILON, power));
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

int GMRFLib_icmp_r(const void *a, const void *b)
{
	return (-GMRFLib_icmp(a, b));
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
	return (-GMRFLib_dcmp(a, b));
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

int GMRFLib_dcmp_abs_r(const void *a, const void *b)
{
	return (-GMRFLib_dcmp_abs(a, b));
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

char *GMRFLib_strdup(const char *ptr)
{
	/*
	 * strdup is not part of ANSI-C 
	 */
#if defined(__STRICT_ANSI__)
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

int GMRFLib_normalize(int n, double *x)
{
	// scale x so the sum is 1

	double sum = 0.0;
	for (int i = 0; i < n; i++) {
		sum += x[i];
	}
	sum = 1.0 / sum;
	for (int i = 0; i < n; i++) {
		x[i] *= sum;
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
#pragma omp critical (Name_bb051132870d1f0b90133946052e91194aa163a5)
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

int GMRFLib_gsl_matrix_count_eq(gsl_matrix * A, double value)
{
	int num = 0;
	for (size_t i = 0; i < A->size1; i++) {
		for (size_t j = 0; j < A->size2; j++) {
			num += (ISNAN(value) ? ISNAN(gsl_matrix_get(A, i, j)) : (gsl_matrix_get(A, i, j) == value));
		}
	}
	return num;
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

forceinline double GMRFLib_logit(double p)
{
	// evaluate log(p/(1-p)) more safe than just log(p/(1-p))
	const double lim = 0.01;

	if (p > lim && p < 1.0 - lim) {
		return log(p / (1.0 - p));
	} else if (p < 0.5) {
		return (log(p) - log1p(-p));
	} else {
		double pp = 1.0 - p;
		return (-log(pp) + log1p(-pp));
	}
}

forceinline double GMRFLib_inv_logit(double x)
{
	// evaluate 1/(1+exp(-x))

	return 1.0 / (2.0 + expm1(-x));
}

const char *GMRFLib_function_name_strip(const char *name)
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

	static map_stri **ddefs = NULL;
	static int *first = NULL;

	if (!ddefs) {
#pragma omp critical (Name_30c48b516c7b1cce1be137af0e429a5e3b52a645)
		{
			if (!ddefs) {
				first = Calloc(GMRFLib_CACHE_LEN, int);
				ddefs = Calloc(GMRFLib_CACHE_LEN, map_stri *);
			}
		}
	}
	int idx = 0;
	GMRFLib_CACHE_SET_ID(idx);

	if (!ddefs[idx]) {
		// format FUN[:N],...
		// prefix's GMRFLib_ and inla_ are removed automatically
		char *def = getenv("INLA_DEBUG");
		int verbose = 0;

		if (def) {
			def = GMRFLib_strdup(def);
		}
		if (verbose) {
			printf("\t\tREAD %s\n", def);
		}

		if (!def) {
			not_defined = 1;
			return 0;
		} else {
			char sep1[] = ",";

			ddefs[idx] = Calloc(1, map_stri);
			map_stri_init_hint(ddefs[idx], 128);
			char *str = def;
			char *s;

			first[idx] = -1;
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
					first[idx] = 2;
				}

				sss = (char *) GMRFLib_function_name_strip((const char *) ss);
				char *nm = NULL;
				if (strlen(ss)) {
					GMRFLib_sprintf(&nm, "%s", sss);
					map_stri_set(ddefs[idx], nm, val);
				}
				if (first[idx] != 2) {
					first[idx] = 0;
				}

				if (verbose) {
					printf("\t\t[%1d] debug init: ADD [%s]=%1d\n", omp_get_thread_num(), sss, val);
				}
			}
		}
	}

	if (!name) {
		return 0;
	} else {
		int *p = map_stri_ptr(ddefs[idx], (char *) (first[idx] == 2 ? "*" : GMRFLib_function_name_strip(name)));
		return (p ? *p : 0);
	}
}

int GMRFLib_trace_functions(const char *name)
{
	static int not_defined = 0;
	if (not_defined) {
		return 0;
	}

	static map_stri **ddefs = NULL;
	static int *first = NULL;

	if (!ddefs) {
#pragma omp critical (Name_3a266edf254a33111bcf4ab49b3acc5833850a29)
		{
			if (!ddefs) {
				first = Calloc(GMRFLib_CACHE_LEN, int);
				ddefs = Calloc(GMRFLib_CACHE_LEN, map_stri *);
			}
		}
	}
	int idx = 0;
	GMRFLib_CACHE_SET_ID(idx);

	if (!ddefs[idx]) {
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
			return 0;
		} else {
			char sep1[] = ",";

			ddefs[idx] = Calloc(1, map_stri);
			map_stri_init_hint(ddefs[idx], 128);
			char *str = def;
			char *s;

			first[idx] = -1;
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
					first[idx] = 2;
				}

				sss = (char *) GMRFLib_function_name_strip((const char *) ss);
				char *nm = NULL;
				if (strlen(ss)) {
					GMRFLib_sprintf(&nm, "%s", sss);
					map_stri_set(ddefs[idx], nm, val);
				}
				if (first[idx] != 2) {
					first[idx] = 0;
				}

				if (verbose) {
					printf("\t\t[%1d] debug init: ADD [%s]=%1d\n", omp_get_thread_num(), sss, val);
				}
			}
		}
	}

	if (!name) {
		return 0;
	} else {
		int *p = map_stri_ptr(ddefs[idx], (char *) (first[idx] == 2 ? "*" : GMRFLib_function_name_strip(name)));
		return (p ? *p : 0);
	}
}

/// ******************************************************************
/// ******************************************************************
/// ******************************************************************

int GMRFLib_idx_create(GMRFLib_idx_tp ** hold)
{
	return GMRFLib_idx_create_x(hold, IDX_ALLOC_INITIAL);
}

int GMRFLib_idx_create_x(GMRFLib_idx_tp ** hold, int len)
{
	len = IMAX(1, len);
	*hold = Calloc(1, GMRFLib_idx_tp);
	(*hold)->idx = Calloc(len, int);
	(*hold)->n_alloc = len;
	(*hold)->n = 0;

	return GMRFLib_SUCCESS;
}

int GMRFLib_str_create(GMRFLib_str_tp ** hold)
{
	return GMRFLib_str_create_x(hold, IDX_ALLOC_INITIAL);
}

int GMRFLib_str_create_x(GMRFLib_str_tp ** hold, int len)
{
	len = IMAX(1, len);
	*hold = Calloc(1, GMRFLib_str_tp);
	(*hold)->str = Calloc(len, char *);
	(*hold)->n_alloc = len;
	(*hold)->n = 0;

	return GMRFLib_SUCCESS;
}

int GMRFLib_idx2_create(GMRFLib_idx2_tp ** hold)
{
	return GMRFLib_idx2_create_x(hold, IDX_ALLOC_INITIAL);
}

int GMRFLib_idx2_create_x(GMRFLib_idx2_tp ** hold, int len)
{
	len = IMAX(1, len);
	*hold = Calloc(1, GMRFLib_idx2_tp);
	(*hold)->idx = Calloc(2, int *);
	(*hold)->idx[0] = Calloc(len, int);
	(*hold)->idx[1] = Calloc(len, int);
	(*hold)->n_alloc = len;
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
	return GMRFLib_idxval_create_x(hold, IDX_ALLOC_INITIAL);
}

int GMRFLib_idxval_create_x(GMRFLib_idxval_tp ** hold, int len)
{
	len = IMAX(1, len);
	*hold = Calloc(1, GMRFLib_idxval_tp);
	(*hold)->idx = Calloc(len, int);
	(*hold)->val = Calloc(len, double);
	(*hold)->n_alloc = len;
	(*hold)->n = 0;
	(*hold)->iaddto = 0;

	return GMRFLib_SUCCESS;
}

GMRFLib_idx_tp **GMRFLib_idx_ncreate(int n)
{
	if (n > 0) {
		GMRFLib_idx_tp **a = Calloc(n, GMRFLib_idx_tp *);
		for (int i = 0; i < n; i++) {
			GMRFLib_idx_create(&(a[i]));
		}
		return a;
	} else {
		return NULL;
	}
}

GMRFLib_idx_tp **GMRFLib_idx_ncreate_x(int n, int len)
{
	if (n > 0) {
		GMRFLib_idx_tp **a = Calloc(n, GMRFLib_idx_tp *);
		for (int i = 0; i < n; i++) {
			GMRFLib_idx_create_x(&(a[i]), len);
		}
		return a;
	} else {
		return NULL;
	}
}

GMRFLib_str_tp **GMRFLib_str_ncreate(int n)
{
	if (n > 0) {
		GMRFLib_str_tp **a = Calloc(n, GMRFLib_str_tp *);
		for (int i = 0; i < n; i++) {
			GMRFLib_str_create(&(a[i]));
		}
		return a;
	} else {
		return NULL;
	}
}

GMRFLib_str_tp **GMRFLib_str_ncreate_x(int n, int len)
{
	if (n > 0) {
		GMRFLib_str_tp **a = Calloc(n, GMRFLib_str_tp *);
		for (int i = 0; i < n; i++) {
			GMRFLib_str_create_x(&(a[i]), len);
		}
		return a;
	} else {
		return NULL;
	}
}

GMRFLib_idx2_tp **GMRFLib_idx2_ncreate(int n)
{
	if (n > 0) {
		GMRFLib_idx2_tp **a = Calloc(n, GMRFLib_idx2_tp *);
		for (int i = 0; i < n; i++) {
			GMRFLib_idx2_create(&(a[i]));
		}
		return a;
	} else {
		return NULL;
	}
}

GMRFLib_idx2_tp **GMRFLib_idx2_ncreate_x(int n, int len)
{
	if (n > 0) {
		GMRFLib_idx2_tp **a = Calloc(n, GMRFLib_idx2_tp *);
		for (int i = 0; i < n; i++) {
			GMRFLib_idx2_create_x(&(a[i]), len);
		}
		return a;
	} else {
		return NULL;
	}
}

GMRFLib_val_tp **GMRFLib_val_ncreate(int n)
{
	if (n > 0) {
		GMRFLib_val_tp **a = Calloc(n, GMRFLib_val_tp *);
		for (int i = 0; i < n; i++) {
			GMRFLib_val_create(&(a[i]));
		}
		return a;
	} else {
		return NULL;
	}
}

GMRFLib_idxval_tp **GMRFLib_idxval_ncreate(int n)
{
	if (n > 0) {
		GMRFLib_idxval_tp **a = Calloc(n, GMRFLib_idxval_tp *);
		for (int i = 0; i < n; i++) {
			GMRFLib_idxval_create(&(a[i]));
		}
		return a;
	} else {
		return NULL;
	}
}

GMRFLib_idxval_tp **GMRFLib_idxval_ncreate_x(int n, int len)
{
	if (n > 0) {
		GMRFLib_idxval_tp **a = Calloc(n, GMRFLib_idxval_tp *);
		for (int i = 0; i < n; i++) {
			GMRFLib_idxval_create_x(&(a[i]), len);
		}
		return a;
	} else {
		return NULL;
	}
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

int GMRFLib_str_printf(FILE * fp, GMRFLib_str_tp * hold, char *msg)
{
	if (hold) {
		fprintf(fp, "[%s] n = %1d  nalloc = %1d\n", msg, hold->n, hold->n_alloc);
		for (int i = 0; i < hold->n; i++) {
			fprintf(fp, "\tidx[%1d] = %s\n", i, hold->str[i]);
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
			fprintf(fp, "\t(idx, val)[%1d] = (%d, %g)\n", i, hold->idx[i], hold->val[i]);
		}
		if (hold->g_i) {
			fprintf(fp, "\tg_n = %1d\n", hold->g_n);
			for (int g = 0; g < hold->g_n; g++) {
				fprintf(fp, "\tgroup %d has length %d and start at index %d (one=%s)\n", g, hold->g_len[g], hold->g_i[g],
					(hold->g_1 && hold->g_1[g] ? "TRUE" : "FALSE"));
				fprintf(fp, "\t\t");
				for (int k = 0; k < IABS(hold->g_len[g]); k++) {
					fprintf(fp, " %1d", hold->idx[hold->g_i[g] + k]);
				}
				fprintf(fp, "\n");
			}
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idxval_info_printf(FILE * fp, GMRFLib_idxval_tp * hold, char *msg)
{
	if (hold) {
		fprintf(fp, "[%s] n=%1d  nalloc=%1d iaddto=%1d\n", (msg ? msg : ""), hold->n, hold->n_alloc, hold->iaddto);
		if (hold->g_i) {
			for (int g = 0; g < hold->g_n; g++) {
				fprintf(fp, "\tgroup %1d/%1d: len=%1d one=%s\n", g, hold->g_n, hold->g_len[g], 
					(hold->g_1 && hold->g_1[g] ? "TRUE" : "FALSE"));
			}
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

int GMRFLib_str_nprune(GMRFLib_str_tp ** a, int n)
{
	if (a) {
		for (int i = 0; i < n; i++) {
			GMRFLib_str_prune(a[i]);
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

	RUN_CODE_BLOCK(nt, 0, 0);
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

	RUN_CODE_BLOCK(nt, 0, 0);
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

	RUN_CODE_BLOCK(nt, 0, 0);
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
			if (hold->idx[j] == hold->idx[i]) {
				if (i > j) {
					hold->val[j] += hold->val[i];
				}
			} else {
				j++;
				hold->idx[j] = hold->idx[i];
				hold->val[j] = hold->val[i];
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

	RUN_CODE_BLOCK(nt, 0, 0);
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

int GMRFLib_str_prune(GMRFLib_str_tp * hold)
{
	if (hold) {
		hold->str = Realloc(hold->str, IMAX(1, hold->n), char *);
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
		int n = IMAX(1, hold->n);
		hold->idx = Realloc(hold->idx, n, int);
		hold->val = Realloc(hold->val, n, double);
		hold->n_alloc = n;
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idxval_sort(GMRFLib_idxval_tp * hold)
{
	return GMRFLib_idxval_nsort(&hold, 1, 1);
}

int GMRFLib_idxval_nsort(GMRFLib_idxval_tp ** hold, int n, int nt)
{
	return GMRFLib_idxval_nsort_x(hold, n, nt, 0);
}

int GMRFLib_idxval_nsort_x(GMRFLib_idxval_tp ** hold, int n, int nt, int prune_zeros)
{
	GMRFLib_ENTER_ROUTINE;

	// prune_zeros, 0=do nothing, >0 remove duplicate with same index within each group, <0 remove all zeros
	int debug = GMRFLib_DEBUG_IF_TRUE();
	
	int nmax = 1;
	for (int i = 0; i < n; i++) {
		GMRFLib_idxval_tp *h = hold[i];
		nmax = IMAX(nmax, h->n);
	}

#define CODE_BLOCK							\
	for (int i = 0; i < n; i++) {					\
		GMRFLib_idxval_tp *h = hold[i];				\
		if (!h) continue;					\
		double tref[5] =  {0, 0, 0, 0, 0};			\
		if (debug) tref[0] -= GMRFLib_cpu();			\
		if (h->n > 1) {						\
			if (0) {					\
				GMRFLib_qsorts((void *) h->idx, (size_t) h->n, sizeof(int), \
					       (void *) h->val, sizeof(double), NULL, (size_t) 0, GMRFLib_icmp); \
			} else {					\
				double *idx_tmp = CODE_BLOCK_WORK_PTR(0); \
				for(int j = 0; j < h->n; j++) {		\
					idx_tmp[j] = (double) h->idx[j]; \
				}					\
				gsl_sort2(idx_tmp, (size_t) 1, h->val, (size_t) 1, (size_t) h->n); \
				for(int j = 0; j < h->n; j++) {		\
					h->idx[j] = (int) idx_tmp[j];	\
				}					\
			}						\
		}							\
		if (debug) tref[0] += GMRFLib_cpu();			\
									\
		/*							\
		 * build basic groups with one group for each sequence and then one for each individual \
		 */							\
		if (debug) tref[1] -= GMRFLib_cpu();			\
		int ng = 1;						\
		for (int j = 1; j < h->n; j++) {			\
			if (h->idx[j] != h->idx[j - 1] + 1) {		\
				ng++;					\
			}						\
		}							\
		int *g_i = Calloc(ng + 1, int);				\
		int *g_len = Calloc(ng + 1, int);			\
									\
		int k = 0;						\
		g_i[0] = 0;						\
		for (int j = 1; j < h->n; j++) {			\
			if (h->idx[j] != h->idx[j - 1] + 1) {		\
				g_len[k] = j - g_i[k];			\
				k++;					\
				g_i[k] = j;				\
			}						\
		}							\
		g_len[ng - 1] = h->n - g_i[ng - 1];			\
		if (debug) tref[1] += GMRFLib_cpu();			\
									\
		if (debug) {						\
			GMRFLib_idxval_printf(stdout, h, "OLD");	\
		}							\
									\
		k = h->n;				       /* so the code below works */ \
		int kg = ng;				       /* so the code below works */ \
									\
		if (debug) tref[2] -= GMRFLib_cpu();			\
		if (prune_zeros) {					\
			/*						\
			 * remove zeros in 'val', either all or all but the first one \
			 */						\
			k = 0;						\
			kg = 0;						\
			int nskip = 0;					\
			for (int g = 0; g < ng; g++) {			\
									\
				int gl = IABS(g_len[g]);		\
				int num_zero = 0;			\
				int idx_zz = -1;			\
				int gn = 0;				\
				int kk = k;				\
									\
				for (int j = 0; j < gl; j++) {		\
					int idx_z = h->idx[g_i[g] + j];	\
									\
					if (prune_zeros < 0) {		\
						if (ISZERO(h->val[g_i[g] + j])) { \
							nskip++;	\
						} else {		\
							h->idx[k] = idx_z; \
							h->val[k] = h->val[g_i[g] + j];	\
							k++;		\
							gn++;		\
						}			\
					} else {			\
						if (idx_z == idx_zz) {	\
							if (ISZERO(h->val[g_i[g] + j])) { \
								num_zero++; \
								if (prune_zeros > 0 && num_zero == 1) {	\
									h->idx[k] = idx_z; \
									h->val[k] = h->val[g_i[g] + j];	\
									k++; \
									gn++; \
								} else { \
									nskip++; \
								}	\
							} else {	\
								h->idx[k] = idx_z; \
								h->val[k] = h->val[g_i[g] + j];	\
								k++;	\
								gn++;	\
							}		\
						} else {		\
							h->idx[k] = idx_z; \
							h->val[k] = h->val[g_i[g] + j];	\
							k++;		\
							gn++;		\
							num_zero = 0;	\
							if (ISZERO(h->val[g_i[g] + j])) { \
								num_zero++; \
							}		\
						}			\
						idx_zz = idx_z;		\
					}				\
				}					\
				if (debug && (g_len[g] != gn)) {	\
					printf("modify group %d from %d to %d, nskip=%1d\n", g, gl, gn, nskip); \
				}					\
				if (gn) {				\
					g_i[kg] = kk;			\
					g_len[kg] = (g_len[g] < 0 ? -gn : gn); \
					kg++;				\
				}					\
			}						\
		}							\
		h->n = k;						\
		h->g_n = kg;						\
		h->g_i = g_i;						\
		h->g_len = g_len;					\
									\
		if (debug) tref[2] += GMRFLib_cpu();			\
									\
		if (debug) {						\
			GMRFLib_idxval_printf(stdout, h, "NEW");	\
		}							\
									\
		/*							\
		 * grow the groups together by merging individual groups, and merging sequential groups \
		 */							\
		if (debug) tref[3] -= GMRFLib_cpu();			\
		int gg = 0;						\
		int irregular = 0;					\
		for (int g = 0; g < h->g_n; g++) {			\
			if (h->g_len[g] == 1) {				\
				if (irregular) {			\
					h->g_len[gg]++;			\
				} else {				\
					h->g_len[gg] = 1;		\
					h->g_i[gg] = h->g_i[g];		\
				}					\
				irregular = 1;				\
			} else if (IABS(h->g_len[g]) > 0) {		\
				if (irregular) {			\
					gg++;				\
				}					\
				irregular = 0;				\
				h->g_len[gg] = -IABS(h->g_len[g]);	\
				h->g_i[gg] = h->g_i[g];			\
				gg++;					\
			} else {					\
				continue;				\
			}						\
		}							\
		h->g_n = gg + irregular;				\
									\
		int fill = 8L;						\
		if (h && h->n > 1) {					\
			int len_extra = 0;				\
			for (int g = 0; g < h->g_n - 1; g++) {		\
				if (h->g_len[g] < 0 && h->g_len[g + 1] < 0) { \
					int last_this = h->idx[h->g_i[g] + IABS(h->g_len[g]) - 1]; \
					int first_next = h->idx[h->g_i[g + 1]];	\
					if (first_next - last_this < fill) { \
						len_extra += first_next - last_this - 1; \
						if (debug) {		\
							printf("\tmerge group %1d and %1d with dist %1d and lengths %1d and %1d\n", \
							       g, g + 1, first_next - last_this, h->g_len[g], h->g_len[g + 1]);	\
						}			\
					}				\
				}					\
			}						\
									\
			if (len_extra > 0) {				\
									\
				int n_new = h->n + len_extra;		\
				int *idx_new = Calloc(n_new, int);	\
				double *val_new = Calloc(n_new, double); \
									\
				int kk = 0;				\
				int gg = 0;				\
									\
				int this_len = h->g_len[gg];		\
				int this_alen = IABS(this_len);		\
									\
				Memcpy(idx_new + kk, h->idx, this_alen * sizeof(int)); \
				Memcpy(val_new + kk, h->val, this_alen * sizeof(double)); \
				kk += this_alen;			\
				h->g_i[gg] = h->g_i[0];			\
				h->g_len[gg] = h->g_len[0];		\
									\
				int gg_len = h->g_len[gg];		\
				int gg_alen = IABS(gg_len);		\
				int gg_first_i = h->g_i[gg];		\
				int gg_last_i = gg_first_i + gg_alen - 1; \
				int gg_last_idx = idx_new[gg_last_i];	\
									\
				int pending = 0;			\
									\
				if (debug) {				\
					printf("n=%1d len_extra=%1d n_new=%1d\n", h->n, len_extra, n_new); \
				}					\
									\
				for (int g = 1; g < h->g_n; g++) {	\
					if (debug) {			\
						printf("\tprocess group=%1d with len=%d kk=%1d\n", g, h->g_len[g], kk);	\
					}				\
									\
					/*				\
					 * merge this group into the current new one or just add it ? \
					 */				\
					int g_len = h->g_len[g];	\
					int g_alen = IABS(g_len);	\
					int g_first_i = h->g_i[g];	\
					int g_first_idx = h->idx[g_first_i]; \
									\
					int idx_diff = g_first_idx - gg_last_idx; \
					if (debug) {			\
						printf("g_len=%1d gg_len=%1d idx_diff=%1d\n", g_len, gg_len, idx_diff);	\
					}				\
									\
					if (g_len < 0 && gg_len < 0 && idx_diff < fill) { \
						if (debug) {		\
							printf("\tmerge group g=%1d into gg=%1d\n", g, gg); \
						}			\
									\
						/*			\
						 * yes, merge it	\
						 */			\
						int len_fill = idx_diff - 1; \
						for (int j = 0; j < len_fill; j++) { \
							idx_new[kk + j] = idx_new[kk - 1] + 1 + j; \
						}			\
						if (debug) {		\
							printf("\t\tlen_fill=%1d\n", len_fill);	\
						}			\
									\
						h->g_len[gg] -= len_fill; \
						kk += len_fill;		\
						Memcpy(idx_new + kk, h->idx + g_first_i, g_alen * sizeof(int));	\
						Memcpy(val_new + kk, h->val + g_first_i, g_alen * sizeof(double)); \
						h->g_len[gg] -= g_alen;	\
						kk += g_alen;		\
									\
						gg_len = h->g_len[gg];	\
						gg_alen = IABS(gg_len);	\
						gg_first_i = h->g_i[gg]; \
						gg_last_i = gg_first_i + gg_alen - 1; \
						gg_last_idx = idx_new[gg_last_i]; \
									\
						if (debug) {		\
							printf("\t\tappend group g=%1d with len %1d\n", g, g_len); \
							printf("\t\tgg_last_idx=%1d kk=%1d\n", gg_last_idx, kk); \
						}			\
						pending = 1;		\
					} else {			\
						/*			\
						 * no, just add a new group \
						 */			\
						gg++;			\
						Memcpy(idx_new + kk, h->idx + g_first_i, g_alen * sizeof(int));	\
						Memcpy(val_new + kk, h->val + g_first_i, g_alen * sizeof(double)); \
						h->g_len[gg] = g_len;	\
						h->g_i[gg] = kk;	\
						kk += g_alen;		\
									\
						if (debug) {		\
							printf("\tmake a new group, gg=%1d with len=%1d\n", gg, g_len);	\
						}			\
									\
						gg_len = h->g_len[gg];	\
						gg_alen = IABS(gg_len);	\
						gg_first_i = h->g_i[gg]; \
						gg_last_i = gg_first_i + gg_alen - 1; \
						gg_last_idx = idx_new[gg_last_i]; \
						pending = (g < h->g_n - 1 ? 0 : 1); \
					}				\
					if (debug) {			\
						printf("gg=%1d kk=%1d\n", gg, kk); \
					}				\
				}					\
									\
				Free(h->idx);				\
				Free(h->val);				\
				h->n = n_new;				\
				h->g_n = gg + pending;			\
				h->idx = idx_new;			\
				h->val = val_new;			\
				h->n_alloc = n_new;			\
			}						\
		}							\
		if (debug) tref[3] += GMRFLib_cpu();			\
									\
		/*							\
		 * add a boolean for all(val[]==1), which makes dot-products into sums \
		 */							\
		if (debug) tref[4] -= GMRFLib_cpu();			\
		int *g_1 = Calloc(h->g_n + 1, int);			\
		for (int g = 0; g < h->g_n; g++) {			\
			int all_1 = 1;					\
			for (int j = 0; j < IABS(h->g_len[g]) && all_1; j++) { \
				all_1 = all_1 && ISEQUAL(h->val[h->g_i[g] + j], 1.0); \
			}						\
			g_1[g] = (all_1 ? 1 : 0);			\
		}							\
		h->g_1 = g_1;						\
		if (debug) tref[4] += GMRFLib_cpu();			\
									\
		if (debug) {						\
			GMRFLib_idxval_info_printf(stdout, h, "\t");	\
		}							\
		if (debug) {						\
			double sum = tref[0] + tref[1] + tref[2] + tref[3] + tref[4]; \
			printf("%.3f %.3f %.3f %.3f %.3f\n", tref[0]/sum, tref[1]/sum, tref[2]/sum, tref[3]/sum, tref[4]/sum); \
		}							\
	}

	RUN_CODE_BLOCK(nt, 1, nmax);						
#undef CODE_BLOCK

	/* 
	 * Add a tag about which is faster, the group or the serial algorithm, for each 'idxval'
	 */
	
	nmax = 1;
	for (int i = 0; i < n; i++) {
		GMRFLib_idxval_tp *h = hold[i];
		if (h->n > 0) {
			nmax = IMAX(nmax, h->idx[h->n - 1]);
		}
	}

	double *x = Calloc(nmax, double);
	for(int i = 0; i < nmax; i++) {
		x[i] = GMRFLib_uniform();
	}
	
	double time_min = 0.0;
	double time_max = 0.0;
	int ntimes = 2;
	
	for(int i = 0; i < n; i++) {
		double tref[2] = {0.0, 0.0};
		double value[2] = {0.0, 0.0};
		
		if (debug) printf("start testing for hold[%1d]...\n", i);
		for(int time = -1; time < ntimes; time++) {
			if (time >= 0) {
				tref[0] -= GMRFLib_cpu();
			}
			DOT_PRODUCT_SERIAL(value[0], hold[i], x);
			if (time >= 0) {
				tref[0] += GMRFLib_cpu();
			}
			
			if (time >= 0) {
				tref[1] -= GMRFLib_cpu();
			}
			DOT_PRODUCT_GROUP(value[1], hold[i], x);
			if (time >= 0) {
				tref[1] += GMRFLib_cpu();
			}
		}

		if ((ABS(value[1] - value[0]) / (1.0 + (ABS(value[0]) + ABS(value[1]))/2.0)) > 1.0E-6) {
			P(value[0]);
			P(value[1]);
			assert(0 == 1);
		}
		if (debug) {
			printf("for h[%1d] with n= %1d : time serial= %.3f  group= %.3f\n", i, hold[i]->n,
			       tref[0] / (tref[0] + tref[1]), tref[1] / (tref[0] + tref[1]));
		}

		hold[i]->preference = (tref[0] < tref[1] ? IDXVAL_SERIAL : IDXVAL_GROUP);

		time_min += DMIN(tref[0], tref[1]) / ntimes;
		time_max += DMAX(tref[0], tref[1]) / ntimes;
	}

	if (debug) {
		printf("idxval opt: saving %.6f seconds/M.eval, %.2f%% improvement\n",
		       (time_max - time_min) * ISQR(1024), 100.0 * (1.0 - time_min / time_max));
	}
	
	Free(x);
									
GMRFLib_LEAVE_ROUTINE;							
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

int GMRFLib_str_free(GMRFLib_str_tp * hold)
{
	if (hold) {
		for (int i = 0; i < hold->n; i++) {
			if (hold->str[i]) {
				Free(hold->str[i]);
			}
		}
		Free(hold->str);
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
		Free(hold->idx);
		Free(hold->val);
		Free(hold->g_i);
		Free(hold->g_len);
		Free(hold->g_1);
		Free(hold);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_idx_add(GMRFLib_idx_tp ** hold, int idx)
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

int GMRFLib_str_add(GMRFLib_str_tp ** hold, char *s)
{
	if (*hold == NULL) {
		GMRFLib_str_create(hold);
	}
	if ((*hold)->n == (*hold)->n_alloc) {
		(*hold)->n_alloc += IMAX(IDX_ALLOC_ADD, (*hold)->n / 4);
		(*hold)->str = Realloc((*hold)->str, (*hold)->n_alloc, char *);
	}
	(*hold)->str[(*hold)->n] = GMRFLib_strdup(s);
	(*hold)->n++;

	return GMRFLib_SUCCESS;
}

int GMRFLib_idx2_add(GMRFLib_idx2_tp ** hold, int idx0, int idx1)
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

int GMRFLib_val_add(GMRFLib_val_tp ** hold, double val)
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

int GMRFLib_idxval_add(GMRFLib_idxval_tp ** hold, int idx, double val)
{
	if (*hold == NULL) {
		GMRFLib_idxval_create(hold);
	}
	if ((*hold)->n == (*hold)->n_alloc) {
		(*hold)->n_alloc += IMAX(IDX_ALLOC_ADD, (*hold)->n / 4);
		(*hold)->idx = Realloc((*hold)->idx, (*hold)->n_alloc, int);
		(*hold)->val = Realloc((*hold)->val, (*hold)->n_alloc, double);
	}
	(*hold)->idx[(*hold)->n] = idx;
	(*hold)->val[(*hold)->n] = val;
	(*hold)->n++;

	return GMRFLib_SUCCESS;
}

int GMRFLib_idxval_addto(GMRFLib_idxval_tp ** hold, int idx, double val)
{
	// if idx exists before, add val to value , otherwise just 'add'.
	// if there are two entries of 'idx', then only the first is used.

	if (*hold == NULL) {
		GMRFLib_idxval_create(hold);
	}

	int i;

	// hopefully this is it.
	i = (*hold)->iaddto;
	if ((*hold)->idx[i] == idx) {
		(*hold)->val[i] += val;
		return GMRFLib_SUCCESS;
	}
	// FIXME: this should be improved in general, but I think for the usage its ok. Since we are likely to add with same or increasing idx,
	// then I added this 'iaddto' which recall the last index, and try to be a little smarter.
	for (i = (*hold)->iaddto + 1; i < (*hold)->n; i++) {
		if ((*hold)->idx[i] == idx) {
			(*hold)->val[i] += val;
			(*hold)->iaddto = i;
			return GMRFLib_SUCCESS;
		}
	}
	for (i = 0; i < (*hold)->iaddto; i++) {
		if ((*hold)->idx[i] == idx) {
			(*hold)->val[i] += val;
			(*hold)->iaddto = i;
			return GMRFLib_SUCCESS;
		}
	}

	if ((*hold)->n == (*hold)->n_alloc) {
		(*hold)->n_alloc += IMAX(IDX_ALLOC_ADD, (*hold)->n / 4);
		(*hold)->idx = Realloc((*hold)->idx, (*hold)->n_alloc, int);
		(*hold)->val = Realloc((*hold)->val, (*hold)->n_alloc, double);
	}
	(*hold)->idx[(*hold)->n] = idx;
	(*hold)->val[(*hold)->n] = val;
	(*hold)->n++;

	return GMRFLib_SUCCESS;
}

int GMRFLib_str_is_member(GMRFLib_str_tp * hold, char *s, int case_sensitive, int *idx_match)
{
	if (hold == NULL) {
		return 0;
	}

	int (*cmp)(const char *, const char *) = (case_sensitive ? strcmp : strcasecmp);
	for (int i = 0; i < hold->n; i++) {
		if (cmp(s, hold->str[i]) == 0) {
			if (idx_match) {
				*idx_match = i;
			}
			return 1;
		}
	}
	return 0;
}





#undef MEMINFO
#undef IDX_ALLOC_ADD
#undef IDX_ALLOC_INITIAL
