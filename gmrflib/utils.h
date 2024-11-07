
/* utils.h
 * 
 * Copyright (C) 2006-2024 Havard Rue
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
  \file utils.h
  \brief Typedefs for \ref utils.c
*/

#ifndef __GMRFLib_UTILS_H__
#define __GMRFLib_UTILS_H__

#include <stdlib.h>
#include <stddef.h>
#include <math.h>

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
#include "GMRFLib/hashP.h"
#include "GMRFLib/graph.h"
#include "GMRFLib/GMRFLibP.h"
#define GMRFLib_printMem(fp_) GMRFLib_printMem_core(fp_, __FILE__, __LINE__)
    typedef struct {
	int nrow;
	map_ivp *vmat;
} GMRFLib_vmatrix_tp;


char *Strdup(const char *s);

char *GMRFLib_memcheck_make_tag(size_t size, const char *file, const char *funcname, int lineno);
char *GMRFLib_rindex(const char *p, int ch);
char *GMRFLib_strtok_r(char *s1, const char *s2, char **lasts);
char *GMRFLib_vec2char(double *arr, int len);
const char *GMRFLib_function_name_strip(const char *name);
double *GMRFLib_vmatrix_get(GMRFLib_vmatrix_tp * vmatrix, int i, int j);
double GMRFLib_cdfnorm(double x);
double GMRFLib_cdfnorm_inv(double p);
double GMRFLib_eps(double power);
double GMRFLib_erf(double x);
double GMRFLib_erf_inv(double x);
double GMRFLib_erfc(double x);
double GMRFLib_erfc_inv(double x);
double GMRFLib_inv_logit(double x);
double GMRFLib_log_apbex(double a, double b);
double GMRFLib_logit(double p);
double GMRFLib_max_value(double *x, int n, int *idx);
double GMRFLib_min_value(double *x, int n, int *idx);
double GMRFLib_signed_pow(double x, double power);
size_t GMRFLib_align(size_t n, size_t size);
size_t GMRFLib_align_simple(size_t n, size_t size);
int GMRFLib_2order_poleq(double *sol1, double *sol2, double a, double b, double c);
int GMRFLib_adjust_vector(double *x, int n);
int GMRFLib_dcmp(const void *a, const void *b);
int GMRFLib_dcmp_abs(const void *a, const void *b);
int GMRFLib_dcmp_abs_r(const void *a, const void *b);
int GMRFLib_dcmp_r(const void *a, const void *b);
int GMRFLib_debug_functions(const char *name);
int GMRFLib_find_nonzero(double *array, int len, int direction);
int GMRFLib_find_value(double *array, int len, int direction, double value);
int GMRFLib_gsl_mat2plain(double **out, gsl_matrix * mat);
int GMRFLib_gsl_matrix_count_eq(gsl_matrix * A, double value);
int GMRFLib_gsl_vec2plain(double **out, gsl_vector * vec);
int GMRFLib_is_zero(double *x, int n);
int GMRFLib_icmp(const void *a, const void *b);
int GMRFLib_icmp_r(const void *a, const void *b);
int GMRFLib_imax_value(int *x, int n, int *idx);
int GMRFLib_imin_value(int *x, int n, int *idx);
int GMRFLib_iamax_value(int *x, int n, int *idx);
int GMRFLib_is_int(char *str, int *value);
int GMRFLib_iuniques(int *nuniques, int **uniques, int *ix, int nx);
int GMRFLib_iwhich_sorted(int val, int *ix, int len);
int GMRFLib_iwhich_sorted_g(int val, int *ix, int len, int *guess_guess);
int GMRFLib_iwhich_sorted_g2(int val, int *ix, int len, int *guess);
int GMRFLib_iwhich_sorted_g_new(int key, int *ix, int len, int *low_guess);
int GMRFLib_memcheck_error(const char *msg, void *p, const char *file, const char *funcname, int lineno);
int GMRFLib_memcheck_printf(FILE * fp);
int GMRFLib_memcheck_register(void *p, size_t size, const char *file, const char *funcname, int lineno);
int GMRFLib_memcheck_remove(void *p, const char *file, const char *funcname, int lineno);
int GMRFLib_normalize(int n, double *x);
int GMRFLib_print_darray(FILE * fp, double *x, int n, const char *desc);
int GMRFLib_print_iarray(FILE * fp, int *x, int n, const char *desc);
int GMRFLib_printf_gsl_matrix(FILE * fp, gsl_matrix * matrix, const char *format);
int GMRFLib_printf_gsl_matrix2(FILE * fp, gsl_matrix * matrix, const char *format, double cutoff);
int GMRFLib_printf_gsl_vector(FILE * fp, gsl_vector * vector, const char *format);
int GMRFLib_printf_matrix(FILE * fp, double *A, int m, int n);
int GMRFLib_scale_vector(double *x, int n);
int GMRFLib_sprintf(char **ptr, const char *fmt, ...);
int GMRFLib_trace_functions(const char *name);
int GMRFLib_unique_additive(int *n, double *x, double eps);
int GMRFLib_unique_additive2(int *n, double *x, double *y, double eps);
int GMRFLib_unique_relative(int *n, double *x, double eps);
int GMRFLib_unique_relative2(int *n, double *x, double *y, double eps);
int GMRFLib_vmatrix_free(GMRFLib_vmatrix_tp * vmatrix, int free_content);
int GMRFLib_vmatrix_init(GMRFLib_vmatrix_tp ** vmatrix, int nrow, GMRFLib_graph_tp * graph);
int GMRFLib_vmatrix_set(GMRFLib_vmatrix_tp * vmatrix, int i, int j, double *vec);
int GMRFLib_which(double val, double *array, int len);
int my_sort2_dd_test_cutoff(int verbose);
int my_sort2_id_test_cutoff(int verbose);
map_id *GMRFLib_duplicate_map_id(map_id * hash);
map_ii *GMRFLib_duplicate_map_ii(map_ii * hash);
mapkit_size_t GMRFLib_nelm_map_id(map_id * hash);
mapkit_size_t GMRFLib_nelm_map_ii(map_ii * hash);
void *GMRFLib_calloc(size_t nmemb, size_t size, const char *file, const char *funcname, int lineno);
void *GMRFLib_malloc(size_t size, const char *file, const char *funcname, int lineno);
void *GMRFLib_memcpy(void *dest, const void *src, size_t n);
void *GMRFLib_realloc(void *old_ptr, size_t size, const char *file, const char *funcname, int lineno);
void GMRFLib_delay(int msec);
void GMRFLib_delay_random(int msec_low, int msec_high);
void GMRFLib_free(void *ptr, const char *file, const char *funcname, int lineno);
void GMRFLib_getMemory(int *currRealMem, int *peakRealMem, int *currVirtMem, int *peakVirtMem);
void GMRFLib_printMem_core(FILE * fp, const char *file, int lineno);
void GMRFLib_malloc_debug_check(void);
void gsl_sort2_dd(double *data1, double *data2, const int n);
void gsl_sort2_id(int *data1, double *data2, const int n);
void gsl_sort2_ii(int *data1, int *data2, const int n);
void my_insertionSort_dd(double *iarr, double *darr, int n);
void my_insertionSort_id(int *iarr, double *darr, int n);
void my_insertionSort_ii(int *iarr, int *darr, int n);
void my_sort2_dd(double *ix, double *x, int n);
void my_sort2_id(int *ix, double *x, int n);
void my_sort2_id_x(int *ix, double *x, int n, void *work);
void my_sort2_ii(int *ix, int *x, int n);


int GMRFLib_is_sorted(void *a, size_t n, size_t size, int (*cmp)(const void *, const void *));
int GMRFLib_is_sorted_ddec(int n, double *a);
int GMRFLib_is_sorted_ddec_plain(int n, double *a);
int GMRFLib_is_sorted_dinc(int n, double *a);
int GMRFLib_is_sorted_dinc_plain(int n, double *a);
int GMRFLib_is_sorted_idec(int n, int *a);
int GMRFLib_is_sorted_idec_plain(int n, int *a);
int GMRFLib_is_sorted_iinc(int n, int *a);
int GMRFLib_is_sorted_iinc_plain(int n, int *a);
void GMRFLib_qsort(void *a, size_t n, size_t size, int (*cmp)(const void *, const void *));
void GMRFLib_qsort2(void *x, size_t nmemb, size_t size_x, void *y, size_t size_y, int (*compar)(const void *, const void *));

int GMRFLib_get_cachelinesize(void);


__END_DECLS
#endif
