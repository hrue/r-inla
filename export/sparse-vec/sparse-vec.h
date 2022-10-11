
/* idxval.h
 * 
 * Copyright (C) 2022-2022 Havard Rue
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

#ifndef __DOT_PRODUCT_H__
#define __DOT_PRODUCT_H__

#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <malloc.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS					       /* empty */
#define __END_DECLS					       /* empty */
#endif

__BEGIN_DECLS typedef enum {
	IDXVAL_UNKNOWN = 0,				       /* do not change */
	IDXVAL_SERIAL,
	IDXVAL_SERIAL_MKL,
	IDXVAL_GROUP,
	IDXVAL_GROUP_MKL
} GMRFLib_idxval_preference_tp;

typedef struct {
	int n;
	int n_alloc;
	int iaddto;
	int *idx;
	GMRFLib_idxval_preference_tp preference;

	int n_n;					       /* new len */
	int *g_idx;
	int g_n;					       /* number of groups with sequential indices */
	int *g_len;					       /* their length */
	int *g_i;					       /* and their starting index */
	int *g_1;					       /* indicator if this group have 'val' all equal to 1.0 */
	int free_g_mem;

	double *val;
	double *g_val;
} GMRFLib_idxval_tp;

GMRFLib_idxval_tp **GMRFLib_idxval_ncreate(int n);
GMRFLib_idxval_tp **GMRFLib_idxval_ncreate_x(int n, int len);

int GMRFLib_idxval_add(GMRFLib_idxval_tp ** hold, int idx, double val);
int GMRFLib_idxval_addto(GMRFLib_idxval_tp ** hold, int idx, double val);
int GMRFLib_idxval_create(GMRFLib_idxval_tp ** hold);
int GMRFLib_idxval_create_x(GMRFLib_idxval_tp ** hold, int len);
int GMRFLib_idxval_create_x(GMRFLib_idxval_tp ** hold, int len);
int GMRFLib_idxval_free(GMRFLib_idxval_tp * hold);
int GMRFLib_idxval_info_printf(FILE * fp, GMRFLib_idxval_tp * hold, const char *msg);
int GMRFLib_idxval_nprune(GMRFLib_idxval_tp ** a, int n, int nt);
int GMRFLib_idxval_nsort(GMRFLib_idxval_tp ** hold, int n, int nt);
int GMRFLib_idxval_nsort_x(GMRFLib_idxval_tp ** hold, int n, int nt, int prune_zeros);
int GMRFLib_idxval_nuniq(GMRFLib_idxval_tp ** a, int n, int nt);
int GMRFLib_idxval_printf(FILE * fp, GMRFLib_idxval_tp * hold, const char *msg);
int GMRFLib_idxval_prune(GMRFLib_idxval_tp * hold);
int GMRFLib_idxval_sort(GMRFLib_idxval_tp * hold);
int GMRFLib_idxval_uniq(GMRFLib_idxval_tp * hold);
int GMRFLib_idxval_prepare(GMRFLib_idxval_tp * hold);

double GMRFLib_ddot_idx_mkl(int n, double *__restrict v, double *__restrict a, int *__restrict idx);
double GMRFLib_ddot(int n, double *__restrict x, double *__restrict y);
double GMRFLib_ddot_idx(int n, double *__restrict v, double *__restrict a, int *__restrict idx);
double GMRFLib_dsum(int n, double *x);
double GMRFLib_dsum_idx(int n, double *__restrict a, int *__restrict idx);

double GMRFLib_dot_product_group(GMRFLib_idxval_tp * __restrict ELM_, double *__restrict ARR_);
double GMRFLib_dot_product_group_mkl(GMRFLib_idxval_tp * __restrict ELM_, double *__restrict ARR_);
double GMRFLib_dot_product_serial(GMRFLib_idxval_tp * __restrict ELM_, double *__restrict ARR_);
double GMRFLib_dot_product_serial_mkl(GMRFLib_idxval_tp * __restrict ELM_, double *__restrict ARR_);
double GMRFLib_dot_product(GMRFLib_idxval_tp * __restrict ELM_, double *__restrict ARR_);

double GMRFLib_min_value(double *x, int n, int *idx);
double GMRFLib_max_value(double *x, int n, int *idx);

// these are better names to export
#define GMRFLib_sparse_vec_tp GMRFLib_idxval_tp
#define GMRFLib_sparse_vec_add(a_,  b_, c_) GMRFLib_idxval_add(a_, b_, c_)
#define GMRFLib_sparse_vec_prepare(a_) GMRFLib_idxval_prepare(a_)
#define GMRFLib_sparse_vec_dot_product(a_, b_) GMRFLib_dot_product(a_, b_)
#define GMRFLib_sparse_vec_free(a_) GMRFLib_idxval_free(a_)

__END_DECLS
#endif
