
/* utils.h
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

/*!
  \file utils.h
  \brief Typedefs for \ref utils.c
*/

#ifndef __GMRFLib_IDXVAL_H__
#define __GMRFLib_IDXVAL_H__

#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif

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
#include "GMRFLib/GMRFLibP.h"
    typedef struct {
	int n;
	int n_alloc;
	int *idx;
} GMRFLib_idx_tp;

typedef struct {
	int n;
	int n_alloc;
	char **str;
} GMRFLib_str_tp;

typedef struct {
	int n;
	int n_alloc;
	int **idx;
} GMRFLib_idx2_tp;

typedef struct {
	int n;
	int n_alloc;
	double *val;
} GMRFLib_val_tp;

typedef enum {
	IDXVAL_UNKNOWN = 0,				       /* do not change */
	IDXVAL_SERIAL,
	IDXVAL_GROUP
} GMRFLib_idxval_preference_tp;

typedef struct {
	int n;
	int n_alloc;
	int iaddto;
	int *idx;
	GMRFLib_idxval_preference_tp preference;

	int g_n;					       /* number of groups with sequential indices */
	int *g_len;					       /* their length */
	int *g_i;					       /* and their starting index */
	int *g_1;					       /* indicator if this group have 'val' all equal to 1.0 */

	double *val;
} GMRFLib_idxval_tp;

#define DDOT(N_, X_, Y_) ddot_(&(N_), X_, &integer_one, Y_, &integer_one)
#define DSUM(N_, X_) my_dsum(N_, X_)
#define DOT_PRODUCT_GROUP(VALUE_, ELM_, ARR_)				\
	if (1) {							\
		double value_ = 0.0;					\
		int integer_one = 1;					\
		for (int g_ = 0; g_ < ELM_->g_n; g_++) {		\
			int istart_ = ELM_->g_i[g_];			\
			int *ii_ = &(ELM_->idx[istart_]);		\
			int len_ = ELM_->g_len[g_];			\
			double *vv_ = &(ELM_->val[istart_]);		\
									\
			if (len_ == 0) continue;			\
			if (len_ > 0) {					\
				double *aa_ = &(ARR_[0]);		\
				if (ELM_->g_1[g_]) {			\
					if (len_ < 8L) {		\
						_Pragma("GCC ivdep")	\
							for (int i_ = 0; i_ < len_; i_++) { \
								value_ += aa_[ii_[i_]]; \
							}		\
					} else {			\
						value_ += my_dsum_idx(len_, aa_, ii_); \
					}				\
				} else {				\
					if (len_ < 8L) {		\
						_Pragma("GCC ivdep")	\
							for (int i_ = 0; i_ < len_; i_++) { \
								value_ += vv_[i_] * aa_[ii_[i_]]; \
							}		\
					} else {			\
						value_ += my_ddot_idx(len_, vv_, aa_, ii_); \
					}				\
				}					\
			} else if (len_ < 0) {				\
				int llen_ = - len_;			\
				double *aa_ = &(ARR_[ii_[0]]);		\
				if (ELM_->g_1[g_]) {			\
					value_ += DSUM(llen_, aa_);	\
				} else {				\
					value_ += DDOT(llen_, vv_, aa_); \
				}					\
			}						\
			if (g_ < ELM_->g_n - 1) __builtin_prefetch(&(ARR_[ELM_->idx[ELM_->g_i[g_ + 1]]])); \
		}							\
		VALUE_ = (typeof(VALUE_)) value_;			\
	}
#define DOT_PRODUCT_SERIAL(VALUE_, ELM_, ARR_)				\
	if (1) {							\
		double value_ = 0.0;					\
		double *vv_ = ELM_->val;				\
		double *aa_ = ARR_;					\
		int *idx_ = ELM_->idx;					\
		if (ELM_->n < 8L) {					\
			_Pragma("GCC ivdep")				\
				for (int i_ = 0; i_ < ELM_->n; i_++) {	\
					value_ += vv_[i_] * aa_[idx_[i_]]; \
				}					\
		} else {						\
			value_ += my_ddot_idx(ELM_->n, vv_, aa_, idx_);	\
		}							\
		VALUE_ = (typeof(VALUE_)) value_;			\
	}
#define DOT_PRODUCT(VALUE_, ELM_, ARR_)					\
	if (1) {							\
		switch(ELM_->preference) {				\
		case IDXVAL_SERIAL:					\
			DOT_PRODUCT_SERIAL(VALUE_, ELM_, ARR_);		\
			break;						\
		case IDXVAL_GROUP:					\
			DOT_PRODUCT_GROUP(VALUE_, ELM_, ARR_);		\
			break;						\
		case IDXVAL_UNKNOWN:					\
		default:						\
			DOT_PRODUCT_GROUP(VALUE_, ELM_, ARR_);		\
			break;						\
		}							\
	}



GMRFLib_idx2_tp **GMRFLib_idx2_ncreate(int n);
GMRFLib_idx2_tp **GMRFLib_idx2_ncreate_x(int n, int len);
GMRFLib_idx_tp **GMRFLib_idx_ncreate(int n);
GMRFLib_idx_tp **GMRFLib_idx_ncreate(int n);
GMRFLib_idx_tp **GMRFLib_idx_ncreate_x(int n, int len);
GMRFLib_idxval_tp **GMRFLib_idxval_ncreate(int n);
GMRFLib_idxval_tp **GMRFLib_idxval_ncreate_x(int n, int len);
GMRFLib_val_tp **GMRFLib_val_ncreate(int n);
int GMRFLib_idx2_add(GMRFLib_idx2_tp ** hold, int idx0, int idx1);
int GMRFLib_idx2_create(GMRFLib_idx2_tp ** hold);
int GMRFLib_idx2_create_x(GMRFLib_idx2_tp ** hold, int len);
int GMRFLib_idx2_free(GMRFLib_idx2_tp * hold);
int GMRFLib_idx2_nprune(GMRFLib_idx2_tp ** a, int n);
int GMRFLib_idx2_printf(FILE * fp, GMRFLib_idx2_tp * hold, char *msg);
int GMRFLib_idx2_prune(GMRFLib_idx2_tp * hold);
int GMRFLib_idx_add(GMRFLib_idx_tp ** hold, int idx);
int GMRFLib_idx_create(GMRFLib_idx_tp ** hold);
int GMRFLib_idx_create_x(GMRFLib_idx_tp ** hold, int len);
int GMRFLib_idx_free(GMRFLib_idx_tp * hold);
int GMRFLib_idx_nprune(GMRFLib_idx_tp ** a, int n);
int GMRFLib_idx_nsort(GMRFLib_idx_tp ** a, int n, int nt);
int GMRFLib_idx_nuniq(GMRFLib_idx_tp ** a, int n, int nt);
int GMRFLib_idx_printf(FILE * fp, GMRFLib_idx_tp * hold, char *msg);
int GMRFLib_idx_prune(GMRFLib_idx_tp * hold);
int GMRFLib_idx_sort(GMRFLib_idx_tp * hold);
int GMRFLib_idx_uniq(GMRFLib_idx_tp * hold);
int GMRFLib_idxval_add(GMRFLib_idxval_tp ** hold, int idx, double val);
int GMRFLib_idxval_addto(GMRFLib_idxval_tp ** hold, int idx, double val);
int GMRFLib_idxval_create(GMRFLib_idxval_tp ** hold);
int GMRFLib_idxval_create_x(GMRFLib_idxval_tp ** hold, int len);
int GMRFLib_idxval_create_x(GMRFLib_idxval_tp ** hold, int len);
int GMRFLib_idxval_free(GMRFLib_idxval_tp * hold);
int GMRFLib_idxval_info_printf(FILE * fp, GMRFLib_idxval_tp * hold, char *msg);
int GMRFLib_idxval_nprune(GMRFLib_idxval_tp ** a, int n, int nt);
int GMRFLib_idxval_nsort(GMRFLib_idxval_tp ** hold, int n, int nt);
int GMRFLib_idxval_nsort_x(GMRFLib_idxval_tp ** hold, int n, int nt, int prune_zeros);
int GMRFLib_idxval_nuniq(GMRFLib_idxval_tp ** a, int n, int nt);
int GMRFLib_idxval_printf(FILE * fp, GMRFLib_idxval_tp * hold, char *msg);
int GMRFLib_idxval_prune(GMRFLib_idxval_tp * hold);
int GMRFLib_idxval_sort(GMRFLib_idxval_tp * hold);
int GMRFLib_idxval_uniq(GMRFLib_idxval_tp * hold);
int GMRFLib_str_add(GMRFLib_str_tp ** hold, char *s);
int GMRFLib_str_create_x(GMRFLib_str_tp ** hold, int len);
int GMRFLib_str_is_member(GMRFLib_str_tp * hold, char *s, int case_sensitive, int *idx_match);
int GMRFLib_str_prune(GMRFLib_str_tp * hold);
int GMRFLib_val_add(GMRFLib_val_tp ** hold, double val);
int GMRFLib_val_create(GMRFLib_val_tp ** hold);
int GMRFLib_val_free(GMRFLib_val_tp * hold);
int GMRFLib_val_nprune(GMRFLib_val_tp ** a, int n);
int GMRFLib_val_printf(FILE * fp, GMRFLib_val_tp * hold, char *msg);
int GMRFLib_val_prune(GMRFLib_val_tp * hold);


__END_DECLS
#endif