
/* cgeneric.h
 * 
 * Copyright (C) 2021-2023 Havard Rue
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
#ifndef __INLA_CGENERIC_H__
#define __INLA_CGENERIC_H__

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
 *
 */
    typedef enum {
	INLA_CGENERIC_VOID = 0,
	INLA_CGENERIC_Q,
	INLA_CGENERIC_GRAPH,
	INLA_CGENERIC_MU,
	INLA_CGENERIC_INITIAL,
	INLA_CGENERIC_LOG_NORM_CONST,
	INLA_CGENERIC_LOG_PRIOR,
	INLA_CGENERIC_QUIT
} inla_cgeneric_cmd_tp;

#define INLA_CGENERIC_CMD_NAME(cmd_) ((cmd_) == INLA_CGENERIC_VOID ? "void" : \
				      ((cmd_) == INLA_CGENERIC_Q ? "Q" : \
				       ((cmd_) == INLA_CGENERIC_GRAPH ? "graph" : \
					((cmd_) == INLA_CGENERIC_MU ? "mu" : \
					 ((cmd_) == INLA_CGENERIC_INITIAL ? "initial" : \
					  ((cmd_) == INLA_CGENERIC_LOG_NORM_CONST ? "log_norm_const" : \
					   ((cmd_) == INLA_CGENERIC_LOG_PRIOR ? "log_prior" : \
					    ((cmd_) == INLA_CGENERIC_QUIT ? "quit" : "(***ERROR***)"))))))))

/*
 *      matrix storage is stored row by row.
 *
 *      In R the (default) storage is column by column, like
 *      > matrix(1:6,2,3)
 *      [,1] [,2] [,3]
 *      [1,]    1    3    5
 *      [2,]    2    4    6
 *
 *      hence
 *      > c(matrix(1:6,2,3))
 *      [1] 1 2 3 4 5 6
 *
 *      while in cgeneric, the matrix elements 'x', is stored as
 *      x[] = { 1, 3, 5, 2, 4, 6 }
 */
typedef struct {
	char *name;
	int nrow;
	int ncol;
	double *x;
} inla_cgeneric_mat_tp;

/* 
 * sparse matrix format, stored used 0-based indices.
 * If the matrix is not a 'dgTMatrix', it is converted using
 * 'inla.as.sparse()'
 *
 * the matrix is stored in the order it appears, like
 *
 *      > A <- inla.as.sparse(matrix(c(1,2,3,0,0,6),2,3))
 *      > A
 *      2 x 3 sparse Matrix of class "dgTMatrix"
 *      [1,] 1 3 .
 *      [2,] 2 . 6
 *      > cbind(i=A@i, j=A@j, x=A@x)
 *            i j x   
 *      [1,]  0 0 1
 *      [2,]  1 0 2
 *      [3,]  0 1 3
 *      [4,]  1 2 6
 *
 *      If you want to have this matrix stored this column by column, 
 *      then swap @i and @j, and rearrange @x, before passing them into
 *      'inla.cgeneric.define'.
 *
 *      For example, if we want to pass only the upper half of a
 *      symmetric sparse matrix, stored by column, then we can do
 *
 *      A <- inla.as.sparse(A)
 *      ii <- A@i
 *      A@i <- A@j
 *      A@j <- ii
 *      idx <- which(A@i <= A@j) 
 *      A@i <- A@i[idx]
 *      A@j <- A@j[idx]
 *      A@x <- A@x[idx]
 *
 *      before passing it on into 'inla.cgeneric.define'
 */
typedef struct {
	char *name;
	int nrow;
	int ncol;
	int n;						       /* number of triplets (i,j,x) */
	int *i;
	int *j;
	double *x;
} inla_cgeneric_smat_tp;

typedef struct {
	char *name;
	int len;
	int *ints;
	double *doubles;
	char *chars;
} inla_cgeneric_vec_tp;

typedef struct {
	// inla .. -t A:B
	// max = maximum number of threads = A*B
	// outer = number of threads in the outer loop (A)
	// inner = number of threads in the inner loop (B)
	int max;
	int outer;
	int inner;
} inla_cgeneric_threads_tp;

typedef struct {
	inla_cgeneric_threads_tp threads;

	int n_ints;
	inla_cgeneric_vec_tp **ints;

	int n_doubles;
	inla_cgeneric_vec_tp **doubles;

	int n_chars;
	inla_cgeneric_vec_tp **chars;

	int n_mats;
	inla_cgeneric_mat_tp **mats;

	int n_smats;
	inla_cgeneric_smat_tp **smats;

	void *cache;
} inla_cgeneric_data_tp;

#if defined(_OPENMP)
// tools useful for creating a cache
#include <omp.h>
#define IMAX_(a_,  b_) ((a_) >= (b_) ? (a_) : (b_))
#define MAX_THREADS(data_) ((data_)->max_threads)
#define CGENERIC_CACHE_LEN(data_) (IMAX_(1, MAX_THREADS(data_)) * (IMAX_(1, MAX_THREADS(data_)) + 1))
#define CGENERIC_CACHE_ASSIGN_IDX(idx_, data_)				\
        if (1) {                                                        \
                int level_ = omp_get_level();                           \
                int tnum_ = omp_get_thread_num();                       \
                if (level_ <= 1)        {                               \
                        idx_ =  tnum_;                                  \
                } else if (level_ == 2) {                               \
                        int level2_ = omp_get_ancestor_thread_num(level_ -1); \
			idx_ = IMAX_(1, 1 + level2_) * IMAX_(1, MAX_THREADS(data_)) + tnum_; \
                } else {                                                \
                        assert(0 == 1);                                 \
                }                                                       \
        }
#endif

typedef double *inla_cgeneric_func_tp(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp * data);

inla_cgeneric_data_tp *inla_cgeneric_read_data(const char *filename, int debug);
inla_cgeneric_func_tp inla_cgeneric_iid_model;
inla_cgeneric_func_tp inla_cgeneric_ar1_model;
inla_cgeneric_func_tp inla_cgeneric_generic0_model;

__END_DECLS
#endif
