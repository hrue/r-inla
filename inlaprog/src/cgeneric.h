
/* cgeneric.h
 * 
 * Copyright (C) 2021-2022 Havard Rue
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
 *      matrix storage, stored column by column, like
 *      > matrix(1:6,2,3)
 *      [,1] [,2] [,3]
 *      [1,]    1    3    5
 *      [2,]    2    4    6
 *      > c(matrix(1:6,2,3))
 *      [1] 1 2 3 4 5 6
 */
typedef struct {
	char *name;
	int nrow;
	int ncol;
	double *x;
} inla_cgeneric_mat_tp;

/* 
 * sparse matrix format, stored used 0-based indices, like
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
	int n_ints;
	inla_cgeneric_vec_tp **ints;

	int n_doubles;
	inla_cgeneric_vec_tp **doubles;

	int n_chars;
	inla_cgeneric_vec_tp **chars;

	int n_mat;
	inla_cgeneric_mat_tp **mats;

	int n_smat;
	inla_cgeneric_smat_tp **smats;
} inla_cgeneric_data_tp;

typedef double *inla_cgeneric_func_tp(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp * data);

inla_cgeneric_data_tp *inla_cgeneric_read_data(const char *filename, int debug);
inla_cgeneric_func_tp inla_cgeneric_iid_model;
inla_cgeneric_func_tp inla_cgeneric_ar1_model;

__END_DECLS
#endif
