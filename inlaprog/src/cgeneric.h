#ifndef __INLA_CGENERIC_H__
#       define __INLA_CGENERIC_H__

#       undef __BEGIN_DECLS
#       undef __END_DECLS
#       ifdef __cplusplus
#              define __BEGIN_DECLS extern "C" {
#              define __END_DECLS }
#       else
#              define __BEGIN_DECLS			       /* empty */
#              define __END_DECLS			       /* empty */
#       endif
__BEGIN_DECLS
#       include <assert.h>
#       include <math.h>
#       include <stdio.h>
#       include <omp.h>
#       include <stdio.h>

/* 
 *
 */
// same definition as in GMRFLibP.h...
#       if !defined(POSSIBLY_UNUSED)
#              ifdef __GNUC__
#                     define POSSIBLY_UNUSED(x) __attribute__((__unused__)) x
#              else
#                     define POSSIBLY_UNUSED(x) x
#              endif
#       endif
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

typedef enum {
	INLA_CLOGLIKE_INITIAL = 1,
	INLA_CLOGLIKE_LOG_PRIOR,
	INLA_CLOGLIKE_LOGLIKE,
	INLA_CLOGLIKE_CDF,
	INLA_CLOGLIKE_QUIT
} inla_cloglike_cmd_tp;

#       define INLA_CGENERIC_CMD_NAME(cmd_) ((cmd_) == INLA_CGENERIC_Q ? "Q" : \
				       ((cmd_) == INLA_CGENERIC_GRAPH ? "graph" : \
					((cmd_) == INLA_CGENERIC_MU ? "mu" : \
					 ((cmd_) == INLA_CGENERIC_INITIAL ? "initial" : \
					  ((cmd_) == INLA_CGENERIC_LOG_NORM_CONST ? "log_norm_const" : \
					   ((cmd_) == INLA_CGENERIC_LOG_PRIOR ? "log_prior" : \
					    ((cmd_) == INLA_CGENERIC_QUIT ? "quit" : "(***ERROR***)")))))))

#       define INLA_CLOGLIKE_CMD_NAME(cmd_) ((cmd_) == INLA_CLOGLIKE_INITIAL ? "initial" : \
				      ((cmd_) == INLA_CLOGLIKE_LOG_PRIOR ? "log_prior" : \
				       ((cmd_) == INLA_CLOGLIKE_LOGLIKE ? "log_like" : \
					((cmd_) == INLA_CLOGLIKE_LOGLIKE ? "CDF" : \
					 ((cmd_) == INLA_CLOGLIKE_QUIT ? "quit" : "(***ERROR***)")))))

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

	int processed;
	void *cache;
} inla_cgeneric_data_tp;

#       if defined(_OPENMP)
// tools useful for creating a cache
#              define IMAX_(a_,  b_) ((a_) >= (b_) ? (a_) : (b_))
#              define MAX_THREADS_(data_) IMAX_(1, ((data_)->threads.max))
#              define CGENERIC_CACHE_LEN(data_) (MAX_THREADS_(data_) * (MAX_THREADS_(data_) + 1))
#              define CGENERIC_CACHE_ASSIGN_IDX(idx_, data_)				\
        {								\
		int level1_ = omp_get_level();				\
		int tnum1_ = omp_get_thread_num();			\
		if (level1_ <= 1) {					\
			idx_ =  tnum1_;					\
		} else if (level1_ == 2) {				\
			int mt_ = MAX_THREADS_(data_);			\
			int tnum2_ = omp_get_ancestor_thread_num(level1_ -1); \
			idx_ = mt_ + tnum1_ + tnum2_ * mt_;		\
		} else {						\
			assert(0 == 1);					\
		}							\
	}
#       endif

typedef double *inla_cgeneric_func_tp(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp * data);
typedef double *inla_cloglike_func_tp(inla_cloglike_cmd_tp cmd, double *theta, inla_cgeneric_data_tp * data,
				      int ny, double *y, int nx, double *x, double *result);

inla_cgeneric_data_tp *inla_cgeneric_read_data(const char *filename, int debug);
inla_cgeneric_func_tp inla_cgeneric_iid_model;
inla_cgeneric_func_tp inla_cgeneric_ar1_model;
inla_cgeneric_func_tp inla_cgeneric_generic0_model;

static void POSSIBLY_UNUSED(inla_cgeneric_data_print) (FILE * fp, inla_cgeneric_data_tp * data) {
	fp = (fp ? fp : stdout);

	fprintf(fp, "\nContent of cgeneric_data\n");
	fprintf(fp, "\tthread.max   = [%1d]\n", data->threads.max);
	fprintf(fp, "\tthread.outer = [%1d]\n", data->threads.outer);
	fprintf(fp, "\tthread.inner = [%1d]\n", data->threads.inner);

	fprintf(fp, "\tnumber of ints = [%1d]\n", data->n_ints);
	for (int i = 0; i < data->n_ints; i++) {
		fprintf(fp, "\t\tname[%s] length[%1d]\n", data->ints[i]->name, data->ints[i]->len);
	}
	fprintf(fp, "\tnumber of doubles = [%1d]\n", data->n_doubles);
	for (int i = 0; i < data->n_doubles; i++) {
		fprintf(fp, "\t\tname[%s] length[%1d]\n", data->doubles[i]->name, data->doubles[i]->len);
	}
	fprintf(fp, "\tnumber of chars = [%1d]\n", data->n_chars);
	for (int i = 0; i < data->n_chars; i++) {
		fprintf(fp, "\t\tname[%s] length[%1d] value[%s]\n", data->chars[i]->name, data->chars[i]->len, data->chars[i]->chars);
	}
	fprintf(fp, "\tnumber of matrices = [%1d]\n", data->n_mats);
	for (int i = 0; i < data->n_mats; i++) {
		fprintf(fp, "\t\tname[%s] dimension[%1d x %1d]\n", data->mats[i]->name, data->mats[i]->nrow, data->mats[i]->ncol);
	}
	fprintf(fp, "\tnumber of sparse matrices = [%1d]\n", data->n_smats);
	for (int i = 0; i < data->n_smats; i++) {
		fprintf(fp, "\t\tname[%s] dimension[%1d x %1d] nelements[%1d]\n", data->smats[i]->name, data->smats[i]->nrow, data->smats[i]->ncol,
			data->smats[i]->n);
	}
	fprintf(fp, "\n");
}

__END_DECLS
#endif
