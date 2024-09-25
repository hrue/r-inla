
/* lapack-interfac.h
 * 
 * Copyright (C) 2001-2024 Havard Rue
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
  \file lapack-interface.h
  \brief Typedefs and defines for \ref lapack-interface.c
*/

#ifndef __GMRFLib_LAPACK_INTERFACE_H__
#define __GMRFLib_LAPACK_INTERFACE_H__

#include <stdlib.h>

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
//
    typedef struct {
	gsl_matrix *U;
	gsl_matrix *D;
	gsl_vector *S;
	gsl_eigen_symmv_workspace *work;
} GMRFLib_gsl_low_rank_store_tp;

typedef struct {
	gsl_matrix *U;
	gsl_matrix *M1;
	gsl_matrix *M2;
	gsl_matrix *AA;
	gsl_vector *S;
	gsl_eigen_symmv_workspace *work;
} GMRFLib_gsl_ensure_spd_store_tp;

typedef struct {
	gsl_matrix *L;
	gsl_vector *S;
} GMRFLib_gsl_spd_solve_store_tp;

typedef struct {
	gsl_matrix *L;
	gsl_vector *xx;
} GMRFLib_gsl_ldnorm_store_tp;

#define BLAS_LEVEL2 2
#define BLAS_LEVEL3 3
int dpbtrf_(const char *, int *, int *, double *, int *, int *, fortran_charlen_t);
int dpbtf2_(const char *, int *, int *, double *, int *, int *, fortran_charlen_t);
int dpotrf_(const char *, int *, double *, int *, int *, fortran_charlen_t);
int dpotf2_(const char *, int *, double *, int *, int *, fortran_charlen_t);
int dtbsv_(const char *, const char *, const char *, int *, int *, double *, int *, double *, int *,
	   fortran_charlen_t, fortran_charlen_t, fortran_charlen_t);
int dpotri_(const char *, int *, double *, int *, int *, fortran_charlen_t);
int dgemm_(const char *, const char *, int *, int *, int *, double *, double *, int *, double *, int
	   *, double *, double *, int *, fortran_charlen_t, fortran_charlen_t);
int dgemv_(const char *, int *, int *, double *, double *, int *, double *, int *, double *, double
	   *, int *, fortran_charlen_t);
int dgemv_failsafe_(const char *, int *, int *, double *, double *, int *, double *, int *, double *, double
		    *, int *, fortran_charlen_t);
double erfi_(double *);
int dpotrs_(const char *, int *, int *, double *, int *, double *, int *, int *, fortran_charlen_t);
int dchdc_(double *, int *, int *, double *, int *, int *, int *, double *);
int dtrmv_(const char *, const char *, const char *, int *, double *, int *, double *, int *,
	   fortran_charlen_t, fortran_charlen_t, fortran_charlen_t);
int idamax_(int *, double *, int *);
void daxpyi_(int *, double *, double *, int *, double *);

int GMRFLib_comp_chol_general(double **chol, double *matrix, int dim, double *logdet, int ecode);
int GMRFLib_comp_chol_semidef(double **chol, int **map, int *rank, double *matrix, int dim, double *logdet, double eps);
int GMRFLib_comp_posdef_inverse(double *matrix, int dim);



GMRFLib_gsl_ensure_spd_store_tp *GMRFLib_gsl_ensure_spd_store_alloc(int n);
int GMRFLib_gsl_ensure_spd_store_free(GMRFLib_gsl_ensure_spd_store_tp * S);
int GMRFLib_gsl_ensure_spd(gsl_matrix * A, double tol, char **msg);
int GMRFLib_gsl_ensure_spd_x(gsl_matrix * A, double tol, char **mgs, GMRFLib_gsl_ensure_spd_store_tp * store);
int GMRFLib_gsl_ensure_spd_inverse(gsl_matrix * A, double tol, char **msg);
int GMRFLib_gsl_ensure_spd_inverse_x(gsl_matrix * A, double tol, char **msg, GMRFLib_gsl_ensure_spd_store_tp * store);
int GMRFLib_gsl_ensure_spd_core(gsl_matrix * A, double tol, int method, char **msg, GMRFLib_gsl_ensure_spd_store_tp * store);
int GMRFLib_ensure_spd(double *A, int dim, double tol, char **msg);
int GMRFLib_ensure_spd_x(double *A, int dim, double tol, char **msg, GMRFLib_gsl_ensure_spd_store_tp * store);
int GMRFLib_solveAxb_posdef(double *sol, double *chol, double *b, int dim, int nrhs);

int daxpy_(int *n, double *alpha, double *x, int *incx, double *y, int *incy);
int dcopy_(int *n, double *x, int *incx, double *y, int *incy);
int dscal_(int *n, double *alpha, double *x, int *inc);
int dscal_(int *n, double *a, double *x, int *incx);

double ddot_(int *len, double *x, int *incx, double *y, int *incy);
double dnrm2_(int *n, double *x, int *inc);
double dasum_(int *n, double *x, int *inc);

int GMRFLib_dscale(int n, double a, double *x);
double GMRFLib_dssqr(int n, double *x);

GMRFLib_gsl_ldnorm_store_tp *GMRFLib_gsl_ldnorm_store_alloc(int n);
int GMRFLib_gsl_ldnorm_store_free(GMRFLib_gsl_ldnorm_store_tp * store);
double GMRFLib_gsl_ldnorm(gsl_vector * x, gsl_vector * mean, gsl_matrix * Q, gsl_matrix * S, int identity);
double GMRFLib_gsl_ldnorm_x(gsl_vector * x, gsl_vector * mean, gsl_matrix * Q, gsl_matrix * S, int identity, GMRFLib_gsl_ldnorm_store_tp * store);

double GMRFLib_gsl_spd_logdet(gsl_matrix * A);
double GMRFLib_gsl_xQx(gsl_vector * x, gsl_matrix * Q);
double GMRFLib_gsl_kld(gsl_vector * m_base, gsl_matrix * Q_base, gsl_vector * m, gsl_matrix * Q, double tol, int *rankdef);
gsl_matrix *GMRFLib_gsl_duplicate_matrix(gsl_matrix * A);
gsl_matrix *GMRFLib_gsl_transpose_matrix(gsl_matrix * A);
gsl_matrix *GMRFLib_gsl_transpose_matrix_x(gsl_matrix * A, gsl_matrix * At);
gsl_vector *GMRFLib_gsl_duplicate_vector(gsl_vector * a);
double GMRFLib_gsl_rms(gsl_vector * a, gsl_vector * b);

GMRFLib_gsl_low_rank_store_tp *GMRFLib_gsl_low_rank_store_alloc(int n);
int GMRFLib_gsl_low_rank_store_free(GMRFLib_gsl_low_rank_store_tp * S);
gsl_matrix *GMRFLib_gsl_low_rank(gsl_matrix * Cov, double tol, gsl_matrix * B);
gsl_matrix *GMRFLib_gsl_low_rank_x(gsl_matrix * Cov, double tol, gsl_matrix * B, GMRFLib_gsl_low_rank_store_tp * store);

int GMRFLib_gsl_gcpo_singular_fix(int *idx_map, size_t idx_node, gsl_matrix * S, double epsilon);
int GMRFLib_gsl_ginv(gsl_matrix * A, double tol, int rankdef);
int GMRFLib_gsl_mgs(gsl_matrix * A);
int GMRFLib_gsl_mv(gsl_matrix * A, gsl_vector * x, gsl_vector * b);
int GMRFLib_gsl_mm(gsl_matrix * A, gsl_matrix * B, gsl_matrix * C);
int GMRFLib_gsl_mmt(gsl_matrix * A, gsl_matrix * B, gsl_matrix * C);
int GMRFLib_gsl_mmm(gsl_matrix * A, gsl_matrix * B, gsl_matrix * C, gsl_matrix * D);
int GMRFLib_gsl_safe_spd_solve(gsl_matrix * A, gsl_vector * b, gsl_vector * x, double tol);
int GMRFLib_gsl_spd_inv(gsl_matrix * A, double tol);
int GMRFLib_gsl_spd_inverse(gsl_matrix * A);

GMRFLib_gsl_spd_solve_store_tp *GMRFLib_gsl_spd_solve_store_alloc(int n);
int GMRFLib_gsl_spd_solve_store_free(GMRFLib_gsl_spd_solve_store_tp * store);
int GMRFLib_gsl_spd_solve(gsl_matrix * A, gsl_vector * b, gsl_vector * x);
int GMRFLib_gsl_spd_solve_x(gsl_matrix * A, gsl_vector * b, gsl_vector * x, GMRFLib_gsl_spd_solve_store_tp * store);

void GMRFLib_daxpb(int n, double a, double *x, double b, double *y);
void GMRFLib_daxpby(int n, double a, double *x, double b, double *y);
void GMRFLib_daxpbyz(int n, double a, double *x, double b, double *y, double *z);
void GMRFLib_daxpbypcz(int n, double a, double *x, double b, double *y, double c, double *z);
void GMRFLib_daxpy(int n, double a, double *x, double *y);
void GMRFLib_fill(int n, double a, double *x);
void GMRFLib_pack(int n, double *a, int *ia, double *y);
void GMRFLib_unpack(int n, double *a, double *y, int *iy);
void GMRFLib_powx(int n, double *x, double a, double *y);

double GMRFLib_dsum(int n, double *x);
double GMRFLib_dsum_idx(int n, double *a, int *idx);
int GMRFLib_isum(int n, int *ix);

int gsl_blas_dgemm_omp(CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB,
		       double alpha, gsl_matrix * A, gsl_matrix * B, double beta, gsl_matrix * C, int nt);
void cblas_dgemm_omp(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA,
		     enum CBLAS_TRANSPOSE TransB, int M, int N,
		     int K, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nt);

__END_DECLS
#endif
