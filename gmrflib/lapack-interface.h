
/* lapack-interfac.h
 * 
 * Copyright (C) 2001-2023 Havard Rue
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

int GMRFLib_comp_chol_general(double **chol, double *matrix, int dim, double *logdet, int ecode);
int GMRFLib_comp_chol_semidef(double **chol, int **map, int *rank, double *matrix, int dim, double *logdet, double eps);
int GMRFLib_comp_posdef_inverse(double *matrix, int dim);
int GMRFLib_gsl_ensure_spd(gsl_matrix * A, double tol, char **msg);
int GMRFLib_gsl_ensure_spd_inverse(gsl_matrix * A, double tol, char **msg);
int GMRFLib_gsl_ensure_spd_core(gsl_matrix * A, double tol, int method, char **msg);
int GMRFLib_ensure_spd(double *A, int dim, double tol, char **msg);
int GMRFLib_solveAxb_posdef(double *sol, double *chol, double *b, int dim, int nrhs);

int daxpy_(int *n, double *alpha, double *x, int *incx, double *y, int *incy);
int dcopy_(int *n, double *x, int *incx, double *y, int *incy);
int dscal_(int *n, double *alpha, double *x, int *inc);
int dscal_(int *n, double *a, double *x, int *incx);

double ddot_(int *len, double *x, int *incx, double *y, int *incy);
double dnrm2_(int *n, double *x, int *inc);
double dasum_(int *n, double *x, int *inc);

int GMRFLib_dscale(int n, double a, double *x);

double GMRFLib_gsl_spd_logdet(gsl_matrix * A);
double GMRFLib_gsl_xQx(gsl_vector * x, gsl_matrix * Q);
double GMRFLib_gsl_log_dnorm(gsl_vector * x, gsl_vector * mean, gsl_matrix * Q, gsl_matrix * S, int identity);
double GMRFLib_gsl_kld(gsl_vector * m_base, gsl_matrix * Q_base, gsl_vector * m, gsl_matrix * Q, double tol, int *rankdef);
gsl_matrix *GMRFLib_gsl_duplicate_matrix(gsl_matrix * A);
gsl_matrix *GMRFLib_gsl_transpose_matrix(gsl_matrix * A);
gsl_matrix *GMRFLib_gsl_low_rank(gsl_matrix * C, double tol);
int GMRFLib_gsl_gcpo_singular_fix(int *idx_map, size_t idx_node, gsl_matrix * S, double epsilon);
int GMRFLib_gsl_ginv(gsl_matrix * A, double tol, int rankdef);
int GMRFLib_gsl_mgs(gsl_matrix * A);
int GMRFLib_gsl_mv(gsl_matrix * A, gsl_vector * x, gsl_vector * b);
int GMRFLib_gsl_mm(gsl_matrix * A, gsl_matrix * B, gsl_matrix * C);
int GMRFLib_gsl_mmm(gsl_matrix * A, gsl_matrix * B, gsl_matrix * C, gsl_matrix * D);
int GMRFLib_gsl_safe_spd_solve(gsl_matrix * A, gsl_vector * b, gsl_vector * x, double tol);
int GMRFLib_gsl_spd_inv(gsl_matrix * A, double tol);
int GMRFLib_gsl_spd_inverse(gsl_matrix * A);

void GMRFLib_daxpb(int n, double a, double *x, double b, double *y);
void GMRFLib_daxpb2(int n, double a, double *x, double b, double *y);
void GMRFLib_daddto(int n, double *x, double *y);

__END_DECLS
#endif
