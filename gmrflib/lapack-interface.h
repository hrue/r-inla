
/* lapack-interfac.h
 * 
 * Copyright (C) 2001-2006 Havard Rue
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

#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
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

/*!
  \brief Define BLAS_LEVEL2
*/
#define BLAS_LEVEL2 2

/*!
  \brief Define BLAS_LEVEL3
*/
#define BLAS_LEVEL3 3
int dpbtrf_(const char *, int *, int *, double *, int *, int *, int);
int dpbtf2_(const char *, int *, int *, double *, int *, int *, int);
int dpotrf_(const char *, int *, double *, int *, int *, int);
int dpotf2_(const char *, int *, double *, int *, int *, int);
int dtbsv_(const char *, const char *, const char *, int *, int *, double *, int *, double *, int *, int, int, int);
int dpotri_(const char *, int *, double *, int *, int *, int);
int dgemm_(const char *, const char *, int *, int *, int *, double *, double *, int *, double *, int
	   *, double *, double *, int *, int, int);
int dgemv_(const char *, int *, int *, double *, double *, int *, double *, int *, double *, double
	   *, int *, int);
int dgemv_failsafe_(const char *, int *, int *, double *, double *, int *, double *, int *, double *, double
		    *, int *, int);
double erfi_(double *);
int dpotrs_(const char *, int *, int *, double *, int *, double *, int *, int *, int);
int dchdc_(double *, int *, int *, double *, int *, int *, int *, double *);
int dtrmv_(const char *, const char *, const char *, int *, double *, int *, double *, int *, int, int, int);

int GMRFLib_comp_chol_semidef(double **chol, int **map, int *rank, double *matrix, int dim, double *logdet, double eps);
int GMRFLib_comp_posdef_inverse(double *matrix, int dim);
int GMRFLib_comp_chol_general(double **chol, double *matrix, int dim, double *logdet, int ecode);
int GMRFLib_solveAxb_posdef(double *sol, double *chol, double *b, int dim, int nrhs);

int daxpy_(int *n, double *alpha, double *x, int *incx, double *y, int *incy);
int dcopy_(int *n, double *x, int *incx, double *y, int *incy);
int dscal_(int *n, double *alpha, double *x, int *inc);
double ddot_(int *n, double *x, int *incx, double *y, int *incy);
double dnrm2_(int *n, double *x, int *inc);

gsl_matrix *GMRFLib_gsl_duplicate_matrix(gsl_matrix * A);
double GMRFLib_gsl_spd_logdet(gsl_matrix * A);
int GMRFLib_gsl_spd_inverse(gsl_matrix * A);
int GMRFLib_gsl_ginv(gsl_matrix * A, double tol, int rankdef);

__END_DECLS
#endif
