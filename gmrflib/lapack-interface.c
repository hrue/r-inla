
/* lapack-interface.c
 * 
 * Copyright (C) 2001-2013 Havard Rue
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

/*!
  \file lapack-interface.c
  \brief The interface towards the LAPACK routines written in fortran.
*/
#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

/* Pre-hg-Id: $Id: lapack-interface.c,v 1.23 2009/05/02 16:54:04 hrue Exp $ */

#if !defined(__FreeBSD__)
#include <malloc.h>
#endif

#include <stdlib.h>
#include <string.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

int GMRFLib_comp_posdef_inverse(double *matrix, int dim)
{
	/*
	 * overwrite a symmetric MATRIX with its inverse 
	 */
	int info = 0, i, j;

	switch (GMRFLib_blas_level) {
	case BLAS_LEVEL2:
		dpotf2_("L", &dim, matrix, &dim, &info, 1);
		break;
	case BLAS_LEVEL3:
		dpotrf_("L", &dim, matrix, &dim, &info, 1);
		break;
	default:
		GMRFLib_ASSERT(1 == 0, GMRFLib_ESNH);
		break;
	}
	if (info)
		GMRFLib_ERROR(GMRFLib_ESINGMAT);

	dpotri_("L", &dim, matrix, &dim, &info, 1);
	if (info)
		GMRFLib_ERROR(GMRFLib_ESINGMAT);

	for (i = 0; i < dim; i++)			       /* fill the U-part */
		for (j = i + 1; j < dim; j++)
			matrix[i + j * dim] = matrix[j + i * dim];

	return GMRFLib_SUCCESS;
}
int GMRFLib_comp_chol_semidef(double **chol, int **map, int *rank, double *matrix, int dim, double *logdet, double eps)
{
	/*
	 * compute the ``cholesky factorisation'' for a positive semidefinite matrix. return malloc'ed factorization in chol
	 * (if !NULL), the malloc'ed mapping in map and the rank in *rank.
	 * 
	 * if logdet, then compute the log determinant of the non-singular part
	 * 
	 * eps is the smalles L[i,i] 
	 */

	double *work = NULL, det, *cchol = NULL;
	int job = 1, info, i;

	cchol = Calloc(ISQR(dim), double);
	*map = Calloc(dim, int);
	work = Calloc(dim, double);

	memcpy(cchol, matrix, ISQR(dim) * sizeof(double));

	dchdc_(cchol, &dim, &dim, work, *map, &job, &info, &eps);
	*rank = info;

	for (i = 0; i < dim; i++) {
		(*map)[i]--;				       /* convert to C index-ing */
	}
	if (logdet) {
		for (det = 0.0, i = 0; i < *rank; i++) {
			det += log(cchol[i + i * dim]);
		}
		*logdet = 2.0 * det;
	}

	Free(work);

	if (chol) {
		*chol = cchol;
	} else {
		Free(cchol);
	}

	return GMRFLib_SUCCESS;
}
int GMRFLib_comp_chol_general(double **chol, double *matrix, int dim, double *logdet, int ecode)
{
	/*
	 * return a malloc'ed cholesky factorisation of MATRIX in *chol and optional the log(determinant). if fail return
	 * `ecode'
	 * 
	 */
	int info = 0, i, j;
	double *a = NULL, det;

	if (0) {
		P(dim);
		for (i = 0; i < dim; i++)
			for (j = 0; j < dim; j++)
				printf("i %d j %d matrix %.12g\n", i, j, matrix[i + dim * j]);
	}

	if (dim == 0) {
		*chol = NULL;
		return GMRFLib_SUCCESS;
	}

	a = Calloc(ISQR(dim), double);
	memcpy(a, matrix, ISQR(dim) * sizeof(double));

	switch (GMRFLib_blas_level) {
	case BLAS_LEVEL2:
		dpotf2_("L", &dim, a, &dim, &info, 1);
		break;
	case BLAS_LEVEL3:
		dpotrf_("L", &dim, a, &dim, &info, 1);
		break;
	default:
		GMRFLib_ASSERT(1 == 0, GMRFLib_ESNH);
		break;
	}

	if (info) {
		Free(a);
		*chol = NULL;

		return ecode;
	}

	if (logdet) {
		for (det = 0.0, i = 0; i < dim; i++) {
			det += log(a[i + i * dim]);
		}
		*logdet = 2.0 * det;
	}

	for (i = 0; i < dim; i++) {			       /* set to zero the upper part */
		for (j = i + 1; j < dim; j++) {
			a[i + j * dim] = 0.0;
		}
	}

	*chol = a;
	return GMRFLib_SUCCESS;
}
int GMRFLib_solveAxb_posdef(double *sol, double *chol, double *b, int dim, int nrhs)
{
	/*
	 * solve Ax=b, where chol is lower Cholesky factor of A. 
	 */
	int info;

	if (sol != b) {
		memcpy(sol, b, dim * nrhs * sizeof(double));
	}
	dpotrs_("L", &dim, &nrhs, chol, &dim, sol, &dim, &info, 1);

	if (info) {
		GMRFLib_ERROR(GMRFLib_EPOSDEF);
	}

	return GMRFLib_SUCCESS;
}
gsl_matrix *GMRFLib_gsl_duplicate_matrix(gsl_matrix * A)
{
	/*
	 * return a new (alloced) copy of matrix A 
	 */
	if (A) {
		gsl_matrix *B = gsl_matrix_alloc(A->size1, A->size2);
		gsl_matrix_memcpy(B, A);

		return B;
	} else {
		return (gsl_matrix *) NULL;
	}
}
double GMRFLib_gsl_spd_logdet(gsl_matrix * A)
{
	/*
	 * compute the log determinant of a SPD matrix A 
	 */
	gsl_matrix *L;
	double logdet = 0.0;
	size_t i;

	L = GMRFLib_gsl_duplicate_matrix(A);
	gsl_linalg_cholesky_decomp(L);

	for (i = 0; i < L->size1; i++) {
		logdet += log(gsl_matrix_get(L, i, i));
	}
	logdet *= 2.0;					       /* |A| = |L|^2 */

	gsl_matrix_free(L);

	return logdet;
}
int GMRFLib_gsl_spd_inverse(gsl_matrix * A)
{
	/*
	 * replace SPD matrix A with its inverse 
	 */
	gsl_matrix *L;
	gsl_vector *x;
	size_t i, n;

	assert(A->size1 == A->size2);
	n = A->size1;

	x = gsl_vector_calloc(n);
	L = GMRFLib_gsl_duplicate_matrix(A);
	gsl_linalg_cholesky_decomp(L);
	for (i = 0; i < n; i++) {
		gsl_vector_set_basis(x, i);
		gsl_linalg_cholesky_svx(L, x);
		gsl_matrix_set_col(A, i, x);
	}

	gsl_vector_free(x);
	gsl_matrix_free(L);

	return GMRFLib_SUCCESS;
}
int GMRFLib_gsl_ginv(gsl_matrix * A, double tol, int rankdef)
{
	/*
	 * replace n x n matrix A with its generlized inverse.  if TOL > 0, use that tolerance. If rankdef is set, use that. If both are set, give an error.
	 */

	assert(A && (A->size1 == A->size2));

	gsl_matrix *U = GMRFLib_gsl_duplicate_matrix(A);
	gsl_matrix *V = gsl_matrix_alloc(A->size1, A->size2);
	gsl_vector *S = gsl_vector_alloc(A->size1);
	gsl_vector *work = gsl_vector_alloc(A->size1);

	gsl_linalg_SV_decomp(U, V, S, work);

	size_t i;
	double one = 1.0, zero = 0.0;
	double s_max = gsl_vector_get(S, 0);
	gsl_matrix *M1 = gsl_matrix_alloc(A->size1, A->size2);
	gsl_matrix *M2 = gsl_matrix_alloc(A->size1, A->size2);

	assert(!(tol > 0.0 && (rankdef >= 0 && rankdef <= (int) A->size1)));
	if (tol > 0.0) {
		for (i = 0; i < A->size1; i++) {
			double s = gsl_vector_get(S, i);

			if (s < tol * s_max) {
				gsl_matrix_set(M2, i, i, 0.0);
			} else {
				gsl_matrix_set(M2, i, i, 1.0 / s);
			}
		}
	} else {
		assert(rankdef >= 0);
		assert(rankdef <= (int) A->size1);

		double first = gsl_vector_get(S, 0);
		double last = gsl_vector_get(S, A->size1 - 1);

		for (i = 0; i < A->size1; i++) {
			double s = gsl_vector_get(S, i);

			if (first > last) {
				// do not use the last 'rdef's
				if (i < (int) A->size1 - rankdef) {
					gsl_matrix_set(M2, i, i, 1.0 / s);
				} else {
					gsl_matrix_set(M2, i, i, 0.0);
				}
			} else {
				// do not use the first 'rdef's
				if (i < rankdef) {
					gsl_matrix_set(M2, i, i, 0.0);
				} else {
					gsl_matrix_set(M2, i, i, 1.0 / s);
				}
			}
		}
	}

	gsl_blas_dgemm(CblasNoTrans, CblasTrans, one, M2, U, zero, M1);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, one, V, M1, zero, M2);
	gsl_matrix_memcpy(A, M2);

	gsl_matrix_free(U);
	gsl_matrix_free(V);
	gsl_matrix_free(M1);
	gsl_matrix_free(M2);
	gsl_vector_free(S);
	gsl_vector_free(work);

	return GMRFLib_SUCCESS;
}
