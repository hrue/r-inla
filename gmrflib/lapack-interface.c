
/* lapack-interface.c
 * 
 * Copyright (C) 2001-2022 Havard Rue
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
#ifndef GITCOMMIT
#define GITCOMMIT
#endif
static const char GitID[] = "file: " __FILE__ "  " GITCOMMIT;

#if !defined(__FreeBSD__)
#include <malloc.h>
#endif

#include <stdlib.h>
#include <string.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

double GMRFLib_gsl_xQx(gsl_vector * x, gsl_matrix * Q)
{
	size_t n = Q->size1;
	assert(n == Q->size2);
	double sqr, sqr2;

	sqr = 0.0;
	for (size_t i = 0; i < n; i++) {
		double xx = gsl_vector_get(x, i);
		sqr += SQR(xx) * gsl_matrix_get(Q, i, i);
		sqr2 = 0.0;
		for (size_t j = i + 1; j < n; j++) {
			sqr2 += gsl_matrix_get(Q, i, j) * gsl_vector_get(x, j);
		}
		sqr += 2.0 * xx * sqr2;
	}
	return sqr;
}

double GMRFLib_gsl_log_dnorm(gsl_vector * x, gsl_vector * mean, gsl_matrix * Q, gsl_matrix * S)
{
	gsl_matrix *L_Q = NULL, *L_S = NULL;
	double log_det_Q = 0.0;
	size_t n;
	if (Q) {
		n = Q->size1;
		L_Q = GMRFLib_gsl_duplicate_matrix(Q);
		gsl_linalg_cholesky_decomp(L_Q);
		for (size_t i = 0; i < n; i++) {
			log_det_Q += log(gsl_matrix_get(L_Q, i, i));
		}
		log_det_Q *= 2.0;
	} else {
		n = S->size1;
		L_S = GMRFLib_gsl_duplicate_matrix(S);
		gsl_linalg_cholesky_decomp(L_S);
		for (size_t i = 0; i < n; i++) {
			log_det_Q += log(gsl_matrix_get(L_S, i, i));
		}
		log_det_Q *= (-2.0);
	}

	gsl_vector *xx = gsl_vector_alloc(n);
	if (x && mean) {
		for (size_t i = 0; i < n; i++) {
			gsl_vector_set(xx, i, gsl_vector_get(x, i) - gsl_vector_get(mean, i));
		}
	} else if (x && !mean) {
		for (size_t i = 0; i < n; i++) {
			gsl_vector_set(xx, i, gsl_vector_get(x, i));
		}
	} else if (!x && mean) {
		for (size_t i = 0; i < n; i++) {
			gsl_vector_set(xx, i, -gsl_vector_get(mean, i));
		}
	} else {
		for (size_t i = 0; i < n; i++) {
			gsl_vector_set(xx, i, 0.0);
		}
	}

	double sqr = 0.0;
	if (Q) {
		sqr = GMRFLib_gsl_xQx(xx, Q);
	} else {
		// copy from randist/mvgauss.c

		/*
		 * compute: work = L^{-1} * (x - mu) 
		 */
		gsl_blas_dtrsv(CblasLower, CblasNoTrans, CblasNonUnit, L_S, xx);
		/*
		 * compute: quadForm = (x - mu)' Sigma^{-1} (x - mu) 
		 */
		gsl_blas_ddot(xx, xx, &sqr);
	}

	if (L_S) {
		gsl_matrix_free(L_S);
	}
	if (L_Q) {
		gsl_matrix_free(L_Q);
	}
	gsl_vector_free(xx);

	return (-(double) n / 2.0 * log(2.0 * M_PI) + 0.5 * log_det_Q - 0.5 * sqr);
}

int GMRFLib_gsl_gcpo_singular_fix(int *idx_map, size_t idx_node, gsl_matrix * S, double epsilon)
{
	// give a covariance matrix S with mid-node 'idx_node', then detect those indices singular with 'idx_node'.
	// idx_map will map contribution from likelihood to which node.
	// S will be overwritten with a 'fixed' version which is non-singular and partly fake.

#define MAT_SYM_SET(M_, i_, j_, val_)		  \
	if (1) {				  \
		gsl_matrix_set(M_, i_, j_, val_); \
		gsl_matrix_set(M_, j_, i_, val_); \
	}

	assert(S->size1 == S->size2);
	if (S->size1 == 0) {
		return GMRFLib_SUCCESS;
	}

	gsl_matrix *C = GMRFLib_gsl_duplicate_matrix(S);
	for (size_t i = 0; i < S->size1; i++) {
		for (size_t j = i + 1; j < S->size1; j++) {
			double val = gsl_matrix_get(C, i, j) / sqrt(gsl_matrix_get(C, i, i) * gsl_matrix_get(C, j, j));
			val = TRUNCATE(ABS(val), 0.0, 1.0);
			val = (ISEQUAL_x(val, 1.0, epsilon) ? 1.0 : 0.0);
			MAT_SYM_SET(C, i, j, val);
		}
	}

	// first connect to idx_node
	for (size_t i = 0; i < S->size1; i++) {
		idx_map[i] = (int) i;
		if (i != idx_node) {
			if (gsl_matrix_get(C, i, idx_node)) {
				// map contribution to 'idx_node' 
				idx_map[i] = (int) idx_node;
				// then make that node independent of the others so it does not do any harm
				for (size_t ii = 0; ii < S->size1; ii++) {
					if (i != ii) {
						MAT_SYM_SET(S, i, ii, 0.0);
						MAT_SYM_SET(C, i, ii, 0.0);
					}
				}
			}
		}
	}

	// then connect any other nodes that are perfectly correlated
	for (size_t i = 0; i < S->size1; i++) {
		for (size_t j = i + 1; j < S->size1; j++) {
			if (!(i == idx_node || j == idx_node)) {
				if (gsl_matrix_get(C, i, j)) {
					idx_map[j] = (int) i;
					for (size_t jj = 0; jj < S->size1; jj++) {
						if (j != jj) {
							MAT_SYM_SET(S, j, jj, 0.0);
							MAT_SYM_SET(C, j, jj, 0.0);
						}
					}
				}
			}
		}
	}

#undef MAT_SYM_SET
	gsl_matrix_free(C);

	return GMRFLib_SUCCESS;
}

int GMRFLib_gsl_mv(gsl_matrix * A, gsl_vector * x, gsl_vector * b)
{
	// b = A x

	gsl_blas_dgemv(CblasNoTrans, 1.0, (const gsl_matrix *) A, (const gsl_vector *) x, 0.0, b);
	return GMRFLib_SUCCESS;
}

int GMRFLib_comp_posdef_inverse(double *matrix, int dim)
{
	/*
	 * overwrite a symmetric MATRIX with its inverse 
	 */
	int info = 0, i, j;

	switch (GMRFLib_blas_level) {
	case BLAS_LEVEL2:
		dpotf2_("L", &dim, matrix, &dim, &info, F_ONE);
		break;
	case BLAS_LEVEL3:
		dpotrf_("L", &dim, matrix, &dim, &info, F_ONE);
		break;
	default:
		GMRFLib_ASSERT(1 == 0, GMRFLib_ESNH);
		break;
	}
	if (info)
		GMRFLib_ERROR(GMRFLib_ESINGMAT);

	dpotri_("L", &dim, matrix, &dim, &info, F_ONE);
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

	Memcpy(cchol, matrix, ISQR(dim) * sizeof(double));

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
	Memcpy(a, matrix, ISQR(dim) * sizeof(double));

	switch (GMRFLib_blas_level) {
	case BLAS_LEVEL2:
		dpotf2_("L", &dim, a, &dim, &info, F_ONE);
		break;
	case BLAS_LEVEL3:
		dpotrf_("L", &dim, a, &dim, &info, F_ONE);
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
	if (sol != b) {
		Memcpy(sol, b, dim * nrhs * sizeof(double));
	}
	int info;
	dpotrs_("L", &dim, &nrhs, chol, &dim, sol, &dim, &info, F_ONE);
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
	assert(L);
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

	gsl_matrix_set_zero(M1);
	gsl_matrix_set_zero(M2);

	// assert(!(tol > 0.0 || (rankdef >= 0 && rankdef <= (int) A->size1)));
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
				if (i < (size_t) A->size1 - rankdef) {
					gsl_matrix_set(M2, i, i, 1.0 / s);
				} else {
					gsl_matrix_set(M2, i, i, 0.0);
				}
			} else {
				// do not use the first 'rdef's
				if (i < (size_t) rankdef) {
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

int GMRFLib_ensure_spd(double *A, int dim, double tol, char **msg)
{
	// this just a plain interface to the GMRFLib_gsl_ensure_spd

	gsl_matrix *AA = gsl_matrix_alloc((size_t) dim, (size_t) dim);
	size_t i, j;

	for (i = 0; i < (size_t) dim; i++) {
		for (j = 0; j <= i; j++) {
			gsl_matrix_set(AA, i, j, A[i + j * dim]);
			gsl_matrix_set(AA, j, i, A[i + j * dim]);
		}
	}
	GMRFLib_gsl_ensure_spd(AA, tol, msg);
	for (i = 0; i < (size_t) dim; i++) {
		for (j = 0; j <= i; j++) {
			A[i + j * dim] = gsl_matrix_get(AA, i, j);
			A[j + i * dim] = gsl_matrix_get(AA, i, j);
		}
	}
	gsl_matrix_free(AA);
	return GMRFLib_SUCCESS;
}

int GMRFLib_gsl_ensure_spd(gsl_matrix * A, double tol, char **msg)
{
	return GMRFLib_gsl_ensure_spd_core(A, tol, 0, msg);
}

int GMRFLib_gsl_ensure_spd_inverse(gsl_matrix * A, double tol, char **msg)
{
	return GMRFLib_gsl_ensure_spd_core(A, tol, 1, msg);
}

int GMRFLib_gsl_ensure_spd_core(gsl_matrix * A, double tol, int method, char **msg)
{
	/*
	 * replace n x n matrix A with its SPD matrix, replacing small eigenvalues with 'tol' * max(|eigenvalue|).
	 *
	 * if tol < 0, then use tol = (min(eigenvalue > 0) / max(eigenvalue)
	 *
	 * method=0: s=max(min,s)
	 * method=1: s=1/max(min,s)
	 */

	assert(A && (A->size1 == A->size2));
	gsl_matrix *U = GMRFLib_gsl_duplicate_matrix(A);
	gsl_vector *S = gsl_vector_alloc(A->size1);
	gsl_eigen_symmv_workspace *work = gsl_eigen_symmv_alloc(A->size1);

	gsl_eigen_symmv(A, S, U, work);

	size_t i;
	double one = 1.0, zero = 0.0, s;

	double s_min, s_max = gsl_vector_max(S);
	if (s_max <= 0.0) {
		s_max = 0.0;				       /* then the whole matrix is zero or INF, as all is wrong... */
	}

	if (tol < 0.0 && s_max > 0.0) {
		int n_neg = 0;
		s_min = s_max;
		for (i = 0; i < A->size1; i++) {
			double s = gsl_vector_get(S, i);
			if (s <= 0.0) {
				n_neg++;
			}
			if (s > 0.0 && s < s_min) {
				s_min = s;
			}
		}
		tol = s_min / s_max;
		if (msg) {
			GMRFLib_sprintf(msg, "set tol=[%.4g]. number of negative eigenvalues=[%1d]", tol, n_neg);
		}
	}

	gsl_matrix *M1 = gsl_matrix_alloc(A->size1, A->size2);
	gsl_matrix *M2 = gsl_matrix_alloc(A->size1, A->size2);

	gsl_matrix_set_zero(M1);
	gsl_matrix_set_zero(M2);

	s_min = tol * s_max;
	if (method == 0) {
		for (i = 0; i < A->size1; i++) {
			s = gsl_vector_get(S, i);
			gsl_matrix_set(M2, i, i, DMAX(s_min, s));
		}
	} else {
		for (i = 0; i < A->size1; i++) {
			s = gsl_vector_get(S, i);
			gsl_matrix_set(M2, i, i, 1.0 / DMAX(s_min, s));
		}
	}

	gsl_blas_dgemm(CblasNoTrans, CblasTrans, one, M2, U, zero, M1);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, one, U, M1, zero, M2);
	gsl_matrix_memcpy(A, M2);

	gsl_matrix_free(U);
	gsl_matrix_free(M1);
	gsl_matrix_free(M2);
	gsl_vector_free(S);
	gsl_eigen_symmv_free(work);

	return GMRFLib_SUCCESS;
}

int GMRFLib_gsl_safe_spd_solve(gsl_matrix * A, gsl_vector * b, gsl_vector * x, double tol)
{
	/*
	 * solve Ax=b, ignoring contributions from eigenvalues < tol*max(eigenval)
	 *
	 * solution is returned in x, while A and b is not changed.
	 *
	 */

	int debug = 0;
	assert(A && (A->size1 == A->size2));
	assert(tol >= 0.0);

	gsl_matrix *U = GMRFLib_gsl_duplicate_matrix(A);
	gsl_vector *S = gsl_vector_alloc(A->size1);

	gsl_eigen_symmv_workspace *work = gsl_eigen_symmv_alloc(A->size1);
	gsl_eigen_symmv(A, S, U, work);

	size_t i;
	double one = 1.0, zero = 0.0, s;
	double s_max = ABS(gsl_vector_max(S));
	gsl_matrix *M1 = gsl_matrix_alloc(A->size1, A->size2);
	gsl_matrix *M2 = gsl_matrix_alloc(A->size1, A->size2);

	gsl_matrix_set_zero(M1);
	gsl_matrix_set_zero(M2);

	double s_min = tol * s_max;
	if (debug && !(s_max > 0.0)) {
		FIXME("s_max > 0 FAILED");
		P(s_max);
		GMRFLib_printf_gsl_matrix(stdout, A, " %g");
	}

	for (i = 0; i < A->size1; i++) {
		s = gsl_vector_get(S, i);
		if (s <= s_min) {
			s = 0.0;
		} else {
			s = 1.0 / s;
		}
		gsl_matrix_set(M2, i, i, s);
	}

	gsl_blas_dgemm(CblasNoTrans, CblasTrans, one, M2, U, zero, M1);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, one, U, M1, zero, M2);

	if (x == b) {
		gsl_vector *xx = gsl_vector_alloc(A->size1);
		gsl_blas_dgemv(CblasNoTrans, one, M2, b, zero, xx);
		gsl_vector_memcpy(x, xx);
		gsl_vector_free(xx);
	} else {
		gsl_blas_dgemv(CblasNoTrans, one, M2, b, zero, x);
	}

	gsl_matrix_free(U);
	gsl_matrix_free(M1);
	gsl_matrix_free(M2);
	gsl_vector_free(S);
	gsl_eigen_symmv_free(work);

	return GMRFLib_SUCCESS;
}

int GMRFLib_gsl_spd_inv(gsl_matrix * A, double tol)
{
	/*
	 * A=inv(A) for symmetric A, ignoring contributions from eigenvalues < tol*max(eigenval)
	 */

	int debug = 0;
	assert(A && (A->size1 == A->size2));
	assert(tol >= 0.0);

	gsl_matrix *U = GMRFLib_gsl_duplicate_matrix(A);
	gsl_vector *S = gsl_vector_alloc(A->size1);

	gsl_eigen_symmv_workspace *work = gsl_eigen_symmv_alloc(A->size1);
	gsl_eigen_symmv(A, S, U, work);

	size_t i;
	double one = 1.0, zero = 0.0, s;
	double s_max = ABS(gsl_vector_max(S));
	gsl_matrix *M1 = gsl_matrix_alloc(A->size1, A->size2);
	gsl_matrix *M2 = gsl_matrix_alloc(A->size1, A->size2);
	gsl_matrix_set_zero(M1);
	gsl_matrix_set_zero(M2);

	double s_min = tol * s_max;
	if (debug && !(s_max > 0.0)) {
		FIXME("s_max > 0 FAILED");
		P(s_max);
		GMRFLib_printf_gsl_matrix(stdout, A, " %g");
	}

	for (i = 0; i < A->size1; i++) {
		s = gsl_vector_get(S, i);
		if (s <= s_min) {
			s = 0.0;
		} else {
			s = 1.0 / s;
		}
		gsl_matrix_set(M2, i, i, s);
	}

	gsl_blas_dgemm(CblasNoTrans, CblasTrans, one, M2, U, zero, M1);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, one, U, M1, zero, A);

	gsl_matrix_free(U);
	gsl_matrix_free(M1);
	gsl_matrix_free(M2);
	gsl_vector_free(S);
	gsl_eigen_symmv_free(work);

	return GMRFLib_SUCCESS;
}

int GMRFLib_gsl_mgs(gsl_matrix * A)
{
	// this is the modified Gram-Schmith ortogonalisation, and it 
	// overwrite matrix A with its ortogonal normalized basis, column by column.
	// the sign of each column is so that max(abs(column)) is positive

	gsl_vector *q = gsl_vector_alloc(A->size1);
	size_t i, j, k, n = A->size1;

	for (i = 0; i < n; i++) {

		double aij_amax = 0.0, r = 0.0;
		for (j = 0; j < n; j++) {
			double elm = gsl_matrix_get(A, j, i);
			r += SQR(elm);
			aij_amax = (ABS(elm) > ABS(aij_amax) ? elm : aij_amax);
		}

		if (aij_amax < 0.0) {			       /* swap the sign of this column */
			for (j = 0; j < n; j++) {
				gsl_matrix_set(A, j, i, -gsl_matrix_get(A, j, i));
			}
		}

		r = sqrt(r);
		for (j = 0; j < n; j++) {
			gsl_vector_set(q, j, gsl_matrix_get(A, j, i) / r);
			gsl_matrix_set(A, j, i, gsl_vector_get(q, j));	// normalize
		}

		for (j = i + 1; j < n; j++) {
			for (r = 0, k = 0; k < n; k++) {
				r += gsl_vector_get(q, k) * gsl_matrix_get(A, k, j);
			}
			for (k = 0; k < n; k++) {
				gsl_matrix_set(A, k, j, gsl_matrix_get(A, k, j) - r * gsl_vector_get(q, k));
			}
		}
	}
	gsl_vector_free(q);

	return (GMRFLib_SUCCESS);
}
