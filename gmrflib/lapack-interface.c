#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <omp.h>

#include "GMRFLib/GMRFLib.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((optimize("O3")))
    __attribute__((target_clones(INLA_CLONE_TARGETS "default")))
void GMRFLib_gsl_dgemm_sym(gsl_matrix *A, gsl_matrix *B, gsl_matrix *C)
{
	// this is a special case.

	// A: n x m, B: m x n, C: n x n and symmetric. all row-major
#define SIZE 128
	int n = (int) A->size1;
	int m = (int) A->size2;
	assert((int) A->tda == m);
	assert((int) B->tda == n);
	assert((int) C->tda == n);
	gsl_matrix_set_zero(C);

	int block_m;
	if (m >= 2 * SIZE) {
		block_m = IMIN(m, SIZE);
	} else if (m > SIZE) {
		block_m = (m + 1) / 2;
	} else {
		block_m = m;
	}

	int block_n;
	if (n >= 2 * SIZE) {
		block_n = IMIN(n, SIZE);
	} else if (n > SIZE) {
		block_n = (n + 1) / 2;
	} else {
		block_n = n;
	}

	for (int ii = 0; ii < n; ii += block_n) {
		int ni = (ii + block_n < n) ? block_n : n - ii;
		for (int jj = ii; jj < n; jj += block_n) {
			int nj = (jj + block_n < n) ? block_n : n - jj;
			for (int kk = 0; kk < m; kk += block_m) {
				int nk = (kk + block_m < m) ? block_m : m - kk;
				cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, ni, nj, nk, 1.0,
					    &A->data[ii * m + kk], m, &B->data[kk * n + jj], n, 1.0, &C->data[ii * n + jj], n);
			}
			// this does not ensure _exact_ symmetry for ii==jj due to numerical sum error, but we can call
			// GMRFLib_gsl_force_symmetric() to do that
			if (ii != jj) {
				for (int p = 0; p < ni; p++) {
					for (int q = 0; q < nj; q++) {
						C->data[(jj + q) * n + (ii + p)] = C->data[(ii + p) * n + (jj + q)];
					}
				}
			}
		}
	}
#undef SIZE
}
#pragma GCC diagnostic pop


#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((optimize("O3")))
    __attribute__((target_clones(INLA_CLONE_TARGETS "default")))
double GMRFLib_gsl_xQx(gsl_vector *x, gsl_matrix *Q)
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
#pragma GCC diagnostic pop

GMRFLib_gsl_ldnorm_store_tp *GMRFLib_gsl_ldnorm_store_alloc(int n)
{
	GMRFLib_gsl_ldnorm_store_tp *S = Calloc(1, GMRFLib_gsl_ldnorm_store_tp);
	S->L = gsl_matrix_alloc(n, n);
	S->xx = gsl_vector_alloc(n);
	return S;
}
int GMRFLib_gsl_ldnorm_store_free(GMRFLib_gsl_ldnorm_store_tp *store)
{
	if (store) {
		gsl_matrix_free(store->L);
		gsl_vector_free(store->xx);
		Free(store);
	}
	return 0;
}

double GMRFLib_gsl_ldnorm(gsl_vector *x, gsl_vector *mean, gsl_matrix *Q, gsl_matrix *S, int identity)
{
	return GMRFLib_gsl_ldnorm_x(x, mean, Q, S, identity, NULL);
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((optimize("O3")))
    __attribute__((target_clones(INLA_CLONE_TARGETS "default")))
double GMRFLib_gsl_ldnorm_x(gsl_vector *x, gsl_vector *mean, gsl_matrix *Q, gsl_matrix *S, int identity, GMRFLib_gsl_ldnorm_store_tp *store)
{
	// 'identity' says that Q=S=I

	gsl_matrix *L_Q = NULL, *L_S = NULL;
	double log_det_Q = 0.0;
	size_t n = 0;

	if (identity) {
		if (Q) {
			n = Q->size1;
		} else if (S) {
			n = S->size1;
		} else if (mean) {
			n = mean->size;
		} else if (x) {
			n = x->size;
		} else {
			assert(0 == 1);
		}
		log_det_Q = 1.0;
	} else {
		if (Q) {
			n = Q->size1;
			if (store) {
				store->L->size1 = n;
				store->L->size2 = n;
				L_Q = store->L;
				gsl_matrix_memcpy(L_Q, Q);
			} else {
				L_Q = GMRFLib_gsl_duplicate_matrix(Q);
			}
			gsl_linalg_cholesky_decomp(L_Q);
			for (size_t i = 0; i < n; i++) {
				log_det_Q += log(gsl_matrix_get(L_Q, i, i));
			}
			log_det_Q *= 2.0;
		} else {
			n = S->size1;
			if (store) {
				store->L->size1 = n;
				store->L->size2 = n;
				L_S = store->L;
				gsl_matrix_memcpy(L_S, S);
			} else {
				L_S = GMRFLib_gsl_duplicate_matrix(S);
			}
			gsl_linalg_cholesky_decomp(L_S);
			for (size_t i = 0; i < n; i++) {
				log_det_Q += log(gsl_matrix_get(L_S, i, i));	/* yes */
			}
			log_det_Q *= (-2.0);		       /* yes */
		}
	}

	gsl_vector *xx = NULL;
	if (store) {
		store->xx->size = n;
		xx = store->xx;
	} else {
		xx = gsl_vector_alloc(n);
	}

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
		gsl_vector_set_zero(xx);
	}

	double sqr = 0.0;
	if (identity) {
		gsl_blas_ddot(xx, xx, &sqr);
	} else {
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
	}

	if (!store) {
		if (L_S) {
			gsl_matrix_free(L_S);
		}
		if (L_Q) {
			gsl_matrix_free(L_Q);
		}
		gsl_vector_free(xx);
	}

	return ((-(double) n * 1.83787706640934548356065947281 + log_det_Q - sqr) * 0.5);
}
#pragma GCC diagnostic pop

int GMRFLib_gsl_gcpo_singular_fix(int *idx_map, size_t idx_node, gsl_matrix *S, double epsilon)
{
	// give a covariance matrix S with mid-node 'idx_node', then detect those indices singular with 'idx_node'.
	// idx_map will map contribution from likelihood to that node.
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

int GMRFLib_gsl_mv(gsl_matrix *A, gsl_vector *x, gsl_vector *b)
{
	// b = A x

	gsl_blas_dgemv(CblasNoTrans, 1.0, (const gsl_matrix *) A, (const gsl_vector *) x, 0.0, b);
	return GMRFLib_SUCCESS;
}

int GMRFLib_gsl_mm(gsl_matrix *A, gsl_matrix *B, gsl_matrix *C)
{
	// C = A B
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, (const gsl_matrix *) A, (const gsl_matrix *) B, 0.0, (gsl_matrix *) C);
	return GMRFLib_SUCCESS;
}

int GMRFLib_gsl_mmt(gsl_matrix *A, gsl_matrix *B, gsl_matrix *C)
{
	// C = A t(B)
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, (const gsl_matrix *) A, (const gsl_matrix *) B, 0.0, (gsl_matrix *) C);
	return GMRFLib_SUCCESS;
}

int GMRFLib_gsl_mmm(gsl_matrix *A, gsl_matrix *B, gsl_matrix *C, gsl_matrix *D)
{
	// D = A B C
	gsl_matrix *T = gsl_matrix_alloc(A->size1, B->size2);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, (const gsl_matrix *) A, (const gsl_matrix *) B, 0.0, (gsl_matrix *) T);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, (const gsl_matrix *) T, (const gsl_matrix *) C, 0.0, (gsl_matrix *) D);
	gsl_matrix_free(T);
	return GMRFLib_SUCCESS;
}

int GMRFLib_comp_posdef_inverse(double *matrix, int dim)
{
	/*
	 * overwrite a symmetric MATRIX with its inverse 
	 */
	int info = 0;

	// dpotf2_("L", &dim, matrix, &dim, &info, F_ONE);
	dpotrf_("L", &dim, matrix, &dim, &info, F_ONE);

	if (info)
		GMRFLib_ERROR(GMRFLib_ESINGMAT);

	dpotri_("L", &dim, matrix, &dim, &info, F_ONE);
	if (info)
		GMRFLib_ERROR(GMRFLib_ESINGMAT);

	for (int i = 0; i < dim; i++)			       /* fill the U-part */
		for (int j = i + 1; j < dim; j++)
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

	double *work = NULL, *cchol = NULL;
	int job = 1, info;

	cchol = Calloc(ISQR(dim), double);
	*map = Calloc(dim, int);
	work = Calloc(dim, double);
	Memcpy(cchol, matrix, ISQR(dim) * sizeof(double));

	dchdc_(cchol, &dim, &dim, work, *map, &job, &info, &eps);
	*rank = info;

	for (int i = 0; i < dim; i++) {
		(*map)[i]--;				       /* convert to C index-ing */
	}
	if (logdet) {
		double ldet = 0.0;
		for (int i = 0; i < *rank; i++) {
			ldet += log(cchol[i + i * dim]);
		}
		*logdet = 2.0 * ldet;
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
	int info = 0;
	double *a = NULL;

	if (0) {
		P(dim);
		for (int i = 0; i < dim; i++)
			for (int j = 0; j < dim; j++)
				printf("i %d j %d matrix %.12g\n", i, j, matrix[i + dim * j]);
	}

	if (dim == 0) {
		*chol = NULL;
		return GMRFLib_SUCCESS;
	}

	a = Calloc(ISQR(dim), double);
	Memcpy(a, matrix, ISQR(dim) * sizeof(double));

	// dpotf2_("L", &dim, a, &dim, &info, F_ONE);
	dpotrf_("L", &dim, a, &dim, &info, F_ONE);

	if (info) {
		Free(a);
		*chol = NULL;

		return ecode;
	}

	if (logdet) {
		double ldet = 0.0;
		for (int i = 0; i < dim; i++) {
			ldet += log(a[i + i * dim]);
		}
		*logdet = 2.0 * ldet;
	}

	for (int i = 0; i < dim; i++) {			       /* set to zero the upper part */
		for (int j = i + 1; j < dim; j++) {
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

gsl_matrix *GMRFLib_gsl_duplicate_matrix(gsl_matrix *A)
{
	/*
	 * return a new (alloced) copy of matrix A 
	 */
	gsl_matrix *B = NULL;
	if (A) {
		B = gsl_matrix_alloc(A->size1, A->size2);
		gsl_matrix_memcpy(B, A);
	}
	return B;
}

gsl_vector *GMRFLib_gsl_duplicate_vector(gsl_vector *a)
{
	/*
	 * return a new (alloced) copy of vector 'a'
	 */
	gsl_vector *b = NULL;
	if (a) {
		b = gsl_vector_alloc(a->size);
		gsl_vector_memcpy(b, a);
	}
	return b;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((optimize("O3")))
    __attribute__((target_clones(INLA_CLONE_TARGETS "default")))
double GMRFLib_gsl_rms(gsl_vector *a, gsl_vector *b)
{
	double rms = 0.0;

	if (!a && !b) {
		return 0.0;
	} else if (!a && b) {
		return GMRFLib_gsl_rms(b, NULL);
	} else if (a && !b) {
		for (size_t i = 0; i < a->size; i++) {
			rms += SQR(gsl_vector_get(a, i));
		}
	} else if (a && b) {
		for (size_t i = 0; i < a->size; i++) {
			rms += SQR(gsl_vector_get(a, i) - gsl_vector_get(b, i));
		}
	} else {
		assert(0 == 1);
	}
	return sqrt(rms / a->size);
}
#pragma GCC diagnostic pop

gsl_matrix *GMRFLib_gsl_transpose_matrix(gsl_matrix *A)
{
	/*
	 * return a new (alloced) t(A)
	 */
	gsl_matrix *At = NULL;

	if (!A) {
		return At;
	}

	if (A->size1 == A->size2) {
		At = GMRFLib_gsl_duplicate_matrix(A);
		gsl_matrix_transpose(At);
	} else {
		At = gsl_matrix_alloc(A->size2, A->size1);
		for (size_t i = 0; i < A->size1; i++) {
			for (size_t j = 0; j < A->size2; j++) {
				gsl_matrix_set(At, j, i, gsl_matrix_get(A, i, j));
			}
		}
	}

	return At;
}

gsl_matrix *GMRFLib_gsl_transpose_matrix_x(gsl_matrix *A, gsl_matrix *At)
{
	// overwrite At with the transpose
	if (!A) {
		return (At ? NULL : At);
	}

	gsl_matrix *AAt = (At ? At : NULL);
	if (!AAt) {
		AAt = gsl_matrix_alloc(A->size2, A->size1);
	} else {
		AAt->size1 = A->size2;
		AAt->size2 = A->size1;
	}
	for (size_t i = 0; i < A->size1; i++) {
		for (size_t j = 0; j < A->size2; j++) {
			gsl_matrix_set(AAt, j, i, gsl_matrix_get(A, i, j));
		}
	}

	return (At ? NULL : At);
}

double GMRFLib_gsl_spd_logdet(gsl_matrix *A)
{
	/*
	 * compute the log determinant of a SPD matrix A 
	 */
	gsl_matrix *L = NULL;
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

int GMRFLib_gsl_force_symmetric(gsl_matrix *A)
{
	// force symmetry A = (A + t(A))/2
	for (size_t i = 0; i < A->size1; i++) {
		for (size_t j = 0; j < i; j++) {
			double a = gsl_matrix_get(A, i, j);
			double aa = gsl_matrix_get(A, j, i);
			double val = (a + aa) / 2.0;
			gsl_matrix_set(A, i, j, val);
			gsl_matrix_set(A, j, i, val);
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_gsl_add_diag(gsl_matrix *A, double value)
{
	// diag(A) = diag(A) + value
	if (!ISZERO(value)) {
		for (size_t i = 0; i < A->size1; i++) {
			double a = gsl_matrix_get(A, i, i) + value;
			gsl_matrix_set(A, i, i, a);
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_gsl_spd_inverse(gsl_matrix *A)
{
	/*
	 * replace SPD matrix A with its inverse 
	 */
	gsl_matrix *L = NULL;
	gsl_vector *x = NULL;
	size_t i, n;

	assert(A->size1 == A->size2);
	n = A->size1;

	L = GMRFLib_gsl_duplicate_matrix(A);
	int ecode = gsl_linalg_cholesky_decomp(L);
	if (ecode != GSL_SUCCESS) {
		gsl_matrix_free(L);
		return !GMRFLib_SUCCESS;
	}

	x = gsl_vector_calloc(n);
	for (i = 0; i < n; i++) {
		gsl_vector_set_basis(x, i);
		gsl_linalg_cholesky_svx(L, x);
		gsl_matrix_set_col(A, i, x);
	}

	gsl_vector_free(x);
	gsl_matrix_free(L);

	return GMRFLib_SUCCESS;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((optimize("O3")))
    __attribute__((target_clones(INLA_CLONE_TARGETS "default")))
int GMRFLib_gsl_ginv(gsl_matrix *A, double tol, int rankdef)
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
	const double one = 1.0, zero = 0.0;
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
#pragma GCC diagnostic pop

int GMRFLib_ensure_spd(double *A, int dim, double tol, char **msg)
{
	return GMRFLib_ensure_spd_x(A, dim, tol, msg, NULL);
}

int GMRFLib_ensure_spd_x(double *A, int dim, double tol, char **msg, GMRFLib_gsl_ensure_spd_store_tp *store)
{
	// this just a plain interface to the GMRFLib_gsl_ensure_spd

	gsl_matrix *AA = NULL;

	if (store) {
		AA = store->AA;
	} else {
		AA = gsl_matrix_alloc((size_t) dim, (size_t) dim);
	}
	size_t i, j;

	for (i = 0; i < (size_t) dim; i++) {
		for (j = 0; j <= i; j++) {
			gsl_matrix_set(AA, i, j, A[i + j * dim]);
			gsl_matrix_set(AA, j, i, A[i + j * dim]);
		}
	}
	GMRFLib_gsl_ensure_spd_x(AA, tol, msg, store);
	for (i = 0; i < (size_t) dim; i++) {
		for (j = 0; j <= i; j++) {
			A[i + j * dim] = gsl_matrix_get(AA, i, j);
			A[j + i * dim] = gsl_matrix_get(AA, i, j);
		}
	}

	if (!store) {
		gsl_matrix_free(AA);
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_gsl_ensure_spd(gsl_matrix *A, double tol, char **msg)
{
	return GMRFLib_gsl_ensure_spd_core(A, tol, 0, msg, NULL);
}

int GMRFLib_gsl_ensure_spd_x(gsl_matrix *A, double tol, char **msg, GMRFLib_gsl_ensure_spd_store_tp *store)
{
	return GMRFLib_gsl_ensure_spd_core(A, tol, 0, msg, store);
}

int GMRFLib_gsl_ensure_spd_inverse(gsl_matrix *A, double tol, char **msg)
{
	return GMRFLib_gsl_ensure_spd_core(A, tol, 1, msg, NULL);
}

int GMRFLib_gsl_ensure_spd_inverse_x(gsl_matrix *A, double tol, char **msg, GMRFLib_gsl_ensure_spd_store_tp *store)
{
	return GMRFLib_gsl_ensure_spd_core(A, tol, 1, msg, store);
}

GMRFLib_gsl_ensure_spd_store_tp *GMRFLib_gsl_ensure_spd_store_alloc(int n)
{
	GMRFLib_gsl_ensure_spd_store_tp *S = Calloc(1, GMRFLib_gsl_ensure_spd_store_tp);
	S->U = gsl_matrix_alloc(n, n);
	S->M1 = gsl_matrix_alloc(n, n);
	S->M2 = gsl_matrix_alloc(n, n);
	S->AA = gsl_matrix_alloc(n, n);
	S->S = gsl_vector_alloc(n);
	S->work = gsl_eigen_symmv_alloc(n);
	return S;
}

int GMRFLib_gsl_ensure_spd_store_free(GMRFLib_gsl_ensure_spd_store_tp *S)
{
	if (S) {
		gsl_matrix_free(S->U);
		gsl_matrix_free(S->M1);
		gsl_matrix_free(S->M2);
		gsl_matrix_free(S->AA);
		gsl_vector_free(S->S);
		gsl_eigen_symmv_free(S->work);
	}
	return 0;
}

int GMRFLib_gsl_ensure_spd_core(gsl_matrix *A, double tol, int method, char **msg, GMRFLib_gsl_ensure_spd_store_tp *store)
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
	gsl_matrix *U = NULL;
	gsl_vector *S = NULL;
	gsl_eigen_symmv_workspace *work = NULL;

	if (store) {
		store->U->size1 = A->size1;
		store->U->size2 = A->size2;
		gsl_matrix_memcpy(store->U, A);
		U = store->U;
		store->S->size = A->size1;
		S = store->S;
		store->work->size = A->size1;
		work = store->work;
	} else {
		U = GMRFLib_gsl_duplicate_matrix(A);
		S = gsl_vector_alloc(A->size1);
		work = gsl_eigen_symmv_alloc(A->size1);
	}

	gsl_eigen_symmv(A, S, U, work);

	size_t i;
	const double one = 1.0, zero = 0.0;
	double s, s_min, s_max = gsl_vector_max(S);
	if (s_max <= 0.0) {
		s_max = 0.0;				       /* then the whole matrix is zero or INF, as all is wrong... */
	}

	if (tol < 0.0 && s_max > 0.0) {
		int n_neg = 0;
		s_min = s_max;
		for (i = 0; i < A->size1; i++) {
			s = gsl_vector_get(S, i);
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

	gsl_matrix *M1 = NULL;
	gsl_matrix *M2 = NULL;

	if (store) {
		store->M1->size1 = A->size1;
		store->M1->size2 = A->size2;
		M1 = store->M1;

		store->M2->size1 = A->size1;
		store->M2->size2 = A->size2;
		M2 = store->M2;
	} else {
		M1 = gsl_matrix_alloc(A->size1, A->size2);
		M2 = gsl_matrix_alloc(A->size1, A->size2);
	}

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

	if (!store) {
		gsl_matrix_free(U);
		gsl_matrix_free(M1);
		gsl_matrix_free(M2);
		gsl_vector_free(S);
		gsl_eigen_symmv_free(work);
	}

	return GMRFLib_SUCCESS;
}

GMRFLib_gsl_spd_solve_store_tp *GMRFLib_gsl_spd_solve_store_alloc(int n)
{
	GMRFLib_gsl_spd_solve_store_tp *S = Calloc(1, GMRFLib_gsl_spd_solve_store_tp);
	S->L = gsl_matrix_alloc(n, n);
	S->S = gsl_vector_alloc(n);
	return S;
}

int GMRFLib_gsl_spd_solve_store_free(GMRFLib_gsl_spd_solve_store_tp *store)
{
	if (store) {
		gsl_matrix_free(store->L);
		gsl_vector_free(store->S);
		Free(store);
	}
	return 0;
}

int GMRFLib_gsl_spd_solve(gsl_matrix *A, gsl_vector *b, gsl_vector *x)
{
	return GMRFLib_gsl_spd_solve_x(A, b, x, NULL);
}

int GMRFLib_gsl_spd_solve_x(gsl_matrix *A, gsl_vector *b, gsl_vector *x, GMRFLib_gsl_spd_solve_store_tp *store)
{
	gsl_matrix *L = NULL;
	gsl_vector *S = NULL;

	if (store) {
		assert(store->L);
		store->L->size1 = A->size1;
		store->L->size2 = A->size2;
		L = store->L;
		gsl_matrix_memcpy(L, A);
		store->S->size = A->size1;
		S = store->S;
	} else {
		L = GMRFLib_gsl_duplicate_matrix(A);
		S = gsl_vector_alloc(A->size1);
	}

	int ecode = gsl_linalg_cholesky_decomp2(L, S);
	if (ecode != GSL_SUCCESS) {
		if (!store) {
			gsl_matrix_free(L);
			gsl_vector_free(S);
		}
		return !GMRFLib_SUCCESS;
	}

	gsl_linalg_cholesky_solve2(L, S, b, x);
	if (!store) {
		gsl_matrix_free(L);
		gsl_vector_free(S);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_gsl_safe_spd_solve(gsl_matrix *A, gsl_vector *b, gsl_vector *x, double tol, int *try_first)
{
	/*
	 * solve Ax=b, ignoring contributions from eigenvalues < tol*max(eigenval)
	 *
	 * solution is returned in x, while A and b is not changed.
	 *
	 * try first a regular one, if fail, the try the svd one
	 */

	assert(A && (A->size1 == A->size2));
	if (try_first && *try_first) {
		int ecode = GMRFLib_gsl_spd_solve_x(A, b, x, NULL);
		if (ecode == GMRFLib_SUCCESS) {
			return GMRFLib_SUCCESS;
		}
	}
	// if try_first did not succeed, then set it to 0
	if (try_first)
		*try_first = 0;

	const int debug = 0;
	assert(tol >= 0.0);

	gsl_matrix *U = GMRFLib_gsl_duplicate_matrix(A);
	gsl_vector *S = gsl_vector_alloc(A->size1);

	gsl_eigen_symmv_workspace *work = gsl_eigen_symmv_alloc(A->size1);
	gsl_eigen_symmv(A, S, U, work);

	size_t i;
	const double one = 1.0, zero = 0.0;
	double s;
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

int GMRFLib_gsl_spd_inv(gsl_matrix *A, double tol, int *try_first)
{
	/*
	 * A=inv(A) for symmetric A, ignoring contributions from eigenvalues < tol*max(eigenval)
	 *
	 * if FORCE, then do not try spd_inverse first
	 */

	assert(A && (A->size1 == A->size2));
	if (try_first && *try_first) {
		int ecode = GMRFLib_gsl_spd_inverse(A);
		if (ecode == GMRFLib_SUCCESS) {
			return GMRFLib_SUCCESS;
		}
	}
	// if try_first did not succeed, then set it to 0
	if (try_first)
		*try_first = 0;

	const int debug = 0;
	assert(tol >= 0.0);

	gsl_matrix *U = GMRFLib_gsl_duplicate_matrix(A);
	gsl_vector *S = gsl_vector_alloc(A->size1);

	gsl_eigen_symmv_workspace *work = gsl_eigen_symmv_alloc(A->size1);
	gsl_eigen_symmv(A, S, U, work);

	size_t i;
	const double one = 1.0, zero = 0.0;
	double s;
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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((optimize("O3")))
    __attribute__((target_clones(INLA_CLONE_TARGETS "default")))
int GMRFLib_gsl_mgs(gsl_matrix *A)
{
	// this is the modified Gram-Schmith ortogonalisation, and it 
	// overwrite matrix A with its ortogonal normalized basis, column by column.
	// the sign of each column is so that max(abs(column)) is positive

	gsl_vector *q = gsl_vector_alloc(A->size1);
	size_t n1 = A->size1;
	size_t n2 = A->size2;

	for (size_t i = 0; i < n2; i++) {
		double aij_amax = 0.0, r = 0.0;

		for (size_t j = 0; j < n1; j++) {
			double elm = gsl_matrix_get(A, j, i);
			r += SQR(elm);
			aij_amax = (ABS(elm) > ABS(aij_amax) ? elm : aij_amax);
		}

		if (aij_amax < 0.0) {			       /* swap the sign of this column */
			for (size_t j = 0; j < n1; j++) {
				gsl_matrix_set(A, j, i, -gsl_matrix_get(A, j, i));
			}
		}

		r = sqrt(r);
		for (size_t j = 0; j < n1; j++) {
			gsl_vector_set(q, j, gsl_matrix_get(A, j, i) / r);
			gsl_matrix_set(A, j, i, gsl_vector_get(q, j));	// normalize
		}

		for (size_t j = i + 1; j < n2; j++) {
			r = 0.0;
			for (size_t k = 0; k < n1; k++) {
				r += gsl_vector_get(q, k) * gsl_matrix_get(A, k, j);
			}
			for (size_t k = 0; k < n1; k++) {
				gsl_matrix_set(A, k, j, gsl_matrix_get(A, k, j) - r * gsl_vector_get(q, k));
			}
		}
	}
	gsl_vector_free(q);

	return (GMRFLib_SUCCESS);
}
#pragma GCC diagnostic pop

GMRFLib_gsl_low_rank_store_tp *GMRFLib_gsl_low_rank_store_alloc(int n)
{
	GMRFLib_gsl_low_rank_store_tp *S = Calloc(1, GMRFLib_gsl_low_rank_store_tp);
	S->U = gsl_matrix_alloc(n, n);
	S->D = gsl_matrix_alloc(n, n);
	S->S = gsl_vector_alloc(n);
	S->work = gsl_eigen_symmv_alloc(n);

	return S;
}

int GMRFLib_gsl_low_rank_store_free(GMRFLib_gsl_low_rank_store_tp *S)
{
	if (S) {
		gsl_matrix_free(S->U);
		gsl_matrix_free(S->D);
		gsl_vector_free(S->S);
		gsl_eigen_symmv_free(S->work);
		Free(S);
	}
	return 0;
}

gsl_matrix *GMRFLib_gsl_low_rank(gsl_matrix *Cov, double tol, gsl_matrix *B)
{
	return GMRFLib_gsl_low_rank_x(Cov, tol, B, NULL);
}

gsl_matrix *GMRFLib_gsl_low_rank_x(gsl_matrix *Cov, double tol, gsl_matrix *B, GMRFLib_gsl_low_rank_store_tp *store)
{
	// if STORE, then assume elements are large enough

	// if B!=NULL then overwrite B (assume its enough space), otherwise return a new B

	// Compute the low-rank representation in terms of x=Bz, where Cov(z)=I, from a given possible singular covariance matrix. We ignore
	// contributions from eigenvalue < tol*max(eigenval)

	assert(Cov && (Cov->size1 == Cov->size2));
	assert(tol >= 0.0 && tol <= 1.0);

	size_t n = Cov->size1;
	gsl_matrix *U = NULL;
	gsl_matrix *D = NULL;
	gsl_vector *S = NULL;
	gsl_eigen_symmv_workspace *work = NULL;

	if (store) {
		store->U->size1 = n;
		store->U->size2 = n;
		gsl_matrix_memcpy(store->U, Cov);

		store->D->size1 = n;
		store->D->size2 = n;
		store->S->size = n;

		D = store->D;
		S = store->S;
		U = store->U;
		store->work->size = n;
		work = store->work;
	} else {
		U = GMRFLib_gsl_duplicate_matrix(Cov);
		D = gsl_matrix_alloc(n, n);
		S = gsl_vector_alloc(n);
		work = gsl_eigen_symmv_alloc(n);
	}

	gsl_eigen_symmv(Cov, S, U, work);
	gsl_eigen_symmv_sort(S, U, GSL_EIGEN_SORT_VAL_DESC);

	double s_max = gsl_vector_max(S);
	assert(s_max > 0.0);
	assert(s_max == gsl_vector_get(S, 0));

	size_t m = 0;
	double s_min = tol * s_max;
	for (size_t i = 0; i < n; i++) {
		if (gsl_vector_get(S, i) >= s_min) {
			m++;
		} else {
			break;
		}
	}

	D->size2 = m;
	gsl_matrix_set_zero(D);
	for (size_t i = 0; i < m; i++) {
		gsl_matrix_set(D, i, i, sqrt(gsl_vector_get(S, i)));
	}

	gsl_matrix *BB = NULL;
	if (B) {
		B->size1 = n;
		B->size2 = m;
		BB = B;
	} else {
		BB = gsl_matrix_alloc(n, m);
	}
	GMRFLib_gsl_mm(U, D, BB);

	if (!store) {
		gsl_matrix_free(U);
		gsl_matrix_free(D);
		gsl_vector_free(S);
		gsl_eigen_symmv_free(work);
	}

	return (B ? NULL : BB);
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((optimize("O3")))
    __attribute__((target_clones(INLA_CLONE_TARGETS "default")))
double GMRFLib_gsl_kld(gsl_vector *m_base, gsl_matrix *Q_base, gsl_vector *m, gsl_matrix *Q, double tol, int *rankdef)
{
	// compute the KLD between two mult-var normals, where either or both matrices can be numerical singular. Make sure that the rank is the
	// same for both, and return the computed rank-deficiency in 'rankdef' if given

	FIXME("gsl_kld is not tested");

	size_t n = Q_base->size1;
	assert(n == Q_base->size2);
	assert(n == Q->size2);
	assert(n == Q->size2);

	gsl_matrix *U_base = GMRFLib_gsl_duplicate_matrix(Q_base);
	gsl_matrix *U = GMRFLib_gsl_duplicate_matrix(Q);
	gsl_vector *S_base = gsl_vector_alloc(n);
	gsl_vector *S = gsl_vector_alloc(n);
	gsl_eigen_symmv_workspace *work_base = gsl_eigen_symmv_alloc(n);
	gsl_eigen_symmv_workspace *work = gsl_eigen_symmv_alloc(n);
	double tolerance_base = tol;
	double tolerance = tol;
	const double one = 1.0, zero = 0.0;

	gsl_vector_set_zero(S_base);
	gsl_vector_set_zero(S);

	gsl_eigen_symmv(Q_base, S_base, U_base, work_base);
	gsl_eigen_symmv(Q, S, U, work);

	size_t i;
	double s = 0.0, s_min_base = 0.0, s_max_base = gsl_vector_max(S_base);
	s_max_base = DMAX(0.0, s_max_base);

	if (s_max_base > 0.0) {
		s_min_base = s_max_base;
		for (i = 0; i < n; i++) {
			s = gsl_vector_get(S_base, i);
			if (s > 0.0 && s < s_min_base) {
				s_min_base = s;
			}
		}
		tolerance_base = s_min_base / s_max_base;
	}

	double s_min = 0.0, s_max = gsl_vector_max(S);
	s_max = DMAX(0.0, s_max);

	if (s_max > 0.0) {
		s_min = s_max;
		for (i = 0; i < n; i++) {
			s = gsl_vector_get(S, i);
			if (s > 0.0 && s < s_min) {
				s_min = s;
			}
		}
		tolerance = s_min / s_max;
	}

	int rdef_base = 0;
	int rdef = 0;

	s_min_base = tolerance_base * s_max_base;
	s_min = tolerance * s_max;

	for (i = 0; i < n; i++) {
		s = gsl_vector_get(S_base, i);
		rdef_base += (s < s_min_base);
		if (s < s_min_base)
			s = 0.0;
		gsl_vector_set(S_base, i, s);

		s = gsl_vector_get(S, i);
		rdef += (s < s_min);
		if (s < s_min)
			s = 0.0;
		gsl_vector_set(S, i, s);
	}

	rdef = IMAX(rdef_base, rdef);

	double ldet_base = 0.0;
	double ldet = 0.0;
	for (i = 0; i < n - rdef; i++) {
		ldet_base += log(gsl_vector_get(S_base, i));
		ldet += log(gsl_vector_get(S, i));
	}

	gsl_matrix *Cov_base = GMRFLib_gsl_duplicate_matrix(Q_base);
	gsl_matrix *M1 = gsl_matrix_alloc(n, n);
	gsl_matrix *M2 = gsl_matrix_alloc(n, n);
	gsl_matrix_set_zero(M1);
	gsl_matrix_set_zero(M2);

	for (i = 0; i < n - rdef; i++) {
		gsl_matrix_set(M2, i, i, 1.0 / gsl_vector_get(S_base, i));
	}
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, one, M2, U_base, zero, M1);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, one, U_base, M1, zero, M2);
	gsl_matrix_memcpy(Cov_base, M2);

	gsl_matrix *Q_corr = GMRFLib_gsl_duplicate_matrix(Q);
	gsl_matrix_set_zero(M1);
	gsl_matrix_set_zero(M2);
	for (i = 0; i < n - rdef; i++) {
		gsl_matrix_set(M2, i, i, gsl_vector_get(S, i));
	}
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, one, M2, U, zero, M1);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, one, U, M1, zero, M2);
	gsl_matrix_memcpy(Q_corr, M2);

	double kld = 0.0;

	// trace-term
	gsl_matrix_set_zero(M1);
	gsl_matrix_set_zero(M2);

	GMRFLib_gsl_mm(Cov_base, Q_corr, M1);
	for (i = 0; i < n; i++) {
		kld += gsl_matrix_get(M1, i, i);
	}
	kld *= (double) (n - rdef) / (double) n;

	// quadratic term
	gsl_vector *v1 = gsl_vector_alloc(n);
	gsl_vector *v2 = gsl_vector_alloc(n);
	gsl_vector_set_zero(v1);
	gsl_vector_set_zero(v2);

	for (i = 0; i < n; i++) {
		gsl_vector_set(v1, i, gsl_vector_get(m_base, i) - gsl_vector_get(m, i));
	}
	GMRFLib_gsl_mv(Cov_base, v1, v2);
	double quadratic = 0.0;
	gsl_blas_ddot(v1, v2, &quadratic);
	kld += quadratic;

	kld = 0.5 * (kld - (double) (n - rdef) + (ldet_base - ldet));
	if (rankdef) {
		*rankdef = rdef;
	}

	gsl_eigen_symmv_free(work_base);
	gsl_eigen_symmv_free(work);
	gsl_matrix_free(M1);
	gsl_matrix_free(M2);
	gsl_matrix_free(Q_corr);
	gsl_matrix_free(U);
	gsl_matrix_free(U_base);
	gsl_vector_free(S);
	gsl_vector_free(S_base);
	gsl_vector_free(v1);
	gsl_vector_free(v2);

	return kld;
}
#pragma GCC diagnostic pop

double GMRFLib_dssqr(int n, double *x)
{
	// sum_i x_i^2
	int one = 1;
	double ret = dnrm2_(&n, x, &one);
	return SQR(ret);
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((optimize("O3")))
    __attribute__((target_clones(INLA_CLONE_TARGETS "default")))
void GMRFLib_dscale(int n, double a, double *x)
{
	// x[i] *= a
	if (n < 32) {
#pragma omp simd
		for (int i = 0; i < n; i++) {
			x[i] *= a;
		}
	} else {
		int one = 1;
		dscal_(&n, &a, x, &one);
	}
}
#pragma GCC diagnostic pop

#define DSCALE2_CORE()				\
	_Pragma("omp simd")			\
	for (int i = 0; i < n; i++) {		\
		y[i] = a * x[i];		\
	}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((optimize("O3")))
    __attribute__((target_clones(INLA_CLONE_TARGETS "default")))
void GMRFLib_dscale2(int n, double a, double *__restrict x, double *__restrict y)
{
	// y[i] = a * x[i]
#if defined(INLA_WITH_SIMDE_AVX512F_) && defined(__AVX512F__)
#       include "intrinsics/simde/dscale2-avx512f.h"
#elif defined(INLA_WITH_SIMDE_AVX2_) && (!defined(__x86_64__) || (defined(__x86_64__) && defined(__AVX2__)))
#       include "intrinsics/simde/dscale2-avx2.h"
#elif defined(INLA_WITH_SIMDE)
#       include "intrinsics/simde/dscale2-sse2.h"
#else
	DSCALE2_CORE();
#endif
}
#pragma GCC diagnostic pop
#undef DSCALE2_CORE

void GMRFLib_daxpby(int n, double a, double *x, double b, double *y)
{
	// y = a * x + b * y
#if defined(INLA_WITH_MKL)
	int inc = 1;
	daxpby_(&n, &a, x, &inc, &b, y, &inc);
#elif defined(INLA_WITH_ARMPL)
	int inc = 1;
	daxpby_(&n, &a, x, &inc, &b, y, &inc);
#else
	GMRFLib_dscale(n, b, y);
	GMRFLib_daxpy(n, a, x, y);
#endif
}

void GMRFLib_daxpbyz(int n, double a, double *x, double b, double *y, double *z)
{
	// z = a * x + b * y
#if defined(INLA_WITH_ARMPL)
	int inc = 1;
	dwaxpby_(&n, &a, x, &inc, &b, y, &inc, z, &inc);
#else
	Memcpy(z, y, n * sizeof(double));
	GMRFLib_daxpby(n, a, x, b, z);
#endif
}

void GMRFLib_daxpbypcz(int n, double a, double *x, double b, double *y, double c, double *z)
{
	// z[i] = a * x[i] + b * y[i] + c;
	GMRFLib_daxpb(n, b, y, c, z);
	GMRFLib_daxpy(n, a, x, z);
}

void GMRFLib_daxpb(int n, double a, double *x, double b, double *y)
{
	// y[i] = a * x[i] + b
	GMRFLib_dfill(n, b, y);
	GMRFLib_daxpy(n, a, x, y);
}

// y = a * x + y
#define DAXPY_CORE(cutoff_)				\
	if (n < cutoff_) {				\
		int limit = n & ~3;			\
		_Pragma("omp simd")			\
		for (int i = 0; i < limit; i += 4) {	\
			y[i] += a * x[i];		\
			y[i+1] += a * x[i+1];		\
			y[i+2] += a * x[i+2];		\
			y[i+3] += a * x[i+3];		\
		}					\
		for (int i = limit; i < n; i++) {	\
			y[i] += a * x[i];		\
		}					\
	} else {					\
		int inc = 1;				\
		daxpy_(&n, &a, x, &inc, y, &inc);	\
	}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((optimize("O3")))
    __attribute__((target_clones(INLA_CLONE_TARGETS "default")))
void GMRFLib_daxpy(int n, double a, double *x, double *y)
{
//      int inc = 1;
//      daxpy_(&n, &a, x, &inc, y, &inc);
	DAXPY_CORE(64);
}
#pragma GCC diagnostic pop

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((optimize("O3")))
    __attribute__((target_clones(INLA_CLONE_TARGETS "default")))
void GMRFLib_daxpy_x(int n, double a, double *x, double *y, int cutoff)
{
	DAXPY_CORE(cutoff);
}
#pragma GCC diagnostic pop
#undef DAXPY_CORE

#define DDOT_CORE(cutoff_)						\
	if (n < cutoff_) {						\
		double sum0 = 0.0, sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;	\
		int limit = n & ~3;					\
		_Pragma("omp simd reduction(+:sum0,sum1,sum2,sum3)")	\
		for (int i = 0; i < limit; i += 4) {			\
			sum0 += x[i] * y[i];				\
			sum1 += x[i+1] * y[i+1];			\
			sum2 += x[i+2] * y[i+2];			\
			sum3 += x[i+3] * y[i+3];			\
		}							\
		for (int i = limit; i < n; i++) {			\
			sum0 += x[i] * y[i];				\
		}							\
		return sum0 + sum1 + sum2 + sum3;			\
	} else {							\
		int one = 1;						\
		return ddot_(&n, x, &one, y, &one);			\
	}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((optimize("O3")))
    __attribute__((target_clones(INLA_CLONE_TARGETS "default")))
double GMRFLib_ddot(int n, double *__restrict x, double *__restrict y)
{
//      int one = 1;
//      return ddot_(&n, x, &one, y, &one);
	DDOT_CORE(32);
}
#pragma GCC diagnostic pop

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((optimize("O3")))
    __attribute__((target_clones(INLA_CLONE_TARGETS "default")))
void GMRFLib_ddot2(double *a, double *b, int n, double *__restrict x, double *__restrict y, double *__restrict z)
{
// simde_mm512_reduce_add_pd is not yet supported so we disable it for the moment (not properly tested!)
//#if defined(INLA_WITH_SIMDE_AVX512F_) && defined(__AVX512F__)
//#       include "intrinsics/simde/ddot2-avx512f.h"
#if defined(INLA_WITH_SIMDE_AVX2_) && (!defined(__x86_64__) || (defined(__x86_64__) && defined(__AVX2__)))
#       include "intrinsics/simde/ddot2-avx2.h"
#elif defined(INLA_WITH_SIMDE)
#       include "intrinsics/simde/ddot2-sse2.h"
#else
	double aa = 0.0, bb = 0.0;
#       pragma omp simd reduction(+: aa, bb)
	for (int i = 0; i < n; i++) {
		aa += x[i] * y[i];
		bb += x[i] * z[i];
	}
	*a = aa;
	*b = bb;
#endif
}
#pragma GCC diagnostic pop

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((optimize("O3")))
    __attribute__((target_clones(INLA_CLONE_TARGETS "default")))
double GMRFLib_ddot_x(int n, double *__restrict x, double *__restrict y, int cutoff)
{
	DDOT_CORE(cutoff);
}
#pragma GCC diagnostic pop
#undef DDOT_CORE

#define SUM_CORE(TYPE_)						\
	TYPE_ r = 0;						\
	_Pragma("omp simd reduction(+: r)")			\
	for (int i = 0; i < n; i++) {				\
		r += x[i];					\
	}							\
	return r + r0

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((optimize("O3")))
    __attribute__((target_clones(INLA_CLONE_TARGETS "default")))
double GMRFLib_dsum(int n, double *x)
{
	double r0 = 0.0;
// simde_mm512_reduce_add_pd is not yet supported so we disable it for the moment (not properly tested!)
//#if defined(INLA_WITH_SIMDE_AVX512F_) && defined(__AVX512F__)
//#       include "intrinsics/simde/dsum-avx512f.h"
#if defined(INLA_WITH_SIMDE_AVX2_) && (!defined(__x86_64__) || (defined(__x86_64__) && defined(__AVX2__)))
#       include "intrinsics/simde/dsum-avx2.h"
#elif defined(INLA_WITH_SIMDE)
#       include "intrinsics/simde/dsum-sse2.h"
#else
	SUM_CORE(double);
#endif
}
#pragma GCC diagnostic pop

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((optimize("O3")))
    __attribute__((target_clones(INLA_CLONE_TARGETS "default")))
int GMRFLib_isum(int n, int *x)
{
	int r0 = 0;
#if defined(INLA_WITH_SIMDE)
#       include "intrinsics/simde/isum-sse2.h"
#else
	SUM_CORE(int);
#endif
}
#undef SUM_CORE
#pragma GCC diagnostic pop

#define SPARSE_DSUM()					\
	double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;	\
	int m = n & ~7;					\
	for (int i = 0; i < m; i += 8) {		\
		s0 += a[idx[i + 0]];			\
		s1 += a[idx[i + 1]];			\
		s2 += a[idx[i + 2]];			\
		s3 += a[idx[i + 3]];			\
		s0 += a[idx[i + 4]];			\
		s1 += a[idx[i + 5]];			\
		s2 += a[idx[i + 6]];			\
		s3 += a[idx[i + 7]];			\
	}						\
	for (int i = m; i < n; i++) {			\
		s0 += a[idx[i]];			\
	}						\
	return s0 + s1 + s2 + s3

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((optimize("O3")))
    __attribute__((target_clones(INLA_CLONE_TARGETS "default")))
double GMRFLib_sparse_dsum(int n, double *__restrict a, int *__restrict idx)
{
	SPARSE_DSUM();
}
#pragma GCC diagnostic pop
#undef SPARSE_DSUM

#define FILL_CORE(TYPE_, LEN_)						\
	if (ISZERO(a)) {						\
		memset((void *) x, 0, (size_t) (n * sizeof(TYPE_)));	\
	} else {							\
		int nn = IMIN(LEN_, n);					\
		_Pragma("omp simd")					\
		for (int i = 0; i < nn; i++) {   			\
			x[i] = a;					\
		}							\
		if (n > LEN_) {						\
			int offset = LEN_;				\
			while (offset < n) {				\
				int len = IMIN(offset, n - offset);	\
				memcpy((void *) (x + offset), (void *) x, (size_t) (len * sizeof(TYPE_))); \
				offset += len;				\
			}						\
		}							\
	}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((optimize("O3")))
    __attribute__((target_clones(INLA_CLONE_TARGETS "default")))
void GMRFLib_dfill(int n, double a, double *x)
{
	FILL_CORE(double, 64);
}
#pragma GCC diagnostic pop

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((optimize("O3")))
    __attribute__((target_clones(INLA_CLONE_TARGETS "default")))
void GMRFLib_ifill(int n, int a, int *x)
{
	FILL_CORE(int, 128);
}
#pragma GCC diagnostic pop

void GMRFLib_bfill(int n, bool a, bool *x)
{
	FILL_CORE(bool, 512);
}
#undef FILL_CORE

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((optimize("O3")))
    __attribute__((target_clones(INLA_CLONE_TARGETS "default")))
void GMRFLib_pack(int n, double *a, int *ia, double *y)
{
	// y[] = a[ia[]]
#if defined(INLA_WITH_MKL)
	vdPackV(n, a, ia, y);
#elif defined(INLA_WITH_SIMDE_AVX512F_) && defined(__AVX512F__)
#       include "intrinsics/simde/pack-avx512f.h"
#elif defined(INLA_WITH_SIMDE_AVX2_) && (!defined(__x86_64__) || (defined(__x86_64__) && defined(__AVX2__)))
#       include "intrinsics/simde/pack-avx2.h"
#elif defined(INLA_WITH_SIMDE)
#       include "intrinsics/simde/pack-sse2.h"
#else
#       pragma omp simd
	for (int i = 0; i < n; i++) {
		y[i] = a[ia[i]];
	}
#endif
}
#pragma GCC diagnostic pop

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((optimize("O3")))
    __attribute__((target_clones(INLA_CLONE_TARGETS "default")))
void GMRFLib_unpack(int n, double *a, double *y, int *iy)
{
	// y[iy[]] = a[]
#if defined(INLA_WITH_MKL)
	vdUnpackV(n, a, y, iy);
#else
#       pragma omp simd
	for (int i = 0; i < n; i++) {
		y[iy[i]] = a[i];
	}
#endif
}
#pragma GCC diagnostic pop

#pragma GCC diagnostic push
__attribute__((optimize("O3")))
void GMRFLib_powx(int n, double *x, double a, double *y)
{
	// y = x^a
#if defined(INLA_WITH_MKL)
	vdPowx(n, x, a, y);
#else
#       pragma omp simd
	for (int i = 0; i < n; i++) {
		y[i] = pow(x[i], a);
	}
#endif
}
#pragma GCC diagnostic pop
