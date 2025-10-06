#include <stdio.h>
#include <stdlib.h>
#include "inla.h"
#include "GMRFLib/GMRFLib.h"

// NOTE: not very efficient but the use is for small matrices only, so its fine

#define SQRT12 3.464101615137755

inla_bm_tp *inla_bm_alloc(int n, int lbw, int ubw)
{
	if (n <= 0 || lbw < 0 || ubw < 0)
		return NULL;
	assert(ubw < n);
	assert(lbw < n);

	inla_bm_tp *A = Calloc(1, inla_bm_tp);
	A->n = n;
	A->ubw = ubw;
	A->lbw = lbw;
	A->ldim = lbw + ubw + 1;
	A->memlen = A->n * A->ldim;
	A->band = Calloc(A->memlen, double);

	return A;
}

void inla_bm_free(inla_bm_tp *A)
{
	if (A) {
		Free(A->band);
		Free(A);
	}
}

void inla_bm_fprintf(FILE *fp, inla_bm_tp *A, char *msg)
{
	fp = (fp ? fp : stdout);

	fprintf(fp, "\n[%s] n= %1d lbw= %1d ubw= %1d ldim= %1d memlen= %1d\n", (msg ? msg : ""), A->n, A->lbw, A->ubw, A->ldim, A->memlen);
	for (int i = 0; i < A->n; i++) {
		for (int j = IMAX(0, i - A->lbw); j <= IMIN(A->n - 1, i + A->ubw); j++) {
			if (BM_VAL(A, i, j)) {
				fprintf(fp, "i %1d j %1d value %.10f\n", i, j, BM_VAL(A, i, j));
			}
		}
	}
}

inla_bm_tp *inla_bm_chol(inla_bm_tp *A, inla_bm_tp *chol)
{
	// return chol (=L) of A. if 'L' overwrite L, otherwise alloc L,

	assert(A->ubw == 0);
	chol = (chol ? chol : inla_bm_alloc(A->n, A->lbw, 0));
	Memcpy(chol->band, A->band, A->memlen * sizeof(double));

	int info = 0;
	dpbtrf_("L", &(chol->n), &(chol->lbw), chol->band, &(chol->ldim), &info, F_ONE);
	assert(info == 0);

	return chol;
}

void inla_bm_solve(inla_bm_tp *chol, double *b, double *x)
{
	// solve Q x = B. if B is NULL then assume b=x. return the solution in x
	assert(chol->ubw == 0);
	int stride = 1;
	if (b) {
		Memcpy(x, b, chol->n * sizeof(double));
	}
	dtbsv_("L", "N", "N", &(chol->n), &(chol->lbw), chol->band, &(chol->ldim), x, &stride, F_ONE, F_ONE, F_ONE);
	dtbsv_("L", "T", "N", &(chol->n), &(chol->lbw), chol->band, &(chol->ldim), x, &stride, F_ONE, F_ONE, F_ONE);
}

void inla_bm_nsolve(inla_bm_tp *chol, double *b, double *x, int nrhs)
{
	for (int i = 0; i < nrhs; i++) {
		int n = chol->n;
		inla_bm_solve(chol, (b ? b + i * n : NULL), x + i * n);
	}
}

inla_bm_tp *inla_bm_partial_inv(inla_bm_tp *chol, inla_bm_tp *inv)
{
	// note that this create a lower-triangular matrix... so better accessed with BM_VAL(inv,IMAX(i,j),IMIN(i,j))

	if (inv) {
		assert(inv->n == chol->n);
		assert(inv->lbw == chol->lbw);
		assert(inv->ubw == chol->ubw);
		GMRFLib_dfill(inv->memlen, 0, inv->band);
	} else {
		inv = inla_bm_alloc(chol->n, chol->lbw, 0);
	}

	for (int i = chol->n - 1; i >= 0; i--) {
		double iLii = 1.0 / BM_VAL(chol, i, i);
		for (int j = IMIN(chol->n - 1, i + chol->lbw); j >= i; j--) {
			double sum = 0.0;
			for (int k = i + 1; k <= IMIN(chol->n - 1, i + chol->lbw); k++) {
				sum += sBM_VAL(chol, k, i) * sBM_VAL(inv, k, j);
			}
			BM_VAL(inv, j, i) = (i == j ? SQR(iLii) : 0.0) - sum * iLii;
		}
	}
	return inv;
}

inla_bm_tp *inla_bm_mm(inla_bm_tp *A, inla_bm_tp *B, inla_bm_tp *AB)
{
	// AB = A %*% B

	assert(A->n == B->n);
	int lbw = IMIN(A->n - 1, A->lbw + B->lbw);
	int ubw = IMIN(A->n - 1, A->ubw + B->ubw);
	if (AB) {
		assert(AB->lbw == lbw && AB->ubw == ubw);
		GMRFLib_dfill(AB->memlen, 0.0, AB->band);
	} else {
		AB = inla_bm_alloc(A->n, lbw, ubw);
	}

	// AB_ij = sum_k A_ik * B_kj
	for (int i = 0; i < AB->n; i++) {
		for (int j = IMAX(0, i - AB->lbw); j <= IMIN(AB->n - 1, i + AB->ubw); j++) {
			double sum = 0.0;
			for (int k = IMAX(0, i - A->lbw); k <= IMIN(A->n - 1, i + A->ubw); k++) {
				if (j - k >= -B->lbw && j - k <= B->ubw) {
					sum += BM_VAL(A, i, k) * BM_VAL(B, k, j);
				}
			}
			BM_VAL(AB, i, j) = sum;
		}
	}
	return AB;
}

inla_bm_tp *inla_bm_mmm(inla_bm_tp *A, inla_bm_tp *B, inla_bm_tp *C, inla_bm_tp *ABC)
{
	// ABC = A %*% B %*% C
	assert(A->n == B->n && A->n == C->n);
	int lbw = IMIN(A->n - 1, A->lbw + B->lbw + C->lbw);
	int ubw = IMIN(A->n - 1, A->ubw + B->ubw + C->ubw);

	if (ABC) {
		assert(ABC->lbw == lbw && ABC->ubw == ubw);
		GMRFLib_dfill(ABC->memlen, 0.0, ABC->band);
	} else {
		ABC = inla_bm_alloc(A->n, lbw, ubw);
	}

	inla_bm_tp *AB = inla_bm_mm(A, B, NULL);
	inla_bm_mm(AB, C, ABC);
	inla_bm_free(AB);

	return ABC;
}

inla_bm_tp *inla_bm_duplicate(inla_bm_tp *A)
{
	inla_bm_tp *Adup = inla_bm_alloc(A->n, A->lbw, A->ubw);
	memcpy(Adup->band, A->band, A->memlen * sizeof(double));
	return Adup;
}

inla_bm_tp *inla_bm_sum(double a, inla_bm_tp *A, double b, inla_bm_tp *B, double c, inla_bm_tp *C, inla_bm_tp *ABC)
{
	// ABC = A->a + b * B + c * C, C is allowed to be NULL

	assert(A->n == B->n);
	if (C) {
		assert(A->n == C->n);
	}

	int lbw = (C ? IMIN(A->n - 1, IMAX(A->lbw, IMAX(B->lbw, C->lbw))) : IMIN(A->n - 1, IMAX(A->lbw, B->lbw)));
	int ubw = (C ? IMIN(A->n - 1, IMAX(A->ubw, IMAX(B->ubw, C->ubw))) : IMIN(A->n - 1, IMAX(A->ubw, B->ubw)));
	if (ABC) {
		ABC->lbw = lbw;
		ABC->ubw = ubw;
	} else {
		ABC = inla_bm_alloc(A->n, lbw, ubw);
	}

	// we could write this better, but...
	for (int i = 0; i < ABC->n; i++) {
		for (int j = IMAX(0, i - ABC->lbw); j <= IMIN(ABC->n - 1, i + ABC->ubw); j++) {
			int ld = i - j, ud = j - i;
			double *ABC_ptr = BM_PTR(ABC, i, j);
			*ABC_ptr = 0.0;

			if (ld <= A->lbw && ud <= A->ubw) {
				*ABC_ptr += a * BM_VAL(A, i, j);
			}
			if (ld <= B->lbw && ud <= B->ubw) {
				*ABC_ptr += b * BM_VAL(B, i, j);
			}
			if (C) {
				if (ld <= C->lbw && ud <= C->ubw) {
					*ABC_ptr += c * BM_VAL(C, i, j);
				}
			}
		}
	}

	return (ABC);
}

inla_bm_tp *inla_bm_sym(inla_bm_tp *A, inla_bm_tp *Asym)
{
	// note that this create a lower-triangular matrix... so better accessed with BM_VAL(inv,IMAX(i,j),IMIN(i,j))
	if (Asym) {
		assert(A->n == Asym->n);
		assert(Asym->ubw == 0);
		assert(A->lbw == Asym->lbw);
		GMRFLib_dfill(Asym->memlen, 0.0, Asym->band);
	} else {
		Asym = inla_bm_alloc(A->n, A->lbw, 0);
	}

	for (int i = 0; i < Asym->n; i++) {
		for (int j = IMAX(0, i - Asym->lbw); j <= i; j++) {
			BM_VAL(Asym, i, j) = BM_VAL(A, i, j);
		}
	}

	return Asym;
}

inla_bm_tp *inla_bm_trans(inla_bm_tp *A, inla_bm_tp *At)
{
	if (At) {
		assert(A->n == At->n);
		assert(A->lbw == At->ubw);
		assert(A->ubw == At->lbw);
		GMRFLib_dfill(At->memlen, 0.0, At->band);
	} else {
		At = inla_bm_alloc(A->n, A->ubw, A->lbw);
	}

	for (int i = 0; i < At->n; i++) {
		for (int j = IMAX(0, i - At->lbw); j <= IMIN(At->n - 1, i + At->ubw); j++) {
			BM_VAL(At, i, j) = BM_VAL(A, j, i);
		}
	}

	return At;
}

void inla_bm_scale(double a, inla_bm_tp *A)
{
	// scale matrix A with 'a'
	GMRFLib_dscale(A->memlen, a, A->band);
}

double inla_prw2_corfunc(double d, double kappa)
{
	double ad = kappa * ABS(d);
	return (1.0 + ad) * exp(-ad);
}

inla_prw2_arg_tp *inla_prw2_create(int n, double *loc)
{
	// create and return prw2_arg
	assert(GMRFLib_OPENMP_IN_SERIAL());

	inla_prw2_arg_tp *arg = Calloc(1, inla_prw2_arg_tp);
	arg->n = n;
	arg->loc = Malloc(n, double);
	Memcpy(arg->loc, loc, n * sizeof(double));

	double *h = Malloc(n - 1, double);
	for (int i = 0; i < n - 1; i++) {
		h[i] = loc[i + 1] - loc[i];
	}
	arg->h = h;

	inla_bm_tp *C = inla_bm_alloc(n, 1, 1);
	for (int i = 0; i < n; i++) {
		if (i == 0) {
			BM_VAL(C, i, i) = h[0] / 3.0;
		} else if (i == n - 1) {
			BM_VAL(C, i, i) = h[n - 2] / 3.0;
		} else {
			BM_VAL(C, i, i) = (h[i] + h[i - 1]) / 3.0;
		}

		int j = i + 1;
		if (j < n) {
			BM_VAL(C, i, j) = BM_VAL(C, j, i) = h[i] / 6.0;
		}
	}
	arg->C = C;

	inla_bm_tp *C_tilde = inla_bm_alloc(n, 0, 0);
	for (int i = 0; i < n; i++) {
		double sum = 0.0;
		for (int j = IMAX(0, i - C->lbw); j <= IMIN(C->n - 1, i + C->ubw); j++) {
			sum += BM_VAL(C, i, j);
		}
		BM_VAL(C_tilde, i, i) = 1.0 / sum;
	}
	arg->C_tilde = C_tilde;

	inla_bm_tp *G = inla_bm_alloc(n, 1, 1);
	for (int i = 0; i < n; i++) {
		if (i == 0) {
			BM_VAL(G, i, i) = 1.0 / h[0];
		} else if (i == n - 1) {
			BM_VAL(G, i, i) = 1.0 / h[n - 2];
		} else {
			BM_VAL(G, i, i) = 1.0 / h[i] + 1.0 / h[i - 1];
		}

		int j = i + 1;
		if (j < n) {
			BM_VAL(G, i, j) = BM_VAL(G, j, i) = -1.0 / h[i];
		}
	}
	arg->G = G;

	inla_bm_tp *M = inla_bm_alloc(n, 1, 1);
	for (int i = 0; i < n; i++) {
		if (i == 0) {
			BM_VAL(M, i, i) = -0.5;
		} else if (i == n - 1) {
			BM_VAL(M, i, i) = 0.5;
		}
		int j = i + 1;
		if (j < n) {
			BM_VAL(M, i, j) = 0.5;
			BM_VAL(M, j, i) = -0.5;
		}
	}
	arg->M = M;

	GMRFLib_graph_mk_linear(&(arg->graph), n, 2, 0);
	arg->cache = Calloc(GMRFLib_CACHE_LEN(), inla_prw2_cache_tp *);
	for (int i = 0; i < GMRFLib_CACHE_LEN(); i++) {
		arg->cache[i] = Calloc(1, inla_prw2_cache_tp);
		arg->cache[i]->range = -1;
		arg->cache[i]->Q = NULL;
	}

	return arg;
}

inla_bm_tp *inla_prw2_build_Q(int thread_id, inla_prw2_arg_tp *arg)
{
	// return Q with prec = 1
	int n = arg->n;
	double range = GMRFLib_SET_RANGE(arg);
	double kappa = sqrt(12.0) / range;
	double tau = 1.0 / (4.0 * POW3(kappa));

	double a = -1.0;				       /* this is still a mystery: why isn't a=1.0 ? */
	double b = 2.0 * kappa;
	double c = SQR(kappa);

	assert(arg->G && arg->M && arg->C && arg->C_tilde);

	inla_bm_tp *H = inla_bm_sum(a, arg->G, b, arg->M, c, arg->C, NULL);
	inla_bm_tp *Ht = inla_bm_trans(H, NULL);
	inla_bm_tp *Q = inla_bm_mmm(Ht, arg->C_tilde, H, NULL);
	inla_bm_scale(tau, Q);

	const int m = 4;
	const int m2 = ISQR(m);
	int mm = n - m;
	int mmn_len = GMRFLib_align_len(mm * n, sizeof(double));
	double *xx = Calloc(2 * mmn_len, double);
	double *yy = xx + mmn_len;
	xx[0] = BM_VAL(Q, 0, 2);
	xx[0 + mm] = BM_VAL(Q, 1, 2);
	xx[1 + mm] = BM_VAL(Q, 1, 3);
	xx[mm - 2 + 2 * mm] = BM_VAL(Q, n - 2, n - 4);
	xx[mm - 1 + 2 * mm] = BM_VAL(Q, n - 2, n - 3);
	xx[mm - 1 + 3 * mm] = BM_VAL(Q, n - 1, n - 3);
	Memcpy(yy, xx, mm * m * sizeof(double));

	inla_bm_tp *QB = inla_bm_alloc(mm, 2, 0);
	for (int i = 0; i < mm; i++) {
		for (int j = IMAX(0, i - QB->lbw); j <= i; j++) {
			int ii = i + 2, jj = j + 2;
			BM_VAL(QB, i, j) = BM_VAL(Q, ii, jj);
		}
	}

	inla_bm_tp *chol_QB = inla_bm_chol(QB, NULL);
	inla_bm_nsolve(chol_QB, NULL, xx, 4);

	int idx[] = { 0, 1, n - 2, n - 1 };
	inla_bm_tp *S = inla_bm_alloc(m, m - 1, 0);
	for (int i = 0; i < m; i++) {
		BM_VAL(S, i, i) = 1.0;
		for (int j = 0; j < i; j++) {
			sBM_VAL(S, i, j) = inla_prw2_corfunc(arg->loc[idx[i]] - arg->loc[idx[j]], kappa);
		}
	}

	inla_bm_tp *chol_S = inla_bm_chol(S, NULL);
	inla_bm_tp *Sinv = inla_bm_partial_inv(chol_S, NULL);

#if 0
	// we might revisit this later
	inla_bm_tp *S0 = inla_bm_alloc(m, m - 1, 0);
	for (int i = 0; i < m; i++) {
		BM_VAL(S0, i, i) = 1.0;
		for (int j = 0; j < i; j++) {
			if ((i == 1 && j == 0) || (i == m - 1 && j == m - 2)) {
				sBM_VAL(S0, i, j) = sBM_VAL(S, i, j);
			}
		}
	}

	inla_bm_tp *chol_S0 = inla_bm_chol(S0, NULL);
	inla_bm_tp *S0inv = inla_bm_partial_inv(chol_S0, NULL);

	P(sBM_VAL(S, 1, 2));
	P(sBM_VAL(Sinv, 0, 0) - sBM_VAL(S0inv, 0, 0));
	inla_bm_fprintf(stdout, Sinv, "Sinv");
	inla_bm_fprintf(stdout, S0inv, "S0inv");

	// maybe use something else? its all about the difference in Sinv and S0inv ? 
	if (0) {
		double small = ABS(sBM_VAL(S, 1, 2));
		if (small > 0.2) {
			fprintf(stderr, "\n *** Warning *** [%s:%1d](%s) correlation between the boundaries is higher than small [%.3g]\n",
				__FILE__, __LINE__, __GMRFLib_FuncName, small);
			fflush(stderr);
		}
	}
#endif

	double correction[m2];
	GMRFLib_dfill(m2, 0.0, correction);
	for (int i = 0; i < m; i++) {
		for (int j = i; j < m; j++) {
			double sum = 0.0;
			for (int k = 0; k < mm; k++) {
				sum += yy[k + i * mm] * xx[k + j * mm];
			}
			correction[i + j * m] = correction[j + i * m] = sum + sBM_VAL(Sinv, i, j);
		}
	}

	BM_VAL(Q, 0, 0) = correction[0 + 0 * m];
	BM_VAL(Q, 1, 0) = correction[1 + 0 * m];
	BM_VAL(Q, 0, 1) = correction[0 + 1 * m];
	BM_VAL(Q, 1, 1) = correction[1 + 1 * m];
	BM_VAL(Q, n - 2, n - 2) = correction[2 + 2 * m];
	BM_VAL(Q, n - 2, n - 1) = correction[2 + 3 * m];
	BM_VAL(Q, n - 1, n - 2) = correction[3 + 2 * m];
	BM_VAL(Q, n - 1, n - 1) = correction[3 + 3 * m];

	inla_bm_tp *Qsym = inla_bm_sym(Q, NULL);

	if (0) {
		// print the correlation function computed from the FEM and the target one
		inla_bm_tp *chol = inla_bm_chol(Qsym, NULL);
		double *x = Calloc(n, double);
		x[n / 2] = 1.0;
		inla_bm_solve(chol, NULL, x);
#pragma omp critical (Name_6445b688a836bdfd82eba490147e8e86c75290f2)
		for (int i = 0; i < n; i++) {
			printf("COR with range = %f i= %1d dist = %f cor %f cor.true %f\n", range, i,
			       ABS(arg->loc[i] - arg->loc[n / 2]), x[i] / x[n / 2], inla_prw2_corfunc(arg->loc[i] - arg->loc[n / 2], kappa));
		}
		inla_bm_free(chol);
		Free(x);
	}

	Free(xx);
	inla_bm_free(H);
	inla_bm_free(Ht);
	inla_bm_free(Q);
	inla_bm_free(QB);
	inla_bm_free(S);
	inla_bm_free(Sinv);
	inla_bm_free(chol_QB);
	inla_bm_free(chol_S);

	return Qsym;
}

double inla_Qfunc_prw2(int thread_id, int i, int j, double *values, void *arg)
{
	int cache_idx;
	GMRFLib_CACHE_SET_IDX(cache_idx);

	inla_prw2_arg_tp *def = (inla_prw2_arg_tp *) arg;
	inla_prw2_cache_tp *cache = def->cache[cache_idx];

	double prec = GMRFLib_SET_PREC(def);
	double range = GMRFLib_SET_RANGE(def);

	if (cache->range != range) {
#pragma omp critical (Name_c8bc8620a5d8dba701c69cf75cf102e4s1d19c421)
		if (cache->range != range) {
			inla_bm_free(cache->Q);
			cache->Q = inla_prw2_build_Q(thread_id, def);
			cache->range = range;
		}
	}

	if (j >= 0) {
		return prec * sBM_VAL(cache->Q, i, j);
	} else {
		int k = 0;
		values[k++] = prec * sBM_VAL(cache->Q, i, i);
		for (int jj = 0; jj < def->graph->lnnbs[i]; jj++) {
			j = def->graph->lnbs[i][jj];
			values[k++] = prec * sBM_VAL(cache->Q, i, j);
		}
		return 0.0;
	}
}

void inla_prw2_pcprior_dist(double *rho, int n, double *d)
{
	// compute d=sqrt(...(rho))

	for (int i = 0; i < n; i++) {
		double lrho = log1p(rho[i] - 1.0);
		d[i] = sqrt(-(3.0 + rho[i]) * POW3(expm1(lrho)) / (1.0 + SQR(rho[i])));
	}
}

void inla_prw2_pcprior_range2rho(double *range, int n, double h_size, double *rho)
{
	for (int i = 0; i < n; i++) {
		rho[i] = exp(-SQRT12 * h_size / range[i]);
	}
}

void inla_prw2_pcprior_cdf_range(double *range, int n, double lambda, double h_size, double *cdf)
{
	// use 'cdf' as a work array
	inla_prw2_pcprior_range2rho(range, n, h_size, cdf);    // rho is now cdf
	inla_prw2_pcprior_dist(cdf, n, cdf);		       // dist is now cdf
	double normc = 1.0 / (1.0 - exp(-lambda * M_SQRT3));
	for (int i = 0; i < n; i++) {
		cdf[i] = exp(-lambda * cdf[i]) * normc;
	}
}

double priorfunc_prw2_pcprior_range(double *x, double *parameters)
{
	double range = exp(x[0]);
	double lambda = parameters[3];			       /* yes */
	double h = 1e-4;
	double r[2] = { range * exp(-h), range * exp(h) };
	double cdf[2];
	double dens;
	double h_size = (parameters[2] <= 0.0 ? 1.0 : parameters[2]);

	inla_prw2_pcprior_cdf_range(r, 2, lambda, h_size, cdf);
	dens = (cdf[1] - cdf[0]) / (2.0 * h) / range;
	return log(dens) + x[0];
}

double priorfunc_prw2_pcprior_range_calibrate_helper(double lambda, double r0, double alpha, double h_size)
{
	double cdf = 0.0;
	inla_prw2_pcprior_cdf_range(&r0, 1, lambda, h_size, &cdf);
	return (1.0 - cdf - alpha);
}

double priorfunc_prw2_pcprior_range_calibrate(double r0, double alpha, double h_size)
{
	// find lambda such that Prob(r > r0) = alpha

	int deriv = 0;
	int verbose = 1;				       /* leave this on */
	int iter_max = 1000;
	int iter = 0;
	double lam0 = 1.0;
	double lam1 = lam0 + 1.0;
	double f0 = priorfunc_prw2_pcprior_range_calibrate_helper(lam0, r0, alpha, h_size);
	double f1 = priorfunc_prw2_pcprior_range_calibrate_helper(lam1, r0, alpha, h_size);
	double fac = 2.0;
	double tref = 0.0;

	assert(h_size > 0);
	if (verbose) {
		tref -= GMRFLib_timer();
		printf("\n\t\tcalibrate prw2_pcprior for range r0 = %.4f alpha = %.4f h_size= %.4f\n", r0, alpha, h_size);
		printf("\t\t\t(lam0 f0) = (%.8f %.8f)\n", lam0, f0);
		printf("\t\t\t(lam1 f1) = (%.8f %.8f)\n", lam1, f1);
	}

	if (f0 * f1 > 0) {
		// same sign. need to do something
		deriv = ((f1 - f0) / (lam1 - lam0) > 0 ? 1.0 : -1.0);
		lam1 = lam0;
		if (f0 > 0) {
			while (f1 > 0) {
				lam1 *= (deriv > 0 ? 1.0 / fac : fac);
				f1 = priorfunc_prw2_pcprior_range_calibrate_helper(lam1, r0, alpha, h_size);
				if (verbose) {
					printf("\t\t\tsearch (lam1 f1) = (%.8f %.8f)\n", lam1, f1);
				}
				if (f1 > 0) {
					f0 = f1;
					lam0 = lam1;
				}
				iter++;
				assert(iter < iter_max);
			}
		} else {
			// f0 < 0
			while (f1 < 0) {
				lam1 *= (deriv > 0 ? fac : 1.0 / fac);
				f1 = priorfunc_prw2_pcprior_range_calibrate_helper(lam1, r0, alpha, h_size);
				if (verbose) {
					printf("\t\t\tsearch (lam1 f1) = (%.8f %.8f)\n", lam1, f1);
				}
				if (f1 < 0) {
					f0 = f1;
					lam0 = lam1;
				}
				iter++;
				assert(iter < iter_max);
			}
		}
	}
	assert(f0 * f1 < 0.0);

	// make sure f0 < 0 < f1. easier coding
	if (f0 > f1) {
		double tmp = f0;
		f0 = f1;
		f1 = tmp;

		tmp = lam0;
		lam0 = lam1;
		lam1 = tmp;
	}

	if (verbose) {
		printf("\t\t\tfound (lam0 f0<0) = (%.8f %.8f)\n", lam0, f0);
		printf("\t\t\tfound (lam1 f1>0) = (%.8f %.8f)\n", lam1, f1);
	}

	double lam_mid = 0.0, f_mid = 0.0, eps = 1e-5;
	while (1) {
		// interpolate in log-scale
		lam_mid = exp((log(lam1) * f0 - log(lam0) * f1) / (f0 - f1));
		f_mid = priorfunc_prw2_pcprior_range_calibrate_helper(lam_mid, r0, alpha, h_size);

		if (verbose) {
			printf("\t\t\t(lam_mid f_mid) = (%.8f %.8f) ERROR %.5f\n", lam_mid, f_mid, ABS(f_mid));
		}

		if (ABS(f_mid) < eps) {
			break;
		}
		if (f_mid < 0) {
			f0 = f_mid;
			lam0 = lam_mid;
		} else {
			f1 = f_mid;
			lam1 = lam_mid;
		}

		iter++;
		assert(iter < iter_max);
	}

	if (verbose) {
		printf("\t\t\tsolution found (lam_mid f_mid) = (%.8f %.8f)\n", lam_mid, f_mid);
		tref += GMRFLib_timer();
		printf("\t\t\ttime used %.3fs\n", tref);
	}

	return lam_mid;
}

#if defined(INLA_WITH_DEVEL)
void inla_bm_test()
{
	const int n = 151;
	double *loc = Calloc(n, double);
	for (int i = 0; i < n; i++) {
		loc[i] = i;
	}
	double rr = loc[n - 1] - loc[0];
	if (0)
		for (int i = 0; i < n; i++) {
			loc[i] = (loc[i] - loc[0]) / rr;
		}

	double **log_range, **log_prec;
	HYPER_NEW(log_range, log(0.2 * rr));
	HYPER_NEW(log_prec, log(1.0));

	P(loc[n - 1] - loc[0]);
	P(exp(log_range[0][0]));

	int thread_id = 0;
	inla_prw2_arg_tp *arg = inla_prw2_create(n, loc);
	arg->log_prec_omp = log_prec;
	arg->log_range_omp = log_range;
	inla_bm_tp *Q = inla_prw2_build_Q(thread_id, arg);
	inla_bm_tp *Qchol = inla_bm_chol(Q, NULL);
	inla_bm_tp *Qinv = inla_bm_partial_inv(Qchol, NULL);

	double *cor = Calloc(n, double);
	cor[n / 2] = 1.0;
	inla_bm_solve(Qchol, NULL, cor);

	double kappa = sqrt(12.0) / exp(log_range[0][0]);
	for (int i = 0; i < n; i++) {
		printf("CORR %f %f %f\n", loc[i], cor[i], inla_prw2_corfunc(loc[n / 2] - loc[i], kappa));
	}

	for (int i = 0; i < n; i++) {
		printf(" %.3g", BM_VAL(Qinv, i, i));
	}
	printf("\n");

	P(exp(**log_range));
	for (int i = 0; i < IMIN(4, n - 2); i++) {
		printf("i %d %g %g %g\n", i,
		       inla_Qfunc_prw2(0, i, i, NULL, (void *) arg),
		       inla_Qfunc_prw2(0, i, i + 1, NULL, (void *) arg), inla_Qfunc_prw2(0, i, i + 2, NULL, (void *) arg));
	}

	**log_range += 0.1;
	P(exp(**log_range));
	for (int i = 0; i < IMIN(4, n - 2); i++) {
		printf("i %d %g %g %g\n", i,
		       inla_Qfunc_prw2(0, i, i, NULL, (void *) arg),
		       inla_Qfunc_prw2(0, i, i + 1, NULL, (void *) arg), inla_Qfunc_prw2(0, i, i + 2, NULL, (void *) arg));
	}

	inla_bm_free(Q);
	inla_bm_free(Qchol);
	inla_bm_free(Qinv);
	Free(loc);
}
void inla_prw2_test(void)
{
	priorfunc_prw2_pcprior_range_calibrate(10.0, 0.9, 1.0);
}

#else
void inla_bm_test()
{
}
void inla_prw2_test(void)
{
}
#endif
