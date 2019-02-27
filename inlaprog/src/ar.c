
/* ar.c
 * 
 * Copyright (C) 2012 Havard Rue
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
#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = HGVERSION;

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

#include "ar.h"
#include "inla.h"

#define PMATRIX(_M, _dim_i, _dim_j, _msg)				\
	if (debug) {							\
		int _i, _j;						\
		printf("\n%s (%1d x %1d)\n", _msg, _dim_i, _dim_j);	\
		for(_i = 0; _i < _dim_i; _i++) {			\
			printf("\t");					\
			for(_j = 0; _j < _dim_j; _j++){			\
				printf(" %10.6f", (_M)[ _i + (_dim_i) * _j]); \
			}						\
			printf("\n");					\
		}							\
		printf("\n");						\
	}


/* 
   functions for the AR(p) model; the pacf2phi and phi2pacf are taken from R's arima.c
 */
int ar_pacf2phi(int p, double *pacf, double *phi)
{
	/*
	 * From arima.c 
	 */
	int j, k;
	double a, *work;

	assert(p > 0);
	work = Calloc(p, double);

	for (j = 0; j < p; j++) {
		work[j] = phi[j] = pacf[j];
	}

	/*
	 * Step two: run the Durbin-Levinson recursions to find phi_{j.}, ( j = 2, ..., p and phi_{p.} are the autoregression coefficients
	 */
	for (j = 1; j < p; j++) {
		a = phi[j];
		for (k = 0; k < j; k++) {
			work[k] -= a * phi[j - k - 1];
		}
		for (k = 0; k < j; k++) {
			phi[k] = work[k];
		}
	}

	Free(work);
	return GMRFLib_SUCCESS;
}

int ar_phi2pacf(int p, double *phi, double *pacf)
{
	/*
	 * From arima.c 
	 */

	int j, k;
	double a, *work;

	assert(p > 0);
	work = Calloc(p, double);

	for (j = 0; j < p; j++) {
		work[j] = pacf[j] = phi[j];
	}

	/*
	 * Run the Durbin-Levinson recursions backwards to find the PACF phi_{j.} from the autoregression coefficients 
	 */
	for (j = p - 1; j > 0; j--) {
		a = pacf[j];
		for (k = 0; k < j; k++) {
			work[k] = (pacf[k] + a * pacf[j - k - 1]) / (1.0 - SQR(a));
		}
		for (k = 0; k < j; k++) {
			pacf[k] = work[k];
		}
	}

	Free(work);
	return GMRFLib_SUCCESS;
}

int ar_marginal_distribution(int p, double *pacf, double *prec, double *Q)
{
	/*
	 * from the set of partial correlation coefficients, PACF, return the marginal precision for a standard innovation process and the p x p -precision matrix
	 * for the first p components
	 */

	size_t i, j, lag, lag_idx, pdim, debug = 0;
	double *phi;

	assert(p > 0);
	pdim = (size_t) p;
	phi = Calloc(pdim, double);
	ar_pacf2phi(pdim, pacf, phi);

	/*
	 * I need to solve the Yule-Walker equations, not for 'phi' but for the correlations 
	 */
	gsl_matrix *A = gsl_matrix_calloc(pdim, pdim);
	gsl_vector *b = gsl_vector_calloc(pdim);
	gsl_vector *x = gsl_vector_calloc(pdim);

	for (i = 0; i < pdim; i++) {
		for (j = 0; j < pdim; j++) {
			if (i == j) {
				gsl_matrix_set(A, i, j, -1.0);
			} else {
				lag = (size_t) IABS((int) i - (int) j);
				lag_idx = lag - 1;	       /* as correlation for lag 1 is x[0], lag 2 is x[1] (lag 0 = 1 and not used) */
				gsl_matrix_set(A, i, lag_idx, phi[j] + gsl_matrix_get(A, i, lag_idx));
			}
		}
		gsl_vector_set(b, i, -phi[i]);
	}

	if (debug) {
		printf("A\n");
		GMRFLib_gsl_matrix_fprintf(stdout, A, NULL);
	}

	int s;
	gsl_permutation *perm = gsl_permutation_alloc(pdim);
	gsl_matrix *Sigma = gsl_matrix_calloc(pdim, pdim);
	gsl_matrix *Sigma_inverse = gsl_matrix_calloc(pdim, pdim);

	gsl_linalg_LU_decomp(A, perm, &s);
	gsl_linalg_LU_solve(A, perm, b, x);

	for (i = 0; i < pdim; i++) {
		for (j = 0; j < pdim; j++) {
			if (i == j) {
				gsl_matrix_set(Sigma, i, j, 1.0);
			} else {
				lag = (size_t) IABS((int) i - (int) j);
				lag_idx = lag - 1;
				gsl_matrix_set(Sigma, i, j, gsl_vector_get(x, lag_idx));
			}
		}
	}

	if (debug) {
		printf("Sigma\n");
		GMRFLib_gsl_matrix_fprintf(stdout, Sigma, NULL);
	}

	gsl_linalg_LU_decomp(Sigma, perm, &s);
	gsl_linalg_LU_invert(Sigma, perm, Sigma_inverse);

	for (i = 0; i < pdim; i++) {
		for (j = 0; j <= i; j++) {
			Q[i + pdim * j] = gsl_matrix_get(Sigma_inverse, i, j);
			Q[j + pdim * i] = gsl_matrix_get(Sigma_inverse, i, j);	/* so its exactly symmetric */
		}
	}

	*prec = 1.0;
	for (i = 0; i < pdim; i++) {
		*prec -= phi[i] * gsl_vector_get(x, i);
	}

	if (debug) {
		P(*prec);
	}

	gsl_permutation_free(perm);
	gsl_vector_free(x);
	gsl_vector_free(b);
	gsl_matrix_free(A);
	gsl_matrix_free(Sigma);
	gsl_matrix_free(Sigma_inverse);
	Free(phi);

	return GMRFLib_SUCCESS;
}

double ar_map_pacf(double arg, map_arg_tp typ, void *param)
{
	/*
	 * the map-function for the PACF
	 */
	switch (typ) {
	case MAP_FORWARD:
		/*
		 * extern = func(local) 
		 */
		return (2.0 * (exp((arg)) / (1.0 + exp((arg)))) - 1.0);
	case MAP_BACKWARD:
		/*
		 * local = func(extern) 
		 */
		return log((1.0 + arg) / (1.0 - arg));
	case MAP_DFORWARD:
		/*
		 * d_extern / d_local 
		 */
		return 2.0 * exp(arg) / ((SQR(1.0 + exp(arg))));
	case MAP_INCREASING:
		/*
		 * return 1.0 if montone increasing and 0.0 otherwise 
		 */
		return 1.0;
	default:
		abort();
	}
	abort();
	return 0.0;
}

double Qfunc_ar(int i, int j, void *arg)
{
	ar_def_tp *def = (ar_def_tp *) arg;

	if (IABS(i - j) > def->p) {
		return 0.0;
	}
	if (def->p == 0) {
		return exp(def->log_prec[GMRFLib_thread_id][0]);
	}

	int debug = 0, ii, jj, eq, dimQ, id;
	assert(def->n >= 2 * def->p);

	dimQ = 2 * def->p + 1;
	id = GMRFLib_thread_id + omp_get_thread_num() * GMRFLib_MAX_THREADS;
	eq = 1;

	for (ii = 0; ii < def->p && eq; ii++) {
		if (def->pacf_intern[ii][GMRFLib_thread_id][0] != def->hold_pacf_intern[id][ii]) {
			eq = 0;
		}
	}

	if (eq) {
		/*
		 * use what we already have
		 */
		int node, nnode;
		double Qmarg_contrib = 0.0, val;

		node = IMIN(i, j);
		nnode = IMAX(i, j);
		if (nnode < def->p) {			       /* recalling (i,j) starts from (0,0) */
			ii = node;
			jj = nnode;
			Qmarg_contrib = def->hold_Qmarg[id][ii + jj * def->p];	/* contribution from the marginal distribution for the first p x's. 
										 */
		} else if (nnode >= def->n - def->p) {
			ii = dimQ - def->p + node - (def->n - def->p);
			jj = dimQ - def->p + nnode - (def->n - def->p);
		} else {
			ii = def->p;
			jj = def->p + (nnode - node);
		}
		assert(LEGAL(ii, dimQ));
		assert(LEGAL(jj, dimQ));

		val = exp(def->log_prec[GMRFLib_thread_id][0]) * (Qmarg_contrib + def->hold_Q[id][ii + jj * dimQ]);
		return (val);
	} else {
		/*
		 * Build the Qmatrix 
		 */
		int k;
		double *phi, *pacf, *L, *Q, *Qmarg, prec;

		phi = Calloc(def->p, double);
		pacf = Calloc(def->p, double);
		L = Calloc(ISQR(dimQ), double);
		Q = Calloc(ISQR(dimQ), double);
		Qmarg = Calloc(ISQR(def->p), double);

		for (ii = 0; ii < def->p; ii++) {
			pacf[ii] = ar_map_pacf(def->pacf_intern[ii][GMRFLib_thread_id][0], MAP_FORWARD, NULL);
		}
		ar_marginal_distribution(def->p, pacf, &prec, Qmarg);
		PMATRIX(Qmarg, def->p, def->p, "Qmarg");
		PMATRIX(&prec, 1, 1, "Prec");
		ar_pacf2phi(def->p, pacf, phi);

		PMATRIX(pacf, def->p, 1, "pacf");
		PMATRIX(phi, def->p, 1, "phi");

		/*
		 * make L, where Lx = z ~ N(0,I) 
		 */
		for (ii = def->p; ii < dimQ; ii++) {
			L[ii + ii * dimQ] = 1.0;
			for (jj = ii - def->p, k = def->p - 1; jj < ii; jj++, k--) {
				L[ii + jj * dimQ] = -phi[k];
			}
		}
		PMATRIX(L, dimQ, dimQ, "Matrix L");

		/*
		 * Q = L' L 
		 */
		for (ii = 0; ii < dimQ; ii++) {
			for (jj = 0; jj < dimQ; jj++) {
				double tmp = 0.0;
				for (k = 0; k < dimQ; k++) {
					tmp += L[k + ii * dimQ] * L[k + jj * dimQ];
				}
				Q[ii + jj * dimQ] = tmp;
			}
		}
		PMATRIX(Q, dimQ, dimQ, "Matrix Q = L' L");

		for (ii = 0; ii < ISQR(dimQ); ii++) {
			Q[ii] /= prec;
		}
		PMATRIX(Q, dimQ, dimQ, "Matrix Q = L' L normalised");

		Free(def->hold_Qmarg[id]);
		Free(def->hold_Q[id]);
		def->hold_Qmarg[id] = Qmarg;
		def->hold_Q[id] = Q;

		for (ii = 0; ii < def->p; ii++) {
			def->hold_pacf_intern[id][ii] = def->pacf_intern[ii][GMRFLib_thread_id][0];
		}

		Free(phi);
		Free(pacf);
		Free(L);

		return Qfunc_ar(i, j, arg);		       /* recursive call */
	}
	assert(0 == 1);

	return 0.0;
}
int ar_test1()
{
	if (1) {
		GMRFLib_graph_tp *g;
		ar_def_tp def;
		double pacf[2] = { 0.5, 0.25 };

		int i, j, k;

		def.n = 10;
		def.p = 2;

		HYPER_NEW(def.log_prec, 0.0);
		def.pacf_intern = Calloc(def.p, double **);
		for (i = 0; i < def.p; i++) {
			double val = pacf[i];
			HYPER_NEW(def.pacf_intern[i], ar_map_pacf(val, MAP_BACKWARD, NULL));
		}

		/*
		 * easier if the storage is setup here 
		 */
		def.hold_pacf_intern = Calloc(ISQR(GMRFLib_MAX_THREADS), double *);
		def.hold_Q = Calloc(ISQR(GMRFLib_MAX_THREADS), double *);
		def.hold_Qmarg = Calloc(ISQR(GMRFLib_MAX_THREADS), double *);
		for (i = 0; i < ISQR(GMRFLib_MAX_THREADS); i++) {
			def.hold_pacf_intern[i] = Calloc(def.p, double);
			for (j = 0; j < def.p; j++) {
				def.hold_pacf_intern[i][j] = GMRFLib_uniform();
			}
		}

		GMRFLib_make_linear_graph(&g, def.n, def.p, 0);
		GMRFLib_print_Qfunc(stdout, g, Qfunc_ar, &def);

		if (0) {
			FILE *fp = fopen("Q.dat", "w");
			for (i = 0; i < def.n; i++)
				for (j = 0; j < def.n; j++)
					fprintf(fp, "%.12g\n", Qfunc_ar(i, j, &def));
			fclose(fp);
		}

		double val = 0.0;
#pragma omp parallel for private(k, i, j)
		for (k = 0; k < 1000; k++) {
			for (i = 0; i < def.n; i++)
				for (j = 0; j < def.n; j++) {
					val += Qfunc_ar(i, j, &def);
				}
		}
		P(val);

		exit(0);
	}

	if (0) {
#define PMAX 10
		int p, i, j;
		double *pacf, *pacf2, *phi;

		pacf = Calloc(PMAX, double);
		pacf2 = Calloc(PMAX, double);
		phi = Calloc(PMAX, double);

		for (p = 1; p <= PMAX; p++) {

			/*
			 * create some pacf's 
			 */
			for (i = 0; i < p; i++) {
				pacf[i] = GMRFLib_uniform();
			}

			/*
			 * compute the phi's 
			 */
			ar_pacf2phi(p, pacf, phi);

			/*
			 * and its inverse 
			 */
			ar_phi2pacf(p, phi, pacf2);


			printf("Result for p = %d\n", p);
			for (i = 0; i < p; i++) {
				j = i + 1;		       /* so phi_j = phi[i], j=1...p */
				printf("j = %2d  \tpacf = %.6g  \tphi %.6g  \tpacf2 %.6g  \tdiff %.6g\n", j, pacf[i], phi[i],
				       pacf2[i], ABS(pacf[i] - pacf2[i]));
			}
			printf("\n");

			double prec;
			double *Q = Calloc(ISQR(p), double);

			ar_marginal_distribution(p, pacf, &prec, Q);

			printf("Marginal distribution\n");
			printf("\tPrec = %g\n", prec);

			for (i = 0; i < p; i++) {
				for (j = 0; j < p; j++) {
					printf(" %10.6g ", Q[i + p * j]);
				}
				printf("\n");
			}

			Free(Q);
		}

		Free(pacf);
		Free(pacf2);
		Free(phi);
	}

	return GMRFLib_SUCCESS;
}
