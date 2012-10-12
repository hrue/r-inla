
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
 *       H{\aa}vard Rue
 *       Department of Mathematical Sciences
 *       The Norwegian University of Science and Technology
 *       N-7491 Trondheim, Norway
 *       Voice: +47-7359-3533    URL  : http://www.math.ntnu.no/~hrue  
 *       Fax  : +47-7359-3524    Email: havard.rue@math.ntnu.no
 *
 */
#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = HGVERSION;

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

#include "ar.h"

/* 
   functions for the AR(p) model; some are taken from R's arima.c
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
			work[k] = (pacf[k] + a * pacf[j - k - 1]) / (1.0 - a * a);
		}
		for (k = 0; k < j; k++) {
			pacf[k] = work[k];
		}
	}

	Free(work);
	return GMRFLib_SUCCESS;
}

int ar_test1()
{
	GMRFLib_graph_tp *g;
	ar_def_tp def;

	def.n = 10;
	def.p = 3;

	GMRFLib_make_linear_graph(&g, def.n, def.p, 0);
	int i, j, jj;

	for(i=0; i<def.n; i++){
		for(j = i;  j < IMIN(def.n, i+def.p); j++){
			Qfunc_ar_core(i, j, &def);
		}
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
				j = i + 1;			       /* so phi_j = phi[i], j=1...p */
				printf("j = %2d  \tpacf = %.6g  \tphi %.6g  \tpacf2 %.6g  \tdiff %.6g\n", j, pacf[i], phi[i], pacf2[i], ABS(pacf[i] - pacf2[i]));
			}
			printf("\n");

			double prec;
			double *Q;

			ar_marginal_distribution(p, pacf, &prec, &Q);

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

int ar_marginal_distribution(int p, double *pacf, double *prec, double **Q)
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
				lag = (size_t) IABS(i - j);
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
				lag = (size_t) IABS(i - j);
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

	*Q = Calloc(ISQR(pdim), double);
	for (i = 0; i < pdim; i++) {
		for (j = 0; j <= i; j++) {
			(*Q)[i + pdim * j] = gsl_matrix_get(Sigma_inverse, i, j);
			(*Q)[j + pdim * i] = gsl_matrix_get(Sigma_inverse, i, j);	/* so its exactly symmetric */
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

double Qfunc_ar_core(int i, int j, void *arg)
{
	char *phi[] = {"-1", "1", "2",  "3"};
	

	ar_def_tp *def = (ar_def_tp *) arg;

	assert(def->n >= 2*def->p);

	if (IABS(i-j) > def->p) {
		return 0.0;
	}
	
	int ii, jj, ki, kj, diff;


	return 0;
}

				
			
			
		
