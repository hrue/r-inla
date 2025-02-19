
/*********************************************************/

/* TAUCS                                                 */

/* Author: Sivan Toledo                                  */

/*********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "taucs.h"

/*********************************************************/

/*                                                       */

/*********************************************************/

#ifndef TAUCS_CORE_GENERAL
int taucs_dtl(ccs_solve_llt) (void *vL, double *x, double *b) {
	taucs_ccs_matrix * L = (taucs_ccs_matrix *) vL;

	int n;
	int i, j;
	int ip, jp;
	double Aij, Ajj, Aii;
	double *y;

	if (!(L->flags & TAUCS_TRIANGULAR)) {
		taucs_printf("taucs_ccs_solve_llt: factor matrix must be triangular\n");
		return -1;
	}
	if (!(L->flags & TAUCS_LOWER)) {
		taucs_printf("taucs_ccs_solve_llt: lower part must be represented\n");
		return -1;
	}

	n = L->n;

	y = (double *) taucs_malloc(n * sizeof(double));
	if (!y)
		return -1;

	for (i = 0; i < n; i++)
		x[i] = b[i];

	/*
	 * Solve L y = b = x 
	 */

	for (j = 0; j < n; j++) {

		/*
		 * we put diagonal elements first 
		 */
		ip = (L->colptr)[j];
		i = (L->rowind)[ip];
		assert(i == j);
		Ajj = L->values[ip];

		/*
		 * for (ip = (L->colptr)[j]; ip < (L->colptr)[j+1]; ip++) { i = (L->rowind)[ip]; if (i==j) { Ajj =
		 * L->values[ip]; break; } } 
		 */

		/*
		 * y[j] = x[j] / Ajj;
		 */
		y[j] = taucs_div(x[j], Ajj);

		for (ip = (L->colptr)[j] + 1; ip < (L->colptr)[j + 1]; ip++) {
			i = (L->rowind)[ip];
			Aij = L->values[ip];
			/*
			 * x[i] -= y[j]*Aij;
			 */
			x[i] = taucs_sub(x[i], taucs_mul(y[j], Aij));
		}

		/*
		 * for (ip = (L->colptr)[j]; ip < (L->colptr)[j+1]; ip++) { i = (L->rowind)[ip]; if (i != j) { Aij =
		 * L->values[ip]; x[i] -= y[j]*Aij; } } 
		 */
	}

	/*
	 * Solve L^T x = y 
	 */

	for (i = n - 1; i >= 0; i--) {

		for (jp = (L->colptr)[i] + 1; jp < (L->colptr)[i + 1]; jp++) {
			j = (L->rowind)[jp];
			Aij = taucs_conj(L->values[jp]);
			/*
			 * y[i] -= x[j]*Aij;
			 */
			y[i] = taucs_sub(y[i], taucs_mul(x[j], Aij));
		}
		/*
		 * for (jp = (L->colptr)[i]; jp < (L->colptr)[i+1]; jp++) { j = (L->rowind)[jp]; if (i != j) { Aij =
		 * L->values[jp]; y[i] -= x[j]*Aij; } } 
		 */

		jp = (L->colptr)[i];
		j = (L->rowind)[jp];
		Aii = L->values[jp];

		/*
		 * for (jp = (L->colptr)[i]; jp < (L->colptr)[i+1]; jp++) { j = (L->rowind)[jp]; if (i==j) { Aii =
		 * L->values[jp]; break; } } 
		 */

		/*
		 * x[i] = y[i] / Aii;
		 */
		x[i] = taucs_div(y[i], Aii);

	}

	taucs_free(y);

	return 0;
}


/*********************************************************/

/* LDL^T solve                                           */

/*********************************************************/

int taucs_dtl(ccs_solve_ldlt) (void *vL, double *x, double *b) {
	taucs_ccs_matrix * L = (taucs_ccs_matrix *) vL;

	int n;
	int i, j;
	int ip, jp;
	double Ajj = taucs_zero_const;			       /* just to suppress the warning */
	double Aij = taucs_zero_const;			       /* just to suppress the warning */
	double *y;

	/*
	 * taucs_printf("taucs_ccs_solve_ldlt: starting\n"); 
	 */

	if (!(L->flags & TAUCS_TRIANGULAR)) {
		taucs_printf("taucs_ccs_solve_ldlt: factor matrix must be triangular\n");
		return -1;
	}
	if (!(L->flags & TAUCS_LOWER)) {
		taucs_printf("taucs_ccs_solve_ldlt: lower part must be represented\n");
		return -1;
	}

	n = L->n;

	y = (double *) taucs_malloc(n * sizeof(double));
	if (!y)
		return -1;

	for (i = 0; i < n; i++)
		x[i] = b[i];

	/*
	 * Solve L y = b = x 
	 */

	/*
	 * taucs_printf("taucs_ccs_solve_ldlt: solving L y = b\n"); 
	 */

	for (j = 0; j < n; j++) {

		y[j] = x[j];
		for (ip = (L->colptr)[j] + 1; ip < (L->colptr)[j + 1]; ip++) {
			i = (L->rowind)[ip];
			Aij = L->values[ip];

			/*
			 * x[i] -= y[j]*Aij;
			 */
			x[i] = taucs_sub(x[i], taucs_mul(y[j], Aij));
		}
	}

	/*
	 * Solve D y = y 
	 */

	for (j = 0; j < n; j++) {

		/*
		 * we put diagonal elements first 
		 */
		ip = (L->colptr)[j];
		i = (L->rowind)[ip];
		assert(i == j);
		Ajj = L->values[ip];

		/*
		 * y[j] = y[j] / Ajj; 
		 */
		y[j] = taucs_div(y[j], Ajj);
	}

	/*
	 * Solve L^T x = y 
	 */

	/*
	 * taucs_printf("taucs_ccs_solve_ldlt: solving L^T x = y\n");
	 */

	for (i = n - 1; i >= 0; i--) {

		for (jp = (L->colptr)[i] + 1; jp < (L->colptr)[i + 1]; jp++) {
			j = (L->rowind)[jp];
			/*
			 * Aij = L->values[jp]; 
			 */
			Aij = taucs_conj(L->values[jp]);
			/*
			 * y[i] -= x[j]*Aij;
			 */
			y[i] = taucs_sub(y[i], taucs_mul(x[j], Aij));
		}
		x[i] = y[i];
	}

	taucs_free(y);

	return 0;
}

#endif							       /* #ifndef TAUCS_CORE_GENERAL */

#ifdef TAUCS_CORE_GENERAL
int taucs_ccs_solve_llt(void *vL, void *x, void *b)
{
	taucs_ccs_matrix *L = (taucs_ccs_matrix *) vL;
	if (L->flags & TAUCS_DOUBLE)
		return taucs_dccs_solve_llt(L, (double *) x, (double *) b);
	assert(0);
	return -1;
}

int taucs_ccs_solve_ldlt(void *vL, void *x, void *b)
{
	taucs_ccs_matrix *L = (taucs_ccs_matrix *) vL;
	if (L->flags & TAUCS_DOUBLE)
		return taucs_dccs_solve_ldlt(L, (double *) x, (double *) b);
	assert(0);
	return -1;

}
#endif							       /* TAUCS_CORE_GENERAL */
