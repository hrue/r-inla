
/*********************************************************/

/* TAUCS                                                 */

/* Author: Sivan Toledo                                  */

/*********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "taucs.h"

#ifndef TAUCS_CORE
#       error "You must define TAUCS_CORE to compile this file"
#endif

/*********************************************************/

/* CCS                                                   */

/*********************************************************/

#ifdef TAUCS_CORE_GENERAL
void taucs_ccs_free(taucs_ccs_matrix *matrix)
{
	taucs_dccs_free(matrix);
}
#endif							       /* TAUCS_CORE_GENERAL */

/* 
   Here the generic and type specific routines are different
   due to a historical accident, which forces us to set the
   flags again in the generic routine.
*/

#ifdef TAUCS_CORE_GENERAL
taucs_ccs_matrix *taucs_ccs_create(int m, int n, int nnz, int flags)
{
	taucs_ccs_matrix *A = NULL;
	if (flags & TAUCS_DOUBLE)
		A = taucs_dccs_create(m, n, nnz);
	if (A) {
		A->flags = flags;
		return A;
	} else {
		assert(0);
	}
	return NULL;
}
#endif							       /* TAUCS_CORE_GENERAL */

#ifndef TAUCS_CORE_GENERAL
taucs_ccs_matrix *taucs_dtl(ccs_create) (int m, int n, int nnz) {
	taucs_ccs_matrix * matrix;

	matrix = (taucs_ccs_matrix *) taucs_malloc(sizeof(taucs_ccs_matrix));
	if (!matrix) {
		taucs_printf("taucs_ccs_create: out of memory\n");
		return NULL;
	}
	matrix->flags = TAUCS_DOUBLE;
	matrix->n = n;
	matrix->m = m;
	matrix->colptr = (int *) taucs_malloc((n + 1) * sizeof(int));
	matrix->rowind = (int *) taucs_malloc(nnz * sizeof(int));
	matrix->values = (double *) taucs_malloc(nnz * sizeof(double));
	if (!(matrix->colptr) || !(matrix->rowind) || !(matrix->values)) {
		taucs_printf("taucs_ccs_create: out of memory (n=%d, nnz=%d)\n", n, nnz);
		taucs_free(matrix->colptr);
		taucs_free(matrix->rowind);
		taucs_free(matrix->values);
		taucs_free(matrix);
		return NULL;
	}

	return matrix;
}

void taucs_dtl(ccs_free) (taucs_ccs_matrix * matrix) {
	if (!matrix)
		return;

	taucs_free(matrix->rowind);
	taucs_free(matrix->colptr);
	taucs_free(matrix->values);
	taucs_free(matrix);
}

#endif							       /* TAUCS_CORE_GENERAL */
