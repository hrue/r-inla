#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

double GMRFLib_stdnormal(void)
{
	/*
	 * use the GSL-version 
	 */
	return gsl_ran_ugaussian(GMRFLib_rng);
}

double GMRFLib_scale_proposal(double F)
{
	if (F <= 1.0) {
		return 1.0;
	} else {
		double len = F - 1 / F;

		if (GMRFLib_uniform() < len / (len + 2 * log(F))) {
			return (1 / F + len * GMRFLib_uniform());
		} else {
			return pow(F, 2.0 * GMRFLib_uniform() - 1.0);
		}
	}
}

double GMRFLib_Wishart_logdens(gsl_matrix *Q, double r, gsl_matrix *R)
{
	/*
	 * this function computes the log-density of the Wishart_p(r, R^-1) density
	 * 
	 * \pi(Q) = ...  |Q|^{(r-p-1)/2} exp( -1/2 trace(Q R))
	 * 
	 * where r >= p+1.
	 * 
	 * This routine recomputes stuff that are constant over repeated evaluations; if this is important to avoid, then write a variant... 
	 */

	GMRFLib_ASSERT_RETVAL(Q, GMRFLib_EINVARG, 0.0);
	GMRFLib_ASSERT_RETVAL(R, GMRFLib_EINVARG, 0.0);
	GMRFLib_ASSERT_RETVAL((Q->size1 == Q->size2) && (R->size1 == R->size2) && (Q->size1 == R->size1), GMRFLib_EINVARG, 0.0);
	GMRFLib_ASSERT_RETVAL(r >= (double) Q->size1 + 1.0, GMRFLib_EINVARG, 0.0);	/* r >= p+1 */

	double logdens, trace, log_c;
	size_t p, i;
	gsl_matrix *C = NULL;
	gsl_matrix *QQ = NULL;

	QQ = GMRFLib_gsl_duplicate_matrix(Q);
	GMRFLib_gsl_ensure_spd(QQ, FLT_EPSILON, NULL);
	p = QQ->size1;
	C = gsl_matrix_calloc(R->size1, R->size2);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, QQ, R, 0.0, C);
	trace = 0;
	for (i = 0; i < C->size1; i++) {
		trace += gsl_matrix_get(C, i, i);
	}
	logdens = -0.5 * trace + (r - (double) p - 1.0) / 2.0 * GMRFLib_gsl_spd_logdet(QQ);

	log_c = r * (double) p / 2.0 * log(2.0) - r / 2.0 * GMRFLib_gsl_spd_logdet(R) + (double) p *((double) p - 1.0) / 4.0 * log(M_PI);
	for (i = 1; i <= p; i++) {
		log_c += gsl_sf_lngamma((r + 1.0 - (double) i) / 2.0);
	}
	logdens -= log_c;

	gsl_matrix_free(C);
	gsl_matrix_free(QQ);

	return logdens;
}
