
/* distributions.c
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
  \file distributions.c
  \brief Functions for sampling from distributions
*/

#include <math.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif

#include <stdlib.h>
#include <stdio.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

/* Pre-hg-Id: $Id: distributions.c,v 1.33 2009/08/26 06:12:46 hrue Exp $ */

double GMRFLib_stdnormal(void)
{
	/*
	 * use the GSL-version 
	 */
	return gsl_ran_ugaussian(GMRFLib_rng);
}

/*!
  \brief Return a sample from a density proportional to \f$1+1/x\f$ in the range \f$[1/F,F]\f$.
*/
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

/*!
  \brief Compute the log-density of a Wishart density \f$\sim W_p(r, R^{-1})\f$
*/
double GMRFLib_Wishart_logdens(gsl_matrix * Q, double r, gsl_matrix * R)
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
	gsl_matrix *C;

	p = Q->size1;
	C = gsl_matrix_calloc(R->size1, R->size2);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Q, R, 0.0, C);
	trace = 0;
	for (i = 0; i < C->size1; i++) {
		trace += gsl_matrix_get(C, i, i);
	}
	logdens = -0.5 * trace + (r - (double) p - 1.0) / 2.0 * GMRFLib_gsl_spd_logdet(Q);

	log_c = r * (double) p / 2.0 * log(2.0) - r / 2.0 * GMRFLib_gsl_spd_logdet(R) + (double) p *((double) p - 1.0) / 4.0 * log(M_PI);
	for (i = 1; i <= p; i++) {
		log_c += gsl_sf_lngamma((r + 1.0 - (double) i) / 2.0);
	}
	logdens -= log_c;

	gsl_matrix_free(C);

	return logdens;
}
