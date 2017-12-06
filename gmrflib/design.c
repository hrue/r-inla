
/* design.c
 * 
 * Copyright (C) 2006 Havard Rue
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
  \file design.c
  \brief Functions to ...

*/

#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

/* Pre-hg-Id: $Id: design.c,v 1.18 2008/04/17 08:24:33 hrue Exp $ */

#include <stddef.h>
#include <float.h>
#include <math.h>
#include <strings.h>
#include <string.h>
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

#include "designP.h"					       /* define the designs */

int GMRFLib_get_design(GMRFLib_design_tp ** design, int nfactors)
{
	/*
	 * return the CCD design with nfactors in design.  the the design computed as described in: Sanchez, S. M. and
	 * P. J. Sanchez, "Very large fractional factorials and central composite designs," ACM Transactions on Modeling and
	 * Computer Simulation 15(4): 362-377.
	 * 
	 * return a scaled design so that the +/-1's are scaled so that z^Tz = 1 
	 */
	int i, j, k, ipos;
	double scale;

	GMRFLib_ASSERT(nfactors >= nfac_from && nfactors <= nfac_to, GMRFLib_EPARAMETER);

	i = nfac_from;
	ipos = 0;
	while (nfactors > i) {
		ipos += nexp[i] * i;
		i++;
	}
	*design = Calloc(1, GMRFLib_design_tp);
	(*design)->nfactors = nfactors;
	(*design)->nexperiments = nexp[nfactors];
	(*design)->experiment = Calloc(nexp[nfactors], double *);
	(*design)->int_weight = Calloc(nexp[nfactors], double);
	(*design)->std_scale = GMRFLib_TRUE;

	for (j = 0; j < nexp[nfactors]; j++) {
		(*design)->experiment[j] = Calloc(nfactors, double);
		(*design)->int_weight[j] = NAN;		       /* meaning that its undefined at this stage */
		scale = 0.0;
		for (k = 0; k < nfactors; k++) {
			(*design)->experiment[j][k] = points[ipos + j * nfactors + k];
			scale += SQR((*design)->experiment[j][k]);
		}
		if (!ISZERO(scale)) {			       /* yes, the origo is within the design points! */
			scale = 1.0 / sqrt(scale);
		}
		for (k = 0; k < nfactors; k++) {
			(*design)->experiment[j][k] *= scale;
		}
	}

	return GMRFLib_SUCCESS;
}
int GMRFLib_read_design(GMRFLib_design_tp ** design, GMRFLib_matrix_tp * D, int std_scale)
{
	/*
	 * read the design from D
	 */
	int j, k, nfac, nexp;
	double sum_w = 0.0;

	*design = Calloc(1, GMRFLib_design_tp);
	(*design)->nfactors = nfac = D->ncol - 1;
	(*design)->nexperiments = nexp = D->nrow;
	assert(nexp > 0);
	(*design)->experiment = Calloc(nexp, double *);
	(*design)->int_weight = Calloc(nexp, double);
	(*design)->std_scale = (std_scale ? GMRFLib_TRUE : GMRFLib_FALSE);

	for (j = 0; j < nexp; j++) {
		(*design)->experiment[j] = Calloc(nfac, double);
		for (k = 0; k < nfac; k++) {
			(*design)->experiment[j][k] = GMRFLib_matrix_get(j, k, D);
		}
		(*design)->int_weight[j] = GMRFLib_matrix_get(j, nfac, D);
		assert((*design)->int_weight[j] >= 0.0);
		sum_w += (*design)->int_weight[j];
	}
	for (j = 0; j < nexp; j++) {
		(*design)->int_weight[j] /= sum_w;
	}

	return GMRFLib_SUCCESS;
}
int GMRFLib_free_design(GMRFLib_design_tp * design)
{
	if (design) {
		int i;

		for (i = 0; i < design->nexperiments; i++) {
			Free(design->experiment[i]);
		}
		Free(design->experiment);
		Free(design->int_weight);
		Free(design);
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_print_design(FILE * fp, GMRFLib_design_tp * design)
{
	int i, j;

	if (!design) {
		return GMRFLib_SUCCESS;
	}
	if (!fp) {
		fp = stdout;
	}

	fprintf(fp, "\tDesign has %d factors and %d experiments, scale=%1s\n", design->nfactors, design->nexperiments,
		(design->std_scale ? "Standardised" : "UserScale"));
	fprintf(fp, "\t\t");
	for (j = 0; j < design->nfactors; j++) {
		fprintf(fp, "     z%1d", j);
	}
	fprintf(fp, " weight\n");

	for (i = 0; i < design->nexperiments; i++) {
		fprintf(fp, "\t\t");
		for (j = 0; j < design->nfactors; j++) {
			fprintf(fp, " %6.3f", design->experiment[i][j]);
		}
		fprintf(fp, " %6.4f\n", design->int_weight[i]);
	}

	return GMRFLib_SUCCESS;
}
