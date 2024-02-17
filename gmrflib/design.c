
/* design.c
 * 
 * Copyright (C) 2006-2024 Havard Rue
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

#include <stddef.h>
#include <float.h>
#include <math.h>
#include <strings.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

#include "designP.h"					       /* define the designs */

int GMRFLib_design_eb(GMRFLib_design_tp **design, int nhyper)
{
	// create an EB design
	*design = Calloc(1, GMRFLib_design_tp);
	(*design)->nfactors = nhyper;
	(*design)->nexperiments = 1;
	(*design)->experiment = Calloc(1, double *);
	(*design)->int_weight = Calloc(1, double);
	(*design)->std_scale = GMRFLib_TRUE;
	(*design)->experiment[0] = Calloc(nhyper, double);
	(*design)->int_weight[0] = 1.0;

	return GMRFLib_SUCCESS;
}

int GMRFLib_design_grid(GMRFLib_design_tp **design, int nhyper)
{
	// this returns the int.strategy="grid" design for the experimental mode, which
	// is different for nhyper=1 and 2, and otherwise equal to the ccd

	// only do dim=1 or 2 here
	if (!(nhyper == 1 || nhyper == 2)) {
		return GMRFLib_design_ccd(design, nhyper);
	}

	double x1[] = {
		-3.5,
		-2.5,
		-1.75,
		-1.0,
		-0.5,
		0.0,
		0.5,
		1.0,
		1.75,
		2.5,
		3.5
	};

	double w1[] = {
		3.187537795,
		1.811358205,
		1.937929918,
		1.431919577,
		1.288639321,
		1.0,
		1.288639321,
		1.431919577,
		1.937929918,
		1.811358205,
		3.187537795
	};

	double x2[] = {
		-2.25, -1.25,
		-2.25, -0.5,
		-2.25, 0.0,
		-2.25, 0.5,
		-2.25, 1.25,
		-1.25, -2.25,
		-1.25, -1.25,
		-1.25, -0.5,
		-1.25, 0.0,
		-1.25, 0.5,
		-1.25, 1.25,
		-1.25, 2.25,
		-0.5, -2.25,
		-0.5, -1.25,
		-0.5, -0.5,
		-0.5, 0.0,
		-0.5, 0.5,
		-0.5, 1.25,
		-0.5, 2.25,
		0.0, -2.25,
		0.0, -1.25,
		0.0, -0.5,
		0.0, 0.0,
		0.0, 0.5,
		0.0, 1.25,
		0.0, 2.25,
		0.5, -2.25,
		0.5, -1.25,
		0.5, -0.5,
		0.5, 0.0,
		0.5, 0.5,
		0.5, 1.25,
		0.5, 2.25,
		1.25, -2.25,
		1.25, -1.25,
		1.25, -0.5,
		1.25, 0.0,
		1.25, 0.5,
		1.25, 1.25,
		1.25, 2.25,
		2.25, -1.25,
		2.25, -0.5,
		2.25, 0.0,
		2.25, 0.5,
		2.25, 1.25
	};

	double w2[] = {
		2.277250821,
		1.248862019,
		1.93160554,
		1.248862019,
		2.277250821,
		2.277250821,
		1.389904145,
		0.762234217,
		1.17894196,
		0.762234217,
		1.389904145,
		2.277250821,
		1.248862019,
		0.762234217,
		0.4180151587,
		0.646540918,
		0.4180151587,
		0.762234217,
		1.248862019,
		1.93160554,
		1.17894196,
		0.646540918,
		1.0,
		0.646540918,
		1.17894196,
		1.93160554,
		1.248862019,
		0.762234217,
		0.4180151587,
		0.646540918,
		0.4180151587,
		0.762234217,
		1.248862019,
		2.277250821,
		1.389904145,
		0.762234217,
		1.17894196,
		0.762234217,
		1.389904145,
		2.277250821,
		2.277250821,
		1.248862019,
		1.93160554,
		1.248862019,
		2.277250821
	};

	int len = (nhyper == 1 ? sizeof(w1) / sizeof(double) : sizeof(w2) / sizeof(double));

	*design = Calloc(1, GMRFLib_design_tp);
	(*design)->nfactors = nhyper;
	(*design)->nexperiments = len;
	(*design)->experiment = Calloc(len, double *);
	(*design)->int_weight = Calloc(len, double);
	(*design)->std_scale = GMRFLib_TRUE;

	switch (nhyper) {
	case 1:
		for (int j = 0; j < len; j++) {
			(*design)->experiment[j] = Calloc(nhyper, double);
			Memcpy((*design)->experiment[j], x1 + j, nhyper * sizeof(double));
			(*design)->int_weight[j] = w1[j];
		}
		break;
	case 2:
		for (int j = 0; j < len; j++) {
			(*design)->experiment[j] = Calloc(nhyper, double);
			Memcpy((*design)->experiment[j], x2 + 2 * j, nhyper * sizeof(double));
			(*design)->int_weight[j] = w2[j];
		}
		break;
	default:
		assert(0 == 1);
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_design_ccd(GMRFLib_design_tp **design, int nfactors)
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
		GMRFLib_dscale(nfactors, scale, (*design)->experiment[j]);
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_design_read(GMRFLib_design_tp **design, GMRFLib_matrix_tp *D, int std_scale)
{
	/*
	 * read the design from D
	 */
	int j, k, nfac, nex;
	double sum_w = 0.0;

	*design = Calloc(1, GMRFLib_design_tp);
	(*design)->nfactors = nfac = D->ncol - 1;
	(*design)->nexperiments = nex = D->nrow;
	assert(nex > 0);
	(*design)->experiment = Calloc(nex, double *);
	(*design)->int_weight = Calloc(nex, double);
	(*design)->std_scale = (std_scale ? GMRFLib_TRUE : GMRFLib_FALSE);

	for (j = 0; j < nex; j++) {
		(*design)->experiment[j] = Calloc(nfac, double);
		for (k = 0; k < nfac; k++) {
			(*design)->experiment[j][k] = GMRFLib_matrix_get(j, k, D);
		}
		(*design)->int_weight[j] = GMRFLib_matrix_get(j, nfac, D);
		assert((*design)->int_weight[j] >= 0.0);
		sum_w += (*design)->int_weight[j];
	}
	for (j = 0; j < nex; j++) {
		(*design)->int_weight[j] /= sum_w;
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_design_free(GMRFLib_design_tp *design)
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

int GMRFLib_design_print(FILE *fp, GMRFLib_design_tp *design)
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
	fprintf(fp, " weight");
	fprintf(fp, "\n");

	for (i = 0; i < design->nexperiments; i++) {
		fprintf(fp, "\t%-.3d\t", i);
		for (j = 0; j < design->nfactors; j++) {
			fprintf(fp, " %6.3f", design->experiment[i][j]);
		}
		fprintf(fp, " %6.4f", design->int_weight[i]);
		fprintf(fp, "\n");
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_design_prune(GMRFLib_design_tp *design, double prob)
{
	if (!design) {
		return GMRFLib_SUCCESS;
	}

	const int debug = 0;
	int *idx = NULL, i = 0;
	double *w = NULL, sumw = 0.0;

	int n = design->nexperiments, m;
	w = Calloc(n, double);
	idx = Calloc(n, int);

	for (i = 0; i < n; i++) {
		idx[i] = i;
		w[i] = design->int_weight[i];
	}

	sumw = GMRFLib_dsum(n, w);
	GMRFLib_dscale(n, 1.0 / sumw, w);
	GMRFLib_qsort2(w, n, sizeof(double), idx, sizeof(int), GMRFLib_dcmp_r);

	sumw = 0.0;
	for (i = 0; i < n; i++) {
		sumw += w[i];
		if (debug) {
			printf("GMRFLib_design_prune: i idx w accw %d %d %f %f\n", i, idx[i], w[i], sumw);
		}
		if (sumw > prob) {
			break;
		}
	}
	m = i + 1;

	double *ww = Calloc(m, double);
	double **ex = Calloc(m, double *);

	for (i = 0; i < m; i++) {
		ex[i] = Calloc(design->nfactors, double);
		Memcpy(ex[i], design->experiment[idx[i]], design->nfactors * sizeof(double));
		ww[i] = design->int_weight[idx[i]];
	}

	sumw = GMRFLib_dsum(m, ww);
	GMRFLib_dscale(m, 1.0 / sumw, ww);

	for (i = 0; i < design->nexperiments; i++) {
		Free(design->experiment[i]);
	}
	Free(design->int_weight);
	design->experiment = ex;
	design->int_weight = ww;
	design->nexperiments = m;

	return GMRFLib_SUCCESS;
}
