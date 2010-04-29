
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
 *       H{\aa}vard Rue
 *       Department of Mathematical Sciences
 *       The Norwegian University of Science and Technology
 *       N-7491 Trondheim, Norway
 *       Voice: +47-7359-3533    URL  : http://www.math.ntnu.no/~hrue  
 *       Fax  : +47-7359-3524    Email: havard.rue@math.ntnu.no
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

	for (j = 0; j < nexp[nfactors]; j++) {
		(*design)->experiment[j] = Calloc(nfactors, double);

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
#undef FNM_DIM
#undef FNM_DAT
}

int GMRFLib_free_design(GMRFLib_design_tp * design)
{
	if (design) {
		int i;

		for (i = 0; i < design->nexperiments; i++) {
			Free(design->experiment[i]);
		}
		Free(design->experiment);
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

	fprintf(fp, "Design has %d factors and %d experiments\n", design->nfactors, design->nexperiments);
	for (i = 0; i < design->nexperiments; i++) {
		fprintf(fp, "\t");
		for (j = 0; j < design->nfactors; j++) {
			fprintf(fp, " %6.3f", design->experiment[i][j]);
		}
		fprintf(fp, "\n");
	}

	return GMRFLib_SUCCESS;
}
