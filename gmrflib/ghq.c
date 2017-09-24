
/* ghq.c
 * 
 * Copyright (C) 2006-2006 Havard Rue
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
  \file ghq.c
  \brief Utilities for Gauss-Hermite quadrature.

*/

#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

/* Pre-hg-Id: $Id: ghq.c,v 1.20 2007/05/27 13:38:44 hrue Exp $ */

#include <math.h>
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"
#include "GMRFLib/hashP.h"

int GMRFLib_ghq__intern(double *x, double *w, int n)
{
	int i, its, j, m;
	double p1, p2, p3, pp, z = 0, z1;
	double eps = 3.0e-14;				       /* shown to work fine for n < 199 */

	/*
	 * this function breaks down for n >= 199 
	 */
	GMRFLib_ASSERT(n < 199, GMRFLib_EINVARG);

	x--;						       /* indices from 1...n */
	if (w) {
		w--;					       /* in this routine */
	}

	m = (n + 1) / 2;
	for (i = 1; i <= m; i++) {
		switch (i) {
		case 1:
			z = sqrt(2.0 * n + 1.0) - 1.85575 * pow(2.0 * n + 1.0, -0.16667);
			break;
		case 2:
			z -= 1.14 * pow((double) n, 0.426) / z;
			break;
		case 3:
			z = 1.86 * z - 0.86 * x[1];
			break;
		case 4:
			z = 1.91 * z - 0.91 * x[2];
			break;
		default:
			z = 2.0 * z - x[i - 2];
		}
		its = 0;
		while (1) {
			its++;

			p1 = 0.7511255444649425;
			p2 = 0.0;
			for (j = 1; j <= n; j++) {
				p3 = p2;
				p2 = p1;
				p1 = z * sqrt(2.0 / j) * p2 - sqrt((j - 1.0) / j) * p3;
			}
			pp = sqrt(2.0 * n) * p2;
			z1 = z;
			z = z1 - p1 / pp;
			if (ABS(z - z1) <= eps)
				break;			       /* converged */
			GMRFLib_ASSERT(its <= 500, GMRFLib_ESNH);
		}

		x[i] = z;
		x[n + 1 - i] = -z;

		if (w) {
			w[i] = 2.0 / (pp * pp);
			w[n + 1 - i] = w[i];
		}
	}
	return 0;
}
int GMRFLib_ghq_abscissas(double **xp, int n)
{
	return GMRFLib_ghq(xp, NULL, n);
}
int GMRFLib_ghq_weights(double **wp, int n)
{
	return GMRFLib_ghq(NULL, wp, n);
}
int GMRFLib_ghq(double **xp, double **wp, int n)
{
	/*
	 * return a pointer to the abscissas and the weights for the Gauss-Hermite quadrature for given `n' with kernel
	 * exp(-x^2/2)/sqrt(2*pi)
	 * 
	 * both abscissas and the weights are optional
	 * 
	 * an internal storage is used to store previously computed abscissas and weights 
	 */

	static map_ivp abscissas;			       /* keep previous computed elements here */

#pragma omp threadprivate(abscissas)
	static map_ivp weights;				       /* keep previous computed elements here */

#pragma omp threadprivate(weights)
	static int first = 1;

#pragma omp threadprivate(first)

	int i;
	double *x, *w;
	void *ptr, *pptr;

	GMRFLib_ASSERT(n > 0, GMRFLib_EINVARG);

	if (first) {
		first = 0;
		map_ivp_init(&abscissas);		       /* init the hash-table */
		map_ivp_init(&weights);			       /* init the hash-table */
	}

	if ((ptr = map_ivp_ptr(&abscissas, n))) {
		/*
		 * use previously computed. note that map_ivp_ptr returns a ptr to the stored ptr. 
		 */
		if (xp) {
			*xp = *((double **) ptr);
		}

		/*
		 * the weights should now be stored as well 
		 */
		if (wp) {
			pptr = map_ivp_ptr(&weights, n);
			GMRFLib_ASSERT(pptr, GMRFLib_ESNH);
			*wp = *((double **) pptr);
		}
	} else {
		/*
		 * compute new ones 
		 */

		x = Calloc(n, double);
		w = Calloc(n, double);

		GMRFLib_EWRAP0(GMRFLib_ghq__intern(x, w, n));

		/*
		 * the Gauss-Hermite is with kernel exp(-x^2), transform to kernel exp(-x^2/2)/sqrt(2*pi)
		 */
		for (i = 0; i < n; i++) {
			x[i] *= M_SQRT2;
			w[i] *= M_SQRT2 / sqrt(2.0 * M_PI);
		}

		/*
		 * reverse the order so its small to large, and sort the weights along 
		 */
		GMRFLib_qsorts((void *) x, (size_t) n, sizeof(double), w, sizeof(double), NULL, 0, GMRFLib_dcmp);

		map_ivp_set(&abscissas, n, (void *) x);
		map_ivp_set(&weights, n, (void *) w);

		if (xp) {
			*xp = x;			       /* return a ptr only */
		}
		if (wp) {
			*wp = w;			       /* return a ptr only */
		}
	}
	return GMRFLib_SUCCESS;
}
