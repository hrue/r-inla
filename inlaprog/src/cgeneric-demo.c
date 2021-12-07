
/* cgeneric-demo.c
 * 
 * Copyright (C) 2021 Havard Rue
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


#include <assert.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <math.h>
#include <stdlib.h>

#include "cgeneric.h"

#if !defined(Calloc)
#define Calloc(n_, type_)  (type_ *)calloc((n_), sizeof(type_))
#endif


#define N 10
double *inla_cgeneric_demo(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp * data)
{
	// this implements a simple IID model for testing purposes.
	// for sparse matrices, return in format (n, len, i, j, Qij), where i<=j.
	// for the graph, then Qij is known to be 1, so its not needed.

	double *ret = NULL, prec = (theta ? exp(theta[0]) : NAN), lprec = (theta ? theta[0] : NAN);

	switch (cmd) {
	case INLA_CGENERIC_VOID:
	{
		assert(!(cmd == INLA_CGENERIC_VOID));
		break;
	}

	case INLA_CGENERIC_GRAPH:
	{
		ret = Calloc(2 + 2 * N, double);
		ret[0] = N;				       /* dimension */
		ret[1] = N;				       /* number of (i <= j) */
		for (int i = 0; i < N; i++) {
			ret[2 + i] = i;
			ret[2 + N + i] = i;
		}
		break;
	}

	case INLA_CGENERIC_Q:
	{
		// optimized format only
		ret = Calloc(2 + N, double);
		ret[0] = -1;				       /* code for optimized output */
		ret[1] = N;				       /* number of (i <= j) */
		for (int i = 0; i < N; i++) {
			ret[2 + i] = prec;
		}
		break;
	}

	case INLA_CGENERIC_MU:
	{
		ret = Calloc(1, double);
		ret[0] = 0;
		break;
	}

	case INLA_CGENERIC_INITIAL:
	{
		ret = Calloc(2, double);
		ret[0] = 1;
		ret[1] = 4.0;
		break;
	}

	case INLA_CGENERIC_LOG_NORM_CONST:
	{
		ret = Calloc(1, double);
		ret[0] = N * (-0.9189385332 + 0.5 * lprec);
		break;
	}

	case INLA_CGENERIC_LOG_PRIOR:
	{
		double u[] = { 1.0, 0.01 }, th, x;
		ret = Calloc(1, double);
		th = -log(u[1]) / u[0];
		x = lprec / 2.0;
		ret[0] = log(th / 2.0) - th * exp(-x) - x;
		break;
	}

	case INLA_CGENERIC_QUIT:
	default:
		break;
	}

	return (ret);
}
