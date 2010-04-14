
/* example-qinv.c
 * 
 * Copyright (C) 2005 Havard Rue
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
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA.
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

#include <stddef.h>
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>

static const char RCSId[] = "$Id: example-qinv.c,v 1.9 2008/10/29 16:32:25 hrue Exp $";

#include "GMRFLib/GMRFLib.h"

typedef struct {
	double phi;					       /* \phi parametere in AR(1) */
	int n;						       /* size of the graph */
} Global_tp;

Global_tp G = { 0.9, 10 };

double Qfunc(int node, int nnode, void *arg)
{
	/*
	 * AR(1) 
	 */
	double scale = 1. / (1 - SQR(G.phi));

	if (node != nnode) {
		return -G.phi * scale;
	}
	if (node == 0 || node == G.n - 1) {
		return scale;
	}
	return (1.0 + SQR(G.phi)) * scale;
}
int main(int argc, char **argv)
{
	/*
	 * the task is compute some elms of the covariance matrix the AR(1) process in the stationary case where we know the truth. 
	 */
	int i, j;
	double *ptr;
	GMRFLib_graph_tp *graph;
	GMRFLib_problem_tp *problem;

	GMRFLib_make_linear_graph(&graph, G.n, 1, 0);
	GMRFLib_init_problem(&problem, NULL, NULL, NULL, NULL, graph, Qfunc, NULL, NULL, NULL, GMRFLib_NEW_PROBLEM);
	GMRFLib_Qinv(problem, GMRFLib_QINV_ALL);

	for (i = 0; i < G.n; i++) {
		if ((ptr = GMRFLib_Qinv_get(problem, i, i))) {
			printf("Var(x[%1d]) = %.12f\n", i, *ptr);
		}
	}
	printf("The covariance matrix\n");
	for (i = 0; i < G.n; i++) {
		for (j = 0; j < G.n; j++) {
			if ((ptr = GMRFLib_Qinv_get(problem, i, j))) {
				printf(" %.12f", *ptr);
			} else {
				printf(" %.12s", "NA");
			}
		}
		printf("\n");
	}

	for (i = 0; i < G.n; i++) {
		if ((ptr = GMRFLib_Qinv_get(problem, i, i))) {
			printf("Stdev(x[%1d]) = %.12f\n", i, sqrt(*ptr));
		}
	}
	return 0;
}
