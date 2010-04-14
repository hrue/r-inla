
/* example-graph2.c
 * 
 * Copyright (C) 2001 Havard Rue
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

#include <stdio.h>

#include "GMRFLib/GMRFLib.h"

int main(int argc, char **argv)
{
	GMRFLib_graph_tp *lgraph;			       /* Graph object, lattice graph */
	GMRFLib_graph_tp *argraph;			       /* Graph object, linear graph */
	GMRFLib_graph_tp *ngraph;			       /* Graph object, folded graph */
	int i, j, nindx, li, lj;
	const int nrow = 6, ncol = 7;			       /* Size of lattice */

	/*
	 * Create a non-cyclic graph on a 6 by 7 lattice, 
	 */
	/*
	 * using a 3 by 3 neighbourhood 
	 */
	printf("\nCreating a lattice graph:\n\n");

	GMRFLib_make_lattice_graph(&lgraph, nrow, ncol, 1, 1, 0);
	GMRFLib_print_graph(stdout, lgraph);

	/*
	 * Extract the ordering of the nodes of the lattice, 
	 */
	/*
	 * and recompute the lattice indices from the node indices 
	 */

	printf("\nExtracting node indices:\n\n");
	for (i = 0; i < nrow; i++) {
		for (j = 0; j < ncol; j++) {
			GMRFLib_lattice2node(&nindx, i, j, nrow, ncol);
			printf("Lattice indices (%d,%d) correspond to node no. %d\n", i, j, nindx);
		}
	}
	printf("\nExtracting lattice indices:\n\n");
	for (i = 0; i < lgraph->n; i++) {
		GMRFLib_node2lattice(i, &li, &lj, nrow, ncol);
		printf("Node no. %d corresponds to lattice indices (%d,%d)\n", i, li, lj);
	}

	/*
	 * Create a cyclic linear graph on an AR(1)-model 
	 */
	printf("\nCreating a linear graph:\n\n");

	GMRFLib_make_linear_graph(&argraph, 10, 1, 1);
	GMRFLib_print_graph(stdout, argraph);

	/*
	 * Expand the neighbourhood of linear graph twice 
	 */
	printf("\nCreating a folded graph:\n\n");
	GMRFLib_nfold_graph(&ngraph, argraph, 2);
	GMRFLib_print_graph(stdout, ngraph);

	GMRFLib_free_graph(lgraph);
	GMRFLib_free_graph(argraph);
	GMRFLib_free_graph(ngraph);

	return 0;
}
