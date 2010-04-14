
/* example-graph1.c
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

#include <stddef.h>
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>

static const char RCSId[] = "$Id: example-graph1.c,v 1.10 2008/10/29 16:26:36 hrue Exp $";

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"


int main(int argc, char **argv)
{
	GMRFLib_graph_tp *graph;			       /* Graph object */
	GMRFLib_graph_tp *sgraph;			       /* Subgraph */
	char gfile[] = "graph1.dat";			       /* Input file */
	char *my_remove;
	int i;

	/*
	 * Read graph from the file gfile, 
	 */
	/*
	 * and print the graph specification to standard output 
	 */
	printf("\nReading a graph from the file '%s':\n\n", gfile);

	GMRFLib_read_graph(&graph, gfile);
	GMRFLib_print_graph(stdout, graph);

	/*
	 * Indicate nodes that are to be removed 
	 */
	my_remove = Calloc(graph->n, char);

	for (i = 0; i < graph->n / 2; i++) {
		my_remove[i] = 1;
	}

	/*
	 * Compute and print subgraph 
	 */
	printf("\nComputing subgraph:\n\n");
	GMRFLib_compute_subgraph(&sgraph, graph, my_remove);
	GMRFLib_print_graph(stdout, sgraph);

	GMRFLib_free_graph(graph);
	GMRFLib_free_graph(sgraph);
	free(my_remove);

	return 0;
}
