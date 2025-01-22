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
