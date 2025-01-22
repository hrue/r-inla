#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>

static const char RCSId[] = "$Id: example-bitmap.c,v 1.5 2007/02/02 19:20:10 hrue Exp $";

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

int main(int argc, char **argv)
{
	/*
	 * in this example we will create bitmaps of the precision matrix, its ordered version and the Cholelsky triangle of a CRW2 model. here 
	 * we create a simple problem to use the GMRFLib_bitmap_problem() function. 
	 */

	int i;
	double prec = 1.0, *c;
	GMRFLib_crwdef_tp *crwdef;
	GMRFLib_problem_tp *problem;
	GMRFLib_graph_tp *graph;

	crwdef = Calloc(1, GMRFLib_crwdef_tp);
	crwdef->n = 50;
	crwdef->order = 2;
	crwdef->prec = &prec;
	crwdef->position = NULL;
	crwdef->layout = GMRFLib_CRW_LAYOUT_PAIRS;

	c = Calloc(2 * crwdef->n, double);		       /* just to make the GMRF proper */

	for (i = 0; i < 2 * crwdef->n; i++) {
		c[i] = 1.0;
	}

	GMRFLib_make_crw_graph(&graph, crwdef);
	GMRFLib_init_problem(&problem, NULL, NULL, c, NULL, graph, GMRFLib_crw, (void *) crwdef, NULL, NULL, GMRFLib_NEW_PROBLEM);

	char *fnm = Strdup("example");

	GMRFLib_bitmap_problem(fnm, problem);

	GMRFLib_free_problem(problem);
	GMRFLib_free_graph(graph);
	Free(crwdef->work);
	Free(crwdef);
	Free(fnm);

	return 0;
}
