
/* example-bitmap.c
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
