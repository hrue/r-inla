
/* example-graph-edit.c
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

static const char RCSId[] = "$Id: example-graph-edit.c,v 1.6 2007/01/20 13:06:02 hrue Exp $";

#include "GMRFLib/GMRFLib.h"

int main(int argc, char **argv)
{
	/*
	 * in this example we will create the graph required for the RW1 model with random effects and a global mean:
	 * 
	 * \mu ~ N(,) x | \mu ~ \mu + RW1(\kappa) z_i | x ~ N(x_i, \tau)
	 * 
	 * so the GMRF is (x,z,\mu), with total length 2*n+1. Here, x has length n, z has length n and \mu has length 1.
	 * 
	 */
	int i;
	GMRFLib_rwdef_tp rwdef;
	GMRFLib_graph_tp *gg, *g;
	GMRFLib_ged_tp *ged;

	rwdef.n = 100;
	rwdef.order = 1;
	rwdef.cyclic = GMRFLib_FALSE;

	/*
	 * create the RW1-graph 
	 */
	GMRFLib_make_rw_graph(&g, &rwdef);

	/*
	 * and create the editable graph-object 
	 */

	GMRFLib_ged_init(&ged, g);

	/*
	 * add the z's: we just add the edge between x_i and z_i, since then the z_i node will be created automatically. 
	 */
	for (i = 0; i < rwdef.n; i++)
		GMRFLib_ged_add(ged, i, i + rwdef.n);

	/*
	 * this just a check that we get it right. it's 2*n nodes in the graph, so the highest number is 2*n-1. 
	 */
	assert(2 * rwdef.n - 1 == GMRFLib_ged_max_node(ged));

	/*
	 * now we add the mean using the same trick as above by just adding the edges between x_i and \mu 
	 */
	for (i = 0; i < rwdef.n; i++)
		GMRFLib_ged_add(ged, i, 2 * rwdef.n);

	/*
	 * done! build the graph, and write it to file 
	 */
	GMRFLib_ged_build(&gg, ged);
	char *fnm = Strdup("graph-edit-example.graph");

	GMRFLib_write_graph(fnm, gg);

	/*
	 * cleanup & exit 
	 */
	GMRFLib_ged_free(ged);
	GMRFLib_free_graph(g);
	GMRFLib_free_graph(gg);
	free(fnm);

	return 0;
}
