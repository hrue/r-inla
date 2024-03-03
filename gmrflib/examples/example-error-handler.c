
/* example-error-handler.c
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

#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"

double Qfunc_invalid(int node, int nnode, void *args);

double Qfunc_invalid(int node, int nnode, void *args)
{
	/*
	 * this is an invalid Qfunction..... 
	 */
	return -1;
}
int main(int argc, char **argv)
{
	/*
	 * this example demonstrate the use of the error-handling system. 
	 */
	int retval;
	char *fnm;
	GMRFLib_error_handler_tp *ehandler;
	GMRFLib_graph_tp *graph;
	GMRFLib_problem_tp *problem;

	GMRFLib_make_linear_graph(&graph, 10, 1, 0);

	/*
	 * turn the error-handling system off so that the error-code it just returned 
	 */
	ehandler = GMRFLib_set_error_handler_off();

	retval = GMRFLib_init_problem(&problem, NULL, NULL, NULL, NULL, graph, Qfunc_invalid, NULL, NULL, NULL, GMRFLib_NEW_PROBLEM);
	printf("retval = %d with meaning [%s]\n", retval, GMRFLib_error_reason(retval));

	fnm = Strdup("this graph does not exists.graph");
	retval = GMRFLib_read_graph(&graph, fnm);
	printf("retval = %d with meaning [%s]\n", retval, GMRFLib_error_reason(retval));
	free(fnm);

	fnm = Strdup("graph-with-error.graph");
	retval = GMRFLib_read_graph(&graph, fnm);
	printf("retval = %d with meaning [%s]\n", retval, GMRFLib_error_reason(retval));
	free(fnm);

	/*
	 * switch back to the previous (or default) error-handler 
	 */
	if (1) {
		GMRFLib_set_error_handler(ehandler);	       /* previous */
	}
	if (0) {
		GMRFLib_set_error_handler(NULL);	       /* default */
	}

	/*
	 * just check that this one fails as `usual'... 
	 */
	GMRFLib_init_problem(&problem, NULL, NULL, NULL, NULL, graph, Qfunc_invalid, NULL, NULL, NULL, GMRFLib_NEW_PROBLEM);

	return 0;
}
