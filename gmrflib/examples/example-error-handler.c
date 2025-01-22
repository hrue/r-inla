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
