#include <stdio.h>
#include <stdlib.h>
#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

#if 0
int main(int argc, char **argv)
{
	if (argc == 1) {
		fprintf(stderr, "Usage: %s  Qij.file\n", argv[0]);
		exit(1);
	}

	double tref;
	GMRFLib_tabulate_Qfunc_tp *tab = NULL;
	GMRFLib_graph_tp *graph = NULL;

	tref = GMRFLib_timer();
	fprintf(stderr, "Read Q-matrix and build graph... ");
	fflush(stderr);
	GMRFLib_tabulate_Qfunc_from_file(&tab, &graph, argv[1], -1, NULL);
	fprintf(stderr, "%gs\n", GMRFLib_timer() - tref);

	fprintf(stderr, "Graph-size %d\n", graph->n);

	tref = GMRFLib_timer();
	fprintf(stderr, "Find a good reordering... ");
	fflush(stderr);
	GMRFLib_optimize_reorder(graph, NULL);
	fprintf(stderr, "[%s] ... %gs\n", GMRFLib_reorder_name(GMRFLib_reorder), GMRFLib_timer() - tref);

	tref = GMRFLib_timer();
	fprintf(stderr, "Factorise Q... ");
	fflush(stderr);
	GMRFLib_problem_tp *problem = NULL;
	GMRFLib_init_problem(&problem, NULL, NULL, NULL, NULL, graph, tab->Qfunc, tab->Qfunc_arg, NULL, NULL, GMRFLib_NEW_PROBLEM);
	fprintf(stderr, "%gs\n", GMRFLib_timer() - tref);

	tref = GMRFLib_timer();
	fprintf(stderr, "Compute Qinv... ");
	fflush(stderr);
	GMRFLib_Qinv(problem, GMRFLib_QINV_ALL);
	fprintf(stderr, "%gs\n", GMRFLib_timer() - tref);

	int i, j, jj;
	for (i = 0; i < graph->n; i++) {
		printf("%d %d %.8g\n", i, i, *GMRFLib_Qinv_get(problem, i, i));
		for (jj = 0; jj < graph->lnnbs[i]; jj++) {
			j = graph->lnbs[i][jj];
			printf("%d %d %.8g\n", i, j, *GMRFLib_Qinv_get(problem, i, j));
		}
	}

	return 0;
}
#endif
