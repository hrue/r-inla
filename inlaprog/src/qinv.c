#include <stdio.h>
#include <stdlib.h>
#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

int main(int argc, char **argv)
{
	if (argc == 1) {
		fprintf(stderr, "Usage: %s  Qij.file\n", argv[0]);
		exit(1);
	}

	double tref;
	GMRFLib_tabulate_Qfunc_tp *tab;
	GMRFLib_graph_tp *graph;

	tref = GMRFLib_cpu();
	fprintf(stderr, "Read Q-matrix and build graph... ");
	fflush(stderr);
	GMRFLib_tabulate_Qfunc_from_file(&tab, &graph, argv[1], NULL, NULL, NULL);
	fprintf(stderr, "%gs\n", GMRFLib_cpu() - tref);

	fprintf(stderr, "Graph-size %d\n", graph->n);

	tref = GMRFLib_cpu();
	fprintf(stderr, "Find a good reordering... ");
	fflush(stderr);
	GMRFLib_optimize_reorder(graph, NULL);
	fprintf(stderr, "[%s] ... %gs\n", GMRFLib_reorder_name(GMRFLib_reorder), GMRFLib_cpu() - tref);

	tref = GMRFLib_cpu();
	fprintf(stderr, "Factorise Q... ");
	fflush(stderr);
	GMRFLib_problem_tp *problem;
	GMRFLib_init_problem(&problem, NULL, NULL, NULL, NULL, graph, tab->Qfunc, tab->Qfunc_arg, NULL, NULL, GMRFLib_NEW_PROBLEM);
	fprintf(stderr, "%gs\n", GMRFLib_cpu() - tref);

	tref = GMRFLib_cpu();
	fprintf(stderr, "Compute Qinv... ");
	fflush(stderr);
	GMRFLib_Qinv(problem, GMRFLib_QINV_ALL);
	fprintf(stderr, "%gs\n", GMRFLib_cpu() - tref);

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
