#include <stddef.h>
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"
static const char RCSId[] = "$Id: example-matern2d.c,v 1.3 2008/10/26 03:46:18 hrue Exp $";

int main(int argc, char **argv, char **env)
{
	int i;
	FILE *fp;
	double prec = 1.0, range = 30.0;
	GMRFLib_graph_tp *g;
	GMRFLib_matern2ddef_tp *def;

	def = Calloc(1, GMRFLib_matern2ddef_tp);

	def->nrow = 100;
	def->ncol = 100;
	def->nu = 1;
	def->cyclic = 0;
	def->prec = &prec;
	def->range = &range;

	GMRFLib_make_matern2d_graph(&g, def);

	GMRFLib_problem_tp *problem;
	GMRFLib_init_problem(&problem, NULL, NULL, NULL, NULL, g, GMRFLib_matern2d, (void *) def, NULL, NULL, GMRFLib_NEW_PROBLEM);
	GMRFLib_sample(problem);
	GMRFLib_Qinv(problem, GMRFLib_QINV_ALL);

	printf("Write a sample to matern2d-sample.dat\n");
	fp = fopen("matern2d-sample.dat", "w");
	for (i = 0; i < g->n; i++) {
		fprintf(fp, "%g\n", problem->sample[i]);
	}
	fclose(fp);

	printf("Write marginal variances to matern2d-variances.dat\n");
	fp = fopen("matern2d-variance.dat", "w");
	for (i = 0; i < g->n; i++) {
		fprintf(fp, "%g\n", *GMRFLib_Qinv_get(problem, i, i));
	}
	fclose(fp);

	// GMRFLib_timer_full_report(NULL);

	return 0;
}
