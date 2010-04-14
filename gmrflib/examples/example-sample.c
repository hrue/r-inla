
/* example-sample.c
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
#include <math.h>
#include <strings.h>
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>

static const char RCSId[] = "$Id: example-sample.c,v 1.16 2008/10/29 16:32:25 hrue Exp $";
#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"


/* The Q-function, returning  
   return Q(i,j), i=node, j=nnode, corresponding to 
   pi(x)  \exp(-\frac{1}{2} \sum_{i\sim j} (x_i-x_j)^2)
*/

double Qfunc(int node, int nnode, void *arg);

double Qfunc(int node, int nnode, void *arg)
{
	GMRFLib_graph_tp *g = (GMRFLib_graph_tp *) arg;

	if (node != nnode) {
		return -1.0;
	} else {
		return g->nnbs[node];
	}
}
int main(int argc, char **argv)
{
	GMRFLib_graph_tp *graph;
	GMRFLib_problem_tp *problem;
	GMRFLib_constr_tp *constr = NULL;

	double *x, *b, *c;
	int i, method, seed;
	void *argQfunc;
	char *fixed = NULL;
	FILE *fp;

	const int nrow = 6, ncol = 6;			       /* Size of lattice graph */

	if (argc < 3) {
		printf("Usage: %s SEED METHOD\n", argv[0]);
		exit(-1);
	}
	seed = atoi(argv[1]);
	method = atoi(argv[2]);

	/*
	 * Initialize random number generator: 
	 */
	GMRFLib_uniform_init((unsigned long int) seed);

	/*
	 * create a lattice graph 
	 */
	GMRFLib_make_lattice_graph(&graph, nrow, ncol, 1, 1, 0);

	x = Calloc(graph->n, double);
	c = Calloc(graph->n, double);
	b = Calloc(graph->n, double);

	/*
	 * Unconditional sampling: 
	 */
	if (method == 1) {
		for (i = 0; i < graph->n; i++) {
			x[i] = 0.0;
			c[i] = 1.0;
		}
		argQfunc = (void *) graph;
		GMRFLib_init_problem(&problem, x, NULL, c, NULL, graph, Qfunc, argQfunc, NULL, NULL, 0);
	}

	/*
	 * Conditional sampling, conditional on fixed values: 
	 */
	if (method == 2) {
		fixed = Calloc(graph->n, char);

		for (i = 0; i < graph->n; i = i + 2) {
			fixed[i] = 1;
		}

		for (i = 0; i < graph->n; i++) {
			x[i] = GMRFLib_uniform();
			c[i] = 0.0;
		}
		argQfunc = (void *) graph;
		GMRFLib_init_problem(&problem, x, NULL, c, NULL, graph, Qfunc, argQfunc, fixed, NULL, 0);
	}

	/*
	 * Conditional sampling, conditional on a linear deterministic constraint: 
	 */
	if (method == 3) {
		fixed = Calloc(graph->n, char);

		for (i = 0; i < graph->n; i = i + 2) {
			fixed[i] = 1;
		}
		for (i = 0; i < graph->n; i++) {
			x[i] = GMRFLib_uniform();
			c[i] = 1.0;
			b[i] = GMRFLib_uniform();
		}

		/*
		 * Create constraint: 
		 */
		GMRFLib_make_empty_constr(&constr);
		constr->nc = 2;
		constr->a_matrix = Calloc(constr->nc * graph->n, double);
		constr->e_vector = Calloc(constr->nc, double);

		constr->e_vector[0] = 0.0;
		constr->e_vector[1] = 0.0;

		constr->a_matrix[1] = 1.0;
		constr->a_matrix[3] = 0.0;
		constr->a_matrix[5] = 2.0;
		for (i = 0; i < graph->n; i++) {
			constr->a_matrix[i * constr->nc] = 1.0;
			if (i > 2) {
				constr->a_matrix[i * constr->nc + 1] = 0.0;
			}
		}
		GMRFLib_prepare_constr(constr, graph, 0);

		argQfunc = (void *) graph;
		GMRFLib_init_problem(&problem, x, b, c, NULL, graph, Qfunc, argQfunc, fixed, constr, 0);
	}

	/*
	 * Conditional sampling, conditional on a linear stochastic constraint: 
	 */

	if (method == 4) {
		fixed = Calloc(graph->n, char);

		for (i = 0; i < graph->n; i++)
			fixed[i] = (GMRFLib_uniform() < 0.5 ? 0 : 1);
		for (i = 0; i < graph->n; i++) {
			x[i] = GMRFLib_uniform();
			c[i] = 1.0;
			b[i] = GMRFLib_uniform();
		}

		/*
		 * Create constraint: 
		 */
		GMRFLib_make_empty_constr(&constr);
		constr->nc = 1;
		constr->a_matrix = Calloc(graph->n, double);

		for (i = 0; i < graph->n; i++) {
			constr->a_matrix[i] = 1.0;
		}
		constr->e_vector = Calloc(1, double);

		constr->e_vector[0] = 0.0;
		constr->errcov_diagonal = Calloc(constr->nc, double);

		constr->errcov_diagonal[0] = 1.0;
		GMRFLib_prepare_constr(constr, graph, 0);

		argQfunc = (void *) graph;
		GMRFLib_init_problem(&problem, x, b, c, NULL, graph, Qfunc, argQfunc, fixed, constr, 0);
	}

	/*
	 * Sample: 
	 */
	GMRFLib_sample(problem);

	/*
	 * Extracting and printing the sample, 
	 */
	/*
	 * the log-density and the mean: 
	 */

	fp = fopen("samples.dat", "w");
	for (i = 0; i < graph->n; i++) {
		fprintf(fp, "Sample %d: %12.6f\n", i, problem->sample[i]);
	}
	fclose(fp);

	printf("Log-density of sample: %.12f\n\n", problem->sub_logdens);
	printf("Sub-index mean (unconstr) mean (constr):\n");
	for (i = 0; i < problem->sub_graph->n; i++)
		printf("%d\t%12.6f\t%12.6f\n", problem->map[i], problem->sub_mean[i], problem->sub_mean_constr[i]);

	printf("Mean (unconstr) mean (constr):\n");
	for (i = 0; i < graph->n; i++) {
		printf("%d\t%12.6f\t%12.6f\n", i, problem->mean[i], problem->mean_constr[i]);
	}

	/*
	 * Free allocated memory: 
	 */
	Free(fixed);
	Free(x);
	Free(c);
	Free(b);

	GMRFLib_free_constr(constr);
	GMRFLib_free_problem(problem);
	GMRFLib_free_graph(graph);

	return 0;
}
