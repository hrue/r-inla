
/* 
 * Copyright (C) 2007 Havard Rue
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
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif

static const char RCSId[] = "$Id: example-approx-1.c,v 1.8 2008/10/29 17:00:34 hrue Exp $";

#include "GMRFLib/GMRFLib.h"

/* 
   the name of the data file
*/
#define DATAFILE "poisson.data"

/*
  length of the data
*/
#define LEN_DATA   300

/* 
   global variable G containg the data and their weights and the graph. 
*/
typedef struct {
	int n;
	double *y;
	double *d;
	GMRFLib_graph_tp *graph;
} Global_tp;
static Global_tp G;

/* 
   structure containing the argumets for the Q function
*/
typedef struct {
	double *log_kappa;
	double *nu;
} Qfunc_args_tp;

double trans_phi(double x)
{
	/*
	 * function transforming nu into phi
	 */
	return (2 * exp(x) / (1 + exp(x)) - 1);
}

int read_data_make_graph(void)
{
	/*
	 * read the data file and make the graph
	 */
	int i;
	FILE *fp;

	G.n = LEN_DATA;

	GMRFLib_make_linear_graph(&G.graph, G.n, 1, 0);	       /* RW of order 1, not cyclic */

	fp = fopen(DATAFILE, "r");
	assert(fp);
	G.y = Calloc(G.n, double);

	for (i = 0; i < G.n; i++)
		fscanf(fp, "%lf\n", &G.y[i]);
	fclose(fp);

	G.d = Calloc(G.graph->n, double);

	for (i = 0; i < (G.graph->n - 1); i++)
		G.d[i] = 1.0;				       /* weights for the data */
	return 0;
}

double Qfunc(int node, int nnode, void *arg)
{
	/*
	 * return Q_ij, i=node, j=nnode, for the prior of the latent field
	 * 
	 * AR1 process
	 */
	Qfunc_args_tp *a;
	double kappa, phi;

	a = (Qfunc_args_tp *) arg;
	kappa = exp(*a->log_kappa);
	phi = trans_phi(*a->nu);

	if (node == nnode) {
		if (node == 0 || node == G.graph->n)
			return kappa;
		else
			return (1 + pow(phi, 2)) * kappa;
	} else
		return -(phi * kappa);
}

int loglik(double *logll, double *x, int m, int idx, double *x_vec, void *arg)
{
	/*
	 * return the log-likelihood for a poisson, y_idx|... ~ Po(E_idx \exp(x_idx))
	 * and its first two dirivatives
	 */
	if (m == 0)
		return GMRFLib_LOGL_COMPUTE_DERIVATIES;	       /* this line is needed since we want exact derivatives */
	logll[0] = G.y[idx] * x[0] - exp(x[0]);
	if (m > 1)
		logll[1] = G.y[idx] - exp(x[1]);
	if (m > 2)
		logll[2] = -exp(x[2]);
	return GMRFLib_SUCCESS;
}

int main(int argc, char **argv, char **env)
{
	Qfunc_args_tp Qarg;
	int i;
	char *compute;
	double log_kappa, nu;

	GMRFLib_ai_param_tp *ai_par;
	GMRFLib_ai_store_tp *ai_store;
	GMRFLib_density_tp **density;

	/*
	 * get graph and data 
	 */
	read_data_make_graph();

	/*
	 *choose the nodes in the hidden field you want for which to compute the density.
	 * Here we chosed the first 10 nodes
	 */
	compute = Calloc(G.graph->n, char);

	for (i = 0; i < 10; i++)
		compute[i] = 1.0;

	/*
	 * Compute approximate posterior density for the chosen nodes in the hidden field for a
	 * value of the hyperparameters kappa and phi
	 */

	/*
	 * allocate space for the densities 
	 */
	density = Calloc(G.graph->n, GMRFLib_density_tp *);
	/*
	 *allocate space for ai_store
	 */
	ai_store = Calloc(1, GMRFLib_ai_store_tp);

	/*
	 * get default options for the approximation-routine
	 * now it's easy to change options if required. 
	 */
	GMRFLib_default_ai_param(&ai_par);
	/*
	 *for example choose a Gaussian strategy instead of the default
	 *simplified Laplace
	 */
	ai_par->strategy = GMRFLib_AI_STRATEGY_GAUSSIAN;

	/*
	 * choose a value for the hyperparameters and setup the argument for the Q function 
	 */
	log_kappa = 10.0;
	nu = 2.0;

	Qarg.log_kappa = &log_kappa;
	Qarg.nu = &nu;

	/*
	 * compute densities for the chosen nodes
	 */

	for (i = 0; i < G.graph->n; i++)
		if (compute[i]) {
			GMRFLib_ai_marginal_hidden(&density[i], NULL,
						   i, NULL, NULL, NULL, NULL, G.d, loglik, NULL, NULL, G.graph, Qfunc, (void *) &Qarg, NULL, ai_par,
						   ai_store);
		}
	/*
	 * print the results 
	 */
	FILE *fp = NULL, *fpp = NULL;

	fp = fopen("lap_results_fix.dat", "w");
	fpp = fopen("gaus_results_fix.dat", "w");
	/*
	 *NB: the densities stored in density[i] are in standardised scale so it is
	 *necessary a transformation to get them in user scale
	 */
	for (i = 0; i < G.graph->n; i++) {
		if (compute[i]) {
			fprintf(fp, "%d ", i);
			double xx, x_real, f_real, f;

			for (xx = -4.0; xx < 4.0; xx += 0.1) {
				GMRFLib_evaluate_density(&f, xx, density[i]);
				x_real = GMRFLib_density_std2user(xx, density[i]);
				f_real = f / density[i]->std_stdev;
				fprintf(fp, " %.10f %.10f ", x_real, f_real);
			}
			fprintf(fp, "\n");
			fprintf(fpp, " %.10f %.10f\n", density[i]->std_mean, density[i]->std_stdev);
		}
	}
	fclose(fpp);
	fclose(fp);
	/*
	 * free memory
	 */
	free(ai_par);
	free(ai_store);
	free(density);

	return 0;
}
