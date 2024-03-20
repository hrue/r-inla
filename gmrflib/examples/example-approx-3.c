
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

static const char RCSId[] = "$Id: example-approx-3.c,v 1.9 2008/10/29 17:11:05 hrue Exp $";

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
   parameters in the gamma-prior for kappa
 */
#define PRIOR_A 0.25
#define PRIOR_B 0.02

/*
  prior precision for nu
*/
#define PREC_PHI_TRANS 0.1

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
	double **log_kappa;
	double **nu;
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
	kappa = exp(a->log_kappa[GMRFLib_thread_id][0]);
	phi = trans_phi(a->nu[GMRFLib_thread_id][0]);

	if (node == nnode) {
		if (node == 0 || node == G.graph->n)
			return kappa;
		else
			return (1 + pow(phi, 2)) * kappa;
	} else
		return -(phi * kappa);
}

double log_prior(double x, double a, double b)
{
	/*
	 * return the unnormalized log-density of a gamma-variable with mean a/b. 
	 */
	return (a - 1.0) * log(x) - x * b;
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

double extra(double *theta, int ntheta, void *args)
{
	/*
	 *returns all terms in GMRF-35 which are constant with respect to x
	 * but depends on the hyperparameters
	 */
	double nu, kappa, val;

	nu = theta[0];
	kappa = exp(theta[1]);

	val = 0.5 * G.n * log(kappa);
	val += log_prior(kappa, PRIOR_A, PRIOR_B);
	val += log(kappa);
	val += -0.5 * PREC_PHI_TRANS * nu * nu;
	return (val);
}

int main(int argc, char **argv, char **env)
{
	Qfunc_args_tp Qarg;
	int i;

	double **log_kappa, **nu;
	double logdens;

	GMRFLib_ai_param_tp *ai_par;
	GMRFLib_ai_store_tp *ai_store;

	/*
	 * get graph and data 
	 */
	read_data_make_graph();

	/*
	 * get default options for the approximation-routine
	 * to compute the marginal for the hyperparameter the default
	 * options are fine 
	 */
	GMRFLib_default_ai_param(&ai_par);

	/*
	 * set initial value for the hyperparameters and setup the argument for the Q function 
	 */
	int tmax = omp_get_max_threads();

	log_kappa = Calloc(tmax, double *);
	nu = Calloc(tmax, double *);
	for (i = 0; i < tmax; i++) {
		log_kappa[i] = Calloc(1, double);
		log_kappa[i][0] = 10.0;
		nu[i] = Calloc(1, double);
		nu[i][0] = 2.0;
	}

	Qarg.nu = nu;
	Qarg.log_kappa = log_kappa;

	double lkappa, nnu;

	for (lkappa = 3; lkappa < 7; lkappa += 0.1)
		for (nnu = 3.5; nnu < 10; nnu += 0.2) {
			ai_store = Calloc(1, GMRFLib_ai_store_tp);

			for (i = 0; i < tmax; i++) {
				nu[i][0] = nnu;
				log_kappa[i][0] = lkappa;
			}
			GMRFLib_ai_marginal_hyperparam(&logdens,
						       NULL, NULL, NULL, NULL, G.d, loglik, NULL, NULL, G.graph, Qfunc, (void *) &Qarg, NULL,
						       ai_par, ai_store);
			double theta[2];
			theta[0] = nnu;
			theta[1] = lkappa;

			logdens += extra(theta, 2, NULL);
			printf("%f %f %f \n", lkappa, nnu, logdens);
			/*
			 * note that it is necessary here to free
			 *  ai_store and reallocate new memory at every
			 * iteration of the loop.
			 */
			free(ai_store);
		}
	return 0;
}
