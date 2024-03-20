
/* 
 * 
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

/* 
   include the required header-file
*/
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
#define PREC_PHI_TRANS 0.01

static const char RCSId[] = "$Id: example-approx-2.c,v 1.16 2010/03/12 12:24:20 hrue Exp $";

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
	 *  return the normalized log density for the  precision kappa
	 */
	return (log(gsl_ran_gamma_pdf(x, a, 1 / b)));
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

	kappa = exp(theta[0]);
	nu = theta[1];

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
	char *compute;
	int n_hyper = 2;				       /* number of hyperparameters */
	double **hyper[n_hyper];			       /* vector of pointers to hyperparameters */
	double **log_kappa, **nu;

	GMRFLib_ai_param_tp *ai_par;
	GMRFLib_ai_store_tp *ai_store;
	GMRFLib_ai_marginal_likelihood_tp *marginal_likelihood;
	GMRFLib_density_tp **density, **gdensity;

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
	 * Compute approximate posterior density for the chosen nodes in the hidden field, hypeparameters 
	 * integrated out
	 */

	/*
	 * allocate memory for the densities 
	 */
	density = Calloc(G.graph->n, GMRFLib_density_tp *);
	gdensity = Calloc(G.graph->n, GMRFLib_density_tp *);

	/*
	 * allocate memory for the marginal_likelihood 
	 */

	marginal_likelihood = Calloc(1, GMRFLib_ai_marginal_likelihood_tp);

	/*
	 *allocate memory for ai_store
	 */
	ai_store = Calloc(1, GMRFLib_ai_store_tp);

	/*
	 * get default options for the approximation-routine
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

	Qarg.log_kappa = log_kappa;
	Qarg.nu = nu;
	/*
	 * set the vector of pointers to the hyperparameters 
	 */
	hyper[0] = log_kappa;
	hyper[1] = nu;

	GMRFLib_ai_INLA(&density, &gdensity, NULL, NULL, NULL, marginal_likelihood, NULL,
			compute,
			hyper, n_hyper, extra, NULL, NULL, NULL, NULL, NULL, G.d, loglik, NULL, NULL, G.graph, Qfunc, (void *) &Qarg, NULL, ai_par,
			ai_store, NULL, NULL, 0, NULL, NULL);
	/*
	 * print the results 
	 */
	FILE *fp1 = NULL, *fp2 = NULL, *fp3 = NULL;

	fp1 = fopen("lap_results.dat", "w");
	fp2 = fopen("gaus_results.dat", "w");
	fp3 = fopen("kl_distance.dat", "w");

	/*
	 *NB: the densities stored in density[i] are in standardised scale so it is
	 *necessary a transformation to get them in user scale
	 */
	for (i = 0; i < G.graph->n; i++) {
		/*
		 * print the simplified laplace approximation 
		 */
		if (density[i]) {

			fprintf(fp1, "%d ", i);
			double xx, x_real, f_real, f;

			for (xx = -4.0; xx < 4.0; xx += 0.1) {
				GMRFLib_evaluate_density(&f, xx, density[i]);
				x_real = GMRFLib_density_std2user(xx, density[i]);
				f_real = f / density[i]->std_stdev;
				fprintf(fp1, " %.10f %.10f ", x_real, f_real);
			}
			fprintf(fp1, "\n");
		}
		/*
		 * print the Gaussian approximation 
		 */
		if (gdensity[i]) {

			fprintf(fp2, "%d ", i);
			double xx, x_real, f_real, f;

			for (xx = -4.0; xx < 4.0; xx += 0.1) {
				GMRFLib_evaluate_density(&f, xx, gdensity[i]);
				x_real = GMRFLib_density_std2user(xx, gdensity[i]);
				f_real = f / gdensity[i]->std_stdev;
				fprintf(fp2, " %.10f %.10f ", x_real, f_real);
			}
			fprintf(fp2, "\n");
		}
		/*
		 *print the Kullback-Leibler distance between the Gaussian approximation
		 * and the simplified Laplace one.
		 */
		if (density[i] && gdensity[i]) {
			double kld;

			GMRFLib_kld_sym(&kld, density[i], gdensity[i]);
			fprintf(fp3, "%d %.10f\n", i, kld);

			GMRFLib_free_density(density[i]);
			GMRFLib_free_density(gdensity[i]);
		}
	}
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);

	/*
	 * free memory
	 */
	free(ai_par);
	free(ai_store);
	free(marginal_likelihood);
	return 0;
}
