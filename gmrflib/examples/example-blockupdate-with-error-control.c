
/* example-blockupdate-with-error-control.c
 * 
 * Copyright (C) 2008 Havard Rue
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

static const char RCSId[] = "$Id: example-blockupdate-with-error-control.c,v 1.8 2010/03/12 12:24:28 hrue Exp $";

/* 
   include the required header-file
 */
#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLib.h"				       /* useful macros */

/* 
   the names of the graph and data
 */
#define GRAPH "germany.graph"
#define DATA  "germany.data"

/* 
   parameters in the gamma-prior for kappa
 */
#define PRIOR_A 0.25
#define PRIOR_B 0.005

/* 
   global variable G containg the data and the graph. 
 */
typedef struct {
	double *y;
	double *E;
	GMRFLib_graph_tp *graph;
} Global;
static Global G;

double eu_gamma(double x, double a, double b);
int read_graph_and_data(void);
int loglik(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
double Qfunc(int node, int nnode, void *arg);

/* 
   This variable determine the sign of the logliklihood contribution, so its 1.0. If we switch the sign, then the blockupdate
   routine fails. We use this to demonstrate how to return the error code instead of calling the default error-handler.
 */
static double loglikelihood_sign = 1.0;

double eu_gamma(double x, double a, double b)
{
	/*
	 * return the unnormalized log-density of a gamma-variable with mean a/b. 
	 */
	return (a - 1.0) * log(x) - x * b;
}
int read_graph_and_data(void)
{
	/*
	 * read the graph and the data 
	 */
	int i;
	double dummy;
	FILE *fp;
	char *fnm;

	fnm = Strdup(GRAPH);
	GMRFLib_read_graph(&G.graph, fnm);
	G.y = Calloc(G.graph->n, double);
	G.E = Calloc(G.graph->n, double);

	fp = fopen(DATA, "r");
	for (i = 0; i < G.graph->n; i++) {
		if (fscanf(fp, "%lf %lf %lf\n", &G.y[i], &G.E[i], &dummy) != 3) {
			abort();
		}
	}
	fclose(fp);
	free(fnm);
	return 0;
}
int loglik(double *logll, double *x, int m, int idx, double *x_vec, void *arg)
{
	/*
	 * return the log-likelihood for a poisson, y_idx|... ~ Po(E_idx \exp(x_idx)) 
	 */
	int i;
	double *y, *E;
	char **args;

	/*
	 * the data y and E come through the arg-pointer 
	 */
	args = (char **) arg;
	y = (double *) args[0];
	E = (double *) args[1];

	/*
	 * compute the log-likelihood 
	 */
	for (i = 0; i < m; i++) {
		logll[i] = loglikelihood_sign * (y[idx] * x[i] - E[idx] * exp(x[i]));
	}

	return GMRFLib_SUCCESS;
}
double Qfunc(int node, int nnode, void *arg)
{
	/*
	 * return Q_ij, i=node, j=nnode, for the prior
	 * 
	 * \pi(x) \sim \exp(-\frac{1}{2} \sum_{i\sim j} (x_i-x_j)^2) 
	 */

	GMRFLib_graph_tp *g;
	void **args;
	double kappa;

	args = (void **) arg;				       /* the required variables come */
	g = (GMRFLib_graph_tp *) args[0];		       /* through the arg-pointer */
	kappa = *((double *) args[1]);

	return kappa * (node == nnode ? g->nnbs[node] : -1.0);
}
int main(int argc, char **argv)
{
	int seed, i, counter = 0;
	double fac, *x, *xnew, *d, kappa, kappa_new, lacc, eprob = 0.0, timeref;
	void *arguments_new[2], *arguments_old[2], *arguments_loglik[2];
	GMRFLib_constr_tp *constr = NULL;
	GMRFLib_blockupdate_param_tp *blockupdate_param;
	GMRFLib_optimize_param_tp *optimize_param;

	printf("Usage: %s SEED\n", *argv);
	assert(argc > 1);

	/*
	 * get the seed and init the random generator in GMRFLib 
	 */
	seed = atoi(argv[1]);
	printf("seed %d\n", seed);
	GMRFLib_uniform_init((unsigned long int) seed);

	/*
	 * get graph and data 
	 */
	read_graph_and_data();

	/*
	 * setup space 
	 */
	x = Calloc(G.graph->n, double);
	xnew = Calloc(G.graph->n, double);
	d = Calloc(G.graph->n, double);

	/*
	 * inital values 
	 */
	kappa = 10.;
	for (i = 0; i < G.graph->n; i++) {
		x[i] = DMIN(2., G.y[i] / G.E[i]);
		d[i] = 1.0;
	}

	/*
	 * the scale-factor for proposing a new value for kappa 
	 */
	fac = 1.5;

	/*
	 * setup arguments to the log-likelihood function, which needs data (y and E), the graph, and kappa. 
	 */
	arguments_loglik[0] = (void *) G.y;
	arguments_loglik[1] = (void *) G.E;
	arguments_new[0] = (void *) G.graph;
	arguments_new[1] = (void *) &kappa_new;
	arguments_old[0] = (void *) G.graph;
	arguments_old[1] = (void *) &kappa;

	/*
	 * choose whether or not the constraint \sum x_i = 0 is required 
	 */
	if (1) {
		GMRFLib_make_empty_constr(&constr);
		constr->nc = 1;
		constr->a_matrix = Calloc(G.graph->n, double);

		for (i = 0; i < G.graph->n; i++) {
			constr->a_matrix[i] = 1.0;
		}
		constr->e_vector = Calloc(1, double);

		constr->e_vector[0] = 0.0;
		GMRFLib_prepare_constr(constr, G.graph, 0);
	} else {
		constr = NULL;
	}

	/*
	 * get default options for the blockupdate-routine and the optimizing-routine which is (optionally) used in the blockupdate-routine.
	 * now it's easy to change options if required. 
	 */
	GMRFLib_default_blockupdate_param(&blockupdate_param);
	GMRFLib_default_optimize_param(&optimize_param);

	/*
	 * find the optimal reordering scheme. this pays off since we are factorising the same matrix many times. 
	 */
	GMRFLib_optimize_reorder(G.graph, NULL);

	/*
	 * start sampling, just go on forever 
	 */
	timeref = GMRFLib_timer();
	while (++counter) {
		/*
		 * propose a new value for kappa 
		 */
		kappa_new = kappa * GMRFLib_scale_proposal(fac);

		/*
		 * conditional on kappa_new, propose a new value for the GMRF x 
		 */

		GMRFLib_error_handler_tp *error_handler = GMRFLib_set_error_handler_off();

		loglikelihood_sign = -1.0;		       /* Force the routine to fail */
		int retval = GMRFLib_blockupdate(&lacc, xnew, x, NULL, NULL,
						 NULL, NULL,
						 NULL, NULL,
						 d, d, loglik, (void *) arguments_loglik, loglik, (void *) arguments_loglik, NULL,
						 G.graph, Qfunc, (void *) arguments_new, Qfunc, (void *) arguments_old, NULL, NULL,
						 NULL, NULL,
						 constr, constr, optimize_param, blockupdate_param);

		if (retval != GMRFLib_SUCCESS) {
			fprintf(stderr, "GMRFLib_blockupdate failed with error_code %d\n", retval);
			fprintf(stderr, "Switch the sign and try again....\n");

			loglikelihood_sign = 1.0;	       /* Make it correct */
			int retval = GMRFLib_blockupdate(&lacc, xnew, x, NULL, NULL,
							 NULL, NULL,
							 NULL, NULL,
							 d, d, loglik, (void *) arguments_loglik, loglik, (void *) arguments_loglik,
							 NULL,
							 G.graph, Qfunc, (void *) arguments_new, Qfunc, (void *) arguments_old,
							 NULL, NULL,
							 NULL, NULL,
							 constr, constr, optimize_param, blockupdate_param);

			if (retval != GMRFLib_SUCCESS) {
				fprintf(stderr, "GMRFLib_blockupdate failed again %d, I give up. \n", retval);
				exit(1);
			}
		}
		/*
		 * reset the error_handler back to the default one. 
		 */
		GMRFLib_set_error_handler(error_handler);

		/*
		 * the blockupdate-routine does not compute the normalization constant terms for the GMRF, nor the prior for kappa 
		 */
		lacc += (G.graph->n - 1.0) / 2.0 * log(kappa_new / kappa)
		    + eu_gamma(kappa_new, PRIOR_A, PRIOR_B)
		    - eu_gamma(kappa, PRIOR_A, PRIOR_B);

		/*
		 * standard accept/reject 
		 */
		if (GMRFLib_uniform() < exp(DMIN(0.0, lacc))) {
			kappa = kappa_new;
			memcpy(x, xnew, G.graph->n * sizeof(double));
		}

		/*
		 * that's it ! 
		 */
		eprob += exp(DMIN(0.0, lacc));
		printf("lacc= %.4f  E(accept_prob)= %.4f iter/src= %.4f  kappa= %.4f\n", lacc, eprob / counter,
		       1. / ((GMRFLib_timer() - timeref) / counter), kappa);
	}
	return 0;
}
