
/* example-auxvar.c
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

static const char RCSId[] = "$Id: example-auxvar.c,v 1.5 2010/03/12 12:24:33 hrue Exp $";

#include "GMRFLib/GMRFLib.h"

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
double Qfunc(int node, int nnode, void *arg);

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
	char *fnm;
	GMRFLib_io_tp *io = NULL;

	fnm = Strdup(GRAPH);
	GMRFLib_read_graph(&G.graph, fnm);
	Free(fnm);

	/*
	 * use the io-interface to read 
	 */
	G.y = Calloc(G.graph->n, double);
	G.E = Calloc(G.graph->n, double);

	GMRFLib_io_open(&io, DATA, "r");
	for (i = 0; i < G.graph->n; i++) {
		GMRFLib_io_read_next(io, &G.y[i], "%lf");
		GMRFLib_io_read_next(io, &G.E[i], "%lf");
		GMRFLib_io_read_next(io, &dummy, "%lf");
	}
	GMRFLib_io_close(io);

	return 0;
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
	int i, counter = 0;
	double fac, *x, *xnew, *d, kappa, kappa_new, lacc, eprob = 0.0, timeref, *b, *c;
	void *arguments_new[2], *arguments_old[2];
	GMRFLib_constr_tp *constr = NULL;

	/*
	 * GMRFLib now sets the seed reading /dev/urandom 
	 */

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
	b = Calloc(G.graph->n, double);
	c = Calloc(G.graph->n, double);

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
	 * setup auxilliary variables. Note that pointers are stored, so a change in 'x' is noted by the _aux-routines.  The vector 'd' is the
	 * same 'd' as would have been used in the _blockupdate()-routine. 
	 */
	GMRFLib_auxs_tp *auxs = NULL;

	GMRFLib_aux_setup_poisson_all(&auxs, G.graph->n, d, G.y, G.E, x);
	GMRFLib_aux_init_all(auxs);

	/*
	 * store intermediate calculations 
	 */
	GMRFLib_store_tp *store = calloc(1, sizeof(GMRFLib_store_tp));

	/*
	 * start sampling, just go on forever 
	 */
	timeref = GMRFLib_timer();
	while (++counter) {

		/*
		 * get the equiv gauss data 
		 */
		GMRFLib_aux_gauss_approx_all(b, c, auxs);

		/*
		 * propose a new value for kappa 
		 */
		kappa_new = kappa * GMRFLib_scale_proposal(fac);

		/*
		 * conditional on kappa_new and the auxilliary variables, propose a new value for the GMRF x
		 *
		 * NOTE: the 'b' and 'c' are not necessarily the 'b' and 'c's in the _blockupdate() routine is the mean is
		 * non-zero. I'm sorry for this `design'-error. //HRue
		 */
		GMRFLib_blockupdate_store(&lacc, xnew, x,
					  b, b,
					  c, c,
					  NULL, NULL,
					  NULL, NULL,
					  NULL, NULL,
					  NULL, NULL,
					  NULL,
					  G.graph, Qfunc, (void *) arguments_new, Qfunc, (void *) arguments_old, NULL, NULL,
					  NULL, NULL, constr, constr, NULL, NULL, store);

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
		 * update the aux-variables 
		 */
		GMRFLib_aux_update_all(auxs);

		/*
		 * that's it ! 
		 */
		eprob += exp(DMIN(0.0, lacc));
		printf("lacc= %.4f  E(accept_prob)= %.4f iter/src= %.4f  kappa= %.4f\n", lacc, eprob / counter, 1. / ((GMRFLib_timer() - timeref) / counter), kappa);
	}
	return 0;
}
