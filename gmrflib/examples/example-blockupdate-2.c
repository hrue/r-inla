
/* example-blockupdate-2.c
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
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"				       /* this is needed */
#include "GMRFLib/GMRFLibP.h"				       /* optional. obtain access some useful macros... */

/* 
   In this example, the model is

   y_ij  ~ Poisson( E_ij * exp( \eta_ij ))

   where

   \eta_ij  = \mu + g0(i) + g1(j) + \epsilon_ij

   We approach this problem using the routines in hgmrfm.c to build the joint graph and Qfunction for

   (\eta, g0, g1, \mu)

   and then apply the blockupdate routine.

   use the R-code in make-data-example-blockupdate-2.r to generate example-data.
*/

static const char RCSId[] = "$Id: example-blockupdate-2.c,v 1.8 2010/03/12 12:24:08 hrue Exp $";

typedef struct {
	double *y;
	double *E;
	int *g0idx;
	int *g1idx;
} Data_tp;

static Data_tp D;

#define NDATA 100
#define NG0 10
#define NG1 10
#define PREC_PRIOR_A 1.0
#define PREC_PRIOR_B 0.01

double log_gamma(double x, double a, double b)
{
	return ((a - 1.0) * log(x) - (x) * b);
}
int read_data(void)
{
	int i;
	FILE *fp = fopen("example-blockupdate-2.dat", "r");
	assert(fp);

	/*
	 * format: E y fidx gidx
	 */
	D.y = Calloc(NDATA, double);
	D.E = Calloc(NDATA, double);
	D.g0idx = Calloc(NDATA, int);
	D.g1idx = Calloc(NDATA, int);

	for (i = 0; i < NDATA; i++) {
		fscanf(fp, "%lf %lf %d %d\n", &D.E[i], &D.y[i], &D.g0idx[i], &D.g1idx[i]);
		// printf("%lf %lf %d %d\n", D.E[i], D.y[i], D.g0idx[i], D.g1idx[i]);
	}
	fclose(fp);

	return 0;
}
int loglikelihood(double *logll, double *x_i, int m, int idx, double *x_vec, void *arg)
{
	if (m == 0) {
		/*
		 * tell that we provide derivaties 
		 */
		return GMRFLib_LOGL_COMPUTE_DERIVATIES;
	}
	assert(idx >= 0 && idx < NDATA);		       /* just a safety check */

	logll[0] = D.y[idx] * x_i[0] - D.E[idx] * exp(x_i[0]);
	if (m > 1) {
		logll[1] = D.y[idx] - D.E[idx] * exp(x_i[1]);
	}
	if (m > 2) {
		logll[2] = -D.E[idx] * exp(x_i[2]);
	}
	return GMRFLib_SUCCESS;
}
int main(int argc, char **argv)
{
	double log_prec_unstruct = -1.0, log_prec_unstruct_new, prec_g0 = 10.0, prec_g0_new, prec_g1 = 10.0, prec_g1_new;

	GMRFLib_rwdef_tp *rw_g0, *rw_g0_new, *rw_g1, *rw_g1_new;
	GMRFLib_graph_tp *graph_g0, *graph_g1;
	int i, N, iter = 0, order = 2;

	double *x, *x_new, *eta, *g0, *g1, *mu;
	GMRFLib_hgmrfm_tp *hgmrfm, *hgmrfm_new;

	char *f_constr;
	int **c;
	GMRFLib_graph_tp **f_graph;
	GMRFLib_Qfunc_tp **f_Qfunc;
	void **f_Qfunc_arg, **f_Qfunc_arg_new;

	double *prec, **covariates;

	read_data();

	rw_g0 = Calloc(1, GMRFLib_rwdef_tp);
	rw_g0_new = Calloc(1, GMRFLib_rwdef_tp);
	rw_g1 = Calloc(1, GMRFLib_rwdef_tp);
	rw_g1_new = Calloc(1, GMRFLib_rwdef_tp);

	rw_g0->n = NG0;
	rw_g0->order = order;
	rw_g0->cyclic = GMRFLib_FALSE;
	rw_g0->prec = &prec_g0;

	rw_g0_new = rw_g0;				       /* identical, but ... */
	rw_g0_new->prec = &prec_g0_new;

	rw_g1->n = NG1;
	rw_g1->order = order;
	rw_g1->cyclic = GMRFLib_FALSE;
	rw_g1->prec = &prec_g1;

	rw_g1_new = rw_g1;
	rw_g1_new->prec = &prec_g1_new;

	GMRFLib_make_rw_graph(&graph_g0, rw_g0);
	GMRFLib_make_rw_graph(&graph_g1, rw_g1);

	/*
	 * build the model 
	 */
	f_Qfunc = Calloc(2, GMRFLib_Qfunc_tp *);
	f_Qfunc[0] = GMRFLib_rw;
	f_Qfunc[1] = GMRFLib_rw;

	f_Qfunc_arg = Calloc(2, void *);
	f_Qfunc_arg[0] = (void *) rw_g0;
	f_Qfunc_arg[1] = (void *) rw_g1;

	f_Qfunc_arg_new = Calloc(2, void *);
	f_Qfunc_arg_new[0] = (void *) rw_g0_new;
	f_Qfunc_arg_new[1] = (void *) rw_g1_new;

	f_graph = Calloc(2, GMRFLib_graph_tp *);
	f_graph[0] = graph_g0;
	f_graph[1] = graph_g1;

	f_constr = Calloc(2, char);
	f_constr[0] = 1;
	f_constr[1] = 1;

	c = Calloc(2, int *);
	c[0] = D.g0idx;
	c[1] = D.g1idx;

	covariates = Calloc(1, double *);
	covariates[0] = Calloc(NDATA, double);

	for (i = 0; i < NDATA; i++) {
		covariates[0][i] = 1.0;
	}
	prec = Calloc(1, double);
	prec[0] = 0.1;

	/*
	 * compute the hierarchical GMRF model 
	 */
	GMRFLib_init_hgmrfm(&hgmrfm, NDATA, &log_prec_unstruct, NULL, 2, c, NULL, f_graph, f_Qfunc, f_Qfunc_arg, f_constr, NULL, NULL, NULL, 1,
			    covariates, prec, 0, NULL, NULL);
	GMRFLib_init_hgmrfm(&hgmrfm_new, NDATA, &log_prec_unstruct_new, NULL, 2, c, NULL, f_graph, f_Qfunc, f_Qfunc_arg_new, f_constr, NULL, NULL,
			    NULL, 1, covariates, prec, 0, NULL, NULL);

	N = hgmrfm->graph->n;
	assert(N == NDATA + NG0 + NG1 + 1);
	x = Calloc(N, double);
	x_new = Calloc(N, double);

	/*
	 * just set the pointers so we know what is what 
	 */
	eta = x;
	g0 = x + NDATA;
	g1 = x + NDATA + NG0;
	mu = x + NDATA + NG0 + NG1;

	/*
	 * singular model, we need proper `priors', add to the diagonal for g0 and g1. 
	 */
	double *diag = Calloc(N, double);
	for (i = NDATA; i < NDATA + NG0 + NG1; i++) {
		diag[i] = 1.0e-09;
	}

	double *d = Calloc(N, double);
	for (i = 0; i < NDATA; i++) {
		d[i] = 1.0;				       /* rest is zero */
	}

	/*
	 * do some more effort to find a good reordering scheme, since the packing so to dense...
	 */
	GMRFLib_optimize_reorder(hgmrfm->graph, NULL);
	printf("Found good reordering scheme[%s]\n", GMRFLib_reorder_name(GMRFLib_reorder));

	double f = 1.1;					       /* step-size in MCMC */
	GMRFLib_store_tp *store = NULL;			       /* use store */
	store = Calloc(1, GMRFLib_store_tp);
	store->store_problems = GMRFLib_TRUE;

	FILE *fp_prec = fopen("trace-prec.dat", "w");
	FILE *fp_g0 = fopen("trace-g0.dat", "w");
	FILE *fp_g1 = fopen("trace-g1.dat", "w");
	FILE *fp_mu = fopen("trace-mu.dat", "w");

	assert(fp_prec);
	assert(fp_g0);
	assert(fp_g1);
	assert(fp_mu);

	double pmean = 0.0;				       /* average accept */

	while (1) {

		double lacc, prec_unstruct, prec_unstruct_new, p;

		/*
		 * this is bit special, as hgmrf_init require the log of the unstructed precision 
		 */
		prec_unstruct = exp(log_prec_unstruct);
		prec_unstruct_new = prec_unstruct * GMRFLib_scale_proposal(f);
		log_prec_unstruct_new = log(prec_unstruct_new);

		/*
		 * these are standard 
		 */
		prec_g0_new = prec_g0 * GMRFLib_scale_proposal(f);
		prec_g1_new = prec_g1 * GMRFLib_scale_proposal(f);

		GMRFLib_blockupdate_store(&lacc, x_new, x, NULL, NULL, diag, diag, NULL, NULL, d, d,
					  loglikelihood, NULL, loglikelihood, NULL, NULL,
					  hgmrfm->graph, hgmrfm_new->Qfunc, hgmrfm_new->Qfunc_arg, hgmrfm->Qfunc, hgmrfm->Qfunc_arg,
					  NULL, NULL, NULL, NULL, hgmrfm->constr, hgmrfm->constr, NULL, NULL, store);

		lacc += (NDATA / 2.0 * log_prec_unstruct_new + log_gamma(prec_unstruct_new, PREC_PRIOR_A, PREC_PRIOR_B)
			 + (NG0 - order) / 2.0 * log(prec_g0_new) + log_gamma(prec_g0_new, PREC_PRIOR_A, PREC_PRIOR_B)
			 + (NG1 - order) / 2.0 * log(prec_g1_new) + log_gamma(prec_g1_new, PREC_PRIOR_A, PREC_PRIOR_B))
		    - (NDATA / 2.0 * log_prec_unstruct + log_gamma(prec_unstruct, PREC_PRIOR_A, PREC_PRIOR_B)
		       + (NG0 - order) / 2.0 * log(prec_g0) + log_gamma(prec_g0, PREC_PRIOR_A, PREC_PRIOR_B)
		       + (NG1 - order) / 2.0 * log(prec_g1) + log_gamma(prec_g1, PREC_PRIOR_A, PREC_PRIOR_B));

		p = exp(DMIN(0.0, lacc));

		pmean = (iter * pmean + p) / (iter + 1.0);

		if (GMRFLib_uniform() < p || iter < 3) {
			memcpy(x, x_new, N * sizeof(double));
			prec_g0 = prec_g0_new;
			prec_g1 = prec_g1_new;
			log_prec_unstruct = log_prec_unstruct_new;
			store->decision = GMRFLib_STORE_ACCEPT;
		} else {
			store->decision = GMRFLib_STORE_REJECT;
		}

		if (!(iter % 10)) {
			// printf("iteration %d prob %f lacc %f\n", iter, p, lacc);
			printf("iteration %d pmean %g\n", iter, pmean);
			fprintf(fp_prec, "%g %g %g\n", exp(log_prec_unstruct), prec_g0, prec_g1);
			for (i = 0; i < NG0; i++) {
				fprintf(fp_g0, "%g ", g0[i]);
			}
			fprintf(fp_g0, "\n");

			for (i = 0; i < NG1; i++) {
				fprintf(fp_g1, "%g ", g1[i]);
			}
			fprintf(fp_g1, "\n");
			fprintf(fp_mu, "%g\n", mu[0]);

			/*
			 * this slow things down, but ensur 'clean' output-files. 
			 */
			fflush(fp_prec);
			fflush(fp_g0);
			fflush(fp_g1);
			fflush(fp_mu);
		}
		iter++;
	}
}
