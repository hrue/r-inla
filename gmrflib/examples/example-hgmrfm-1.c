
/* example-hgmrf-1.c
 * 
 * Copyright (C) 2007-08 Havard Rue
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
   In this example, we will demonstrate how GMRFLib_init_hgmrfm() can be used to model a Log-Gaussian Cox-process, using
   the model
   
   \eta | ... ~ N( \beta_0 + \beta_1 * z_1 + \beta_2 * z_2 + spatial-smooth-term, 1/ prec_unstruct)
   
   This an example similar to the one in Section~5 in the INLA report by Rue, Martino and Chopin, from 2007.
*/

#define DATA_FILE        "data_small.dat"
#define NROW             20
#define NCOL             10
#define AREA             18.29

#define LEN_DATA         (NROW*NCOL)

#define PREC_PRIOR_A 1.0
#define PREC_PRIOR_B 0.01

static const char RCSId[] = "$Id: example-hgmrfm-1.c,v 1.17 2010/03/12 12:24:15 hrue Exp $";

typedef struct {
	int ndata;					       /* sizeof graph & data */
	double *data;
	double *z1;
	double *z2;
} Data_tp;

static Data_tp D;

double log_gamma(double x, double a, double b)
{
	return ((a - 1.0) * log(x) - x * b);
}
int read_data(void)
{
	int i;
	FILE *fp = fopen(DATA_FILE, "r");

	assert(fp);
	/*
	 * format: count z1 z2, where z1 and z2 are covariates 
	 */

	D.data = Calloc(LEN_DATA, double);
	D.z1 = Calloc(LEN_DATA, double);
	D.z2 = Calloc(LEN_DATA, double);

	for (i = 0; i < LEN_DATA; i++) {
		fscanf(fp, "%lf %lf %lf\n", &D.data[i], &D.z1[i], &D.z2[i]);
	}
	fclose(fp);

	return 0;
}
int loglikelihood(double *logll, double *x_i, int m, int idx, double *x_vec, void *arg)
{
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_DERIVATIES;
	}
	logll[0] = D.data[idx] * x_i[0] - AREA * exp(x_i[0]);
	assert(idx >= 0 && idx < LEN_DATA);
	if (m > 1) {
		logll[1] = D.data[idx] - AREA * exp(x_i[1]);
	}
	if (m > 2) {
		logll[2] = -AREA * exp(x_i[2]);
	}
	return GMRFLib_SUCCESS;
}
double extra(double *theta, int ntheta, void *args)
{
	double prec_spatial, prec_unstruct;

	prec_unstruct = exp(theta[0]);
	prec_spatial = exp(theta[1]);

	return 0.5 * (LEN_DATA - 3.0) * log(prec_spatial) + 0.5 * LEN_DATA * log(prec_unstruct) +
	    log_gamma(prec_spatial, PREC_PRIOR_A, PREC_PRIOR_B) + log_gamma(prec_unstruct, PREC_PRIOR_A, PREC_PRIOR_B)
	    + log(prec_spatial) + log(prec_unstruct);
}
int main(int argc, char **argv, char **env)
{
	double **log_prec_spatial, **log_prec_unstruct, **covariates, prec[3] = { 0.01, 0.01, 0.01 }, **hyperparam[2], *d;
	int i, **c;

	FILE *fp, *sfp;
	GMRFLib_Qfunc_tp **f_Qfunc;
	GMRFLib_ai_param_tp *ai_par;
	GMRFLib_ai_store_tp *ai_store;
	GMRFLib_density_tp **density;
	GMRFLib_graph_tp *rw2d_graph, **f_graph;
	GMRFLib_hgmrfm_tp *hgmrfm;
	char *f_sumzero, *compute;
	void **f_Qfunc_arg;
	int tmax = omp_get_max_threads();

	read_data();					       /* get the data and the covariates */

	/*
	 * define the RW2D model 
	 */
	GMRFLib_rw2ddef_tp *rw2d = Calloc(1, GMRFLib_rw2ddef_tp);

	rw2d->nrow = NROW;
	rw2d->ncol = NCOL;
	rw2d->cyclic = 0;
	rw2d->order = 2;
	log_prec_spatial = Calloc(tmax, double *);
	for (i = 0; i < tmax; i++) {
		log_prec_spatial[i] = Calloc(1, double);
		log_prec_spatial[i][0] = 1.0;
	}
	rw2d->log_prec_omp = log_prec_spatial;

	f_Qfunc_arg = Calloc(1, void *);

	f_Qfunc_arg[0] = (void *) rw2d;
	f_Qfunc = Calloc(1, GMRFLib_Qfunc_tp *);
	f_Qfunc[0] = GMRFLib_rw2d;

	GMRFLib_make_rw2d_graph(&rw2d_graph, rw2d);
	f_graph = Calloc(1, GMRFLib_graph_tp *);
	f_graph[0] = rw2d_graph;

	f_sumzero = Calloc(1, char);

	f_sumzero[0] = 1;				       /* yes please */

	c = Calloc(1, int *);
	c[0] = Calloc(LEN_DATA, int);

	for (i = 0; i < LEN_DATA; i++) {
		c[0][i] = i;
	}

	/*
	 * define the covariates 
	 */
	covariates = Calloc(3, double *);
	covariates[0] = Calloc(LEN_DATA, double);

	for (i = 0; i < LEN_DATA; i++) {
		covariates[0][i] = 1.0;			       /* constant term */
	}
	covariates[1] = D.z1;
	covariates[2] = D.z2;

	log_prec_unstruct = Calloc(tmax, double *);
	for (i = 0; i < tmax; i++) {
		log_prec_unstruct[i] = Calloc(1, double);
		log_prec_unstruct[i][0] = 1.0;
	}

	/*
	 * compute the hierarchical GMRF model 
	 */
	GMRFLib_init_hgmrfm(&hgmrfm, LEN_DATA, NULL, log_prec_unstruct,
			    1, c, NULL, f_graph, f_Qfunc, f_Qfunc_arg, f_sumzero, NULL, NULL, NULL, 3, covariates, prec, 0, NULL, NULL);

	/*
	 * ...and now we feed this into the _INLA() routine! 
	 */
	compute = Calloc(hgmrfm->graph->n, char);

	for (i = 0; i < hgmrfm->graph->n; i++) {
		compute[i] = 1;
	}

	ai_store = Calloc(1, GMRFLib_ai_store_tp);

	GMRFLib_default_ai_param(&ai_par);		       /* use default options */
	ai_par->strategy = GMRFLib_AI_STRATEGY_MEANSKEWCORRECTED_GAUSSIAN;
	ai_par->hessian_finite_difference_step_len = ai_par->gradient_finite_difference_step_len = 0.005;

	hyperparam[0] = log_prec_unstruct;
	hyperparam[1] = log_prec_spatial;
	d = Calloc(hgmrfm->graph->n, double);

	for (i = 0; i < LEN_DATA; i++) {
		d[i] = 1.0;				       /* the first LEN_DATA are \eta */
	}

	GMRFLib_ai_INLA(&density, NULL, NULL, NULL, NULL, NULL, NULL, compute, hyperparam, 2, extra, NULL, NULL, NULL, NULL, NULL,
			d, loglikelihood, NULL, NULL, hgmrfm->graph, hgmrfm->Qfunc, hgmrfm->Qfunc_arg, hgmrfm->constr, ai_par, ai_store,
			NULL, NULL, 0, NULL, NULL);

	/*
	 * write out the results 
	 */
	fp = fopen("lgcp-marginal-densities.dat", "w");
	sfp = fopen("lgcp-summary.dat", "w");

	for (i = 0; i < hgmrfm->graph->n; i++) {
		if (density[i]) {
			fprintf(fp, "%d ", i);
			double xx, x_real, f_real, f;

			for (xx = -4.0; xx < 4.0; xx += 0.1) {
				GMRFLib_evaluate_density(&f, xx, density[i]);

				x_real = xx * density[i]->std_stdev + density[i]->std_mean;
				f_real = f / density[i]->std_stdev;

				fprintf(fp, " %.10f %.10f ", x_real, f_real);
			}
			fprintf(fp, "\n");
			fprintf(sfp, "%d %f %f\n", i, density[i]->user_mean, density[i]->user_stdev);
		}
	}
	fclose(fp);
	fclose(sfp);

	/*
	 * cleanup 
	 */
	if (density) {
		for (i = 0; i < hgmrfm->graph->n; i++) {
			GMRFLib_free_density(density[i]);
		}
		Free(density);
	}
	Free(d);
	Free(ai_par);
	Free(compute);
	Free(c[0]);
	Free(c);
	GMRFLib_free_graph(rw2d_graph);
	Free(rw2d);
	Free(D.data);
	Free(D.z1);
	Free(D.z2);
	Free(f_Qfunc_arg);
	Free(f_Qfunc);
	Free(f_graph);
	GMRFLib_free_hgmrfm(hgmrfm);

	GMRFLib_timer_full_report(NULL);

	return 0;
}
