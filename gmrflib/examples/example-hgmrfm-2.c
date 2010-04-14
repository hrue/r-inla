
/* example-hgmrf-2.c
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

static const char RCSId[] = "$Id: example-hgmrfm-2.c,v 1.15 2010/03/12 12:24:12 hrue Exp $";

#include "GMRFLib/GMRFLib.h"

/* 
   This is a similar example to example-hgmrfm-1.c, but now we replace the effect of covariates z1 and z2 with smooth curve
   (RW2).
*/

#define DATA_FILE        "data_small.dat"
#define NROW             20
#define NCOL             10
#define AREA             18.29

#define LEN_DATA         (NROW*NCOL)

#define PREC_PRIOR_A 1.0
#define PREC_PRIOR_B 0.01

#define LEN_F 100


typedef struct {
	int ndata;					       /* sizeof graph & data */
	double *data;
	double *z1;
	double *z2;
} Data_tp;

static Data_tp D;

double log_gamma(double x, double a, double b)
{
	return ((a - 1.0) * log(x) - (x) * b);
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
	double prec_spatial, prec_unstruct, prec_z1, prec_z2;

	prec_unstruct = exp(theta[0]);
	prec_spatial = exp(theta[1]);
	prec_z1 = exp(theta[2]);
	prec_z2 = exp(theta[3]);

	return 0.5 * (LEN_DATA - 3.0) * log(prec_spatial)
	    + 0.5 * LEN_DATA * log(prec_unstruct) + +0.5 * (LEN_F - 2.0) * log(prec_z1)
	    + 0.5 * (LEN_F - 2.0) * log(prec_z2)
	    + log_gamma(prec_spatial, PREC_PRIOR_A, PREC_PRIOR_B)
	    + log_gamma(prec_unstruct, PREC_PRIOR_A, PREC_PRIOR_B)
	    + log_gamma(prec_z1, PREC_PRIOR_A, PREC_PRIOR_B)
	    + log_gamma(prec_z2, PREC_PRIOR_A, PREC_PRIOR_B)
	    + log(prec_spatial)
	    + log(prec_unstruct)
	    + log(prec_z1)
	    + log(prec_z2);
}
double GMRFLib_rw_extra(int i, int j, void *arg)
{
	/*
	 * we add a small term on the diagonal; otherwise the precision matrix is singular. 
	 */
	return GMRFLib_rw(i, j, arg) + (i == j ? 1.0e-6 : 0.0);
}
int main(int argc, char **argv, char **env)
{
	double **log_prec_spatial, **log_prec_unstruct, **log_prec_z1, **log_prec_z2, **covariates, prec[1] = { 0.01 }, **hyperparam[4], *d;
	int i;

	FILE *fp, *sfp;
	GMRFLib_ai_param_tp *ai_par;
	GMRFLib_ai_store_tp *ai_store;
	GMRFLib_density_tp **density;
	GMRFLib_hgmrfm_tp *hgmrfm;
	char *f_constr, *compute;

	int **c;
	GMRFLib_rwdef_tp *rw_z1, *rw_z2;
	GMRFLib_graph_tp *graph_z1, *graph_z2, **f_graph, *rw2d_graph;
	GMRFLib_Qfunc_tp **f_Qfunc;
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
		log_prec_spatial[i][0] = 2.0;
	}
	rw2d->log_prec_omp = log_prec_spatial;

	f_Qfunc_arg = Calloc(3, void *);

	f_Qfunc_arg[0] = (void *) rw2d;
	f_Qfunc = Calloc(3, GMRFLib_Qfunc_tp *);
	f_Qfunc[0] = GMRFLib_rw2d;

	GMRFLib_make_rw2d_graph(&rw2d_graph, rw2d);
	f_graph = Calloc(3, GMRFLib_graph_tp *);
	f_graph[0] = rw2d_graph;

	f_constr = Calloc(3, char);

	f_constr[0] = 1;				       /* yes please */

	c = Calloc(3, int *);
	c[0] = Calloc(LEN_DATA, int);

	for (i = 0; i < LEN_DATA; i++) {
		c[0][i] = i;
	}

	log_prec_z1 = Calloc(tmax, double *);
	for (i = 0; i < tmax; i++) {
		log_prec_z1[i] = Calloc(1, double);
		log_prec_z1[i][0] = 2.0;
	}
	log_prec_z2 = Calloc(tmax, double *);
	for (i = 0; i < tmax; i++) {
		log_prec_z2[i] = Calloc(1, double);
		log_prec_z2[i][0] = 2.0;
	}

	/*
	 * smooth functions to account for covariates 
	 */
	rw_z1 = Calloc(1, GMRFLib_rwdef_tp);
	rw_z1->n = LEN_F;
	rw_z1->order = 2;
	rw_z1->cyclic = 0;
	rw_z1->cyclic = 0;
	rw_z1->log_prec_omp = log_prec_z1;

	rw_z2 = Calloc(1, GMRFLib_rwdef_tp);
	rw_z2->n = LEN_F;
	rw_z2->order = 2;
	rw_z2->cyclic = 0;
	rw_z2->cyclic = 0;
	rw_z2->log_prec_omp = log_prec_z2;

	GMRFLib_make_rw_graph(&graph_z1, rw_z1);
	GMRFLib_make_rw_graph(&graph_z2, rw_z2);

	f_graph[1] = graph_z1;
	f_graph[2] = graph_z2;
	f_Qfunc[1] = GMRFLib_rw_extra;
	f_Qfunc[2] = GMRFLib_rw_extra;
	f_Qfunc_arg[1] = (void *) rw_z1;
	f_Qfunc_arg[2] = (void *) rw_z2;

	c[1] = Calloc(LEN_DATA, int);
	c[2] = Calloc(LEN_DATA, int);

	for (i = 0; i < LEN_DATA; i++) {
		c[1][i] = (int) (D.z1[i] * 100);	       /* z1 is between 0 and 1; use 100 different values */
		c[2][i] = (int) (D.z2[i] * 100);	       /* z2 is between 0 and 1; use 100 different values */
	}

	f_constr[1] = 1;				       /* yes please */
	f_constr[2] = 1;				       /* yes please */

	/*
	 * define the covariates 
	 */
	covariates = Calloc(1, double *);
	covariates[0] = Calloc(LEN_DATA, double);

	for (i = 0; i < LEN_DATA; i++) {
		covariates[0][i] = 1.0;			       /* constant term */
	}

	log_prec_unstruct = Calloc(tmax, double *);
	for (i = 0; i < tmax; i++) {
		log_prec_unstruct[i] = Calloc(1, double);
		log_prec_unstruct[i][0] = 2.0;
	}

	/*
	 * compute the hierarchical GMRF model 
	 */
	GMRFLib_init_hgmrfm(&hgmrfm, LEN_DATA, NULL, log_prec_unstruct,
			    3, c, NULL, f_graph, f_Qfunc, f_Qfunc_arg, f_constr, NULL, NULL, NULL, 1, covariates, prec, 0, NULL, NULL);

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
	ai_par->int_strategy = GMRFLib_AI_INT_STRATEGY_CCD;

	hyperparam[0] = log_prec_unstruct;
	hyperparam[1] = log_prec_spatial;
	hyperparam[2] = log_prec_z1;
	hyperparam[3] = log_prec_z2;
	d = Calloc(hgmrfm->graph->n, double);

	for (i = 0; i < LEN_DATA; i++) {
		d[i] = 1.0;				       /* the first LEN_DATA are \eta */
	}

	GMRFLib_ai_INLA(&density, NULL, NULL, NULL, NULL, NULL, NULL, compute, hyperparam, 4, extra, NULL, NULL, NULL, NULL, NULL,
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
	Free(covariates[0]);
	Free(covariates);
	GMRFLib_free_graph(rw2d_graph);
	Free(rw2d);
	Free(D.data);
	Free(D.z1);
	Free(D.z2);
	Free(f_Qfunc_arg);
	Free(f_Qfunc);
	Free(f_graph);
	GMRFLib_free_hgmrfm(hgmrfm);

	Free(rw_z1);
	Free(rw_z2);
	GMRFLib_free_graph(graph_z1);
	GMRFLib_free_graph(graph_z2);
	Free(c[0]);
	Free(c[1]);
	Free(c[2]);
	Free(c);
	Free(f_constr);

	GMRFLib_timer_full_report(NULL);

	return 0;
}
