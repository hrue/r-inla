
/* fgn.c
 * 
 * Copyright (C) 2016-17 Havard Rue
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
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * The author's contact information:
 *
 *        Haavard Rue
 *        CEMSE Division
 *        King Abdullah University of Science and Technology
 *        Thuwal 23955-6900, Saudi Arabia
 *        Email: haavard.rue@kaust.edu.sa
 *        Office: +966 (0)12 808 0640
 *
 */
#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = HGVERSION;

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

#include "fgn.h"
#include "interpol.h"

int inla_make_fgn_graph(GMRFLib_graph_tp ** graph, inla_fgn_arg_tp * def)
{
	int i, j;
	GMRFLib_graph_tp *g_ar1 = NULL, *g_I = 0;
	GMRFLib_ged_tp *ged = NULL;

	GMRFLib_make_linear_graph(&g_ar1, def->n, 1, 0);
	GMRFLib_make_linear_graph(&g_I, def->n, 0, 0);
	GMRFLib_ged_init(&ged, NULL);

	GMRFLib_ged_insert_graph2(ged, g_I, 0, 0);
	for (i = 0; i < def->k; i++) {
		GMRFLib_ged_insert_graph2(ged, g_I, 0, (1 + i) * def->n);
	}
	for (i = 0; i < def->k; i++) {
		GMRFLib_ged_insert_graph2(ged, g_ar1, (i + 1) * def->n, (i + 1) * def->n);
		for (j = 1; j < def->k - i; j++) {
			GMRFLib_ged_insert_graph2(ged, g_I, (1 + i) * def->n, (1 + i + j) * def->n);
		}
	}
	assert(GMRFLib_ged_max_node(ged) == def->N - 1);
	GMRFLib_ged_build(graph, ged);

	GMRFLib_ged_free(ged);
	GMRFLib_free_graph(g_ar1);
	GMRFLib_free_graph(g_I);

	return (GMRFLib_SUCCESS);
}

int inla_make_fgn2_graph(GMRFLib_graph_tp ** graph, inla_fgn2_arg_tp * def)
{
	int i;
	GMRFLib_graph_tp *g_ar1 = NULL;
	GMRFLib_ged_tp *ged = NULL;

	GMRFLib_make_linear_graph(&g_ar1, def->n, 1, 0);
	GMRFLib_ged_init(&ged, NULL);

	for (i = 0; i < def->k; i++) {
		GMRFLib_ged_insert_graph2(ged, g_ar1, i * def->n, i * def->n);
	}
	for (i = 0; i < def->k - 1; i++) {
		GMRFLib_ged_insert_graph2(ged, g_ar1, i * def->n, (i + 1) * def->n);
	}
	assert(GMRFLib_ged_max_node(ged) == def->N - 1);
	GMRFLib_ged_build(graph, ged);
	GMRFLib_free_graph(g_ar1);

	return (GMRFLib_SUCCESS);
}



double Qfunc_fgn(int i, int j, void *arg)
{
	// the model (z,x1,x2,x3,...), where z = 1/\sqrt{prec} * \sum_i \sqrt{w_i} x_i + tiny.noise,
	// where each x is standard AR1

	int debug = 0;
	static double **phi_cache = NULL, **w_cache = NULL, *H_intern_cache = NULL;

	if (!arg) {
		assert(i < 0 && j < 0);			       /* safety check */
		if (phi_cache == NULL) {
#pragma omp critical 
			{
				if (phi_cache == NULL) {
					phi_cache = Calloc(ISQR(GMRFLib_MAX_THREADS), double *);
					w_cache = Calloc(ISQR(GMRFLib_MAX_THREADS), double *);
					H_intern_cache = Calloc(ISQR(GMRFLib_MAX_THREADS), double);
					
					for (int j = 0; j < ISQR(GMRFLib_MAX_THREADS); j++) {
						phi_cache[j] = Calloc(2 * FGN_KMAX - 1, double);
						w_cache[j] = Calloc(2 * FGN_KMAX - 1, double);
					}
					if (debug) {
						printf("Qfunc_fgn: initialize cache\n");
					}
				}
			}
		}
		return NAN;				       /* so it will break if used wrong */
	}

	inla_fgn_arg_tp *a = (inla_fgn_arg_tp *) arg;
	double H_intern, prec, val = 0.0, *phi, *w, kappa;
	int id = omp_get_thread_num() * GMRFLib_MAX_THREADS + GMRFLib_thread_id;

	phi = phi_cache[id];
	w = w_cache[id];

	H_intern = a->H_intern[GMRFLib_thread_id][0];
	prec = map_precision(a->log_prec[GMRFLib_thread_id][0], MAP_FORWARD, NULL);
	kappa = a->prec_eps * prec;

	if (!ISEQUAL(H_intern, H_intern_cache[id])) {
		if (debug) {
			printf("Qfunc_fgn: update cache H_intern[%1d]= %f\n", id, H_intern);
		}
		inla_fgn_get(phi, w, H_intern, a->k);
		H_intern_cache[id] = H_intern;
		if (debug) {
			for (int k = 0; k < a->k; k++)
				printf("\tphi[%1d]= %f   w[%1d]= %f\n", k, phi[k], k, w[k]);
		}
	}

	div_t ii, jj;
	ii = div(IMIN(i, j), a->n);
	jj = div(IMAX(i, j), a->n);

	// the x^i's are the scaled AR1's, and the FGN is then just the sum of the
	// components.
	if (ii.quot == 0) {
		val = kappa * (jj.quot == 0 ? 1.0 : -1.0);
	} else {
		if (ii.quot == jj.quot) {
			// this is the AR1
			double prec_cond = 1.0 / (1.0 - SQR(phi[ii.quot - 1L]));
			double scale = prec / w[ii.quot - 1L];
			if (ii.rem != jj.rem) {
				// off-diagonal
				val = -scale * prec_cond * phi[ii.quot - 1L];
			} else {
				// diagonal
				val = scale * prec_cond * ((ii.rem == 0 || ii.rem == a->n - 1L) ? 1.0 : (1.0 + SQR(phi[ii.quot - 1L])));
				val += kappa;
			}
		} else {
			val = kappa;
		}
	}

	return val;
}

double inla_fgn2_helper(int i, int j, int n, double phi)
{
	double phi2 = SQR(phi);
	return ((i != j ? -phi : ((i == 0 || i == n - 1) ? 1.0 : (1.0 + phi2))) / (1.0 - phi2));

	if (0) {
		double val;

		if (i != j) {
			val = -phi;
		} else {
			if (i == 0 || i == n - 1) {
				val = 1.0;
			} else {
				val = 1.0 + phi2;
			}
		}
		return val / (1.0 - phi2);		       /* correction to get unit variance */
	}
}

double Qfunc_fgn2(int i, int j, void *arg)
{
	// the x^i's are the scaled AR1's, and FGN is the cummulative sum of the components.

	int debug = 0;
	static double **phi_cache = NULL, **w_cache = NULL, *H_intern_cache = NULL;

	if (!arg) {
		assert(phi_cache == NULL);		       /* do not initialize twice */
		phi_cache = Calloc(ISQR(GMRFLib_MAX_THREADS), double *);
		w_cache = Calloc(ISQR(GMRFLib_MAX_THREADS), double *);
		H_intern_cache = Calloc(ISQR(GMRFLib_MAX_THREADS), double);

		for (int j = 0; j < ISQR(GMRFLib_MAX_THREADS); j++) {
			phi_cache[j] = Calloc(2 * FGN_KMAX - 1, double);
			w_cache[j] = Calloc(2 * FGN_KMAX - 1, double);
		}
		if (debug) {
			printf("Qfunc_fgn2: initialize cache\n");
		}
		return NAN;
	}

	inla_fgn2_arg_tp *a = (inla_fgn2_arg_tp *) arg;
	double H_intern, prec, val = 0.0, *phi, *w;
	int id = omp_get_thread_num() * GMRFLib_MAX_THREADS + GMRFLib_thread_id;

	phi = phi_cache[id];
	w = w_cache[id];

	H_intern = a->H_intern[GMRFLib_thread_id][0];
	prec = map_precision(a->log_prec[GMRFLib_thread_id][0], MAP_FORWARD, NULL);

	if (!ISEQUAL(H_intern, H_intern_cache[id])) {
		if (debug) {
			printf("Qfunc_fgn2: update cache H_intern[%1d]= %f\n", id, H_intern);
		}
		inla_fgn_get(phi, w, H_intern, a->k);
		H_intern_cache[id] = H_intern;
		if (debug) {
			for (int k = 0; k < a->k; k++)
				printf("\tphi[%1d]= %f   w[%1d]= %f\n", k, phi[k], k, w[k]);
		}
	}

	div_t ii = div(IMIN(i, j), a->n);
	div_t jj = div(IMAX(i, j), a->n);
	int k, kk;
	double qval, qqval;

	// this is to get the slowest components last, otherwise use k = ii.quot and kk=k-1
	k = a->k - ii.quot - 1;
	kk = k + 1;

	qval = inla_fgn2_helper(ii.rem, jj.rem, a->n, phi[k]) / w[k];
	qqval = ((ii.quot > 0) ? inla_fgn2_helper(ii.rem, jj.rem, a->n, phi[kk]) / w[kk] : 0.0);
	val = prec * ((ii.quot == jj.quot) ? qval + qqval : -qval);

	return val;
}


int inla_fgn_get(double *phi, double *w, double H_intern, int k)
{
	// fill in the weights and the phis for a given H_intern
	if (k == 3) {
#include "fgn-tables-3.h"
#include "fgn-code.h"
	} else if (k == 4) {
#include "fgn-tables-4.h"
#include "fgn-code.h"
	} else {
		GMRFLib_ASSERT(k == 3 || k == 4, GMRFLib_EPARAMETER);
	}
	abort();
}

double priorfunc_fgn_priorH(double *H_intern, double *param)
{
	// return the log-prior for H_intern
	double lprior;
#include "fgn-prior-tables.h"
#pragma omp critical
	{
		static GMRFLib_spline_tp *dist_spline = NULL;
		if (!dist_spline) {
			dist_spline = inla_spline_create(H_int, Dist, sizeof(H_int) / sizeof(double));
		}

		double U_intern, lambda;
		U_intern = map_H(param[0], MAP_BACKWARD, NULL);
		lambda = -log(param[1]) / inla_spline_eval(U_intern, dist_spline);
		lprior = log(lambda) - lambda * inla_spline_eval(*H_intern, dist_spline) +
		    log(fabs(inla_spline_eval_deriv(*H_intern, dist_spline)));

		if (0) {
			P(*H_intern);
			P(lambda);
			P(inla_spline_eval(*H_intern, dist_spline));
			P(inla_spline_eval_deriv(*H_intern, dist_spline));
			P(lprior);
		}
	}
	return lprior;
}
