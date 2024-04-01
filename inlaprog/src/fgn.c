
/* fgn.c
 * 
 * Copyright (C) 2016-2024 Havard Rue
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

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

#include "fgn.h"

int inla_make_fgn_graph(GMRFLib_graph_tp **graph, inla_fgn_arg_tp *def)
{
	int i, j;
	GMRFLib_graph_tp *g_ar1 = NULL, *g_I = 0;
	GMRFLib_ged_tp *ged = NULL;

	GMRFLib_graph_mk_linear(&g_ar1, def->n, 1, 0);
	GMRFLib_graph_mk_linear(&g_I, def->n, 0, 0);
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
	assert(ged->n == def->N);
	GMRFLib_ged_build(graph, ged);

	GMRFLib_ged_free(ged);
	GMRFLib_graph_free(g_ar1);
	GMRFLib_graph_free(g_I);

	return (GMRFLib_SUCCESS);
}

int inla_make_fgn2_graph(GMRFLib_graph_tp **graph, inla_fgn2_arg_tp *def)
{
	int i;
	GMRFLib_graph_tp *g_ar1 = NULL;
	GMRFLib_ged_tp *ged = NULL;

	GMRFLib_graph_mk_linear(&g_ar1, def->n, 1, 0);
	GMRFLib_ged_init(&ged, NULL);

	for (i = 0; i < def->k; i++) {
		GMRFLib_ged_insert_graph2(ged, g_ar1, i * def->n, i * def->n);
	}
	for (i = 0; i < def->k - 1; i++) {
		GMRFLib_ged_insert_graph2(ged, g_ar1, i * def->n, (i + 1) * def->n);
	}
	assert(ged->n == def->N);
	GMRFLib_ged_build(graph, ged);
	GMRFLib_graph_free(g_ar1);

	return (GMRFLib_SUCCESS);
}

double Qfunc_fgn(int thread_id, int i, int j, double *UNUSED(values), void *arg)
{
	if (j < 0) {
		return NAN;
	}
	// the model (z,x1,x2,x3,...), where z = 1/\sqrt{prec} * \sum_i \sqrt{w_i} x_i + tiny.noise,
	// where each x is standard AR1

	const int debug = 0;
	static double **phi_cache = NULL, **w_cache = NULL, *H_intern_cache = NULL;

	if (phi_cache == NULL) {
#pragma omp critical (Name_6cee800e55124771d0e7fd552ae7e48a27e4f94e)
		{
			if (phi_cache == NULL) {
				double **cache = Calloc(GMRFLib_CACHE_LEN(), double *);
				w_cache = Calloc(GMRFLib_CACHE_LEN(), double *);
				H_intern_cache = Calloc(GMRFLib_CACHE_LEN(), double);

				for (int jj = 0; jj < GMRFLib_CACHE_LEN(); jj++) {
					cache[jj] = Calloc(2 * FGN_KMAX - 1, double);
					w_cache[jj] = Calloc(2 * FGN_KMAX - 1, double);
				}
				phi_cache = cache;
			}
		}
	}

	inla_fgn_arg_tp *a = (inla_fgn_arg_tp *) arg;
	double H_intern, prec, val = 0.0, *phi = NULL, *w = NULL, kappa;
	int id = 0;

	GMRFLib_CACHE_SET_ID(id);
	phi = phi_cache[id];
	w = w_cache[id];

	H_intern = a->H_intern[thread_id][0];
	prec = map_precision(a->log_prec[thread_id][0], MAP_FORWARD, NULL);
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
	assert(0 == 1);
	return 0.0;
}

double Qfunc_fgn2(int thread_id, int i, int j, double *UNUSED(values), void *arg)
{
	if (j < 0) {
		return NAN;
	}
	// the x^i's are the scaled AR1's, and FGN is the cummulative sum of the components.

	const int debug = 0;
	static double **phi_cache = NULL, **w_cache = NULL, *H_intern_cache = NULL;

	if (!phi_cache) {
#pragma omp critical (Name_31036ca2cfd217477a399b276d2192bbc39a5fb7)
		if (!phi_cache) {
			double **cache = Calloc(GMRFLib_CACHE_LEN(), double *);
			w_cache = Calloc(GMRFLib_CACHE_LEN(), double *);
			H_intern_cache = Calloc(GMRFLib_CACHE_LEN(), double);

			for (int jj = 0; jj < GMRFLib_CACHE_LEN(); jj++) {
				cache[jj] = Calloc(2 * FGN_KMAX - 1, double);
				w_cache[jj] = Calloc(2 * FGN_KMAX - 1, double);
			}
			phi_cache = cache;
		}
	}

	inla_fgn2_arg_tp *a = (inla_fgn2_arg_tp *) arg;
	double H_intern, prec, val = 0.0, *phi = NULL, *w = NULL;
	int id = 0;

	GMRFLib_CACHE_SET_ID(id);
	phi = phi_cache[id];
	w = w_cache[id];

	H_intern = a->H_intern[thread_id][0];
	prec = map_precision(a->log_prec[thread_id][0], MAP_FORWARD, NULL);

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

	static GMRFLib_spline_tp *dist_spline = NULL;
#pragma omp critical (Name_f88269b9720b21345f72723d8de2fc329de96a39)
	{
		if (!dist_spline) {
			dist_spline = GMRFLib_spline_create(H_int, Dist, sizeof(H_int) / sizeof(double));
		}
	}

	double U_intern, lambda;
	U_intern = map_H(param[0], MAP_BACKWARD, NULL);
	lambda = -log(param[1]) / GMRFLib_spline_eval(U_intern, dist_spline);
	lprior = log(lambda) - lambda * GMRFLib_spline_eval(*H_intern, dist_spline) + log(fabs(GMRFLib_spline_eval_deriv(*H_intern, dist_spline)));

	return lprior;
}

void priorfunc_fgn_priorH_extract(void)
{
	// extract the prior so we can get it into R
#include "fgn-prior-tables.h"

	static GMRFLib_spline_tp *dist_spline = NULL;
#pragma omp critical (Name_1083cc49be9497f3e0b14820ba227f6584988f41)
	{
		if (!dist_spline) {
			dist_spline = GMRFLib_spline_create(H_int, Dist, sizeof(H_int) / sizeof(double));
		}
	}

	for (double H_intern = -10.0, dH = 0.05; H_intern <= 10.0 + dH / 2.0; H_intern += dH) {
		double dist = GMRFLib_spline_eval(H_intern, dist_spline);
		printf("%.12f %.12f\n", H_intern, dist);
	}

	double param[] = { 0.9, 0.1 };
	printf("\n\n## check:  param = 0.9 0.1\n");

	double theta;

	theta = -1.1;
	printf("## check:  log.prior(%.8f) =  %.8f\n", theta, priorfunc_fgn_priorH(&theta, param));
	theta = 0.2;
	printf("## check:  log.prior(%.8f) =  %.8f\n", theta, priorfunc_fgn_priorH(&theta, param));
	theta = 1.3;
	printf("## check:  log.prior(%.8f) =  %.8f\n", theta, priorfunc_fgn_priorH(&theta, param));
}
