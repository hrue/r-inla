
/* fgn.c
 * 
 * Copyright (C) 2016 Havard Rue
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
 *       H{\aa}vard Rue
 *       Department of Mathematical Sciences
 *       The Norwegian University of Science and Technology
 *       N-7491 Trondheim, Norway
 *       Voice: +47-7359-3533    URL  : http://www.math.ntnu.no/~hrue  
 *       Fax  : +47-7359-3524    Email: havard.rue@math.ntnu.no
 *
 */
#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = HGVERSION;

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

#include "fgn.h"

#define KMAX (6L)

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

double Qfunc_fgn(int i, int j, void *arg)
{
	// the model (z,x1,x2,x3,...), where z = 1/\sqrt{prec} * \sum_i \sqrt{w_i} x_i + tiny.noise,
	// where each x is standard AR1

	static double phi[2L * KMAX - 1L], w[2L * KMAX - 1L], H_ref = -1.0;
#pragma omp threadprivate(phi, w, H_ref)

	inla_fgn_arg_tp *a = (inla_fgn_arg_tp *) arg;
	double H = 0.7888, prec = 1.0, val = 0.0;
	int debug = 0;

	if (!ISEQUAL(H, H_ref)) {
#pragma omp critical
		{
			if (!ISEQUAL(H, H_ref)) {
				inla_fng_get(phi, w, H, a->k);
				H_ref = H;
			}
		}
	}

	if (debug) {
#pragma omp critical
		{
			for (int k = 0; k < a->k; k++)
				printf("\t phi[%1d] = %f   w[%1d] = %f\n", k, phi[k], k, w[k]);
		}
	}

	div_t ii, jj;
	ii = div(IMIN(i, j), a->n);
	jj = div(IMAX(i, j), a->n);

	if (ii.quot == 0) {
		val = (jj.quot == 0 ? a->prec_eps : -sqrt(w[jj.quot - 1L] / prec) * a->prec_eps);
	} else {
		if (ii.quot == jj.quot) {
			// this is the AR1
			double prec_cond = 1.0 / (1.0 - SQR(phi[ii.quot - 1L]));
			if (ii.rem != jj.rem) {
				// off-diagonal
				val = -prec_cond * phi[ii.quot - 1L];
			} else {
				// diagonal
				val = prec_cond * ((ii.rem == 0 || ii.rem == a->n - 1L) ? 1.0 : (1.0 + SQR(phi[ii.quot - 1L])));
				val += a->prec_eps * w[ii.quot - 1L] / prec;
			}
		} else {
			val = sqrt(w[ii.quot - 1L] * w[jj.quot - 1L]) * a->prec_eps / prec;
		}
	}

	return val;
}

int inla_fng_get(double *phi, double *w, double H, int k)
{
	// fill in the weights and the phis for a given H
#include "fgn-tables.h"

	int idx, i, len_par;
	double weight, tmp;

	assert(k == K);
	assert(H >= H_start && H <= H_end);
	idx = (int) floor((H - H_start) / H_by);	       /* idx is the block-index */
	weight = (H - (H_start + idx * H_by)) / H_by;
	len_par = 2 * K - 1;
	idx *= len_par;					       /* and now the index in the table */

	double *fit_par = Calloc(len_par, double);
	for (i = 0; i < len_par; i++) {
		fit_par[i] = (1.0 - weight) * param[idx + i] + weight * param[idx + len_par + i];
	}

	// the first K are phi
	for (i = 0, tmp = 0.0; i < K; i++) {
		tmp += exp(-fit_par[i]);
		phi[i] = 1.0 / (1.0 + tmp);
	}

	// the remaining K-1 are the weights
	double psum, *par = Calloc(len_par, double);
	par[0] = psum = 1;
	for (i = 1; i < K; i++) {
		par[i] = exp(fit_par[K + (i - 1)]);
		psum += par[i];
	}
	for (i = 0; i < K; i++) {
		w[i] = par[i] / psum;
	}

	Free(fit_par);
	Free(par);

	return GMRFLib_SUCCESS;
}
