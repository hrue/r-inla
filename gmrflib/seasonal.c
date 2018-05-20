
/* seasonal.c
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
static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

/* Pre-hg-Id: $Id: seasonal.c,v 1.13 2008/10/26 03:21:48 hrue Exp $ */

#include <time.h>
#include <strings.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>
#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

/*!
 *  \file seasonal.c
 *  \brief This file contains support for a seasonal component in a time-series model.
 */


/*!
 * \brief This function returns element Q(i,j), with i=node and j=nnode, of the precision matrix for the seasonal component
 * defined in \c def.

 \param[in] node First node
 \param[in] nnode Second node
 \param[in] def The definition of the seasonal component

 \sa \ref GMRFLib_seasonaldef_tp, \ref GMRFLib_make_seasonal_graph
 */

double GMRFLib_seasonal(int node, int nnode, void *def)
{
	int imax, imin, diff;
	double prec, val = 0.0;
	GMRFLib_seasonaldef_tp *sdef = (GMRFLib_seasonaldef_tp *) def;

	if (sdef->n < sdef->s) {
		return 0.0;
	}
	if (sdef->cyclic) {
		diff = IMIN(IABS(node - nnode), IABS(IMIN(node, nnode) - IMAX(node, nnode) + sdef->n));
		if (diff >= sdef->s) {
			return 0.0;
		} else {
			val = (double) (sdef->s - diff);
		}
	} else {
		imin = IMIN(node, nnode);
		imax = IMAX(node, nnode);
		diff = imax - imin;

		if (diff >= sdef->s) {
			return 0.0;
		} else if (imin <= (sdef->s - diff - 1)) {
			val = (imin + 1.0);
		} else if ((imin > (sdef->s - diff - 1)) && (imin < (sdef->n - sdef->s))) {
			val = (double) (sdef->s - diff);
		} else if (imin >= (sdef->n - sdef->s)) {
			val = (double) (sdef->n - diff - imin);
		} else {
			return 0.0;
		}
	}
	prec = GMRFLib_SET_PREC(sdef);
	prec *= (sdef->prec_scale ? sdef->prec_scale[0] : 1.0);

	return val * prec;
};


int GMRFLib_seasonal_scale(GMRFLib_seasonaldef_tp * def)
{
	GMRFLib_seasonaldef_tp *sdef = Calloc(1, GMRFLib_seasonaldef_tp);

	if (def->s == 1) {
		def->prec_scale = Calloc(1, double);
		def->prec_scale[0] = 1.0;
		return GMRFLib_SUCCESS;
	}

	int i, ii, j, k, n, m, s, nc;

	n = sdef->n = (def->n / def->s) * def->s;	       // make sure its a multiplum of 's'
	s = sdef->s = def->s;
	m = n / s;
	assert(m > 0);

	GMRFLib_graph_tp *graph = NULL;
	GMRFLib_make_seasonal_graph(&graph, sdef);

	GMRFLib_constr_tp *constr = NULL;
	GMRFLib_make_empty_constr(&constr);

	nc = constr->nc = s - 1;
	constr->a_matrix = Calloc(nc * n, double);

	for (ii = 0; ii < m; ii++) {
		k = ii * s;
		for (j = 0; j < nc; j++) {
			i = k + j;
			constr->a_matrix[i * nc + j] = 1.0;
		}
		i = k + nc;
		for (j = 0; j < nc; j++) {
			constr->a_matrix[i * nc + j] = -1.0;
		}
	}

	constr->e_vector = Calloc(nc, double);
	GMRFLib_prepare_constr(constr, graph, GMRFLib_TRUE);

	double *c = Calloc(n, double), eps = GMRFLib_eps(.75);
	GMRFLib_problem_tp *problem;

	for (i = 0; i < n; i++) {
		c[i] = eps;
	}

	int retval = GMRFLib_SUCCESS, ok = 0, num_try = 0, num_try_max = 100;
	GMRFLib_error_handler_tp *old_handler = GMRFLib_set_error_handler_off();

	while (!ok) {
		retval = GMRFLib_init_problem(&problem, NULL, NULL, c, NULL, graph,
					      GMRFLib_seasonal, (void *) sdef, NULL, constr, GMRFLib_NEW_PROBLEM);
		switch (retval) {
		case GMRFLib_EPOSDEF:
			for (i = 0; i < n; i++) {
				c[i] *= 10.0;
			}
			problem = NULL;
			break;
		case GMRFLib_SUCCESS:
			ok = 1;
			break;
		default:
			GMRFLib_set_error_handler(old_handler);
			GMRFLib_ERROR(retval);
			abort();
			break;
		}

		if (++num_try >= num_try_max) {
			FIXME("This should not happen. Contact developers...");
			abort();
		}
	}
	GMRFLib_set_error_handler(old_handler);
	GMRFLib_Qinv(problem, GMRFLib_QINV_DIAG);

	double sum = 0.0;
	for (i = 0; i < n; i++) {
		sum += log(*(GMRFLib_Qinv_get(problem, i, i)));
	}
	def->prec_scale = Calloc(1, double);
	def->prec_scale[0] = exp(sum / n);

	Free(c);
	Free(sdef);
	GMRFLib_free_constr(constr);
	GMRFLib_free_graph(graph);
	GMRFLib_free_problem(problem);

	return GMRFLib_SUCCESS;
}


/*!
  \brief Make the graph suitable to the seasonal model defined in \c def.

  \param[out] graph  The graph for the seasonal model.

  \param[in] def The definition of the seasonal model

  \sa GMRFLib_seasonal(), GMRFLib_seasonaldef_tp
*/
int GMRFLib_make_seasonal_graph(GMRFLib_graph_tp ** graph, GMRFLib_seasonaldef_tp * def)
{
	GMRFLib_make_linear_graph(graph, def->n, def->s - 1, def->cyclic);
	return GMRFLib_SUCCESS;
}
