
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
 *
 * This code was initially contributed by Sara Martino (NTNU), and adapted to GMRFLib by H. Rue.
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

	return val * prec;
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
