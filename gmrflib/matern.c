
/* matern.c
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

/*!
  \file matern.c
  \brief Functions to define Matern models
*/

#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

/* Pre-hg-Id: $Id: matern.c,v 1.10 2008/11/11 18:46:46 hrue Exp $ */

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

/*!
  \brief This function returns element Q(i,j), with i=node and j=nnode, of the precision matrix for the Matern defined in \c def, which is of type \c
  GMRFLib_matern2ddef_tp.

  \param[in] node First node
  \param[in] nnode Second node
  \param[in] def The definition of the Matern-model (2d and regular lattice)

  \sa \ref GMRFLib_matern2ddef_tp, \ref GMRFLib_make_matern2d_graph
*/
double GMRFLib_matern2d(int node, int nnode, void *def)
{
	double prec, range, kappa, std_variance, a, val = 0;
	GMRFLib_matern2ddef_tp *arg = NULL;

	arg = (GMRFLib_matern2ddef_tp *) def;
	prec = GMRFLib_SET_PREC(arg);
	range = GMRFLib_SET_RANGE(arg);
	kappa = 2.0 * sqrt(2.0 * arg->nu) / range;
	a = 4.0 + SQR(kappa);

	if (arg->nu > 0) {
		std_variance = pow(kappa, -2.0 * arg->nu) / (4.0 * M_PI * arg->nu);
	} else {
		FIXME1("nu = 0 is not well-defined for the moment; sorry about that...");
		GMRFLib_ASSERT_RETVAL(0 == 1, GMRFLib_EPARAMETER, 0.0);
	}

	if (node == nnode) {
		switch (arg->nu) {
		case 0:
			val = a;
			break;
		case 1:
			val = 4.0 + SQR(a);
			break;
		case 2:
			val = a * (SQR(a) + 12.0);
			break;
		case 3:
			val = SQR(SQR(a) + 6.0) + 12 * SQR(a);
			break;
		default:
			GMRFLib_ASSERT_RETVAL(arg->nu >= 0 && arg->nu <= 3, GMRFLib_EPARAMETER, 0.0);
		}
		return val * prec * std_variance;
	}

	int irow, icol, jrow, jcol, drow, dcol, dmin, dmax;

	GMRFLib_node2lattice(node, &irow, &icol, arg->nrow, arg->ncol);
	GMRFLib_node2lattice(nnode, &jrow, &jcol, arg->nrow, arg->ncol);

	drow = IABS(irow - jrow);
	dcol = IABS(icol - jcol);
	if (arg->cyclic) {
		drow = IMIN(drow, arg->nrow - drow);
		dcol = IMIN(dcol, arg->ncol - dcol);
	}

	dmax = IMAX(drow, dcol);
	dmin = IMIN(drow, dcol);

	switch (arg->nu) {
	case 0:
		if (dmin == 0 && dmax == 1) {
			val = -1.0;
		}
		break;
	case 1:
		switch (dmin) {
		case 0:
			switch (dmax) {
			case 1:
				val = -2 * a;
				break;
			case 2:
				val = 1.0;
				break;
			default:
				val = 0.0;
			}
			break;
		case 1:
			val = (dmax == 1 ? 2.0 : 0.0);
		}
		break;
	case 2:
		switch (dmin) {
		case 0:
			switch (dmax) {
			case 1:
				val = -3.0 * (SQR(a) + 3.0);
				break;
			case 2:
				val = 3.0 * a;
				break;
			case 3:
				val = -1.0;
			}
			break;
		case 1:
			switch (dmax) {
			case 1:
				val = 6 * a;
				break;
			case 2:
				val = -3.0;
			}
			break;
		}
		break;
	case 3:
		switch (dmin) {
		case 0:
			switch (dmax) {
			case 1:
				val = -4.0 * a * (SQR(a) + 9.0);
				break;
			case 2:
				val = 2 * (3 * SQR(a) + 8.0);
				break;
			case 3:
				val = -4 * a;
				break;
			case 4:
				val = 1.0;
			}
			break;
		case 1:
			switch (dmax) {
			case 1:
				val = 12.0 * (SQR(a) + 2.0);
				break;
			case 2:
				val = -12 * a;
				break;
			case 3:
				val = 4.0;
			}
			break;
		case 2:
			val = (dmax == 2 ? 6.0 : 0.0);
			break;
		}
		break;
	default:
		GMRFLib_ASSERT_RETVAL(arg->nu >= 0 && arg->nu <= 3, GMRFLib_EPARAMETER, 0.0);
	}
	return val * prec * std_variance;
}

/*!
  \brief Creates the graph for the Martern 2d model in \c def

  \param[out] graph   The ppointer to the graph
  \param[in] def  The definition of the Matern-model (2d and regular lattice)

  \sa  \ref GMRFLib_matern2d, \ref GMRFLib_matern2ddef_tp
*/
int GMRFLib_make_matern2d_graph(GMRFLib_graph_tp ** graph, GMRFLib_matern2ddef_tp * def)
{
	GMRFLib_graph_tp *g = NULL;

	GMRFLib_make_lattice_graph(&g, def->nrow, def->ncol, def->nu + 1, def->nu + 1, def->cyclic);
	GMRFLib_prune_graph(graph, g, (GMRFLib_Qfunc_tp *) GMRFLib_matern2d, (void *) def);
	GMRFLib_free_graph(g);

	return GMRFLib_SUCCESS;
}

/*
  Example for manual
 */

/*!

  \page ex_matern A simple example of a Matern-2d model
  
\par Program code:

\verbinclude example-doxygen-matern2d.txt

*/
