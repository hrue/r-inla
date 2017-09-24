
/* sphere.c
 * 
 * Copyright (C) 2005-2006 Havard Rue
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
  \file sphere.c
  \brief Functions for defining spherical IGMRFs

*/

#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

/* Pre-hg-Id: $Id: sphere.c,v 1.29 2008/09/10 09:07:08 hrue Exp $ */

#include <string.h>
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

/*!

\brief This function defines a spherical IGMRF of a given resoultion
   
This function defines a spherical IGMRF of a given \c resoultion and scaling \c scale. The
resolution must be 0, 1, 2, 3, 4, 5 or 6, giving a IGMRF with size 12, 42, 162, 642, 2562, 10242 or
40962, respectively. The argument \c scale, gives a pointer to a possible scaling of the elements in
the precision matrix. If \c scale is \c NULL, then the scaling is 1.

\param[in,out] sphere At output, <em>(*sphere)</em> is a pointer to a \c GMRFLib_sphere_tp -object,
  initialised and defined according to the problem specification.
  
\param[in] resolution The resolution of the spherical IGMRF: either 0, 1, 2, 3, 4, 5 or 6.

\param[in] scale A pointer to the scaling of the precision matrix. If \c scale is \c NULL, then \c log_scale is used. 

\param[in] log_scale A pointer to the log of the scaling of the precision matrix. If \c log_scale is \c NULL, then \c log_scale_omp is used.

\param[in] log_scale_omp A ppointer to the log of the scaling of the precision matrix. If \c scale and \c log_scale is \c NULL, then \c
log_scale_omp[GMRFLib_thread_id] is used, otherwise the \c scale is set to 1.

  \par Example

  Here is an example using spherical IGMRFs using \c GMRFLib_make_spherical_igmrf()
  
  \par Program code:

  \verbinclude example-doxygen-sphere.txt

*/
int GMRFLib_make_spherical_igmrf(GMRFLib_sphere_tp ** sphere, int resolution, double *scale, double *log_scale, double **log_scale_omp)
{
	char *fnm, *p;
	int i, n;
	GMRFLib_tabulate_Qfunc_tp *tab = NULL;
	GMRFLib_io_tp *io = NULL;

	GMRFLib_ASSERT(resolution >= 0 && resolution <= 6, GMRFLib_EPARAMETER);
	*sphere = Calloc(1, GMRFLib_sphere_tp);

	GMRFLib_EWRAP0(GMRFLib_sprintf(&fnm, "sphere-%1d.graph", resolution));
	GMRFLib_io_find_file_in_path(&p, fnm, 0);
	if (!p) {
		Free(fnm);
		GMRFLib_EWRAP0(GMRFLib_sprintf(&fnm, "sphere-%1d.graph.gz", resolution));
		GMRFLib_EWRAP0(GMRFLib_io_find_file_in_path(&p, fnm, 1));
	}
	GMRFLib_EWRAP0(GMRFLib_read_graph(&((*sphere)->graph), p));
	Free(fnm);
	Free(p);

	n = (*sphere)->graph->n;

	GMRFLib_EWRAP0(GMRFLib_sprintf(&fnm, "sphere-%1d.Q", resolution));
	GMRFLib_io_find_file_in_path(&p, fnm, 0);
	if (!p) {
		Free(fnm);
		GMRFLib_EWRAP0(GMRFLib_sprintf(&fnm, "sphere-%1d.Q.gz", resolution));
		GMRFLib_EWRAP0(GMRFLib_io_find_file_in_path(&p, fnm, 1));
	}

	/*
	 * TODO: we can use the new implenetation of _Qfunc_from_file() which build the graph simultanously. 
	 */
	GMRFLib_EWRAP0(GMRFLib_tabulate_Qfunc_from_file_OLD(&tab, (*sphere)->graph, p, scale, log_scale, log_scale_omp));
	Free(fnm);
	Free(p);

	(*sphere)->Qfunc = tab->Qfunc;
	(*sphere)->Qfunc_arg = tab->Qfunc_arg;
	Free(tab);					       /* yes */

	GMRFLib_EWRAP0(GMRFLib_sprintf(&fnm, "sphere-%1d.area", resolution));
	GMRFLib_io_find_file_in_path(&p, fnm, 0);
	if (!p) {
		Free(fnm);
		GMRFLib_EWRAP0(GMRFLib_sprintf(&fnm, "sphere-%1d.area.gz", resolution));
		GMRFLib_EWRAP0(GMRFLib_io_find_file_in_path(&p, fnm, 1));
	}
	(*sphere)->area = Calloc(n, double);

	GMRFLib_EWRAP0(GMRFLib_io_open(&io, p, "r"));
	for (i = 0; i < n; i++) {
		GMRFLib_EWRAP0(GMRFLib_io_read_next(io, &((*sphere)->area[i]), "%lf"));
	}
	GMRFLib_EWRAP0(GMRFLib_io_close(io));
	Free(fnm);
	Free(p);

	GMRFLib_EWRAP0(GMRFLib_sprintf(&fnm, "sphere-%1d.coord", resolution));
	GMRFLib_io_find_file_in_path(&p, fnm, 0);
	if (!p) {
		Free(fnm);
		GMRFLib_EWRAP0(GMRFLib_sprintf(&fnm, "sphere-%1d.coord.gz", resolution));
		GMRFLib_EWRAP0(GMRFLib_io_find_file_in_path(&p, fnm, 1));
	}
	GMRFLib_EWRAP0(GMRFLib_io_open(&io, p, "r"));

	(*sphere)->coord = Calloc(n, GMRFLib_sphere_coord_tp);
	for (i = 0; i < n; i++) {
		GMRFLib_EWRAP0(GMRFLib_io_read_next(io, &((*sphere)->coord[i].x), "%lf"));
		GMRFLib_EWRAP0(GMRFLib_io_read_next(io, &((*sphere)->coord[i].y), "%lf"));
		GMRFLib_EWRAP0(GMRFLib_io_read_next(io, &((*sphere)->coord[i].z), "%lf"));
	}
	GMRFLib_EWRAP0(GMRFLib_io_close(io));
	Free(fnm);
	Free(p);

	return GMRFLib_SUCCESS;
}

/*!
  \brief Free a spherical iGMRF

  This function will free an \c GMRFLib_sphere_tp object, created by \c
  GMRFLib_make_spherical_igmrf()
*/
int GMRFLib_free_spherical_igmrf(GMRFLib_sphere_tp * sphere)
{
	GMRFLib_tabulate_Qfunc_tp *tab;

	if (!sphere) {
		return GMRFLib_SUCCESS;
	}

	GMRFLib_free_graph(sphere->graph);
	Free(sphere->area);
	Free(sphere->coord);

	tab = Calloc(1, GMRFLib_tabulate_Qfunc_tp);
	tab->Qfunc = sphere->Qfunc;
	tab->Qfunc_arg = sphere->Qfunc_arg;
	GMRFLib_free_tabulate_Qfunc(tab);

	Free(sphere);

	return GMRFLib_SUCCESS;
}
