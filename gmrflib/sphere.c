
/* sphere.c
 * 
 * Copyright (C) 2005-2021 Havard Rue
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

#ifndef GITCOMMIT
#define GITCOMMIT
#endif
static const char GitID[] = "file: " __FILE__ "  " GITCOMMIT;

#include <string.h>
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

int GMRFLib_make_spherical_igmrf(GMRFLib_sphere_tp ** sphere, int resolution, double *UNUSED(scale), double *UNUSED(log_scale),
				 double **UNUSED(log_scale_omp))
{
	char *fnm = NULL, *p;
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
	GMRFLib_EWRAP0(GMRFLib_graph_read(&((*sphere)->graph), p));
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
	FIXME("FIX THIS");
	abort();

	// GMRFLib_EWRAP0(GMRFLib_tabulate_Qfunc_from_file_OLD(&tab, (*sphere)->graph, p, scale, log_scale, log_scale_omp));
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

int GMRFLib_free_spherical_igmrf(GMRFLib_sphere_tp * sphere)
{
	GMRFLib_tabulate_Qfunc_tp *tab;

	if (!sphere) {
		return GMRFLib_SUCCESS;
	}

	GMRFLib_graph_free(sphere->graph);
	Free(sphere->area);
	Free(sphere->coord);

	tab = Calloc(1, GMRFLib_tabulate_Qfunc_tp);
	tab->Qfunc = sphere->Qfunc;
	tab->Qfunc_arg = sphere->Qfunc_arg;
	GMRFLib_free_tabulate_Qfunc(tab);

	Free(sphere);

	return GMRFLib_SUCCESS;
}
