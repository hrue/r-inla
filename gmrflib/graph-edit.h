
/* graph-edit.h
 * 
 * Copyright (C) 2006-2024 Havard Rue
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
 *
 */

/*!
  \file graph-edit.h
  \brief Typedefs used to edit graphs
*/

#ifndef __GMRFLib_GRAPH_EDIT_H__
#define __GMRFLib_GRAPH_EDIT_H__

#include <stdlib.h>
#include <stddef.h>
#include <math.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS					       /* empty */
#define __END_DECLS					       /* empty */
#endif

__BEGIN_DECLS

/*!
   \brief The various values for a tag
*/

/*
 */
    typedef struct {
	int n;
	int n_alloc;
	map_ii *Q;
} GMRFLib_ged_tp;

int GMRFLib_ged_add(GMRFLib_ged_tp * ged, int node, int nnode);
int GMRFLib_ged_append_graph(GMRFLib_ged_tp * ged, GMRFLib_graph_tp * graph);
int GMRFLib_ged_build(GMRFLib_graph_tp ** graph, GMRFLib_ged_tp * ged);
int GMRFLib_ged_free(GMRFLib_ged_tp * ged);
int GMRFLib_ged_init(GMRFLib_ged_tp ** ged, GMRFLib_graph_tp * graph);
int GMRFLib_ged_init2(GMRFLib_ged_tp ** ged, int max_node);
int GMRFLib_ged_insert_graph(GMRFLib_ged_tp * ged, GMRFLib_graph_tp * graph, int at_node);
int GMRFLib_ged_insert_graph2(GMRFLib_ged_tp * ged, GMRFLib_graph_tp * graph, int at_i_node, int at_j_node);

__END_DECLS
#endif
