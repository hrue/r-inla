
/* graph-edit.h
 * 
 * Copyright (C) 2006 Havard Rue
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
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif

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
    typedef enum {

	/**
	 *  Normal node tag 
	 */
	GMRFLib_GED_TAG_NORMAL,

	/**
	 *  Global node tag 
	 */
	GMRFLib_GED_TAG_GLOBAL,

	/**
	 *  Independent node tag 
	 */
	GMRFLib_GED_TAG_INDEP
} GMRFLib_ged_tag_tp;

typedef struct {
	int max_node;					       /* keep track of the max{all nodes} */
	spmatrix Q;					       /* the edges with value 1 for present and 0 for absent */
	map_ii tags;					       /* a list of tagged nodes and their tags */
} GMRFLib_ged_tp;

int GMRFLib_ged_add(GMRFLib_ged_tp * ged, int node, int nnode);
int GMRFLib_ged_append_graph(GMRFLib_ged_tp * ged, GMRFLib_graph_tp * graph);
int GMRFLib_ged_append_node(GMRFLib_ged_tp * ged, GMRFLib_ged_tag_tp);
int GMRFLib_ged_build(GMRFLib_graph_tp ** graph, GMRFLib_ged_tp * ged);
int GMRFLib_ged_free(GMRFLib_ged_tp * ged);
int GMRFLib_ged_init(GMRFLib_ged_tp ** ged, GMRFLib_graph_tp * graph);
int GMRFLib_ged_insert_graph(GMRFLib_ged_tp * ged, GMRFLib_graph_tp * graph, int at_node);
int GMRFLib_ged_insert_graph2(GMRFLib_ged_tp * ged, GMRFLib_graph_tp * graph, int at_i_node, int at_j_node);
int GMRFLib_ged_max_node(GMRFLib_ged_tp * ged);
int GMRFLib_ged_print__intern(FILE * fp, GMRFLib_ged_tp * ged);
int GMRFLib_ged_remove(GMRFLib_ged_tp * ged, int node, int nnode);
int GMRFLib_ged_tag(GMRFLib_ged_tp * ged, int node, GMRFLib_ged_tag_tp);

__END_DECLS
#endif
