
/* graph-edit.c
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
 */

/**
  \file graph-edit.c
  \brief Functions to edit graphs

  This set of functions makes it easy (or at least easier) to edit graphs in GMRFLib, for example by adding or removing nodes in
  a graph. The available functions are

  - GMRFLib_ged_init()  which from an (optional) graph, create an editable graph-object.

  This editable graph-object, can then be edited by by issuing a sequence of the following commands

  - GMRFLib_ged_add() Add an arbitrary new node or edge. By convention, if an edge between \a i and \a j is
    added, the node \a i and \a j is created at the same time.
  - GMRFLib_ged_append_graph() Append a graph to the end of the editable graph-object.
  - GMRFLib_ged_insert_graph() Insert a graph into the editable graph-object.
  - GMRFLib_ged_insert_graph2() Insert a graph into the editable graph-object.
  - GMRFLib_ged_append_node() Append a new node to the end of the editable graph-object and (optionally) tag it; See \c
    GMRFLib_ged_tag() 
  - GMRFLib_ged_remove()  Remove a node from the editable graph-object

  - GMRFLib_ged_tag() Tag a node with a special property. Current properties are GMRFLib_ged_tag_tp:GMRFLib_GED_TAG_NORMAL which
  is the default, GMRFLib_ged_tag_tp::GMRFLib_GED_TAG_GLOBAL which denote that this node is a neighbour of all other nodes that
  are present, and GMRFLib_ged_tag_tp::GMRFLib_GED_TAG_INDEP which says that this node is has no neighbours.  By convention,
  GMRFLib_ged_tag_tp::GMRFLib_GED_TAG_INDEP take precedence over GMRFLib_ged_tag_tp::GMRFLib_GED_TAG_GLOBAL.

  After the series of editing command, we can build a new \c GMRFLib_graph_tp -object, by issuing

  - GMRFLib_ged_build()

  If some nodes are removed, then the mapping between the nodes in the create graph and the nodes in the editable graph-object, is
  given in GMRFLib_graph_tp::mothergraph_idx, as usual.

  We can free the editable graph by

  - GMRFLib_ged_free()

  \par Example:

  In this small example, we add 2 global nodes to an simple linear graph

  \verbatim
    GMRFLib_graph_tp *new_graph, *g;
    GMRFLib_ged_tp *ged;

    GMRFLib_make_linear_graph(&g, n, 1, 0); // Create a linear graph

    GMRFLib_ged_init(&ged, g);  
    GMRFLib_ged_append_node(ged, GMRFLib_GED_TAG_GLOBAL);
    GMRFLib_ged_append_node(ged, GMRFLib_GED_TAG_GLOBAL);
    GMRFLib_ged_build(&new_graph, ged);    //  new_graph is the new graph
    
    GMRFLib_ged_free(ged);
  \endverbatim
  
  Using the \a _append_ -functions adds the node to the editable graph. Alternatively, we can specify node-numbers directly as
  in the following example

  \verbatim
    GMRFLib_graph_tp *new_graph, *g;
    GMRFLib_ged_tp *ged;

    GMRFLib_make_linear_graph(&g, n, 1, 0); // Create a linear graph

    GMRFLib_ged_init(&ged, g);  
    GMRFLib_ged_add(ged, n, n);
    GMRFLib_ged_add(ged, n+1, n+1);
    GMRFlib_ged_tag(ged, n, GMRFLib_GED_TAG_GLOBAL);
    GMRFlib_ged_tag(ged, n+1, GMRFLib_GED_TAG_GLOBAL);
    GMRFLib_ged_build(&new_graph, ged);    //  new_graph is the new graph
    
    GMRFLib_ged_free(ged);
  \endverbatim
  Only nodes that are added to the graph will be used, so we achieve the same new-graph by using
  \verbatim
    GMRFLib_graph_tp *new_graph, *g;
    GMRFLib_ged_tp *ged;

    GMRFLib_make_linear_graph(&g, n, 1, 0); // Create a linear graph

    GMRFLib_ged_init(&ged, g);  
    GMRFLib_ged_add(ged, n+1, n+1);
    GMRFLib_ged_add(ged, n+2, n+2);
    GMRFlib_ged_tag(ged, n+1, GMRFLib_GED_TAG_GLOBAL);
    GMRFlib_ged_tag(ged, n+2, GMRFLib_GED_TAG_GLOBAL);
    GMRFLib_ged_build(&new_graph, ged);    //  new_graph is the new graph
    
    GMRFLib_ged_free(ged);
  \endverbatim
  as node \a n will be empty  and removed when `building' the new graph using GMRFLib_ged_build().

  \par Yet another example:
  \verbinclude example-doxygen-graph-edit.txt

*/

#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

/* Pre-hg-Id: $Id: graph-edit.c,v 1.41 2009/11/04 18:24:31 hrue Exp $ */

#include <float.h>
#include <math.h>
#include <strings.h>
#include <string.h>
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"
#include "GMRFLib/hashP.h"

/**
  \brief Create an editable graph starting with (optionally) a \a graph.

  \param[out] ged  A editable graph (of type GMRFLib_ged_tp) is returned in \a *ged
  \param[in] graph An optional graph to initialise the editable graph object with

  \sa GMRFLib_ged_build(), GMRFLib_ged_free()
*/
int GMRFLib_ged_init(GMRFLib_ged_tp ** ged, GMRFLib_graph_tp * graph)
{
	/*
	 * initialise the 'ged' object starting with graph. graph can of'course be NULL. 
	 */
	*ged = Calloc(1, GMRFLib_ged_tp);
	spmatrix_init_hint(&((*ged)->Q), (mapkit_size_t) (graph ? (2 * graph->n) : 1000));
	map_ii_init(&((*ged)->tags));
	(*ged)->max_node = -1;				       /* yes, this is correct */

	if (graph) {
		GMRFLib_ged_append_graph(*ged, graph);
	}

	return GMRFLib_SUCCESS;
}

/**
  \brief Remove a node or an edge from an editable graph-object

  \param[in,out] ged The editable graph-object.
  \param[in] node First node
  \param[in] nnode Second node

  If \a node is different from \a nnode, then remove the edge between \a node and \a nnode. If \a node equals \a nnode, then
  remove node \a node itself. Note that if node \a node is removed, then so are all edges where \a node is a part of. 

  \sa GMRFLib_ged_add(), GMRFLib_ged_append_graph(), GMRFLib_ged_append_node()

*/
int GMRFLib_ged_remove(GMRFLib_ged_tp * ged, int node, int nnode)
{
	/*
	 * mark the edge between node and nnode as 'removed', or the node itself if they're equal 
	 */
	spmatrix_set(&(ged->Q), IMIN(node, nnode), IMAX(node, nnode), 0.0);

	return GMRFLib_SUCCESS;
}

/**
  \brief Add a node or an edge from an editable graph-object

  \param[in,out] ged The editable graph-object.
  \param[in] node First node
  \param[in] nnode Second node

  If \a node is different from \a nnode, then add the edge between \a node and \a nnode, and create node \a node and node \a
  nnode if they do not exists.  If \a node equals \a nnode, then add node \a node itself.

  \sa GMRFLib_ged_remove(), GMRFLib_ged_append_node(), GMRFLib_ged_append_graph()

*/
int GMRFLib_ged_add(GMRFLib_ged_tp * ged, int node, int nnode)
{
	/*
	 * add edge between node and nnode. add 'node' or 'nnode' to the set if not already present 
	 */
	if (node == nnode) {
		spmatrix_set(&(ged->Q), node, node, 1.0);
		/*
		 * workaround for internal ``bug'' in hash.c 
		 */
		spmatrix_value(&(ged->Q), node, node);
	} else {
		spmatrix_set(&(ged->Q), IMIN(node, nnode), IMAX(node, nnode), 1.0);
		spmatrix_set(&(ged->Q), node, node, 1.0);
		spmatrix_set(&(ged->Q), nnode, nnode, 1.0);
		/*
		 * workaround for internal ``bug'' in hash.c 
		 */
		spmatrix_value(&(ged->Q), IMIN(node, nnode), IMAX(node, nnode));
		spmatrix_value(&(ged->Q), node, node);
		spmatrix_value(&(ged->Q), nnode, nnode);
	}
	ged->max_node = IMAX(ged->max_node, IMAX(node, nnode));

	return GMRFLib_SUCCESS;
}

/**
  \brief Append a graph to an editable graph-object

  \param[in,out] ged the editable graph-object
  \param[in] graph The graph to be appended to the editable graph-object.

  The first node in \a graph will have node-number GMRFLib_ged_tp::max_node+1 where \a max_node is the maximum node number in
  the editable graph-object.

  \sa GMRFLib_ged_max_node(), GMRFLib_ged_append_node()
*/
int GMRFLib_ged_append_graph(GMRFLib_ged_tp * ged, GMRFLib_graph_tp * graph)
{
	GMRFLib_ged_insert_graph(ged, graph, ged->max_node + 1);

	return GMRFLib_SUCCESS;
}

/**
  \brief Insert a graph into an editable graph-object: version 1

  \param[in,out] ged the editable graph-object
  \param[in] graph The graph to be appended to the editable graph-object.
  \param[in] at_node The node where the graph will be inserted
  
  \sa GMRFLib_ged_append_graph(), GMRFLib_ged_insert_graph2()
*/
int GMRFLib_ged_insert_graph(GMRFLib_ged_tp * ged, GMRFLib_graph_tp * graph, int at_node)
{
	return GMRFLib_ged_insert_graph2(ged, graph, at_node, at_node);
}

/**
  \brief Insert a graph into an editable graph-object: version 2

  \param[in,out] ged the editable graph-object
  \param[in] graph The graph to be appended to the editable graph-object.
  \param[in] at_i_node The first node where the graph will be inserted
  \param[in] at_j_node The second node where the graph will be inserted
  
  \sa GMRFLib_ged_append_graph(), GMRFLib_ged_insert_graph()
*/
int GMRFLib_ged_insert_graph2(GMRFLib_ged_tp * ged, GMRFLib_graph_tp * graph, int at_i_node, int at_j_node)
{
	/*
	 * append graph to 'ged' 
	 */
	if (graph) {
		int i, j, jj;

		for (i = 0; i < graph->n; i++) {
			GMRFLib_ged_add(ged, i + at_i_node, i + at_j_node);
			for (jj = 0; jj < graph->nnbs[i]; jj++) {
				j = graph->nbs[i][jj];
				GMRFLib_ged_add(ged, i + at_i_node, j + at_j_node);
			}
		}
	}

	return GMRFLib_SUCCESS;
}

/**
  \brief Tag a node in the editable graph-object

  \param[in,out] ged The editable graph-object.
  \param[in] node The node number
  \param[in] tag The tag (one of GMRFLib_ged_tag_tp)

  Tag a node in the editable graph-object as either
  - GMRFLib_ged_tag_tp::GMRFLib_GED_TAG_NORMAL (this is default and usually not needed)
  - GMRFLib_ged_tag_tp::GMRFLib_GED_TAG_GLOBAL Tag node \a node as a global node, i.e. it is neigbour to all other nodes in the
  graph. 
  - GMRFLib_ged_tag_tp::GMRFLib_GED_TAG_INDEP  Tag node \a node as a singelton node, i.e. it is not neigbour to any other nodes
  in the graph. 
  
  By convention, GMRFLib_ged_tag_tp::GMRFLib_GED_TAG_INDEP take precedence over GMRFLib_ged_tag_tp::GMRFLib_GED_TAG_GLOBAL

  \sa GMRFLib_ged_node_add()

*/
int GMRFLib_ged_tag(GMRFLib_ged_tp * ged, int node, GMRFLib_ged_tag_tp tag)
{
	/*
	 * tag a node. if node does not exists, add it as well. 
	 */
	GMRFLib_ged_add(ged, node, node);
	if (tag != GMRFLib_GED_TAG_NORMAL) {
		map_ii_set(&(ged->tags), node, (int) tag);
	}

	return GMRFLib_SUCCESS;
}

/**
  \brief  Append a new node to the editable graph-object

  \param[in,out] ged The editable graph-object
  \param[in] tag The tag (one of GMRFLib_ged_tag_tp)

  Append a new node to the editable graph-object. The new node will have node-number GMRFLib_ged_tp::max_node + 1.

  \sa GMRFLib_ged_max_node(), GMRFLib_ged_add()
*/
int GMRFLib_ged_append_node(GMRFLib_ged_tp * ged, GMRFLib_ged_tag_tp tag)
{
	/*
	 * append a node to 'ged' 
	 */
	return GMRFLib_ged_tag(ged, ged->max_node + 1, tag);
}

/**
  \brief Return the maximum node-number in the editable graph-object

  \param[in] ged The editable graph-object

  This function returns the current maximum node-number in the editable graph-object, for example

  \par Example:
  \verbatim
    GMRFLib_make_linear_graph(&g, 3, 1, 0);
    GMRFLib_ged_init(&ged, g);
    GMRFLib_ged_remove(ged, 2, 2);
    printf("max_node before append %d\n", GMRFLib_ged_max_node(ged));
    GMRFLib_ged_append_node(ged, GMRFLib_GED_TAG_GLOBAL);
    printf("max_node after append %d\n", GMRFLib_ged_max_node(ged));
 \endverbatim
 will output
 \verbatim
    max_node before append 2
    max_node after append 3
 \endverbatim
 Note that 'removing a node' just tag is as 'to be removed'. The actual removal is done when building the graph issuing
  GMRFLib_ged_build(). 
*/
int GMRFLib_ged_max_node(GMRFLib_ged_tp * ged)
{
	return ged->max_node;
}

/**
  \brief Build a graph of type GMRFLib_graph_tp from the editable graph-object

  \param[out] graph The graph of type GMRFLib_graph_tp is returned as \a *graph
  \param[in]  ged   The editable graph-object

  \sa GMRFLib_graph_tp
*/
int GMRFLib_ged_build(GMRFLib_graph_tp ** graph, GMRFLib_ged_tp * ged)
{
#define NOMAP (-1)
	/*
	 * build the graph 
	 */
	GMRFLib_graph_tp *g;
	int i, j, jj, n, *nnbs, node = -1, nnode = -1, **nbs, *map, *imap, n_new;
	unsigned char *node_in_use, *indep;
	map_ii **hash;
	spmatrix_storage *sptr;
	map_ii_storage *iptr;

	/*
	 * number of (possible) nodes 
	 */
	n = ged->max_node + 1;

	/*
	 * flag nodes in use 
	 */
	node_in_use = Calloc(n, unsigned char);

	for (sptr = NULL; (sptr = spmatrix_nextptr(&(ged->Q), sptr)) != NULL;) {
		if (sptr->value != 0.0 && sptr->key.key1 == sptr->key.key2) {
			node_in_use[sptr->key.key1] = 1;
		}
	}

	/*
	 * add global nodes (if any) 
	 */
	for (j = 0, iptr = NULL; (iptr = map_ii_nextptr(&(ged->tags), iptr)) != NULL;) {
		if (iptr->value == GMRFLib_GED_TAG_GLOBAL) {
			i = iptr->key;
			if (node_in_use[i]) {
				for (jj = 0; jj < n; jj++) {
					if (node_in_use[jj] && jj != i) {
						GMRFLib_ged_add(ged, i, jj);
					}
				}
			}
		}
	}

	/*
	 * flag indep nodes (if any) 
	 */
	indep = Calloc(n, unsigned char);

	for (j = 0, iptr = NULL; (iptr = map_ii_nextptr(&(ged->tags), iptr)) != NULL;) {
		if (iptr->value == GMRFLib_GED_TAG_INDEP) {
			indep[iptr->key] = 1;
		}
	}

	imap = Calloc(n, int);
	map = Calloc(n, int);
	nbs = Calloc(n, int *);
	nnbs = Calloc(n, int);

	hash = Calloc(n, map_ii *);
	for (i = 0; i < n; i++) {
		hash[i] = Calloc(1, map_ii);
		map_ii_init(hash[i]);
	}

	/*
	 * find the mapping 
	 */
	for (i = 0; i < n; i++) {
		map[i] = imap[i] = NOMAP;
	}
	for (i = j = n_new = 0; i < n; i++) {
		if (node_in_use[i]) {
			map[i] = j;
			imap[j] = i;
			j++;
			n_new++;
		}
	}

	/*
	 * find the neighbours 
	 */
	for (sptr = NULL; (sptr = spmatrix_nextptr(&(ged->Q), sptr)) != NULL;) {
		if (sptr->value != 0.0 && sptr->key.key1 != sptr->key.key2) {

			node = IMIN(sptr->key.key1, sptr->key.key2);
			nnode = IMAX(sptr->key.key1, sptr->key.key2);
			assert(LEGAL(node, n) && LEGAL(nnode, n));

			if (node_in_use[node] && node_in_use[nnode] && !indep[node] && !indep[nnode]) {
				map_ii_set(hash[node], nnode, 1);
				map_ii_set(hash[nnode], node, 1);
				nnbs[node]++;
				nnbs[nnode]++;
			}
		}
	}

	/*
	 * make the nbs-array 
	 */
	for (i = 0; i < n; i++) {
		if (nnbs[i]) {
			nbs[i] = Calloc(nnbs[i], int);

			for (j = 0, iptr = NULL; (iptr = map_ii_nextptr(hash[i], iptr)) != NULL;) {
				nbs[i][j++] = iptr->key;
			}
		}
	}

	/*
	 * make the graph-object. copy it into a new one to make the correct memory-layout 
	 */
	for (i = 0; i < n_new; i++) {
		nnbs[i] = nnbs[imap[i]];		       /* ok, as imap[i] >= i */
		nbs[i] = nbs[imap[i]];			       /* ok, as imap[i] >= i */
	}

	for (i = 0; i < n_new; i++) {
		for (jj = 0; jj < nnbs[i]; jj++) {
			nbs[i][jj] = map[nbs[i][jj]];	       /* map to the new nodes */
		}
	}

	GMRFLib_make_empty_graph(&g);
	g->n = n_new;
	g->nnbs = nnbs;
	g->nbs = nbs;
	g->mothergraph_idx = imap;			       /* preserve the mapping */

	GMRFLib_copy_graph(graph, g);
	GMRFLib_prepare_graph(*graph);

	if (0) {
		/*
		 * validate the graph? this should not be needed, as this function will ensurethe graph should be always be consistent. 
		 * In any case, enable this for the moment until these tools are sufficiently validated. 
		 */
		GMRFLib_EWRAP0(GMRFLib_validate_graph(stderr, *graph));
	}

	/*
	 * cleanup 
	 */
	for (i = 0; i < n_new; i++) {
		Free(nbs[i]);
	}
	Free(nbs);
	Free(nnbs);
	Free(node_in_use);
	Free(map);
	Free(imap);
	Free(g);
	Free(indep);

	for (i = 0; i < n; i++) {
		map_ii_free(hash[i]);
		Free(hash[i]);
	}
	Free(hash);

	return GMRFLib_SUCCESS;
#undef NOMAP
}

/**
  \brief Free a editable graph-object

  \param[in] ged The editable graph-object to be free'd.

  \sa GMRFLib_ged_init()
 */
int GMRFLib_ged_free(GMRFLib_ged_tp * ged)
{
	spmatrix_free(&(ged->Q));
	map_ii_free(&(ged->tags));
	Free(ged);

	return GMRFLib_SUCCESS;
}

int GMRFLib_ged_print__intern(FILE * fp, GMRFLib_ged_tp * ged)
{
	spmatrix_storage *sptr;
	map_ii_storage *iptr;
	FILE *fpp;

	fpp = (fp ? fp : stdout);

	fprintf(fpp, "Contents of ged=0x%" PRIxPTR "\n", (uintptr_t) ged);
	fprintf(fpp, "\tmax_node = %1d\n", ged->max_node);
	for (sptr = NULL; (sptr = spmatrix_nextptr(&(ged->Q), sptr)) != NULL;) {
		fprintf(fpp, "\tQ[%1d, %1d] = %.1f\n", sptr->key.key1, sptr->key.key2, sptr->value);
	}

	for (iptr = NULL; (iptr = map_ii_nextptr(&(ged->tags), iptr)) != NULL;) {
		fprintf(fpp, "\ttags[%1d] = %s\n", iptr->key,
			(iptr->value == GMRFLib_GED_TAG_INDEP ? "Independent" : (iptr->value == GMRFLib_GED_TAG_GLOBAL ? "Global" : "Normal")));
	}

	fflush(fpp);

	return GMRFLib_SUCCESS;
}
