
/* graph-edit.c
 * 
 * Copyright (C) 2006-2021 Havard Rue
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

int GMRFLib_ged_remove(GMRFLib_ged_tp * ged, int node, int nnode)
{
	/*
	 * mark the edge between node and nnode as 'removed', or the node itself if they're equal 
	 */
	spmatrix_set(&(ged->Q), IMIN(node, nnode), IMAX(node, nnode), 0.0);

	return GMRFLib_SUCCESS;
}

int GMRFLib_ged_add(GMRFLib_ged_tp * ged, int node, int nnode)
{
	/*
	 * add edge between node and nnode. add 'node' or 'nnode' to the set if not already present 
	 */
	if (node == nnode) {
		spmatrix_set(&(ged->Q), node, node, 1.0);
	} else {
		if (!spmatrix_value(&(ged->Q), IMIN(node, nnode), IMAX(node, nnode))) {
			spmatrix_set(&(ged->Q), IMIN(node, nnode), IMAX(node, nnode), 1.0);
			spmatrix_set(&(ged->Q), node, node, 1.0);
			spmatrix_set(&(ged->Q), nnode, nnode, 1.0);
		}
	}
	ged->max_node = IMAX(ged->max_node, IMAX(node, nnode));

	return GMRFLib_SUCCESS;
}

int GMRFLib_ged_append_graph(GMRFLib_ged_tp * ged, GMRFLib_graph_tp * graph)
{
	GMRFLib_ged_insert_graph(ged, graph, ged->max_node + 1);

	return GMRFLib_SUCCESS;
}

int GMRFLib_ged_insert_graph(GMRFLib_ged_tp * ged, GMRFLib_graph_tp * graph, int at_node)
{
	return GMRFLib_ged_insert_graph2(ged, graph, at_node, at_node);
}

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

int GMRFLib_ged_append_node(GMRFLib_ged_tp * ged, GMRFLib_ged_tag_tp tag)
{
	/*
	 * append a node to 'ged' 
	 */
	return GMRFLib_ged_tag(ged, ged->max_node + 1, tag);
}

int GMRFLib_ged_max_node(GMRFLib_ged_tp * ged)
{
	return ged->max_node;
}

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

	GMRFLib_graph_mk_empty(&g);
	g->n = n_new;
	g->nnbs = nnbs;
	g->nbs = nbs;
	g->mothergraph_idx = imap;			       /* preserve the mapping */

	GMRFLib_graph_duplicate(graph, g);

	if (0) {
		/*
		 * validate the graph? this should not be needed, as this function will ensurethe graph should be always be consistent. 
		 * In any case, enable this for the moment until these tools are sufficiently validated. 
		 */
		GMRFLib_EWRAP0(GMRFLib_graph_validate(stderr, *graph));
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
