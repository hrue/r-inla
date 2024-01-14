
/* graph-edit.c
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
 */

#include <float.h>
#include <math.h>
#include <strings.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"
#include "GMRFLib/hashP.h"

#define GED_INIT 64
#define GED_GROW 1024

int GMRFLib_ged_init2(GMRFLib_ged_tp **ged, int n)
{
	*ged = Calloc(1, GMRFLib_ged_tp);
	(*ged)->n = 0;
	(*ged)->n_alloc = (IMAX(0, n) / GED_GROW + 1) * GED_GROW;
	(*ged)->Q = Calloc((*ged)->n_alloc, map_ii);
	for (int i = 0; i < (*ged)->n_alloc; i++) {
		map_ii_init_hint(&((*ged)->Q[i]), (size_t) GED_INIT);
	}
	if (n > 0) {
		GMRFLib_ged_add(*ged, n - 1, n - 1);
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_ged_init(GMRFLib_ged_tp **ged, GMRFLib_graph_tp *graph)
{
	GMRFLib_ged_init2(ged, (graph ? graph->n : 0));
	if (graph) {
		GMRFLib_ged_insert_graph(*ged, graph, 0);
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_ged_add(GMRFLib_ged_tp *ged, int node, int nnode)
{
	int imax = IMAX(node, nnode);
	int imin = IMIN(node, nnode);

	if (imax >= ged->n_alloc) {
		int np = ged->n_alloc;
		ged->n_alloc += (IMAX(1, imax + 1 - ged->n_alloc) / GED_GROW + 1) * GED_GROW;
		ged->Q = Realloc(ged->Q, ged->n_alloc, map_ii);
		for (int i = np; i < ged->n_alloc; i++) {
			map_ii_init_hint(&(ged->Q[i]), (size_t) GED_INIT);
		}
	}

	if (imin != imax) {
		map_ii_set(&(ged->Q[imin]), imax, 1);
	}
	ged->n = IMAX(ged->n, imax + 1);

	return GMRFLib_SUCCESS;
}

int GMRFLib_ged_append_graph(GMRFLib_ged_tp *ged, GMRFLib_graph_tp *graph)
{
	GMRFLib_ged_insert_graph(ged, graph, ged->n);
	return GMRFLib_SUCCESS;
}

int GMRFLib_ged_insert_graph(GMRFLib_ged_tp *ged, GMRFLib_graph_tp *graph, int at_node)
{
	return GMRFLib_ged_insert_graph2(ged, graph, at_node, at_node);
}

int GMRFLib_ged_insert_graph2(GMRFLib_ged_tp *ged, GMRFLib_graph_tp *graph, int at_i_node, int at_j_node)
{
	if (graph) {
		for (int i = 0; i < graph->n; i++) {
			int ii = i + at_i_node;
			GMRFLib_ged_add(ged, i + at_i_node, i + at_j_node);
			for (int jj = 0; jj < graph->nnbs[i]; jj++) {
				int j = graph->nbs[i][jj];
				int jjj = j + at_j_node;
				if (jjj > ii) {
					GMRFLib_ged_add(ged, ii, jjj);
				}
			}
		}
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_ged_build(GMRFLib_graph_tp **graph, GMRFLib_ged_tp *ged)
{
	GMRFLib_graph_tp *g = NULL;
	int n, *nnbs = NULL, **nbs = NULL;

	n = ged->n;
	nbs = Calloc(n, int *);
	nnbs = Calloc(n, int);

	for (int i = 0; i < n; i++) {
		for (map_ii_storage * p = NULL; (p = map_ii_nextptr(&(ged->Q[i]), p)) != NULL;) {
			if (p->key > i) {
				map_ii_set(&(ged->Q[p->key]), i, 1);
			}
		}
	}

#define CODE_BLOCK							\
	for (int i = 0; i < n; i++) {					\
		int j;							\
		map_ii_storage *p;					\
		for (j = 0, p = NULL; (p = map_ii_nextptr(&(ged->Q[i]), p)) != NULL;) {	\
			j++;						\
		}							\
		nnbs[i] = j;						\
		if (nnbs[i]) {						\
			nbs[i] = Calloc(nnbs[i], int);			\
			for (j = 0, p = NULL; (p = map_ii_nextptr(&(ged->Q[i]), p)) != NULL; j++) { \
				nbs[i][j] = p->key;			\
			}						\
			assert(j == nnbs[i]);				\
		} else {						\
			nbs[i] = NULL;					\
		}							\
	}
	RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 0, 0);
#undef CODE_BLOCK

	GMRFLib_graph_mk_empty(&g);
	g->n = n;
	g->nnbs = nnbs;
	g->nbs = nbs;

	GMRFLib_graph_prepare(g);
	GMRFLib_graph_duplicate(graph, g);

	for (int i = 0; i < n; i++) {
		Free(nbs[i]);
	}
	Free(nbs);
	Free(nnbs);
	Free(g);

	return GMRFLib_SUCCESS;
}

int GMRFLib_ged_free(GMRFLib_ged_tp *ged)
{
	for (int i = 0; i < ged->n_alloc; i++) {
		map_ii_free(&(ged->Q[i]));
	}
	Free(ged->Q);
	Free(ged);

	return GMRFLib_SUCCESS;
}
