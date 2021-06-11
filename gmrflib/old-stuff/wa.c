
/* wa.c
 * 
 * Copyright (C) 2001-2021 Havard Rue
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

#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"
#include "GMRFLib/hashP.h"

#ifndef GITCOMMIT
#define GITCOMMIT
#endif
static const char GitID[] = "file: " __FILE__ "  " GITCOMMIT;

int GMRFLib_free_wa_problem(GMRFLib_wa_problem_tp * wa_problem)
{
	if (wa_problem) {
		int i;
		GMRFLib_waQfunc_arg_tp *arg = NULL;

		arg = (GMRFLib_waQfunc_arg_tp *) wa_problem->Qfunc_arg;

		if (arg->neigb_info) {
			for (i = 0; i < arg->n_neigb_info; i++) {
				if (arg->neigb_info[i]) {
					Free(arg->neigb_info[i]->node_list);
					Free(arg->neigb_info[i]);
				}
			}
			Free(arg->neigb_info);
			spmatrix_free(&(arg->neigb_idx_hash));
		}

		GMRFLib_graph_free(wa_problem->graph);
		GMRFLib_graph_free(arg->waQgraph);
		GMRFLib_graph_free(arg->wagraph);

		Free(wa_problem->Qfunc_arg);
		Free(wa_problem);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_init_wa_problem(GMRFLib_wa_problem_tp ** wa_problem, GMRFLib_graph_tp * wagraph, GMRFLib_Qfunc_tp * wafunc, void *wafunc_arg)
{
	/*
	 * wagraph is the graph defining
	 * 
	 * (*) \sum_i (w_ii x_i - \sum_{j~i} w_ij x_j)^2 = x^T Q x
	 * 
	 * where 'wafunc(i,j)' returns w_ij with arguments wafunc_arg. this routine compute the corresponding graph defining Q, 
	 * a function returning Q_ij and its arguments. 
	 */

	int i, j, k, kk, jj, jjj, kkk, n, *memsiz = NULL, nnb, *hold = NULL, indx;
	GMRFLib_graph_tp *graph = NULL, *ngraph = NULL;
	GMRFLib_waQfunc_arg_tp *wa_arg = NULL;

	GMRFLib_graph_mk_empty(&graph);
	n = wagraph->n;
	graph->n = n;
	graph->nnbs = Calloc(n, int);
	graph->nbs = Calloc(n, int *);
	memsiz = Calloc(n, int);

	for (i = 0; i < n; i++) {
		graph->nnbs[i] = wagraph->nnbs[i];
		memsiz[i] = IMAX(1, 2 * graph->nnbs[i]);
		graph->nbs[i] = Calloc(memsiz[i], int);

		if (graph->nnbs[i])
			Memcpy(graph->nbs[i], wagraph->nbs[i], graph->nnbs[i] * sizeof(int));
	}

	/*
	 * build the new graph 
	 */
	if (0) {
		/*
		 * old code, slow algorithm..... 
		 */
		for (i = 0; i < n; i++)
			for (j = 0; j < wagraph->nnbs[i]; j++)
				for (jj = j + 1; jj < wagraph->nnbs[i]; jj++) {
					k = wagraph->nbs[i][j];
					kk = wagraph->nbs[i][jj];

					if (!GMRFLib_graph_is_nb(k, kk, graph)) {
						/*
						 * add, must do it symmetrically! must also sort the neighbours, since the
						 * 'is_neig' function assume they are sorted. 
						 */
						graph->nnbs[k]++;
						if (graph->nnbs[k] > memsiz[k]) {
							memsiz[k] *= 2;
							graph->nbs[k] = Realloc(graph->nbs[k], memsiz[k], int);
						}
						graph->nbs[k][graph->nnbs[k] - 1] = kk;
						qsort(graph->nbs[k], (size_t) graph->nnbs[k], sizeof(int), GMRFLib_icmp);

						graph->nnbs[kk]++;
						if (graph->nnbs[kk] > memsiz[kk]) {
							memsiz[kk] *= 2;
							graph->nbs[kk] = Realloc(graph->nbs[kk], memsiz[kk], int);
						}
						graph->nbs[kk][graph->nnbs[kk] - 1] = k;
						qsort(graph->nbs[kk], (size_t) graph->nnbs[kk], sizeof(int), GMRFLib_icmp);
					}
				}
	} else {
		/*
		 * new version, runs faster, use slightly more memory, but not that much 
		 */
		if (0)
			printf("\n\n%s:%1d:NEW CODE HERE, NOT PROPERLY TESTED!!!\n\n", __FILE__, __LINE__);
		for (i = 0; i < n; i++)
			for (j = 0; j < wagraph->nnbs[i]; j++) {
				k = wagraph->nbs[i][j];
				for (jj = j + 1; jj < wagraph->nnbs[i]; jj++) {
					kk = wagraph->nbs[i][jj];

					if (graph->nnbs[k] >= memsiz[k]) {
						qsort(graph->nbs[k], (size_t) graph->nnbs[k], sizeof(int), GMRFLib_icmp);
						for (jjj = 1, kkk = 0; jjj < graph->nnbs[k]; jjj++)
							if (graph->nbs[k][jjj] != graph->nbs[k][kkk])
								graph->nbs[k][++kkk] = graph->nbs[k][jjj];
						graph->nnbs[k] = kkk + 1;
					}
					if (graph->nnbs[k] >= memsiz[k]) {
						memsiz[k] *= 2;
						graph->nbs[k] = Realloc(graph->nbs[k], memsiz[k], int);
					}
					graph->nbs[k][graph->nnbs[k]++] = kk;

					if (graph->nnbs[kk] >= memsiz[kk]) {
						qsort(graph->nbs[kk], (size_t) graph->nnbs[kk], sizeof(int), GMRFLib_icmp);
						for (jjj = 1, kkk = 0; jjj < graph->nnbs[kk]; jjj++)
							if (graph->nbs[kk][jjj] != graph->nbs[kk][kkk])
								graph->nbs[kk][++kkk] = graph->nbs[kk][jjj];
						graph->nnbs[kk] = kkk + 1;
					}
					if (graph->nnbs[kk] >= memsiz[kk]) {
						memsiz[kk] *= 2;
						graph->nbs[kk] = Realloc(graph->nbs[kk], memsiz[kk], int);
					}
					graph->nbs[kk][graph->nnbs[kk]++] = k;
				}
			}
		for (i = 0; i < n; i++)			       /* this step is needed */
			if (graph->nnbs[i]) {
				qsort(graph->nbs[i], (size_t) graph->nnbs[i], sizeof(int), GMRFLib_icmp);
				for (j = 1, k = 0; j < graph->nnbs[i]; j++)
					if (graph->nbs[i][j] != graph->nbs[i][k])
						graph->nbs[i][++k] = graph->nbs[i][j];
				graph->nnbs[i] = k + 1;
			}
	}

	for (i = 0, nnb = 0; i < n; i++)
		nnb += graph->nnbs[i];
	if (nnb) {
		hold = Calloc(nnb, int);		       /* use a linear storage */
	}
	for (i = 0, indx = 0; i < n; i++) {
		if (graph->nnbs[i]) {
			Memcpy(&hold[indx], graph->nbs[i], graph->nnbs[i] * sizeof(int));
			Free(graph->nbs[i]);
			graph->nbs[i] = &hold[indx];
		} else
			Free(graph->nbs[i]);

		indx += graph->nnbs[i];
	}

	Free(memsiz);
	GMRFLib_graph_prepare(graph);

	/*
	 * setup the new types 
	 */
	wa_arg = Calloc(1, GMRFLib_waQfunc_arg_tp);
	wa_arg->waQgraph = graph;
	wa_arg->waQfunc = wafunc;
	wa_arg->waQfunc_arg = wafunc_arg;
	GMRFLib_graph_duplicate(&(wa_arg->wagraph), wagraph);  /* yes, make a copy! */

	if (GMRFLib_use_wa_table_lookup) {
		/*
		 * build hash table for idx = (i,j) -> neigb_info[idx] 
		 */
		double idx = 0.0;
		int nelm;

		GMRFLib_graph_nnodes(&nelm, graph);
		spmatrix_init_hint(&(wa_arg->neigb_idx_hash), (mapkit_size_t) nelm);	/* give a hint of the size */
		wa_arg->neigb_idx_hash.alwaysdefault = 0;
		for (i = 0; i < graph->n; i++) {
			spmatrix_set(&(wa_arg->neigb_idx_hash), i, i, idx);
			idx++;
			for (j = 0; j < graph->nnbs[i]; j++) {
				spmatrix_set(&(wa_arg->neigb_idx_hash), i, graph->nbs[i][j], idx);
				idx++;
			}
		}
		wa_arg->n_neigb_info = (int) idx;
		/*
		 * hold list of neighbors and common neigbhbors, to gain some real speedup.... 
		 */
		wa_arg->neigb_info = Calloc(wa_arg->n_neigb_info, GMRFLib_node_list_tp *);
	} else
		wa_arg->neigb_info = NULL;

	*wa_problem = Calloc(1, GMRFLib_wa_problem_tp);
	(*wa_problem)->graph = graph;
	(*wa_problem)->Qfunc = GMRFLib_waQfunc;
	(*wa_problem)->Qfunc_arg = (void *) wa_arg;

	/*
	 * prune graph. note that waQgraph is the original unpruned graph which is needed for correct computing of Q_ij after
	 * pruning. 
	 */
	GMRFLib_graph_prune(&ngraph, (*wa_problem)->graph, (*wa_problem)->Qfunc, (*wa_problem)->Qfunc_arg);
	(*wa_problem)->graph = ngraph;

	return GMRFLib_SUCCESS;
}

double GMRFLib_waQfunc(int node, int nnode, double *UNUSED(values), void *arg)
{
	if (node >= 0 && nnode < 0) {
		return NAN;
	}

	int j, jj, k, count, idx;
	double val, tmp, *ptr = NULL;
	GMRFLib_waQfunc_arg_tp *args = NULL;

	args = (GMRFLib_waQfunc_arg_tp *) arg;

	if (GMRFLib_use_wa_table_lookup) {		       /* use internal tables */
		if ((ptr = spmatrix_ptr(&(args->neigb_idx_hash), node, nnode))) {
			idx = (int) *ptr;
		} else {
			GMRFLib_ASSERT(ptr != NULL, GMRFLib_ESNH);
		}

		if (node == nnode) {
			if (!args->neigb_info[idx]) {
				/*
				 * then build up information 
				 */
				for (j = 0, count = 0; j < args->waQgraph->nnbs[node]; j++) {
					k = args->waQgraph->nbs[node][j];
					if (GMRFLib_graph_is_nb(k, node, args->wagraph))
						count++;
				}
				args->neigb_info[idx] = Calloc(1, GMRFLib_node_list_tp);
				if (count)
					args->neigb_info[idx]->node_list = Calloc(count, int);

				args->neigb_info[idx]->n_nodes = count;
				for (j = 0, jj = 0; j < args->waQgraph->nnbs[node]; j++) {
					k = args->waQgraph->nbs[node][j];
					if (GMRFLib_graph_is_nb(k, node, args->wagraph))
						args->neigb_info[idx]->node_list[jj++] = k;
				}
			}

			tmp = (*args->waQfunc) (node, node, NULL, args->waQfunc_arg);
			val = SQR(tmp);
			for (j = 0; j < args->neigb_info[idx]->n_nodes; j++) {
				k = args->neigb_info[idx]->node_list[j];
				tmp = (*args->waQfunc) (k, node, NULL, args->waQfunc_arg);
				val += SQR(tmp);
			}
			return val;
		} else {
			if (!args->neigb_info[idx]) {
				/*
				 * then build up information 
				 */

				args->neigb_info[idx] = Calloc(1, GMRFLib_node_list_tp);

				args->neigb_info[idx]->node_nnode = GMRFLib_graph_is_nb(node, nnode, args->wagraph);
				args->neigb_info[idx]->nnode_node = GMRFLib_graph_is_nb(nnode, node, args->wagraph);

				for (j = 0, count = 0; j < args->waQgraph->nnbs[node]; j++) {
					k = args->waQgraph->nbs[node][j];
					if (GMRFLib_graph_is_nb(k, node, args->wagraph)
					    && GMRFLib_graph_is_nb(k, nnode, args->wagraph))
						count++;
				}
				if (count)
					args->neigb_info[idx]->node_list = Calloc(count, int);

				args->neigb_info[idx]->n_nodes = count;
				for (j = 0, jj = 0; j < args->waQgraph->nnbs[node]; j++) {
					k = args->waQgraph->nbs[node][j];
					if (GMRFLib_graph_is_nb(k, node, args->wagraph)
					    && GMRFLib_graph_is_nb(k, nnode, args->wagraph))
						args->neigb_info[idx]->node_list[jj++] = k;
				}
			}

			val = 0.0;
			if (args->neigb_info[idx]->node_nnode)
				val -= (*args->waQfunc) (node, node, NULL, args->waQfunc_arg)
				    * (*args->waQfunc) (node, nnode, NULL, args->waQfunc_arg);
			if (args->neigb_info[idx]->nnode_node)
				val -= (*args->waQfunc) (nnode, nnode, NULL, args->waQfunc_arg)
				    * (*args->waQfunc) (nnode, node, NULL, args->waQfunc_arg);

			for (j = 0; j < args->neigb_info[idx]->n_nodes; j++) {
				k = args->neigb_info[idx]->node_list[j];
				val += (*args->waQfunc) (k, node, NULL, args->waQfunc_arg)
				    * (*args->waQfunc) (k, nnode, NULL, args->waQfunc_arg);
			}
			return val;
		}
	} else {					       /* plain version */

		if (node == nnode) {
			tmp = (*args->waQfunc) (node, node, NULL, args->waQfunc_arg);
			val = SQR(tmp);

			for (j = 0; j < args->waQgraph->nnbs[node]; j++) {
				k = args->waQgraph->nbs[node][j];
				if (GMRFLib_graph_is_nb(node, k, args->wagraph)) {
					tmp = (*args->waQfunc) (k, node, NULL, args->waQfunc_arg);
					val += SQR(tmp);
				}
			}
			return val;
		} else {
			val = 0.0;
			if (GMRFLib_graph_is_nb(node, nnode, args->wagraph))
				val -= (*args->waQfunc) (node, node, NULL, args->waQfunc_arg)
				    * (*args->waQfunc) (node, nnode, NULL, args->waQfunc_arg);
			if (GMRFLib_graph_is_nb(nnode, node, args->wagraph))
				val -= (*args->waQfunc) (nnode, nnode, NULL, args->waQfunc_arg)
				    * (*args->waQfunc) (nnode, node, NULL, args->waQfunc_arg);

			for (j = 0; j < args->waQgraph->nnbs[node]; j++) {
				k = args->waQgraph->nbs[node][j];
				if (GMRFLib_graph_is_nb(k, node, args->wagraph) && GMRFLib_graph_is_nb(k, nnode, args->wagraph))
					val += (*args->waQfunc) (k, node, NULL, args->waQfunc_arg)
					    * (*args->waQfunc) (k, nnode, NULL, args->waQfunc_arg);
			}
			return val;
		}
	}

}

int GMRFLib_init_nwa_problem(GMRFLib_nwa_problem_tp ** nwa_problem,
			     int n_wa, GMRFLib_graph_tp ** wagraph, GMRFLib_Qfunc_tp ** wafunc, void **wafunc_arg,
			     int n_g, GMRFLib_graph_tp ** graph, GMRFLib_Qfunc_tp ** Qfunc, void **Qfunc_arg)
{
	/*
	 * init this problem:
	 * 
	 * \sum_s \sum_i (w_ii x_i - \sum_{j~i} w_ij x_j)^2 + \sum_k x^T Q x = x^T Q x
	 * 
	 * where all w_.. and ~ depends on 's = 0...n_wa-1'.
	 * 
	 * this routine setup a sequence of wa_problems and ordinary graphs/Qfuncs and return the graph and the Qfunction for
	 * the sum.
	 * 
	 */

	int k;
	GMRFLib_nwa_problem_tp **nwa = NULL;
	GMRFLib_graph_tp **nwa_graphs;

	GMRFLib_ASSERT(n_wa >= 0 && n_g >= 0, GMRFLib_EINVARG);

	if (n_wa + n_g <= 0) {
		*nwa_problem = NULL;
		return GMRFLib_SUCCESS;
	}

	nwa = Calloc(n_wa + n_g + 1, GMRFLib_nwa_problem_tp *);	/* yes! use a NULL-terminated list */
	nwa_graphs = Calloc(n_wa + n_g, GMRFLib_graph_tp *);

	for (k = 0; k < n_wa; k++) {			       /* do the wa_problems */
		GMRFLib_wa_problem_tp *wa = NULL;

		GMRFLib_init_wa_problem(&wa, wagraph[k], wafunc[k], (wafunc_arg ? wafunc_arg[k] : NULL));

		nwa[k] = Calloc(1, GMRFLib_nwa_problem_tp);
		nwa[k]->graph = wa->graph;
		nwa[k]->Qfunc = wa->Qfunc;
		nwa[k]->Qfunc_arg = wa->Qfunc_arg;
		nwa[k]->wa_problem = 1;			       /* must set this flag for the FREE to be ok */
		nwa_graphs[k] = wa->graph;		       /* use it to build the union-graph */
		Free(wa);
	}
	for (k = 0; k < n_g; k++) {			       /* then add ordinary Qfunc/graph pairs */
		nwa[k + n_wa] = Calloc(1, GMRFLib_nwa_problem_tp);
		nwa[k + n_wa]->graph = graph[k];	       /* do not copy graph */
		nwa[k + n_wa]->Qfunc = Qfunc[k];
		nwa[k + n_wa]->Qfunc_arg = (Qfunc_arg ? Qfunc_arg[k] : NULL);	/* no copy here either */
		nwa[k + n_wa]->wa_problem = 0;		       /* do not FREE as a wa_problem */
		nwa_graphs[k + n_wa] = graph[k];	       /* this is just to compute the union graph */
	}

	*nwa_problem = Calloc(1, GMRFLib_nwa_problem_tp);

	GMRFLib_graph_union(&((*nwa_problem)->graph), nwa_graphs, n_wa + n_g);
	Free(nwa_graphs);

	(*nwa_problem)->Qfunc = GMRFLib_nwaQfunc;
	(*nwa_problem)->Qfunc_arg = (void *) nwa;
	(*nwa_problem)->wa_problem = 0;

	/*
	 * prune graph is not needed, as each wa_graph[k] is pruned. 
	 */

	return GMRFLib_SUCCESS;
}

int GMRFLib_free_nwa_problem(GMRFLib_nwa_problem_tp * nwa_problem)
{
	/*
	 * free a nwa_problem 
	 */

	int k;
	GMRFLib_nwa_problem_tp **nwa = NULL;

	if (!nwa_problem)
		return GMRFLib_SUCCESS;

	nwa = (GMRFLib_nwa_problem_tp **) (nwa_problem->Qfunc_arg);

	for (k = 0; nwa[k]; k++) {
		/*
		 * check if this is a wa_problem or not. matters for how this is free'ed 
		 */
		if (nwa[k]->wa_problem) {
			GMRFLib_wa_problem_tp *wa = NULL;

			wa = Calloc(1, GMRFLib_wa_problem_tp);
			wa->graph = nwa[k]->graph;
			wa->Qfunc = nwa[k]->Qfunc;
			wa->Qfunc_arg = nwa[k]->Qfunc_arg;

			GMRFLib_free_wa_problem(wa);
		} else {
			/*
			 * nothing special to do here 
			 */
		}
		Free(nwa[k]);
	}
	GMRFLib_graph_free(nwa_problem->graph);
	Free(nwa_problem->Qfunc_arg);
	Free(nwa_problem);

	return GMRFLib_SUCCESS;
}

double GMRFLib_nwaQfunc(int node, int nnode, double *UNUSED(values), void *arg)
{
	if (node >= 0 && nnode < 0) {
		return NAN;
	}

	int k, equal;
	double val = 0.0;
	GMRFLib_nwa_problem_tp **nwa = NULL;

	nwa = (GMRFLib_nwa_problem_tp **) arg;

	for (k = 0, equal = (node == nnode); nwa[k]; k++)
		if (equal || GMRFLib_graph_is_nb(node, nnode, nwa[k]->graph))
			val += (*(nwa[k]->Qfunc)) (node, nnode, NULL, nwa[k]->Qfunc_arg);

	return val;
}
