
/* wa.c
 * 
 * Copyright (C) 2001-2006 Havard Rue
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
  \file wa.c
  \brief Setup for the generation of graphs for a weighted average model.

  <em>The wa-graph</em>: The neighbourhood-structure of the graph of on which a GMRF <em>\b x</em>
  is defined, is specified indirectly by specifying the weights, \f$ w_{ij} \f$, and the graph
  corresponding to these weights, and then \em GMRFLib will compute the graph of <em>\b x</em>. The
  graph is specified in three steps:

  - Create a function \em wafunc, say, returning the weights \f$ w_{ij};\; i=1,\ldots,n;\;
  j=1,\ldots,n \f$ for each pair of nodes \em i and \em j. The function should be of the same format
  as the function computing the elements of <em>\b Q</em> of a general graph, that is, of type \c
  GMRFLib_Qfunc_tp().
  
  - Specify the graph of the weights \f$ w_{ij} \f$, creating a \c GMRFLib_graph_tp -object.

  - Call \c GMRFLib_init_wa_problem(), creating the graph and the <em>\b Q</em> -function of the
  GMRF <em>\b x</em> corresponding to the user-defined weight function. This function returns an
  object of type \c GMRFLib_wa_problem_tp, holding the graph and the <em>\b Q</em> -function needed
  for sampling and evaluation.\n The members of the data structure can then be used just as the
  graph and the <em>\b Q</em> -function and it's arguments for a general graph.\n\n \c
  GMRFLib_init_nwa_problem() returns an object of type \c GMRFLib_nwa_problem_tp, holding the graph
  and the <em>\b Q</em> -function needed for sampling and evaluation. \n The members of the data
  structure can then be used just as the graph and the <em>\b Q</em> -function and it's arguments
  for a general graph.

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

#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

/* Pre-hg-Id: $Id: wa.c,v 1.23 2008/08/26 07:07:13 hrue Exp $ */

/*!
  \brief Deallocates a \c GMRFLib_wa_problem_tp -object and it's members, allocated by 
  \c GMRFLib_init_wa_problem.

  \param[in,out] wa_problem A \c GMRFLib_wa_problem_tp -object, allocated by \c
  GMRFLib_init_wa_problem(). At output, \a wa_problem and it's members have been deallocated.

  \note DO NOT USE \c GMRFLib_free_wa_problem() to deallocate a wa-problem created by \c
  GMRFLib_init_nwa_problem(), as the internal data structures are different.

  \sa GMRFLib_init_wa_problem
 */
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

		GMRFLib_free_graph(wa_problem->graph);
		GMRFLib_free_graph(arg->waQgraph);
		GMRFLib_free_graph(arg->wagraph);

		Free(wa_problem->Qfunc_arg);
		Free(wa_problem);
	}
	return GMRFLib_SUCCESS;
}

/*!

  \brief Computes the graph corresponding to the precision matrix <em>\b Q</em> based on the weight
  function.

  This function converts a graph definition on the weights \f$ w_{ij} \f$ of a weighted average GMRF
  model, as defined by <b>(GMRF-8)</b> in \ref specification, to the graph defining the precision
  matrix <em>\b Q</em> of the GMRF <em>\b x</em>.

  \param[out] wa_problem At output, \a wa_problem is allocated as a pointer to a \c
  GMRFLib_wa_problem_tp, holding the graph of the GMRF <em>\b x</em>, the function defining the
  <em>\b Q</em> -matrix and a \em void -pointer holding the address of the variable or data
  structure holding it's arguments.

  \param[in] wagraph The graph of the weights \f$ {w_{ij}} \f$ in the density function of <em>\b
  x</em>, of the form <b>(GMRF-8)</b>.

  \param[in] wafunc A pointer to a function returning the weights \f$ w_{ij} \f$. The function is to
  be of the same format as the function defining the elements of the <em>\b Q</em> -matrix of a
  general GMRF, such that the argument list and return value should be the same as for the template
  \c GMRFLib_Qfunc_tp().  The arguments should be the indices \em i and \em j and a \em void
  -pointer referring to additional arguments to the function. If \f$ w_{ij}=0 \f$ initially, it is
  required to be kept unchanged.

  \param[in] wafunc_arg A \em void -pointer holding the address of a variable or data structure
  defining additional arguments to the function \a wafunc.

  \remarks Based on the weight function, defining the weights in <b>(GMRF-8)</b>, the routine
  computes the graph corresponding to the precision matrix <em>\b Q</em> of the GMRF \em x. Also,
  the function defining the elements of <em>\b Q</em> and the set of arguments (in addition to the
  indices \em i and \em j) are generated. These are returned as members of the \c
  GMRFLib_wa_problem_tp -object <em>(*problem)</em>.

  \note The weights \f$ w_{ij} \f$ that are initially defined to be 0, should be kept fixed. The
  function \c GMRFLib_prune_graph() is called on the resulting graph of the GMRF <em>\b x</em>,
  removing elements of the graph corresponding to <em>\b Q (i,j) = 0</em>.  Re-setting the zero
  weights, for example by letting the weights depend on parameters that might change while running
  the programs, will invalidate this graph reduction.\n This routine will in most cases compute
  <em>\b Q (i,j)</em> less efficiently using more memory than a tailored implementation, but may
  save you for a lot of work!!! \n There is a global variable \c GMRFLib_use_wa_table_lookup() which
  controls the internal behaviour: if it is \c #GMRFLib_TRUE (default value) then internal
  lookup-tables are build that (can really) speed up the computation, and if it is \c
  #GMRFLib_FALSE, then it does not use internal lookup-tables.  The storage requirement is \f$
  {\mathcal O}(n^{3/2}) \f$.

  \par Example:
  See \ref ex_wa
  
  \sa GMRFLib_free_wa_problem, GMRFLib_prune_graph.
 */
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

	GMRFLib_make_empty_graph(&graph);
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
			memcpy(graph->nbs[i], wagraph->nbs[i], graph->nnbs[i] * sizeof(int));
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

					if (!GMRFLib_is_neighb(k, kk, graph)) {
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
			memcpy(&hold[indx], graph->nbs[i], graph->nnbs[i] * sizeof(int));
			Free(graph->nbs[i]);
			graph->nbs[i] = &hold[indx];
		} else
			Free(graph->nbs[i]);

		indx += graph->nnbs[i];
	}

	Free(memsiz);
	GMRFLib_prepare_graph(graph);

	/*
	 * setup the new types 
	 */
	wa_arg = Calloc(1, GMRFLib_waQfunc_arg_tp);
	wa_arg->waQgraph = graph;
	wa_arg->waQfunc = wafunc;
	wa_arg->waQfunc_arg = wafunc_arg;
	GMRFLib_copy_graph(&(wa_arg->wagraph), wagraph);       /* yes, make a copy! */

	if (GMRFLib_use_wa_table_lookup) {
		/*
		 * build hash table for idx = (i,j) -> neigb_info[idx] 
		 */
		double idx = 0.0;
		int nelm;

		GMRFLib_nQelm(&nelm, graph);
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
	GMRFLib_prune_graph(&ngraph, (*wa_problem)->graph, (*wa_problem)->Qfunc, (*wa_problem)->Qfunc_arg);
	(*wa_problem)->graph = ngraph;

	return GMRFLib_SUCCESS;
}

/*
  NOT DOCUMENTED...
  
  \brief Returns the value Q(node,nnode), that is, the element of the precision matrix corresponding to
  nodes node and nnode.

  Starting with the expression \f[ (*) \sum_i (w_{ii} x_i - \sum_{j~i} w_{ij} x_j)^2, \f] for given
  weights \f$ w_{ij} \f$, this function returns \f$ Q_{ij} \f$, <em>i = node, j = nnode</em>, so
  that \f[ (*) = x^T Q x \f] for *graph*, and where the function <em>wa_func(i,j)</em> return the
  weight \f$ w_{ij} \f$. Note that <em>wa_func(i,j) != wa_func(j,i)</em> in general and the graph
  defining <em>\b Q</em> is \em graph and not \em wagraph. \em graph can be found from \em wagraph
  using \c GMRFLib_wagraph2graph().
  
  \param[in] node,nnode The nodes of the graph for which to compute the precision.  \param[in] arg
  Must be of type \c GMRFLib_waQfunc_arg_tp, containing \em graph, \em wagraph and \em wa_func.

  \note There is potential great speedup in this routine by storing all neighbour-relations to
  reduce all the lookups. This is now implemented and controlled by the global varable
  GMRFLib_use_wa_table_lookup.
*/
double GMRFLib_waQfunc(int node, int nnode, void *arg)
{
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
					if (GMRFLib_is_neighb(k, node, args->wagraph))
						count++;
				}
				args->neigb_info[idx] = Calloc(1, GMRFLib_node_list_tp);
				if (count)
					args->neigb_info[idx]->node_list = Calloc(count, int);

				args->neigb_info[idx]->n_nodes = count;
				for (j = 0, jj = 0; j < args->waQgraph->nnbs[node]; j++) {
					k = args->waQgraph->nbs[node][j];
					if (GMRFLib_is_neighb(k, node, args->wagraph))
						args->neigb_info[idx]->node_list[jj++] = k;
				}
			}

			tmp = (*args->waQfunc) (node, node, args->waQfunc_arg);
			val = SQR(tmp);
			for (j = 0; j < args->neigb_info[idx]->n_nodes; j++) {
				k = args->neigb_info[idx]->node_list[j];
				tmp = (*args->waQfunc) (k, node, args->waQfunc_arg);
				val += SQR(tmp);
			}
			return val;
		} else {
			if (!args->neigb_info[idx]) {
				/*
				 * then build up information 
				 */

				args->neigb_info[idx] = Calloc(1, GMRFLib_node_list_tp);

				args->neigb_info[idx]->node_nnode = GMRFLib_is_neighb(node, nnode, args->wagraph);
				args->neigb_info[idx]->nnode_node = GMRFLib_is_neighb(nnode, node, args->wagraph);

				for (j = 0, count = 0; j < args->waQgraph->nnbs[node]; j++) {
					k = args->waQgraph->nbs[node][j];
					if (GMRFLib_is_neighb(k, node, args->wagraph)
					    && GMRFLib_is_neighb(k, nnode, args->wagraph))
						count++;
				}
				if (count)
					args->neigb_info[idx]->node_list = Calloc(count, int);

				args->neigb_info[idx]->n_nodes = count;
				for (j = 0, jj = 0; j < args->waQgraph->nnbs[node]; j++) {
					k = args->waQgraph->nbs[node][j];
					if (GMRFLib_is_neighb(k, node, args->wagraph)
					    && GMRFLib_is_neighb(k, nnode, args->wagraph))
						args->neigb_info[idx]->node_list[jj++] = k;
				}
			}

			val = 0.0;
			if (args->neigb_info[idx]->node_nnode)
				val -= (*args->waQfunc) (node, node, args->waQfunc_arg)
				    * (*args->waQfunc) (node, nnode, args->waQfunc_arg);
			if (args->neigb_info[idx]->nnode_node)
				val -= (*args->waQfunc) (nnode, nnode, args->waQfunc_arg)
				    * (*args->waQfunc) (nnode, node, args->waQfunc_arg);

			for (j = 0; j < args->neigb_info[idx]->n_nodes; j++) {
				k = args->neigb_info[idx]->node_list[j];
				val += (*args->waQfunc) (k, node, args->waQfunc_arg)
				    * (*args->waQfunc) (k, nnode, args->waQfunc_arg);
			}
			return val;
		}
	} else {					       /* plain version */

		if (node == nnode) {
			tmp = (*args->waQfunc) (node, node, args->waQfunc_arg);
			val = SQR(tmp);

			for (j = 0; j < args->waQgraph->nnbs[node]; j++) {
				k = args->waQgraph->nbs[node][j];
				if (GMRFLib_is_neighb(node, k, args->wagraph)) {
					tmp = (*args->waQfunc) (k, node, args->waQfunc_arg);
					val += SQR(tmp);
				}
			}
			return val;
		} else {
			val = 0.0;
			if (GMRFLib_is_neighb(node, nnode, args->wagraph))
				val -= (*args->waQfunc) (node, node, args->waQfunc_arg)
				    * (*args->waQfunc) (node, nnode, args->waQfunc_arg);
			if (GMRFLib_is_neighb(nnode, node, args->wagraph))
				val -= (*args->waQfunc) (nnode, nnode, args->waQfunc_arg)
				    * (*args->waQfunc) (nnode, node, args->waQfunc_arg);

			for (j = 0; j < args->waQgraph->nnbs[node]; j++) {
				k = args->waQgraph->nbs[node][j];
				if (GMRFLib_is_neighb(k, node, args->wagraph) && GMRFLib_is_neighb(k, nnode, args->wagraph))
					val += (*args->waQfunc) (k, node, args->waQfunc_arg)
					    * (*args->waQfunc) (k, nnode, args->waQfunc_arg);
			}
			return val;
		}
	}

}

/*!
  \brief This function is a generalization of \c GMRFLib_init_wa_problem().

  This function is a generalization of GMRFLib_init_wa_problem() and compuate the equivalent graph
  and Q-function for
  \f[ \pi(\mbox{\boldmath $x$}) \propto \mbox{exp}\left(-\frac{1}{2}\sum_{s=0}^{n_{wa}-1}
  \sum_{i}\left(w_{ii,s} x_i - \sum_{j\sim_s i} w_{ij,s} x_{j}\right)^2\right) \times
  \mbox{exp}\left(-\frac{1}{2}\sum_{s=0}^{n_{g}-1} \mbox{\boldmath $x$}^T\mbox{\boldmath
  $Q$}_s\mbox{\boldmath $x$}\right) \hspace{1cm} (*) \f] allowing for a sum of weighted avarages and
  sum of "ordinary" terms \f$ \mbox{\boldmath $x$}^T\mbox{\boldmath $Q$}\mbox{\boldmath $x$} \f$.

  \param[out] nwa_problem At output, \a nwa_problem is allocated as a pointer to a \c
  GMRFLib_nwa_problem_tp, holding the graph of the GMRF <em>\b x</em>, the function defining the
  <em>\b Q</em> -matrix and a \em void -pointer holding the address of the variable or data
  structure holding it's arguments.

  \param[in] n_wa Number of weighted terms (\f$ n_{wa} \f$ in <b>(*)</b>).

  \param[in] wagraph An array with length \a n_wa of pointers to graphs, where <em>wa_graph[s]</em>
  corresponds to term \em s in <b>(*)</b>.

  \param[in] wafunc An array with length \a n_wa of pointers to a function returning the weights \f$
  w_{ij,s} \f$. The function is to be of the same format as the function defining the elements of
  the <em>\b Q</em> -matrix of a general GMRF, such that the argument list and return value should
  be the same as for the template \c GMRFLib_Qfunc_tp(). The arguments should be the indices \em i
  and \em j and a \em void -pointer referring to additional arguments to the function. If \f$
  w_{ij,s}=0 \f$ initially, it is required to be kept unchanged.

  \param[in] wafunc_arg An array with length \a n_wa of \em void-pointers holding the address
  of a variable or data structure defining additional arguments to the function \a wafunc.

  \param[in] n_g Number of ordinary terms (\f$ n_{wa} \f$ in <b>(*)</b>).

  \param[in] graph An array with length \a n_g of pointers to graphs, where <em>graph[s]</em>
  corresponds to term \em s in <b>(*)</b>).

  \param[in] Qfunc An array with length \a n_g of pointers to a function returning \f$ Q_{ij,s} \f$.

  \param[in] Qfunc_arg An array with length \a n_g of \em void-pointers holding the address of
  a variable or data structure defining additional arguments to the function \a Qfunc.
  
  \remarks Based on the input arguments the routine computes the graph corresponding to the
  precision matrix <em>\b Q</em> of the GMRF <em>\b x</em>. Also, the function defining the elements
  of <em>\b Q</em> and the set of arguments (in addition to the indices \em i and \em j ) are
  generated. These are returned as members of the \c GMRFLib_nwa_problem_tp -object
  <em>(*problem)</em>.

  \note The weights \f$ w_{ij,s} \f$ that are initially defined to be 0, should be kept fixed. The
  function \c GMRFLib_prune_graph() is called on the resulting graph of the GMRF <em>\b x</em>,
  removing elements of the graph corresponding to <em>\b Q (i,j)=0</em>. Re-setting the zero
  weights, for example by letting the weights depend on parameters that might change while running
  the programs, will invalidate this graph reduction. \n This routine will in most cases compute
  <em>\b Q (i,j)</em> less efficiently than a tailored implementation, but may save you for a lot of
  work!

  \sa GMRFLib_free_nwa_problem, GMRFLib_prune_graph.
 */
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

	GMRFLib_union_graph(&((*nwa_problem)->graph), nwa_graphs, n_wa + n_g);
	Free(nwa_graphs);

	(*nwa_problem)->Qfunc = GMRFLib_nwaQfunc;
	(*nwa_problem)->Qfunc_arg = (void *) nwa;
	(*nwa_problem)->wa_problem = 0;

	/*
	 * prune graph is not needed, as each wa_graph[k] is pruned. 
	 */

	return GMRFLib_SUCCESS;
}

/*!  \brief Deallocates a \c GMRFLib_nwa_problem_tp -object and it's members, allocated by \c
  GMRFLib_init_nwa_problem().

  \param[in,out] nwa_problem A \c GMRFLib_nwa_problem_tp -object, allocated by \c
  GMRFLib_init_nwa_problem(). At output, \a nwa_problem and it's members have been deallocated.

  \note DO NOT USE \c GMRFLib_free_nwa_problem() to deallocate a wa-problem created by \c
  GMRFLib_init_wa_problem(), as the internal data structures are different.

  \sa GMRFLib_init_nwa_problem
 */
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
	GMRFLib_free_graph(nwa_problem->graph);
	Free(nwa_problem->Qfunc_arg);
	Free(nwa_problem);

	return GMRFLib_SUCCESS;
}
double GMRFLib_nwaQfunc(int node, int nnode, void *arg)
{
	int k, equal;
	double val = 0.0;
	GMRFLib_nwa_problem_tp **nwa = NULL;

	nwa = (GMRFLib_nwa_problem_tp **) arg;

	for (k = 0, equal = (node == nnode); nwa[k]; k++)
		if (equal || GMRFLib_is_neighb(node, nnode, nwa[k]->graph))
			val += (*(nwa[k]->Qfunc)) (node, nnode, nwa[k]->Qfunc_arg);

	return val;
}

/*
  Example for manual
 */

/*! \page ex_wa A worked out example smoothing a time-series data, using the wa_problem feature
  
  In this example, we will illustrate the \em wa_problem feature, which is quite useful in situations like the following.
  Consider the time-series data <em>\b y</em> in Fig..., which we assume to be observed according to the model \f[ y_i = x_i +
  \epsilon_i,\f] where <em>\b x</em> is a hidden smooth curve and \f$ \bf\epsilon \f$ is \em iid Gaussian noise with zero mean
  and variance \f$ 0.25^2 \f$ .

  As a prior, or a smoothing-prior, for <em>\b x</em>, we will consider two models: a first order model \f[\pi(\bf{x}) \propto
  \exp\left(-\frac{h}{2} \sum_{i=2}^n (x_i-x_{i-1})^2\right) \hspace{2cm} (wa-1) \f] and a second order model \f[ \pi(\bf{x})
  \propto \exp\left(-\frac{h}{2} \sum_{i=3}^n (x_i-2x_{i-1}+x_{i-2})^2\right). \hspace{2cm} (wa-2) \f] The quantities of
  interest are the posterior mean for <em>\b x|y</em> and samples from the posterior \f[ \pi(\bf{x}|\bf{y}) \propto
  \exp\left(-\frac{h}{2} \sum_{i=2}^n (x_i-x_{i-1})^2+ \sum_{i=1}^n (y_ix_i/\sigma^{2}+x_i^2/\sigma^{2}) \right). \hspace{2cm}
  (wa-3) \f] A similar expression can be obtained for the second order model.

  To implement Equation <b>(wa-3)</b>, we observe that the prior is of the \em wa_problem form <b>(GMRF-8)</b> in \ref
  specification, for which we use the function \c GMRFLib_init_wa_problem() to convert the prior into the standard <em>\b
  x'Qx</em>-form. The data will define the <em>\b b</em> and <em>\b c</em> vector in \c GMRFLib_init_problem().

  \par Program code:

  \verbinclude example-doxygen-wa.txt

  The program provides the output in the figures below, for order=1 and <em>h=50</em>, and order=2
  and <em>h=100</em>.

  \htmlonly
  <table>
  <td><img src="figs/order1.gif" width="400"> </td>
  <td><img src="figs/order2.gif" width="400"> </td>
  </table> 
  The observed and smoothed time-series, order=1 and <em>h=50</em>
  to the left, and order=2 and <em>h=100</em> to the right. The figures
  display the observed data, the posterior mean and 10 iid
  samples from the posterior.
  \endhtmlonly
*/
