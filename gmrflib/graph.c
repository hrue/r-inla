
/* graph.c
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
  \file graph.c
  \brief Functions for creating and handling general and regular graphs
  
  There are two ways of specifying a general graph.

  - <em> Reading from a file.</em>
  The simplest way to specify a general
  graph is to read the graph specifications from a file, using the
  function \c GMRFLib_read_graph_ascii(). This function allocates
  and initializes the members based on the information given on the
  input file, and customizes the graph for use in other library
  functions.
  - <em> Explicitly creating a \c GMRFLib_graph_tp - object.</em> 
  The user can create a general graph by allocating and
  initializing a variable of type \c GMRFLib_graph_tp, and then 
  (IMPORTANT!) calling the function \c GMRFLib_prepare_graph() to customize 
  the graph and check that the graph is consistently defined.  The user 
  should specify the the members \em n, \em nnbs and \em nbs in 
  \c GMRFLib_graph_tp. The member \em mothergraph_idx is initialized, 
  and needed, only if the graph is generated as a subgraph of another graph, 
  using \c GMRFLib_compute_subgraph().

  In addition to the functions for specifying general graphs, the library
  provides two separate functions for the specification of
  two-dimensional lattice graphs and one-dimensional linear graphs, and
  a separate setup for the generation of graphs for a weighted average
  model (see wa.c).

  - <em> A lattice graph:</em> A two-dimensional graph on a lattice is
  most easily specified using the function \c GMRFLib_make_lattice_graph(), 
  specifying the grid sizes \f$ n_{row} \f$  and
  \f$ n_{col} \f$ and the parameters \f$ m_{row} \f$ and \f$ m_{col} \f$, 
  defining the neighbourhood  \f$ (2 m_{row} + 1)\times (2 m_{col} + 1) \f$.
  - <em> A linear graph:</em> To specify a one-dimensional linear graph,
  i.e. an autoregressive AR(\em p)-model, use the function 
  \c GMRFLib_make_linear_graph(), specifying the number of nodes
  and the size of the neighbourhood.
*/

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

/* Pre-hg-Id: $Id: graph.c,v 1.102 2010/02/15 08:26:37 hrue Exp $ */

/*!
  \brief Creates an empty graph
  \param[in,out] graph A pointer to a \c GMRFLib_graph_tp pointer. 
  At output, \em graph points to an empty graph.
  \par Example:
  <tt> main(int argc, char* argv[])\n
  {\n
  ::::::::\n
  GMRFLib_graph_tp *graph;\n
  GMRFLib_make_empty_graph(&graph);\n
  ::::::::\n
  }
  </tt>
*/
int GMRFLib_make_empty_graph(GMRFLib_graph_tp ** graph)
{
	/*
	 * this function creates an empty graph 
	 */

	*graph = Calloc(1, GMRFLib_graph_tp);

	/*
	 * user variables 
	 */
	(*graph)->n = 0;
	(*graph)->nbs = NULL;
	(*graph)->nnbs = NULL;

	/*
	 * private variables 
	 */
	(*graph)->mothergraph_idx = NULL;

	return GMRFLib_SUCCESS;
}

/*!
  \brief Reads a graph from a file

  \param[in,out] graph At output, <em>(*graph)</em> has been initialized 
  by using the information on the file \em filename. 

  \param[in] filename The name of the file, formatted as described 
  below, containing the specification of the graph.

  \remarks The file is assumed to be of the format
  \verbinclude doxygen_file_format.txt
  Here, \em n, \em nbs and \em nnbs refer to the members of 
  \c GMRFLib_graph_tp, and <em>nn[i]</em> is the node number
  of the node with index \em i, running from \em 0 to <em>n-1</em>.  
  In words, the first row contains the number of nodes, \em n, and the
  successive rows, rows <em>i = 2,...,n+1</em>, list the node number,
  <em>nn[i]</em>, the number of neighbours, <em>nnbs[i]</em>, and
  the node numbers, <em>nbs[i][j]</em>, for the neighbours of each
  node \em i of the graph. The number of neighbours <em>nnbs[i]</em>
  might be zero.\n\n
  The function calls \c GMRFLib_prepare_graph() to
  customize the graph for further computations.

  \par Example:
  See \ref ex_graph

  \sa GMRFLib_prepare_graph, GMRFLib_print_graph
*/
int GMRFLib_read_graph(GMRFLib_graph_tp ** graph, const char *filename)
{
	GMRFLib_read_graph_binary(graph, filename);
	if (*graph != NULL) {
		return GMRFLib_SUCCESS;
	} else {
		return GMRFLib_read_graph_ascii(graph, filename);
	}
}
int GMRFLib_read_graph_ascii(GMRFLib_graph_tp ** graph, const char *filename)
{
#define TO_INT(_ix, _x) \
	if (1) {							\
		_ix = (int) (_x);					\
		if (!ISEQUAL((double) (_ix), (_x))) {			\
			char *msg;					\
			GMRFLib_sprintf(&msg, "Error reading graph. This is not an integer [%g]", _x); \
			GMRFLib_ERROR_MSG(GMRFLib_EGRAPH, msg);		\
		}							\
	}

	/*
	 * read a graph in the following format
	 * 
	 * N node[0] nnbs[0] nbs[node[0]][0] nbs[node[0]][1] ... nbs[node[0]][nnbs[0]-1] node[1] nnbs[1] nbs[node[1]][0]
	 * nbs[node[1]][1] ... nbs[node[1]][nnbs[1]-1] : node[N-1] nnbs[N-1] nbs[node[N-1]][0] nbs[node[N-1]][1] ...
	 * nbs[node[N-1]][nnbs[N-1]-1] 
	 */

	/*
	 * dec 09: add the posibility to give a 1-based graph. if min(node)=1 and max(node)=n, then this defines a 1-based graph.
	 * otherwise, its a 0-based graph.
	 */

	int *storage = NULL, n_neig_tot = 0, storage_indx, itmp;
	int i, j, tnode, min_node, max_node;
	double dnode, tmp;
	GMRFLib_io_tp *io = NULL;

	GMRFLib_EWRAP0(GMRFLib_io_open(&io, filename, "r"));
	GMRFLib_EWRAP0(GMRFLib_make_empty_graph(graph));

	GMRFLib_EWRAP0(GMRFLib_io_read_next(io, &tmp, "%lf"));
	TO_INT((*graph)->n, tmp);
	GMRFLib_ASSERT((*graph)->n >= 0, GMRFLib_EPARAMETER);

	(*graph)->nnbs = Calloc((*graph)->n + 1, int);
	(*graph)->nbs = Calloc((*graph)->n + 1, int *);

	min_node = INT_MAX;
	max_node = INT_MIN;

	for (i = 0; i < (*graph)->n; i++) {
		GMRFLib_EWRAP0(GMRFLib_io_read_next(io, &dnode, "%lf"));
		TO_INT(tnode, dnode);
		if (tnode < 0 || tnode > (*graph)->n) {
			GMRFLib_io_error(io, GMRFLib_IO_ERR_READLINE);
			return GMRFLib_EREADFILE;
		}
		min_node = IMIN(min_node, tnode);
		max_node = IMAX(max_node, tnode);

		GMRFLib_EWRAP0(GMRFLib_io_read_next(io, &tmp, "%lf"));
		TO_INT(itmp, tmp);
		if (itmp < 0 || itmp > (*graph)->n) {
			GMRFLib_io_error(io, GMRFLib_IO_ERR_READLINE);
			return GMRFLib_EREADFILE;
		}
		(*graph)->nnbs[tnode] = itmp;
		n_neig_tot += (*graph)->nnbs[tnode];

		if ((*graph)->nnbs[tnode]) {
			(*graph)->nbs[tnode] = Calloc((*graph)->nnbs[tnode], int);

			for (j = 0; j < (*graph)->nnbs[tnode]; j++) {
				GMRFLib_EWRAP0(GMRFLib_io_read_next(io, &tmp, "%lf"));
				TO_INT(itmp, tmp);
				if (itmp < 0 || itmp > (*graph)->n) {
					GMRFLib_io_error(io, GMRFLib_IO_ERR_READLINE);
					return GMRFLib_EREADFILE;
				}
				(*graph)->nbs[tnode][j] = itmp;
				min_node = IMIN(min_node, itmp);
				max_node = IMAX(max_node, itmp);
			}
		} else {
			(*graph)->nbs[tnode] = NULL;
		}
	}

	if (min_node == 1 && max_node == (*graph)->n) {
		/*
		 * This is a 1-based graph; convert it to a zero-based graph 
		 */
		int im;

		for (i = 1; i <= (*graph)->n; i++) {	       /* YES! */
			im = i - 1;
			(*graph)->nnbs[im] = (*graph)->nnbs[i];
			(*graph)->nbs[im] = (*graph)->nbs[i];

			if ((*graph)->nnbs[im]) {
				for (j = 0; j < (*graph)->nnbs[im]; j++) {
					(*graph)->nbs[im][j]--;
				}
			}
		}
	} else {
		if (!(min_node == 0 && max_node == (*graph)->n - 1)) {
			/*
			 * then there is something wrong; write a message and say so 
			 */
			fprintf(stderr, "\n\n\n");
			fprintf(stderr, " *** Error in graph: n=%1d but min_node=%1d and max_node=%1d\n", (*graph)->n, min_node, max_node);
			fprintf(stderr, " *** The position in reported in the file is not correct...\n");
			GMRFLib_io_error(io, GMRFLib_IO_ERR_READLINE);
			return GMRFLib_EREADFILE;
		}
	}

	GMRFLib_EWRAP0(GMRFLib_io_close(io));


	/*
	 * map the graph to a more computational convenient memory layout! use just one long vector to store all the neighbors. 
	 */
	if (n_neig_tot) {
		storage = Calloc(n_neig_tot, int);

		storage_indx = 0;
		for (i = 0; i < (*graph)->n; i++) {
			if ((*graph)->nnbs[i]) {
				memcpy(&storage[storage_indx], (*graph)->nbs[i], (size_t) (sizeof(int) * (*graph)->nnbs[i]));
				Free((*graph)->nbs[i]);
				(*graph)->nbs[i] = &storage[storage_indx];
				storage_indx += (*graph)->nnbs[i];
			} else {
				(*graph)->nbs[i] = NULL;
			}
		}
	}
	if (GMRFLib_verify_graph_read_from_disc) {
		GMRFLib_EWRAP0(GMRFLib_validate_graph(stderr, *graph));
	}

	GMRFLib_EWRAP0(GMRFLib_prepare_graph(*graph));	       /* prepare the graph for computations */
#undef TO_INT
	return GMRFLib_SUCCESS;
}

/*!
  \brief Prints the specification of a graph to standard
    output or a file
  \param[out] fp The \em FILE* on which to print the graph.
  \param[in] graph The graph to be printed.
  \sa GMRFLib_read_graph
 */
int GMRFLib_print_graph(FILE * fp, GMRFLib_graph_tp * graph)
{
	int i, j;
	FILE *fpp = NULL;

	fpp = (fp ? fp : stdout);

	fprintf(fpp, "graph has %1d nodes\n", graph->n);
	for (i = 0; i < graph->n; i++) {
		fprintf(fpp, "node %1d has %1d neighbours:", i, graph->nnbs[i]);
		for (j = 0; j < graph->nnbs[i]; j++) {
			fprintf(fpp, " %1d", graph->nbs[i][j]);
		}
		fprintf(fpp, "\n");
	}
	if (graph->mothergraph_idx) {
		for (i = 0; i < graph->n; i++) {
			fprintf(fpp, "node %1d corresponds to mother-node %1d\n", i, graph->mothergraph_idx[i]);
		}
	}

	return GMRFLib_SUCCESS;
}

/*!
  \brief Write a graph to file in ascii format
  \param[in] filename The name of the file to store the graph in
    ascii format.
  \param[in] graph The graph to be written.
  \sa GMRFLib_read_graph.
 */
int GMRFLib_write_graph(const char *filename, GMRFLib_graph_tp * graph)
{
	/*
	 * write graph to file filename in the format so it can be read by 'read_graph' 
	 */

	FILE *fp = NULL;

	if (!filename) {
		return GMRFLib_SUCCESS;
	}
	if (!graph) {
		return GMRFLib_SUCCESS;
	}

	fp = fopen(filename, "w");
	if (!fp) {
		GMRFLib_ERROR(GMRFLib_EOPENFILE);
	}

	GMRFLib_write_graph_2(fp, graph);
	fclose(fp);

	return GMRFLib_SUCCESS;
}
int GMRFLib_write_graph_2(FILE * fp, GMRFLib_graph_tp * graph)
{
	int i, j;

	fprintf(fp, "%1d\n", graph->n);
	for (i = 0; i < graph->n; i++) {
		fprintf(fp, "%1d %1d", i, graph->nnbs[i]);
		for (j = 0; j < graph->nnbs[i]; j++) {
			fprintf(fp, " %1d", graph->nbs[i][j]);
		}
		fprintf(fp, "\n");
	}
	return GMRFLib_SUCCESS;
}

/*!
  \brief Write a graph in binary format
  \param[in] filename The name of the file
  \param[in] graph The graph to be written
  \sa GMRFLib_read_graph_binary()
 */
int GMRFLib_write_graph_binary(const char *filename, GMRFLib_graph_tp * graph)
{
	/*
	 * use base 1: so the nodes are 1...n, not 0..n-1. this makes the connection to R and R-inla easier. However, the read_graph routines will autodetect if
	 * the nodes are 0..n-1 or 1...n. 
	 */

	int i, tag = GMRFLib_BINARY_GRAPH_FILE_MAGIC, offset = 1, idx, iidx, j;
	GMRFLib_io_tp *io = NULL;

	if (!filename || !graph) {
		return GMRFLib_SUCCESS;
	}

	GMRFLib_EWRAP0(GMRFLib_io_open(&io, filename, "wb"));
	GMRFLib_io_write(io, (const void *) &tag, sizeof(int));	/* so we can detect that this is of correct format */
	GMRFLib_io_write(io, (const void *) &(graph->n), sizeof(int));
	for (i = 0; i < graph->n; i++) {
		idx = i + offset;
		GMRFLib_io_write(io, (const void *) &idx, sizeof(int));
		GMRFLib_io_write(io, (const void *) &(graph->nnbs[i]), sizeof(int));
		if (graph->nnbs[i]) {
			if (offset == 0) {
				GMRFLib_io_write(io, (const void *) (graph->nbs[i]), (unsigned int) (graph->nnbs[i] * sizeof(int)));
			} else {
				for (j = 0; j < graph->nnbs[i]; j++) {
					iidx = graph->nbs[i][j] + offset;
					GMRFLib_io_write(io, (const void *) &iidx, (unsigned int) sizeof(int));
				}
			}
		}
	}
	GMRFLib_io_close(io);

	return GMRFLib_SUCCESS;
}

/*!
  \brief Read a graph from file written by \c GMRFLib_write_graph_binary()
  \param[in] filename The name of the file
  \param[in] graph The graph 
  \sa GMRFLib_write_graph_binary()
 */
int GMRFLib_read_graph_binary(GMRFLib_graph_tp ** graph, const char *filename)
{
	int i, j, ii, idx, tag;
	GMRFLib_io_tp *io = NULL;
	GMRFLib_graph_tp *g = NULL;

	GMRFLib_EWRAP0(GMRFLib_io_open(&io, filename, "rb"));
	GMRFLib_EWRAP0(GMRFLib_io_read(io, (void *) &tag, sizeof(int)));
	if (tag != GMRFLib_BINARY_GRAPH_FILE_MAGIC) {
		/*
		 * this is not a binary graph file 
		 */
		GMRFLib_io_close(io);
		*graph = NULL;

		return !GMRFLib_SUCCESS;
	}

	GMRFLib_make_empty_graph(&g);
	GMRFLib_EWRAP0(GMRFLib_io_read(io, (void *) &(g->n), sizeof(int)));
	GMRFLib_ASSERT(g->n >= 0, GMRFLib_EPARAMETER);

	g->nnbs = Calloc(g->n + 1, int);		       /* yes. */
	g->nbs = Calloc(g->n + 1, int *);		       /* yes. */

	for (i = 0; i < g->n; i++) {
		GMRFLib_io_read(io, (void *) &ii, sizeof(int));
		GMRFLib_io_read(io, (void *) &(g->nnbs[ii]), sizeof(int));
		if (g->nnbs[ii]) {
			g->nbs[ii] = Calloc(g->nnbs[ii], int);
			GMRFLib_io_read(io, (void *) (g->nbs[ii]), (unsigned int) (g->nnbs[ii] * sizeof(int)));
		} else {
			g->nbs[ii] = Calloc(1, int);	       /* so we can check if this is read */
		}
	}
	GMRFLib_io_close(io);

	/*
	 * now we check if this graph is 1-based or 0-based. 
	 */
	int min_node = g->n + 1, max_node = -1;

	for (i = 0; i < g->n + 1; i++) {
		if (g->nbs[i]) {			       /* yes, we check the ptr if this is read or not */
			min_node = IMIN(min_node, i);
			max_node = IMAX(max_node, i);

			for (j = 0; j < g->nnbs[i]; j++) {
				idx = g->nbs[i][j];
				min_node = IMIN(min_node, idx);
				max_node = IMAX(max_node, idx);
			}
		}
	}

	if (min_node == 1 && max_node == g->n) {
		/*
		 * convert to 0-based graph 
		 */
		Free(g->nbs[0]);
		for (i = 0; i < g->n; i++) {
			g->nnbs[i] = g->nnbs[i + 1];
			g->nbs[i] = g->nbs[i + 1];
			g->nbs[i + 1] = NULL;
			for (j = 0; j < g->nnbs[i]; j++)
				g->nbs[i][j]--;
		}
	} else if (min_node == 0 && max_node == g->n - 1) {
		/*
		 * ok 
		 */
	} else {
		fprintf(stderr, "\n\nmin_node = %1d max_node = %1d. this should not happen.\n", min_node, max_node);
		GMRFLib_ERROR(GMRFLib_ESNH);
	}

	GMRFLib_EWRAP0(GMRFLib_prepare_graph(g));
	GMRFLib_EWRAP0(GMRFLib_copy_graph(graph, g));

	for (i = 0; i < g->n + 1; i++) {		       /* yes, its +1 */
		Free(g->nbs[i]);
	}
	Free(g->nbs);
	Free(g->nnbs);
	Free(g);

	return GMRFLib_SUCCESS;
}

/*!
  \brief Free a graph read or generated by GMRFLib

  Frees the memory held by a \c GMRFLib_graph_tp -variable.  NOTE: To ensure safe memory handling, this function should only be
  applied to graphs read or generated by the graph-generating routines of the library, and NOT to user-specified graphs.

  \param[in,out] graph  A graph generated by a graph-generating
  library routine. At output, the graph and it's array members are
  all deallocated.
  \sa GMRFLib_read_graph, GMRFLib_make_linear_graph, 
  GMRFLib_make_lattice_graph, GMRFLib_prepare_graph, GMRFLib_read_graph_binary, etc...
 */
int GMRFLib_free_graph(GMRFLib_graph_tp * graph)
{
	/*
	 * free a graph build with ``GMRFLib_read_graph'' 
	 */
	int i;

	if (!graph) {
		return GMRFLib_SUCCESS;
	}

	for (i = 0; i < graph->n; i++) {
		if (graph->nnbs[i]) {
			Free(graph->nbs[i]);
			break;				       /* new memory layout, only `free' the first!!! */
		}
	}
	Free(graph->nbs);
	Free(graph->nnbs);
	Free(graph->mothergraph_idx);
	Free(graph);

	return GMRFLib_SUCCESS;
}

/*!
  \brief Compute the number of non-zero elements in Q.
  \param[out] nelm Return the number of non-zero elements in Q in \a *nelm
  \param[in] graph The graph.
 */
int GMRFLib_nQelm(int *nelm, GMRFLib_graph_tp * graph)
{
	/*
	 * Return the number of non-zero elements in Q 
	 */
	int nn, i;

	for (i = 0, nn = graph->n; i < graph->n; i++) {
		nn += graph->nnbs[i];
	}

	*nelm = nn;
	return GMRFLib_SUCCESS;
}

/*!
  \brief Return bit-number BITNO, bitno = 0, 1, 2, ..., 7.
 */
int GMRFLib_getbit(GMRFLib_uchar c, unsigned int bitno)
{
	/*
	 * return bit-number BITNO, bitno = 0, 1, 2, ..., 7 
	 */
	unsigned int zero = 0;

	return (int) ((c >> bitno) & ~(~zero << 1));
}

/*!
  \brief Set bitno 0, 1, 2, ..., or 7, to TRUE.
 */
int GMRFLib_setbit(GMRFLib_uchar * c, unsigned int bitno)
{
	/*
	 * set bitno 0, 1, 2, ..., 7, to TRUE 
	 */
	unsigned int zero = 0;

	*c = *c | ((~(~zero << (bitno + 1)) & (~zero << bitno)));

	return GMRFLib_SUCCESS;
}

/*!
  \brief Print all the bits in \c c to \c fp
 */
int GMRFLib_printbits(FILE * fp, GMRFLib_uchar c)
{
	/*
	 * just print all the bits in C to FP 
	 */

	int j, nn = 8 * sizeof(GMRFLib_uchar);

	fprintf(fp, "int=[%u] : ", c);

	for (j = 0; j < nn; j++) {
		fprintf(fp, "%1d", (int) GMRFLib_getbit((GMRFLib_uchar) c, (unsigned int) (nn - j - 1)));
	}
	fprintf(fp, "\n");

	return GMRFLib_SUCCESS;
}

/*!
  \brief Checks whether two nodes \em i and \em j are neighbours.
  \param[in] node, nnode The nodes \em i and \em j to be checked.
  \param[in] graph The graph to be checked.
  \return Returns \c GMRFLib_TRUE if the nodes \em node and \em nnode are neighbours, 
          otherwise it returns \c GMRFLib_FALSE.
  \par Example:
  \verbinclude doxygen_is_neighbour.txt
 */
int GMRFLib_is_neighb(int node, int nnode, GMRFLib_graph_tp * graph)
{
	/*
	 * plain version.  return 1 if nnode is a neighbour of node, otherwise 0. assume that the nodes are sorted. (if node == 
	 * nnode, then they are not neighbours.)
	 * 
	 * Sat Nov 22 12:43:09 CET 2003 : I compared this version agains a hash-table variant, both using spmatrix and map_ii
	 * for each i. the last version was slower by a factor 2, the other even slower. it seems like this plain version is
	 * quite fast though... 
	 */

	int j, m, k;

	/*
	 * make this extention to ease its use 
	 */
	if (node < 0 || node >= graph->n || nnode < 0 || nnode >= graph->n)
		return GMRFLib_FALSE;

	m = graph->nnbs[node];

	if (!m || nnode < graph->nbs[node][0] || nnode > graph->nbs[node][m - 1]) {
		return GMRFLib_FALSE;
	}

	for (j = 0; j < m; j++) {
		k = graph->nbs[node][j];
		if (k > nnode) {
			return GMRFLib_FALSE;
		}
		if (k == nnode) {
			return GMRFLib_TRUE;
		}
	}
	return GMRFLib_FALSE;
}

/*!
  \brief Prepare the graph by sort the vertices in increasing orders.

  Post-processes a user-specified \c GMRFLib_graph_tp object, such that 
  the resulting graph is of the format required by the library functions 
  operating on graphs.

  \param[in,out] graph  A graph explicitly specified by the user
  by creating and initializing a \c GMRFLib_graph_tp -object.  
  At output, the graph is customized.

  \note If the user spesify its own graph, then this function MUST used to prepare internal
  structures and ensure that they are correct, otherwise, strange errors can occure.

  \sa GMRFLib_read_graph
 */
int GMRFLib_prepare_graph(GMRFLib_graph_tp * graph)
{
	/*
	 * prepare the graph by sort the vertices in increasing orders 
	 */
	GMRFLib_EWRAP0(GMRFLib_sort_nodes(graph));
	GMRFLib_EWRAP0(GMRFLib_make_nodes_unique(graph));

	return GMRFLib_SUCCESS;
}

int GMRFLib_make_nodes_unique(GMRFLib_graph_tp * graph)
{
	/*
	 * ensure the neigbours are unique. the neigbours must be sorted. 
	 */

	if (!graph) {
		return GMRFLib_SUCCESS;
	}

	int i;

	for (i = 0; i < graph->n; i++) {
		if (graph->nnbs[i]) {
			int k = 0, j;

			for (j = 1; j < graph->nnbs[i]; j++) {
				if (graph->nbs[i][k] != graph->nbs[i][j]) {
					graph->nbs[i][++k] = graph->nbs[i][j];
				}
			}
			graph->nnbs[i] = k + 1;
		}
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_sort_nodes(GMRFLib_graph_tp * graph)
{
	/*
	 * sort the vertices in increasing order 
	 */

	int i;

	if (!graph) {
		return GMRFLib_SUCCESS;
	}

	for (i = 0; i < graph->n; i++) {
		if (graph->nnbs[i]) {
			qsort(graph->nbs[i], (size_t) graph->nnbs[i], sizeof(int), GMRFLib_icmp);
		}
	}

	return GMRFLib_SUCCESS;
}
int GMRFLib_compute_bandwidth(int *bandwidth, GMRFLib_graph_tp * graph, int *remap)
{
	int bw = 0, i, j, node;

	for (i = 0; i < graph->n; i++) {
		node = remap[i];
		for (j = 0; j < graph->nnbs[i]; j++) {
			bw = IMAX(bw, node - remap[graph->nbs[i][j]]);
		}
	}
	*bandwidth = bw;

	return GMRFLib_SUCCESS;
}
int GMRFLib_find_idx(int *idx, int n, int *iarray, int value)
{
	int i;

	for (i = 0; i < n; i++)
		if (iarray[i] == value) {
			if (idx) {
				*idx = i;
			}
			return GMRFLib_SUCCESS;
		}
	return GMRFLib_EINDEX;
}

/*!
  \brief Validates the graph by checking that the members of a \c GMRFLib_graph_tp -object are defined consistently.

  \param[out] fp If \c !NULL, the function will write error messages stating the type of validation
  error on this file, in addition to the default standard output error message issued by
  GMRFLib_error_handler().  If \c NULL, only the default error message is issued.

  \param[in] graph The graph to be validated.

  \sa GMRFLib_error_handler
 */
int GMRFLib_validate_graph(FILE * fp, GMRFLib_graph_tp * graph)
{

	int i, j, jj, error = 0;

	if (graph->n == 0) {
		return GMRFLib_SUCCESS;
	}
	if (graph->n < 0) {
		if (fp) {
			fprintf(fp, "%s: error: graph->n = %1d < 0\n", __GMRFLib_FuncName, graph->n);
		}
		GMRFLib_ERROR(GMRFLib_EPARAMETER);
	}

	/*
	 * check nbs 
	 */
	for (i = 0; i < graph->n; i++)
		for (j = 0; j < graph->nnbs[i]; j++) {
			if (graph->nbs[i][j] < 0 || graph->nbs[i][j] >= graph->n) {
				error++;
				if (fp) {
					fprintf(fp, "\n\n%s: error: nbs[%1d][%1d]=[%1d] is out of range\n", __GMRFLib_FuncName, i, j,
						graph->nbs[i][j]);
				}
			}
			if (graph->nbs[i][j] == i) {
				error++;
				if (fp) {
					fprintf(fp, "\n\n%s: error: nbs[%1d][%1d] = node[%1d]. not allowed!\n", __GMRFLib_FuncName, i, j, i);
				}
			}
		}
	if (error)
		return GMRFLib_EGRAPH;

	/*
	 * check if i~j then j~i 
	 */

	for (i = 0; i < graph->n; i++)
		for (j = 0; j < graph->nnbs[i]; j++) {
			jj = graph->nbs[i][j];
			if (GMRFLib_find_idx(NULL, graph->nnbs[jj], graph->nbs[jj], i) != GMRFLib_SUCCESS) {
				if (fp) {
					fprintf(fp, "\n\n%s: error: node[%1d] has neighbor node[%1d], but not oposite\n", __GMRFLib_FuncName, i,
						jj);
				}
				error++;
			}
		}
	if (error) {
		return GMRFLib_EGRAPH;
	}

	/*
	 * ok 
	 */
	if (0 && fp) {
		fprintf(fp, "%s: graph is OK.\n", __GMRFLib_FuncName);
	}

	return GMRFLib_SUCCESS;
}
int GMRFLib_remap_graph(GMRFLib_graph_tp ** ngraph, GMRFLib_graph_tp * graph, int *remap)
{
	/*
	 * return the remapped graph based on 'graph'. the returned graph has the identity mapping. 
	 */

	int i, j, k, nnb, indx, *hold = NULL;

	if (!graph) {
		*ngraph = (GMRFLib_graph_tp *) NULL;
		return GMRFLib_SUCCESS;
	}

	GMRFLib_make_empty_graph(ngraph);
	(*ngraph)->n = graph->n;
	(*ngraph)->nnbs = Calloc((*ngraph)->n, int);
	(*ngraph)->nbs = Calloc((*ngraph)->n, int *);

	for (i = 0; i < graph->n; i++) {
		k = remap[i];
		(*ngraph)->nnbs[k] = graph->nnbs[i];
		if ((*ngraph)->nnbs[k]) {
			(*ngraph)->nbs[k] = Calloc((*ngraph)->nnbs[k], int);

			for (j = 0; j < (*ngraph)->nnbs[k]; j++) {
				(*ngraph)->nbs[k][j] = remap[graph->nbs[i][j]];
			}
		} else {
			(*ngraph)->nbs[k] = NULL;
		}
	}

	/*
	 * rearrange into linear storage and free temporary storage 
	 */
	for (i = 0, nnb = 0; i < (*ngraph)->n; i++) {
		nnb += (*ngraph)->nnbs[i];
	}
	if (nnb) {
		hold = Calloc(nnb, int);
	} else {
		hold = NULL;
	}

	for (i = 0, indx = 0; i < (*ngraph)->n; i++) {
		if ((*ngraph)->nnbs[i]) {
			memcpy(&hold[indx], (*ngraph)->nbs[i], (size_t) ((*ngraph)->nnbs[i] * sizeof(int)));
			Free((*ngraph)->nbs[i]);
			(*ngraph)->nbs[i] = &hold[indx];
		} else {
			Free((*ngraph)->nbs[i]);
		}
		indx += (*ngraph)->nnbs[i];
	}
	GMRFLib_EWRAP0(GMRFLib_prepare_graph(*ngraph));

	return GMRFLib_SUCCESS;
}

/*!
  \brief Returns a copy of a graph

  \param[out] graph_new A \c GMRFLib_graph_tp -object. At output, it contains a copy of \em graph.
  \param[in] graph_old A \c GMRFLib_graph_tp -object.
 */
int GMRFLib_copy_graph(GMRFLib_graph_tp ** graph_new, GMRFLib_graph_tp * graph_old)
{
	/*
	 * there is no need to do call _prepare_graph is the old graph is assumed to be ok. 
	 */
	int m, i, n, *hold = NULL, hold_idx;
	GMRFLib_graph_tp *g = NULL;

	GMRFLib_ENTER_ROUTINE;
	if (!graph_old) {
		*graph_new = NULL;
		GMRFLib_LEAVE_ROUTINE;
		return GMRFLib_SUCCESS;
	}

	if (0) {
		/*
		 * old slow code 
		 */
		n = graph_old->n;
		GMRFLib_compute_subgraph(&g, graph_old, NULL);
		if (graph_old->mothergraph_idx) {
			if (!g->mothergraph_idx)
				g->mothergraph_idx = Calloc(n, int);
			memcpy(g->mothergraph_idx, graph_old->mothergraph_idx, (size_t) (n * sizeof(int)));
		}
		*graph_new = g;
	} else {
		/*
		 * new and better code 
		 */
		GMRFLib_make_empty_graph(&g);
		g->n = n = graph_old->n;
		g->nnbs = Calloc(n, int);
		memcpy(g->nnbs, graph_old->nnbs, (size_t) (n * sizeof(int)));

		for (i = m = 0; i < n; i++) {
			m += g->nnbs[i];
		}
		if (m) {
			hold = Calloc(m, int);
			g->nbs = Calloc(n, int *);

			for (i = hold_idx = 0; i < n; i++) {
				if (g->nnbs[i]) {
					g->nbs[i] = &hold[hold_idx];
					memcpy(g->nbs[i], graph_old->nbs[i], (size_t) (g->nnbs[i] * sizeof(int)));
					hold_idx += g->nnbs[i];
				}
			}
		}
		if (graph_old->mothergraph_idx) {
			if (!g->mothergraph_idx) {
				g->mothergraph_idx = Calloc(n, int);
			}
			memcpy(g->mothergraph_idx, graph_old->mothergraph_idx, (size_t) (n * sizeof(int)));
		}
		*graph_new = g;
	}

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

GMRFLib_sizeof_tp GMRFLib_sizeof_graph(GMRFLib_graph_tp * graph)
{
	/*
	 * return, approximately, the sizeof GRAPH 
	 */

	if (!graph) {
		return 0;
	}

	GMRFLib_sizeof_tp siz = 0;
	int i, m, n;

	n = graph->n;
	for (i = m = 0; i < n; i++) {
		m += graph->nnbs[i];
	}

	siz += sizeof(int) + m * sizeof(int) + n * sizeof(int) + n * sizeof(int *);

	if (graph->mothergraph_idx) {
		siz += n * sizeof(int);
	}

	return siz;
}

/*!
  \brief Computes a subgraph of a graph by removing the nodes indicated by the \em remove_flag argument.

  \param[out] subgraph The subgraph generated by removing the nodes for which \em remove_flag is 1
    from the graph \em graph.
  \param[in] graph The graph from which to compute the subgraph.
  \param[in] remove_flag An array of length \em n, the number of nodes in \em graph, specifying
    which nodes to be removed creating the subgraphs. The nodes for which \em remove_flags is \em 0,
    IS included in the subgraph.

  \par Example:
  See \ref ex_graph
 */
int GMRFLib_compute_subgraph(GMRFLib_graph_tp ** subgraph, GMRFLib_graph_tp * graph, char *remove_flag)
{

	/*
	 * return a subgraph of graph, by ruling out those nodes for which remove_flag[i] is true, keeping those which
	 * remove_flag[i] false. 
	 */
	int i, j, nneig, nn, k, n_neig_tot, storage_indx, *sg_iidx, *storage = NULL, free_remove_flag = 0;

	GMRFLib_ENTER_ROUTINE;

	if (!graph) {
		GMRFLib_LEAVE_ROUTINE;
		return GMRFLib_SUCCESS;
	}

	/*
	 * if the graph is empty, then just return an empty graph 
	 */
	if (graph->n == 0) {
		*subgraph = Calloc(1, GMRFLib_graph_tp);
		(*subgraph)->n = 0;
		GMRFLib_LEAVE_ROUTINE;

		return GMRFLib_SUCCESS;
	}

	/*
	 * to ease the interface: remove_flag = NULL is ok. then this just do a (slow) copy of the graph. 
	 */
	if (!remove_flag) {
		remove_flag = Calloc(graph->n, char);

		free_remove_flag = 1;
	}

	GMRFLib_make_empty_graph(subgraph);

	for (i = 0, nn = 0; i < graph->n; i++) {
		nn += (!remove_flag[i]);
	}
	(*subgraph)->n = nn;
	if (!((*subgraph)->n)) {
		GMRFLib_LEAVE_ROUTINE;
		return GMRFLib_SUCCESS;
	}

	/*
	 * create space 
	 */
	(*subgraph)->nnbs = Calloc((*subgraph)->n, int);
	(*subgraph)->nbs = Calloc((*subgraph)->n, int *);

	(*subgraph)->mothergraph_idx = Calloc((*subgraph)->n, int);
	sg_iidx = Calloc(graph->n, int);

	/*
	 * make the mapping of nodes 
	 */
	for (i = 0, k = 0; i < graph->n; i++) {
		if (!remove_flag[i]) {
			(*subgraph)->mothergraph_idx[k] = i;
			sg_iidx[i] = k;
			k++;
		} else {
			sg_iidx[i] = -1;		       /* to force a failure if used wrong */
		}
	}

	/*
	 * parse the graph and collect nodes not to be removed. 
	 */
	for (i = 0, k = 0, n_neig_tot = 0; i < graph->n; i++) {
		if (!remove_flag[i]) {
			for (j = 0, nneig = 0; j < graph->nnbs[i]; j++) {
				nneig += (!remove_flag[graph->nbs[i][j]]);
			}
			n_neig_tot += nneig;
			(*subgraph)->nnbs[k] = nneig;

			if (nneig > 0) {
				(*subgraph)->nbs[k] = Calloc(nneig, int);

				for (j = 0, nneig = 0; j < graph->nnbs[i]; j++) {
					if (!remove_flag[graph->nbs[i][j]]) {
						(*subgraph)->nbs[k][nneig] = sg_iidx[graph->nbs[i][j]];
						nneig++;
					}
				}
			} else {
				(*subgraph)->nbs[k] = NULL;
			}

			k++;
		}
	}

	Free(sg_iidx);

	/*
	 * map the graph to a more computational convenient memory layout! use just one long vector to store all the neighbours. 
	 */
	if (n_neig_tot) {
		storage = Calloc(n_neig_tot, int);

		storage_indx = 0;
		for (i = 0; i < (*subgraph)->n; i++) {
			if ((*subgraph)->nnbs[i]) {
				memcpy(&storage[storage_indx], (*subgraph)->nbs[i], (size_t) (sizeof(int) * (*subgraph)->nnbs[i]));
				Free((*subgraph)->nbs[i]);
				(*subgraph)->nbs[i] = &storage[storage_indx];
				storage_indx += (*subgraph)->nnbs[i];
			} else {
				(*subgraph)->nbs[i] = NULL;
			}
		}
	}
	GMRFLib_EWRAP0(GMRFLib_prepare_graph(*subgraph));

	if (free_remove_flag) {
		Free(remove_flag);			       /* if we have used our own */
	}

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}
int GMRFLib_convert_to_mapped(double *destination, double *source, GMRFLib_graph_tp * graph, int *remap)
{
	/*
	 * convert from the real-world to the mapped world. source might be NULL. 
	 */
	int i;

	if ((destination && source) && (destination != source)) {
		for (i = 0; i < graph->n; i++) {
			destination[remap[i]] = source[i];
		}
	} else {
		double *work = Malloc(graph->n, double);
		memcpy(work, destination, graph->n * sizeof(double));
		for (i = 0; i < graph->n; i++) {
			destination[remap[i]] = work[i];
		}
		Free(work);
	}
	return GMRFLib_SUCCESS;
}
int GMRFLib_convert_from_mapped(double *destination, double *source, GMRFLib_graph_tp * graph, int *remap)
{
	/*
	 * convert from the mapped-world to the real world. source might be NULL. 
	 */

	int i;

	if ((destination && source) && (destination != source)) {
		for (i = 0; i < graph->n; i++) {
			destination[i] = source[remap[i]];
		}
	} else {
		double *work = Malloc(graph->n, double);
		memcpy(work, destination, graph->n * sizeof(double));
		for (i = 0; i < graph->n; i++) {
			destination[i] = work[remap[i]];
		}
		Free(work);
	}
	return GMRFLib_SUCCESS;
}

/*!
  \brief  Compute \f$ \mbox{\boldmath $Q$}\mbox{\boldmath $x$} \f$.
  
  Computes \f$ \mbox{\boldmath $Q$} \mbox{\boldmath $x$} \f$, the matrix-vector product of the
  precision matrix <em>\b Q</em> and an array <em>\b x</em>. The elements of <em>\b x</em> should be
  in the original ordering of the nodes of the graph.

  \param[out] result At output, the length \em n array \em result contains 
  the value of the product  \f$ \mbox{\boldmath $Q$}\times \mbox{\boldmath $x$} \f$.
  \param[in] x The array <em>\b x</em>, as a length \em n array.
  \param[in] graph The graph on which <em>\b x</em> is defined.
  \param[in] Qfunc The function defining the precision matrix <em>\b Q</em>.
  \param[in] Qfunc_arg A \em void -pointer defining the
    address of a variable or data structure holding the arguments to
    the function \em Qfunc.

  \sa GMRFLib_Qfunc_tp, GMRFLib_xQx
 */
int GMRFLib_Qx(double *result, double *x, GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg)
{
	/*
	 * compute RESULT = Q*x, (RESULT is a vector).
	 */
	int i, id;

	id = GMRFLib_thread_id;
	memset(result, 0, graph->n * sizeof(double));

#pragma omp parallel for private(i)
	for (i = 0; i < graph->n; i++) {
		int j, jj;

		GMRFLib_thread_id = id;
		result[i] += Qfunc(i, i, Qfunc_arg) * x[i];
		for (j = 0; j < graph->nnbs[i]; j++) {
			jj = graph->nbs[i][j];
			result[i] += Qfunc(i, jj, Qfunc_arg) * x[jj];
		}
	}
	GMRFLib_thread_id = id;

	return GMRFLib_SUCCESS;
}
int GMRFLib_print_Qfunc(FILE * fp, GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg)
{
	int i, j, jj;

	if (!fp) {
		fp = stdout;
	}
	for (i = 0; i < graph->n; i++) {
		fprintf(fp, "Q[ %1d , %1d ] = %.12f\n", i, i, Qfunc(i, i, Qfunc_arg));
		for (j = 0; j < graph->nnbs[i]; j++) {
			jj = graph->nbs[i][j];
			fprintf(fp, "\tQ[ %1d , %1d ] = %.12f\n", i, jj, Qfunc(i, jj, Qfunc_arg));
		}
	}
	return GMRFLib_SUCCESS;
}

/*!
  \brief  Compute \f$ \mbox{\boldmath $x^TQx$} \f$.
  
  Compute \f$ \mbox{\boldmath $x^TQx$} \f$: the scalar product of <em>\b x</em>, and the matrix-vector product of the precision
  matrix <em>\b Q</em> and an array <em>\b x</em>. The elements of <em>\b x</em> should be in the original ordering of the nodes
  of the graph.

  \param[out] result The result is returned as \a *result.
  \param[in] x The array <em>\b x</em>, as a length \em n array.
  \param[in] graph The graph on which <em>\b x</em> is defined.
  \param[in] Qfunc The function defining the precision matrix <em>\b Q</em>.
  \param[in] Qfunc_arg A \em void -pointer defining the
    address of a variable or data structure holding the arguments to
    the function \em Qfunc.

  \sa GMRFLib_xQx
 */
int GMRFLib_xQx(double *result, double *x, GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg)
{
	int i;
	double *y = NULL, res;

	y = Calloc(graph->n, double);

	GMRFLib_Qx(y, x, graph, Qfunc, Qfunc_arg);
	for (i = 0, res = 0.0; i < graph->n; i++) {
		res += y[i] * x[i];
	}
	Free(y);

	*result = res;
	return GMRFLib_SUCCESS;
}

/*!  \brief Creates a graph on a lattice, customizing the graph for use in other library functions.

  \param[out] graph At output, <em>(*graph)</em> contains
  the specification of the lattice graph, customized for use in
  other library routines.
  \param[in] nrow, ncol The number of rows ( \f$ n_{row} \f$ ) and columns
  ( \f$ n_{col} \f$ ) in the \f$ n_{row}\times n_{col} \f$ lattice.
  \param[in] nb_row, nb_col Specification of the parameters \f$ nb_{row} \f$ 
  and \f$ nb_{col} \f$ of the \f$ (2 nb_{row}+1)\times(2 nb_{col}+1) \f$
  -neighbourhood of the lattice.
  \param[in] cyclic_flag If this argument is <em>>0</em>, the graph is made
  cyclic.

  \par Example:
  See \ref ex_graph

  \sa GMRFLib_lattice2node, GMRFLib_node2lattice
 */
int GMRFLib_make_lattice_graph(GMRFLib_graph_tp ** graph, int nrow, int ncol, int nb_row, int nb_col, int cyclic_flag)
{
	/*
	 * make an lattice graph, with nodes from 0...n-1, and a (2 x nb_row +1) x (2 x nb_col + 1) neighborhood. if
	 * `cyclic_flag`, then make the graph cyclic. the prune_graph function might be useful to derive non-square
	 * neighborhoods. 
	 */
	int *hold = NULL, n, nnb, ir, ic, node, nnode, irow, icol;

	nb_row = IMIN(nrow - 1, nb_row);		       /* otherwise (i+nrow)%nrow will fail if cyclic */
	nb_col = IMIN(ncol - 1, nb_col);		       /* same here */
	n = ncol * nrow;
	nnb = (2 * nb_row + 1) * (2 * nb_col + 1);

	GMRFLib_make_empty_graph(graph);
	(*graph)->n = n;
	(*graph)->nnbs = Calloc(n, int);
	(*graph)->nbs = Calloc(n, int *);

	/*
	 * if nh is large compared to n, the this graph may contain double (or more) occurences of one neigbor, but the
	 * `prepare_graph' function fix this. 
	 */
	if (nnb) {
		hold = Calloc(n * nnb, int);		       /* use a linear storage */

		for (node = 0; node < n; node++) {
			(*graph)->nbs[node] = &hold[node * nnb];	/* set pointers to it */
		}

		for (node = 0; node < n; node++) {
			GMRFLib_node2lattice(node, &irow, &icol, nrow, ncol);
			if (cyclic_flag) {
				for (ir = irow - nb_row; ir <= irow + nb_row; ir++)
					for (ic = icol - nb_col; ic <= icol + nb_col; ic++) {
						GMRFLib_lattice2node(&nnode, (ir + nrow) % nrow, (ic + ncol) % ncol, nrow, ncol);
						if (nnode != node) {
							(*graph)->nbs[node][(*graph)->nnbs[node]++] = nnode;
						}
					}
			} else {
				for (ir = IMAX(0, irow - nb_row); ir <= IMIN(nrow - 1, irow + nb_row); ir++) {
					for (ic = IMAX(0, icol - nb_col); ic <= IMIN(ncol - 1, icol + nb_col); ic++) {
						GMRFLib_lattice2node(&nnode, ir, ic, nrow, ncol);
						if (nnode != node) {
							(*graph)->nbs[node][(*graph)->nnbs[node]++] = nnode;
						}
					}
				}
			}
		}
	} else {
		for (node = 0; node < n; node++) {
			(*graph)->nbs[node] = NULL;	       /* if zero neighbours, the pointer should be zero */
		}
	}

	GMRFLib_EWRAP0(GMRFLib_prepare_graph(*graph));

	return GMRFLib_SUCCESS;
}

/*!
  \brief  Return the node number of a point in a lattice.

  In the case of a rectangular neighbourhood, the ordering of the
  nodes will be the one minimizing the bandwidth of the precision
  matrix of a GMRF on the graph.
  
  \param[out] node The node index \f$ (\in [0,n_{row}\times n_{col}-1]) \f$
    corresponding to the row and column indices of the graph.
  \param[in] irow Row index \f$ (\in [0,n_{row}-1]) \f$.
  \param[in] icol Column index \f$ (\in [0,n_{col}-1]) \f$.
  \param[in] nrow The number of rows of the lattice.
  \param[in] ncol The number of columns of the lattice.

  \remarks The node index \c node is computed as follows: <em>node = icol + irow*ncol</em>, so column-wise storage is used.

  \par Example:
  See \ref ex_graph
  \sa GMRFLib_node2lattice, GMRFLib_make_lattice_graph
 */
int GMRFLib_lattice2node(int *node, int irow, int icol, int nrow, int ncol)
{
	// *node = icol + irow * ncol;
	*node = irow + icol * nrow;
	return GMRFLib_SUCCESS;
}

/*!
  \brief  Compute the lattice indices of a node
  
  \param[in] node The node index \f$ (\in [0,n_{row}\times n_{col}-1]) \f$
    corresponding to the row and column indices of the graph.
  \param[out] irow Row index \f$ (\in [0,n_{row}-1]) \f$ corresponding
  to node number \em node in the graph.
  \param[out] icol Column index \f$ (\in [0,n_{col}-1]) \f$ corresponding
  to node number \em node in the graph.
  \param[in] nrow The number of rows of the lattice.
  \param[in] ncol The number of columns of the lattice.

  \par Example:
  See \ref ex_graph
  \sa GMRFLib_lattice2node, GMRFLib_make_lattice_graph
 */
int GMRFLib_node2lattice(int node, int *irow, int *icol, int nrow, int ncol)
{
	// *irow = node / ncol;
	// *icol = node - (*irow) * ncol;

	*icol = node / nrow;
	*irow = node - (*icol) * nrow;

	return GMRFLib_SUCCESS;
}

/*!

  \brief Returns a new graph which is a copy of graph but where entries where the elements of the
  precision matrix <em>Q(i,j)=0, i~j,</em> are removed.

  \param[out] new_graph A \c GMRFLib_graph_tp -object. At
  output, it contains a copy of \em graph but all entries
  where <em>Q(i,j)=0</em> are removed.
  \param[in] graph At output, <em>(*graph)</em> has been initialized 
  by using the information on the file \em filename. 
  \param[in] Qfunc A function defining the elements of the \em \b Q -matrix.
  \param[in] Qfunc_arg A \em void -pointer defining the
  address of a variable or data structure holding the arguments to
  the function \em Qfunc.

  \sa GMRFLib_init_wa_problem, GMRFLib_init_nwa_problem.
 */
int GMRFLib_prune_graph(GMRFLib_graph_tp ** new_graph, GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg)
{

	/*
	 * return a modified graph by removing entries where Q(i,j)=0, i~j. 
	 */
	int i, j, k, ii, *free_ptr = NULL, found;

	GMRFLib_copy_graph(new_graph, graph);

	/*
	 * this is a bit tricky. as GMRFLib_free_graph free's the first ptr where nnbs[i]>0, then we must make sure that this
	 * is the same ptr after pruning, as we modify `new_graph'. so, we need to store that ptr and make sure its ok after
	 * pruning. 
	 */
	for (i = 0; i < (*new_graph)->n; i++) {
		if ((*new_graph)->nnbs[i]) {
			free_ptr = (*new_graph)->nbs[i];
			break;
		}
	}

	for (i = 0; i < graph->n; i++) {
		if (graph->nnbs[i]) {
			for (j = 0, k = 0; j < graph->nnbs[i]; j++) {
				ii = graph->nbs[i][j];
				if (Qfunc(i, ii, Qfunc_arg) != 0.0) {
					(*new_graph)->nbs[i][k] = (*new_graph)->nbs[i][j];
					k++;
				}
			}
			(*new_graph)->nnbs[i] = k;
			if ((*new_graph)->nnbs[i] == 0)
				(*new_graph)->nbs[i] = NULL;
		}
	}
	GMRFLib_EWRAP0(GMRFLib_prepare_graph(*new_graph));

	for (i = 0, found = 0; i < (*new_graph)->n; i++) {
		if ((*new_graph)->nnbs[i]) {
			if (graph->nbs[i] != free_ptr) {
				/*
				 * 
				 * *pos = &((*new_graph)->nnbs[i]);
				 * 
				 */
				if (0) {
					GMRFLib_msg(stdout, "\n\nNEW CODE HERE, HOPE ITS OK\n\n");
				}

				for (j = 0; j < (*new_graph)->nnbs[i]; j++) {
					free_ptr[j] = (*new_graph)->nbs[i][j];
				}
				(*new_graph)->nbs[i] = free_ptr;
			}
			found = 1;
			break;
		}
	}
	if (!found) {
		/*
		 * then all neighbours are gone, and we have to free free_ptr ourself. 
		 */
		Free(free_ptr);
	}

	return GMRFLib_SUCCESS;
}

/*!
  \brief Creates a linear graph, corresponding to an autoregressive AR(\em p)-model

  Make an linear graph, with nodes from <em>0...n-1</em>, and \em bw neighbours in each
  direction. hence, <em>bw = p</em>, makes the graph for an auto-regressive process or order \em
  p. if `cyclic_flag`, then make the graph cyclic.

  \param graph At output, (*graph) contains the specifications of the 
  linear graph, customized for use in other library routines.
  \param n The number of nodes in the graph.
  \param bw Specifies the bandwidth of the graph, such that each node 
  has \em bw neighbours in each direction. A bandwidth of \em bw 
  corresponds to an AR(\em bw)-model.
  \param cyclic_flag If this argument is > 0, the graph is made cyclic.

  \par Example:
  See \ref ex_graph
 */
int GMRFLib_make_linear_graph(GMRFLib_graph_tp ** graph, int n, int bw, int cyclic_flag)
{
	int i, j, k, *hold = NULL;

	bw = IMIN(n - 1, IMAX(0, bw));
	GMRFLib_make_empty_graph(graph);
	(*graph)->n = n;
	(*graph)->nnbs = Calloc(n, int);
	(*graph)->nbs = Calloc(n, int *);

	if (bw) {
		hold = Calloc(n * 2 * bw, int);		       /* use a linear storage */

		for (i = 0; i < n; i++) {
			(*graph)->nbs[i] = &hold[i * 2 * bw];  /* set pointers to it */
		}

		if (cyclic_flag) {
			for (i = 0; i < n; i++) {
				(*graph)->nnbs[i] = 2 * bw;
				for (j = i - bw, k = 0; j <= i + bw; j++) {
					if (j != i) {
						(*graph)->nbs[i][k++] = MOD(j, n);
					}
				}
			}
		} else {
			for (i = 0; i < n; i++) {
				(*graph)->nnbs[i] = (i - IMAX(i - bw, 0)) + (IMIN(n - 1, i + bw) - i);
				for (j = i - bw, k = 0; j <= i + bw; j++) {
					if (j != i && LEGAL(j, n)) {
						(*graph)->nbs[i][k++] = j;
					}
				}
				if (!k) {
					(*graph)->nbs[i] = NULL;	/* if zero neighbours, the pointer should be zero */
				}
			}
		}
	} else {
		for (i = 0; i < n; i++) {
			(*graph)->nbs[i] = NULL;	       /* if zero neighbours, the pointer should be zero */
		}
	}

	GMRFLib_EWRAP0(GMRFLib_prepare_graph(*graph));
	return GMRFLib_SUCCESS;
}


/* NOT DOCUMENTED
   
  \brief Internal function in \c GMRFLib_nfold_graph().

  Returns a new graph created by expanding one graph \em g with another 
  graph \em  gg. The two graphs must have the same number of nodes.

  \param[out] ng The new expanded graph. The pointer <em>(*ng)</em> will be 
  allocated within the routine.
  \param[in] g, gg The graphs to be combined.

  \remarks Given two graphs \em g and \em gg, the new graph \em ng is 
  created by expanding the neighbourhood of \em g by the neighbourhood 
  of \em gg. That is, for each neighbour \f$ j_g \f$ of a node \f$ i_g \f$ 
  of \em g, a node \f$ k_{gg} \f$ in \em gg is included in the 
  neighbourhood of node \f$ i_g \f$ in \em g if \f$ k_{gg} \f$ and 
  \f$ j_g \f$ are neighbours in \em gg. The new graph \em ng is equal to the
  expanded version of of \em g.

  \sa GMRFLib_nfold_graph
 */
int GMRFLib_fold_graph(GMRFLib_graph_tp ** ng, GMRFLib_graph_tp * g, GMRFLib_graph_tp * gg)
{
	/*
	 * return ng = g * gg, meaning that ng is g expanded by gg 
	 */
	int i, ii, j, jj, k, kk, nneig, nnb, *hold, indx;
	GMRFLib_graph_tp *newg = NULL;

	/*
	 * first the easy cases 
	 */
	if (!g && !gg) {
		*ng = NULL;
		return GMRFLib_SUCCESS;
	}
	if ((g && !gg) || (gg && !g)) {
		GMRFLib_copy_graph(ng, (g ? g : gg));
		return GMRFLib_SUCCESS;
	}

	GMRFLib_ASSERT(g->n == gg->n, GMRFLib_EPARAMETER);     /* of same size? */

	GMRFLib_make_empty_graph(&newg);
	newg->n = g->n;
	newg->nnbs = Calloc(g->n, int);
	newg->nbs = Calloc(g->n, int *);

	for (i = 0; i < newg->n; i++) {
		/*
		 * count number of neighbours neighbours 
		 */
		for (j = 0, nneig = g->nnbs[i]; j < g->nnbs[i]; j++) {
			nneig += gg->nnbs[g->nbs[i][j]];
		}
		if (nneig) {
			newg->nbs[i] = Calloc(nneig, int);

			newg->nnbs[i] = 0;
			for (j = 0, k = 0; j < g->nnbs[i]; j++) {
				jj = newg->nbs[i][k++] = g->nbs[i][j];
				for (ii = 0; ii < gg->nnbs[jj]; ii++) {
					kk = gg->nbs[jj][ii];
					if (kk != i) {
						newg->nbs[i][k++] = kk;
					}
				}
			}
			newg->nnbs[i] = k;

			/*
			 * make them unique 
			 */
			qsort((void *) newg->nbs[i], (size_t) newg->nnbs[i], sizeof(int), GMRFLib_icmp);
			for (j = k = 0; j < newg->nnbs[i]; j++) {
				if (j == 0 || newg->nbs[i][j] != newg->nbs[i][k - 1]) {
					newg->nbs[i][k++] = newg->nbs[i][j];
				}
			}
			newg->nnbs[i] = k;
		} else {
			newg->nnbs[i] = 0;
			newg->nbs[i] = NULL;
		}
	}

	/*
	 * rearrange into linear storage and free temporary storage 
	 */
	for (i = 0, nnb = 0; i < newg->n; i++) {
		nnb += newg->nnbs[i];
	}

	if (nnb) {
		hold = Calloc(nnb, int);
	} else
		hold = NULL;

	for (i = 0, indx = 0; i < newg->n; i++) {
		if (newg->nnbs[i]) {
			memcpy(&hold[indx], newg->nbs[i], newg->nnbs[i] * sizeof(int));
			Free(newg->nbs[i]);
			newg->nbs[i] = &hold[indx];
		} else {
			Free(newg->nbs[i]);
		}
		indx += newg->nnbs[i];
	}
	GMRFLib_EWRAP0(GMRFLib_prepare_graph(newg));
	*ng = newg;

	return GMRFLib_SUCCESS;
}

/*!
  \brief Returns a new graph created by multiple expanding of one graph by itself.

  <tt>nfold = 0: ng = I \n
  nfold = 1: ng = og \n
  nfold = 2: ng = og*og \n
  nfold = 3: ng = og*og*og \n
  etc </tt>

  \param[out] ng The new expanded graph. The pointer <em>(*ng)</em> will be 
  allocated within the routine.
  \param[in] og The graph to be expanded.
  \param[in] nfold The number of times to expand the graph.

  \remarks Given the old graph \em og, the new graph \em ng is created 
  by expanding the neighbours of \em og by it's own neighbourhood 
  <em>nfold-1</em> times. The function GMRFLib_fold_graph is called 
  successively <em>nfold-1</em> times, the first time using 
  <em>g = gg = og</em>, and in successive runs letting <em>g = gprev</em> 
  and <em>gg = og</em>, where \em gprev is the resulting graph of 
  the previous run.  For completeness, <em>nfold = 0</em> is allowed for, 
  meaning that a graph with \em n independent nodes is returned as 
  \em ng, where \em n is the number of nodes in the old graph \em og. If
  <em>nfold = 1</em>, the new graph is a copy of the old graph.

  \par Example
  See \ref ex_graph
 */
int GMRFLib_nfold_graph(GMRFLib_graph_tp ** ng, GMRFLib_graph_tp * og, int nfold)
{
	/*
	 * make new graph, 'ng', that is 'nfold' of 'og'
	 * 
	 * nfold = 0: ng = I nfold = 1: ng = og nfold = 2: ng = og*og nfold = 3: ng = og*og*og etc
	 * 
	 */
	int i;
	GMRFLib_graph_tp *newg = NULL, *oldg = NULL;

	if (nfold == 0) {
		GMRFLib_make_empty_graph(&newg);
		newg->n = og->n;
		newg->nnbs = Calloc(newg->n, int);
		newg->nbs = Calloc(newg->n, int *);

		GMRFLib_EWRAP0(GMRFLib_prepare_graph(newg));
	} else if (nfold == 1) {
		GMRFLib_copy_graph(&newg, og);
	} else {
		for (i = 0, oldg = newg = NULL; i < nfold; i++) {
			GMRFLib_fold_graph(&newg, oldg, og);
			GMRFLib_free_graph(oldg);
			if (i < nfold - 1) {
				oldg = newg;
				newg = NULL;
			}
		}
	}

	*ng = newg;
	return GMRFLib_SUCCESS;
}

/*! 
  \brief Return a new graph which is the union of n_graphs graphs.

  <em>i~j</em> in union_graph, if <em>i~j</em> in
  <em>graph_array[0]...graph_array[n_graphs-1]</em> \n\n
*/
int GMRFLib_union_graph(GMRFLib_graph_tp ** union_graph, GMRFLib_graph_tp ** graph_array, int n_graphs)
{
	/*
	 * return a new graph which is the union of n_graphs graphs: i~j in union_graph, if i~j in
	 * graph_array[0]...graph_array[n_graphs-1] 
	 */

	int i, k, node, nnbs, idx, *hold = NULL, hold_idx;
	GMRFLib_graph_tp *tmp_graph;

	if (!graph_array || n_graphs <= 0) {
		*union_graph = NULL;
		return GMRFLib_SUCCESS;
	}
	for (k = 1; k < n_graphs; k++)
		GMRFLib_ASSERT(graph_array[0]->n == graph_array[k]->n, GMRFLib_EPARAMETER);

	GMRFLib_make_empty_graph(union_graph);
	(*union_graph)->n = graph_array[0]->n;
	(*union_graph)->nnbs = Calloc((*union_graph)->n, int);
	(*union_graph)->nbs = Calloc((*union_graph)->n, int *);

	for (node = 0, nnbs = 0; node < (*union_graph)->n; node++) {
		for (k = 0; k < n_graphs; k++) {
			nnbs += graph_array[k]->nnbs[node];
		}
	}

	if (nnbs) {
		hold = Calloc(nnbs, int);
	} else {
		hold = NULL;
	}

	for (node = 0, hold_idx = 0; node < (*union_graph)->n; node++) {
		for (k = 0, nnbs = 0; k < n_graphs; k++) {
			nnbs += graph_array[k]->nnbs[node];
		}
		if (nnbs) {
			(*union_graph)->nnbs[node] = nnbs;     /* will include multiple counts */
			(*union_graph)->nbs[node] = &hold[hold_idx];

			for (k = 0, idx = 0; k < n_graphs; k++) {
				for (i = 0; i < graph_array[k]->nnbs[node]; i++) {
					(*union_graph)->nbs[node][idx++] = graph_array[k]->nbs[node][i];
				}
			}

			qsort((*union_graph)->nbs[node], (size_t) (*union_graph)->nnbs[node], sizeof(int), GMRFLib_icmp);
			for (i = 1, nnbs = 0; i < (*union_graph)->nnbs[node]; i++) {
				if ((*union_graph)->nbs[node][i] != (*union_graph)->nbs[node][nnbs]) {
					(*union_graph)->nbs[node][++nnbs] = (*union_graph)->nbs[node][i];
				}
			}
			(*union_graph)->nnbs[node] = nnbs + 1;
			hold_idx += (*union_graph)->nnbs[node];
		} else {
			(*union_graph)->nnbs[node] = 0;
			(*union_graph)->nbs[node] = NULL;
		}
	}
	GMRFLib_EWRAP0(GMRFLib_prepare_graph(*union_graph));   /* this is required */

	/*
	 * the union_graph is now (probably) to large as it acounts for multiple counts. the easiest way out of this, is to
	 * make (and use) a new copy, then free the current one. 
	 */
	GMRFLib_copy_graph(&tmp_graph, *union_graph);
	GMRFLib_free_graph(*union_graph);
	*union_graph = tmp_graph;

	return GMRFLib_SUCCESS;
}

/*!
  \brief Add missing neighbours to a uncomplete or malformed graph; if <em>i~j</em> but <em>j!~i</em>, then add <em>j~i</em>.

  \param[out] n_graph A pointer to a new completed graph, buildt
    from \em graph.
  \param[in] graph The incomplete graph.

  \note If \em graph is ok in the sense that <em>i~j</em>
    then <em>j~i</em>, then this routine just returns a copy of the
    graph, allthough it will be a bit slow.
 */
int GMRFLib_complete_graph(GMRFLib_graph_tp ** n_graph, GMRFLib_graph_tp * graph)
{
	/*
	 * return a new graph that is complete: if i~j but j!~i, then add j~i 
	 */

	int i, ii, j, k, *neigh_size, n_neig_tot, *storage = NULL, storage_idx;

	if (!graph) {
		*n_graph = NULL;
		return GMRFLib_SUCCESS;
	}

	/*
	 * setup new graph 
	 */
	GMRFLib_make_empty_graph(n_graph);
	neigh_size = Calloc(graph->n, int);
	(*n_graph)->n = graph->n;
	(*n_graph)->nbs = Calloc(graph->n, int *);
	(*n_graph)->nnbs = Calloc(graph->n, int);

	for (i = 0; i < (*n_graph)->n; i++) {
		(*n_graph)->nnbs[i] = graph->nnbs[i];
		neigh_size[i] = IMAX(16, 2 * (*n_graph)->nnbs[i]);
		(*n_graph)->nbs[i] = Calloc(neigh_size[i], int);

		if ((*n_graph)->nnbs[i]) {
			memcpy((*n_graph)->nbs[i], graph->nbs[i], graph->nnbs[i] * sizeof(int));
		}
	}

	/*
	 * if i ~ j, then add j ~ i 
	 */
	for (i = 0; i < (*n_graph)->n; i++) {
		for (j = 0; j < (*n_graph)->nnbs[i]; j++) {
			ii = (*n_graph)->nbs[i][j];
			if ((*n_graph)->nnbs[ii] == neigh_size[ii]) {
				neigh_size[ii] *= 2;
				(*n_graph)->nbs[ii] = Realloc((*n_graph)->nbs[ii], neigh_size[ii], int);
			}
			(*n_graph)->nbs[ii][(*n_graph)->nnbs[ii]++] = i;
		}
	}

	/*
	 * some may appear more than once, count number of unique neighbours 
	 */
	for (i = 0, n_neig_tot = 0; i < (*n_graph)->n; i++) {
		if ((*n_graph)->nnbs[i]) {
			qsort((*n_graph)->nbs[i], (size_t) (*n_graph)->nnbs[i], sizeof(int), GMRFLib_icmp);
			for (j = 1, k = 0; j < (*n_graph)->nnbs[i]; j++) {
				if ((*n_graph)->nbs[i][j] != (*n_graph)->nbs[i][k]) {
					(*n_graph)->nbs[i][++k] = (*n_graph)->nbs[i][j];
				}
			}
			(*n_graph)->nnbs[i] = k + 1;
			n_neig_tot += (*n_graph)->nnbs[i];
		}
	}

	if (n_neig_tot) {
		storage = Calloc(n_neig_tot, int);

		for (i = 0, storage_idx = 0; i < (*n_graph)->n; i++) {
			if ((*n_graph)->nnbs[i]) {
				memcpy(&storage[storage_idx], (*n_graph)->nbs[i], sizeof(int) * (*n_graph)->nnbs[i]);
				Free((*n_graph)->nbs[i]);
				(*n_graph)->nbs[i] = &storage[storage_idx];
				storage_idx += (*n_graph)->nnbs[i];
			} else {
				Free((*n_graph)->nbs[i]);
			}
		}
	}
	GMRFLib_EWRAP0(GMRFLib_prepare_graph(*n_graph));

	Free(neigh_size);

	return GMRFLib_SUCCESS;
}
int GMRFLib_offset_graph(GMRFLib_graph_tp ** new_graph, int n_new, int offset, GMRFLib_graph_tp * graph)
{
	/*
	 * insert graph into a larger graph with n_new nodes such that node `i' in graph corresponds to node `offset +i' in the 
	 * larger graph. n_new must be larger than graph->n, and `offset' such that `0<= offset', and `offset + graph->n <
	 * new_graph->n' 
	 */

	int i, ii, j, n_neig, *hold = NULL, hold_idx;
	GMRFLib_graph_tp *g = NULL;

	for (i = 0, n_neig = 0; i < graph->n; i++) {
		n_neig += graph->nnbs[i];
	}

	GMRFLib_make_empty_graph(&g);
	g->n = n_new;
	g->nnbs = Calloc(n_new, int);
	g->nbs = Calloc(n_new, int *);

	if (n_neig) {
		hold = Calloc(n_neig, int);
	} else {
		hold = NULL;
	}

	for (i = 0, hold_idx = 0; i < graph->n; i++) {
		ii = i + offset;
		g->nnbs[ii] = graph->nnbs[i];
		g->nbs[ii] = &hold[hold_idx];
		hold_idx += graph->nnbs[i];

		for (j = 0; j < graph->nnbs[i]; j++) {
			g->nbs[ii][j] = graph->nbs[i][j] + offset;
		}
	}
	GMRFLib_EWRAP0(GMRFLib_prepare_graph(g));

	*new_graph = g;

	return GMRFLib_SUCCESS;
}
double GMRFLib_offset_Qfunc(int node, int nnode, void *arg)
{
	GMRFLib_offset_arg_tp *a = NULL;

	a = (GMRFLib_offset_arg_tp *) arg;

	if (IMIN(node, nnode) < a->offset || IMAX(node, nnode) >= a->offset + a->n) {
		return 0.0;
	}

	return (*a->Qfunc) (node - a->offset, nnode - a->offset, a->Qfunc_arg);
}
int GMRFLib_offset(GMRFLib_offset_tp ** off, int n_new, int offset, GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg)
{
	/*
	 * make an object containg the new graph, Qfunc and Qfunc_arg when an graph if inserted into a graph. 
	 */

	GMRFLib_offset_tp *val = NULL;
	GMRFLib_offset_arg_tp *offset_arg;

	offset_arg = Calloc(1, GMRFLib_offset_arg_tp);
	offset_arg->Qfunc = Qfunc;
	offset_arg->Qfunc_arg = Qfunc_arg;
	offset_arg->offset = offset;
	offset_arg->n = graph->n;

	val = Calloc(1, GMRFLib_offset_tp);
	GMRFLib_offset_graph(&(val->graph), n_new, offset, graph);

	val->Qfunc = GMRFLib_offset_Qfunc;
	val->Qfunc_arg = (void *) offset_arg;

	*off = val;

	return GMRFLib_SUCCESS;
}
int *GMRFLib_connected_components(GMRFLib_graph_tp * g)
{
	/*
	 * return a vector of length n, indicating which connecting component each node belongs to 
	 */

	if (g == NULL || g->n == 0) {
		return NULL;
	}

	int i, n, *cc, ccc;
	char *visited;

	n = g->n;
	cc = Calloc(n, int);
	ccc = -1;					       /* the counter. yes, start at -1 */
	visited = Calloc(n, char);

	for (i = 0; i < n; i++) {
		if (!visited[i]) {
			ccc++;
			GMRFLib_connected_components_do(i, g, cc, visited, &ccc);
		}
	}

	Free(visited);

	return cc;
}
int GMRFLib_connected_components_do(int node, GMRFLib_graph_tp * g, int *cc, char *visited, int *ccc)
{
	if (visited[node]) {				       /* I don't need this but include it for clarity */
		return GMRFLib_SUCCESS;
	}

	visited[node] = 1;
	cc[node] = *ccc;

	int i, nnode;
	for (i = 0; i < g->nnbs[node]; i++) {
		nnode = g->nbs[node][i];
		if (!visited[nnode]) {			       /* faster to do a check here than to doit inside the funcall */
			GMRFLib_connected_components_do(nnode, g, cc, visited, ccc);
		}
	}

	return GMRFLib_SUCCESS;
}


/*
  Example for manual
 */

/*! \page ex_graph Graph specification and handling. 
  This page includes the examples
  \ref ex_graph_sec1 and \ref ex_graph_sec2 

  \section ex_graph_sec1 Reading and printing graphs

  \par Description:

  This program performs the following tasks:
  
  - Reads a graph from the file \c graph1.dat, given below,
  of the form required by the function \c GMRFLib_read_graph(). 
  The specified graph has 10 nodes, with each of them having between 1 and 5 neighbours. 
  \verbinclude doxygen_graph_2.txt
  
  - Prints the graph specification to standard output.

  - Computes a subgraph, by specifying nodes to be removed from
  the graph. More specifically, the first half of the nodes,
  that is node 0, ..., 4, are specified to be removed.  To
  compute the subgraph, the function \c GMRFLib_compute_subgraph() is called.
  
  - Frees allocated memory.

  \par Program code:

  \verbinclude example-doxygen-graph1.txt

  \par Output:

  \verbinclude doxygen_graph_3.txt

  \section ex_graph_sec2 Creating lattice graphs, linear graphs and folded graphs 

  \par Description:

  This program describes how to create a lattice graph and a linear
  graph, and how to expand the neighbourhood of a graph. The lattice
  graph is defined on a 6 \c x 7 lattice, such that \c nr = 6 and
  \c nc = 7. The neighbourhood is defined to be a 3 \c x
  3-neighbourhood, such that \c mr = \c mc = 1 in (ref 1).  The
  node numbers corresponding to the lattice indices are extracted, and
  the reverse operation is also illustrated.  The linear graph is
  defined by a cyclic AR(1)-model with 10 elements. A third graph is
  generated by expanding the neighbourhood of the linear graph twice.

  \par Program code:

  \verbinclude example-doxygen-graph1.txt

  \par Output:

  \verbinclude doxygen_graph_5.txt

*/
