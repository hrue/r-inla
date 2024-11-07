
/* graph.c
 * 
 * Copyright (C) 2001-2024 Havard Rue
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

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

static int graph_store_use = 1;
static map_strvp graph_store;
static int graph_store_must_init = 1;
static int graph_store_debug = 0;

#define NUM_THREADS_GRAPH(graph_) ((graph_)->n > 1024 ? 2 : 1)

int GMRFLib_graph_init_store(void)
{
	GMRFLib_ENTER_ROUTINE;
	graph_store_debug = GMRFLib_DEBUG_IF_TRUE();

	if (graph_store_use) {
		if (graph_store_must_init) {
			map_strvp_init_hint(&graph_store, 128);
			graph_store.alwaysdefault = 1;
			graph_store_must_init = 0;
			if (graph_store_debug) {
				printf("\tgraph_store: init storage\n");
			}
		}
	}
	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_graph_mk_empty(GMRFLib_graph_tp **graph)
{
	/*
	 * this function creates an empty graph 
	 */

	*graph = Calloc(1, GMRFLib_graph_tp);
	(*graph)->n = 0;
	(*graph)->nbs = NULL;
	(*graph)->lnbs = NULL;
	(*graph)->snbs = NULL;
	(*graph)->nnbs = NULL;
	(*graph)->lnnbs = NULL;
	(*graph)->snnbs = NULL;
	(*graph)->sha = NULL;

	return GMRFLib_SUCCESS;
}

int GMRFLib_graph_read(GMRFLib_graph_tp **graph, const char *filename)
{
	GMRFLib_graph_read_binary(graph, filename);
	if (*graph != NULL) {
		return GMRFLib_SUCCESS;
	} else {
		return GMRFLib_graph_read_ascii(graph, filename);
	}
}

int GMRFLib_graph_read_ascii(GMRFLib_graph_tp **graph, const char *filename)
{
#define TO_INT(_ix, _x) \
	if (1) {							\
		_ix = (int) (_x);					\
		if (!ISEQUAL((double) (_ix), (_x))) {			\
			char *msg = NULL;				\
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
	GMRFLib_graph_mk_empty(graph);
	assert(graph);

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
				Memcpy(&storage[storage_indx], (*graph)->nbs[i], (size_t) (sizeof(int) * (*graph)->nnbs[i]));
				Free((*graph)->nbs[i]);
				(*graph)->nbs[i] = &storage[storage_indx];
				storage_indx += (*graph)->nnbs[i];
			} else {
				(*graph)->nbs[i] = NULL;
			}
		}
	}
	GMRFLib_EWRAP0(GMRFLib_graph_validate(stderr, *graph));

	GMRFLib_graph_prepare(*graph);			       /* prepare the graph for computations */
#undef TO_INT
	return GMRFLib_SUCCESS;
}

int GMRFLib_printf_graph(FILE *fp, GMRFLib_graph_tp *graph)
{
	int i, j;
	FILE *fpp = NULL;

	fpp = (fp ? fp : stdout);

	fprintf(fpp, "graph has n=%1d nodes and nnz=%1d non-zero off-diagonals\n", graph->n, graph->nnz);
	fprintf(fpp, "\tlnnz=%1d snnz=%1d n_ptr=%1d n_idx=%1d\n", graph->lnnz, graph->snnz, graph->n_ptr, graph->n_idx);
	for (i = 0; i < graph->n; i++) {
		fprintf(fpp, "node %1d has %1d neighbors, %1d lneighbors and %1d sneighbors:", i, graph->nnbs[i], graph->lnnbs[i], graph->snnbs[i]);
		for (j = 0; j < graph->nnbs[i]; j++) {
			fprintf(fpp, " %1d", graph->nbs[i][j]);
		}
		fprintf(fpp, "\n");
	}
	for (i = 0; i < graph->n; i++) {
		fprintf(fpp, "node %1d has %1d lneighbors:", i, graph->lnnbs[i]);
		for (j = 0; j < graph->lnnbs[i]; j++) {
			fprintf(fpp, " %1d", graph->lnbs[i][j]);
		}
		fprintf(fpp, "\n");
	}
	for (i = 0; i < graph->n; i++) {
		fprintf(fpp, "node %1d has %1d sneighbors:", i, graph->snnbs[i]);
		for (j = 0; j < graph->snnbs[i]; j++) {
			fprintf(fpp, " %1d", graph->snbs[i][j]);
		}
		fprintf(fpp, "\n");
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_graph_write(const char *filename, GMRFLib_graph_tp *graph)
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

	GMRFLib_graph_write2(fp, graph);
	fclose(fp);

	return GMRFLib_SUCCESS;
}

int GMRFLib_graph_write2(FILE *fp, GMRFLib_graph_tp *graph)
{
	fprintf(fp, "%1d\n", graph->n);
	for (int i = 0; i < graph->n; i++) {
		fprintf(fp, "%1d %1d", i, graph->nnbs[i]);
		for (int j = 0; j < graph->nnbs[i]; j++) {
			fprintf(fp, " %1d", graph->nbs[i][j]);
		}
		fprintf(fp, "\n");
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_graph_write_b(const char *filename, GMRFLib_graph_tp *graph)
{
	/*
	 * use base 1: so the nodes are 1...n, not 0..n-1. this makes the connection to R and R-inla easier. However, the read_graph routines will autodetect if
	 * the nodes are 0..n-1 or 1...n. 
	 */

	int tag = GMRFLib_BINARY_GRAPH_FILE_MAGIC, offset = 1, idx, iidx;
	GMRFLib_io_tp *io = NULL;

	if (!filename || !graph) {
		return GMRFLib_SUCCESS;
	}

	GMRFLib_EWRAP0(GMRFLib_io_open(&io, filename, "wb"));
	GMRFLib_io_write(io, (const void *) &tag, sizeof(int));	/* so we can detect that this is of correct format */
	GMRFLib_io_write(io, (const void *) &(graph->n), sizeof(int));
	for (int i = 0; i < graph->n; i++) {
		idx = i + offset;
		GMRFLib_io_write(io, (const void *) &idx, sizeof(int));
		GMRFLib_io_write(io, (const void *) &(graph->nnbs[i]), sizeof(int));
		if (graph->nnbs[i]) {
			if (offset == 0) {
				GMRFLib_io_write(io, (const void *) (graph->nbs[i]), (unsigned int) (graph->nnbs[i] * sizeof(int)));
			} else {
				for (int j = 0; j < graph->nnbs[i]; j++) {
					iidx = graph->nbs[i][j] + offset;
					GMRFLib_io_write(io, (const void *) &iidx, (unsigned int) sizeof(int));
				}
			}
		}
	}
	GMRFLib_io_close(io);

	return GMRFLib_SUCCESS;
}

int GMRFLib_graph_read_binary(GMRFLib_graph_tp **graph, const char *filename)
{
	int tag = 0;
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

	GMRFLib_graph_mk_empty(&g);
	GMRFLib_EWRAP0(GMRFLib_io_read(io, (void *) &(g->n), sizeof(int)));
	GMRFLib_ASSERT(g->n >= 0, GMRFLib_EPARAMETER);

	g->nnbs = Calloc(g->n + 1, int);		       /* yes. */
	g->nbs = Calloc(g->n + 1, int *);		       /* yes. */

	for (int i = 0; i < g->n; i++) {
		int ii = 0;
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
	for (int i = 0; i < g->n + 1; i++) {
		if (g->nbs[i]) {			       /* yes, we check the ptr if this is read or not */
			min_node = IMIN(min_node, i);
			max_node = IMAX(max_node, i);
			for (int j = 0; j < g->nnbs[i]; j++) {
				int idx = g->nbs[i][j];
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
		for (int i = 0; i < g->n; i++) {
			g->nnbs[i] = g->nnbs[i + 1];
			g->nbs[i] = g->nbs[i + 1];
			g->nbs[i + 1] = NULL;
			for (int j = 0; j < g->nnbs[i]; j++) {
				g->nbs[i][j]--;
			}
		}
	} else if (min_node == 0 && max_node == g->n - 1) {
		/*
		 * ok 
		 */
	} else {
		fprintf(stderr, "\n\nmin_node = %1d max_node = %1d. this should not happen.\n", min_node, max_node);
		GMRFLib_ERROR(GMRFLib_ESNH);
	}

	GMRFLib_graph_prepare(g);
	GMRFLib_graph_duplicate(graph, g);

	for (int i = 0; i < g->n + 1; i++) {		       /* yes, its +1 */
		Free(g->nbs[i]);
	}
	Free(g->nbs);
	Free(g->nnbs);
	Free(g);

	return GMRFLib_SUCCESS;
}

int GMRFLib_graph_free(GMRFLib_graph_tp *graph)
{
	/*
	 * free a graph build with ``GMRFLib_graph_read'' 
	 */
	int i;
	GMRFLib_DEBUG_INIT();

	if (!graph) {
		return GMRFLib_SUCCESS;
	}

	if (graph_store_use && graph->sha) {
		void *p = NULL;
		p = map_strvp_ptr(&graph_store, (char *) graph->sha);
		if (graph_store_debug) {
			if (p) {
				printf("\t[%1d] graph_store: graph is found in store: do not free.\n", omp_get_thread_num());
			} else {
				printf("\t[%1d] graph_store: graph is not found in store: free.\n", omp_get_thread_num());
			}
		}
		if (p) {
			return GMRFLib_SUCCESS;
		}
	}

	for (i = 0; i < graph->n; i++) {
		if (graph->nnbs[i]) {
			Free(graph->nbs[i]);
			break;				       /* new memory layout, only `free' the first!!! */
		}
	}
	Free(graph->nbs);
	Free(graph->nnbs);
	Free(graph->lnbs);
	Free(graph->snbs);
	Free(graph->lnnbs);
	Free(graph->snnbs);
	Free(graph->sha);
	Free(graph->rowptr);
	Free(graph->colptr);
	Free(graph->rowidx);
	Free(graph->colidx);
	Free(graph);

	return GMRFLib_SUCCESS;
}

int GMRFLib_getbit(GMRFLib_uchar c, unsigned int bitno)
{
	/*
	 * return bit-number BITNO, bitno = 0, 1, 2, ..., 7 
	 */
	unsigned int zero = 0;

	return (int) ((c >> bitno) & ~(~zero << 1));
}

int GMRFLib_setbit(GMRFLib_uchar *c, unsigned int bitno)
{
	/*
	 * set bitno 0, 1, 2, ..., 7, to TRUE 
	 */
	unsigned int zero = 0;

	*c = *c | ((~(~zero << (bitno + 1)) & (~zero << bitno)));

	return GMRFLib_SUCCESS;
}

int GMRFLib_printbits(FILE *fp, GMRFLib_uchar c)
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

void *GMRFLib_bsearch2(int key, int n, int *array, int *guess)
{
	int mid, top, val, *piv = NULL, *base = array;
	int low = 0;

	if (array[guess[0]] <= key) {
		low = guess[0];
	}

	base += low;
	mid = top = n - low;

	while (mid) {
		mid = top / 2;
		piv = base + mid;
		val = key - *piv;
		if (val == 0) {
			guess[0] = piv - array;
			return piv;
		}
		if (val > 0) {
			base = piv;
		}
		top -= mid;
	}

	return NULL;
}

void *GMRFLib_bsearch(int key, int n, int *array)
{
	int mid, top, val, *piv = NULL, *base = array;
	mid = top = n;

	while (mid) {
		mid = top / 2;
		piv = base + mid;
		val = key - *piv;
		if (val == 0) {
			return piv;
		}
		if (val > 0) {
			base = piv;
		}
		top -= mid;
	}

	return NULL;
}

int GMRFLib_graph_is_nb(int node, int nnode, GMRFLib_graph_tp *graph)
{
	int imin, imax;
	if (node < nnode) {
		imin = node;
		imax = nnode;
	} else {
		assert(node != nnode);
		imin = nnode;
		imax = node;
	}

	int m = graph->lnnbs[imin];
	if (m) {
		int *nb = graph->lnbs[imin];
		if (nnode <= nb[m - 1]) {
			return (GMRFLib_bsearch(imax, m, nb) != NULL);
		}
	}
	return 0;
}

int GMRFLib_graph_is_nb_g(int node, int nnode, GMRFLib_graph_tp *graph, int *g)
{
	/*
	 * return 1 if nnode is a neighbour of node, otherwise 0. assume that the nodes are sorted. note that if node == nnode,
	 * then they are not neighbours.
	 */

	assert(node < nnode);

	int m = graph->lnnbs[node];
	if (m) {
		int *nb = graph->lnbs[node];
		if (nnode <= nb[m - 1]) {
			return (GMRFLib_bsearch2(nnode, m, nb, g) != NULL);
		}
	}
	return 0;
}

int GMRFLib_graph_add_crs_crc(GMRFLib_graph_tp *graph)
{

	if (!graph) {
		return GMRFLib_SUCCESS;
	}

	int n = graph->n;
	int N = graph->n + graph->nnz / 2;

	// TAUCS
	int *colptr = Calloc(graph->n + 1, int);
	int *rowidx = Calloc(N, int);

	// PARDISO
	int *rowptr = Calloc(graph->n + 1, int);
	int *colidx = Calloc(N, int);

	colptr[0] = 0;
	for (int i = 0, k = 0; i < n; i++) {
		rowidx[k] = i;
		k++;
		Memcpy(&(rowidx[k]), graph->snbs[i], graph->snnbs[i] * sizeof(int));
		k += graph->snnbs[i];
		colptr[i + 1] = colptr[i] + graph->snnbs[i] + 1;
	}

	rowptr[0] = 0;
	for (int i = 0, k = 0; i < n; i++) {
		colidx[k] = i;
		k++;
		Memcpy(&(colidx[k]), graph->lnbs[i], graph->lnnbs[i] * sizeof(int));
		k += graph->lnnbs[i];
		rowptr[i + 1] = rowptr[i] + graph->lnnbs[i] + 1;
	}

	graph->n_ptr = graph->n + 1;
	graph->n_idx = N;
	graph->rowptr = rowptr;
	graph->colidx = colidx;
	graph->colptr = colptr;
	graph->rowidx = rowidx;

	return GMRFLib_SUCCESS;
}


int GMRFLib_graph_add_row2col(GMRFLib_graph_tp *graph)
{
	// mapping between the CSR and CSC; see build_TAUCS.
	// only used for TAUCS

	if (!graph) {
		return GMRFLib_SUCCESS;
	}

	int n = graph->n;
	int N = graph->n + graph->nnz / 2;
	int *row2col = Calloc(N, int);
	assert(N > 0);

	if (0) {
		int *row = Calloc(graph->n + 1, int);
		row[0] = 0;
		for (int i = 1; i < n + 1; i++) {
			row[i] = row[i - 1] + 1 + graph->lnnbs[i - 1];
		}
	}
#define Q(i_, j_, kk_) (graph->rowptr[IMIN(i_, j_)] + kk_)
	for (int i = 0, k = 0; i < n; i++) {
		row2col[k++] = Q(i, i, 0);
		for (int jj = 0; jj < graph->snnbs[i]; jj++) {
			int j = graph->snbs[i][jj];
			int kk = 1 + GMRFLib_iwhich_sorted(i, graph->lnbs[j], graph->lnnbs[j]);
			row2col[k++] = Q(i, j, kk);
		}
	}
	graph->row2col = row2col;

#undef Q
	return GMRFLib_SUCCESS;
}

int GMRFLib_graph_prepare(GMRFLib_graph_tp *graph)
{
	/*
	 * prepare the graph by sort the vertices in increasing orders 
	 */
	graph->nnz = GMRFLib_isum(graph->n, graph->nnbs);

	GMRFLib_graph_sort(graph);			       /* must be before lnbs */
	GMRFLib_graph_add_lnbs_info(graph);		       /* must be before sha */
	// need this check as graph is also used in the non-symmetric case for matrix
	if (graph->lnnz == graph->snnz) {
		GMRFLib_graph_add_crs_crc(graph);
		GMRFLib_graph_add_row2col(graph);	       /* needs to come after crs_crc */
	}
	GMRFLib_graph_add_sha(graph);

	graph->max_nnbs = -1;
	graph->max_nnbs = GMRFLib_graph_max_nnbs(graph);

	return GMRFLib_SUCCESS;
}

int GMRFLib_graph_add_lnbs_info(GMRFLib_graph_tp *graph)
{
	// these nodes are sorted

	if (!graph) {
		return GMRFLib_SUCCESS;
	}

	int n = graph->n;
	graph->lnnbs = Calloc(n, int);
	graph->lnbs = Calloc(n, int *);
	graph->snnbs = Calloc(n, int);
	graph->snbs = Calloc(n, int *);

#define CODE_BLOCK							\
	for (int i = 0; i < n; i++) {					\
		int k = graph->nnbs[i];					\
		for (int jj = 0; jj < graph->nnbs[i]; jj++) {		\
			int j = graph->nbs[i][jj];			\
			if (j > i) {					\
				k = jj;					\
				graph->lnbs[i] = &(graph->nbs[i][jj]);	\
				break;					\
			}						\
		}							\
		graph->lnnbs[i] = graph->nnbs[i] - k;			\
									\
		k = 0;							\
		for (int jj = 0; jj < graph->nnbs[i]; jj++) {		\
			int j = graph->nbs[i][jj];			\
			if (j < i) {					\
				k++;					\
			} else {					\
				break;					\
			}						\
		}							\
		graph->snnbs[i] = k;					\
		graph->snbs[i] = (k > 0 ? graph->nbs[i] : NULL);	\
		assert(graph->snnbs[i] + graph->lnnbs[i] == graph->nnbs[i]); \
	}

	RUN_CODE_BLOCK(NUM_THREADS_GRAPH(graph), 0, 0);
#undef CODE_BLOCK

	graph->lnnz = GMRFLib_isum(graph->n, graph->lnnbs);
	graph->snnz = GMRFLib_isum(graph->n, graph->snnbs);

	return GMRFLib_SUCCESS;
}

int GMRFLib_graph_mk_unique(GMRFLib_graph_tp *graph)
{
	/*
	 * ensure the neigbours are unique. the neigbours must be sorted. 
	 */

	if (!graph) {
		return GMRFLib_SUCCESS;
	}
#define CODE_BLOCK							\
	for (int i = 0; i < graph->n; i++) {				\
		if (graph->nnbs[i]) {					\
			int k = 0;					\
			for (int j = 1; j < graph->nnbs[i]; j++) {	\
				if (graph->nbs[i][k] != graph->nbs[i][j]) { \
					graph->nbs[i][++k] = graph->nbs[i][j]; \
				}					\
			}						\
			graph->nnbs[i] = k + 1;				\
		}							\
	}

	RUN_CODE_BLOCK(NUM_THREADS_GRAPH(graph), 0, 0);
#undef CODE_BLOCK


	return GMRFLib_SUCCESS;
}

int GMRFLib_graph_sort(GMRFLib_graph_tp *graph)
{
	/*
	 * sort the vertices in increasing order 
	 */

	if (!graph) {
		return GMRFLib_SUCCESS;
	}
#define CODE_BLOCK							\
	for (int i = 0; i < graph->n; i++) {				\
		if (graph->nnbs[i]) {					\
			GMRFLib_qsort(graph->nbs[i], (size_t) graph->nnbs[i], sizeof(int), GMRFLib_icmp); \
		}							\
	}

	RUN_CODE_BLOCK(NUM_THREADS_GRAPH(graph), 0, 0);
#undef CODE_BLOCK

	return GMRFLib_SUCCESS;
}

int GMRFLib_graph_comp_bw(int *bandwidth, GMRFLib_graph_tp *graph, int *remap)
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

int GMRFLib_graph_validate(FILE *fp, GMRFLib_graph_tp *graph)
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

int GMRFLib_graph_remap(GMRFLib_graph_tp **ngraph, GMRFLib_graph_tp *graph, int *remap)
{
	/*
	 * return the remapped graph based on 'graph'. the returned graph has the identity mapping. 
	 */

	int i, j, k, nnb, indx, *hold = NULL;

	if (!graph) {
		*ngraph = (GMRFLib_graph_tp *) NULL;
		return GMRFLib_SUCCESS;
	}

	GMRFLib_graph_mk_empty(ngraph);
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
	nnb = GMRFLib_isum((*ngraph)->n, (*ngraph)->nnbs);
	hold = Calloc(IMAX(1, nnb), int);

	for (i = 0, indx = 0; i < (*ngraph)->n; i++) {
		if ((*ngraph)->nnbs[i]) {
			Memcpy(&hold[indx], (*ngraph)->nbs[i], (size_t) ((*ngraph)->nnbs[i] * sizeof(int)));
			Free((*ngraph)->nbs[i]);
			(*ngraph)->nbs[i] = &hold[indx];
		} else {
			Free((*ngraph)->nbs[i]);
		}
		indx += (*ngraph)->nnbs[i];
	}
	GMRFLib_graph_prepare(*ngraph);

	// doit like this, to prevent compiler warings
	if (nnb == 0) {
		Free(hold);
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_graph_duplicate(GMRFLib_graph_tp **graph_new, GMRFLib_graph_tp *graph_old)
{
	int m, i, n, *hold = NULL, hold_idx;
	GMRFLib_graph_tp *g = NULL;

	GMRFLib_ENTER_ROUTINE;
	if (!graph_old) {
		*graph_new = NULL;
		GMRFLib_LEAVE_ROUTINE;
		return GMRFLib_SUCCESS;
	}

	if (graph_store_use && graph_old->sha) {
		void **p = NULL;
		p = map_strvp_ptr(&graph_store, (char *) graph_old->sha);
		if (graph_store_debug) {
			if (p) {
				printf("\t[%1d] graph_store: graph is found in store: do not duplicate.\n", omp_get_thread_num());
			} else {
				printf("\t[%1d] graph_store: graph is not found in store: duplicate.\n", omp_get_thread_num());
			}
		}
		if (p) {
			*graph_new = (GMRFLib_graph_tp *) * p;
			GMRFLib_LEAVE_ROUTINE;
			return GMRFLib_SUCCESS;
		}
	}

	GMRFLib_graph_mk_empty(&g);
	g->n = n = graph_old->n;
	g->nnbs = Calloc(n, int);
	Memcpy(g->nnbs, graph_old->nnbs, (size_t) (n * sizeof(int)));

	m = graph_old->nnz;
	hold = Calloc(IMAX(1, m), int);
	g->nbs = Calloc(n, int *);

	for (i = hold_idx = 0; i < n; i++) {
		if (g->nnbs[i]) {
			g->nbs[i] = &hold[hold_idx];
			Memcpy(g->nbs[i], graph_old->nbs[i], (size_t) (g->nnbs[i] * sizeof(int)));
			hold_idx += g->nnbs[i];
		}
	}

	*graph_new = g;
	GMRFLib_graph_prepare(g);

	if (graph_store_use && graph_old->sha) {
		if (graph_store_debug) {
			printf("\t[%1d] graph_store: store graph 0x%p\n", omp_get_thread_num(), (void *) g);
		}
#pragma omp critical (Name_c524502943d363cb45e15d587b32804a133415b2)
		{
			map_strvp_set(&graph_store, (char *) g->sha, (void *) g);
		}
	}

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

size_t GMRFLib_graph_sizeof(GMRFLib_graph_tp *graph)
{
	/*
	 * return, approximately, the sizeof GRAPH 
	 */

	if (!graph) {
		return 0;
	}

	size_t siz = 0;
	int m, n;

	n = graph->n;
	m = GMRFLib_isum(n, graph->nnbs);
	siz += sizeof(int) + m * sizeof(int) + 2 * n * sizeof(int) + 2 * n * sizeof(int *);

	return siz;
}

int GMRFLib_graph_comp_subgraph(GMRFLib_graph_tp **subgraph, GMRFLib_graph_tp *graph, char *remove_flag, int **node_map)
{
	if (!remove_flag) {
		if (node_map) {
			assert(graph->n > 0);
			*node_map = Calloc(graph->n, int);
#pragma omp simd
			for (int i = 0; i < graph->n; i++) {
				(*node_map)[i] = i;
			}
		}
		return GMRFLib_graph_duplicate(subgraph, graph);
	} else {

		/*
		 * return a subgraph of graph, by ruling out those nodes for which remove_flag[i] is true, keeping those which
		 * remove_flag[i] false. 
		 */
		int nneig, nn = 0, n_neig_tot = 0, storage_indx, *nmap = NULL, *sg_iidx = NULL, *storage = NULL, free_remove_flag = 0;

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
			if (node_map) {
				*node_map = nmap;
			}
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

		GMRFLib_graph_mk_empty(subgraph);

#pragma omp simd reduction(+: nn)
		for (int i = 0; i < graph->n; i++) {
			nn += (!remove_flag[i]);
		}
		(*subgraph)->n = nn;

		if (!((*subgraph)->n)) {
			GMRFLib_LEAVE_ROUTINE;
			if (node_map) {
				*node_map = nmap;
			}
			return GMRFLib_SUCCESS;
		}

		/*
		 * create space 
		 */
		(*subgraph)->nnbs = Calloc((*subgraph)->n, int);
		(*subgraph)->nbs = Calloc((*subgraph)->n, int *);

		/*
		 * make the mapping of nodes 
		 */
		sg_iidx = Calloc(graph->n, int);
		nmap = Calloc(nn, int);
		int k = 0;
		for (int i = 0; i < graph->n; i++) {
			if (!remove_flag[i]) {
				sg_iidx[i] = k;
				nmap[k] = i;
				k++;
			} else {
				sg_iidx[i] = -1;	       /* to force a failure if used wrong */
			}
		}
		assert(k == nn);
		if (node_map) {
			*node_map = nmap;
		}

		/*
		 * parse the graph and collect nodes not to be removed. 
		 */
		k = 0;
		n_neig_tot = 0;
		for (int i = 0; i < graph->n; i++) {
			if (!remove_flag[i]) {
				nneig = 0;
#pragma omp simd reduction(+: nneig)
				for (int j = 0; j < graph->nnbs[i]; j++) {
					nneig += (!remove_flag[graph->nbs[i][j]]);
				}
				n_neig_tot += nneig;
				(*subgraph)->nnbs[k] = nneig;

				if (nneig > 0) {
					(*subgraph)->nbs[k] = Calloc(nneig, int);
					nneig = 0;
					for (int j = 0; j < graph->nnbs[i]; j++) {
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

		/*
		 * map the graph to a more computational convenient memory layout! use just one long vector to store all the neighbours. 
		 */
		if (n_neig_tot) {
			storage = Calloc(n_neig_tot, int);
			storage_indx = 0;
			for (int i = 0; i < (*subgraph)->n; i++) {
				if ((*subgraph)->nnbs[i]) {
					Memcpy(&storage[storage_indx], (*subgraph)->nbs[i], (size_t) (sizeof(int) * (*subgraph)->nnbs[i]));
					Free((*subgraph)->nbs[i]);
					(*subgraph)->nbs[i] = &storage[storage_indx];
					storage_indx += (*subgraph)->nnbs[i];
				} else {
					(*subgraph)->nbs[i] = NULL;
				}
			}
		}
		GMRFLib_graph_prepare(*subgraph);

		if (free_remove_flag) {
			Free(remove_flag);		       /* if we have used our own */
		}
		Free(sg_iidx);

		GMRFLib_LEAVE_ROUTINE;
		return GMRFLib_SUCCESS;
	}
}

int GMRFLib_convert_to_mapped(double *destination, double *source, GMRFLib_graph_tp *graph, int *remap)
{
	/*
	 * convert from the real-world to the mapped world. source might be NULL. 
	 */
	if ((destination && source) && !OVERLAP(destination, source, graph->n)) {
		// for (int i = 0; i < graph->n; i++) destination[remap[i]] = source[i];
		GMRFLib_unpack(graph->n, source, destination, remap);
	} else {
		static double **wwork = NULL;
		static int *wwork_len = NULL;
		double *work = NULL;
		if (!wwork) {
#pragma omp critical (Name_3fd48b64af36f6b6464b44ae9bf0fc94062b77fc)
			{
				if (!wwork) {
					wwork_len = Calloc(GMRFLib_CACHE_LEN(), int);
					wwork = Calloc(GMRFLib_CACHE_LEN(), double *);
				}
			}
		}

		int cache_idx = 0;
		GMRFLib_CACHE_SET_ID(cache_idx);
		if (graph->n > wwork_len[cache_idx]) {
			Free(wwork[cache_idx]);
			wwork_len[cache_idx] = graph->n;
			wwork[cache_idx] = Calloc(wwork_len[cache_idx], double);
		}
		work = wwork[cache_idx];
		Memcpy(work, destination, graph->n * sizeof(double));
		// for (int i = 0; i < graph->n; i++) destination[remap[i]] = work[i];
		GMRFLib_unpack(graph->n, work, destination, remap);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_convert_from_mapped(double *destination, double *source, GMRFLib_graph_tp *graph, int *remap)
{
	/*
	 * convert from the mapped-world to the real world. source might be NULL. 
	 */
	if ((destination && source) && !OVERLAP(destination, source, graph->n)) {
		// for (int i = 0; i < graph->n; i++) destination[i] = source[remap[i]];
		GMRFLib_pack(graph->n, source, remap, destination);
	} else {
		static double **wwork = NULL;
		static int *wwork_len = NULL;
		double *work = NULL;
		if (!wwork) {
#pragma omp critical (Name_1a7225070120a13086515e5c250e27a43a27bdd8)
			{
				if (!wwork) {
					wwork_len = Calloc(GMRFLib_CACHE_LEN(), int);
					wwork = Calloc(GMRFLib_CACHE_LEN(), double *);
				}
			}
		}

		int cache_idx = 0;
		GMRFLib_CACHE_SET_ID(cache_idx);
		if (graph->n > wwork_len[cache_idx]) {
			Free(wwork[cache_idx]);
			wwork_len[cache_idx] = graph->n;
			wwork[cache_idx] = Calloc(wwork_len[cache_idx], double);
		}
		work = wwork[cache_idx];

		Memcpy(work, destination, graph->n * sizeof(double));
		// for (int i = 0; i < graph->n; i++) destination[i] = work[remap[i]];
		GMRFLib_pack(graph->n, work, remap, destination);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_graph_max_nnbs(GMRFLib_graph_tp *graph)
{
	return (graph->max_nnbs >= 0 ? graph->max_nnbs : GMRFLib_imax_value(graph->nnbs, graph->n, NULL));
}

int GMRFLib_graph_max_lnnbs(GMRFLib_graph_tp *graph)
{
	return (GMRFLib_imax_value(graph->lnnbs, graph->n, NULL));
}

int GMRFLib_graph_max_snnbs(GMRFLib_graph_tp *graph)
{
	return (GMRFLib_imax_value(graph->snnbs, graph->n, NULL));
}

int GMRFLib_Qx(int thread_id, double *result, double *x, GMRFLib_graph_tp *graph, GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_arg)
{
	return (GMRFLib_Qx2(thread_id, result, x, graph, Qfunc, Qfunc_arg, NULL));
}

int GMRFLib_Qx2(int thread_id, double *result, double *x, GMRFLib_graph_tp *graph, GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_arg, double *diag)
{
	GMRFLib_ENTER_ROUTINE;

	const int debug = 0;
	int run_parallel = (GMRFLib_Qx_strategy != 0);
	int max_t;
	double *values = NULL, res;

	max_t = IMAX(GMRFLib_openmp->max_threads_inner, GMRFLib_openmp->max_threads_outer);
	assert(result);
	GMRFLib_fill(graph->n, 0.0, result);

	int m = GMRFLib_align(1 + GMRFLib_graph_max_nnbs(graph), sizeof(double));
	Calloc_init(m + (!diag ? graph->n : 0), 2);
	values = Calloc_get(m);
	assert(values);
	if (!diag) {
		diag = Calloc_get(graph->n);
	}
	res = Qfunc(thread_id, 0, -1, values, Qfunc_arg);

	if (ISNAN(res)) {
		if (run_parallel) {
			if (debug) {
				FIXME("Qx2: run parallel");
			}
#define CODE_BLOCK							\
			for (int i = 0; i < graph->n; i++) {		\
				double sum = (Qfunc(thread_id, i, i, NULL, Qfunc_arg) + diag[i]) * x[i]; \
				int *j_a = graph->nbs[i];		\
				for (int jj = 0; jj < graph->nnbs[i]; jj++) { \
					int j = j_a[jj];		\
					sum += Qfunc(thread_id, i, j, NULL, Qfunc_arg) * x[j]; \
				}					\
				result[i] += sum;			\
			}

			RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 0, 0);
#undef CODE_BLOCK
		} else {
			if (debug) {
				FIXME("Qx2: run serial");
			}
			for (int i = 0; i < graph->n; i++) {
				double sum = 0.0, xi = x[i];
				int *j_a = graph->lnbs[i];

				result[i] += (Qfunc(thread_id, i, i, NULL, Qfunc_arg) + diag[i]) * xi;
				for (int jj = 0; jj < graph->lnnbs[i]; jj++) {
					int j = j_a[jj];
					double qij = Qfunc(thread_id, i, j, NULL, Qfunc_arg);
					sum += qij * x[j];
					result[j] += qij * xi;
				}
				result[i] += sum;
			}
		}
	} else {
		if (run_parallel) {
			if (debug) {
				FIXME("Qx2: run block parallel");
			}
			int n1 = GMRFLib_align(graph->n, sizeof(int));
			double *local_result = Calloc(max_t * n1, double);
			char *used = Calloc(max_t, char);
#define CODE_BLOCK							\
			for (int i = 0; i < graph->n; i++) {		\
				double *r = NULL, *local_values = NULL, xi = x[i]; \
				int *j_a = graph->lnbs[i];		\
				/* may run in serial */			\
				int tnum = (nt__ > 1 ? omp_get_thread_num() : 0); \
				used[tnum] = 1;				\
				r = local_result + tnum * n1;		\
				local_values = CODE_BLOCK_WORK_PTR_x(0, tnum);	\
				Qfunc(thread_id, i, -1, local_values, Qfunc_arg); \
				double sum = (local_values[0] + diag[i]) * xi; \
				for (int jj = 0; jj < graph->lnnbs[i]; jj++) { \
					int j = j_a[jj];		\
					double lval = local_values[jj+1]; \
					sum += lval * x[j];		\
					r[j] += lval * xi;		\
				}					\
				r[i] += sum;				\
			}

			RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 1, m);
#undef CODE_BLOCK

			for (int j = 0; j < max_t; j++) {
				if (used[j]) {
					int offset = j * n1;
					double *r = local_result + offset;
					GMRFLib_daddto(graph->n, r, result);
				}
			}

			Free(used);
			Free(local_result);
		} else {
			if (debug) {
				FIXME("Qx2: run block serial");
			}
			for (int i = 0; i < graph->n; i++) {
				double sum = 0.0, xi = x[i];
				int *j_a = graph->lnbs[i];

				res = Qfunc(thread_id, i, -1, values, Qfunc_arg);
				result[i] += (values[0] + diag[i]) * xi;
				for (int jj = 0; jj < graph->lnnbs[i]; jj++) {
					double val = values[jj + 1];
					int j = j_a[jj];

					sum += val * x[j];
					result[j] += val * xi;
				}
				result[i] += sum;
			}
		}
	}

	Calloc_free();
	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_get_Qrow(int thread_id, int row, int *nelm, int *idx, double *vals, GMRFLib_graph_tp *graph, GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_arg)
{
	// return the row of Q in terms of (idx,valus) with length nelm
	// idx and vals must pre-exists with enough storage

	int ii = 0;
	for (int j = 0; j < graph->snnbs[row]; j++) {
		int col = graph->snbs[row][j];
		vals[ii++] = Qfunc(thread_id, row, col, NULL, Qfunc_arg);
	}

	double res = Qfunc(thread_id, row, -1, vals + ii, Qfunc_arg);
	if (ISNAN(res)) {
		vals[ii++] = Qfunc(thread_id, row, row, NULL, Qfunc_arg);
		for (int j = 0; j < graph->lnnbs[row]; j++) {
			int col = graph->lnbs[row][j];
			vals[ii++] = Qfunc(thread_id, row, col, NULL, Qfunc_arg);
		}
	}

	Memcpy(idx, graph->snbs[row], graph->snnbs[row] * sizeof(int));
	Memcpy(idx + graph->snnbs[row], &row, sizeof(int));
	Memcpy(idx + graph->snnbs[row] + 1, graph->lnbs[row], graph->lnnbs[row] * sizeof(int));
	*nelm = graph->nnbs[row] + 1;

	return GMRFLib_SUCCESS;
}

int GMRFLib_QM(int thread_id, gsl_matrix *result, gsl_matrix *x, GMRFLib_graph_tp *graph, GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_arg)
{
	GMRFLib_ENTER_ROUTINE;

	// taken from GMRFLibP.h
	int nt = ((GMRFLib_OPENMP_IN_PARALLEL_ONE_THREAD() || GMRFLib_OPENMP_IN_SERIAL())?
		  IMAX(GMRFLib_openmp->max_threads_inner, GMRFLib_openmp->max_threads_outer) : GMRFLib_openmp->max_threads_inner);

	int ncol = result->size2;
	int len = GMRFLib_align(1 + GMRFLib_graph_max_nnbs(graph), sizeof(double));
	double *values = Calloc(len * nt, double);
	double *dot_values = Calloc(ncol * nt, double);
	int *indx = Calloc(len * nt, int);

	gsl_matrix_set_zero(result);
	if (nt > 1) {
#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_inner)
		for (int i = 0; i < graph->n; i++) {
			int tnum = omp_get_thread_num();
			double *val = values + tnum * len;
			double *dval = dot_values + tnum * ncol;
			int *id = indx + tnum * len;

			int nelm = 0;
			GMRFLib_get_Qrow(thread_id, i, &nelm, id, val, graph, Qfunc, Qfunc_arg);
#pragma omp simd
			for (int kk = 0; kk < nelm; kk++) {
				id[kk] *= x->tda;
			}
			GMRFLib_fill(ncol, 0.0, dval);
			for (int kk = 0; kk < nelm; kk++) {
				double v = val[kk];
				double *pp = x->data + id[kk];
				GMRFLib_daxpy(ncol, v, pp, dval);
			}
			for (int k = 0; k < ncol; k++) {
				gsl_matrix_set(result, i, k, dval[k]);
			}
		}
	} else {
		double *p1 = NULL, *p2 = NULL, *p3 = NULL, *p4 = NULL;
		for (int i = 0; i < graph->n; i++) {
			double res = Qfunc(thread_id, i, -1, values, Qfunc_arg);
			if (ISNAN(res)) {
				int ii = 0;
				values[ii++] = Qfunc(thread_id, i, i, NULL, Qfunc_arg);
				for (int jj = 0; jj < graph->lnnbs[i]; jj++) {
					int j = graph->lnbs[i][jj];
					values[ii++] = Qfunc(thread_id, i, j, NULL, Qfunc_arg);
				}
			}
			p1 = gsl_matrix_ptr(result, i, 0);
			p3 = gsl_matrix_ptr(x, i, 0);
			// for (int k = 0; k < ncol; k++) p1[k] += p3[k] * values[0];
			GMRFLib_daxpy(ncol, values[0], p3, p1);
			int *j_a = graph->lnbs[i];
			for (int jj = 0; jj < graph->lnnbs[i]; jj++) {
				int j = j_a[jj];
				double qij = values[1 + jj];
				p2 = gsl_matrix_ptr(result, j, 0);
				p4 = gsl_matrix_ptr(x, j, 0);
				// for (int k = 0; k < ncol; k++) {
				// p1[k] += qij * p4[k];
				// p2[k] += qij * p3[k];
				GMRFLib_daxpy(ncol, qij, p4, p1);
				GMRFLib_daxpy(ncol, qij, p3, p2);
			}
		}
	}

	Free(values);
	Free(indx);
	Free(dot_values);

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_QM_ORIG(int thread_id, gsl_matrix *result, gsl_matrix *x, GMRFLib_graph_tp *graph, GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_arg)
{
	GMRFLib_ENTER_ROUTINE;

#define ADDTO(M_, i_, j_, val_) gsl_matrix_set(M_, i_, j_, gsl_matrix_get(M_, i_, j_) + val_)

	int ncol = result->size2;
	double res, *values = NULL;

	int len_values = GMRFLib_align(1 + GMRFLib_graph_max_nnbs(graph), sizeof(double));
	if (GMRFLib_OPENMP_IN_PARALLEL() && GMRFLib_openmp->max_threads_inner > 1) {
		values = Calloc(len_values * GMRFLib_openmp->max_threads_inner, double);
	} else {
		values = Calloc(len_values, double);
	}

	gsl_matrix_set_zero(result);
	res = Qfunc(thread_id, 0, -1, values, Qfunc_arg);
	if (ISNAN(res)) {
		if (0 &&				       // TURN OFF THIS AS THE SERIAL IS JUST SO MUCH BETTER for the moment
		    GMRFLib_OPENMP_IN_PARALLEL() && GMRFLib_openmp->max_threads_inner > 1) {
#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_inner)
			for (int k = 0; k < ncol; k++) {
				for (int i = 0; i < graph->n; i++) {
					double qij = Qfunc(thread_id, i, i, NULL, Qfunc_arg);
					ADDTO(result, i, k, gsl_matrix_get(x, i, k) * qij);
					int *j_a = graph->lnbs[i];
					for (int jj = 0; jj < graph->lnnbs[i]; jj++) {
						int j = j_a[jj];
						qij = Qfunc(thread_id, i, j, NULL, Qfunc_arg);
						ADDTO(result, i, k, qij * gsl_matrix_get(x, j, k));
						ADDTO(result, j, k, qij * gsl_matrix_get(x, i, k));
					}
				}
			}
		} else {
			double *p1 = NULL, *p2 = NULL, *p3 = NULL, *p4 = NULL, qij;
			for (int i = 0; i < graph->n; i++) {
				qij = Qfunc(thread_id, i, i, NULL, Qfunc_arg);
				p1 = gsl_matrix_ptr(result, i, 0);
				p3 = gsl_matrix_ptr(x, i, 0);
				// for (int k = 0; k < ncol; k++) p1[k] += p3[k] * qij;
				GMRFLib_daxpy(ncol, qij, p3, p1);
				int *j_a = graph->lnbs[i];
				for (int jj = 0; jj < graph->lnnbs[i]; jj++) {
					int j = j_a[jj];
					qij = Qfunc(thread_id, i, j, NULL, Qfunc_arg);
					p2 = gsl_matrix_ptr(result, j, 0);
					p4 = gsl_matrix_ptr(x, j, 0);
					// for (int k = 0; k < ncol; k++) {
					// p1[k] += qij * p4[k];
					// p2[k] += qij * p3[k];
					GMRFLib_daxpy(ncol, qij, p4, p1);
					GMRFLib_daxpy(ncol, qij, p3, p2);
				}
			}
		}
	} else {
		if (0 &&				       // TURN OFF THIS AS THE SERIAL IS JUST SO MUCH BETTER for the moment
		    GMRFLib_OPENMP_IN_PARALLEL() && GMRFLib_openmp->max_threads_inner > 1) {
			// I think is less good as it index the matrices in the wrong direction
#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_inner)
			for (int k = 0; k < ncol; k++) {
				double *val = values + omp_get_thread_num() * len_values;
				assert(omp_get_thread_num() < GMRFLib_openmp->max_threads_inner);
				for (int i = 0; i < graph->n; i++) {
					Qfunc(thread_id, i, -1, val, Qfunc_arg);
					ADDTO(result, i, k, gsl_matrix_get(x, i, k) * val[0]);
					int *j_a = graph->lnbs[i];
					for (int jj = 0; jj < graph->lnnbs[i]; jj++) {
						int j = j_a[jj];
						double qij = val[1 + jj];
						ADDTO(result, i, k, qij * gsl_matrix_get(x, j, k));
						ADDTO(result, j, k, qij * gsl_matrix_get(x, i, k));
					}
				}
			}
		} else {
			// better one, as the 'k' loop is sequential with inc=1
			double *p1 = NULL, *p2 = NULL, *p3 = NULL, *p4 = NULL;
			for (int i = 0; i < graph->n; i++) {
				Qfunc(thread_id, i, -1, values, Qfunc_arg);
				p1 = gsl_matrix_ptr(result, i, 0);
				p3 = gsl_matrix_ptr(x, i, 0);
				// for (int k = 0; k < ncol; k++) p1[k] += p3[k] * values[0];
				GMRFLib_daxpy(ncol, values[0], p3, p1);
				int *j_a = graph->lnbs[i];
				for (int jj = 0; jj < graph->lnnbs[i]; jj++) {
					int j = j_a[jj];
					double qij = values[1 + jj];
					p2 = gsl_matrix_ptr(result, j, 0);
					p4 = gsl_matrix_ptr(x, j, 0);
					// for (int k = 0; k < ncol; k++) {
					// p1[k] += qij * p4[k];
					// p2[k] += qij * p3[k];
					GMRFLib_daxpy(ncol, qij, p4, p1);
					GMRFLib_daxpy(ncol, qij, p3, p2);
				}
			}

			if (0) {
				assert(0 == 1);
				// the old code which is more readable.
				for (int i = 0; i < graph->n; i++) {
					Qfunc(thread_id, i, -1, values, Qfunc_arg);
#pragma GCC ivdep
					for (int k = 0; k < ncol; k++) {
						ADDTO(result, i, k, gsl_matrix_get(x, i, k) * values[0]);
					}
					int *j_a = graph->lnbs[i];
					for (int jj = 0; jj < graph->lnnbs[i]; jj++) {
						int j = j_a[jj];
						double qij = values[1 + jj];
#pragma GCC ivdep
						for (int k = 0; k < ncol; k++) {
							ADDTO(result, i, k, qij * gsl_matrix_get(x, j, k));
							ADDTO(result, j, k, qij * gsl_matrix_get(x, i, k));
						}
					}
				}
			}

		}
	}

	Free(values);
#undef ADDTO

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_printf_Qfunc2(int thread_id, FILE *fp, GMRFLib_graph_tp *graph, GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_arg)
{
	// print in sparse matrix style. only for small graphs...
	double value;

	for (int i = 0; i < graph->n; i++) {
		for (int j = 0; j < graph->n; j++) {
			if (i == j || GMRFLib_graph_is_nb(i, j, graph)) {
				value = Qfunc(thread_id, i, j, NULL, Qfunc_arg);
			} else {
				value = 0.0;
			}

			if (ISZERO(value)) {
				fprintf(fp, "  .   ");
			} else {
				fprintf(fp, "%5.2f ", value);
			}
		}
		fprintf(fp, "\n");
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_printf_Qfunc(int thread_id, FILE *fp, GMRFLib_graph_tp *graph, GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_arg)
{
	if (!fp) {
		fp = stdout;
	}
	for (int i = 0; i < graph->n; i++) {
		fprintf(fp, "Q[ %1d , %1d ] = %.10f\n", i, i, Qfunc(thread_id, i, i, NULL, Qfunc_arg));
		for (int j = 0; j < graph->nnbs[i]; j++) {
			int jj = graph->nbs[i][j];
			fprintf(fp, "\tQ[ %1d , %1d ] = %.10f\n", i, jj, Qfunc(thread_id, i, jj, NULL, Qfunc_arg));
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_xQx(int thread_id, double *result, double *x, GMRFLib_graph_tp *graph, GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_arg)
{
	return (GMRFLib_xQx2(thread_id, result, x, graph, Qfunc, Qfunc_arg, NULL));
}

int GMRFLib_xQx2(int thread_id, double *result, double *x, GMRFLib_graph_tp *graph, GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_arg, double *diag)
{
	double *y = NULL;
	Calloc_init(graph->n, 1);
	y = Calloc_get(graph->n);

	GMRFLib_Qx2(thread_id, y, x, graph, Qfunc, Qfunc_arg, diag);
	*result = GMRFLib_ddot(graph->n, x, y);

	Calloc_free();
	return GMRFLib_SUCCESS;
}

int GMRFLib_graph_mk_lattice(GMRFLib_graph_tp **graph, int nrow, int ncol, int nb_row, int nb_col, int cyclic_flag)
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

	GMRFLib_graph_mk_empty(graph);
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

	GMRFLib_graph_prepare(*graph);

	return GMRFLib_SUCCESS;
}

int GMRFLib_lattice2node(int *node, int irow, int icol, int nrow, int UNUSED(ncol))
{
	// *node = icol + irow * ncol;
	*node = irow + icol * nrow;
	return GMRFLib_SUCCESS;
}

int GMRFLib_node2lattice(int node, int *irow, int *icol, int nrow, int UNUSED(ncol))
{
	// *irow = node / ncol;
	// *icol = node - (*irow) * ncol;

	*icol = node / nrow;
	*irow = node - (*icol) * nrow;

	return GMRFLib_SUCCESS;
}

int GMRFLib_graph_mk_linear(GMRFLib_graph_tp **graph, int n, int bw, int cyclic_flag)
{
	int *hold = NULL;

	assert(n > 0);
	bw = IMIN(n - 1, IMAX(0, bw));
	bw = IMAX(0, bw);
	GMRFLib_graph_mk_empty(graph);
	(*graph)->n = n;
	(*graph)->nnbs = Calloc(n, int);
	(*graph)->nbs = Calloc(n, int *);

	if (bw > 0) {
		hold = Calloc(n * 2 * bw, int);		       /* use a linear storage */
		for (int i = 0; i < n; i++) {
			(*graph)->nbs[i] = &hold[i * 2 * bw];  /* set pointers to it */
		}

		if (cyclic_flag) {
			for (int i = 0; i < n; i++) {
				(*graph)->nnbs[i] = 2 * bw;
				int k = 0;
				for (int j = i - bw; j <= i + bw; j++) {
					if (j != i) {
						(*graph)->nbs[i][k++] = MOD(j, n);
					}
				}
			}
		} else {
			for (int i = 0; i < n; i++) {
				(*graph)->nnbs[i] = (i - IMAX(i - bw, 0)) + (IMIN(n - 1, i + bw) - i);
				int k = 0;
				for (int j = i - bw; j <= i + bw; j++) {
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
		for (int i = 0; i < n; i++) {
			(*graph)->nbs[i] = NULL;	       /* if zero neighbours, the pointer should be zero */
		}
	}

	GMRFLib_graph_prepare(*graph);
	return GMRFLib_SUCCESS;
}

int GMRFLib_graph_fold(GMRFLib_graph_tp **ng, GMRFLib_graph_tp *g, GMRFLib_graph_tp *gg)
{
	/*
	 * return ng = g * gg, meaning that ng is g expanded by gg 
	 */

	if (!g && !gg) {
		*ng = NULL;
		return GMRFLib_SUCCESS;
	}
	if ((g && !gg) || (gg && !g)) {
		GMRFLib_graph_duplicate(ng, (g ? g : gg));
		return GMRFLib_SUCCESS;
	}

	GMRFLib_ASSERT(g->n == gg->n, GMRFLib_EPARAMETER);     /* of same size? */

	GMRFLib_graph_tp *newg = NULL;
	GMRFLib_graph_mk_empty(&newg);
	newg->n = g->n;
	newg->nnbs = Calloc(g->n, int);
	newg->nbs = Calloc(g->n, int *);
	assert(newg->nnbs && newg->nbs);

	for (int i = 0; i < newg->n; i++) {
		/*
		 * count number of neighbours neighbours 
		 */
		int nneig = g->nnbs[i];
		for (int j = 0; j < g->nnbs[i]; j++) {
			nneig += gg->nnbs[g->nbs[i][j]];
		}
		if (nneig) {
			newg->nbs[i] = Calloc(nneig, int);
			assert(newg->nbs[i]);

			newg->nnbs[i] = 0;
			int k = 0;
			for (int j = 0; j < g->nnbs[i]; j++) {
				int jj = newg->nbs[i][k++] = g->nbs[i][j];
				for (int ii = 0; ii < gg->nnbs[jj]; ii++) {
					int kk = gg->nbs[jj][ii];
					if (kk != i) {
						newg->nbs[i][k++] = kk;
					}
				}
			}
			newg->nnbs[i] = k;

			/*
			 * make them unique 
			 */
			QSORT_FUN((void *) newg->nbs[i], (size_t) newg->nnbs[i], sizeof(int), GMRFLib_icmp);
			k = 0;
			for (int j = 0; j < newg->nnbs[i]; j++) {
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
	int nnb = GMRFLib_isum(newg->n, newg->nnbs);
	int *hold = Calloc(IMAX(1, nnb), int);
	assert(hold);

	for (int i = 0, indx = 0; i < newg->n; i++) {
		if (newg->nnbs[i]) {
			assert(newg->nbs[i]);
			Memcpy(&hold[indx], newg->nbs[i], newg->nnbs[i] * sizeof(int));
			Free(newg->nbs[i]);
			newg->nbs[i] = &hold[indx];
		} else {
			Free(newg->nbs[i]);
		}
		indx += newg->nnbs[i];
	}
	GMRFLib_graph_prepare(newg);
	*ng = newg;
	// doit like this to prevent compiler warnings
	if (nnb == 0) {
		Free(hold);
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_graph_nfold(GMRFLib_graph_tp **ng, GMRFLib_graph_tp *og, int nfold)
{
	/*
	 * make new graph, 'ng', that is 'nfold' of 'og'
	 * 
	 * nfold = 0: ng = I nfold = 1: ng = og nfold = 2: ng = og*og nfold = 3: ng = og*og*og etc
	 * 
	 */
	GMRFLib_graph_tp *newg = NULL, *oldg = NULL;

	if (nfold == 0) {
		GMRFLib_graph_mk_empty(&newg);
		newg->n = og->n;
		newg->nnbs = Calloc(newg->n, int);
		newg->nbs = Calloc(newg->n, int *);

		GMRFLib_graph_prepare(newg);
	} else if (nfold == 1) {
		GMRFLib_graph_duplicate(&newg, og);
	} else {
		for (int i = 0; i < nfold; i++) {
			GMRFLib_graph_fold(&newg, oldg, og);
			GMRFLib_graph_free(oldg);
			if (i < nfold - 1) {
				oldg = newg;
				newg = NULL;
			}
		}
	}

	*ng = newg;
	return GMRFLib_SUCCESS;
}

int GMRFLib_graph_union(GMRFLib_graph_tp **union_graph, GMRFLib_graph_tp **graph_array, int n_graphs)
{
	GMRFLib_ged_tp *ged = NULL;
	GMRFLib_ged_init(&ged, NULL);
	for (int i = 0; i < n_graphs; i++) {
		GMRFLib_ged_insert_graph(ged, graph_array[i], 0);
	}
	GMRFLib_ged_build(union_graph, ged);
	GMRFLib_graph_prepare(*union_graph);
	GMRFLib_ged_free(ged);

	return GMRFLib_SUCCESS;
}

int GMRFLib_graph_complete(GMRFLib_graph_tp **n_graph, GMRFLib_graph_tp *graph)
{
	/*
	 * return a new graph that is complete: if i~j but j!~i, then add j~i 
	 */

	if (!graph) {
		*n_graph = NULL;
		return GMRFLib_SUCCESS;
	}

	/*
	 * setup new graph 
	 */
	GMRFLib_graph_mk_empty(n_graph);
	int *neigh_size = Calloc(graph->n, int);
	(*n_graph)->n = graph->n;
	(*n_graph)->nbs = Calloc(graph->n, int *);
	(*n_graph)->nnbs = Calloc(graph->n, int);

	for (int i = 0; i < (*n_graph)->n; i++) {
		(*n_graph)->nnbs[i] = graph->nnbs[i];
		neigh_size[i] = IMAX(16, 2 * (*n_graph)->nnbs[i]);
		(*n_graph)->nbs[i] = Calloc(neigh_size[i], int);

		if ((*n_graph)->nnbs[i]) {
			Memcpy((*n_graph)->nbs[i], graph->nbs[i], graph->nnbs[i] * sizeof(int));
		}
	}

	/*
	 * if i ~ j, then add j ~ i 
	 */
	for (int i = 0; i < (*n_graph)->n; i++) {
		for (int j = 0; j < (*n_graph)->nnbs[i]; j++) {
			int ii = (*n_graph)->nbs[i][j];
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
	int n_neig_tot = 0;
	for (int i = 0; i < (*n_graph)->n; i++) {
		if ((*n_graph)->nnbs[i]) {
			QSORT_FUN((*n_graph)->nbs[i], (size_t) (*n_graph)->nnbs[i], sizeof(int), GMRFLib_icmp);
			int k = 0;
			for (int j = 1; j < (*n_graph)->nnbs[i]; j++) {
				if ((*n_graph)->nbs[i][j] != (*n_graph)->nbs[i][k]) {
					(*n_graph)->nbs[i][++k] = (*n_graph)->nbs[i][j];
				}
			}
			(*n_graph)->nnbs[i] = k + 1;
			n_neig_tot += (*n_graph)->nnbs[i];
		}
	}

	if (n_neig_tot) {
		int *storage = Calloc(n_neig_tot, int);
		int storage_idx = 0;
		for (int i = 0; i < (*n_graph)->n; i++) {
			if ((*n_graph)->nnbs[i]) {
				Memcpy(&storage[storage_idx], (*n_graph)->nbs[i], sizeof(int) * (*n_graph)->nnbs[i]);
				Free((*n_graph)->nbs[i]);
				(*n_graph)->nbs[i] = &storage[storage_idx];
				storage_idx += (*n_graph)->nnbs[i];
			} else {
				Free((*n_graph)->nbs[i]);
			}
		}
	}
	GMRFLib_graph_prepare(*n_graph);
	Free(neigh_size);

	return GMRFLib_SUCCESS;
}

int GMRFLib_graph_insert(GMRFLib_graph_tp **new_graph, int n_new, int offset, GMRFLib_graph_tp *graph)
{
	/*
	 * insert graph into a larger graph with n_new nodes such that node `i' in graph corresponds to node `offset +i' in the 
	 * larger graph. n_new must be larger than graph->n, and `offset' such that `0<= offset', and `offset + graph->n <
	 * new_graph->n' 
	 */

	int n_neig, *hold = NULL;
	GMRFLib_graph_tp *g = NULL;

	n_neig = GMRFLib_isum(graph->n, graph->nnbs);
	GMRFLib_graph_mk_empty(&g);
	g->n = n_new;
	g->nnbs = Calloc(n_new, int);
	g->nbs = Calloc(n_new, int *);

	if (n_neig) {
		hold = Calloc(n_neig, int);
	} else {
		hold = NULL;
	}

	for (int i = 0, hold_idx = 0; i < graph->n; i++) {
		int ii = i + offset;
		g->nnbs[ii] = graph->nnbs[i];
		g->nbs[ii] = &hold[hold_idx];
		hold_idx += graph->nnbs[i];

		for (int j = 0; j < graph->nnbs[i]; j++) {
			g->nbs[ii][j] = graph->nbs[i][j] + offset;
		}
	}
	GMRFLib_graph_prepare(g);

	*new_graph = g;

	return GMRFLib_SUCCESS;
}

double GMRFLib_offset_Qfunc(int thread_id, int node, int nnode, double *UNUSED(values), void *arg)
{
	if (nnode < 0) {
		return NAN;
	}

	GMRFLib_offset_arg_tp *a = (GMRFLib_offset_arg_tp *) arg;
	if (IMIN(node, nnode) < a->offset || IMAX(node, nnode) >= a->offset + a->n) {
		return 0.0;
	}

	return (*a->Qfunc) (thread_id, node - a->offset, nnode - a->offset, NULL, a->Qfunc_arg);
}

int GMRFLib_offset(GMRFLib_offset_tp **off, int n_new, int offset, GMRFLib_graph_tp *graph, GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_arg)
{
	/*
	 * make an object containg the new graph, Qfunc and Qfunc_arg when an graph if inserted into a graph. 
	 */

	GMRFLib_offset_tp *val = NULL;
	GMRFLib_offset_arg_tp *offset_arg = NULL;

	offset_arg = Calloc(1, GMRFLib_offset_arg_tp);
	offset_arg->Qfunc = Qfunc;
	offset_arg->Qfunc_arg = Qfunc_arg;
	offset_arg->offset = offset;
	offset_arg->n = graph->n;

	val = Calloc(1, GMRFLib_offset_tp);
	GMRFLib_graph_insert(&(val->graph), n_new, offset, graph);

	val->Qfunc = GMRFLib_offset_Qfunc;
	val->Qfunc_arg = (void *) offset_arg;

	*off = val;

	return GMRFLib_SUCCESS;
}

int *GMRFLib_graph_cc(GMRFLib_graph_tp *g)
{
	/*
	 * return a vector of length n, indicating which connecting component each node belongs to 
	 */

	if (g == NULL || g->n == 0) {
		return NULL;
	}

	int n, *cc = NULL, ccc;
	char *visited = NULL;

	n = g->n;
	cc = Calloc(n, int);
	ccc = -1;					       /* the counter. yes, start at -1 */
	visited = Calloc(n, char);
	assert(visited);

	for (int i = 0; i < n; i++) {
		if (!visited[i]) {
			ccc++;
			GMRFLib_graph_cc_do(i, g, cc, visited, &ccc);
		}
	}
	Free(visited);

	return cc;
}

int GMRFLib_graph_cc_do(int node, GMRFLib_graph_tp *g, int *cc, char *visited, int *ccc)
{
	if (visited[node]) {				       /* I don't need this but include it for clarity */
		return GMRFLib_SUCCESS;
	}

	visited[node] = 1;
	cc[node] = *ccc;

	for (int i = 0; i < g->nnbs[node]; i++) {
		int nnode = g->nbs[node][i];
		if (!visited[nnode]) {			       /* faster to do a check here than to doit inside the funcall */
			GMRFLib_graph_cc_do(nnode, g, cc, visited, ccc);
		}
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_graph_add_sha(GMRFLib_graph_tp *g)
{
	if (g->n == 0) {
		g->sha = NULL;
		return GMRFLib_SUCCESS;
	}

	GMRFLib_SHA_TP c;
	unsigned char *md = Calloc(GMRFLib_SHA_DIGEST_LEN + 1, unsigned char);

	Memset(md, 0, GMRFLib_SHA_DIGEST_LEN + 1);
	GMRFLib_SHA_Init(&c);

	GMRFLib_SHA_IUPDATE(&(g->n), 1);
	GMRFLib_SHA_IUPDATE(&(g->nnz), 1);
	GMRFLib_SHA_IUPDATE(g->nnbs, g->n);
	for (int i = 0; i < g->n; i++) {
		GMRFLib_SHA_IUPDATE(g->nbs[i], g->nnbs[i]);
	}

	GMRFLib_SHA_Final(md, &c);
	md[GMRFLib_SHA_DIGEST_LEN] = '\0';
	g->sha = md;

	return GMRFLib_SUCCESS;
}
