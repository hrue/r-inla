
/* experimental.c
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

/*!
  \file experimental.c
  \brief Various functions that may be moved into GMRFLib.

*/

#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

/* Pre-hg-Id: $Id: experimental.c,v 1.16 2009/01/01 10:40:00 hrue Exp $ */

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

/*!
  \brief Reads a graph from a file, binary format

  \param[in,out] graph At output, <em>(*graph)</em> has been initialized 
  by using the information on the file \em filename.

  \param[in] filename The name of the file, binary formatted similar to
  \c GMRFLib_read_graph() (without newlines), containing the specification
  of the graph.

  \sa GMRFLib_read_graph, GMRFLib_write_graph_binary
 */
int GMRFLib_read_graph_binary_EXPERIMENTAL(GMRFLib_graph_tp ** graph, const char *filename)
{

	/*
	 * read a graph in the following format
	 * 
	 * N node[0] nnbs[0] nbs[node[0]][0] nbs[node[0]][1] ... nbs[node[0]][nnbs[0]-1] node[1] nnbs[1] nbs[node[1]][0]
	 * nbs[node[1]][1] ... nbs[node[1]][nnbs[1]-1] : node[N-1] nnbs[N-1] nbs[node[N-1]][0] nbs[node[N-1]][1] ...
	 * nbs[node[N-1]][nnbs[N-1]-1] 
	 */

#define READ_ERROR() do { \
    if (1){\
        fprintf(stderr,"\n\n\t%s: error: file [%s]:\n",__GMRFLib_FuncName,filename);\
	fprintf(stderr,"\t\tfail to read [%1lu] bytes from byte [%1lu]\n", (unsigned long)(nelm*sizeof(int)), (unsigned long)byte); \
    }    \
    if (fp) fclose(fp);\
    (*graph) = NULL;\
    GMRFLib_ERROR(GMRFLib_EREADFILE);} while(0)

	int *storage, n_neig_tot = 0, storage_indx;
	int i, byte, tnode;
	size_t nelm;
	FILE *fp;

	if (!filename) {
		return GMRFLib_SUCCESS;
	}
	if (!(fp = fopen(filename, "r"))) {
		GMRFLib_ERROR(GMRFLib_EOPENFILE);
	}

	GMRFLib_make_empty_graph(graph);

	byte = 0;
	nelm = 1;
	if (fread(&((*graph)->n), sizeof(int), nelm, fp) != 1) {
		READ_ERROR();
	}
	byte += sizeof(int);

	(*graph)->nnbs = Calloc((*graph)->n, int);
	(*graph)->nbs = Calloc((*graph)->n, int *);

	for (i = 0; i < (*graph)->n; i++) {
		nelm = 1;
		if (fread(&tnode, sizeof(int), nelm, fp) != 1) {
			READ_ERROR();			       /* target node */
		}
		byte += sizeof(int);
		if (tnode < 0 || tnode >= (*graph)->n) {
			fprintf(stderr, "\n\n\t%s: error: file [%s]: byte[%1d]\n", __GMRFLib_FuncName, filename, byte);
			fprintf(stderr, "\t\tnode-number[%1d] is not in the range[0:%1d]\n", tnode, (*graph)->n);
			if (fp) {
				fclose(fp);
			}
			GMRFLib_ERROR(GMRFLib_EPARAMETER);
		}

		nelm = 1;
		if (fread(&((*graph)->nnbs[tnode]), sizeof(int), nelm, fp) != 1) {
			READ_ERROR();
		}
		byte += sizeof(int);
		n_neig_tot += (*graph)->nnbs[tnode];
		if ((*graph)->nnbs[tnode]) {
			(*graph)->nbs[tnode] = Calloc((*graph)->nnbs[tnode], int);

			nelm = (*graph)->nnbs[tnode];
			if (fread((*graph)->nbs[tnode], sizeof(int), nelm, fp) != nelm) {
				READ_ERROR();
			}
			byte += sizeof(int) * ((*graph)->nnbs[tnode]);
		} else {
			(*graph)->nbs[tnode] = NULL;
		}
	}
	fclose(fp);

	/*
	 * map the graph to a more computational convenient memory layout! use just one long vector to store all the neighbors. 
	 */
	storage = Calloc(n_neig_tot, int);

	storage_indx = 0;
	for (i = 0; i < (*graph)->n; i++)
		if ((*graph)->nnbs[i]) {
			memcpy(&storage[storage_indx], (*graph)->nbs[i], sizeof(int) * (*graph)->nnbs[i]);
			Free((*graph)->nbs[i]);
			(*graph)->nbs[i] = &storage[storage_indx];
			storage_indx += (*graph)->nnbs[i];
		} else {
			(*graph)->nbs[i] = NULL;
		}

	if (GMRFLib_verify_graph_read_from_disc) {
		GMRFLib_EWRAP0(GMRFLib_validate_graph(stderr, *graph));
	}

	GMRFLib_EWRAP0(GMRFLib_prepare_graph(*graph));	       /* prepare the graph for computations */
	return GMRFLib_SUCCESS;
#undef READ_ERROR
}

/*!
  \brief Write a graph to file, binary format
  \param[in] filename The name of the file to store the graph in
    binary format. Can be read with \c GMRFLib_read_graph_binary().
  \param[in] graph The graph
  \sa GMRFLib_read_graph_binary, GMRFLib_write_graph.
 */
int GMRFLib_write_graph_binary_EXPERIMENTAL(const char *filename, GMRFLib_graph_tp * graph)
{
	/*
	 * write graph to file filename in binary format so it can be read by 'read_graph_binary' 
	 */

	int i;
	FILE *fp;
	size_t ret;

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
	ret = fwrite(&(graph->n), sizeof(int), 1, fp);
	for (i = 0; i < graph->n; i++) {
		ret = fwrite(&i, sizeof(int), 1, fp);
		ret = fwrite(&(graph->nnbs[i]), sizeof(int), 1, fp);
		if (graph->nnbs[i]) {
			ret = fwrite(graph->nbs[i], sizeof(int), (unsigned int) graph->nnbs[i], fp);
		}
	}
	fclose(fp);

	return GMRFLib_SUCCESS;
}
