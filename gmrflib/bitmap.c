
/* bitmap.c
 * 
 * Copyright (C) 2004-2006 Havard Rue
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
  \file bitmap.c
  \brief Functions to create bitmaps of precision matrices and Cholesky triangles.

  This file contains functions to create bitmaps (in the  pbm - portable bitmap file format) of
  - the precision matrix (or graph),
  - the reordered precision matrix (or graph)
  - and the Cholesky triangle.

  The main function is \c  GMRFLib_bitmap_problem() which provides a simple interface to create
  a bitmap of the precision matrix, the reordered precision matrix and the Cholesky triangle for
  a given <em> problem </em>.

  \sa GMRFLib_init_problem

  \par Example:
  \verbinclude example-doxygen-bitmap.txt

  which produce the following bitmaps
  \htmlonly
  <table>
  <tr><td align="center">
    <table>
    <tr>
    <td width="250" align="center"><img src="figs/example.gif" width="200" height="200"></td>
    <td width="250" align="center"><img src="figs/example-reordered.gif" width="200" height="200"></td>
    <td width="250" align="center"><img src="figs/example_L.gif" width="200" height="200"></td>
    </tr>
    </table></td></tr>
  </table>
  \endhtmlonly
  showing the precision matrix, the reordered precision matrix and the Cholesky triangle of the reordered precision matrix.

*/

#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

/* Pre-hg-Id: $Id: bitmap.c,v 1.23 2008/12/28 19:43:43 hrue Exp $ */


int GMRFLib_bitmap_image(const char *filename, GMRFLib_uchar * image, int nx, int ny)
{
	/*
	 * Create a PNB file of image, x-based storage 
	 */

	FILE *fp;
	int i, j, k, bit;
	GMRFLib_uchar c = 0, *iptr = image;

	fp = fopen(filename, "wb");
	if (fp) {
		fprintf(fp, "P4 %1d %1d ", nx, ny);
		for (j = 0; j < ny; j++) {
			for (i = 0; i < nx; i++) {
				k = 7 - i % 8;
				bit = (*iptr ? 1 : 0);
				c = c | (bit << k);
				iptr++;
				if ((i + 1) % 8 == 0 || i == (nx - 1)) {
					fputc(c, fp);
					c = 0;
				}

			}
		}
		fclose(fp);
	} else {
		GMRFLib_ERROR(GMRFLib_EOPENFILE);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_bitmap_graph__intern(GMRFLib_graph_tp * graph, const char *filename, int *mapping)
{
#define ROUND(_i) ((int) ((_i) * reduce_factor))
#define SET(_i, _j) bitmap[ROUND(_i) + ROUND(_j) * N] = 1

	int i, j, n, N, err, im, jm;
	double reduce_factor;
	GMRFLib_uchar *bitmap = NULL;

	n = graph->n;

	if (GMRFLib_bitmap_max_dimension > 0 && n > GMRFLib_bitmap_max_dimension) {
		N = GMRFLib_bitmap_max_dimension;
		reduce_factor = (double) N / (double) n;
	} else {
		N = n;
		reduce_factor = 1.0;
	}

	bitmap = Calloc(ISQR(N), GMRFLib_uchar);
	for (i = 0; i < graph->n; i++) {
		im = mapping[i];
		SET(im, im);
		for (j = 0; j < graph->nnbs[i]; j++) {
			jm = mapping[graph->nbs[i][j]];
			SET(im, jm);
		}
	}

	err = GMRFLib_bitmap_image(filename, bitmap, N, N);
	Free(bitmap);

#undef SET
#undef ROUND
	return err;
}

/*!
  \brief Create bitmaps of the original and reordered precision matrix (or the graph) in the pbm
  (portable bitmap) file format.

  \param[in] filename_body The body of the filename to store the original and reordered precision
  matrix (or the graph). The files are <em>filename_body</em>.pbm and
  <em>filename_body</em>-reordered.pbm

  \param[in] remap The reordering of the indices, so the new position of node <em>i</em> is
    <em>remap[i]</em>

  \param[in] graph The graph
    
  \remarks If <em>filename_body</em> is \em NULL, then the name \em graph is used.
    
  \remarks Normal usage is to use the \c GMRFLib_bitmap_problem() instead of \c
  GMRFLib_bitmap_graph()
    
  \sa GMRFLib_bitmap_problem
*/
int GMRFLib_bitmap_graph(const char *filename_body, int *remap, GMRFLib_graph_tp * graph)
{
	int i, *mapping;
	char *filename;

	mapping = Calloc(graph->n, int);
	for (i = 0; i < graph->n; i++) {
		mapping[i] = i;
	}

	GMRFLib_EWRAP0(GMRFLib_sprintf(&filename, "%s.pbm", (filename_body ? filename_body : "graph")));
	GMRFLib_EWRAP0(GMRFLib_bitmap_graph__intern(graph, filename, mapping));
	Free(mapping);
	Free(filename);

	if (remap) {
		GMRFLib_EWRAP0(GMRFLib_sprintf(&filename, "%s-reordered.pbm", (filename_body ? filename_body : "graph")));
		GMRFLib_EWRAP0(GMRFLib_bitmap_graph__intern(graph, filename, remap));
		Free(filename);
	}

	return GMRFLib_SUCCESS;
}

/*!  \brief Create bitmaps of the precision matrix (or graph), the reordered precision matrix (or
  graph) and the Cholesky triangle, for a <em>problem</em>

  \param[in] filename_body The body of the filename to store the orignal and reordered precision
  matrix (or the graph). The files are <em>filename_body</em>.pbm,
  <em>filename_body</em>-reordered.pbm and <em>filename_body</em>_L.pbm
  
  \param[in] problem A pointer to the \em problem
  
  \par Example:
  
  The call <tt>GMRFLib_bitmap_problem("Q", problem);</tt> will create the bitmaps <em>Q.pbm</em>,
  <em>Q-reordered.pbm</em> and <em>Q_L.pbm</em>
  
  \sa GMRFLib_init_problem
*/
int GMRFLib_bitmap_problem(const char *filename_body, GMRFLib_problem_tp * problem)
{
	GMRFLib_EWRAP0(GMRFLib_bitmap_graph(filename_body, problem->sub_sm_fact.remap, problem->sub_graph));
	GMRFLib_EWRAP0(GMRFLib_bitmap_factorisation(filename_body, &(problem->sub_sm_fact), problem->sub_graph));

	return GMRFLib_SUCCESS;
}
