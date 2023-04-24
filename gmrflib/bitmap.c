
/* bitmap.c
 * 
 * Copyright (C) 2004-2023 Havard Rue
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

#include <stdio.h>
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

int GMRFLib_bitmap_image(const char *filename, GMRFLib_uchar *image, int nx, int ny)
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

int GMRFLib_bitmap_graph__intern(GMRFLib_graph_tp *graph, const char *filename, int *mapping)
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

int GMRFLib_bitmap_graph(const char *filename_body, int *remap, GMRFLib_graph_tp *graph)
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

int GMRFLib_bitmap_problem(const char *filename_body, GMRFLib_problem_tp *problem)
{
	GMRFLib_EWRAP0(GMRFLib_bitmap_graph(filename_body, problem->sub_sm_fact.remap, problem->sub_graph));
	GMRFLib_EWRAP0(GMRFLib_bitmap_factorisation(filename_body, &(problem->sub_sm_fact), problem->sub_graph));

	return GMRFLib_SUCCESS;
}
