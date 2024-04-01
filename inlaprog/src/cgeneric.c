
/* cgeneric.c
 * 
 * Copyright (C) 2021-2024 Havard Rue
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
#include <math.h>
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"
#include "cgeneric.h"

#if !defined(Calloc)
#define Calloc(n_, type_)  (type_ *)calloc((n_), sizeof(type_))
#endif

inla_cgeneric_data_tp *inla_cgeneric_read_data(const char *filename, int debug)
{
#define READ_NAME(nm_) \
	if (1) {							\
		int j_;							\
		nread = fread((void *) &j_, sizeof(int), (size_t) 1, fp); assert(nread == (size_t) 1); \
		nm_ = Calloc(j_ + 1L, char);				\
		nread = fread((void *) nm_, sizeof(char), (size_t) (j_ + 1L), fp); assert(nread == (size_t) (j_ + 1L)); \
	}

	FILE *fp = NULL;
	size_t nread;
	inla_cgeneric_data_tp *data = Calloc(1, inla_cgeneric_data_tp);
	int i, j, k, len;

	data->threads.max = GMRFLib_MAX_THREADS();
	data->threads.outer = GMRFLib_openmp->max_threads_nested[0];
	data->threads.inner = GMRFLib_openmp->max_threads_nested[1];

	fp = fopen(filename, "rb");
	assert(fp);

	nread = fread((void *) &len, sizeof(int), (size_t) 1, fp);
	assert(nread == (size_t) 1);
	if (debug) {
		printf("\tNumber of ints %1d\n", len);
	}
	data->n_ints = len;
	data->ints = Calloc(len, inla_cgeneric_vec_tp *);
	for (k = 0; k < len; k++) {
		data->ints[k] = Calloc(1, inla_cgeneric_vec_tp);
		READ_NAME(data->ints[k]->name);
		nread = fread((void *) &j, sizeof(int), (size_t) 1, fp);
		assert(nread == (size_t) 1);
		data->ints[k]->len = j;
		data->ints[k]->ints = Calloc(j, int);
		nread = fread((void *) data->ints[k]->ints, sizeof(int), (size_t) j, fp);
		assert(nread == (size_t) j);
		if (debug) {
			printf("\t\tname[%s] length[%1d] type[%s]\n", data->ints[k]->name, data->ints[k]->len, "INT");
			for (i = 0; i < j; i++) {
				printf("\t\t\tidx=%1d %d\n", i, data->ints[k]->ints[i]);
			}
		}
	}

	nread = fread((void *) &len, sizeof(int), (size_t) 1, fp);
	assert(nread == (size_t) 1);
	if (debug) {
		printf("\tNumber of doubles %1d\n", len);
	}
	data->n_doubles = len;
	data->doubles = Calloc(len, inla_cgeneric_vec_tp *);
	for (k = 0; k < len; k++) {
		data->doubles[k] = Calloc(1, inla_cgeneric_vec_tp);
		READ_NAME(data->doubles[k]->name);
		nread = fread((void *) &j, sizeof(int), (size_t) 1, fp);
		assert(nread == (size_t) 1);

		data->doubles[k]->len = j;
		data->doubles[k]->doubles = Calloc(j, double);
		nread = fread((void *) data->doubles[k]->doubles, sizeof(double), (size_t) j, fp);
		assert(nread == (size_t) j);
		if (debug) {
			printf("\t\tname[%s] length[%1d] type[%s]\n", data->doubles[k]->name, data->doubles[k]->len, "DOUBLE");
			for (i = 0; i < j; i++) {
				printf("\t\t\tidx=%1d %g\n", i, data->doubles[k]->doubles[i]);
			}
		}
	}

	nread = fread((void *) &len, sizeof(int), (size_t) 1, fp);
	assert(nread == (size_t) 1);
	if (debug) {
		printf("\tNumber of chars %1d\n", len);
	}
	data->n_chars = len;
	data->chars = Calloc(len, inla_cgeneric_vec_tp *);
	for (k = 0; k < len; k++) {
		data->chars[k] = Calloc(1, inla_cgeneric_vec_tp);
		READ_NAME(data->chars[k]->name);
		nread = fread((void *) &j, sizeof(int), (size_t) 1, fp);
		assert(nread == (size_t) 1);

		data->chars[k]->len = j;
		data->chars[k]->chars = Calloc(j + 1L, char);
		nread = fread((void *) data->chars[k]->chars, sizeof(char), (size_t) (j + 1L), fp);
		assert(nread == (size_t) (j + 1L));
		if (debug) {
			printf("\t\tname[%s] length[%1d] type[%s]\n", data->chars[k]->name, data->chars[k]->len, "CHAR");
			printf("\t\t\t%s\n", data->chars[k]->chars);
		}
	}

	nread = fread((void *) &len, sizeof(int), (size_t) 1, fp);
	assert(nread == (size_t) 1);
	if (debug) {
		printf("\tNumber of matrices %1d\n", len);
	}
	data->n_mats = len;
	data->mats = Calloc(len, inla_cgeneric_mat_tp *);
	for (k = 0; k < len; k++) {
		data->mats[k] = Calloc(1, inla_cgeneric_mat_tp);
		READ_NAME(data->mats[k]->name);

		int dim[2], nn;
		nread = fread((void *) dim, sizeof(int), (size_t) 2, fp);
		assert(nread == (size_t) 2);
		data->mats[k]->nrow = dim[0];
		data->mats[k]->ncol = dim[1];
		nn = dim[0] * dim[1];

		data->mats[k]->nrow = dim[0];
		data->mats[k]->ncol = dim[1];
		data->mats[k]->x = Calloc(nn, double);
		nread = fread((void *) data->mats[k]->x, sizeof(double), (size_t) nn, fp);
		assert(nread == (size_t) nn);
		if (debug) {
			printf("\t\tname[%s] nrow[%1d] ncol[%1d] type[%s]\n", data->mats[k]->name, data->mats[k]->nrow, data->mats[k]->ncol,
			       "MATRIX");
			for (i = 0; i < nn; i++) {
				printf("\t\t\tidx=%1d %g\n", i, data->mats[k]->x[i]);
			}
		}
	}

	nread = fread((void *) &len, sizeof(int), (size_t) 1, fp);
	assert(nread == (size_t) 1);
	if (debug) {
		printf("\tNumber of sparse matrices %1d\n", len);
	}
	data->n_smats = len;
	data->smats = Calloc(len, inla_cgeneric_smat_tp *);
	for (k = 0; k < len; k++) {
		data->smats[k] = Calloc(1, inla_cgeneric_smat_tp);
		READ_NAME(data->smats[k]->name);
		int dim[3], n;
		nread = fread((void *) dim, sizeof(int), (size_t) 3, fp);
		assert(nread == (size_t) 3);
		data->smats[k]->nrow = dim[0];
		data->smats[k]->ncol = dim[1];
		data->smats[k]->n = n = dim[2];
		data->smats[k]->i = Calloc(n, int);
		data->smats[k]->j = Calloc(n, int);
		data->smats[k]->x = Calloc(n, double);
		nread = fread((void *) data->smats[k]->i, sizeof(int), (size_t) n, fp);
		assert(nread == (size_t) n);
		nread = fread((void *) data->smats[k]->j, sizeof(int), (size_t) n, fp);
		assert(nread == (size_t) n);
		nread = fread((void *) data->smats[k]->x, sizeof(double), (size_t) n, fp);
		assert(nread == (size_t) n);
		if (debug) {
			printf("\t\tname[%s] nrow[%1d] ncol[%1d] n[%1d] type[%s]\n",
			       data->smats[k]->name, data->smats[k]->nrow, data->smats[k]->ncol, data->smats[k]->n, "SPARSE.MATRIX");
			for (i = 0; i < n; i++) {
				printf("\t\t\tidx=%1d (i=%1d j=%1d x=%g)\n", i, data->smats[k]->i[i], data->smats[k]->j[i], data->smats[k]->x[i]);
			}
		}
	}
	fclose(fp);
#undef READ_NAME

	return data;
}
