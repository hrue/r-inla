
/* cgeneric.c
 * 
 * Copyright (C) 2021 Havard Rue
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
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
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
#define iDEBUG(msg_, i_) if (debug) printf("\tread_data: %s %d\n", msg_, i_)
#define dDEBUG(msg_, x_) if (debug) printf("\tread_data: %s %g\n", msg_, x_)
#define idDEBUG(msg_, idx_, x_) if (debug) printf("\tread_data: %s[%1d] %g\n", msg_, idx_, x_)
#define cDEBUG(msg_, c_) if (debug) printf("\tread_data: %s %s\n", msg_, c_)
#define ijxDEBUG(msg_, idx_, i_, j_, x_) if (debug) printf("\tread_data: %s[%1d] (%d, %d, %g)\n", msg_, idx_, i_, j_, x_)

#define READ_NAME(nm_) if (1) {						\
		int j;							\
		nread = fread((void *) &j, sizeof(int), (size_t) 1, fp); assert(nread == (size_t) 1); \
		data->name_ ## nm_[k] = Calloc(j + 1L, char);		\
		if (0) data->name_ ## nm_[k][j+1L] = '\0';		\
		nread = fread((void *) data->name_ ## nm_[k], sizeof(char), (size_t) (j + 1L), fp); assert(nread == (size_t) (j + 1L)); \
		cDEBUG("name", data->name_ ## nm_[k]);			\
	}

	FILE *fp;
	size_t nread;
	inla_cgeneric_data_tp *data = Calloc(1, inla_cgeneric_data_tp);
	int i, j, k, len;

	fp = fopen(filename, "rb");
	assert(fp);

	nread = fread((void *) &len, sizeof(int), (size_t) 1, fp);
	assert(nread == (size_t) 1);
	iDEBUG("Number of ints", len);
	data->n_ints = len;
	data->name_ints = Calloc(len, char *);
	data->ints = Calloc(len, int *);
	for (k = 0; k < len; k++) {
		READ_NAME(ints);
		nread = fread((void *) &j, sizeof(int), (size_t) 1, fp);
		assert(nread == (size_t) 1);
		iDEBUG("lenght", j);

		data->ints[k] = Calloc(j, int);
		nread = fread((void *) data->ints[k], sizeof(int), (size_t) j, fp);
		assert(nread == (size_t) j);
		for (i = 0; i < j; i++) {
			iDEBUG("contents", data->ints[k][i]);
		}
	}

	nread = fread((void *) &len, sizeof(int), (size_t) 1, fp);
	assert(nread == (size_t) 1);
	iDEBUG("Number of doubles", len);
	data->n_doubles = len;
	data->name_doubles = Calloc(len, char *);
	data->doubles = Calloc(len, double *);
	for (k = 0; k < len; k++) {
		READ_NAME(doubles);
		nread = fread((void *) &j, sizeof(int), (size_t) 1, fp);
		assert(nread == (size_t) 1);
		iDEBUG("lenght", j);

		data->doubles[k] = Calloc(j, double);
		nread = fread((void *) data->doubles[k], sizeof(double), (size_t) j, fp);
		assert(nread == (size_t) j);
		for (i = 0; i < j; i++) {
			dDEBUG("contents", data->doubles[k][i]);
		}
	}

	nread = fread((void *) &len, sizeof(int), (size_t) 1, fp);
	assert(nread == (size_t) 1);
	iDEBUG("Number of chars", len);
	data->n_chars = len;
	data->name_chars = Calloc(len, char *);
	data->chars = Calloc(len, char *);
	for (k = 0; k < len; k++) {
		READ_NAME(chars);
		nread = fread((void *) &j, sizeof(int), (size_t) 1, fp);
		assert(nread == (size_t) 1);
		// iDEBUG("lenght", j);

		data->chars[k] = Calloc(j + 1L, char);
		nread = fread((void *) data->chars[k], sizeof(char), (size_t) (j + 1L), fp);
		assert(nread == (size_t) (j + 1L));
		cDEBUG("contents", data->chars[k]);
	}

	nread = fread((void *) &len, sizeof(int), (size_t) 1, fp);
	assert(nread == (size_t) 1);
	iDEBUG("Number of mat", len);
	data->n_mat = len;
	data->name_mat = Calloc(len, char *);
	data->mat = Calloc(len, inla_cgeneric_mat_tp *);
	for (k = 0; k < len; k++) {
		READ_NAME(chars);
		data->mat[k] = Calloc(1, inla_cgeneric_mat_tp);

		int dim[2], nn;
		nread = fread((void *) dim, sizeof(int), (size_t) 2, fp);
		assert(nread == (size_t) 2);
		data->mat[k]->nrow = dim[0];
		data->mat[k]->ncol = dim[1];
		nn = dim[0] * dim[1];

		iDEBUG("nrow", data->mat[k]->nrow);
		iDEBUG("ncol", data->mat[k]->ncol);

		data->mat[k]->x = Calloc(nn, double);
		nread = fread((void *) data->mat[k]->x, sizeof(double), (size_t) nn, fp);
		assert(nread == (size_t) nn);
		for (i = 0; i < nn; i++) {
			idDEBUG("\tx", i, data->mat[k]->x[i]);
		}
	}

	nread = fread((void *) &len, sizeof(int), (size_t) 1, fp);
	assert(nread == (size_t) 1);
	iDEBUG("Number of smat", len);
	data->n_smat = len;
	data->name_smat = Calloc(len, char *);
	data->smat = Calloc(len, inla_cgeneric_smat_tp *);
	for (k = 0; k < len; k++) {
		READ_NAME(chars);
		data->smat[k] = Calloc(1, inla_cgeneric_smat_tp);

		int dim[3], n;
		nread = fread((void *) dim, sizeof(int), (size_t) 3, fp);
		assert(nread == (size_t) 3);
		data->smat[k]->nrow = dim[0];
		data->smat[k]->ncol = dim[1];
		data->smat[k]->n = n = dim[2];

		iDEBUG("nrow", data->smat[k]->nrow);
		iDEBUG("ncol", data->smat[k]->ncol);
		iDEBUG("n", data->smat[k]->n);

		data->smat[k]->i = Calloc(n, int);
		data->smat[k]->j = Calloc(n, int);
		data->smat[k]->x = Calloc(n, double);
		nread = fread((void *) data->smat[k]->i, sizeof(int), (size_t) n, fp);
		assert(nread == (size_t) n);
		nread = fread((void *) data->smat[k]->j, sizeof(int), (size_t) n, fp);
		assert(nread == (size_t) n);
		nread = fread((void *) data->smat[k]->x, sizeof(double), (size_t) n, fp);
		assert(nread == (size_t) n);

		for (i = 0; i < data->smat[k]->n; i++) {
			ijxDEBUG("\tx", i, data->smat[k]->i[i], data->smat[k]->j[i], data->smat[k]->x[i]);
		}
	}
	fclose(fp);

	return data;
}
