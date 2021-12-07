
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

inla_cgeneric_data_tp * inla_cgeneric_read_data(const char *filename) 
{
#define iDEBUG(msg_, i_) if (debug) printf("\tread_data: %s %d\n", msg_, i_)
#define dDEBUG(msg_, x_) if (debug) printf("\tread_data: %s %g\n", msg_, x_)
#define idDEBUG(msg_, idx_, x_) if (debug) printf("\tread_data: %s[%1d] %g\n", msg_, idx_, x_)
#define cDEBUG(msg_, c_) if (debug) printf("\tread_data: %s %s\n", msg_, c_)
#define ijxDEBUG(msg_, idx_, i_, j_, x_) if (debug) printf("\tread_data: %s[%1d] (%d, %d, %g)\n", msg_, idx_, i_, j_, x_)
	
#define READ_NAME(nm_) if (1) {						\
		int j;							\
		fread((void *) &j, sizeof(int), (size_t) 1, fp);	\
		data->name_ ## nm_[k] = Calloc(j + 1L, char);		\
		if(0)data->name_ ## nm_[k][j+1L] = '\0';		\
		fread((void *) data->name_ ## nm_[k], sizeof(char), (size_t) (j + 1L), fp); \
		cDEBUG("name", data->name_ ## nm_[k]);			\
	}

	FILE *fp;
	inla_cgeneric_data_tp * data = Calloc(1, inla_cgeneric_data_tp);
	int i, j, k, len, debug = 1;

	fp = fopen(filename, "rb");
	assert(fp);

	fread((void *) &len, sizeof(int), (size_t) 1, fp);
	iDEBUG("Number of ints", len);
	data->n_ints = len;
	data->name_ints = Calloc(len, char *);
	data->ints = Calloc(len, int *);
	for(k = 0; k < len; k++) {
		READ_NAME(ints);
		fread((void *) &j, sizeof(int), (size_t) 1, fp);
		iDEBUG("lenght", j);

		data->ints[k] = Calloc(j, int);
		fread((void *) data->ints[k], sizeof(int), (size_t) j, fp);
		for(i = 0; i < j; i++) {
			iDEBUG("contents", data->ints[k][i]);
		}
	}

	fread((void *) &len, sizeof(int), (size_t) 1, fp);
	iDEBUG("Number of doubles", len);
	data->n_doubles = len;
	data->name_doubles = Calloc(len, char *);
	data->doubles = Calloc(len, double *);
	for(k = 0; k < len; k++) {
		READ_NAME(doubles);
		fread((void *) &j, sizeof(int), (size_t) 1, fp);
		iDEBUG("lenght", j);

		data->doubles[k] = Calloc(j, double);
		fread((void *) data->doubles[k], sizeof(double), (size_t) j, fp);
		for(i = 0; i < j; i++) {
			dDEBUG("contents", data->doubles[k][i]);
		}
	}

	fread((void *) &len, sizeof(int), (size_t) 1, fp);
	iDEBUG("Number of chars", len);
	data->n_chars = len;
	data->name_chars = Calloc(len, char *);
	data->chars = Calloc(len, char *);
	for(k = 0; k < len; k++) {
		READ_NAME(chars);
		fread((void *) &j, sizeof(int), (size_t) 1, fp);
		//iDEBUG("lenght", j);

		data->chars[k] = Calloc(j+1L, char);
		fread((void *) data->chars[k], sizeof(char), (size_t) (j+1L), fp);
		cDEBUG("contents", data->chars[k]);
	}

	fread((void *) &len, sizeof(int), (size_t) 1, fp);
	iDEBUG("Number of matrices", len);
	data->n_matrices = len;
	data->name_matrices = Calloc(len, char *);
	data->matrices = Calloc(len, inla_cgeneric_matrix_tp *);
	for(k = 0; k < len; k++) {
		READ_NAME(chars);
		data->matrices[k] = Calloc(1, inla_cgeneric_matrix_tp);

		int dim[2], nn;
		fread((void *) dim, sizeof(int), (size_t) 2, fp);
		data->matrices[k]->nrow = dim[0];
		data->matrices[k]->ncol = dim[1];
		nn = dim[0] * dim[1];

		iDEBUG("nrow", data->matrices[k]->nrow);
		iDEBUG("ncol", data->matrices[k]->ncol);

		data->matrices[k]->x = Calloc(nn, double);
		fread((void *) data->matrices[k]->x, sizeof(double), (size_t) nn, fp);
		for(i = 0; i < nn; i++) {
			idDEBUG("\tx", i, data->matrices[k]->x[i]);
		}
	}

	fread((void *) &len, sizeof(int), (size_t) 1, fp);
	iDEBUG("Number of smatrices", len);
	data->n_smatrices = len;
	data->name_smatrices = Calloc(len, char *);
	data->smatrices = Calloc(len, inla_cgeneric_smatrix_tp *);
	for(k = 0; k < len; k++) {
		READ_NAME(chars);
		data->smatrices[k] = Calloc(1, inla_cgeneric_smatrix_tp);

		int dim[3], n;
		fread((void *) dim, sizeof(int), (size_t) 3, fp);
		data->smatrices[k]->nrow = dim[0];
		data->smatrices[k]->ncol = dim[1];
		data->smatrices[k]->n = n = dim[2];

		iDEBUG("nrow", data->smatrices[k]->nrow);
		iDEBUG("ncol", data->smatrices[k]->ncol);
		iDEBUG("n", data->smatrices[k]->n);

		data->smatrices[k]->i = Calloc(n, int);
		data->smatrices[k]->j = Calloc(n, int);
		data->smatrices[k]->x = Calloc(n, double);
		fread((void *) data->smatrices[k]->i, sizeof(int), (size_t) n, fp);
		fread((void *) data->smatrices[k]->j, sizeof(int), (size_t) n, fp);
		fread((void *) data->smatrices[k]->x, sizeof(double), (size_t) n, fp);

		for(i = 0; i < data->smatrices[k]->n; i++) {
			ijxDEBUG("\tx", i, data->smatrices[k]->i[i], data->smatrices[k]->j[i], data->smatrices[k]->x[i]);
		}
	}
	fclose(fp);

	return data;
}
