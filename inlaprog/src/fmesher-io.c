
/* fmesher-io.c
 * 
 * Copyright (C) 2010 Havard Rue
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
 *       H{\aa}vard Rue
 *       Department of Mathematical Sciences
 *       The Norwegian University of Science and Technology
 *       N-7491 Trondheim, Norway
 *       Voice: +47-7359-3533    URL  : http://www.math.ntnu.no/~hrue  
 *       Fax  : +47-7359-3524    Email: havard.rue@math.ntnu.no
 *
 */
#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"
#include "fmesher-io.h"

inla_matrix_tp *inla_read_fmesher_file(const char *filename)
{
#define ERROR(msg)							\
	{								\
		fprintf(stderr, "\n\n%s:%1d: *** ERROR *** \n\t%s\n\n", __FILE__,  __LINE__,  msg); \
		exit(EXIT_FAILURE);					\
		return (inla_matrix_tp *)NULL;				\
	}
	
#define READ(ptr, n, type)						\
	{								\
		long int position;					\
		size_t nread;						\
									\
		position = ftell(fp);					\
		nread = fread((void *)ptr, sizeof(type), (size_t) n, (FILE *) fp); \
		if (nread != (size_t) n){				\
			char *m;					\
			GMRFLib_sprintf(&m, "Fail to read [%1u] elems of size [%1u] from file [%s], at position %ld\n", \
					n, sizeof(type), filename, position); \
			ERROR(m);					\
		}							\
	}
	
	FILE *fp = NULL;
	char *msg = NULL; 
	int *header =  NULL;
	int len_header = 0;
	int verbose = 1, debug = 0, i, j, k;
	inla_matrix_tp *M = NULL;

	if (debug)
		verbose = 1;
	
	fp = fopen(filename,  "rb");
	if (!fp){
		GMRFLib_sprintf(&msg, "Fail to open file [%s]", filename);
		ERROR(msg);
	}
	if (verbose)
		printf("Open file [%s]\n", filename);
	
	READ(&len_header, 1, int);
	if (len_header < 8){
		GMRFLib_sprintf(&msg, "Header in file [%s] is only %1d (< 8) ints long.", filename, len_header);
		ERROR(msg);
	}

	header = Calloc(len_header, int);
	READ(header, len_header, int);
	
	M = Calloc(1, inla_matrix_tp);

	int version = header[0];
	int elems = header[1];
	int nrow = header[2];
	int ncol = header[3];
	int datatype = header[4];
	int valuetype = header[5];
	int matrixtype = header[6];
	int storagetype = header[7];

	if (verbose){
		printf("\tversion     \t%d\n", version);
		printf("\telems       \t%d\n", elems);
		printf("\tnrow        \t%d\n", nrow);
		printf("\tncol        \t%d\n", ncol);
		printf("\tdatatype    \t%d\n", datatype);
		printf("\tvaluetype   \t%d\n", valuetype);
		printf("\tmatrixtype  \t%d\n", matrixtype);
		printf("\tstoragetype \t%d\n", storagetype);
	}

	M->nrow = nrow;
	M->ncol = ncol;
	M->elems = elems;

	int rowmajor = (storagetype == 0);
	int integer = (valuetype == 0);
	int dense = (datatype == 0);
	int general = (matrixtype == 0);
	int symmetric = (matrixtype == 1);
	int diagonal = (matrixtype == 2);
	
	if (dense) {
		if (symmetric || diagonal){
			ERROR(" (dense && (symmetric || diagonal)) is not yet implemented.");
		}
		assert(general);
		
		M->values = Calloc(elems, double);
		if (integer){
			int *ivalues = NULL;

			ivalues = Calloc(elems, int);
			READ(ivalues, elems, int);
			for(i = 0; i<elems; i++){
				M->values[i] = (double) ivalues[i];
			}
			Free(ivalues);
		} else {
			READ(M->values, elems, double);
		}
		if (rowmajor){
			/* 
			   swap
			*/
			
			double *swapped_values = Calloc(elems, double);

			for(i=0; i<nrow; i++) {
				for(j=0; j<ncol; j++) {
					int idx = i + j*nrow;
					int idxx = j + i*ncol;

					swapped_values[idxx] = M->values[idx];
				}
			}
			Free(M->values);
			M->values = swapped_values;
		}

		if (debug) {
			printf("\n\t%d x %d\n", M->nrow,  M->ncol);
			for(i=0; i<nrow; i++){
				printf("\t");
				for(j=0; j<ncol; j++){
					printf("%.3f ", M->values[i + j * nrow]);
				}
				printf("\n");
			}
		}
		
	} else {
		/* 
		   sparse
		*/
		
		M->i = Calloc(elems, int);
		M->j = Calloc(elems, int);
		M->values = Calloc(elems, double);

		if (rowmajor) {
			for(k=0; k<elems; k++){
				READ(&(M->i[k]), 1, int);
				READ(&(M->j[k]), 1, int);
				if (integer){
					int itmp;
					READ(&itmp, 1, int);
					M->values[k] = (double) itmp;
				} else {
					READ(&(M->values[k]), 1, double);
				}
			}
		} else {
			READ(M->i, elems, int);
			READ(M->j, elems, int);
			if (integer){
				int *ivalues = Calloc(elems, int);
				READ(ivalues, elems, int);
				for(k=0; k<elems; k++){
					M->values[k] = (double) ivalues[k];
				}
				Free(ivalues);
			} else {
				READ(M->values, elems, double);
			}
		}
		if (symmetric) {
			int nneq = 0;

			for(k=0; k<elems; k++){
				if (M->i[k] != M->j[k])
					nneq++;
			}
			if (nneq > 0){
				M->i = Realloc(M->i, elems + nneq, int);
				M->j = Realloc(M->j, elems + nneq, int);
				M->values = Realloc(M->values, elems + nneq, double);

				int kk = elems;
				for(k=0; k<elems; k++) {
					if (M->i[k] != M->j[k]){
						M->i[kk] = M->j[k]; /* yes */
						M->j[kk] = M->i[k]; /* yes */
						M->values[kk] = M->values[k];
						kk++;
					}
					assert(kk == elems + nneq);
				}
				M->elems += nneq;
			}
		}

		if (debug) {
			double *A = Calloc(M->nrow * M->nrow, double);
			double fix= GMRFLib_uniform();

			if (0){
				for(k=0; k<elems; k++)
					printf("%d: i j values %d %d %f\n", k, M->i[k], M->j[k], M->values[k]);
			}

			for(k=0; k<M->nrow*M->ncol; k++)
				A[k] = fix;
			for(k=0; k<M->elems; k++){
				A[ M->i[k] + M->j[k] * M->nrow ] = M->values[k];
			}
				
			printf("\n\t%d x %d\n", M->nrow,  M->ncol);
			for(i=0; i<M->nrow; i++){
				printf("\t");
				for(j=0; j<M->ncol; j++){
					int idx = i+j*M->nrow;

					if(A[idx] != fix)
						printf("%.3f ", A[idx]);
					else
						printf("      ");
				}
				printf("\n");
			}
		}
	}

#undef READ
#undef ERROR

	return (M);
}

#ifdef TESTME
int main(int argc, char **argv)
{
	int i;
	for(i=1; i< argc; i++)
		inla_read_fmesher_file(argv[i]);
	return 0;
}
#endif
