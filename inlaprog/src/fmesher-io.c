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
	/* 
	   read a fmesher_file. 
	 */

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
		if (nread != (size_t) n) {				\
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
	int verbose = 0, debug = 0, i, j, k;
	inla_matrix_tp *M = NULL;

	if (debug)
		verbose = 1;
	
	fp = fopen(filename,  "rb");
	if (!fp) {
		GMRFLib_sprintf(&msg, "Fail to open file [%s]", filename);
		ERROR(msg);
	}
	if (verbose)
		printf("Open file [%s]\n", filename);
	
	READ(&len_header, 1, int);
	if (len_header < 8) {
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

	Free(header);

	if (verbose) {
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
		if (symmetric || diagonal) {
			ERROR(" (dense && (symmetric || diagonal)) is not yet implemented.");
		}
		assert(general);
		assert(elems == nrow*ncol);
		
		M->A = Calloc(elems, double);
		if (integer) {
			M->iA = Calloc(elems, int);
			READ(M->iA, elems, int);
			for(i = 0; i<elems; i++) {
				M->A[i] = (double) M->iA[i];
			}
		} else {
			READ(M->A, elems, double);
		}
		if (rowmajor) {
			/* 
			   swap
			*/
			double *swapped_values = Calloc(elems, double);

			for(i=0; i<nrow; i++) {
				for(j=0; j<ncol; j++) {
					swapped_values[j + i*ncol] = M->A[i + j*nrow];
				}
			}
			Free(M->A);
			M->A = swapped_values;

			if (integer) {
				for(k=0; k<elems; k++) {
					M->iA[k] = (int) M->A[k];
				}
			}
		}

		if (debug) {
			printf("\n\t%d x %d %s\n", M->nrow,  M->ncol,  (integer ? "integer" :  "double"));
			for(i=0; i<nrow; i++) {
				printf("\t");
				if (integer) {
					for(j=0; j<ncol; j++) {
						printf("%3d ", M->iA[i + j * nrow]);
					}
				} else {
					for(j=0; j<ncol; j++) {
						printf("%.3f ", M->A[i + j * nrow]);
					}
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
		if (integer) {
			M->ivalues = Calloc(elems, int);
		}
		
		if (rowmajor) {
			for(k=0; k<elems; k++) {
				READ(&(M->i[k]), 1, int);
				READ(&(M->j[k]), 1, int);
				if (integer) {
					READ(&(M->ivalues[k]), 1, int);
					M->values[k] = (double) M->ivalues[k];
				} else {
					READ(&(M->values[k]), 1, double);
				}
			}
		} else {
			READ(M->i, elems, int);
			READ(M->j, elems, int);
			if (integer) {
				READ(M->ivalues, elems, int);
				for(k=0; k<elems; k++) {
					M->values[k] = (double) M->ivalues[k];
				}
			} else {
				READ(M->values, elems, double);
			}
		}

		if (symmetric) {
			/* 
			   In the symmetric case, I need to check if either all M_ij, i>=j, or all M_ij, i<=j, is given. What to do if both sides of the matrix is
			   given is not spesificed, or, we need to define what we mean in this case. So let's do that if that case appear... he he.
			*/
			
			int all_i_st_j = 1,  all_j_st_i = 1;

			for(k=0; k<elems; k++) {
				int ij[2];

				all_i_st_j &= (M->i[k] <= M->j[k]);
				all_j_st_i &= (M->j[k] <= M->i[k]);
				
				ij[0] = M->i[k];
				ij[1] = M->j[k];

				M->i[k] = IMIN(ij[0], ij[1]);
				M->j[k] = IMAX(ij[0], ij[1]);
			}
			if ((all_i_st_j && all_j_st_i) || (!all_i_st_j && !all_j_st_i)) {
				ERROR("Not all entries satisfy all i >= j, or all j <= i, in the symmetric case. Do not know what to do...");
			}

			int nneq = 0;
			for(k=0; k<elems; k++) {
				if (M->i[k] != M->j[k]) {
					nneq++;
				}
			}
			if (nneq > 0) {
				M->i = Realloc(M->i, elems + nneq, int);
				M->j = Realloc(M->j, elems + nneq, int);
				M->values = Realloc(M->values, elems + nneq, double);
				if (integer) {
					M->ivalues = Realloc(M->ivalues, elems + nneq, int);
				}

				int kk = elems;
				for(k=0; k<elems; k++) {
					if (M->i[k] != M->j[k]) {
						M->i[kk] = M->j[k]; /* yes */
						M->j[kk] = M->i[k]; /* yes */
						M->values[kk] = M->values[k];
						if (integer) {
							M->ivalues[kk] = M->ivalues[k];
						}
						kk++;
					}
				}
				assert(kk == elems + nneq);
				M->elems += nneq;
			}
		}

		if (debug) {
			double *A = Calloc(M->nrow * M->nrow, double);
			int *iA = Calloc(M->nrow * M->nrow, int);
			double fix= GMRFLib_uniform();

			if (1) {
				if (integer) {
					for(k=0; k<elems; k++)
						printf("\t%d: i j values %d %d %d\n", k, M->i[k], M->j[k], M->ivalues[k]);
				} else {
					for(k=0; k<elems; k++)
						printf("\t%d: i j values %d %d %f\n", k, M->i[k], M->j[k], M->values[k]);
				}
			}

			for(k=0; k<M->nrow*M->ncol; k++) {
				A[k] = fix;
			}
			if (integer) {
				for(k=0; k<M->elems; k++) {
					iA[ M->i[k] + M->j[k] * M->nrow ] = M->ivalues[k];
				}
			} 
			for(k=0; k<M->elems; k++) {
				A[ M->i[k] + M->j[k] * M->nrow ] = M->values[k];
			}
				
			printf("\n\t%d x %d %s\n", M->nrow,  M->ncol, (integer ? "integer" : "double"));
			for(i=0; i<M->nrow; i++) {
				printf("\t");
				for(j=0; j<M->ncol; j++) {
					int idx = i+j*M->nrow;

					if(A[idx] != fix) {
						if (integer) {
							printf("%4d ", iA[idx]);
						} else {
							printf("%.3f ", A[idx]);
						}
					} else {
						printf(" .    ");
					}
				}
				printf("\n");
			}
			Free(A);
			Free(iA);
		}
	}
#undef READ
#undef ERROR

	return (M);
}
int inla_write_fmesher_file(inla_matrix_tp *M, const char *filename)
{
	/* 
	   write fmesher-file
	*/

#define ERROR(msg)							\
	{								\
		fprintf(stderr, "\n\n%s:%1d: *** ERROR *** \n\t%s\n\n", __FILE__,  __LINE__,  msg); \
		exit(EXIT_FAILURE);					\
		return 1;						\
	}
	
#define WRITE(ptr, n, type)						\
	{								\
		size_t nwrite;						\
									\
		nwrite = fwrite((const void *)ptr, sizeof(type), (size_t) n, (FILE *) fp); \
		if (nwrite != (size_t) n) {				\
			char *m;					\
			GMRFLib_sprintf(&m, "Fail to write [%1u] elems of size [%1u] from file [%s], at position %ld\n", \
					n, sizeof(type), filename, ftell(fp)); \
			ERROR(m);					\
		}							\
	}
	
	FILE *fp = NULL;
	char *msg = NULL; 
	int *header =  NULL;
	int len_header;
	int verbose = 0, dense, integer;
	int i;
	
	if (!M)
		return 0;
	
	fp = fopen(filename,  "wb");
	if (!fp) {
		GMRFLib_sprintf(&msg, "Fail to open file [%s]", filename);
		ERROR(msg);
	}
	if (verbose)
		printf("Open file [%s]\n", filename);
	
	len_header = 8;
	header = Calloc(len_header, int);

	if (M->A || M->iA) {
		dense = 1;
		integer = (M->iA ? 1 : 0);
	} else {
		dense = 0;
		integer = (M->ivalues ? 1 : 0);
	}

	header[0] = 0;					       /* version */
	header[1] = M->elems;
	header[2] = M->nrow;
	header[3] = M->ncol;
	header[4] = (dense ? 0 : 1);
	header[5] = (integer ? 0 : 1);
	header[6] = 0;					       /* general */
	header[7] = 1;					       /* columnmajor */

	if (verbose) {
		for(i=0; i<len_header; i++)
			printf("\theader[%1d] = %1d\n", i, header[i]);
	}
	
	WRITE(&len_header, 1, int);
	WRITE(header, len_header, int);

	if (dense) {
		if (integer) {
			WRITE(M->iA, M->elems, int);
		} else {
			WRITE(M->A, M->elems, double);
		}
	} else {
		WRITE(M->i, M->elems, int);
		WRITE(M->j, M->elems, int);
		if (integer) {
			WRITE(M->ivalues, M->elems, int);
		} else {
			WRITE(M->values, M->elems, double);
		}
	}

	fclose(fp);
	Free(header);
	
#undef ERROR
#undef WRITE	
	return (0);
}
double *inla_matrix_get_diagonal(inla_matrix_tp *M)
{
	/* 
	   return the diagonal of the matrix as a new and alloced double vector.
	*/

	double *diag = NULL;
	int i, j, k;
	
	if (M) {
		if (M->nrow != M->ncol) {
			fprintf(stderr, "*** %s:%1d ***  Not a diagonal matrix: %1d != %1d\n",
				__FILE__,  __LINE__,  M->nrow,  M->ncol);
			exit(1);
		}

		if (M->nrow) {
			diag = Calloc(M->nrow, double);
			if (M->A){
				for(i=0; i<M->nrow; i++)
					diag[i] = M->A[ i + j * M->nrow ];
			} else {
				for(k=0; k<M->elems; k++) {
					if (M->i[k] == M->j[k])
						diag[ M->i[k] ] = M->values[k];
				}
			}
		}
	}
	return diag;
}
int inla_matrix_free(inla_matrix_tp *M)
{
	if (M) {
		Free(M->i);
		Free(M->j);
		Free(M->values);
		Free(M->ivalues);
		Free(M->A);
		Free(M->iA);
		Free(M);
	}
	return (0);
}

inla_matrix_tp *inla_matrix_1(int n)
{
	/* 
	   return a 1-matrix with given dimension
	*/

	if (n > 0){
		inla_matrix_tp *M = Calloc(1, inla_matrix_tp);

		M->nrow = M->elems = n;
		M->ncol = 1;
		M->A = Calloc(n, double);

		int i;
		for(i=0; i<n; i++)
			M->A[i] = 1.0;

		return M;
	} else {
		return NULL;
	}
}
	
int inla_file_check(const char *filename, const char *mode)
{
	/* 
	   open file with given mode and return 0 if ok, and 1 if failure
	 */

	FILE *fp = fopen(filename, mode);
	if (fp) {
		fclose(fp);
		return 0;
	} else {
		return 1;
	}
}

#ifdef TESTME
int main(int argc, char **argv)
{
	int i;
	inla_matrix_tp *M;
	
	for(i=1; i< argc; i++) {
		printf("\n\n\n\n *** Check file %s\n\n\n",  argv[i]);
		M = inla_read_fmesher_file(argv[i]);
		printf("\n\n");
		inla_write_fmesher_file(M,  "testmatrix.dat");
		printf("\n\n");
		M = inla_read_fmesher_file("testmatrix.dat");
	}
	return 0;
}
#endif
