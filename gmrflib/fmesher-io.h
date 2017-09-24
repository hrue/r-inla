
/* fmesher-io.h
 * 
 * Copyright (C) 2010-2011 Havard Rue
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
#ifndef __FMESHER_IO_H__
#define __FMESHER_IO_H__
#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS					       /* empty */
#define __END_DECLS					       /* empty */
#endif
__BEGIN_DECLS

/* 
   this follows std (i,j,values) storage. in case of dense matrix, then the values are stored in A, columnwise. 
 */
    typedef struct {
	int nrow;
	int ncol;

	/*
	 * sparse 
	 */
	int elems;					       /* number of elements */
	int *i;
	int *j;


	/*
	 * these are only defined if the matrix is sparse 
	 */
	int htable_column_order;			       /* stored by column or row? default is row which is = 0 */
	GMRFLib_graph_tp *graph;			       /* (possibly nonsymmetric) graph */
	map_id **htable;				       /* array of hash-tables for the values */

	/*
	 * on reading; 'values' are always set. on writing; only one of these can be set. 
	 */
	int *ivalues;
	double *values;

	/*
	 * dense, columnwise. on reading `A' is always set. on writing; only one of these can be set. 
	 */
	int *iA;
	double *A;

	/*
	 * if read from file, add the fileinfo here 
	 */
	char *filename;					       /* filename if any, where this file is read from */
	long int offset;				       /* offset in the file */
	int whence;					       /* whence of the file */
	long int tell;					       /* the position where this matrix ended */
} GMRFLib_matrix_tp;


GMRFLib_matrix_tp *GMRFLib_matrix_1(int n);
GMRFLib_matrix_tp *GMRFLib_read_fmesher_file(const char *filename, long int offset, int whence);
GMRFLib_matrix_tp *GMRFLib_matrix_transpose(GMRFLib_matrix_tp * M);
double *GMRFLib_matrix_get_diagonal(GMRFLib_matrix_tp * M);
double GMRFLib_matrix_get(int i, int j, GMRFLib_matrix_tp * M);
int GMRFLib_file_exists(const char *filename, const char *mode);
int GMRFLib_is_fmesher_file(const char *filename, long int offset, int whence);
int GMRFLib_matrix_free(GMRFLib_matrix_tp * M);
int GMRFLib_write_fmesher_file(GMRFLib_matrix_tp * M, const char *filename, long int offset, int whence);
int GMRFLib_matrix_add_graph_and_hash(GMRFLib_matrix_tp * M);
int GMRFLib_matrix_get_row(double *values, int i, GMRFLib_matrix_tp * M);

__END_DECLS
#endif
