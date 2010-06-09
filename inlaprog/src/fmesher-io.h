/* fmesher-io.h
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
#ifndef __FMESHER_IO_H__
#define __FMESHER_IO_H__
#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS					       /* empty */
# define __END_DECLS					       /* empty */
#endif
__BEGIN_DECLS


/* 
   this follows std (i,j,values) storage. in case of dense matrix, then the values are stored in A, columnwise. 
 */
typedef struct
{
	int nrow;
	int ncol;

	/* 
	   sparse
	*/
	int elems;					       /* number of elements */
	int *i;
	int *j;

	/* 
	   on reading; 'values' are always set. on writing; only one of these can be set. 
	*/
	int *ivalues;
	double *values;

	/* 
	   dense, columnwise. on reading `A' is always set. on writing; only one of these can be set. ls
	*/
	int *iA;
	double *A;
}
	inla_matrix_tp;


double *inla_matrix_get_diagonal(inla_matrix_tp *M);
inla_matrix_tp *inla_read_fmesher_file(const char *filename);
int inla_matrix_free(inla_matrix_tp *M);
int inla_write_fmesher_file(inla_matrix_tp *M, const char *filename);
int inla_file_check(const char *filename, const char *mode);
inla_matrix_tp *inla_matrix_1(int n);

__END_DECLS

#endif
