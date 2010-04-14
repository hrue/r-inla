// GMRFLib-Tutorial: Copyright (C) 2010 Bradley M. Bell
// Distributed under the terms of the GNU Public License version 2 or higher.

/*!
\file multiply_matrix.c
\brief Dense and simple matrix multiplication.
*/

/// identify who made the last change to this file and when
static const char RCSId[] = "$Id: multiply_matrix.c,v 1.6 2010/03/16 13:25:02 bradbell Exp $"; 

// check include file prototype is correct
# include "multiply_matrix.h"


/*!
\brief
A simple implementation of dense (not sparse) matrix multiplication.

We are given two matrices 
\f$ A \in {\bf R}^{\ell \times m} \f$ and \f$ B \in {\bf R}^{m \times n} \f$ 
and which to compute \f$ C \in {\bf R}^{\ell \times n} \f$ given by
\f[
	C_{i,j} = \sum_{k=0}^{m-1} A_{i,k} B_{k,j}
\f]

\param[in] ell
is the row size for \a A and the row size for \a C.

\param[in] m
is the column size for \a A and the row size for \a B.

\param[in] n
is the column size for \a B and the column size for \a C.

\param[in] A
For 
<tt>i=0,...,ell-1, k=0,...,m-1,  A[i + k*ell]</tt> 
is equal to \f$ A_{i,k} \f$.

\param[in] B
For
<tt>k=0,...,m-1, j=0,...,n-1,  B[i + j*m]</tt> 
is equal to \f$ B_{k,j} \f$.

\param[out] C
The input value of the elements of \a C does not matter.
On output,
<tt>i=0,...,ell-1, j=0,...,n-1,  C[i + j*ell]</tt> 
is equal to \f$ C_{i,j} \f$.

\par Assumptions
None of the memory corresponding to the matrix \a C is the same
as memory for the matrices \a A or \a B.
*/
void multiply_matrix(int ell, int m, int n, double* A, double* B, double* C)
{	int i, j, k;

	for (i = 0; i < ell; i++)
	{	for(j = 0; j < n; j++ )
		{	C[i + j*ell] = 0.;
			for(k = 0; k < m; k++)
				C[i + j*ell] += A[i + k*ell] * B[k + j*m];
		}
	}
	return;
}
