// GMRFLib-Tutorial: Copyright (C) 2010 Bradley M. Bell
// Distributed under the terms of the GNU Public License version 2 or higher.

/*!
\file check_matrix.c
\brief Check that one matrix is close to another matrix.
*/
/// identify who made the last change to this file and when
static const char RCSId[] = "$Id: check_matrix.c,v 1.4 2010/03/16 13:25:01 bradbell Exp $";

/*!
\brief
Check that one matrix is close to another matrix.

We are given two matrices \f$ A \f$ and \f$ B \f$ 
in \f$ {\bf R}^{m \times n} \f$ and a value \f$ \varepsilon \f$.
We say that \f$ A \f$ is \f$ \varepsilon \f$ equal to \f$ B \f$ if
\f[
	\| A - B \|_\infty \leq \varepsilon
\f]
where for any matrix \f$ C \in {\bf R}^{m \times n} \f$,
\f[
\| C \|_\infty = \max \left\{ \left. \sum_{j=1}^n | C_{i,j} | \; \right| \;
	i = 0 , \cdots m-1 \right\}
\f]

\param[in] m
is the row size of the matrices \f$ A \f$ and \f$ B \f$.

\param[in] n
is the column size of the matrices \f$ A \f$ and \f$ B \f$.

\param[in] A
For 
<tt>i=0,...,m-1, j=0,...,n-1,  A[i + j*m]</tt> 
is equal to \f$ A_{i,j} \f$.

\param[in] B
For
<tt>i=0,...,m-1, j=0,...,n-1,  B[i + j*m]</tt> 
is equal to \f$ B_{i,j} \f$.

\param[in] epsilon
is the value for \f$ \varepsilon \f$.

\return
If \f$ A \f$ is \f$ \varepsilon \f$ equal to \f$ B \f$,
the return value is one,
otherwise it is zero.
*/
# include <math.h>

// check include file prototype is correct
# include "check_matrix.h"

int check_matrix(int m, int n, double *A, double *B, double epsilon)
{	int i, j;
	int all_pass = 1;

	for (i = 0; i < m; i++)
	{	double sum_abs_diff = 0.;
		for(j = 0; j < n; j++ )
			sum_abs_diff += fabs( A[i + j*m] - B[i + j*m] );
 		all_pass &= sum_abs_diff <=  epsilon;
	}
	return all_pass;
}
