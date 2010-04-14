# ifndef GMRF_TUTORIAL_POS_INV_H
# define GMRF_TUTORIAL_POS_INV_H
// $Id: pos_inv_.h,v 1.4 2010/03/16 13:32:34 bradbell Exp $
//
// GMRFLib-Tutorial: Copyright (C) 2010 Bradley M. Bell
// Distributed under the terms of the GNU Public License version 2 or higher.

/*!
\file pos_inv_.h
\brief 
Prototype for routine that computes inverse of a positive definate equations.
*/
 
/*!
Inverse of a symmetric positive definate dense (not sparse) matrix A.

\param[in] n
<tt>*n</tt> is the row and column dimension for the matrix A.
This value is not changed.

\param[in] A
Contains the values for the positive definate matrix \a A. To be specific,
for \f$ i = 0 , \ldots , n-1 \f$
and \f$ j = 0 , \ldots , n-1 \f$,
<tt>A[i + j*n]</tt> is the value of \f$ A_{i,j} \f$.
These values are not changed.

\param work
Is a vector of length <tt>n * n</tt>.
The input values do not matter and the output values
are undefined.

\param[out] Ainv
Is a vector of length <tt>n * n</tt>.
The input values do not matter.
On output,
for \f$ i = 0 , \ldots , n-1 \f$
and \f$ j = 0 , \ldots , n-1 \f$,
<tt>Ainv[i + j*n]</tt> is the value of \f$ A_{i,j}^{-1} \f$.

\param[out] logdet
The input value of \a logdet does not matter.
On output it is the log of the determinant of the matrix \a A.

\param[out] info
\li if <tt>*info == 0</tt>, the problem was solved.

\li if <tt>*info < 0</tt>, a program error occurred in \c pos_inv_.

\li if <tt>*info > 0</tt>, the leading square minor of \a A with size \c info
is not positive definate. 

*/
extern void pos_inv_(
int* n, double* A, double* work, double* Ainv, double *logdet, int* info
);

# endif
