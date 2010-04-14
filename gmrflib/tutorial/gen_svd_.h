# ifndef GMRF_TUTORIAL_GEN_SVD_H
# define GMRF_TUTORIAL_GEN_SVD_H
// $Id: gen_svd_.h,v 1.3 2010/03/16 13:32:34 bradbell Exp $
//
// GMRFLib-Tutorial: Copyright (C) 2010 Bradley M. Bell
// Distributed under the terms of the GNU Public License version 2 or higher.

/*!
\file gen_svd_.h
\brief 
Prototype for singular value decomposition for a general matrix.
*/
 
/*!
Singular value decomposition for a dense (not sparse) general matrix.

The singular value decomposition of \f$ A \in {\bf R}^{m \times n}\f$, 
is a factorization in the form
\f[
	A = U * \Sigma * V^{\rm T}
\f]
where \f$ U \in {\bf R}^{m \times m} \f$ is an orthogonal matrix,
\f$ V \in {\bf R}^{n \times n} \f$ is an orthogonal matrix,
and \f$ \Sigma \in {\bf R}^{m \times n} \f$ is a diagonal matrix; i.e.,
all of the elements of \f$ \Sigma \f$ are zero except for 
\f$ \min ( m , n ) \f$ elements along its diagonal.


\param[in] m
<tt>*m</tt> is the row dimension for the matrix A.
This value is not changed.

\param[in] n
<tt>*n</tt> is the column dimension for the matrix A.
This value is not changed.

\param[in] A
Contains the values for the matrix \a A. To be specific,
for \f$ i = 0 , \ldots , m-1 \f$
and \f$ j = 0 , \ldots , n-1 \f$,
<tt>A[i + j*m]</tt> is the value of \f$ A_{i,j} \f$.
These values are not changed.

\param temp
Is a vector of length <tt>m * n</tt>.
The input values do not matter and the output values
are undefined.

\param work
Is a vector of length <tt>5 * min(m, n) + max(m, n) + 1</tt>.
The input values do not matter and the output values
are undefined.

\param[out] U
Is a vector of length <tt>m * m</tt>.
The input values do not matter.
On output,
for \f$ i = 0 , \ldots , m-1 \f$
and \f$ j = 0 , \ldots , m-1 \f$,
<tt>U[i + j*n]</tt> is the value of \f$ U_{i,j} \f$.

\param[out] S
Is a vector of length <tt>min(m, n)</tt>.
The input values do not matter.
On output,
for \f$ i = 0 , \ldots , min(m,n)-1 \f$,
<tt>S[i]</tt> is the value of \f$ \Sigma_{i,i} \f$.

\param[out] Vt
Is a vector of length <tt>n * n</tt>.
The input values do not matter.
On output,
for \f$ i = 0 , \ldots , n-1 \f$
and \f$ j = 0 , \ldots , n-1 \f$,
<tt>Vt[i + j*n]</tt> is the value of \f$ V_{i,j}^{\rm T} = V_{j,i} \f$.

\param[out] info
\li if <tt>*info == 0</tt>, the problem was solved.

\li if <tt>*info < 0</tt>, a program error occurred in \c gen_svd_.

\li if <tt>*info > 0</tt>, it is the number of superdiagonals that 
did not converge.
*/
extern void gen_svd_(
	int*     m     , 
	int*     n     , 
	double*  A     , 
	double*  temp  , 
	double*  work  , 
	double*  U     , 
	double*  S     ,
	double*  Vt    ,
	int*     info
);

# endif
