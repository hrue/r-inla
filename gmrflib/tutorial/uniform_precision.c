// GMRFLib-Tutorial: Copyright (C) 2010 Bradley M. Bell
// Distributed under the terms of the GNU Public License version 2 or higher.

/*!
\file uniform_precision.c
\brief Example precision matrix where the off diagonal entries are all equal.
*/

/// identify who made the last change to this file and when
static const char RCSId[] = "$Id: uniform_precision.c,v 1.12 2010/03/17 15:31:47 bradbell Exp $"; 

#include "GMRFLib/GMRFLib.h"

// check include file prototype is correct
# include "uniform_precision.h"

/*!
\brief 
Example precision matrix where the off diagonal entries are all equal.

We use \f$ n \f$ to denote the number of nodes in the graph.
For \f$ i = 0 , \ldots , n-1 \f$,
\f$ m(i) \f$ denotes the number of neighbours of node \f$ i \f$.
For \f$ k = 0 , \ldots , m(i) - 1 \f$, 
\f$ J_k \f$ denotes the \f$ k \f$-th neighbour of node \f$ i \f$.
This example precision matrix is given by
\f[
Q_{i,j} = \left\{ \begin{array}{cl}
	\kappa m(i) & {\rm if} \; i = j \\
	- \kappa  & {\rm if} \; j \in \{ J_0 , \ldots , J_{m(i)-1} \} \\
	0   & {\rm otherwise}
\end{array} \right.
\f]

\param[in] i
is the row index in the precision matrix
(\c i and \c j must be graph neighbors or equal).

\param[in] j
is the column index in the precision matrix
(\c i and \c j must be graph neighbors or equal).

\param[in] kappa
is the absolute value of the off diagonal elements.

\param[in] g
The graph for this GMRF is
<tt>g = (GMRFLib_graph_tp *) v</tt>.
\n
\f$ n \f$: The number of nodes in the graph is <tt>g->n</tt>.
\n
\f$ m(i) \f$: The number of neighbours corresponding to node \a i
is <tt>g->nnbs[i]</tt>.
\n
\f$ \{ J_0 , \ldots , J_{m(i)-1} \} \f$:
The set of neighours corresponding to node \a i is
{ <tt>g->nbs[i][0]</tt>, \f$ \cdots \f$ , <tt>g->nbs[i][g->nnbs[i]-1]</tt> }.

\return
The return value is \f$ Q_{i,j} \f$.

\par Assumptions
\li <tt>0 <= i < g->n</tt>
\li <tt>0 <= j < g->n</tt>
\li either <tt>i == j</tt> or \c i and \c j are neighbours
    (not both; i.e., a node is not its own neighbor).
\li Suppose  
	<tt>0 <= ell < g->nnbs[i]</tt>,
	<tt>0 <= k   < g->nnbs[i]</tt>, 
	and <tt>ell != k</tt>.
    It follows that
    <tt>g->nbs[i][ell] != <tt>g->nbs[i][k]</tt></tt>.
*/
 
double uniform_precision(int i, int j, double kappa, GMRFLib_graph_tp* g)
{
	// check assumptions
# ifndef NDEBUG
	static int call_number = 0;
	++call_number;
	// check these assumptions once every 20 calls
	if( call_number % 20 == 1 ) 
	{	int k;
		int neighbours_or_equal;
		int nnbs = g->nnbs[i];
		int* nbs = g->nbs[i];

		neighbours_or_equal = (i == j);
		for(k = 0; k < nnbs; k++) 
		{	int ell;
			assert ( i != nbs[k] );
			for(ell = 0; ell < nnbs; ell++)
				assert( ell == k || nbs[ell] != nbs[k] ); 
			neighbours_or_equal |= (j == nbs[k]);
		}
		assert( neighbours_or_equal );
	}
	assert( 0 <= i && i < g->n );
	assert( 0 <= j && j < g->n );
# endif

	if (i != j) {
		// each off diagonal case occurs once
		return - kappa;
	} else {
		// diagonal case (number of neighbouring points)
		return kappa * g->nnbs[i];
	}
}
