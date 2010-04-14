// GMRFLib-Tutorial: Copyright (C) 2010 Bradley M. Bell
// Distributed under the terms of the GNU Public License version 2 or higher.

/*!
\file ln_gauss_density.c
\brief Evaluate the log-density function for a Gaussian.
*/

/// identify who made the last change to this file and when
static const char RCSId[] = "$Id: ln_gauss_density.c,v 1.3 2010/03/21 14:21:29 bradbell Exp $"; 

# include <math.h>
# include "ln_gauss_density.h"

// routine to invert a symmetric positive definate matrix
# include "pos_inv_.h"

/*!
\brief
Evaluate the log-density function for a Gaussian distribution.

Given a mean \f$ \mu \in {\bf R}^n \f$, 
a variance \f$ V \in {\bf R}^{n \times n} \f$,
and a sample \f$ x \in {\bf R}^n \f$,
evaluate the log-density for the corresponding Gaussian distribution; 
i.e.
\f[
	\log [ {\bf p} ( x | \mu , V ) ] 
	=
	- \frac{n}{2} \log( 2 \pi )
	- \frac{1}{2} \log \det ( V ) 
	- \frac{1}{2} ( x - \mu )^{\rm T} V^{-1} ( x - \mu )
\f]

\param[in] n
is the number of components in the random variable.

\param[in] x
is a vector of length \a n specifying the
point at which to evaluate the log-density.

\param[in] mu
is a vector of length \a n containing the mean
for the distribution; i.e., \f$ \mu \f$. 

\param[in] inv_V
is a vector of length <tt>n * n</tt> containing
the inverse of the varaince for the distribution; i.e. \f$ V^{-1} \f$.
To be specific,
<tt>i=0,...,n-1, j=0,...,n-1,  inv_V[i + j*m]</tt> 
is equal to \f$ V_{i,j}^{-1} \f$.

\param[in] logdet
is the log-determinant for \f$ V^{-1} \f$
(which is the negative of the log-determinant for \f$ V \f$. 

\return
is the log-density at the point \a x,
for the  Gaussian distribution with mean \a mu  and variance \a V; i.e.,
\f$ \log [ {\bf p} ( x | \mu , V ) ]  \f$.
*/

double ln_gauss_density(
	int n, double *x, double *mu, double *inv_V, double logdet
)
{	double log_2pi;
	double sum;
	int i, j;

	sum = 0.;
	for(i = 0; i < n; i++)
		for(j = 0; j < n; j++)
			sum -= (x[i]-mu[i]) * inv_V[i + j*n] * (x[j]-mu[j]); 

	sum    +=  logdet;
	log_2pi = log( 8. * atan(1) );
	sum    -= n * log_2pi;

	return sum / 2.;
}
