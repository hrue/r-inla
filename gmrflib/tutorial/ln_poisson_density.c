// GMRFLib-Tutorial: Copyright (C) 2010 Bradley M. Bell
// Distributed under the terms of the GNU Public License version 2 or higher.

/*!
\file ln_poisson_density.c
\brief Evaluate the log-density function for a Possion distribution.
*/

/// identify who made the last change to this file and when
static const char RCSId[] = "$Id: ln_poisson_density.c,v 1.2 2010/03/21 14:21:30 bradbell Exp $"; 

# include <math.h>
# include "ln_poisson_density.h"
# include "ln_gamma.h"


/*!
\brief
Evaluate the log-density for a Possion distribution.

The Possion distribution has the following density function
with respect to \f$ y \f$:
\f[
	{\bf p} (y | \lambda ) = \frac{ \lambda^y } { y ! } \exp( - \lambda  )
\f]
where \f$ y \f$ is the number of occurnaces,
and \f$ \lambda  \f$ is a parameter that is equal
to the expected number of occurances in the time interval
(during which the occurances are counted).
This distribution has mean 
\f$ \lambda  \f$ and variance \f$ \lambda  \f$.

\param[in] y
is the number of occurances that occured in the time interval
(must be non-negative).

\param[in] \lambda 
is the the expected number of occurances in the time interval.

\return
The return value is \f$ \log [ {\bf p} ( y | \lambda  ) ] \f$.
*/
double ln_poisson_density(int y, double lambda )
{	return y * log(lambda ) - ln_gamma(y + 1.0) - lambda ; }

/*!
Example and test of \ref ln_poisson_density .

\returns
one if the test passes and zero otherwise.
*/
# include <math.h>
int ln_poisson_density_test(void)
{	double lambda = 4.;
	int    y      = 2.;
	// lambda^y * exp(-lambda) / y !
	double poisson = pow(lambda, (double) y) * exp(-lambda) / 2.; 

	int ok = fabs( ln_poisson_density(y, lambda) - log(poisson) ) < 1e-8; 
	return ok;
}
