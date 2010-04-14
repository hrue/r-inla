//
// GMRFLib-Tutorial: Copyright (C) 2010 Bradley M. Bell
// Distributed under the terms of the GNU Public License version 2 or higher.

# include "ln_gamma.h"

/*!
\file ln_gamma.c
\brief 
Prototype for log of the Gamma function.
*/
/// identify who made the last change to this file and when
static const char RCSId[] = "$Id: ln_gamma.c,v 1.2 2010/03/17 04:49:28 bradbell Exp $";

/*!
Log of the Gamma function.

The Gamma function \f$ \Gamma(a) \f$ is defined by
\f[
	\Gamma(a) = \int_0^\infty x^{a-1} \exp( - x ) \; {\bf d} x
\f]
Note that for \f$ a \f$ a positive integer,
\f$ Gamma(a) = a ! \f$.

\param[in] a
is the point at which we are computing the log of the Gamma function
(\f$ 0 < a \f$).

\return
the return value is \f$ \log [ \Gamma (a) ] \f$.
*/
double ln_gamma(double a)
{	extern double gsl_sf_lngamma(double a);
	return gsl_sf_lngamma(a);
}

/*!
Example and test of \ref ln_gamma .

\returns
one if the test passes and zero otherwise.
*/
# include <math.h>
int ln_gamma_test(void)
{	int  ok = 1;
	double factorial = 1.;
	double i, lg;
	for(i = 1.; i < 5.; i += 1.)
	{	factorial *= i;
		lg  = ln_gamma(i+1.);
		ok &= fabs( lg - log(factorial) ) < 1e-8;
	}
	return ok;
}
