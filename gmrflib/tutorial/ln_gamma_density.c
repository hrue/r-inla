// GMRFLib-Tutorial: Copyright (C) 2010 Bradley M. Bell
// Distributed under the terms of the GNU Public License version 2 or higher.

/*!
\file ln_gamma_density.c
\brief Evaluate the log-density function for a Gamma distribution.
*/

/// identify who made the last change to this file and when
static const char RCSId[] = "$Id: ln_gamma_density.c,v 1.3 2010/03/21 14:21:29 bradbell Exp $"; 

# include <math.h>
# include "ln_gamma.h"
# include "ln_gamma_density.h"


/*!
\brief
Evaluate the log-density for a Gamma distribution.

The Gamma distribution has the following density function
with respect to \f$ x \f$:
\f[
	{\bf p} ( x | a, b) = \frac{b^a}{ \Gamma (a) } x^{a-1} \exp ( - b x )
\f]
where \f$ a \f$ and \f$ b \f$ are parameters and
\f$ \Gamma (a) \f$ is the \ref ln_gamma "gamma function".
This distribution has mean \f$ a / b \f$ and variance \f$ a / b^2 \f$.

\param[in] x
is the value of the Gamma distributed random variate.

\param[in] a
is the shape parameter in the Gamma distribution

\param[in] b
is the reciprical of the rate parameter in the Gamma distribution
(\f$ 1 / b \f$ is called the scale parameter).

\return
The return value is \f$ \log [ {\bf p} (x | a, b) ] \f$.

*/
double ln_gamma_density(double x, double a, double b)
{	return a*log(b) - ln_gamma(a) + (a-1.0)*log(x) - b*x; }

/*!
Example and test of \ref ln_gamma_density .

\returns
one if the test passes and zero otherwise.
*/
# include <math.h>
int ln_gamma_density_test(void)
{	double a = 4.;
	double b = 3.;
	double x = 2.;
	// x^(a-1) * exp(-b*x) * b^a / (a-1)!
	double gamma = pow(x, a-1) * exp(-b*x) * pow(b, a) / 6.; 

	int ok = fabs( ln_gamma_density(x, a, b) - log(gamma) ) < 1e-8; 
	return ok;
}
