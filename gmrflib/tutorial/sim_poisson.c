// GMRFLib-Tutorial: Copyright (C) 2010 Bradley M. Bell
// Distributed under the terms of the GNU Public License version 2 or higher.

/*!
\file sim_poisson.c
\brief Generate one realization from a Possion distribution
*/

/// identify who made the last change to this file and when
static const char RCSId[] = "$Id: sim_poisson.c,v 1.2 2010/03/19 14:53:01 bradbell Exp $"; 

# include <math.h>
# include <assert.h>
# include <gsl/gsl_rng.h>
# include <gsl/gsl_randist.h>
# include "ln_poisson_density.h"
# include "ln_gamma.h"

/// A GSL random number generator used by \ref sim_poisson
static gsl_rng *rng_;

/// Initialize \ref sim_poisson 
void sim_poisson_init(void)
{	// just pick on generator and do not worry about seed
	rng_  = gsl_rng_alloc(gsl_rng_taus);
	assert( rng_ != NULL );	
}

/// Free memory allocated for \ref sim_poisson
void sim_poisson_free(void)
{	gsl_rng_free(rng_); }


/*!
Generate one realization from a Possion distribution

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

\param[in] \lambda 
is the the expected number of occurances in the time interval.

\return
is the simulated number of occurances that occured in the time interval
(must be non-negative).

\par Assumptions
- The function \ref sim_poisson_init must be called before the
first call to \ref sim_poisson.

- The function \c sim_poisson_free should be called after the last
call to \c sim_poisson.
*/
unsigned int sim_poisson(double lambda )
{	return gsl_ran_poisson(rng_, lambda); }

/*!
Example and test of \ref sim_poisson.

\returns
one if the test passes and zero otherwise.
*/
# include <math.h>
int sim_poisson_test(void)
{	double  lambda  = 4.;
	int     n_sim   = 50000;

	double sum, sumsq, sample_mean, sample_var;
	int i, y, ok = 1;

	sim_poisson_init();
	sum   = 0.;
	sumsq = 0.;
	for(i = 0; i < n_sim; i++)
	{	y      = sim_poisson(lambda); 
		sum   += (double) y;
		sumsq += (double) (y * y);
	} 

	sample_mean = sum / (double) n_sim;
	sample_var  = sumsq / (double) n_sim - sample_mean * sample_mean;

	ok &= fabs( sample_mean - lambda ) < 10. / sqrt( (double) n_sim ); 
	ok &= fabs( sample_var  - lambda ) < 10. / sqrt( (double) n_sim );

	return ok;
}
