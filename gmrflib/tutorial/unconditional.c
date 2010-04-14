// GMRFLib-Tutorial: Copyright (C) 2010 Bradley M. Bell
// Distributed under the terms of the GNU Public License version 2 or higher.

/*!
\file unconditional.c
\brief Unconditional Sampling of a GMRF
*/

/// identify who made the last change to this file and when
static const char RCSId[] = "$Id: unconditional.c,v 1.14 2010/03/17 15:31:47 bradbell Exp $"; 

#include "GMRFLib/GMRFLib.h"

// evaluate the log-density for a Gaussian
# include "ln_gauss_density.h"

// invert a dense (not sparse) symmetric positive definate matrix 
# include "pos_inv_.h"

// simple (not sparse or fast) matrix multiply
# include "multiply_matrix.h"

// computes entries in the uniform precision matrix example
# include "uniform_precision.h"

// checks if one matrix close to another
# include "check_matrix.h"

/*!
\brief 
Use uniform precision matrix, with \f$ \kappa = 1 \f$, for this example.

\param[in] i
is the row index in the precision matrix.

\param[in] j
is the column index in the precision matrix.

\param[in] arg
The graph for this GMRF is <tt>g = (GMRFLib_graph_tp *) arg</tt>.
See \ref uniform_precision for details.
*/

static double Qfunc(int i, int j, void *arg)
{
	// convert arg to type expected by uniform_precision
	GMRFLib_graph_tp *g = (GMRFLib_graph_tp *) arg;

	double kappa = 1.0;
	return uniform_precision(i, j, kappa, g);
}

/*!
\brief 
Unconditional sampling of a GMRF.

\section unconditional_Prior_For_x Prior For x
Let \f$ Q \in {\bf R}^{n \times n} \f$ be the matrix corresponding
to the function \c Qfunc in this file.
Given a vector \f$ \mu \in {\bf R}^n \f$
a vector \f$ c \in {\bf R}^n \f$, and
a vector \f$ b \in {\bf R}^n \f$,
the prior probability density for the value of the nodes in the graph
\f$ x \in {\bf R}^n \f$ is specified by
\f[
{\bf p} ( x ) \propto \exp \left\{ 
	b^{\rm T} x
	- \frac{1}{2} ( x - \mu )^{\rm T}
		[ Q + {\rm diag} (c) ] 
			( x - \mu ) 
\right\}
\f]
We use \f$ \bar{x} \f$ for the mean corresponding to \f$ {\bf p} (x) \f$
defined above.
It follows that the gradient \f$ \nabla \log {\bf p} ( \bar{x} ) = 0 \f$
and hence
\f[
\bar{x} = {\bf E} (x) = \mu + [ Q + {\rm diag} (c) ]^{-1} b 
\f]
The corresponding variance for \f$ x \f$ is given by 
\f[
	{\bf V} ( x , x )
	=
	[ \nabla^2 \log {\bf p} ( \bar{x} ) ]^{-1}
	=
	[ Q + {\rm diag} (c) ]^{-1}
\f]
The corresponding negative log-likelihood is given by
\f[
\log [ {\bf p} (x) ]
=
- \frac{n}{2} \log( 2 \pi ) 
+ \frac{1}{2} \log \det [ Q + {\rm diag} (c)  ]
- \frac{1}{2} ( x - \bar{x} )^{\rm T} [ Q + {\rm diag} (c)  ] ( x - \bar{x} )
\f]
\return
1 if test passes and 0 otherwise.
*/
int unconditional(void)
{
	// problem parameters  (you can change these)
	const int nsim        = 50000; // number of simulations
	const int seed        = 1234;  // random number seed
	const int nrow        = 4;     // number of rows in lattice
	const int ncol        = 4;     // number of column in lattice
	const int nb_row      = 1;     // 2*nb_row is # neighbouring rows
	const int nb_col      = 1;     // 2*nb_col is # neighbouring colums
	const int cyclic_flag = 0;     // is graph cyclic (on neighbourhoods)

	// a convergence criteria for sum of nsim independent samples
	const double epsilon = 4. / sqrt( (double) nsim );

	// initialize return flag
	int ok = 1;

	// the true mean for the field
	double *x_bar;

	// the sample mean for the field
	double *sample_mean;

	// the sample covariance for the field
	double *sample_cov;

	// the information matrix; i.e., Q + diag(c)
	double *information;

	// the covariance matrix; i.e., [Q + diag(c)]^{-1}
	double *covariance;

	// vectors that appear in the definition of p(x)
	double *x, *b, *c, *mu;

	// pointer to the beginning of memory allocated with Calloc
	double *ptr_calloc;

	// number of doubles allocated
	int   n_calloc;

	// log-determinant returned by pos_inv_
	double logdet;

	// flag returned by pos_inv_
	int info;

	// a work vector 
	double *work;

	// graph defining the connectivity for the random field
	GMRFLib_graph_tp   *graph;

	// structure defining the problem to be solved
	GMRFLib_problem_tp *problem;

	// structure defining the constraints
	GMRFLib_constr_tp   *constr;

	// temporary indices
	int i, j, k, ell;

	// See below for the meaning of these variables 
	char        *fixed;
	void        *argQfunc;
	int          n;
	unsigned int keep;

	// Initialize random number generator: 
	GMRFLib_uniform_init((unsigned long int) seed);

	// create a lattice graph 
	GMRFLib_make_lattice_graph(
		&graph, 
		nrow, 
		ncol, 
		nb_row, 
		nb_col, 
		cyclic_flag
	);
	// number of nodes in the graph
	n = graph->n;

	// number of doubles to allocate
	n_calloc = 6 * n + 4 * n * n;
	// allocate the memory
	ptr_calloc = (double*) Calloc(n_calloc, double);
	// split the memory up into individual vectors
	i  = 0;
	x           = ptr_calloc + i; i += n;   // x           has length n
	b           = ptr_calloc + i; i += n;   // b           has length n
	c           = ptr_calloc + i; i += n;   // c           has length n
	mu          = ptr_calloc + i; i += n;   // mu          has length n
	x_bar       = ptr_calloc + i; i += n;   // x_bar       has length n
	sample_mean = ptr_calloc + i; i += n;   // sample_mean has length n
	sample_cov  = ptr_calloc + i; i += n*n; // sample_cov  has length n*n
	information = ptr_calloc + i; i += n*n; // information has length n*n
	covariance  = ptr_calloc + i; i += n*n; // covariance  has length n*n
	work        = ptr_calloc + i; i += n*n; // work        has length n*n
	assert( i == n_calloc );

	// Create unconditional sampling problem
	argQfunc = (void *) graph;     // extra arguments to Qfunc
	fixed    = NULL;               // no components of x are fixed
	constr   = NULL;               // no constraints are present
	keep     = 0;                  // problem has not yet been created
	for (k = 0; k < n; k++)
	{	x[k]    = 0.0; // values not used becasue fixed == NULL
		b[k]    = 0.5; // value of b 
		c[k]    = 1.0; // value of c
		mu[k]   = 1.5; // value of \mu
	}
	GMRFLib_init_problem(
		&problem, 
		x, 
		b, 
		c, 
		mu, 
		graph, 
		Qfunc, 
		argQfunc, 
		fixed, 
		constr, 
		keep
	);


	// compute information matrix
	for(i = 0; i < n; i++)
	{	int nnbs = graph->nnbs[i];
		int *nbs = graph->nbs[i];
		for(j = 0; j < n; j++ )
			information[i + j*n] = 0.;
		information[i + i*n] = c[i] + Qfunc(i, i, argQfunc);
		for(k = 0; k < nnbs; k++)
		{	j = nbs[k];
			information[i + j*n] = Qfunc(i, j, argQfunc);
		}
	}
	// compute the covariance and log-determinant of information matrix
	pos_inv_(&n, information, work, covariance, &logdet, &info);
	assert( info == 0 );

	// x_bar = mu + [Q + diag(c)]^{-1} b
	multiply_matrix(n, n, 1, covariance, b, x_bar);
	for(i = 0; i < n; i++)
		x_bar[i] += mu[i]; 

	// compute sample mean and covariance
	for(i = 0; i < n; i++)
	{	sample_mean[i] = 0.;
		for(j = 0; j < n; j++)
			sample_cov[i + j*n] = 0.;
	}
	for(ell = 0; ell < nsim; ell++)
	{	GMRFLib_sample(problem);
		for(i = 0; i < n; i++)
		{	sample_mean[i] += problem->sample[i] / nsim;
			for(j = 0; j < n; j++)
			{	double si = problem->sample[i] - x_bar[i];
				double sj = problem->sample[j] - x_bar[j];
				sample_cov[i + j*n] += si * sj / nsim;
			}
		}
		// Example and check of evaluation of the log-likelihood
		if( ell < 10 )
		{	// check the log-density evaluation by GMRFLib
			double delta = 5. * n * DBL_EPSILON;
			double logdens = ln_gauss_density(
				n, problem->sample, x_bar, information, logdet
			);  
			GMRFLib_evaluate(problem);
			ok &= fabs(problem->sub_logdens - logdens) <= delta;
		}
	}
	// check the sample mean
	for(i = 0; i < n; i++)
		ok &= fabs(x_bar[i] - sample_mean[i]) <=  epsilon;

	// check if covariance matrix and sample covariance are near inverses
	ok &= check_matrix(n, n, covariance, sample_cov, n * epsilon);

	/*
	 * Free allocated memory: 
	 */
	Free(ptr_calloc);

	GMRFLib_free_constr(constr);
	GMRFLib_free_problem(problem);
	GMRFLib_free_graph(graph);

	return ok;
}
