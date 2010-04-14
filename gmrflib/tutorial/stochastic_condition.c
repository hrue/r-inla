// GMRFLib-Tutorial: Copyright (C) 2010 Bradley M. Bell
// Distributed under the terms of the GNU Public License version 2 or higher.

/*!
\file stochastic_condition.c
\brief Sampling a GMRF conditioned on a linear stochastic constraint.
*/

/// identify who made the last change to this file and when
static const char RCSId[] = "$Id: stochastic_condition.c,v 1.13 2010/03/17 15:31:47 bradbell Exp $"; 

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
Sampling a GMRF conditioned on a stochastic linear constraint.

\section stochastic_condition_Prior_For_x Prior For x
The \ref unconditional_Prior_For_x
is given as Gaussian with mean and variance given by
\f[
\begin{array}{rcl}
{\bf E} (x) & = & \mu + [ Q + {\rm diag} (c) ]^{-1} b 
\\
{\bf V} ( x , x ) & = & [ Q + {\rm diag} (c) ]^{-1}
\end{array}
\f]

\section stochastic_condition_linear_measurement linear measurement 
Given a matrix \f$ A \in {\bf R}^{m \times n} \f$,
a symmetric positive definate matrix 
\f$ \Sigma \in {\bf R}^{m \times m} \f$,
we model the measurement vector \f$ e \in {\bf R}^m \f$ by
\f[
e = A x + \varepsilon
\f]
where \f$ \varepsilon \sim {\bf N}( 0 , \Sigma ) \f$.
The probability density for the measurement \f$ e \f$,
conditioned on the field \f$ x \f$, is given by
\f[
{\bf p} ( e | x ) \propto \exp \left[ 
- \frac{1}{2} ( A x - e )^{\rm T} \Sigma^{-1} ( A x - e )
\right]
\f]

\section stochastic_condition_conditional_for_x conditional for x
The conditional distribution for \f$ x \f$ is given by
\f[
{\bf p} ( x | e ) \propto {\bf p} (e | x) {\bf p} (x)
\f]
We use \f$ \hat{x} \f$ for the mean of this distribution.
It follows from the fact that the conditional for \f$ x \f$ is Gaussian that 
\f[
\begin{array}{rcl}
0 & = & 
\nabla_x \log [ {\bf p} (e | x ) {\bf p} ( x ) ]_{x = \hat{x}}
\\ & = &
b 
- [ Q + {\rm diag} (c) ] ( \hat{x} - \mu ) 
- A^{\rm T} \Sigma^{-1} (A \hat{x} - e )
\\ & = &
- \left[ 
	Q + {\rm diag} (c) + A^{\rm T} \Sigma^{-1} A \hat{x} 
\right] \hat{x} 
b 
+ [ Q + {\rm diag} (c) ] \mu 
+ A^{\rm T} \Sigma^{-1} e
\\
\hat{x} & = &
\left[ 
	Q + {\rm diag} (c) + A^{\rm T} \Sigma^{-1} A \hat{x} 
\right]^{-1} 
\left[ 
	Q \mu + {\rm diag} (c) \mu + b + A^{\rm T} \Sigma^{-1} e
\right]
\end{array}
\f]
It also follows that the variance of the conditional is given by
\f[
\left[
	\nabla_{xx}^2 \log [ {\bf p} (e | x ) {\bf p} ( x ) ]_{x = \hat{x}}
\right]^{-1}
=
\left[ 
	Q + {\rm diag} (c) + A^{\rm T} \Sigma^{-1} A \hat{x} 
\right]^{-1} 
\f]
Hence, since the conditional distribution is Gaussian,
its negative log-likelihood is given by
\f[
\log [ {\bf p} ( x | e ) ] 
=
- \frac{n}{2} \log( 2 \pi ) 
+ \frac{1}{2} \log \det \left[ 
	Q + {\rm diag} (c) + A^{\rm T} \Sigma^{-1} A \hat{x} 
\right]
- \frac{1}{2} ( x - \hat{x} )^{\rm T} \left[ 
	Q + {\rm diag} (c) + A^{\rm T} \Sigma^{-1} A \hat{x} 
\right] ( x - \hat{x} )
\f]

\return
1 if test passes and 0 otherwise.
*/
int stochastic_condition(void)
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

	// mean for the conditional distribution
	double *x_hat;

	// rhs = b + [Q  + diag(c)] * mu + A^T * Sigma^{-1} * e
	double *rhs;

	// the sample mean for the field
	double *sample_mean;

	// the sample covariance for the field
	double *sample_cov;

	// the information matrix
	double *information;

	// the covariance matrix
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

	// work vector used by pos_inv_
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
	int          nc;
	double      *e;
	double      *A;
	double      *Sigma;

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
	n_calloc = 7 * n + 4 * n * n;
	// allocate the memory
	ptr_calloc = (double*) Calloc(n_calloc, double);
	// split the memory up into individual vectors
	i  = 0;
	x           = ptr_calloc + i; i += n;   // x           has length n
	b           = ptr_calloc + i; i += n;   // b           has length n
	c           = ptr_calloc + i; i += n;   // c           has length n
	mu          = ptr_calloc + i; i += n;   // mu          has length n
	x_hat       = ptr_calloc + i; i += n;   // x_hat       has length n
	rhs         = ptr_calloc + i; i += n;   // rhs         has length n
	sample_mean = ptr_calloc + i; i += n;   // sample_mean has length n
	sample_cov  = ptr_calloc + i; i += n*n; // sample_cov  has length n*n
	information = ptr_calloc + i; i += n*n; // information has length n*n
	covariance  = ptr_calloc + i; i += n*n; // covariance  has length n*n
	work        = ptr_calloc + i; i += n*n; // work        has length n*n
	assert( i == n_calloc );

	// -------------------------------------------------------------------
	// Create the stochastic constraint 
	GMRFLib_make_empty_constr(&constr);
	// # of columns in the matrix A
	nc         = 2;
	constr->nc = nc;  
	// Define the matrix A with first row [0, 1] , second row [1, 0],
	// third row [0, 1], fourth row [1, 0], ...
	A                = Calloc(nc * n, double); 
	constr->a_matrix = A;
	for(i = 0; i < n; i++)
	{	for(k = 0; k < nc; k++)
			A[k + i*nc] = (double) (((i + k) % nc) == 0); 
	}
	// Define the matrix Sigma as diagonal with .1 along the diagonal.
	Sigma                   = Calloc(nc, double);
	constr->errcov_diagonal = Sigma;
	for(k = 0; k < nc; k++)
		Sigma[k] = .25;
	// Define the vector e.  We could choose a true value x_true, 
	// simulate N(0, .1 * I) measurement noise epsilon, then set 
	// 	e = A * x_true + epsilon.
	// but that would complicate this example.
	e                = Calloc(nc, double);
	constr->e_vector = e;
	for(k = 0; k < nc; k++)
		e[k] = 0.;
	// Compute (and store in constr) unspecified internal information that 
	// is needed by the sampling routines.
	{	int scale_constr = 0;
		GMRFLib_prepare_constr(constr, graph, scale_constr);
	}
	// -------------------------------------------------------------------

	// Create stochastic conditional sampling problem
	argQfunc = (void *) graph;     // extra arguments to Qfunc
	fixed    = NULL;               // no components of x are fixed
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


	// initialize the information matrix as [ Q + diag(c) ]
	for(i = 0; i < n; i++)
	{	int nnbs = graph->nnbs[i];
		int *nbs = graph->nbs[i];
		for(j = 0; j < n; j++ )
			information[i + j*n] = 0.;
		information[i +  i*n] = c[i] + Qfunc(i, i, argQfunc);
		for(k = 0; k < nnbs; k++)
		{	j = nbs[k];
			information[i + j*n] = Qfunc(i, j, argQfunc);
		}
	}
	// rhs = b + [Q  + diag(c)] * mu + A^T * Sigma^{-1} * e 
	for(i = 0; i < n; i++)
	{	rhs[i] = b[i];
		for(j = 0; j < n; j++)
			rhs[i] += information[i + j*n] * mu[j];
		for(k = 0; k < nc; k++)
			rhs[i] += A[k + i*nc] * e[k] / Sigma[k];
	}
			
	// add A^T * Sigma^{-1} * A to the information matrix
	for(i = 0; i < n; i++)
	{	for(j = 0; j < n; j++)
		{	for(k = 0; k < nc; k++)
			{	double term = 	A[k + i*nc] * A[k + j*nc] 
				            / Sigma[k];
				information[i + j*n] += term;
			}
		}
	} 
	// compute the covaraince matrix as
	// [ Q + diag(c) + A^T * Sigma^{-1} * A ]^{-1}
	pos_inv_(&n, information, work, covariance, &logdet, &info);
	assert( info == 0 );

	// compute x_hat = 
	// covariance * { b + [Q  + diag(c)] * mu + A^T * Sigma^{-1} * e }
	multiply_matrix(n, n, 1, covariance, rhs, x_hat);

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
			{	double si = problem->sample[i] - x_hat[i];
				double sj = problem->sample[j] - x_hat[j];
				sample_cov[i + j*n] += si * sj / nsim;
			}
		}
		// Example evaluation of the log-likelihood
		if( ell < 10 )
		{	// check the log-density evaluation by GMRFLib
			double delta = 100. * n * DBL_EPSILON;
			double logdens = ln_gauss_density(
				n, problem->sample, x_hat, information, logdet
			);  
			GMRFLib_evaluate(problem);
			ok &= fabs(problem->sub_logdens - logdens) <= delta;
		}
	}
	// check the sample mean
	for(i = 0; i < n; i++)
		ok &= fabs(x_hat[i] - sample_mean[i]) <=  epsilon;


	// check if information matrix and sample covariance are near inverses
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
