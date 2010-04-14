// GMRFLib-Tutorial: Copyright (C) 2010 Bradley M. Bell
// Distributed under the terms of the GNU Public License version 2 or higher.

/*!
\file block_update.c
\brief MCMC sampling of a GMRF and parameters in its conditional distribution.
*/

/// identify who made the last change to this file and when
static const char RCSId[] = "$Id: block_update.c,v 1.6 2010/03/21 14:33:46 bradbell Exp $"; 

# include <assert.h>
# include "GMRFLib/GMRFLib.h"

# include "ln_poisson_density.h"
# include "uniform_precision.h"
# include "sim_poisson.h"

/*!
Unormalized log-likelihood of the data given the field.

\param[out] loglik
a vector of lenght \c m containing the unormalized posterior log likelihood 
for the field at index \c idx and at the field values specified by \c x.

\param[in] x
a vector of length \c m containing values for the field at index \c idx.

\param[in] m
is the length of the vectors \c loglik and \c x
(must be zero, one, or three).

\param[in] idx
the field index (also the graph index) at which to evaluate
the posterior log likelihood.

\param[in] x_vec
not used.

\param[in] arg (see below)

\return
The value zero, indicating no error, is returned. 
Note that \c assert is used to check the correctness of the input arguments.

\par arg

- ((void **) arg )[0] 
is a <tt>const GMRFLib_graph_tp *</tt> pointing to the graph
(see \ref uniform_precision for details).

- ((void **) arg )[1] 
is a <tt>const double *</tt> pointing to a vector containing the 
measured number of counts for each element of the graph.

- ((void **) arg)[2]
is a <tt>const double *</tt> pointing to the value of \f$ \kappa \f$
to use for precision matrix of the GMRF.
*/
static int data_loglikelihood(
	double *loglik         ,
	double *x              , 
	int     m              , 
	int     idx            , 	
	double *x_vec          , 
	void   *arg            )
{	int i;

	GMRFLib_graph_tp *g = (GMRFLib_graph_tp *) ((void **) arg)[0];
	double *y           = (double *) ((void **) arg)[1];

	assert( 0 <= idx && (size_t) idx < (g->n) ); 
	assert( 0 == m || m == 1 || m == 3 );
	
	for(i = 0; i < m; i++)
		loglik[i] = ln_poisson_density(y[idx], exp( x[i] ) );

	// return no error
	return 0;
}

/*!
\brief 
Use uniform precision matrix where \f$ \kappa \f$ is an argument.

\param[in] i
is the row index in the precision matrix
(\c i and \c j must be graph neighbors or equal).

\param[in] j
is the column index in the precision matrix
(\c i and \c j must be graph neighbors or equal).

\param[in] arg (see below)

\return
The return value is the precision matrix 
(inverse covaraince matrix) value at index \f$ (i, j) \f$.

\par arg

- ((void **) arg )[0] 
is a <tt>const GMRFLib_graph_tp *</tt> pointing to the graph
(see \ref uniform_precision for details).

- ((void **) arg )[1] 
is a <tt>const double *</tt> pointing to a vector containing the 
measured number of counts for each element of the graph.

- ((void **) arg)[2]
is a <tt>const double *</tt> pointing to the value of \f$ \kappa \f$
to use for precision matrix of the GMRF.
*/

static double Qfunc(int i, int j, void *arg)
{
	GMRFLib_graph_tp *g = (GMRFLib_graph_tp *) ((void **) arg)[0];
	double *kappa       = (double *) ((void **) arg)[2];

	return uniform_precision(i, j, *kappa, g);
}

/*!
MCMC block updating to approximate the posterior for a GMRF and its parameters.

\par Graph
The following one dimensional graph is used for this example:
\verbatim
	(0) - (1) - ... - (n-1)
\endverbatim
For \f$ i = 0 , \ldots , n-1 \f$,
we use \f$ x_i \f$ to denote the value of the field at the corresponding node.

\par Precision Matrix
The precision matrix for this example is given by
\f[
Q_{i,j} ( \kappa ) 
= \left\{ \begin{array}{ll}
	\kappa   & {\rm if} \; i = j = 0 \; {\rm or} \; i = j = n-1 \\
	2 \kappa & {\rm if} \; i = j \; {\rm and} \; 0 < i < n-1  \\
	- \kappa & {\rm if} \; | i - j | = 1  \\
	0        * {\rm otherwise}
\end{array} \right.
\f]
where \f$ \kappa > 0 \f$ is a fixed parameter.

\par GMRF Distribution
The GMRF \f$ x \f$ is distributed with mean \f$ \log( \lambda ) \f$
and varaince \f$ [ Q ( \kappa ) + {\rm I} ]^{-1} \f$ 
where \f$ Q ( \kappa ) \f$ is the precision matrix defined above
and \f$ {\rm I} \f$ is the identity matrix.
In otherwords, 
\f[
\log [ {\bf p} ( x | \kappa ) ]
= 
- \frac{n}{2} \log( 2 \pi ) 
+ \frac{1}{2} \log \det [ Q ( \kappa ) + {\rm diag} (c)  ]
- \frac{1}{2} ( x - \mu )^{\rm T} 
	[ Q ( \kappa) + {\rm diag} (c)  ] 
		( x - \mu )
\f]
where for \f$ i = 0 , \ldots , n-1 \f$,
\f$ c_i = 1 \f$ and \f$ \mu_i = \log( \lambda ) \f$


\par Measurement Distribution
Let \f$ n \f$ be the number of components in the field (and the graph).
Given the value for the field \f$ x \in {\bf R}^n \f$
(and the field parameter \f$ \kappa \f$),
for \f$ i = 0 , \ldots , n-1 \f$
the measurement \f$ y_i \f$ has a Poission distribution with
rate \f$ \exp ( x_i ) \f$.
In otherwords,
\f[
{\bf p} ( y_i | x , \kappa ) 
= 
\frac{ \exp( y_i x_i ) } { y_i ! } \exp[ - \exp( x_i ) ]
\f]
Note that the \f$ \exp( y_i x_i ) \f$ 
is equal to \f$ \exp( x_i ) \f$ to the power \f$ y_i \f$.

\par Parameter Prior Distribution
The parameter \f$ \lambda \f$ is known during the MCMC simulation
and hence no prior for it is necessary.
The parameter \f$ \kappa \f$ is unknown. 
Hence we need to specify a prior and proposal distribution 
for \f$ \kappa \f$
The prior distribution for \f$ \kappa \f$ is 
\f[
{\bf p} ( \kappa ) = \left\{ \begin{array}{ll}
1 & {\rm if} \; \kappa \in [ \kappa_{\rm sim} / 10 \; , \; 10 \kappa_{\rm sim} ]
\\
0 & {\rm otherwise}
\end{array} \right.
\f]
where \f$ \kappa_{\rm sim} \f$ 
is the value used for \f$ \kappa \f$ during the simualtions.

\par Parameter Proposal Distribution
The proposal density for a transition from \f$ \kappa_{\rm old} \f$
to \f$ \kappa_{\rm new} \f$ is the same as the prior; i.e.,
\f[
{\bf p} ( \kappa_{\rm new} | \kappa_{\rm old} ) 
=
\left\{ \begin{array}{ll}
1 & {\rm if} \; \kappa \in [ \kappa_{\rm sim} / 10 \; , \; 10 \kappa_{\rm sim} ]
\\
0 & {\rm otherwise}
\end{array} \right.
\f]

*/
int block_update(void)
{
	// problem parameters (you can change these) -------------------------
	const int n_sim         = 100;   // number of MCMC samples
	const int n_row         = 1;     // number of rows in the lattice
	const int n_col         = 100;   // number of columns in the lattice
	const int seed          = 123;	 // the seed
	const int nb_row        = 0;     // 2*nb_row is # neighbouring rows
	const int nb_col        = 1;     // 2*nb_col is # neighbouring columns
	const int cyclic_flag   = 0;     // is the graph cyclic

	const double kappa_sim  = 10;    // value of kappa during simulations
	const double lambda_sim = 20.;   // rate during data simulation
	// --------------------------------------------------------------------
	// define other simple variables in alphabetical order
	double              accept;
	void              **arg,        **arg_new;
	double             *b;
	double             *c;
	double             *d;
	double              diff;
	GMRFLib_constr_tp  *constr;
	char               *fixed;
	unsigned int        i;
	unsigned int        k;
	unsigned int        keep;
	double              kappa,       kappa_new;
	double              laccept;
	double             *mu;
	int                 n;
	unsigned int        n_calloc;
	int                 ok;
	double             *ptr_calloc;
	double              sample_mean;
	double              sample_var;
	double              sum;
	double              sumsq;
	double             *x,         *x_new;
	double             *y;
	// --------------------------------------------------------------------

	// options for the optimizer
	GMRFLib_optimize_param_tp *optimize_param;

	// some specifications for the block-sampler
	GMRFLib_blockupdate_param_tp *blockupdate_param;

	// do not need to Initialize random number generator, it does so automatically.
	GMRFLib_uniform_init((unsigned long int) seed);

	// structure defining the problem used to simulate the data
	GMRFLib_problem_tp *problem;

	// graph defining the connectivity for the random field
	GMRFLib_graph_tp   *graph;

	// create the graph ----------------------------------------------------
	GMRFLib_make_lattice_graph(
		&graph,
		n_row,
		n_col,
		nb_row,
		nb_col,
		cyclic_flag
	);
	// number of nodes in the graph
	n = graph->n;
	// ---------------------------------------------------------------------
	// Allocate double memory
	n_calloc = 6 * n;
	ptr_calloc = (double *) Calloc(n_calloc, double);
	i          = 0;
	c          = ptr_calloc + i; i += n;  // c[n]
	b          = NULL;
	d          = ptr_calloc + i; i += n;  // d[n]
	mu         = ptr_calloc + i; i += n;  // mu[n]
	x          = ptr_calloc + i; i += n;  // x[n]
	x_new      = ptr_calloc + i; i += n;  // x_new[n]
	y          = ptr_calloc + i; i += n;  // x_new[n]
	assert( i == n_calloc );
	// Allocate other memory
	arg      = (void **) Calloc(3, void*);
	arg_new  = (void **) Calloc(3, void*);
	// ----------------------------------------------------------------------
	// value of extra argument vector 
	arg[0]     = (void *) graph;
	arg[1]     = (void *) y;
	arg[2]     = (void *) &kappa;
	arg_new[0] = (void *) graph;
	arg_new[1] = (void *) y;
	arg_new[2] = (void *) &kappa_new;

	// Sample from the prior distribution for the GMRF x
	fixed    = NULL;               // no components of x are fixed
	constr   = NULL;               // no constraints are present
	keep     = 0;                     // problem has not yet been created
	for (i = 0; i < n; i++)
	{	x[i]    = 0.0;               // values not used becasue fixed == NULL
		c[i]    = 1.0e-8;	     // add something small to Q to make it proper
		d[i]    = 1.0;               // multiplier for loglik
		mu[i]   = log( lambda_sim ); // so that exp( mean of x ) = lambda_sim
	}
	kappa = kappa_sim;
	GMRFLib_init_problem(
		&problem, 
		x, 
		b, 
		c, 
		mu, 
		graph, 
		Qfunc, 
		(void *) arg, 
		fixed, 
		constr, 
		keep
	);
	GMRFLib_sample(problem);

	// recenter the sample so its sums to zero (could be done using `constr' above)
	sum = 0.0;
	for(i=0; i<n; i++)
		sum += problem->sample[i];
	sum /= n;
	for(i=0; i<n; i++)
		problem->sample[i] -=  sum - log(lambda_sim);

	// now simulate data for each point on the graph
	sim_poisson_init();
	for(i = 0; i < n; i++){
		//This is wrong, the sample is not 'x' but problem->sample
		y[i] = sim_poisson( exp( problem->sample[i] ) );
		//printf("%d %g %d\n", i, problem->sample[i],  (int) y[i]);
	}
	sim_poisson_free();

	// options for the blockupdate-routine and the optimizing-routine 
	GMRFLib_default_blockupdate_param(&blockupdate_param);
	GMRFLib_default_optimize_param(&optimize_param);

	// initial valuse for the MCMC
	kappa = kappa_sim;        // This initialization is overly accurate.
	for(i = 0; i < n; i++)    // Once this works, we plan to try a less 
		x[i] = mu[i];        // accurate initialization.

	// now run the Markov Chain Monte-Carlo
	sum   = 0.;
	sumsq = 0.;
	for(k = 0; k < n_sim; k++)
	{
		// next value for kappa, there is a builtin functin for this
		kappa_new = kappa_sim * GMRFLib_scale_proposal(2.0);

		// conditonial on kappa_new, propose a new value for the GMRF x
		GMRFLib_blockupdate(
			&laccept              ,
			x_new                 ,
			x                     ,
			b                     ,
			b                     ,
			c                     ,
			c                     ,
			mu                    ,
			mu                    ,
			d                     ,
			d                     ,
			data_loglikelihood    ,
			(void *) arg          ,
			data_loglikelihood    ,
			(void *) arg          ,
			fixed                 ,
			graph                 ,
			Qfunc                 ,
			(void *) arg_new      ,
			Qfunc                 ,
			(void *) arg          ,
			NULL                  ,
			NULL                  ,
			NULL                  ,
			NULL                  ,
			constr                ,
			constr                ,
			optimize_param        ,
			blockupdate_param
			);
		// We need to correct for the normalizing constant in
		//	p( x_new | kappa_new ) / p( x | kappa )
		// which is equal to the determinant ratio
		//	det( Q(kappa_new) ) / det( Q(kappa) )
		laccept += log(kappa_new / kappa) * (n - 1) / 2.0; 
		// Our MCMC proposal density is such that
		//	p( kappa_new | kappa ) / p( kappa | kappa_new ) = 1.
		accept   = exp(DMIN(0.0, laccept));
		//printf("accept %g\n",  accept);
		// accpet new or keep old (note unifrom is almost everywhere < 1)
		if( GMRFLib_uniform() < accept || k == 0)
		{	kappa = kappa_new;
			for(i = 0; i < n; i++)
				x[i] = x_new[i];
		}		
		sum   += kappa;
		sumsq += kappa * kappa;
	}
	sample_mean = sum / (double) n_sim;
	sample_var  = sumsq / (double) n_sim - sample_mean * sample_mean;

	diff = kappa_sim - sample_mean;
	ok   = 1;
	ok  &= fabs(diff) < kappa_sim/sqrt(n_sim);
	ok  &= sample_var >=  diff * diff / sqrt(n_sim);

	if (ok)
		printf("OK:\tblock_update\n");
	else
		printf(
			"block_update() is not yet working. Here are some of the details:\n"
			"\tThe simulation value for kappa is %f.\n"
			"\tThe block_update sample mean and variance for kappa are:\n"
			"\tsample_mean = %f, sample_var = %f.\n", 
			kappa_sim, sample_mean, sample_var
			);

	// Free allocated memory
	Free(arg_new);
	Free(arg);
	Free(ptr_calloc);

	GMRFLib_free_constr(constr);
	GMRFLib_free_problem(problem);
	GMRFLib_free_graph(graph);

	return ok;
}
