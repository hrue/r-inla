// GMRFLib-Tutorial: Copyright (C) 2010 Bradley M. Bell
// Distributed under the terms of the GNU Public License version 2 or higher.

/*!
\file fixed_condition.c
\brief Sampling a GMRF conditioned on a fixed subset of the field.
*/

/// identify who made the last change to this file and when
static const char RCSId[] = "$Id: fixed_condition.c,v 1.9 2010/03/17 15:31:47 bradbell Exp $"; 

#include "GMRFLib/GMRFLib.h"

// evaluate the log-density for a Gaussian
# include "ln_gauss_density.h"

// invert a dense (not sparse) symmetric positive definate matrix 
# include "pos_inv_.h"

// simple (not sparse or fast) matrix multiply
# include "multiply_matrix.h"

// computes entries in the uniform precision matrix example
# include "uniform_precision.h"

// checks if one matrix is close to another
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
Sampling a GMRF conditioned on a fixed subset of the field.

\section fixed_condition_Prior_For_x Prior For x
The \ref unconditional_Prior_For_x
is given as Gaussian with mean and variance given by
\f[
\begin{array}{rcl}
{\bf E} (x) & = & \mu + [ Q + {\rm diag} (c) ]^{-1} b 
\\
{\bf V} ( x , x ) & = & [ Q + {\rm diag} (c) ]^{-1}
\end{array}
\f]

\section Conditional_Gaussian_Theorem Conditional Gaussian Theorem
Suppose that \f$ u \in {\bf R}^m \f$ and \f$ w \in {\bf R}^\ell \f$ 
are Gaussian random varables and define
\f[
w_\perp = w - {\bf V} (w , u) {\bf V} (u, u)^{-1} [ u - {\bf E} (u) ]
\f]
It follows that \f$ w_\perp \f$ is independent of \f$ u \f$ and 
\f[
\begin{array}{rcl}
{\bf E} ( w | u ) & = & 
{\bf E} ( w ) + {\bf V} (w , u) {\bf V} (u, u)^{-1} [ u - {\bf E} (u) ]
\\
{\bf V} ( w , w | u ) 
& = &
{\bf V} ( w , w ) 
-
{\bf V} (w , u) {\bf V} (u, u)^{-1} {\bf V} (u , w ) 
\\
& = & 
{\bf V} ( w_\perp , w_\perp )
\end{array}
\f]
\b Proof:
\f[
\begin{array}{rcl}
{\bf V} ( w_\perp , u )
& = &
{\bf V} ( w , w ) 
-
{\bf V} (w , u) {\bf V} (u, u)^{-1} 
	{\bf E} \{ [ u - {\bf E} (u) ] [ u - {\bf E} (u) ]^{\rm T} \}
\\
& = &
{\bf V} ( w , w ) 
-
{\bf V} (w , u) {\bf V} (u, u)^{-1} {\bf V} (u, u)
\\
& = &
0
\end{array}
\f]
This establishes the first assertion in the theorem; i.e.,
that \f$ w_\perp \f$ and \f$ u \f$ are independent.
We also have
\f[
\begin{array}{rcl}
{\bf E} ( w | u )
& = &
{\bf E} ( w_\perp | u ) 
+ 
{\bf V} (w , u) {\bf V} (u, u)^{-1} [ u - {\bf E} (u) ]
\\
& = &
{\bf E} ( w_\perp ) 
+ 
{\bf V} (w , u) {\bf V} (u, u)^{-1} [ u - {\bf E} (u) ]
\\
& = &
{\bf E} ( w ) 
+ 
{\bf V} (w , u) {\bf V} (u, u)^{-1} [ u - {\bf E} (u) ]
\end{array}
\f]
This establishes the equation for \f$ {\bf E} ( w | u ) \f$ in the theorem.
We also have
\f[
\begin{array}{rcl}
{\bf V} ( w , w | u )
& = &
{\bf V} ( w_\perp , w_\perp | u )
\\
& = &
{\bf V} ( w_\perp, w_\perp )
\\
& = &
{\bf V} ( w , w ) 
- 
2 * {\bf V} ( w , u ) {\bf V} (u, u)^{-1}  {\bf V} ( u , w )
+
{\bf V} (w , u) {\bf V} (u, u)^{-1} 
	{\bf V} (u, u) 
	{\bf V} (u, u)^{-1}  {\bf V} ( u , w )
\\
& = &
{\bf V} ( w , w ) 
- 
{\bf V} ( w , u ) {\bf V} (u, u)^{-1}  {\bf V} ( u , w )
\end{array}
\f]
The establishes the equation for 
\f$ {\bf V} ( w , w | u ) \f$
and thereby completes the proof of the theorem.

\return
1 if test passes and 0 otherwise.
*/
int fixed_condition(void)
{
	// problem parameters  (you can change these)
	const int nsim        = 50000; // number of simulations
	const int seed        = 1234;  // random number seed
	const int nrow        = 4;     // number of rows in lattice
	const int ncol        = 4;     // number of column in lattice
	const int nb_row      = 1;     // 2*nb_row is # neighbouring rows
	const int nb_col      = 1;     // 2*nb_col is # neighbouring colums
	const int cyclic_flag = 0;     // is graph cyclic (on neighbourhoods)
	const int n_fixed     = nrow * ncol / 2; // number fixed components

	// a convergence criteria for sum of nsim independent samples
	const double epsilon = 4. / sqrt( (double) nsim );

	// initialize return flag
	int ok = 1;

	// the sample mean for the field
	double *sample_mean;

	// the sample covariance for the field
	double *sample_cov;

	// temporary matrices
	double *work, *temp;

	// the corresponding variance matrices
	double *V_xx, *V_uu, *V_uw, *V_wu, *V_ww, *inv_V_ww, *inv_V_uu;

	// x_bar = mu + [Q  + diag(c)]^{-1} * b = mu + V_xx * b
	double *x_bar;

	// expecteced value of w given u
	double *w_hat;

	// covaraince of w given u
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

	// graph defining the connectivity for the random field
	GMRFLib_graph_tp   *graph;

	// structure defining the problem to be solved
	GMRFLib_problem_tp *problem;

	// structure defining the constraints
	GMRFLib_constr_tp   *constr;

	// number of nodes in the graph; i.e., number of components in x
	int          n;

	// number of fixed nodes; i.e., number of components in u
	int          nu;

	// number of free nodes; i.e., number of components in w
	int          nw;

	// temporary indices
	int i, j, k, ell;

	// See below for the meaning of these variables 
	char        *fixed;
	void        *argQfunc;
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
	// number of fixed components; i.e., size of u
	nu = n_fixed;
	assert( nu < n );
	// number of free components; i.e., size of w
	nw = n - nu;

	// number of doubles to allocate
	n_calloc = 5*n + 2*nw + 3*n*n + 4*nw*nw + 2*nw*nu + 2*nu*nu;
	// allocate the memory
	ptr_calloc = (double*) Calloc(n_calloc, double);
	// split the memory up into individual vectors
	i  = 0;
	x           = ptr_calloc + i; i += n;     // x           length n
	b           = ptr_calloc + i; i += n;     // b           length n
	c           = ptr_calloc + i; i += n;     // c           length n
	mu          = ptr_calloc + i; i += n;     // mu          length n
	x_bar       = ptr_calloc + i; i += n;     // x_bar       length n
	sample_mean = ptr_calloc + i; i += nw;    // sample_mean length nw
	w_hat       = ptr_calloc + i; i += nw;    // w_hat       length nw
	V_xx        = ptr_calloc + i; i += n*n;   // V_xx        length n*n
	work        = ptr_calloc + i; i += n*n;   // work        lenght n*n
	temp        = ptr_calloc + i; i += n*n;   // temp        lenght n*n
	sample_cov  = ptr_calloc + i; i += nw*nw; // sample_cov  length nw*nw
	covariance  = ptr_calloc + i; i += nw*nw; // covariance  length nw*nw
	V_ww        = ptr_calloc + i; i += nw*nw; // V_ww        length nw*nw
	inv_V_ww    = ptr_calloc + i; i += nw*nw; // inv_V_ww    length nw*nw
	V_wu        = ptr_calloc + i; i += nw*nu; // V_wu        length nw*nu
	V_uw        = ptr_calloc + i; i += nu*nw; // V_uw        length nu*nw
	V_uu        = ptr_calloc + i; i += nu*nu; // V_uu        length nu*nu
	inv_V_uu    = ptr_calloc + i; i += nu*nu; // inv_V_uu    length nu*nu
	assert( i == n_calloc );

	// Specify the fixed subset for this problem;
	fixed = Calloc(n, char);     // allocate memory for this flag
	for(i = 0; i < nu; i++)     // indices in x corresponding to u
		fixed[i] = 1;
	for(i = nu; i < n; i++)     // indices in x corresponding to w
		fixed[i] = 0;
	

	// Create fixed conditional sampling problem
	argQfunc = (void *) graph;     // extra arguments to Qfunc
	constr   = NULL;               // no constraints are present
	keep     = 0;                  // problem has not yet been created
	for (k = 0; k < n; k++)
	{	x[k]    = 0.0; // only fixed values matter
		b[k]    = 0.5; // value of b 
		c[k]    = 1.0; // value of c
		mu[k]   = 1.5; // value of mu
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

	// set the V_xx = [ Q + diag(c) ]^{-1} = variance for x
	for(i = 0; i < n; i++)
	{	int nnbs = graph->nnbs[i];
		int *nbs = graph->nbs[i];
		for(j = 0; j < n; j++ )
			temp[i + j*n] = 0.;
		temp[i + i*n] = c[i] + Qfunc(i, i, argQfunc);
		for(k = 0; k < nnbs; k++)
		{	j = nbs[k];
			temp[i + j*n] = Qfunc(i, j, argQfunc);
		}
	}
	pos_inv_(&n, temp, work, V_xx, &logdet, &info);
	assert( info == 0 );

	// extract V_uu from V_xx
	for(i = 0; i < nu; i++)
		for(j = 0; j < nu; j++)
			V_uu[i + j*nu] = V_xx[i + j*n];

	// extract V_wu from V_xx
	for(i = 0; i < nw; i++)
		for(j = 0; j < nu; j++)
			V_wu[i + j*nw] = V_xx[(i+nu) +j*n];

	// extract V_uw from V_xx
	for(i = 0; i < nu; i++)
		for(j = 0; j < nw; j++)
			V_uw[i + j*nu] = V_xx[i + (j+nu)*n];

	// extract V_ww from V_xx
	for(i = 0; i < nw; i++)
		for(j = 0; j < nw; j++)
			V_ww[i + j*nw] = V_xx[(i+nu) + (j+nu)*n];

	// compute inv_V_uu = V_uu^{-1}
	pos_inv_(&nu, V_uu, work, inv_V_uu, &logdet, &info);
	
	// compute x_bar = mu + V_xx * b = E[x] = (E[u] , E[w])
	multiply_matrix(n, n, 1, V_xx, b, x_bar);
	for(i = 0; i < n; i++)
		x_bar[i] += mu[i];

	// compute w_hat =  E[w] + V(w,u)*V(u,u)^{-1} * [ u - E(u) ]
	// ---------------------------------------------------------
	// work = u - E(u)
	for(i = 0; i < nu; i++)
		work[i] = x[i] - x_bar[i];
	// temp = V(u,u)^{-1} * [ u - E(u) ]
	multiply_matrix(nu, nu, 1, inv_V_uu, work, temp);
	// work = V(w,u)*V(u,u)^{-1} * [ u - E(u) ]
	multiply_matrix(nw, nu, 1, V_wu, temp, work);
	// w_hat = E[w] + V(w,u)*V(u,u)^{-1} * [ u - E(u) ]
	for(i = 0; i < nw; i++)
		w_hat[i] = x_bar[i + nu] + work[i];
	// ---------------------------------------------------------

	// compute covariance = V(w,w) - V(w,u) * V(u,u)^{-1} * V(u,w)
	// ---------------------------------------------------------
	// work =  V(u,u)^{-1} * V(u,w)
	multiply_matrix(nu, nu, nw, inv_V_uu, V_uw, work);
	// temp =  V(w,u) * V(u,u)^{-1} * V(u,w)
	multiply_matrix(nw, nu, nw, V_wu, work, temp);
	// covariance = V(w,w) - V(w,u) * V(u,u)^{-1} * V(u,w)
	for(i = 0; i < nw * nw; i++)
		covariance[i] = V_ww[i] - temp[i];
	// ---------------------------------------------------------
	// inv_V_ww = inverse of covariance of w given u
	pos_inv_(&nw, covariance, work, inv_V_ww, &logdet, &info);
	assert( info == 0 );
	logdet = - logdet;  // log determinant of inv_V_ww
	  
	for(ell = 0; ell < nsim; ell++)
	{	GMRFLib_sample(problem);
		for(i = 0; i < nw; i++)
		{	sample_mean[i] += problem->sample[i + nu] / nsim;
			for(j = 0; j < nw; j++)
			{	double si = problem->sample[i + nu] - w_hat[i];
				double sj = problem->sample[j + nu] - w_hat[j];
				sample_cov[i + j*nw] += si * sj / nsim;
			}
		}
		// Example and check of evaluation of the log-likelihood
		if( ell < 10 )
		{	// check the log-density evaluation by GMRFLib
			double delta = 10. * n * DBL_EPSILON;
			double logdens = ln_gauss_density(
			nw, problem->sample + nu, w_hat, inv_V_ww, logdet
			);  
			GMRFLib_evaluate(problem);
			ok &= fabs(problem->sub_logdens - logdens) <= delta;
		}
	}
	// check the sample mean
	for(i = 0; i < nw; i++)
		ok &= fabs(w_hat[i] - sample_mean[i]) <=  epsilon;

	// check if the sample covariance 
	ok &= check_matrix(nw, nw, sample_cov, covariance, n * epsilon);

	/*
	 * Free allocated memory: 
	 */
	Free(fixed);
	Free(ptr_calloc);

	GMRFLib_free_constr(constr);
	GMRFLib_free_problem(problem);
	GMRFLib_free_graph(graph);

	return ok;
}
