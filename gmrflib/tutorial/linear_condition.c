// GMRFLib-Tutorial: Copyright (C) 2010 Bradley M. Bell
// Distributed under the terms of the GNU Public License version 2 or higher.

/*!
\file linear_condition.c
\brief Sampling of a GMRF conditioned on a linear deterministic constraint.
*/

/// identify who made the last change to this file and when
static const char RCSId[] = "$Id: linear_condition.c,v 1.9 2010/03/17 15:31:47 bradbell Exp $"; 

#include "GMRFLib/GMRFLib.h"

// evaluate the log-density for a Gaussian
# include "ln_gauss_density.h"

// compute singular decomposition for a general dense (not sparse) matrix
# include "gen_svd_.h"

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

\param i
is the row index in the precision matrix.

\param j
is the column index in the precision matrix.

\param arg
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
Sampling a GMRF conditioned on a determinstic linear constraint.

\section linear_condition_Prior_For_x Prior For x
The \ref unconditional_Prior_For_x
is given as Gaussian with mean and variance given by
\f[
\begin{array}{rcl}
{\bf E} (x) & = & \mu + [ Q + {\rm diag} (c) ]^{-1} b 
\\
{\bf V} ( x , x ) & = & [ Q + {\rm diag} (c) ]^{-1}
\end{array}
\f]

\section linear_condition_conditional_for_x conditional for x
We are given a matrix \f$ A \in {\bf R}^{m \times n} \f$,
and a value \f$ e \in {\bf R}^m \f$.
We consider the conditional distribution of \f$ x \f$ given \f$ A x = e \f$.
The matrix \f$ A \f$ has a singular decomposition
\f[
	A = U \Sigma W
\f]
where \f$ U \f$ is orthogonal (the inverse of \f$ U \f$ is \f$ U^{\rm T} \f$),
\f$ \Sigma \f$ is diagonal,
and \f$ W \f$ is orthogonal.
We define the random variables \f$ y \in {\bf R}^m \f$ by \f$ y = W x \f$.
It follows that the prior distribution for \f$ y \f$ 
is Gaussian with 
\f[
\begin{array}{rcl}
{\bf E} (y) & = & W \left\{ \mu + [ Q + {\rm diag} (c) ]^{-1} b \right\} 
\\
{\bf V} ( y , y ) 
& = & 
\left\{ W [ Q + {\rm diag} (c) ] W^{\rm T} \right\}^{-1}
=   
W [ Q + {\rm diag} (c) ]^{-1} W^{\rm T}
\end{array}
\f]
Furthermore (since \f$ m < n \f$)
the first \f$ {\bf R}^m \f$ components of \f$ y \f$ are fixed by
the equation \f$ U \Sigma y = e \f$
and the rest of the components are unconstrained.
If we use \f$ \Sigma_{u,u} \in {\bf R}^{m \times m} \f$
to denote the first \f$ m \f$ columns of \f$ \Sigma \f$,
the first \f$ m \f$ elements of \f$ y \f$ are constrained to have value
\f$ \Sigma_{u,u}^{-1} U^{\rm T} e \f$.
In summary
\f[
y  
= W x
= \left( \begin{array}{c} 
	u \\ v 
\end{array} \right)  
= \left( \begin{array}{c} 
	\Sigma_{u,u}^{-1} U^{\rm T} e \\ v 
\end{array} \right)  
\f]
where \f$ u \f$ is fixed.
Thus, we have reduced this case to the \ref fixed_condition case
with a modified prior distribution.
 
\return
1 if test passes and 0 otherwise.
*/
int linear_condition(void)
{
	// problem parameters  (you can change these)
	const int nsim        = 50000; // number of simulations
	const int seed        = 1234;  // random number seed
	const int nrow        = 4;     // number of rows in lattice
	const int ncol        = 4;     // number of column in lattice
	const int nb_row      = 1;     // 2*nb_row is # neighbouring rows
	const int nb_col      = 1;     // 2*nb_col is # neighbouring colums
	const int cyclic_flag = 0;     // is graph cyclic (on neighbourhoods)
	const int nconstr     = 2;     // number of constraints

	// a convergence criteria for sum of nsim independent samples
	const double epsilon = 4. / sqrt( (double) nsim );

	// initialize return flag
	int ok = 1;

	// the sample mean for w given u
	double *sample_mean;

	// the sample covariance for w given u
	double *sample_cov;

	// temporary matrices
	double *work, *temp;

	// corresponding variance matrices
	double *V_yy, *V_uu, *V_uw, *V_wu, *V_ww, *inv_V_ww, *inv_V_uu;

	// matrices in the SVD decomposition of A
	double *U, *S, *W, *Wt;

	// Qcinv = [Q  + diag(c)]^{-1} 
	double *Qcinv;

	// y = W * x
	double *y;

	// y_bar = W * ( mu + [Q  + diag(c)]^{-1} * b }
	double *y_bar;
	
	// u_fixed = Sigma_{u,u}^{-1} * U^T * e
	double *u_fixed;

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

	// number of free nodes; i.e., number of components in v
	int          nw;

	// number of fixed nodes; i.e., number of components in u
	int          nu;

	// size of the work array
	int          lwork;

	// temporary indices
	int i, j, k, ell;

	// See below for the meaning of these variables 
	double      *A;
	double      *e;
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
	n    = graph->n;
	// number of components in u (constraints)
	nu    = nconstr;
	assert( nu < n );
	// number of components in v
	nw    = n - nu;
	// make sure work is large enough for an n*n array and for gen_svd
	lwork = 5 * nu + n * n;

	// number of doubles to allocate
	n_calloc  = 6*n + 5*n*n + 2*nu + 2*nw + 3*nu*nu + 2*nu*nw + 4*nw*nw;
	n_calloc += lwork;
	// allocate the memory
	ptr_calloc = (double*) Calloc(n_calloc, double);
	// split the memory up into individual vectors
	i  = 0;
	x           = ptr_calloc + i; i += n;     // x           length n
	b           = ptr_calloc + i; i += n;     // b           length n
	c           = ptr_calloc + i; i += n;     // c           length n
	mu          = ptr_calloc + i; i += n;     // mu          length n
	y           = ptr_calloc + i; i += n;     // y           length n
	y_bar       = ptr_calloc + i; i += n;     // y_bar       length n
	V_yy        = ptr_calloc + i; i += n*n;   // V_yy        length n*n
	Qcinv       = ptr_calloc + i; i += n*n;   // Qcinv       length n*n
	temp        = ptr_calloc + i; i += n*n;   // temp        lenght n*n
	W           = ptr_calloc + i; i += n*n;   // W           lenght n*n
	Wt          = ptr_calloc + i; i += n*n;   // Wt          lenght n*n
	S           = ptr_calloc + i; i += nu;    // S           lenght nu 
	u_fixed     = ptr_calloc + i; i += nu;    // u_fixed     lenght nu 
	sample_mean = ptr_calloc + i; i += nw;    // sample_mean length nw
	w_hat       = ptr_calloc + i; i += nw;    // w_hat       length nw
	U           = ptr_calloc + i; i += nu*nu; // U           lenght nu*nu
	V_uu        = ptr_calloc + i; i += nu*nu; // V_uu        length nu*nu
	inv_V_uu    = ptr_calloc + i; i += nu*nu; // inv_V_uu    length nu*nu
	V_uw        = ptr_calloc + i; i += nu*nw; // V_uw        length nu*nw
	V_wu        = ptr_calloc + i; i += nw*nu; // V_uw        length nw*nu
	V_ww        = ptr_calloc + i; i += nw*nw; // V_uu        length nw*nw
	inv_V_ww    = ptr_calloc + i; i += nw*nw; // inv_V_ww    length nw*nw
	sample_cov  = ptr_calloc + i; i += nw*nw; // sample_cov  length nw*nw
	covariance  = ptr_calloc + i; i += nw*nw; // covariance  length nw*nw
	work        = ptr_calloc + i; i += lwork; // work        lenght lwork;
	assert( i == n_calloc );

	// -------------------------------------------------------------------
	// Create the linear deterministic constraint 
	GMRFLib_make_empty_constr(&constr);
	// # of columns in the matrix A
	constr->nc = nu;  
	// Define the matrix A with first row [1, 0] , second row [1, 1],
	// third row [1, 0], fourth row [1, 1], ...
	// Note that the order of the elements of A required by GMRFLib is
	// the transpose of the standard order used in C.
	A                = Calloc(nu * n, double); 
	constr->a_matrix = A;
	for(i = 0; i < n; i++)
	{	A[0 + i*nu] = 1.;	
		for(k = 1; k < nu; k++)
			A [k + i*nu]  = (double) (((i + k) % nu) == 0); 
	}
	// Define the vector e.  
	e                = Calloc(nu, double);
	constr->e_vector = e;
	for(k = 0; k < nu; k++)
		e[k] = 2. / (double) (1 + k);
	// Compute (and store in constr) unspecified internal information that 
	// is needed by the sampling routines.
	{	int scale_constr = 0;
		GMRFLib_prepare_constr(constr, graph, scale_constr);
	}
	// -------------------------------------------------------------------

	// Create fixed conditional sampling problem
	argQfunc = (void *) graph;     // extra arguments to Qfunc
	fixed    = NULL;
	keep     = 0;                  // problem has not yet been created
	for (k = 0; k < n; k++)
	{	x[k]    = 0.0; // values not used becasue fixed == NULL
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

	// compute the singular decomposition of A = U * diag(S) * W
	gen_svd_(&nu, &n, A, temp, work, U, S, W, &info);
	assert( info == 0 );
	// ------------------------------------------------------------
	// check the singular value computation
	// work = diag(S)
	for(i = 0; i < nu; i++)
	{	for(j = 0; j < n; j++)
			work[i + j*nu] = 0.; 
		work[i + i*nu] = S[i];
	}
	// temp = diag(S) * W
	multiply_matrix(nu, n, n, work, W, temp);
	// work = U * diag(S) * W
	multiply_matrix(nu, nu, n, U, temp, work);
	// compare with A
	ok &= check_matrix(nu, n, A, work, 10. * n * DBL_EPSILON);
	// ------------------------------------------------------------

	// set Wt =  W^T 
	for(i = 0; i < n; i++)
		for(j = 0; j < n; j++)
			Wt[i + j*n] = W[j + i*n];

	// compute Qcinv = [ Q + diag(c) ]^{-1} 
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
	pos_inv_(&n, temp, work, Qcinv, &logdet, &info);
	assert( info == 0 );

	// compute V_yy = W * [ Q + diag(c) ]^{-1} * W^T
	multiply_matrix(n, n, n, Qcinv, Wt, temp);
	multiply_matrix(n, n, n, W, temp, V_yy);

	// extract V_uu from V_yy
	for(i = 0; i < nu; i++)
		for(j = 0; j < nu; j++)
			V_uu[i + j*nu] = V_yy[i + j*n];

	// extract V_uw from V_yy
	for(i = 0; i < nu; i++)
		for(j = 0; j < nw; j++)
			V_uw[i + j*nu] = V_yy[i + (j+nu)*n];

	// extract V_wu from V_yy
	for(i = 0; i < nw; i++)
		for(j = 0; j < nu; j++)
			V_wu[i + j*nw] = V_yy[(i+nu) +j*n];

	// extract V_ww from V_yy
	for(i = 0; i < nw; i++)
		for(j = 0; j < nw; j++)
			V_ww[i + j*nw] = V_yy[(i+nu) + (j+nu)*n];

	// compute inv_V_uu = V_uu^{-1}
	pos_inv_(&nu, V_uu, work, inv_V_uu, &logdet, &info);
	
	// compute y_bar = W { mu + [ Q + diag(c) ]^{-1} b }
	multiply_matrix(n, n, 1, Qcinv, b, temp);
	for(i = 0; i < n; i++)
		temp[i] += mu[i];
	multiply_matrix(n, n, 1, W, temp, y_bar);

	// compute u_fixed = Sigma_{u,u}^{-1} * U^T * e
	for(i = 0; i < nu; i++)
	{	double sum = 0;
		for(k = 0; k < nu; k++)
			sum += U[k + i*nu] * e[k];
		u_fixed[i] = sum / S[i];
	}  

	// compute w_hat =  E[w] + V(w,u)*V(u,u)^{-1} * [ u_fixed - E(u) ]
	// ---------------------------------------------------------
	// work = u_fixed - E(u)
	for(i = 0; i < nu; i++)
		work[i] = u_fixed[i] - y_bar[i];
	// temp = V(u,u)^{-1} * [ u_fixed - E(u) ]
	multiply_matrix(nu, nu, 1, inv_V_uu, work, temp);
	// work = V(w,u)*V(u,u)^{-1} * [ u_fixed - E(u) ]
	multiply_matrix(nw, nu, 1, V_wu, temp, work);
	// w_hat = E[w] + V(w,u)*V(u,u)^{-1} * [ u_fixed - E(u) ]
	for(i = 0; i < nw; i++)
		w_hat[i] = y_bar[i + nu] + work[i];
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
	logdet = - logdet;  // log( det( inv_V_ww ) )
	  
	for(ell = 0; ell < nsim; ell++)
	{	GMRFLib_sample(problem);
		// y = W * problem->sample
		multiply_matrix(n, n, 1, W, problem->sample, y);
		for(i = 0; i < nw; i++)
		{	sample_mean[i] += y[i + nu] / nsim;
			for(j = 0; j < nw; j++)
			{	double si = y[i + nu] - w_hat[i];
				double sj = y[j + nu] - w_hat[j];
				sample_cov[i + j*nw] += si * sj / nsim;
			}
		}
		// Example and check of evaluation of the log-likelihood
		if( ell < 10 )
		{	// check the log-density evaluation by GMRFLib
			double delta = 100. * n * DBL_EPSILON;
			double logdens = ln_gauss_density(
				nw, y + nu, w_hat, inv_V_ww, logdet
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
