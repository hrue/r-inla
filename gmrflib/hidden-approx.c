
/* hidden-approx.c
 * 
 * Copyright (C) 2001-2006 Havard Rue
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * The author's contact information:
 *
 *        Haavard Rue
 *        CEMSE Division
 *        King Abdullah University of Science and Technology
 *        Thuwal 23955-6900, Saudi Arabia
 *        Email: haavard.rue@kaust.edu.sa
 *        Office: +966 (0)12 808 0640
 *
 */

/*!
  \file hidden-approx.c
  \brief Sampling from a hidden GMRF

  \section hidden Approximating hidden GMRFs
  Approximations of a non-normal density of the form

  \f[ \pi(\mbox{\boldmath $x|y$}) \propto \exp\left( -\frac{1}{2} (\mbox{\boldmath
  $x$}-\mbox{\boldmath $\mu$})^T (\mbox{\boldmath $Q$} + \mbox{diag}(\mbox{\boldmath $c$}))
  (\mbox{\boldmath $x$}-\mbox{\boldmath $\mu$}) + \mbox{\boldmath $b$}^T\mbox{\boldmath $x$} +\sum
  d_i f_i(x_i,y_i) \right). \hspace{2cm} (HGMRF-1) \f]

  We assume that the terms \f$ f_i(x_i,y_i) \f$ contain also non-quadratic terms of \f$ x_i \f$, so
  that <em>\b x|y</em> is not normal. \n If <em>\b x</em> is a GMRF wrt the graph \em G which is
  partially observed through <em>\b y</em>, then <em>\b x|y</em> is called a hidden GMRF (HGMRF).
  Note that (HGMRF-1) defines <em>\b x</em> as a Markov random field wrt (the same graph) \em G, but
  it is not Gaussian.\n

  See also \ref hidden_GMRF
*/

#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

/* Pre-hg-Id: $Id: hidden-approx.c,v 1.48 2010/01/16 12:31:52 hrue Exp $ */

#include <assert.h>
#include <float.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif

#include <gsl/gsl_statistics.h>
#include <gsl/gsl_cdf.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

/*!
  \brief 
  \param hidden_par A pointer to a \c GMRFLib_hidden_param_tp pointer. 
  At output the \c GMRFLib_hidden_param_tp -object contains the default values

  \par Description of elements in \c GMRFLib_hidden_param_tp -object:
  \em range: Range of the spline approximation \n
  <b>Default value: 6.0 (standard deviations) </b> \n
  \em step_len: Step-length for Taylor expansion \n
  <b>Default value: 1.0e-4 </b> \n
  \em neighpar: How many neigbhour neighbors to include when making the neighbour-graph \n
  <b>Default value: 1 </b> \n
  \em neightype: Use the neighbours (#GMRFLib_NEIGHTYPE_GRAPH) or just indices with
  lower index (#GMRFLib_NEIGHTYPE_LINEAR) \n
  <b>Default value: \c GMRFLib_NEIGHTYPE_GRAPH</b> \n
  \em norder: Order in the fitted log-spline: 1 or 2 \n
  <b>Default value: 2 </b> \n
  \em nresolution: Number of points in the spline \n
  <b>Default value: 12 </b> \n
  \em nsample: How many samples to produce when estimating \f$\hat I\f$ \n
  <b>Default value: 4 </b> \n
  \em modeoption: Exapand around the current value (#GMRFLib_MODEOPTION_CURRENT) or
  the mode (#GMRFLib_MODEOPTION_MODE) \n
  <b>Default value: #GMRFLib_MODEOPTION_MODE </b> \n
  \em gaussapprox: If #GMRFLib_TRUE, use the Gaussian approximation instead of computing the non-Gaussian
  approximation, overriding all other options. \n
  <b>Default value: 0 </b> \n
 */
int GMRFLib_default_hidden_par(GMRFLib_hidden_param_tp ** hidden_par)
{
	*hidden_par = Calloc(1, GMRFLib_hidden_param_tp);

	(*hidden_par)->neightype = GMRFLib_NEIGHTYPE_LINEAR;
	(*hidden_par)->neightype = GMRFLib_NEIGHTYPE_GRAPH;
	(*hidden_par)->neighpar = 1;			       /* how to make the neighbor-graph */
	(*hidden_par)->nresolution = 12;		       /* number of points in the spline */
	(*hidden_par)->norder = 2;			       /* order in the fitted spline: 1 or 2 */
	(*hidden_par)->nsample = 4;			       /* number of basic samples */
	(*hidden_par)->nantithetic = 1;			       /* # of antithetic samples is GMRFLib_NUM_ANTITHETIC (=4) times this: +- x
							        * (u,1-u)scale */
	(*hidden_par)->range = 6.0;			       /* + = range * stdev */
	(*hidden_par)->cmeanmode = GMRFLib_COND_MODE;
	(*hidden_par)->cmeanmode = GMRFLib_COND_MEAN;
	(*hidden_par)->modeoption = GMRFLib_MODEOPTION_CURRENT;
	(*hidden_par)->modeoption = GMRFLib_MODEOPTION_MODE;
	(*hidden_par)->gaussapprox = 0;			       /* if true, use the gaussian approximation instead */

	(*hidden_par)->ngraph = NULL;			       /* this will be computed later on */
	(*hidden_par)->step_len = GMRFLib_eps(0.25);	       /* step-length for Taylor expantion */
	(*hidden_par)->stencil = 5;			       /* 3,5,7 */

	/*
	 * the expert-options are always default OFF! 
	 */
	(*hidden_par)->remap = NULL;			       /* use default remap */
	(*hidden_par)->stop_idx = 0;			       /* full sweeep */
	(*hidden_par)->marginal_density = NULL;		       /* no marginal-density [last idx] */
	(*hidden_par)->marginal_mean = NULL;		       /* no marginal-mean [last idx] */
	(*hidden_par)->marginal_stdev = NULL;		       /* no marginal-stdev [last idx] */

	return GMRFLib_SUCCESS;
}
int GMRFLib_loglFunc_wrapper(double *logll, double *x, int m, int idx, double *x_vec, double *y_cdf, void *arg)
{
	/*
	 * adjust for index, the second order contrib and zero mean 
	 */

	int i;
	GMRFLib_logl_arg_tp *args;

	static double *ff = NULL, *xx = NULL, a, b, c, mean;
	static int len = 0;

	if (m > len) {					       /* static storage, never free'ed */
		Free(ff);
		Free(xx);
		len = m;
		ff = Calloc(len, double);
		xx = Calloc(len, double);
	}

	args = (GMRFLib_logl_arg_tp *) arg;
	mean = args->sub_mean[idx];
	for (i = 0; i < m; i++)
		xx[i] = x[i] + mean;
	(*args->loglFunc) (ff, xx, m, args->map[idx], x_vec, NULL, args->loglFunc_arg);

	a = args->acoof[idx];
	b = args->bcoof[idx];
	c = -0.5 * args->ccoof[idx];

	for (i = 0; i < m; i++)
		logll[i] = ff[i] - (a + (b + c * xx[i]) * xx[i]);

	if (0)						       /* FIXME */
		for (i = 0; i < m; i++)
			logll[i] = DMIN(logll[i], 0.0);

	return GMRFLib_SUCCESS;
}

/*!
  \brief Initializes and specifies a \c GMRFLib_hidden_problem_tp -object holding all 
  information needed for sampling from a HGMRF.
  
  \param[in,out] hidden_problem At output \a hidden_problem is a pointer to a 
  \c GMRFLib_hidden_problem_tp -object, initialized and defined according to
  the problem specification.
  \param[in] x A length \em n array, where \n is the number of nodes in the 
  graph, of initial values of the HGMRF. If \a x = \c NULL then all elements
  are taken to be zero.  If \a fixed_value \f$ \neq \f$ \c NULL, the elements of 
  <em>\b x</em> corresponding to \a fixed_value=1 are the fixed values in a
  conditional simulation. The remaining elements of <em>\b x</em> are not used, 
  and can take arbitrary values. If \a fixed_value = \c NULL, all values can be
  arbitrary.
  \param[in] b If <tt>!NULL</tt>, a length \em n array holding the elements of 
  the vector <em>\b b</em> in the general expression of the density if a HGMRF, 
  as given in <b>(HGMRF-1)</b> in \ref hidden.
  \param[in] c If <tt>!NULL</tt>, a length \em n array of elements to add to 
  the diagonal of the precision matrix defined by the function \a Qfunc. The 
  argument \a c should hold the elements of the vector <em>\b c</em> in 
  <b>(HGMRF-1)</b>.
  \param[in] mean If <tt>!NULL</tt>, a length \em n array holding the elements 
  of the vector \f$ \mbox{\boldmath $\mu$} \f$ in <b>(HGMRF-1)</b>.  
  \param[in] graph The graph on which the HGMRF is defined.
  \param[in] Qfunc A pointer to a user-defined function defining the 
  precision matrix <em>\b Q</em> of a HGMRF
  \param[in] Qfunc_args The arguments to the function \a Qfunc defining 
  the precision matrix.
  \param[in] fixed_value If \c !NULL, the sampling is done conditionally 
  on fixed values.
  \param[in] d If <tt>!NULL</tt>, a length \em n array holding the elements of 
  the vector <em>\b d</em> in <b>(HGMRF-1)</b>.
  \param[in] loglFunc A function of type \c GMRFLib_logl_tp(), returning the value
  of the function \f$ f_i(x_i,y_i) \f$ in (HGMRF-1) in \ref hidden, in many applications
  equal to the log-likelihood of the problem.
  \param[in] loglFunc_arg A \em void pointer holding the address of a variable
  or data structure defining additional arguments to the function \a loglFunc
  \param[in] optpar The options to the optimizer. If \c NULL, default options are used.
  \param[in] hidden_par If \c NULL, created by \c GMRFLib_default_hidden_par()

  \note If you want to use the approximation \f$ \hat I = 0\f$,
  put \a neighpar = 0 in the \c GMRFLib_hidden_param_tp -object \a hidden_par.
 */
int GMRFLib_init_problem_hidden(GMRFLib_hidden_problem_tp ** hidden_problem,
				double *x, double *b, double *c, double *mean,
				GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_args,
				char *fixed_value,
				double *d, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, GMRFLib_optimize_param_tp * optpar,
				GMRFLib_hidden_param_tp * hidden_par)
{
	int retval;

	GMRFLib_ENTER_ROUTINE;
	retval = GMRFLib_init_problem_hidden_store(hidden_problem, x, b, c, mean, graph, Qfunc, Qfunc_args, fixed_value,
						   d, loglFunc, loglFunc_arg, optpar, hidden_par, NULL);
	GMRFLib_LEAVE_ROUTINE;
	return retval;
}
int GMRFLib_init_problem_hidden_store(GMRFLib_hidden_problem_tp ** hidden_problem,
				      double *x, double *b, double *c, double *mean,
				      GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_args,
				      char *fixed_value,
				      double *d, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg,
				      GMRFLib_optimize_param_tp * optpar, GMRFLib_hidden_param_tp * hidden_par, GMRFLib_store_tp * store)
{
	int sub_n, i, j, node, nnode, free_fixed_value = 0, free_x = 0, *inv_remap;
	double *mode, *cc;
	GMRFLib_sm_fact_tp tmp_fact;
	GMRFLib_graph_tp *tmp_graph = NULL;
	GMRFLib_logl_arg_tp *args;

	GMRFLib_ASSERT(graph, GMRFLib_EINVARG);
	GMRFLib_ASSERT(Qfunc, GMRFLib_EINVARG);

	GMRFLib_ENTER_ROUTINE;

	if (store)
		FIXME("store-functionality is not yet implemented.");

	if (!fixed_value) {
		fixed_value = Calloc(graph->n, char);

		free_fixed_value = 1;
	}
	if (!x) {
		x = Calloc(graph->n, double);

		free_x = 1;
	}

	(*hidden_problem) = Calloc(1, GMRFLib_hidden_problem_tp);
	if (hidden_par) {
		/*
		 * copy remap (if !NULL) later, as i need to know the size of sub_graph 
		 */
		(*hidden_problem)->hidden_par = Calloc(1, GMRFLib_hidden_param_tp);
		memcpy((*hidden_problem)->hidden_par, hidden_par, sizeof(GMRFLib_hidden_param_tp));
		(*hidden_problem)->hidden_par->remap = NULL;   /* yes, fix this later */
	} else
		GMRFLib_default_hidden_par(&((*hidden_problem)->hidden_par));

	/*
	 * first locate the mode, or similar, to approximate around 
	 */
	mode = Calloc(graph->n, double);
	memcpy(mode, x, graph->n * sizeof(double));
	if ((*hidden_problem)->hidden_par->modeoption == GMRFLib_MODEOPTION_MODE) {
		GMRFLib_EWRAP1(GMRFLib_optimize(mode, b, c, mean, graph, Qfunc, Qfunc_args, fixed_value, NULL, d, loglFunc, loglFunc_arg, optpar));
	}

	/*
	 * convert the variables so we get rid of the fixed_values and then write the function in the canonical form
	 * 
	 * large parts here is adapted from `problem-setup.c'!!! 
	 */

	/*
	 * define the new graph (ok if fixed_value = NULL) 
	 */
	GMRFLib_EWRAP1(GMRFLib_compute_subgraph(&tmp_graph, graph, fixed_value));
	sub_n = tmp_graph->n;
	if (sub_n == 0) {				       /* fast return if there is nothing todo */
		Free((*hidden_problem));
		GMRFLib_free_graph(tmp_graph);

		if (free_fixed_value)
			Free(fixed_value);
		if (free_x)
			Free(x);

		GMRFLib_LEAVE_ROUTINE;
		return GMRFLib_SUCCESS;
	}

	/*
	 * check the hidden_par for a vital parameter regarding the marginal 
	 */
	GMRFLib_ASSERT((*hidden_problem)->hidden_par->stop_idx >= 0, GMRFLib_EINVARG);
	GMRFLib_ASSERT((*hidden_problem)->hidden_par->stop_idx < sub_n, GMRFLib_EINVARG);

	/*
	 * compute reordering, either given or to be computed 
	 */
	tmp_fact.smtp = GMRFLib_smtp;
	if (hidden_par && hidden_par->remap) {
		tmp_fact.remap = Calloc(sub_n, int);
		memcpy(tmp_fact.remap, hidden_par->remap, sub_n * sizeof(int));
		GMRFLib_EWRAP1(GMRFLib_compute_bandwidth(&tmp_fact.bandwidth, tmp_graph, tmp_fact.remap));

		(*hidden_problem)->hidden_par->remap = Calloc(sub_n, int);
		memcpy((*hidden_problem)->hidden_par->remap, hidden_par->remap, sub_n * sizeof(int));

		if (0)
			printf("use users remap with bandw %d\n", tmp_fact.bandwidth);
	} else {
		GMRFLib_EWRAP1(GMRFLib_compute_reordering(&tmp_fact, tmp_graph, NULL));
	}

	/*
	 * map the indices so the mapping is identity
	 * 
	 * j in remapped graph is mothergraph_idx[inv_remap[j]] in real world
	 * 
	 */
	(*hidden_problem)->sub_sm_fact.smtp = tmp_fact.smtp;
	(*hidden_problem)->sub_sm_fact.bandwidth = tmp_fact.bandwidth;
	(*hidden_problem)->sub_sm_fact.remap = Calloc(sub_n, int);

	for (i = 0; i < sub_n; i++)
		(*hidden_problem)->sub_sm_fact.remap[i] = i;
	inv_remap = Calloc(sub_n, int);
	(*hidden_problem)->map = Calloc(sub_n, int);

	for (i = 0; i < sub_n; i++)
		inv_remap[tmp_fact.remap[i]] = i;
	for (i = 0; i < sub_n; i++)
		(*hidden_problem)->map[i] = tmp_graph->mothergraph_idx[inv_remap[i]];
	GMRFLib_EWRAP1(GMRFLib_remap_graph(&((*hidden_problem)->sub_graph), tmp_graph, tmp_fact.remap));
	(*hidden_problem)->sub_graph->mothergraph_idx = (*hidden_problem)->map;	/* dont think this is needed? */
	(*hidden_problem)->sub_graph->mothergraph_idx = NULL;  /* ...therefore i put it NULL to fail if it is */

	Free(inv_remap);
	GMRFLib_free_graph(tmp_graph);
	GMRFLib_free_reordering(&tmp_fact);

	/*
	 * setup space & misc. sub_mean is also sub_mode as we're generating a gaussian approximation centered at the mode 
	 */
	(*hidden_problem)->sample = Calloc(graph->n, double);
	(*hidden_problem)->sample_star = Calloc(graph->n, double);
	(*hidden_problem)->sub_sample = Calloc(graph->n, double);
	(*hidden_problem)->sub_mean = Calloc(sub_n, double);
	(*hidden_problem)->sub_b = Calloc(sub_n, double);
	(*hidden_problem)->sub_d = Calloc(sub_n, double);
	cc = Calloc(sub_n, double);

	memcpy((*hidden_problem)->sample, x, graph->n * sizeof(double));
	memcpy((*hidden_problem)->sample_star, mode, graph->n * sizeof(double));
	for (i = 0; i < sub_n; i++)
		(*hidden_problem)->sub_mean[i] = mode[(*hidden_problem)->map[i]];
	if (b)
		for (i = 0; i < sub_n; i++)
			(*hidden_problem)->sub_b[i] = b[(*hidden_problem)->map[i]];
	if (d)
		for (i = 0; i < sub_n; i++)
			(*hidden_problem)->sub_d[i] = d[(*hidden_problem)->map[i]];
	if (c)
		for (i = 0; i < sub_n; i++)
			cc[i] = c[(*hidden_problem)->map[i]];

	(*hidden_problem)->x_vec = Calloc(graph->n, double);
	memcpy((*hidden_problem)->x_vec, mode, graph->n * sizeof(double));

	(*hidden_problem)->ran_approx_state = GMRFLib_uniform_getstate(NULL);

	/*
	 * make the arguments to the wrapper function 
	 */
	(*hidden_problem)->sub_Qfunc = GMRFLib_Qfunc_wrapper;
	(*hidden_problem)->sub_Qfunc_arg = Calloc(1, GMRFLib_Qfunc_arg_tp);
	(*hidden_problem)->sub_Qfunc_arg->map = (*hidden_problem)->map;	/* yes, this ptr is needed */
	(*hidden_problem)->sub_Qfunc_arg->user_Qfunc = Qfunc;
	(*hidden_problem)->sub_Qfunc_arg->user_Qfunc_args = Qfunc_args;
	(*hidden_problem)->sub_Qfunc_arg->diagonal_adds = cc;

	GMRFLib_ASSERT(((*hidden_problem)->hidden_par->neightype == GMRFLib_NEIGHTYPE_LINEAR ||
			(*hidden_problem)->hidden_par->neightype == GMRFLib_NEIGHTYPE_GRAPH), GMRFLib_EPARAMETER);
	GMRFLib_ASSERT((*hidden_problem)->hidden_par->cmeanmode == GMRFLib_COND_MEAN ||
		       (*hidden_problem)->hidden_par->cmeanmode == GMRFLib_COND_MODE, GMRFLib_EPARAMETER);
	GMRFLib_ASSERT((*hidden_problem)->hidden_par->norder == 1 || (*hidden_problem)->hidden_par->norder == 2, GMRFLib_EPARAMETER);
	GMRFLib_ASSERT((*hidden_problem)->hidden_par->neighpar >= 0, GMRFLib_EPARAMETER);
	GMRFLib_ASSERT((*hidden_problem)->hidden_par->nresolution >= 0, GMRFLib_EPARAMETER);
	GMRFLib_ASSERT((*hidden_problem)->hidden_par->nsample >= 0, GMRFLib_EPARAMETER);
	GMRFLib_ASSERT((*hidden_problem)->hidden_par->nantithetic >= 0, GMRFLib_EPARAMETER);
	GMRFLib_ASSERT((*hidden_problem)->hidden_par->range >= 0., GMRFLib_EPARAMETER);

	/*
	 * now compute the new 'effective' b, and then the mean. recall to add the 'c' term manually, since we're using the
	 * original Qfunc.
	 * 
	 * x=(x1,x2), then x1|x2 has b = Q11 \mu1 - Q12(x2-\mu2) 
	 */
	for (i = 0; i < sub_n; i++) {			       /* loop over all sub_nodes */
		node = (*hidden_problem)->map[i];
		if (mean)
			(*hidden_problem)->sub_b[i] += ((*Qfunc) (node, node, Qfunc_args) + cc[i]) * mean[node];	/* add diagonal-term */

		for (j = 0; j < graph->nnbs[node]; j++) {      /* then over all neighbors */
			double qvalue;

			nnode = graph->nbs[node][j];
			qvalue = (*Qfunc) (node, nnode, Qfunc_args);

			if (fixed_value[nnode]) {
				/*
				 * nnode is fixed 
				 */
				if (mean)
					(*hidden_problem)->sub_b[i] -= qvalue * (mode[nnode] - mean[nnode]);
				else
					(*hidden_problem)->sub_b[i] -= qvalue * mode[nnode];
			} else {
				/*
				 * nnone is not fixed 
				 */
				if (mean)
					(*hidden_problem)->sub_b[i] += qvalue * mean[nnode];
			}
		}
	}

	/*
	 * setup wrapper for the loglFunc, this involves BOTH the 'index' and correction for the quadratic term 
	 */
	(*hidden_problem)->loglFunc = GMRFLib_loglFunc_wrapper;
	args = Calloc(1, GMRFLib_logl_arg_tp);
	args->loglFunc = loglFunc;
	args->loglFunc_arg = loglFunc_arg;
	args->map = (*hidden_problem)->map;
	args->sub_mean = (*hidden_problem)->sub_mean;
	args->acoof = Calloc(sub_n, double);
	args->bcoof = Calloc(sub_n, double);
	args->ccoof = Calloc(sub_n, double);

	(*hidden_problem)->loglFunc_arg = (void *) args;

	/*
	 * add term from the loglikelihoodterm 
	 */
	for (i = 0; i < sub_n; i++)
		(*hidden_problem)->x_vec[(*hidden_problem)->map[i]] = mode[i];
	for (i = 0; i < sub_n; i++) {
		if ((*hidden_problem)->sub_d[i]) {
			GMRFLib_2order_approx(&args->acoof[i], &args->bcoof[i], &args->ccoof[i],
					      (*hidden_problem)->sub_d[i], (*hidden_problem)->sub_mean[i],
					      (*hidden_problem)->map[i], (*hidden_problem)->x_vec,
					      loglFunc, loglFunc_arg, &((*hidden_problem)->hidden_par->step_len),
					      &((*hidden_problem)->hidden_par->stencil));
			args->ccoof[i] = DMAX(0.0, args->ccoof[i]);
		} else
			args->acoof[i] = args->bcoof[i] = args->ccoof[i] = 0.0;

		(*hidden_problem)->sub_b[i] += args->bcoof[i];
		(*hidden_problem)->sub_Qfunc_arg->diagonal_adds[i] += args->ccoof[i];
	}

	/*
	 * now compute the Cholesky-factorisation. I solve one additonal time to get the mean, even though I know it already as 
	 * the mode [cost is nothing...] 
	 */

	GMRFLib_EWRAP1(GMRFLib_build_sparse_matrix(&((*hidden_problem)->sub_sm_fact), (*hidden_problem)->sub_Qfunc,
						   (char *) ((*hidden_problem)->sub_Qfunc_arg), (*hidden_problem)->sub_graph));
	GMRFLib_EWRAP1(GMRFLib_factorise_sparse_matrix(&((*hidden_problem)->sub_sm_fact), (*hidden_problem)->sub_graph));

	memcpy((*hidden_problem)->sub_mean, (*hidden_problem)->sub_b, sub_n * sizeof(double));
	GMRFLib_EWRAP1(GMRFLib_solve_llt_sparse_matrix
		       ((*hidden_problem)->sub_mean, 1, &((*hidden_problem)->sub_sm_fact), (*hidden_problem)->sub_graph));

	/*
	 * determine the neighbor-graph 
	 */
	switch ((*hidden_problem)->hidden_par->neightype) {
	case GMRFLib_NEIGHTYPE_LINEAR:
		GMRFLib_EWRAP1(GMRFLib_make_linear_graph(&((*hidden_problem)->hidden_par->ngraph),
							 (*hidden_problem)->sub_graph->n, (*hidden_problem)->hidden_par->neighpar, 0));
		break;
	case GMRFLib_NEIGHTYPE_GRAPH:
		GMRFLib_EWRAP1(GMRFLib_nfold_graph(&((*hidden_problem)->hidden_par->ngraph),
						   (*hidden_problem)->sub_graph, (*hidden_problem)->hidden_par->neighpar));
		break;
	}

	Free(mode);
	if (free_fixed_value)
		Free(fixed_value);
	if (free_x)
		Free(x);

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

/*!
  \brief Deallocates all allocated arrays in a \c GMRFLib_hidden_problem_tp -object.
*/
int GMRFLib_free_hidden(GMRFLib_hidden_problem_tp * hidden_problem)
{
	/*
	 * free a hidden_problem 
	 */
	GMRFLib_hidden_problem_tp *h = hidden_problem;

	if (!h)
		return GMRFLib_SUCCESS;

	Free(h->sample);
	Free(h->sample_star);
	Free(h->sub_sample);
	Free(h->ran_approx_state);
	Free(h->ran_current_state);
	Free(h->map);
	Free(h->sub_mean);
	Free(h->sub_b);
	Free(h->sub_d);
	GMRFLib_free_fact_sparse_matrix(&(h->sub_sm_fact));
	GMRFLib_free_reordering(&(h->sub_sm_fact));
	Free(h->x_vec);
	GMRFLib_free_graph(h->sub_graph);

	if (h->sub_Qfunc_arg) {
		Free(h->sub_Qfunc_arg->diagonal_adds);
		Free(h->sub_Qfunc_arg);
	}
	if (h->loglFunc_arg) {
		Free(((GMRFLib_logl_arg_tp *) h->loglFunc_arg)->acoof);
		Free(((GMRFLib_logl_arg_tp *) h->loglFunc_arg)->bcoof);
		Free(((GMRFLib_logl_arg_tp *) h->loglFunc_arg)->ccoof);
		Free(h->loglFunc_arg);
	}
	GMRFLib_free_graph(h->hidden_par->ngraph);
	Free(h->hidden_par->remap);
	/*
	 * to free this is the responsibility of the user! don't free it here. 
	 */
	if (0)
		if (h->hidden_par->marginal_density)
			GMRFLib_gdens_Free(*(h->hidden_par->marginal_density));

	Free(h->hidden_par);
	Free(h);

	return GMRFLib_SUCCESS;
}

/*!
  \brief Samples one realization of the elements of a HGMRF.

  \param[in,out] hidden_problem At input \a hidden_problem should contain the problem specification as defined by a call to \c
  GMRFLib_init_problem_hidden().  At output, a sample of the HGMRF has been generated, and is stored in the \em sample member of
  the data structure \a hidden_problem.  Also, the log-density of the sample is computed, and stored in the \a sub_logdens
  member. If \a hidden_problem = \c NULL, nothing is done, and the routine returns 0.
  
  \remark To sample one realization from a HGMRF <em>\b x|y</em>, first invoke \c GMRFLib_init_problem_hidden() generating the
  problem specification by initializing a \c GMRFLib_hidden_problem_tp -object \a hidden_problem, and then call \c
  GMRFLib_hidden_sample() using \a hidden_problem as an argument. Repeated samples are generated by repeated calls to \c
  GMRFLib_hidden_sample(), using the same problem specification object. There is no need to call \c
  GMRFLib_init_problem_hidden() more than once for each sampling problem. \n Prior to sampling, the built-in random number
  generator has to be initialized.  This is done by calling the function (<tt>*GMRFLib_uniform_init</tt>), taking one integer
  argument, specifying the seed of the generator. The user might replace the default random generator by another, see \c
  globals.h for details. \n At output, the sampling results are stored in the \c GMRFLib_hidden_problem_tp -object. No functions
  are provided for extracting the results of the sampling, so the results are retrieved by explicitly accessing the members of
  \a hidden_problem after running the routine. The sampled values are stored in <tt>hidden_problem->sample</tt> and the
  corresponding log-density in <tt>hidden_problem->sub_logdens</tt>. The fixed values are copied into the corresponding elements
  in \a sample.

  \sa GMRFLib_init_problem_hidden, GMRFLib_evaluate_hidden.
 */
int GMRFLib_sample_hidden(GMRFLib_hidden_problem_tp * hidden_problem)
{
	GMRFLib_ENTER_ROUTINE;
	GMRFLib_EWRAP1(GMRFLib_doit_hidden(hidden_problem, 1));
	GMRFLib_EWRAP1(GMRFLib_doit_hidden(hidden_problem, 0));
	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

/*!
  \brief Evaluates the log-density of a sample.

  \param[in,out] hidden_problem At input \a hidden_problem should contain
  the problem specification and a sample.
  At output, the log-density of the current sample has been computed, and is stored
  in the \em sub_logdens member of the data structure \a hidden_problem.
  
  \sa GMRFLib_init_problem_hidden, GMRFLib_sample_hidden.
*/
int GMRFLib_evaluate_hidden(GMRFLib_hidden_problem_tp * hidden_problem)
{
	GMRFLib_ENTER_ROUTINE;
	GMRFLib_EWRAP1(GMRFLib_doit_hidden(hidden_problem, 0));
	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

/*!
  \brief Internal routine for \c GMRFLib_sample_hidden() and \c GMRFLib_evaluate_hidden().

  \param[in,out] hidden_problem At input \a hidden_problem should contain the problem specification
  as defined by a call to \c GMRFLib_init_problem_hidden().  At output \a hidden_problem contains a sample
  or the log-density of the sample depending of the value of \a sample_flag.

  \param sample_flag If #GMRFLib_TRUE sample a hidden problem. If #GMRFLib_FALSE evaluate a hidden problem.
*/
int GMRFLib_doit_hidden(GMRFLib_hidden_problem_tp * hidden_problem, int sample_flag)
{
	/*
	 * sample (sample_flag=TRUE) or evaluate (sample_flag=FALSE) a hidden problem 
	 */
	int sub_n, i;
	GMRFLib_hidden_problem_tp *h = hidden_problem;

	if (!h)
		return GMRFLib_SUCCESS;

	/*
	 * save the current state, and switch to the state used for the construction 
	 */
	h->ran_current_state = GMRFLib_uniform_getstate(NULL);
	GMRFLib_EWRAP0(GMRFLib_uniform_setstate(h->ran_approx_state));

	sub_n = h->sub_graph->n;
	if (!sub_n)
		return GMRFLib_SUCCESS;

	for (i = h->hidden_par->stop_idx; i < sub_n; i++)
		h->sub_sample[i] = h->sample[h->map[i]] - h->sub_mean[i];

	for (i = sub_n - 1, h->sub_logdens = 0.0; i >= h->hidden_par->stop_idx; i--) {
		double ldens;

		GMRFLib_EWRAP0(GMRFLib_doit_hidden_i(h, sample_flag, i, &ldens));
		h->sub_logdens += ldens;
	}

	if (sample_flag)
		for (i = h->hidden_par->stop_idx; i < sub_n; i++)
			h->sample[h->map[i]] = h->sub_sample[i] + h->sub_mean[i];

	/*
	 * switch back to the users state 
	 */
	GMRFLib_EWRAP0(GMRFLib_uniform_setstate(h->ran_current_state));
	Free(h->ran_current_state);

	return GMRFLib_SUCCESS;
}
int GMRFLib_doit_hidden_i(GMRFLib_hidden_problem_tp * h, int sample_flag, int idx, double *ldens)
{
	int i, j, ii, jj, debug = 0, sub_n, bw, fidx;
	double cm, csd, x_low, x_high, hold;
	GMRFLib_gdens_tp *gdens = NULL;

	static GMRFLib_ConditionalFunc_arg_tp *args = NULL;
	static double **samples = NULL, *cmean0 = NULL, *cmean1 = NULL, *cmean = NULL, *cmode = NULL;
	static SampleScale_tp *sample_scale = NULL;
	static int num_sample = 0, len_sample = 0, n_antithetic = 0;

	sub_n = h->sub_graph->n;
	bw = h->sub_sm_fact.bandwidth;

	if (h->hidden_par->ngraph->nnbs[idx]) {
		if (h->hidden_par->ngraph->nnbs[idx] && h->hidden_par->ngraph->nbs[idx][0] < idx)
			fidx = h->hidden_par->ngraph->nbs[idx][0];	/* sorted increasingly */
		else
			fidx = idx;
		GMRFLib_ASSERT(fidx >= 0, GMRFLib_EPARAMETER);
	} else
		fidx = idx;

	if (debug)
		printf("%s: enter with idx=%1d fidx=%1d\n", __GMRFLib_FuncName, idx, fidx);

	/*
	 * setup workspace which is static & never free'ed FIXME: free this when calling ...free_hidden FIXME: dont need so
	 * large storage, only a subset will do. 
	 */
	if (h->hidden_par->nsample > num_sample || sub_n > len_sample || h->hidden_par->nantithetic > n_antithetic) {
		if (samples)
			for (i = 0; i < num_sample; i++)
				Free(samples[i]);
		if (sample_scale)
			for (i = 0; i < num_sample; i++)
				Free(sample_scale[i].s);
		Free(samples);
		Free(sample_scale);
		Free(cmean0);
		Free(cmean1);
		Free(cmean);
		Free(cmode);

		num_sample = h->hidden_par->nsample;
		len_sample = sub_n;
		n_antithetic = h->hidden_par->nantithetic;

		if (debug)
			printf("%s: malloc static storage for the samples\n", __GMRFLib_FuncName);

		if (num_sample) {
			samples = Calloc(num_sample, double *);

			for (i = 0; i < num_sample; i++) {
				samples[i] = Calloc(len_sample, double);
			}
		} else {
			samples = NULL;
		}

		cmean0 = Calloc(len_sample, double);
		cmean1 = Calloc(len_sample, double);
		cmean = Calloc(len_sample, double);
		cmode = Calloc(len_sample, double);

		if (num_sample) {
			sample_scale = Calloc(num_sample, SampleScale_tp);
			for (i = 0; i < num_sample; i++) {
				if (h->hidden_par->nantithetic > 0)
					sample_scale[i].n_sample_scale = h->hidden_par->nantithetic * GMRFLib_NUM_ANTITHETIC;
				else
					sample_scale[i].n_sample_scale = 1;
				sample_scale[i].s = Calloc(sample_scale[i].n_sample_scale, double);
			}
		} else {
			sample_scale = NULL;
		}
	}

	/*
	 * first we need to sample the Gaussians with zero conditional mean. cmean is the conditonal mean for a fixed value of
	 * x[idx], any choice will do. then we subtract cond.mean from the samples to obtain zero cond.mean ones 
	 */
	for (i = idx + 1; i < IMIN(idx + bw + 1, sub_n); i++) {
		cmean0[i] = cmean1[i] = cmode[i] = h->sub_sample[i];
	}
	for (i = idx - 1; i >= fidx; i--) {
		cmean0[i] = cmean1[i] = cmode[i] = 0.0;
	}
	cmean0[idx] = 0.0;				       /* fixed value */
	cmean1[idx] = 1.0;				       /* fixed value */
	GMRFLib_EWRAP0(GMRFLib_solve_lt_sparse_matrix_special(cmean0, &(h->sub_sm_fact), h->sub_graph, idx, fidx, 1));	/* no remap * * as remap=I */
	GMRFLib_EWRAP0(GMRFLib_solve_lt_sparse_matrix_special(cmean1, &(h->sub_sm_fact), h->sub_graph, idx, fidx, 1));

	/*
	 * if we're using the gaussian approximation, then we dont need any samples to estimate the correction factor. 
	 */
	if (!h->hidden_par->gaussapprox) {
		if (1 ||				       /* do this always, se below why */
		    h->hidden_par->neightype == GMRFLib_NEIGHTYPE_LINEAR || h->hidden_par->cmeanmode == GMRFLib_COND_MODE) {	/* FIXME: is this
																 * what I whant? */
			/*
			 * sample every index 
			 */
			for (i = 0; i < h->hidden_par->nsample; i++) {
				double nsqr, fac;

				for (j = idx + 1; j < IMIN(idx + bw + 1, sub_n); j++) {
					samples[i][j] = cmean0[j];
				}
				samples[i][idx] = 0.0;	       /* same value as above */

				if (h->hidden_par->nsample == 1 && h->hidden_par->nantithetic == 0) {
					/*
					 * this is special, then we use the conditional mean 
					 */
					for (j = idx - 1; j >= fidx; j--)
						samples[i][j] = 0.0;
				} else {
					/*
					 * this is the ordinary case, here we do as 'usual' 
					 */
					for (j = idx - 1, nsqr = 0.0; j >= fidx; j--) {
						samples[i][j] = GMRFLib_stdnormal();
						nsqr += SQR(samples[i][j]);
					}
					if (nsqr > 0) {
						fac = 1.0 / sqrt(nsqr);
						for (j = idx - 1, nsqr = 0.0; j >= fidx; j--)
							samples[i][j] *= fac;
					}
				}

				if (idx - fidx > 0) {
					double u = GMRFLib_uniform(), v;
					int kk;

					if (h->hidden_par->nantithetic > 0) {
						for (jj = 0, kk = 0; jj < h->hidden_par->nantithetic; jj++) {
							v = u + (double) jj / (double) h->hidden_par->nantithetic;
							if (v > 1.0)
								v -= 1.0;
							sample_scale[i].s[kk] = sqrt(gsl_cdf_chisq_Pinv(v, (double) (idx - fidx)));
							sample_scale[i].s[kk + 1] = -sample_scale[i].s[kk];
							sample_scale[i].s[kk + 2] = sqrt(gsl_cdf_chisq_Pinv(1.0 - v, (double) (idx - fidx)));
							sample_scale[i].s[kk + 3] = -sample_scale[i].s[kk + 2];

							kk += GMRFLib_NUM_ANTITHETIC;	/* has to be 4 */
						}
					} else
						sample_scale[i].s[0] = sqrt(gsl_cdf_chisq_Pinv(GMRFLib_uniform(), (double) (idx - fidx)));
				} else
					memset(sample_scale[i].s, 0, sample_scale[i].n_sample_scale * sizeof(double));

				GMRFLib_EWRAP0(GMRFLib_solve_lt_sparse_matrix_special(samples[i], &(h->sub_sm_fact), h->sub_graph, idx, fidx, 1));
				for (j = idx; j >= fidx; j--)
					samples[i][j] -= cmean0[j];
			}
		} else {				       /* == GMRFLib_NEIGHTYPE_GRAPH */

			/*
			 * sample just those we need 
			 */
			double *cov = NULL, *ccov = NULL, *sol = NULL, *chol_ccov = NULL;
			int nneig, ndim, cndim, *indexs = NULL;

			/*
			 * dette er ikke helt rett ettersom vi maa betinge mhp paa hele x i bandwidth, og ikke bare x[idx]!!!
			 * men vi kan finne L (hvor Lz=x) ved aa se paa sampler med z=0 og z_i=1 deretter. da finner vi vel
			 * kolonne 'i' i L eller noe.
			 * 
			 * profile'ing viser at ikke saa mye effort gaar med paa dtbsvspecial_, saa er det egentlig saa mye
			 * vits? dette er et realistisk eksempel:
			 * 
			 * Each sample counts as 0.01 seconds. % cumulative self self total time seconds seconds calls ms/call
			 * ms/call name 33.35 56.49 56.49 505138953 0.00 0.00 GMRFLib_loglFunc_wrapper 29.30 106.12 49.63
			 * 187392 0.26 0.81 GMRFLib_ConditionalFunc 27.03 151.90 45.78 505245459 0.00 0.00 loglik 3.22 157.36
			 * 5.46 2962186 0.00 0.00 gratio_ 1.91 160.60 3.24 40992 0.08 0.08 dtbsvspecial_
			 * 
			 */
			FIXME1("ENTER ILLEGAL CODE!!!");
			abort();

			for (j = 0, nneig = 0; j < h->hidden_par->ngraph->nnbs[idx]; j++)
				if (h->hidden_par->ngraph->nbs[idx][j] < idx)
					nneig++;
			ndim = nneig + 1;
			cndim = nneig;

			if (cndim > 0) {
				sol = Calloc(sub_n, double);
				cov = Calloc(ISQR(ndim), double);
				ccov = Calloc(ISQR(cndim), double);
				indexs = Calloc(ndim, int);

				indexs[0] = idx;
				for (j = 0, i = 1; j < h->hidden_par->ngraph->nnbs[idx]; j++)
					if (h->hidden_par->ngraph->nbs[idx][j] < idx)
						indexs[i++] = h->hidden_par->ngraph->nbs[idx][j];

				for (j = 0; j < ndim; j++) {
					memset(sol, 0, sub_n * sizeof(double));
					sol[indexs[j]] = 1.0;
					GMRFLib_EWRAP0(GMRFLib_solve_llt_sparse_matrix(sol, 1, &(h->sub_sm_fact), h->sub_graph));
					for (i = 0; i < ndim; i++)
						cov[j + i * ndim] = cov[i + j * ndim] = sol[indexs[i]];
				}
				for (j = 0; j < cndim; j++)
					for (i = 0; i < cndim; i++)
						ccov[i + cndim * j] = ccov[j + cndim * i] =
						    cov[(i + 1) + ndim * (j + 1)] - cov[i + 1] * cov[j + 1] / cov[0];

				GMRFLib_EWRAP0(GMRFLib_comp_chol_general(&chol_ccov, ccov, cndim, NULL, GMRFLib_ESINGMAT));

				for (i = 0; i < h->hidden_par->nsample; i++) {
					for (j = 0; j < cndim; j++)
						sol[j] = GMRFLib_stdnormal();
					for (ii = 0; ii < cndim; ii++) {
						for (jj = 0, hold = 0.0; jj <= ii; jj++)
							hold += chol_ccov[ii + jj * cndim] * sol[jj];
						samples[i][indexs[ii + 1]] = hold;
					}
				}

				Free(sol);
				Free(cov);
				Free(ccov);
				Free(indexs);
				Free(chol_ccov);
			}
		}
	}

	if (!args)
		args = Calloc(1, GMRFLib_ConditionalFunc_arg_tp);

	args->cmean = cmean;
	args->cmean0 = cmean0;
	args->cmean1 = cmean1;
	args->cmode = cmode;
	args->fidx = fidx;
	args->h = h;
	args->idx = idx;
	args->samples = samples;
	args->sample_scale = sample_scale;

	if (h->hidden_par->gaussapprox || (h->hidden_par->nsample == 0 && h->sub_d[idx] == 0.0)) {
		GMRFLib_EWRAP0(GMRFLib_comp_cond_meansd(&cm, &csd, idx, cmode, 1, &(h->sub_sm_fact), h->sub_graph));
	} else {
		if (h->hidden_par->cmeanmode == GMRFLib_COND_MODE) {
			GMRFLib_EWRAP0(GMRFLib_locate_cmode(cmode, idx, h, &cm, &csd));
		} else {
			GMRFLib_EWRAP0(GMRFLib_comp_cond_meansd(&cm, &csd, idx, cmode, 1, &(h->sub_sm_fact), h->sub_graph));
			cmode[idx] = cm;
		}
		x_low = cmode[idx] - h->hidden_par->range * csd;	/* low end of approx */
		x_high = cmode[idx] + h->hidden_par->range * csd;	/* high end of approx */

		if (1) {
			/*
			 * old 
			 */
			gdens = GMRFLib_gdens_Init(x_low, x_high, h->hidden_par->nresolution, h->hidden_par->nresolution,
						   h->hidden_par->norder, GMRFLib_ConditionalFunc, (void *) args);
		} else {
			/*
			 * new experimental code. note that the arguments are different!!! 
			 */
			gdens =
			    GMRFLib_gdens_InitNew(h->hidden_par->range, cmode[idx], csd, h->hidden_par->nresolution, GMRFLib_ConditionalFunc,
						  (void *) args);
		}

	}

	if (sample_flag) {
		void *state;

		/*
		 * ...oops, here we need to be a bit carful about the random stream. if we want to sample, we doit from the
		 * 'current' state, ie we need to switch to that one just here, and then return to the 'approx_state'
		 * afterwards. 
		 */
		state = GMRFLib_uniform_getstate(NULL);
		GMRFLib_EWRAP0(GMRFLib_uniform_setstate(h->ran_current_state));
		Free(h->ran_current_state);

		if (h->hidden_par->gaussapprox || (h->hidden_par->nsample == 0 && h->sub_d[idx] == 0.0)) {
			double zz;

			zz = GMRFLib_stdnormal();
			h->sub_sample[idx] = zz * csd + cm;
			*ldens = -0.5 * SQR(zz) - log(csd) - 0.5 * log(2.0 * M_PI);
		} else {
			GMRFLib_gdens_Sample(gdens, &(h->sub_sample[idx]), ldens);
		}

		/*
		 * then switch back. 
		 */
		h->ran_current_state = GMRFLib_uniform_getstate(NULL);
		GMRFLib_EWRAP0(GMRFLib_uniform_setstate(state));
		Free(state);
	} else {
		if (h->hidden_par->gaussapprox || (h->hidden_par->nsample == 0 && h->sub_d[idx] == 0.0)) {
			double zz;

			zz = (h->sub_sample[idx] - cm) / csd;
			*ldens = -0.5 * SQR(zz) - log(csd) - 0.5 * log(2.0 * M_PI);
		} else {
			GMRFLib_EWRAP0(GMRFLib_gdens_LDens(gdens, &(h->sub_sample[idx]), ldens));
		}
	}

	if (0) {
		/*
		 * for checking... 
		 */
		if (!(h->hidden_par->gaussapprox)) {
			double xp, ld;

			// for(xp=x_low;xp<x_high;xp += (x_high-x_low)/500.)
			for (xp = -0.5; xp <= 0.5; xp += 0.05) {
				GMRFLib_gdens_LDens(gdens, &xp, &ld);
				if (exp(ld) > 0.000001)
					printf("MARG %.16d %.18g %.18g %.18g %f %f\n", idx, xp, exp(ld),
					       exp(-0.5 * SQR((xp - cm) / csd)) / csd / pow(2.0 * M_PI, 0.5), cm, csd);
			}
			printf("%.6d\n", idx);
		}
	}

	/*
	 * special xpert-option for idx=n-1; copy the marginal density and the marginal mean to the user. 
	 */
	if (h->hidden_par->marginal_density && idx == h->sub_graph->n - 1) {
		*(h->hidden_par->marginal_density) = gdens;
		if (h->hidden_par->marginal_mean)
			*(h->hidden_par->marginal_mean) = h->sub_mean[idx];
		if (h->hidden_par->marginal_stdev)
			*(h->hidden_par->marginal_stdev) = csd;
	} else
		GMRFLib_gdens_Free(gdens);

	return GMRFLib_SUCCESS;
}
int GMRFLib_locate_cmode(double *cmode, int idx, GMRFLib_hidden_problem_tp * h, double *cond_mean, double *cond_stdev)
{
	/*
	 * locate the conditonal mode for 'idx' taking the data also into account.  if cond_mean then return the conditonal
	 * mean and similar with cond_stdev 
	 */
	double xopt, hivar, cm, csd, fac;
	int iter_max = 50;

	GMRFLib_EWRAP0(GMRFLib_comp_cond_meansd(&cm, &csd, idx, cmode, 1, &(h->sub_sm_fact), h->sub_graph));
	if (cond_mean)
		*cond_mean = cm;
	if (cond_stdev)
		*cond_stdev = csd;

	hivar = 0.5 / SQR(csd);
	xopt = cm;

	fac = 0.25;
	while (iter_max--) {
		double step = 1.0e-4, eps = 1.0e-4, f[3], xx[3], fl, fc, fr, deriv, dderiv;

		xx[0] = xopt - step;
		xx[1] = xopt;
		xx[2] = xopt + step;
		if (h->sub_d[idx])
			(*(h->loglFunc)) (f, xx, 3, idx, NULL, NULL, h->loglFunc_arg);
		else
			f[0] = f[1] = f[2] = 0.0;

		fl = -SQR(xx[0] - cm) * hivar + h->sub_d[idx] * f[0];
		fc = -SQR(xx[1] - cm) * hivar + h->sub_d[idx] * f[1];
		fr = -SQR(xx[2] - cm) * hivar + h->sub_d[idx] * f[2];

		deriv = (fr - fl) / 2.0;		       /* ...ignore 1/step */
		dderiv = (fr - 2.0 * fc + fl) / step;	       /* ...ignore 1/step */
		if (!ISZERO(dderiv)) {
			xopt -= DMIN(fac, 1.0) * deriv / dderiv;	/* ...because 1/step's cancels here */
			if (ABS(deriv / dderiv) < eps)
				break;
		}
		fac += 0.25;
	}

	if (iter_max > 0)
		cmode[idx] = xopt;
	else
		cmode[idx] = cm;
	if (0)
		GMRFLib_ASSERT(iter_max > 0, GMRFLib_EOPTNR);

	return GMRFLib_SUCCESS;
}
double GMRFLib_ConditionalFunc(double xn, void *arg)
{
	GMRFLib_ConditionalFunc_arg_tp *a;
	double ldens = 0.0, fac, cm, sum, tmp;
	int node, i, j, jj, k, kk, ns;

	static double *cond_stdev = NULL, *ldens_f = NULL;
	static int len_cond_stdev = 0, len_ldens_f = 0;
	static double *ff = NULL, *xx = NULL;
	static int len_xx = 0;

	a = (GMRFLib_ConditionalFunc_arg_tp *) arg;

	if (a->h->sub_graph->n > len_cond_stdev) {
		len_cond_stdev = a->h->sub_graph->n;	       /* FIXME: dont need that large storage */
		Free(cond_stdev);			       /* FIXME: how to free the static storage? */
		cond_stdev = Calloc(len_cond_stdev, double);
	}

	if (a->sample_scale)
		ns = a->sample_scale[0].n_sample_scale;	       /* they are all the same, pr.def */
	else
		ns = 0;

	if (a->h->hidden_par->nsample > len_ldens_f) {
		len_ldens_f = a->h->hidden_par->nsample * ns;
		Free(ldens_f);
		ldens_f = Calloc(len_ldens_f, double);
	}

	if (ns * a->h->hidden_par->nsample > len_xx) {
		Free(ff);
		Free(xx);
		len_xx = ns * a->h->hidden_par->nsample;
		ff = Calloc(len_xx, double);
		xx = Calloc(len_xx, double);
	}

	/*
	 * get cm and csd at idx. 
	 */
	GMRFLib_comp_cond_meansd(&cm, &cond_stdev[a->idx], a->idx, a->cmode, 1, &(a->h->sub_sm_fact), a->h->sub_graph);

	/*
	 * compute the conditonal mode or mean 
	 */

	if (a->h->hidden_par->cmeanmode == GMRFLib_COND_MODE) {
		a->cmode[a->idx] = xn;
		for (j = a->idx - 1; j >= a->fidx; j--)
			GMRFLib_locate_cmode(a->cmode, j, a->h, NULL, &cond_stdev[j]);

		if (0) {
			double ccmm;

			a->cmode[a->idx] = xn;
			for (j = a->idx - 1; j >= a->fidx; j--) {
				GMRFLib_locate_cmode(a->cmode, j, a->h, &ccmm, &cond_stdev[j]);
				printf("idx %d j %d cmean %f cmode %f dev %f\n", a->idx, j, ccmm, a->cmode[j],
				       (ccmm - a->cmode[j]) / cond_stdev[j]);
			}
		}
	}

	fac = (xn - cm) / cond_stdev[a->idx];
	if ((a->h->hidden_par->neightype == GMRFLib_NEIGHTYPE_LINEAR) || (a->h->hidden_par->cmeanmode == GMRFLib_COND_MODE)) {	/* need all of them 
																 */
		for (j = a->idx - 1; j >= a->fidx; j--)
			a->cmean[j] = a->cmean0[j] + fac * (a->cmean1[j] - a->cmean0[j]);
	} else {					       /* do not need all cmean's */

		if (0)
			for (j = a->idx + 1; j >= a->fidx; j--)
				a->cmean[j] = DBL_MAX;	       /* just a test... */
		for (j = 0; j < a->h->hidden_par->ngraph->nnbs[a->idx]; j++) {
			node = a->h->hidden_par->ngraph->nbs[a->idx][j];
			if (node < a->idx)
				a->cmean[node] = a->cmean0[node] + fac * (a->cmean1[node] - a->cmean0[node]);
		}
	}

	if (a->h->hidden_par->nsample * ns > 0) {
		memset(ldens_f, 0, a->h->hidden_par->nsample * ns * sizeof(double));

		if (a->h->hidden_par->cmeanmode == GMRFLib_COND_MODE) {	/* only needed then */
			for (k = 0; k < a->h->hidden_par->nsample; k++) {
				double factor, bb;

				kk = k * ns;
				for (j = a->idx - 1; j >= a->fidx; j--) {
					factor = -0.5 / SQR(cond_stdev[j]);
					bb = a->cmode[j] - a->cmean[j];

					for (jj = 0; jj < ns; jj++)
						ldens_f[kk + jj] += factor * ((2.0 * a->sample_scale[k].s[jj] * a->samples[k][j] + bb) * bb);
				}
			}
			if (0)
				for (k = 0; k < a->h->hidden_par->nsample; k++)
					for (jj = 0; jj < ns; jj++)
						printf("correction for cmode/mean %d %f\n", kk + jj, ldens_f[kk + jj]);
		}
		for (j = 0; j < a->h->hidden_par->ngraph->nnbs[a->idx]; j++) {
			double c_mode, c_mean, x_s, d_s;

			node = a->h->hidden_par->ngraph->nbs[a->idx][j];
			d_s = a->h->sub_d[node];

			if (node < a->idx && d_s) {
				c_mode = a->cmode[node];
				c_mean = a->cmean[node];

				/*
				 * WARNING: don't change the loop order of k and jj !!!!! 
				 */
				if (a->h->hidden_par->cmeanmode == GMRFLib_COND_MODE) {
					for (k = 0, kk = 0; k < a->h->hidden_par->nsample; k++) {
						x_s = a->samples[k][node];
						for (jj = 0; jj < ns; jj++)
							xx[kk++] = a->sample_scale[k].s[jj] * x_s + c_mode;	/* use cmode */
					}
				} else {
					for (k = 0, kk = 0; k < a->h->hidden_par->nsample; k++) {
						x_s = a->samples[k][node];
						for (jj = 0; jj < ns; jj++)
							xx[kk++] = a->sample_scale[k].s[jj] * x_s + c_mean;	/* use cmean */
					}
				}

				(*(a->h->loglFunc)) (ff, xx, a->h->hidden_par->nsample * ns, node, NULL, NULL, a->h->loglFunc_arg);
				for (k = 0, kk = 0; k < a->h->hidden_par->nsample; k++)
					for (jj = 0; jj < ns; jj++) {
						ldens_f[kk] += d_s * ff[kk];
						if (0)
							printf("ff[%1d] %f ds %f\n", kk, ldens_f[kk], d_s);
						kk++;
					}
			}
		}

		for (i = 0, sum = 0.0; i < a->h->hidden_par->nsample * ns; i++)
			sum += exp(ldens_f[i]);
		sum /= (ns * a->h->hidden_par->nsample);
		if (sum <= FLT_MIN)
			ldens = -FLT_MAX;
		else
			ldens = log(sum);
	} else {
		ldens = 0.0;
	}

	/*
	 * then at last, add the contribs at idx, if there is any 
	 */
	if (a->h->sub_d[a->idx])
		(*(a->h->loglFunc)) (&tmp, &xn, 1, a->idx, NULL, NULL, a->h->loglFunc_arg);
	else
		tmp = 0.0;

	ldens += -0.5 * SQR(xn - cm) / SQR(cond_stdev[a->idx]) + a->h->sub_d[a->idx] * tmp;

	if (0)
		printf("xn %f cm %f csd %f ldens %f \n", xn, cm, cond_stdev[a->idx], ldens);

	return ldens;
}
