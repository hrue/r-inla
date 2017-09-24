
/* hidden-approx.h
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
 *
 */

/*!
  \file hidden-approx.h
  \brief Typedefs and defines for \ref hidden-approx.c
*/

#ifndef __GMRFLib_HIDDEN_APPROX_H__
#define __GMRFLib_HIDDEN_APPROX_H__

#include <math.h>
#include <strings.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS					       /* empty */
#define __END_DECLS					       /* empty */
#endif

__BEGIN_DECLS

/**
 *  \brief One antithetic sample consists of GMRFLib_NUM_ANTITHETIC ones, which is
 * switching the sign, times, using \f$u\f$ and \f$1-u\f$ to determine the scale.
*/
#define GMRFLib_NUM_ANTITHETIC 4

/**
 * \brief Use the linear graph to determine the neigbours.
*/
#define GMRFLib_NEIGHTYPE_LINEAR 0

/**
 * \brief Use the folded graph to determine the neigbours.
*/
#define GMRFLib_NEIGHTYPE_GRAPH 1

/*
  Undocumented features
*/
#define GMRFLib_COND_MEAN 2
#define GMRFLib_COND_MODE 3

/*!
  \struct GMRFLib_hidden_param_tp hidden-approx.h
  \brief Holding information on how to set up the hidden approximation problem.

  \sa GMRFLib_default_hidden_par
 */

/*
 */
    typedef struct {

	/**
	 *  \brief Hold the neigbour graph. Computed internally. Must be initialised to \c NULL. 
	 */
	GMRFLib_graph_tp *ngraph;

	/**
	 *  \brief Range of the spline approximation. 
	 */
	double range;

	/**
	 *  \brief Step-length for Taylor expansion. 
	 */
	double step_len;
	int stencil;

	/**
	 *  \brief How many neigbhour neighbors to include when making the neighbour-graph. 
	 */
	int neighpar;

	/**
	 *  \brief Use the neighbours (see #GMRFLib_NEIGHTYPE_GRAPH) or just indices with lower index (see
	 * #GMRFLib_NEIGHTYPE_LINEAR). 
	 */
	int neightype;

	/**
	 *  \brief Order in the fitted log-spline: 1 or 2. 
	 */
	int norder;

	/**
	 *  \brief Number of regions in the spline. 
	 */
	int nresolution;

	/**
	 *  \brief How many samples to produce when estimating \f$\hat I\f$ . 
	 */
	int nsample;

	/**
	 *  \brief How many antithetic samples shall be constructed for each sample.
	 * 
	 * \sa #GMRFLib_NUM_ANTITHETIC 
	 */
	int nantithetic;

	/**
	 *  \brief Expand around the current value (see #GMRFLib_MODEOPTION_CURRENT) or the mode (see
	 * #GMRFLib_MODEOPTION_MODE). 
	 */
	int modeoption;

	/**
	 *  \brief If \c TRUE, use the Gaussian approximation instead of computing the non-Gaussian approximation. This choice 
	 * will override all other options. 
	 */
	int gaussapprox;

	/*
	 * Not documented.... 
	 */
	int cmeanmode;

	/*
	 * Undocumented xpert-options...
	 * 
	 * if we control the remap, and only sample from n-1 .... stop_idx, then we sample from a marginal of the full
	 * conditional. this have useful applications. these options are intented used togther, but does not. note that `remap' 
	 * is default NULL and stop_idx is default 0, ie sample from the full conditional. 
	 */
	int *remap;					       /* force to use this remap n */
	int stop_idx;					       /* sample/evalutate the marginal x_from n-1 to stop_idx, in the permuted world. */
	GMRFLib_gdens_tp **marginal_density;		       /* if not NULL, then store the LAST marginal as (*marginal) instead of FREE'ing it. */
	double *marginal_mean, *marginal_stdev;		       /* the density stored in `marginal' is adjusted for the marginal mean. hence to evaluate the density 
							        * in 'x' then marginal density must be evaluated at 'x-marginal_mean'. the marginal_stdev is
							        * returned [optionally] for convenice. */
} GMRFLib_hidden_param_tp;

/*!
  \struct GMRFLib_hidden_problem_tp hidden-approx.h
  \brief Specification of the hidden sampling problem.
 */
typedef struct {

	/**
	 *  \brief Mapping to the real world 
	 */
	int *map;

	/**
	 *  \brief The sample produced or the configuration to be evaluated. 
	 */
	double *sample;

	/**
	 *  \brief The point we expand around in real world, current configuration or the mode 
	 */
	double *sample_star;

	/**
	 *  \brief The sample produced on the sub_graph and with zero mean (Internal use only) 
	 */
	double *sub_sample;

	/**
	 *  \brief The mean on the sub_graph. 
	 */
	double *sub_mean;

	/**
	 *  \brief Internal use only 
	 */
	double *sub_b;

	/**
	 *  \brief Internal use only 
	 */
	double *sub_d;

	/**
	 *  \brief Internal use only 
	 */
	double *x_vec;

	/**
	 *  \brief The log-density of the sample produced or the configuration to be evaluated 
	 */
	double sub_logdens;

	/**
	 *  \brief Internal use only 
	 */
	void *ran_approx_state;

	/**
	 *  \brief Internal use only 
	 */
	void *ran_current_state;

	/**
	 *  \brief The Cholesky triangle (internal use only) 
	 */
	GMRFLib_sm_fact_tp sub_sm_fact;

	/**
	 *  \brief The subgraph (internal use only) 
	 */
	GMRFLib_graph_tp *sub_graph;

	/**
	 *  \brief Internal use only 
	 */
	GMRFLib_Qfunc_tp *sub_Qfunc;

	/**
	 *  \brief Internal use only 
	 */
	GMRFLib_Qfunc_arg_tp *sub_Qfunc_arg;

	/**
	 *  \brief Internal use only. 
	 */
	GMRFLib_logl_tp *loglFunc;

	/**
	 *  \brief Internal use only 
	 */
	void *loglFunc_arg;

	/**
	 *  \brief Parameters which control the approximation \sa GMRFLib_default_hidden_par 
	 */
	GMRFLib_hidden_param_tp *hidden_par;
} GMRFLib_hidden_problem_tp;

typedef struct {
	int n_sample_scale;
	double *s;
} SampleScale_tp;

typedef struct {
	int idx;
	int fidx;
	double *cmode;
	double *cmean;
	double *cmean0;
	double *cmean1;
	double **samples;
	SampleScale_tp *sample_scale;
	GMRFLib_hidden_problem_tp *h;
} GMRFLib_ConditionalFunc_arg_tp;

typedef struct {
	GMRFLib_logl_tp *loglFunc;
	void *loglFunc_arg;
	int *map;
	double *sub_mean;
	double *acoof;
	double *bcoof;
	double *ccoof;
} GMRFLib_logl_arg_tp;

int GMRFLib_default_hidden_par(GMRFLib_hidden_param_tp ** hidden_par);
int GMRFLib_loglFunc_wrapper(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int GMRFLib_init_problem_hidden(GMRFLib_hidden_problem_tp ** hidden_problem,
				double *x, double *b, double *c, double *mean,
				GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_args,
				char *fixed_value,
				double *d, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg,
				GMRFLib_optimize_param_tp * optpar, GMRFLib_hidden_param_tp * hidden_par);
int GMRFLib_init_problem_hidden_store(GMRFLib_hidden_problem_tp ** hidden_problem,
				      double *x, double *b, double *c, double *mean,
				      GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_args,
				      char *fixed_value,
				      double *d, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg,
				      GMRFLib_optimize_param_tp * optpar, GMRFLib_hidden_param_tp * hidden_par, GMRFLib_store_tp * store);
int GMRFLib_free_hidden(GMRFLib_hidden_problem_tp * hidden_problem);
int GMRFLib_sample_hidden(GMRFLib_hidden_problem_tp * hidden_problem);
int GMRFLib_evaluate_hidden(GMRFLib_hidden_problem_tp * hidden_problem);
int GMRFLib_doit_hidden(GMRFLib_hidden_problem_tp * hidden_problem, int flag);
int GMRFLib_doit_hidden_i(GMRFLib_hidden_problem_tp * h, int sample, int idx, double *ldens);
int GMRFLib_locate_cmode(double *cmode, int idx, GMRFLib_hidden_problem_tp * h, double *cond_mean, double *cond_stdev);
double GMRFLib_ConditionalFunc(double xn, void *arg);

__END_DECLS
#endif
