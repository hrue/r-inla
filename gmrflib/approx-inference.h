
/* approx-inference.h
 * 
 * Copyright (C) 2006-2022 Havard Rue
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
  \file approx-inference.h
  \brief Typedefs used to do approximative inference
*/

#ifndef __GMRFLib_APPROX_INFERENCE_H__
#define __GMRFLib_APPROX_INFERENCE_H__

#include <math.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>
#include "GMRFLib/GMRFLibP.h"

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

/*
 *
 */
typedef double GMRFLib_mfunc_tp(int thread_id, int node, void *argument);

typedef struct {
	GMRFLib_graph_tp *graph;
	GMRFLib_Qfunc_tp *Qfunc;
	void *Qfunc_arg;
	double diagonal;

	GMRFLib_mfunc_tp *mfunc;
	void *mfunc_arg;

	int n;
	int nreplicate;
	int ngroup;
} GMRFLib_bfunc2_tp;

typedef struct {
	int idx;
	GMRFLib_bfunc2_tp *bdef;			       /* doit like this, as many has the same 'bdef' */
} GMRFLib_bfunc_tp;

/** 
 * Available strategies
 */

typedef enum {

	/**
	 * \brief Use the GMRF-approximation
	 */
	GMRFLib_AI_STRATEGY_GAUSSIAN = 0,

	/**
	 * \brief Use the mean-corrected Gaussian approximation
	 */
	GMRFLib_AI_STRATEGY_MEANCORRECTED_GAUSSIAN,

	/**
	 * \brief Use the mean and skewness-corrected Gaussian approximation
	 */
	GMRFLib_AI_STRATEGY_MEANSKEWCORRECTED_GAUSSIAN,

	/**
	 *  \brief Fit a Spline-corrected Gaussian 
	 */
	GMRFLib_AI_STRATEGY_FIT_SCGAUSSIAN,

	/**
	 * \brief Adaptive strategy
	 */
	GMRFLib_AI_STRATEGY_ADAPTIVE
} GMRFLib_ai_strategy_tp;

typedef enum {
	GMRFLib_AI_VB_INVALID = 0,
	GMRFLib_AI_VB_MEAN,
	GMRFLib_AI_VB_VARIANCE
} GMRFLib_ai_vb_strategy_tp;

typedef enum {
	GMRFLib_VB_HESSIAN_STRATEGY_INVALID = 0,
	GMRFLib_VB_HESSIAN_STRATEGY_FULL,
	GMRFLib_VB_HESSIAN_STRATEGY_PARTIAL,
	GMRFLib_VB_HESSIAN_STRATEGY_DIAGONAL
} GMRFLib_ai_vb_hessian_strategy_tp;

#define VB_HESSIAN_STRATEGY_NAME(v_) ((v_) == GMRFLib_VB_HESSIAN_STRATEGY_FULL ? "full" : \
				      ((v_) == GMRFLib_VB_HESSIAN_STRATEGY_PARTIAL ? "partial" : \
				       ((v_) == GMRFLib_VB_HESSIAN_STRATEGY_DIAGONAL ? "diagonal" : "invalid")))

typedef enum {

	/**
	 * \brief Use  a grid strategy for integration
	 */
	GMRFLib_AI_INT_STRATEGY_GRID,

	/**
	 * \brief Use a CCD design for integration
	 */
	GMRFLib_AI_INT_STRATEGY_CCD,

	/**
	 * \brief Use an Empirical Bayes approach
	 */
	GMRFLib_AI_INT_STRATEGY_EMPIRICAL_BAYES,

	/**
	 * \brief Auto 
	 */
	GMRFLib_AI_INT_STRATEGY_AUTO,

	/**
	 * \brief USER (real scale)
	 */
	GMRFLib_AI_INT_STRATEGY_USER,

	/**
	 * \brief USER_STD (std scale)
	 */
	GMRFLib_AI_INT_STRATEGY_USER_STD,

	/**
	 * \brief USER_EXPERT (expert option: the density and the integration weights in the real scale)
	 */
	GMRFLib_AI_INT_STRATEGY_USER_EXPERT,

	/**
	 * \brief USER_PART2 (expert option: this is part2 of 'twostage' for which the weights are the corrections to 'log_dens')
	 */
	GMRFLib_AI_INT_STRATEGY_USER_PART2
} GMRFLib_ai_int_strategy_tp;

/** 
 * Types of linear approximations: \f$\log\pi(x_{-i}|x_i,\theta,y) \approx a x_i + \mbox{constant}\f$ 
 */
typedef enum {

	/**
	 * \brief Do not use a linear approximation, compute \f$\log\pi(x_{-i}|x_i,\theta,y)\f$ exact (default)
	 */
	GMRFLib_AI_LINEAR_CORRECTION_OFF = 0,		       /* MUST BE ZERO! */

	/**
	 * \brief Compute the derivative by correcting unconditional variances. This is the fast option.
	 */
	GMRFLib_AI_LINEAR_CORRECTION_FAST,

	/**
	 * \brief Compute the (large scale) derivative using a central difference approximation (large h), involves two
	 * factorisations.
	 */
	GMRFLib_AI_LINEAR_CORRECTION_CENTRAL_DIFFERENCE
} GMRFLib_ai_linear_correction_tp;

/**
 * Types of interpolators for computing marginals for each hyperparamter
 */
typedef enum {

	/**
	 * \brief Chose interpolation type based on the integration strategy
	 */
	GMRFLib_AI_INTERPOLATOR_AUTO = 0,

	/**
	 * \brief Use the nearest point only
	 */
	GMRFLib_AI_INTERPOLATOR_NEAREST,

	/**
	 * \brief Linear interpolation using the dim+1 nearest points
	 */
	GMRFLib_AI_INTERPOLATOR_LINEAR,

	/**
	 * \brief Quadratic interpolation using the dim+1 nearest points
	 */
	GMRFLib_AI_INTERPOLATOR_QUADRATIC,

	/**
	 * \brief Linear interpolation using weighted distance
	 */
	GMRFLib_AI_INTERPOLATOR_WEIGHTED_DISTANCE,

	/**
	 * \brief Special interpolation for the CCD integration scheme (analytic approximation)
	 */
	GMRFLib_AI_INTERPOLATOR_CCD,

	/**
	 * \brief Special interpolation for the CCD integration scheme (numerical integration)
	 */
	GMRFLib_AI_INTERPOLATOR_CCD_INTEGRATE,

	/**
	 * \brief Special interpolation for the GRID integration scheme and when hessian_force_diagonal is TRUE
	 */
	GMRFLib_AI_INTERPOLATOR_GRIDSUM,

	/**
	 * \brief Use the plain Gaussian
	 */
	GMRFLib_AI_INTERPOLATOR_GAUSSIAN
} GMRFLib_ai_interpolator_tp;

#define INTERPOLATOR_NAME(interp) ((interp) == GMRFLib_AI_INTERPOLATOR_GAUSSIAN ? "Gaussian" : \
				   ((interp) == GMRFLib_AI_INTERPOLATOR_WEIGHTED_DISTANCE ? "Weigted distance" : \
				    ((interp) == GMRFLib_AI_INTERPOLATOR_NEAREST ? "Nearest" : \
				     ((interp) == GMRFLib_AI_INTERPOLATOR_LINEAR ? "Linear" : \
				      ((interp) == GMRFLib_AI_INTERPOLATOR_QUADRATIC ? "Quadratic" : \
				       ((interp) == GMRFLib_AI_INTERPOLATOR_CCD ? "CCD" : \
					((interp) == GMRFLib_AI_INTERPOLATOR_CCD_INTEGRATE ? "CCDIntegrate" : \
					 ((interp) == GMRFLib_AI_INTERPOLATOR_GRIDSUM ? "GRIDSUM" : \
					  ((interp) == GMRFLib_AI_INTERPOLATOR_AUTO ? "Auto" : "Unknown: ERROR")))))))))

typedef enum {

	/**
	 * \brief BFGS(2) implementation in GSL
	 */
	GMRFLib_AI_OPTIMISER_GSL = 2,

	/**
	 * \brief The default choice
	 */
	GMRFLib_AI_OPTIMISER_DEFAULT
} GMRFLib_ai_optimiser_tp;

#define GMRFLib_AI_OPTIMISER_NAME(opt) \
	((opt) == GMRFLib_AI_OPTIMISER_GSL ? "GSL-BFGS2" : \
	 ((opt) == GMRFLib_AI_OPTIMISER_DEFAULT ? "DEFAULT METHOD" : "unknown!!!"))

/**
 * \brief Parameters for doing approximate inference
 */
typedef struct {

	int preopt_mode;

	/**
	 * \brief The stategy used to compute marginals 
	 */
	GMRFLib_ai_strategy_tp strategy;

	int adapt_max;
	int adapt_len;
	GMRFLib_ai_strategy_tp *adapt_strategy;

	/**
	 *  \brief Fast mode? If TRUE, compromise accurancy for the sake of speed. Usually OK.
	 */
	int fast;

	/**
	 * \brief Gaussian data?  If TRUE, then we know that all data is Gaussian
	 */
	int gaussian_data;

	/**
	 * \brief Type of linear correction by approximating \f$\log\pi(x_{-i}|x_i,\theta,y)\f$ linearly wrt \f$x_i\f$.
	 *
	 * This option is default OFF.
	 */
	GMRFLib_ai_linear_correction_tp linear_correction;

	/**
	 *  \brief The number of points evaluated when computing the Laplace-approximation 
	 */
	int n_points;

	/**
	 * \brief The step-length \f$h\f$ in discrete Taylor-approximation
	 *
	 * Step length in the (discrete) computation of a Taylor expansion or second order approximation of the log-likelihood
	 * around a point <em> \b x_0</em> (CG and NR).  \n\n
	 */
	double step_len;

	/**
	 * \brief Number of points in the stencil: 3, 5 or 7, to compute numerical derivaties
	 *
	 */
	int stencil;

	/**
	 *  \brief Ignore nodes where the derivative in the conditional mean is less than \a eps. 
	 */
	double cutoff;

	/**
	 * \brief Write log when integrating the hyperparameter
	 */
	FILE *fp_log;

	/**
	 * \brief Write info about the marginal for the hyperparameters
	 */
	FILE *fp_hyperparam;

	/**
	 * \brief  The integration strategy.
	 */
	GMRFLib_ai_int_strategy_tp int_strategy;

	/**
	 * \brief  The design, if strategy is _USER or _USER_STD
	 */
	GMRFLib_design_tp *int_design;

	/**
	 * \brief The scaling for \c GMRFLib_AI_INT_STRATEGY_CCD, must be > 1.
	 *
	 * The scaling used is \f$f = f_0 \sqrt{d}\f$ where \f$d\f$ is the number of hyperparameters.
	 */
	double f0;

	/**
	 * \brief The steplength when integrating the hyperparameters in N(0,1) scale. (Only for \c GMRFLib_AI_INT_STRATEGY_GRID)
	 */
	double dz;

	/** 
	 * \brief The maximal difference in log-density before skipping. (Only for \c GMRFLib_AI_INT_STRATEGY_GRID)
	 */
	double diff_log_dens;

	/**
	 * \brief Skip configurations based on ordering? (Only for \c GMRFLib_AI_INT_STRATEGY_GRID)
	 */
	int skip_configurations;

	/**
	 * \brief Adjust the integration weights? (Only for \c GMRFLib_AI_INT_STRATEGY_GRID)
	 */
	int adjust_weights;

	/**
	 * \brief Use forward finite difference to compute the gradient?
	 *
	 * If TRUE then use finite difference to compute the gradient. If FALSE, then use finite central difference which is
	 * more accurate but more expensive.
	 */
	int gradient_forward_finite_difference;

	/**
	 * \brief Use forward finite difference to compute the Hessian?
	 *
	 * If TRUE then use finite difference to compute the Hessian. If FALSE, then use finite central difference which is more
	 * accurate but more expensive.
	 */
	int hessian_forward_finite_difference;

	/**
	 * \brief The step-length in the finite difference calculations for the gradient.
	 */
	double gradient_finite_difference_step_len;

	/**
	 * \brief The step-length in the finite difference calculations for the Hessian.
	 */
	double hessian_finite_difference_step_len;

	/**
	 * \brief Setting this option to \c GMRFLib_TRUE forces the Hessian matrix to be diagonal
	 */
	int hessian_force_diagonal;

	/**
	 * \brief Compute and display effective number of parameters? (Only if \c fp_log is non-NULL)
	 */
	int compute_nparam_eff;

	/**
	 * \brief Do a MC-test for the error?  (Only if \c fp_log is non-NULL)
	 *
	 * If \c do_MC_error_check is positive, then use the defaul number of samples in the MC test. If \c do_MC_error_check is
	 * negative, use then use this many (removing the sign) samples in the MC test.
	 */
	int do_MC_error_check;

	/**
	 * \brief Type of interpolator used to compute the marginal of each hyperparameter, one of \c GMRFLib_ai_interpolator_tp
	 */
	GMRFLib_ai_interpolator_tp interpolator;

	/**
	 * \brief Type of optimiser to use
	 */
	GMRFLib_ai_optimiser_tp optimiser;

	/**
	 * \brief Run the optimiser twice by restarting the optmiser at the first found solution
	 */
	int restart;

	/**
	 * \brief GSL parameter tolerance (something with linesearch, recommended 0.1)
	 */
	double gsl_tol;

	/**
	 * \brief GSL parameter epsg. Stopping parameter for |grad f|.
	 */
	double gsl_epsg;

	/**
	 * \brief GSL parameter epsf. Stopping parameter for |f|.
	 */
	double gsl_epsf;

	/**
	 * \brief GSL parameter epsg. Stopping parameter for |x|.
	 */
	double gsl_epsx;

	/**
	 * \brief GSL parameter step_size. Size of the initial step in the line-search
	 */
	double gsl_step_size;

	/**
	 * \brief Gaussian approximation optmiser parameter: abserr_func
	 */
	double optpar_abserr_func;

	/**
	 * \brief Gaussian approximation optimiser parameter: abserr_func
	 */
	double optpar_abserr_step;

	/**
	 * \brief Gaussian approximation optimiser parameter: fp
	 */
	FILE *optpar_fp;

	/**
	 * \brief Step factor for the Newton Raphson algorithm: nr_step_factor
	 */
	double optpar_nr_step_factor;

	/**
	 * \brief A flag to say that the initial values of the hyperparameters, are the known mode
	 */
	int mode_known;

	/**
	 * \brief Accepted limit for computing the CPO-density
	 *
	 * How large difference in log-scale of the CPO-density, between the mode and the value at the border for bing classified as ok.
	 */
	double cpo_req_diff_logdens;

	/**
	 * \brief Use the stupid search algorithm
	 */
	int stupid_search_mode;

	/**
	 * \brief Maximum number of iterations when for the stupid search algorithm
	 */
	int stupid_search_max_iter;

	/**
	 * \brief The scale-factor for the stupid-search algorithm
	 */
	double stupid_search_factor;

	/**
	 * \brief Flag for manual-cpo calculation. (Expert use only.)
	 */
	int cpo_manual;

	/**
	 * \brief Maximum function evaluations for numerical integration (hyperparameters).
	 */
	int numint_max_fn_eval;

	/**
	 * \brief Relative error for numerical integration (hyperparameters).
	 */
	double numint_rel_err;

	/**
	 * \brief Absolute error for numerical integration (hyperparameters).
	 */
	double numint_abs_err;

	/**
	 * \brief Minimum value for the hessian from the likelihood used in the optimisation.
	 */
	double cmin;

	/**
	 * \brief The strategy for the linear term, either 0 (skip) or 1 (keep).
	 */
	int b_strategy;

	/**
	 * \brief Enable VB-correct?
	 */
	int vb_enable;

	/**
	 * \brief VB correction (0=mean)
	 */
	int vb_strategy;

	/**
	 * \brief Verbose while VB-correct?
	 */
	int vb_verbose;

	/**
	 * \brief Number of max iterations (>= 0)
	 */
	int vb_iter_max;

	/**
	 * \brief Emergency option
	 */
	double vb_emergency;

	/**
	 * \brief update_hessian
	 */
	int vb_hessian_update;

	/**
	 * \brief strategy for hessian
	 */
	GMRFLib_ai_vb_hessian_strategy_tp vb_hessian_strategy;

	/**
	 * \brief N-limit when to enable a f() component
	 */
	int vb_f_enable_limit_mean;
	int vb_f_enable_limit_variance;

	/**
	 * \brief List of nodes to VB correct for, if any.
	 */
	char *vb_nodes_mean;
	char *vb_nodes_variance;

	/** 
	 * Try to be smart when optimizing in INLA?   
	 */
	int optimise_smart;

	/**
         * Use directions when optimising?
	 */
	int optimise_use_directions;
	gsl_matrix *optimise_use_directions_m;

	// use the experimental version?
	int improved_simplified_laplace;

	int parallel_linesearch;
	int fixed_mode;
} GMRFLib_ai_param_tp;

/**
 * \brief Store computations to improve the speed
 */
typedef struct {

	/**
	 *\brief Store computations for the marginal of the hyperparameters
	 */
	GMRFLib_store_tp *store;

	/**
	 *\brief Store the previously computed mode
	 */
	double *mode;

	/**
	 *\brief Store the Gaussian approximation for the latent field
	 */
	GMRFLib_problem_tp *problem;
	int nc_orig;
	double *aa;
	double *bb;					       /* the 'bb' and 'cc' array */
	double *cc;

	int Npred;

	/**
	 *\brief Store indices for the datapoints
	 */
	int nd;

	/**
	 *\brief Store indices for the datapoints
	 */
	int *d_idx;

	/**
	 *\brief Store standard deviations for the Gaussian approximation
	 */
	double *stdev;

	/**
	 *\brief Store the correction terms if \a use_linear_correction = TRUE
	 */
	double *correction_term;			       /* hold the correction term for the use_linear_correction */

	/**
	 *\brief Store the third (and optionally fourth) derivative of the log likelihood
	 */
	double *derivative3;
	double *derivative4;
	int *correction_idx;				       /* hold the idx's for the corrections, ie those with d[i] != 0 */
	int nidx;					       /* number of those */
} GMRFLib_ai_store_tp;

#define GMRFLib_AI_STORE_NEFF_NOT_COMPUTED (-1.23456789)

/**
 * \brief Store the integrated likelihood for the model
 */
typedef struct {

	/**
	   \brief integrated likehood computed via numerical integration
	*/
	double marginal_likelihood_integration;

	/**
	   \brief integrated likehood computed from a Gaussian approximation
	*/
	double marginal_likelihood_gaussian_approx;
} GMRFLib_ai_marginal_likelihood_tp;

/**
 * \brief Results of the DIC computations
 */
typedef struct {

	/**
	 * \brief Mean of the devianace  
	 */
	double mean_of_deviance;

	/**
	 * \brief Mean of the devianace  (saturated)
	 */
	double mean_of_deviance_sat;

	/**
	 * \brief Devianace of the mean  
	 */
	double deviance_of_mean;

	/**
	 * \brief Devianace of the mean  (saturated)
	 */
	double deviance_of_mean_sat;

	/**
	 * \brief The effective number of parameters
	 */
	double p;

	/**
	 * \brief The DIC value
	 */
	double dic;

	/**
	 * \brief The DIC value (saturated)
	 */
	double dic_sat;

	/**
	 * \brief The length of the E(Deviance) contribution
	 */
	int n_deviance;

	/**
	 * \brief The E(deviance) contribution
	 */
	double *e_deviance;

	/**
	 * \brief The E(deviance) contribution (saturated)
	 */
	double *e_deviance_sat;

	/**
	 * \brief The deviance(E) contribution
	 */
	double *deviance_e;

	/**
	 * \brief The deviance(E) contribution (saturated)
	 */
	double *deviance_e_sat;

} GMRFLib_ai_dic_tp;

/**
 *\brief A template computing all terms in <b>(GMRF-35)</b> which are constant with respect to \f$\mbox{\boldmath$x$}\f$ but
 * depend on \f$ \mbox{\boldmath$\theta$}\f$.
 */
typedef double GMRFLib_ai_log_extra_tp(int thread_id, double *hyperparam, int nhyper, void *arg);

typedef struct {
	int nhyper;
	int idx;
	int hyper_count;
	int return_log;
	double *hyper_z;
	double *hyper_ldens;
	double theta_fixed;
	double *theta_mode;
	gsl_vector *sqrt_eigen_values;
	gsl_matrix *eigen_vectors;
	double dz;
	double *stdev_corr_pos;
	double *stdev_corr_neg;

	/*
	 * for internal use 
	 */
	double *z, *theta;
	GMRFLib_ai_interpolator_tp interpolator;
} GMRFLib_ai_integrator_arg_tp;

/**
 *   \brief The type of the cpo-object returned by \c GMRFLib_INLA().
 */
typedef struct {

	/**
	 * \brief Total number of nodes
	 */
	int n;

	/**
	 * The CPO-value, if cpo_value[i] is non-NULL, then cpo_value[i] points to the CPO value.
	 */
	double **value;

	/**
	 * The arithmetic mean of the cpo-values (\f$\sum_i p_i / n\f$)
	 */
	double mean_value;

	/**
	 * The geometrical mean of the cpo-values (\f$ exp(\sum_i \log(p_i)/n) \f$)
	 */
	double gmean_value;

	/**
	 * The PIT-value, if pit_value[i] is non-NULL, then pit_value[i] points to the PIT value.  PIT = Prob(y.NEW < y.OBS |
	 * all data but not y.OBS), for each i.  These probabilities can be used to detect surprising observations.
	 */
	double **pit_value;

	/**
	 * Failure indicator; if 0 then all seems ok, if 1 then all is very bad
	 */
	double **failure;
} GMRFLib_ai_cpo_tp;

typedef enum {
	GMRFLib_GCPO_BUILD_STRATEGY_POSTERIOR = 0,
	GMRFLib_GCPO_BUILD_STRATEGY_PRIOR = 1
} GMRFLib_gcpo_build_strategy_tp;

#define GMRFLib_GCPO_BUILD_STRATEGY_NAME(arg_) ((arg_) == GMRFLib_GCPO_BUILD_STRATEGY_POSTERIOR ? "Posterior" : \
						 ((arg_) == GMRFLib_GCPO_BUILD_STRATEGY_PRIOR ? "Prior" : \
						  "UNKNOWN"))

typedef struct {
	int num_level_sets;
	int size_max;
	int ngroups;
	int verbose;
	int correct_hyperpar;
	double epsilon;
	double prior_diagonal;
	GMRFLib_gcpo_build_strategy_tp build_strategy;	       /* 0=posterior, 1=prior (see above) */
	GMRFLib_idxval_tp **groups;
	GMRFLib_idx_tp *selection;

	int remove_fixed;
	GMRFLib_str_tp *remove;				       /* only one of these can be !NULL */
	GMRFLib_str_tp *keep;

	// this is to build the strategy="prior" corrections. these are just copy of ptr's, do not free!
	int idx_tot;
	int *idx_start;
	int *idx_n;
	char **idx_tag;
} GMRFLib_gcpo_param_tp;

typedef struct {
	int n;
	double *value;
	double *kld;
	double *mean;
	double *sd;
	GMRFLib_idxval_tp **groups;
} GMRFLib_gcpo_tp;

typedef struct {
	int Npred;
	GMRFLib_idxval_tp **groups;
	GMRFLib_idx2_tp **missing;
	GMRFLib_idxsubmat_vector_tp **missing2;
} GMRFLib_gcpo_groups_tp;

typedef struct {
	GMRFLib_idxval_tp *idxs;			       /* list of nodes in the matrix, sorted */
	gsl_matrix *cov_mat;				       /* the covariance matrix */
	double value;
	double lpred_mean;
	double lpred_sd;
	double kld;
	double marg_theta_correction;
	int node_min;					       /* min(nodes) */
	int node_max;					       /* max(nodes) */
	int idx_node;					       /* the index for the central node, 'i' in 'gcpo[i]', within cov_mat */
} GMRFLib_gcpo_elm_tp;

/**
 *   \brief The type of the cpo-object returned by \c GMRFLib_INLA().
 */
typedef struct {

	/**
	 * \brief Total number of nodes
	 */
	int n;

	/**
	 * The PO-value, if po_value[i] is non-NULL, then po_value[i] points to the PO value.
	 */
	double **value;

} GMRFLib_ai_po_tp;

typedef struct {
	int dens_count;
	double log_posterior;				       /* May have been adjusted for integration weight */
	double log_posterior_orig;			       /* Not adjusted for integration weight */
	double *theta;					       /* */
	double *mean;					       /* mean */
	double *improved_mean;				       /* improved mean */
	double *skewness;				       /* the skewness in the skew-normal=E[((x-mu)/sd)^3] */
	double *Q;					       /* the Q_ij-values */
	double *Qinv;					       /* the Qinv_ij-values */
	double *Qprior;
} GMRFLib_store_config_tp;

typedef struct {
	int n;						       /* size of the graph */
	int nz;						       /* size of the number of unique elements of Q */
	int ntheta;					       /* */
	int *i;						       /* i-values in Qij */
	int *j;						       /* j-values in Qij */
	int nconfig;					       /* number of configurations */
	GMRFLib_constr_tp *constr;			       /* */
	GMRFLib_store_config_tp **config;		       /* the configurations */
} GMRFLib_store_configs_tp;

typedef struct {
	double log_posterior;				       /* May have been adjusted for integration weight */
	double log_posterior_orig;			       /* Not adjusted for integration weight */
	double *theta;					       /* */
	double *mean;					       /* mean */
	double *improved_mean;				       /* improved mean */
	double *skewness;				       /* the skewness in the skew-normal=E[((x-mu)/sd)^3] */
	double *Q;					       /* the Q_ij-values */
	double *Qinv;					       /* the Qinv_ij-values */
	double *Qprior;
	double *cpodens_moments;
	double *gcpodens_moments;
} GMRFLib_store_config_preopt_tp;

typedef struct {
	int mnpred;
	int Npred;
	int n;						       /* size of the graph */
	int nz;						       /* size of the number of unique elements of Q */
	int prior_nz;					       /* size of the number of unique elements of Q */
	int ntheta;					       /* */
	int *i;						       /* i-values in Qij */
	int *j;						       /* j-values in Qij */
	int *iprior;					       /* i-values in Qij */
	int *jprior;					       /* j-values in Qij */
	int nconfig;					       /* number of configurations */
	GMRFLib_constr_tp *constr;			       /* */
	GMRFLib_store_config_preopt_tp **config;	       /* the configurations */
	GMRFLib_matrix_tp *A;
	GMRFLib_matrix_tp *pA;
} GMRFLib_store_configs_preopt_tp;


typedef struct {
	int nhyper;
	double *cov_m;

	double *eigenvalues;				       /* Need also the eigen-stuff as the corrections depends on the sign of the
							        * eigenvectors. */
	double *eigenvectors;
	double *stdev_corr_pos;
	double *stdev_corr_neg;

	double log_posterior_mode;

	/*
	 * [0] is the preparation in INLA
	 * [1] is the optimisation
	 * [2] is the integration
	 * [3] is the post-processing part including computing
	 *     marginals for the hyperparameters.
	 */
	double wall_clock_time_used[4];

	int len_reordering;
	int *reordering;

	int compute_corr_lin;
	double *corr_lin;				       /* correlation of the lincombs (derived only) */
	double *cov_lin;				       /* covariance of the lincombs (derived only) */

	int mode_status;				       /* 0 for ok, 1 not ok. */
	int nfunc;					       /* number of function evaluations */

	GMRFLib_matrix_tp *opt_directions;

	GMRFLib_store_configs_tp **configs;		       /* configs[id][...] */
	GMRFLib_store_configs_preopt_tp **configs_preopt;      /* configs[id][...] */
} GMRFLib_ai_misc_output_tp;

typedef struct {
	size_t nhyper;
	size_t nconfig;
	size_t *idx_mapping;
	size_t idx_next;
	int *configurations;

	char all_out;
	char *out;
	double diff_log_dens;
} GMRFLib_ai_pool_tp;

typedef struct {
	int first_nonzero;				       /* first nonzero idx = min(idx) */
	int last_nonzero;				       /* last nonzero idx = max(idx) */
	int first_nonzero_mapped;			       /* first nonzero idx of L^-1 a = min(remap(idx)). automatically added */
	int last_nonzero_mapped;			       /* last nonzero idx of L^-1 a. automatically added */
} GMRFLib_lc_tinfo_tp;

typedef struct {
	int n;						       /* length */
	int *idx;					       /* list of indices */
	double *weight;
	GMRFLib_lc_tinfo_tp *tinfo;			       /* thread-info */
} GMRFLib_lc_tp;

typedef struct {
	int i;
	int j;
} GMRFLib_lc_ij_tp;

typedef struct {
	double *stdev_corr_neg;
	double *stdev_corr_pos;
	gsl_vector *sqrt_eigen_values;
	gsl_matrix *eigen_vectors;
} GMRFLib_userfunc2_arg_tp;

typedef struct {
	double *stdev_corr_neg;
	double *stdev_corr_pos;
	gsl_vector *sqrt_eigen_values;
	gsl_matrix *eigen_vectors;
} GMRFLib_userfunc3_arg_tp;

typedef enum {
	GMRFLib_TRANSFORM_FORWARD = 1,			       /* same as in inla.h */
	GMRFLib_TRANSFORM_BACKWARD = 2,			       /* same as in inla.h */
	GMRFLib_TRANSFORM_DFORWARD = 3,			       /* same as in inla.h */
	GMRFLib_TRANSFORM_INCREASING = 4		       /* same as in inla.h */
} GMRFLib_transform_func_arg_tp;

typedef double GMRFLib_transform_func_tp(int thread_id, double arg, GMRFLib_transform_func_arg_tp typ, void *param, double *cov);

typedef struct {
	GMRFLib_transform_func_tp *func;
	void *arg;
	double *cov;
} GMRFLib_transform_array_func_tp;

typedef struct {
	// c[0] + c[1]*x + 1/2*c[2]*x^2
	double coofs[3];
} GMRFLib_vb_coofs_tp;

#define GMRFLib_AI_POOL_GET 1
#define GMRFLib_AI_POOL_SET 2

#include "GMRFLib/pre-opt.h"

int GMRFLib_ai_pool_free(GMRFLib_ai_pool_tp * pool);
int GMRFLib_ai_pool_get(GMRFLib_ai_pool_tp * pool, int *iz, size_t *idx);
int GMRFLib_ai_pool_init(GMRFLib_ai_pool_tp ** pool, GMRFLib_ai_param_tp * ai_par, int nhyper);
int GMRFLib_ai_pool_intern(GMRFLib_ai_pool_tp * pool, int *iz, size_t *idx, double logdens, int action);
int GMRFLib_ai_pool_set(GMRFLib_ai_pool_tp * pool, size_t idx, double logdens);
int GMRFLib_pool_cmp(const void *a, const void *b);
int GMRFLib_pool_cmp1(const void *a, const void *b);

int GMRFLib_ai_marginal_one_hyperparamter(GMRFLib_density_tp ** density, int idx, int nhyper, int hyper_count, double *hyper_z,
					  double *hyper_ldens, double *theta_mode, gsl_vector * sqrt_eigen_values,
					  gsl_matrix * eigen_vectors, double *std_stdev_theta, double dz,
					  double *stdev_corr_pos, double *stdev_corr_neg, GMRFLib_ai_interpolator_tp interpolator,
					  GMRFLib_ai_param_tp * ai_par, double *covmat);
double GMRFLib_ai_integrator_func(unsigned ndim, const double *x, void *arg);
double GMRFLib_interpolator_linear(int ndim, int nobs, double *x, double *xobs, double *yobs, void *arg);
double GMRFLib_interpolator_quadratic(int ndim, int nobs, double *x, double *xobs, double *yobs, void *arg);
double GMRFLib_interpolator_wdistance(int ndim, int nobs, double *x, double *xobs, double *yobs, void *arg);
double GMRFLib_interpolator_ccd(int ndim, int nobs, double *x, double *xobs, double *yobs, void *arg);
double GMRFLib_interpolator_distance(int ndim, double *x, double *xx);
double GMRFLib_interpolator_distance2(int ndim, double *x, double *xx);

int GMRFLib_default_ai_param(GMRFLib_ai_param_tp ** aipar);
int GMRFLib_print_ai_param(FILE * fp, GMRFLib_ai_param_tp * ai_par);
int GMRFLib_ai_marginal_hyperparam(int thread_id, double *logdens,
				   double *x, double *b, double *c, double *mean, double *d,
				   GMRFLib_logl_tp * loglFunc, void *loglFunc_arg,
				   GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg,
				   GMRFLib_constr_tp * constr, GMRFLib_ai_param_tp * aipar, GMRFLib_ai_store_tp * store,
				   GMRFLib_preopt_tp * preopt);
int GMRFLib_ai_log_posterior(int thread_id, double *logdens,
			     double *x, double *b, double *c, double *mean, double *d,
			     GMRFLib_logl_tp * loglFunc, void *loglFunc_arg,
			     GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg, GMRFLib_constr_tp * constr);
int GMRFLib_ai_log_posterior_restricted(int thread_id, double *logdens, double *x, double *x_mode, double *x_gradient, double delta, double *b,
					double *c, double *mean, double *d, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg,
					GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg,
					GMRFLib_constr_tp * constr, GMRFLib_graph_tp * subgraph, GMRFLib_ai_store_tp * ai_store,
					int *node_map, double *ql);
int GMRFLib_ai_marginal_hidden(int thread_id, GMRFLib_density_tp ** density, GMRFLib_density_tp ** cpo_density,
			       int lookup_tables,
			       int idx, double *x, double *b, double *c, double *mean, double *d,
			       GMRFLib_logl_tp * loglFunc, void *loglFunc_arg,
			       GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg,
			       GMRFLib_constr_tp * constr, GMRFLib_ai_param_tp * ai_par, GMRFLib_ai_store_tp * ai_store,
			       GMRFLib_preopt_tp * preopt);
int GMRFLib_ai_update_conditional_mean(int thread_id, GMRFLib_problem_tp * pproblem, double *x, double *mean,
				       GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_args,
				       GMRFLib_constr_tp * constr, double *bbb, double *ccc, double **covariances, int idx);
int GMRFLib_ai_update_conditional_mean2(double *cond_mean, GMRFLib_problem_tp * problem, int idx, double evalue, double **covariances);
int GMRFLib_init_GMRF_approximation_store__intern(int thread_id, GMRFLib_problem_tp ** problem, double *x, double *b, double *c, double *mean,
						  double *d, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg,
						  GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg,
						  GMRFLib_constr_tp * constr, GMRFLib_optimize_param_tp * optpar,
						  GMRFLib_blockupdate_param_tp * blockupdate_par, GMRFLib_store_tp * store, double *aa, double *bb,
						  double *cc, int gaussian_data, double c_min, int b_strategy, int nested,
						  GMRFLib_preopt_tp * preopt);
int GMRFLib_free_ai_store(GMRFLib_ai_store_tp * ai_store);

int GMRFLib_ai_INLA(GMRFLib_density_tp *** density,
		    GMRFLib_density_tp *** density_transform, GMRFLib_transform_array_func_tp ** tfunc,
		    GMRFLib_density_tp *** density_hyper,
		    GMRFLib_ai_cpo_tp ** cpo, GMRFLib_ai_po_tp ** po, GMRFLib_ai_dic_tp * dic,
		    GMRFLib_ai_marginal_likelihood_tp * marginal_likelihood,
		    char *compute, double ***hyperparam, int nhyper, GMRFLib_ai_log_extra_tp * log_extra, void *log_extra_arg,
		    double *x, double *b, double *c, double *mean, GMRFLib_bfunc_tp ** bfunc, double *d,
		    GMRFLib_logl_tp * loglFunc, void *loglFunc_arg,
		    GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg,
		    GMRFLib_constr_tp * constr, GMRFLib_ai_param_tp * ai_par, GMRFLib_ai_store_tp * ai_store,
		    int nlin, GMRFLib_lc_tp ** Alin, GMRFLib_density_tp *** dlin, GMRFLib_ai_misc_output_tp * misc_output,
		    GMRFLib_preopt_tp * preopt, GMRFLib_preopt_res_tp * rpreopt);
int GMRFLib_ai_INLA_experimental(GMRFLib_density_tp *** density,
				 GMRFLib_density_tp *** density_transform, GMRFLib_transform_array_func_tp ** tfunc,
				 GMRFLib_density_tp *** density_hyper,
				 GMRFLib_gcpo_tp ** gcpo, GMRFLib_gcpo_param_tp * gcpo_param,
				 GMRFLib_ai_cpo_tp ** cpo, GMRFLib_ai_po_tp ** po, GMRFLib_ai_dic_tp * dic,
				 GMRFLib_ai_marginal_likelihood_tp * marginal_likelihood,
				 double ***hyperparam, int nhyper,
				 GMRFLib_ai_log_extra_tp * log_extra, void *log_extra_arg,
				 double *x, double *b, double *c, double *mean,
				 GMRFLib_bfunc_tp ** bfunc, double *d,
				 GMRFLib_logl_tp * loglFunc, void *loglFunc_arg,
				 GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg,
				 GMRFLib_constr_tp * constr, GMRFLib_ai_param_tp * ai_par, GMRFLib_ai_store_tp * ai_store,
				 int nlin, GMRFLib_lc_tp ** Alin, GMRFLib_density_tp *** dlin, GMRFLib_ai_misc_output_tp * misc_output,
				 GMRFLib_preopt_tp * preopt);
int GMRFLib_ai_store_config(int thread_id, GMRFLib_ai_misc_output_tp * mo,
			    int ntheta, double *theta, double log_posterior, double log_posterior_orig,
			    double *improved_mean, double *skewness, GMRFLib_problem_tp * gmrf_approx,
			    GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg, double *c, int dens_count);
int GMRFLib_ai_store_config_preopt(int thread_id, GMRFLib_ai_misc_output_tp * mo, int ntheta, double *theta, double log_posterior,
				   double log_posterior_orig, GMRFLib_problem_tp * problem, double *mean_corrected,
				   GMRFLib_preopt_tp * preopt, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg,
				   double *cpodens_moments, double *gcpodens_moments);
int GMRFLib_compute_cpodens(int thread_id, GMRFLib_density_tp ** cpo_density, GMRFLib_density_tp * density,
			    int idx, double d, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, GMRFLib_ai_param_tp * ai_par);
int GMRFLib_ai_compute_lincomb(GMRFLib_density_tp *** lindens, double **cross, int nlin, GMRFLib_lc_tp ** Alin, GMRFLib_ai_store_tp * ai_store,
			       double *improved_mean, int lookup_tables);
GMRFLib_ai_store_tp *GMRFLib_duplicate_ai_store(GMRFLib_ai_store_tp * ai_store, int skeleton, int copy_ptr, int copy_pardiso_ptr);
GMRFLib_ai_store_tp *GMRFLib_assign_ai_store(GMRFLib_ai_store_tp * to, GMRFLib_ai_store_tp * from);
double GMRFLib_ai_vb_mEll(int idx, GMRFLib_density_tp * density, double d, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg,
			  double *x_vec, double mean_add, double var_scale);
int GMRFLib_ai_vb_prepare(int thread_id, GMRFLib_vb_coofs_tp * coofs, int idx, GMRFLib_density_tp * density, double d, GMRFLib_logl_tp * loglFunc,
			  void *loglFunc_arg, double *x_vec);
int GMRFLib_ai_vb_prepare_mean(int thread_id,
			       GMRFLib_vb_coofs_tp * coofs, int idx, double d, GMRFLib_logl_tp * loglFunc,
			       void *loglFunc_arg, double *x_vec, double mean, double sd);
int GMRFLib_ai_vb_prepare_variance(int thread_id,
				   GMRFLib_vb_coofs_tp * coofs, int idx, double d, GMRFLib_logl_tp * loglFunc,
				   void *loglFunc_arg, double *x_vec, double mean, double sd);
char *GMRFLib_ai_tag(int *iz, int len);
double GMRFLib_ai_cpopit_integrate(int thread_id, double *cpo, double *pit, int idx, GMRFLib_density_tp * cpo_density, double d,
				   GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, double *x_vec);
double GMRFLib_ai_dic_integrate(int thread_id, int idx, GMRFLib_density_tp * density, double d, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg,
				double *x_vec);
double GMRFLib_ai_po_integrate(int thread_id, double *po, double *po2, double *po3, int idx, GMRFLib_density_tp * po_density, double d,
			       GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, double *x_vec);
double GMRFLib_interpolator_nearest(int ndim, int nobs, double *x, double *xobs, double *yobs, void *arg);
int GMRFLib_ai_add_Qinv_to_ai_store(GMRFLib_ai_store_tp * ai_store);
int GMRFLib_ai_adjust_integration_weights(double *adj_weights, double *weights, double **izs, int n, int nhyper, double dz);
int GMRFLib_ai_correct_cpodens(double *dens, double *x, int *n, GMRFLib_ai_param_tp * ai_par);
int GMRFLib_ai_cpo_free(GMRFLib_ai_cpo_tp * cpo);
int GMRFLib_ai_param_duplicate(GMRFLib_ai_param_tp ** ai_par_new, GMRFLib_ai_param_tp * ai_par);
int GMRFLib_ai_param_free(GMRFLib_ai_param_tp * ai_par);
int GMRFLib_ai_po_free(GMRFLib_ai_po_tp * po);
int GMRFLib_ai_skip_configurations(map_strd * hash_table, int k, int *iz, int *izz, int *len, int *k_max, int len_length, int nhyper);
int GMRFLib_ai_theta2z(double *z, int nhyper, double *theta_mode, double *theta, gsl_vector * sqrt_eigen_values, gsl_matrix * eigen_vectors);
int GMRFLib_ai_validate_cpodens(GMRFLib_density_tp * cpo_density);
int GMRFLib_ai_z2theta(double *theta, int nhyper, double *theta_mode, double *z, gsl_vector * sqrt_eigen_values, gsl_matrix * eigen_vectors);

int GMRFLib_ai_vb_correct_mean(int thread_id, GMRFLib_density_tp *** density,	// need two types
			       int dens_count,
			       GMRFLib_density_tp ** dens_local,
			       double *c,
			       double *d,
			       GMRFLib_ai_param_tp * ai_par,
			       GMRFLib_ai_store_tp * ai_store,
			       GMRFLib_graph_tp * graph,
			       GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg,
			       GMRFLib_preopt_tp * preopt);
int GMRFLib_ai_vb_correct_mean_std(int thread_id, GMRFLib_density_tp *** density,	// need two types
				   int dens_count,
				   GMRFLib_density_tp ** dens_local,
				   double *c,
				   double *d,
				   GMRFLib_ai_param_tp * ai_par,
				   GMRFLib_ai_store_tp * ai_store,
				   GMRFLib_graph_tp * graph,
				   GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg);
int GMRFLib_ai_vb_correct_mean_preopt(int thread_id, GMRFLib_density_tp *** density,
				      int dens_count,
				      double *c,
				      double *d,
				      GMRFLib_ai_param_tp * ai_par,
				      GMRFLib_ai_store_tp * ai_store,
				      GMRFLib_graph_tp * graph,
				      GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg,
				      GMRFLib_preopt_tp * preopt);
int GMRFLib_ai_vb_correct_variance_preopt(int thread_id,
					  GMRFLib_density_tp *** density,
					  int dens_count,
					  double *UNUSED(c),
					  double *d,
					  GMRFLib_ai_param_tp * ai_par,
					  GMRFLib_ai_store_tp * ai_store,
					  GMRFLib_graph_tp * graph,
					  GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg,
					  GMRFLib_preopt_tp * preopt);

double GMRFLib_bfunc_eval(int thread_id, double *con, GMRFLib_bfunc_tp * bfunc);
int GMRFLib_bnew(int thread_id, double **bnew, double *constant, int n, double *b, GMRFLib_bfunc_tp ** bfunc);
int GMRFLib_transform_density(GMRFLib_density_tp ** tdensity, GMRFLib_density_tp * density, GMRFLib_transform_array_func_tp * func);
int GMRFLib_vb_fit_gaussian(int n, double *x, double *ld, double *mean, double *sd);

GMRFLib_gcpo_elm_tp **GMRFLib_gcpo(int thread_id, GMRFLib_ai_store_tp * ai_store_id, double *lpred_mean, double *lpred_mode,
				   double *lpred_variance, GMRFLib_preopt_tp * preopt,
				   GMRFLib_gcpo_groups_tp * groups, double *d, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg,
				   GMRFLib_ai_param_tp * ai_par, GMRFLib_gcpo_param_tp * gcpo_param, double *gcpodens_moments);
GMRFLib_gcpo_groups_tp *GMRFLib_gcpo_build(int thread_id, GMRFLib_ai_store_tp * ai_store, GMRFLib_preopt_tp * preopt,
					   GMRFLib_gcpo_param_tp * gcpo_param);

double inla_compute_saturated_loglik(int, int, GMRFLib_logl_tp *, double *, void *);

__END_DECLS
#endif
