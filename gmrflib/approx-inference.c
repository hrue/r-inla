
/* approx-inference.c
 * 
 * Copyright (C) 2006-2008 Havard Rue
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
 *       H{\aa}vard Rue
 *       Department of Mathematical Sciences
 *       The Norwegian University of Science and Technology
 *       N-7491 Trondheim, Norway
 *       Voice: +47-7359-3533    URL  : http://www.math.ntnu.no/~hrue  
 *       Fax  : +47-7359-3524    Email: havard.rue@math.ntnu.no
 *
 */

/*!
  \file approx-inference.c
  \brief Functions for approximate inference using Integrated Nested Laplace Approximation (INLA). (UNDER CONSTRUCTION)

  These functions perform approximate inference for the posterior
  where 
  \f[ \pi(\mbox{\boldmath$x$}, \mbox{\boldmath$\theta$} \mid \mbox{\boldmath$y$}) \propto
  \pi(\mbox{\boldmath$\theta$})
  \pi(\mbox{\boldmath$x$} \mid \mbox{\boldmath$\theta$}) \prod_{i} \pi(y_i \mid x_i, \mbox{\boldmath$\theta$})
  \f]
  and \f$ \pi(\mbox{\boldmath$x$} \mid \mbox{\boldmath$\theta$}) \f$ is a GMRF.

  The routines compute (approximately)
  - the posterior marginal of \f$\mbox{\boldmath$\theta$}\f$,
  \f[ \pi(\mbox{\boldmath$\theta$} \mid \mbox{\boldmath$y$}) \f] 
  - the posterior marginal of \f$x_i\f$ for fixed value of \f$ \mbox{\boldmath$\theta$} \f$,
  \f[ \pi(x_i \mid \mbox{\boldmath$\theta$}, \mbox{\boldmath$y$}) \f]
  - the posterior marginal of \f$x_i\f$
  \f[ \pi(x_i \mid  \mbox{\boldmath$y$}) \f]
  - the integrated likelihood for the model
  \f[ \pi(\mbox{\boldmath$y$}) \f]
  See also \ref INLA for some more details, and \ref ex_ai for an example.\n
*/

#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

/* Pre-hg-Id: $Id: approx-inference.c,v 1.713 2010/04/10 20:07:01 hrue Exp $ */

#include <assert.h>
#include <time.h>
#include <stddef.h>
#include <float.h>
#include <math.h>
#include <strings.h>
#include <string.h>
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"
#include "GMRFLib/hashP.h"

int domin_(void);
int domin_get_results_(double *theta, double *log_dens, int *ierr);
int domin_seteps_(double *epsx, double *epsf, double *epsg);

static int pool_nhyper = -1;



/*!
  \brief Create a \c GMRFLib_ai_param_tp -object holding the default values.
  
  \param[out] ai_par A pointer to a \c GMRFLib_ai_param_tp pointer. At output the \c GMRFLib_ai_param_tp -object contains the
  default values.

  \par Description of elements in \c GMRFLib_ai_param_tp -object:

  \em GMRFLib_ai_param_tp::strategy : The strategy used to compute approximations to the posterior marginals \f$ \pi(x_i \mid
  \mbox{\boldmath$\theta$}, \mbox{\boldmath$y$}) \f$. One of GMRFLib_ai_strategy_tp. \n Default value:
  #GMRFLib_AI_STRATEGY_MEANCORRECTED_GAUSSIAN\n

  \em GMRFLib_ai_param_tp::fast : Use the conditinal mean as an approximation to the conditional mode, and approximate the log
  of the denominator in the Laplace approximation by a linear term.\n <b>Default value: #GMRFLib_TRUE </b>\n

  \em GMRFLib_ai_param_tp::linear_correction One of GMRFLib_ai_linear_correction_tp. Various ways to approximate the log of the
  denominator in the Laplace approximation by a linear term. If GMRFLib_ai_param_tp::fast, and the meancorrected or
  skewmeancorrected strategy is selected, then GMRFLib_ai_param_tp::linear_correction is #GMRFLib_AI_LINEAR_CORRECTION_FAST,
  which computes the linear term exact.

  \em GMRFLib_ai_param_tp::int_strategy : Use a grid strategy or a central composite design to integrate out the
  hyperparameters.\n <b>Default value: #GMRFLib_AI_INT_STRATEGY_GRID </b>\n
  
  \em GMRFLib_ai_param_tp::f0 : Only used if \a int_strategy = \a #GMRFLib_AI_INT_STRATEGY_CCD. \n Parameter used in the CCD
  integration.\n <b>Default value: 1.1 </b>\n

  \em GMRFLib_ai_param_tp::dz : Step length for the integration procedure.\n <b>Default value: 1 </b>\n

  \em GMRFLib_ai_param_tp::adjust_weights : Adjust the weights in the integration procedure so that a Gaussian density is
  integrated with no error.\n <b>Default value: #GMRFLib_TRUE </b>\n
  
  \em GMRFLib_ai_param_tp::diff_logdens : Only used if \a int_strategy = #GMRFLib_AI_INT_STRATEGY_GRID. \n Threshhold for
  accepting a configuration.\n <b>Default value: 2.5 </b>\n
  
  \em GMRFLib_ai_param_tp::skip_configurations : Only used if \a int_strategy = #GMRFLib_AI_INT_STRATEGY_GRID. \n Skip fill-in
  configuration larger than a non-accepted one.\n <b>Default value: #GMRFLib_TRUE </b>\n
  
  \em GMRFLib_ai_param_tp::gradient_forward_finite_difference : Use forward finite difference to compute the gradient.\n
  <b>Default value: #GMRFLib_TRUE </b>\n
  
  \em GMRFLib_ai_param_tp::gradient_finite_difference_step_len : Step length to compute the gradient.\n <b>Default value: 1.0e-4
  </b>\n
  
  \em GMRFLib_ai_param_tp::hessian_forward_finite_difference : Use forward finite difference to compute the hessian at the
  mode.\n <b>Default value: #GMRFLib_TRUE </b>\n
  
  \em GMRFLib_ai_param_tp::hessian_finite_difference_step_len : Step length to compute the hessian at the mode.\n <b>Default
  value: 1.0e-4 </b>\n

  \em GMRFLib_ai_param_tp::hessian_force_diagonal : Force the hessian to be diagonal.\n <b>Default value: #GMRFLib_FALSE </b>\n
  
*/
int GMRFLib_default_ai_param(GMRFLib_ai_param_tp ** ai_par)
{
	/*
	 * if the MEANCORRECTED_GAUSSIAN strategy is selected, then some kind of linear_correction need to be selected. if
	 * linear_correction is OFF then its selected on-site depending on the value of 'fast'. 
	 */

	*ai_par = Calloc(1, GMRFLib_ai_param_tp);

	(*ai_par)->strategy = GMRFLib_AI_STRATEGY_GAUSSIAN;
	(*ai_par)->strategy = GMRFLib_AI_STRATEGY_FIT_SCGAUSSIAN;
	(*ai_par)->strategy = GMRFLib_AI_STRATEGY_MEANCORRECTED_GAUSSIAN;
	(*ai_par)->strategy = GMRFLib_AI_STRATEGY_MEANSKEWCORRECTED_GAUSSIAN;

	(*ai_par)->fast = GMRFLib_FALSE;		       /* compute conditional mode */
	(*ai_par)->fast = GMRFLib_TRUE;			       /* use mode = conditional mean */

	(*ai_par)->gaussian_data = GMRFLib_TRUE;
	(*ai_par)->gaussian_data = GMRFLib_FALSE;

	(*ai_par)->linear_correction = GMRFLib_AI_LINEAR_CORRECTION_CENTRAL_DIFFERENCE;
	(*ai_par)->linear_correction = GMRFLib_AI_LINEAR_CORRECTION_OFF;
	(*ai_par)->linear_correction = GMRFLib_AI_LINEAR_CORRECTION_FAST;	/* to match the ->fast mode above */

	(*ai_par)->si_directory = NULL;

	/*
	 * none of these are used, but they are the defaults if the user wants improved approximations 
	 */
	(*ai_par)->n_points = 9;			       /* how many points to evaluate */
	(*ai_par)->step_len = GMRFLib_eps(0.25);	       /* If the derivaties has to be computed numerically */
	(*ai_par)->cutoff = 0.0;			       /* the cutoff for the gradient in the (Gaussian) conditional mean */

	/*
	 * defaults for the integration itself 
	 */
	(*ai_par)->fp_log = NULL;
	(*ai_par)->fp_log = stdout;
	(*ai_par)->fp_hyperparam = NULL;
	(*ai_par)->int_strategy = GMRFLib_AI_INT_STRATEGY_EMPIRICAL_BAYES;
	(*ai_par)->int_strategy = GMRFLib_AI_INT_STRATEGY_CCD;
	(*ai_par)->int_strategy = GMRFLib_AI_INT_STRATEGY_GRID;
	(*ai_par)->f0 = 1.1;
	(*ai_par)->dz = 1.0;
	(*ai_par)->adjust_weights = GMRFLib_FALSE;
	(*ai_par)->adjust_weights = GMRFLib_TRUE;
	(*ai_par)->diff_log_dens = 2.5;
	(*ai_par)->skip_configurations = GMRFLib_FALSE;
	(*ai_par)->skip_configurations = GMRFLib_TRUE;

	/*
	 * use forward differences for the gradient and central differences for the (final) hessian 
	 */
	(*ai_par)->gradient_forward_finite_difference = GMRFLib_FALSE;	/* use central difference */
	(*ai_par)->gradient_forward_finite_difference = GMRFLib_TRUE;	/* use forward difference */
	(*ai_par)->gradient_finite_difference_step_len = 0.01;

	(*ai_par)->hessian_forward_finite_difference = GMRFLib_TRUE;	/* use forward difference */
	(*ai_par)->hessian_forward_finite_difference = GMRFLib_FALSE;	/* use central difference */
	(*ai_par)->hessian_finite_difference_step_len = sqrt((*ai_par)->gradient_finite_difference_step_len);

	(*ai_par)->hessian_force_diagonal = GMRFLib_TRUE;
	(*ai_par)->hessian_force_diagonal = GMRFLib_FALSE;

	/*
	 * these are only valid if fp_log is !NULL. 
	 */
	(*ai_par)->compute_nparam_eff = GMRFLib_TRUE;
	(*ai_par)->do_MC_error_check = GMRFLib_FALSE;

	/*
	 * choice of interpolator 
	 */
	(*ai_par)->interpolator = GMRFLib_AI_INTERPOLATOR_GAUSSIAN;
	(*ai_par)->interpolator = GMRFLib_AI_INTERPOLATOR_NEAREST;
	(*ai_par)->interpolator = GMRFLib_AI_INTERPOLATOR_LINEAR;
	(*ai_par)->interpolator = GMRFLib_AI_INTERPOLATOR_QUADRATIC;
	(*ai_par)->interpolator = GMRFLib_AI_INTERPOLATOR_CCD;
	(*ai_par)->interpolator = GMRFLib_AI_INTERPOLATOR_WEIGHTED_DISTANCE;
	(*ai_par)->interpolator = GMRFLib_AI_INTERPOLATOR_AUTO;	/* automatic choice */


	/*
	 * Type of optimiser to use: GMRFLib_AI_OPTIMISER_DOMIN, GMRFLib_AI_OPTIMISER_GSL
	 */
	(*ai_par)->optimiser = GMRFLib_AI_OPTIMISER_DEFAULT;
	(*ai_par)->restart = 0;
	(*ai_par)->domin_epsx = 0.005;
	(*ai_par)->domin_epsf = 0.005;
	(*ai_par)->domin_epsg = 0.005;
	(*ai_par)->gsl_tol = 0.1;
	(*ai_par)->gsl_epsg = 0.005;
	(*ai_par)->gsl_step_size = 1.0;
	(*ai_par)->mode_known = 0;

	/*
	 * parameters for the Gaussian approximations 
	 */
	(*ai_par)->optpar_abserr_func = 0.01;
	(*ai_par)->optpar_abserr_step = 0.01;
	(*ai_par)->optpar_fp = NULL;

	(*ai_par)->cpo_req_diff_logdens = 3.0;

	(*ai_par)->adaptive_hessian_mode = GMRFLib_FALSE;
	(*ai_par)->adaptive_hessian_mode = GMRFLib_TRUE;
	(*ai_par)->adaptive_hessian_max_trials = 1000;
	(*ai_par)->adaptive_hessian_scale = 1.01;

	(*ai_par)->huge = GMRFLib_FALSE;
	(*ai_par)->cpo_manual = GMRFLib_FALSE;

	return GMRFLib_SUCCESS;
}

/*!
  \brief Print the values of a \c GMRFLib_ai_param_tp -object

  \param[out] fp The *FILE  on which to print the output
  
  \param[in] ai_par The \c GMRFLib_ai_param_tp -object to be printed
*/
int GMRFLib_print_ai_param(FILE * fp, GMRFLib_ai_param_tp * ai_par)
{
	int show_expert_options = 1, i;

	if (!ai_par) {
		return GMRFLib_SUCCESS;
	}
	fp = (fp ? fp : stdout);

	fprintf(fp, "Contents of ai_param %p\n", (void *) ai_par);

	fprintf(fp, "\tOptimiser: %s\n", GMRFLib_AI_OPTIMISER_NAME(ai_par->optimiser));
	fprintf(fp, "\t\tOption for %s: epsx = %.6g\n", GMRFLib_AI_OPTIMISER_NAME(GMRFLib_AI_OPTIMISER_DOMIN), ai_par->domin_epsx);
	fprintf(fp, "\t\tOption for %s: epsf = %.6g\n", GMRFLib_AI_OPTIMISER_NAME(GMRFLib_AI_OPTIMISER_DOMIN), ai_par->domin_epsf);
	fprintf(fp, "\t\tOption for %s: epsg = %.6g\n", GMRFLib_AI_OPTIMISER_NAME(GMRFLib_AI_OPTIMISER_DOMIN), ai_par->domin_epsg);
	fprintf(fp, "\t\tOption for %s: tol  = %.6g\n", GMRFLib_AI_OPTIMISER_NAME(GMRFLib_AI_OPTIMISER_GSL), ai_par->gsl_tol);
	fprintf(fp, "\t\tOption for %s: epsg = %.6g\n", GMRFLib_AI_OPTIMISER_NAME(GMRFLib_AI_OPTIMISER_GSL), ai_par->gsl_epsg);
	fprintf(fp, "\t\tOption for %s: step_size = %.6g\n", GMRFLib_AI_OPTIMISER_NAME(GMRFLib_AI_OPTIMISER_GSL), ai_par->gsl_step_size);
	fprintf(fp, "\t\tRestart: %1d\n", ai_par->restart);
	fprintf(fp, "\t\tMode known: %s\n", (ai_par->mode_known ? "Yes" : "No"));

	fprintf(fp, "\tGaussian approximation:\n");
	fprintf(fp, "\t\tabserr_func = %.6g\n", ai_par->optpar_abserr_func);
	fprintf(fp, "\t\tabserr_step = %.6g\n", ai_par->optpar_abserr_step);
	fprintf(fp, "\t\toptpar_fp = %lx\n", (long unsigned int) (ai_par->optpar_fp));

	fprintf(fp, "\tGaussian data: %s\n", (ai_par->gaussian_data ? "Yes" : "No"));

	fprintf(fp, "\tStrategy: \t");
	if (ai_par->strategy == GMRFLib_AI_STRATEGY_GAUSSIAN) {
		fprintf(fp, "Use the Gaussian approximation\n");
	}
	if (ai_par->strategy == GMRFLib_AI_STRATEGY_FIT_SCGAUSSIAN) {
		fprintf(fp, "Fit a spline-corrected Gaussian\n");
	}
	if (ai_par->strategy == GMRFLib_AI_STRATEGY_MEANCORRECTED_GAUSSIAN) {
		fprintf(fp, "Use a mean-corrected Gaussian\n");
	}
	if (ai_par->strategy == GMRFLib_AI_STRATEGY_MEANSKEWCORRECTED_GAUSSIAN) {
		fprintf(fp, "Use a mean-skew corrected Gaussian by fitting a Skew-Normal\n");
	}

	fprintf(fp, "\tFast mode: \t%s\n", (ai_par->fast == GMRFLib_FALSE ? "Off" : "On"));

	fprintf(fp, "\tUse linear approximation to log(|Q +c|)? %s\n", (ai_par->linear_correction == GMRFLib_AI_LINEAR_CORRECTION_OFF ? "No" : "Yes"));
	if (ai_par->linear_correction != GMRFLib_AI_LINEAR_CORRECTION_OFF) {
		fprintf(fp, "\t\tMethod:\t ");
		if (ai_par->linear_correction == GMRFLib_AI_LINEAR_CORRECTION_CENTRAL_DIFFERENCE) {
			fprintf(fp, "Estimate the derivative using central difference\n");
		}
		if (ai_par->linear_correction == GMRFLib_AI_LINEAR_CORRECTION_FAST) {
			fprintf(fp, "Compute the derivative exact\n");
		}
	}
	fprintf(fp, "\tSI directory: \t%s\n", (ai_par->si_directory == NULL ? "<NONE>" : ai_par->si_directory));

	fprintf(fp, "\tParameters for improved approximations\n");
	fprintf(fp, "\t\tNumber of points evaluate:\t %d\n", ai_par->n_points);
	fprintf(fp, "\t\tStep length to compute derivatives numerically:\t %f\n", ai_par->step_len);
	fprintf(fp, "\t\tCutoff value to construct local neigborhood:\t %f\n", ai_par->cutoff);

	fprintf(fp, "\tLog calculations:\t %s\n", (ai_par->fp_log ? "On" : "Off"));
	fprintf(fp, "\tLog calculated marginal for the hyperparameters:\t %s\n", (ai_par->fp_hyperparam ? "On" : "Off"));

	fprintf(fp, "\tIntegration strategy:\t ");
	if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_CCD) {
		fprintf(fp, "Use points from Central Composite Design (CCD)\n");
	}
	if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_GRID) {
		fprintf(fp, "Use adaptive grid-approach (GRID)\n");
	}
	if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_EMPIRICAL_BAYES) {
		fprintf(fp, "Use only the modal configuration (EMPIRICAL_BAYES)\n");
	}
	fprintf(fp, "\t\tf0 (CCD only):\t %f\n", ai_par->f0);
	fprintf(fp, "\t\tdz (GRID only):\t %f\n", ai_par->dz);
	fprintf(fp, "\t\tAdjust weights (GRID only):\t %s\n", (ai_par->adjust_weights == GMRFLib_FALSE ? "Off" : "On"));
	fprintf(fp, "\t\tDifference in log-density limit (GRID only):\t %f\n", ai_par->diff_log_dens);
	fprintf(fp, "\t\tSkip configurations with (presumed) small density (GRID only):\t %s\n", (ai_par->skip_configurations == GMRFLib_FALSE ? "Off" : "On"));

	fprintf(fp, "\tGradient is computed using %s with step-length %f\n",
		(ai_par->gradient_forward_finite_difference == GMRFLib_TRUE ? "Forward difference" : "Central difference"),
		ai_par->gradient_finite_difference_step_len);
	fprintf(fp, "\tHessian is computed using %s with step-length %f\n",
		(ai_par->hessian_forward_finite_difference == GMRFLib_TRUE ? "Forward difference" : "Central difference"),
		ai_par->hessian_finite_difference_step_len);
	fprintf(fp, "\tHessian matrix is forced to be a diagonal matrix? [%s]\n", (ai_par->hessian_force_diagonal ? "Yes" : "No"));

	if (ai_par->fp_log) {
		fprintf(fp, "\tCompute effective number of parameters? [%s]\n", (ai_par->compute_nparam_eff ? "Yes" : "No"));
		fprintf(fp, "\tPerform a Monte Carlo error-test? [%s]\n", (ai_par->do_MC_error_check ? "Yes" : "No"));
		if (ai_par->do_MC_error_check) {
			if (ai_par->do_MC_error_check > 0) {
				fprintf(fp, "\t\tUsing [default] number of samples in the test\n");
			} else {
				fprintf(fp, "\t\tUse [%1d] number of samples in the test\n", -ai_par->do_MC_error_check);
			}
		}
	}

	fprintf(fp, "\tInterpolator [%s]\n", INTERPOLATOR_NAME(ai_par->interpolator));
	fprintf(fp, "\tCPO required diff in log-density [%g]\n", ai_par->cpo_req_diff_logdens);

	fprintf(fp, "\tAdaptive estimation of the Hessian:\n");
	fprintf(fp, "\t\tStatus     [%s]\n", (ai_par->adaptive_hessian_mode ? "On" : "Off"));
	fprintf(fp, "\t\tMax trials [%d]\n", ai_par->adaptive_hessian_max_trials);
	fprintf(fp, "\t\tScale      [%g]\n", ai_par->adaptive_hessian_scale);

	fprintf(fp, "\tHuge model [%s]\n", (ai_par->huge ? "Yes" : "No"));

	if (show_expert_options) {
		/*
		 * expert options goes here 
		 */
		fprintf(fp, "\tCPO manual calculation[%s]\n", (ai_par->cpo_manual ? "Yes" : "No"));
	}

	if (ai_par->si_idx) {
		fprintf(fp, "\tDump information\n");
		fprintf(fp, "\t\tnd = %1d\n", ai_par->si_idx->nd);
		for (i = 0; i < ai_par->si_idx->nd; i++) {
			fprintf(fp, "\t\t block %1d: start=%1d len=%1d tag=[%s]\n", i, ai_par->si_idx->start[i], ai_par->si_idx->len[i], ai_par->si_idx->tag[i]);
		}
	}
	fprintf(fp, "\n");

	return GMRFLib_SUCCESS;
}

/*!

  \brief Returns the unnormalized (approximate) value of \f[ \pi(\mbox{\boldmath$\theta$} \mid \mbox{\boldmath$y$}) \f]
  
  \param[out] logdens The log of the unnormalized posterior marginal for \f$\mbox{\boldmath$\theta$}\f$
  
  \param[in] x A length \em n array, where \n is the number of nodes in the graph. If \a x = \c NULL then all elements are taken
  to be zero.  If \a fixed_value \f$ \neq \f$ \c NULL, the elements of <em>\b x</em> corresponding to \a fixed_value=1 are the
  fixed values in a conditional simulation. The remaining elements of <em>\b x</em> are not used, and can take arbitrary
  values. If \a fixed_value = \c NULL, all values can be arbitrary.

  \param[in] b If <tt>!NULL</tt>, a length \em n array holding the elements of the vector <em>\b b</em> in <b>(GMRF-34)</b> in
  \ref INLA.  If \c NULL, <em>\b b</em> is assumed to be 0.

  \param[in] c If <tt>!NULL</tt>, a length \em n array holding the elements of the vector \f$ \mbox{\boldmath$\mu$} \f$ in
  <b>(GMRF-34)</b> in \ref INLA.  If \c NULL, <em>\b c</em> is assumed to be 0.

  \param[in] mean If <tt>!NULL</tt>, a length \em n array holding the elements of the vector <em>\b </em> in <b>(GMRF-34)</b> in
  \ref INLA. If \c NULL, \f$ \mbox{\boldmath$\mu$} \f$ is assumed to be 0.

  \param[in] d If <tt>!NULL</tt>, a length \em n array holding the elements of the vector <em>\b d</em> in <b>(GMRF-34)</b> in
  \ref INLA.

  \param[in] loglFunc A function of type \c GMRFLib_logl_tp(), returning the value of the function \f$ f_i(x_i,y_i) \f$ in
  (GMRF-34) in \ref INLA, in many applications equal to the log-likelihood of the problem. If the function \f$ f_i(x_i,y_i) \f$
  depends on unknown hyper-parameters, these are passed throug in \a loglFunc_arg

  \param[in] loglFunc_arg \em void pointer holding the address of a variable or data structure defining additional arguments to
  the function \a loglFunc.

  \param[in] fixed_value The elements of this array, of length \c graph->n, should take the value 0 or 1. The conditional mean
  and covariance matrix of \f$ \{x_i:\mbox{\small\tt fixed\_value}[i]=0\} \f$ given \f$ \{x_i:\mbox{\small\tt
  fixed\_value}[i]=1\} \f$ are computed, and elements of the GMRF <em>\b x</em> corresponding to <tt>fixed_value = 1</tt>, are
  kept fixed at the corresponding values specified in the argument \a x. It is allowed that all values of <tt>fixed_value</tt>
  equals 1.

  \param[in] graph The graph on which the GMRF <em>\b x</em> is defined.

  \param[in] Qfunc A function of type \c GMRFLib_Qfunc_tp(), computing the values of the precision matrix <em>\b Q</em>. If the
  function depends on unknown hyper-parameters, these are passed through \a Qfunc_arg

  \param[in] Qfunc_arg A \em void pointer holding the address of a variable or data structure defining additional arguments to
  the function \a Qfunc.

  \param[in] constr A pointer to a \c GMRFLib_constr_tp -object holding the information on the type of linear constraint, in the
  case of constrained problem.

  \param[in] ai_par Some specifications for the approximation, as a \c GMRFLib_ai_param_tp -object.

  \param[in] ai_store A pointer to a \c GMRFLib_ai_store_tp object where temporary calculations will be stored and retrieved to
  increase the speed.

*/
int GMRFLib_ai_marginal_hyperparam(double *logdens,
				   double *x, double *b, double *c, double *mean, double *d,
				   GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, char *fixed_value,
				   GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg,
				   GMRFLib_constr_tp * constr, GMRFLib_ai_param_tp * ai_par, GMRFLib_ai_store_tp * ai_store)
{
	/*
	 * return the unnormalised log marginal density for the hyperparamers in `logdens' 
	 */

	double ldens;
	int n, free_ai_par = 0;

	GMRFLib_problem_tp *problem = NULL;
	GMRFLib_optimize_param_tp *optpar = NULL;
	GMRFLib_blockupdate_param_tp *blockpar = NULL;

	GMRFLib_ENTER_ROUTINE;

	GMRFLib_ASSERT(ai_store, GMRFLib_EPARAMETER);	       /* this is required */

	GMRFLib_default_optimize_param(&optpar);
	GMRFLib_default_blockupdate_param(&blockpar);
	if (!ai_par) {
		GMRFLib_default_ai_param(&ai_par);
		free_ai_par = 1;
	}
	optpar->step_len = ai_par->step_len;
	optpar->abserr_func = ai_par->optpar_abserr_func;
	optpar->abserr_step = ai_par->optpar_abserr_func;
	optpar->abserr_step = ai_par->optpar_abserr_func;
	optpar->fp = ai_par->optpar_fp;
	blockpar->step_len = ai_par->step_len;
	blockpar->modeoption = GMRFLib_MODEOPTION_MODE;
	blockpar->fp = ai_par->optpar_fp;

	if (!ai_store->store) {
		ai_store->store = Calloc(1, GMRFLib_store_tp);
	}

	n = graph->n;
	Free(ai_store->correction_term);
	Free(ai_store->correction_idx);
	Free(ai_store->derivative3);
	Free(ai_store->bb);
	Free(ai_store->cc);
	ai_store->bb = Calloc(n, double);
	ai_store->cc = Calloc(n, double);


	/*
	 * tabulate the Qfunction, which speedup... 
	 */
	GMRFLib_tabulate_Qfunc_tp *tabQfunc = NULL;
	int use_tabulated_Qfunc = 1;

	GMRFLib_tabulate_Qfunc(&tabQfunc, graph, Qfunc, Qfunc_arg, NULL, NULL, NULL);

	/*
	 * first compute the GMRF-approximation 
	 */
	if (use_tabulated_Qfunc) {
		/*
		 * Use the tabulated Qfunc here 
		 */
		if (ai_store->mode) {
			GMRFLib_EWRAP1(GMRFLib_init_GMRF_approximation_store__intern
				       (&problem, ai_store->mode, b, c, mean, d, loglFunc, loglFunc_arg, fixed_value, graph,
					tabQfunc->Qfunc, tabQfunc->Qfunc_arg,
					constr, optpar, blockpar, ai_store->store, ai_store->bb, ai_store->cc, ai_par->gaussian_data));
		} else {
			GMRFLib_EWRAP1(GMRFLib_init_GMRF_approximation_store__intern
				       (&problem, x, b, c, mean, d, loglFunc, loglFunc_arg, fixed_value, graph,
					tabQfunc->Qfunc, tabQfunc->Qfunc_arg,
					constr, optpar, blockpar, ai_store->store, ai_store->bb, ai_store->cc, ai_par->gaussian_data));
		}
	} else {

		/*
		 * old version! 
		 */
		if (ai_store->mode) {
			GMRFLib_EWRAP1(GMRFLib_init_GMRF_approximation_store__intern
				       (&problem, ai_store->mode, b, c, mean, d, loglFunc, loglFunc_arg, fixed_value, graph, Qfunc,
					Qfunc_arg, constr, optpar, blockpar, ai_store->store, ai_store->bb, ai_store->cc, ai_par->gaussian_data));
		} else {
			GMRFLib_EWRAP1(GMRFLib_init_GMRF_approximation_store__intern
				       (&problem, x, b, c, mean, d, loglFunc, loglFunc_arg, fixed_value, graph, Qfunc, Qfunc_arg,
					constr, optpar, blockpar, ai_store->store, ai_store->bb, ai_store->cc, ai_par->gaussian_data));
		}
	}

	/*
	 * if store, then store the mode to use as the initial point at later calls 
	 */
	Free(ai_store->mode);
	ai_store->mode = Calloc(n, double);
	memcpy(ai_store->mode, problem->mean_constr, n * sizeof(double));

	/*
	 * then evaluate it in the mode (ie the mean) to get the log-density 
	 */
	memcpy(problem->sample, problem->mean_constr, n * sizeof(double));
	GMRFLib_EWRAP1(GMRFLib_evaluate(problem));

	if (use_tabulated_Qfunc) {
		/*
		 * Use the tabulated Qfunction 
		 */
		GMRFLib_ai_log_posterior(&ldens, problem->mean_constr, b, c, mean, d, loglFunc, loglFunc_arg, fixed_value, graph,
					 tabQfunc->Qfunc, tabQfunc->Qfunc_arg, constr);
	} else {
		/*
		 * old code 
		 */
		GMRFLib_ai_log_posterior(&ldens, problem->mean_constr, b, c, mean, d, loglFunc, loglFunc_arg, fixed_value, graph, Qfunc, Qfunc_arg, constr);
	}

	// printf("ai_marginal_hyper thread_id %d ldens %.12f sub_logdens %.12f\n", GMRFLib_thread_id, ldens, problem->sub_logdens);

	/*
	 * End of tabulate_Qfunc 
	 */
	GMRFLib_free_tabulate_Qfunc(tabQfunc);
	tabQfunc = NULL;

	/*
	 **
	 */

	*logdens = ldens - problem->sub_logdens;

	/*
	 * store the GMRF-approximation in store 
	 */
	GMRFLib_free_problem(ai_store->problem);
	ai_store->problem = problem;

	/*
	 * cleanup 
	 */
	Free(optpar);
	Free(blockpar);
	if (free_ai_par)
		Free(ai_par);

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}
int GMRFLib_ai_log_posterior(double *logdens,
			     double *x, double *b, double *c, double *mean, double *d,
			     GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, char *fixed_value,
			     GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg, GMRFLib_constr_tp * constr)
{
	/*
	 * compute the log posterior of configuration 'x' up to an additive constant and return the value in 'logdens'. 
	 */

	int i, n, id, run_with_omp;
	double *xx = NULL, val, logll, sqr_term, result;

	GMRFLib_ENTER_ROUTINE;

	run_with_omp = (omp_get_max_threads() > 1 ? 1 : 0);
	n = graph->n;
	xx = Calloc(n, double);				       /* xx = x - mean */

	if (mean) {
		for (i = 0; i < n; i++) {
			xx[i] = x[i] - mean[i];
		}
	} else {
		memcpy(xx, x, n * sizeof(double));
	}

	id = GMRFLib_thread_id;

	if (run_with_omp) {
		double tmp0 = 0.0, tmp1 = 0.0, tmp2 = 0.0, tmp3 = 0.0, tmp4 = 0.0;

		{
			{
				GMRFLib_xQx(&result, xx, graph, Qfunc, Qfunc_arg);
				tmp4 = -0.5 * result;
			}

			{
				/*
				 * add the diagonal if 'c' 
				 */
				if (c) {
					int ii;

					for (ii = 0; ii < n; ii++) {
						tmp0 += c[ii] * SQR(xx[ii]);
					}
					tmp0 = -0.5 * tmp0;
				}
			}
			{
				/*
				 * add the linear term 
				 */
				if (b) {
					int ii;

					for (ii = 0; ii < n; ii++) {
						tmp1 += b[ii] * x[ii];
					}
				}
			}
			if (d) {
				/*
				 * do not include fixed points 
				 */
				int ii;

				if (1) {
					/*
					 * new code; better for omp 
					 */
					int *idxs = NULL, nidx, iii;

					idxs = Calloc(n, int);
					nidx = 0;
					if (fixed_value) {
						for (ii = 0; ii < n; ii++) {
							if (d[ii] && !fixed_value[ii]) {
								idxs[nidx++] = ii;
							}
						}
#pragma omp parallel for private(ii, iii, logll) reduction(+: tmp2) schedule(static)
						for (iii = 0; iii < nidx; iii++) {
							GMRFLib_thread_id = id;
							ii = idxs[iii];
							loglFunc(&logll, &x[ii], 1, ii, x, loglFunc_arg);
							tmp2 += d[ii] * logll;
						}
						GMRFLib_thread_id = id;
					} else {
						for (ii = 0; ii < n; ii++) {
							if (d[ii]) {
								idxs[nidx++] = ii;
							}
						}
#pragma omp parallel for private(ii, iii, logll) reduction(+: tmp2) schedule(static)
						for (iii = 0; iii < nidx; iii++) {
							GMRFLib_thread_id = id;
							ii = idxs[iii];
							loglFunc(&logll, &x[ii], 1, ii, x, loglFunc_arg);
							tmp2 += d[ii] * logll;
						}
						GMRFLib_thread_id = id;
					}
					Free(idxs);
				} else {
					/*
					 * old version 
					 */
					GMRFLib_thread_id = id;
					if (fixed_value) {
						for (ii = 0; ii < n; ii++) {
							if (d[ii] && !fixed_value[ii]) {
								loglFunc(&logll, &x[ii], 1, ii, x, loglFunc_arg);
								tmp2 += d[ii] * logll;
							}
						}
					} else {
						for (ii = 0; ii < n; ii++) {
							if (d[ii]) {
								loglFunc(&logll, &x[ii], 1, ii, x, loglFunc_arg);
								tmp2 += d[ii] * logll;
							}
						}
					}
				}
			}
			{
				/*
				 * adjust if stochastic constraint 
				 */

				if (STOCHASTIC_CONSTR(constr)) {
					GMRFLib_eval_constr(NULL, &sqr_term, x, constr, graph);
					tmp3 = -0.5 * sqr_term;
				}
			}
		}
		val = tmp0 + tmp1 + tmp2 + tmp3 + tmp4;
	} else {
		GMRFLib_EWRAP1(GMRFLib_xQx(&result, xx, graph, Qfunc, Qfunc_arg));
		val = -0.5 * result;

		double tmp;

		/*
		 * add the diagonal if 'c' 
		 */
		if (c) {
			tmp = 0.0;
			for (i = 0; i < n; i++) {
				tmp += c[i] * SQR(xx[i]);
			}
			val += -0.5 * tmp;
		}

		/*
		 * add the linear term 
		 */
		if (b) {
			tmp = 0.0;
			for (i = 0; i < n; i++) {
				tmp += b[i] * x[i];
			}
			val += tmp;
		}
		if (d) {
			/*
			 * do not include fixed points 
			 */
			tmp = 0.0;
			if (fixed_value) {
				for (i = 0; i < n; i++) {
					if (d[i] && !fixed_value[i]) {
						loglFunc(&logll, &x[i], 1, i, x, loglFunc_arg);
						tmp += d[i] * logll;
					}
				}
			} else {
				for (i = 0; i < n; i++) {
					if (d[i]) {
						loglFunc(&logll, &x[i], 1, i, x, loglFunc_arg);
						tmp += d[i] * logll;
					}
				}
			}
			val += tmp;
		}
		/*
		 * adjust if stochastic constraint 
		 */
		if (STOCHASTIC_CONSTR(constr)) {
			GMRFLib_EWRAP1(GMRFLib_eval_constr(NULL, &sqr_term, x, constr, graph));
			val += -0.5 * sqr_term;
		}
	}

	Free(xx);

	*logdens = val;

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}
int GMRFLib_ai_log_posterior_restricted(double *logdens, double *x, double *x_mode, double *x_gradient, double delta,
					double *b, double *c, double *mean, double *d,
					GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, char *fixed_value,
					GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg,
					GMRFLib_constr_tp * constr, GMRFLib_graph_tp * subgraph, GMRFLib_ai_store_tp * ai_store)
{
	/*
	 * this is the same function as GMRFLib_ai_log_posterior, BUT we only include those terms where at least one component
	 * is NOT FIXED.
	 * 
	 * the added last argument, subgraph, is added as an argument since its fixed for many repeated calls to this function.
	 * 
	 * if logdens==NULL, then the routine is initialised and the linear and quadratic term are computed. these terms are
	 * used for later successive calls.
	 * 
	 * if subgraph == NULL, then we use the GMRFLib_ai_log_posterior()-function
	 *
	 */

	/*
	 * ASSUME NOW THAT SUBGRAPH IS A COPY OF GRAPH!!! ie FIXED = NULL 
	 */

	int i, j, ii, jj, ns;
	double xx, *f = NULL, *g = NULL, val, tmp, logll = 0.0, q_value;

	static double quadratic_term = 0.0, linear_term = 0.0; /* compute those if logdens == NULL */
#pragma omp threadprivate(quadratic_term, linear_term)

	/*
	 * if subgraph is not available, use the default routine 
	 */
	if (!subgraph) {
		GMRFLib_EWRAP0(GMRFLib_ai_log_posterior(logdens, x, b, c, mean, d, loglFunc, loglFunc_arg, fixed_value, graph, Qfunc, Qfunc_arg, constr));
		return GMRFLib_SUCCESS;
	}
	GMRFLib_ENTER_ROUTINE;

	ns = subgraph->n;

	if (!logdens) {
		/*
		 * compute the quadratic and linear term 
		 */
		f = Calloc(ns, double);
		g = Calloc(ns, double);


		for (ii = 0; ii < ns; ii++) {
			// i = subgraph->mothergraph_idx[ii];
			i = ii;
			xx = (mean ? x_mode[i] - mean[i] : x_mode[i]);

			if (c) {
				q_value = (Qfunc(i, i, Qfunc_arg) + c[i]);
				f[ii] += x_gradient[i] * q_value;
				g[ii] += xx * q_value;
			} else {
				q_value = Qfunc(i, i, Qfunc_arg);
				f[ii] += x_gradient[i] * q_value;
				g[ii] += xx * q_value;
			}
			if (mean) {
				for (jj = 0; jj < graph->nnbs[i]; jj++) {
					j = graph->nbs[i][jj];
					q_value = Qfunc(i, j, Qfunc_arg);
					f[ii] += x_gradient[j] * q_value;
					g[ii] += (x_mode[j] - mean[j]) * q_value;
				}
			} else {
				for (jj = 0; jj < graph->nnbs[i]; jj++) {
					j = graph->nbs[i][jj];
					q_value = Qfunc(i, j, Qfunc_arg);
					f[ii] += x_gradient[j] * q_value;
					g[ii] += x_mode[j] * q_value;
				}
			}
		}

		linear_term = 0.0;
		quadratic_term = 0.0;
		for (ii = 0; ii < ns; ii++) {
			// i = subgraph->mothergraph_idx[ii];
			i = ii;
			linear_term -= g[ii] * x_gradient[i];
			quadratic_term += f[ii] * x_gradient[i];
		}
		if (b) {
			for (ii = 0; ii < ns; ii++) {
				// i = subgraph->mothergraph_idx[ii];
				i = ii;
				linear_term += x_gradient[i] * b[i];
			}
		}
		Free(f);
		Free(g);
	} else {
		val = -0.5 * SQR(delta) * quadratic_term + delta * linear_term;
		if (ai_store) {
			if (d) {
				tmp = 0.0;
				for (ii = 0; ii < ai_store->nd; ii++) {
					i = ai_store->d_idx[ii];
					loglFunc(&logll, &x[i], 1, i, x, loglFunc_arg);
					tmp += d[i] * logll;
				}
				val += tmp;
			}
		} else {
			if (d) {
				tmp = 0.0;
				for (ii = 0; ii < ns; ii++) {
					// i = subgraph->mothergraph_idx[ii];
					i = ii;
					if (d[i]) {
						loglFunc(&logll, &x[i], 1, i, x, loglFunc_arg);
						tmp += d[i] * logll;
					}
				}
				val += tmp;
			}
		}

		/*
		 * adjust if stochastic constraint 
		 */
		if (STOCHASTIC_CONSTR(constr)) {
			double sqr_term = 0.0;

			GMRFLib_EWRAP0(GMRFLib_eval_constr(NULL, &sqr_term, x, constr, graph));
			val += -0.5 * sqr_term;
		}
		*logdens = val;
	}

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}
int GMRFLib_ai_log_posterior_restricted_ORIG(double *logdens,
					     double *x, double *b, double *c, double *mean, double *d,
					     GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, char *fixed_value,
					     GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg,
					     GMRFLib_constr_tp * constr, GMRFLib_graph_tp * subgraph)
{
	/*
	 * this is the same function as GMRFLib_ai_log_posterior, BUT we only include those terms where at least one component
	 * is NOT FIXED.
	 * 
	 * the added last argument, subgraph, is added as an argument since its fixed for many repeated calls to this function.
	 * 
	 * if subgraph == NULL, then we use the GMRFLib_ai_log_posterior()-function 
	 */

	int i, j, ii, jj, ns;
	double *xx = NULL, val, tmp, logll = 0.0, sqr_term = 0.0;

	/*
	 * if subgraph is not available, use the default routine 
	 */
	if (!subgraph) {
		GMRFLib_EWRAP0(GMRFLib_ai_log_posterior(logdens, x, b, c, mean, d, loglFunc, loglFunc_arg, fixed_value, graph, Qfunc, Qfunc_arg, constr));
		return GMRFLib_SUCCESS;
	}
	GMRFLib_ENTER_ROUTINE;

	ns = subgraph->n;
	xx = Calloc(ns, double);			       /* xx = x - mean */
	if (mean) {
		for (ii = 0; ii < ns; ii++) {
			i = subgraph->mothergraph_idx[ii];
			xx[ii] = x[i] - mean[i];
		}
	} else {
		for (ii = 0; ii < ns; ii++) {
			i = subgraph->mothergraph_idx[ii];
			xx[ii] = x[i];
		}
	}

	sqr_term = 0.0;
	for (ii = 0; ii < ns; ii++) {
		i = subgraph->mothergraph_idx[ii];
		if (c) {
			sqr_term += SQR(xx[ii]) * (Qfunc(i, i, Qfunc_arg) + c[i]);
		} else {
			sqr_term += SQR(xx[ii]) * Qfunc(i, i, Qfunc_arg);
		}

		/*
		 * we have to compute the quadratic term like this, at least one of the components can be !fixed.
		 * 
		 * we need a term `2.0' for some cross-terms because we only loop over those `i' which are not fixed. 
		 */
		if (mean) {
			for (jj = 0; jj < graph->nnbs[i]; jj++) {
				j = graph->nbs[i][jj];
				if (fixed_value[j]) {
					sqr_term += 2.0 * xx[ii] * (x[j] - mean[j]) * Qfunc(i, j, Qfunc_arg);
				} else {
					sqr_term += xx[ii] * (x[j] - mean[j]) * Qfunc(i, j, Qfunc_arg);
				}
			}
		} else {
			for (jj = 0; jj < graph->nnbs[i]; jj++) {
				j = graph->nbs[i][jj];
				if (fixed_value[j]) {
					sqr_term += 2.0 * xx[ii] * x[j] * Qfunc(i, j, Qfunc_arg);
				} else {
					sqr_term += xx[ii] * x[j] * Qfunc(i, j, Qfunc_arg);
				}
			}
		}
	}
	val = -0.5 * sqr_term;				       /* val = -1/2 * x^T (Q + diag(c)) x */

	if (b) {
		tmp = 0.0;
		for (ii = 0; ii < ns; ii++) {
			i = subgraph->mothergraph_idx[ii];
			tmp += b[i] * x[i];
		}
		val += tmp;
	}
	if (d) {
		tmp = 0.0;
		for (ii = 0; ii < ns; ii++) {
			i = subgraph->mothergraph_idx[ii];
			if (d[i]) {
				loglFunc(&logll, &x[i], 1, i, x, loglFunc_arg);
				tmp += d[i] * logll;
			}
		}
		val += tmp;
	}

	/*
	 * adjust if stochastic constraint 
	 */
	if (STOCHASTIC_CONSTR(constr)) {
		GMRFLib_EWRAP0(GMRFLib_eval_constr(NULL, &sqr_term, x, constr, graph));
		val += -0.5 * sqr_term;
	}

	*logdens = val;

	Free(xx);
	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}
int GMRFLib_ai_nparam_eff(double *nparam_eff, double *nparam_eff_rel, GMRFLib_problem_tp * problem, double *c, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg)
{
	/*
	 * compute the effective number of paramerts using the DIC formula, p_D = p - trace( -P * V), where P is the Hessian of the log prior
	 * and V is the posterior covariance matrix.
	 * 
	 * the values returned are the effective number of paramters, which is computed either absolute (nparam_eff) or relative
	 * (nparam_eff_rel) 
	 */

	double correction = 0.0, *cov = NULL, n;
	int i, j, ii, jj, sub_n;

	sub_n = problem->sub_graph->n;

	for (i = 0; i < sub_n; i++) {

		ii = problem->map[i];
		cov = GMRFLib_Qinv_get(problem, ii, ii);
		GMRFLib_ASSERT(cov, GMRFLib_ESNH);
		correction += (Qfunc(ii, ii, Qfunc_arg) + (c ? c[ii] : 0.0)) * *cov;

		for (j = 0; j < problem->sub_graph->nnbs[i]; j++) {
			/*
			 * if the graph and Q match, then all covariances needed are computed, so the check ``if (cov)'' is for odd cases
			 * really. 
			 */
			jj = problem->map[problem->sub_graph->nbs[i][j]];
			cov = GMRFLib_Qinv_get(problem, ii, jj);
			if (cov) {
				correction += Qfunc(ii, jj, Qfunc_arg) * *cov;
			} else {
				/*
				 * do nothing 
				 */
			}
		}
	}

	n = sub_n - (problem->sub_constr ? problem->sub_constr->nc : 0.0);
	if (nparam_eff) {
		*nparam_eff = n - correction;
	}
	if (nparam_eff_rel) {
		*nparam_eff_rel = (n - correction) / n;
	}

	return GMRFLib_SUCCESS;
}
int GMRFLib_ai_do_MC_error_check(double *err, GMRFLib_problem_tp * problem, double *d, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, int nsamp)
{
	/*
	 * return the results in 'err', which is (at least) of dimension 5. err[0] : mean err[1] : (adjusted) stdev err[2] :
	 * lower 2.5% quantile err[3] : upper 2.5% quantile, err[4] : empirical mean of r*() not scaled with n.
	 */

	int i, ii, k, nsamples = 1000, sub_n, ndata = 0;
	double s[2], *a = NULL, *b = NULL, *c = NULL, loglikelihood, xval, diff, *samples = NULL;

	if (nsamp < 0) {
		nsamples = -nsamp;
	}

	a = Calloc(problem->n, double);
	b = Calloc(problem->n, double);
	c = Calloc(problem->n, double);

	sub_n = problem->sub_graph->n;
	for (ii = 0; ii < sub_n; ii++) {
		i = problem->map[ii];
		if (d && d[i]) {
			GMRFLib_2order_approx(&a[i], &b[i], &c[i], d[i], problem->mean_constr[i], i, problem->mean_constr, loglFunc, loglFunc_arg, NULL);
		}
	}

	/*
	 * count the number of observations 
	 */
	for (ii = 0; ii < sub_n; ii++) {
		i = problem->map[ii];
		if (d && d[i]) {
			ndata++;
		}
	}

	printf("ndata %d\n", ndata);
	if (ndata == 0) {
		err[0] = err[1] = err[2] = err[3] = err[4] = 0.0;
	} else {
		samples = Calloc(nsamples, double);

		s[0] = s[1] = 0.0;
		for (k = 0; k < nsamples; k++) {
			GMRFLib_sample(problem);
			diff = 0.0;
			for (ii = 0; ii < sub_n; ii++) {
				i = problem->map[ii];
				if (d && d[i]) {
					xval = problem->sample[i];
					loglFunc(&loglikelihood, &xval, 1, i, problem->sample, loglFunc_arg);
					diff += a[i] + (b[i] - 0.5 * c[i] * xval) * xval - d[i] * loglikelihood;
				}
			}
			diff /= ndata;
			samples[k] = diff;
			s[0] += diff;
			s[1] += SQR(diff);
		}

		s[0] /= nsamples;
		s[1] /= nsamples;

		/*
		 * correct the stdev estimate accoring to the student-t, so we can treat is as a `gaussian case' 
		 */
		double alpha = 0.025;

		err[0] = s[0];
		err[1] = sqrt(s[1] - SQR(s[0]));
		err[1] *= gsl_cdf_tdist_Pinv(alpha, nsamples - 1.0) / gsl_cdf_ugaussian_Pinv(alpha);
		err[4] = s[0] * ndata;

		/*
		 * compute the empirical quantiles 
		 */
		qsort((void *) samples, (size_t) nsamples, sizeof(double), GMRFLib_dcmp);
		err[2] = samples[(int) (alpha * (nsamples - 1.0))];
		err[3] = samples[(int) ((1.0 - alpha) * (nsamples - 1.0))];
	}

	Free(a);
	Free(b);
	Free(c);
	Free(samples);

	return GMRFLib_SUCCESS;
}
int GMRFLib_ai_marginal_hidden(GMRFLib_density_tp ** density, GMRFLib_density_tp ** cpo_density,
			       int idx, double *x, double *b, double *c, double *mean, double *d,
			       GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, char *fixed_value,
			       GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg,
			       GMRFLib_constr_tp * constr, GMRFLib_ai_param_tp * ai_par, GMRFLib_ai_store_tp * ai_store)
{
	/*
	 * compute the approximation to the marginal for the hidden field at index 'idx' and return the density in *density. if
	 * (cpo_density), then compute also the density when y_idx is removed.
	 */

	char *fix = NULL, *fixx = NULL;
	int i, j, k, nd = -1, n = -1, sub_n = -1, count, free_ai_par = 0, debug = 0, n_points, ns = -1, ii, free_ai_store = 0, *i_idx, *j_idx, one = 1;
	double *x_points = NULL, x_sd, x_mean, *cond_mode = NULL, *fixed_mode = NULL, *log_density = NULL,
	    log_dens_cond, deriv_log_dens_cond = 0.0, a, *derivative = NULL, *mean_and_variance = NULL, lc0, lc1, ld0, ld1, c0, c1, deldif =
	    GMRFLib_eps(1.0 / 3.0), h2 = 0.0, inv_stdev, *cov = NULL, corr, corr_term, *covariances = NULL, alpha;

	GMRFLib_graph_tp *subgraph = NULL;
	GMRFLib_Qinv_tp *store_Qinv = NULL;
	GMRFLib_optimize_param_tp *optpar = NULL;
	GMRFLib_blockupdate_param_tp *blockpar = NULL;
	GMRFLib_problem_tp *newp = NULL;
	GMRFLib_store_tp *store = NULL;

#define COMPUTE_CPO_DENSITY						\
	if (cpo_density) {						\
		if (d[idx]) {						\
			double *xp = NULL, *ld = NULL, *logcor = NULL, *x_user = NULL, *work = NULL, _alpha=-1.0; \
			int np = 35, _one = 1, _debug = 0;		\
									\
			work = Calloc(4*np, double); /* storage */	\
			xp = &work[0];					\
			ld = &work[np];					\
			logcor = &work[2*np];				\
			x_user = &work[3*np];				\
			GMRFLib_ghq_abscissas(&xp, np);			\
			GMRFLib_evaluate_nlogdensity(ld, xp, np, *density); \
			GMRFLib_density_std2user_n(x_user, xp, np, *density); \
			loglFunc(logcor, x_user, np, idx, fixed_mode, loglFunc_arg); \
			daxpy_(&np, &_alpha, logcor, &_one, ld, &_one); /* ld = ld - logcor */ \
			if (_debug && np && idx == 0) {			\
				int _i;					\
				for(_i = 0; _i < np; _i++)		\
					printf("CPO: %d %g %g\n", idx, xp[_i], ld[_i]);	\
				printf("CPO: \n");			\
			}						\
			GMRFLib_ai_correct_cpodens(ld, xp, &np, ai_par); \
			if (np > 4) {					\
				GMRFLib_density_create(cpo_density, GMRFLib_DENSITY_TYPE_SCGAUSSIAN, np, xp, ld, \
						       (*density)->std_mean, (*density)->std_stdev, GMRFLib_FALSE); \
				if (_debug) {				\
					P((*density)->std_mean);	\
					P((*density)->std_stdev);	\
				}					\
			} else {					\
				*cpo_density = NULL;			\
			}						\
			Free(work);					\
		} else {						\
			if (ai_par->cpo_manual){			\
				GMRFLib_density_duplicate(cpo_density, *density); \
			} else {					\
				*cpo_density = NULL;			\
			}						\
		}							\
	}

	GMRFLib_ENTER_ROUTINE;


	if (0) {
		FILE *fffp = fopen("constr.dat", "w");
		GMRFLib_print_constr(fffp, constr, graph);
		exit(0);
	}

	if (fixed_value && fixed_value[idx]) {
		*density = NULL;
		GMRFLib_LEAVE_ROUTINE;
		return GMRFLib_SUCCESS;			       /* no need request the density for a fixed value */
	}

	if (!ai_store) {				       /* use a temporary storage? */
		free_ai_store = 1;
		ai_store = Calloc(1, GMRFLib_ai_store_tp);
	}

	n = graph->n;

	/*
	 * put array where there is data into ai_store 
	 */
	if (d) {
		if (!ai_store->d_idx) {
#pragma omp critical
			if (!ai_store->d_idx) {
				for (i = nd = 0; i < n; i++) {
					if (d[i]) {
						nd++;
					}
				}
				ai_store->nd = nd;
				ai_store->d_idx = Calloc(nd, int);
				for (i = j = 0; i < n; i++) {
					if (d[i]) {
						ai_store->d_idx[j++] = i;
					}
				}
				assert(j == nd);
			}
		}
	} else {
		ai_store->nd = 0;
		ai_store->d_idx = NULL;
	}

	GMRFLib_default_optimize_param(&optpar);
	GMRFLib_default_blockupdate_param(&blockpar);
	if (!ai_par) {
		free_ai_par = 1;
		GMRFLib_default_ai_param(&ai_par);
	}

	/*
	 * if we use the ...CORRECTED_GAUSSIAN strategy, we have to use the linear correction. if the user have not selected an
	 * option for this, then do so depending if 'fast' option is set or not.
	 */
	if ((ai_par->strategy == GMRFLib_AI_STRATEGY_MEANCORRECTED_GAUSSIAN ||
	     ai_par->strategy == GMRFLib_AI_STRATEGY_MEANSKEWCORRECTED_GAUSSIAN) && ai_par->linear_correction == GMRFLib_AI_LINEAR_CORRECTION_OFF) {
		if (ai_par->fast) {
#pragma omp critical
			{
				if (ai_par->fast) {
					ai_par->linear_correction = GMRFLib_AI_LINEAR_CORRECTION_FAST;
				}
			}
		}
	}
	optpar->step_len = ai_par->step_len;
	blockpar->step_len = ai_par->step_len;
	blockpar->modeoption = GMRFLib_MODEOPTION_MODE;

	if (!(ai_store->problem)) {
		/*
		 * we need to store these two arrays for later usage 
		 */
		Free(ai_store->bb);
		Free(ai_store->cc);
		ai_store->bb = Calloc(n, double);
		ai_store->cc = Calloc(n, double);

		GMRFLib_EWRAP1(GMRFLib_init_GMRF_approximation_store__intern(&ai_store->problem,
									     (ai_store->mode ? ai_store->mode : x),
									     b, c, mean, d, loglFunc, loglFunc_arg, fixed_value,
									     graph, Qfunc, Qfunc_arg, constr, optpar, blockpar,
									     ai_store->store, ai_store->bb, ai_store->cc, ai_par->gaussian_data));
		GMRFLib_ai_add_Qinv_to_ai_store(ai_store);
		if (ai_par->compute_nparam_eff) {
			GMRFLib_ai_nparam_eff(&(ai_store->neff), NULL, ai_store->problem, c, Qfunc, Qfunc_arg);
		} else {
			ai_store->neff = GMRFLib_AI_STORE_NEFF_NOT_COMPUTED;
		}

		/*
		 * store the mode for later usage 
		 */
		Free(ai_store->mode);
		ai_store->mode = Calloc(n, double);
		memcpy(ai_store->mode, ai_store->problem->mean_constr, n * sizeof(double));
	} else {
		GMRFLib_ai_add_Qinv_to_ai_store(ai_store);
		if (ai_par->compute_nparam_eff && ai_store->neff == GMRFLib_AI_STORE_NEFF_NOT_COMPUTED) {
			GMRFLib_ai_nparam_eff(&(ai_store->neff), NULL, ai_store->problem, c, Qfunc, Qfunc_arg);
		}
	}

	/*
	 * ensure that the Qinv is present to be sure.
	 */
	GMRFLib_ai_add_Qinv_to_ai_store(ai_store);

	/*
	 * for internal use 
	 */
	x_mean = ai_store->problem->mean_constr[idx];
	x_sd = ai_store->stdev[idx];

	/*
	 * if just the Gaussian is needed, exit here 
	 */
	if (ai_par->strategy == GMRFLib_AI_STRATEGY_GAUSSIAN) {
		/*
		 * build the density-object 
		 */
		GMRFLib_density_create_normal(density, 0.0, 1.0, x_mean, x_sd);
		fixed_mode = ai_store->problem->mean_constr;

		COMPUTE_CPO_DENSITY;

		if (free_ai_store) {
			GMRFLib_free_ai_store(ai_store);
		}

		Free(blockpar);
		Free(optpar);

		GMRFLib_LEAVE_ROUTINE;
		return GMRFLib_SUCCESS;
	}

	if (fixed_value) {
		FIXME("Is this correct with FIXED_VALUE? FIXME.");
		abort();
	}
	if (!ai_store) {
		FIXME("I assume ai_store is !NULL.");
		abort();
	}

	sub_n = ai_store->problem->sub_graph->n;
	n_points = ai_par->n_points;
	GMRFLib_ghq_abscissas(&x_points, n_points);	       /* get the x-points */
	qsort(x_points, (size_t) n_points, sizeof(double), GMRFLib_dcmp_abs);	/* sort them using ABS() */

	log_density = Calloc(n_points, double);		       /* values of the log_density */
	cond_mode = Calloc(n, double);
	fixed_mode = Calloc(n, double);
	memcpy(cond_mode, ai_store->problem->mean_constr, n * sizeof(double));
	memcpy(fixed_mode, ai_store->problem->mean_constr, n * sizeof(double));

	fix = Calloc(n, char);
	fixx = Calloc(n, char);
	derivative = Calloc(n, double);

	/*
	 * the conditonal_mean is returned through 'derivative' 
	 */
	GMRFLib_ai_update_conditional_mean2(derivative, ai_store->problem, idx, x_mean + 1.0, &covariances);

	alpha = -1.0;
	daxpy_(&n, &alpha, fixed_mode, &one, derivative, &one);	/* derivative = derivative - fixed_mode */

	/*
	 * if we do not use the meancorrected gaussian and the fast-option, then locate local neigb. set the derivative to zero
	 * for those sites that are not in the local neigb.
	 */
	if ((ai_par->strategy == GMRFLib_AI_STRATEGY_MEANCORRECTED_GAUSSIAN || ai_par->strategy == GMRFLib_AI_STRATEGY_MEANSKEWCORRECTED_GAUSSIAN)
	    && ai_par->linear_correction == GMRFLib_AI_LINEAR_CORRECTION_FAST) {
		if (fixed_value) {
			memcpy(fix, fixed_value, n * sizeof(char));
			memcpy(fixx, fixed_value, n * sizeof(char));
		}
	} else {
		/*
		 * here we need the 'neigbours'.... 
		 */

		if (!ISZERO(ai_par->cutoff)) {
			if (fixed_value) {
				for (i = 0; i < n; i++) {
					if (!fixed_value[i]) {
						a = x_sd * derivative[i] / ai_store->stdev[i];
						if (ABS(a) < ai_par->cutoff) {
							fix[i] = fixx[i] = 1;
							derivative[i] = 0.0;
						}
					}
				}
			} else {
				for (i = 0; i < n; i++) {
					a = x_sd * derivative[i] / ai_store->stdev[i];
					if (ABS(a) < ai_par->cutoff) {
						fix[i] = fixx[i] = 1;
						derivative[i] = 0.0;
					}
				}
			}
		}
		if (fixed_value) {			       /* if there are fixed values already: add these */
			for (i = 0; i < n; i++) {
				if (fixed_value[i]) {
					fix[i] = fixx[i] = 1;
					derivative[i] = 0.0;
				}
			}
		}

		if (debug) {
			for (i = count = 0; i < n; i++)
				if (fix[i]) {
					count++;
				}
		}

		/*
		 * note that idx is included in subgraph 
		 */
		GMRFLib_EWRAP1(GMRFLib_compute_subgraph(&subgraph, graph, fixx));
		ns = subgraph->n;
	}

	fixx[idx] = 0;					       /* this is how 'fix' and 'fixx' differ */
	fix[idx] = 1;

	optpar->opt_type = GMRFLib_OPTTYPE_NR;		       /* force this option */
	store = Calloc(1, GMRFLib_store_tp);		       /* this can be used ;-) */

	if ((ai_par->linear_correction == GMRFLib_AI_LINEAR_CORRECTION_FAST) && !(ai_store->correction_term)) {
		double s = 1.0 / (2.0 * deldif);
		ai_store->correction_term = Calloc(n, double); /* compute this */
		ai_store->derivative3 = Calloc(n, double);     /* and this */
		ai_store->correction_idx = Calloc(n, int);     /* and this one */

		/*
		 * the same code splitted 
		 */
		ai_store->nidx = 0;
		if (fixed_value) {
			/*
			 * NOT IN USE 
			 */
			for (i = 0; i < n; i++) {
				if (d[i] && !fixed_value[i]) {
					ai_store->correction_idx[ai_store->nidx++] = i;
					GMRFLib_2order_approx(NULL, NULL, &c0, d[i], fixed_mode[i] - deldif, i,
							      fixed_mode, loglFunc, loglFunc_arg, &(ai_par->step_len));
					GMRFLib_2order_approx(NULL, NULL, &c1, d[i], fixed_mode[i] + deldif, i,
							      fixed_mode, loglFunc, loglFunc_arg, &(ai_par->step_len));
					ai_store->derivative3[i] = -(c1 - c0) * s;	/* `-' since c is negative 2.deriv */
					ai_store->correction_term[i] = -SQR(ai_store->stdev[i]) * ai_store->derivative3[i];
				}
			}
		} else {
			for (ii = 0; ii < ai_store->nd; ii++) {
				i = ai_store->d_idx[ii];
				ai_store->correction_idx[ai_store->nidx++] = i;
				GMRFLib_2order_approx(NULL, NULL, &c0, d[i], fixed_mode[i] - deldif, i, fixed_mode, loglFunc, loglFunc_arg, &(ai_par->step_len));
				GMRFLib_2order_approx(NULL, NULL, &c1, d[i], fixed_mode[i] + deldif, i, fixed_mode, loglFunc, loglFunc_arg, &(ai_par->step_len));
				ai_store->derivative3[i] = -(c1 - c0) * s;	/* `-' since c is negative 2.deriv */
				ai_store->correction_term[i] = -SQR(ai_store->stdev[i]) * ai_store->derivative3[i];
			}
		}
	}

	switch (ai_par->linear_correction) {
	case GMRFLib_AI_LINEAR_CORRECTION_OFF:
		break;					       /* including this; warnings in gcc */

	case GMRFLib_AI_LINEAR_CORRECTION_FAST:
		/*
		 * here, the term for i=idx will be zero, so we do not need to test for it.  it is a small hack to get the
		 * Qinv-values. Qinv_get only use the sub_inverse points, so we make a fake one and just retrive the elements we
		 * need. the version here, is somewhat faster than the old version below.
		 */
		deriv_log_dens_cond = 0.0;
		inv_stdev = 1.0 / x_sd;
		if (covariances) {
			for (j = 0; j < ai_store->nidx; j++) {
				i = ai_store->correction_idx[j];
				corr = covariances[i] * inv_stdev / ai_store->stdev[i];
				corr_term = 1.0 - SQR(corr);
				deriv_log_dens_cond += ai_store->correction_term[i] * corr_term * derivative[i];
				// printf("Using covariances Cov[%1d, %1d] = %.12f\n", idx, i, covariances[i]);
			}
			Free(covariances);
		} else {
			i_idx = map_ii_ptr(store_Qinv->mapping, idx);
			for (j = 0; j < ai_store->nidx; j++) {
				i = ai_store->correction_idx[j];
				j_idx = map_ii_ptr(store_Qinv->mapping, i);
				cov = map_id_ptr(store_Qinv->Qinv[IMIN(*i_idx, *j_idx)], IMAX(*i_idx, *j_idx));
				if (cov) {
					corr = *cov * inv_stdev / ai_store->stdev[i];
					corr_term = 1.0 - SQR(corr);
				} else {
					// corr_term = 0.9975; /* assume correlation = 0.05 for the rest */
					corr_term = 0.9375;    /* assume correlation = 0.25 for the rest */
				}
				deriv_log_dens_cond += ai_store->correction_term[i] * corr_term * derivative[i];
			}
		}
		deriv_log_dens_cond *= x_sd / 2.0;
		break;

	case GMRFLib_AI_LINEAR_CORRECTION_CENTRAL_DIFFERENCE:
		/*
		 * use central difference for the log-density. this involved computing the determinant of Q twice. use a large h
		 * to capture the large(r) scale behaviour.
		 */
		if (ai_par->fast) {
			blockpar->modeoption = GMRFLib_MODEOPTION_CURRENT;
		}
		h2 = 0.5 * x_sd;
		for (ii = 0; ii < ns; ii++) {
			i = subgraph->mothergraph_idx[ii];
			cond_mode[i] = fixed_mode[i] - h2 * derivative[i];
		}
		GMRFLib_EWRAP1(GMRFLib_init_GMRF_approximation_store
			       (&newp, cond_mode, b, c, mean, d, loglFunc, loglFunc_arg, fix, graph, Qfunc, Qfunc_arg, constr, optpar, blockpar, store));
		if (newp) {
			memcpy(newp->sample, newp->mean_constr, n * sizeof(double));
			GMRFLib_EWRAP1(GMRFLib_evaluate(newp));
			ld0 = newp->sub_logdens;
			lc0 = newp->log_normc;
			GMRFLib_free_problem(newp);
		} else {
			ld0 = lc0 = 0.0;
		}

		for (ii = 0; ii < ns; ii++) {
			i = subgraph->mothergraph_idx[ii];
			cond_mode[i] = fixed_mode[i] + h2 * derivative[i];
		}
		GMRFLib_EWRAP1(GMRFLib_init_GMRF_approximation_store
			       (&newp, cond_mode, b, c, mean, d, loglFunc, loglFunc_arg, fix, graph, Qfunc, Qfunc_arg, constr, optpar, blockpar, store));
		if (newp) {
			memcpy(newp->sample, newp->mean_constr, n * sizeof(double));
			GMRFLib_EWRAP1(GMRFLib_evaluate(newp));
			ld1 = newp->sub_logdens;
			lc1 = newp->log_normc;
			GMRFLib_free_problem(newp);
		} else {
			ld1 = lc1 = 0.0;
		}
		deriv_log_dens_cond = x_sd * (ld1 - ld0) / (2.0 * h2);
		// printf("DERIV %f %f ", deriv_log_dens_cond, x_sd*(lc1-lc0)/(2.0*h2));

		break;
	}

	if (ai_par->strategy != GMRFLib_AI_STRATEGY_MEANCORRECTED_GAUSSIAN && ai_par->strategy != GMRFLib_AI_STRATEGY_MEANSKEWCORRECTED_GAUSSIAN) {
		/*
		 * this is the main loop, where we evaluate different values for x_i. 
		 */

		int debug_me = 0;

		for (k = 0; k < n_points; k++) {
			/*
			 * find first the conditional mode: compute the initial value `cond_mode' 
			 */
			for (ii = 0; ii < ns; ii++) {
				cond_mode[ii] = fixed_mode[ii] + x_points[k] * x_sd * derivative[ii];
			}

			/*
			 * then find the conditional mode 
			 */
			if (ai_par->linear_correction != GMRFLib_AI_LINEAR_CORRECTION_OFF) {
				/*
				 * this option just use a linear correction ... 
				 */
				log_dens_cond = deriv_log_dens_cond * x_points[k];
				if (debug_me)
					printf("linear_correction x %f log_dens_cond %f\n", x_points[k], log_dens_cond);
			} else {
				/*
				 * ... while the default action is to compute the correction-term 
				 */
				if (ai_par->fast) {
					blockpar->modeoption = GMRFLib_MODEOPTION_CURRENT;
				}
				GMRFLib_EWRAP1(GMRFLib_init_GMRF_approximation_store
					       (&newp, cond_mode, b, c, mean, d, loglFunc, loglFunc_arg, fix, graph, Qfunc,
						Qfunc_arg, constr, optpar, blockpar, store));
				if (newp) {
					/*
					 * ....which is returned in newp->mean_constr 
					 */
					memcpy(newp->sample, newp->mean_constr, n * sizeof(double));
					memcpy(cond_mode, newp->mean_constr, n * sizeof(double));

					/*
					 * with density 
					 */
					GMRFLib_EWRAP1(GMRFLib_evaluate(newp));
					log_dens_cond = newp->sub_logdens;

					/*
					 * free this object 
					 */
					GMRFLib_free_problem(newp);
				} else {
					/*
					 * we end here if all nodes are fixed 
					 */
					log_dens_cond = 0.0;
				}

				if (debug_me)
					printf("true x %f log_dens_cond %f\n", x_points[k], log_dens_cond);
			}

			if (0) {
				/*
				 * this is the old version which is slow(er) 
				 */

				GMRFLib_EWRAP1(GMRFLib_ai_log_posterior_restricted_ORIG
					       (&log_density[k], cond_mode, b, c, mean, d, loglFunc, loglFunc_arg, fixx, graph,
						Qfunc, Qfunc_arg, constr, subgraph));
			} else {
				/*
				 * this is the fast version that take into account that x = x_mode + delta * gradient for the
				 * quadratic term
				 * 
				 * we first initialise the routine computing the linear and quadratic term, and then we can get
				 * the speedup for successive calls
				 *
				 * TODO: we should be able to speedup this, as we know what the linear and quadratic term is for
				 * the full conditional where the likelihood-term is included as well.
				 */
				if (k == 0) {
					GMRFLib_EWRAP1(GMRFLib_ai_log_posterior_restricted(NULL,
											   fixed_mode, fixed_mode, derivative,
											   0.0, b, c, mean, d, loglFunc,
											   loglFunc_arg, fixx, graph, Qfunc, Qfunc_arg, constr,
											   subgraph, ai_store));
				}
				GMRFLib_EWRAP1(GMRFLib_ai_log_posterior_restricted(&log_density[k],
										   cond_mode, fixed_mode, derivative,
										   x_points[k] * x_sd, b, c, mean, d, loglFunc,
										   loglFunc_arg, fixx, graph, Qfunc, Qfunc_arg, constr, subgraph, ai_store));
			}
			log_density[k] -= log_dens_cond;
		}
	}

	GMRFLib_free_store(store);

	if (ai_par->strategy == GMRFLib_AI_STRATEGY_MEANCORRECTED_GAUSSIAN) {
		GMRFLib_density_create_normal(density, -deriv_log_dens_cond, 1.0, x_mean, x_sd);
	} else if (ai_par->strategy == GMRFLib_AI_STRATEGY_MEANSKEWCORRECTED_GAUSSIAN) {
		int np = 11, err, fail = 0;
		double *ld = NULL, *xp = NULL, xx, low, high, third_order_derivative, a_sigma, cc, sol1, aa, tmp;
		GMRFLib_sn_param_tp snp;

		third_order_derivative = 0.0;
		for (j = 0; j < ai_store->nidx; j++) {
			i = ai_store->correction_idx[j];
			third_order_derivative += ai_store->derivative3[i] * gsl_pow_3(derivative[i]);
		}
		third_order_derivative *= gsl_pow_3(x_sd);

		/*
		 * match the third order derivative at the "mode" and then match the mean and variance with the meancorrected one. 
		 */
		if (ISZERO(third_order_derivative)) {
			/*
			 * in this case we are back to the Gaussian 
			 */
			snp.alpha = 0.0;
			snp.omega = 1.0;
			snp.xi = -deriv_log_dens_cond;
		} else {
			/*
			 * this require the Skew-Normal 
			 */
			a_sigma = GMRFLib_signed_pow(third_order_derivative / 0.2180136141449902, 1. / 3.);
			cc = 1.0 / a_sigma;
			err = GMRFLib_2order_poleq(&sol1, NULL, SQR(cc) * (1.0 - 2.0 / M_PI), SQR(cc) - 1.0, -1.0);
			if (err == GMRFLib_SUCCESS) {
				snp.alpha = aa = sqrt(sol1) * (third_order_derivative > 0.0 ? 1.0 : -1.0);
				tmp = 1.0 / (1.0 - (2.0 / M_PI) * SQR(aa) / (1.0 + SQR(aa)));
				if (tmp > 0.0) {
					snp.omega = sqrt(tmp);
					snp.xi = -deriv_log_dens_cond - snp.omega * sqrt(2.0 / M_PI) * aa / sqrt(1.0 + SQR(aa));
				} else {
					fail = 1;
				}
			} else {
				fail = 1;
			}
		}

		if (fail) {
			low = -deriv_log_dens_cond - 1.0;
			high = -deriv_log_dens_cond + 1.0;
			ld = Calloc(2 * np, double);	       /* xp = Calloc(np,double) */

			xp = &ld[np];

			for (k = 0; k < np; k++) {
				xp[k] = xx = low + k * (high - low) / (np - 1.0);
				ld[k] = -0.5 * SQR(xx) - deriv_log_dens_cond * xx + (third_order_derivative / 6.0) * gsl_pow_3(xx + deriv_log_dens_cond);
			}
			GMRFLib_EWRAP1(GMRFLib_sn_fit(&snp, NULL, xp, ld, np));
			Free(ld);
		}
		GMRFLib_density_create_sn(density, snp, x_mean, x_sd, GMRFLib_FALSE);
	} else {
		GMRFLib_density_create(density, GMRFLib_DENSITY_TYPE_SCGAUSSIAN, n_points, x_points, log_density, x_mean, x_sd, GMRFLib_FALSE);
	}

	COMPUTE_CPO_DENSITY;

	Free(fixed_mode);
	Free(fix);
	Free(fixx);
	Free(cond_mode);
	Free(derivative);
	Free(log_density);
	Free(optpar);
	Free(blockpar);
	if (free_ai_par) {
		Free(ai_par);
	}
	Free(mean_and_variance);
	GMRFLib_free_graph(subgraph);
	if (free_ai_store) {
		GMRFLib_free_ai_store(ai_store);
	}
#undef COMPUTE_CPO_DENSITY
	GMRFLib_LEAVE_ROUTINE;

	return GMRFLib_SUCCESS;
}
int GMRFLib_ai_update_conditional_mean(GMRFLib_problem_tp * pproblem, double *x, double *mean, char *fixed_value,
				       GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_args,
				       GMRFLib_constr_tp * constr, double *bbb, double *ccc, double **covariances, int idx)
{
	GMRFLib_ENTER_ROUTINE;

	if (0) {
		/*
		 * this is the fail-safe alternative, slower... 
		 */
		GMRFLib_constr_tp *c = NULL;

		GMRFLib_duplicate_constr(&c, constr, graph);
		GMRFLib_EWRAP1(GMRFLib_init_problem(&pproblem, x, bbb, ccc, mean, graph, Qfunc, Qfunc_args, fixed_value, c,
						    GMRFLib_KEEP_graph | GMRFLib_KEEP_chol | GMRFLib_UPDATE_constr));
		if (covariances) {
			*covariances = NULL;
		}
	} else {
		/*
		 * this is extracted from the problem_setup routine compute the new conditional mean using a rank-1 update assuming the old
		 * problem is the same. 
		 */

		GMRFLib_problem_tp **problem = &pproblem;
		int i, j, sub_n, nc, k, kk, inc = 1, n, use_old_code = 0;
		double alpha, beta, *aqat_m, *tmp_vector, *qi_at_m_store = NULL, *t_vector, *constr_m;

		assert((*problem)->sub_constr == constr);

		n = graph->n;
		sub_n = (*problem)->sub_graph->n;
		constr_m = (*problem)->constr_m;
		(*problem)->constr_m = NULL;
		Free((*problem)->l_aqat_m);

		nc = (*problem)->sub_constr->nc;	       /* shortname */
		if (nc == 1) {
			Free((*problem)->qi_at_m);
		} else {
			if ((*problem)->qi_at_m) {
				qi_at_m_store = (*problem)->qi_at_m;
				(*problem)->qi_at_m = NULL;
			}
		}

		(*problem)->qi_at_m = Calloc(nc * sub_n, double);
		memcpy((*problem)->qi_at_m, qi_at_m_store, (nc - 1) * sub_n * sizeof(double));
		Free(qi_at_m_store);

		/*
		 * this solves the equation for the last constraint only... 
		 */
		k = nc - 1;
		kk = k * sub_n;
		for (i = 0; i < sub_n; i++) {
			(*problem)->qi_at_m[i + kk] = (*problem)->sub_constr->a_matrix[k + nc * i];
		}

		if (1) {
			/*
			 * use the special solver which is faster just for this particular case where we know what the rhs is (its zero except
			 * for a one at index idx). 
			 */
			GMRFLib_EWRAP1(GMRFLib_solve_llt_sparse_matrix_special(&((*problem)->qi_at_m[kk]), &((*problem)->sub_sm_fact), (*problem)->sub_graph, idx));
		} else {
			/*
			 * or solve as usual 
			 */
			GMRFLib_EWRAP1(GMRFLib_solve_llt_sparse_matrix(&((*problem)->qi_at_m[kk]), &((*problem)->sub_sm_fact), (*problem)->sub_graph));
		}

		if (covariances) {
			int ii, idx_map;

			kk = (nc - 1) * sub_n;
			*covariances = Calloc(n, double);

			for (i = 0; i < sub_n; i++) {
				ii = (*problem)->map[i];
				(*covariances)[ii] = (*problem)->qi_at_m[kk + i];
			}

			if ((*problem)->sub_constr->nc > 1) {
				for (i = 0; i < sub_n; i++) {
					if (idx == (*problem)->map[i]) {
						break;
					}
				}
				idx_map = i;
				for (i = 0; i < sub_n; i++) {
					ii = (*problem)->map[i];
					for (kk = 0; kk < (*problem)->sub_constr->nc - 1; kk++) {
						(*covariances)[ii] -= constr_m[i + kk * sub_n] * (*problem)->qi_at_m[idx_map + kk * sub_n];
					}
				}
			}
		}

		/*
		 * compute l_aqat_m = chol(AQ^{-1}A^T)^{-1}) = chol(A qi_at_m)^{-1}, size = nc x nc 
		 */
		aqat_m = Calloc(nc * nc, double);

		alpha = 1.0;
		beta = 0.0;
		dgemm_("N", "N", &nc, &nc, &sub_n, &alpha, (*problem)->sub_constr->a_matrix, &nc, (*problem)->qi_at_m, &sub_n, &beta, aqat_m, &nc, 1, 1);
		if (STOCHASTIC_CONSTR((*problem)->sub_constr)) {
			/*
			 * add the covariance matrix, AQ^-1A^t + \Sigma 
			 */
			if ((*problem)->sub_constr->errcov_diagonal) {
				for (i = 0; i < nc; i++) {
					aqat_m[i + i * nc] += (*problem)->sub_constr->errcov_diagonal[i];
				}
			} else {
				for (i = 0; i < ISQR(nc); i++) {
					aqat_m[i] += (*problem)->sub_constr->errcov_general[i];
				}
			}
		}

		/*
		 * compute chol(aqat_m), recall that GMRFLib_comp_chol_general returns a new malloced L 
		 */
		GMRFLib_EWRAP1(GMRFLib_comp_chol_general(&((*problem)->l_aqat_m), aqat_m, nc, &((*problem)->logdet_aqat), GMRFLib_ESINGCONSTR));
		Free(aqat_m);

		/*
		 * ...and the constr-matrix Q^-1A^T inv(AQ^{-1}A^T + Sigma) 
		 */
		(*problem)->constr_m = Calloc(sub_n * nc, double);
		tmp_vector = Calloc(sub_n * nc, double);

		if (use_old_code) {
			for (i = 0, k = 0; i < sub_n; i++) {
				for (j = 0; j < nc; j++) {
					tmp_vector[k++] = (*problem)->qi_at_m[i + j * sub_n];
				}
			}
		} else {
			switch (nc) {
			case 1:
				for (i = 0, k = 0; i < sub_n; i++, k++) {
					tmp_vector[k] = (*problem)->qi_at_m[i];
				}
				break;
			case 2:
				for (i = 0, k = 0; i < sub_n; i++, k += 2) {
					tmp_vector[k] = (*problem)->qi_at_m[i];
					tmp_vector[k + 1] = (*problem)->qi_at_m[i + sub_n];
				}
				break;
			case 3:
				for (i = 0, k = 0; i < sub_n; i++, k += 3) {
					tmp_vector[k] = (*problem)->qi_at_m[i];
					tmp_vector[k + 1] = (*problem)->qi_at_m[i + sub_n];
					tmp_vector[k + 2] = (*problem)->qi_at_m[i + 2 * sub_n];
				}
				break;
			case 4:
				for (i = 0, k = 0; i < sub_n; i++, k += 4) {
					tmp_vector[k] = (*problem)->qi_at_m[i];
					tmp_vector[k + 1] = (*problem)->qi_at_m[i + sub_n];
					tmp_vector[k + 2] = (*problem)->qi_at_m[i + 2 * sub_n];
					tmp_vector[k + 3] = (*problem)->qi_at_m[i + 3 * sub_n];
				}
				break;
			case 5:
				for (i = 0, k = 0; i < sub_n; i++, k += 5) {
					tmp_vector[k] = (*problem)->qi_at_m[i];
					tmp_vector[k + 1] = (*problem)->qi_at_m[i + sub_n];
					tmp_vector[k + 2] = (*problem)->qi_at_m[i + 2 * sub_n];
					tmp_vector[k + 3] = (*problem)->qi_at_m[i + 3 * sub_n];
					tmp_vector[k + 4] = (*problem)->qi_at_m[i + 4 * sub_n];
				}
				break;
			case 6:
				for (i = 0, k = 0; i < sub_n; i++, k += 6) {
					tmp_vector[k] = (*problem)->qi_at_m[i];
					tmp_vector[k + 1] = (*problem)->qi_at_m[i + sub_n];
					tmp_vector[k + 2] = (*problem)->qi_at_m[i + 2 * sub_n];
					tmp_vector[k + 3] = (*problem)->qi_at_m[i + 3 * sub_n];
					tmp_vector[k + 4] = (*problem)->qi_at_m[i + 4 * sub_n];
					tmp_vector[k + 5] = (*problem)->qi_at_m[i + 5 * sub_n];
				}
				break;
			case 7:
				for (i = 0, k = 0; i < sub_n; i++, k += 7) {
					tmp_vector[k] = (*problem)->qi_at_m[i];
					tmp_vector[k + 1] = (*problem)->qi_at_m[i + sub_n];
					tmp_vector[k + 2] = (*problem)->qi_at_m[i + 2 * sub_n];
					tmp_vector[k + 3] = (*problem)->qi_at_m[i + 3 * sub_n];
					tmp_vector[k + 4] = (*problem)->qi_at_m[i + 4 * sub_n];
					tmp_vector[k + 5] = (*problem)->qi_at_m[i + 5 * sub_n];
					tmp_vector[k + 6] = (*problem)->qi_at_m[i + 6 * sub_n];
				}
				break;
			case 8:
				for (i = 0, k = 0; i < sub_n; i++, k += 8) {
					tmp_vector[k] = (*problem)->qi_at_m[i];
					tmp_vector[k + 1] = (*problem)->qi_at_m[i + sub_n];
					tmp_vector[k + 2] = (*problem)->qi_at_m[i + 2 * sub_n];
					tmp_vector[k + 3] = (*problem)->qi_at_m[i + 3 * sub_n];
					tmp_vector[k + 4] = (*problem)->qi_at_m[i + 4 * sub_n];
					tmp_vector[k + 5] = (*problem)->qi_at_m[i + 5 * sub_n];
					tmp_vector[k + 6] = (*problem)->qi_at_m[i + 6 * sub_n];
					tmp_vector[k + 7] = (*problem)->qi_at_m[i + 7 * sub_n];
				}
				break;
			default:
				for (i = 0, k = 0; i < sub_n; i++) {
					for (j = 0; j < nc; j++) {
						tmp_vector[k++] = (*problem)->qi_at_m[i + j * sub_n];
					}
				}
			}
		}

		GMRFLib_EWRAP1(GMRFLib_solveAxb_posdef(tmp_vector, (*problem)->l_aqat_m, tmp_vector, nc, sub_n));

		if (use_old_code) {
			for (i = 0, k = 0; i < sub_n; i++) {
				for (j = 0; j < nc; j++) {
					(*problem)->constr_m[i + j * sub_n] = tmp_vector[k++];
				}
			}
		} else {
			switch (nc) {
			case 1:
				for (i = 0, k = 0; i < sub_n; i++, k++) {
					(*problem)->constr_m[i] = tmp_vector[k];
				}
				break;
			case 2:
				for (i = 0, k = 0; i < sub_n; i++, k += 2) {
					(*problem)->constr_m[i] = tmp_vector[k];
					(*problem)->constr_m[i + sub_n] = tmp_vector[k + 1];
				}
				break;
			case 3:
				for (i = 0, k = 0; i < sub_n; i++, k += 3) {
					(*problem)->constr_m[i] = tmp_vector[k];
					(*problem)->constr_m[i + sub_n] = tmp_vector[k + 1];
					(*problem)->constr_m[i + 2 * sub_n] = tmp_vector[k + 2];
				}
				break;
			case 4:
				for (i = 0, k = 0; i < sub_n; i++, k += 4) {
					(*problem)->constr_m[i] = tmp_vector[k];
					(*problem)->constr_m[i + sub_n] = tmp_vector[k + 1];
					(*problem)->constr_m[i + 2 * sub_n] = tmp_vector[k + 2];
					(*problem)->constr_m[i + 3 * sub_n] = tmp_vector[k + 3];
				}
				break;
			case 5:
				for (i = 0, k = 0; i < sub_n; i++, k += 5) {
					(*problem)->constr_m[i] = tmp_vector[k];
					(*problem)->constr_m[i + sub_n] = tmp_vector[k + 1];
					(*problem)->constr_m[i + 2 * sub_n] = tmp_vector[k + 2];
					(*problem)->constr_m[i + 3 * sub_n] = tmp_vector[k + 3];
					(*problem)->constr_m[i + 4 * sub_n] = tmp_vector[k + 4];
				}
				break;
			case 6:
				for (i = 0, k = 0; i < sub_n; i++, k += 6) {
					(*problem)->constr_m[i] = tmp_vector[k];
					(*problem)->constr_m[i + sub_n] = tmp_vector[k + 1];
					(*problem)->constr_m[i + 2 * sub_n] = tmp_vector[k + 2];
					(*problem)->constr_m[i + 3 * sub_n] = tmp_vector[k + 3];
					(*problem)->constr_m[i + 4 * sub_n] = tmp_vector[k + 4];
					(*problem)->constr_m[i + 5 * sub_n] = tmp_vector[k + 5];
				}
				break;
			case 7:
				for (i = 0, k = 0; i < sub_n; i++, k += 7) {
					(*problem)->constr_m[i] = tmp_vector[k];
					(*problem)->constr_m[i + sub_n] = tmp_vector[k + 1];
					(*problem)->constr_m[i + 2 * sub_n] = tmp_vector[k + 2];
					(*problem)->constr_m[i + 3 * sub_n] = tmp_vector[k + 3];
					(*problem)->constr_m[i + 4 * sub_n] = tmp_vector[k + 4];
					(*problem)->constr_m[i + 5 * sub_n] = tmp_vector[k + 5];
					(*problem)->constr_m[i + 6 * sub_n] = tmp_vector[k + 6];
				}
				break;
			case 8:
				for (i = 0, k = 0; i < sub_n; i++, k += 8) {
					(*problem)->constr_m[i] = tmp_vector[k];
					(*problem)->constr_m[i + sub_n] = tmp_vector[k + 1];
					(*problem)->constr_m[i + 2 * sub_n] = tmp_vector[k + 2];
					(*problem)->constr_m[i + 3 * sub_n] = tmp_vector[k + 3];
					(*problem)->constr_m[i + 4 * sub_n] = tmp_vector[k + 4];
					(*problem)->constr_m[i + 5 * sub_n] = tmp_vector[k + 5];
					(*problem)->constr_m[i + 6 * sub_n] = tmp_vector[k + 6];
					(*problem)->constr_m[i + 7 * sub_n] = tmp_vector[k + 7];
				}
				break;
			default:
				for (i = 0, k = 0; i < sub_n; i++) {
					for (j = 0; j < nc; j++) {
						(*problem)->constr_m[i + j * sub_n] = tmp_vector[k++];
					}
				}
			}
		}

		Free(tmp_vector);
		nc = constr->nc;
		memcpy((*problem)->sub_mean_constr, (*problem)->sub_mean, sub_n * sizeof(double));
		t_vector = Calloc(nc, double);

		GMRFLib_EWRAP1(GMRFLib_eval_constr(t_vector, NULL, (*problem)->sub_mean, (*problem)->sub_constr, (*problem)->sub_graph));

		alpha = -1.0;
		beta = 1.0;				       /* mean_constr = mean - cond_m*t_vector */
		dgemv_("N", &sub_n, &nc, &alpha, (*problem)->constr_m, &sub_n, t_vector, &inc, &beta, (*problem)->sub_mean_constr, &inc, 1);
		Free(t_vector);

		for (i = 0; i < sub_n; i++) {
			j = (*problem)->map[i];
			(*problem)->mean_constr[j] = (*problem)->sub_mean_constr[i];
		}

		Free(constr_m);
	}

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}
int GMRFLib_ai_update_conditional_mean2(double *cond_mean, GMRFLib_problem_tp * problem, int idx, double evalue, double **covariances)
{
	/*
	 * this is version 2 (or 3) of the same routine, which use the recursive formula for computing the new condition mean without altering the contents of
	 * PROBLEM.
	 * 
	 * assume FIXED == NULL
	 * 
	 */

/* 
   this one is FIXED by design and such that A[IDX(i,j,n)] = A_ij, i=0...n-1, j = 0..k-1 for n x k matrix A 
 */
#define IDX(i, j, n) ((i) + (j)*(n))
#define WORK(n) &work[work_p]; work_p += (n)

	int i, k, n, nc, ncc, one = 1, work_p, work_size, ndiv;
	double *c = NULL, *v = NULL, *w = NULL, *z = NULL, alpha = 0.0, beta = 0.0, b22 = 0.0, *constr_m_new = NULL, *t_vec = NULL, *work = NULL, *tmp_m =
	    NULL, val;

	GMRFLib_ENTER_ROUTINE;

	n = problem->n;
	assert(n == problem->sub_graph->n);
	nc = (problem->sub_constr ? problem->sub_constr->nc : 0);
	ncc = nc + 1;

	/*
	 * setup workspace for small-mem's for the hole routine here. 
	 */
	work_size = n + ncc + (nc ? n + nc + nc + ISQR(nc) : 0);
	work = Calloc(work_size, double);
	work_p = 0;
	c = WORK(n);
	t_vec = WORK(ncc);
	if (nc) {
		z = WORK(n);
		v = WORK(nc);
		w = WORK(nc);
		tmp_m = WORK(ISQR(nc));
	}
	assert(work_p == work_size);

	c[idx] = 1.0;
	GMRFLib_solve_llt_sparse_matrix_special(c, &(problem->sub_sm_fact), problem->sub_graph, idx);

	if (covariances) {
		*covariances = Calloc(n, double);
		memcpy(*covariances, c, n * sizeof(double));

		if (nc) {
			for (k = 0; k < nc; k++) {
				/*
				 * replacement for
				 * 
				 * val = problem->qi_at_m[IDX(idx, k, n)]; for (i = 0; i < n; i++) (*covariances)[i] -= problem->constr_m[IDX(i, k, n)] * val; 
				 */
				alpha = -problem->qi_at_m[IDX(idx, k, n)];
				daxpy_(&n, &alpha, &(problem->constr_m[IDX(0, k, n)]), &one, &((*covariances)[0]), &one);
			}
		}
	}

	constr_m_new = Calloc(n * ncc, double);
	if (nc) {
		/*
		 * add inv(A Q^-1 A^t) if it does not exists. Be careful, as we need to add this and that this routine can be called threaded with the same
		 * ai_store. 
		 */
		if (!problem->inv_aqat_m) {
#pragma omp critical
			{
				if (!problem->inv_aqat_m) {
					double *m;

					m = Calloc(ISQR(nc), double);
					alpha = 1.0;
					beta = 0.0;
					dgemm_("N", "N", &nc, &nc, &n, &alpha, problem->sub_constr->a_matrix, &nc, problem->qi_at_m, &n, &beta, m, &nc, 1, 1);
					GMRFLib_comp_posdef_inverse(m, nc);
					problem->inv_aqat_m = m;
				}
			}
		}

		/*
		 * v = A c.  w = inv(AQ^-1A^T) v.  z = (Q^-1 A^T) w. 
		 */
		alpha = 1.0;
		beta = 0.0;
		dgemv_("N", &nc, &n, &alpha, problem->sub_constr->a_matrix, &nc, c, &one, &beta, v, &one, 1);
		dgemv_("N", &nc, &nc, &alpha, problem->inv_aqat_m, &nc, v, &one, &beta, w, &one, 1);
		dgemv_("N", &n, &nc, &alpha, problem->qi_at_m, &n, w, &one, &beta, z, &one, 1);


		/*
		 * replacement for:: b22 = c[idx]; for (i = 0; i < nc; i++) b22 -= v[i] * w[i]; b22 = 1.0 / b22; 
		 */

		b22 = 1.0 / (c[idx] - ddot_(&nc, v, &one, w, &one));
		assert(b22 > 0);

		for (k = 0; k < nc; k++) {
			val = b22 * w[k];
			daxpy_(&nc, &val, v, &one, tmp_m, &one);
			tmp_m[IDX(k, k, nc)] += 1.0;
		}

		alpha = 1.0;
		beta = 0.0;
		dgemm_("N", "N", &n, &nc, &nc, &alpha, problem->constr_m, &n, tmp_m, &nc, &beta, constr_m_new, &n, 1, 1);

		for (k = 0; k < nc; k++) {
			val = -b22 * w[k];
			daxpy_(&n, &val, c, &one, &constr_m_new[IDX(0, k, n)], &one);
		}

		k = (ncc - 1) * n;
		ndiv = 4 * (n / 4);
		for (i = 0; i < ndiv; i += 4) {
			constr_m_new[i + k] = b22 * (c[i] - z[i]);
			constr_m_new[i + k + 1] = b22 * (c[i + 1] - z[i + 1]);
			constr_m_new[i + k + 2] = b22 * (c[i + 2] - z[i + 2]);
			constr_m_new[i + k + 3] = b22 * (c[i + 3] - z[i + 3]);
		}
		for (i = ndiv; i < n; i++) {
			constr_m_new[i + k] = b22 * (c[i] - z[i]);
		}

		// No longer needed: GMRFLib_eval_constr(t_vec, NULL, problem->sub_mean, problem->sub_constr, problem->sub_graph);
		assert(problem->sub_constr_value);
		memcpy(t_vec, problem->sub_constr_value, nc * sizeof(double));
	} else {
		double cinv = 1.0 / c[idx];

		// for(i=0; i<n; i++) constr_m_new[k + i] = cinv*c[i];
		daxpy_(&n, &cinv, c, &one, &constr_m_new[nc * n], &one);
	}

	t_vec[ncc - 1] = problem->sub_mean[idx] - evalue;

	alpha = -1.0;
	beta = 1.0;					       /* mean_constr = mean - cond_m*t_vec */
	memcpy(cond_mean, problem->sub_mean, n * sizeof(double));
	dgemv_("N", &n, &ncc, &alpha, constr_m_new, &n, t_vec, &one, &beta, cond_mean, &one, 1);

	Free(work);
	Free(constr_m_new);
	GMRFLib_LEAVE_ROUTINE;

#undef IDX
#undef WORK
	return GMRFLib_SUCCESS;
}
int GMRFLib_free_ai_store(GMRFLib_ai_store_tp * ai_store)
{
	/*
	 * free that part of 'store' that is used in this module. 
	 */

	if (ai_store) {
		GMRFLib_free_store(ai_store->store);
		GMRFLib_free_problem(ai_store->problem);
		Free(ai_store->mode);
		Free(ai_store->bb);
		Free(ai_store->cc);
		Free(ai_store->stdev);
		Free(ai_store->d_idx);
		Free(ai_store->correction_term);
		Free(ai_store->correction_idx);
		Free(ai_store->derivative3);
		Free(ai_store);
	}
	return GMRFLib_SUCCESS;
}
int GMRFLib_ai_z2theta(double *theta, int nhyper, double *theta_mode, double *z, gsl_vector * eigen_values, gsl_matrix * eigen_vectors)
{
	/*
	 * compute new theta-values for given vector of z (which is N(0,I)), using the relationship
	 * 
	 * theta = theta_mode + eigen_vectors * diag(1/sqrt(eigen_values)) * z 
	 */

	size_t i, j;
	double tmp, v_ij, lambda, *u = NULL;

	u = Calloc(nhyper, double);

	for (i = 0; i < (size_t) nhyper; i++) {
		lambda = gsl_vector_get(eigen_values, i);
		u[i] = z[i] / sqrt(lambda);
	}

	for (i = 0; i < (size_t) nhyper; i++) {
		for (j = 0, tmp = 0.0; j < (size_t) nhyper; j++) {
			v_ij = gsl_matrix_get(eigen_vectors, i, j);
			tmp += v_ij * u[j];
		}
		theta[i] = theta_mode[i] + tmp;
	}
	Free(u);

	return GMRFLib_SUCCESS;
}
int GMRFLib_ai_theta2z(double *z, int nhyper, double *theta_mode, double *theta, gsl_vector * eigen_values, gsl_matrix * eigen_vectors)
{
	/*
	 * compute z-values for given vector of theta, using the relationship
	 * 
	 * theta = theta_mode + eigen_vectors * diag(1/sqrt(eigen_values)) * z 
	 */

	size_t i, j;
	double tmp, v_ji, lambda, *u = NULL;

	u = Calloc(nhyper, double);

	for (i = 0; i < (size_t) nhyper; i++) {
		u[i] = theta[i] - theta_mode[i];
	}

	for (i = 0; i < (size_t) nhyper; i++) {
		for (j = 0, tmp = 0.0; j < (size_t) nhyper; j++) {
			v_ji = gsl_matrix_get(eigen_vectors, j, i);
			tmp += v_ji * u[j];
		}
		lambda = gsl_vector_get(eigen_values, i);
		z[i] = tmp * sqrt(lambda);
	}
	Free(u);

	return GMRFLib_SUCCESS;
}
int GMRFLib_init_GMRF_approximation_store__intern(GMRFLib_problem_tp ** problem, double *x, double *b, double *c, double *mean,
						  double *d, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, char *fixed_value,
						  GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg,
						  GMRFLib_constr_tp * constr, GMRFLib_optimize_param_tp * optpar,
						  GMRFLib_blockupdate_param_tp * blockupdate_par, GMRFLib_store_tp * store, double *bb, double *cc,
						  int gaussian_data)
{
	/*
	 * This is copy of the original routine but with optional two last arguments 
	 */

	int i, free_x = 0, free_b = 0, free_c = 0, free_mean = 0, free_d = 0, free_blockpar = 0, free_bb = 0, free_cc = 0, n, id,
	    *idxs = NULL, nidx = 0, use_old_code = 0;
	double *mode = NULL;
	static int new_idea = 0;

	if (new_idea == 99) {
		if (getenv("INLA_NEW_IDEA")) {
			new_idea = 1;
		} else {
			new_idea = 0;
		}
		FIXME("READ INLA_NEW_IDEA: ");
		P(new_idea);
	}
#define FREE_ALL if (1) { if (free_x) Free(x); if (free_b) Free(b); if (free_c) Free(c); if (free_d) Free(d); \
		if (free_mean) Free(mean); if (free_blockpar) Free(blockupdate_par); if (free_bb) Free(bb); \
		if (free_cc) Free(cc); Free(mode); Free(idxs); }

	GMRFLib_ENTER_ROUTINE;

	id = GMRFLib_thread_id;
	n = graph->n;
	if (n == 0) {
		*problem = NULL;
		return GMRFLib_SUCCESS;
	}

	if (!x) {
		free_x = 1;
		x = Calloc(n, double);
	}
	if (!b) {
		free_b = 1;
		b = Calloc(n, double);
	}
	if (!c) {
		free_c = 1;
		c = Calloc(n, double);
	}
	if (!d) {
		free_d = 1;
		d = Calloc(n, double);
	}
	if (!mean) {
		free_mean = 1;
		mean = Calloc(n, double);
	}
	if (!bb) {
		free_bb = 1;
		bb = Calloc(n, double);
	}
	if (!cc) {
		free_cc = 1;
		cc = Calloc(n, double);
	}
	mode = Calloc(n, double);
	memcpy(mode, x, n * sizeof(double));

	if (use_old_code || (optpar && optpar->opt_type != GMRFLib_OPTTYPE_NR)) {
		/*
		 * this is the old version 
		 */

		FIXME("OLD VERSION");

		if (!blockupdate_par) {
			GMRFLib_default_blockupdate_param(&blockupdate_par);
			free_blockpar = 1;
		}

		if (blockupdate_par->modeoption == GMRFLib_MODEOPTION_MODE && d) {
			GMRFLib_EWRAP1(GMRFLib_optimize_store
				       (mode, b, c, mean, graph, Qfunc, Qfunc_arg, fixed_value, constr, d, loglFunc, loglFunc_arg, optpar, store));
		}

		/*
		 * compute the terms from loglFunc 
		 */
		if (d) {
#pragma omp parallel for private(i) schedule(static)
			for (i = 0; i < n; i++) {
				if (d[i] && (!fixed_value || !fixed_value[i])) {
					GMRFLib_thread_id = id;
					GMRFLib_2order_approx(NULL, &bb[i], &cc[i], d[i], mode[i], i, mode, loglFunc, loglFunc_arg, &(blockupdate_par->step_len));
					cc[i] = DMAX(0.0, cc[i]);	/* do not want negative terms on the diagonal */
				}
			}
			GMRFLib_thread_id = id;
		}

		if (b) {
			if (mean) {
				for (i = 0; i < n; i++) {
					bb[i] += b[i] - cc[i] * mean[i];
				}
			} else {
				for (i = 0; i < n; i++) {
					bb[i] += b[i];
				}
			}
		} else {
			if (mean) {
				for (i = 0; i < n; i++) {
					bb[i] -= cc[i] * mean[i];
				}
			}
		}
		if (c) {
			for (i = 0; i < n; i++) {
				cc[i] += c[i];
			}
		}
		GMRFLib_EWRAP1(GMRFLib_init_problem_store(problem, x, bb, cc, mean, graph, Qfunc, Qfunc_arg, fixed_value, constr, GMRFLib_NEW_PROBLEM, store));
	} else {
		/*
		 * the NEW implementation which do optimisation and GMRF_approximation in the same operation. this is tailored for INLA of'course, and only ment
		 * for that. 
		 */

		if (optpar && optpar->fp)
			fprintf(optpar->fp, "\nComputing GMRF approximation\n------------------------------\n");
		nidx = 0;
		idxs = Calloc(n, int);
		for (i = 0; i < n; i++) {
			if (d[i]) {
				idxs[nidx++] = i;
			}
		}

		int iter, itmax = optpar->max_iter;

		for (iter = 0; iter < itmax; iter++) {

			memcpy(bb, b, n * sizeof(double));
			memcpy(cc, c, n * sizeof(double));

#pragma omp parallel for private(i) schedule(static)
			for (i = 0; i < nidx; i++) {
				int idx;
				double bcoof, ccoof;

				GMRFLib_thread_id = id;
				idx = idxs[i];
				GMRFLib_2order_approx(NULL, &bcoof, &ccoof, d[idx], mode[idx], idx, mode, loglFunc, loglFunc_arg, &(optpar->step_len));
				bb[idx] += bcoof;
				cc[idx] += DMAX(0.0, ccoof);
			}
			GMRFLib_thread_id = id;

			for (i = 0; i < n; i++) {
				bb[i] += -c[i] * mean[i];
				cc[i] += c[i];
			}
			GMRFLib_EWRAP1(GMRFLib_init_problem_store(problem, x, bb, cc, mean, graph, Qfunc, Qfunc_arg, fixed_value, constr,
								  GMRFLib_NEW_PROBLEM, store));

			double err = 0.0, f;
			f = DMIN(1.0, (iter + 1.0) * optpar->nr_step_factor);
			for (i = 0; i < n; i++) {
				err += SQR((*problem)->mean_constr[i] - mode[i]);
				mode[i] += f * ((*problem)->mean_constr[i] - mode[i]);
			}
			err = sqrt(err / n);
			if (optpar && optpar->fp)
				fprintf(optpar->fp, "iteration %d error %.6g\n", iter, err);

			if (err < optpar->abserr_step || gaussian_data) {
				/*
				 * we're done! 
				 */
				break;
			}

			if (0) {
				if (err < 1.0) {
					/*
					 * NOT PROPERLY TESTED!!!! do one step more step without touching Q and its factorisation.
					 */
					memcpy(bb, b, n * sizeof(double));
#pragma omp parallel for private(i) schedule(static)
					for (i = 0; i < nidx; i++) {
						int idx;
						double bcoof;

						GMRFLib_thread_id = id;
						idx = idxs[i];
						GMRFLib_2order_approx(NULL, &bcoof, NULL, d[idx], (*problem)->mean_constr[idx], idx,
								      (*problem)->mean_constr, loglFunc, loglFunc_arg, &(optpar->step_len));
						bb[idx] += bcoof;
					}
					GMRFLib_thread_id = id;
					for (i = 0; i < n; i++) {
						bb[i] += -c[i] * mean[i];
					}
					GMRFLib_EWRAP1(GMRFLib_init_problem_store(problem, x, bb, cc, mean, graph, Qfunc, Qfunc_arg, fixed_value, constr,
										  GMRFLib_KEEP_graph | GMRFLib_KEEP_constr | GMRFLib_KEEP_chol, NULL));
					for (i = 0; i < n; i++) {
						mode[i] = (*problem)->mean_constr[i];
					}
				}
			}
			GMRFLib_free_problem(*problem);
			*problem = NULL;

			if (gsl_isnan(err))
				break;
		}
		if (!*problem) {
			/*
			 * fail to converge. restart with a reduced step_factor. 
			 */
			FREE_ALL;
			optpar->nr_step_factor /= 10.0;
			optpar->max_iter *= 2;
			if (optpar && optpar->fp) {
				fprintf(optpar->fp, "\n\n%s: Optimisation fail to converge.\n\t\t\tRetry with a new optpar->nr_step_factor = %f\n",
					__GMRFLib_FuncName, optpar->nr_step_factor);
			}
			if (optpar->nr_step_factor < 0.00001) {
				return GMRFLib_EOPTNR;
			} else {
				return GMRFLib_init_GMRF_approximation_store__intern(problem, x, b, c, mean, d,
										     loglFunc, loglFunc_arg, fixed_value, graph, Qfunc, Qfunc_arg,
										     constr, optpar, blockupdate_par, store, bb, cc, gaussian_data);
			}
		}
		if (*problem && new_idea) {
			double *sd = Calloc(n, double), fac = 3, step_len, err = 0.0, f, itmax_local = 40;

			for (iter = 0; iter < itmax; iter++) {

				GMRFLib_Qinv(*problem, GMRFLib_QINV_ALL);

				for (i = 0; i < n; i++) {
					double *p = GMRFLib_Qinv_get(*problem, i, i);
					sd[i] = (p ? sqrt(*p) : 0.0);
				}

				memcpy(bb, b, n * sizeof(double));
				memcpy(cc, c, n * sizeof(double));

#pragma omp parallel for private(i) schedule(static)
				for (i = 0; i < nidx; i++) {
					int idx;
					double bcoof, ccoof;

					GMRFLib_thread_id = id;
					idx = idxs[i];
					step_len = -(fac * sd[idx]);	/* yes! */
					GMRFLib_2order_approx(NULL, &bcoof, &ccoof, d[idx], mode[idx], idx, mode, loglFunc, loglFunc_arg, &step_len);
					bb[idx] += bcoof;
					cc[idx] += DMAX(0.0, ccoof);
				}
				GMRFLib_thread_id = id;

				for (i = 0; i < n; i++) {
					bb[i] += -c[i] * mean[i];
					cc[i] += c[i];
				}
				GMRFLib_EWRAP1(GMRFLib_init_problem_store(problem, x, bb, cc, mean, graph, Qfunc, Qfunc_arg, fixed_value, constr,
									  GMRFLib_NEW_PROBLEM, store));

				f = DMIN(1.0, (iter + 1.0) * optpar->nr_step_factor);
				for (i = 0; i < n; i++) {
					err += SQR((*problem)->mean_constr[i] - mode[i]);
					mode[i] += f * ((*problem)->mean_constr[i] - mode[i]);
				}
				err = sqrt(err / n);
				printf("PART 2: iter %d err %g\n", iter, err);

				if (err < optpar->abserr_step || iter >= itmax_local) {
					/*
					 * we're done! 
					 */
					break;
				}
				GMRFLib_free_problem(*problem);
				*problem = NULL;
			}
			if (!*problem) {
				/*
				 * fail to converge. restart with a reduced step_factor. 
				 */
				FREE_ALL;
				optpar->nr_step_factor /= 10.0;
				optpar->max_iter *= 2;
				if (optpar && optpar->fp) {
					fprintf(stderr, "\n\n%s: Optimisation fail to converge in PART 2.\n\t\t\tRetry with a new optpar->nr_step_factor = %f\n",
						__GMRFLib_FuncName, optpar->nr_step_factor);
				}
				if (optpar->nr_step_factor < 0.00001) {
					return GMRFLib_EOPTNR;
				} else {
					return GMRFLib_init_GMRF_approximation_store__intern(problem, x, b, c, mean, d,
											     loglFunc, loglFunc_arg, fixed_value, graph, Qfunc, Qfunc_arg,
											     constr, optpar, blockupdate_par, store, bb, cc, gaussian_data);
				}
			}
		}
	}

	FREE_ALL;
	GMRFLib_LEAVE_ROUTINE;

	return GMRFLib_SUCCESS;

#undef FREE_ALL
}
char *GMRFLib_ai_tag(int *iz, int len)
{
	/*
	 * return a tag for iz 
	 */
	char *tag = NULL;
	int i, blen = 6;
	size_t len_alloc = len * (blen + 1) + 1;

	tag = Calloc(len_alloc, char);

	tag[0] = '\0';
	for (i = 0; i < len; i++) {
		sprintf(&tag[strlen(tag)], " %2d", iz[i]);
	}
	return tag;
}
int GMRFLib_ai_skip_configurations(map_strd * hash_table, int k, int *iz, int *izz, int *len, int *k_max, int len_length, int nhyper)
{
	/*
	 * mark all configurations >= than 'iz' as to be skipped. 
	 */
	int *iz_local, *izz_local, kk, i, larger;
	char *tag = NULL;

	iz_local = Calloc(nhyper, int);
	izz_local = Calloc(nhyper, int);
	memcpy(iz_local, iz, nhyper * sizeof(int));
	memcpy(izz_local, izz, nhyper * sizeof(int));

	for (kk = k; kk < len_length; kk++) {
		/*
		 * compute the next configuration 
		 */
		for (i = nhyper - 1; i >= 0; i--) {
			if ((izz_local[i] = (izz_local[i] + 1) % len[i])) {
				break;
			}
		}
		for (i = 0; i < nhyper; i++) {
			iz_local[i] = (izz_local[i] <= k_max[i] ? izz_local[i] : k_max[i] - izz_local[i]);
		}

		for (i = 0, larger = 1; i < nhyper && larger; i++) {
			if ((iz[i] < 0 && iz_local[i] > iz[i]) || (iz[i] > 0 && iz_local[i] < iz[i])) {
				larger = 0;
			}
		}
		if (larger) {
			tag = GMRFLib_ai_tag(iz_local, nhyper);
			if (map_strd_ptr(hash_table, tag)) {
				Free(tag);
			} else {
				map_strd_set(hash_table, tag, 1.0);
			}
		}
	}
	Free(iz_local);
	Free(izz_local);

	return GMRFLib_SUCCESS;
}

/*!

  \brief Computes an approximation to the posterior marginal of \f$x_i\f$ \f$\mbox{\boldmath$\theta$}\f$, \f$ \pi(x_i \mid
  \mbox{\boldmath$y$}) \f$ at those indices where compute[i] = #GMRFLib_TRUE . It also computes the integrated likelihood \f$
  \pi(\mbox{\boldmath$y$}) \f$ as the normalizing constant for \f$ \pi(\mbox{\boldmath$\theta$} \mid \mbox{\boldmath$y$}) \f$.

  \param[out] density A pppointer to a \c GMRFLib_density_tp structure. If not \c NULL, then on output it will hold the
  approximated marginal posteriors for those indeces where compute[i] = #GMRFLib_TRUE based on the approximation type specified
  in \a ai_par

  \param[out] gdensity A pppointer to a \c GMRFLib_density_tp structure. If not \c NULL, then on output it will hold the
  approximated marginal posteriors for those indeces where compute[i] = #GMRFLib_TRUE based on the Gaussian mixture
  
  \param[out] density_hyper A pppointer to a \c GMRFLib_density_tp structure. If not \c NULL, then on output it will hold the
  approximated marginal densities for all the hyperparameters.  This feature verified for \c GMRFLib_ai_param_tp::strategy =
  #GMRFLib_AI_INT_STRATEGY_GRID, but can be used at the users own risk for the case \c GMRFLib_ai_param_tp::strategy =
  #GMRFLib_AI_INT_STRATEGY_CCD. Accurate marginals for each hyperparameter requires a more dense and larger grid of
  \f$\theta\f$-values. Spesifically, \c GMRFLib_ai_param_tp::diff_log_dens should be increased, \c
  GMRFLib_ai_param_tp::skip_configurations could be turned off and \c GMRFLib_ai_param_tp::dz could perhaps be descreased.

  \param[out] cpo A ppointer to a \c GMRFLib_ai_cpo_tp stucture. If not \c NULL, then on output it will hold a alloced \c
  GMRFLib_ai_cpo_tp object, containing the values of \f$\log \pi(y_i | y_{-i})\f$ for each i for which the density is computed
  and observed. Note that \c cpo is not created unless either \c density or \c gdensity is requested.  The created \c cpo
  -object, can be free'd, by using the routine \c GMRFLib_ai_cpo_free(). Note that computing the correct \c cpo require that the
  likelihood function is correctly normalised, hence it must be implemented with its normalising constants as well.
  
  \param[out] dic A pointer to a \c GMRFLib_ai_dic_tp - object. If not \c NULL, then on output, it contains the compute DIC
  (deviance information criteria) assuming the latent Gaussian field is in focus.

  \param[out] marginal_likelihood A pointer to a \c GMRFLib_ai_marginal_likelihood_tp - object. If not \c NULL, then on output,
  it contains the log of the integrated likelihood \f$ \pi(\mbox{\boldmath$y$}) \f$.

  \param[out] neffp A pointer to a \c GMRFLib_ai_neffp_tp - object. If not \c NULL, then on output, its contents are computed.
  
  \param[in] compute A length \em n, the non-zero element of the array indentify the index for which the posterior marginals are
  computed.

  \param[in] hyperparam A pointer to a vector of pointers to the hyperparameters of the model. Initially they should contain a
  reasonable guess for the modal configuration.

  \param[in] nhyper The number of  hypeparameters

  \param[in] log_extra A pointer to a user defined function of type \a GMRFLib_ai_log_extra_tp computing all terms in
  <b>(GMRF-35)</b> which are constant with respect to \f$\mbox{\boldmath$x$}\f$ but depend on \f$ \mbox{\boldmath$\theta$}\f$.
  
  \param[in] log_extra_arg A \em void pointer holding the address of a variable or data structure defining additional arguments
  to the function \a logl_extra.
  
  \param[in] x A length \em n array, where \n is the number of nodes in the graph. If \a x = \c NULL then all elements are taken
  to be zero.  If \a fixed_value \f$ \neq \f$ \c NULL, the elements of <em>\b x</em> corresponding to \a fixed_value=1 are the
  fixed values in a conditional simulation. The remaining elements of <em>\b x</em> are not used, and can take arbitrary
  values. If \a fixed_value = \c NULL, all values can be arbitrary.

  \param[in] b If <tt>!NULL</tt>, a length \em n array holding the elements of the vector <em>\b b</em> in <b>(GMRF-34)</b> in
  \ref INLA.  If \c NULL, <em>\b b</em> is assumed to be 0.

  \param[in] c If <tt>!NULL</tt>, a length \em n array holding the elements of the vector \f$ \mbox{\boldmath$\mu$} \f$ in
  <b>(GMRF-34)</b> in \ref INLA.  If \c NULL, <em>\b c</em> is assumed to be 0.

  \param[in] mean If <tt>!NULL</tt>, a length \em n array holding the elements of the vector <em>\b </em> in <b>(GMRF-34)</b> in
  \ref INLA. If \c NULL, \f$ \mbox{\boldmath$\mu$} \f$ is assumed to be 0.

  \param[in] d If <tt>!NULL</tt>, a length \em n array holding the elements of the vector <em>\b d</em> in <b>(GMRF-34)</b> in
  \ref INLA.

  \param[in] loglFunc A function of type \c GMRFLib_logl_tp(), returning the value of the function \f$ f_i(x_i,y_i) \f$ in
  (GMRF-34) in \ref INLA, in many applications equal to the log-likelihood of the problem. If the function \f$ f_i(x_i,y_i) \f$
  depends on unknown hyper-parameters, these are passed throug in \a loglFunc_arg

  \param[in] loglFunc_arg  A \em void pointer holding the address of a variableq or data structure defining additional
  arguments to the function \a loglFunc.

  \param[in] fixed_value The elements of this array, of length \c graph->n, should take the value 0 or 1. The conditional mean
  and covariance matrix of \f$ \{x_i:\mbox{\small\tt fixed\_value}[i]=0\} \f$ given \f$ \{x_i:\mbox{\small\tt
  fixed\_value}[i]=1\} \f$ are computed, and elements of the GMRF <em>\b x</em> corresponding to <tt>fixed_value = 1</tt>, are
  kept fixed at the corresponding values specified in the argument \a x. It is allowed that all values of <tt>fixed_value</tt>
  equals 1.

  \param[in] graph The graph on which the GMRF <em>\b x</em> is defined.

  \param[in] Qfunc A function of type \c GMRFLib_Qfunc_tp(), computing the values of the precision matrix <em>\b Q</em>. If the
  function depends on unknown hyper-parameters, these are passed through \a Qfunc_arg

  \param[in] Qfunc_arg A \em void pointer holding the address of a variable or data structure defining additional arguments to
  the function \a Qfunc.

  \param[in] constr A pointer to a \c GMRFLib_constr_tp -object holding the information on the type of linear constraint, in the
  case of constrained problem.

  \param[in] ai_par Some specifications for the approximation, as a \c GMRFLib_ai_param_tp -object.

  \param[in] ai_store A pointer to a \c GMRFLib_ai_store_tp object where temporary calculations will be stored and retrieved to
  increase the speed.

  \param[in] linear_term_func A pointer to a function returning linear terms. Set to \c NULL if you do not have this.

  \param[in] linear_term_func_arg A pointer to a function-arguments of the function returning linear terms. Set to \c NULL if you do not have this.
 */
int GMRFLib_ai_INLA(GMRFLib_density_tp *** density, GMRFLib_density_tp *** gdensity, GMRFLib_density_tp *** density_hyper,
		    GMRFLib_ai_cpo_tp ** cpo, GMRFLib_ai_dic_tp * dic,
		    GMRFLib_ai_marginal_likelihood_tp * marginal_likelihood, GMRFLib_ai_neffp_tp * neffp,
		    char *compute, double ***hyperparam, int nhyper,
		    GMRFLib_ai_log_extra_tp * log_extra, void *log_extra_arg,
		    double *x, double *b, double *c, double *mean, double *d,
		    GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, char *fixed_value,
		    GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg,
		    GMRFLib_constr_tp * constr, GMRFLib_ai_param_tp * ai_par, GMRFLib_ai_store_tp * ai_store,
		    GMRFLib_linear_term_func_tp * linear_term_func, void *linear_term_func_arg, int nlin, double *Alin, GMRFLib_density_tp *** dlin,
		    GMRFLib_ai_misc_output_tp ** misc_output)
{
	/*
	 * 
	 * compute integrated marginals for the hidden gmrf at those indices where compute[i] = true. the hyperparameters is
	 * specified by a vector of pointers to the hyperparameters (hyperparam) and number of hyperparameters (nhyper). they
	 * contain initially, a reasonalbe guess for the mode.
	 *
	 * output the densities in DENSITY, the density based on the Gaussian mixture in GDENSITY, and the approximate marginals
	 * for the hyperparameters in DENSITY_HYPER.
	 * 
	 */

#define COMPUTE_LINDENS(_store)						\
	if (nlin) {							\
		int _i;							\
		double *improved_mean = Calloc(graph->n, double);	\
		for(_i = 0; _i<graph->n; _i++) {				\
			if (dens[_i] && dens[_i][dens_count]){			\
				improved_mean[_i] = dens[_i][dens_count]->user_mean; \
			} else {					\
				improved_mean[_i] = ai_store->problem->mean_constr[_i]; \
			}						\
		}							\
		lin_dens[dens_count] = GMRFLib_ai_compute_lincomb(nlin, Alin, _store, improved_mean); \
		Free(improved_mean);					\
	}

#define CHECK_HYPER_STORAGE_FORCE(num_) CHECK_HYPER_STORAGE_INTERN(num_, 1)
#define CHECK_HYPER_STORAGE CHECK_HYPER_STORAGE_INTERN(10, 0)
#define CHECK_HYPER_STORAGE_INTERN(num_, force_)			\
	if ((hyper_count >= hyper_len) || (force_)) {			\
		int old_hyper_len = hyper_len;				\
		int ii_;						\
		hyper_len += num_;					\
		hyper_z = Realloc(hyper_z, hyper_len * nhyper, double);	\
		hyper_ldens = Realloc(hyper_ldens, hyper_len, double);	\
		if (nlin > 0) {						\
			lin_dens = Realloc(lin_dens, hyper_len, GMRFLib_density_tp **); \
			for(ii_ = old_hyper_len; ii_ < hyper_len; ii_++) \
				lin_dens[ii_] = NULL;			\
		}							\
	}

#define CHECK_DENS_STORAGE_FORCE(num_)  CHECK_DENS_STORAGE_INTERN(num_, 1)
#define CHECK_DENS_STORAGE CHECK_DENS_STORAGE_INTERN(10, 0)
#define CHECK_DENS_STORAGE_INTERN(num_, force_)				\
	if ((dens_count >= dens_max) || (force_)) {			\
		int ii_, jj_;						\
		int old_dens_max = dens_max;				\
		dens_max += num_;					\
		weights = Realloc(weights, dens_max, double);		\
		if (GMRFLib_ai_INLA_userfunc0) {				\
			userfunc_values = Realloc(userfunc_values, dens_max, double *); \
		}							\
		izs = Realloc(izs, dens_max, double *);			\
		memset(&(izs[old_dens_max]), 0, (num_) * sizeof(double *)); \
		neff = Realloc(neff, dens_max, double);			\
		for (ii_ = 0; ii_ < graph->n; ii_++) {			\
			if (dens[ii_]){					\
				dens[ii_] = Realloc(dens[ii_], dens_max, GMRFLib_density_tp *); \
				for(jj_ = old_dens_max; jj_ < dens_max; jj_++) \
					dens[ii_][jj_] = NULL;		\
			}						\
		}							\
		if (cpo) {						\
			for (ii_ = 0; ii_ < compute_n; ii_++) {		\
				jj_ = compute_idx[ii_];			\
				if (d[jj_] || ai_par->cpo_manual){	\
					cpo_theta[jj_] = Realloc(cpo_theta[jj_], dens_max, double); \
					pit_theta[jj_] = Realloc(pit_theta[jj_], dens_max, double); \
					failure_theta[jj_] = Realloc(failure_theta[jj_], dens_max, double); \
				}					\
			}						\
		}							\
		if (dic) {						\
			for (ii_ = 0; ii_ < compute_n; ii_++) {		\
				jj_ = compute_idx[ii_];			\
				if (d[jj_]){				\
					deviance_theta[jj_] = Realloc(deviance_theta[jj_], dens_max, double); \
				}					\
			}						\
		}							\
	}

#define SET_THETA_MODE							\
	if (theta_mode) {						\
		int i_, j_;						\
		for(j_=0; j_ < tmax; j_++) {				\
			for(i_ = 0; i_ < nhyper; i_++){			\
				hyperparam[i_][j_][0] = theta_mode[i_]; \
			}						\
		}							\
	}


#define COMPUTE_CPO_AND_DIC						\
	if (d[ii] || ai_par->cpo_manual) {				\
		if (cpo || ai_par->cpo_manual) {			\
			failure_theta[ii][dens_count] = GMRFLib_ai_cpopit_integrate(&cpo_theta[ii][dens_count], \
										    &pit_theta[ii][dens_count], ii, cpodens, \
										    loglFunc, loglFunc_arg, xx_mode); \
		}							\
		if (dic) {						\
			deviance_theta[ii][dens_count] = GMRFLib_ai_dic_integrate(ii, dens[ii][dens_count], \
										  loglFunc, loglFunc_arg, xx_mode); \
		}							\
	}

#define COMPUTE_CPO_AND_DIC_LOCAL					\
	if (d[ii] || ai_par->cpo_manual) {				\
		if (cpo || ai_par->cpo_manual) {			\
			failure_theta_local[ii] +=			\
				GMRFLib_ai_cpopit_integrate(&cpo_theta_local[ii], &pit_theta_local[ii],	\
							    ii, cpodens, loglFunc, loglFunc_arg, xx_mode); \
		}							\
		if (0) printf("COMPUTE cpo for %d %g %g\n", ii,  cpo_theta_local[ii], pit_theta_local[ii]); \
		if (dic) {						\
			deviance_theta_local[ii] =			\
				GMRFLib_ai_dic_integrate(ii, dens_local[ii], loglFunc, loglFunc_arg, xx_mode); \
		}							\
	}


/* 
   since all threads compute the same quantity, this is it well defined
*/
#define COMPUTE_NEFF					\
	if (run_with_omp ) {				\
		neff[dens_count] = ai_store_id->neff;	\
	} else {					\
		neff[dens_count] = ai_store->neff;	\
	}

#define COMPUTE_NEFF_LOCAL neff_local = ai_store_id->neff

#define ADD_LINEAR_TERM							\
	if (linear_term_func) {						\
		double mu_add = linear_term_func(ii, linear_term_func_arg); \
		dens[ii][dens_count]->user_mean += mu_add;		\
		dens[ii][dens_count]->std_mean += mu_add;		\
		assert(dens[ii][dens_count]->spline_P == NULL);		\
	}

#define ADD_LINEAR_TERM_LOCAL						\
	if (linear_term_func) {						\
		double mu_add = linear_term_func(ii, linear_term_func_arg); \
		dens_local[ii]->user_mean += mu_add;			\
		dens_local[ii]->std_mean += mu_add;			\
		assert(dens_local[ii]->spline_P == NULL);		\
	}

#define COMPUTE       COMPUTE_NEFF;       COMPUTE_CPO_AND_DIC;       ADD_LINEAR_TERM
#define COMPUTE_LOCAL COMPUTE_NEFF_LOCAL; COMPUTE_CPO_AND_DIC_LOCAL; ADD_LINEAR_TERM_LOCAL

	int i, j, k, *k_max = NULL, *k_min = NULL, *k_maxx = NULL, *k_minn = NULL, ierr, *iz = NULL, *izz = NULL, *len =
	    NULL, *iz_axes = NULL, skip, dir, len_length, free_ai_par = 0, config_count = 0, free_compute = 0, dens_count =
	    0, dens_max, hyper_len = 0, hyper_count = 0, *compute_idx = NULL, compute_n, tmax, run_with_omp, need_Qinv = 1;
	double *hessian = NULL, *theta = NULL, *theta_mode = NULL, *x_mode = NULL, log_dens_mode, log_dens, *z = NULL, **izs =
	    NULL, *stdev_corr_pos = NULL, *stdev_corr_neg = NULL, f, w, w_origo, tref, tu, *weights = NULL, *adj_weights =
	    NULL, *hyper_z = NULL, *hyper_ldens = NULL, **userfunc_values = NULL, *inverse_hessian = NULL,
	    **cpo_theta = NULL, **pit_theta = NULL, **deviance_theta = NULL, *neff = NULL, **failure_theta = NULL;
	char *tag = NULL;
	gsl_matrix *H = NULL, *eigen_vectors = NULL;
	gsl_eigen_symmv_workspace *work = NULL;
	gsl_vector *eigen_values = NULL;
	map_strd hash_table;
	GMRFLib_density_tp ***dens = NULL;
	GMRFLib_density_tp ***lin_dens = NULL;
	GMRFLib_ai_store_tp **ais = NULL;

	if (fixed_value) {
		FIXME("\n\n\n\nGMRFLib_INLA() do not longer work with FIXED_VALUE; please write a wrapper.\n");
		abort();
	}

	tmax = omp_get_max_threads();
	run_with_omp = (omp_get_max_threads() > 1 ? 1 : 0);

	if (!ai_par) {
		GMRFLib_default_ai_param(&ai_par);
		free_ai_par = 1;
	}
	/*
	 * otherwise, it might go very wrong below 
	 */
	GMRFLib_ASSERT(ai_par && (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_GRID ||
				  ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_EMPIRICAL_BAYES ||
				  ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_CCD), GMRFLib_EPARAMETER);

	GMRFLib_ENTER_ROUTINE;

	if (misc_output) {
		*misc_output = Calloc(1, GMRFLib_ai_misc_output_tp);
	}

	/*
	 * this has to be true, I think... 
	 */
	if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_CCD) {
		ai_par->dz = 1.0;
	}

	if (!compute) {
		free_compute = 1;
		compute = Calloc(graph->n, char);
	}

	nhyper = IMAX(0, nhyper);
	dens_max = IMIN(20, IMAX((int) gsl_pow_int(2.0, nhyper), 1));
	dens = Calloc(graph->n, GMRFLib_density_tp **);
	weights = Calloc(dens_max, double);
	izs = Calloc(dens_max, double *);
	neff = Calloc(dens_max, double);

	map_strd_init_hint(&hash_table, dens_max);
	hash_table.alwaysdefault = 0;

	if ((density || gdensity) && cpo) {
		(*cpo) = Calloc(1, GMRFLib_ai_cpo_tp);
		(*cpo)->n = graph->n;
		(*cpo)->value = Calloc(graph->n, double *);
		(*cpo)->pit_value = Calloc(graph->n, double *);
		(*cpo)->failure = Calloc(graph->n, double *);
	} else {
		cpo = NULL;
	}
	if (GMRFLib_ai_INLA_userfunc0) {
		userfunc_values = Calloc(dens_max, double *);
	}
	/*
	 * make a list of those idxs we will compute 
	 */
	compute_idx = Calloc(graph->n, int);

	compute_n = 0;
	for (i = 0; i < graph->n; i++) {
		if (compute[i]) {
			compute_idx[compute_n++] = i;
		}
	}

	/*
	 * only one of the marginal are computed with cpo_manual is TRUE. The code depends on the this assumption I think. 
	 */
	if (ai_par->cpo_manual) {
		GMRFLib_ASSERT(compute_n == 1, GMRFLib_ESNH);
		GMRFLib_ASSERT(d[compute_idx[0]] == 0.0, GMRFLib_ESNH);
		if (dic) {				       /* meaningless to compute the DIC in this case */
			dic = NULL;
		}
	}

	need_Qinv = ((compute_n || ai_par->compute_nparam_eff) ? 1 : 0);

	for (i = 0; i < compute_n; i++) {
		j = compute_idx[i];
		dens[j] = Calloc(dens_max, GMRFLib_density_tp *);	/* storage for the marginals */
	}

	if (cpo) {
		cpo_theta = Calloc(graph->n, double *);	       /* cpo-value conditioned on theta */
		pit_theta = Calloc(graph->n, double *);	       /* pit-value conditioned on theta */
		failure_theta = Calloc(graph->n, double *);    /* failure indicator on theta */
		for (i = 0; i < compute_n; i++) {
			j = compute_idx[i];
			if (d[j] || ai_par->cpo_manual) {
				cpo_theta[j] = Calloc(dens_max, double);
				pit_theta[j] = Calloc(dens_max, double);
				failure_theta[j] = Calloc(dens_max, double);
			}
		}
	}
	if (dic) {
		deviance_theta = Calloc(graph->n, double *);   /* mean of deviance conditioned on theta */
		for (i = 0; i < compute_n; i++) {
			j = compute_idx[i];
			if (d[j]) {
				deviance_theta[j] = Calloc(dens_max, double);
			}
		}
	}

	if (nhyper > 0) {
		/*
		 * the first step is to locate the mode of \pi(\theta | y). here we use the domin-optimiser routine.  NOTE that this
		 * '_setup' ensure that ai_store is changed for each call to _domin_f. this is a bit dirty programming, but there is no
		 * good way to get around it for the moment.
		 */

		GMRFLib_domin_setup(hyperparam, nhyper, log_extra, log_extra_arg, compute, x, b, c, mean, d, loglFunc, loglFunc_arg,
				    fixed_value, graph, Qfunc, Qfunc_arg, constr, ai_par, ai_store);
		domin_seteps_(&(ai_par->domin_epsx), &(ai_par->domin_epsf), &(ai_par->domin_epsg));
		/*
		 * the optimizer runs most smoothly when #threads is about nhyper+1, which is the number of `natural' threads for
		 * computing the gradient.
		 */
		int tmax_local;


		theta = Calloc(nhyper, double);		       /* theta is the hyperparameters */
		theta_mode = Calloc(nhyper, double);
		x_mode = Calloc(graph->n, double);
		z = Calloc(nhyper, double);
		tmax_local = IMIN(nhyper + 1, tmax);
		if (tmax_local < tmax && !ai_par->huge) {
			/*
			 * only if there is a change we do not have a HUGE model
			 */
			omp_set_num_threads(tmax_local);
		}

		/*
		 * if not set to be known, then optimise 
		 */
		if (!(ai_par->mode_known)) {
			
			if (ai_par->fp_log) {
				fprintf(ai_par->fp_log, "Optimise using %s\n", GMRFLib_AI_OPTIMISER_NAME(ai_par->optimiser));
			}

			switch (ai_par->optimiser) {
			case GMRFLib_AI_OPTIMISER_DOMIN:
				domin_();		       /* this is the optimizer */
				if (ai_par->restart) {
					for (k = 0; k < IMAX(0, ai_par->restart); k++)
						domin_();      /* restart */
				}
				domin_get_results_(theta_mode, &log_dens_mode, &ierr);
				break;

			case GMRFLib_AI_OPTIMISER_GSL:
			case GMRFLib_AI_OPTIMISER_DEFAULT:
				GMRFLib_gsl_optimize(ai_par);
				if (ai_par->restart) {
					for (k = 0; k < IMAX(0, ai_par->restart); k++)
						GMRFLib_gsl_optimize(ai_par);	/* restart */
				}
				GMRFLib_gsl_get_results(theta_mode, &log_dens_mode);
				break;

			default:
				GMRFLib_ASSERT((ai_par->optimiser == GMRFLib_AI_OPTIMISER_DOMIN) ||
					       (ai_par->optimiser == GMRFLib_AI_OPTIMISER_GSL) ||
					       (ai_par->optimiser == GMRFLib_AI_OPTIMISER_DEFAULT), GMRFLib_EPARAMETER);
				break;
			}

			if (ai_par->fp_log) {
				fprintf(ai_par->fp_log, "Number of function evaluations = %1d\n", GMRFLib_domin_get_f_count());
			}
		} else {
			/*
			 * use the initial values only 
			 */
			for (i = 0; i < nhyper; i++) {
				theta_mode[i] = hyperparam[i][0][0];
			}
			if (ai_par->fp_log) {
				fprintf(ai_par->fp_log, "Using known modal configuration = [");
				for (i = 0; i < nhyper; i++) {
					fprintf(ai_par->fp_log, " %6.3f", theta_mode[i]);
				}
				fprintf(ai_par->fp_log, "]\n");
			}
			GMRFLib_domin_f(theta_mode, &log_dens_mode, &ierr);
			log_dens_mode *= -1.0;
			if (ai_par->fp_log) {
				fprintf(ai_par->fp_log, "Compute mode: %10.4f\n", log_dens_mode);
			}
		}


		if (tmax_local < tmax && !ai_par->huge) {
			/*
			 * set it back 
			 */
			omp_set_num_threads(tmax);
		}

		SET_THETA_MODE;

		if (ai_par->fp_log) {
			fprintf(ai_par->fp_log, "Compute the Hessian using %s differences and step_size[%g]. Matrix-type [%s]\n",
				(ai_par->hessian_forward_finite_difference ? "forward" : "central"),
				ai_par->hessian_finite_difference_step_len, (ai_par->hessian_force_diagonal ? "diagonal" : "dense"));
		}


		/*
		 * The parameters for the adaptive hessian estimation is set in ai_par (hence G.ai_par in domin-interface.c).
		 */
		int hess_count = 0;
		double log_dens_mode_save = log_dens_mode;
		int count_log_dens_mode_fail = 0;
		int count_log_dens_mode_fail_max = ai_par->adaptive_hessian_max_trials;

		hessian = Calloc(ISQR(nhyper), double);
		while (GMRFLib_domin_estimate_hessian(hessian, theta_mode, &log_dens_mode, hess_count) != GMRFLib_SUCCESS) {
			if (!hess_count) {
				if (ai_par->fp_log)
					fprintf(ai_par->fp_log, "Mode not sufficient accurate; switch to a stupid local search strategy.\n");
			}
			hess_count++;

			count_log_dens_mode_fail += (log_dens_mode_save > log_dens_mode ? 1 : 0);
			if (log_dens_mode_save > log_dens_mode && count_log_dens_mode_fail > count_log_dens_mode_fail_max) {
				if (ai_par->fp_log) {
					fprintf(stderr, "\n\n*** Mode is not accurate yet but we have reached the rounding error level. Break.\n\n");
				}
				break;
			}
			// printf("%.12g %.12g\n", log_dens_mode_save, log_dens_mode);
			log_dens_mode_save = log_dens_mode;

			if (GMRFLib_request_optimiser_to_stop) {
				fprintf(stderr, "\n\n*** Optimiser requested to stop; stop local search..\n");
				break;
			}
			if (hess_count >= ai_par->adaptive_hessian_max_trials) {
				fprintf(stderr, "\n\n*** Mode not found using the stupid local search strategy; I give up.\n");
				fprintf(stderr, "*** Try to modify the initial values.\n");
				GMRFLib_ASSERT(hess_count < ai_par->adaptive_hessian_max_trials, GMRFLib_EMISC);
			}
		}

		/*
		 * do this again to get the ai_store set correctly.
		 */
		SET_THETA_MODE;
		if (hess_count) {
			if (ai_par->fp_log) {
				fprintf(ai_par->fp_log, "Restart optimisation at the best mode found so far and reesteimate the Hessian\n");
			}
			switch (ai_par->optimiser) {
			case GMRFLib_AI_OPTIMISER_DOMIN:
				domin_();		       /* this is the optimizer */
				domin_get_results_(theta_mode, &log_dens_mode, &ierr);
				break;

			case GMRFLib_AI_OPTIMISER_GSL:
			case GMRFLib_AI_OPTIMISER_DEFAULT:
				GMRFLib_gsl_optimize(ai_par);
				GMRFLib_gsl_get_results(theta_mode, &log_dens_mode);
				break;

			default:
				GMRFLib_ASSERT((ai_par->optimiser == GMRFLib_AI_OPTIMISER_DOMIN) ||
					       (ai_par->optimiser == GMRFLib_AI_OPTIMISER_GSL) ||
					       (ai_par->optimiser == GMRFLib_AI_OPTIMISER_DEFAULT), GMRFLib_EPARAMETER);
				break;
			}
			GMRFLib_domin_estimate_hessian(hessian, theta_mode, &log_dens_mode, 0);
			SET_THETA_MODE;
		}

		if (ai_par->fp_log) {
			for (i = 0; i < nhyper; i++) {
				for (j = 0; j < nhyper; j++) {
					fprintf(ai_par->fp_log, " %12.6f", hessian[i + j * nhyper]);
				}
				fprintf(ai_par->fp_log, "\n");
			}
		}

		H = gsl_matrix_calloc((size_t) nhyper, (size_t) nhyper);
		for (i = 0; i < nhyper; i++) {
			for (j = 0; j < nhyper; j++) {
				gsl_matrix_set(H, (size_t) i, (size_t) j, hessian[i + nhyper * j]);
			}
		}
		work = gsl_eigen_symmv_alloc((size_t) nhyper);
		eigen_vectors = gsl_matrix_calloc((size_t) nhyper, (size_t) nhyper);
		eigen_values = gsl_vector_calloc((size_t) nhyper);
		gsl_eigen_symmv(H, eigen_values, eigen_vectors, work);
		gsl_eigen_symmv_free(work);

		if (ai_par->fp_log) {
			fprintf(ai_par->fp_log, "Eigenvectors of the Hessian\n");
			GMRFLib_gsl_matrix_fprintf(ai_par->fp_log, eigen_vectors, "\t%f");
			fprintf(ai_par->fp_log, "Eigenvalues of the Hessian\n");
			gsl_vector_fprintf(ai_par->fp_log, eigen_values, "\t%f");
		}

		/*
		 * check that the hessian is positive definite 
		 */

		// OLD CODE:
		// GMRFLib_ASSERT(gsl_vector_get(eigen_values, (unsigned int) i) > 0.0, GMRFLib_EPOSDEF);

		// NEW:
		double min_pos_eigenvalue = DBL_MAX;
		for (i = 0; i < nhyper; i++) {
			double eigv = gsl_vector_get(eigen_values, (unsigned int) i);

			if (eigv > 0.0 && eigv < min_pos_eigenvalue)
				min_pos_eigenvalue = eigv;
		}
		if (min_pos_eigenvalue == DBL_MAX) {
			min_pos_eigenvalue = 1.0;	       /* if all are negative */
		} else {
			min_pos_eigenvalue /= 100.0;	       /* JUST A CHOICE */
		}
		int a_change = 0;

		for (i = 0; i < nhyper; i++) {
			double eigv = gsl_vector_get(eigen_values, (unsigned int) i);

			if (eigv < 0.0) {
				fprintf(stderr, "\n");
				fprintf(stderr, "\t*** WARNING *** Eigenvalue %1d of the Hessian is %.6g < 0\n", i, eigv);
				fprintf(stderr, "\t*** WARNING *** Set this eigenvalue to %.6g\n", min_pos_eigenvalue);
				fprintf(stderr, "\t*** WARNING *** This might have consequence for the accurancy of\n");
				fprintf(stderr, "\t*** WARNING *** the approximations; please check!\n");
				fprintf(stderr, "\t*** WARNING *** R-inla: Use option inla(..., control.inla = list(h = h.value), ...) \n");
				fprintf(stderr, "\t*** WARNING *** R-inla: to chose a different  `h.value'.\n");
				fprintf(stderr, "\n");

				gsl_vector_set(eigen_values, (unsigned int) i, min_pos_eigenvalue);
				a_change++;
			}
		}

		if (a_change) {
			/*
			 * rebuild the Hessian using the new eigenvalues. I should have used matrix-multiplication routines, but I had this code already from
			 * af-program.c ;-) In any case, the matrix is small... 
			 */

			for (i = 0; i < nhyper; i++) {
				for (j = i; j < nhyper; j++) {
					double sum = 0.0;
					for (k = 0; k < nhyper; k++) {
						sum += gsl_matrix_get(eigen_vectors, i, k) * gsl_matrix_get(eigen_vectors, j, k)
						    * gsl_vector_get(eigen_values, k);
					}
					hessian[i + j * nhyper] = hessian[j + i * nhyper] = sum;
				}
			}
		}


		/*
		 * compute the inverse hessian, for scaling purposes 
		 */
		inverse_hessian = Calloc(ISQR(nhyper), double);
		memcpy(inverse_hessian, hessian, ISQR(nhyper) * sizeof(double));
		GMRFLib_comp_posdef_inverse(inverse_hessian, nhyper);

		if (misc_output) {
			(*misc_output)->nhyper = nhyper;
			(*misc_output)->cov_m = Calloc(ISQR(nhyper), double);
			memcpy((*misc_output)->cov_m, inverse_hessian, ISQR(nhyper) * sizeof(double));
		}

		if (ai_par->fp_log) {
			/*
			 * print the stdev/correlation matrix: stdevs on the diagonal and the correlations on the off-diagonal.
			 */
			int ii, jj;
			double val;

			fprintf(ai_par->fp_log, "StDev/Correlation matrix (scaled inverse Hessian)\n");
			for (ii = 0; ii < nhyper; ii++) {
				for (jj = 0; jj < nhyper; jj++) {
					if (jj >= ii) {
						if (ii == jj) {
							val = sqrt(inverse_hessian[ii + jj * nhyper]);
						} else {
							val = inverse_hessian[ii + jj * nhyper] /
							    sqrt(inverse_hessian[ii + ii * nhyper] * inverse_hessian[jj + jj * nhyper]);
						}
						fprintf(ai_par->fp_log, " %12.6f", val);
					} else {
						fprintf(ai_par->fp_log, " %12s", "");
					}
				}
				fprintf(ai_par->fp_log, "\n");
			}
		}

		/*
		 * setup space for storage; used for openmp 
		 */
		ais = Calloc(tmax, GMRFLib_ai_store_tp *);

		/*
		 * search the space. first, look at main directions and stop if the density differ more than dlog_dens from the value
		 * at the mode, log_dens_mode. outside the main directions, only add the point if the corresponding values at the main
		 * directions is in. 
		 *
		 * include a check of user-arguments. if a node a fixed, we cannot compute its marginal density...
		 */
		if (fixed_value) {
			for (i = 0; i < graph->n; i++) {
				if (compute[i] && fixed_value[i]) {
					compute[i] = 0;
				}
			}
		}

		iz = Calloc(nhyper, int);
		memset(iz, 0, nhyper * sizeof(int));

		hyper_len = dens_max;
		hyper_count = 0;
		hyper_z = Calloc(hyper_len * nhyper, double);
		hyper_ldens = Calloc(hyper_len, double);
		if (nlin > 0) {
			lin_dens = Calloc(hyper_len, GMRFLib_density_tp **);
		} else {
			nlin = 0;
		}

		/*
		 * compute the corrected scalings/stdevs, if required. 
		 */
		if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_CCD
		    || (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_GRID && density_hyper && ai_par->interpolator == GMRFLib_AI_INTERPOLATOR_CCD)) {

			stdev_corr_pos = Calloc(nhyper, double);
			stdev_corr_neg = Calloc(nhyper, double);

			/*
			 * two versions: 1. a nhyper loop, 2. a 2*nhyper loop. 
			 */
			if (tmax > nhyper) {
#pragma omp parallel for private(k) schedule(static)
				for (k = 0; k < 2 * nhyper; k++) {
					double f0, *zz = NULL, *ttheta = NULL, llog_dens;
					int kk, opt;
					GMRFLib_ai_store_tp *s = NULL;

					if (k < nhyper) {
						kk = k;
						opt = 0;
					} else {
						kk = k - nhyper;
						opt = 1;
					}
					zz = Calloc(nhyper, double);
					ttheta = Calloc(nhyper, double);
					memset(zz, 0, nhyper * sizeof(double));
					GMRFLib_thread_id = omp_get_thread_num();

					if (omp_in_parallel()) {
						if (!ais[GMRFLib_thread_id]) {
							ais[GMRFLib_thread_id] = GMRFLib_duplicate_ai_store(ai_store);
						}
						s = ais[GMRFLib_thread_id];
					} else {
						s = ai_store;  /* the common one */
					}

					if (opt == 0) {
						zz[kk] = 2.0;
						GMRFLib_ai_z2theta(ttheta, nhyper, theta_mode, zz, eigen_values, eigen_vectors);
						GMRFLib_domin_f_intern(ttheta, &llog_dens, &ierr, s);
						llog_dens *= -1.0;
						f0 = log_dens_mode - llog_dens;
						stdev_corr_pos[kk] = (f0 > 0.0 ? sqrt(2.0 / f0) : 1.0);
					} else {
						zz[kk] = -2.0;
						GMRFLib_ai_z2theta(ttheta, nhyper, theta_mode, zz, eigen_values, eigen_vectors);
						GMRFLib_domin_f_intern(ttheta, &llog_dens, &ierr, s);
						llog_dens *= -1.0;
						f0 = log_dens_mode - llog_dens;
						stdev_corr_neg[kk] = (f0 > 0.0 ? sqrt(2.0 / f0) : 1.0);
					}

					Free(zz);
					Free(ttheta);
				}
			} else {
#pragma omp parallel for private(k) schedule(static)
				for (k = 0; k < nhyper; k++) {
					double f0, *zz = NULL, *ttheta = NULL, llog_dens;
					GMRFLib_ai_store_tp *s = NULL;

					zz = Calloc(nhyper, double);
					ttheta = Calloc(nhyper, double);
					memset(zz, 0, nhyper * sizeof(double));
					GMRFLib_thread_id = omp_get_thread_num();

					if (omp_in_parallel()) {
						if (!ais[GMRFLib_thread_id]) {
							ais[GMRFLib_thread_id] = GMRFLib_duplicate_ai_store(ai_store);
						}
						s = ais[GMRFLib_thread_id];
					} else {
						s = ai_store;  /* the common one */
					}

					zz[k] = 2.0;
					GMRFLib_ai_z2theta(ttheta, nhyper, theta_mode, zz, eigen_values, eigen_vectors);
					GMRFLib_domin_f_intern(ttheta, &llog_dens, &ierr, s);
					llog_dens *= -1.0;
					f0 = log_dens_mode - llog_dens;
					stdev_corr_pos[k] = (f0 > 0.0 ? sqrt(2.0 / f0) : 1.0);

					zz[k] = -2.0;
					GMRFLib_ai_z2theta(ttheta, nhyper, theta_mode, zz, eigen_values, eigen_vectors);
					GMRFLib_domin_f_intern(ttheta, &llog_dens, &ierr, s);
					llog_dens *= -1.0;
					f0 = log_dens_mode - llog_dens;
					stdev_corr_neg[k] = (f0 > 0.0 ? sqrt(2.0 / f0) : 1.0);

					Free(zz);
					Free(ttheta);
				}
			}

			for (k = 0; k < nhyper; k++) {
				if (ai_par->fp_log) {
					fprintf(ai_par->fp_log,
						"Compute corrected stdev for theta[%1d]: negative %f  positive %f\n", k, stdev_corr_neg[k], stdev_corr_pos[k]);
				}
			}
		}

		if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_EMPIRICAL_BAYES) {
			if (need_Qinv) {
				GMRFLib_ai_add_Qinv_to_ai_store(ai_store);	/* add Qinv if required */
				GMRFLib_ai_si(ai_par, 0.0, theta_mode, nhyper, graph, ai_store);
			}
			ai_store->neff = GMRFLib_AI_STORE_NEFF_NOT_COMPUTED;

			if (run_with_omp) {
#pragma omp parallel
				{
					GMRFLib_ai_store_tp *ai_store_id = GMRFLib_duplicate_ai_store(ai_store);
#pragma omp for private(i) schedule(static) nowait
					for (i = 0; i < compute_n; i++) {
						int ii = compute_idx[i];
						GMRFLib_density_tp *cpodens = NULL;

						GMRFLib_thread_id = 0;

						GMRFLib_ai_marginal_hidden(&dens[ii][dens_count],
									   (cpo && (d[ii] || ai_par->cpo_manual) ? &cpodens : NULL), ii, x, b, c, mean, d,
									   loglFunc, loglFunc_arg, fixed_value, graph, Qfunc,
									   Qfunc_arg, constr, ai_par, ai_store_id);
						double *xx_mode = ai_store_id->mode;

						COMPUTE;
						GMRFLib_free_density(cpodens);
					}
					GMRFLib_free_ai_store(ai_store_id);
				}
			} else {
				GMRFLib_ai_store_tp *ai_store_id = NULL;
				for (i = 0; i < compute_n; i++) {
					int ii = compute_idx[i];
					GMRFLib_density_tp *cpodens = NULL;

					GMRFLib_thread_id = 0;

					GMRFLib_ai_marginal_hidden(&dens[ii][dens_count], (cpo && (d[ii] || ai_par->cpo_manual) ? &cpodens : NULL),
								   ii, x, b, c, mean, d,
								   loglFunc, loglFunc_arg, fixed_value, graph, Qfunc, Qfunc_arg, constr, ai_par, ai_store);
					double *xx_mode = ai_store->mode;
					COMPUTE;
					GMRFLib_free_density(cpodens);
				}
			}
			if (GMRFLib_ai_INLA_userfunc0) {
				userfunc_values[dens_count] = GMRFLib_ai_INLA_userfunc0(ai_store->problem, theta, nhyper);
			}
			COMPUTE_LINDENS(ai_store);

			izs[dens_count] = Calloc(nhyper, double);
			for (i = 0; i < nhyper; i++) {
				izs[dens_count][i] = 0;
			}
			weights[dens_count] = 0.0;
			dens_count++;

			/*
			 * END OF GMRFLib_AI_INT_STRATEGY_EMPIRICAL_BAYES 
			 */
		} else if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_CCD) {
			/*
			 * use points from the ccd-design to do the integration
			 */

			GMRFLib_design_tp *design = NULL;
			GMRFLib_get_design(&design, nhyper);

			f = DMAX(ai_par->f0, 1.0) * sqrt((double) nhyper);
			w = 1.0 / ((design->nexperiments - 1.0) * (1.0 + exp(-0.5 * SQR(f)) * (SQR(f) / nhyper - 1.0)));
			w_origo = 1.0 - (design->nexperiments - 1.0) * w;
			if (run_with_omp) {
				/*
				 * this new code parallelise over each configuration, and not within each configuration. 
				 */
				GMRFLib_ASSERT(dens_count == 0, GMRFLib_ESNH);
				GMRFLib_ASSERT(hyper_count == 0, GMRFLib_ESNH);

				/*
				 * ensure enough storage 
				 */
				CHECK_DENS_STORAGE_FORCE(design->nexperiments);
				CHECK_HYPER_STORAGE_FORCE(design->nexperiments);

#pragma omp parallel for private(k, i, log_dens, dens_count, hyper_count, tref, tu, ierr) schedule(static)
				for (k = 0; k < design->nexperiments; k++) {

					double *z_local, *theta_local, log_dens_orig;
					GMRFLib_ai_store_tp *ai_store_id = NULL;

					dens_count = k;
					hyper_count = k;
					GMRFLib_thread_id = omp_get_thread_num();

					if (omp_in_parallel()) {
						if (!ais[GMRFLib_thread_id]) {
							ais[GMRFLib_thread_id] = GMRFLib_duplicate_ai_store(ai_store);
						}
						ai_store_id = ais[GMRFLib_thread_id];
					} else {
						ai_store_id = ai_store;	/* the common one */
					}

					z_local = Calloc(nhyper, double);
					theta_local = Calloc(nhyper, double);

					for (i = 0; i < nhyper; i++) {
						z_local[i] = f * design->experiment[k][i]
						    * (design->experiment[k][i] > 0.0 ? stdev_corr_pos[i] : stdev_corr_neg[i]);
					}
					GMRFLib_ai_z2theta(theta_local, nhyper, theta_mode, z_local, eigen_values, eigen_vectors);
					GMRFLib_domin_f_intern(theta_local, &log_dens, &ierr, ai_store_id);
					log_dens *= -1.0;
					log_dens_orig = log_dens;

					/*
					 * correct the log_dens due to the integration weights which is special for the CCD integration:
					 * double the weights for the points not in the center
					 */
					if (nhyper > 1) {
						/*
						 * the weight formula is only valid for nhyper > 1. 
						 */
						int origo = 1;

						for (i = 0; i < nhyper; i++) {
							origo = (origo && ISZERO(z_local[i]));
						}
						if (origo) {
							if (ISZERO(w_origo)) {
								log_dens += -DBL_MAX;
							} else {
								log_dens += log(w_origo);
							}
						} else {
							log_dens += log(w);
						}
					}

					/*
					 * register the density for the marginal of the hyperparameters computations. first check space. 
					 */
					for (i = 0; i < nhyper; i++) {
						hyper_z[hyper_count * nhyper + i] = z_local[i];
					}
					hyper_ldens[hyper_count] = log_dens - log_dens_mode;

					/*
					 * compute the marginals for this point. check storage 
					 */
					if (nhyper > 0) {
						weights[dens_count] = log_dens;
					} else {
						weights[dens_count] = 0.0;
					}
					izs[dens_count] = Calloc(nhyper, double);

					for (i = 0; i < nhyper; i++) {
						izs[dens_count][i] = z_local[i];
					}

					tref = GMRFLib_cpu();

					if (need_Qinv) {
						GMRFLib_ai_add_Qinv_to_ai_store(ai_store_id);	/* add Qinv */
						GMRFLib_ai_si(ai_par, log_dens, theta_local, nhyper, graph, ai_store_id);
					}
					ai_store_id->neff = GMRFLib_AI_STORE_NEFF_NOT_COMPUTED;
					for (i = 0; i < compute_n; i++) {
						int ii = compute_idx[i];
						GMRFLib_density_tp *cpodens = NULL;

						GMRFLib_ai_marginal_hidden(&dens[ii][dens_count],
									   (cpo && (d[ii] || ai_par->cpo_manual) ? &cpodens : NULL),
									   ii, x, b, c, mean, d,
									   loglFunc, loglFunc_arg, fixed_value,
									   graph, Qfunc, Qfunc_arg, constr, ai_par, ai_store_id);
						double *xx_mode = ai_store_id->mode;
						COMPUTE;
						GMRFLib_free_density(cpodens);
					}
					if (GMRFLib_ai_INLA_userfunc0) {
						userfunc_values[dens_count] = GMRFLib_ai_INLA_userfunc0(ai_store_id->problem, theta_local, nhyper);
					}
					COMPUTE_LINDENS(ai_store_id);
					tu = GMRFLib_cpu() - tref;
					if (ai_par->fp_log) {
#pragma omp critical
						{
							fprintf(ai_par->fp_log, "\tconfig %2d/%1d=[", config_count++, design->nexperiments);
							for (i = 0; i < nhyper; i++) {
								fprintf(ai_par->fp_log, " %5.2f", z_local[i]);
							}
							/*
							 * we need to use the log_dens_orig as the other one is also included the integration weights. 
							 */
							fprintf(ai_par->fp_log, "] log(rel.dens)=%5.2f, [%1d] accept, compute,", log_dens_orig - log_dens_mode,
								omp_get_thread_num());
							fprintf(ai_par->fp_log, " %.2fs\n", tu);
						}
					}

					Free(z_local);
					Free(theta_local);
				}

				/*
				 * set the values back 
				 */
				dens_count = design->nexperiments;
				hyper_count = design->nexperiments;

			} else {
				/*
				 * old code, which parallelise within each configuration only. 
				 */
				for (k = 0; k < design->nexperiments; k++) {
					for (i = 0; i < nhyper; i++) {
						z[i] = f * design->experiment[k][i]
						    * (design->experiment[k][i] > 0.0 ? stdev_corr_pos[i] : stdev_corr_neg[i]);
					}
					GMRFLib_ai_z2theta(theta, nhyper, theta_mode, z, eigen_values, eigen_vectors);
					GMRFLib_domin_f(theta, &log_dens, &ierr);
					log_dens *= -1.0;

					if (ai_par->fp_log) {
						fprintf(ai_par->fp_log, "\tconfig %2d=[", config_count++);
						for (i = 0; i < nhyper; i++) {
							fprintf(ai_par->fp_log, " %5.2f", z[i]);
						}
						fprintf(ai_par->fp_log, "] log(rel.dens)=%5.2f, [%1d] accept, compute,", log_dens - log_dens_mode,
							omp_get_thread_num());
					}

					/*
					 * correct the log_dens due to the integration weights which is special for the CCD integration:
					 * double the weights for the points not in the center
					 */
					if (nhyper > 1) {
						/*
						 * the weight formula is only valid for nhyper > 1. 
						 */
						int origo = 1;

						for (i = 0; i < nhyper; i++) {
							origo = (origo && ISZERO(z[i]));
						}
						if (origo) {
							if (ISZERO(w_origo)) {
								log_dens += -DBL_MAX;
							} else {
								log_dens += log(w_origo);
							}
						} else {
							log_dens += log(w);
						}
					}

					/*
					 * register the density for the marginal of the hyperparameters computations. first check space. 
					 */
					CHECK_HYPER_STORAGE;
					for (i = 0; i < nhyper; i++) {
						hyper_z[hyper_count * nhyper + i] = z[i];
					}
					hyper_ldens[hyper_count] = log_dens - log_dens_mode;
					hyper_count++;

					/*
					 * compute the marginals for this point. check storage 
					 */
					CHECK_DENS_STORAGE;
					if (nhyper > 0) {
						weights[dens_count] = log_dens;
					} else {
						weights[dens_count] = 0.0;
					}
					izs[dens_count] = Calloc(nhyper, double);

					for (i = 0; i < nhyper; i++) {
						izs[dens_count][i] = z[i];
					}

					tref = GMRFLib_cpu();

					if (need_Qinv) {
						GMRFLib_ai_add_Qinv_to_ai_store(ai_store);	/* add Qinv */
						GMRFLib_ai_si(ai_par, log_dens, theta, nhyper, graph, ai_store);
					}
					ai_store->neff = GMRFLib_AI_STORE_NEFF_NOT_COMPUTED;

					if (run_with_omp) {
#pragma omp parallel
						{
							GMRFLib_ai_store_tp *ai_store_id = GMRFLib_duplicate_ai_store(ai_store);

#pragma omp for private(i) schedule(static) nowait
							for (i = 0; i < compute_n; i++) {
								int ii = compute_idx[i];
								GMRFLib_density_tp *cpodens = NULL;

								GMRFLib_thread_id = 0;

								GMRFLib_ai_marginal_hidden(&dens[ii][dens_count],
											   (cpo
											    && (d[ii] || ai_par->cpo_manual) ? &cpodens : NULL), ii, x, b,
											   c, mean, d, loglFunc, loglFunc_arg,
											   fixed_value, graph, Qfunc, Qfunc_arg, constr, ai_par, ai_store_id);
								double *xx_mode = ai_store_id->mode;
								COMPUTE;
								GMRFLib_free_density(cpodens);
							}
							GMRFLib_free_ai_store(ai_store_id);
						}
					} else {
						GMRFLib_ai_store_tp *ai_store_id = NULL;
						for (i = 0; i < compute_n; i++) {
							int ii = compute_idx[i];
							GMRFLib_density_tp *cpodens = NULL;

							GMRFLib_thread_id = 0;

							GMRFLib_ai_marginal_hidden(&dens[ii][dens_count],
										   (cpo
										    && (d[ii] || ai_par->cpo_manual) ? &cpodens : NULL), ii, x, b, c, mean,
										   d, loglFunc, loglFunc_arg, fixed_value, graph,
										   Qfunc, Qfunc_arg, constr, ai_par, ai_store);

							double *xx_mode = ai_store->mode;
							COMPUTE;
							GMRFLib_free_density(cpodens);
						}
					}
					if (GMRFLib_ai_INLA_userfunc0) {
						userfunc_values[dens_count] = GMRFLib_ai_INLA_userfunc0(ai_store->problem, theta, nhyper);
					}
					COMPUTE_LINDENS(ai_store);
					tu = GMRFLib_cpu() - tref;
					if (ai_par->fp_log) {
						fprintf(ai_par->fp_log, " %.2fs\n", tu);
					}
					dens_count++;
				}
			}

			/*
			 * END OF GMRFLib_AI_INT_STRATEGY_CCD 
			 */
		} else {
			/*
			 * integrate using GRID 
			 */
			if (run_with_omp) {
				/*
				 * new code which parallise over configurations 
				 */
				GMRFLib_ai_pool_tp *pool = NULL;
				unsigned int kk;

				GMRFLib_ai_pool_init(&pool, ai_par, nhyper);

				GMRFLib_ASSERT(dens_count == 0, GMRFLib_ESNH);
				GMRFLib_ASSERT(hyper_count == 0, GMRFLib_ESNH);

#pragma omp parallel for private(i, log_dens, tref, tu, ierr) schedule(static)
				for (kk = 0; kk < pool->nconfig; kk++) {
					GMRFLib_ai_store_tp *ai_store_id = NULL;
					GMRFLib_density_tp **dens_local = NULL;
					double *z_local = NULL, *theta_local = NULL, *cpo_theta_local = NULL, *pit_theta_local = NULL, *failure_theta_local = NULL,
					    *userfunc_values_local = NULL, weights_local, val, *deviance_theta_local = NULL, neff_local = 0.0;
					int err, *iz_local = NULL;
					size_t idx;

					iz_local = Calloc(nhyper, int);
					err = GMRFLib_ai_pool_get(pool, iz_local, &idx);
					/*
					 * if we get a new config, then go on, otherwise, do nothing 
					 */
					if (err == GMRFLib_SUCCESS) {
						tref = GMRFLib_cpu();
						GMRFLib_thread_id = omp_get_thread_num();

						if (omp_in_parallel()) {
							if (!ais[GMRFLib_thread_id]) {
								ais[GMRFLib_thread_id] = GMRFLib_duplicate_ai_store(ai_store);
							}
							ai_store_id = ais[GMRFLib_thread_id];
						} else {
							ai_store_id = ai_store;	/* the common one */
						}

						z_local = Calloc(nhyper, double);
						theta_local = Calloc(nhyper, double);
						for (i = 0; i < nhyper; i++) {
							z_local[i] = iz_local[i] * ai_par->dz;
						}
						GMRFLib_ai_z2theta(theta_local, nhyper, theta_mode, z_local, eigen_values, eigen_vectors);
						GMRFLib_domin_f_intern(theta_local, &log_dens, &ierr, ai_store_id);
						log_dens *= -1.0;

						val = log_dens - log_dens_mode;
						if (-val > ai_par->diff_log_dens) {
							GMRFLib_ai_pool_set(pool, idx, val);
							if (ai_par->fp_log) {
#pragma omp critical
								{
									fprintf(ai_par->fp_log, "\tconfig %2d=[", config_count++);
									for (i = 0; i < nhyper; i++) {
										fprintf(ai_par->fp_log, " %5.2f", z_local[i]);
									}
									fprintf(ai_par->fp_log, "] log(rel.dens)=%5.2f, reject, %.2fs\n", val,
										GMRFLib_cpu() - tref);
								}
							}
						} else {
							/*
							 * compute the marginals for this point. check storage 
							 */
							weights_local = log_dens;
							if (need_Qinv) {
								GMRFLib_ai_add_Qinv_to_ai_store(ai_store_id);	/* add Qinv */
								GMRFLib_ai_si(ai_par, log_dens, theta_local, nhyper, graph, ai_store_id);
							}
							ai_store_id->neff = GMRFLib_AI_STORE_NEFF_NOT_COMPUTED;
							dens_local = Calloc(graph->n, GMRFLib_density_tp *);
							if (cpo) {
								cpo_theta_local = Calloc(graph->n, double);
								pit_theta_local = Calloc(graph->n, double);
								failure_theta_local = Calloc(graph->n, double);
							}
							if (dic) {
								deviance_theta_local = Calloc(graph->n, double);
							}
							for (i = 0; i < compute_n; i++) {
								GMRFLib_density_tp *cpodens = NULL;
								int ii;
								double *xx_mode = NULL;

								ii = compute_idx[i];

								GMRFLib_ai_marginal_hidden(&dens_local[ii],
											   (cpo && (d[ii] || ai_par->cpo_manual) ? &cpodens : NULL), ii, x, b,
											   c, mean, d, loglFunc, loglFunc_arg,
											   fixed_value, graph, Qfunc, Qfunc_arg, constr, ai_par, ai_store_id);
								xx_mode = ai_store_id->mode;

								COMPUTE_LOCAL;
								GMRFLib_free_density(cpodens);
							}
							if (GMRFLib_ai_INLA_userfunc0) {
								userfunc_values_local = GMRFLib_ai_INLA_userfunc0(ai_store_id->problem, theta_local, nhyper);
							}
							tu = GMRFLib_cpu() - tref;

#pragma omp critical
							{
								int ii;
								if (ai_par->fp_log) {
									{
										fprintf(ai_par->fp_log, "\tconfig %2d=[", config_count++);
										for (i = 0; i < nhyper; i++) {
											fprintf(ai_par->fp_log, " %5.2f", z_local[i]);
										}
										fprintf(ai_par->fp_log, "] log(rel.dens)=%5.2f, [%1d] accept, compute,", val,
											omp_get_thread_num());
										fprintf(ai_par->fp_log, " %.2fs\n", tu);
									}
								}
								CHECK_DENS_STORAGE;
								CHECK_HYPER_STORAGE;

								weights[dens_count] = weights_local;
								memcpy(&hyper_z[hyper_count * nhyper], z_local, nhyper * sizeof(double));
								hyper_ldens[hyper_count] = log_dens - log_dens_mode;
								izs[dens_count] = Calloc(nhyper, double);
								memcpy(izs[dens_count], z_local, nhyper * sizeof(double));
								neff[dens_count] = neff_local;

								for (i = 0; i < compute_n; i++) {
									ii = compute_idx[i];
									dens[ii][dens_count] = dens_local[ii];
								}
								COMPUTE_LINDENS(ai_store_id);
								if (cpo) {
									for (i = 0; i < compute_n; i++) {
										ii = compute_idx[i];
										if (d[ii] || ai_par->cpo_manual) {
											cpo_theta[ii][dens_count] = cpo_theta_local[ii];
											pit_theta[ii][dens_count] = pit_theta_local[ii];
											failure_theta[ii][dens_count] = failure_theta_local[ii];
										}
									}
								}
								if (dic) {
									for (i = 0; i < compute_n; i++) {
										ii = compute_idx[i];
										if (d[ii]) {
											deviance_theta[ii][dens_count] = deviance_theta_local[ii];
										}
									}
								}
								if (GMRFLib_ai_INLA_userfunc0) {
									userfunc_values[dens_count] = userfunc_values_local;
								}

								hyper_count++;
								dens_count++;
							}
						}
					}
					Free(cpo_theta_local);
					Free(dens_local);
					Free(deviance_theta_local);
					Free(iz_local);
					Free(pit_theta_local);
					Free(failure_theta_local);
					Free(theta_local);
					Free(z_local);
				}
			} else {
				k_max = Calloc(nhyper, int);
				k_maxx = Calloc(nhyper, int);
				k_min = Calloc(nhyper, int);
				k_minn = Calloc(nhyper, int);

				for (k = 0; k < nhyper; k++) {
					for (dir = -1; dir < 2; dir += 2) {

						if (ai_par->fp_log) {
							fprintf(ai_par->fp_log, "Search: coordinate %d direction %d\n", k, dir);
						}

						memset(z, 0, nhyper * sizeof(double));
						memset(iz, 0, nhyper * sizeof(int));
						while (1) {
							z[k] += dir * ai_par->dz;
							iz[k] += dir;

							/*
							 * note that _domin_f stores calculations in the ai_store even though its not in the argument
							 * list 
							 */
							GMRFLib_ai_z2theta(theta, nhyper, theta_mode, z, eigen_values, eigen_vectors);
							GMRFLib_domin_f(theta, &log_dens, &ierr);
							log_dens *= -1.0;
							tag = GMRFLib_ai_tag(iz, nhyper);

							/*
							 * register the density for the marginal of the hyperparameters computations. first check space. 
							 */
							CHECK_HYPER_STORAGE;
							for (i = 0; i < nhyper; i++) {
								hyper_z[hyper_count * nhyper + i] = 0.0;
							}
							hyper_z[hyper_count * nhyper + k] = z[k];
							hyper_ldens[hyper_count] = log_dens - log_dens_mode;
							hyper_count++;

							k_maxx[k] = IMAX(k_maxx[k], iz[k]);	/* count all points... */
							k_minn[k] = IMIN(k_minn[k], iz[k]);

							/*
							 * print out and check if we should stop or not 
							 */
							if (ai_par->fp_log) {
								fprintf(ai_par->fp_log, "\tconfig %2d=[%s] log(rel.dens)=%6.2f,",
									config_count++, tag, log_dens - log_dens_mode);
							}
							if (log_dens_mode - log_dens > ai_par->diff_log_dens) {
								if (ai_par->fp_log) {
									fprintf(ai_par->fp_log, " diff to large, stop searching\n");
								}
								Free(tag);
								break;
							} else {
								if (ai_par->fp_log) {
									fprintf(ai_par->fp_log, " accept, compute,");
								}
							}
							map_strd_set(&hash_table, tag, log_dens);

							k_max[k] = IMAX(k_max[k], iz[k]);	/* count 'accepted' points... */
							k_min[k] = IMIN(k_min[k], iz[k]);

							/*
							 * compute the marginals for this point. check first that the storage is ok 
							 */
							CHECK_DENS_STORAGE;

							weights[dens_count] = log_dens;
							izs[dens_count] = Calloc(nhyper, double);
							for (i = 0; i < nhyper; i++) {
								izs[dens_count][i] = (double) iz[i];
							}

							tref = GMRFLib_cpu();
							if (need_Qinv) {
								GMRFLib_ai_add_Qinv_to_ai_store(ai_store);	/* add Qinv */
								GMRFLib_ai_si(ai_par, log_dens, theta, nhyper, graph, ai_store);
							}
							ai_store->neff = GMRFLib_AI_STORE_NEFF_NOT_COMPUTED;

							if (run_with_omp) {
#pragma omp parallel
								{
									/*
									 * let each thread gets its own (temporary) copy of ai_store; then just split the indices 
									 */
									GMRFLib_ai_store_tp *ai_store_id = GMRFLib_duplicate_ai_store(ai_store);

#pragma omp for private(i) schedule(static) nowait
									for (i = 0; i < compute_n; i++) {
										int ii = compute_idx[i];
										GMRFLib_density_tp *cpodens = NULL;

										GMRFLib_thread_id = 0;

										GMRFLib_ai_marginal_hidden(&dens[ii][dens_count],
													   (cpo
													    && (d[ii] || ai_par->cpo_manual) ? &cpodens :
													    NULL), ii, x, b, c,
													   mean, d, loglFunc,
													   loglFunc_arg,
													   fixed_value, graph,
													   Qfunc, Qfunc_arg, constr, ai_par, ai_store_id);
										double *xx_mode = ai_store_id->mode;
										COMPUTE;
										GMRFLib_free_density(cpodens);
									}
									GMRFLib_free_ai_store(ai_store_id);
								}
							} else {
								GMRFLib_ai_store_tp *ai_store_id = NULL;
								for (i = 0; i < compute_n; i++) {
									int ii = compute_idx[i];
									GMRFLib_density_tp *cpodens = NULL;

									GMRFLib_thread_id = 0;

									GMRFLib_ai_marginal_hidden(&dens[ii][dens_count],
												   (cpo && (d[ii] || ai_par->cpo_manual) ? &cpodens : NULL),
												   ii, x, b, c, mean, d,
												   loglFunc, loglFunc_arg,
												   fixed_value, graph, Qfunc, Qfunc_arg, constr, ai_par, ai_store);
									double *xx_mode = ai_store->mode;
									COMPUTE;
									GMRFLib_free_density(cpodens);
								}
							}
							if (GMRFLib_ai_INLA_userfunc0) {
								userfunc_values[dens_count] = GMRFLib_ai_INLA_userfunc0(ai_store->problem, theta, nhyper);
							}
							COMPUTE_LINDENS(ai_store);
							tu = GMRFLib_cpu() - tref;
							if (ai_par->fp_log) {
								fprintf(ai_par->fp_log, " %.2fs\n", tu);
							}
							dens_count++;
						}
					}
				}

				len = Calloc(nhyper, int);

				for (i = 0, len_length = 1; i < nhyper; i++) {
					len[i] = 1 + k_max[i] - k_min[i];
					len_length *= len[i];
				}

				iz_axes = Calloc(nhyper, int);
				izz = Calloc(nhyper, int);
				memset(izz, 0, nhyper * sizeof(int));

				if (ai_par->fp_log) {
					fprintf(ai_par->fp_log, "Fill-in computations\n");
				}
				for (k = 0; k < len_length; k++) {
					/*
					 * this is bit special as i want for k_max=2, the ordering 0,1,2,-1,.. 
					 */
					for (i = 0; i < nhyper; i++) {
						iz[i] = (izz[i] <= k_max[i] ? izz[i] : k_max[i] - izz[i]);
						z[i] = iz[i] * ai_par->dz;
					}
					tag = GMRFLib_ai_tag(iz, nhyper);

					/*
					 * check that we're not at the coordinate axis itself 
					 */
					skip = (map_strd_ptr(&hash_table, tag) != NULL ? 1 : 0);
					if (!skip) {
						GMRFLib_ai_z2theta(theta, nhyper, theta_mode, z, eigen_values, eigen_vectors);
						GMRFLib_domin_f(theta, &log_dens, &ierr);
						log_dens *= -1.0;
						if (ai_par->fp_log) {
							fprintf(ai_par->fp_log, "\tconfig %2d=[%s] log(rel.dens)=%6.2f,",
								config_count++, tag, log_dens - log_dens_mode);
						}

						/*
						 * register the density for the marginal of the hyperparameters computations. first check space. 
						 */
						CHECK_HYPER_STORAGE;
						for (i = 0; i < nhyper; i++) {
							hyper_z[hyper_count * nhyper + i] = z[i];
						}
						hyper_ldens[hyper_count] = log_dens - log_dens_mode;
						hyper_count++;

						if (log_dens_mode - log_dens > ai_par->diff_log_dens) {
							if (ai_par->fp_log) {
								if (ai_par->skip_configurations) {
									fprintf(ai_par->fp_log, " diff to large, stop and skip skip more distant configs\n");
								} else {
									fprintf(ai_par->fp_log, " diff to large, stop\n");
								}
							}
							skip = 1;

							if (ai_par->skip_configurations) {
								/*
								 * mark all those configurations >= than this as to be skipped as well 
								 */
								GMRFLib_ai_skip_configurations(&hash_table, k, iz, izz, len, k_max, len_length, nhyper);
							}
						} else {
							if (ai_par->fp_log) {
								fprintf(ai_par->fp_log, " accept, compute,");
							}
						}
						if (!skip) {
							/*
							 * compute the marginals for this point. check first that the storage is ok 
							 */
							CHECK_DENS_STORAGE;

							weights[dens_count] = log_dens;
							izs[dens_count] = Calloc(nhyper, double);
							for (i = 0; i < nhyper; i++) {
								izs[dens_count][i] = (double) iz[i];
							}
							tref = GMRFLib_cpu();
							if (need_Qinv) {
								GMRFLib_ai_add_Qinv_to_ai_store(ai_store);	/* add Qinv */
								GMRFLib_ai_si(ai_par, log_dens, theta, nhyper, graph, ai_store);
							}
							ai_store->neff = GMRFLib_AI_STORE_NEFF_NOT_COMPUTED;

							if (run_with_omp) {
#pragma omp parallel
								{
									/*
									 * let each thread gets its own (temporary) copy of ai_store; then just split the indices 
									 */
									GMRFLib_ai_store_tp *ai_store_id = GMRFLib_duplicate_ai_store(ai_store);

#pragma omp for private(i) schedule(static) nowait
									for (i = 0; i < compute_n; i++) {
										int ii = compute_idx[i];
										GMRFLib_density_tp *cpodens = NULL;

										GMRFLib_thread_id = 0;

										GMRFLib_ai_marginal_hidden(&dens[ii][dens_count],
													   (cpo
													    && (d[ii] || ai_par->cpo_manual) ? &cpodens :
													    NULL), ii, x, b, c,
													   mean, d, loglFunc,
													   loglFunc_arg,
													   fixed_value, graph,
													   Qfunc, Qfunc_arg, constr, ai_par, ai_store_id);
										double *xx_mode = ai_store_id->mode;
										COMPUTE;
										GMRFLib_free_density(cpodens);
									}
									GMRFLib_free_ai_store(ai_store_id);
								}
							} else {
								GMRFLib_ai_store_tp *ai_store_id = NULL;
								for (i = 0; i < compute_n; i++) {
									int ii = compute_idx[i];
									GMRFLib_density_tp *cpodens = NULL;

									GMRFLib_thread_id = 0;

									GMRFLib_ai_marginal_hidden(&dens[ii][dens_count],
												   (cpo && (d[ii] || ai_par->cpo_manual) ? &cpodens : NULL),
												   ii, x, b, c, mean, d,
												   loglFunc, loglFunc_arg,
												   fixed_value, graph, Qfunc, Qfunc_arg, constr, ai_par, ai_store);

									double *xx_mode = ai_store->mode;
									COMPUTE;
									GMRFLib_free_density(cpodens);
								}
							}
							if (GMRFLib_ai_INLA_userfunc0) {
								userfunc_values[dens_count] = GMRFLib_ai_INLA_userfunc0(ai_store->problem, theta, nhyper);
							}
							COMPUTE_LINDENS(ai_store);
							tu = GMRFLib_cpu() - tref;
							if (ai_par->fp_log) {
								fprintf(ai_par->fp_log, " %.2fs\n", tu);
							}
							dens_count++;
						}
					}
					Free(tag);

					/*
					 * compute the next configuration 
					 */
					for (i = nhyper - 1; i >= 0; i--) {
						if ((izz[i] = (izz[i] + 1) % len[i])) {
							break;
						}
					}
				}
			}
			/*
			 * END OF GMRFLib_AI_INT_STRATEGY_GRID 
			 */
		}
		/*
		 * END OF nhyper>0 
		 */
	} else {

		/*
		 * this is the case for nhyper = 0 
		 */

		if (nlin) {
			lin_dens = Calloc(1, GMRFLib_density_tp **);
		}
		if (need_Qinv) {
			/*
			 * In this case the contents of ai_store is NULL 
			 */
			GMRFLib_ai_add_Qinv_to_ai_store(ai_store);	/* add Qinv if required */
			GMRFLib_ai_si(ai_par, 0.0, NULL, 0, graph, ai_store);
		}
		ai_store->neff = GMRFLib_AI_STORE_NEFF_NOT_COMPUTED;

		if (run_with_omp) {
#pragma omp parallel
			{
				GMRFLib_ai_store_tp *ai_store_id = GMRFLib_duplicate_ai_store(ai_store);

#pragma omp for private(i) schedule(static) nowait
				for (i = 0; i < compute_n; i++) {
					int ii = compute_idx[i];
					GMRFLib_density_tp *cpodens = NULL;

					GMRFLib_thread_id = 0;

					GMRFLib_ai_marginal_hidden(&dens[ii][dens_count],
								   (cpo && (d[ii] || ai_par->cpo_manual) ? &cpodens : NULL), ii, x, b, c, mean, d,
								   loglFunc, loglFunc_arg, fixed_value, graph, Qfunc, Qfunc_arg, constr, ai_par, ai_store_id);

					double *xx_mode = ai_store_id->mode;
					COMPUTE;
					GMRFLib_free_density(cpodens);
				}
				GMRFLib_free_ai_store(ai_store_id);
			}
		} else {
			GMRFLib_ai_store_tp *ai_store_id = NULL;
			for (i = 0; i < compute_n; i++) {
				int ii = compute_idx[i];
				GMRFLib_density_tp *cpodens = NULL;

				GMRFLib_thread_id = 0;

				GMRFLib_ai_marginal_hidden(&dens[ii][dens_count], (cpo && (d[ii] || ai_par->cpo_manual) ? &cpodens : NULL),
							   ii, x, b, c, mean, d,
							   loglFunc, loglFunc_arg, fixed_value, graph, Qfunc, Qfunc_arg, constr, ai_par, ai_store);

				double *xx_mode = ai_store->mode;

				COMPUTE;
				GMRFLib_free_density(cpodens);
			}
		}
		if (GMRFLib_ai_INLA_userfunc0) {
			userfunc_values[dens_count] = GMRFLib_ai_INLA_userfunc0(ai_store->problem, theta, nhyper);
		}
		COMPUTE_LINDENS(ai_store);
		weights[dens_count] = 0.0;
		dens_count++;

		/*
		 * END OF nhyper == 0 
		 */
	}

	/*
	 * collect terms and build the output density 
	 */
	if (density) {
		*density = Calloc(graph->n, GMRFLib_density_tp *);
	}
	if (gdensity) {
		*gdensity = Calloc(graph->n, GMRFLib_density_tp *);
	}
	if (density_hyper) {
		*density_hyper = Calloc(nhyper, GMRFLib_density_tp *);
	}
	if (dlin && nlin) {
		*dlin = Calloc(nlin, GMRFLib_density_tp *);
	}

	/*
	 * if ai_par->adj_weights is false, then adj_weights and weights are the same. 
	 */
	GMRFLib_adjust_vector(weights, dens_count);
	for (j = 0; j < dens_count; j++) {
		weights[j] = exp(weights[j]);
	}
	adj_weights = Calloc(dens_count, double);

	if (ai_par->adjust_weights && (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_GRID
				       || (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_CCD && nhyper == 1))) {
		/*
		 * we may chose to adjust the weights whe int_strategy == GMRFLib_AI_INT_STRATEGY_GRID. the
		 * GMRFLib_AI_INT_STRATEGY_CCD case, this is already done above, except for the case where nhyper=1.
		 */
		GMRFLib_ai_adjust_integration_weights(adj_weights, weights, izs, dens_count, nhyper, ai_par->dz);
	} else {
		memcpy(adj_weights, weights, dens_count * sizeof(double));
	}

	if (ai_par->fp_log) {
		fprintf(ai_par->fp_log, "Combine the densities with relative weights:\n");
		for (j = 0; j < dens_count; j++) {
			fprintf(ai_par->fp_log, "\tconfig %2d/%2d=[", j, dens_count);
			for (k = 0; k < nhyper; k++) {
				fprintf(ai_par->fp_log, " %5.2f", izs[j][k]);
			}
			fprintf(ai_par->fp_log, "] weight = %.3f", weights[j]);
			if (ai_par->adjust_weights && (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_GRID)) {
				fprintf(ai_par->fp_log, " adjusted weight = %.3f", adj_weights[j]);
			}
			if (ai_par->compute_nparam_eff) {
				fprintf(ai_par->fp_log, "  neff = %.2f\n", neff[j]);
			} else {
				fprintf(ai_par->fp_log, "\n");
			}
		}
	}

	/*
	 * normalise the weights so they sum to 1.
	 */
	{
		double sum = 0.0;

		for (j = 0; j < dens_count; j++) {
			sum += adj_weights[j];
		}
		sum = 1.0 / sum;
		for (j = 0; j < dens_count; j++) {
			adj_weights[j] *= sum;
		}
	}

	if (density || gdensity) {
#pragma omp parallel for private(j) schedule(static)
		for (j = 0; j < compute_n; j++) {
			int ii = compute_idx[j];
			GMRFLib_density_tp *dens_combine, *gdens_combine;

			GMRFLib_density_combine((density ? &dens_combine : NULL), (gdensity ? &gdens_combine : NULL), dens_count, dens[ii], adj_weights);
			if (density) {
				(*density)[ii] = dens_combine;
			}
			if (gdensity) {
				(*gdensity)[ii] = gdens_combine;
			}
		}
	}
	if (dlin && nlin) {
		GMRFLib_density_tp **dtmp, *dcombine;

		assert(lin_dens);
		dtmp = Calloc(dens_count, GMRFLib_density_tp *);
		for (j = 0; j < nlin; j++) {
			/*
			 * I need to do this as the storage is wrong.... 
			 */
			for (k = 0; k < dens_count; k++) {
				dtmp[k] = lin_dens[k][j];
			}
			GMRFLib_density_combine(&dcombine, NULL, dens_count, dtmp, adj_weights);
			(*dlin)[j] = dcombine;
		}
		Free(dtmp);
	}
	if (ai_par->fp_log) {
		fprintf(ai_par->fp_log, "Done.\n");
	}

	if (ai_par->compute_nparam_eff) {
		/*
		 * compute the effective number of parameters 
		 */
		double neff_eval = 0.0, neff2_eval = 0.0, neff_sd, ndata = 0.0;

		for (j = 0; j < dens_count; j++) {
			neff_eval += adj_weights[j] * neff[j];
			neff2_eval += adj_weights[j] * SQR(neff[j]);
		}
		neff_sd = sqrt(DMAX(0.0, neff2_eval - SQR(neff_eval)));

		if (d) {
			for (i = 0; i < graph->n; i++) {
				if (d[i]) {
					if (!fixed_value || (fixed_value && !fixed_value[i])) {
						ndata++;
					}
				}
			}
		}
		if (ai_par->fp_log) {
			if (ndata && !ISZERO(neff_eval)) {
				printf("Expected effective number of parameters: %.3f(%.3f),  eqv.#replicates: %.3f\n", neff_eval, neff_sd, ndata / neff_eval);
			} else {
				printf("Expected effective number of parameters: %.3f(%.3f)\n", neff_eval, neff_sd);
			}
		}

		if (neffp) {
			neffp->mean = neff_eval;
			neffp->stdev = neff_sd;
			neffp->nrep = ndata / neff_eval;
		}
	}

	if (cpo) {
		/*
		 * first we need to compute the normalising constants for \pi(theta|y_i) for each i. This we call Z[i].
		 * 
		 * Note that \pi(theta_j | y_{-i}) = adj_weights[j] / cpo_theta[i][j] / Z[i]; 
		 */
		double *Z = Calloc(graph->n, double);

		for (j = 0; j < compute_n; j++) {
			int jj, ii;

			ii = compute_idx[j];
			if (cpo_theta[ii]) {
				for (jj = 0; jj < dens_count; jj++) {
					if (cpo_theta[ii][jj]) /* we ignore those that have failed */
						Z[ii] += adj_weights[jj] / cpo_theta[ii][jj];
				}
			}
		}

		for (j = 0; j < compute_n; j++) {
			int ii, jj;
			double evalue, evalue2, evalue_one;

			ii = compute_idx[j];
			if (cpo_theta[ii]) {
				(*cpo)->value[ii] = Calloc(1, double);

				for (jj = 0, evalue = evalue_one = 0.0; jj < dens_count; jj++) {
					if (cpo_theta[ii][jj]) {
						evalue += cpo_theta[ii][jj] * adj_weights[jj] / cpo_theta[ii][jj] / Z[ii];
						evalue_one += adj_weights[jj] / cpo_theta[ii][jj] / Z[ii];
					}
				}
				if (evalue_one) {
					(*cpo)->value[ii][0] = evalue / evalue_one;
				} else {
					(*cpo)->value[ii][0] = 0.0;
				}
			} else {
				(*cpo)->value[ii] = NULL;
			}

			if (cpo_theta[ii]) {
				(*cpo)->pit_value[ii] = Calloc(1, double);
				(*cpo)->failure[ii] = Calloc(1, double);
				for (jj = 0, evalue = evalue2 = evalue_one = 0.0; jj < dens_count; jj++) {
					if (cpo_theta[ii][jj]) {
						evalue += pit_theta[ii][jj] * adj_weights[jj] / cpo_theta[ii][jj] / Z[ii];
						evalue_one += adj_weights[jj] / cpo_theta[ii][jj] / Z[ii];
					}
					/*
					 * this is defined over the unadjusted weights 
					 */
					evalue2 += failure_theta[ii][jj] * adj_weights[jj];
				}
				if (evalue_one) {
					evalue = TRUNCATE(evalue / evalue_one, 0.0, 1.0);
				} else {
					evalue = 0.0;
				}
				(*cpo)->pit_value[ii][0] = evalue;
				(*cpo)->failure[ii][0] = evalue2;
			} else {
				(*cpo)->pit_value[ii] = NULL;
				(*cpo)->failure[ii] = NULL;
			}
		}
		(*cpo)->mean_value = (*cpo)->gmean_value = 0.0;
		if (compute_n) {
			int count = 0;

			for (j = 0; j < compute_n; j++) {
				int ii = compute_idx[j];

				if (cpo_theta[ii]) {
					(*cpo)->mean_value += *((*cpo)->value[ii]);
					(*cpo)->gmean_value += log(*((*cpo)->value[ii]));
					count++;
				}
			}
			if (count) {
				(*cpo)->mean_value /= (double) count;
				(*cpo)->gmean_value = exp((*cpo)->gmean_value / (double) count);
			} else {
				(*cpo)->mean_value = (*cpo)->gmean_value = 0.0;
			}
		}
		Free(Z);
	}

	if (dic) {
		double mean_deviance = 0.0, deviance_mean = 0.0, *x_vec = NULL;

		SET_THETA_MODE;

		/*
		 * need this for loglFunc() we need that compute is TRUE for all indices that enters loglFunc. There is no way to check this here. 
		 */
		x_vec = Calloc(graph->n, double);
		for (j = 0; j < compute_n; j++) {
			int ii = compute_idx[j];
			x_vec[ii] = (*density)[ii]->user_mean;
		}

		for (j = 0; j < compute_n; j++) {
			double md, dm;
			int ii = compute_idx[j];

			if (d[ii]) {
				int jj;
				double evalue, sum, logll;

				for (jj = 0, evalue = sum = 0.0; jj < dens_count; jj++) {
					evalue += deviance_theta[ii][jj] * adj_weights[jj];
					sum += adj_weights[jj];
				}
				md = evalue / sum;

				if (!(density && (*density)[ii])) {
					fprintf(stderr, "\n\n\nFIXME FIXME!!!!!!!!\n\n\n");
					abort();
				}
				loglFunc(&logll, &((*density)[ii]->user_mean), 1, ii, x_vec, loglFunc_arg);
				dm = -2.0 * logll;
			} else {
				dm = md = 0.0;
			}

			deviance_mean += dm;
			mean_deviance += md;
		}
		Free(x_vec);

		dic->mean_of_deviance = mean_deviance;
		dic->deviance_of_mean = deviance_mean;
		dic->p = mean_deviance - deviance_mean;
		dic->dic = dic->p + mean_deviance;

		if (ai_par->fp_log) {
			fprintf(ai_par->fp_log, "DIC:\n");
			fprintf(ai_par->fp_log, "\tMean of Deviance................. %g\n", dic->mean_of_deviance);
			fprintf(ai_par->fp_log, "\tDeviance at Mean................. %g\n", dic->deviance_of_mean);
			fprintf(ai_par->fp_log, "\tEffective number of parameters... %g\n", dic->p);
			fprintf(ai_par->fp_log, "\tDIC.............................. %g\n", dic->dic);
		}
	}
	if (GMRFLib_ai_INLA_userfunc0 && GMRFLib_ai_INLA_userfunc0_dim > 0) {

		int dim = GMRFLib_ai_INLA_userfunc0_dim;
		GMRFLib_ai_INLA_userfunc0_density = Calloc(dim, GMRFLib_density_tp *);

		for (j = 0; j < dim; j++) {
			double val = 0.0, wsum = 0.0, val2 = 0.0, mmean, vvar, ssd;

			for (i = 0; i < dens_count; i++) {
				wsum += adj_weights[i];
				val += adj_weights[i] * userfunc_values[i][j];
				val2 += adj_weights[i] * SQR(userfunc_values[i][j]);
			}

			mmean = val / wsum;
			vvar = DMAX(DBL_EPSILON, val2 / wsum - SQR(mmean));
			ssd = sqrt(vvar);
			GMRFLib_density_create_normal(&(GMRFLib_ai_INLA_userfunc0_density[j]), 0.0, 1.0, mmean, ssd);

			// The densities are to ``unstable'' to fit...
			// GMRFLib_density_create(&(GMRFLib_ai_INLA_userfunc0_density[j]), GMRFLib_DENSITY_TYPE_SCGAUSSIAN,
			// dens_count, values, ldens, mmean, ssd, GMRFLib_TRUE);

			if (0) {
				if (ai_par->fp_log) {
					fprintf(ai_par->fp_log, "User-defined function0[%d] = %.12f (%.12f)\n", j, mmean, ssd);
				}
			}
		}
	}

	/*
	 * Compute the marginal likelihood; compute both the Gaussian approximatin and a non-parametric one. The marginal likelhood is the
	 * normalising constant for the posterior marginal for \theta. 
	 */
	if (marginal_likelihood && nhyper) {
		marginal_likelihood->marginal_likelihood_gaussian_approx = 0.5 * nhyper * log(2.0 * M_PI) + log_dens_mode;
		for (i = 0; i < nhyper; i++) {
			marginal_likelihood->marginal_likelihood_gaussian_approx -= 0.5 * log(gsl_vector_get(eigen_values, (unsigned int) i));
		}

		if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_CCD) {
			marginal_likelihood->marginal_likelihood_integration = marginal_likelihood->marginal_likelihood_gaussian_approx;
		} else {
			double integral = 0.0, log_jacobian = 0.0;

			for (j = 0; j < dens_count; j++) {
				integral += weights[j];
			}
			integral *= ai_par->dz;
			for (i = 0; i < nhyper; i++) {
				log_jacobian -= 0.5 * log(gsl_vector_get(eigen_values, (unsigned int) i));
			}
			marginal_likelihood->marginal_likelihood_integration = log(integral) + log_jacobian + log_dens_mode;
		}
		if (ai_par->fp_log) {
			fprintf(ai_par->fp_log, "Marginal likelihood: Integration %f Gaussian-approx %f\n",
				marginal_likelihood->marginal_likelihood_integration, marginal_likelihood->marginal_likelihood_gaussian_approx);
		}
	}

	/*
	 * compute the posterior marginals for each hyperparameter, if possible 
	 */
	if (hyper_z && density_hyper && nhyper) {
		if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_EMPIRICAL_BAYES || ai_par->interpolator == GMRFLib_AI_INTERPOLATOR_GAUSSIAN) {
			/*
			 * Just use the modal values and the stdev's found from the Hessian. 
			 */
			for (k = 0; k < nhyper; k++) {
				GMRFLib_density_create_normal(&((*density_hyper)[k]), 0.0, 1.0, theta_mode[k], sqrt(inverse_hessian[k + nhyper * k]));
			}
		} else {
			if (ai_par->fp_log) {
				fprintf(ai_par->fp_log, "Compute the marginal for each of the %1d hyperparameters\n", nhyper);
				fprintf(ai_par->fp_log, "Interpolation method: %s\n", INTERPOLATOR_NAME(ai_par->interpolator));
			}
			/*
			 * add points one step outwards to put a guard-zone around. we do that just by looping over all possible
			 * configurations, and then check if we're on the boundary. the process may be performed twice to pindown the
			 * density. Note: only for strategy == GRID.
			 */
			int ntimes;

			GMRFLib_ai_interpolator_tp interpol = (ai_par->interpolator == GMRFLib_AI_INTERPOLATOR_AUTO ?
							       /*
							        *  ...AUTO is requested. unless we *can* use the GRIDSUM
							        *  approach, we use default the CCD
							        */
							       (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_GRID
								&& ai_par->hessian_force_diagonal ? GMRFLib_AI_INTERPOLATOR_GRIDSUM : GMRFLib_AI_INTERPOLATOR_CCD)
							       /*
							        *  ...if not AUTO, then use whatever is requested
							        */
							       : ai_par->interpolator);

			if (interpol != GMRFLib_AI_INTERPOLATOR_CCD && interpol != GMRFLib_AI_INTERPOLATOR_GRIDSUM) {
				double bvalue = GMRFLib_min_value(hyper_ldens, hyper_count);

				for (ntimes = 0; ntimes < 2; ntimes++) {
					int guard_count = 0;

					for (i = 0, len_length = 1; i < nhyper; i++) {
						k_maxx[i]++;   /* add one to construct the guard-zone */
						k_minn[i]--;   /* subtract one ... */
						len[i] = 1 + k_maxx[i] - k_minn[i];
						len_length *= len[i];
					}
					memset(izz, 0, nhyper * sizeof(int));
					for (k = 0; k < len_length; k++) {
						for (i = 0; i < nhyper; i++) {
							iz[i] = (izz[i] <= k_maxx[i] ? izz[i] : k_maxx[i] - izz[i]);
							z[i] = iz[i] * ai_par->dz;
						}

						/*
						 * are we ON the boundary? 
						 */
						int test1 = 0;

						for (i = 0; i < nhyper; i++) {
							test1 += (iz[i] == k_maxx[i] || iz[i] == k_minn[i] ? 1 : 0);
						}
						if (test1 >= IMAX(1, nhyper - 1)) {
							/*
							 * yes, we are; add this configuration 
							 */
							guard_count++;
							CHECK_HYPER_STORAGE;
							for (i = 0; i < nhyper; i++) {
								hyper_z[hyper_count * nhyper + i] = z[i];
							}
							hyper_ldens[hyper_count] = (ntimes == 0 ? 1.5 : 3.0) * bvalue;
							hyper_count++;
						}
						/*
						 * compute the next configuration 
						 */
						for (i = nhyper - 1; i >= 0; i--) {
							if ((izz[i] = (izz[i] + 1) % len[i])) {
								break;
							}
						}
					}
					if (ai_par->fp_log) {
						fprintf(ai_par->fp_log,
							"\tLoop %1d: Added %1d guard-points with log-dens-diff equal to %f\n",
							ntimes, guard_count, hyper_ldens[hyper_count - 1]);
					}
				}
			}
			double *std_stdev_theta = Calloc(nhyper, double);

			for (k = 0; k < nhyper; k++) {
				std_stdev_theta[k] = sqrt(inverse_hessian[k + nhyper * k]);
			}

			if (ai_par->fp_hyperparam) {
				/*
				 * write out the hole set 
				 */
				double *theta_tmp = Calloc(nhyper, double), log_jacobian = 0.0;

				if (eigen_values) {
					for (k = 0; k < nhyper; k++) {
						log_jacobian -= 0.5 * log(gsl_vector_get(eigen_values, (unsigned int) k));
					}
				}
				for (k = 0; k < hyper_count; k++) {
					int kk;

					GMRFLib_ai_z2theta(theta_tmp, nhyper, theta_mode, &(hyper_z[k * nhyper]), eigen_values, eigen_vectors);
					// fprintf(ai_par->fp_hyperparam, "%s: ", __GMRFLib_FuncName);
					for (kk = 0; kk < nhyper; kk++) {
						fprintf(ai_par->fp_hyperparam, " %.10g", theta_tmp[kk]);
					}
					fprintf(ai_par->fp_hyperparam, " %.10g\n", hyper_ldens[k] + log_dens_mode + log_jacobian);
				}
				fflush(ai_par->fp_hyperparam);
				Free(theta_tmp);
			}
			if (run_with_omp) {
				if (ai_par->fp_log) {
					fprintf(ai_par->fp_log,
						"\tCompute the marginal for theta[%1d] to theta[%1d] using numerical integration...\n", 0, nhyper - 1);
				}
			}
#pragma omp parallel for private(k) schedule(static)
			for (k = 0; k < nhyper; k++) {
				if (!run_with_omp) {
					if (ai_par->fp_log) {
						fprintf(ai_par->fp_log, "\tCompute the marginal for theta[%1d] using numerical integration\n", k);
					}
				}
				GMRFLib_ai_marginal_for_one_hyperparamter(&((*density_hyper)[k]), k, nhyper, hyper_count, hyper_z,
									  hyper_ldens, theta_mode, eigen_values, eigen_vectors,
									  std_stdev_theta, ai_par->dz, stdev_corr_pos, stdev_corr_neg, interpol);
			}
			if (run_with_omp) {
				if (ai_par->fp_log) {
					fprintf(ai_par->fp_log,
						"\tCompute the marginal for theta[%1d] to theta[%1d] using numerical integration... Done.\n", 0, nhyper - 1);
				}
			}
			Free(std_stdev_theta);

			if (ai_par->fp_log) {
				fprintf(ai_par->fp_log, "Compute the marginal for the hyperparameters... done.\n");
			}
		}
	}

	/*
	 * return the mode in hyperparam and in 'x'
	 */
	SET_THETA_MODE;
	if (x && x_mode) {
		memcpy(x, x_mode, graph->n * sizeof(double));
	}


	/*
	 * userfunction1 
	 */
	if (GMRFLib_ai_INLA_userfunc1)
		GMRFLib_ai_INLA_userfunc1(theta_mode, nhyper, inverse_hessian);

	/*
	 * cleanup 
	 */
	if (izs) {
		for (j = 0; j < dens_count; j++) {
			Free(izs[j]);
		}
		Free(izs);
	}
	if (lin_dens && nlin) {
		if (dens_count) {
			for (j = 0; j < dens_count; j++) {
				for (i = 0; i < nlin; i++)
					GMRFLib_free_density(lin_dens[j][i]);
				Free(lin_dens[j]);
			}
			Free(lin_dens);
		} else {
			for (i = 0; i < nlin; i++)
				GMRFLib_free_density(lin_dens[0][i]);
			Free(lin_dens[0]);
			Free(lin_dens);
		}
	}

	Free(adj_weights);
	Free(hessian);
	Free(inverse_hessian);
	Free(iz);
	Free(iz_axes);
	Free(izz);
	Free(k_max);
	Free(k_maxx);
	Free(k_min);
	Free(k_minn);
	Free(len);
	Free(stdev_corr_neg);
	Free(stdev_corr_pos);
	Free(theta);
	Free(theta_mode);
	Free(userfunc_values);
	Free(weights);
	Free(z);
	Free(neff);
	if (cpo_theta) {
		for (i = 0; i < compute_n; i++) {
			j = compute_idx[i];
			if (d[j] || ai_par->cpo_manual) {
				Free(cpo_theta[j]);
			}
		}
		Free(cpo_theta);
	}
	if (pit_theta) {
		for (i = 0; i < compute_n; i++) {
			j = compute_idx[i];
			if (d[j]) {
				Free(pit_theta[j]);
			}
		}
		Free(pit_theta);
	}
	if (failure_theta) {
		for (i = 0; i < compute_n; i++) {
			j = compute_idx[i];
			if (d[j]) {
				Free(failure_theta[j]);
			}
		}
		Free(failure_theta);
	}
	if (deviance_theta) {
		for (i = 0; i < compute_n; i++) {
			j = compute_idx[i];
			if (d[j]) {
				Free(deviance_theta[j]);
			}
		}
		Free(deviance_theta);
	}
	Free(compute_idx);
	if (free_ai_par) {
		Free(ai_par);
	}
	if (free_compute) {
		Free(compute);
	}
	if (H) {
		gsl_matrix_free(H);
	}
	if (eigen_vectors) {
		gsl_matrix_free(eigen_vectors);
	}
	if (eigen_values) {
		gsl_vector_free(eigen_values);
	}
	for (k = -1; (k = map_strd_next(&hash_table, k)) != -1;) {
		Free(hash_table.contents[k].key);	       /* the keys are alloced... */
	}

	map_strd_free(&hash_table);

	Free(hyper_z);
	Free(hyper_ldens);

	if (compute) {
		for (i = 0; i < graph->n; i++) {
			if (compute[i]) {
				for (j = 0; j < dens_count; j++) {
					GMRFLib_free_density(dens[i][j]);
				}
				Free(dens[i]);
			}
		}
		Free(dens);
	}

	if (ais) {
		for (k = 0; k < tmax; k++) {
			if (ais[k]) {
				GMRFLib_free_ai_store(ais[k]);
			}
		}
		Free(ais);
	}

	if (nhyper) {
		GMRFLib_domin_exit();
	}

	GMRFLib_LEAVE_ROUTINE;
#undef CHECK_HYPER_STORAGE
#undef CHECK_DENS_STORAGE
#undef COMPUTE_CPO_AND_DIC
#undef COMPUTE_CPO_AND_DIC_LOCAL
#undef COMPUTE
#undef COMPUTE_LOCAL
#undef COMPUTE_NEFF
#undef COMPUTE_NEFF_LOCAL
#undef ADD_LINEAR_TERM
#undef ADD_LINEAR_TERM_LOCAL

	return GMRFLib_SUCCESS;
}
GMRFLib_density_tp **GMRFLib_ai_compute_lincomb(int nlin, double *Alin, GMRFLib_ai_store_tp * ai_store, double *improved_mean)
{
	/*
	 * Compute the marginals for the linear combinations using just the Gaussians 
	 */
	GMRFLib_problem_tp *problem = ai_store->problem;
	int i, j, k, n, debug = 0;
	double *v, *vv, *AA, var, mean, imean;
	GMRFLib_density_tp **d;

	assert(problem != NULL);
	if (nlin <= 0)
		return NULL;

	n = problem->n;
	v = Calloc(n * nlin, double);

	for (i = 0; i < nlin; i++) {
		j = i * n;
		if (debug) {
			for (k = 0; k < n; k++)
				printf("Alin %d %g\n", k, Alin[j + k]);
		}
		memcpy(v + j, Alin + j, n * sizeof(double));
		GMRFLib_solve_l_sparse_matrix(v + j, &(problem->sub_sm_fact), problem->sub_graph);
	}

	/*
	 * the correction matrix due to linear constraints 
	 */
	int nc = (problem->sub_constr ? problem->sub_constr->nc : 0);
	double *var_corr = Calloc(nlin, double);
	if (nc) {
		int one = 1;
		double *w = Calloc(nc * nlin, double);
		double *ww = Calloc(nc * nlin, double);

		if (0) {
			/*
			 * SAFE CODE 
			 */

			/*
			 * w = AA^T CONSTR_M 
			 */
			for (i = 0; i < nlin; i++) {
				for (j = 0; j < nc; j++) {
					for (k = 0; k < n; k++) {
						w[i + j * nlin] += Alin[k + i * n] * problem->constr_m[k + j * n];
					}
				}
			}

			/*
			 * ww = AA * QI_AT 
			 */
			for (i = 0; i < nlin; i++) {
				for (j = 0; j < nc; j++) {
					for (k = 0; k < n; k++) {
						ww[i + j * nlin] += Alin[k + i * n] * problem->qi_at_m[k + j * n];
					}
				}
			}

			for (i = 0; i < nlin; i++) {
				for (j = 0; j < nc; j++)
					var_corr[i] += w[i + j * nlin] * ww[i + j * nlin];
			}
		} else {
			/*
			 * MORE OPTIMISED CODE. TODO: add proper gdemm() calls (I'm to confused to do this now...) 
			 */
			double *wp, *Alinp, *cp, *wwp;

			/*
			 * w = AA^T CONSTR_M 
			 */
			for (i = 0; i < nlin; i++) {
				Alinp = Alin + i * n;
				for (j = 0; j < nc; j++) {
					wp = w + i + j * nlin;
					cp = problem->constr_m + j * n;
					wp[0] = ddot_(&n, Alinp, &one, cp, &one);
				}
			}

			/*
			 * ww = AA * QI_AT 
			 */
			for (i = 0; i < nlin; i++) {
				Alinp = Alin + i * n;
				for (j = 0; j < nc; j++) {
					wwp = ww + i + j * nlin;
					cp = problem->qi_at_m + j * n;
					wwp[0] = ddot_(&n, Alinp, &one, cp, &one);
				}
			}

			for (i = 0; i < nlin; i++) {
				var_corr[i] = ddot_(&nc, w + i, &nlin, ww + i, &nlin);
			}
		}

		Free(w);
		Free(ww);
	}

	d = Calloc(nlin, GMRFLib_density_tp *);
	for (i = 0; i < nlin; i++) {
		vv = v + i * n;
		AA = Alin + i * n;
		var = mean = imean = 0.0;
		for (j = 0; j < n; j++) {
			var += vv[j] * vv[j];
			mean += AA[j] * problem->mean_constr[j];
			imean += AA[j] * improved_mean[j];
		}
		if (debug) {
			printf("BEFORE sd[%d] = %g\t mean = %g imean = %g\n", i, sqrt(var), mean, imean);
			printf("AFTER sd[%d] = %g\t mean = %g imean = %g\n", i, sqrt(DMAX(0.0, var - var_corr[i])), mean, imean);
		}
		var = DMAX(0.0, var - var_corr[i]);
		GMRFLib_density_create_normal(&d[i], (imean - mean) / sqrt(var), 1.0, mean, sqrt(var));
	}

	Free(var_corr);
	Free(v);

	return d;
}
int GMRFLib_ai_si(GMRFLib_ai_param_tp * ai_par, double logdens, double *theta, int nhyper, GMRFLib_graph_tp * graph, GMRFLib_ai_store_tp * ai_store)
{
	/*
	 * SI: write the covariances, marginals, theta to file. There are two output-formats: standard txt-format, and R-output 
	 */

	if (!ai_par || ai_par->si_directory == NULL)
		return GMRFLib_SUCCESS;

	GMRFLib_ai_si_tp *d = ai_par->si_idx;
	double *p = NULL, *pi = NULL, *pj = NULL;
	int i, j = 0, jj, k, c, use_R_format = 1, use_txt_format = 0;
	static int config[2] = { 1, 1 };
	char *fnm;
	FILE *fp;

	if (!d || d->nd <= 0) {
		return GMRFLib_SUCCESS;
	}

	if (use_txt_format) {
		FIXME("use_txt_format is not yet rewritten for si_idx");
		exit(1);
#pragma omp critical
		{
			GMRFLib_sprintf(&fnm, "%s/%s-%.4d.txt", ai_par->si_directory, "configuration", config[0]);
			c = config[0];
			config[0]++;
		}

		fp = fopen(fnm, "w");
		GMRFLib_ASSERT(fp != NULL, GMRFLib_EOPENFILE);

		/*
		 * ordinary output 
		 */
		fprintf(fp, "Configuration:\n\t%d\n", c);
		fprintf(fp, "Log-density:\n\t%.6g\n", logdens);
		fprintf(fp, "Theta:\n\t");
		for (i = 0; i < nhyper; i++)
			fprintf(fp, " %.6g", theta[i]);
		fprintf(fp, "\n");
		fprintf(fp, "Marginals (i,mean,stdev):\n");
		for (i = 0; i < graph->n; i++) {
			p = GMRFLib_Qinv_get(ai_store->problem, i, i);
			fprintf(fp, "\t%d %.6g %.6g\n", i, ai_store->problem->mean_constr[i], (p ? sqrt(*p) : 0.0));
		}
		fprintf(fp, "Covariances (i,j,Cov(i,j)):\n");
		for (i = 0; i < graph->n; i++) {
			for (jj = 0; jj < graph->nnbs[i]; jj++) {
				j = graph->nbs[i][jj];
				if (i < j) {
					p = GMRFLib_Qinv_get(ai_store->problem, i, j);
					fprintf(fp, "\t%d %d %.6g\n", i, j, (p ? *p : 0.0));
				}
			}
		}
		fclose(fp);
	}

	if (use_R_format) {
#pragma omp critical
		{
			GMRFLib_sprintf(&fnm, "%s/%s-%.4d.R", ai_par->si_directory, "configuration", config[1]);
			c = config[1];
			config[1]++;
		}

		fp = fopen(fnm, "w");
		GMRFLib_ASSERT(fp != NULL, GMRFLib_EOPENFILE);

		/*
		 * output as an R-list 
		 */

		fprintf(fp, "if (exists(\"inla.si.configuration\") && !is.list(inla.si.configuration)) inla.si.configuration = list()\n");
		fprintf(fp, "if (!exists(\"inla.si.configuration\")) inla.si.configuration = list()\n");
		fprintf(fp, "inla.si.configuration[[%1d]] = list(\n", c);
		fprintf(fp, "log.dens=c(%.6g),\n", logdens);
		fprintf(fp, "theta=c(");
		for (i = 0; i < nhyper - 1; i++)
			fprintf(fp, " %.6g,", theta[i]);
		fprintf(fp, " %.6g),\n", theta[nhyper - 1]);

		fprintf(fp, "tags=c(");
		for (k = 0; k < d->nd; k++) {
			fprintf(fp, "\"%s\",", d->tag[k]);
		}
		fseek(fp, (long) -1, SEEK_CUR);
		fprintf(fp, "),\n");

		fprintf(fp, "start=c(");
		for (k = 0; k < d->nd; k++) {
			fprintf(fp, "%1d,", d->start[k]);
		}
		fseek(fp, (long) -1, SEEK_CUR);
		fprintf(fp, "),\n");

		fprintf(fp, "len=c(");
		for (k = 0; k < d->nd; k++) {
			fprintf(fp, "%1d,", d->len[k]);
		}
		fseek(fp, (long) -1, SEEK_CUR);
		fprintf(fp, "),\n");

		fprintf(fp, "mean=list(");
		for (k = 0; k < d->nd; k++) {
			fprintf(fp, "\"%s\" = c(", d->tag[k]);
			for (i = d->start[k]; i < d->start[k] + d->len[k]; i++) {
				fprintf(fp, " %.6g,", ai_store->problem->mean_constr[i]);
			}
			fseek(fp, (long) -1, SEEK_CUR);
			fprintf(fp, "),\n");
		}
		fseek(fp, (long) -2, SEEK_CUR);		       /* yes, need also the newline */
		fprintf(fp, "),\n");

		fprintf(fp, "sd=list(");
		for (k = 0; k < d->nd; k++) {
			fprintf(fp, "\"%s\" = c(", d->tag[k]);
			for (i = d->start[k]; i < d->start[k] + d->len[k]; i++) {
				p = GMRFLib_Qinv_get(ai_store->problem, i, i);
				fprintf(fp, " %.6g,", (p ? sqrt(*p) : 0.0));
			}
			fseek(fp, (long) -1, SEEK_CUR);
			fprintf(fp, "),\n");
		}
		fseek(fp, (long) -2, SEEK_CUR);		       /* yes, need also the newline */
		fprintf(fp, "),\n");

		fprintf(fp, "cor=list(");
		for (k = 0; k < d->nd; k++) {
			fprintf(fp, "\"%s\" = c(", d->tag[k]);
			for (i = d->start[k]; i < d->start[k] + d->len[k]; i++) {
				pi = GMRFLib_Qinv_get(ai_store->problem, i, i);
				for (j = d->start[k]; j < d->start[k] + d->len[k]; j++) {
					p = GMRFLib_Qinv_get(ai_store->problem, i, j);
					pj = GMRFLib_Qinv_get(ai_store->problem, j, j);
					fprintf(fp, "%.6g,", (p ? *p / sqrt(*pi * *pj) : 0.0));
				}
			}
			fseek(fp, (long) -1, SEEK_CUR);
			fprintf(fp, "),\n");
		}
		fseek(fp, (long) -2, SEEK_CUR);		       /* yes, need also the newline */
		fprintf(fp, "))\n");
		fclose(fp);
	}

	return GMRFLib_SUCCESS;
}
int GMRFLib_ai_correct_cpodens(double *logdens, double *x, int *n, GMRFLib_ai_param_tp * ai_par)
{
	/*
	 * remove 'local' maxima at the extremes and return !GMRFLib_SUCCESS if the density is monotone otherwise GMRFLib_SUCCESS. If maximum, or nearly, is on the
	 * extremes, then flag an error.
	 */
	int idx, i;
	char *code = Calloc(*n, char);
	double mode;

	idx = *n;
	do {
		idx--;
	} while (idx > 0 && (logdens[idx] > logdens[idx - 1]));

	if (idx < *n)
		memset(&code[idx], 1, *n - idx);

	idx = -1;
	do {
		idx++;
	} while (idx < *n && (logdens[idx] > logdens[idx + 1]));

	if (idx > 0)
		memset(&code[0], 1, idx);

	for (idx = i = 0; i < *n; i++) {
		if (!code[i]) {
			logdens[idx] = logdens[i];
			x[idx] = x[i];
			idx++;
		}
	}
	*n = idx;
	Free(code);

	if (*n) {
		/*
		 * a final check 
		 */
		mode = GMRFLib_max_value(logdens, *n);
		if (mode - DMAX(logdens[0], logdens[*n - 1]) < ai_par->cpo_req_diff_logdens) {
			/*
			 * this is no good... 
			 */
			*n = 0;
		}
	}

	return GMRFLib_SUCCESS;
}
double GMRFLib_ai_cpopit_integrate(double *cpo, double *pit,
				   int idx, GMRFLib_density_tp * cpo_density, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, double *x_vec)
{
	/*
	 * cpo_density is the marginal for x_idx without y_idx, density: is the marginal for x_idx with y_idx.
	 */
	int retval, compute_cpo = 1, i, k, np = GMRFLib_faster_integration_np;
	double low, dx, dxi, *xp = NULL, *xpi = NULL, *dens = NULL, *prob = NULL, *work = NULL,
	    integral = 0.0, integral2 = 0.0, w[2] = { 4.0, 2.0 }, integral_one, *loglik = NULL, fail = 0.0;



	if (!cpo_density) {
		if (cpo) {
			*cpo = 0.0;
		}
		if (pit) {
			*pit = 1.0;
		}
		fail = 1.0;

		return fail;
	}

	retval = loglFunc(NULL, NULL, 0, idx, x_vec, loglFunc_arg);
	if (!(retval == GMRFLib_LOGL_COMPUTE_CDF || retval == GMRFLib_LOGL_COMPUTE_DERIVATIES_AND_CDF)) {
		compute_cpo = 0;
	}

	GMRFLib_ASSERT_RETVAL(np > 3, GMRFLib_ESNH, 0.0);

	work = Calloc(5 * np, double);
	xp = work;
	xpi = work + np;
	dens = work + 2 * np;
	prob = work + 3 * np;
	loglik = work + 4 * np;

	dxi = (cpo_density->x_max - cpo_density->x_min) / (np - 1.0);
	low = GMRFLib_density_std2user(cpo_density->x_min, cpo_density);
	dx = (GMRFLib_density_std2user(cpo_density->x_max, cpo_density) - low) / (np - 1.0);

	xp[0] = low;
	xpi[0] = cpo_density->x_min;
	for (i = 1; i < np; i++) {
		xp[i] = xp[0] + i * dx;
		xpi[i] = xpi[0] + i * dxi;
	}
	GMRFLib_evaluate_ndensity(dens, xpi, np, cpo_density);

	if (compute_cpo) {
		loglFunc(prob, xp, -np, idx, x_vec, loglFunc_arg);
	} else {
		memset(prob, 0, np * sizeof(double));
	}
	loglFunc(loglik, xp, np, idx, x_vec, loglFunc_arg);

	if (0) {
		if (idx == 0) {
			FIXME("idx = 0; write cpo_density");
			FILE *fp = fopen("cpo-density.dat", "w");
			for (i = 0; i < np; i++) {
				fprintf(fp, "%g %g %g\n", xp[i], dens[i], exp(loglik[i]));
			}
			fclose(fp);
			FIXME("info for cpo_dens for idx=0");
			GMRFLib_density_printf(stdout, cpo_density);
		}
	}

	for (i = 0; i < np; i++) {
		xp[i] = prob[i] * dens[i];		       /* resuse and redefine xp! */
		xpi[i] = exp(loglik[i]) * dens[i];	       /* resuse and redefine xpi! */
	}

	integral = xp[0] + xp[np - 1];
	integral2 = xpi[0] + xpi[np - 1];
	integral_one = dens[0] + dens[np - 1];
	for (i = 1, k = 0; i < np - 1; i++, k = (k + 1) % 2) {
		integral += w[k] * xp[i];
		integral2 += w[k] * xpi[i];
		integral_one += w[k] * dens[i];
	}
	if (ISZERO(integral_one)) {
		fail = 1.0;
		integral = integral2 = 0.0;
	} else {
		integral /= integral_one;
		integral2 /= integral_one;
	}
	integral = TRUNCATE(integral, 0.0, 1.0);
	integral2 = DMAX(DBL_MIN, integral2);

	if (cpo) {
		*cpo = integral2;
	}
	if (pit) {
		*pit = integral;
	}

	Free(work);
	return fail;
}
double GMRFLib_ai_dic_integrate(int idx, GMRFLib_density_tp * density, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, double *x_vec)
{
	/*
	 * compute the integral of -2*loglikelihood * density(x), wrt x
	 */
	int i, k, np = GMRFLib_faster_integration_np;
	double low, dx, dxi, *xp = NULL, *xpi = NULL, *dens = NULL, *loglik = NULL, *work = NULL, integral = 0.0, w[2] = { 4.0, 2.0 }, integral_one;

	GMRFLib_ASSERT_RETVAL(np > 3, GMRFLib_ESNH, 0.0);

	work = Calloc(4 * np, double);

	xp = work;
	xpi = work + np;
	dens = work + 2 * np;
	loglik = work + 3 * np;

	dxi = (density->x_max - density->x_min) / (np - 1.0);
	low = GMRFLib_density_std2user(density->x_min, density);
	dx = (GMRFLib_density_std2user(density->x_max, density) - low) / (np - 1.0);

	xp[0] = low;
	xpi[0] = density->x_min;
	for (i = 1; i < np; i++) {
		xp[i] = xp[0] + i * dx;
		xpi[i] = xpi[0] + i * dxi;
	}
	GMRFLib_evaluate_ndensity(dens, xpi, np, density);
	loglFunc(loglik, xp, np, idx, x_vec, loglFunc_arg);

	integral = loglik[0] * dens[0] + loglik[np - 1] * dens[np - 1];
	integral_one = dens[0] + dens[np - 1];
	for (i = 1, k = 0; i < np - 1; i++, k = (k + 1) % 2) {
		integral += w[k] * loglik[i] * dens[i];
		integral_one += w[k] * dens[i];
	}
	integral = -2.0 * integral / integral_one;

	Free(work);
	return integral;
}

/**
 *   \brief Free an \c GMRFLib_ai_cpo_tp -object created by \c GMRFLib_INLA()
 */
int GMRFLib_ai_cpo_free(GMRFLib_ai_cpo_tp * cpo)
{
	int i;

	if (!cpo) {
		return GMRFLib_SUCCESS;
	}
	if (cpo->value) {
		for (i = 0; i < cpo->n; i++) {
			if (cpo->value[i]) {
				Free(cpo->value[i]);
			}
			if (cpo->pit_value[i]) {
				Free(cpo->pit_value[i]);
			}
		}
		Free(cpo->value);
		Free(cpo->pit_value);
	}
	Free(cpo);

	return GMRFLib_SUCCESS;
}

int GMRFLib_ai_add_Qinv_to_ai_store(GMRFLib_ai_store_tp * ai_store)
{
	if (!ai_store || !(ai_store->problem)) {
		return GMRFLib_SUCCESS;
	}
	if (!ai_store->problem->sub_inverse) {
		int i, n;

		GMRFLib_Qinv(ai_store->problem, GMRFLib_QINV_NEIGB);
		Free(ai_store->stdev);
		n = ai_store->problem->n;
		ai_store->stdev = Calloc(n, double);
		for (i = 0; i < n; i++) {
			double *var = GMRFLib_Qinv_get(ai_store->problem, i, i);

			if (var) {
				ai_store->stdev[i] = sqrt(*var);
			}
		}
	}
	return GMRFLib_SUCCESS;
}
int GMRFLib_ai_adjust_integration_weights(double *adj_weights, double *weights, double **izs, int n, int nhyper, double dz)
{
	/*
	 * a simple adjustment of the weights: compare f_midpoints*dz with the integral of the k-dim iid gaussian in the dz box.
	 *
	 * this is for the case: ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_GRID
	 */

	int i, k, fail = 0;
	double w, ww, correction, x, f = 1 / sqrt(2.0 * M_PI);

	for (i = 0; i < n && !fail; i++) {
		correction = 1.0;
		for (k = 0; k < nhyper; k++) {
			x = izs[i][k];
			correction *= (gsl_cdf_ugaussian_P(x + dz / 2.0) - gsl_cdf_ugaussian_P(x - dz / 2.0)) / (dz * f * exp(-0.5 * SQR(x)));
		}
		adj_weights[i] = correction * weights[i];
		if (ISZERO(correction)) {
			fail = 1;
		}
	}

	if (!fail) {
		/*
		 * scale the adjusted weigts so them sum to the same as the origial weights 
		 */
		w = ww = 0.0;
		for (i = 0; i < n; i++) {
			w += weights[i];
			ww += adj_weights[i];
		}
		w = w / ww;
		for (i = 0; i < n; i++) {
			adj_weights[i] *= w;
		}
	} else {
		w = 0.0;
		for (i = 0; i < n; i++) {
			w += weights[i];
		}
		for (i = 0; i < n; i++) {
			adj_weights[i] = weights[i] / w;
		}
	}
	return GMRFLib_SUCCESS;
}
int GMRFLib_ai_marginal_for_one_hyperparamter(GMRFLib_density_tp ** density, int idx, int nhyper, int hyper_count, double *hyper_z,
					      double *hyper_ldens, double *theta_mode, gsl_vector * eigen_values,
					      gsl_matrix * eigen_vectors, double *std_stdev_theta, double dz,
					      double *stdev_corr_pos, double *stdev_corr_neg, GMRFLib_ai_interpolator_tp interpolator)
{
#define NEXTRA 11
	int i, j;
	double *points = NULL, *ldens_values, *theta_max, *theta_min, sd;
	double extra_points[NEXTRA] = { -3.0, -2.0, -1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 2.0, 3.0 };
	int npoints;

	GMRFLib_ENTER_ROUTINE;

	// 
	// Do not use this option, better off just doing what's done below. This is especially bad if the int_strategy = CCD.
	// 
	// if (idx == 0 && nhyper == 1) 
	// sd = std_stdev_theta[idx];
	// GMRFLib_density_create(density, GMRFLib_DENSITY_TYPE_SCGAUSSIAN, hyper_count, hyper_z, hyper_ldens, theta_mode[idx], sd, GMRFLib_TRUE);
	// 

	if (interpolator == GMRFLib_AI_INTERPOLATOR_GRIDSUM) {
		/*
		 * this require that hessian_force_diagonal == TRUE, but this is verified in the call, unless the use wants it... 
		 */
		double *dens = NULL, *theta_tmp = NULL, *work = NULL, *npoints_j = NULL;

		npoints = 0;
		work = Calloc(5 * hyper_count + nhyper, double);

		points = work;
		dens = work + hyper_count;
		ldens_values = work + 2 * hyper_count;
		theta_tmp = work + 3 * hyper_count;
		npoints_j = work + 4 * hyper_count;

		for (i = 0; i < hyper_count; i++) {
			GMRFLib_ai_z2theta(theta_tmp, nhyper, theta_mode, &(hyper_z[i * nhyper]), eigen_values, eigen_vectors);
			j = GMRFLib_which(theta_tmp[idx], points, npoints);
			if (j >= 0) {
				/*
				 * point we have already 
				 */
				dens[j] += exp(hyper_ldens[i]);
				npoints_j[j]++;
			} else {
				/*
				 * new point 
				 */
				points[npoints] = theta_tmp[idx];
				dens[npoints] = exp(hyper_ldens[i]);
				npoints_j[npoints] = 1.0;
				npoints++;
			}
		}

		if (nhyper > 1) {
			int new_npoints = 0;

			for (i = j = 0; i < npoints; i++) {
				if (npoints_j[i] > 1.0) {
					points[j] = points[i];
					dens[j] = dens[i];

					new_npoints++;
					j++;
				}
			}
			npoints = new_npoints;
		}

		sd = std_stdev_theta[idx];

		for (i = 0; i < npoints; i++) {
			ldens_values[i] = log(dens[i]);
			points[i] = (points[i] - theta_mode[idx]) / sd;
			// printf("points[%1d] = %f ldens %f \n", i, points[i], ldens_values[i]);
		}
		GMRFLib_density_create(density, GMRFLib_DENSITY_TYPE_SCGAUSSIAN, npoints, points, ldens_values,
				       theta_mode[idx], std_stdev_theta[idx], GMRFLib_TRUE);
		Free(work);
	} else {
		/*
		 * compute the marginal for the idx'th hyperparameter 
		 */
		npoints = 31;
		points = Calloc(npoints + NEXTRA, double);
		ldens_values = Calloc(npoints + NEXTRA, double);

		sd = std_stdev_theta[idx];
		theta_max = Calloc(nhyper - 1, double);
		theta_min = Calloc(nhyper - 1, double);

		double *theta_tmp = Calloc(nhyper, double);
		double *theta_min_all = Calloc(nhyper, double);
		double *theta_max_all = Calloc(nhyper, double);

		for (i = 0; i < hyper_count; i++) {
			GMRFLib_ai_z2theta(theta_tmp, nhyper, theta_mode, &(hyper_z[i * nhyper]), eigen_values, eigen_vectors);
			if (i == 0) {
				for (j = 0; j < nhyper; j++) {
					theta_max_all[j] = theta_min_all[j] = theta_tmp[j];
				}
			} else {
				for (j = 0; j < nhyper; j++) {
					theta_max_all[j] = DMAX(theta_max_all[j], theta_tmp[j]);
					theta_min_all[j] = DMIN(theta_min_all[j], theta_tmp[j]);
				}
			}
		}
		for (i = j = 0; i < nhyper; i++) {
			if (i != idx) {
				theta_min[j] = theta_min_all[i];
				theta_max[j] = theta_max_all[i];
				j++;
			}
		}
		GMRFLib_ai_integrator_arg_tp *arg = Calloc(1, GMRFLib_ai_integrator_arg_tp);

		arg->nhyper = nhyper;
		arg->idx = idx;
		arg->hyper_count = hyper_count;
		arg->hyper_z = hyper_z;
		arg->hyper_ldens = hyper_ldens;
		arg->theta_mode = theta_mode;
		arg->eigen_values = eigen_values;
		arg->eigen_vectors = eigen_vectors;
		arg->z = Calloc(nhyper, double);
		arg->theta = Calloc(nhyper, double);

		arg->stdev_corr_pos = stdev_corr_pos;
		arg->stdev_corr_neg = stdev_corr_neg;
		arg->dz = dz;
		arg->interpolator = interpolator;

		/*
		 * we need to bound the maximum function evaluations, otherwise it can just go on forever, especially for _linear
		 * and _quadratic interpolation. seems like they produce to `rough' integrands... 
		 */
		unsigned int max_eval = (unsigned int) gsl_pow_int(400.0, ITRUNCATE(nhyper - 1, 1, 3));	/* 0 for none */
		double abs_err = 0.0001, rel_err = 0.0001, value, err;
		int retval;

		for (i = 0; i < npoints; i++) {
			arg->theta_fixed = theta_min_all[idx] + i * (theta_max_all[idx] - theta_min_all[idx]) / (npoints - 1.0);
			points[i] = (arg->theta_fixed - theta_mode[idx]) / sd;
			retval =
			    adapt_integrate(GMRFLib_ai_integrator_func, arg, (unsigned int) nhyper - 1, (const double *) theta_min,
					    (const double *) theta_max, max_eval, abs_err, rel_err, &value, &err);
			value = DMAX(DBL_MIN, value);
			ldens_values[i] = log(value);
			if (retval) {
				fprintf(stderr, "\n\tGMRFLib_ai_marginal_for_one_hyperparamter: warning:\n");
				fprintf(stderr, "\t\tMaximum number of function evaluations is reached\n");
				fprintf(stderr, "\t\tPart 1, i=%1d, point=%g\n", i, points[i]);
			}
		}
		for (i = 0; i < NEXTRA; i++) {
			arg->theta_fixed = theta_mode[idx] + extra_points[i] * sd;
			points[i + npoints] = extra_points[i];
			retval =
			    adapt_integrate(GMRFLib_ai_integrator_func, arg, (unsigned int) nhyper - 1, (const double *) theta_min,
					    (const double *) theta_max, max_eval, abs_err, rel_err, &value, &err);
			value = DMAX(DBL_MIN, value);
			ldens_values[i + npoints] = log(value);
			if (retval) {
				fprintf(stderr, "\n\tGMRFLib_ai_marginal_for_one_hyperparamter: warning:\n");
				fprintf(stderr, "\t\tMaximum number of function evaluations is reached\n");
				fprintf(stderr, "\t\tPart 2, i=%1d, point=%g\n", i, points[i]);
			}
		}

		GMRFLib_density_create(density, GMRFLib_DENSITY_TYPE_SCGAUSSIAN, npoints + NEXTRA, points, ldens_values, theta_mode[idx], sd, GMRFLib_TRUE);

		Free(arg->z);
		Free(arg->theta);
		Free(arg);
		Free(theta_max_all);
		Free(theta_min_all);
		Free(theta_tmp);
		Free(theta_min);
		Free(theta_max);
		Free(ldens_values);
		Free(points);
	}

#undef NEXTRA
	GMRFLib_LEAVE_ROUTINE;

	return GMRFLib_SUCCESS;
}
double GMRFLib_ai_integrator_func(unsigned ndim, const double *x, void *arg)
{
	/*
	 * x is theta but for ndim = a->nhyper - 1. (this function also works for ndim == a->nhyper).
	 */

	GMRFLib_ai_integrator_arg_tp *a = (GMRFLib_ai_integrator_arg_tp *) arg;
	double val;
	int i, j;

	if (a->idx >= 0) {				       /* ie, ndim == nhyper-1 */
		for (i = j = 0; i < a->nhyper; i++) {
			if (i != a->idx) {
				a->theta[i] = x[j];
				j++;
			}
		}
		a->theta[a->idx] = a->theta_fixed;
	} else {
		memcpy(a->theta, x, a->nhyper * sizeof(double));	/* ndim == nhypera */
	}
	GMRFLib_ai_theta2z(a->z, a->nhyper, a->theta_mode, a->theta, a->eigen_values, a->eigen_vectors);
	switch (a->interpolator) {
	case GMRFLib_AI_INTERPOLATOR_NEAREST:
		val = GMRFLib_interpolator_nearest(a->nhyper, a->hyper_count, a->z, a->hyper_z, a->hyper_ldens, (void *) &(a->dz));
		break;
	case GMRFLib_AI_INTERPOLATOR_LINEAR:
		val = GMRFLib_interpolator_linear(a->nhyper, a->hyper_count, a->z, a->hyper_z, a->hyper_ldens, (void *) &(a->dz));
		break;
	case GMRFLib_AI_INTERPOLATOR_QUADRATIC:
		val = GMRFLib_interpolator_quadratic(a->nhyper, a->hyper_count, a->z, a->hyper_z, a->hyper_ldens, (void *) &(a->dz));
		break;
	case GMRFLib_AI_INTERPOLATOR_WEIGHTED_DISTANCE:
		val = GMRFLib_interpolator_wdistance(a->nhyper, a->hyper_count, a->z, a->hyper_z, a->hyper_ldens, (void *) &(a->dz));
		break;
	case GMRFLib_AI_INTERPOLATOR_CCD:
		val = GMRFLib_interpolator_ccd(a->nhyper, a->hyper_count, a->z, a->hyper_z, a->hyper_ldens, (void *) a);
		break;
	default:
		GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
	}

	if (0) {
		for (i = 0; i < (int) ndim; i++) {
			printf(" %.10g %.10g", a->z[0], a->z[1]);
		}
		printf(" Lin %.10g Quad %.10g Wdist %.10g",
		       GMRFLib_interpolator_linear(a->nhyper, a->hyper_count, a->z, a->hyper_z, a->hyper_ldens, (void *) &(a->dz)),
		       GMRFLib_interpolator_quadratic(a->nhyper, a->hyper_count, a->z, a->hyper_z, a->hyper_ldens,
						      (void *) &(a->dz)), GMRFLib_interpolator_wdistance(a->nhyper, a->hyper_count,
													 a->z, a->hyper_z, a->hyper_ldens, (void *) &(a->dz)));
		if (a->stdev_corr_pos && a->stdev_corr_neg) {
			printf(" CCD %.10g\n", GMRFLib_interpolator_ccd(a->nhyper, a->hyper_count, a->z, a->hyper_z, a->hyper_ldens, (void *) a));
		} else {
			printf("\n");
		}
	}

	return exp(val);
}
double GMRFLib_interpolator_distance2(int ndim, double *x, double *xx)
{
	/*
	 * return the squared Eucledian distance between x and xx 
	 */
	int i;
	double d = 0.0;

	for (i = 0; i < ndim; i++) {
		d += SQR(x[i] - xx[i]);
	}

	return d;
}
double GMRFLib_interpolator_distance(int ndim, double *x, double *xx)
{
	return sqrt(GMRFLib_interpolator_distance2(ndim, x, xx));
}
double GMRFLib_interpolator_nearest(int ndim, int nobs, double *x, double *xobs, double *yobs, void *arg)
{
	/*
	 * Just use the nearest point
	 */
	int i, imin = 0;
	double dist = 0.0, dtmp, value;

	/*
	 * compute the distances and find the indices for the smallest distances. 
	 */
	for (i = 0; i < nobs; i++) {
		dtmp = GMRFLib_interpolator_distance2(ndim, x, &(xobs[i * ndim]));
		if (dtmp < dist || i == 0) {
			dist = dtmp;
			imin = i;
		}
	}
	value = yobs[imin];

	return value;
}
double GMRFLib_interpolator_linear(int ndim, int nobs, double *x, double *xobs, double *yobs, void *arg)
{
	/*
	 * Compute the interpolated value at x for nobs observations: xobs, yobs. dimension of xobs is ndim*nobs and stored
	 * C-wise.  This routine just implements linear interpolation from the ndim+1 nearest points.
	 */
	size_t i, ii, j;
	double value, ymin, ymax;

	double *dd = Calloc(nobs, double);
	size_t *idxs = Calloc(nobs, size_t);
	size_t m = (size_t) ndim + 1;
	gsl_matrix *A = gsl_matrix_calloc(m, m);
	gsl_vector *b = gsl_vector_calloc(m);
	gsl_vector *sol = gsl_vector_calloc(m);
	gsl_permutation *p = gsl_permutation_alloc(m);

	/*
	 * compute the distances and find the indices for the smallest `m' distances. 
	 */
	for (i = 0; i < (size_t) nobs; i++) {
		dd[i] = GMRFLib_interpolator_distance2(ndim, x, &(xobs[i * ndim]));
		idxs[i] = i;
	}
	gsl_sort_smallest_index(idxs, m, dd, 1, (size_t) nobs);

	/*
	 * make the linear interpolator 
	 */
	for (i = 0; i < m; i++) {
		ii = idxs[i];
		for (j = 0; j < m - 1; j++) {
			gsl_matrix_set(A, i, j, xobs[ii * ndim + j]);
		}
		gsl_matrix_set(A, i, m - 1, 1.0);
		/*
		 * make sure the matrix is not singular. but then it strictly not an interpolator, but... 
		 */
		gsl_matrix_set(A, i, i, gsl_matrix_get(A, i, i) + FLT_EPSILON);
		gsl_vector_set(b, i, yobs[ii]);
	}
	int s;

	gsl_linalg_LU_decomp(A, p, &s);
	gsl_linalg_LU_solve(A, p, b, sol);

	/*
	 * compute the interpolated values 
	 */
	for (i = 0, value = 0.0; i < m - 1; i++) {
		value += gsl_vector_get(sol, i) * x[i];
	}
	value += gsl_vector_get(sol, m - 1);

	/*
	 * prevent odd cases where the point `x' is not in the 'triangle'... 
	 */
	ymax = gsl_vector_max(b);
	ymin = gsl_vector_min(b);
	value = TRUNCATE(value, ymin, ymax);

	/*
	 * cleanup 
	 */
	gsl_permutation_free(p);
	gsl_matrix_free(A);
	gsl_vector_free(b);
	gsl_vector_free(sol);
	Free(idxs);
	Free(dd);

	return value;
}
double GMRFLib_interpolator_quadratic(int ndim, int nobs, double *x, double *xobs, double *yobs, void *arg)
{
	/*
	 * Compute the interpolated value at x for nobs observations: xobs, yobs. dimension of xobs is ndim*nobs and stored
	 * C-wise.  This routine just implements quadratic interpolation from the ndim+1 nearest points, assuming the curvature
	 * is 1, as for the standard-normal.
	 */

	int s;
	size_t i, ii, j;
	double value, correction, dz, ymax = 0, ymin = 0;

	double *dd = Calloc(nobs, double);
	size_t *idxs = Calloc(nobs, size_t);
	size_t m = (size_t) ndim + 1;
	gsl_matrix *A = gsl_matrix_calloc(m, m);
	gsl_vector *b = gsl_vector_calloc(m);
	gsl_vector *sol = gsl_vector_calloc(m);
	gsl_permutation *p = gsl_permutation_alloc(m);

	dz = (arg ? *((double *) arg) : 1.0);
	/*
	 * compute the distances and find the indices for the smallest `m' distances. 
	 */
	for (i = 0; i < (size_t) nobs; i++) {
		dd[i] = GMRFLib_interpolator_distance2(ndim, x, &(xobs[i * ndim]));
		idxs[i] = i;
	}
	gsl_sort_smallest_index(idxs, m, dd, 1, (size_t) nobs);

	/*
	 * make the linear interpolator 
	 */
	for (i = 0; i < m; i++) {
		ii = idxs[i];
		correction = 0.0;
		for (j = 0; j < m - 1; j++) {
			gsl_matrix_set(A, i, j, xobs[ii * ndim + j]);
			correction += 0.5 * SQR(xobs[ii * ndim + j] * dz);
		}
		gsl_matrix_set(A, i, m - 1, 1.0);
		/*
		 * make sure the matrix is not singular. but then it strictly not an interpolator, but... 
		 */
		gsl_matrix_set(A, i, i, gsl_matrix_get(A, i, i) + FLT_EPSILON);
		gsl_vector_set(b, i, yobs[ii] + correction);

		if (i == 0) {
			ymax = ymin = yobs[ii];
		} else {
			ymax = DMAX(ymax, yobs[ii]);
			ymin = DMIN(ymin, yobs[ii]);
		}
	}

	gsl_linalg_LU_decomp(A, p, &s);
	gsl_linalg_LU_solve(A, p, b, sol);

	/*
	 * compute the interpolated values 
	 */
	for (i = 0, value = 0.0; i < m - 1; i++) {
		value += gsl_vector_get(sol, i) * x[i] - 0.5 * SQR(x[i] * dz);
	}
	value += gsl_vector_get(sol, m - 1);
	value = TRUNCATE(value, ymin, ymax);

	/*
	 * cleanup 
	 */
	gsl_permutation_free(p);
	gsl_matrix_free(A);
	gsl_vector_free(b);
	gsl_vector_free(sol);
	Free(idxs);
	Free(dd);

	return value;
}
double GMRFLib_interpolator_wdistance(int ndim, int nobs, double *x, double *xobs, double *yobs, void *arg)
{
	/*
	 * Compute the interpolated value at x for nobs observations: xobs, yobs. dimension of xobs is ndim*nobs and stored
	 * C-wise.  This routine just implements the easiest choice: the Shephard method.
	 * 
	 * WARNING: This functions is tuned for the problem where xobs are standarised like the standard-normal.
	 * 
	 * WARNING: is it not at all good in extrapolating values outside the 'region of interest'...
	 */

	int i;
	double w, wsum = 0.0, value = 0.0, d2, dz, idz2;

	dz = (arg ? *((double *) arg) : 1.0);
	idz2 = 1.0 / SQR(dz);

	for (i = 0; i < nobs; i++) {
		d2 = GMRFLib_interpolator_distance2(ndim, x, &(xobs[i * ndim]));

		/*
		 * this produce ``spikes'' not good for numerical integration.... so I turn this off for the moment 
		 */
		if (d2 < 0.0) {
			return yobs[i];
		}

		/*
		 * use either weights like exp(-d) or exp(-d^2). I prefer the `-d' as the 'd^2' seem to produce to smooth
		 * interpolators.
		 */
		// w = exp(-1.38629 * d2 * idz2);
		// w = exp(-1.38629 * sqrt(d2 * idz2));
		w = exp(-2.77258 * sqrt(d2 * idz2));

		wsum += w;
		value += w * yobs[i];
	}
	return value / wsum;
}
double GMRFLib_interpolator_ccd(int ndim, int nobs, double *x, double *xobs, double *yobs, void *arg)
{
	/*
	 * This is rather special; we use explicitely that the posterior is approximately Gaussian, but we use the corrected
	 * stdev's for positive and negative. (If they are not computed, then assume the corrected stdevs are all 1.)
	 */

	int i;
	double value = 0.0, sd;
	GMRFLib_ai_integrator_arg_tp *a = (GMRFLib_ai_integrator_arg_tp *) arg;

	if (a->stdev_corr_pos && a->stdev_corr_neg) {
		for (i = 0; i < ndim; i++) {
			sd = (x[i] > 0 ? a->stdev_corr_pos[i] : a->stdev_corr_neg[i]);
			value += -0.5 * SQR(x[i] / sd);
		}
	} else {
		for (i = 0; i < ndim; i++) {
			value += -0.5 * SQR(x[i]);
		}
	}

	return value;
}
GMRFLib_sizeof_tp GMRFLib_sizeof_ai_store(GMRFLib_ai_store_tp * ai_store)
{
	/*
	 * return, approximately, the size in bytes of ai_store 
	 */
	if (!ai_store)
		return 0;

	GMRFLib_sizeof_tp siz = 0;
	int n = (ai_store->problem ? ai_store->problem->n : 0);

	siz += sizeof(GMRFLib_ai_store_tp);
	siz += GMRFLib_sizeof_store(ai_store->store);
	siz += sizeof(double) * n;
	siz += GMRFLib_sizeof_problem(ai_store->problem);
	siz += sizeof(double) * n * 6;

	return siz;
}
GMRFLib_ai_store_tp *GMRFLib_duplicate_ai_store(GMRFLib_ai_store_tp * ai_store)
{
	/*
	 * duplicate AI_STORE 
	 */

#define DUPLICATE(name, len, tp)   if (1) {				\
		if (ai_store->name && len){				\
			new_ai_store->name = Calloc(len, tp);		\
			memcpy(new_ai_store->name, ai_store->name, len*sizeof(tp)); \
		} else {						\
			new_ai_store->name = NULL;			\
		}							\
	}
#define COPY(name) new_ai_store->name = ai_store->name

	GMRFLib_ENTER_ROUTINE;
	if (!ai_store) {
		GMRFLib_LEAVE_ROUTINE;
		return NULL;
	}

	GMRFLib_ai_store_tp *new_ai_store = Calloc(1, GMRFLib_ai_store_tp);
	int n = (ai_store->problem ? ai_store->problem->n : 0);
	int nd = ai_store->nd;

	new_ai_store->store = GMRFLib_duplicate_store(ai_store->store);
	DUPLICATE(mode, n, double);

	new_ai_store->problem = GMRFLib_duplicate_problem(ai_store->problem);

	COPY(nidx);
	COPY(neff);
	COPY(nd);

	DUPLICATE(bb, n, double);
	DUPLICATE(cc, n, double);
	DUPLICATE(stdev, n, double);
	DUPLICATE(correction_term, n, double);
	DUPLICATE(derivative3, n, double);
	DUPLICATE(correction_idx, n, int);
	DUPLICATE(d_idx, nd, int);

	GMRFLib_LEAVE_ROUTINE;
	return new_ai_store;

#undef DUPLICATE
#undef COPY
}
int GMRFLib_ai_pool_init(GMRFLib_ai_pool_tp ** pool, GMRFLib_ai_param_tp * ai_par, int nhyper)
{
	size_t i, j, k, debug = 0;
	int len, half_len;
	GMRFLib_ai_pool_tp *p;
	int *iz, *izz;

	GMRFLib_ASSERT(nhyper > 0, GMRFLib_EPARAMETER);

	*pool = Calloc(1, GMRFLib_ai_pool_tp);
	p = *pool;					       /* shorter name... */
	p->nhyper = nhyper;
	p->diff_log_dens = ai_par->diff_log_dens;
	p->all_out = 0;
	half_len = (int) ((sqrt(2.0 * p->diff_log_dens) + 1.0) / ai_par->dz + 2.0);
	len = 2 * half_len + 1;
	p->nconfig = (size_t) pow((double) len, (double) p->nhyper);
	p->configurations = Calloc((size_t) (p->nconfig * p->nhyper), GMRFLib_int8);
	p->idx_mapping = Calloc(p->nconfig, size_t);
	p->out = Calloc(p->nconfig, char);
	p->idx_next = 0;
	iz = Calloc(p->nhyper, int);
	izz = Calloc(p->nhyper, int);

	/*
	 * iz[i] goes from 0... len-1 izz[i] goes from -half_len ... half_len 
	 */
	for (i = k = 0; i < p->nconfig; i++) {

		for (j = 0; j < p->nhyper; j++) {
			izz[j] = (iz[j] <= half_len ? iz[j] : half_len - iz[j]);
		}
		if (debug) {
			printf("configuration %zu = [ ", i);
			for (j = 0; j < p->nhyper; j++) {
				printf("  %1d", izz[j]);
			}
			printf(" ]\n");
		}
		for (j = 0; j < p->nhyper; j++) {
			p->configurations[k + j] = (GMRFLib_int8) izz[j];
		}
		k += p->nhyper;

		int jj;
		for (jj = (int) p->nhyper - 1; jj >= 0; jj--) {
			if ((iz[jj] = (iz[jj] + 1) % len)) {
				break;
			}
		}
	}

	/*
	 * now we need to ``sort'' the configurations... recall to set pool_hyper which is required.
	 */
	pool_nhyper = (int) p->nhyper;
	if (omp_get_max_threads() > 1) {
		qsort(p->configurations, p->nconfig, p->nhyper * sizeof(GMRFLib_int8), GMRFLib_pool_cmp);
	} else {
		/*
		 * alternative sorting: _pool_cmp1: seems like _pool_cmp runs faster (better wrt 'reject') 
		 */
		qsort(p->configurations, p->nconfig, p->nhyper * sizeof(GMRFLib_int8), GMRFLib_pool_cmp);
	}
	pool_nhyper = -1;

	if (debug) {
		for (i = k = 0; i < p->nconfig; i++) {
			printf("sorted configuration %zu = [ ", i);
			for (j = 0; j < p->nhyper; j++) {
				printf("  %1d", p->configurations[k + j]);
			}
			printf(" ]\n");
			k += p->nhyper;
		}
	}

	Free(iz);
	Free(izz);

	return GMRFLib_SUCCESS;
}
int GMRFLib_pool_cmp(const void *a, const void *b)
{
	/*
	 * sort by Eucledian distance, otherwise, use the skip_config ordering.
	 * 
	 * Note that pool_hyper must be set properly. 
	 */
	const GMRFLib_int8 *ia, *ib;
	int i, larger, dist_a, dist_b;

	GMRFLib_ASSERT(pool_nhyper > 0, GMRFLib_ESNH);

	ia = (const GMRFLib_int8 *) a;
	ib = (const GMRFLib_int8 *) b;
	dist_a = dist_b = 0.0;
	for (i = 0; i < pool_nhyper; i++) {
		dist_a += ISQR((int) ia[i]);
		dist_b += ISQR((int) ib[i]);
	}
	if (!ISZERO(dist_a - dist_b)) {
		return (dist_a > dist_b ? 1 : -1);
	} else {
		larger = 1;
		for (i = 0; i < pool_nhyper; i++) {
			if (ia[i] > 0) {
				larger = larger && (ib[i] >= ia[i]);
			}
			if (ia[i] < 0) {
				larger = larger && (ib[i] <= ia[i]);
			}
		}
		return (larger ? -1 : 1);
	}
	return 0;
}
int GMRFLib_pool_cmp1(const void *a, const void *b)
{
	/*
	 * sort by skip_config ordering.
	 * 
	 * Note that pool_hyper must be set properly. 
	 */
	const GMRFLib_int8 *ia, *ib;
	int i, larger, eq;

	GMRFLib_ASSERT(pool_nhyper > 0, GMRFLib_ESNH);

	ia = (const GMRFLib_int8 *) a;
	ib = (const GMRFLib_int8 *) b;
	larger = 1;
	eq = 1;
	for (i = 0; i < pool_nhyper; i++) {
		eq = (eq && (ib[i] == ia[i]));
		if (ia[i]) {
			if (ia[i] > 0) {
				larger = larger && (ib[i] >= ia[i]);
			} else {
				larger = larger && (ib[i] <= ia[i]);
			}
		}
	}
	if (eq) {
		return 0;
	} else {
		return (larger ? -1 : 1);
	}
}
int GMRFLib_ai_pool_free(GMRFLib_ai_pool_tp * pool)
{
	if (pool) {
		Free(pool->configurations);
		Free(pool->idx_mapping);
		Free(pool->out);
		Free(pool);
	}
	return GMRFLib_SUCCESS;
}
int GMRFLib_ai_pool_get(GMRFLib_ai_pool_tp * pool, int *iz, size_t * idx)
{
	return GMRFLib_ai_pool_intern(pool, iz, idx, 0.0, GMRFLib_AI_POOL_GET);
}
int GMRFLib_ai_pool_set(GMRFLib_ai_pool_tp * pool, size_t idx, double logdens)
{
	return GMRFLib_ai_pool_intern(pool, NULL, &idx, logdens, GMRFLib_AI_POOL_SET);
}
int GMRFLib_ai_pool_intern(GMRFLib_ai_pool_tp * pool, int *iz, size_t * idx, double logdens, int action)
{
	int debug = 0;
	int retval = 0;

	GMRFLib_ASSERT(idx, GMRFLib_EPARAMETER);
#pragma omp critical					       /* yes, only one at the time */
	{
		size_t i, j, k;

		if (action == GMRFLib_AI_POOL_GET) {
			/*
			 * find the first free configuration. return the entry number in IDX and keeping the mapping between IDX
			 * and the internal index, in pool->idx_mapping.
			 */
			if (pool->all_out) {		       /* fast exit? */
				retval = !GMRFLib_SUCCESS;
			} else {
				int found = 0;
				for (i = 0; i < pool->nconfig && !found; i++) {
					if (!(pool->out[i])) {
						pool->out[i] = 1;
						pool->idx_mapping[pool->idx_next] = i;
						*idx = pool->idx_next;
						pool->idx_next++;
						for (j = 0; j < pool->nhyper; j++) {
							iz[j] = (int) pool->configurations[i * pool->nhyper + j];
						}
						if (debug) {
							printf("pool get idx %1zu i %zu: ", *idx, i);
							for (j = 0; j < pool->nhyper; j++) {
								printf(" %1d", iz[j]);
							}
							printf("\n");
						}
						found = 1;
						retval = GMRFLib_SUCCESS;
					}
				}
				if (!found) {
					pool->all_out = 1;     /* so we can make a fast return next time */
					retval = !GMRFLib_SUCCESS;
				}
			}
		} else if (action == GMRFLib_AI_POOL_SET) {
			if (-logdens > pool->diff_log_dens) {
				int *izz = NULL, *izz_local = NULL, larger;

				izz = Calloc(pool->nhyper, int);
				izz_local = Calloc(pool->nhyper, int);
				k = pool->idx_mapping[*idx];
				for (j = 0; j < pool->nhyper; j++) {
					izz[j] = (int) pool->configurations[k * pool->nhyper + j];
				}
				if (debug) {
					printf("pool set idx %1zu k %1zu: ", *idx, k);
					for (j = 0; j < pool->nhyper; j++) {
						printf(" %1d", izz[j]);
					}
					printf("\n");
				}


				for (i = 0; i < pool->nconfig; i++) {
					if (!pool->out[i]) {
						for (j = 0; j < pool->nhyper; j++) {
							izz_local[j] = (int) pool->configurations[i * pool->nhyper + j];
						}

						larger = 1;
						for (j = 0; j < pool->nhyper && larger; j++) {
							if (izz[j] > 0) {
								larger = larger && (izz_local[j] >= izz[j]);
							}
							if (izz[j] < 0) {
								larger = larger && (izz_local[j] <= izz[j]);
							}
						}
						if (larger) {
							pool->out[i] = 1;
							if (debug) {
								printf("\t\tpool set idx %1zu to OUT : ", i);
								for (j = 0; j < pool->nhyper; j++) {
									printf(" %1d", izz_local[j]);
								}
								printf("\n");
							}
						}
					}
				}
				Free(izz_local);
				Free(izz);
			}
			retval = GMRFLib_SUCCESS;
		}
	}
	return retval;
}

/*
  Example for manual
*/

/** \page ex_ai A worked out example in approximate inference using INLA
  
We consider here a model based on a first order auto-regressive latent field with unknown precision \f$ \tau \f$ and
autocorrelation parameter \f$ \phi \f$: \f[ x_t \mid x_1,\dots,x_{t-1},\phi,\kappa\sim {\mathcal N}(\phi x_{t-1},1/\tau).  \f]
Our observations are Poisson distributed \f[ y_t \mid x_t \sim \mbox{Po}(\exp(x_t)).  \f] We transform the hyperparameter
\f$\tau\f$ and \f$ \phi \f$ so that \f$\mbox{\boldmath$\theta$}\in\mathcal{R}^2\f$ \f[ \theta_1 = \log \kappa, \quad \theta_2 =
\mbox{logit}\frac{\phi + 1}{2}.  \f] Finally we choose priors for \f$\mbox{\boldmath$\theta$}\f$ \f[
\theta_1\sim\mbox{Gamma}(a,b), \quad \theta_2\sim{\mathcal N}(0,\tau_{\theta_2}).  \f]

- Our first task it to compute an approximation to the posterior marginal for the nodes in the hidden field given a fixed value
for the hyperparameter \f$\mbox{\boldmath$\theta$}\f$

\par Program code:

\verbinclude example-doxygen-approx-1.txt

- The next task it to compute an approximation to the posterior marginal for the nodes in the hidden field when we integrate out
the hyperparameters

\par Program code:

\verbinclude example-doxygen-approx-2.txt

- Finally, we compute an approximate density for the posterior of the hyperparameters \f[
\widetilde{\pi}(\mbox{\boldmath$\theta$} \mid \mbox{\boldmath$y$}) \propto
\frac{\pi(\mbox{\boldmath$x$},\mbox{\boldmath$\theta$}|\mbox{\boldmath$y$})}{
\widetilde{\pi}_{G}(\mbox{\boldmath$x$}|\mbox{\boldmath$\theta$}, \mbox{\boldmath$y$})}\Bigg|_{\mbox{\boldmath$x$} =
\mbox{\boldmath$x$}^{\star}(\mbox{\boldmath$\theta$})}\propto\mbox{ } \frac{\pi(\mbox{\boldmath$x$}\mid\mbox{\boldmath$\theta$})
\pi(\mbox{\boldmath$y$}\mid\mbox{\boldmath$x$})} {\widetilde{\pi}_{G}(\mbox{\boldmath$x$}|\mbox{\boldmath$\theta$},
\mbox{\boldmath$y$})}\pi(\mbox{\boldmath$\theta$}).\hspace{2cm} (INLA-1) \f] The routine \c GMRFLib_ai_marginal_hyperparam()
computes the first term in <b>(INLA-1)</b>, <em>except</em> for any term wich is constant with respect to \f$\mbox{\boldmath$x$}
\f$ but depends on \f$\mbox{\boldmath$\theta$} \f$. These terms, together with the prior for the hyperparameters
\f$\pi(\mbox{\boldmath$\theta$}) \f$ have to be provided through a \a GMRFLib_ai_log_extra_tp function.\n

\par Program code:

\verbinclude example-doxygen-approx-3.txt

*/

