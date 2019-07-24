
/* approx-inference.c
 * 
 * Copyright (C) 2006-2018 Havard Rue
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

#include <stdint.h>
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
double inla_compute_saturated_loglik(int, GMRFLib_logl_tp *, double *, void *);

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

	(*ai_par)->adapt_max = 5;
	(*ai_par)->adapt_len = 0;
	(*ai_par)->adapt_strategy = NULL;

	(*ai_par)->fast = GMRFLib_FALSE;		       /* compute conditional mode */
	(*ai_par)->fast = GMRFLib_TRUE;			       /* use mode = conditional mean */

	(*ai_par)->gaussian_data = GMRFLib_TRUE;
	(*ai_par)->gaussian_data = GMRFLib_FALSE;

	(*ai_par)->linear_correction = GMRFLib_AI_LINEAR_CORRECTION_CENTRAL_DIFFERENCE;
	(*ai_par)->linear_correction = GMRFLib_AI_LINEAR_CORRECTION_OFF;
	(*ai_par)->linear_correction = GMRFLib_AI_LINEAR_CORRECTION_FAST;	/* to match the ->fast mode above */

	/*
	 * none of these are used, but they are the defaults if the user wants improved approximations 
	 */
	(*ai_par)->n_points = 9;			       /* how many points to evaluate */
	(*ai_par)->step_len = GMRFLib_eps(0.25);	       /* If the derivaties has to be computed numerically */
	(*ai_par)->stencil = 5;				       /* number of points to use */
	(*ai_par)->cutoff = 0.0;			       /* the cutoff for the gradient in the (Gaussian) conditional mean */

	/*
	 * defaults for the integration itself 
	 */
	(*ai_par)->fp_log = NULL;
	(*ai_par)->fp_log = stdout;
	(*ai_par)->fp_hyperparam = NULL;
	(*ai_par)->int_strategy = GMRFLib_AI_INT_STRATEGY_EMPIRICAL_BAYES;
	(*ai_par)->int_strategy = GMRFLib_AI_INT_STRATEGY_CCD;
	(*ai_par)->int_strategy = GMRFLib_AI_INT_STRATEGY_USER;
	(*ai_par)->int_strategy = GMRFLib_AI_INT_STRATEGY_USER_STD;
	(*ai_par)->int_strategy = GMRFLib_AI_INT_STRATEGY_GRID;
	(*ai_par)->int_strategy = GMRFLib_AI_INT_STRATEGY_AUTO;
	(*ai_par)->int_design = NULL;
	(*ai_par)->f0 = 1.1;
	(*ai_par)->dz = 0.75;
	(*ai_par)->adjust_weights = GMRFLib_FALSE;
	(*ai_par)->adjust_weights = GMRFLib_TRUE;
	(*ai_par)->diff_log_dens = 4.0;
	(*ai_par)->skip_configurations = GMRFLib_FALSE;
	(*ai_par)->skip_configurations = GMRFLib_TRUE;

	/*
	 * use forward differences for the gradient and central differences for the (final) hessian 
	 */
	(*ai_par)->gradient_forward_finite_difference = GMRFLib_TRUE;	/* use forward difference */
	(*ai_par)->gradient_forward_finite_difference = GMRFLib_FALSE;	/* use central difference */
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
	(*ai_par)->interpolator = GMRFLib_AI_INTERPOLATOR_CCD_INTEGRATE;
	(*ai_par)->interpolator = GMRFLib_AI_INTERPOLATOR_CCD;
	(*ai_par)->interpolator = GMRFLib_AI_INTERPOLATOR_WEIGHTED_DISTANCE;
	(*ai_par)->interpolator = GMRFLib_AI_INTERPOLATOR_AUTO;	/* automatic choice */

	(*ai_par)->optimiser = GMRFLib_AI_OPTIMISER_GSL;
	(*ai_par)->restart = 0;
	(*ai_par)->gsl_tol = 0.1;
	(*ai_par)->gsl_epsg = 0.005;
	(*ai_par)->gsl_epsf = pow(0.005, 1.5);		       /* this is the default relationship used in R-INLA */
	(*ai_par)->gsl_epsx = 0.005;
	(*ai_par)->gsl_step_size = 1.0;
	(*ai_par)->mode_known = 0;

	/*
	 * parameters for the Gaussian approximations 
	 */
	(*ai_par)->optpar_abserr_func = 0.0005;
	(*ai_par)->optpar_abserr_step = 0.0005;
	(*ai_par)->optpar_fp = NULL;
	(*ai_par)->optpar_nr_step_factor = 1.0;

	(*ai_par)->cpo_req_diff_logdens = 3.0;

	(*ai_par)->stupid_search_mode = GMRFLib_FALSE;
	(*ai_par)->stupid_search_mode = GMRFLib_TRUE;
	(*ai_par)->stupid_search_max_iter = 1000;
	(*ai_par)->stupid_search_factor = 1.01;

	(*ai_par)->cpo_manual = GMRFLib_FALSE;

	/*
	 * for numerical integration 
	 */
	(*ai_par)->numint_max_fn_eval = 10000;
	(*ai_par)->numint_rel_err = 1e-3;
	(*ai_par)->numint_abs_err = 1e-4;

	/*
	 * for numerical optimisation
	 */
	(*ai_par)->cmin = 0.0;

	/*
	 * default is no correction
	 */
	(*ai_par)->correct = NULL;
	(*ai_par)->correct_factor = 1.0;		       /* set but is default not used */
	(*ai_par)->correct_strategy = GMRFLib_AI_STRATEGY_MEANCORRECTED_GAUSSIAN;
	(*ai_par)->correct_verbose = GMRFLib_FALSE;

	return GMRFLib_SUCCESS;
}

int GMRFLib_ai_param_duplicate(GMRFLib_ai_param_tp ** ai_par_new, GMRFLib_ai_param_tp * ai_par)
{
	if (!ai_par) {
		*ai_par_new = NULL;
		return GMRFLib_SUCCESS;
	}

	*ai_par_new = Calloc(1, GMRFLib_ai_param_tp);
	memcpy((void *) *ai_par_new, (void *) ai_par, sizeof(GMRFLib_ai_param_tp));

	if (ai_par->adapt_strategy) {
		(*ai_par_new)->adapt_strategy = Calloc(ai_par->adapt_len, GMRFLib_ai_strategy_tp);
		memcpy((void *) (*ai_par_new)->adapt_strategy, (void *) ai_par->adapt_strategy, sizeof(GMRFLib_ai_strategy_tp) * ai_par->adapt_len);
	}

	return GMRFLib_SUCCESS;
}
int GMRFLib_ai_param_free(GMRFLib_ai_param_tp * ai_par)
{
	if (ai_par) {
		if (ai_par->adapt_len > 0 && ai_par->adapt_strategy) {
			Free(ai_par->adapt_strategy);
		}
		Free(ai_par);
	}

	return GMRFLib_SUCCESS;
}



/*!
  \brief Print the values of a \c GMRFLib_ai_param_tp -object

  \param[out] fp The *FILE  on which to print the output
  
  \param[in] ai_par The \c GMRFLib_ai_param_tp -object to be printed
*/
int GMRFLib_print_ai_param(FILE * fp, GMRFLib_ai_param_tp * ai_par)
{
	int show_expert_options = 1;

	if (!ai_par) {
		return GMRFLib_SUCCESS;
	}
	fp = (fp ? fp : stdout);

	fprintf(fp, "Contents of ai_param %p\n", (void *) ai_par);

	fprintf(fp, "\tOptimiser: %s\n", GMRFLib_AI_OPTIMISER_NAME(ai_par->optimiser));
	fprintf(fp, "\t\tOption for %s: tol  = %.6g\n", GMRFLib_AI_OPTIMISER_NAME(GMRFLib_AI_OPTIMISER_GSL), ai_par->gsl_tol);
	fprintf(fp, "\t\tOption for %s: step_size = %.6g\n", GMRFLib_AI_OPTIMISER_NAME(GMRFLib_AI_OPTIMISER_GSL), ai_par->gsl_step_size);
	fprintf(fp, "\t\tOption for %s: epsx = %.6g\n", GMRFLib_AI_OPTIMISER_NAME(GMRFLib_AI_OPTIMISER_GSL), ai_par->gsl_epsx);
	fprintf(fp, "\t\tOption for %s: epsf = %.6g\n", GMRFLib_AI_OPTIMISER_NAME(GMRFLib_AI_OPTIMISER_GSL), ai_par->gsl_epsf);
	fprintf(fp, "\t\tOption for %s: epsg = %.6g\n", GMRFLib_AI_OPTIMISER_NAME(GMRFLib_AI_OPTIMISER_GSL), ai_par->gsl_epsg);
	fprintf(fp, "\t\tRestart: %1d\n", ai_par->restart);
	fprintf(fp, "\t\tMode known: %s\n", (ai_par->mode_known ? "Yes" : "No"));

	fprintf(fp, "\tGaussian approximation:\n");
	fprintf(fp, "\t\ttolerance_func = %.6g\n", ai_par->optpar_abserr_func);
	fprintf(fp, "\t\ttolerance_step = %.6g\n", ai_par->optpar_abserr_step);
	fprintf(fp, "\t\toptpar_fp = %" PRIxPTR "\n", (uintptr_t) (ai_par->optpar_fp));
	fprintf(fp, "\t\toptpar_nr_step_factor = %.6g\n", ai_par->optpar_nr_step_factor);

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
	if (ai_par->strategy == GMRFLib_AI_STRATEGY_ADAPTIVE) {
		fprintf(fp, "Use an adaptive strategy (max=%1d)\n", ai_par->adapt_max);
	}

	fprintf(fp, "\tFast mode: \t%s\n", (ai_par->fast == GMRFLib_FALSE ? "Off" : "On"));

	fprintf(fp, "\tUse linear approximation to log(|Q +c|)? %s\n",
		(ai_par->linear_correction == GMRFLib_AI_LINEAR_CORRECTION_OFF ? "No" : "Yes"));
	if (ai_par->linear_correction != GMRFLib_AI_LINEAR_CORRECTION_OFF) {
		fprintf(fp, "\t\tMethod:\t ");
		if (ai_par->linear_correction == GMRFLib_AI_LINEAR_CORRECTION_CENTRAL_DIFFERENCE) {
			fprintf(fp, "Estimate the derivative using central difference\n");
		}
		if (ai_par->linear_correction == GMRFLib_AI_LINEAR_CORRECTION_FAST) {
			fprintf(fp, "Compute the derivative exact\n");
		}
	}
	fprintf(fp, "\tParameters for improved approximations\n");
	fprintf(fp, "\t\tNumber of points evaluate:\t %d\n", ai_par->n_points);
	fprintf(fp, "\t\tStep length to compute derivatives numerically:\t %g\n", ai_par->step_len);
	fprintf(fp, "\t\tStencil to compute derivatives numerically:\t %d\n", ai_par->stencil);
	fprintf(fp, "\t\tCutoff value to construct local neigborhood:\t %g\n", ai_par->cutoff);

	fprintf(fp, "\tLog calculations:\t %s\n", (ai_par->fp_log ? "On" : "Off"));
	fprintf(fp, "\tLog calculated marginal for the hyperparameters:\t %s\n", (ai_par->fp_hyperparam ? "On" : "Off"));

	fprintf(fp, "\tIntegration strategy:\t ");
	if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_AUTO) {
		fprintf(fp, "Automatic (GRID for dim(theta)=1 and 2 and otherwise CCD)\n");
	}
	if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_CCD) {
		fprintf(fp, "Use points from Central Composite Design (CCD)\n");
	}
	if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_GRID) {
		fprintf(fp, "Use adaptive grid-approach (GRID)\n");
	}
	if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_EMPIRICAL_BAYES) {
		fprintf(fp, "Use only the modal configuration (EMPIRICAL_BAYES)\n");
	}
	if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_USER) {
		fprintf(fp, "Use user-defined integration points\n");
	}
	if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_USER_STD) {
		fprintf(fp, "Use user-defined standardized integration points\n");
	}
	if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_USER_EXPERT) {
		fprintf(fp, "Use user-defined expert integration points and weights\n");
	}
	GMRFLib_print_design(fp, ai_par->int_design);

	fprintf(fp, "\t\tf0 (CCD only):\t %f\n", ai_par->f0);
	fprintf(fp, "\t\tdz (GRID only):\t %f\n", ai_par->dz);
	fprintf(fp, "\t\tAdjust weights (GRID only):\t %s\n", (ai_par->adjust_weights == GMRFLib_FALSE ? "Off" : "On"));
	fprintf(fp, "\t\tDifference in log-density limit (GRID only):\t %f\n", ai_par->diff_log_dens);
	fprintf(fp, "\t\tSkip configurations with (presumed) small density (GRID only):\t %s\n",
		(ai_par->skip_configurations == GMRFLib_FALSE ? "Off" : "On"));

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

	fprintf(fp, "\tStupid search mode:\n");
	fprintf(fp, "\t\tStatus     [%s]\n", (ai_par->stupid_search_mode ? "On" : "Off"));
	fprintf(fp, "\t\tMax iter   [%d]\n", ai_par->stupid_search_max_iter);
	fprintf(fp, "\t\tFactor     [%g]\n", ai_par->stupid_search_factor);

	fprintf(fp, "\tNumerical integration of hyperparameters:\n");
	fprintf(fp, "\t\tMaximum number of function evaluations [%1d]\n", ai_par->numint_max_fn_eval);
	fprintf(fp, "\t\tRelative error ....................... [%g]\n", ai_par->numint_rel_err);
	fprintf(fp, "\t\tAbsolute error ....................... [%g]\n", ai_par->numint_abs_err);

	fprintf(fp, "\tTo stabilise the numerical optimisation:\n");
	fprintf(fp, "\t\tMinimum value of the -Hessian [%g]\n", ai_par->cmin);

	if (show_expert_options) {
		/*
		 * expert options goes here 
		 */
		fprintf(fp, "\tCPO manual calculation[%s]\n", (ai_par->cpo_manual ? "Yes" : "No"));
	}

	if (ai_par->correct) {
		fprintf(fp, "\tLaplace-correction is Enabled with correction factor[%.4f]\n", ai_par->correct_factor);
		if (ai_par->correct_strategy == GMRFLib_AI_STRATEGY_MEANCORRECTED_GAUSSIAN ||
		    ai_par->correct_strategy == GMRFLib_AI_STRATEGY_MEANSKEWCORRECTED_GAUSSIAN)
			fprintf(fp, "\t\tstrategy = [simplified.laplace]\n");
		if (ai_par->correct_strategy == GMRFLib_AI_STRATEGY_FIT_SCGAUSSIAN)
			fprintf(fp, "\t\tstrategy = [laplace]\n");
		fprintf(fp, "\t\tverbose = [%s]\n", (ai_par->correct_verbose ? "TRUE" : "FALSE"));
	} else {
		fprintf(fp, "\tLaplace-correction is Disabled.\n");
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
	 * return the unnormalised log marginal density for the hyperparamers in `logdens'. note that we assume that 'b' and 'mean' are constants. when this is not
	 * the case, the missing term is added elsewhere. 
	 */

	double ldens;
	int n, free_ai_par = 0;

	/*
	 * this is a special option for _INLA(), so it works only when calling it the first time 
	 */
	static int nr_step_factor_first_time_only = 1;
#pragma omp threadprivate(nr_step_factor_first_time_only)

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
	optpar->stencil = ai_par->stencil;
	optpar->abserr_func = ai_par->optpar_abserr_func;
	optpar->abserr_step = ai_par->optpar_abserr_func;
	optpar->fp = ai_par->optpar_fp;

	if (ai_par->optpar_nr_step_factor < 0) {
		/*
		 *  then do this only the first time for each thread. This would improve initial values
		 */
		if (nr_step_factor_first_time_only) {
			optpar->nr_step_factor = -ai_par->optpar_nr_step_factor;
			nr_step_factor_first_time_only = 0;
		} else {
			optpar->nr_step_factor = 1.0;
		}
	} else {
		optpar->nr_step_factor = ai_par->optpar_nr_step_factor;
	}

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
	Free(ai_store->aa);
	Free(ai_store->bb);
	Free(ai_store->cc);
	ai_store->aa = Calloc(n, double);
	ai_store->bb = Calloc(n, double);
	ai_store->cc = Calloc(n, double);

	/*
	 * first compute the GMRF-approximation 
	 */
	if (ai_store->mode) {
		GMRFLib_EWRAP1(GMRFLib_init_GMRF_approximation_store__intern
			       (&problem, ai_store->mode, b, c, mean, d, loglFunc, loglFunc_arg, fixed_value, graph, Qfunc,
				Qfunc_arg, constr, optpar, blockpar, ai_store->store, ai_store->aa, ai_store->bb, ai_store->cc,
				ai_par->gaussian_data, ai_par->cmin, 0));
	} else {
		GMRFLib_EWRAP1(GMRFLib_init_GMRF_approximation_store__intern
			       (&problem, x, b, c, mean, d, loglFunc, loglFunc_arg, fixed_value, graph, Qfunc, Qfunc_arg,
				constr, optpar, blockpar, ai_store->store, ai_store->aa, ai_store->bb, ai_store->cc,
				ai_par->gaussian_data, ai_par->cmin, 0));
	}

	GMRFLib_ASSERT(problem, GMRFLib_EOPTNR);

	/*
	 * if store, then store the mode to use as the initial point at later calls 
	 */
	Free(ai_store->mode);
	ai_store->mode = Calloc(n, double);
	memcpy(ai_store->mode, problem->mean_constr, n * sizeof(double));

	if (mean == NULL && 1) {
		/*
		 * Here we use the joint expression and take advantage of that we have already evaluated the log-likelihood in the mode.
		 * 
		 * the expression below is not correct if mean != 0, so then we need a term -1/2 mu*(Q+C)mu as well. so therefore the if (mean==NULL) test. 
		 */
		memset(problem->sample, 0, n * sizeof(double));
		GMRFLib_evaluate(problem);

		double A = 0;
		int i;
		for (i = 0; i < n; i++) {
			A += ai_store->aa[i];
		}
		*logdens = A - problem->sub_logdens;
	} else {
		/*
		 * then evaluate it in the mode (ie the mean) to get the log-density 
		 */
		memcpy(problem->sample, problem->mean_constr, n * sizeof(double));
		GMRFLib_EWRAP1(GMRFLib_evaluate(problem));
		GMRFLib_ai_log_posterior(&ldens, problem->mean_constr, b, c, mean, d, loglFunc, loglFunc_arg, fixed_value, graph,
					 Qfunc, Qfunc_arg, constr);

		// printf("ai_marginal_hyper thread_id %d ldens %.12f sub_logdens %.12f\n", GMRFLib_thread_id, ldens,
		// problem->sub_logdens);

		*logdens = ldens - problem->sub_logdens;
	}

	/*
	 * store the GMRF-approximation in store 
	 */
	GMRFLib_free_problem(ai_store->problem);
	ai_store->problem = problem;

	if (ai_par->correct) {
		/*
		 * compute the correction to the LA
		 */
		int i, j, compute_n = 0, *compute_idx = NULL;
		for (i = 0; i < n; i++) {
			compute_n += (int) ai_par->correct[i];
		}
		assert(compute_n > 0);

		compute_idx = Calloc(compute_n, int);
		for (i = j = 0; i < n; i++) {
			if (ai_par->correct[i]) {
				compute_idx[j++] = i;
			}
		}
		assert(j == compute_n);

		GMRFLib_marginal_hidden_store_tp *marginal_hidden_store = Calloc(1, GMRFLib_marginal_hidden_store_tp);
		GMRFLib_density_tp **dens = Calloc(compute_n, GMRFLib_density_tp *);
		GMRFLib_ai_param_tp *ai_par_local = Calloc(1, GMRFLib_ai_param_tp);
		marginal_hidden_store->n = graph->n;
		marginal_hidden_store->subgraphs = Calloc(graph->n, GMRFLib_graph_tp *);

		memcpy(ai_par_local, ai_par, sizeof(GMRFLib_ai_param_tp));
		ai_par_local->strategy = ai_par->correct_strategy;
		GMRFLib_ai_add_Qinv_to_ai_store(ai_store);

		if (ai_par->correct_verbose) {
			printf("\tCorrect: Compute marginals for %d nodes\n", compute_n);
		}
		for (i = 0; i < compute_n; i++) {
			GMRFLib_ai_marginal_hidden(&dens[i],
						   NULL, compute_idx[i], x, b, c, mean, d,
						   loglFunc, loglFunc_arg, fixed_value, graph, Qfunc,
						   Qfunc_arg, constr, ai_par_local, ai_store, marginal_hidden_store);
		}
		Free(ai_par_local);

		GMRFLib_lc_tp **Alin = Calloc(compute_n, GMRFLib_lc_tp *);
		for (i = 0; i < compute_n; i++) {
			Alin[i] = Calloc(1, GMRFLib_lc_tp);
			Alin[i]->n = 1;
			Alin[i]->idx = Calloc(1, int);
			Alin[i]->idx[0] = compute_idx[i];
			Alin[i]->weight = Calloc(1, float);
			Alin[i]->weight[0] = 1.0;
			Alin[i]->tinfo = Calloc(ISQR(GMRFLib_MAX_THREADS), GMRFLib_lc_tinfo_tp);
			for (j = 0; j < ISQR(GMRFLib_MAX_THREADS); j++) {
				Alin[i]->tinfo[j].first_nonzero = -1;
				Alin[i]->tinfo[j].last_nonzero = -1;
				Alin[i]->tinfo[j].first_nonzero_mapped = -1;
				Alin[i]->tinfo[j].last_nonzero_mapped = -1;
			}
		}

		double *improved_mean = Calloc(n, double);     /* just with zeros as we do not care about the mean... */
		GMRFLib_density_tp **lin_dens = Calloc(compute_n, GMRFLib_density_tp *);
		double *cov = NULL;

		if (ai_par->correct_verbose) {
			printf("\tCorrect: Compute the covariance for the %d nodes\n", compute_n);
		}
		GMRFLib_ai_compute_lincomb(&lin_dens, &cov, compute_n, Alin, ai_store, improved_mean);
		if (ai_par->correct_verbose) {
			printf("\t\tCovariance matrix:\n");
			for (i = 0; i < compute_n; i++) {
				printf("\t\t");
				for (j = 0; j < compute_n; j++) {
					printf("%10.5f ", cov[i + j * compute_n]);
				}
				printf("\n");
			}
			printf("\t\tDifference in the mean:\n");
			for (i = 0; i < compute_n; i++) {
				printf("\t\t%10.5f\n", dens[i]->std_mean - dens[i]->user_mean);
			}
		}

		double corr = 0.0, *icov = cov;
		GMRFLib_comp_posdef_inverse(cov, compute_n);
		for (i = 0; i < compute_n; i++) {
			for (j = 0; j < compute_n; j++) {
				corr += (dens[i]->std_mean - dens[i]->user_mean)
				    * icov[i + j * compute_n]
				    * (dens[j]->std_mean - dens[j]->user_mean);
			}
		}

#define FUNCORR(_x)  (2.0/(1+exp(-2.0 * (_x))) -1.0)	       // makes the derivative in 0 eq to 1
		double upper = compute_n * ai_par->correct_factor;
		*logdens += 0.5 * upper * FUNCORR(corr / upper);
		if (ai_par->correct_verbose) {
			printf("\t\tCorrect: correction: raw = %.6f adjusted = %.6f\n", 0.5 * corr, 0.5 * upper * FUNCORR(corr / upper));
		}
#undef FUNCORR

		GMRFLib_free_marginal_hidden_store(marginal_hidden_store);
		Free(improved_mean);
		for (i = 0; i < compute_n; i++) {
			GMRFLib_free_density(lin_dens[i]);
			GMRFLib_free_density(dens[i]);
			Free(Alin[i]->idx);
			Free(Alin[i]->weight);
			Free(Alin[i]->tinfo);
		}
		Free(lin_dens);
		Free(dens);
		Free(Alin);
		Free(cov);
		Free(compute_idx);
	}


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

	run_with_omp = (GMRFLib_openmp->max_threads_inner > 1 ? 1 : 0);

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

						for (iii = 0; iii < nidx; iii++) {
							GMRFLib_thread_id = id;
							ii = idxs[iii];
							loglFunc(&logll, &x[ii], 1, ii, x, NULL, loglFunc_arg);
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
							loglFunc(&logll, &x[ii], 1, ii, x, NULL, loglFunc_arg);
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
								loglFunc(&logll, &x[ii], 1, ii, x, NULL, loglFunc_arg);
								tmp2 += d[ii] * logll;
							}
						}
					} else {
						for (ii = 0; ii < n; ii++) {
							if (d[ii]) {
								loglFunc(&logll, &x[ii], 1, ii, x, NULL, loglFunc_arg);
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
						loglFunc(&logll, &x[i], 1, i, x, NULL, loglFunc_arg);
						tmp += d[i] * logll;
					}
				}
			} else {
				for (i = 0; i < n; i++) {
					if (d[i]) {
						loglFunc(&logll, &x[i], 1, i, x, NULL, loglFunc_arg);
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
int GMRFLib_ai_log_posterior_restricted_OLD(double *logdens, double *x, double *x_mode, double *x_gradient, double delta,
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
		GMRFLib_EWRAP0(GMRFLib_ai_log_posterior
			       (logdens, x, b, c, mean, d, loglFunc, loglFunc_arg, fixed_value, graph, Qfunc, Qfunc_arg, constr));
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
					loglFunc(&logll, &x[i], 1, i, x, NULL, loglFunc_arg);
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
						loglFunc(&logll, &x[i], 1, i, x, NULL, loglFunc_arg);
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
	 */

	int i, j, ii, jj, ns;
	double xx, *f = NULL, *g = NULL, val, tmp, logll = 0.0, q_value;

	static double quadratic_term = 0.0, linear_term = 0.0; /* compute those if logdens == NULL */
#pragma omp threadprivate(quadratic_term, linear_term)

	GMRFLib_ENTER_ROUTINE;

	assert(subgraph);
	ns = subgraph->n;

	if (!logdens) {
		/*
		 * compute the quadratic and linear term 
		 */
		f = Calloc(ns, double);
		g = Calloc(ns, double);


		for (ii = 0; ii < ns; ii++) {
			i = subgraph->mothergraph_idx[ii];
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
			i = subgraph->mothergraph_idx[ii];
			linear_term -= g[ii] * x_gradient[i];
			quadratic_term += f[ii] * x_gradient[i];
		}
		if (b) {
			for (ii = 0; ii < ns; ii++) {
				i = subgraph->mothergraph_idx[ii];
				linear_term += x_gradient[i] * b[i];
			}
		}
		Free(f);
		Free(g);
	} else {
		val = -0.5 * SQR(delta) * quadratic_term + delta * linear_term;
		if (d) {
			tmp = 0.0;
			for (ii = 0; ii < ns; ii++) {
				i = subgraph->mothergraph_idx[ii];
				if (d[i]) {
					loglFunc(&logll, &x[i], 1, i, x, NULL, loglFunc_arg);
					tmp += d[i] * logll;
				}
			}
			val += tmp;
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
int GMRFLib_ai_nparam_eff(double *nparam_eff, double *nparam_eff_rel, GMRFLib_problem_tp * problem, double *c,
			  GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg)
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
			GMRFLib_2order_approx(&a[i], &b[i], &c[i], d[i], problem->mean_constr[i], i, problem->mean_constr, loglFunc,
					      loglFunc_arg, NULL, NULL);
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
					loglFunc(&loglikelihood, &xval, 1, i, problem->sample, NULL, loglFunc_arg);
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
			       GMRFLib_constr_tp * constr, GMRFLib_ai_param_tp * ai_par, GMRFLib_ai_store_tp * ai_store,
			       GMRFLib_marginal_hidden_store_tp * marginal_hidden_store)
{
	/*
	 * compute the approximation to the marginal for the hidden field at index 'idx' and return the density in *density. if
	 * (cpo_density), then compute also the density when y_idx is removed.
	 */

	char *fix = NULL, *fixx = NULL;
	int i, j, k, nd = -1, n = -1, free_ai_par = 0, n_points, ns = -1, ii, free_ai_store = 0, *i_idx, *j_idx, one = 1;
	double *x_points = NULL, x_sd, x_mean, *cond_mode = NULL, *fixed_mode = NULL, *log_density = NULL,
	    log_dens_cond, deriv_log_dens_cond = 0.0, a, *derivative = NULL, *mean_and_variance = NULL, ld0, ld1, c0, c1, deldif =
	    GMRFLib_eps(1.0 / 6.0), h2 = 0.0, inv_stdev, *cov = NULL, corr, corr_term, *covariances = NULL, alpha;

	GMRFLib_graph_tp *subgraph = NULL;
	GMRFLib_Qinv_tp *store_Qinv = NULL;
	GMRFLib_optimize_param_tp *optpar = NULL;
	GMRFLib_blockupdate_param_tp *blockpar = NULL;
	GMRFLib_problem_tp *newp = NULL;
	GMRFLib_store_tp *store = NULL;
	GMRFLib_ai_strategy_tp strategy;

#define COMPUTE_CPO_DENSITY						\
	if (cpo_density) {						\
		if (d[idx]) {						\
			double *xp = NULL, *xp_tmp = NULL,		\
				*ld = NULL, *logcor = NULL, *x_user = NULL, *work = NULL, _alpha=-1.0; \
			int itry, flag, np, np_orig = 51, _debug = 0, _one = 1, _i, npx = 8, itmp, np_new = np_orig + 2*npx; \
			double cor_eps = GMRFLib_eps(0.75), cor_max, range;	\
									\
			work = Calloc(4*np_new, double);		\
			for(itry = 0; itry < 2;	itry++)			\
			{						\
				np = np_orig;				\
				flag = 0;				\
				ld = &work[0];				\
				logcor = &work[np_new];			\
				x_user = &work[2*np_new];		\
				xp = &work[3*np_new];			\
				GMRFLib_ghq_abscissas(&xp_tmp, np);	\
				range = xp_tmp[np-1];			\
				memcpy(xp+npx, xp_tmp, np*sizeof(double)); \
				for(itmp = 0; itmp < npx; itmp++) {	\
					xp[itmp] = xp[npx] - range * (npx - itmp)/(double)npx; \
					xp[np + npx + itmp] = xp[npx + np - 1]  + range * (itmp + 1.0)/(double)npx; \
				}					\
				np = np_new;				\
				if (_debug) {				\
					if (0) \
						for(itmp = 0; itmp < np; itmp++) \
							printf("xp[%1d] = %f\n", itmp, xp[itmp]); \
					GMRFLib_density_printf(stdout, *density); \
				}					\
				GMRFLib_evaluate_nlogdensity(ld, xp, np, *density); \
				GMRFLib_density_std2user_n(x_user, xp, np, *density); \
				loglFunc(logcor, x_user, np, idx, fixed_mode, NULL, loglFunc_arg); \
				for(_i=0; _i < np; _i++) {		\
					logcor[_i] *= d[idx];		\
				}					\
				if (_debug && np) {			\
					for(_i = 0; _i < np; _i++)	\
						printf("CPO: %d BEFORE x_user %g xp %g ld %g logcor %g ld-logcor %g\n", idx,\
						       x_user[_i], xp[_i], ld[_i], logcor[_i], ld[_i]-logcor[_i]); \
				}					\
				if (itry == 1 && cor_eps > 0.0) {	\
					flag = 1;			\
					cor_max = exp(log(cor_eps) + GMRFLib_max_value(logcor, np, NULL)); \
					for(_i=0; _i < np; _i++) {	\
						ld[_i] = ld[_i] + logcor[_i] - 2.0*GMRFLib_log_apbex(cor_max, logcor[_i]); \
					}				\
				} else {				\
					daxpy_(&np, &_alpha, logcor, &_one, ld, &_one);  /* ld = ld + logcor */ \
				}					\
				GMRFLib_ai_correct_cpodens(ld, xp, &np, ai_par); \
				if (_debug && np) {			\
					for(_i = 0; _i < np; _i++)	\
						printf("CPO AFTER: %d %g %g\n", idx, xp[_i], ld[_i]);	\
				}					\
				if (np > 4) {				\
					GMRFLib_density_create(cpo_density, GMRFLib_DENSITY_TYPE_SCGAUSSIAN, np, xp, ld, \
							       (*density)->std_mean, (*density)->std_stdev, GMRFLib_FALSE); \
					if (flag && cpo_density) GMRFLib_setbit(&((*cpo_density)->flags), DENSITY_FLAGS_FAILURE); \
				} else {				\
					*cpo_density = NULL;		\
				}					\
				if (*cpo_density || itry == 1)		\
					break;				\
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
	if (ai_par->strategy == GMRFLib_AI_STRATEGY_ADAPTIVE) {
		assert(ai_par->adapt_len > 0);
		assert(ai_par->adapt_strategy != NULL);
		strategy = ai_par->adapt_strategy[idx];
	} else {
		strategy = ai_par->strategy;
	}

	/*
	 * if we use the ...CORRECTED_GAUSSIAN strategy, we have to use the linear correction. if the user have not selected an
	 * option for this, then do so depending if 'fast' option is set or not.
	 */
	if ((strategy == GMRFLib_AI_STRATEGY_MEANCORRECTED_GAUSSIAN || strategy == GMRFLib_AI_STRATEGY_MEANSKEWCORRECTED_GAUSSIAN)
	    && ai_par->linear_correction == GMRFLib_AI_LINEAR_CORRECTION_OFF) {
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
		Free(ai_store->aa);
		Free(ai_store->bb);
		Free(ai_store->cc);
		ai_store->aa = Calloc(n, double);
		ai_store->bb = Calloc(n, double);
		ai_store->cc = Calloc(n, double);

		GMRFLib_EWRAP1(GMRFLib_init_GMRF_approximation_store__intern(&(ai_store->problem),
									     (ai_store->mode ? ai_store->mode : x),
									     b, c, mean, d, loglFunc, loglFunc_arg, fixed_value,
									     graph, Qfunc, Qfunc_arg, constr, optpar, blockpar,
									     ai_store->store, ai_store->aa, ai_store->bb,
									     ai_store->cc, ai_par->gaussian_data, ai_par->cmin, 0));
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
	if (strategy == GMRFLib_AI_STRATEGY_GAUSSIAN) {
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

	if (0) {
		P(idx);
		P(x_mean);
		for (i = 0; i < n; i++)
			printf("derivative %d %.12g  fixed_mode %.12g  covariances %.12g\n", i, derivative[i], fixed_mode[i], covariances[i]);
		// exit(0);
	}


	/*
	 * if we do not use the meancorrected gaussian and the fast-option, then locate local neigb. set the derivative to zero
	 * for those sites that are not in the local neigb.
	 */
	if ((strategy == GMRFLib_AI_STRATEGY_MEANCORRECTED_GAUSSIAN || strategy == GMRFLib_AI_STRATEGY_MEANSKEWCORRECTED_GAUSSIAN)
	    && ai_par->linear_correction == GMRFLib_AI_LINEAR_CORRECTION_FAST) {
		if (fixed_value) {
			memcpy(fix, fixed_value, n * sizeof(char));
			memcpy(fixx, fixed_value, n * sizeof(char));
		}
	} else {

		if (marginal_hidden_store->subgraphs && marginal_hidden_store->subgraphs[idx]) {
			subgraph = marginal_hidden_store->subgraphs[idx];
			ns = subgraph->n;
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
			if (fixed_value) {		       /* if there are fixed values already: add these */
				for (i = 0; i < n; i++) {
					if (fixed_value[i]) {
						fix[i] = fixx[i] = 1;
						derivative[i] = 0.0;
					}
				}
			}

			/*
			 * note that idx is included in subgraph 
			 */
			GMRFLib_EWRAP1(GMRFLib_compute_subgraph(&subgraph, graph, fixx));

			/*
			 * store it for lated usage
			 */
			marginal_hidden_store->subgraphs[idx] = subgraph;
			ns = subgraph->n;
		}
	}

	fixx[idx] = 0;					       /* this is how 'fix' and 'fixx' differ */
	fix[idx] = 1;

	optpar->opt_type = GMRFLib_OPTTYPE_NR;		       /* force this option */
	store = Calloc(1, GMRFLib_store_tp);		       /* this can be used ;-) */

	if ((ai_par->linear_correction == GMRFLib_AI_LINEAR_CORRECTION_FAST) && !(ai_store->correction_term)) {
#pragma omp critical
		{
			if (!(ai_store->correction_term)) {
				double s = 1.0 / (2.0 * deldif);
				ai_store->correction_term = Calloc(n, double);	/* compute this */
				ai_store->derivative3 = Calloc(n, double);	/* and this */
				ai_store->correction_idx = Calloc(n, int);	/* and this one */

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
									      fixed_mode, loglFunc, loglFunc_arg,
									      &(ai_par->step_len), &(ai_par->stencil));
							GMRFLib_2order_approx(NULL, NULL, &c1, d[i], fixed_mode[i] + deldif, i,
									      fixed_mode, loglFunc, loglFunc_arg,
									      &(ai_par->step_len), &(ai_par->stencil));
							ai_store->derivative3[i] = -(c1 - c0) * s;	/* `-' since c is negative 2.deriv */
							ai_store->correction_term[i] = -SQR(ai_store->stdev[i]) * ai_store->derivative3[i];
						}
					}
				} else {
					// printf("RECOMPUTE derivative3 for thread %d and idx %d\n", omp_get_thread_num(), idx);
					for (ii = 0; ii < ai_store->nd; ii++) {
						i = ai_store->d_idx[ii];
						ai_store->correction_idx[ai_store->nidx++] = i;
						GMRFLib_2order_approx(NULL, NULL, &c0, d[i], fixed_mode[i] - deldif, i, fixed_mode,
								      loglFunc, loglFunc_arg, &(ai_par->step_len), &(ai_par->stencil));
						GMRFLib_2order_approx(NULL, NULL, &c1, d[i], fixed_mode[i] + deldif, i, fixed_mode,
								      loglFunc, loglFunc_arg, &(ai_par->step_len), &(ai_par->stencil));
						ai_store->derivative3[i] = -(c1 - c0) * s;	/* `-' since c is negative 2.deriv */
						ai_store->correction_term[i] = -SQR(ai_store->stdev[i]) * ai_store->derivative3[i];
					}
				}
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
			       (&newp, cond_mode, b, c, mean, d, loglFunc, loglFunc_arg, fix, graph, Qfunc, Qfunc_arg, constr,
				optpar, blockpar, store));
		if (newp) {
			memcpy(newp->sample, newp->mean_constr, n * sizeof(double));
			GMRFLib_EWRAP1(GMRFLib_evaluate(newp));
			ld0 = newp->sub_logdens;
			GMRFLib_free_problem(newp);
		} else {
			ld0 = 0.0;
		}

		for (ii = 0; ii < ns; ii++) {
			i = subgraph->mothergraph_idx[ii];
			cond_mode[i] = fixed_mode[i] + h2 * derivative[i];
		}
		GMRFLib_EWRAP1(GMRFLib_init_GMRF_approximation_store
			       (&newp, cond_mode, b, c, mean, d, loglFunc, loglFunc_arg, fix, graph, Qfunc, Qfunc_arg, constr,
				optpar, blockpar, store));
		if (newp) {
			memcpy(newp->sample, newp->mean_constr, n * sizeof(double));
			GMRFLib_EWRAP1(GMRFLib_evaluate(newp));
			ld1 = newp->sub_logdens;
			GMRFLib_free_problem(newp);
		} else {
			ld1 = 0.0;
		}
		deriv_log_dens_cond = x_sd * (ld1 - ld0) / (2.0 * h2);
		break;
	}

	if (strategy != GMRFLib_AI_STRATEGY_MEANCORRECTED_GAUSSIAN && strategy != GMRFLib_AI_STRATEGY_MEANSKEWCORRECTED_GAUSSIAN) {
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

			/*
			 * this is the fast version that take into account that x = x_mode + delta * gradient for the
			 * quadratic term
			 * 
			 * we first initialise the routine computing the linear and quadratic term, and then we can get
			 * the speedup for successive calls
			 */
			assert(subgraph);
			if (k == 0) {
				GMRFLib_ai_log_posterior_restricted(NULL,
								    fixed_mode, fixed_mode, derivative,
								    0.0, b, c, mean, d, loglFunc,
								    loglFunc_arg, fixx, graph, Qfunc, Qfunc_arg, constr, subgraph, ai_store);
			}
			GMRFLib_ai_log_posterior_restricted(&log_density[k],
							    cond_mode, fixed_mode, derivative,
							    x_points[k] * x_sd, b, c, mean, d, loglFunc,
							    loglFunc_arg, fixx, graph, Qfunc, Qfunc_arg, constr, subgraph, ai_store);
			log_density[k] -= log_dens_cond;
		}
	}

	GMRFLib_free_store(store);

	if (strategy == GMRFLib_AI_STRATEGY_MEANCORRECTED_GAUSSIAN) {
		GMRFLib_density_create_normal(density, -deriv_log_dens_cond, 1.0, x_mean, x_sd);
	} else if (strategy == GMRFLib_AI_STRATEGY_MEANSKEWCORRECTED_GAUSSIAN) {
		int np = 11, err, fail = 0;
		double *ld = NULL, *xp = NULL, xx, low, high, third_order_derivative, a_sigma, cc, sol1, aa, tmp;
		GMRFLib_sn_param_tp snp;

		int iii, jjj;
		third_order_derivative = 0.0;
		for (jjj = 0; jjj < ai_store->nidx; jjj++) {
			iii = ai_store->correction_idx[jjj];
			third_order_derivative += ai_store->derivative3[iii] * gsl_pow_3(derivative[iii]);
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
				ld[k] =
				    -0.5 * SQR(xx) - deriv_log_dens_cond * xx + (third_order_derivative / 6.0) * gsl_pow_3(xx +
															   deriv_log_dens_cond);
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
		GMRFLib_ai_param_free(ai_par);
	}
	Free(mean_and_variance);
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

		// this solves the equation for the last constraint only... 
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
			GMRFLib_EWRAP1(GMRFLib_solve_llt_sparse_matrix_special
				       (&((*problem)->qi_at_m[kk]), &((*problem)->sub_sm_fact), (*problem)->sub_graph, idx));
		} else {
			/*
			 * or solve as usual 
			 */
			GMRFLib_EWRAP1(GMRFLib_solve_llt_sparse_matrix
				       (&((*problem)->qi_at_m[kk]), 1, &((*problem)->sub_sm_fact), (*problem)->sub_graph));
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
		dgemm_("N", "N", &nc, &nc, &sub_n, &alpha, (*problem)->sub_constr->a_matrix, &nc, (*problem)->qi_at_m, &sub_n,
		       &beta, aqat_m, &nc, 1, 1);
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
#define WORK(n) &(work[work_p]); work_p += (n)

	int i, k, n, nc, ncc, one = 1, work_p, work_size, ndiv;
	double *c = NULL, *v = NULL, *w = NULL, *z = NULL, alpha = 0.0, beta = 0.0, b22 = 0.0, *constr_m_new = NULL, *t_vec =
	    NULL, *work = NULL, *tmp_m = NULL, val;

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
					dgemm_("N", "N", &nc, &nc, &n, &alpha, problem->sub_constr->a_matrix, &nc, problem->qi_at_m,
					       &n, &beta, m, &nc, 1, 1);
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

		if (0) {
			FIXME("print problem->inv_aqat_m");
			GMRFLib_matrix_fprintf(stdout, problem->inv_aqat_m, nc, nc);
			printf("c[%1d] = %g  n=%d\n", idx, c[idx], n);
			double sum = 0.0;
			for (i = 0; i < nc; i++) {
				sum += v[i] * w[i];
				printf("idx %d c %f i %d v %f w %f\n", idx, c[idx], i, v[i], w[i]);
			}
		}


		/*
		 * replacement for:: b22 = c[idx]; for (i = 0; i < nc; i++) b22 -= v[i] * w[i]; b22 = 1.0 / b22; 
		 */

		// P(ddot_(&nc, v, &one, w, &one));

		b22 = 1.0 / (c[idx] - ddot_(&nc, v, &one, w, &one));
		if (b22 <= 0.0) {
			fprintf(stderr, "\n\n*** Warning *** Numerical error gives b22 = %g <= 0.0 for idx=%1d. setting b22=%g\n\n",
				b22, idx, GSL_SQRT_DBL_MIN);
			b22 = GSL_SQRT_DBL_MIN;
		}

		for (k = 0; k < nc; k++) {
			val = b22 * w[k];
			daxpy_(&nc, &val, v, &one, &tmp_m[k * nc], &one);
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
		Free(ai_store->aa);
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
int GMRFLib_ai_z2theta(double *theta, int nhyper, double *theta_mode, double *z, gsl_vector * sqrt_eigen_values, gsl_matrix * eigen_vectors)
{
	/*
	 * compute new theta-values for given vector of z (which is N(0,I)), using the relationship
	 * 
	 * theta = theta_mode + eigen_vectors * diag(1/sqrt(eigen_values)) * z 
	 */

	size_t i, j;
	double tmp, v_ij, *u = NULL;

	u = Calloc(nhyper, double);

	for (i = 0; i < (size_t) nhyper; i++) {
		u[i] = z[i] / gsl_vector_get(sqrt_eigen_values, i);
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
int GMRFLib_ai_theta2z(double *z, int nhyper, double *theta_mode, double *theta, gsl_vector * sqrt_eigen_values, gsl_matrix * eigen_vectors)
{
	/*
	 * compute z-values for given vector of theta, using the relationship
	 * 
	 * theta = theta_mode + eigen_vectors * diag(1/sqrt_eigen_values) * z 
	 */

	size_t i, j;
	double tmp, *u = NULL;

	u = Calloc(nhyper, double);

	for (i = 0; i < (size_t) nhyper; i++) {
		u[i] = theta[i] - theta_mode[i];
	}

	for (i = 0; i < (size_t) nhyper; i++) {
		for (j = 0, tmp = 0.0; j < (size_t) nhyper; j++) {
			tmp += gsl_matrix_get(eigen_vectors, j, i) * u[j];
		}
		z[i] = tmp * gsl_vector_get(sqrt_eigen_values, i);
	}
	Free(u);

	return GMRFLib_SUCCESS;
}
int GMRFLib_init_GMRF_approximation_store__intern(GMRFLib_problem_tp ** problem, double *x, double *b, double *c, double *mean,
						  double *d, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, char *fixed_value,
						  GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg,
						  GMRFLib_constr_tp * constr, GMRFLib_optimize_param_tp * optpar,
						  GMRFLib_blockupdate_param_tp * blockupdate_par, GMRFLib_store_tp * store,
						  double *aa, double *bb, double *cc, int gaussian_data, double cmin, int nested)
{
	/*
	 * This is copy of the original routine but with optional arguments 
	 */

	int i, free_x = 0, free_b = 0, free_c = 0, free_mean = 0, free_d = 0, free_blockpar = 0, free_aa = 0, free_bb = 0, free_cc =
	    0, n, id, *idxs = NULL, nidx = 0;
	double *mode = NULL;

#define FREE_ALL if (1) { if (free_x) Free(x); if (free_b) Free(b); if (free_c) Free(c); if (free_d) Free(d); \
		if (free_mean) Free(mean); if (free_blockpar) Free(blockupdate_par); if (free_aa) Free(aa); if (free_bb) Free(bb); \
		if (free_cc) Free(cc); Free(idxs); }

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
	if (!aa) {
		free_aa = 1;
		aa = Calloc(n, double);
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

	/*
	 * these tricks are currently disabled 
	 */
	int cc_positive = 1;
	int cc_is_negative = 0;
	double cc_factor = 0.1;
	double cc_factor_mult = 1.2;
	int catch_error = 0;

	GMRFLib_problem_tp *lproblem = NULL;
	double *mode_initial = Calloc(n, double);
	double err_previous = 0;

	memcpy(mode_initial, mode, n * sizeof(double));	       /* store the starting value */

	for (iter = 0; iter < itmax; iter++) {

		memset(aa, 0, n * sizeof(double));
		memcpy(bb, b, n * sizeof(double));
		memcpy(cc, c, n * sizeof(double));

		cc_is_negative = 0;
#pragma omp parallel for private(i) schedule(static)
		for (i = 0; i < nidx; i++) {
			int idx;
			double bcoof, ccoof;

			GMRFLib_thread_id = id;
			idx = idxs[i];
			GMRFLib_2order_approx(&(aa[idx]), &bcoof, &ccoof, d[idx], mode[idx], idx, mode, loglFunc, loglFunc_arg,
					      &(optpar->step_len), &(optpar->stencil));
			cc_is_negative = (cc_is_negative || ccoof < 0.0);	/* this line IS OK! also for multithread.. */
			//if (ccoof < 0.0) printf("idx %d ccoof %.12g\n", idx, ccoof);
			if (cc_positive) {
				if (cmin == 0.0) {
					if (ccoof <= cmin) {
						FIXME1("TRY NEW STRATEGY FOR bcoof");
						// do nothing. set ccoof=0 and bcoof=0.
					} else {
						bb[idx] += bcoof;
						cc[idx] += DMAX(cmin, ccoof);
					}
				} else {
					bb[idx] += bcoof;
					cc[idx] += DMAX(cmin, ccoof);
				}
			} else {
				if (ccoof > 0.0) {
					bb[idx] += bcoof;
					cc[idx] += ccoof;
				} else {
					bb[idx] += bcoof;
					// bb[idx] += cc_factor*bcoof; /* what to use?? if any...*/
					cc[idx] += cc_factor * ccoof;
				}
			}
		}
		GMRFLib_thread_id = id;
		if (!cc_positive) {
			cc_factor = DMIN(1.0, cc_factor * cc_factor_mult);
		}

		for (i = 0; i < n; i++) {
			bb[i] += -c[i] * mean[i];
		}

		/*
		 * I thought this was quicker without store, as there is just reuse and no copy... but not.  I free lproblem below and set it to NULL, so
		 * it will always be lproblem = NULL 
		 */
		int idum;
#pragma omp parallel for private(idum) num_threads(1)
		for (idum = 0; idum < 1; idum++) {
			if (!lproblem) {
				if (GMRFLib_catch_error_for_inla) {
					int ret;
					ret = GMRFLib_init_problem_store(&lproblem, x, bb, cc, mean, graph, Qfunc, Qfunc_arg, fixed_value,
									 constr, GMRFLib_NEW_PROBLEM, store);
					if (ret != GMRFLib_SUCCESS) {
						catch_error = 1;
					}
				} else {
					GMRFLib_init_problem_store
					    (&lproblem, x, bb, cc, mean, graph, Qfunc, Qfunc_arg, fixed_value, constr, GMRFLib_NEW_PROBLEM, store);
				}
			} else {
				/*
				 * store could be NULL here I presume...? 
				 */
				GMRFLib_init_problem_store
				    (&lproblem, x, bb, cc, mean, graph, Qfunc, Qfunc_arg, fixed_value, constr, GMRFLib_KEEP_graph, store);
			}
		}

		if (catch_error) {
			lproblem = NULL;
			break;
		}

		int flag_cycle_behaviour = 0;
		double err = 0.0, f;
		if (gaussian_data) {
			f = 1.0;
		} else {
			f = DMIN(1.0, (iter + 1.0) * optpar->nr_step_factor);
		}

		// if (f != 1.0) printf("%d:%d: f = %f\n", omp_get_thread_num(), GMRFLib_thread_id, f);

		for (i = 0; i < n; i++) {
			err += SQR((lproblem)->mean_constr[i] - mode[i]);
			mode[i] += f * ((lproblem)->mean_constr[i] - mode[i]);
		}
		err = sqrt(err / n);

		if (iter == 0) {
			err_previous = err;
		} else {
			if ((float) (10.0 * err) == (float) (10.0 * err_previous)) {
				/*
				 * we're down to some rounding error and cannot get any further. this weird situation has happend. 
				 */
				flag_cycle_behaviour = 1;
			}
			if (err > 4.0 * err_previous) {
				iter += itmax;		       /* so we can restart... */
			}
			err_previous = err;
		}

		if (optpar && optpar->fp)
			fprintf(optpar->fp, "[%1d] iteration %d error %.12g\n", GMRFLib_thread_id, iter, err);

		if (gaussian_data) {
			/*
			 * I need to update 'aa' as this is not evaluated in the mode! The sum of the a's are used later
			 */
#pragma omp parallel for private(i) schedule(static)
			for (i = 0; i < nidx; i++) {
				int idx = idxs[i];
				GMRFLib_thread_id = id;
				GMRFLib_2order_approx(&(aa[idx]), NULL, NULL, d[idx], mode[idx], idx, mode, loglFunc, loglFunc_arg,
						      &(optpar->step_len), &(optpar->stencil));
			}
			GMRFLib_thread_id = id;
		}

		if (err < optpar->abserr_step || gaussian_data || flag_cycle_behaviour) {
			/*
			 * we're done!  unless we have negative elements on the diagonal...
			 */

			/*
			 * disable these cc_positive tricks now... 
			 */
			break;

			if (cc_is_negative && !cc_positive && (cc_factor < 1.0)) {
				/*
				 * do nothing 
				 */
			} else if (cc_is_negative && cc_positive) {
				FIXME("switch to cc_positive = 0");
				cc_positive = 0;
			} else {
				break;
			}
		}

		if (gsl_isnan(err))
			break;

		GMRFLib_free_problem(lproblem);
		lproblem = NULL;
	}

	if (iter < itmax) {
		*problem = lproblem;
	} else {
		*problem = NULL;
		GMRFLib_free_problem(lproblem);
	}

	if (!*problem) {
		if (nested == 1) {
			GMRFLib_ASSERT(*problem, GMRFLib_EOPTNR);
			return GMRFLib_EOPTNR;
		} else if (nested == 2) {
			return GMRFLib_EOPTNR;
		} else {
			/*
			 * fail to converge. restart with a reduced step_factor. 
			 */
			memcpy(mode, mode_initial, n * sizeof(double));	/* store the starting value */
			Free(mode_initial);
			FREE_ALL;
			GMRFLib_optimize_param_tp new_optpar;
			memcpy(&new_optpar, optpar, sizeof(GMRFLib_optimize_param_tp));
			new_optpar.nr_step_factor *= 0.5;
			new_optpar.max_iter *= 2;
			if (new_optpar.fp) {
				fprintf(new_optpar.fp,
					"\n\n%s: Optimisation fail to converge.\n\t\t\tRetry with optpar->nr_step_factor = %g and add trust-region\n",
					__GMRFLib_FuncName, new_optpar.nr_step_factor);
			}
			if (new_optpar.nr_step_factor < 1e-3) {
				return GMRFLib_EOPTNR;
			} else {
				/*
				 * add trust region; try to find the smallest 'lambda' that work fine. well, approximatly only...
				 */
				int retval, kk, ntimes = 1000, stop = 0;
				double lambda = 10000.0,       /* first value for lambda */
				    lambda_fac = 0.1,	       /* decrease it with this ammount for each iteration */
				    lambda_lim = 1e-6;	       /* value of lambda where we exit the loop */
				double *c_new = Calloc(graph->n, double);

				for (kk = 0; kk < ntimes; kk++) {

					for (i = 0; i < graph->n; i++) {
						c_new[i] = lambda * Qfunc(i, i, Qfunc_arg) + (1.0 + lambda) * c[i];
						if (ISNAN(x[i]) || ISINF(x[i])) {
							if (!mode) FIXME("MODE is NULL");
							x[i] = mode[i];
						}
					}
					retval = GMRFLib_init_GMRF_approximation_store__intern(problem, x, b, c_new, mean, d,
											       loglFunc, loglFunc_arg, fixed_value,
											       graph, Qfunc, Qfunc_arg, constr,
											       &new_optpar, blockupdate_par, store,
											       aa, bb, cc, gaussian_data, cmin, 2);
					if (stop && retval == GMRFLib_SUCCESS) {
						break;
					}
					GMRFLib_ASSERT(lambda < 1.0 / lambda_lim, GMRFLib_EOPTNR);	/* exit if lambda is to large */

					if (retval == GMRFLib_SUCCESS && lambda <= lambda_lim) {
						/*
						 * lambda is small enough; we're done
						 */
						break;
					}

					if (retval != GMRFLib_SUCCESS) {
						/*
						 * it means that previous lambda worked fine, but not this one, so lets retry the previous one and then exit.
						 */
						stop = 1;
						lambda /= lambda_fac;
					} else {
						/*
						 * we're ok, decrease lambda
						 */
						lambda *= lambda_fac;
					}

					if (retval == GMRFLib_SUCCESS) {
						/*
						 *  we're ok, restart with the obtained mode
						 */
						memcpy(x, (*problem)->mean_constr, graph->n * sizeof(double));
						GMRFLib_free_problem(*problem);
					} else {
						*problem = NULL;
						if (!stop) {
							return retval;
						}
					}
				}
				if (stop) {
					assert(stop && retval == GMRFLib_SUCCESS);
				}
				Free(c_new);
			}
		}
	}

	Free(mode_initial);
	Free(mode);					       /* not part of 'FREE_ALL' */

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
  
  \param[out] po A ppointer to a \c GMRFLib_ai_po_tp stucture. If not \c NULL, then on output it will hold a alloced \c
  GMRFLib_ai_po_tp object, containing the values of \f$\log \pi(y_i | y)\f$ for each i for which the density is computed
  and observed. Note that \c po is not created unless either \c density or \c gdensity is requested.  The created \c po
  -object, can be free'd, by using the routine \c GMRFLib_ai_po_free(). Note that computing the correct \c po require that the
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
int GMRFLib_ai_INLA(GMRFLib_density_tp *** density, GMRFLib_density_tp *** gdensity,
		    GMRFLib_density_tp *** density_transform, GMRFLib_transform_array_func_tp ** tfunc,
		    GMRFLib_density_tp *** density_hyper,
		    GMRFLib_ai_cpo_tp ** cpo, GMRFLib_ai_po_tp ** po, GMRFLib_ai_dic_tp * dic,
		    GMRFLib_ai_marginal_likelihood_tp * marginal_likelihood, GMRFLib_ai_neffp_tp * neffp,
		    char *compute, double ***hyperparam, int nhyper,
		    GMRFLib_ai_log_extra_tp * log_extra, void *log_extra_arg,
		    double *x, double *b, double *c, double *mean,
		    GMRFLib_bfunc_tp ** bfunc, double *d,
		    GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, char *fixed_value,
		    GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg,
		    GMRFLib_constr_tp * constr, GMRFLib_ai_param_tp * ai_par, GMRFLib_ai_store_tp * ai_store,
		    int nlin, GMRFLib_lc_tp ** Alin, GMRFLib_density_tp *** dlin, GMRFLib_ai_misc_output_tp * misc_output)
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
		double *_improved_mean = Calloc(graph->n, double);	\
		for(_i = 0; _i<graph->n; _i++) {			\
			if (dens[_i] && dens[_i][dens_count]){		\
				_improved_mean[_i] = dens[_i][dens_count]->user_mean; \
			} else {					\
				_improved_mean[_i] = ai_store->problem->mean_constr[_i]; \
			}						\
		}							\
		GMRFLib_ai_compute_lincomb(&(lin_dens[dens_count]), (lin_cross ? &(lin_cross[dens_count]) : NULL), nlin, Alin, _store, _improved_mean); \
		Free(_improved_mean);					\
	}

#define ADD_CONFIG(_store, _theta, _log_posterior, _log_posterior_orig)	\
	if (1) {							\
		int _i;							\
		double *_improved_mean = Calloc(graph->n, double);	\
		double *_skewness = Calloc(graph->n, double);		\
		for(_i = 0; _i<graph->n; _i++) {			\
			_skewness[_i] = NAN;				\
			if (dens[_i] && dens[_i][dens_count]){		\
				_improved_mean[_i] = dens[_i][dens_count]->user_mean; \
				_skewness[_i] = dens[_i][dens_count]->skewness;	\
			} else {					\
				_improved_mean[_i] = (_store)->problem->mean_constr[_i]; \
			}						\
		}							\
		GMRFLib_ai_store_config(misc_output, nhyper, _theta, _log_posterior, _log_posterior_orig, _improved_mean, _skewness, (_store)->problem); \
		Free(_improved_mean);					\
		Free(_skewness);					\
	}


#define CHECK_HYPER_STORAGE_FORCE(num_) CHECK_HYPER_STORAGE_INTERN(num_, 4)
#define CHECK_HYPER_STORAGE CHECK_HYPER_STORAGE_INTERN(1, 0)
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
			if (misc_output && misc_output->compute_corr_lin){ \
				lin_cross = Realloc(lin_cross, hyper_len, double *); \
			}						\
		}							\
	}

#define CHECK_DENS_STORAGE_FORCE(num_)  CHECK_DENS_STORAGE_INTERN(num_, 1)
#define CHECK_DENS_STORAGE CHECK_DENS_STORAGE_INTERN(1, 0)
#define CHECK_DENS_STORAGE_INTERN(num_, force_)				\
	if ((dens_count >= dens_max) || (force_)) {			\
		int ii_, jj_, kk_;					\
		int old_dens_max = dens_max;				\
		dens_max += num_;					\
		weights = Realloc(weights, dens_max, double);		\
		if (GMRFLib_ai_INLA_userfunc0) {			\
			userfunc_values = Realloc(userfunc_values, dens_max, double *); \
		}							\
		izs = Realloc(izs, dens_max, double *);			\
		memset(&(izs[old_dens_max]), 0, (num_) * sizeof(double *)); \
		neff = Realloc(neff, dens_max, double);			\
		for (kk_ = 0; kk_ < compute_n; kk_++) {			\
			ii_ = compute_idx[kk_];				\
			if (dens[ii_]){					\
				dens[ii_] = Realloc(dens[ii_], dens_max, GMRFLib_density_tp *); \
				for(jj_ = old_dens_max; jj_ < dens_max; jj_++) \
					dens[ii_][jj_] = NULL;		\
			}						\
			if (dens_transform[ii_]){			\
				dens_transform[ii_] = Realloc(dens_transform[ii_], dens_max, GMRFLib_density_tp *); \
				for(jj_ = old_dens_max; jj_ < dens_max; jj_++) \
					dens_transform[ii_][jj_] = NULL; \
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
		if (po) {						\
			for (ii_ = 0; ii_ < compute_n; ii_++) {		\
				jj_ = compute_idx[ii_];			\
				if (d[jj_]){	\
					po_theta[jj_] = Realloc(po_theta[jj_], dens_max, double); \
					po2_theta[jj_] = Realloc(po2_theta[jj_], dens_max, double); \
					po3_theta[jj_] = Realloc(po3_theta[jj_], dens_max, double); \
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


/* 
 * if cpo_manual, then by definition, d[ii] = 0, but the observation is still there, so we have set, temporary, d[ii] = 1.
 */
#define COMPUTE_CPO_AND_DIC						\
	if (d[ii] || ai_par->cpo_manual) {				\
		if (cpo || ai_par->cpo_manual) {			\
			failure_theta[ii][dens_count] = GMRFLib_ai_cpopit_integrate(&cpo_theta[ii][dens_count], \
										    &pit_theta[ii][dens_count], ii, cpodens, \
										    (ai_par->cpo_manual ? 1.0 : d[ii]), loglFunc, loglFunc_arg, xx_mode); \
			if (cpodens && GMRFLib_getbit(cpodens->flags, DENSITY_FLAGS_FAILURE)) { \
				failure_theta[ii][dens_count] = 1.0;	\
			}						\
		}							\
		if (dic) {						\
			deviance_theta[ii][dens_count] = GMRFLib_ai_dic_integrate(ii, dens[ii][dens_count], \
										  (ai_par->cpo_manual ? 1.0 : d[ii]), loglFunc, loglFunc_arg, xx_mode); \
}									\
	}

#define COMPUTE_PO							\
	if (d[ii]) {							\
		if (po) {						\
			GMRFLib_ai_po_integrate(&po_theta[ii][dens_count], &po2_theta[ii][dens_count], &po3_theta[ii][dens_count], \
						ii, dens[ii][dens_count], d[ii], loglFunc, loglFunc_arg, xx_mode); \
		}							\
	}

#define COMPUTE_CPO_AND_DIC_LOCAL					\
	if (d[ii] || ai_par->cpo_manual) {				\
		if (cpo || ai_par->cpo_manual) {			\
			failure_theta_local[ii] +=			\
				GMRFLib_ai_cpopit_integrate(&cpo_theta_local[ii], &pit_theta_local[ii],	\
							    ii, cpodens, (ai_par->cpo_manual ? 1.0 : d[ii]), loglFunc, loglFunc_arg, xx_mode); \
			if (cpodens && GMRFLib_getbit(cpodens->flags, DENSITY_FLAGS_FAILURE)) { \
				failure_theta_local[ii] += 1.0; \
			}						\
		}							\
		if (dic) {						\
			deviance_theta_local[ii] =			\
				GMRFLib_ai_dic_integrate(ii, dens_local[ii], (ai_par->cpo_manual ? 1.0 : d[ii]), loglFunc, loglFunc_arg, xx_mode); \
		}							\
	}

#define COMPUTE_PO_LOCAL						\
	if (d[ii]) {							\
		if (po) {						\
			GMRFLib_ai_po_integrate(&po_theta_local[ii], &po2_theta_local[ii], &po3_theta_local[ii], \
						ii, dens_local[ii], d[ii], loglFunc, loglFunc_arg, xx_mode); \
		}							\
	}


/* 
   since all threads compute the same quantity, this is it well defined
*/
#define COMPUTE_NEFF					\
	if (run_with_omp) {				\
		neff[dens_count] = ai_store_id->neff;   \
	} else {					\
		neff[dens_count] = ai_store->neff;	\
	}

#define COMPUTE_NEFF2						\
	if (run_with_omp) {					\
		neff[dens_count] = ai_store_id[id]->neff;	\
	} else {						\
		neff[dens_count] = ai_store->neff;		\
	}

#define COMPUTE_NEFF_LOCAL neff_local = ai_store_id->neff

#define COMPUTE       COMPUTE_NEFF;       COMPUTE_CPO_AND_DIC; COMPUTE_PO;
#define COMPUTE2      COMPUTE_NEFF2;      COMPUTE_CPO_AND_DIC; COMPUTE_PO;
#define COMPUTE_LOCAL COMPUTE_NEFF_LOCAL; COMPUTE_CPO_AND_DIC_LOCAL; COMPUTE_PO_LOCAL;

	int i, j, k, *k_max = NULL, *k_min = NULL, *k_maxx = NULL, *k_minn = NULL, ierr, *iz = NULL, *izz = NULL, *len =
	    NULL, *iz_axes = NULL, skip, dir, len_length, free_ai_par = 0, config_count = 0, free_compute = 0, dens_count =
	    0, dens_max, hyper_len = 0, hyper_count = 0, *compute_idx = NULL, compute_n, tmax, run_with_omp, need_Qinv = 1;

	double *hessian = NULL, *theta = NULL, *theta_mode = NULL, *x_mode = NULL, log_dens_mode, log_dens, *z = NULL, **izs =
	    NULL, *stdev_corr_pos = NULL, *stdev_corr_neg = NULL, f, w, w_origo, tref, tu, *weights = NULL, *adj_weights =
	    NULL, *hyper_z = NULL, *hyper_ldens = NULL, **userfunc_values = NULL, *inverse_hessian = NULL, *neff = NULL, *timer;
	double **cpo_theta = NULL, **po_theta = NULL, **po2_theta = NULL, **po3_theta = NULL, **pit_theta = NULL, **deviance_theta =
	    NULL, **failure_theta = NULL;
	char *tag = NULL;
	gsl_matrix *H = NULL, *eigen_vectors = NULL;
	gsl_eigen_symmv_workspace *work = NULL;
	gsl_vector *eigen_values = NULL;
	gsl_vector *sqrt_eigen_values = NULL;
	map_strd hash_table;
	GMRFLib_density_tp ***dens = NULL;
	GMRFLib_density_tp ***dens_transform = NULL;
	GMRFLib_density_tp ***lin_dens = NULL;
	GMRFLib_ai_store_tp **ais = NULL;
	double **lin_cross = NULL;
	GMRFLib_marginal_hidden_store_tp *marginal_hidden_store = NULL;

	if (fixed_value) {
		FIXME("\n\n\n\nGMRFLib_INLA() do not longer work with FIXED_VALUE; please write a wrapper.\n");
		abort();
	}

	tmax = GMRFLib_MAX_THREADS;
	run_with_omp = (IMAX(GMRFLib_openmp->max_threads_outer, GMRFLib_openmp->max_threads_inner) > 1 ? 1 : 0);

	if (!ai_par) {
		GMRFLib_default_ai_param(&ai_par);
		free_ai_par = 1;
	}
	/*
	 * otherwise, it might go very wrong below 
	 */
	GMRFLib_ASSERT(ai_par && (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_GRID ||
				  ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_AUTO ||
				  ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_EMPIRICAL_BAYES ||
				  ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_USER ||
				  ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_USER_STD ||
				  ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_USER_EXPERT ||
				  ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_CCD), GMRFLib_EPARAMETER);
	/*
	 * Simply chose int-strategy here
	 */
	if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_AUTO) {
		ai_par->int_strategy = (nhyper <= 2 ? GMRFLib_AI_INT_STRATEGY_GRID : GMRFLib_AI_INT_STRATEGY_CCD);
	}

	GMRFLib_ENTER_ROUTINE;

	if (misc_output) {
		timer = misc_output->wall_clock_time_used;
		misc_output->mode_status = 0;
	} else {
		timer = NULL;
	}

	if (timer) {
		timer[0] = GMRFLib_cpu();
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
	dens_max = 1;
	dens = Calloc(graph->n, GMRFLib_density_tp **);
	dens_transform = Calloc(graph->n, GMRFLib_density_tp **);
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
	if ((density || gdensity) && po) {
		(*po) = Calloc(1, GMRFLib_ai_po_tp);
		(*po)->n = graph->n;
		(*po)->value = Calloc(graph->n, double *);
	} else {
		po = NULL;
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
		GMRFLib_ASSERT(compute_n > 0, GMRFLib_ESNH);
		for (i = 0; i < compute_n; i++) {
			GMRFLib_ASSERT(d[compute_idx[i]] == 0.0, GMRFLib_ESNH);
		}
		if (dic) {				       /* meaningless to compute the DIC in this case */
			dic = NULL;
		}
	}

	need_Qinv = ((compute_n || ai_par->compute_nparam_eff) ? 1 : 0);

	for (i = 0; i < compute_n; i++) {
		j = compute_idx[i];
		dens[j] = Calloc(dens_max, GMRFLib_density_tp *);	/* storage for the marginals */
		if (tfunc && tfunc[j]) {
			dens_transform[j] = Calloc(dens_max, GMRFLib_density_tp *);
		}
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
	if (po) {
		po_theta = Calloc(graph->n, double *);	       /* po-value conditioned on theta */
		po2_theta = Calloc(graph->n, double *);	       /* po-value conditioned on theta */
		po3_theta = Calloc(graph->n, double *);	       /* po-value conditioned on theta */
		for (i = 0; i < compute_n; i++) {
			j = compute_idx[i];
			if (d[j]) {
				po_theta[j] = Calloc(dens_max, double);
				po2_theta[j] = Calloc(dens_max, double);
				po3_theta[j] = Calloc(dens_max, double);
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

	x_mode = Calloc(graph->n, double);

	marginal_hidden_store = Calloc(1, GMRFLib_marginal_hidden_store_tp);
	if (ai_par->strategy == GMRFLib_AI_STRATEGY_FIT_SCGAUSSIAN || ai_par->strategy == GMRFLib_AI_STRATEGY_ADAPTIVE) {
		marginal_hidden_store->n = graph->n;
		marginal_hidden_store->subgraphs = Calloc(graph->n, GMRFLib_graph_tp *);
	} else {
		marginal_hidden_store->n = 0;
		marginal_hidden_store->subgraphs = NULL;
	}

	if (timer) {
		/*
		 * preparations 
		 */
		timer[0] = GMRFLib_cpu() - timer[0];
		timer[1] = GMRFLib_cpu();
	}

	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_OPTIMIZE, (void *) &nhyper, NULL);
	if (nhyper > 0) {
		/*
		 * the first step is to locate the mode of \pi(\theta | y). here we use the domin-optimiser routine.  NOTE that this
		 * '_setup' ensure that ai_store is changed for each call to _domin_f. this is a bit dirty programming, but there is no
		 * good way to get around it for the moment.
		 */
		GMRFLib_domin_setup(hyperparam, nhyper, log_extra, log_extra_arg, compute, x, b, c, mean, bfunc, d, loglFunc,
				    loglFunc_arg, fixed_value, graph, Qfunc, Qfunc_arg, constr, ai_par, ai_store);
		/*
		 * the optimizer runs most smoothly when #threads is about nhyper+1, which is the number of `natural' threads for
		 * computing the gradient.
		 */
		theta = Calloc(nhyper, double);		       /* theta is the hyperparameters */
		theta_mode = Calloc(nhyper, double);
		z = Calloc(nhyper, double);

		/*
		 * if not set to be known, then optimise 
		 */
		if (!(ai_par->mode_known)) {

			if (ai_par->fp_log) {
				fprintf(ai_par->fp_log, "Optimise using %s\n", GMRFLib_AI_OPTIMISER_NAME(ai_par->optimiser));
			}

			switch (ai_par->optimiser) {
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
				GMRFLib_ASSERT(0 == 1, GMRFLib_EPARAMETER);
				break;
			}

			if (ai_par->fp_log) {
				fprintf(ai_par->fp_log, "Optim: Number of function evaluations = %1d\n", GMRFLib_domin_get_f_count());
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
			GMRFLib_domin_f(theta_mode, &log_dens_mode, &ierr, NULL, NULL);
			log_dens_mode *= -1.0;
			if (ai_par->fp_log) {
				fprintf(ai_par->fp_log, "Compute mode: %10.4f\n", log_dens_mode);
			}
		}

		SET_THETA_MODE;

		GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_HESSIAN, (void *) &nhyper, NULL);

		if (ai_par->fp_log) {
			fprintf(ai_par->fp_log, "Compute the Hessian using %s differences and step_size[%g]. Matrix-type [%s]\n",
				(ai_par->hessian_forward_finite_difference ? "forward" : "central"),
				ai_par->hessian_finite_difference_step_len, (ai_par->hessian_force_diagonal ? "diagonal" : "dense"));
		}

		/*
		 * The parameters for the adaptive hessian estimation is set in ai_par (hence G.ai_par in domin-interface.c).
		 */
		double log_dens_mode_save = log_dens_mode;
		int stupid_mode_iter = 0;

		hessian = Calloc(ISQR(nhyper), double);
		while (GMRFLib_domin_estimate_hessian(hessian, theta_mode, &log_dens_mode, stupid_mode_iter) != GMRFLib_SUCCESS) {
			if (!stupid_mode_iter) {
				if (ai_par->fp_log)
					fprintf(ai_par->fp_log, "Mode not sufficient accurate; switch to a stupid local search strategy.\n");
			}
			stupid_mode_iter++;

			if (log_dens_mode_save > log_dens_mode && stupid_mode_iter > ai_par->stupid_search_max_iter) {
				if (ai_par->fp_log) {
					fprintf(stderr,
						"\n\n*** Mode is not accurate yet but we have reached the rounding error level. Break.\n\n");
				}
				break;
			}
			// printf("%.12g %.12g\n", log_dens_mode_save, log_dens_mode);
			log_dens_mode_save = log_dens_mode;

			if (GMRFLib_request_optimiser_to_stop) {
				fprintf(stderr, "\n\n*** Optimiser requested to stop; stop local search..\n");
				break;
			}
			if (stupid_mode_iter >= ai_par->stupid_search_max_iter) {
				fprintf(stderr, "\n\n");
				fprintf(stderr, "***\n");
				fprintf(stderr, "*** WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING\n");
				fprintf(stderr, "***\n");
				fprintf(stderr, "*** Mode not found using the stupid local search strategy; I give up.\n");
				fprintf(stderr,
					"*** I continue with best mode found and the correspondingly Hessian-matrix (can be diagonal only).\n");
				fprintf(stderr, "*** Please rerun with possible improved initial values or do other changes!!!\n");
				fprintf(stderr, "***\n");
				fprintf(stderr, "\n\n");
				break;
				// GMRFLib_ASSERT(stupid_mode_iter < ai_par->stupid_search_max_iter, GMRFLib_EMISC);
			}
		}

		/*
		 * do this again to get the ai_store set correctly.
		 */
		SET_THETA_MODE;
		if (x_mode) {
			memcpy(x_mode, ai_store->mode, graph->n * sizeof(double));
		}

		if (stupid_mode_iter) {
			// FIXME("------------> do one function call");
			for (i = 0; i < nhyper; i++) {
				theta_mode[i] = hyperparam[i][0][0];
			}
			GMRFLib_domin_f(theta_mode, &log_dens_mode, &ierr, NULL, NULL);
			log_dens_mode *= -1.0;
			SET_THETA_MODE;
			if (x_mode) {
				memcpy(x_mode, ai_store->mode, graph->n * sizeof(double));
			}
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

		double min_pos_eigenvalue = DBL_MAX;
		for (i = 0; i < nhyper; i++) {
			double eigv = gsl_vector_get(eigen_values, (unsigned int) i);

			if (eigv > 0.0) {
				min_pos_eigenvalue = DMIN(min_pos_eigenvalue, eigv);
			}
		}
		if (min_pos_eigenvalue == DBL_MAX) {
			min_pos_eigenvalue = 1.0;	       /* if all are negative, zero included */
		}
		int a_change = 0, all_negative = 1;

		for (i = 0; i < nhyper; i++) {
			double eigv = gsl_vector_get(eigen_values, (unsigned int) i);

			all_negative = (all_negative && (eigv <= 0.0 || ISZERO(eigv)));
			if (eigv < 0.0) {
				fprintf(stderr, "\n");
				fprintf(stderr, "\t*** WARNING *** Eigenvalue %1d of the Hessian is %.6g < 0\n", i, eigv);
				fprintf(stderr, "\t*** WARNING *** Set this eigenvalue to %.6g\n", min_pos_eigenvalue);
				fprintf(stderr, "\t*** WARNING *** This have consequence for the accurancy of\n");
				fprintf(stderr, "\t*** WARNING *** the approximations; please check!!!\n");
				fprintf(stderr, "\t*** WARNING *** R-inla: Use option inla(..., control.inla = list(h = h.value), ...) \n");
				fprintf(stderr, "\t*** WARNING *** R-inla: to chose a different  `h.value'.\n");
				fprintf(stderr, "\n");

				gsl_vector_set(eigen_values, (unsigned int) i, min_pos_eigenvalue);
				a_change += 1000;
			}
		}

		if (a_change) {
			if (misc_output) {
				misc_output->mode_status += a_change;	/* not a 'good mode'... */
			}
		}

		sqrt_eigen_values = gsl_vector_alloc((unsigned int) nhyper);
		for (i = 0; i < nhyper; i++) {
			gsl_vector_set(sqrt_eigen_values, (unsigned int) i, sqrt(gsl_vector_get(eigen_values, (unsigned int) i)));
		}

		if (a_change) {
			/*
			 * rebuild the Hessian using the new eigenvalues. I should have used matrix-multiplication routines, but I had this code already from
			 * af-program.c ;-) In any case, the matrix is small...
			 */
			if (all_negative) {
				/*
				 * if all eigenvalues are negative, just set the Hessian to a diagonal matrix, and go on... 
				 */

				fprintf(stderr,
					"\n\t*** WARNING *** R-inla: All eigenvalues of the Hessian are negative. Go on with Hessian = Identity\n\n");
				memset(hessian, 0, ISQR(nhyper) * sizeof(double));
				for (i = 0; i < nhyper; i++)
					hessian[i + i * nhyper] = 1.0;
			} else {
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
		}

		/*
		 * compute the inverse hessian, for scaling purposes 
		 */
		inverse_hessian = Calloc(ISQR(nhyper), double);
		memcpy(inverse_hessian, hessian, ISQR(nhyper) * sizeof(double));
		GMRFLib_comp_posdef_inverse(inverse_hessian, nhyper);

		if (misc_output) {
			misc_output->nhyper = nhyper;
			misc_output->cov_m = Calloc(ISQR(nhyper), double);
			memcpy(misc_output->cov_m, inverse_hessian, ISQR(nhyper) * sizeof(double));
			misc_output->log_posterior_mode = log_dens_mode;

			/*
			 * I need these as well, as the correction terms needs it (and we need also the sign of the eigenvectors...). 
			 */
			misc_output->eigenvalues = Calloc(nhyper, double);
			for (i = 0; i < nhyper; i++) {
				misc_output->eigenvalues[i] = 1.0 / gsl_vector_get(eigen_values, i);	/* need the eigenvalues of the cov.mat not
													 * hessian */
			}
			misc_output->eigenvectors = Calloc(ISQR(nhyper), double);
			for (i = 0; i < nhyper; i++) {
				for (j = 0; j < nhyper; j++) {
					misc_output->eigenvectors[i + j * nhyper] = gsl_matrix_get(eigen_vectors, i, j);
				}
			}
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
			if (misc_output && misc_output->compute_corr_lin) {
				lin_cross = Calloc(hyper_len, double *);
			}
		} else {
			nlin = 0;
		}

		/*
		 * compute the corrected scalings/stdevs, if required. 
		 */
		if ((ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_CCD)
		    || (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_GRID && density_hyper &&
			(ai_par->interpolator == GMRFLib_AI_INTERPOLATOR_CCD || ai_par->interpolator == GMRFLib_AI_INTERPOLATOR_CCD_INTEGRATE))
		    // as the scalings are used for the inla.sample.hyper() function... and they do not take much time in any case
		    || 1) {
			GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_HESSIAN_SCALE, (void *) &nhyper, NULL);

			stdev_corr_pos = Calloc(nhyper, double);
			stdev_corr_neg = Calloc(nhyper, double);

			/*
			 * two versions: 1. a nhyper loop, 2. a 2*nhyper loop. 
			 */
			if (omp_get_max_threads() > nhyper) {
#pragma omp parallel for private(k) num_threads(GMRFLib_openmp->max_threads_outer)
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
							ais[GMRFLib_thread_id] =
							    GMRFLib_duplicate_ai_store(ai_store, GMRFLib_TRUE, GMRFLib_TRUE, GMRFLib_FALSE);
						}
						s = ais[GMRFLib_thread_id];
					} else {
						s = ai_store;  /* the common one */
					}

					if (opt == 0) {
						zz[kk] = 2.0;
						GMRFLib_ai_z2theta(ttheta, nhyper, theta_mode, zz, sqrt_eigen_values, eigen_vectors);
						GMRFLib_domin_f_intern(ttheta, &llog_dens, &ierr, s, NULL, NULL);
						llog_dens *= -1.0;
						f0 = log_dens_mode - llog_dens;
						stdev_corr_pos[kk] = (f0 > 0.0 ? sqrt(2.0 / f0) : 1.0);
					} else {
						zz[kk] = -2.0;
						GMRFLib_ai_z2theta(ttheta, nhyper, theta_mode, zz, sqrt_eigen_values, eigen_vectors);
						GMRFLib_domin_f_intern(ttheta, &llog_dens, &ierr, s, NULL, NULL);
						llog_dens *= -1.0;
						f0 = log_dens_mode - llog_dens;
						stdev_corr_neg[kk] = (f0 > 0.0 ? sqrt(2.0 / f0) : 1.0);
					}

					Free(zz);
					Free(ttheta);
				}
			} else {
#pragma omp parallel for private(k) num_threads(GMRFLib_openmp->max_threads_outer)
				for (k = 0; k < nhyper; k++) {
					double f0, *zz = NULL, *ttheta = NULL, llog_dens;
					GMRFLib_ai_store_tp *s = NULL;

					zz = Calloc(nhyper, double);
					ttheta = Calloc(nhyper, double);
					memset(zz, 0, nhyper * sizeof(double));
					GMRFLib_thread_id = omp_get_thread_num();

					if (omp_in_parallel()) {
						if (!ais[GMRFLib_thread_id]) {
							ais[GMRFLib_thread_id] =
							    GMRFLib_duplicate_ai_store(ai_store, GMRFLib_TRUE, GMRFLib_TRUE, GMRFLib_FALSE);
						}
						s = ais[GMRFLib_thread_id];
					} else {
						s = ai_store;  /* the common one */
					}

					zz[k] = 2.0;
					GMRFLib_ai_z2theta(ttheta, nhyper, theta_mode, zz, sqrt_eigen_values, eigen_vectors);
					GMRFLib_domin_f_intern(ttheta, &llog_dens, &ierr, s, NULL, NULL);
					llog_dens *= -1.0;
					f0 = log_dens_mode - llog_dens;
					stdev_corr_pos[k] = (f0 > 0.0 ? sqrt(2.0 / f0) : 1.0);

					zz[k] = -2.0;
					GMRFLib_ai_z2theta(ttheta, nhyper, theta_mode, zz, sqrt_eigen_values, eigen_vectors);
					GMRFLib_domin_f_intern(ttheta, &llog_dens, &ierr, s, NULL, NULL);
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
						"Compute corrected stdev for theta[%1d]: negative %f  positive %f\n", k,
						stdev_corr_neg[k], stdev_corr_pos[k]);
				}
			}

			if (misc_output) {
				misc_output->stdev_corr_pos = Calloc(nhyper, double);
				memcpy(misc_output->stdev_corr_pos, stdev_corr_pos, nhyper * sizeof(double));
				misc_output->stdev_corr_neg = Calloc(nhyper, double);
				memcpy(misc_output->stdev_corr_neg, stdev_corr_neg, nhyper * sizeof(double));
			}
		} else {
			// just fill with 1's
			if (misc_output) {
				// these are now computed, hence we use the Gaussian approximation
				misc_output->stdev_corr_pos = Calloc(nhyper, double);
				misc_output->stdev_corr_neg = Calloc(nhyper, double);
				stdev_corr_pos = Calloc(nhyper, double);
				stdev_corr_neg = Calloc(nhyper, double);
				for (k = 0; k < nhyper; k++) {
					stdev_corr_pos[k] = misc_output->stdev_corr_pos[k] = stdev_corr_neg[k] =
					    misc_output->stdev_corr_neg[k] = 1.0;
				}
			} else {
				stdev_corr_pos = Calloc(nhyper, double);
				stdev_corr_neg = Calloc(nhyper, double);
				for (k = 0; k < nhyper; k++) {
					stdev_corr_pos[k] = stdev_corr_neg[k] = 1.0;
				}
			}
		}

		if (misc_output) {
			for (k = 0; k < nhyper; k++) {
				if (ISEQUAL(misc_output->stdev_corr_pos[k], 1.0)) {
					misc_output->mode_status++;
				}
				if (ISEQUAL(misc_output->stdev_corr_neg[k], 1.0)) {
					misc_output->mode_status++;
				}
			}
		}
		// need to reset this, as ai_store is not set correctly
		if (x_mode) {
			memcpy(x_mode, ai_store->mode, graph->n * sizeof(double));
		}
		GMRFLib_domin_f(theta_mode, &log_dens_mode, &ierr, NULL, NULL);
		log_dens_mode *= -1.0;
		SET_THETA_MODE;
		if (x_mode) {
			memcpy(x_mode, ai_store->mode, graph->n * sizeof(double));
		}

		if (timer) {
			timer[1] = GMRFLib_cpu() - timer[1];
			timer[2] = GMRFLib_cpu();
		}

		if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_EMPIRICAL_BAYES) {
			if (need_Qinv) {
				GMRFLib_ai_add_Qinv_to_ai_store(ai_store);	/* add Qinv if required */
			}
			ai_store->neff = GMRFLib_AI_STORE_NEFF_NOT_COMPUTED;
			// GMRFLib_ai_store_config(misc_output, nhyper, theta_mode, 0.0, ai_store->problem);
			if (run_with_omp) {
				GMRFLib_ai_store_tp **ai_store_id = Calloc(GMRFLib_MAX_THREADS, GMRFLib_ai_store_tp *);
				GMRFLib_pardiso_thread_safe = GMRFLib_FALSE;
#pragma omp parallel for private(i) num_threads(GMRFLib_openmp->max_threads_outer)
				for (i = 0; i < compute_n; i++) {
					int ii = compute_idx[i];
					int id = omp_get_thread_num();
					GMRFLib_density_tp *cpodens = NULL;
					GMRFLib_thread_id = 0;
					if (!ai_store_id[id]) {
						ai_store_id[id] = GMRFLib_duplicate_ai_store(ai_store, GMRFLib_FALSE, GMRFLib_TRUE, GMRFLib_TRUE);
					}
					GMRFLib_ai_marginal_hidden(&dens[ii][dens_count], (cpo && (d[ii]
												   || ai_par->cpo_manual) ? &cpodens : NULL),
								   ii, x, b, c, mean, d, loglFunc, loglFunc_arg, fixed_value, graph,
								   Qfunc, Qfunc_arg, constr, ai_par, ai_store_id[id], marginal_hidden_store);
					if (tfunc && tfunc[ii]) {
						GMRFLib_transform_density(&dens_transform[ii][dens_count], dens[ii][dens_count], tfunc[ii]);
					}
					double *xx_mode = ai_store_id[id]->mode;
					COMPUTE2;
					GMRFLib_free_density(cpodens);
				}
				GMRFLib_pardiso_thread_safe = GMRFLib_TRUE;
				for (i = 0; i < GMRFLib_MAX_THREADS; i++) {
					if (!ai_store_id[i]) {
						GMRFLib_free_ai_store(ai_store_id[i]);
					}
				}
				Free(ai_store_id);

			} else {
				GMRFLib_ai_store_tp *ai_store_id = NULL;
				for (i = 0; i < compute_n; i++) {
					int ii = compute_idx[i];
					GMRFLib_density_tp *cpodens = NULL;

					GMRFLib_thread_id = 0;
					GMRFLib_ai_marginal_hidden(&dens[ii][dens_count], (cpo && (d[ii]
												   || ai_par->cpo_manual) ? &cpodens : NULL),
								   ii, x, b, c, mean, d, loglFunc, loglFunc_arg, fixed_value, graph,
								   Qfunc, Qfunc_arg, constr, ai_par, ai_store, marginal_hidden_store);
					if (tfunc && tfunc[ii]) {
						GMRFLib_transform_density(&dens_transform[ii][dens_count], dens[ii][dens_count], tfunc[ii]);
					}
					double *xx_mode = ai_store->mode;
					COMPUTE;
					GMRFLib_free_density(cpodens);
				}
			}
			if (GMRFLib_ai_INLA_userfunc0) {
				userfunc_values[dens_count] = GMRFLib_ai_INLA_userfunc0(ai_store->problem, theta, nhyper);
			}
			COMPUTE_LINDENS(ai_store);
			ADD_CONFIG(ai_store, theta_mode, 0.0, 0.0);

			izs[dens_count] = Calloc(nhyper, double);
			for (i = 0; i < nhyper; i++) {
				izs[dens_count][i] = 0;
			}
			weights[dens_count] = 0.0;
			dens_count++;

			/*
			 * END OF GMRFLib_AI_INT_STRATEGY_EMPIRICAL_BAYES 
			 */
		} else if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_CCD ||
			   ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_USER ||
			   ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_USER_STD ||
			   ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_USER_EXPERT) {
			/*
			 * use points from the ccd-design to do the integration. This includes also the deterministic
			 * integration points.
			 */
			GMRFLib_design_tp *design = NULL;

			if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_CCD) {
				GMRFLib_get_design(&design, nhyper);
			} else {
				design = ai_par->int_design;
			}

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

#pragma omp parallel for private(k, i, log_dens, dens_count, hyper_count, tref, tu, ierr) num_threads(GMRFLib_openmp->max_threads_outer)
				for (k = 0; k < design->nexperiments; k++) {

					double *z_local, *theta_local, log_dens_orig;
					GMRFLib_ai_store_tp *ai_store_id = NULL;
					GMRFLib_tabulate_Qfunc_tp *tabQfunc = NULL;
					double *bnew = NULL;

					dens_count = k;
					hyper_count = k;
					GMRFLib_thread_id = omp_get_thread_num();

					if (omp_in_parallel()) {
						if (!ais[GMRFLib_thread_id]) {
							ais[GMRFLib_thread_id] =
							    GMRFLib_duplicate_ai_store(ai_store, GMRFLib_FALSE, GMRFLib_TRUE, GMRFLib_FALSE);
						}
						ai_store_id = ais[GMRFLib_thread_id];
					} else {
						ai_store_id = ai_store;	/* the common one */
					}

					z_local = Calloc(nhyper, double);
					theta_local = Calloc(nhyper, double);

					if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_CCD) {
						for (i = 0; i < nhyper; i++) {
							z_local[i] = f * design->experiment[k][i]
							    * (design->experiment[k][i] > 0.0 ? stdev_corr_pos[i] : stdev_corr_neg[i]);
						}
					} else if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_USER_STD) {
						for (i = 0; i < nhyper; i++) {
							z_local[i] = design->experiment[k][i]
							    * (design->experiment[k][i] > 0.0 ? stdev_corr_pos[i] : stdev_corr_neg[i]);
						}
					} else if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_USER ||
						   ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_USER_EXPERT) {
						for (i = 0; i < nhyper; i++) {
							z_local[i] = design->experiment[k][i];
						}
					} else {
						assert(0 == 1);
					}

					if (design->std_scale) {
						// convert to theta_local
						GMRFLib_ai_z2theta(theta_local, nhyper, theta_mode, z_local, sqrt_eigen_values, eigen_vectors);
					} else {
						// theta_local is, by request, the same as z_local
						memcpy(theta_local, z_local, nhyper * sizeof(double));
					}
					GMRFLib_domin_f_intern(theta_local, &log_dens, &ierr, ai_store_id, &tabQfunc, &bnew);
					log_dens *= -1.0;
					log_dens_orig = log_dens;

					// make sure z_local's are aligned with theta_local's, for later usage.
					GMRFLib_ai_theta2z(z_local, nhyper, theta_mode, theta_local, sqrt_eigen_values, eigen_vectors);

					/*
					 * correct the log_dens due to the integration weights which is special for the CCD
					 * integration and the deterministic integration points
					 * 
					 */
					if (ISNAN(design->int_weight[k])) {
						// integration weights are undefined. use these for the CCD design (as it _IS_
						// the CCD design in this case
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
					} else {
						// integration weights are _given_. this is the deterministic integration points
						if (ai_par->int_strategy != GMRFLib_AI_INT_STRATEGY_USER_EXPERT) {
							log_dens += log(design->int_weight[k]);
						} else {
							// we do that later
						}
					}

					/*
					 * register the density for the marginal of the hyperparameters computations. first check space. 
					 */
					for (i = 0; i < nhyper; i++) {
						hyper_z[hyper_count * nhyper + i] = z_local[i];
					}
					hyper_ldens[hyper_count] = log_dens_orig - log_dens_mode;

					/*
					 * compute the marginals for this point. check storage 
					 */
					if (nhyper > 0) {
						if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_USER_EXPERT) {
							// In this case, the int_weights INCLUDE the log_dens
							weights[dens_count] = log(design->int_weight[k]);
						} else {
							weights[dens_count] = log_dens;
						}
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
					}
					// GMRFLib_ai_store_config(misc_output, nhyper, theta_local, log_dens,
					// ai_store_id->problem);
					ai_store_id->neff = GMRFLib_AI_STORE_NEFF_NOT_COMPUTED;
					for (i = 0; i < compute_n; i++) {
						int ii = compute_idx[i];
						GMRFLib_density_tp *cpodens = NULL;

						GMRFLib_ai_marginal_hidden(&dens[ii][dens_count],
									   (cpo && (d[ii] || ai_par->cpo_manual) ? &cpodens : NULL),
									   ii, x, bnew, c, mean, d,
									   loglFunc, loglFunc_arg, fixed_value,
									   graph, tabQfunc->Qfunc, tabQfunc->Qfunc_arg,
									   constr, ai_par, ai_store_id, marginal_hidden_store);
						if (tfunc && tfunc[ii]) {
							GMRFLib_transform_density(&dens_transform[ii][dens_count], dens[ii][dens_count], tfunc[ii]);
						}
						double *xx_mode = ai_store_id->mode;
						COMPUTE;
						GMRFLib_free_density(cpodens);
					}
					if (GMRFLib_ai_INLA_userfunc0) {
						userfunc_values[dens_count] = GMRFLib_ai_INLA_userfunc0(ai_store_id->problem, theta_local, nhyper);
					}
					COMPUTE_LINDENS(ai_store_id);
					ADD_CONFIG(ai_store_id, theta_local, log_dens, log_dens_orig);
					tu = GMRFLib_cpu() - tref;
					if (ai_par->fp_log) {
#pragma omp critical
						{
							fprintf(ai_par->fp_log, "\tconfig %2d/%1d=[", config_count++, design->nexperiments);
							for (i = 0; i < nhyper; i++) {
								fprintf(ai_par->fp_log, " %5.3f", z_local[i]);
							}
							/*
							 * we need to use the log_dens_orig as the other one is also included the integration weights. 
							 */
							fprintf(ai_par->fp_log, "] log(rel.dens)=%5.3f, [%1d] accept, compute,",
								log_dens_orig - log_dens_mode, omp_get_thread_num());
							fprintf(ai_par->fp_log, " %.2fs\n", tu);
						}
					}

					GMRFLib_free_tabulate_Qfunc(tabQfunc);
					Free(bnew);
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

					GMRFLib_tabulate_Qfunc_tp *tabQfunc = NULL;
					double *bnew = NULL;
					double log_dens_orig;

					if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_CCD) {
						for (i = 0; i < nhyper; i++) {
							z[i] = f * design->experiment[k][i]
							    * (design->experiment[k][i] > 0.0 ? stdev_corr_pos[i] : stdev_corr_neg[i]);
						}
					} else if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_USER_STD) {
						for (i = 0; i < nhyper; i++) {
							z[i] = design->experiment[k][i]
							    * (design->experiment[k][i] > 0.0 ? stdev_corr_pos[i] : stdev_corr_neg[i]);
						}
					} else if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_USER ||
						   ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_USER_EXPERT) {
						for (i = 0; i < nhyper; i++) {
							z[i] = design->experiment[k][i];
						}
					} else {
						assert(0 == 1);
					}

					if (design->std_scale) {
						// convert to theta
						GMRFLib_ai_z2theta(theta, nhyper, theta_mode, z, sqrt_eigen_values, eigen_vectors);
					} else {
						// theta is, by request, the same as z
						memcpy(theta, z, nhyper * sizeof(double));
					}
					// make sure z's are aligned with theta's, for later usage.
					GMRFLib_ai_theta2z(z, nhyper, theta_mode, theta, sqrt_eigen_values, eigen_vectors);

					GMRFLib_domin_f(theta, &log_dens, &ierr, &tabQfunc, &bnew);
					log_dens *= -1.0;
					log_dens_orig = log_dens;

					if (ai_par->fp_log) {
						fprintf(ai_par->fp_log, "\tconfig %2d=[", config_count++);
						for (i = 0; i < nhyper; i++) {
							fprintf(ai_par->fp_log, " %5.3f", z[i]);
						}
						fprintf(ai_par->fp_log, "] log(rel.dens)=%5.3f, [%1d] accept, compute,",
							log_dens - log_dens_mode, omp_get_thread_num());
					}

					/*
					 * correct the log_dens due to the integration weights which is special for the CCD integration:
					 * double the weights for the points not in the center
					 */
					if (ISNAN(design->int_weight[k])) {
						// integration weights are undefined. use these for the CCD design (as it _IS_
						// the CCD design in this case
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
					} else {
						// integration weights are _given_. this is the deterministic integration points
						log_dens += log(design->int_weight[k]);
					}

					/*
					 * register the density for the marginal of the hyperparameters computations. first check space. 
					 */
					CHECK_HYPER_STORAGE;
					for (i = 0; i < nhyper; i++) {
						hyper_z[hyper_count * nhyper + i] = z[i];
					}
					hyper_ldens[hyper_count] = log_dens_orig - log_dens_mode;
					hyper_count++;

					/*
					 * compute the marginals for this point. check storage 
					 */
					CHECK_DENS_STORAGE;
					if (nhyper > 0) {
						if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_USER_EXPERT) {
							// In this case, the int_weights INCLUDE the log_dens
							weights[dens_count] = log(design->int_weight[k]);
						} else {
							weights[dens_count] = log_dens;
						}
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
					}
					ai_store->neff = GMRFLib_AI_STORE_NEFF_NOT_COMPUTED;

					// GMRFLib_ai_store_config(misc_output, nhyper, theta, log_dens, ai_store->problem);
					if (run_with_omp) {
						GMRFLib_ai_store_tp **ai_store_id = Calloc(GMRFLib_MAX_THREADS, GMRFLib_ai_store_tp *);
						GMRFLib_pardiso_thread_safe = GMRFLib_FALSE;
#pragma omp parallel for private(i) num_threads(GMRFLib_openmp->max_threads_outer)
						for (i = 0; i < compute_n; i++) {
							int ii = compute_idx[i];
							GMRFLib_density_tp *cpodens = NULL;
							int id = omp_get_thread_num();
							if (!ai_store_id[id]) {
								ai_store_id[id] =
								    GMRFLib_duplicate_ai_store(ai_store, GMRFLib_FALSE, GMRFLib_TRUE, GMRFLib_TRUE);
							}
							GMRFLib_thread_id = 0;
							GMRFLib_ai_marginal_hidden(&dens[ii][dens_count], (cpo && (d[ii]
														   || ai_par->cpo_manual) ?
													   &cpodens : NULL), ii, x,
										   bnew, c, mean, d, loglFunc, loglFunc_arg,
										   fixed_value, graph, tabQfunc->Qfunc,
										   tabQfunc->Qfunc_arg, constr, ai_par, ai_store_id[id],
										   marginal_hidden_store);
							if (tfunc && tfunc[ii]) {
								GMRFLib_transform_density(&dens_transform[ii][dens_count],
											  dens[ii][dens_count], tfunc[ii]);
							}
							double *xx_mode = ai_store_id[id]->mode;
							COMPUTE2;
							GMRFLib_free_density(cpodens);
						}
						GMRFLib_pardiso_thread_safe = GMRFLib_TRUE;
						for (i = 0; i < GMRFLib_MAX_THREADS; i++) {
							if (!ai_store_id[i]) {
								GMRFLib_free_ai_store(ai_store_id[i]);
							}
						}
						Free(ai_store_id);

					} else {

						for (i = 0; i < compute_n; i++) {
							int ii = compute_idx[i];
							GMRFLib_density_tp *cpodens = NULL;

							GMRFLib_thread_id = 0;
							GMRFLib_ai_marginal_hidden(&dens[ii][dens_count], (cpo && (d[ii]
														   || ai_par->cpo_manual) ?
													   &cpodens : NULL), ii, x,
										   bnew, c, mean, d, loglFunc, loglFunc_arg,
										   fixed_value, graph, tabQfunc->Qfunc,
										   tabQfunc->Qfunc_arg, constr, ai_par, ai_store,
										   marginal_hidden_store);
							if (tfunc && tfunc[ii]) {
								GMRFLib_transform_density(&dens_transform[ii][dens_count],
											  dens[ii][dens_count], tfunc[ii]);
							}
							double *xx_mode = ai_store->mode;
							GMRFLib_ai_store_tp *ai_store_id = ai_store;
							COMPUTE;
							GMRFLib_free_density(cpodens);
						}
					}

					if (GMRFLib_ai_INLA_userfunc0) {
						userfunc_values[dens_count] = GMRFLib_ai_INLA_userfunc0(ai_store->problem, theta, nhyper);
					}
					COMPUTE_LINDENS(ai_store);
					ADD_CONFIG(ai_store, theta, log_dens, log_dens_orig);
					tu = GMRFLib_cpu() - tref;
					if (ai_par->fp_log) {
						fprintf(ai_par->fp_log, " %.2fs\n", tu);
					}
					GMRFLib_free_tabulate_Qfunc(tabQfunc);
					Free(bnew);
					dens_count++;
				}
			}

			/*
			 * END OF GMRFLib_AI_INT_STRATEGY_CCD / USER / USER_STD
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

#pragma omp parallel for private(i, log_dens, tref, tu, ierr, kk) num_threads(GMRFLib_openmp->max_threads_outer)
				for (kk = 0; kk < pool->nconfig; kk++) {
					GMRFLib_ai_store_tp *ai_store_id = NULL;
					GMRFLib_density_tp **dens_local = NULL;
					GMRFLib_density_tp **dens_local_transform = NULL;
					double *z_local = NULL, *theta_local = NULL, *userfunc_values_local =
					    NULL, weights_local, val, neff_local = 0.0;
					double *cpo_theta_local = NULL, *po_theta_local = NULL, *po2_theta_local =
					    NULL, *po3_theta_local = NULL, *pit_theta_local = NULL, *failure_theta_local =
					    NULL, *deviance_theta_local = NULL;
					int err, *iz_local = NULL;
					size_t idx;
					GMRFLib_tabulate_Qfunc_tp *tabQfunc = NULL;
					double *bnew = NULL;

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
								ais[GMRFLib_thread_id] =
								    GMRFLib_duplicate_ai_store(ai_store, GMRFLib_FALSE, GMRFLib_TRUE,
											       GMRFLib_FALSE);
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
						GMRFLib_ai_z2theta(theta_local, nhyper, theta_mode, z_local, sqrt_eigen_values, eigen_vectors);
						GMRFLib_domin_f_intern(theta_local, &log_dens, &ierr, ai_store_id, &tabQfunc, &bnew);
						log_dens *= -1.0;

						val = log_dens - log_dens_mode;
						if ((ISINF(val) || ISNAN(val)) || -val > ai_par->diff_log_dens) {
							GMRFLib_ai_pool_set(pool, idx, val);
							if (ai_par->fp_log) {
#pragma omp critical
								{
									fprintf(ai_par->fp_log, "\tconfig %2d=[", config_count++);
									for (i = 0; i < nhyper; i++) {
										fprintf(ai_par->fp_log, " %5.3f", z_local[i]);
									}
									fprintf(ai_par->fp_log,
										"] log(rel.dens)=%5.3f, reject, %.2fs\n", val,
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
							}
							ai_store_id->neff = GMRFLib_AI_STORE_NEFF_NOT_COMPUTED;
							dens_local = Calloc(graph->n, GMRFLib_density_tp *);
							dens_local_transform = Calloc(graph->n, GMRFLib_density_tp *);
							if (cpo) {
								cpo_theta_local = Calloc(graph->n, double);
								pit_theta_local = Calloc(graph->n, double);
								failure_theta_local = Calloc(graph->n, double);
							}
							if (po) {
								po_theta_local = Calloc(graph->n, double);
								po2_theta_local = Calloc(graph->n, double);
								po3_theta_local = Calloc(graph->n, double);
							}
							if (dic) {
								deviance_theta_local = Calloc(graph->n, double);
							}
							for (i = 0; i < compute_n; i++) {
								GMRFLib_density_tp *cpodens = NULL;
								int ii;
								double *xx_mode = NULL;

								ii = compute_idx[i];
								GMRFLib_ai_marginal_hidden(&dens_local[ii], (cpo && (d[ii]
														     || ai_par->cpo_manual)
													     ? &cpodens : NULL), ii,
											   x, bnew, c, mean, d, loglFunc,
											   loglFunc_arg, fixed_value, graph,
											   tabQfunc->Qfunc, tabQfunc->Qfunc_arg,
											   constr, ai_par, ai_store_id, marginal_hidden_store);
								if (tfunc && tfunc[ii]) {
									GMRFLib_transform_density(&dens_local_transform[ii],
												  dens_local[ii], tfunc[ii]);
								}
								xx_mode = ai_store_id->mode;

								COMPUTE_LOCAL;
								GMRFLib_free_density(cpodens);
							}
							if (GMRFLib_ai_INLA_userfunc0) {
								userfunc_values_local =
								    GMRFLib_ai_INLA_userfunc0(ai_store_id->problem, theta_local, nhyper);
							}
							tu = GMRFLib_cpu() - tref;

#pragma omp critical
							{
								int ii;
								if (ai_par->fp_log) {
									{
										fprintf(ai_par->fp_log, "\tconfig %2d=[", config_count++);
										for (i = 0; i < nhyper; i++) {
											fprintf(ai_par->fp_log, " %5.3f", z_local[i]);
										}
										fprintf(ai_par->fp_log,
											"] log(rel.dens)=%5.3f, [%1d] accept, compute,",
											val, omp_get_thread_num());
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
									if (tfunc && tfunc[ii]) {
										dens_transform[ii][dens_count] = dens_local_transform[ii];
									}
								}
								COMPUTE_LINDENS(ai_store_id);
								ADD_CONFIG(ai_store_id, theta_local, log_dens, log_dens);
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
								if (po) {
									for (i = 0; i < compute_n; i++) {
										ii = compute_idx[i];
										if (d[ii]) {
											po_theta[ii][dens_count] = po_theta_local[ii];
											po2_theta[ii][dens_count] = po2_theta_local[ii];
											po3_theta[ii][dens_count] = po3_theta_local[ii];
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
					Free(po_theta_local);
					Free(po2_theta_local);
					Free(po3_theta_local);
					Free(dens_local);
					Free(dens_local_transform);
					Free(deviance_theta_local);
					Free(iz_local);
					Free(pit_theta_local);
					Free(failure_theta_local);
					Free(theta_local);
					Free(z_local);
					GMRFLib_free_tabulate_Qfunc(tabQfunc);
					Free(bnew);
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
							GMRFLib_tabulate_Qfunc_tp *tabQfunc = NULL;
							double *bnew = NULL;

							z[k] += dir * ai_par->dz;
							iz[k] += dir;

							/*
							 * note that _domin_f stores calculations in the ai_store even though its not in the argument
							 * list 
							 */
							GMRFLib_ai_z2theta(theta, nhyper, theta_mode, z, sqrt_eigen_values, eigen_vectors);
							GMRFLib_domin_f(theta, &log_dens, &ierr, &tabQfunc, &bnew);
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
							}
							ai_store->neff = GMRFLib_AI_STORE_NEFF_NOT_COMPUTED;

							// GMRFLib_ai_store_config(misc_output, nhyper, theta, log_dens,
							// ai_store->problem);
							if (run_with_omp) {
								GMRFLib_ai_store_tp **ai_store_id =
								    Calloc(GMRFLib_MAX_THREADS, GMRFLib_ai_store_tp *);
								GMRFLib_pardiso_thread_safe = GMRFLib_FALSE;
#pragma omp parallel for private(i) num_threads(GMRFLib_openmp->max_threads_outer)
								for (i = 0; i < compute_n; i++) {
									int ii = compute_idx[i];
									GMRFLib_density_tp *cpodens = NULL;
									int id = omp_get_thread_num();
									if (!ai_store_id[id]) {
										ai_store_id[id] =
										    GMRFLib_duplicate_ai_store(ai_store, GMRFLib_FALSE,
													       GMRFLib_TRUE, GMRFLib_TRUE);
									}

									GMRFLib_thread_id = 0;
									GMRFLib_ai_marginal_hidden(&dens[ii][dens_count], (cpo && (d[ii]
																   ||
																   ai_par->cpo_manual)
															   ? &cpodens : NULL), ii,
												   x, bnew, c, mean, d, loglFunc, loglFunc_arg,
												   fixed_value, graph, tabQfunc->Qfunc,
												   tabQfunc->Qfunc_arg, constr, ai_par,
												   ai_store_id[id], marginal_hidden_store);
									if (tfunc && tfunc[ii]) {
										GMRFLib_transform_density(&dens_transform[ii]
													  [dens_count],
													  dens[ii][dens_count], tfunc[ii]);
									}
									double *xx_mode = ai_store_id[id]->mode;
									COMPUTE2;
									GMRFLib_free_density(cpodens);
								}
								GMRFLib_pardiso_thread_safe = GMRFLib_TRUE;
								for (i = 0; i < GMRFLib_MAX_THREADS; i++) {
									if (!ai_store_id[i]) {
										GMRFLib_free_ai_store(ai_store_id[i]);
									}
								}
								Free(ai_store_id);

							} else {

								for (i = 0; i < compute_n; i++) {
									int ii = compute_idx[i];
									GMRFLib_density_tp *cpodens = NULL;

									GMRFLib_thread_id = 0;
									GMRFLib_ai_marginal_hidden(&dens[ii][dens_count], (cpo && (d[ii]
																   ||
																   ai_par->cpo_manual)
															   ? &cpodens : NULL), ii,
												   x, bnew, c, mean, d, loglFunc, loglFunc_arg,
												   fixed_value, graph, tabQfunc->Qfunc,
												   tabQfunc->Qfunc_arg, constr, ai_par, ai_store,
												   marginal_hidden_store);
									if (tfunc && tfunc[ii]) {
										GMRFLib_transform_density(&dens_transform[ii]
													  [dens_count],
													  dens[ii][dens_count], tfunc[ii]);
									}
									double *xx_mode = ai_store->mode;
									GMRFLib_ai_store_tp *ai_store_id = ai_store;
									COMPUTE;
									GMRFLib_free_density(cpodens);
								}
							}
							if (GMRFLib_ai_INLA_userfunc0) {
								userfunc_values[dens_count] =
								    GMRFLib_ai_INLA_userfunc0(ai_store->problem, theta, nhyper);
							}
							COMPUTE_LINDENS(ai_store);
							ADD_CONFIG(ai_store, theta, log_dens, log_dens);
							tu = GMRFLib_cpu() - tref;
							if (ai_par->fp_log) {
								fprintf(ai_par->fp_log, " %.2fs\n", tu);
							}

							GMRFLib_free_tabulate_Qfunc(tabQfunc);
							Free(bnew);
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
						GMRFLib_tabulate_Qfunc_tp *tabQfunc = NULL;
						double *bnew = NULL;

						GMRFLib_ai_z2theta(theta, nhyper, theta_mode, z, sqrt_eigen_values, eigen_vectors);
						GMRFLib_domin_f(theta, &log_dens, &ierr, &tabQfunc, &bnew);
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
									fprintf(ai_par->fp_log,
										" diff to large, stop and skip skip more distant configs\n");
								} else {
									fprintf(ai_par->fp_log, " diff to large, stop\n");
								}
							}
							skip = 1;

							if (ai_par->skip_configurations) {
								/*
								 * mark all those configurations >= than this as to be skipped as well 
								 */
								GMRFLib_ai_skip_configurations(&hash_table, k, iz, izz, len, k_max,
											       len_length, nhyper);
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
							}
							ai_store->neff = GMRFLib_AI_STORE_NEFF_NOT_COMPUTED;

							// GMRFLib_ai_store_config(misc_output, nhyper, theta, log_dens,
							// ai_store->problem);
							if (run_with_omp) {
								GMRFLib_ai_store_tp **ai_store_id =
								    Calloc(GMRFLib_MAX_THREADS, GMRFLib_ai_store_tp *);
								GMRFLib_pardiso_thread_safe = GMRFLib_FALSE;
#pragma omp parallel for private(i) num_threads(GMRFLib_openmp->max_threads_outer)
								for (i = 0; i < compute_n; i++) {
									int ii = compute_idx[i];
									GMRFLib_density_tp *cpodens = NULL;
									int id = omp_get_thread_num();
									if (!ai_store_id[id]) {
										ai_store_id[id] =
										    GMRFLib_duplicate_ai_store(ai_store, GMRFLib_FALSE,
													       GMRFLib_TRUE, GMRFLib_TRUE);
									}

									GMRFLib_thread_id = 0;
									GMRFLib_ai_marginal_hidden(&dens[ii][dens_count], (cpo && (d[ii]
																   ||
																   ai_par->cpo_manual)
															   ? &cpodens : NULL), ii,
												   x, bnew, c, mean, d, loglFunc, loglFunc_arg,
												   fixed_value, graph, tabQfunc->Qfunc,
												   tabQfunc->Qfunc_arg, constr, ai_par,
												   ai_store_id[id], marginal_hidden_store);
									if (tfunc && tfunc[ii]) {
										GMRFLib_transform_density(&dens_transform[ii]
													  [dens_count],
													  dens[ii][dens_count], tfunc[ii]);
									}
									double *xx_mode = ai_store_id[id]->mode;
									COMPUTE2;
									GMRFLib_free_density(cpodens);
								}
								GMRFLib_pardiso_thread_safe = GMRFLib_TRUE;
								for (i = 0; i < GMRFLib_MAX_THREADS; i++) {
									if (!ai_store_id[i]) {
										GMRFLib_free_ai_store(ai_store_id[i]);
									}
								}
								Free(ai_store_id);

							} else {

								GMRFLib_ai_store_tp *ai_store_id = NULL;
								for (i = 0; i < compute_n; i++) {
									int ii = compute_idx[i];
									GMRFLib_density_tp *cpodens = NULL;

									GMRFLib_thread_id = 0;
									GMRFLib_ai_marginal_hidden(&dens[ii][dens_count], (cpo && (d[ii]
																   ||
																   ai_par->cpo_manual)
															   ? &cpodens : NULL), ii,
												   x, bnew, c, mean, d, loglFunc, loglFunc_arg,
												   fixed_value, graph, tabQfunc->Qfunc,
												   tabQfunc->Qfunc_arg, constr, ai_par, ai_store,
												   marginal_hidden_store);
									if (tfunc && tfunc[ii]) {
										GMRFLib_transform_density(&dens_transform[ii]
													  [dens_count],
													  dens[ii][dens_count], tfunc[ii]);
									}
									double *xx_mode = ai_store->mode;
									COMPUTE;
									GMRFLib_free_density(cpodens);
								}
							}
							if (GMRFLib_ai_INLA_userfunc0) {
								userfunc_values[dens_count] =
								    GMRFLib_ai_INLA_userfunc0(ai_store->problem, theta, nhyper);
							}
							COMPUTE_LINDENS(ai_store);
							ADD_CONFIG(ai_store, theta, log_dens, log_dens);
							tu = GMRFLib_cpu() - tref;
							if (ai_par->fp_log) {
								fprintf(ai_par->fp_log, " %.2fs\n", tu);
							}
							dens_count++;
						}
						GMRFLib_free_tabulate_Qfunc(tabQfunc);
						Free(bnew);
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
		if (timer) {
			timer[1] = 0.0;
			timer[2] = GMRFLib_cpu();
		}

		if (nlin) {
			lin_dens = Calloc(1, GMRFLib_density_tp **);
			if (misc_output && misc_output->compute_corr_lin) {
				lin_cross = Calloc(1, double *);
			}
		}
		GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_OPTIMIZE, NULL, NULL);

		int idum;
#pragma omp parallel for private(idum) num_threads(GMRFLib_openmp->max_threads_outer)
		for (idum = 0; idum < 1; idum++) {	       /* YES YES YES, otherwise PARDISO go nuts! */
			/*
			 * In this case the contents of ai_store is NULL, we need to recompute the Gaussian approximation since
			 * the contents of ai_store is NULL in this case.
			 */
			double tmp_logdens;
			double *bnew = NULL, con = 0.0;
			GMRFLib_bnew(&bnew, &con, graph->n, b, bfunc);
			GMRFLib_ai_marginal_hyperparam(&tmp_logdens, x, bnew, c, mean, d,
						       loglFunc, loglFunc_arg, fixed_value, graph, Qfunc, Qfunc_arg, constr, ai_par, ai_store);
			log_dens_mode = tmp_logdens + con + log_extra(NULL, nhyper, log_extra_arg);	/* nhyper=0, so theta=NULL is ok */
			GMRFLib_ai_add_Qinv_to_ai_store(ai_store);
			Free(bnew);
		}
		ai_store->neff = GMRFLib_AI_STORE_NEFF_NOT_COMPUTED;

		if (run_with_omp) {
			GMRFLib_ai_store_tp **ai_store_id = Calloc(GMRFLib_MAX_THREADS, GMRFLib_ai_store_tp *);
			double *bnew = NULL, con = 0.0;
			GMRFLib_bnew(&bnew, &con, graph->n, b, bfunc);
			GMRFLib_pardiso_thread_safe = GMRFLib_FALSE;
#pragma omp parallel for private(i) num_threads(GMRFLib_openmp->max_threads_outer)
			for (i = 0; i < compute_n; i++) {
				int id = omp_get_thread_num();
				int ii = compute_idx[i];
				GMRFLib_density_tp *cpodens = NULL;

				if (!ai_store_id[id]) {
					ai_store_id[id] = GMRFLib_duplicate_ai_store(ai_store, GMRFLib_FALSE, GMRFLib_TRUE, GMRFLib_TRUE);
				}
				GMRFLib_thread_id = 0;
				GMRFLib_ai_marginal_hidden(&dens[ii][dens_count], (cpo && (d[ii]
											   || ai_par->cpo_manual) ? &cpodens : NULL), ii, x,
							   bnew, c, mean, d, loglFunc, loglFunc_arg, fixed_value, graph, Qfunc,
							   Qfunc_arg, constr, ai_par, ai_store_id[id], marginal_hidden_store);
				if (tfunc && tfunc[ii]) {
					GMRFLib_transform_density(&dens_transform[ii][dens_count], dens[ii][dens_count], tfunc[ii]);
				}
				double *xx_mode = ai_store_id[id]->mode;
				COMPUTE2;
				GMRFLib_free_density(cpodens);
			}
			GMRFLib_pardiso_thread_safe = GMRFLib_TRUE;
			for (i = 0; i < GMRFLib_MAX_THREADS; i++) {
				if (!ai_store_id[i]) {
					GMRFLib_free_ai_store(ai_store_id[i]);
				}
			}
			Free(ai_store_id);
			Free(bnew)
		} else {
			GMRFLib_ai_store_tp *ai_store_id = NULL;
			double *bnew = NULL, con = 0.0;
			GMRFLib_bnew(&bnew, &con, graph->n, b, bfunc);

			for (i = 0; i < compute_n; i++) {
				int ii = compute_idx[i];
				GMRFLib_density_tp *cpodens = NULL;

				GMRFLib_thread_id = 0;

				GMRFLib_ai_marginal_hidden(&dens[ii][dens_count], (cpo && (d[ii]
											   || ai_par->cpo_manual) ? &cpodens : NULL), ii, x,
							   bnew, c, mean, d, loglFunc, loglFunc_arg, fixed_value, graph, Qfunc,
							   Qfunc_arg, constr, ai_par, ai_store, marginal_hidden_store);
				if (tfunc && tfunc[ii]) {
					GMRFLib_transform_density(&dens_transform[ii][dens_count], dens[ii][dens_count], tfunc[ii]);
				}
				double *xx_mode = ai_store->mode;
				COMPUTE;
				GMRFLib_free_density(cpodens);
			}
			memcpy(x_mode, ai_store->mode, graph->n * sizeof(double));
			Free(bnew);
		}
		if (GMRFLib_ai_INLA_userfunc0) {
			userfunc_values[dens_count] = GMRFLib_ai_INLA_userfunc0(ai_store->problem, theta, nhyper);
		}
		COMPUTE_LINDENS(ai_store);
		ADD_CONFIG(ai_store, NULL, log_dens_mode, log_dens_mode);
		weights[dens_count] = 0.0;
		dens_count++;

		/*
		 * END OF nhyper == 0 
		 */
	}

	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_DEFAULT, NULL, NULL);
	if (timer) {
		timer[2] = GMRFLib_cpu() - timer[2];
		timer[3] = GMRFLib_cpu();
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
	if (density_transform) {
		*density_transform = Calloc(graph->n, GMRFLib_density_tp *);
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
				fprintf(ai_par->fp_log, " %5.3f", izs[j][k]);
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
		/*
		 * need a separate strategy here, as this might take time if the number of points are large, and we essentially want to do this loop with max
		 * num_threads.
		 */
		GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_COMBINE, NULL, NULL);
#pragma omp parallel for private(j) num_threads(GMRFLib_openmp->max_threads_outer)
		for (j = 0; j < compute_n; j++) {
			int ii = compute_idx[j];
			// fprintf(stderr, "thead %d ii %d\n", omp_get_thread_num(), ii);
			GMRFLib_density_tp *dens_combine, *gdens_combine;
			GMRFLib_density_combine((density ? &dens_combine : NULL), (gdensity ? &gdens_combine : NULL), dens_count,
						dens[ii], adj_weights);
			if (density) {
				(*density)[ii] = dens_combine;
				// printf("ii %d mean before %g mean after %g\n", ii, dens[ii][0]->user_mean,
				// dens_combine->user_mean);
			}
			if (gdensity) {
				(*gdensity)[ii] = gdens_combine;
			}

			if (tfunc && tfunc[ii]) {
				GMRFLib_density_tp *dens_c;
				GMRFLib_density_combine((density_transform ? &dens_c : NULL), NULL, dens_count, dens_transform[ii], adj_weights);
				(*density_transform)[ii] = dens_c;
			}
		}
		GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_DEFAULT, NULL, NULL);
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

		if (misc_output && misc_output->compute_corr_lin) {
			double *ptmp;
			misc_output->corr_lin = ptmp = Calloc(ISQR(nlin), double);
			misc_output->cov_lin = Calloc(ISQR(nlin), double);

			for (i = 0; i < nlin; i++) {
				for (j = i; j < nlin; j++) {
					for (k = 0; k < dens_count; k++) {
						ptmp[i + j * nlin] += adj_weights[k] * lin_cross[k][i + j * nlin];
					}
					ptmp[j + i * nlin] = ptmp[i + j * nlin];
				}
			}
			memcpy(misc_output->cov_lin, ptmp, ISQR(nlin) * sizeof(double));

			double *ptmp_scale = Calloc(ISQR(nlin), double);
			for (i = 0; i < nlin; i++) {
				ptmp_scale[i + i * nlin] = 1.0 / sqrt(ptmp[i + i * nlin]);
			}

			for (i = 0; i < nlin; i++) {
				for (j = i + 1; j < nlin; j++) {
					ptmp[i + j * nlin] = ptmp[i + j * nlin] * ptmp_scale[i + i * nlin] * ptmp_scale[j + j * nlin];
					ptmp[j + i * nlin] = ptmp[i + j * nlin];
				}
			}

			for (i = 0; i < nlin; i++) {
				ptmp[i + i * nlin] = 1.0;
			}

			Free(ptmp_scale);
		}

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
				printf("Expected effective number of parameters: %.3f(%.3f),  eqv.#replicates: %.3f\n", neff_eval,
				       neff_sd, ndata / neff_eval);
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
		if (!(ai_par->cpo_manual)) {
			/*
			 * In this case, the \pi(theta|y) is computed with all the data. then we need to correct
			 *
			 * first we need to compute the normalising constants for \pi(theta|y_i) for each i. This we call Z[i].
			 * 
			 * Note that \pi(theta_j | y_{-i}) = adj_weights[j] / cpo_theta[i][j] / Z[i];
			 *
			 * We do not correct for int.strategy = user.expert, for which the weights are given.
			 */
			double *Z = Calloc(graph->n, double);

			for (j = 0; j < compute_n; j++) {
				int jj, ii;

				ii = compute_idx[j];
				if (cpo_theta[ii]) {
					for (jj = 0; jj < dens_count; jj++) {
						if (!ISNAN(cpo_theta[ii][jj]))	/* we ignore those that have failed */
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
					
					if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_USER_EXPERT) {
						for (jj = 0, evalue = evalue_one = 0.0; jj < dens_count; jj++) {
							if (!ISNAN(cpo_theta[ii][jj])) {
								evalue += cpo_theta[ii][jj] * adj_weights[jj];
								evalue_one += adj_weights[jj]; 
							}
						}
					} else {
						// here, we correct for adjusting pi(theta_j|y_{-i})
						for (jj = 0, evalue = evalue_one = 0.0; jj < dens_count; jj++) {
							if (!ISNAN(cpo_theta[ii][jj])) {
								evalue += cpo_theta[ii][jj] * adj_weights[jj] / cpo_theta[ii][jj] / Z[ii];
								evalue_one += adj_weights[jj] / cpo_theta[ii][jj] / Z[ii];
							}
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
						if (!ISNAN(cpo_theta[ii][jj])) {
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

			Free(Z);
		} else {
			/*
			 * cpo_manual. In this case, \pi(theta|y) is compute with the cpo-data removed, hence we do not need to correct
			 */
			for (j = 0; j < compute_n; j++) {
				int ii, jj;
				double evalue, evalue2;

				ii = compute_idx[j];
				if (cpo_theta[ii]) {
					(*cpo)->value[ii] = Calloc(1, double);
					for (jj = 0, evalue = 0.0; jj < dens_count; jj++) {
						evalue += cpo_theta[ii][jj] * adj_weights[jj];
					}
					(*cpo)->value[ii][0] = evalue;
				} else {
					(*cpo)->value[ii] = NULL;
				}

				if (cpo_theta[ii]) {
					(*cpo)->pit_value[ii] = Calloc(1, double);
					(*cpo)->failure[ii] = Calloc(1, double);
					for (jj = 0, evalue = evalue2 = 0.0; jj < dens_count; jj++) {
						evalue += pit_theta[ii][jj] * adj_weights[jj];
						evalue2 += failure_theta[ii][jj] * adj_weights[jj];
					}
					evalue = TRUNCATE(evalue, 0.0, 1.0);

					(*cpo)->pit_value[ii][0] = evalue;
					(*cpo)->failure[ii][0] = evalue2;
				} else {
					(*cpo)->pit_value[ii] = NULL;
					(*cpo)->failure[ii] = NULL;
				}
			}
		}
		(*cpo)->mean_value = (*cpo)->gmean_value = 0.0;
		if (compute_n) {
			int count = 0;
			int gmean_inf = 0;

			for (j = 0; j < compute_n; j++) {
				int ii = compute_idx[j];

				if (cpo_theta[ii]) {
					(*cpo)->mean_value += *((*cpo)->value[ii]);
					if (*((*cpo)->value[ii]) > 0.0) {
						(*cpo)->gmean_value += log(*((*cpo)->value[ii]));
					} else {
						/*
						 * flag the case cpo=0
						 */
						(*cpo)->gmean_value = 0.0;
						gmean_inf = 1;
					}
					count++;
				}
			}
			if (count) {
				(*cpo)->mean_value /= (double) count;
				if (!gmean_inf) {
					(*cpo)->gmean_value = exp((*cpo)->gmean_value / (double) count);
				} else {
					/*
					 * cpo=0, hence the geometric mean is -Inf, which is more or less -DBL_MAX
					 */
					(*cpo)->gmean_value = -DBL_MAX;
				}
			} else {
				(*cpo)->mean_value = (*cpo)->gmean_value = 0.0;
			}
		}
	}

	if (po) {
		// including waic
		for (j = 0; j < compute_n; j++) {
			int ii, jj;
			double evalue, evalue2, evalue3, evalue_one;

			ii = compute_idx[j];
			if (po_theta[ii]) {
				(*po)->value[ii] = Calloc(2, double);

				for (jj = 0, evalue = evalue2 = evalue3 = evalue_one = 0.0; jj < dens_count; jj++) {
					if (po_theta[ii][jj]) {
						evalue += po_theta[ii][jj] * adj_weights[jj];
						evalue2 += po2_theta[ii][jj] * adj_weights[jj];
						evalue3 += po3_theta[ii][jj] * adj_weights[jj];
						evalue_one += adj_weights[jj];
					}
				}
				if (evalue_one) {
					(*po)->value[ii][0] = evalue / evalue_one;
					(*po)->value[ii][1] = DMAX(0.0, evalue3 / evalue_one - SQR(evalue2 / evalue_one));
				} else {
					(*po)->value[ii][0] = 0.0;
					(*po)->value[ii][1] = 0.0;
				}
			} else {
				(*po)->value[ii] = NULL;
			}
		}
	}

	if (dic) {
		double mean_deviance = 0.0, mean_deviance_sat = 0.0, deviance_mean = 0.0, deviance_mean_sat = 0.0, *x_vec = NULL;

		SET_THETA_MODE;

		/*
		 * need this for loglFunc() we need that compute is TRUE for all indices that enters loglFunc. There is no way to check this here. 
		 */
		x_vec = Calloc(graph->n, double);
		for (j = 0; j < compute_n; j++) {
			int ii = compute_idx[j];
			x_vec[ii] = (*density)[ii]->user_mean;
		}


		/*
		 * find the min length of the data contribution that cover all data points 
		 */
		int ndev = 0;
		for (j = 0; j < compute_n; j++) {
			int ii = compute_idx[j];
			if (d[ii]) {
				ndev = IMAX(ndev, ii);
			}
		}
		ndev++;

		double *e_deviance = Calloc(ndev, double), *e_deviance_sat = Calloc(ndev, double),
		    *deviance_e = Calloc(ndev, double), *deviance_e_sat = Calloc(ndev, double);

		for (j = 0; j < ndev; j++) {
			e_deviance[j] = e_deviance_sat[j] = deviance_e[j] = deviance_e_sat[j] = NAN;
		}

		for (j = 0; j < compute_n; j++) {
			double md, md_sat, dm, dm_sat, logl_saturated;
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
				double double_tmp = (double) ((*density)[ii]->user_mean);
				loglFunc(&logll, &double_tmp, 1, ii, x_vec, NULL, loglFunc_arg);
				logll *= d[ii];
				logl_saturated = d[ii] * inla_compute_saturated_loglik(ii, loglFunc, x_vec, loglFunc_arg);
				dm = -2.0 * logll;
				dm_sat = -2.0 * (logll - logl_saturated);
				md_sat = md + 2.0 * logl_saturated;
				e_deviance[ii] = md;
				e_deviance_sat[ii] = md_sat;
				deviance_e[ii] = dm;
				deviance_e_sat[ii] = dm_sat;
			} else {
				dm = md = dm_sat = md_sat = 0.0;
			}

			deviance_mean += dm;
			deviance_mean_sat += dm_sat;
			mean_deviance += md;
			mean_deviance_sat += md_sat;
		}
		Free(x_vec);

		dic->mean_of_deviance = mean_deviance;
		dic->mean_of_deviance_sat = mean_deviance_sat;
		dic->deviance_of_mean = deviance_mean;
		dic->deviance_of_mean_sat = deviance_mean_sat;
		dic->p = mean_deviance - deviance_mean;
		dic->dic = dic->p + mean_deviance;
		dic->dic_sat = dic->p + mean_deviance_sat;
		dic->n_deviance = ndev;
		dic->e_deviance = e_deviance;
		dic->e_deviance_sat = e_deviance_sat;
		dic->deviance_e = deviance_e;
		dic->deviance_e_sat = deviance_e_sat;

		if (ai_par->fp_log) {
			fprintf(ai_par->fp_log, "DIC:\n");
			fprintf(ai_par->fp_log, "\tMean of Deviance ................. %g\n", dic->mean_of_deviance);
			fprintf(ai_par->fp_log, "\tDeviance at Mean ................. %g\n", dic->deviance_of_mean);
			fprintf(ai_par->fp_log, "\tEffective number of parameters ... %g\n", dic->p);
			fprintf(ai_par->fp_log, "\tDIC .............................. %g\n", dic->dic);
			fprintf(ai_par->fp_log, "DIC (Saturated):\n");
			fprintf(ai_par->fp_log, "\tMean of Deviance ................. %g\n", dic->mean_of_deviance_sat);
			fprintf(ai_par->fp_log, "\tDeviance at Mean ................. %g\n", dic->deviance_of_mean_sat);
			fprintf(ai_par->fp_log, "\tEffective number of parameters ... %g\n", dic->p);
			fprintf(ai_par->fp_log, "\tDIC .............................. %g\n", dic->dic_sat);
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
	if (marginal_likelihood) {
		if (nhyper > 0) {
			marginal_likelihood->marginal_likelihood_gaussian_approx = 0.5 * nhyper * log(2.0 * M_PI) + log_dens_mode;
			for (i = 0; i < nhyper; i++) {
				marginal_likelihood->marginal_likelihood_gaussian_approx -=
				    0.5 * log(gsl_vector_get(eigen_values, (unsigned int) i));
			}

			if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_CCD) {
				/*
				 * in this case we integrate the 'ccd' approximation; the normal with stdev corrections. 
				 */
				marginal_likelihood->marginal_likelihood_integration = 0.5 * nhyper * log(2.0 * M_PI) + log_dens_mode;
				for (i = 0; i < nhyper; i++) {
					marginal_likelihood->marginal_likelihood_integration -=
					    0.5 * (log(gsl_vector_get(eigen_values, (unsigned int) i)) +
						   0.5 * (log(SQR(stdev_corr_pos[i])) + log(SQR(stdev_corr_neg[i]))));
				}
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
					marginal_likelihood->marginal_likelihood_integration,
					marginal_likelihood->marginal_likelihood_gaussian_approx);
			}
		} else {
			/*
			 * nhyper = 0 
			 */
			marginal_likelihood->marginal_likelihood_gaussian_approx = log_dens_mode;
			marginal_likelihood->marginal_likelihood_integration = log_dens_mode;
		}
	}

	/*
	 * compute the posterior marginals for each hyperparameter, if possible 
	 */
	if (hyper_z && density_hyper && nhyper) {

		GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_INTEGRATE_HYPERPAR, (void *) &nhyper, NULL);

		if (ai_par->interpolator == GMRFLib_AI_INTERPOLATOR_GAUSSIAN) {
			/*
			 * Just use the modal values and the stdev's found from the Hessian. 
			 */
			for (k = 0; k < nhyper; k++) {
				GMRFLib_density_create_normal(&((*density_hyper)[k]), 0.0, 1.0, theta_mode[k],
							      sqrt(inverse_hessian[k + nhyper * k]));
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
								&& (ai_par->hessian_force_diagonal || nhyper == 1)
								? GMRFLib_AI_INTERPOLATOR_GRIDSUM : GMRFLib_AI_INTERPOLATOR_CCD)
							       /*
							        *  ...if not AUTO, then use whatever is requested
							        */
							       : ai_par->interpolator);

			if (interpol != GMRFLib_AI_INTERPOLATOR_CCD && interpol != GMRFLib_AI_INTERPOLATOR_CCD_INTEGRATE &&
			    interpol != GMRFLib_AI_INTERPOLATOR_GRIDSUM) {
				double bvalue = GMRFLib_min_value(hyper_ldens, hyper_count, NULL);

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

					GMRFLib_ai_z2theta(theta_tmp, nhyper, theta_mode, &(hyper_z[k * nhyper]), sqrt_eigen_values, eigen_vectors);
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
						"\tCompute the marginal for theta[%1d] to theta[%1d] using numerical integration...\n",
						0, nhyper - 1);
				}
			}
#pragma omp parallel for private(k) num_threads(GMRFLib_openmp->max_threads_outer)
			for (k = 0; k < nhyper; k++) {
				if (!run_with_omp) {
					if (ai_par->fp_log) {
						fprintf(ai_par->fp_log, "\tCompute the marginal for theta[%1d] using numerical integration\n", k);
					}
				}
				GMRFLib_ai_marginal_for_one_hyperparamter(&((*density_hyper)[k]), k, nhyper, hyper_count, hyper_z,
									  hyper_ldens, theta_mode, sqrt_eigen_values, eigen_vectors,
									  std_stdev_theta, ai_par->dz, stdev_corr_pos,
									  stdev_corr_neg, interpol, ai_par, inverse_hessian);
			}

			if (run_with_omp) {
				if (ai_par->fp_log) {
					fprintf(ai_par->fp_log,
						"\tCompute the marginal for theta[%1d] to theta[%1d] using numerical integration... Done.\n",
						0, nhyper - 1);
				}
			}
			Free(std_stdev_theta);

			if (ai_par->fp_log) {
				fprintf(ai_par->fp_log, "Compute the marginal for the hyperparameters... done.\n");
			}
		}
	}

	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_DEFAULT, NULL, NULL);

	/*
	 * return the mode in hyperparam and in 'x'
	 */
	SET_THETA_MODE;
	if (x && x_mode) {
		memcpy(x, x_mode, graph->n * sizeof(double));
	}

	if (misc_output) {
		/*
		 * store the reordering as well. 
		 */
		if (ai_store) {
			if (ai_store->problem->sub_sm_fact.remap != NULL) {
				misc_output->len_reordering = ai_store->problem->sub_graph->n;
				misc_output->nfunc = GMRFLib_domin_get_f_count();
				misc_output->reordering = Calloc(misc_output->len_reordering, int);
				memcpy(misc_output->reordering, ai_store->problem->sub_sm_fact.remap, misc_output->len_reordering * sizeof(int));
			}
		}
	}

	/*
	 * userfunction1 
	 */
	if (GMRFLib_ai_INLA_userfunc1) {
		GMRFLib_ai_INLA_userfunc1(theta_mode, nhyper, inverse_hessian);
	}

	if (GMRFLib_ai_INLA_userfunc2) {
		/*
		 * OOPS! This loop CANNOT be run in parallel!!! 
		 */
		GMRFLib_userfunc2_arg_tp *arg = Calloc(1, GMRFLib_userfunc2_arg_tp);

		arg->stdev_corr_neg = stdev_corr_neg;
		arg->stdev_corr_pos = stdev_corr_pos;
		arg->sqrt_eigen_values = sqrt_eigen_values;
		arg->eigen_vectors = eigen_vectors;

		for (i = 0; i < GMRFLib_ai_INLA_userfunc2_n; i++) {
			GMRFLib_ai_INLA_userfunc2[i] (i, theta_mode, nhyper, inverse_hessian, (void *) arg);
		}
		Free(arg);
	}

	if (GMRFLib_ai_INLA_userfunc3) {
		/*
		 * OOPS! This loop CANNOT be run in parallel!!! 
		 */
		GMRFLib_userfunc3_arg_tp *arg = Calloc(1, GMRFLib_userfunc3_arg_tp);

		arg->stdev_corr_neg = stdev_corr_neg;
		arg->stdev_corr_pos = stdev_corr_pos;
		arg->sqrt_eigen_values = sqrt_eigen_values;
		arg->eigen_vectors = eigen_vectors;

		for (i = 0; i < GMRFLib_ai_INLA_userfunc3_n; i++) {
			GMRFLib_ai_INLA_userfunc3[i] (i, theta_mode, nhyper, inverse_hessian, (void *) arg);
		}
		Free(arg);
	}

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

		if (lin_cross) {
			for (i = 0; i < dens_count; i++) {
				Free(lin_cross[i]);
			}
			Free(lin_cross);
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
	GMRFLib_free_marginal_hidden_store(marginal_hidden_store);
	if (cpo_theta) {
		for (i = 0; i < compute_n; i++) {
			j = compute_idx[i];
			if (d[j] || ai_par->cpo_manual) {
				Free(cpo_theta[j]);
			}
		}
		Free(cpo_theta);
	}
	if (po_theta) {
		for (i = 0; i < compute_n; i++) {
			j = compute_idx[i];
			if (d[j]) {
				Free(po_theta[j]);
				Free(po2_theta[j]);
				Free(po3_theta[j]);
			}
		}
		Free(po_theta);
		Free(po2_theta);
		Free(po3_theta);
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
	if (sqrt_eigen_values) {
		gsl_vector_free(sqrt_eigen_values);
	}
	for (k = -1; (k = (int) map_strd_next(&hash_table, k)) != -1;) {
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
	}
	Free(dens);

	if (tfunc) {
		for (i = 0; i < graph->n; i++) {
			if (tfunc[i]) {
				for (j = 0; j < dens_count; j++) {
					GMRFLib_free_density(dens_transform[i][j]);
				}
				Free(dens_transform[i]);
			}
		}
	}
	Free(dens_transform);

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

	if (timer) {
		timer[3] = GMRFLib_cpu() - timer[3];
	}

	GMRFLib_LEAVE_ROUTINE;
#undef CHECK_HYPER_STORAGE
#undef CHECK_DENS_STORAGE
#undef COMPUTE_CPO_AND_DIC
#undef COMPUTE_CPO_AND_DIC_LOCAL
#undef COMPUTE
#undef COMPUTE2
#undef COMPUTE_LOCAL
#undef COMPUTE_NEFF
#undef COMPUTE_NEFF2
#undef COMPUTE_NEFF_LOCAL
#undef ADD_LINEAR_TERM
#undef ADD_LINEAR_TERM_LOCAL

	return GMRFLib_SUCCESS;
}
int GMRFLib_transform_density(GMRFLib_density_tp ** tdensity, GMRFLib_density_tp * density, GMRFLib_transform_array_func_tp * func)
{
	fprintf(stderr,
		"\n\n\nDISABLE THIS FEATURE FOR NOW, DO NOT KNOW HOW TO DO THIS WELL AT THE MOMENT. SINCE THE SCALE OF THE Xs CAN BE SO DIFFERENT, WE WILL NEED A NEW APPROACH OF HOW TO REPRESENT AND COMPUTE MIXTURES OF THESE DENSITIES.  SEE INFO ABOUT THIS AT: ISSUES\n\n.");
	GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);		       /* never enter this function */
	abort();

	/*
	 * OLD CODE THAT CAN ''FAIL'' FOR DESIGN REASONS 
	 */

	double *x, x_user, *ld, m, s;
	int len_x, i;

	if (!(density->P) || !(density->Pinv)) {
		GMRFLib_init_density(density, GMRFLib_TRUE);
	}
	GMRFLib_density_layout_x(&x, &len_x, density);	       /* in the standarized scale */

	ld = Calloc(len_x, double);
	GMRFLib_evaluate_nlogdensity(ld, x, len_x, density);

	m = func->func(density->std_mean, GMRFLib_TRANSFORM_FORWARD, func->arg, func->cov);
	s = ABS(func->func(density->std_mean, GMRFLib_TRANSFORM_DFORWARD, func->arg, func->cov)) * density->std_stdev;
	for (i = 0; i < len_x; i++) {
		x_user = GMRFLib_density_std2user(x[i], density);
		ld[i] -= log(ABS(func->func(x_user, GMRFLib_TRANSFORM_DFORWARD, func->arg, func->cov)));
		x[i] = (func->func(x_user, GMRFLib_TRANSFORM_FORWARD, func->arg, func->cov) - m) / s;
	}

	GMRFLib_density_create(tdensity, GMRFLib_DENSITY_TYPE_SCGAUSSIAN, len_x, x, ld, m, s, GMRFLib_FALSE);

	if (0) {
		printf("OLD\n");
		GMRFLib_density_printf(stdout, density);
		printf("NEW\n");
		GMRFLib_density_printf(stdout, *tdensity);
	}

	Free(x);
	Free(ld);

	return GMRFLib_SUCCESS;
}
int GMRFLib_ai_store_config(GMRFLib_ai_misc_output_tp * mo, int ntheta, double *theta, double log_posterior,
			    double log_posterior_orig, double *improved_mean, double *skewness, GMRFLib_problem_tp * gmrf_approx)
{
	if (!mo || !(mo->configs)) {
		return GMRFLib_SUCCESS;
	}

	int debug = 0;
	int id = omp_get_thread_num();

	if (!(mo->configs[id])) {
		mo->configs[id] = Calloc(1, GMRFLib_store_configs_tp);
		GMRFLib_graph_tp *g;
		GMRFLib_copy_graph(&g, gmrf_approx->sub_graph);
		mo->configs[id]->graph = g;
		if (debug) {
			printf("remapped graph\n");
			GMRFLib_print_graph(stdout, g);
		}

		int nelm;				       /* number of elements in Q; double conting */
		GMRFLib_nQelm(&nelm, gmrf_approx->sub_graph);
		mo->configs[id]->n = gmrf_approx->sub_graph->n;
		mo->configs[id]->nz = (nelm - mo->configs[id]->n) / 2 + mo->configs[id]->n;
		mo->configs[id]->ntheta = ntheta;

		if (debug) {
			printf("n nz ntheta %d %d %d\n", mo->configs[id]->n, mo->configs[id]->nz, mo->configs[id]->ntheta);
		}

		GMRFLib_constr_tp *c;
		GMRFLib_duplicate_constr(&c, gmrf_approx->sub_constr, gmrf_approx->sub_graph);	/* might or might not be mapped ???? */
		mo->configs[id]->constr = c;

		if (debug) {
			printf("constraints\n");
			GMRFLib_print_constr(stdout, c, gmrf_approx->sub_graph);
		}

		int *i, *j, ii, jj, k, kk;
		i = Calloc(mo->configs[id]->nz, int);
		j = Calloc(mo->configs[id]->nz, int);

		for (ii = k = 0; ii < g->n; ii++) {
			i[k] = ii;
			j[k] = ii;
			k++;
			for (kk = 0; kk < g->nnbs[ii]; kk++) {
				jj = g->nbs[ii][kk];
				if (ii < jj) {
					i[k] = ii;
					j[k] = jj;
					k++;
				}
			}
		}
		mo->configs[id]->i = i;
		mo->configs[id]->j = j;
		mo->configs[id]->nconfig = 0;
		mo->configs[id]->config = NULL;
	}


	mo->configs[id]->config = Realloc(mo->configs[id]->config, mo->configs[id]->nconfig + 1, GMRFLib_store_config_tp *);
	mo->configs[id]->config[mo->configs[id]->nconfig] = Calloc(1, GMRFLib_store_config_tp);

	int ii, jj, k, kk;
	double *Qinv, *Q, *mean, *imean, *skew;
	GMRFLib_graph_tp *g = gmrf_approx->sub_graph;

	Q = Calloc(mo->configs[id]->nz, double);
	for (ii = k = 0; ii < g->n; ii++) {
		Q[k++] = gmrf_approx->tab->Qfunc(ii, ii, gmrf_approx->tab->Qfunc_arg);
		for (kk = 0; kk < g->nnbs[ii]; kk++) {
			jj = g->nbs[ii][kk];
			if (ii < jj) {
				Q[k++] = gmrf_approx->tab->Qfunc(ii, jj, gmrf_approx->tab->Qfunc_arg);
			}
		}
	}

	mean = Calloc(g->n, double);
	imean = Calloc(g->n, double);
	skew = Calloc(g->n, double);
	memcpy(mean, gmrf_approx->mean_constr, g->n * sizeof(double));
	memcpy(imean, improved_mean, g->n * sizeof(double));
	memcpy(skew, skewness, g->n * sizeof(double));

	Qinv = Calloc(mo->configs[id]->nz, double);
	for (k = 0; k < mo->configs[id]->nz; k++) {
		double *tmp = GMRFLib_Qinv_get(gmrf_approx, mo->configs[id]->i[k], mo->configs[id]->j[k]);
		Qinv[k] = (tmp ? *tmp : NAN);
	}


	if (debug) {
		printf("i mean\n");
		for (k = 0; k < g->n; k++) {
			printf("%d\t %.12g\n", k, mean[k]);
		}
		printf("i\tj\tQij\n");
		for (k = 0; k < mo->configs[id]->nz; k++) {
			printf("%d\t %d\t %.12g\n", mo->configs[id]->i[k], mo->configs[id]->j[k], Q[k]);
		}
	}

	mo->configs[id]->config[mo->configs[id]->nconfig]->Q = Q;
	mo->configs[id]->config[mo->configs[id]->nconfig]->Qinv = Qinv;
	mo->configs[id]->config[mo->configs[id]->nconfig]->mean = mean;
	mo->configs[id]->config[mo->configs[id]->nconfig]->improved_mean = imean;
	mo->configs[id]->config[mo->configs[id]->nconfig]->skewness = skew;
	mo->configs[id]->config[mo->configs[id]->nconfig]->log_posterior = log_posterior;	/* may include integration weights */
	mo->configs[id]->config[mo->configs[id]->nconfig]->log_posterior_orig = log_posterior_orig;	/* do NOT include integration weights */
	if (mo->configs[id]->ntheta) {
		mo->configs[id]->config[mo->configs[id]->nconfig]->theta = Calloc(mo->configs[id]->ntheta, double);
		memcpy(mo->configs[id]->config[mo->configs[id]->nconfig]->theta, theta, mo->configs[id]->ntheta * sizeof(double));
	} else {
		mo->configs[id]->config[mo->configs[id]->nconfig]->theta = NULL;
	}
	mo->configs[id]->nconfig++;

	return GMRFLib_SUCCESS;
}
int GMRFLib_bnew(double **bnew, double *constant, int n, double *b, GMRFLib_bfunc_tp ** bfunc)
{
	/*
	 * bnew is a new alloced ptr for the new b. constant is the missing constant to be added due to b=b(theta).
	 */

	int i;
	double *bb = Calloc(n, double);
	double con = 0.0, con_add = 0.0;

	memcpy((void *) bb, (void *) b, n * sizeof(double));
	if (bfunc) {
		for (i = 0; i < n; i++) {
			if (bfunc[i]) {
				bb[i] += GMRFLib_bfunc_eval(&con_add, bfunc[i]);
				con += con_add;
			}
		}
	}

	*bnew = bb;
	*constant = -con / 2.0;

	return GMRFLib_SUCCESS;
}
int GMRFLib_free_marginal_hidden_store(GMRFLib_marginal_hidden_store_tp * m)
{
	int i;

	if (m) {
		if (m->subgraphs) {
			for (i = 0; i < m->n; i++)
				GMRFLib_free_graph(m->subgraphs[i]);
			Free(m->subgraphs);
		}
		Free(m);
	}
	return GMRFLib_SUCCESS;
}
int GMRFLib_ai_compute_lincomb(GMRFLib_density_tp *** lindens, double **cross, int nlin, GMRFLib_lc_tp ** Alin,
			       GMRFLib_ai_store_tp * ai_store, double *improved_mean)
{
	/*
	 * Compute the marginals for the linear combinations using just the Gaussians. The computations gets a bit messy since we will try to avoid dependency of n,
	 * in both a, in a^T x, and in remap(a). We only need the range of non-zero terms 'remap(a)' and the non-zero terms in 'a'.
	 */

	GMRFLib_problem_tp *problem = ai_store->problem;
	int *remap = problem->sub_sm_fact.remap;
	int i, j, k, n, nc = 0, one = 1, id;
	GMRFLib_density_tp **d;

	// I disable optimatisation as there is something going on with pardiso, in _some_ cases.
	int disable_opt = 1;


	typedef struct {
		double *v;
		int from_idx;
		int to_idx;
	} cross_tp;
	cross_tp *cross_store = NULL;

	if (GMRFLib_smtp == GMRFLib_SMTP_TAUCS || GMRFLib_smtp == GMRFLib_SMTP_BAND) {
		remap = problem->sub_sm_fact.remap;
	} else {
		// pardiso has a dynamic permutation...
		remap = problem->sub_sm_fact.PARDISO_fact->pstore->perm;
	}

	// id = GMRFLib_thread_id;
	id = GMRFLib_thread_id + omp_get_thread_num() * GMRFLib_MAX_THREADS;
	assert(problem != NULL);
	if (nlin <= 0) {
		return !GMRFLib_SUCCESS;
	}

	n = problem->n;
	d = Calloc(nlin, GMRFLib_density_tp *);
	nc = (problem->sub_constr ? problem->sub_constr->nc : 0);

	if (cross) {
		cross_store = Calloc(nlin, cross_tp);
	}

	// there is some bad designed code which require some code calling this to be run as
	// a critical region, which mess up if run with pardiso and this loop in parallel.
	// so I disable it for the moment.
	int use_pardiso = (GMRFLib_openmp->strategy == GMRFLib_OPENMP_STRATEGY_PARDISO_SERIAL ||
	                   GMRFLib_openmp->strategy == GMRFLib_OPENMP_STRATEGY_PARDISO_PARALLEL);
#pragma omp parallel for private(i, j, k) num_threads(GMRFLib_openmp->max_threads_inner) if (!use_pardiso)
	for (i = 0; i < nlin; i++) {

		int from_idx, to_idx, len, from_idx_a, to_idx_a, len_a, jj;
		double var, mean, imean, *a = NULL, *b = NULL, *v = NULL, *vv = NULL, var_corr, weight;

		if (Alin[i]->tinfo[id].first_nonzero < 0) {
			/*
			 * we know that the idx's are sorted, so its easier to find the first and last non-zero 
			 */
			// Alin[i]->tinfo[id].first_nonzero = GMRFLib_imin_value(Alin[i]->idx, Alin[i]->n);
			Alin[i]->tinfo[id].first_nonzero = Alin[i]->idx[0];
		}
		if (Alin[i]->tinfo[id].last_nonzero < 0) {
			/*
			 * we know that the idx's are sorted, so its easier to find the first and last non-zero 
			 */
			// Alin[i]->tinfo[id].last_nonzero = GMRFLib_imax_value(Alin[i]->idx, Alin[i]->n);
			Alin[i]->tinfo[id].last_nonzero = Alin[i]->idx[Alin[i]->n - 1];
		}

		if (disable_opt) {
			// disable optimization
			Alin[i]->tinfo[id].first_nonzero = 0;
			Alin[i]->tinfo[id].last_nonzero = n - 1;
		}

		from_idx_a = Alin[i]->tinfo[id].first_nonzero;
		to_idx_a = Alin[i]->tinfo[id].last_nonzero;
		assert(LEGAL(from_idx_a, n));
		assert(LEGAL(to_idx_a, n));
		len_a = to_idx_a - from_idx_a + 1;

		a = Calloc(len_a, double);
		for (j = 0; j < Alin[i]->n; j++) {
			a[Alin[i]->idx[j] - from_idx_a] = (double) Alin[i]->weight[j];
		}

		/*
		 * compute the first non-zero index (mapped) if not already there
		 */
		if (Alin[i]->tinfo[id].first_nonzero_mapped < 0) {
			int findx = n;

			for (j = 0; j < Alin[i]->n; j++) {
				k = remap[Alin[i]->idx[j]];
				findx = IMIN(findx, k);
			}
			Alin[i]->tinfo[id].first_nonzero_mapped = findx;
			Alin[i]->tinfo[id].last_nonzero_mapped = -1;
		}

		if (disable_opt) {
			// disable optimization
			Alin[i]->tinfo[id].first_nonzero_mapped = 0;
			Alin[i]->tinfo[id].last_nonzero_mapped = n - 1;
		}

		from_idx = Alin[i]->tinfo[id].first_nonzero_mapped;
		to_idx = (Alin[i]->tinfo[id].last_nonzero_mapped < 0 ? n - 1 : Alin[i]->tinfo[id].last_nonzero_mapped);

		len = to_idx - from_idx + 1;
		vv = Calloc(n, double);
		for (j = 0; j < Alin[i]->n; j++) {
			vv[remap[Alin[i]->idx[j]]] = (double) Alin[i]->weight[j];
		}
		GMRFLib_solve_l_sparse_matrix_special(vv, &(problem->sub_sm_fact), problem->sub_graph, from_idx, from_idx + len - 1, 1);
		v = Calloc(len, double);
		memcpy(v, vv + from_idx, len * sizeof(double));
		Free(vv);

		/*
		 * compute the last non-zero index (mapped) if not already there
		 */
		if (Alin[i]->tinfo[id].last_nonzero_mapped < 0) {
			Alin[i]->tinfo[id].last_nonzero_mapped = GMRFLib_find_nonzero(v, len, -1) + from_idx;
		}

		/*
		 * we do not need to map back since the innerproduct is the same in any case. It is here possible to store only
		 * the non-zero indices of v, and do the inner product over those, as they will remain the same for various
		 * calls to this function, but I don't think this is worth it.
		 */
		var = ddot_(&len, v, &one, v, &one);
		if (cross) {
			cross_store[i].from_idx = from_idx;
			cross_store[i].to_idx = to_idx;
			cross_store[i].v = v;
			v = NULL;
		} else {
			Free(v);
		}
		Free(b);

		/*
		 * the correction matrix due to linear constraints 
		 */
		var_corr = 0.0;
		if (nc) {
			for (j = 0; j < nc; j++) {
				/*
				 * w = AA^T CONSTR_M 
				 */
				double *p, *pp, w, ww;

				w = ww = 0.0;
				p = &(problem->constr_m[j * n]);
				pp = &(problem->qi_at_m[j * n]);

				for (jj = 0; jj < Alin[i]->n; jj++) {
					k = Alin[i]->idx[jj];
					weight = (double) Alin[i]->weight[jj];

					w += weight * p[k];
					ww += weight * pp[k];
				}
				var_corr += w * ww;
			}
		}

		mean = imean = 0.0;
		for (j = 0; j < Alin[i]->n; j++) {
			k = Alin[i]->idx[j];
			weight = (double) Alin[i]->weight[j];

			mean += weight * problem->mean_constr[k];
			imean += weight * improved_mean[k];
		}
		var = DMAX(DBL_EPSILON, var - var_corr);
		GMRFLib_density_create_normal(&d[i], (imean - mean) / sqrt(var), 1.0, mean, sqrt(var));
		Free(a);
	}

	if (cross) {
		/*
		 * do calculations for the E(xi*x_j) 
		 */

		*cross = Calloc(ISQR(nlin), double);

		/*
		 * need this table for OPENMP 
		 */
		int klen = (ISQR(nlin) + nlin) / 2;
		GMRFLib_lc_ij_tp *arr = Calloc(klen, GMRFLib_lc_ij_tp);
		for (k = 0, i = 0; i < nlin; i++) {
			for (j = i; j < nlin; j++) {
				arr[k].i = i;
				arr[k].j = j;
				k++;
			}
		}
		assert(k == klen);

		/*
		 * this loop is quick in any case, so no need to make do it in parallel unless we have constraints ? 
		 */
#pragma omp parallel for private(i, j, k) if(nc) num_threads(GMRFLib_openmp->max_threads_inner)
		for (k = 0; k < klen; k++) {
			i = arr[k].i;
			j = arr[k].j;

			int ij_from, ij_to;

			ij_from = IMAX(cross_store[i].from_idx, cross_store[j].from_idx);
			ij_to = IMIN(cross_store[i].to_idx, cross_store[j].to_idx);

			if (ij_from <= ij_to) {
				double *v_i = NULL, *v_j = NULL;
				int ij_len = 0;

				v_i = &(cross_store[i].v[ij_from - cross_store[i].from_idx]);
				v_j = &(cross_store[j].v[ij_from - cross_store[j].from_idx]);
				ij_len = ij_to - ij_from + 1;
				(*cross)[i + j * nlin] = ddot_(&ij_len, v_i, &one, v_j, &one);
			} else {
				(*cross)[i + j * nlin] = 0.0;
			}

			if (nc) {
				/*
				 * the correction matrix due to linear constraints 
				 */
				double correction = 0.0;
				int kk;

				for (kk = 0; kk < nc; kk++) {
					/*
					 * w = AA^T CONSTR_M 
					 */
					double *p = NULL, *pp = NULL, w, ww, weight;
					int jj, idx;

					w = ww = 0.0;
					p = &(problem->constr_m[kk * n]);
					pp = &(problem->qi_at_m[kk * n]);

					for (jj = 0; jj < Alin[i]->n; jj++) {
						idx = Alin[i]->idx[jj];
						weight = (double) Alin[i]->weight[jj];
						w += weight * p[idx];
					}
					for (jj = 0; jj < Alin[j]->n; jj++) {
						idx = Alin[j]->idx[jj];
						weight = (double) Alin[j]->weight[jj];
						ww += weight * pp[idx];
					}
					correction += w * ww;
				}
				(*cross)[i + j * nlin] += -correction;
			}
			(*cross)[j + i * nlin] = (*cross)[i + j * nlin];
		}

		Free(arr);
		for (i = 0; i < nlin; i++) {
			Free(cross_store[i].v);
		}
	}

	*lindens = d;

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
		mode = GMRFLib_max_value(logdens, *n, NULL);
		if (mode - DMAX(logdens[0], logdens[*n - 1]) < ai_par->cpo_req_diff_logdens) {
			/*
			 * this is no good... 
			 */
			*n = 0;
		}
	}

	return GMRFLib_SUCCESS;
}
double GMRFLib_ai_cpopit_integrate(double *cpo, double *pit, int idx, GMRFLib_density_tp * cpo_density, double d,
				   GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, double *x_vec)
{
	/*
	 * cpo_density is the marginal for x_idx without y_idx, density: is the marginal for x_idx with y_idx.
	 */
	int retval, compute_cpo = 1, i, k, np = GMRFLib_faster_integration_np;
	double low, dx, dxi, *xp = NULL, *xpi = NULL, *dens = NULL, *prob = NULL, *work = NULL,
	    integral = 0.0, integral2 = 0.0, w[2] = { 4.0, 2.0 }, integral_one, *loglik = NULL;
	double fail = 0.0;
	if (!cpo_density) {
		if (cpo) {
			*cpo = NAN;
		}
		if (pit) {
			*pit = NAN;
		}
		fail = 1.0;

		return fail;
	}

	retval = loglFunc(NULL, NULL, 0, idx, x_vec, NULL, loglFunc_arg);
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
		loglFunc(prob, xp, -np, idx, x_vec, NULL, loglFunc_arg);	/* no correction for 'd' here; should we? */
	} else {
		memset(prob, 0, np * sizeof(double));
	}
	loglFunc(loglik, xp, np, idx, x_vec, NULL, loglFunc_arg);
	for (i = 0; i < np; i++) {
		loglik[i] *= d;
	}

	if (0) {
		P(idx);
		FIXME("write cpo_density");
		char *ff;
		GMRFLib_sprintf(&ff, "cpo-density-%1d.dat", idx);

		FILE *fp = fopen(ff, "w");
		for (i = 0; i < np; i++) {
			fprintf(fp, "%g %g %g %g\n", xp[i], dens[i], exp(loglik[i]), prob[i]);
		}
		fclose(fp);
		printf("write file %s\n", ff);
		Free(ff);
		FIXME("info for cpo_dens");
		GMRFLib_density_printf(stdout, cpo_density);
	}

	for (i = 0; i < np; i++) {
		xp[i] = prob[i] * dens[i];		       /* reuse and redefine xp! */
		xpi[i] = exp(loglik[i]) * dens[i];	       /* reuse and redefine xpi! */
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
double GMRFLib_ai_po_integrate(double *po, double *po2, double *po3, int idx, GMRFLib_density_tp * po_density,
			       double d, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, double *x_vec)
{
	int i, k, np = GMRFLib_faster_integration_np;
	double low, dx, dxi, *xp = NULL, *xpi = NULL, *xpi3 = NULL, *xpi4 = NULL,
	    *dens = NULL, *work = NULL, integral2 = 0.0, integral3 = 0.0, integral4 = 0.0, w[2] = { 4.0, 2.0 }, integral_one, *loglik = NULL;
	double fail = 0.0;
	if (!po_density) {
		if (po) {
			*po = NAN;
		}
		if (po2) {
			*po2 = NAN;
		}
		if (po3) {
			*po3 = NAN;
		}
		fail = 1.0;

		return fail;
	}

	GMRFLib_ASSERT_RETVAL(np > 3, GMRFLib_ESNH, 0.0);
	work = Calloc(6 * np, double);
	xp = work;
	xpi = work + 1 * np;
	dens = work + 2 * np;
	loglik = work + 3 * np;
	xpi3 = work + 4 * np;
	xpi4 = work + 5 * np;

	dxi = (po_density->x_max - po_density->x_min) / (np - 1.0);
	low = GMRFLib_density_std2user(po_density->x_min, po_density);
	dx = (GMRFLib_density_std2user(po_density->x_max, po_density) - low) / (np - 1.0);

	xp[0] = low;
	xpi[0] = po_density->x_min;
	for (i = 1; i < np; i++) {
		xp[i] = xp[0] + i * dx;
		xpi[i] = xpi[0] + i * dxi;
	}
	GMRFLib_evaluate_ndensity(dens, xpi, np, po_density);

	loglFunc(loglik, xp, np, idx, x_vec, NULL, loglFunc_arg);
	for (i = 0; i < np; i++) {
		loglik[i] *= d;
	}
	for (i = 0; i < np; i++) {
		xpi[i] = exp(loglik[i]) * dens[i];	       /* reuse and redefine xpi! */
		xpi3[i] = loglik[i] * dens[i];		       /* yes, first moment */
		xpi4[i] = SQR(loglik[i]) * dens[i];	       /* yes, second moment */
	}

	integral2 = xpi[0] + xpi[np - 1];
	integral3 = xpi3[0] + xpi3[np - 1];
	integral4 = xpi4[0] + xpi4[np - 1];
	integral_one = dens[0] + dens[np - 1];
	for (i = 1, k = 0; i < np - 1; i++, k = (k + 1) % 2) {
		integral2 += w[k] * xpi[i];
		integral3 += w[k] * xpi3[i];
		integral4 += w[k] * xpi4[i];
		integral_one += w[k] * dens[i];
	}
	if (ISZERO(integral_one)) {
		fail = 1.0;
		integral2 = integral3 = integral4 = 0.0;
	} else {
		integral2 /= integral_one;
		integral3 /= integral_one;
		integral4 /= integral_one;
	}
	integral2 = DMAX(DBL_MIN, integral2);

	if (po) {
		*po = integral2;
	}
	if (po2) {
		*po2 = integral3;
	}
	if (po3) {
		*po3 = integral4;
	}
	Free(work);

	return fail;
}
double GMRFLib_ai_dic_integrate(int idx, GMRFLib_density_tp * density, double d, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, double *x_vec)
{
	/*
	 * compute the integral of -2*loglikelihood * density(x), wrt x
	 */
	double inla_compute_saturated_loglik();

	int i, k, np = GMRFLib_faster_integration_np;
	double low, dx, dxi, *xp = NULL, *xpi = NULL, *dens = NULL, *loglik = NULL, *work = NULL, integral = 0.0, w[2] =
	    { 4.0, 2.0 }, integral_one, logl_saturated;

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
	loglFunc(loglik, xp, np, idx, x_vec, NULL, loglFunc_arg);
	for (i = 0; i < np; i++) {
		loglik[i] *= d;
	}
	// logl_saturated = d * inla_compute_saturated_loglik(idx, loglFunc, x_vec, NULL, loglFunc_arg);
	logl_saturated = 0.0;

	integral = loglik[0] * dens[0] + loglik[np - 1] * dens[np - 1];
	integral_one = dens[0] + dens[np - 1];
	for (i = 1, k = 0; i < np - 1; i++, k = (k + 1) % 2) {
		integral += w[k] * loglik[i] * dens[i];
		integral_one += w[k] * dens[i];
	}
	integral = -2.0 * (integral / integral_one - logl_saturated);

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

/**
 *   \brief Free an \c GMRFLib_ai_po_tp -object created by \c GMRFLib_INLA()
 */
int GMRFLib_ai_po_free(GMRFLib_ai_po_tp * po)
{
	int i;

	if (!po) {
		return GMRFLib_SUCCESS;
	}
	if (po->value) {
		for (i = 0; i < po->n; i++) {
			if (po->value[i]) {
				Free(po->value[i]);
			}
		}
		Free(po->value);
	}
	Free(po);

	return GMRFLib_SUCCESS;
}

int GMRFLib_ai_add_Qinv_to_ai_store(GMRFLib_ai_store_tp * ai_store)
{
	if (!ai_store || !(ai_store->problem)) {
		return GMRFLib_SUCCESS;
	}

	if (!ai_store->problem->sub_inverse) {
		int i, n, idum;

#pragma omp parallel for private(idum) num_threads(GMRFLib_openmp->max_threads_outer)
		for (idum = 0; idum < 1; idum++) {
			GMRFLib_Qinv(ai_store->problem, GMRFLib_QINV_NEIGB);
		}
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
					      double *hyper_ldens, double *theta_mode, gsl_vector * sqrt_eigen_values,
					      gsl_matrix * eigen_vectors, double *std_stdev_theta, double dz,
					      double *stdev_corr_pos, double *stdev_corr_neg,
					      GMRFLib_ai_interpolator_tp interpolator, GMRFLib_ai_param_tp * ai_par, double *covmat)
{
#define COV(i, j)  covmat[ (i) + (j)*nhyper ]
#define NEXTRA 19
	int i, j;
	double *points = NULL, *ldens_values, *theta_max, *theta_min, sd;
	double extra_points[] = { -15.0, -10.0, -7.0, -5.0, -3.0, -2.0, -1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0 };
	int npoints;

	GMRFLib_ENTER_ROUTINE;
	assert(sizeof(extra_points) / sizeof(double) == NEXTRA);

	// 
	// Do not use this option, better off just doing what's done below. This is especially bad if the int_strategy = CCD.
	// 
	// if (idx == 0 && nhyper == 1) 
	// sd = std_stdev_theta[idx];
	// GMRFLib_density_create(density, GMRFLib_DENSITY_TYPE_SCGAUSSIAN, hyper_count, hyper_z, hyper_ldens, theta_mode[idx], sd, 
	// 
	// GMRFLib_TRUE);
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
			GMRFLib_ai_z2theta(theta_tmp, nhyper, theta_mode, &(hyper_z[i * nhyper]), sqrt_eigen_values, eigen_vectors);
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
	} else if (interpolator == GMRFLib_AI_INTERPOLATOR_CCD) {
		/*
		 * this is the new analytic approximation 
		 */

		GMRFLib_ai_integrator_arg_tp *arg = Calloc(1, GMRFLib_ai_integrator_arg_tp);

		arg->nhyper = nhyper;
		arg->idx = -1;
		arg->return_log = GMRFLib_TRUE;		       /* return log(density), yes. */
		arg->hyper_count = hyper_count;
		arg->hyper_z = hyper_z;
		arg->hyper_ldens = hyper_ldens;
		arg->theta_mode = theta_mode;
		arg->sqrt_eigen_values = sqrt_eigen_values;
		arg->eigen_vectors = eigen_vectors;
		arg->z = Calloc(nhyper, double);
		arg->theta = Calloc(nhyper, double);
		arg->stdev_corr_pos = stdev_corr_pos;
		arg->stdev_corr_neg = stdev_corr_neg;
		arg->dz = dz;
		arg->interpolator = interpolator;

		npoints = 71;
		double theta_fixed, *x = NULL, *xx = NULL, *xxx = NULL;

		GMRFLib_ghq_abscissas(&xx, npoints);
		xxx = Calloc(npoints + NEXTRA, double);
		memcpy(xxx, xx, npoints * sizeof(double));
		memcpy(xxx + npoints, extra_points, NEXTRA * sizeof(double));

		npoints += NEXTRA;
		qsort((void *) xxx, (size_t) npoints, sizeof(double), GMRFLib_dcmp);
		xxx[0] = DMIN(xxx[0], -GMRFLib_DENSITY_INTEGRATION_LIMIT);
		xxx[npoints - 1] = DMAX(xxx[npoints - 1], GMRFLib_DENSITY_INTEGRATION_LIMIT);
		x = Calloc(nhyper, double);
		ldens_values = Calloc(npoints, double);

		for (i = 0; i < npoints; i++) {
			theta_fixed = theta_mode[idx] + std_stdev_theta[idx] * xxx[i];
			for (j = 0; j < nhyper; j++) {
				if (j != idx) {
					x[j] = theta_mode[j] + (COV(idx, j) / COV(idx, idx)) * (theta_fixed - theta_mode[idx]);
				} else {
					x[j] = theta_fixed;
				}
			}
			ldens_values[i] = GMRFLib_ai_integrator_func(nhyper, x, arg);
			// printf("i %d x %g ldens %g\n", i, x[idx], ldens_values[i]);
		}

		GMRFLib_density_create(density, GMRFLib_DENSITY_TYPE_SCGAUSSIAN, npoints, xxx, ldens_values, theta_mode[idx],
				       std_stdev_theta[idx], GMRFLib_TRUE);

		Free(x);
		Free(xxx);
		Free(ldens_values);
		Free(points);
		Free(arg->z);
		Free(arg->theta);
		Free(arg);
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

		double limit = 6.0;			       /* defines the integration limit */
		for (i = j = 0; i < nhyper; i++) {
			if (i != idx) {
				theta_min[j] = theta_mode[i] - limit * stdev_corr_neg[i] * std_stdev_theta[i];
				theta_max[j] = theta_mode[i] + limit * stdev_corr_pos[i] * std_stdev_theta[i];
				j++;
			}
			theta_min_all[i] = theta_mode[i] - limit * stdev_corr_neg[i] * std_stdev_theta[i];
			theta_max_all[i] = theta_mode[i] + limit * stdev_corr_pos[i] * std_stdev_theta[i];
		}

		int tmax = GMRFLib_MAX_THREADS;
		GMRFLib_ai_integrator_arg_tp **arg = Calloc(tmax, GMRFLib_ai_integrator_arg_tp *);

		for (i = 0; i < tmax; i++) {
			arg[i] = Calloc(tmax, GMRFLib_ai_integrator_arg_tp);

			arg[i]->nhyper = nhyper;
			arg[i]->idx = idx;
			arg[i]->return_log = GMRFLib_FALSE;
			arg[i]->hyper_count = hyper_count;
			arg[i]->hyper_z = hyper_z;
			arg[i]->hyper_ldens = hyper_ldens;
			arg[i]->theta_mode = theta_mode;
			arg[i]->sqrt_eigen_values = sqrt_eigen_values;
			arg[i]->eigen_vectors = eigen_vectors;
			arg[i]->z = Calloc(nhyper, double);
			arg[i]->theta = Calloc(nhyper, double);

			arg[i]->stdev_corr_pos = stdev_corr_pos;
			arg[i]->stdev_corr_neg = stdev_corr_neg;
			arg[i]->dz = dz;
			arg[i]->interpolator = interpolator;
		}

		/*
		 * we need to bound the maximum function evaluations, otherwise it can just go on forever, especially for _linear
		 * and _quadratic interpolation. seems like they produce to `rough' integrands... 
		 */
		unsigned int max_eval = (unsigned int) ai_par->numint_max_fn_eval;

#pragma omp parallel for private(i) num_threads(GMRFLib_openmp->max_threads_outer)
		for (i = 0; i < npoints + NEXTRA; i++) {
			int retval;
			double abs_err = ai_par->numint_abs_err, rel_err = ai_par->numint_rel_err, value, err;
			int thread = omp_get_thread_num();

			if (i < npoints) {
				arg[thread]->theta_fixed = theta_min_all[idx] + i * (theta_max_all[idx] - theta_min_all[idx]) / (npoints - 1.0);
				points[i] = (arg[thread]->theta_fixed - theta_mode[idx]) / sd;
			} else {
				int ii = i - npoints;
				arg[thread]->theta_fixed = theta_mode[idx] + extra_points[ii] * sd;
				points[i] = extra_points[ii];
			}

			retval = adapt_integrate(GMRFLib_ai_integrator_func, arg[thread], (unsigned int) (nhyper - 1),
						 (const double *) theta_min, (const double *) theta_max, max_eval, abs_err, rel_err, &value, &err);
			value = DMAX(DBL_MIN, value);
			ldens_values[i] = log(value);
			if (retval) {
#pragma omp critical
				{
					fprintf(stderr, "\n\tGMRFLib_ai_marginal_for_one_hyperparamter: warning:\n");
					fprintf(stderr, "\t\tMaximum number of function evaluations is reached\n");
					fprintf(stderr, "\t\ti=%1d, point=%g, thread=%1d\n", i, points[i], thread);
				}
			}
			// printf("eval points %.12g ldens %.12g REF %.12g\n", points[i], ldens_values[i], -0.5*SQR(points[i]));
		}

		GMRFLib_density_create(density, GMRFLib_DENSITY_TYPE_SCGAUSSIAN, npoints + NEXTRA, points, ldens_values,
				       theta_mode[idx], sd, GMRFLib_TRUE);

		for (i = 0; i < tmax; i++) {
			Free(arg[i]->z);
			Free(arg[i]->theta);
			Free(arg[i]);
		}
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
#undef COV
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
		memcpy(a->theta, x, a->nhyper * sizeof(double));	/* ndim == nhyper */
	}

	GMRFLib_ai_theta2z(a->z, a->nhyper, a->theta_mode, a->theta, a->sqrt_eigen_values, a->eigen_vectors);
	switch (a->interpolator) {
	case GMRFLib_AI_INTERPOLATOR_CCD:
	case GMRFLib_AI_INTERPOLATOR_CCD_INTEGRATE:
		val = GMRFLib_interpolator_ccd(a->nhyper, a->hyper_count, a->z, a->hyper_z, a->hyper_ldens, (void *) a);
		break;
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
													 a->z, a->hyper_z,
													 a->hyper_ldens, (void *) &(a->dz)));
		if (a->stdev_corr_pos && a->stdev_corr_neg) {
			printf(" CCD %.10g\n", GMRFLib_interpolator_ccd(a->nhyper, a->hyper_count, a->z, a->hyper_z, a->hyper_ldens, (void *) a));
		} else {
			printf("\n");
		}
	}

	return (a->return_log ? val : exp(val));
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
GMRFLib_ai_store_tp *GMRFLib_duplicate_ai_store(GMRFLib_ai_store_tp * ai_store, int skeleton, int copy_ptr, int copy_pardiso_ptr)
{
	/*
	 * duplicate AI_STORE. 'skeleton' only duplicate 'required' features. 'copy_ptr' only copies pointers to some objects known to be 'read only'
	 */
	int id = 0;
	GMRFLib_meminfo_thread_id = id;

#define DUPLICATE(name, len, tp, skeleton_)				\
	if (1) {							\
		if (ai_store->name && len && !skeleton_){		\
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

//#pragma omp parallel sections
	{
//#pragma omp section
		{
			GMRFLib_meminfo_thread_id = id;
			new_ai_store->store = GMRFLib_duplicate_store(ai_store->store, skeleton, copy_ptr, copy_pardiso_ptr);
		}
//#pragma omp section
		{
			GMRFLib_meminfo_thread_id = id;
			new_ai_store->problem = GMRFLib_duplicate_problem(ai_store->problem, skeleton, copy_ptr, copy_pardiso_ptr);
			COPY(nidx);
			COPY(neff);
			COPY(nd);

			DUPLICATE(mode, n, double, 0);
			DUPLICATE(aa, n, double, skeleton);
			DUPLICATE(bb, n, double, skeleton);
			DUPLICATE(cc, n, double, skeleton);
			DUPLICATE(stdev, n, double, skeleton);
			DUPLICATE(correction_term, n, double, skeleton);
			DUPLICATE(derivative3, n, double, skeleton);
			DUPLICATE(correction_idx, n, int, skeleton);
			DUPLICATE(d_idx, nd, int, 0);
		}
	}


	GMRFLib_meminfo_thread_id *= -1;
	char *tmp = Calloc(1, char);
	Free(tmp);

	GMRFLib_LEAVE_ROUTINE;
	return new_ai_store;

#undef DUPLICATE
#undef COPY
}
GMRFLib_ai_store_tp *GMRFLib_assign_ai_store(GMRFLib_ai_store_tp * to, GMRFLib_ai_store_tp * from)
{
	/*
	 * set contents of TO = FROM
	 */

#define ASSIGN(name) to->name = from->name

	GMRFLib_ENTER_ROUTINE;
	if (!to || !from) {
		GMRFLib_LEAVE_ROUTINE;
		return NULL;
	}

	ASSIGN(aa);
	ASSIGN(bb);
	ASSIGN(cc);
	ASSIGN(correction_idx);
	ASSIGN(correction_term);
	ASSIGN(d_idx);
	ASSIGN(derivative3);
	ASSIGN(mode);
	ASSIGN(nc_orig);
	ASSIGN(nd);
	ASSIGN(neff);
	ASSIGN(nidx);
	ASSIGN(problem);
	ASSIGN(stdev);
	ASSIGN(store);

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
#undef ASSIGN
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
	half_len = (int) (((sqrt(2.0 * p->diff_log_dens) + 1.0) / ai_par->dz + 2.0) *
			  (nhyper == 1 ? 6.0 : (nhyper == 2 ? 4.0 : (nhyper == 3 ? 2.0 : 1.0))));
	len = (2 * half_len + 1);
	p->nconfig = (size_t) pow((double) len, (double) p->nhyper);
	p->configurations = Calloc((size_t) (p->nconfig * p->nhyper), GMRFLib_short_int);
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
			p->configurations[k + j] = (GMRFLib_short_int) izz[j];
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
	if (GMRFLib_MAX_THREADS > 1) {
		qsort(p->configurations, p->nconfig, p->nhyper * sizeof(GMRFLib_short_int), GMRFLib_pool_cmp);
	} else {
		/*
		 * alternative sorting: _pool_cmp1: seems like _pool_cmp runs faster (better wrt 'reject') 
		 */
		qsort(p->configurations, p->nconfig, p->nhyper * sizeof(GMRFLib_short_int), GMRFLib_pool_cmp);
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
	const GMRFLib_short_int *ia, *ib;
	int i, larger, dist_a, dist_b;

	GMRFLib_ASSERT(pool_nhyper > 0, GMRFLib_ESNH);

	ia = (const GMRFLib_short_int *) a;
	ib = (const GMRFLib_short_int *) b;
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
	const GMRFLib_short_int *ia, *ib;
	int i, larger, eq;

	GMRFLib_ASSERT(pool_nhyper > 0, GMRFLib_ESNH);

	ia = (const GMRFLib_short_int *) a;
	ib = (const GMRFLib_short_int *) b;
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
			if ((ISNAN(logdens) || ISINF(logdens)) || -logdens > pool->diff_log_dens) {
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
double GMRFLib_bfunc_eval(double *constant, GMRFLib_bfunc_tp * bfunc)
{
	// evaluate bfunc: b[idx] = sum_i Q[idx,i]*mean[i]. 'con' is the contribution to the constant m'(Q+c)m.

#define MAPIDX(_idx, _d) MOD(MOD(_idx, (_d)->n * (_d)->ngroup), (_d)->n)

	if (bfunc == NULL || bfunc->bdef == NULL || bfunc->idx < 0) {
		return 0.0;
	}

	double b = 0.0;
	int i, j, idx = bfunc->idx;
	GMRFLib_bfunc2_tp *d = bfunc->bdef;

	// fprintf(stderr, "idx %d mapidx %d n %d nr %d ng %d\n", idx, MAPIDX(idx, d), d->n, d->nreplicate, d->ngroup);

	b = (d->diagonal + d->Qfunc(idx, idx, d->Qfunc_arg)) * d->mfunc(MAPIDX(idx, d), d->mfunc_arg);
	for (i = 0; i < d->graph->nnbs[idx]; i++) {
		j = d->graph->nbs[idx][i];
		b += d->Qfunc(idx, j, d->Qfunc_arg) * d->mfunc(MAPIDX(j, d), d->mfunc_arg);
	}

	*constant = b * d->mfunc(MAPIDX(idx, d), d->mfunc_arg);

#undef MAPIDX
	return b;
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
