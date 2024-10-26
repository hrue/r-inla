
/* approx-inference.c
 *
 * Copyright (C) 2006-2024 Havard Rue
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

#include <stdint.h>
#include <assert.h>
#include <time.h>
#include <stddef.h>
#include <float.h>
#include <math.h>
#include <strings.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"
#include "GMRFLib/hashP.h"

int GMRFLib_default_ai_param(GMRFLib_ai_param_tp **ai_par)
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

	(*ai_par)->linear_correction = GMRFLib_AI_LINEAR_CORRECTION_CENTRAL_DIFFERENCE;
	(*ai_par)->linear_correction = GMRFLib_AI_LINEAR_CORRECTION_OFF;
	(*ai_par)->linear_correction = GMRFLib_AI_LINEAR_CORRECTION_FAST;	/* to match the ->fast mode above */

	/*
	 * none of these are used, but they are the defaults if the user wants improved approximations 
	 */
	(*ai_par)->n_points = 9;			       /* how many points to evaluate */
	(*ai_par)->step_len = GSL_ROOT4_DBL_EPSILON;	       /* If the derivaties has to be computed numerically */
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

	(*ai_par)->optimise_smart = GMRFLib_FALSE;
	(*ai_par)->optimise_use_directions = GMRFLib_FALSE;
	(*ai_par)->optimise_use_directions_m = NULL;
	(*ai_par)->optimiser = GMRFLib_AI_OPTIMISER_GSL;
	(*ai_par)->restart = 0;
	(*ai_par)->gsl_tol = 0.1;
	(*ai_par)->gsl_epsg = 0.005;
	(*ai_par)->gsl_epsf = pow(0.005, 1.5);		       /* this is the default relationship used in R-INLA */
	(*ai_par)->gsl_epsx = 0.005;
	(*ai_par)->gsl_step_size = 1.0;
	(*ai_par)->mode_restart = 1;
	(*ai_par)->mode_fixed = 0;
	(*ai_par)->mode_use_mode = 0;
	(*ai_par)->parallel_linesearch = 0;

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
	(*ai_par)->b_strategy = 1;			       /* keep */

	/*
	 * default is no correction of any kind
	 */
	(*ai_par)->vb_enable = 0;
	(*ai_par)->vb_strategy = GMRFLib_AI_VB_MEAN;
	(*ai_par)->vb_verbose = 0;
	(*ai_par)->vb_iter_max = 5;
	(*ai_par)->vb_f_enable_limit_mean = 20;
	(*ai_par)->vb_f_enable_limit_variance = 10;
	(*ai_par)->vb_f_enable_limit_mean_max = 1024;
	(*ai_par)->vb_f_enable_limit_variance_max = 768;
	(*ai_par)->vb_nodes_mean = NULL;
	(*ai_par)->vb_nodes_variance = NULL;

	(*ai_par)->hessian_correct_skewness_only = 0;

	return GMRFLib_SUCCESS;
}

int GMRFLib_ai_param_duplicate(GMRFLib_ai_param_tp **ai_par_new, GMRFLib_ai_param_tp *ai_par)
{
	if (!ai_par) {
		*ai_par_new = NULL;
		return GMRFLib_SUCCESS;
	}

	*ai_par_new = Calloc(1, GMRFLib_ai_param_tp);
	Memcpy((void *) *ai_par_new, (void *) ai_par, sizeof(GMRFLib_ai_param_tp));

	if (ai_par->optimise_use_directions_m) {
		(*ai_par_new)->optimise_use_directions_m = GMRFLib_gsl_duplicate_matrix(ai_par->optimise_use_directions_m);
	}

	if (ai_par->adapt_strategy) {
		(*ai_par_new)->adapt_strategy = Calloc(ai_par->adapt_len, GMRFLib_ai_strategy_tp);
		Memcpy((void *) (*ai_par_new)->adapt_strategy, (void *) ai_par->adapt_strategy, sizeof(GMRFLib_ai_strategy_tp) * ai_par->adapt_len);
	}

	return GMRFLib_SUCCESS;
}
int GMRFLib_ai_param_free(GMRFLib_ai_param_tp *ai_par)
{
	if (ai_par) {
		if (ai_par->adapt_len > 0 && ai_par->adapt_strategy) {
			Free(ai_par->adapt_strategy);
		}
		Free(ai_par);
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_print_ai_param(FILE *fp, GMRFLib_ai_param_tp *ai_par)
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
	fprintf(fp, "\t\tOptimise: try to be smart: %s\n", (ai_par->optimise_smart ? "Yes" : "No"));
	fprintf(fp, "\t\tOptimise: use directions: %s\n", (ai_par->optimise_use_directions ? "Yes" : "No"));
	fprintf(fp, "\t\tMode restart: %s\n", (ai_par->mode_restart ? "Yes" : "No"));
	fprintf(fp, "\t\tMode fixed: %s\n", (ai_par->mode_fixed ? "Yes" : "No"));
	fprintf(fp, "\t\tMode use_mode: %s\n", (ai_par->mode_use_mode ? "Yes" : "No"));
	fprintf(fp, "\t\tparallel linesearch [%1d]\n", ai_par->parallel_linesearch);

	gsl_matrix *m = ai_par->optimise_use_directions_m;
	if (m) {
		fprintf(fp, "\tDirection matrix:\n");
		for (size_t i = 0; i < m->size1; i++) {
			fprintf(fp, "\t\t");
			for (size_t j = 0; j < m->size2; j++) {
				fprintf(fp, " %6.3f", gsl_matrix_get(m, i, j));
			}
			fprintf(fp, "\n");
		}
	}

	fprintf(fp, "\tGaussian approximation:\n");
	fprintf(fp, "\t\ttolerance_func = %.6g\n", ai_par->optpar_abserr_func);
	fprintf(fp, "\t\ttolerance_step = %.6g\n", ai_par->optpar_abserr_step);
	fprintf(fp, "\t\toptpar_fp = %" PRIxPTR "\n", (uintptr_t) (ai_par->optpar_fp));
	fprintf(fp, "\t\toptpar_nr_step_factor = %.6g\n", ai_par->optpar_nr_step_factor);

	fprintf(fp, "\tGaussian data: %s\n", (GMRFLib_gaussian_data ? "Yes" : "No"));

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
	if (ai_par->improved_simplified_laplace) {
		fprintf(fp, "\t\t\tUse an experimental improved simplified.laplace\n");
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
	GMRFLib_design_print(fp, ai_par->int_design);

	fprintf(fp, "\t\tf0 (CCD only):\t %.3f\n", ai_par->f0);
	fprintf(fp, "\t\tdz (GRID only):\t %.3f\n", ai_par->dz);
	fprintf(fp, "\t\tAdjust weights (GRID only):\t %s\n", (ai_par->adjust_weights == GMRFLib_FALSE ? "Off" : "On"));
	fprintf(fp, "\t\tDifference in log-density limit (GRID only):\t %.3f\n", ai_par->diff_log_dens);
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
	fprintf(fp, "\t\tStrategy for the linear term [%s]\n", (ai_par->b_strategy == 0 ? "Skip" : "Keep"));

	if (show_expert_options) {
		/*
		 * expert options goes here 
		 */
		fprintf(fp, "\tCPO manual calculation[%s]\n", (ai_par->cpo_manual ? "Yes" : "No"));
	}

	if (ai_par->vb_enable) {
		fprintf(fp, "\tVB correction is [Enabled]\n");
		fprintf(fp, "\t\tstrategy                    = [%s]\n", (ai_par->vb_strategy == GMRFLib_AI_VB_MEAN ? "mean" :
									 (ai_par->vb_strategy ==
									  GMRFLib_AI_VB_VARIANCE ? "mean and variance" : "UNKNOWN")));
		fprintf(fp, "\t\tverbose                     = [%s]\n", (ai_par->vb_verbose ? "Yes" : "No"));
		fprintf(fp, "\t\tf_enable_limit_mean         = [%1d]\n", ai_par->vb_f_enable_limit_mean);
		fprintf(fp, "\t\tf_enable_limit_var          = [%1d]\n", ai_par->vb_f_enable_limit_variance);
		fprintf(fp, "\t\tf_enable_limit_mean_max     = [%1d]\n", ai_par->vb_f_enable_limit_mean_max);
		fprintf(fp, "\t\tf_enable_limit_variance_max = [%1d]\n", ai_par->vb_f_enable_limit_variance_max);
		fprintf(fp, "\t\titer_max                    = [%1d]\n", ai_par->vb_iter_max);
		fprintf(fp, "\t\temergency                   = [%.2f]\n", ai_par->vb_emergency);
		fprintf(fp, "\t\thessian_update              = [%1d]\n", ai_par->vb_hessian_update);
		fprintf(fp, "\t\thessian_strategy            = [%s]\n", VB_HESSIAN_STRATEGY_NAME(ai_par->vb_hessian_strategy));
	} else {
		fprintf(fp, "\tVB-correction is [Disabled]\n");
	}

	fprintf(fp, "\tMisc options: \n");
	fprintf(fp, "\t\tHessian correct skewness only [%1d]\n", ai_par->hessian_correct_skewness_only);

	return GMRFLib_SUCCESS;
}

int GMRFLib_ai_marginal_hyperparam(int thread_id,
				   double *logdens,
				   double *x, double *b, double *c, double *mean, double *d, int *fl,
				   GMRFLib_logl_tp *loglFunc, void *loglFunc_arg,
				   GMRFLib_graph_tp *graph, GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_arg,
				   GMRFLib_constr_tp *constr, GMRFLib_ai_param_tp *ai_par, GMRFLib_ai_store_tp *ai_store,
				   GMRFLib_preopt_tp *preopt, GMRFLib_idx_tp *d_idx)
{
	/*
	 * return the unnormalised log marginal density for the hyperparamers in `logdens'. note that we assume that 'b' and 'mean' are constants. when this is not
	 * the case, the missing term is added elsewhere. 
	 */

	double ldens = 0.0;
	int n, free_ai_par = 0;
	int Npred = (preopt ? preopt->Npred : graph->n);
	ai_store->Npred = Npred;

	static int *nnr_step_factor_first_time_only = NULL;
	if (!nnr_step_factor_first_time_only) {
#pragma omp critical (Name_4afa98d4e0e7cc3aec97acb922c2fa7fb65a660f)
		{
			if (!nnr_step_factor_first_time_only) {
				nnr_step_factor_first_time_only = Calloc(GMRFLib_CACHE_LEN(), int);
				for (int i = 0; i < GMRFLib_CACHE_LEN(); i++) {
					nnr_step_factor_first_time_only[i] = 1;
				}
			}
		}
	}
	int idx = 0;
	GMRFLib_CACHE_SET_ID(idx);

	int *nr_step_factor_first_time_only = &(nnr_step_factor_first_time_only[idx]);

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
	optpar->abserr_step = ai_par->optpar_abserr_step;
	optpar->fp = ai_par->optpar_fp;

	if (ai_par->optpar_nr_step_factor < 0) {
		/*
		 *  then do this only the first time for each thread. This would improve initial values
		 */
		if (*nr_step_factor_first_time_only) {
			optpar->nr_step_factor = -ai_par->optpar_nr_step_factor;
			*nr_step_factor_first_time_only = 0;
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
	Free(ai_store->derivative4);
	Free(ai_store->aa);
	Free(ai_store->bb);
	Free(ai_store->cc);
	ai_store->aa = Calloc(Npred, double);
	ai_store->bb = Calloc(Npred, double);
	ai_store->cc = Calloc(Npred, double);

	/*
	 * first compute the GMRF-approximation 
	 */
	if (ai_store->mode) {
		GMRFLib_EWRAP1(GMRFLib_init_GMRF_approximation_store__intern
			       (thread_id, &problem, ai_store->mode, b, c, mean, d, fl, loglFunc, loglFunc_arg, graph, Qfunc,
				Qfunc_arg, constr, optpar, blockpar, ai_store->store, ai_store->aa, ai_store->bb, ai_store->cc,
				ai_par->cmin, ai_par->b_strategy, 0, preopt, d_idx));
	} else {
		GMRFLib_EWRAP1(GMRFLib_init_GMRF_approximation_store__intern
			       (thread_id, &problem, x, b, c, mean, d, fl, loglFunc, loglFunc_arg, graph, Qfunc, Qfunc_arg,
				constr, optpar, blockpar, ai_store->store, ai_store->aa, ai_store->bb, ai_store->cc,
				ai_par->cmin, ai_par->b_strategy, 0, preopt, d_idx));
	}
	GMRFLib_ASSERT(problem, GMRFLib_EOPTNR);

	/*
	 * if store, then store the mode to use as the initial point at later calls 
	 */
	Free(ai_store->mode);
	ai_store->mode = Calloc(n, double);
	Memcpy(ai_store->mode, problem->mean_constr, n * sizeof(double));

	if (mean == NULL) {
		/*
		 * Here we use the joint expression and take advantage of that we have already evaluated the log-likelihood in the mode.
		 * 
		 * the expression below is not correct if mean != 0, so then we need a term -1/2 mu*(Q+C)mu as well. so therefore the if (mean==NULL) test. 
		 */
		Memset(problem->sample, 0, n * sizeof(double));
		GMRFLib_evaluate(problem);

		double A = GMRFLib_dsum(Npred, ai_store->aa);
		*logdens = A - problem->sub_logdens;
	} else {
		/*
		 * then evaluate it in the mode (ie the mean) to get the log-density 
		 */
		Memcpy(problem->sample, problem->mean_constr, n * sizeof(double));
		GMRFLib_EWRAP1(GMRFLib_evaluate(problem));
		GMRFLib_ai_log_posterior(thread_id, &ldens, problem->mean_constr, b, c, mean, d, loglFunc, loglFunc_arg, graph, Qfunc, Qfunc_arg,
					 constr);

		*logdens = ldens - problem->sub_logdens;
	}

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

int GMRFLib_free_ai_store(GMRFLib_ai_store_tp *ai_store)
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
		GMRFLib_idx_free(ai_store->d_idx);
		Free(ai_store->correction_term);
		Free(ai_store->correction_idx);
		Free(ai_store->derivative3);
		Free(ai_store->derivative4);
		Free(ai_store);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_ai_z2theta(double *theta, int nhyper, double *theta_mode, double *z, gsl_vector *sqrt_eigen_values, gsl_matrix *eigen_vectors)
{
	/*
	 * compute new theta-values for given vector of z (which is N(0,I)), using the relationship
	 * 
	 * theta = theta_mode + eigen_vectors * diag(1/sqrt(eigen_values)) * z 
	 */

	size_t i, j;
	double tmp, v_ij, *u = NULL;

	assert(z);
	assert(nhyper > 0);
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

int GMRFLib_ai_theta2z(double *z, int nhyper, double *theta_mode, double *theta, gsl_vector *sqrt_eigen_values, gsl_matrix *eigen_vectors)
{
	/*
	 * compute z-values for given vector of theta, using the relationship
	 * 
	 * theta = theta_mode + eigen_vectors * diag(1/sqrt_eigen_values) * z 
	 */

	size_t i, j;
	double tmp, *u = NULL;

	assert(nhyper > 0);
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

int GMRFLib_init_GMRF_approximation_store__intern(int thread_id,
						  GMRFLib_problem_tp **problem, double *x, double *b, double *c, double *mean,
						  double *d, int *fl, GMRFLib_logl_tp *loglFunc, void *loglFunc_arg,
						  GMRFLib_graph_tp *graph, GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_arg,
						  GMRFLib_constr_tp *constr, GMRFLib_optimize_param_tp *optpar,
						  GMRFLib_blockupdate_param_tp *blockupdate_par, GMRFLib_store_tp *store,
						  double *aa, double *bb, double *cc, double cmin, int b_strategy, int nested,
						  GMRFLib_preopt_tp *preopt, GMRFLib_idx_tp *d_idx)
{
	// this is implicitely assumed
	assert(mean == NULL);
	assert(x);

	if (0) {
		FIXME("Enter with x=");
		for (int i = 0; i < graph->n; i++)
			printf(" %.4f", x[i]);
		printf("\n");
	}

	/*
	 * This is copy of the original routine but with optional arguments 
	 */

	int free_b = 0, free_c = 0, free_mean = 0, free_d = 0, free_blockpar = 0, free_aa = 0, free_bb = 0, free_cc = 0, n;
	int Npred = (GMRFLib_inla_mode == GMRFLib_MODE_COMPACT ? preopt->Npred : graph->n);
	double *mode = NULL;

#define FREE_ALL if (1) { if (free_b) Free(b); if (free_c) Free(c); if (free_d) Free(d); \
		if (free_mean) Free(mean); if (free_blockpar) Free(blockupdate_par); if (free_aa) Free(aa); if (free_bb) Free(bb); \
		if (free_cc) Free(cc); if (free_d_idx) GMRFLib_idx_free(d_idx);  }

	n = graph->n;
	if (n == 0) {
		*problem = NULL;
		return GMRFLib_SUCCESS;
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
		aa = Calloc(Npred, double);
	}
	if (!bb) {
		free_bb = 1;
		bb = Calloc(Npred, double);
	}
	if (!cc) {
		free_cc = 1;
		cc = Calloc(Npred, double);
	}
	mode = Calloc(n, double);

	Memcpy(mode, x, n * sizeof(double));

	int free_d_idx = 0;
	if (d_idx == NULL) {
		GMRFLib_idx_create(&d_idx);
		for (int i = 0; i < (preopt ? preopt->Npred : graph->n); i++) {
			if (d[i]) {
				GMRFLib_idx_add(&d_idx, i);
			}
		}
		free_d_idx = 1;
	}

	/*
	 * the NEW implementation which do optimisation and GMRF_approximation in the same operation. this is tailored for INLA of'course, and only ment
	 * for that. 
	 */

	GMRFLib_ENTER_ROUTINE;

	if (optpar && optpar->fp)
		fprintf(optpar->fp, "\n[%1d] Computing GMRF approximation\n------------------------------\n", omp_get_thread_num());
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

	Memcpy(mode_initial, mode, n * sizeof(double));	       /* store the starting value */

	double *bcoof = NULL, *ccoof = NULL, *linear_predictor = NULL;
	int free_linear_predictor = 0;

	bcoof = Calloc(Npred, double);
	ccoof = Calloc(Npred, double);

	for (iter = 0; iter < itmax; iter++) {

		if (0) {
			for (int i = 0; i < n; i++) {
				printf("i mode %d %f \n", i, mode[i]);
			}
		}

		Memset(aa, 0, Npred * sizeof(double));
		Memset(bb, 0, Npred * sizeof(double));
		Memset(cc, 0, Npred * sizeof(double));

		if (GMRFLib_inla_mode == GMRFLib_MODE_COMPACT) {
			if (!free_linear_predictor) {
				linear_predictor = Calloc(Npred, double);
				free_linear_predictor = 1;
			}
			GMRFLib_preopt_predictor(linear_predictor, mode, preopt);

			if (0) {
				for (int i = 0; i < preopt->Npred; i++) {
					printf("i predictor %d %f\n", i, linear_predictor[i]);
				}
			}

		} else {
			linear_predictor = mode;
		}

		cc_is_negative = 0;
		Memset(bcoof, 0, Npred * sizeof(double));
		Memset(ccoof, 0, Npred * sizeof(double));

#define CODE_BLOCK							\
		for (int i_ = 0; i_ < d_idx->n; i_++) {			\
			int idx = d_idx->idx[i_];			\
			double ccmin = cmin;				\
			double step_len = DMAX(FLT_EPSILON, optpar->step_len); \
			if (fl && fl[idx]) {				\
				/* then do nothing fancy and no checks, cmin = NULL*/ \
				GMRFLib_2order_approx(thread_id, &(aa[idx]), &(bcoof[idx]), &(ccoof[idx]), NULL, d[idx], \
						      linear_predictor[idx], idx, mode, loglFunc, loglFunc_arg, \
						      &step_len, &(optpar->stencil), NULL); \
			} else {					\
				if (ISINF(ccmin) == 1) {		\
					/* Enter adaptive mode. Try with increasing step_len until ccoof >0 */ \
					/* If not successful then fall back to default step_len with cmin=0 */ \
					while(1) {			\
						GMRFLib_2order_approx(thread_id, &(aa[idx]), &(bcoof[idx]), &(ccoof[idx]), NULL, d[idx], \
								      linear_predictor[idx], idx, mode, loglFunc, loglFunc_arg, \
								      &step_len, &(optpar->stencil), NULL); \
						/* if ok, we are done */ \
						if (ccoof[idx] > 0.0) break; \
						/* otherwise, increase the step_len and retry */ \
						step_len *= 10.0;	\
						/* unless we have gone to far... */ \
						if (step_len > 1.0) {	\
							ccmin = DBL_EPSILON; \
							GMRFLib_2order_approx(thread_id, &(aa[idx]), &(bcoof[idx]), &(ccoof[idx]), NULL, d[idx], \
									      linear_predictor[idx], idx, mode, loglFunc, loglFunc_arg, \
									      &(optpar->step_len), &(optpar->stencil), &ccmin); \
							break;		\
						}			\
					}				\
				} else {				\
					GMRFLib_2order_approx(thread_id, &(aa[idx]), &(bcoof[idx]), &(ccoof[idx]), NULL, d[idx], \
							      linear_predictor[idx], idx, mode, loglFunc, loglFunc_arg, \
							      &(optpar->step_len), &(optpar->stencil), &cmin); \
				}					\
				cc_is_negative = (cc_is_negative || ccoof[idx] < 0.0); \
				if (ccoof[idx] == cmin && b_strategy == INLA_B_STRATEGY_SKIP) { \
					bcoof[idx] = 0.0;		\
				}					\
			}						\
			bb[idx] += bcoof[idx];				\
			cc[idx] += ccoof[idx];				\
		}

		if (GMRFLib_openmp->adaptive && omp_get_level() == 0) {
			// this is the exception of the rule, as we want to run this in parallel if we are in adaptive-mode and level=0.
			int nt = GMRFLib_PARDISO_MAX_NUM_THREADS();
			omp_set_num_threads(nt);
			RUN_CODE_BLOCK(nt, 0, 0);
			omp_set_num_threads(GMRFLib_openmp->max_threads_outer);
		} else {
			RUN_CODE_BLOCK(GMRFLib_openmp->max_threads_inner, 0, 0);
		}
#undef CODE_BLOCK

		double *bb_use = NULL, *cc_use = NULL;
		if (preopt) {
			GMRFLib_preopt_update(thread_id, preopt, bb, cc);
			bb_use = preopt->total_b[thread_id];
			GMRFLib_daddto(n, b, bb_use);
			cc_use = c;			       /* that what is there from before */
		} else {
			assert(Npred == n);

			GMRFLib_daddto(n, b, bb);
			GMRFLib_daddto(n, c, cc);
			bb_use = bb;
			cc_use = cc;
		}

		if (!cc_positive) {
			cc_factor = DMIN(1.0, cc_factor * cc_factor_mult);
		}

		/*
		 * I thought this was quicker without store, as there is just reuse and no copy... but not.  I free lproblem below and set it to NULL, so
		 * it will always be lproblem = NULL 
		 */

		if (!lproblem) {
			if (preopt) {
				assert((void *) preopt == (void *) Qfunc_arg);
			}
			// GMRFLib_printf_Qfunc2(stdout, graph, Qfunc, Qfunc_arg);

			int ret;
			ret = GMRFLib_init_problem_store(thread_id, &lproblem, x, bb_use, cc_use, mean, graph, Qfunc, Qfunc_arg, constr, store);
			if (ret != GMRFLib_SUCCESS) {
				catch_error = 1;
			}
		} else {
			GMRFLib_init_problem_store(thread_id, &lproblem, x, bb_use, cc_use, mean, graph, Qfunc, Qfunc_arg, constr, store);
		}

		if (catch_error) {
			lproblem = NULL;
			break;
		}

		int flag_cycle_behaviour = 0;
		double f;
		if (GMRFLib_gaussian_data) {
			f = 1.0;
		} else {
			f = DMIN(1.0, (iter + 1.0) * optpar->nr_step_factor);
		}

		double err = 0.0;
#pragma omp simd reduction(+: err)
		for (int i = 0; i < n; i++) {
			err += SQR((lproblem)->mean_constr[i] - mode[i]);
		}
		err = sqrt(err / n);

#pragma omp simd
		for (int i = 0; i < n; i++) {
			mode[i] += f * ((lproblem)->mean_constr[i] - mode[i]);
		}

		if (iter > 0) {
			if ((float) err == (float) err_previous) {
				// we're down to some rounding error and cannot get any further. this weird situation has happend. 
				// flag_cycle_behaviour = 1;

				// I'll disable this for now. It is better to just run out of iterations. [Apr'22]
				flag_cycle_behaviour = 0;
			}

			if (err > 4.0 * err_previous) {
				iter += itmax;		       /* so we can restart... */
			}
		}

		if (optpar && optpar->fp)
			fprintf(optpar->fp, "[%1d] iteration %d error %.4g\n", thread_id, iter, err);

		if (GMRFLib_gaussian_data) {
			/*
			 * I need to update 'aa' as this is not evaluated in the mode! The sum of the a's are/might-be used later
			 */

			if (GMRFLib_inla_mode == GMRFLib_MODE_COMPACT) {
				GMRFLib_preopt_predictor(linear_predictor, mode, preopt);
			} else {
				linear_predictor = mode;
			}
#define CODE_BLOCK							\
			for (int i = 0; i < d_idx->n; i++) {		\
				int idx = d_idx->idx[i];		\
				double local_bb, local_cc;		\
				GMRFLib_2order_approx(thread_id,	\
						      &(aa[idx]), &local_bb, &local_cc, NULL, d[idx], linear_predictor[idx], idx, mode, loglFunc, \
						      loglFunc_arg, &(optpar->step_len), &(optpar->stencil), NULL); \
			}
			RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 0, 0);
#undef CODE_BLOCK
		}

		if (ISNAN(err)) {
			break;
		}

		if (GMRFLib_gaussian_data || err < optpar->abserr_step || flag_cycle_behaviour) {
			break;
		}

		GMRFLib_free_problem(lproblem);
		lproblem = NULL;

		err_previous = err;
	}

	Free(ccoof);
	Free(bcoof);

	if (free_linear_predictor)
		Free(linear_predictor);

	if (iter < itmax) {
		*problem = lproblem;
	} else {
		*problem = NULL;
		GMRFLib_free_problem(lproblem);
	}

	if (0) {
		if (*problem) {
			FIXME("Leave with x=");
			for (int i = 0; i < graph->n; i++)
				printf(" %.4f", (*problem)->sub_mean_constr[i]);
			printf("\n");
		}
	}

	if (!*problem) {
		if (nested == 1) {
			GMRFLib_ASSERT(*problem, GMRFLib_EOPTNR);
			GMRFLib_LEAVE_ROUTINE;
			return GMRFLib_EOPTNR;
		} else if (nested == 2) {
			GMRFLib_LEAVE_ROUTINE;
			return GMRFLib_EOPTNR;
		} else {
			/*
			 * fail to converge. restart with a reduced step_factor. 
			 */
			Memcpy(mode, mode_initial, n * sizeof(double));	/* store the starting value */
			Free(mode_initial);
			FREE_ALL;
			GMRFLib_optimize_param_tp new_optpar;
			Memcpy(&new_optpar, optpar, sizeof(GMRFLib_optimize_param_tp));
			new_optpar.nr_step_factor *= 0.5;
			new_optpar.max_iter *= 2;
			if (new_optpar.fp) {
				fprintf(new_optpar.fp,
					"\n\n%s: Optimisation fail to converge.\n\t\t\tRetry with optpar->nr_step_factor = %g and add trust-region\n",
					__GMRFLib_FuncName, new_optpar.nr_step_factor);
			}
			if (new_optpar.nr_step_factor < 1e-3) {
				GMRFLib_LEAVE_ROUTINE;
				return GMRFLib_EOPTNR;
			} else {
				/*
				 * add trust region; try to find the smallest 'lambda' that work fine. well, approximatly only...
				 */
				int retval, ntimes = 1000, stop = 0;
				double lambda = 10000.0,       /* first value for lambda */
				    lambda_fac = 0.1,	       /* decrease it with this ammount for each iteration */
				    lambda_lim = 1e-6;	       /* value of lambda where we exit the loop */
				double *c_new = Calloc(graph->n, double);

				for (int kk = 0; kk < ntimes; kk++) {
					for (int i = 0; i < graph->n; i++) {
						c_new[i] = lambda * Qfunc(thread_id, i, i, NULL, Qfunc_arg) + lambda * c[i] + lambda;
						if (x && (ISNAN(x[i]) || ISINF(x[i]))) {
							x[i] = mode[i];
						}
						if (x && (ISNAN(x[i]) || ISINF(x[i]))) {
							x[i] = 0.0;
						}
					}
					retval = GMRFLib_init_GMRF_approximation_store__intern(thread_id, problem, x, b, c_new, mean, d, fl,
											       loglFunc, loglFunc_arg,
											       graph, Qfunc, Qfunc_arg, constr,
											       &new_optpar, blockupdate_par, store,
											       aa, bb, cc, cmin, b_strategy, 2, preopt, d_idx);
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
						Memcpy(x, (*problem)->mean_constr, graph->n * sizeof(double));
						GMRFLib_free_problem(*problem);
					} else {
						*problem = NULL;
						if (!stop) {
							GMRFLib_LEAVE_ROUTINE;
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

int GMRFLib_ai_INLA_experimental(GMRFLib_density_tp ***density,
				 GMRFLib_density_tp ***density_transform, GMRFLib_transform_array_func_tp **tfunc,
				 GMRFLib_density_tp ***density_hyper,
				 GMRFLib_gcpo_tp **gcpo, GMRFLib_gcpo_param_tp *gcpo_param,
				 GMRFLib_ai_cpo_tp **cpo, GMRFLib_ai_po_tp **po, GMRFLib_ai_dic_tp *dic,
				 GMRFLib_ai_marginal_likelihood_tp *marginal_likelihood,
				 double ***hyperparam, int nhyper,
				 GMRFLib_ai_log_extra_tp *log_extra, void *log_extra_arg,
				 double *x, double *b, double *c, double *mean,
				 GMRFLib_bfunc_tp **bfunc, GMRFLib_prior_mean_tp **prior_mean,
				 double *d, int *fl,
				 GMRFLib_logl_tp *loglFunc, void *loglFunc_arg,
				 GMRFLib_graph_tp *graph, GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_arg,
				 GMRFLib_constr_tp *constr, GMRFLib_ai_param_tp *ai_par, GMRFLib_ai_store_tp *ai_store,
				 int nlin, GMRFLib_lc_tp **Alin, GMRFLib_density_tp ***dlin, GMRFLib_ai_misc_output_tp *misc_output,
				 GMRFLib_preopt_tp *preopt)
{
#define SET_MODE							\
	GMRFLib_opt_get_latent(x_mode);					\
	if (theta_mode) {						\
		GMRFLib_opt_get_hyper(theta_mode);			\
		for(int j_=0; j_ < GMRFLib_MAX_THREADS(); j_++) {	\
			for(int i_ = 0; i_ < nhyper; i_++){		\
				hyperparam[i_][j_][0] = theta_mode[i_]; \
			}						\
		}							\
	}

	int *k_max = NULL, *k_min = NULL, *k_maxx = NULL, *k_minn = NULL, ierr, *iz = NULL, *izz = NULL, *len =
	    NULL, *iz_axes = NULL, free_ai_par = 0, config_count = 0, dens_count = 0, dens_max, hyper_len = 0;

	double *hessian = NULL, *theta = NULL, *theta_mode = NULL, *x_mode = NULL, log_dens_mode = 0, log_dens, *z = NULL, **izs =
	    NULL, *stdev_corr_pos = NULL, *stdev_corr_neg = NULL, f, w, w_origo, tref, tu, *weights = NULL, *adj_weights =
	    NULL, *hyper_z = NULL, *hyper_ldens = NULL, **userfunc_values = NULL, *inverse_hessian = NULL, *timer = NULL,
	    **cpo_theta = NULL, **po_theta = NULL, **po2_theta = NULL, **po3_theta = NULL, **pit_theta = NULL, ***deviance_theta =
	    NULL, **failure_theta = NULL;

	GMRFLib_gcpo_elm_tp ***gcpo_theta = NULL;
	gsl_matrix *H = NULL, *eigen_vectors = NULL;
	gsl_eigen_symmv_workspace *work = NULL;
	gsl_vector *eigen_values = NULL;
	gsl_vector *sqrt_eigen_values = NULL;
	GMRFLib_density_tp ***dens = NULL;
	GMRFLib_density_tp ***lpred = NULL;
	GMRFLib_density_tp ***dens_transform = NULL;
	GMRFLib_density_tp ***lin_dens = NULL;
	GMRFLib_ai_store_tp **ais = NULL;
	double **lin_cross = NULL;
	GMRFLib_gcpo_groups_tp *gcpo_groups = NULL;

	GMRFLib_idx_tp *d_idx = NULL;
	GMRFLib_idx_create_x(&d_idx, preopt->Npred);
	for (int i = 0; i < preopt->Npred; i++) {
		if (d[i]) {
			GMRFLib_idx_add(&d_idx, i);
		}
	}

	assert(GMRFLib_inla_mode == GMRFLib_MODE_COMPACT);
	assert(preopt);

	if (!ai_par) {
		GMRFLib_default_ai_param(&ai_par);
		free_ai_par = 1;
	}

	GMRFLib_ASSERT(ai_par && (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_AUTO ||
				  ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_EMPIRICAL_BAYES ||
				  ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_USER ||
				  ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_USER_STD ||
				  ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_USER_EXPERT ||
				  ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_GRID ||
				  ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_CCD), GMRFLib_EPARAMETER);

	if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_AUTO) {
		ai_par->int_strategy = GMRFLib_AI_INT_STRATEGY_GRID;
	}

	GMRFLib_ENTER_ROUTINE;

	int n_warnings = 0;
	if (misc_output) {
		timer = misc_output->wall_clock_time_used;
		misc_output->mode_status = 0;
		misc_output->warnings = Calloc(1, char *);
	} else {
		timer = NULL;
	}

	if (timer) {
		timer[0] = GMRFLib_timer();
	}

	/*
	 * this has to be true, I think... 
	 */
	if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_CCD) {
		ai_par->dz = 1.0;
	}

	nhyper = IMAX(0, nhyper);
	ais = Calloc(GMRFLib_MAX_THREADS(), GMRFLib_ai_store_tp *);

	// need to determine dens_max
	GMRFLib_design_tp *tdesign = NULL;

	if (ai_par->mode_fixed) {
		// easier to override the design here
		GMRFLib_design_eb(&tdesign, nhyper);
	} else if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_CCD && nhyper > 0) {
		GMRFLib_design_ccd(&tdesign, nhyper);
	} else if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_GRID && nhyper > 0) {
		GMRFLib_design_grid(&tdesign, nhyper);
	} else if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_EMPIRICAL_BAYES || nhyper == 0) {
		GMRFLib_design_eb(&tdesign, nhyper);
	} else {
		tdesign = ai_par->int_design;
	}
	dens_max = tdesign->nexperiments;
	weights = Calloc(dens_max, double);
	izs = Calloc(dens_max, double *);
	x_mode = Calloc(graph->n, double);

	if (gcpo) {
		(*gcpo) = Calloc(1, GMRFLib_gcpo_tp);
		(*gcpo)->n = preopt->Npred;
		(*gcpo)->value = Calloc(preopt->Npred, double);
		(*gcpo)->kld = Calloc(preopt->Npred, double);
		(*gcpo)->mean = Calloc(preopt->Npred, double);
		(*gcpo)->sd = Calloc(preopt->Npred, double);
		(*gcpo)->groups = NULL;			       /* to be completed later */
	}

	if (cpo) {
		(*cpo) = Calloc(1, GMRFLib_ai_cpo_tp);
		(*cpo)->n = preopt->Npred;
		(*cpo)->value = Calloc(preopt->Npred, double *);
		(*cpo)->pit_value = Calloc(preopt->Npred, double *);
		(*cpo)->failure = Calloc(preopt->Npred, double *);
	}
	if (po) {
		(*po) = Calloc(1, GMRFLib_ai_po_tp);
		(*po)->n = preopt->Npred;
		(*po)->value = Calloc(preopt->Npred, double *);
	}

	if (GMRFLib_ai_INLA_userfunc0) {
		userfunc_values = Calloc(dens_max, double *);
	}

	dens = Calloc(graph->n, GMRFLib_density_tp **);
	dens_transform = Calloc(graph->n, GMRFLib_density_tp **);
	for (int j = 0; j < graph->n; j++) {
		dens[j] = Calloc(dens_max, GMRFLib_density_tp *);	/* storage for the marginals */
		if (tfunc && tfunc[j]) {
			dens_transform[j] = Calloc(dens_max, GMRFLib_density_tp *);
		}
	}

	lpred = Calloc(preopt->mnpred, GMRFLib_density_tp **);
	for (int j = 0; j < preopt->mnpred; j++) {
		lpred[j] = Calloc(dens_max, GMRFLib_density_tp *);	/* storage for the marginals */
	}

	if (gcpo) {
		gcpo_theta = Calloc(dens_max, GMRFLib_gcpo_elm_tp **);
	}

	if (cpo) {
		cpo_theta = Calloc(preopt->Npred, double *);   /* cpo-value conditioned on theta */
		pit_theta = Calloc(preopt->Npred, double *);   /* pit-value conditioned on theta */
		failure_theta = Calloc(preopt->Npred, double *);	/* failure indicator on theta */

		for (int jj = 0; jj < d_idx->n; jj++) {
			int j = d_idx->idx[jj];
			if (fl[j]) {
				continue;
			}
			cpo_theta[j] = Calloc(dens_max, double);
			pit_theta[j] = Calloc(dens_max, double);
			failure_theta[j] = Calloc(dens_max, double);
		}
	}
	if (po) {
		po_theta = Calloc(preopt->Npred, double *);    /* po-value conditioned on theta */
		po2_theta = Calloc(preopt->Npred, double *);   /* po-value conditioned on theta */
		po3_theta = Calloc(preopt->Npred, double *);   /* po-value conditioned on theta */
		for (int jj = 0; jj < d_idx->n; jj++) {
			int j = d_idx->idx[jj];
			if (fl[j]) {
				continue;
			}
			po_theta[j] = Calloc(dens_max, double);
			po2_theta[j] = Calloc(dens_max, double);
			po3_theta[j] = Calloc(dens_max, double);
		}
	}
	if (dic) {
		deviance_theta = Calloc(preopt->Npred, double **);	/* mean of deviance conditioned on theta */
		for (int jj = 0; jj < d_idx->n; jj++) {
			int j = d_idx->idx[jj];
			if (fl[j]) {
				continue;
			}
			deviance_theta[j] = Calloc(dens_max, double *);
		}
	}

	if (timer) {
		timer[0] = GMRFLib_timer() - timer[0];	       /* preparation */
		timer[1] = GMRFLib_timer();
	}

	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_OPTIMIZE, (void *) &nhyper, NULL);
	GMRFLib_opt_setup(hyperparam, nhyper, log_extra, log_extra_arg, NULL, x, b, c, mean, bfunc, d, fl, loglFunc,
			  loglFunc_arg, graph, Qfunc, Qfunc_arg, constr, ai_par, ai_store, preopt, d_idx);

	if (nhyper > 0) {
		/*
		 * the first step is to locate the mode of \pi(\theta | y). here we use the opt-optimiser routine.  NOTE that this
		 * '_setup' ensure that ai_store is changed for each call to _opt_f. this is a bit dirty programming, but there is no
		 * good way to get around it for the moment.
		 */
		GMRFLib_opt_setup(hyperparam, nhyper, log_extra, log_extra_arg, NULL, x, b, c, mean, bfunc, d, fl, loglFunc,
				  loglFunc_arg, graph, Qfunc, Qfunc_arg, constr, ai_par, ai_store, preopt, d_idx);

		/*
		 * the optimizer runs most smoothly when #threads is about nhyper+1, which is the number of `natural' threads for
		 * computing the gradient.
		 */
		theta = Calloc(nhyper, double);		       /* theta is the hyperparameters */
		theta_mode = Calloc(nhyper, double);
		z = Calloc(nhyper, double);

		if (ai_par->mode_restart) {

			if (ai_par->fp_log) {
				fprintf(ai_par->fp_log, "Optimise using %s\n", GMRFLib_AI_OPTIMISER_NAME(ai_par->optimiser));
			}

			switch (ai_par->optimiser) {
			case GMRFLib_AI_OPTIMISER_GSL:
			case GMRFLib_AI_OPTIMISER_DEFAULT:
			{
				int retval;
				int fd_save = ai_par->gradient_forward_finite_difference;
				if (ai_par->optimise_smart) {
					// I'm not sure this is a good idea, as its f()'s are slower with more threads unless for huge models.
					// So I turn this off (and revert back to old behavious) until I know more.
					if (1 || GMRFLib_openmp->max_threads_outer < nhyper * 2) {
						ai_par->gradient_forward_finite_difference = GMRFLib_TRUE;
						if (ai_par->fp_log) {
							fprintf(ai_par->fp_log,
								"Smart optimise part I: estimate gradient using forward differences\n");
						}
					} else {
						if (ai_par->fp_log) {
							fprintf(ai_par->fp_log,
								"Smart optimise part I: estimate gradient using central differences\n");
						}
					}
				}

				double gsl_epsg = ai_par->gsl_epsg;
				double gsl_epsf = ai_par->gsl_epsf;
				double gsl_epsx = ai_par->gsl_epsx;

				ai_par->gsl_epsg *= 5.0;
				ai_par->gsl_epsf *= 5.0;
				ai_par->gsl_epsx *= 5.0;

				retval = GMRFLib_gsl_optimize(ai_par);
				if (retval != GMRFLib_SUCCESS) {
					retval = GMRFLib_gsl_optimize(ai_par);
				}

				ai_par->gsl_epsg = gsl_epsg;
				ai_par->gsl_epsf = gsl_epsf;
				ai_par->gsl_epsx = gsl_epsx;

				if (ai_par->optimise_smart) {
					ai_par->restart = IMAX(1, ai_par->restart);
					ai_par->gradient_forward_finite_difference = GMRFLib_FALSE;
					if (ai_par->fp_log) {
						fprintf(ai_par->fp_log, "Smart optimise part II: estimate gradient using central differences\n");
						fprintf(ai_par->fp_log, "Smart optimise part II: restart optimiser\n");
					}
				}

				if (ai_par->restart) {
					for (int k = 0; k < IMAX(0, ai_par->restart); k++) {
						retval = GMRFLib_gsl_optimize(ai_par);	/* restart */
						if (retval != GMRFLib_SUCCESS) {
							retval = GMRFLib_gsl_optimize(ai_par);
						}
					}
				}

				if (ai_par->parallel_linesearch) {
					if (ai_par->fp_log) {
						fprintf(ai_par->fp_log, "Smart optimise part II: turn off parallel linesearch\n");
					}
					GMRFLib_opt_turn_off_parallel_linesearch();
					GMRFLib_gsl_optimize(ai_par);	/* restart */
				}

				GMRFLib_gsl_get_results(theta_mode, &log_dens_mode);
				ai_par->gradient_forward_finite_difference = fd_save;
			}
				break;

			default:
				GMRFLib_ASSERT(0 == 1, GMRFLib_EPARAMETER);
				break;
			}

			if (ai_par->fp_log) {
				fprintf(ai_par->fp_log, "Optim: Number of function evaluations = %1d\n", GMRFLib_opt_get_f_count());
			}
		} else {
			/*
			 * use the initial values only 
			 */
			for (int i = 0; i < nhyper; i++) {
				theta_mode[i] = hyperparam[i][0][0];
			}
			if (ai_par->fp_log) {
				fprintf(ai_par->fp_log, "Using known modal configuration = [");
				for (int i = 0; i < nhyper; i++) {
					fprintf(ai_par->fp_log, " %6.3f", theta_mode[i]);
				}
				fprintf(ai_par->fp_log, " ]\n");
			}
			// this is not needed as we do that below. I am not quite sure if we need this in general, but...
			int thread_id = 0;
			assert(omp_get_thread_num() == 0);
			GMRFLib_opt_f(thread_id, theta_mode, &log_dens_mode, &ierr, NULL, NULL);
			log_dens_mode *= -1.0;
			if (ai_par->fp_log) {
				fprintf(ai_par->fp_log, "Compute mode: %10.3f\n", log_dens_mode);
			}
		}

		SET_MODE;
		GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_HESSIAN, (void *) &nhyper, NULL);

		if (ai_par->fp_log) {
			fprintf(ai_par->fp_log, "Compute the Hessian using %s differences and step_size[%g]. Matrix-type [%s]\n",
				(ai_par->hessian_forward_finite_difference ? "forward" : "central"),
				ai_par->hessian_finite_difference_step_len, (ai_par->hessian_force_diagonal ? "diagonal" : "dense"));
		}

		/*
		 * The parameters for the adaptive hessian estimation is set in ai_par (hence G.ai_par in opt-interface.c).
		 */
		double log_dens_mode_save = log_dens_mode;
		int stupid_mode_iter = 0, smart_success = 0;
		int fd_save = ai_par->hessian_forward_finite_difference;

		hessian = Calloc(ISQR(nhyper), double);

		if (ai_par->mode_fixed) {
			if (ai_par->fp_log) {
				fprintf(ai_par->fp_log, "fixed_mode=1 so, artificially, Hessian=diag(1)\n");
			}
			for (int i = 0; i < nhyper; i++) {
				hessian[i + nhyper * i] = 1.0;
			}
		} else {
			if (!(ai_par->optimise_smart) || !smart_success) {

				if (ai_par->optimise_smart) {
					ai_par->hessian_forward_finite_difference = GMRFLib_FALSE;
					if (ai_par->fp_log) {
						fprintf(ai_par->fp_log, "Smart optimise part IV: estimate Hessian using central differences\n");
					}
				}

				while (GMRFLib_opt_estimate_hessian(hessian, theta_mode, &log_dens_mode, stupid_mode_iter) != GMRFLib_SUCCESS) {
					if (!stupid_mode_iter) {
#pragma omp critical (Name_4edc3808412cc155a1bc15b9d4cd319a85e133b9)
						{
							if (ai_par->fp_log) {
								fprintf(ai_par->fp_log,
									"Mode not sufficient accurate; switch to a stupid local search strategy.\n");
							}
							n_warnings++;
							misc_output->warnings = Realloc(misc_output->warnings, n_warnings + 1, char *);
							misc_output->warnings[n_warnings - 1] =
							    Strdup
							    ("Stupid local search strategy used: This is usually a sign of a ill-defined model and/or non-informative data.");
							misc_output->warnings[n_warnings] = NULL;
						}
						stupid_mode_iter++;
					}

					if (log_dens_mode_save > log_dens_mode && stupid_mode_iter > ai_par->stupid_search_max_iter) {
						if (ai_par->fp_log) {
							fprintf(stderr,
								"\n\n*** Mode is not accurate yet but we have reached the rounding error level. Break.\n\n");
						}
						break;
					}
					// printf("%.12g %.12g\n", log_dens_mode_save, log_dens_mode);
					log_dens_mode_save = log_dens_mode;

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
			}
		}

		ai_par->hessian_forward_finite_difference = fd_save;

		/*
		 * do this again to get the ai_store set correctly.
		 */
		SET_MODE;
		if (stupid_mode_iter) {
			int thread_id = 0;
			assert(omp_get_thread_num() == 0);
			GMRFLib_opt_f(thread_id, theta_mode, &log_dens_mode, &ierr, NULL, NULL);
			log_dens_mode *= -1.0;
			SET_MODE;
		}

		if (ai_par->fp_log) {
			for (int i = 0; i < nhyper; i++) {
				for (int j = 0; j < nhyper; j++) {
					fprintf(ai_par->fp_log, " %10.3f", hessian[i + j * nhyper]);
				}
				fprintf(ai_par->fp_log, "\n");
			}
		}

		H = gsl_matrix_calloc((size_t) nhyper, (size_t) nhyper);
		for (int i = 0; i < nhyper; i++) {
			for (int j = 0; j < nhyper; j++) {
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
			GMRFLib_printf_gsl_matrix(ai_par->fp_log, eigen_vectors, "\t%7.3f");
			fprintf(ai_par->fp_log, "Eigenvalues of the Hessian\n");
			gsl_vector_fprintf(ai_par->fp_log, eigen_values, "\t%12.3f");
		}

		/*
		 * check that the hessian is positive definite 
		 */
		double min_pos_eigenvalue = DBL_MAX;
		for (int i = 0; i < nhyper; i++) {
			double eigv = gsl_vector_get(eigen_values, (unsigned int) i);
			if (eigv > 0.0) {
				min_pos_eigenvalue = DMIN(min_pos_eigenvalue, eigv);
			}
		}
		if (min_pos_eigenvalue == DBL_MAX) {
			min_pos_eigenvalue = 1.0;	       /* if all are negative, zero included */
		}
		int a_change = 0, all_negative = 1;

		for (int i = 0; i < nhyper; i++) {
			double eigv = gsl_vector_get(eigen_values, (unsigned int) i);

			all_negative = (all_negative && (eigv <= 0.0 || ISZERO(eigv)));
			if (eigv < 0.0) {
#pragma omp critical (Name_5de6e658f0ecd07258a12caddb9702c5fdf878dc)
				{
					fprintf(stderr, "\n");
					fprintf(stderr, "\t*** WARNING *** Eigenvalue %1d of the Hessian is %.6g < 0\n", i, eigv);
					fprintf(stderr, "\t*** WARNING *** This have consequence for the accurancy of the hyperpar\n");
					fprintf(stderr, "\t*** WARNING *** Continue with a diagonal Hessian.\n");
					fprintf(stderr, "\n");
					char *msg = NULL;
					GMRFLib_sprintf(&msg,
							"Hessian.eigen.value[%1d] = %.3f < 0. This is usually a sign of a ill-defined model and/or non-informative data.",
							i, eigv);
					n_warnings++;
					misc_output->warnings = Realloc(misc_output->warnings, n_warnings + 1, char *);
					misc_output->warnings[n_warnings - 1] = msg;
					misc_output->warnings[n_warnings] = NULL;

					gsl_vector_set(eigen_values, (unsigned int) i, min_pos_eigenvalue);
					a_change += 1000;
				}
			}
		}

		if (a_change) {
			if (misc_output) {
				misc_output->mode_status += a_change;	/* not a 'good mode'... */
			}
		}

		sqrt_eigen_values = gsl_vector_alloc((unsigned int) nhyper);
		for (int i = 0; i < nhyper; i++) {
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
					"\n\t*** WARNING *** R-inla: All eigenvalues of the Hessian are negative. Move on with Hessian = Identity\n\n");
				Memset(hessian, 0, ISQR(nhyper) * sizeof(double));
				for (int i = 0; i < nhyper; i++)
					hessian[i + i * nhyper] = 1.0;
			} else {
				// new. revert back to a diagonal hessian
				for (int i = 0; i < nhyper; i++) {
					hessian[i + i * nhyper] = DMAX(DBL_EPSILON, hessian[i + i * nhyper]);
					for (int j = i + 1; j < nhyper; j++) {
						hessian[i + j * nhyper] = hessian[j + i * nhyper] = 0.0;
					}
				}
				// need the new eigenvalues/vectors for futher calculations. its easy, we just compute
				// them again.
				for (int i = 0; i < nhyper; i++) {
					for (int j = 0; j < nhyper; j++) {
						gsl_matrix_set(H, (size_t) i, (size_t) j, hessian[i + nhyper * j]);
					}
				}
				work = gsl_eigen_symmv_alloc((size_t) nhyper);
				gsl_eigen_symmv(H, eigen_values, eigen_vectors, work);
				gsl_eigen_symmv_free(work);
			}
		}

		/*
		 * compute the inverse hessian, for scaling purposes 
		 */
		inverse_hessian = Calloc(ISQR(nhyper), double);
		Memcpy(inverse_hessian, hessian, ISQR(nhyper) * sizeof(double));
		GMRFLib_comp_posdef_inverse(inverse_hessian, nhyper);

		if (misc_output) {
			misc_output->nhyper = nhyper;
			misc_output->cov_m = Calloc(ISQR(nhyper), double);
			Memcpy(misc_output->cov_m, inverse_hessian, ISQR(nhyper) * sizeof(double));
			misc_output->log_posterior_mode = log_dens_mode;

			/*
			 * I need these as well, as the correction terms needs it (and we need also the sign of the eigenvectors...). 
			 */
			misc_output->eigenvalues = Calloc(nhyper, double);
			for (int i = 0; i < nhyper; i++) {
				misc_output->eigenvalues[i] = 1.0 / gsl_vector_get(eigen_values, i);	/* need the eigenvalues of the cov.mat not
													 * hessian */
			}
			GMRFLib_gsl_mat2plain(&(misc_output->eigenvectors), eigen_vectors);
		}

		if (ai_par->fp_log) {
			/*
			 * print the stdev/correlation matrix: stdevs on the diagonal and the correlations on the off-diagonal.
			 */
			double val;

			fprintf(ai_par->fp_log, "StDev/Correlation matrix (scaled inverse Hessian)\n");
			for (int ii = 0; ii < nhyper; ii++) {
				for (int jj = 0; jj < nhyper; jj++) {
					if (jj >= ii) {
						if (ii == jj) {
							val = sqrt(inverse_hessian[ii + jj * nhyper]);
						} else {
							val = inverse_hessian[ii + jj * nhyper] /
							    sqrt(inverse_hessian[ii + ii * nhyper] * inverse_hessian[jj + jj * nhyper]);
						}
						fprintf(ai_par->fp_log, " %7.3f", val);
					} else {
						fprintf(ai_par->fp_log, " %7s", "");
					}
				}
				fprintf(ai_par->fp_log, "\n");
			}
		}

		/*
		 * search the space. first, look at main directions and stop if the density differ more than dlog_dens from the value
		 * at the mode, log_dens_mode. outside the main directions, only add the point if the corresponding values at the main
		 * directions is in. 
		 *
		 * include a check of user-arguments. if a node a fixed, we cannot compute its marginal density...
		 */

		iz = Calloc(nhyper, int);
		Memset(iz, 0, nhyper * sizeof(int));

		hyper_len = dens_max;
		hyper_z = Calloc(hyper_len * nhyper, double);
		hyper_ldens = Calloc(hyper_len, double);

		GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_HESSIAN_SCALE, (void *) &nhyper, NULL);

		stdev_corr_pos = Calloc(nhyper, double);
		stdev_corr_neg = Calloc(nhyper, double);

		if (ai_par->mode_fixed) {
			if (ai_par->fp_log) {
				fprintf(ai_par->fp_log, "mode_fixed=1, so artificially, scaling of sd is set to 1.0\n");
			}
			for (int i = 0; i < nhyper; i++) {
				stdev_corr_neg[i] = 1.0;
				stdev_corr_pos[i] = 1.0;
			}
		} else {
#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_outer)
			for (int k = 0; k < nhyper; k++) {
				int thread_id = omp_get_thread_num();
				double step = M_SQRT2;
				double f0, *zz = NULL, *ttheta = NULL, llog_dens;
				GMRFLib_ai_store_tp *s = NULL;

				zz = Calloc(nhyper, double);
				ttheta = Calloc(nhyper, double);
				Memset(zz, 0, nhyper * sizeof(double));

				if (GMRFLib_OPENMP_IN_PARALLEL_ONEPLUS_THREAD()) {
					if (!ais[thread_id]) {
						ais[thread_id] = GMRFLib_duplicate_ai_store(ai_store, GMRFLib_TRUE, GMRFLib_TRUE, GMRFLib_FALSE);
					}
					s = ais[thread_id];
				} else {
					s = ai_store;	       /* the common one */
				}

				zz[k] = step;
				GMRFLib_ai_z2theta(ttheta, nhyper, theta_mode, zz, sqrt_eigen_values, eigen_vectors);
				GMRFLib_opt_f_intern(thread_id, ttheta, &llog_dens, &ierr, s, NULL, NULL);
				llog_dens *= -1.0;
				f0 = log_dens_mode - llog_dens;
				stdev_corr_pos[k] = (f0 > 0.0 ? sqrt(SQR(step) / (2.0 * f0)) : 1.0);

				zz[k] = -step;
				GMRFLib_ai_z2theta(ttheta, nhyper, theta_mode, zz, sqrt_eigen_values, eigen_vectors);
				GMRFLib_opt_f_intern(thread_id, ttheta, &llog_dens, &ierr, s, NULL, NULL);
				llog_dens *= -1.0;
				f0 = log_dens_mode - llog_dens;
				stdev_corr_neg[k] = (f0 > 0.0 ? sqrt(SQR(step) / (2.0 * f0)) : 1.0);

				double gmean = sqrt(stdev_corr_pos[k] * stdev_corr_neg[k]);
				double lim = 0.3;
				if ((stdev_corr_neg[k] < lim && stdev_corr_pos[k] < lim) ||
				    (stdev_corr_neg[k] > 1.0 / lim && stdev_corr_pos[k] > 1.0 / lim)) {
#pragma omp critical (Name_4267c78b945bb634ad14f44c9e0a55dc1213e726)
					{
						char *w1 = NULL;
						char *w2 = NULL;
						GMRFLib_sprintf(&w1,
								"Skewness correction for transf.hyperpar[%1d] is to high/low: gmean = %.2g, corr=(%.2f,%.2f).",
								k, gmean, stdev_corr_neg[k], stdev_corr_pos[k]);
						GMRFLib_sprintf(&w2, "%s.",
								(ai_par->hessian_correct_skewness_only ?
								 "This IS corrected for, but is usually a sign of a ill-defined model and/or issues with the fit"
								 :
								 "This IS NOT corrected for and is usually a sign of a ill-defined model and/or issues with the fit"));
						if (ai_par->fp_log) {
							fprintf(ai_par->fp_log, "\n*** Warning *** %s\n                %s\n\n", w1, w2);
						}
						if (misc_output) {
							// yes, we have warnings[n_warnings] be the NULL-ptr so we do not need
							// to pass 'n_warnings'
							n_warnings++;
							misc_output->warnings = Realloc(misc_output->warnings, n_warnings + 1, char *);
							char *w12 = NULL;
							GMRFLib_sprintf(&w12, "%s %s", w1, w2);
							misc_output->warnings[n_warnings - 1] = w12;
							misc_output->warnings[n_warnings] = NULL;
						}
						Free(w1);
						Free(w2);
					}
				}

				if (ai_par->hessian_correct_skewness_only) {
					stdev_corr_pos[k] /= gmean;
					stdev_corr_neg[k] /= gmean;
				}

				Free(zz);
				Free(ttheta);
			}
		}

		if (misc_output) {
			misc_output->stdev_corr_pos = Calloc(nhyper, double);
			Memcpy(misc_output->stdev_corr_pos, stdev_corr_pos, nhyper * sizeof(double));
			misc_output->stdev_corr_neg = Calloc(nhyper, double);
			Memcpy(misc_output->stdev_corr_neg, stdev_corr_neg, nhyper * sizeof(double));
		}

		SET_MODE;
	} else {
		// just fill with 1's
		if (misc_output) {
			// these are now computed, hence we use the Gaussian approximation
			misc_output->stdev_corr_pos = Calloc(nhyper, double);
			misc_output->stdev_corr_neg = Calloc(nhyper, double);
			stdev_corr_pos = Calloc(nhyper, double);
			stdev_corr_neg = Calloc(nhyper, double);
			for (int k = 0; k < nhyper; k++) {
				stdev_corr_pos[k] = misc_output->stdev_corr_pos[k] = stdev_corr_neg[k] = misc_output->stdev_corr_neg[k] = 1.0;
			}
		} else {
			stdev_corr_pos = Calloc(nhyper, double);
			stdev_corr_neg = Calloc(nhyper, double);
			for (int k = 0; k < nhyper; k++) {
				stdev_corr_pos[k] = stdev_corr_neg[k] = 1.0;
			}
		}

		int thread_id = 0;
		assert(omp_get_thread_num() == 0);
		GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_OPTIMIZE, (void *) &nhyper, NULL);
		GMRFLib_opt_f(thread_id, theta_mode, &log_dens_mode, &ierr, NULL, NULL);
		log_dens_mode *= -1.0;
		SET_MODE;
	}

	misc_output->opt_trace = (nhyper ? GMRFLib_opt_trace_get() : NULL);
	if (nlin > 0) {
		lin_dens = Calloc(dens_max, GMRFLib_density_tp **);
		if (misc_output && misc_output->compute_corr_lin) {
			lin_cross = Calloc(dens_max, double *);
		}
	} else {
		nlin = 0;
	}

	for (int k = 0; k < nhyper; k++) {
		if (ai_par->fp_log) {
			fprintf(ai_par->fp_log,
				"%s corrected stdev for theta[%1d]: negative %.3f  positive %.3f\n",
				"Compute", k, stdev_corr_neg[k], stdev_corr_pos[k]);
		}
	}

	if (misc_output) {
		for (int k = 0; k < nhyper; k++) {
			if (ISEQUAL(misc_output->stdev_corr_pos[k], 1.0)) {
				misc_output->mode_status++;
			}
			if (ISEQUAL(misc_output->stdev_corr_neg[k], 1.0)) {
				misc_output->mode_status++;
			}
		}
	}

	if (timer) {
		timer[1] = GMRFLib_timer() - timer[1];
		timer[2] = GMRFLib_timer();
	}

	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_INTEGRATE_HYPERPAR, NULL, NULL);

	// if fixed_mode=1, then we need to use EB
	if (ai_par->mode_fixed) {
		if (ai_par->int_strategy != GMRFLib_AI_INT_STRATEGY_EMPIRICAL_BAYES) {
			ai_par->int_strategy = GMRFLib_AI_INT_STRATEGY_EMPIRICAL_BAYES;
			if (ai_par->fp_log) {
				fprintf(ai_par->fp_log, "int.strategy is set to EB, since mode_fixed=1\n");
			}
		}
	}

	GMRFLib_design_tp *design = NULL;
	if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_CCD && nhyper > 0) {
		GMRFLib_design_ccd(&design, nhyper);
	} else if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_GRID && nhyper > 0) {
		GMRFLib_design_grid(&design, nhyper);
	} else if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_EMPIRICAL_BAYES || nhyper == 0) {
		// collect these two case into one
		GMRFLib_design_eb(&design, nhyper);
	} else {
		design = ai_par->int_design;
	}

	if (design->nexperiments > 1 && (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_GRID)) {
		f = DMAX(ai_par->f0, 1.0) * sqrt((double) nhyper);
		w = 1.0 / ((design->nexperiments - 1.0) * (1.0 + exp(-0.5 * SQR(f)) * (SQR(f) / nhyper - 1.0)));
		w_origo = 1.0 - (design->nexperiments - 1.0) * w;
	} else {
		f = w = w_origo = 1.0;
	}

	GMRFLib_ASSERT(dens_count == 0, GMRFLib_ESNH);

	if (nlin) {
		// we have to shift indices for the preopt format. Any indices belonging to Predictor and/or APredictor will then become
		// negative, and we do not support that. revert this change after the next parallel loop
		for (int i = 0; i < nlin; i++) {
			for (int j = 0; j < Alin[i]->n; j++) {
				Alin[i]->idx[j] -= preopt->mnpred;
				if (Alin[i]->idx[j] < 0) {
					char *s = Strdup("Using Predictor and/or Apredictor in lincomb in 'experimental-mode' is not supported.");
					GMRFLib_ERROR_MSG_NO_RETURN(GMRFLib_EPARAMETER, s);
					exit(1);
				}
			}
		}
	}

	if (gcpo) {
		// build the gcpo-groups, but change the openmp-strategy, temporary, as _gcpo_build() parallelise in one level only
		int thread_id = 0;
		assert(omp_get_thread_num() == 0);
		GMRFLib_openmp_place_tp place = GMRFLib_openmp->place;
		GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_GCPO_BUILD, NULL, NULL);
		gcpo_groups = GMRFLib_gcpo_build(thread_id, ai_store, preopt, gcpo_param, fl);
		GMRFLib_openmp_implement_strategy(place, NULL, NULL);
	}
	// if we have to many threads in outer we can move them to the inner level. Note that this will not increase the number of threads for
	// PARDISO:chol/Qinv/reorder, but will do for PARDISO:solve.
	GMRFLib_openmp_place_tp place_save = GMRFLib_OPENMP_PLACES_DEFAULT;

	int nt = IMAX(1, IMIN(design->nexperiments, GMRFLib_openmp->max_threads_outer));
	if (GMRFLib_openmp->max_threads_outer > design->nexperiments) {
		int outer = design->nexperiments;
		int inner = IMAX(1, GMRFLib_MAX_THREADS() / outer);
		place_save = GMRFLib_openmp->place;
		GMRFLib_openmp_implement_strategy_special(outer, inner);
	}
#pragma omp parallel for private(log_dens, dens_count, tref, tu, ierr) num_threads(nt)
	for (int k = 0; k < design->nexperiments; k++) {
		int thread_id = omp_get_thread_num();

		double *z_local = NULL, *theta_local = NULL, log_dens_orig;
		GMRFLib_ai_store_tp *ai_store_id = NULL;
		GMRFLib_tabulate_Qfunc_tp *tabQfunc = NULL;
		double *bnew = NULL;

		dens_count = k;

		if (GMRFLib_OPENMP_IN_PARALLEL_ONEPLUS_THREAD()) {
			if (!ais[thread_id]) {
				ais[thread_id] = GMRFLib_duplicate_ai_store(ai_store, GMRFLib_FALSE, GMRFLib_TRUE, GMRFLib_FALSE);
			}
			ai_store_id = ais[thread_id];
		} else {
			ai_store_id = ai_store;		       /* the common one */
		}
		assert(ai_store_id);

		z_local = Calloc(nhyper, double);
		theta_local = Calloc(nhyper, double);

		if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_CCD || ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_GRID) {
			for (int i = 0; i < nhyper; i++) {
				z_local[i] = f * design->experiment[k][i]
				    * (design->experiment[k][i] > 0.0 ? stdev_corr_pos[i] : stdev_corr_neg[i]);
			}
		} else if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_USER_STD) {
			for (int i = 0; i < nhyper; i++) {
				z_local[i] = design->experiment[k][i]
				    * (design->experiment[k][i] > 0.0 ? stdev_corr_pos[i] : stdev_corr_neg[i]);
			}
		} else if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_USER || ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_USER_EXPERT) {
			for (int i = 0; i < nhyper; i++) {
				z_local[i] = design->experiment[k][i];
			}
		} else {
			// nothing
		}

		if ((nhyper > 0 || GMRFLib_OPENMP_IN_PARALLEL_ONEPLUS_THREAD()) && !ai_par->mode_fixed) {
			if (design->std_scale) {
				// convert to theta_local
				GMRFLib_ai_z2theta(theta_local, nhyper, theta_mode, z_local, sqrt_eigen_values, eigen_vectors);
			} else {
				// theta_local is, by request, the same as z_local
				Memcpy(theta_local, z_local, nhyper * sizeof(double));
			}

			GMRFLib_opt_f_intern(thread_id, theta_local, &log_dens, &ierr, ai_store_id, &tabQfunc, &bnew);
			log_dens *= -1.0;
			log_dens_orig = log_dens;

			// make sure z_local's are aligned with theta_local's, for later usage.
			GMRFLib_ai_theta2z(z_local, nhyper, theta_mode, theta_local, sqrt_eigen_values, eigen_vectors);
		} else {
			log_dens = log_dens_orig = log_dens_mode;
		}

		/*
		 * correct the log_dens due to the integration weights which is special for the CCD
		 * integration and the deterministic integration points
		 * 
		 */
		if (ISNAN(design->int_weight[k]) || ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_CCD) {
			// integration weights are undefined. use these for the CCD design (as it _IS_
			// the CCD design in this case
			if (nhyper > 1) {
				/*
				 * the weight formula is only valid for nhyper > 1. 
				 */
				int origo = 1;

				for (int i = 0; i < nhyper; i++) {
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
			// integration weights are _given_. this is the deterministic integration points (which includes the _GRID)
			if (ai_par->int_strategy != GMRFLib_AI_INT_STRATEGY_USER_EXPERT) {
				log_dens += log(design->int_weight[k]);
			} else {
				// we do that later
			}
		}

		for (int i = 0; i < nhyper; i++) {
			hyper_z[dens_count * nhyper + i] = z_local[i];
		}
		if (nhyper) {
			hyper_ldens[dens_count] = log_dens_orig - log_dens_mode;
		}

		if (nhyper > 0) {
			if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_USER_EXPERT) {
				// In this case, the int_weights INCLUDE the log_dens
				weights[dens_count] = log(design->int_weight[k]);
			} else {
				weights[dens_count] = log_dens;
			}
		} else {
			weights[dens_count] = 1.0;
		}
		izs[dens_count] = Calloc(nhyper, double);

		for (int i = 0; i < nhyper; i++) {
			izs[dens_count][i] = z_local[i];
		}

		tref = GMRFLib_timer();
		GMRFLib_ai_add_Qinv_to_ai_store(ai_store_id);  /* add Qinv if its not there already */

#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_inner)
		for (int i = 0; i < graph->n; i++) {
			GMRFLib_density_create_normal(&dens[i][dens_count], 0.0, 1.0, ai_store_id->mode[i], ai_store_id->stdev[i], 0);
			if (tfunc && tfunc[i]) {
				GMRFLib_transform_density(&dens_transform[i][dens_count], dens[i][dens_count], tfunc[i]);
			}
		}

		if (ai_par->vb_enable && (ai_par->vb_strategy == GMRFLib_AI_VB_MEAN || ai_par->vb_strategy == GMRFLib_AI_VB_VARIANCE)) {
			GMRFLib_ai_vb_correct_mean_preopt(thread_id, dens, dens_count,
							  c, d, prior_mean, ai_par, ai_store_id, graph,
							  (tabQfunc ? tabQfunc->Qfunc : Qfunc), (tabQfunc ? tabQfunc->Qfunc_arg : Qfunc_arg),
							  loglFunc, loglFunc_arg, preopt, d_idx);
		}

		double *c_corrected = NULL;
		if (ai_par->vb_enable && (ai_par->vb_strategy == GMRFLib_AI_VB_VARIANCE)) {
			c_corrected = Calloc(graph->n, double);
			GMRFLib_ai_vb_correct_variance_preopt(thread_id, dens, dens_count,
							      c, d, ai_par, ai_store_id, graph,
							      (tabQfunc ? tabQfunc->Qfunc : Qfunc), (tabQfunc ? tabQfunc->Qfunc_arg : Qfunc_arg),
							      loglFunc, loglFunc_arg, preopt, c_corrected, d_idx);
		}

		double *mean_corrected = Calloc(graph->n, double);
		double *lpred_mean = Calloc(preopt->mnpred, double);
		double *lpred_mode = Calloc(preopt->mnpred, double);
		double *lpred_variance = Calloc(preopt->mnpred, double);

		for (int i = 0; i < graph->n; i++) {
			mean_corrected[i] = dens[i][dens_count]->user_mean;
		}
		GMRFLib_preopt_predictor_moments(lpred_mean, lpred_variance, preopt, ai_store_id->problem, mean_corrected);
		GMRFLib_preopt_predictor_moments(lpred_mode, NULL, preopt, ai_store_id->problem, NULL);

#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_inner)
		for (int i = 0; i < preopt->mnpred; i++) {
			GMRFLib_density_create_normal(&lpred[i][dens_count], 0.0, 1.0, lpred_mean[i], sqrt(lpred_variance[i]), 0);
		}

		double *gcpodens_moments = NULL;
		if (gcpo) {
			if (misc_output->configs_preopt) {
				gcpodens_moments = Calloc(3 * preopt->Npred, double);
				GMRFLib_fill(3 * preopt->Npred, NAN, gcpodens_moments);
			}
			gcpo_theta[dens_count] = GMRFLib_gcpo(thread_id, ai_store_id, lpred_mean, lpred_mode, lpred_variance, preopt, gcpo_groups,
							      d, loglFunc, loglFunc_arg, ai_par, gcpo_param, gcpodens_moments, d_idx);
		}

		if (GMRFLib_ai_INLA_userfunc0) {
			userfunc_values[dens_count] = GMRFLib_ai_INLA_userfunc0(thread_id, ai_store_id->problem, theta_local, nhyper);
		}

		if (nlin) {
			GMRFLib_ai_compute_lincomb(&(lin_dens[dens_count]), (lin_cross ? &(lin_cross[dens_count]) : NULL),
						   nlin, Alin, ai_store_id, mean_corrected, 0);
		}

		double *cpodens_moments = NULL;
		if (misc_output->configs_preopt && cpo) {
			cpodens_moments = Calloc(3 * preopt->Npred, double);
			GMRFLib_fill(3 * preopt->Npred, NAN, cpodens_moments);
		}

		if (cpo || dic || po) {
#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_inner)
			for (int ii = 0; ii < d_idx->n; ii++) {
				int i = d_idx->idx[ii];
				if (fl[i]) {
					continue;
				}
				GMRFLib_density_tp *cpodens = NULL;
				if (cpo) {
					GMRFLib_compute_cpodens(thread_id, &cpodens, lpred[i][dens_count], i, d[i], loglFunc, loglFunc_arg, ai_par);
					if (cpodens_moments) {
						if (cpodens) {
							cpodens_moments[3 * i + 0] = cpodens->user_mean;
							cpodens_moments[3 * i + 1] = SQR(cpodens->user_stdev);
							cpodens_moments[3 * i + 2] = cpodens->skewness;
						} else {
							cpodens_moments[3 * i + 0] = NAN;
							cpodens_moments[3 * i + 1] = NAN;
							cpodens_moments[3 * i + 2] = NAN;
						}
					}
					failure_theta[i][dens_count] = GMRFLib_ai_cpopit_integrate(thread_id,
												   &cpo_theta[i][dens_count],
												   &pit_theta[i][dens_count], i, cpodens,
												   d[i], loglFunc, loglFunc_arg, lpred_mean);
					if (cpodens && GMRFLib_getbit(cpodens->flags, DENSITY_FLAGS_FAILURE)) {
						failure_theta[i][dens_count] = 1.0;
					}
					GMRFLib_free_density(cpodens);
				}
				if (dic) {
					deviance_theta[i][dens_count] = GMRFLib_ai_dic_integrate(thread_id, i, lpred[i][dens_count], d[i],
												 loglFunc, loglFunc_arg, lpred_mean);
				}
				if (po) {
					GMRFLib_ai_po_integrate(thread_id, &po_theta[i][dens_count], &po2_theta[i][dens_count],
								&po3_theta[i][dens_count], i, lpred[i][dens_count], d[i], loglFunc, loglFunc_arg,
								lpred_mean);
				}
			}
		}

		char **arg_str = NULL;
		if (misc_output->likelihood_info) {
			assert(misc_output->configs_preopt);
			arg_str = Calloc(preopt->Npred, char *);
			for (int jj = 0; jj < d_idx->n; jj++) {
				int j = d_idx->idx[jj];
				double dummy;
				loglFunc(thread_id, &dummy, &(lpred_mean[j]), 1, j, NULL, NULL, loglFunc_arg, &(arg_str[j]));
			}
		}

		double *ll_info = NULL;
		if (misc_output->configs_preopt) {
			ll_info = Calloc(3 * preopt->Npred, double);
#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_inner)
			for (int j = 0; j < preopt->Npred; j++) {
				int jj = 3 * j;
				double local_aa;
				if (d[j]) {
					GMRFLib_2order_taylor(thread_id, &local_aa, &(ll_info[jj]), &(ll_info[jj + 1]), &(ll_info[jj + 2]), d[j],
							      lpred_mode[j], j, lpred_mode, loglFunc, loglFunc_arg, &ai_par->step_len,
							      &ai_par->stencil);
				} else {
					ll_info[jj] = ll_info[jj + 1] = ll_info[jj + 2] = NAN;
				}
			}
		}

		int free_if_not_configs = 1;
		if (misc_output->configs_preopt) {
			GMRFLib_ai_store_config_preopt(thread_id, misc_output, nhyper, theta_local, log_dens, log_dens_orig, ai_store_id->problem,
						       mean_corrected, preopt, Qfunc, Qfunc_arg, cpodens_moments, gcpodens_moments, arg_str,
						       ll_info, lpred_mean, lpred_variance, c_corrected);
			free_if_not_configs = 0;
		}

		tu = GMRFLib_timer() - tref;
		if (ai_par->fp_log) {
#pragma omp critical (Name_8a7254c4a570078955ae0e221dd0594e23386e57)
			{
				fprintf(ai_par->fp_log, "config %2d/%1d=[", config_count++, design->nexperiments);
				for (int i = 0; i < nhyper; i++) {
					fprintf(ai_par->fp_log, " %6.3f", z_local[i]);
				}
				/*
				 * we need to use the log_dens_orig as the other one is also included the integration weights. 
				 */
				fprintf(ai_par->fp_log, " ] log(rel.dens)= %6.3f, [%1d] accept, compute,",
					log_dens_orig - log_dens_mode, omp_get_thread_num());
				fprintf(ai_par->fp_log, " %.2fs\n", tu);
			}
		}

		GMRFLib_free_tabulate_Qfunc(tabQfunc);
		Free(bnew);
		Free(z_local);
		Free(theta_local);
		if (free_if_not_configs) {
			// if configs_preopt, then these vectors are store there hence only Free if we do not have configs=TRUE
			Free(lpred_mean);
			Free(lpred_variance);
			Free(cpodens_moments);
			Free(gcpodens_moments);
		}
		Free(mean_corrected);
		Free(c_corrected);
	}

	if (place_save) {
		GMRFLib_openmp_implement_strategy(place_save, NULL, NULL);
	}

	if (nlin) {
		// revert back as it was before
		for (int i = 0; i < nlin; i++) {
			for (int j = 0; j < Alin[i]->n; j++) {
				Alin[i]->idx[j] += preopt->mnpred;
			}
		}
	}
	// save (x, theta) adding the predictors
	preopt->mode_theta = Calloc(nhyper, double);
	Memcpy(preopt->mode_theta, theta_mode, nhyper * sizeof(double));
	preopt->mode_x = Calloc(preopt->mnpred + preopt->n, double);
	// GMRFLib_opt_get_latent(&(preopt->mode_x[preopt->mnpred]));
	Memcpy(&(preopt->mode_x[preopt->mnpred]), x_mode, preopt->n * sizeof(double));
	GMRFLib_preopt_full_predictor(preopt->mode_x, &(preopt->mode_x[preopt->mnpred]), preopt);

	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_DEFAULT, NULL, NULL);
	if (timer) {
		timer[2] = GMRFLib_timer() - timer[2];
		timer[3] = GMRFLib_timer();
	}

	*density = Calloc(preopt->n + preopt->mnpred, GMRFLib_density_tp *);
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
	GMRFLib_adjust_vector(weights, dens_max);
	for (int j = 0; j < dens_max; j++) {
		weights[j] = exp(weights[j]);
	}
	adj_weights = Calloc(dens_max, double);

	if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_CCD && nhyper == 1) {
		/*
		 * the CCD case, this is already done above, except for the case where nhyper=1.
		 */
		GMRFLib_ai_adjust_integration_weights(adj_weights, weights, izs, dens_max, nhyper, ai_par->dz);
	} else {
		Memcpy(adj_weights, weights, dens_max * sizeof(double));
	}

	if (ai_par->fp_log) {
		fprintf(ai_par->fp_log, "\nCombine the densities with relative weights:\n");
		for (int j = 0; j < dens_max; j++) {
			fprintf(ai_par->fp_log, "config %2d/%2d=[", j, dens_max);
			for (int k = 0; k < nhyper; k++) {
				fprintf(ai_par->fp_log, " %6.3f", izs[j][k]);
			}
			fprintf(ai_par->fp_log, " ] weight= %6.3f", weights[j]);
			if (ai_par->adjust_weights && (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_GRID)) {
				fprintf(ai_par->fp_log, " adjusted weight= %6.3f", adj_weights[j]);
			}
			fprintf(ai_par->fp_log, "\n");
		}
	}

	if (ai_par->fp_log) {
		GMRFLib_printMem(ai_par->fp_log);
	}

	/*
	 * to keep the same code as before 
	 */
	GMRFLib_idxval_tp *probs = GMRFLib_density_prune_weights(adj_weights, dens_max, GMRFLib_weight_prob_one);
	GMRFLib_idxval_tp *probs_combine = GMRFLib_density_prune_weights(adj_weights, dens_max, GMRFLib_weight_prob);

	// merge the two loops into one larger one for better omp
	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_COMBINE, NULL, NULL);

#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_outer)
	for (int ii = 0; ii < preopt->mnpred + graph->n; ii++) {
		int i;
		if (ii < preopt->mnpred) {
			i = ii;
			GMRFLib_density_tp *dens_combine = NULL;
			if (GMRFLib_save_memory) {
				// if skewness is to large then it will switch to the default...
				GMRFLib_density_combine_x(&dens_combine, lpred[i], probs_combine, GMRFLib_DENSITY_TYPE_SKEWNORMAL);
			} else {
				GMRFLib_density_combine(&dens_combine, lpred[i], probs_combine);
			}
			(*density)[i] = dens_combine;

			for (int k = 0; k < probs_combine->n; k++) {
				GMRFLib_free_density(lpred[i][k]);
				lpred[i][k] = NULL;
			}
			Free(lpred[i]);
		} else {
			i = ii - preopt->mnpred;
			GMRFLib_density_tp *dens_combine = NULL;
			if (GMRFLib_save_memory) {
				// if skewness is to large then it will switch to the default...
				GMRFLib_density_combine_x(&dens_combine, dens[i], probs, GMRFLib_DENSITY_TYPE_SKEWNORMAL);
			} else {
				GMRFLib_density_combine(&dens_combine, dens[i], probs);
			}
			(*density)[ii] = dens_combine;	       /* yes, its 'ii' */

			for (int k = 0; k < probs_combine->n; k++) {
				GMRFLib_free_density(dens[i][k]);
				dens[i][k] = NULL;
			}
			Free(dens[i]);

			if (tfunc && tfunc[i]) {
				GMRFLib_density_tp *dens_c = NULL;
				GMRFLib_density_combine(&dens_c, dens_transform[i], probs_combine);
				(*density_transform)[i] = dens_c;

				for (int k = 0; k < probs_combine->n; k++) {
					GMRFLib_free_density(dens_transform[i][k]);
					dens_transform[i][k] = NULL;
				}
				Free(dens_transform[i]);
			}
		}
	}

	if (ai_par->fp_log) {
		GMRFLib_printMem(ai_par->fp_log);
	}
	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_DEFAULT, NULL, NULL);

	if (dlin && nlin) {
		GMRFLib_density_tp **dtmp = NULL, *dcombine = NULL;

		assert(lin_dens);
		dtmp = Calloc(dens_max, GMRFLib_density_tp *);
		for (int j = 0; j < nlin; j++) {
			/*
			 * I need to do this as the storage is wrong.... 
			 */
			for (int k = 0; k < dens_max; k++) {
				dtmp[k] = lin_dens[k][j];
			}
			GMRFLib_density_combine(&dcombine, dtmp, probs_combine);
			(*dlin)[j] = dcombine;
		}
		Free(dtmp);

		if (misc_output && misc_output->compute_corr_lin) {
			double *ptmp = NULL;
			misc_output->corr_lin = ptmp = Calloc(ISQR(nlin), double);
			misc_output->cov_lin = Calloc(ISQR(nlin), double);

			for (int i = 0; i < nlin; i++) {
				for (int j = i; j < nlin; j++) {
					for (int kk = 0; kk < probs->n; kk++) {
						int k = probs->idx[kk];
						ptmp[i + j * nlin] += probs->val[kk] * lin_cross[k][i + j * nlin];
					}
					ptmp[j + i * nlin] = ptmp[i + j * nlin];
				}
			}
			Memcpy(misc_output->cov_lin, ptmp, ISQR(nlin) * sizeof(double));

			double *ptmp_scale = Calloc(nlin, double);
			for (int i = 0; i < nlin; i++) {
				ptmp_scale[i] = 1.0 / sqrt(ptmp[i + i * nlin]);
			}

			for (int i = 0; i < nlin; i++) {
				for (int j = i + 1; j < nlin; j++) {
					ptmp[i + j * nlin] = ptmp[i + j * nlin] * ptmp_scale[i] * ptmp_scale[j];
					ptmp[j + i * nlin] = ptmp[i + j * nlin];
				}
			}

			for (int i = 0; i < nlin; i++) {
				ptmp[i + i * nlin] = 1.0;
			}

			Free(ptmp_scale);
		}
	}
	if (ai_par->fp_log) {
		fprintf(ai_par->fp_log, "\n");
	}

	if (gcpo) {
		(*gcpo)->n = preopt->Npred;
		(*gcpo)->groups = gcpo_groups->groups;

		// if theta_correction is turned off, then all correction terms are 0
#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_outer)
		for (int j = 0; j < preopt->Npred; j++) {
			double lcorr_max = gcpo_theta[0][j]->marg_theta_correction;
			for (int jjj = 1; jjj < dens_max; jjj++) {
				lcorr_max = DMAX(lcorr_max, gcpo_theta[jjj][j]->marg_theta_correction);
			}
			for (int jjj = 0; jjj < dens_max; jjj++) {
				gcpo_theta[jjj][j]->marg_theta_correction -= lcorr_max;
				// prevent the corrections to be to large for robustness. exp(-15)=3.1E-07..
				// gcpo_theta[jjj][j]->marg_theta_correction = DMAX(-15.0, gcpo_theta[jjj][j]->marg_theta_correction);
				// P(exp(gcpo_theta[jjj][j]->marg_theta_correction));
			}

			double evalue = 0.0, evalue_one = 0.0;
			for (int jjj = 0; jjj < probs->n; jjj++) {
				int jj = probs->idx[jjj];
				evalue += gcpo_theta[jj][j]->value * probs->val[jjj] * exp(gcpo_theta[jj][j]->marg_theta_correction);
				evalue_one += probs->val[jjj] * exp(gcpo_theta[jj][j]->marg_theta_correction);
			}
			(*gcpo)->value[j] = evalue / evalue_one;

			evalue = 0.0;
			for (int jjj = 0; jjj < probs->n; jjj++) {
				int jj = probs->idx[jjj];
				evalue += gcpo_theta[jj][j]->kld * probs->val[jjj] * exp(gcpo_theta[jj][j]->marg_theta_correction);
			}
			(*gcpo)->kld[j] = evalue / evalue_one;

			evalue = 0.0;
			for (int jjj = 0; jjj < probs->n; jjj++) {
				int jj = probs->idx[jjj];
				evalue += gcpo_theta[jj][j]->lpred_mean * probs->val[jjj] * exp(gcpo_theta[jj][j]->marg_theta_correction);
			}
			(*gcpo)->mean[j] = evalue / evalue_one;

			evalue = 0.0;
			for (int jjj = 0; jjj < probs->n; jjj++) {
				int jj = probs->idx[jjj];
				evalue += (SQR(gcpo_theta[jj][j]->lpred_sd) + SQR(gcpo_theta[jj][j]->lpred_mean)) *
				    probs->val[jjj] * exp(gcpo_theta[jj][j]->marg_theta_correction);
			}
			(*gcpo)->sd[j] = sqrt(DMAX(0.0, evalue / evalue_one - SQR((*gcpo)->mean[j])));
		}
	}

	if (cpo) {
		double *Z = Calloc(preopt->Npred, double);

#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_outer)
		for (int j = 0; j < preopt->Npred; j++) {
			double evalue, evalue2, evalue_one = 1.0;
			int ii = j;

			if (cpo_theta[ii]) {
				for (int jjj = 0; jjj < probs->n; jjj++) {
					int jj = probs->idx[jjj];
					if (!ISNAN(cpo_theta[ii][jj]))	/* we ignore those that have failed */
						Z[ii] += probs->val[jjj] / cpo_theta[ii][jj];
				}
				(*cpo)->value[ii] = Calloc(1, double);
				if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_USER_EXPERT) {
					evalue = 0.0;
					for (int jjj = 0; jjj < probs->n; jjj++) {
						int jj = probs->idx[jjj];
						if (!ISNAN(cpo_theta[ii][jj])) {
							evalue += cpo_theta[ii][jj] * probs->val[jjj];
						}
					}
				} else {
					// here, we correct for adjusting pi(theta_j|y_{-i})
					evalue = evalue_one = 0.0;
					for (int jjj = 0; jjj < probs->n; jjj++) {
						int jj = probs->idx[jjj];
						if (!ISNAN(cpo_theta[ii][jj])) {
							evalue += cpo_theta[ii][jj] * probs->val[jjj] / cpo_theta[ii][jj] / Z[ii];
							evalue_one += probs->val[jjj] / cpo_theta[ii][jj] / Z[ii];
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
				evalue = evalue2 = evalue_one = 0.0;
				for (int jjj = 0; jjj < probs->n; jjj++) {
					int jj = probs->idx[jjj];
					if (!ISNAN(cpo_theta[ii][jj])) {
						evalue += pit_theta[ii][jj] * probs->val[jjj] / cpo_theta[ii][jj] / Z[ii];
						evalue_one += probs->val[jjj] / cpo_theta[ii][jj] / Z[ii];
					}
					/*
					 * this is defined over the unadjusted weights 
					 */
					evalue2 += failure_theta[ii][jj] * probs->val[jjj];
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

		(*cpo)->mean_value = (*cpo)->gmean_value = 0.0;
		if (preopt->Npred) {
			int count = 0;
			int gmean_inf = 0;

			for (int j = 0; j < preopt->Npred; j++) {
				int ii = j;

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
		SET_MODE;
#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_outer)
		for (int j = 0; j < preopt->Npred; j++) {
			double evalue, evalue2, evalue3, evalue_one;
			int ii = j;
			if (po_theta[ii]) {
				(*po)->value[ii] = Calloc(2, double);
				evalue_one = 1.0;
				evalue = evalue2 = evalue3 = 0.0;
				for (int jjj = 0; jjj < probs->n; jjj++) {
					int jj = probs->idx[jjj];
					if (po_theta[ii][jj]) {
						evalue += po_theta[ii][jj] * probs->val[jjj];
						evalue2 += po2_theta[ii][jj] * probs->val[jjj];
						evalue3 += po3_theta[ii][jj] * probs->val[jjj];
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
		SET_MODE;

		double mean_deviance = 0.0, mean_deviance_sat = 0.0, deviance_mean = 0.0, deviance_mean_sat = 0.0, *x_vec = NULL;
		int thread_id = 0;
		assert(omp_get_thread_num() == 0);

		/*
		 * need this for loglFunc() we need that compute is TRUE for all indices that enters loglFunc. There is no way to check this here. 
		 */
		x_vec = Calloc(preopt->mnpred, double);
		for (int j = 0; j < preopt->mnpred; j++) {
			if ((*density)[j]) {
				x_vec[j] = (*density)[j]->user_mean;
			} else {
				x_vec[j] = NAN;
			}
		}

		/*
		 * find the min length of the data contribution that cover all data points 
		 */
		int ndev = preopt->Npred;
		double *e_deviance = Calloc(ndev, double);
		double *e_deviance_sat = Calloc(ndev, double);
		double *deviance_e = Calloc(ndev, double);
		double *deviance_e_sat = Calloc(ndev, double);
		double *sign = Calloc(ndev, double);

		GMRFLib_fill(ndev, NAN, e_deviance);
		GMRFLib_fill(ndev, NAN, e_deviance_sat);
		GMRFLib_fill(ndev, NAN, deviance_e);
		GMRFLib_fill(ndev, NAN, deviance_e_sat);
		GMRFLib_fill(ndev, NAN, sign);

#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_outer) reduction(+ : deviance_mean,  deviance_mean_sat, mean_deviance, mean_deviance_sat)
		for (int j = 0; j < d_idx->n; j++) {
			double md = 0.0, md_sat = 0.0, dm = 0.0, dm_sat = 0.0, logl_sat = 0.0;
			int ii = d_idx->idx[j];

			if (fl[ii]) {
				continue;
			}

			double evalue = 0.0, evalue_sat = 0.0;
			for (int jjj = 0; jjj < probs->n; jjj++) {
				int jj = probs->idx[jjj];
				evalue += deviance_theta[ii][jj][0] * probs->val[jjj];
				evalue_sat += deviance_theta[ii][jj][1] * probs->val[jjj];
			}
			md = evalue;
			md_sat = evalue_sat;

			double x_tmp = (double) ((*density)[ii]->user_mean);
			double logl = 0.0;

			loglFunc(thread_id, &logl, &x_tmp, 1, ii, x_vec, NULL, loglFunc_arg, NULL);
			logl_sat = inla_compute_saturated_loglik(thread_id, ii, loglFunc, x_vec, loglFunc_arg);

			dm = -2.0 * d[ii] * logl;
			dm_sat = -2.0 * d[ii] * (logl - logl_sat);
			e_deviance[ii] = md;
			e_deviance_sat[ii] = md_sat;
			deviance_e[ii] = dm;
			deviance_e_sat[ii] = dm_sat;

			deviance_mean += dm;
			deviance_mean_sat += dm_sat;
			mean_deviance += md;
			mean_deviance_sat += md_sat;

			// neither of these options are fail-safe. I cannot see how to do this fail-safe without really mapping to the
			// real data doing the comparison there. But this information is not available at this level
			double sig = 0.0;
			if (loglFunc(0, NULL, NULL, 0, ii, NULL, NULL, loglFunc_arg, NULL) == GMRFLib_LOGL_COMPUTE_CDF) {
				loglFunc(0, &sig, &x_tmp, -1, ii, NULL, NULL, loglFunc_arg, NULL);
				sig = (sig <= 0.5 ? -1.0 : 1.0);
			} else {
				double xx[2], ld[2] = { 0.0, 0.0 };
				xx[0] = x_tmp;
				xx[1] = xx[0] + 0.01 * (*density)[ii]->user_stdev;
				loglFunc(0, ld, xx, 2, ii, NULL, NULL, loglFunc_arg, NULL);
				sig = (ld[1] > ld[0] ? 1.0 : -1.0);
			}
			sign[ii] = sig;
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
		dic->sign = sign;

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

		for (int j = 0; j < dim; j++) {
			double val = 0.0, wsum = 0.0, val2 = 0.0, mmean, vvar, ssd;

			for (int ii = 0; ii < probs->n; ii++) {
				int i = probs->idx[ii];
				wsum += probs->val[ii];
				val += probs->val[ii] * userfunc_values[i][j];
				val2 += probs->val[ii] * SQR(userfunc_values[i][j]);
			}

			mmean = val / wsum;
			vvar = DMAX(DBL_EPSILON, val2 / wsum - SQR(mmean));
			ssd = sqrt(vvar);
			GMRFLib_density_create_normal(&(GMRFLib_ai_INLA_userfunc0_density[j]), 0.0, 1.0, mmean, ssd, GMRFLib_TRUE);
		}
	}

	/*
	 * Compute the marginal likelihood; compute both the Gaussian approximatin and a non-parametric one. The marginal likelhood is the
	 * normalising constant for the posterior marginal for \theta. 
	 */
	if (marginal_likelihood) {
		if (nhyper > 0) {
			marginal_likelihood->marginal_likelihood_gaussian_approx = 0.5 * nhyper * log(2.0 * M_PI) + log_dens_mode;
			for (int i = 0; i < nhyper; i++) {
				marginal_likelihood->marginal_likelihood_gaussian_approx -=
				    0.5 * log(gsl_vector_get(eigen_values, (unsigned int) i));
			}

			if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_CCD) {
				/*
				 * in this case we integrate the 'ccd' approximation; the normal with stdev corrections. 
				 */
				marginal_likelihood->marginal_likelihood_integration = 0.5 * nhyper * log(2.0 * M_PI) + log_dens_mode;
				for (int i = 0; i < nhyper; i++) {
					marginal_likelihood->marginal_likelihood_integration -=
					    0.5 * (log(gsl_vector_get(eigen_values, (unsigned int) i)) +
						   0.5 * (log(SQR(stdev_corr_pos[i])) + log(SQR(stdev_corr_neg[i]))));
				}
			} else {
				double integral = 0.0, log_jacobian = 0.0;

				for (int j = 0; j < dens_max; j++) {
					integral += weights[j];
				}
				integral *= ai_par->dz;
				for (int i = 0; i < nhyper; i++) {
					log_jacobian -= 0.5 * log(gsl_vector_get(eigen_values, (unsigned int) i));
				}
				marginal_likelihood->marginal_likelihood_integration = log(integral) + log_jacobian + log_dens_mode;
			}
			if (ai_par->fp_log) {
				fprintf(ai_par->fp_log, "Marginal likelihood: Integration %.3f Gaussian-approx %.3f\n",
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
		if (ai_par->interpolator == GMRFLib_AI_INTERPOLATOR_GAUSSIAN) {
			/*
			 * Just use the modal values and the stdev's found from the Hessian. 
			 */
			for (int k = 0; k < nhyper; k++) {
				GMRFLib_density_create_normal(&((*density_hyper)[k]), 0.0, 1.0, theta_mode[k],
							      sqrt(inverse_hessian[k + nhyper * k]), GMRFLib_TRUE);
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
			double *std_stdev_theta = Calloc(nhyper, double);
			for (int k = 0; k < nhyper; k++) {
				std_stdev_theta[k] = sqrt(inverse_hessian[k + nhyper * k]);
			}

			/*
			 * write out the hole set 
			 */
			double *theta_tmp = NULL, log_jacobian = 0.0;

			theta_tmp = Calloc((int) IMAX(0, nhyper), double);
			if (eigen_values) {
				for (int k = 0; k < nhyper; k++) {
					log_jacobian -= 0.5 * log(gsl_vector_get(eigen_values, (unsigned int) k));
				}
			}

			if (nhyper > 0) {
				for (int kkk = 0; kkk < probs->n; kkk++) {
					int k = probs->idx[kkk];
					GMRFLib_ai_z2theta(theta_tmp, nhyper, theta_mode, &(hyper_z[k * nhyper]), sqrt_eigen_values, eigen_vectors);
					if (ai_par->fp_hyperparam) {
						for (int kk = 0; kk < nhyper; kk++) {
							fprintf(ai_par->fp_hyperparam, " %.10g", theta_tmp[kk]);
						}
						fprintf(ai_par->fp_hyperparam, " %.10g %.10g\n", hyper_ldens[k] + log_dens_mode + log_jacobian,
							probs->val[kkk]);
					}
				}
			}
			if (ai_par->fp_hyperparam) {
				fflush(ai_par->fp_hyperparam);
			}
			Free(theta_tmp);

			if (ai_par->fp_log) {
				fprintf(ai_par->fp_log,
					"\tCompute the marginal for theta[%1d] to theta[%1d] using numerical integration...\n", 0, nhyper - 1);
			}
			for (int k = 0; k < nhyper; k++) {
				GMRFLib_ai_marginal_one_hyperparamter(&((*density_hyper)[k]), k, nhyper, dens_max, hyper_z,
								      hyper_ldens, theta_mode, sqrt_eigen_values, eigen_vectors,
								      std_stdev_theta, ai_par->dz, stdev_corr_pos,
								      stdev_corr_neg, interpol, ai_par, inverse_hessian);
			}

			if (ai_par->fp_log) {
				fprintf(ai_par->fp_log,
					"\tCompute the marginal for theta[%1d] to theta[%1d] using numerical integration... Done.\n",
					0, nhyper - 1);
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
	SET_MODE;
	Memcpy(x, x_mode, graph->n * sizeof(double));

	if (misc_output) {
		/*
		 * store the reordering as well. 
		 */
		if (ai_store) {
			if (ai_store->problem->sub_sm_fact.remap != NULL) {
				misc_output->len_reordering = ai_store->problem->sub_graph->n;
				misc_output->nfunc = GMRFLib_opt_get_f_count();
				misc_output->opt_directions = GMRFLib_opt_get_directions();
				misc_output->reordering = Calloc(misc_output->len_reordering, int);
				Memcpy(misc_output->reordering, ai_store->problem->sub_sm_fact.remap, misc_output->len_reordering * sizeof(int));
			}
		}
	}

	if (GMRFLib_ai_INLA_userfunc1) {
		int thread_id = 0;
		assert(omp_get_thread_num() == 0);
		GMRFLib_ai_INLA_userfunc1(thread_id, theta_mode, nhyper, inverse_hessian);
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

		for (int i = 0; i < GMRFLib_ai_INLA_userfunc2_n; i++) {
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

		for (int i = 0; i < GMRFLib_ai_INLA_userfunc3_n; i++) {
			GMRFLib_ai_INLA_userfunc3[i] (i, theta_mode, nhyper, inverse_hessian, (void *) arg);
		}
		Free(arg);
	}

	if (ai_par->fp_log) {
		GMRFLib_printMem(ai_par->fp_log);
	}

	/*
	 * cleanup 
	 */
	if (izs) {
		for (int j = 0; j < dens_max; j++) {
			Free(izs[j]);
		}
		Free(izs);
	}
	if (lin_dens && nlin) {
		if (dens_max) {
			for (int j = 0; j < dens_max; j++) {
				for (int i = 0; i < nlin; i++)
					GMRFLib_free_density(lin_dens[j][i]);
				Free(lin_dens[j]);
			}
			Free(lin_dens);
		} else {
			for (int i = 0; i < nlin; i++)
				GMRFLib_free_density(lin_dens[0][i]);
			Free(lin_dens[0]);
			Free(lin_dens);
		}

		if (lin_cross) {
			for (int i = 0; i < dens_max; i++) {
				Free(lin_cross[i]);
			}
			Free(lin_cross);
		}
	}

	Free(stdev_corr_neg);
	Free(stdev_corr_pos);
	Free(hessian);
	Free(inverse_hessian);
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

	Free(adj_weights);
	Free(iz);
	Free(iz_axes);
	Free(izz);
	Free(k_max);
	Free(k_maxx);
	Free(k_min);
	Free(k_minn);
	Free(len);
	Free(theta);
	Free(theta_mode);
	Free(userfunc_values);
	Free(weights);
	Free(z);
	if (cpo_theta) {
		for (int i = 0; i < preopt->Npred; i++) {
			int j = i;
			if (d[j]) {
				Free(cpo_theta[j]);
			}
		}
		Free(cpo_theta);
	}
	if (po_theta) {
		for (int i = 0; i < preopt->Npred; i++) {
			int j = i;
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
		for (int i = 0; i < preopt->Npred; i++) {
			int j = i;
			if (d[j]) {
				Free(pit_theta[j]);
			}
		}
		Free(pit_theta);
	}
	if (failure_theta) {
		for (int i = 0; i < preopt->Npred; i++) {
			int j = i;
			if (d[j]) {
				Free(failure_theta[j]);
			}
		}
		Free(failure_theta);
	}
	if (deviance_theta) {
		for (int i = 0; i < preopt->Npred; i++) {
			int j = i;
			if (d[j]) {
				Free(deviance_theta[j][0]);
				Free(deviance_theta[j]);
			}
		}
		Free(deviance_theta);
	}
	if (free_ai_par) {
		Free(ai_par);
	}

	Free(hyper_z);
	Free(hyper_ldens);

	for (int i = 0; i < graph->n; i++) {
		if (dens[i]) {
			for (int j = 0; j < dens_max; j++) {
				GMRFLib_free_density(dens[i][j]);
			}
			Free(dens[i]);
		}
	}
	Free(dens);

	if (tfunc) {
		for (int i = 0; i < graph->n; i++) {
			if (tfunc[i]) {
				for (int j = 0; j < dens_max; j++) {
					GMRFLib_free_density(dens_transform[i][j]);
				}
				Free(dens_transform[i]);
			}
		}
	}
	Free(dens_transform);

	if (ais) {
		for (int k = 0; k < GMRFLib_MAX_THREADS(); k++) {
			if (ais[k]) {
				GMRFLib_free_ai_store(ais[k]);
			}
		}
		Free(ais);
	}

	if (nhyper) {
		GMRFLib_opt_exit();
	}

	if (timer) {
		timer[3] = GMRFLib_timer() - timer[3];
	}

	GMRFLib_idx_free(d_idx);

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
#undef SET_MODE

	return GMRFLib_SUCCESS;
}

int GMRFLib_equal_cor(double c1, double c2, GMRFLib_gcpo_param_tp *param)
{
#define COR2INTERN(c_) log((1.0 + (c_))/(1.0 - (c_)))

	// quick check
	if (ABS(c1 - c2) > param->sqrt_epsilon) {
		return 0;
	}

	double clim = 1.0 - FLT_EPSILON;
	double tc1 = TRUNCATE(c1, -clim, clim);
	double tc2 = TRUNCATE(c2, -clim, clim);

	if (ABS(COR2INTERN(tc1) - COR2INTERN(tc2)) < param->epsilon) {
		return 1;
	}

	return 0;
}

GMRFLib_gcpo_groups_tp *GMRFLib_gcpo_build(int thread_id, GMRFLib_ai_store_tp *ai_store, GMRFLib_preopt_tp *preopt,
					   GMRFLib_gcpo_param_tp *gcpo_param, int *UNUSED(fl))
{
#define A_idx(node_) (preopt->pAA_idxval ? preopt->pAA_idxval[node_] : preopt->A_idxval[node_])
#define W(node_) (gcpo_param->weights[node_])
#define LEGAL_TO_ADD(node_) (!(gcpo_param->group_selection) ? 1 :	\
			     GMRFLib_iwhich_sorted(node_,		\
						   gcpo_param->group_selection->idx, \
						   gcpo_param->group_selection->n) >= 0)

	GMRFLib_ENTER_ROUTINE;

	TIMER_INIT(0, 10);

	int detailed_output = GMRFLib_DEBUG_IF();
	int Npred = preopt->Npred;
	int mnpred = preopt->mnpred;
	int N = IMAX(preopt->n, Npred);
	GMRFLib_idxval_tp **groups = NULL;

	if (!(gcpo_param->weights) || (gcpo_param->weights && gcpo_param->len_weights < Npred)) {
		double *w = Calloc(Npred, double);
		GMRFLib_fill(Npred, 1.0, w);
		if (gcpo_param->weights) {
			// use those who already are defined
			Memcpy(w, gcpo_param->weights, gcpo_param->len_weights * sizeof(double));
			Free(gcpo_param->weights);
		}
		gcpo_param->weights = w;
		gcpo_param->len_weights = Npred;
	}

	TIMER_CHECK;

	if (!(gcpo_param->groups)) {
		if (gcpo_param->verbose || detailed_output) {
			printf("%s[%1d]: Build groups, strategy[%s]\n", __GMRFLib_FuncName, omp_get_thread_num(),
			       GMRFLib_GCPO_BUILD_STRATEGY_NAME(gcpo_param->build_strategy));
		}

		GMRFLib_ai_store_tp *build_ai_store = NULL;
		GMRFLib_ai_store_tp *local_ai_store = NULL;

		if (gcpo_param->build_strategy == GMRFLib_GCPO_BUILD_STRATEGY_PRIOR) {

			int nn = ai_store->problem->n;
			double *mask = Calloc(nn, double);
			double *diag = Calloc(nn, double);

			int n_offset = 0;
			int idx_offset = 0;
			if (!strcmp(gcpo_param->idx_tag[0], "APredictor")) {
				assert(!strcmp(gcpo_param->idx_tag[1], "Predictor"));
				n_offset = gcpo_param->idx_n[0] + gcpo_param->idx_n[1];
				idx_offset = 2;
			} else {
				assert(!strcmp(gcpo_param->idx_tag[0], "Predictor"));
				n_offset = gcpo_param->idx_n[0];
				idx_offset = 1;
			}

			double big = 1.0 / GSL_DBL_EPSILON;

			GMRFLib_fill(nn, 1.0, mask);
			GMRFLib_fill(nn, 0.0, diag);

			if (gcpo_param->remove) {
				int *visited = Calloc(gcpo_param->remove->n, int);
				for (int j = idx_offset; j < gcpo_param->idx_tot; j++) {
					char *tag = gcpo_param->idx_tag[j];
					int idx_found = -1;
					if (GMRFLib_str_is_member(gcpo_param->remove, tag, 0, &idx_found)) {
						for (int i = 0; i < gcpo_param->idx_n[j]; i++) {
							int k = gcpo_param->idx_start[j] - n_offset + i;
							mask[k] = 0.0;
							diag[k] = big;
						}
						visited[idx_found] = 1;
					}
				}

				int err = 0;
				for (int j = 0; j < gcpo_param->remove->n; j++) {
					if (!visited[j]) {
						err++;
						printf("\n[%1d] %s:%1d: *** error *** gcpo_param->remove[%1d]=[%s] is not found, abort!\n\n",
						       omp_get_thread_num(), __GMRFLib_FuncName, __LINE__, j, gcpo_param->remove->str[j]);
					}
				}
				assert(err == 0);
				Free(visited);
			}

			if (gcpo_param->keep) {
				int *visited = Calloc(gcpo_param->keep->n, int);
				// create a new default
#pragma omp simd
				for (int i = 0; i < nn; i++) {
					diag[i] = big;
				}
				Memset(mask, 0, nn * sizeof(double));

				for (int j = idx_offset; j < gcpo_param->idx_tot; j++) {
					char *tag = gcpo_param->idx_tag[j];
					int idx_found = -1;
					if (GMRFLib_str_is_member(gcpo_param->keep, tag, 0, &idx_found)) {
						for (int i = 0; i < gcpo_param->idx_n[j]; i++) {
							int k = gcpo_param->idx_start[j] - n_offset + i;
							mask[k] = 1.0;
							diag[k] = 0.0;
						}
						visited[idx_found] = 1;
					}
				}

				int err = 0;
				for (int j = 0; j < gcpo_param->keep->n; j++) {
					if (!visited[j]) {
						err++;
						printf("\n[%1d] %s:%1d: *** error *** gcpo_param->keep[%1d]=[%s] is not found, abort!\n\n",
						       omp_get_thread_num(), __GMRFLib_FuncName, __LINE__, j, gcpo_param->keep->str[j]);
					}
				}
				assert(err == 0);
				Free(visited);
			}

			if (gcpo_param->remove_fixed) {
				for (int j = idx_offset; j < gcpo_param->idx_tot; j++) {
					if (gcpo_param->idx_n[j] == 1) {
						int k = gcpo_param->idx_start[j] - n_offset;
						mask[k] = 0.0;
						diag[k] = big;
					}
				}
			}

			if (0) {
				for (int i = 0; i < gcpo_param->idx_tot; i++) {
					printf("\t\t%-30s %10d %10d\n", gcpo_param->idx_tag[i], gcpo_param->idx_start[i], gcpo_param->idx_n[i]);
				}
				for (int i = 0; i < nn; i++) {
					printf("i %d mask %g diag %g\n", i, mask[i], diag[i]);
				}
			}

			// pass args to gcpo_Qfunc
			preopt->gcpo_mask = mask;
			preopt->gcpo_diag = diag;

			GMRFLib_problem_tp *problem = NULL;
			double *c = Calloc(nn, double);
#pragma omp simd
			for (int i = 0; i < nn; i++) {
				c[i] = gcpo_param->prior_diagonal;
			}
			GMRFLib_init_problem(thread_id, &problem, NULL, NULL, c, NULL,
					     preopt->preopt_graph, preopt->gcpo_Qfunc, preopt->preopt_Qfunc_arg, preopt->latent_constr);
			local_ai_store = Calloc(1, GMRFLib_ai_store_tp);
			local_ai_store->problem = problem;
			preopt->gcpo_mask = NULL;	       /* set it back */
			preopt->gcpo_diag = NULL;	       /* set it back */

			Free(c);
			Free(mask);
			Free(diag);

			build_ai_store = local_ai_store;
		} else {
			assert(gcpo_param->build_strategy == GMRFLib_GCPO_BUILD_STRATEGY_POSTERIOR);
			build_ai_store = ai_store;
		}
		assert(build_ai_store != NULL);

		double *isd = Calloc(mnpred, double);
		groups = GMRFLib_idxval_ncreate_x(Npred, 1 + IABS((int) gcpo_param->num_level_sets));
		GMRFLib_ai_add_Qinv_to_ai_store(ai_store);
		GMRFLib_ai_add_Qinv_to_ai_store(build_ai_store);

		GMRFLib_preopt_predictor_moments(NULL, isd, preopt, build_ai_store->problem, NULL);
		double min_sd = sqrt(GMRFLib_min_value(isd, Npred, NULL));
#pragma omp simd
		for (int i = 0; i < Npred; i++) {
			isd[i] = 1.0 / sqrt(isd[i]);
		}

		GMRFLib_idx_tp *selection = NULL;
		if (!(gcpo_param->selection)) {
			GMRFLib_idx_create_x(&selection, Npred);
			for (int i = 0; i < Npred; i++) {
				GMRFLib_idx_add(&selection, i);
			}
		} else {
			GMRFLib_idx_create_x(&selection, gcpo_param->selection->n);
			for (int i = 0; i < gcpo_param->selection->n; i++) {
				int j = gcpo_param->selection->idx[i];
				GMRFLib_idx_add(&selection, j);
			}
			// selection = gcpo_param->selection;
			assert(GMRFLib_imax_value(selection->idx, selection->n, NULL) < Npred);
			assert(GMRFLib_imin_value(selection->idx, selection->n, NULL) >= 0);
		}
		if (gcpo_param->verbose || detailed_output) {
			assert(selection);
			printf("%s[%1d]: Use selection of %1d indices and num.level.sets %1g\n", __GMRFLib_FuncName,
			       omp_get_thread_num(), selection->n, gcpo_param->num_level_sets);
		}
		assert(selection);

#define CODE_BLOCK							\
		for(int ii = 0;	ii < selection->n; ii++) {		\
			CODE_BLOCK_ALL_WORK_ZERO();			\
			int node = selection->idx[ii];			\
									\
			if (gcpo_param->num_level_sets == -1.0) {	\
				GMRFLib_idxval_add(&(groups[node]), node, 1.0); \
				continue;				\
			}						\
									\
			GMRFLib_idxval_tp *v = A_idx(node);		\
			double *cor = CODE_BLOCK_WORK_PTR(0);		\
			double *a = CODE_BLOCK_WORK_PTR(1);		\
			double *Sa = CODE_BLOCK_WORK_PTR(2);		\
			double *cor_abs = CODE_BLOCK_WORK_PTR(3);	\
			size_t *largest = (size_t *) CODE_BLOCK_WORK_PTR(4); \
									\
			GMRFLib_unpack(v->n, v->val, a, v->idx);	\
			GMRFLib_Qsolve(Sa, a, build_ai_store->problem, -1); \
									\
			double eps = 1.0E-4 * min_sd / isd[node];	\
			_Pragma("omp simd")				\
				for (int iii = 0; iii < build_ai_store->problem->sub_graph->n; iii++) { \
					if (ABS(Sa[iii]) < eps) Sa[iii] = 0.0; \
				}					\
									\
			for (int nnode = 0; nnode < Npred; nnode++) {	\
				GMRFLib_idxval_tp *vv = A_idx(nnode);	\
				double sum = 0.0;			\
				GMRFLib_dot_product_INLINE(sum, vv, Sa); \
				sum *= isd[node] * isd[nnode];		\
				cor[nnode] = TRUNCATE(sum, -1.0, 1.0);	\
				cor_abs[nnode] = ABS(cor[nnode]);	\
			}						\
									\
			cor[node] = cor_abs[node] = 1.0;		\
			int levels_ok = 0;				\
			double levels_magnify = 1.0;			\
									\
			while (!levels_ok) {				\
				groups[node]->n = 0;			\
				int siz_g = IMIN(Npred, (int) (levels_magnify * (ABS(gcpo_param->num_level_sets) + 4L))); \
				levels_magnify *= 4.0;			\
				GMRFLib_DEBUG_idddd("node siz_g Npred num_level_sets levels_magnify", node, (double) siz_g, (double) Npred, \
						    gcpo_param->num_level_sets, levels_magnify); \
				gsl_sort_largest_index(largest, (size_t) siz_g, cor_abs, (size_t) 1, (size_t) Npred); \
									\
				double sumw = W(node);			\
				int i_prev = (int) largest[0];		\
				double cor_abs_prev = 1.0;		\
				GMRFLib_idxval_add(&(groups[node]), i_prev, cor_abs_prev); \
				for (int i = 1; i < siz_g && !levels_ok; i++) {	\
					int i_new = (int) largest[i];	\
					double cor_abs_new = cor_abs[i_new]; \
					if (LEGAL_TO_ADD(i_new)) {	\
						/* we have to go to one more before we stop as we need to add all equal ones first */ \
						if (!GMRFLib_equal_cor(cor_abs_new, cor_abs_prev, gcpo_param)) { \
							if ((sumw >= gcpo_param->num_level_sets)) { \
								/* then we will go over if adding, then skip */ \
								levels_ok = 1; \
							} else {	\
								sumw += W(i_new); \
								i_prev = i_new;	\
								cor_abs_prev = cor_abs_new; \
								GMRFLib_DEBUG_id("add new level  i_new cor_abs_new", i_new, cor_abs_new); \
								GMRFLib_idxval_add(&(groups[node]), i_new , cor[i_new]); \
							}		\
						} else {		\
							cor_abs[i_new] = cor_abs_prev; \
							cor[i_new] = DSIGN(cor[i_new]) * cor_abs_prev; \
							GMRFLib_idxval_add(&(groups[node]), i_new, cor[i_new]); \
							GMRFLib_DEBUG_id("add to old level  i_new cor_abs_prev", i_new, cor_abs_prev); \
							/* use the maximum weight when they are equal */ \
							if (W(i_new) > W(i_prev)) { \
								/* correct sumw, reset i_prev to point to the max weight one */ \
								sumw += W(i_new) - W(i_prev); \
								i_prev =  i_new; \
							}		\
						}			\
					}				\
					if (!levels_ok) {		\
						if ((sumw > gcpo_param->num_level_sets) || \
						    (gcpo_param->size_max > 0 && groups[node]->n >= gcpo_param->size_max)) { \
							levels_ok = 1;	\
						}			\
					}				\
					if (groups[node]->n >= Npred) levels_ok = 1; /* emergency option */ \
				}					\
				if (levels_ok) {			\
					if (gcpo_param->verbose || detailed_output) { \
						printf("%s[%1d]: for node=%1d : sumw %g, num.nodes %1d\n", \
						       __GMRFLib_FuncName, omp_get_thread_num(), node,  sumw, groups[node]->n); \
						printf("%s[%1d]: stop because there are no more levels or size.max is reached.\n", \
						       __GMRFLib_FuncName, omp_get_thread_num()); } \
				}					\
				GMRFLib_DEBUG_d("found group with sum of weights", sumw); \
				GMRFLib_DEBUG_i("levels_ok", levels_ok); \
			}						\
			if (gcpo_param->friends) {			\
				int group_n = groups[node]->n;		\
				for (int i = 0; i < group_n; i++) {	\
					int new_node = groups[node]->idx[i]; \
					if (new_node < gcpo_param->friends_n) { \
						for (int j = 0; j < gcpo_param->friends[new_node]->n; j++) { \
							int new_node2 = gcpo_param->friends[new_node]->idx[j]; \
							if (LEGAL(new_node2, Npred)) { \
								GMRFLib_idxval_add(&(groups[node]), new_node2, cor[new_node2]); \
							}		\
						}			\
					}				\
				}					\
			}						\
			/* no prepare or accumulate */			\
			GMRFLib_idxval_nsort_x(&(groups[node]), 1, 1, 0, 0); \
			/* this can happen: ensure node is part of its own group, as it might have been thrown out */ \
			/* due to max_size is reached or there are to many with |correlation| = 1 */ \
			if (GMRFLib_iwhich_sorted(node, groups[node]->idx, groups[node]->n) < 0) { \
				GMRFLib_idxval_add(&(groups[node]), node, 1.0);	\
				GMRFLib_idxval_nsort_x(&(groups[node]), 1, 1, 0, 0); \
			}						\
			if (0) P(node);					\
			if (0) GMRFLib_idxval_printf(stdout, groups[node], "after adding friends"); \
		}

		RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 5, N);
#undef CODE_BLOCK

		GMRFLib_idx_free(selection);
		Free(isd);
		GMRFLib_free_ai_store(local_ai_store);
	} else {
		if (gcpo_param->verbose || detailed_output) {
			printf("%s[%1d]: Use user-defined groups\n", __GMRFLib_FuncName, omp_get_thread_num());
		}

		// if the number of given groups are to short compared to Npred, then pad with empty groups. if its longer, then that is an error.
		if (gcpo_param->ngroups < Npred) {
			gcpo_param->groups = Realloc(gcpo_param->groups, Npred, GMRFLib_idxval_tp *);
			for (int i = gcpo_param->ngroups; i < Npred; i++) {
				GMRFLib_idxval_create_x(&(gcpo_param->groups[i]), 1 + IABS((int) gcpo_param->num_level_sets));
			}
		} else if (gcpo_param->ngroups > Npred) {
			assert(gcpo_param->ngroups > Npred);
		}
		groups = gcpo_param->groups;
	}


	TIMER_CHECK;

	// add first off-diagonals
	GMRFLib_idx2_tp **missing = GMRFLib_idx2_ncreate_x(Npred, 1 + IABS((int) gcpo_param->num_level_sets));
	for (int node = 0; node < Npred; node++) {
		if (groups[node]->n > 1) {
			for (int i = 0; i < groups[node]->n; i++) {
				int ii = groups[node]->idx[i];
				for (int j = i + 1; j < groups[node]->n; j++) {
					int jj = groups[node]->idx[j];
					GMRFLib_idx2_add(&(missing[IMIN(ii, jj)]), IMAX(ii, jj), node);
				}
			}
		}
	}

	// then diagonals
	for (int node = 0; node < Npred; node++) {
		if (groups[node]->n > 0) {
			GMRFLib_idx2_add(&(missing[node]), node, node);
		}
	}

	// build what to return
	GMRFLib_gcpo_groups_tp *ggroups = Calloc(1, GMRFLib_gcpo_groups_tp);
	ggroups->Npred = Npred;
	ggroups->groups = groups;
	ggroups->missing = missing;

	if (detailed_output) {
#pragma omp critical (Name_0c006e103a84c0a6e6169eed5e739b8065a95b95)
		{
			for (int node = 0; node < Npred; node++) {
				char *msg = NULL;
				GMRFLib_sprintf(&msg, "%s[%1d]: node %d", __GMRFLib_FuncName, omp_get_thread_num(), node);
				if (groups[node]->n > 0) {
					GMRFLib_idxval_printf(stdout, groups[node], msg);
				}
				if (missing[node]->n > 0) {
					GMRFLib_idx2_printf(stdout, missing[node], msg);
				}
				Free(msg);
			}
		}
	}

	TIMER_SUMMARY;

#undef A_idx
#undef W
#undef LEGAL_TO_ADD
	GMRFLib_LEAVE_ROUTINE;
	return ggroups;
}

GMRFLib_gcpo_elm_tp **GMRFLib_gcpo(int thread_id, GMRFLib_ai_store_tp *ai_store_id, double *lpred_mean, double *lpred_mode,
				   double *lpred_variance, GMRFLib_preopt_tp *preopt,
				   GMRFLib_gcpo_groups_tp *groups, double *d, GMRFLib_logl_tp *loglFunc, void *loglFunc_arg,
				   GMRFLib_ai_param_tp *ai_par, GMRFLib_gcpo_param_tp *gcpo_param, double *gcpodens_moments, GMRFLib_idx_tp *d_idx)
{
#define A_idx(node_) (preopt->pAA_idxval ? preopt->pAA_idxval[node_] : preopt->A_idxval[node_])

	GMRFLib_ENTER_ROUTINE;

	TIMER_INIT(0, 10);

	int detailed_output = GMRFLib_DEBUG_IF();
	int Npred = preopt->Npred;
	int mnpred = preopt->mnpred;
	int nn = preopt->n;
	int N = IMAX(nn, Npred);
	int max_ng = -1;
	int corr_hyper = gcpo_param->correct_hyperpar;
	const int np = GMRFLib_INT_GHQ_POINTS;
	double zero = 0.0;
	// double spd_eps = GSL_SQRT_DBL_EPSILON;
	double diag_eps = GSL_ROOT4_DBL_EPSILON;
	double diag_scale = 1.0 + diag_eps;

	if (gcpo_param->verbose || detailed_output) {
		printf("enter _gcpo with...\n");
		for (int i = 0; i < Npred; i++) {
			printf("\ti %d lpred_mean %f lpred_mode %f lpred_variance %f\n", i, lpred_mean[i], lpred_mode[i], lpred_variance[i]);
		}
	}

	TIMER_CHECK;

	Calloc_init(mnpred + 3 * Npred, 4);
	double *sd = Calloc_get(mnpred);
	GMRFLib_gcpo_elm_tp **gcpo = Calloc(Npred, GMRFLib_gcpo_elm_tp *);
	for (int i = 0; i < Npred; i++) {
		sd[i] = sqrt(lpred_variance[i]);
		gcpo[i] = Calloc(1, GMRFLib_gcpo_elm_tp);
		gcpo[i]->idxs = groups->groups[i];	       // just a copy!
		max_ng = IMAX(max_ng, gcpo[i]->idxs->n);
		if (gcpo[i]->idxs->n) {
			gcpo[i]->cov_mat = gsl_matrix_calloc((size_t) gcpo[i]->idxs->n, (size_t) gcpo[i]->idxs->n);
			gsl_matrix_set_all(gcpo[i]->cov_mat, NAN);
		} else {
			gcpo[i]->cov_mat = NULL;
		}

		// this matrix is alone, just set it here
		if (gcpo[i]->idxs->n == 1 && gcpo[i]->idxs->idx[0] == i) {
			gcpo[i]->node_min = i;
			gcpo[i]->node_max = i;
			gsl_matrix_set(gcpo[i]->cov_mat, 0, 0, lpred_variance[i]);
		}
	}

	GMRFLib_idx_tp *node_idx = NULL;
	unsigned char *skip = Calloc(Npred, unsigned char);
	assert(skip);

	TIMER_CHECK;

	for (int node = 0; node < Npred; node++) {
		// this case does not need to be computed
		if (groups->missing[node]->n == 1 && groups->missing[node]->idx[0][0] == node && groups->missing[node]->idx[1][0] == node) {
			if (gcpo_param->verbose || detailed_output) {
				printf("%s[%1d]: node %d is singleton, skip solve\n", __GMRFLib_FuncName, omp_get_thread_num(), node);
			}
			skip[node] = 1;
		}
		if (groups->missing[node]->n > 0) {
			GMRFLib_idx_add(&node_idx, node);
		}
	}

	TIMER_CHECK;

#define CODE_BLOCK							\
	for (int inode = 0; inode < node_idx->n; inode++) {		\
		int node = node_idx->idx[inode];			\
		double *a = CODE_BLOCK_WORK_PTR(0);			\
		double *Sa = CODE_BLOCK_WORK_PTR(1);			\
		CODE_BLOCK_ALL_WORK_ZERO();				\
									\
		if (gcpo_param->verbose || detailed_output) {		\
			if (skip[node]) {				\
				printf("%s[%1d]: Skip solve for node %d\n", __GMRFLib_FuncName, omp_get_thread_num(), node); \
			} else {					\
				printf("%s[%1d]: Solve for node %d\n", __GMRFLib_FuncName, omp_get_thread_num(), node); \
			}						\
		}							\
		gcpo[node]->node_min = gcpo[node]->idxs->idx[0];	\
		gcpo[node]->node_max = gcpo[node]->idxs->idx[IMAX(0, gcpo[node]->idxs->n - 1)]; \
		gcpo[node]->idx_node = GMRFLib_iwhich_sorted(node, (int *) (gcpo[node]->idxs->idx), gcpo[node]->idxs->n); \
									\
		if (gcpo[node]->idxs->n > 0) {				\
			if (gcpo[node]->idx_node < 0) {			\
				P(inode);				\
				P(node);				\
				P(gcpo[node]->idxs->n);			\
				P(gcpo[node]->idx_node);		\
				GMRFLib_idxval_printf(stdout, gcpo[node]->idxs, "gcpo[node]->idxs"); \
			}						\
			assert(gcpo[node]->idx_node >= 0);		\
		}							\
									\
		int need_Sa = 1;					\
		for(int k = 0; k < groups->missing[node]->n; k++) {	\
			int nnode = groups->missing[node]->idx[0][k];	\
			int cm_idx = groups->missing[node]->idx[1][k];	\
			gsl_matrix *mat = gcpo[cm_idx]->cov_mat;	\
			int ii = GMRFLib_iwhich_sorted(node, (int *) gcpo[cm_idx]->idxs->idx, gcpo[cm_idx]->idxs->n); \
			int jj = GMRFLib_iwhich_sorted(nnode, (int *) gcpo[cm_idx]->idxs->idx, gcpo[cm_idx]->idxs->n); \
			assert(ii >= 0 && jj >= 0);			\
			gsl_matrix_set(mat, ii, ii, lpred_variance[node]); \
			if (jj != ii) {					\
				if (need_Sa) {				\
					assert(!skip[node]);		\
					GMRFLib_idxval_tp *v = A_idx(node); \
					GMRFLib_unpack(v->n, v->val, a, v->idx); \
					GMRFLib_Qsolve(Sa, a, ai_store_id->problem, -1); \
					need_Sa = 0;			\
				}					\
				GMRFLib_idxval_tp *v = A_idx(nnode);	\
				double sum = 0.0;			\
				GMRFLib_dot_product_INLINE(sum, v, Sa); \
				double f = sd[node] * sd[nnode];	\
				sum /= f;				\
				double cov = TRUNCATE(sum, -1.0, 1.0) * f; \
				gsl_matrix_set(mat, jj, jj, lpred_variance[nnode]); \
				gsl_matrix_set(mat, ii, jj, cov);	\
				gsl_matrix_set(mat, jj, ii, cov);	\
			}						\
		}							\
	}

	if (node_idx) {
		RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 2, N);
	}
#undef CODE_BLOCK

	GMRFLib_idx_free(node_idx);
	Free(skip);

	TIMER_CHECK;

	if (detailed_output) {
#pragma omp critical (Name_0139eb204165e8e82ee3aaaaff59eab1d5b3cc14)
		{
			for (int node = 0; node < Npred; node++) {
				if (gcpo[node]->cov_mat && gcpo[node]->cov_mat->size1 > 0) {
					printf("\ncov_mat for node=%d size=%d idx_node=%d\n", node, (int) gcpo[node]->cov_mat->size1,
					       (int) gcpo[node]->idx_node);
					GMRFLib_idxval_printf(stdout, gcpo[node]->idxs, Strdup("nodes in this group"));
					GMRFLib_printf_gsl_matrix(stdout, gcpo[node]->cov_mat, " %.8f");
				}
			}
		}
	}

	int num_error = 0;
	for (int node = 0; node < Npred; node++) {
		if (gcpo[node]->cov_mat) {
			for (int i = 0; i < (int) gcpo[node]->cov_mat->size1; i++) {
				for (int j = i; j < (int) gcpo[node]->cov_mat->size2; j++) {
					if (ISNAN(gsl_matrix_get(gcpo[node]->cov_mat, i, j))) {
						num_error++;
						printf("%s[%1d]: ERROR: covmat for node %d, element i %d j %d, is NAN\n",
						       __GMRFLib_FuncName, omp_get_thread_num(), node, i, j);
					}
				}
			}
		}
	}
	assert(num_error == 0);

	TIMER_CHECK;

	// precompute all the logll terms to avoid to compute them to many times
	double *local_aa = Calloc_get(Npred);
	double *local_bb = Calloc_get(Npred);
	double *local_cc = Calloc_get(Npred);

#define CODE_BLOCK							\
	for(int i = 0; i < d_idx->n; i++) {				\
		int nnode = d_idx->idx[i];				\
		double xx = lpred_mode[nnode];				\
		GMRFLib_2order_approx(thread_id, &(local_aa[nnode]), &(local_bb[nnode]), &(local_cc[nnode]), NULL, d[nnode], xx, nnode, \
				      lpred_mode, loglFunc, loglFunc_arg, &ai_par->step_len, &ai_par->stencil, &zero); \
	}

	RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 0, 0);
#undef CODE_BLOCK

	GMRFLib_idx_tp *node_idx2 = NULL;
	for (int node = 0; node < Npred; node++) {
		if (gcpo[node]->cov_mat) {
			GMRFLib_idx_add(&node_idx2, node);
		} else {
			gcpo[node]->value = NAN;
		}
	}

	typedef struct {
		gsl_matrix *B;
		gsl_matrix *Bt;
		gsl_matrix *BtHB;
		gsl_matrix *H;
		gsl_matrix *QQ;
		gsl_matrix *SS;
		gsl_matrix *S;
		gsl_vector *Bzmean;
		gsl_vector *mean;
		gsl_vector *mean_old;
		gsl_vector *zb;
		gsl_vector *zmean;
		gsl_vector *ztmp;

		gsl_vector *zstar;
		gsl_vector *zstarp;
		gsl_vector *xstar;
		gsl_vector *b;
		gsl_vector *vec;
		gsl_vector *vec2;
		gsl_matrix *C;
		gsl_matrix *Qstar;

		GMRFLib_gsl_low_rank_store_tp *low_rank_store;
		GMRFLib_gsl_ensure_spd_store_tp *ensure_spd_store;
		GMRFLib_gsl_spd_solve_store_tp *spd_solve_store;
		GMRFLib_gsl_ldnorm_store_tp *log_dnorm_store;
	} local_storage_tp;

// need this for the CODE_BLOCK..._x
#define CODE_BLOCK_WORK_TP_FREE(p_)			\
	if (1) {					\
		gsl_matrix_free((p_)->B);		\
		gsl_matrix_free((p_)->Bt);		\
		gsl_matrix_free((p_)->BtHB);		\
		gsl_matrix_free((p_)->H);		\
		gsl_matrix_free((p_)->QQ);		\
		gsl_matrix_free((p_)->SS);		\
		gsl_matrix_free((p_)->S);		\
		gsl_vector_free((p_)->mean);		\
		gsl_vector_free((p_)->mean_old);	\
		gsl_vector_free((p_)->zb);		\
		gsl_vector_free((p_)->zmean);		\
		gsl_vector_free((p_)->ztmp);		\
		gsl_vector_free((p_)->zstar);		\
		gsl_vector_free((p_)->zstarp);		\
		gsl_vector_free((p_)->xstar);		\
		gsl_vector_free((p_)->b);		\
		gsl_vector_free((p_)->vec);		\
		gsl_vector_free((p_)->vec2);		\
		gsl_matrix_free((p_)->C);		\
		gsl_matrix_free((p_)->Qstar);		\
		GMRFLib_gsl_low_rank_store_free((p_)->low_rank_store);	\
		GMRFLib_gsl_ensure_spd_store_free((p_)->ensure_spd_store); \
		GMRFLib_gsl_spd_solve_store_free((p_)->spd_solve_store); \
		GMRFLib_gsl_ldnorm_store_free((p_)->log_dnorm_store); \
	}

#define CODE_BLOCK							\
	for(int inode = 0; inode < node_idx2->n; inode++) {		\
		int node = node_idx2->idx[inode];			\
		CODE_BLOCK_ALL_WORK_ZERO();				\
		int ng = gcpo[node]->idxs->n;				\
		local_storage_tp *lstore = CODE_BLOCK_WORK_TP_PTR();	\
		if (lstore->B == NULL) {				\
			lstore->B = gsl_matrix_calloc(max_ng, max_ng);	\
			lstore->Bt = gsl_matrix_calloc(max_ng, max_ng);	\
			lstore->BtHB = gsl_matrix_calloc(max_ng, max_ng); \
			lstore->Bzmean = gsl_vector_alloc(max_ng);	\
			lstore->H = gsl_matrix_alloc(max_ng, max_ng);	\
			lstore->QQ = gsl_matrix_alloc(max_ng, max_ng);	\
			lstore->SS = gsl_matrix_alloc(max_ng, max_ng);	\
			lstore->S = gsl_matrix_alloc(max_ng, max_ng);	\
			lstore->mean = gsl_vector_calloc(max_ng);	\
			lstore->mean_old = gsl_vector_calloc(max_ng);	\
			lstore->zb = gsl_vector_alloc(max_ng);		\
			lstore->zmean = gsl_vector_alloc(max_ng);	\
			lstore->ztmp = gsl_vector_alloc(max_ng);	\
			lstore->low_rank_store = GMRFLib_gsl_low_rank_store_alloc(max_ng); \
			lstore->ensure_spd_store = GMRFLib_gsl_ensure_spd_store_alloc(max_ng); \
									\
			lstore->zstar = gsl_vector_alloc(max_ng);	\
			lstore->zstarp = gsl_vector_alloc(max_ng);	\
			lstore->xstar = gsl_vector_alloc(max_ng);	\
			lstore->b = gsl_vector_alloc(max_ng);		\
			lstore->vec = gsl_vector_alloc(max_ng);		\
			lstore->vec2 = gsl_vector_alloc(max_ng);	\
			lstore->C = gsl_matrix_alloc(max_ng, max_ng);	\
			lstore->Qstar = gsl_matrix_alloc(max_ng, max_ng); \
			lstore->spd_solve_store = GMRFLib_gsl_spd_solve_store_alloc(max_ng); \
			lstore->log_dnorm_store = GMRFLib_gsl_ldnorm_store_alloc(max_ng); \
		}							\
		lstore->B->size1 = lstore->B->size2 = (size_t) ng;	\
		lstore->S->size1 = gcpo[node]->cov_mat->size1;		\
		lstore->S->size2 = gcpo[node]->cov_mat->size2;		\
		gsl_matrix_memcpy(lstore->S, gcpo[node]->cov_mat);	\
		lstore->Bt->size1 = lstore->Bt->size2 = (size_t) ng;	\
		lstore->BtHB->size1 = lstore->BtHB->size2 = (size_t) ng; \
		lstore->Bzmean->size = (size_t) ng;			\
		lstore->H->size1 = lstore->H->size2 = (size_t) ng;	\
		lstore->H->size1 = lstore->H->size2 = (size_t) ng;	\
		lstore->QQ->size1 = lstore->QQ->size2 = (size_t) ng;	\
		lstore->mean->size = (size_t) ng;			\
		lstore->mean_old->size = (size_t) ng;			\
		lstore->zb->size = (size_t) ng;				\
		lstore->zmean->size = (size_t) ng;			\
		lstore->ztmp->size = (size_t) ng;			\
									\
		int *idxs = gcpo[node]->idxs->idx;			\
		size_t idx_node = gcpo[node]->idx_node;			\
		double bb_idx_node = NAN, cc_idx_node = NAN;		\
		double *bb = CODE_BLOCK_WORK_PTR(0);			\
		double *cc = CODE_BLOCK_WORK_PTR(1);			\
									\
		/* use local names */					\
		gsl_matrix *B = lstore->B;				\
		gsl_matrix *Bt = lstore->Bt;				\
		gsl_matrix *BtHB = lstore->BtHB;			\
		gsl_matrix *H = lstore->H;				\
		gsl_matrix *S = lstore->S;				\
		gsl_matrix *QQ = lstore->QQ;				\
		gsl_matrix *SS = lstore->SS;				\
		gsl_vector *Bzmean = lstore->Bzmean;			\
		gsl_vector *mean = lstore->mean;			\
		gsl_vector *mean_old = lstore->mean_old;		\
		gsl_vector *zb = lstore->zb;				\
		gsl_vector *zmean = lstore->zmean;			\
		gsl_vector *ztmp = lstore->ztmp;			\
									\
		if (detailed_output) {					\
			printf("node %d, idx_node %zu,cov mat\n", node, idx_node); \
			GMRFLib_printf_gsl_matrix(stdout, S, " %.8f ");	\
		}							\
									\
		gcpo[node]->marg_theta_correction = 0.0;		\
		for(int i = 0; i < ng; i++) {				\
			int nnode = idxs[i];				\
			gsl_vector_set(mean_old, (size_t) i, lpred_mode[nnode]); \
			bb[i] += local_bb[nnode];			\
			cc[i] += local_cc[nnode];			\
			if (i == (int) idx_node) {			\
				bb_idx_node = local_bb[nnode];		\
				cc_idx_node = local_cc[nnode];		\
			}						\
		}							\
									\
		if (detailed_output) {					\
			printf("node %d, cov.mat and mean\n", node);	\
			GMRFLib_printf_gsl_matrix(stdout, S, " %.8f ");	\
			GMRFLib_printf_gsl_vector(stdout, mean_old, " %.8f "); \
		}							\
									\
		double low_rank_eps = GSL_ROOT3_DBL_EPSILON;		\
		/* x = B_{n x m} z + mean.x */				\
		GMRFLib_gsl_low_rank_x(S, low_rank_eps, B, lstore->low_rank_store); \
		GMRFLib_gsl_transpose_matrix_x(B, Bt);			\
									\
		size_t n = (size_t) ng;					\
		size_t m = B->size2;					\
									\
		gsl_matrix_set_zero(H);					\
		for(size_t i = 0; i < n; i++) {				\
			gsl_matrix_set(H, i, i, cc[i]);			\
		}							\
		GMRFLib_gsl_mv(H, mean_old, ztmp);			\
		for(size_t i = 0; i < n; i++) {				\
			gsl_vector_set(ztmp, i, gsl_vector_get(ztmp, i) - bb[i]); \
		}							\
									\
		zb->size = m;						\
		GMRFLib_gsl_mv(Bt, ztmp, zb);				\
		BtHB->size1 = BtHB->size2 = m;				\
		GMRFLib_gsl_mmm(Bt, H, B, BtHB);			\
		QQ->size1 = QQ->size2 = m;				\
		SS->size1 = SS->size2 = m;				\
		gsl_matrix_set_zero(QQ);				\
		for(size_t i = 0; i < m; i++) {				\
			gsl_matrix_set(QQ, i, i, DMAX(0, 1.0 - gsl_matrix_get(BtHB, i, i))); \
			if (0) {					\
				for(size_t j = 0; j < i; j++) {		\
					/* the val should be zero since H is diagonal  as Bt %*% B = I_m */ \
					double val = - gsl_matrix_get(BtHB, i, j); \
					gsl_matrix_set(QQ, i, j, val);	\
					gsl_matrix_set(QQ, j, i, val);	\
				}					\
			}						\
		}							\
		gsl_matrix_memcpy(SS, QQ);				\
		GMRFLib_gsl_ensure_spd_inverse_x(SS, low_rank_eps, NULL, lstore->ensure_spd_store); \
		zmean->size = m;					\
		GMRFLib_gsl_mv(SS, zb, zmean);				\
		if (detailed_output) {					\
			printf("node %d, cov.mat and mean after low-rank transformation\n", node); \
			GMRFLib_printf_gsl_matrix(stdout, SS, " %.8f "); \
			GMRFLib_printf_gsl_vector(stdout, zmean, " %.8f "); \
		}							\
		GMRFLib_gsl_mv(B, zmean, Bzmean);			\
		for(size_t i = 0; i < n; i++) {				\
			gsl_vector_set(mean, i, gsl_vector_get(mean_old, i) + gsl_vector_get(Bzmean, i)); \
		}							\
		GMRFLib_gsl_mmm(B, SS, Bt, S);				\
									\
		double mean_idx_node = gsl_vector_get(mean, idx_node) + (lpred_mean[node] - lpred_mode[node]); \
		gcpo[node]->lpred_mean = mean_idx_node;			\
		gcpo[node]->lpred_sd = sqrt(DMAX(DBL_EPSILON, gsl_matrix_get(S, idx_node, idx_node) / (1 + 0.0 * diag_scale))); \
		gcpo[node]->kld =  0.5 * (SQR(gcpo[node]->lpred_sd) / lpred_variance[node] - 1.0 + \
					  SQR(gcpo[node]->lpred_mean - lpred_mean[node]) / lpred_variance[node] + \
					  log(lpred_variance[node] / SQR(gcpo[node]->lpred_sd))); \
									\
		if (d[node] && (!(gcpo_param->type) || gcpo_param->type[node] == 0)) { \
			/* do the integral by approximating phi(x)*exp(ll(x)) with a normal, and */ \
			/* then do the GHQ with respect to that normal as the kernel, with the */ \
			/* ``errors'' as the function */		\
			double *weights = NULL, *xx = NULL;		\
			GMRFLib_ghq(&xx, &weights, np);			\
			double *xp = CODE_BLOCK_WORK_PTR(2);		\
			double *loglik = CODE_BLOCK_WORK_PTR(3);	\
									\
			double val = 0.0;				\
			double loc_prec, loc_mean, loc_sd;		\
			double ll_prec = cc_idx_node;			\
			double ll_mean = bb_idx_node / cc_idx_node;	\
			double lp_prec = 1.0 / SQR(gcpo[node]->lpred_sd); \
			double lp_mean = gcpo[node]->lpred_mean;	\
									\
			loc_prec = lp_prec + ll_prec;			\
			loc_sd = 1.0 / sqrt(loc_prec);			\
			loc_mean = (lp_prec * lp_mean + ll_prec * ll_mean) / loc_prec; \
			for (int i = 0; i < np; i++) {			\
				xp[i] = loc_mean + loc_sd * xx[i];	\
			}						\
			loglFunc(thread_id, loglik, xp, np, node, lpred_mean, NULL, loglFunc_arg, NULL); \
									\
			double d_tmp = d[node];				\
			for (int i = 0; i < np; i++) {			\
				val += exp(d_tmp * loglik[i]		\
					   - 0.5 * lp_prec * SQR(xp[i] - lp_mean) \
					   + 0.5 * SQR(xx[i]))		\
					* weights[i];			\
			}						\
			gcpo[node]->value = val * sqrt(lp_prec/loc_prec); \
			if (corr_hyper) {				\
				gcpo[node]->marg_theta_correction = - log(gcpo[node]->value); \
			}						\
		} else if (d[node] && gcpo_param->type && gcpo_param->type[node] != 0) { \
			/* x = B z + mean_old */			\
			/* zmean */					\
			/* SS = cov(z), QQ=prec(z) */			\
			/* the nodes in the group */			\
			/* for(int i = 0; i < ng; i++) nnode = idxs[i];	*/ \
			double eps = 1.0E-6;				\
			int max_iter = 100;				\
									\
			lstore->zstar->size = m;			\
			gsl_vector *zstar = lstore->zstar;		\
									\
			lstore->zstarp->size = m;			\
			gsl_vector *zstarp = lstore->zstarp;		\
									\
			lstore->xstar->size = ng;			\
			gsl_vector *xstar = lstore->xstar;		\
									\
			lstore->b->size = ng;				\
			gsl_vector *b = lstore->b;			\
									\
			lstore->vec->size = m;				\
			gsl_vector *vec = lstore->vec;			\
									\
			lstore->vec2->size = m;				\
			gsl_vector *vec2 = lstore->vec2;		\
									\
			lstore->C->size1 = lstore->C->size2 = ng;	\
			gsl_matrix *C = lstore->C;			\
			gsl_matrix_set_zero(C);				\
									\
			lstore->Qstar->size1 = lstore->Qstar->size2 = m; \
			gsl_matrix *Qstar = lstore->Qstar;		\
			double lla = 0.0;				\
									\
			for (int iter = 0; iter < max_iter; iter++) {	\
				lla = 0.0;				\
				if (iter == 0) {			\
					/* This is the initial value corresponding to full-data */ \
					gsl_vector_set_zero(zstar);	\
					gsl_vector_memcpy(xstar, mean_old); \
					for (int i = 0; i < ng; i++) {	\
						int nnode = idxs[i];	\
						double xx = gsl_vector_get(xstar, i); \
						gsl_vector_set(b, i, local_bb[nnode]); \
						gsl_matrix_set(C, i, i, local_cc[nnode]); \
						lla += local_aa[nnode] + xx * (local_bb[nnode] - 0.5 * xx * local_cc[nnode]); \
					}				\
				} else {				\
					gsl_vector_memcpy(zstarp, zstar); \
					GMRFLib_gsl_mv(B, zstar, xstar); \
					for (int i = 0; i < ng; i++) {	\
						gsl_vector_set(xstar, i, gsl_vector_get(xstar, i) + gsl_vector_get(mean_old, i)); \
					}				\
					for (int i = 0; i < ng; i++) {	\
						int nnode = idxs[i];	\
						double ll_a = 0.0, ll_b = 0.0, ll_c = 0.0; \
						double xx = gsl_vector_get(xstar, i); \
						GMRFLib_2order_approx(thread_id, &ll_a, &ll_b, &ll_c, NULL, d[nnode], xx, nnode, \
								      lpred_mode, loglFunc, loglFunc_arg, &ai_par->step_len, &ai_par->stencil, &zero); \
						gsl_vector_set(b, i, ll_b); \
						gsl_matrix_set(C, i, i, ll_c); \
						lla += ll_a + xx * (ll_b - 0.5 * xx * ll_c); \
					}				\
				}					\
				GMRFLib_gsl_mmm(Bt, C, B, Qstar);	\
				for (size_t i = 0; i < m; i++) {	\
					gsl_matrix_set(Qstar, i, i, gsl_matrix_get(QQ, i, i) + gsl_matrix_get(Qstar, i, i)); \
					for (size_t j=0; j < i ; j++) {	\
						double val = gsl_matrix_get(QQ, i, j) + gsl_matrix_get(Qstar, i, j); \
						gsl_matrix_set(Qstar, i, j, val); \
						gsl_matrix_set(Qstar, j, i, val); \
					}				\
				}					\
				for (int i = 0; i < ng; i++) {		\
					gsl_vector_set(b, i, gsl_vector_get(b, i) - gsl_matrix_get(C, i, i) * gsl_vector_get(mean_old, i)); \
				}					\
				GMRFLib_gsl_mv(Bt, b, vec);		\
				GMRFLib_gsl_mv(QQ, zmean, vec2);	\
				for (size_t i = 0; i < m; i++) {	\
					gsl_vector_set(vec, i, gsl_vector_get(vec, i) + gsl_vector_get(vec2, i)); \
				}					\
				GMRFLib_gsl_spd_solve_x(Qstar, vec, zstar, lstore->spd_solve_store); \
				double rms = GMRFLib_gsl_rms(zstar, zstarp); \
									\
				if (gcpo_param->verbose || detailed_output) { \
					printf("%s[%1d]: iter %1d for node %d rms %.8g\n", __GMRFLib_FuncName, \
					       omp_get_thread_num(), iter, node, GMRFLib_gsl_rms(zstar, zstarp)); \
				}					\
				if (detailed_output) {			\
					GMRFLib_printf_gsl_vector(stdout, zstar, "\t%.6g"); \
				}					\
				if (rms < eps) {			\
					break;				\
				}					\
			}						\
			/* need to divide the call to _dnorm to avoid potential issue with common store */ \
			lla += GMRFLib_gsl_ldnorm_x(zstar, zmean, QQ, NULL, 0, lstore->log_dnorm_store); \
			gcpo[node]->value = lla - GMRFLib_gsl_ldnorm_x(NULL, NULL, Qstar, NULL, 0, lstore->log_dnorm_store); \
			if (corr_hyper) {				\
				/* ->value is in log-scale at this point */ \
				gcpo[node]->marg_theta_correction = - gcpo[node]->value; \
			}						\
			gcpo[node]->value = exp(gcpo[node]->value);	\
		} else {						\
			gcpo[node]->value = NAN;			\
		}							\
									\
		if (gcpodens_moments) {					\
			gcpodens_moments[node * 3 + 0] = gcpo[node]->lpred_mean; \
			gcpodens_moments[node * 3 + 1] = SQR(gcpo[node]->lpred_sd); \
			gcpodens_moments[node * 3 + 2] = gcpo[node]->marg_theta_correction; \
		}							\
									\
		if (gcpo_param->verbose || detailed_output) {		\
			printf("%s[%1d]: node %d lpred_mean %f lpred_sd %f kld %f value %f\n", \
			       __GMRFLib_FuncName, omp_get_thread_num(), \
			       node, gcpo[node]->lpred_mean, gcpo[node]->lpred_sd, gcpo[node]->kld, gcpo[node]->value); \
		}							\
	}

	RUN_CODE_BLOCK_X(GMRFLib_MAX_THREADS(), 4, IMAX(np, max_ng), local_storage_tp);
#undef CODE_BLOCK

	GMRFLib_idx_free(node_idx2);
	Calloc_free();
#undef A_idx

	TIMER_SUMMARY;

	GMRFLib_LEAVE_ROUTINE;
	return gcpo;
}

int GMRFLib_compute_cpodens(int thread_id, GMRFLib_density_tp **cpo_density, GMRFLib_density_tp *density,
			    int idx, double d, GMRFLib_logl_tp *loglFunc, void *loglFunc_arg, GMRFLib_ai_param_tp *ai_par)
{
	if (!d) {
		*cpo_density = NULL;
		return GMRFLib_SUCCESS;
	}

	const int debug = 0;
	int itry, flag, np, np_orig = GMRFLib_INT_GHQ_POINTS + 4, npx = 8, itmp, np_new = np_orig + 2 * npx;
	double *xp = NULL, *xp_tmp = NULL, *ld = NULL, *logcor = NULL, *x_user = NULL;
	double cor_eps = (GSL_SQRT_DBL_EPSILON * GSL_ROOT4_DBL_EPSILON), cor_max, range;

	Calloc_init(4 * np_new, 4);
	ld = Calloc_get(np_new);
	logcor = Calloc_get(np_new);
	x_user = Calloc_get(np_new);
	xp = Calloc_get(np_new);

	for (itry = 0; itry < 2; itry++) {
		np = np_orig;
		flag = 0;
		GMRFLib_ghq_abscissas(&xp_tmp, np);
		range = xp_tmp[np - 1];
		Memcpy(xp + npx, xp_tmp, np * sizeof(double));
		for (itmp = 0; itmp < npx; itmp++) {
			xp[itmp] = xp[npx] - range * (npx - itmp) / (double) npx;
			xp[np + npx + itmp] = xp[npx + np - 1] + range * (itmp + 1.0) / (double) npx;
		}
		np = np_new;
		if (debug) {
#pragma omp critical (Name_6bbd084d591e27d9aacdb520c3e898945c242003)
			{
				for (itmp = 0; itmp < np; itmp++)
					printf("xp[%1d] = %.3f\n", itmp, xp[itmp]);
				GMRFLib_density_printf(stdout, density);
			}
		}
		GMRFLib_evaluate_nlogdensity(ld, xp, np, density);
		GMRFLib_density_std2user_n(x_user, xp, np, density);
		loglFunc(thread_id, logcor, x_user, np, idx, NULL, NULL, loglFunc_arg, NULL);
		GMRFLib_dscale(np, d, logcor);

		if (debug && np) {
#pragma omp critical (Name_45542d32821a8fbfd2cec71e8219d7eeb4b423f2)
			{
				for (int i = 0; i < np; i++)
					printf("CPO: %d BEFORE x_user %g xp %g ld %g logcor %g ld-logcor %g\n", idx,
					       x_user[i], xp[i], ld[i], logcor[i], ld[i] - logcor[i]);
			}
		}
		cor_max = exp(log(cor_eps) + GMRFLib_max_value(logcor, np, NULL));
#pragma omp simd
		for (int i = 0; i < np; i++) {
			ld[i] += logcor[i] - 2.0 * GMRFLib_log_apbex(cor_max, logcor[i]);
		}
		int npp = np;
		GMRFLib_ai_correct_cpodens(ld, xp, &np, ai_par);
		flag = (npp > np);
		if (debug && np) {
#pragma omp critical (Name_c6e59ebf504f17645e98f57731cc4de48bd2748a)
			{
				for (int i = 0; i < np; i++)
					printf("CPO AFTER: %d %g %g\n", idx, xp[i], ld[i]);
			}
		}
		if (np > 4) {
			GMRFLib_density_create(cpo_density, GMRFLib_DENSITY_TYPE_SCGAUSSIAN, np, xp, ld,
					       density->std_mean, density->std_stdev, GMRFLib_FALSE);
			if (flag && cpo_density)
				GMRFLib_setbit(&((*cpo_density)->flags), DENSITY_FLAGS_FAILURE);
		} else {
			*cpo_density = NULL;
		}
		if (*cpo_density || itry == 1)
			break;
	}

	Calloc_free();
	return GMRFLib_SUCCESS;
}

int GMRFLib_ai_vb_prepare_mean(int thread_id,
			       GMRFLib_vb_coofs_tp *coofs, int idx, double d, GMRFLib_logl_tp *loglFunc,
			       void *loglFunc_arg, double *x_vec, double mean, double sd, double *workspace)
{
	// compute the Taylor-expansion of integral of -loglikelihood * density(x), around the mean of x.
	// optional workspace: size >= 2 * GMRFLib_INT_GHQ_ALLOC_LEN 

	// Normal kernel: deriv: ... * (x-m)/s^2
	// dderiv: ... * ((x-m)^2 - s^2)/s^4
	// GMRFLib_density_type_tp type;

	if (ISZERO(d)) {
		coofs->coofs[0] = coofs->coofs[1] = coofs->coofs[2] = 0.0;
		return GMRFLib_SUCCESS;
	}

	static double *xp = NULL;
	static double *wp = NULL;
	static double *wxp = NULL;
	static double *wxp2 = NULL;

	if (!wp) {
#pragma omp critical (Name_00c5c0bab9ee4213c2351e3b2275ded2f8b87d22)
		{
			if (!wp) {
				double *wtmp = NULL;
				GMRFLib_ghq(&xp, &wtmp, GMRFLib_INT_GHQ_POINTS);	/* just give ptr to storage */
				wxp = Calloc(GMRFLib_INT_GHQ_POINTS * 2, double);
				wxp2 = wxp + GMRFLib_INT_GHQ_POINTS;
				for (int i = 0; i < GMRFLib_INT_GHQ_POINTS; i++) {
					wxp[i] = wtmp[i] * xp[i];
					wxp2[i] = wtmp[i] * (SQR(xp[i]) - 1.0);
				}
				wp = wtmp;
			}
		}
	}

	double *x_user = NULL, *loglik = NULL;
	if (workspace) {
		x_user = workspace;
	} else {
		x_user = Calloc(2 * GMRFLib_INT_GHQ_ALLOC_LEN, double);
	}
	loglik = x_user + GMRFLib_INT_GHQ_ALLOC_LEN;
	GMRFLib_daxpb(GMRFLib_INT_GHQ_POINTS, sd, xp, mean, x_user);
	loglFunc(thread_id, loglik, x_user, GMRFLib_INT_GHQ_POINTS, idx, x_vec, NULL, loglFunc_arg, NULL);

	double s_inv = 1.0 / sd, s2_inv = SQR(s_inv);
	double B = GMRFLib_ddot(GMRFLib_INT_GHQ_POINTS, loglik, wxp);
	double C = GMRFLib_ddot(GMRFLib_INT_GHQ_POINTS, loglik, wxp2);

	// coofs->coofs[0] = -d * A;
	coofs->coofs[0] = NAN;
	coofs->coofs[1] = -d * B * s_inv;
	coofs->coofs[2] = -d * C * s2_inv;

	if (!workspace) {
		Free(x_user);
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_ai_vb_prepare_variance(int thread_id, GMRFLib_vb_coofs_tp *coofs, int idx, double d,
				   GMRFLib_logl_tp *loglFunc, void *loglFunc_arg, double *x_vec, double mean, double sd, double *workspace)
{
	// compute the Taylor-expansion of the integral of -loglikelihood * density(x), for the variance of x (assuming its N())
	// optional workspace: size >= 2 * GMRFLib_INT_GHQ_ALLOC_LEN 

	if (ISZERO(d)) {
		coofs->coofs[0] = coofs->coofs[1] = coofs->coofs[2] = 0.0;
		return GMRFLib_SUCCESS;
	}

	static double *wp = NULL;
	static double *xp = NULL;
	static double *wxp2 = NULL;
	static double *wxp3 = NULL;

	if (!wp) {
#pragma omp critical (Name_0713ff01bf46f0328663d7242f8e788872085a66)
		{
			if (!wp) {
				double *wtmp = NULL;
				GMRFLib_ghq(&xp, &wtmp, GMRFLib_INT_GHQ_POINTS);	/* just give ptr to storage */
				wxp2 = Calloc(2 * GMRFLib_INT_GHQ_ALLOC_LEN, double);
				wxp3 = wxp2 + GMRFLib_INT_GHQ_ALLOC_LEN;
				for (int i = 0; i < GMRFLib_INT_GHQ_POINTS; i++) {
					double z2 = SQR(xp[i]);
					wxp2[i] = wtmp[i] * (z2 - 1.0);
					wxp3[i] = wtmp[i] * (3.0 - 6.0 * z2 + SQR(z2));
				}
				wp = wtmp;
			}
		}
	}

	double *x_user = NULL, *loglik = NULL;
	double s2_inv = 1.0 / SQR(sd);

	if (workspace) {
		x_user = workspace;
	} else {
		x_user = Calloc(2 * GMRFLib_INT_GHQ_ALLOC_LEN, double);
	}
	loglik = x_user + GMRFLib_INT_GHQ_ALLOC_LEN;

	GMRFLib_daxpb(GMRFLib_INT_GHQ_POINTS, sd, xp, mean, x_user);
	loglFunc(thread_id, loglik, x_user, GMRFLib_INT_GHQ_POINTS, idx, x_vec, NULL, loglFunc_arg, NULL);

	double B = GMRFLib_ddot(GMRFLib_INT_GHQ_POINTS, wxp2, loglik);
	double C = GMRFLib_ddot(GMRFLib_INT_GHQ_POINTS, wxp3, loglik);

	// coofs->coofs[0] = -d * A;
	coofs->coofs[0] = NAN;
	coofs->coofs[1] = -d * B * 0.5 * s2_inv;
	coofs->coofs[2] = -d * C * 0.25 * SQR(s2_inv);

	if (!workspace) {
		Free(x_user);
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_ai_vb_correct_mean(int thread_id, GMRFLib_density_tp ***density,	// need two types
			       int dens_count,
			       GMRFLib_density_tp **dens_local,
			       double *c,
			       double *d,
			       GMRFLib_ai_param_tp *ai_par,
			       GMRFLib_ai_store_tp *ai_store,
			       GMRFLib_graph_tp *graph,
			       GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_arg, GMRFLib_logl_tp *loglFunc, void *loglFunc_arg,
			       GMRFLib_preopt_tp *UNUSED(preopt))
{
	if (GMRFLib_inla_mode == GMRFLib_MODE_COMPACT) {
		// nothing to do here
		return GMRFLib_SUCCESS;
	} else {
		return GMRFLib_ai_vb_correct_mean_std(thread_id, density, dens_count, dens_local,
						      c, d, ai_par, ai_store, graph, Qfunc, Qfunc_arg, loglFunc, loglFunc_arg);
	}
}

int GMRFLib_ai_vb_correct_mean_preopt(int thread_id,
				      GMRFLib_density_tp ***density,
				      int dens_count,
				      double *UNUSED(c),
				      double *d,
				      GMRFLib_prior_mean_tp **prior_mean,
				      GMRFLib_ai_param_tp *ai_par,
				      GMRFLib_ai_store_tp *ai_store,
				      GMRFLib_graph_tp *graph,
				      GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_arg, GMRFLib_logl_tp *loglFunc, void *loglFunc_arg,
				      GMRFLib_preopt_tp *preopt, GMRFLib_idx_tp *d_idx)
{
	GMRFLib_ENTER_ROUTINE;
	Calloc_init(7 * graph->n + 2 * preopt->mnpred + 4 * preopt->Npred, 13);
	FILE *fp = (ai_par->fp_log ? ai_par->fp_log : stdout);
	int verbose = ai_par->vb_verbose && ai_par->fp_log;

#define SHOW_TIME(msg_)							\
	if (debug) {							\
		fprintf(fp, "[%1d] vb_preopt: %s %.3f\n", omp_get_thread_num(), msg_, GMRFLib_timer()-_tref); \
		_tref = GMRFLib_timer();				\
	}

	assert(GMRFLib_inla_mode == GMRFLib_MODE_COMPACT);

	// save time: only compute MM the first time, and keep MM and its factorisation fixed during the iterations. the motivation is that the
	// 2nd order properties will hardly change while the 1st order properties, ie the mean, will.
	int keep_MM = 1;

	int niter = ai_par->vb_iter_max;
	int debug = GMRFLib_DEBUG_IF();
	int emergency = 0;
	double one = 1.0, mone = -1.0, zero = 0.0;
	double max_correct = 5.001;
	double tref = GMRFLib_timer();
	GMRFLib_tabulate_Qfunc_tp *tabQ = NULL;

	if (!(ai_par->vb_enable && ai_par->vb_nodes_mean)) {
		GMRFLib_LEAVE_ROUTINE;
		return GMRFLib_SUCCESS;
	}

	// need the idx's for the vb correction and the data locations
	GMRFLib_idx_tp *vb_idx = NULL;

	for (int i = 0; i < graph->n; i++) {
		if (ai_par->vb_nodes_mean[i]) {
			GMRFLib_idx_add(&vb_idx, i);
		}
	}
	if (!vb_idx) {
		GMRFLib_LEAVE_ROUTINE;
		return GMRFLib_SUCCESS;
	}

	GMRFLib_ai_add_Qinv_to_ai_store(ai_store);	       /* add Qinv if required */

	double *sd = Calloc_get(graph->n);
	for (int i = 0; i < graph->n; i++) {
		double *var = GMRFLib_Qinv_get(ai_store->problem, i, i);
		sd[i] = (var ? sqrt(*var) : NAN);
	}

	double *x_mean = Calloc_get(graph->n);
	double *x_mean_orig = Calloc_get(graph->n);
	double *dx = Calloc_get(graph->n);
	double *pmean = Calloc_get(preopt->mnpred);
	double *pvar = Calloc_get(preopt->mnpred);

	for (int i = 0; i < graph->n; i++) {
		if (density[i][dens_count]) {
			x_mean[i] = density[i][dens_count]->user_mean;
		} else {
			x_mean[i] = ai_store->problem->mean_constr[i];
		}
	}
	Memcpy(x_mean_orig, x_mean, graph->n * sizeof(double));
	if (debug) {
		for (int i = 0; i < graph->n; i++) {
			printf("[%1d] x_mean[%1d] = %.12g\n", omp_get_thread_num(), i, x_mean[i]);
		}
	}

	double *tmp = Calloc_get(graph->n);
	// matrix with Cov(), alloc later as we can 'save one' large alloc
	gsl_matrix *M = NULL;				       // gsl_matrix_alloc(graph->n, vb_idx->n); 
	gsl_matrix *QM = gsl_matrix_alloc(graph->n, vb_idx->n);
	gsl_vector *B = gsl_vector_alloc(graph->n);
	gsl_matrix *MM = gsl_matrix_alloc(vb_idx->n, vb_idx->n);
	gsl_permutation *perm = gsl_permutation_alloc(vb_idx->n);
	gsl_vector *MB = gsl_vector_alloc(vb_idx->n);
	gsl_vector *delta = gsl_vector_alloc(vb_idx->n);
	gsl_vector *delta_prev = NULL;
	gsl_vector *delta_mu = gsl_vector_alloc(graph->n);

	double *BB = Calloc_get(preopt->Npred);
	double *CC = Calloc_get(preopt->Npred);

	double *prior_mean_tmp = Calloc_get(graph->n);
	double *prior_mean_fix = Calloc_get(graph->n);
	GMRFLib_prior_mean_get(thread_id, prior_mean_fix, graph->n, prior_mean);

	GMRFLib_preopt_tp *a = (GMRFLib_preopt_tp *) Qfunc_arg;
	GMRFLib_tabulate_Qfunc_tp *prior = NULL;
	assert(a == preopt && Qfunc == GMRFLib_preopt_Qfunc);

	// need to store these values as we change them, and so we need to set them back afterwards
	double *like_b_save = Calloc_get(preopt->Npred);
	double *like_c_save = Calloc_get(preopt->Npred);
	Memcpy(like_b_save, preopt->like_b[thread_id], preopt->Npred * sizeof(double));
	Memcpy(like_c_save, preopt->like_c[thread_id], preopt->Npred * sizeof(double));

	GMRFLib_tabulate_Qfunc_core(thread_id, &tabQ, graph, Qfunc, Qfunc_arg, NULL, 1);
	GMRFLib_tabulate_Qfunc_core(thread_id, &prior, preopt->latent_graph, GMRFLib_preopt_Qfunc_prior, Qfunc_arg, NULL, 1);
	gsl_matrix_set_zero(QM);

	int hessian_update = ai_par->vb_hessian_update;
	assert(hessian_update > 0);

	M = gsl_matrix_alloc(graph->n, vb_idx->n);
#define CODE_BLOCK							\
	for (int jj = 0; jj < vb_idx->n; jj++) {			\
		CODE_BLOCK_WORK_ZERO(0);				\
		int j = vb_idx->idx[jj];				\
		double *b = CODE_BLOCK_WORK_PTR(0);			\
		double *cov = CODE_BLOCK_WORK_PTR(1);			\
		b[j] = 1.0;						\
		GMRFLib_Qsolve(cov, b, ai_store->problem, j);		\
		for (int i = 0; i < graph->n; i++) {			\
			gsl_matrix_set(M, i, jj, cov[i]);		\
		}							\
	}
	RUN_CODE_BLOCK(IMIN(vb_idx->n, GMRFLib_MAX_THREADS()), 2, graph->n);
#undef CODE_BLOCK

	if (0) {
		// runs slower...
		double *b = Calloc(vb_idx->n * graph->n, double);
		double *cov = Calloc(vb_idx->n * graph->n, double);

		for (int jj = 0; jj < vb_idx->n; jj++) {
			int j = vb_idx->idx[jj];
			double *b_ = b + jj * graph->n;
			b_[j] = 1.0;
		}
		GMRFLib_Qsolves(cov, b, vb_idx->n, ai_store->problem);
		Free(b);

		M = gsl_matrix_alloc(graph->n, vb_idx->n);
		for (int jj = 0; jj < vb_idx->n; jj++) {
			double *cov_ = cov + jj * graph->n;
			for (int i = 0; i < graph->n; i++) {
				gsl_matrix_set(M, i, jj, cov_[i]);
			}
		}
		Free(cov);
	}

	double dxs[niter];
	GMRFLib_fill(niter, 0.0, dxs);

	for (int iter = 0; iter < niter + 1; iter++) {
		int update_MM = ((iter + 1 <= hessian_update) || (iter >= 2 && (dxs[iter - 1] > dxs[iter - 2])) || !keep_MM);
		double err_dx = 0.0;

		dxs[iter] = 0.0;
		gsl_vector_set_zero(B);
		gsl_vector_set_zero(MB);
		if (iter > 0) {
			if (delta_prev) {
				gsl_vector_free(delta_prev);
			}
			delta_prev = GMRFLib_gsl_duplicate_vector(delta);
		}
		gsl_vector_set_zero(delta);
		gsl_vector_set_zero(delta_mu);
		if (update_MM) {
			gsl_matrix_set_zero(MM);
		}

		// no need to compute the variance more than once since we're doing just the mean correction
		if (iter == 0) {
			GMRFLib_preopt_predictor_moments(pmean, pvar, preopt, ai_store->problem, x_mean);
		}
		// I know I compute the mean twice for iter=0, but then the timing gets right
		GMRFLib_preopt_predictor_moments(pmean, NULL, preopt, ai_store->problem, x_mean);

#define CODE_BLOCK							\
		for (int ii = 0; ii < d_idx->n; ii++) {			\
			int i = d_idx->idx[ii];				\
			GMRFLib_vb_coofs_tp vb_coof = {.coofs = {NAN, NAN, NAN}}; \
			GMRFLib_ai_vb_prepare_mean(thread_id, &vb_coof, i, d[i], loglFunc, loglFunc_arg, x_mean, pmean[i], sqrt(pvar[i]), CODE_BLOCK_WORK_PTR(0)); \
			if (debug) {					\
				fprintf(fp, "[%1d] i %d (mean,sd) = %.6f %.6f (A,B,C) = %.6f %.6f %.6f\n", omp_get_thread_num(), i, \
					pmean[i], sqrt(pvar[i]), vb_coof.coofs[0], vb_coof.coofs[1], vb_coof.coofs[2]); \
			}						\
			BB[i] = vb_coof.coofs[1];			\
			CC[i] = vb_coof.coofs[2];			\
			if (ISNAN(CC[i]) || ISNAN(BB[i])) {		\
				if (0) printf("idx %d CC <= 0, or BB or CC is NAN\n", i); \
				BB[i] = CC[i] = 0.0;			\
			}						\
		}

		RUN_CODE_BLOCK(IMIN(d_idx->n, GMRFLib_MAX_THREADS()), 1, 2 * GMRFLib_INT_GHQ_ALLOC_LEN);
#undef CODE_BLOCK

		GMRFLib_preopt_update(thread_id, preopt, BB, CC);
#pragma omp simd
		for (int ii = 0; ii < graph->n; ii++) {
			prior_mean_tmp[ii] = x_mean[ii] - prior_mean_fix[ii];
		}
		GMRFLib_Qx(thread_id, tmp, prior_mean_tmp, preopt->latent_graph, prior->Qfunc, prior->Qfunc_arg);

		for (int i = 0; i < graph->n; i++) {
			tmp[i] += preopt->total_b[thread_id][i];
			gsl_vector_set(B, i, tmp[i]);
		}

		if (update_MM) {
#define CODE_BLOCK	{						\
				GMRFLib_QM(thread_id, QM, M, graph, tabQ->Qfunc, tabQ->Qfunc_arg); \
				if (nt__ > 1) {				\
					gsl_blas_dgemm_omp(CblasTrans, CblasNoTrans, one, M, QM, zero, MM, nt__); \
				} else {				\
					gsl_blas_dgemm(CblasTrans, CblasNoTrans, one, M, QM, zero, MM);	\
				}					\
			}
			RUN_CODE_BLOCK_PLAIN(GMRFLib_MAX_THREADS(), 0, 0);
#undef CODE_BLOCK
		}

		gsl_blas_dgemv(CblasTrans, mone, M, B, zero, MB);
		if (debug) {
			GMRFLib_printf_gsl_matrix(stdout, M, "%.6f ");
			GMRFLib_printf_gsl_vector(stdout, B, "%.6f ");
		}

		// the system can be singular, like with intrinsic model components. its safe to invert the non-singular part only
		if (keep_MM) {
			// in this case, keep the inv of MM through the iterations
			if (update_MM) {
				GMRFLib_gsl_spd_inv(MM, GSL_ROOT3_DBL_EPSILON);
			}
			if (debug) {
				printf("MM\n");
				GMRFLib_printf_gsl_matrix(stdout, MM, "%.6f ");
				printf("MB\n");
				GMRFLib_printf_gsl_vector(stdout, MB, "%.6f ");
			}

			// hence just do inv(MM) %*% MB
			gsl_blas_dgemv(CblasNoTrans, one, MM, MB, zero, delta);
		} else {
			// solve MM %*% delta = MB
			GMRFLib_gsl_safe_spd_solve(MM, MB, delta, GSL_ROOT3_DBL_EPSILON);
		}

		int delta_is_NAN = 0;
		for (int i = 0; i < (int) delta->size; i++) {
			if (ISNAN(gsl_vector_get(delta, i))) {
				delta_is_NAN = i + 1;
				gsl_vector_set_zero(delta);
				break;
			}
		}

		int flag_cyclic = 0;
		if (iter > 0) {
			double mean_delta = 0.0;
			for (int i = 0; i < (int) delta->size; i++) {
				mean_delta += SQR(gsl_vector_get(delta, i) + gsl_vector_get(delta_prev, i));
			}
			mean_delta = sqrt(mean_delta / delta->size);
			if (mean_delta < FLT_EPSILON) {
				// take the half and then exit later
				flag_cyclic = 1;
				for (int i = 0; i < (int) delta->size; i++) {
					gsl_vector_set(delta, i, 0.5 * gsl_vector_get(delta, i));
				}
			}
		}

		gsl_blas_dgemv(CblasNoTrans, one, M, delta, zero, delta_mu);
		err_dx = 0.0;
		for (int i = 0; i < graph->n; i++) {
			dx[i] = gsl_vector_get(delta_mu, i);
			double adx = ABS(dx[i] / sd[i]);
			err_dx = DMAX(err_dx, adx);
			// truncate individual components
			if (adx > max_correct) {
				dx[i] = max_correct * sd[i] * DSIGN(dx[i]);
			}
		}
		dxs[iter] = err_dx;

		int diverge = 0;
		if (iter >= 3) {
			// 'diverge' is defined as increasing three times in a row
			diverge = ((dxs[iter - 0] > dxs[iter - 1]) && (dxs[iter - 1] > dxs[iter - 2]) && (dxs[iter - 2] > dxs[iter - 3]));
		}

		GMRFLib_daddto(graph->n, dx, x_mean);
		double max_correction = 0.0;
#pragma omp simd
		for (int i = 0; i < graph->n; i++) {
			max_correction = DMAX(max_correction, ABS(x_mean[i] - x_mean_orig[i]) / sd[i]);
		}

		int max_corr_flag = (max_correction >= ai_par->vb_emergency);

		if (max_corr_flag || delta_is_NAN || diverge) {
#pragma omp critical (Name_1169f76e685daed4d69fb5a745f9e95b4f5f633b)
			{
				if (delta_is_NAN) {
					fprintf(stderr, "\n\n\t***[%1d] warning *** delta[%1d] is NAN, 'vb.correction' is aborted\n", thread_id,
						delta_is_NAN - 1);
				}
				if (diverge) {
					fprintf(stderr,
						"\n\n\t***[%1d] warning *** iterative process seems to diverge, 'vb.correction' is aborted\n",
						thread_id);
				}
				if (max_corr_flag) {
					fprintf(stderr, "\n\n\t***[%1d] warning *** max_correction = %.2f >= %.2f, 'vb.correction' is aborted\n",
						thread_id, max_correction, ai_par->vb_emergency);
					fprintf(stderr, "\t*** You can change the emergency value (current value=%.2f) by \n",
						ai_par->vb_emergency);
					fprintf(stderr, "\t*** \t'control.inla=list(control.vb=list(emergency=...))'\n\n");
				}
				fprintf(stderr, "\t*** Please (re-)consider your model, priors, confounding, etc.\n");

				if (fp != stderr) {
					if (delta_is_NAN) {
						fprintf(fp, "\n\n\t***[%1d] warning *** delta[%1d] is NAN, 'vb.correction' is aborted\n",
							thread_id, delta_is_NAN - 1);
					}
					if (diverge) {
						fprintf(fp,
							"\n\n\t***[%1d] warning *** iterative process seems to diverge, 'vb.correction' is aborted\n",
							thread_id);
					}
					if (max_corr_flag) {
						fprintf(fp, "\n\n\t*** warning *** max_correction = %.2f >= %.2f, 'vb.correction' is aborted\n",
							max_correction, ai_par->vb_emergency);
						fprintf(fp, "\t*** You can change the emergency value (current value=%.2f) by \n",
							ai_par->vb_emergency);
						fprintf(fp, "\t*** \t'control.inla=list(control.vb=list(emergency=...))'\n\n");
					}
					fprintf(fp, "\t*** Please (re-)consider your model, priors, confounding, etc.\n");
				}
			}
			emergency = 1;
			break;
		}

		// this is RMS standardized change between the iterations, otherwise, just run the max iterations.
		// test this here so we can adjust the verbose output below
		int do_break = 0;
		if (err_dx < 0.01 || iter == niter - 1 || flag_cyclic) {
			do_break = 1;
		}

		if (verbose) {
#pragma omp critical (Name_d9343cf5e9cd69d222c869579102b5231d628874)
			{
				fprintf(fp, "\t[%1d]Iter [%1d/%1d] VB correct[MEAN] in total[%.2f sec/iter] cyclic[%s]\n",
					omp_get_thread_num(), iter, niter, (GMRFLib_timer() - tref) / (iter + 1.0),
					(flag_cyclic ? strdup("Yes") : strdup("No")));
				fprintf(fp, "\t\tNumber of nodes corrected for [%1d] max(dx/sd)[%.4f]\n", (int) delta->size, err_dx);
				if (do_break) {
					for (int jj = 0; jj < vb_idx->n; jj++) {
						int j = vb_idx->idx[jj];
						fprintf(fp, "\t\tNode[%1d] delta[%.3f] dx/sd[%.3f] (x-mode)/sd[%.3f]\n", j,
							gsl_vector_get(delta_mu, j), dx[j] / sd[j],
							(x_mean[j] - ai_store->problem->mean_constr[j]) / sd[j]);
					}
					fprintf(fp, "\t\tImplied correction for [%1d] nodes\n", preopt->mnpred + graph->n - vb_idx->n);
				}
			}
		}

		if (iter == niter) {
			// no covergence after 'niter' iterations, then we skip it
			emergency = 1;
		}

		if (do_break) {
			break;
		}
	}

	// we need to update those in any case
	GMRFLib_preopt_update(thread_id, preopt, like_b_save, like_c_save);

	// update the mean unless we're in an emergency
	if (!emergency) {
		for (int i = 0; i < graph->n; i++) {
			if (density[i][dens_count]) {
				GMRFLib_density_new_user_mean(density[i][dens_count], x_mean[i]);
			}
		}
	}

	GMRFLib_free_tabulate_Qfunc(tabQ);
	GMRFLib_free_tabulate_Qfunc(prior);
	gsl_matrix_free(M);
	gsl_matrix_free(MM);
	gsl_matrix_free(QM);
	gsl_permutation_free(perm);
	gsl_vector_free(B);
	gsl_vector_free(MB);
	gsl_vector_free(delta);
	gsl_vector_free(delta_mu);
	if (delta_prev) {
		gsl_vector_free(delta_prev);
	}

	GMRFLib_idx_free(vb_idx);
	Calloc_free();

	GMRFLib_LEAVE_ROUTINE;
#undef SHOW_TIME
	return GMRFLib_SUCCESS;
}

int GMRFLib_ai_vb_correct_variance_preopt(int thread_id,
					  GMRFLib_density_tp ***density,
					  int dens_count,
					  double *c,
					  double *d,
					  GMRFLib_ai_param_tp *ai_par,
					  GMRFLib_ai_store_tp *ai_store,
					  GMRFLib_graph_tp *graph,
					  GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_arg, GMRFLib_logl_tp *loglFunc, void *loglFunc_arg,
					  GMRFLib_preopt_tp *preopt, double *c_corrected, GMRFLib_idx_tp *d_idx)
{
	GMRFLib_ENTER_ROUTINE;
	assert(GMRFLib_inla_mode == GMRFLib_MODE_COMPACT);

	static double tref_a[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	static double count_tref_a = 0.0;		       /* only for testing */
	const int enable_tref_a = 0;			       /* only for testing */

	int hessian_diagonal = (ai_par->vb_hessian_strategy == GMRFLib_VB_HESSIAN_STRATEGY_DIAGONAL);
	int hessian_partial = (ai_par->vb_hessian_strategy == GMRFLib_VB_HESSIAN_STRATEGY_PARTIAL);
	int hessian_full = (ai_par->vb_hessian_strategy == GMRFLib_VB_HESSIAN_STRATEGY_FULL);
	int hessian_update = ai_par->vb_hessian_update;
	assert(hessian_update > 0);

	double tref = GMRFLib_timer();
	double grad_err = 0.0;
	int niter = ai_par->vb_iter_max;
	int debug = GMRFLib_DEBUG_IF();
	int tn = omp_get_thread_num();
	FILE *fp = (ai_par->fp_log ? ai_par->fp_log : stdout);
	int verbose = ai_par->vb_verbose && ai_par->fp_log;

	if (!(ai_par->vb_enable && ai_par->vb_nodes_variance)) {
		GMRFLib_LEAVE_ROUTINE;
		return GMRFLib_SUCCESS;
	}

	// need the idx's for the vb correction and the data locations
	GMRFLib_idx_tp *vb_idx = NULL;

	for (int i = 0; i < graph->n; i++) {
		if (ai_par->vb_nodes_variance[i]) {
			GMRFLib_idx_add(&vb_idx, i);
		}
	}
	if (!vb_idx) {
		GMRFLib_LEAVE_ROUTINE;
		return GMRFLib_SUCCESS;
	}

	assert(graph->n > 0);
	Calloc_init(7 * graph->n + 2 * preopt->mnpred + 4 * preopt->Npred + 1 * vb_idx->n, 13);
	double *x_mean = Calloc_get(graph->n);
	double *mean_constr = Calloc_get(graph->n);
	double *pmean = Calloc_get(preopt->mnpred);
	double *pvar = Calloc_get(preopt->mnpred);

#pragma omp simd
	for (int i = 0; i < graph->n; i++) {
		x_mean[i] = (density[i][dens_count] ? density[i][dens_count]->user_mean : ai_store->problem->mean_constr[i]);
	}
	Memcpy(mean_constr, ai_store->problem->mean_constr, graph->n * sizeof(double));

	double *BB = Calloc_get(preopt->Npred);
	double *CC = Calloc_get(preopt->Npred);

	GMRFLib_preopt_tp *a = (GMRFLib_preopt_tp *) Qfunc_arg;
	GMRFLib_tabulate_Qfunc_tp *prior = NULL;
	assert(a == preopt && Qfunc == GMRFLib_preopt_Qfunc);

	// need to store these values as we change them, and so we need to set them back afterwards
	double *like_b_save = Calloc_get(preopt->Npred);
	double *like_c_save = Calloc_get(preopt->Npred);
	Memcpy(like_b_save, preopt->like_b[thread_id], preopt->Npred * sizeof(double));
	Memcpy(like_c_save, preopt->like_c[thread_id], preopt->Npred * sizeof(double));

	GMRFLib_tabulate_Qfunc_core(thread_id, &prior, preopt->latent_graph, GMRFLib_preopt_Qfunc_prior, Qfunc_arg, NULL, 1);

	double *theta = Calloc_get(vb_idx->n);
	double *c_like = Calloc_get(graph->n);
	double *c_add = Calloc_get(graph->n);
	double *sd_prev = Calloc_get(graph->n);
	double *sd_orig = Calloc_get(graph->n);

	gsl_vector *delta = gsl_vector_alloc(vb_idx->n);
	gsl_vector *gradient = gsl_vector_alloc(vb_idx->n);
	gsl_matrix *hessian = gsl_matrix_alloc(vb_idx->n, vb_idx->n);
	gsl_vector_set_zero(delta);
	gsl_vector_set_zero(gradient);
	gsl_matrix_set_zero(hessian);

#define FUN(theta_) (theta_)
#define FUN_DERIV(theta_) (1.0)

	for (int i = 0; i < graph->n; i++) {
		sd_orig[i] = sd_prev[i] = sqrt(*GMRFLib_Qinv_get(ai_store->problem, i, i));
	}

	// the gradient of variance wrt theta (additive) is -var()^2, which sugguests c_like = 1/var()^2, but it seems a little to much. we add
	// a scaleing depending on the marginal stdev and normalize them to unit geometric mean.
	double csum = 0.0;
	for (int jj = 0; jj < vb_idx->n; jj++) {
		int j = vb_idx->idx[jj];
		c_like[j] = 1.0 / sd_orig[j];
		csum += log(c_like[j]);
	}
	csum = 1.0 / exp(csum / vb_idx->n);		       /* normalize the scalings */
	for (int jj = 0; jj < vb_idx->n; jj++) {
		int j = vb_idx->idx[jj];
		c_like[j] *= csum;
	}

	// can define shorter storage, but it should be fine really, I hope... this is properly yet to be verified
#if 0
#define STORAGE_TP float
#define TYPE_CAST (double)
#else
#define STORAGE_TP double
#define STORAGE_TYPE_CAST
#endif

	int storage_is_double = (sizeof(STORAGE_TP) == sizeof(double));
	static int print_once_only = 1;
	if (verbose && print_once_only) {
		print_once_only = 0;
		fprintf(fp, "\t[%1d]Cache cov(eta, latent)[%.1fMb] and cov(latent)[%.1fMb]\n", tn,
			d_idx->n * vb_idx->n * sizeof(STORAGE_TP) / SQR(1024.0), graph->n * vb_idx->n * sizeof(STORAGE_TP) / SQR(1024.0));
	}

	STORAGE_TP **cov_eta_latent_store = NULL;
	cov_eta_latent_store = Calloc(vb_idx->n, STORAGE_TP *);
	for (int ii = 0; ii < vb_idx->n; ii++) {
		cov_eta_latent_store[ii] = Calloc(d_idx->n, STORAGE_TP);
	}

	STORAGE_TP **cov_latent_store = Calloc(vb_idx->n, STORAGE_TP *);
	for (int ii = 0; ii < vb_idx->n; ii++) {
		cov_latent_store[ii] = Calloc(graph->n, STORAGE_TP);
	}

	double diff_sigma_max = 0.0;
	double diff_sigma_max_limit = 0.005;
	GMRFLib_problem_tp *problem = NULL;

	// main loop
	int iter = 0;
	for (iter = 0; iter < niter + 1; iter++) {	       /* yes, it should be '+1': we 'break' at the top */

		if (enable_tref_a) {
			tref_a[0] -= GMRFLib_timer();
		}
		// first we need to define a new 'problem' with the current values of thetas
		Memcpy(c_add, c, graph->n * sizeof(double));
		for (int jj = 0; jj < vb_idx->n; jj++) {
			int j = vb_idx->idx[jj];
			c_add[j] += c_like[j] * FUN(theta[jj]);
		}
		if (iter == 0) {
			// then we can use this one that is already computed
			problem = ai_store->problem;
		} else {
			// with iter==1 we know we used ai_store->problem, which should not be free'd, so we need to wait until iter>1
			if (iter > 1) {
				GMRFLib_free_problem(problem);
			}
			GMRFLib_error_handler_tp *old_handler = GMRFLib_set_error_handler_off();
			int retval = GMRFLib_init_problem(thread_id, &problem, NULL, NULL, c_add, x_mean, graph, Qfunc, Qfunc_arg,
							  ai_store->problem->sub_constr);
			if (retval != GMRFLib_SUCCESS) {
				P(retval);
				problem = NULL;
				break;
			}
			GMRFLib_set_error_handler(old_handler);
			GMRFLib_Qinv(problem);
		}
		if (enable_tref_a) {
			tref_a[0] += GMRFLib_timer();
		}

		if (iter > 0) {
			diff_sigma_max = 0.0;
			for (int i = 0; i < graph->n; i++) {
				double sd_new = sqrt(*GMRFLib_Qinv_get(problem, i, i));
				// KLD
				double kld = 0.5 * (SQR(sd_new / sd_prev[i]) - 1.0 + 2.0 * log(sd_prev[i] / sd_new));
				diff_sigma_max = DMAX(diff_sigma_max, kld);
				sd_prev[i] = sd_new;
			}
			diff_sigma_max = sqrt(diff_sigma_max);
			int pass = (diff_sigma_max < diff_sigma_max_limit || iter == niter);

			if (verbose) {
#pragma omp critical (Name_665f8b032e79706ed21a89a78c139f0eac1f454a)
				{
					fprintf(fp, "\t[%1d]Iter [%1d/%1d] VB correct with strategy [VARIANCE] in total[%.3fsec]\n",
						tn, iter, niter, GMRFLib_timer() - tref);
					fprintf(fp, "\t\tNumber of nodes corrected for [%1d] |gradient|[%.4g] max.step.sigma[%.4g](%s)\n",
						(int) vb_idx->n, grad_err, diff_sigma_max, (pass ? "PASS" : "FAIL"));
					if (pass) {
						for (int jj = 0; jj < vb_idx->n; jj++) {
							int j = vb_idx->idx[jj];
							fprintf(fp, "\t\tNode[%1d] delta[%.4f] sd/sd.orig[%.4f] sd[%.4f]\n", j,
								c_like[j] * FUN(theta[jj]), sd_prev[j] / sd_orig[j], sd_prev[j]);
						}
						fprintf(fp, "\t\tImplied correction for [%1d] nodes\n", preopt->mnpred + graph->n - vb_idx->n);
					}
				}
			}
			if (pass) {
				break;
			}
		}
		GMRFLib_preopt_predictor_moments(pmean, pvar, preopt, problem, x_mean);

		if (enable_tref_a) {
			tref_a[1] -= GMRFLib_timer();
		}

#define CODE_BLOCK							\
		for (int ii = 0; ii < d_idx->n; ii++) {			\
			int i = d_idx->idx[ii];				\
			GMRFLib_vb_coofs_tp vb_coof = {.coofs = {NAN, NAN, NAN}}; \
			GMRFLib_ai_vb_prepare_variance(thread_id, &vb_coof, i, d[i], loglFunc, loglFunc_arg, x_mean, pmean[i], sqrt(pvar[i]), CODE_BLOCK_WORK_PTR(0)); \
			if (debug && 0) {				\
				fprintf(fp, "\t[%1d] i %d (mean,sd) = %.6f %.6f (A,B,C) = %.6f %.6f %.6f\n", tn, i, \
				       pmean[i], sqrt(pvar[i]), vb_coof.coofs[0], vb_coof.coofs[1], vb_coof.coofs[2]); \
			}						\
			BB[i] = vb_coof.coofs[1];			\
			CC[i] = DMAX(0.0, vb_coof.coofs[2]);		\
		}

		RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 1, 2 * GMRFLib_INT_GHQ_ALLOC_LEN);
#undef CODE_BLOCK

		if (enable_tref_a) {
			tref_a[1] += GMRFLib_timer();
		}

		GMRFLib_idxval_tp **A = NULL;
		if (preopt->pA_idxval) {
			A = preopt->pAA_idxval;
		} else {
			A = preopt->A_idxval;
		}

		int UNUSED(integer_one) = 1;

//#define COV_ETA_LATENT(value_, k_, cov_latent_) (value_) = GMRFLib_dot_product(A[k_], cov_latent_)
#define COV_ETA_LATENT(value_, k_, cov_latent_) GMRFLib_dot_product_INLINE(value_, A[k_], cov_latent_)

#define COMPUTE_COV_LATENT(cov_latent_, j_, b_)				\
		if (1) {						\
			int _j = j_;					\
			Memset(b_, 0, graph->n * sizeof(double));	\
			b_[_j] = 1.0;					\
			GMRFLib_Qsolve(cov_latent_, b_, problem, _j);	\
		}

		if (enable_tref_a) {
			tref_a[2] -= GMRFLib_timer();
		}

#define CODE_BLOCK							\
		for (int ii = 0; ii < vb_idx->n; ii++) {		\
			int i = vb_idx->idx[ii];			\
			double *cov_latent_i = CODE_BLOCK_WORK_PTR(0);	\
			double *b_i = CODE_BLOCK_WORK_PTR(1);		\
			COMPUTE_COV_LATENT(cov_latent_i, i, b_i);	\
									\
			if (storage_is_double) {			\
				Memcpy(cov_latent_store[ii], cov_latent_i, graph->n * sizeof(double)); \
			} else {					\
				for(int k = 0; k < graph->n; k++) {	\
					cov_latent_store[ii][k] = (typeof(STORAGE_TP)) cov_latent_i[k];	\
				}					\
			}						\
									\
			STORAGE_TP *cov_eta_latent_i = cov_eta_latent_store[ii]; \
			for (int kk = 0; kk < d_idx->n; kk++) {		\
				int k = d_idx->idx[kk];			\
				COV_ETA_LATENT(cov_eta_latent_i[kk], k, cov_latent_i); \
			}						\
		}

		RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 2, graph->n);
#undef CODE_BLOCK

		if (enable_tref_a) {
			tref_a[2] += GMRFLib_timer();
		}

#define CODE_BLOCK							\
		for (int ii = 0; ii < vb_idx->n; ii++) {		\
			if (enable_tref_a) {				\
				tref_a[3] -= GMRFLib_timer();		\
			}						\
			int i = vb_idx->idx[ii];			\
			STORAGE_TP *cov_latent_i = cov_latent_store[ii]; \
			STORAGE_TP *cov_latent_j = NULL;		\
			double *values= CODE_BLOCK_WORK_PTR(0);		\
			double param_correction_i = c_like[i] * FUN_DERIV(theta[ii]); \
			double mell = 0.0;				\
			STORAGE_TP *cov_eta_latent_i = cov_eta_latent_store[ii]; \
			_Pragma("GCC ivdep")				\
				for (int kk = 0; kk < d_idx->n; kk++) {	\
					int k = d_idx->idx[kk];		\
					double S_ki = STORAGE_TYPE_CAST cov_eta_latent_i[kk]; \
					mell -= BB[k] * SQR(S_ki);	\
				}					\
			double ldet = STORAGE_TYPE_CAST cov_latent_i[i]; \
			double trace0 = 0.0;				\
			double trace1 = 0.0;				\
			for (int j = 0; j < preopt->latent_graph->n; j++) { \
				prior->Qfunc(thread_id, j, -1, values, prior->Qfunc_arg); \
				trace0 += values[0] * (-SQR(STORAGE_TYPE_CAST cov_latent_i[j])); \
				double trace_tmp = 0.0;			\
				double *valuesp1 = values + 1;		\
				int *jj_a = preopt->latent_graph->lnbs[j]; \
				for (int jjj = 0; jjj < preopt->latent_graph->lnnbs[j]; jjj++) { \
					int jj = jj_a[jjj];		\
					trace_tmp -= valuesp1[jjj] * STORAGE_TYPE_CAST cov_latent_i[jj]; \
				}					\
				trace1 += 2.0 * STORAGE_TYPE_CAST cov_latent_i[j] * trace_tmp; \
			}						\
			gsl_vector_set(gradient, (size_t) ii, (mell + 0.5 * (trace0 + trace1 + ldet)) * param_correction_i); \
			if (enable_tref_a) {				\
				tref_a[3] += GMRFLib_timer();		\
			}						\
			if (debug) {					\
				fprintf(fp, "\t[%1d]Iter [%1d/%1d] gradient[%1d] = %f\n", tn, iter, niter, ii, gsl_vector_get(gradient, (size_t) ii)); \
			}						\
									\
			if (enable_tref_a) {				\
				tref_a[4] -= GMRFLib_timer();		\
			}						\
			if (iter < hessian_update) {			\
				int jj_upper = (hessian_diagonal ? ii + 1 : vb_idx->n); \
				for (int jj = ii; jj < jj_upper; jj++) { \
					int j = vb_idx->idx[jj];	\
					STORAGE_TP *cov_eta_latent_j = cov_eta_latent_store[jj]; \
					double twoC1 = 2.0 * STORAGE_TYPE_CAST cov_latent_i[j]; \
					cov_latent_j = cov_latent_store[jj]; \
					double param_correction_j = c_like[j] * FUN_DERIV(theta[jj]); \
					mell = 0.0;			\
					/* in the PARTIAL strategy, just use the diagonal from the likelihood term */ \
					if ((hessian_partial && (ii == jj)) || hessian_full) { \
						_Pragma("GCC ivdep")	\
							for (int kk = 0; kk < d_idx->n; kk++) { \
								int k = d_idx->idx[kk];	\
								double S_kikj = STORAGE_TYPE_CAST cov_eta_latent_i[kk] * STORAGE_TYPE_CAST cov_eta_latent_j[kk]; \
								mell += S_kikj * (CC[k] * S_kikj  + BB[k] * twoC1); \
							}		\
					}				\
					ldet = -SQR(STORAGE_TYPE_CAST cov_latent_i[j]); \
					trace0 = 0.0;			\
					trace1 = 0.0;			\
					for (int k = 0; k < preopt->latent_graph->n; k++) { \
						double C2 = STORAGE_TYPE_CAST cov_latent_i[k]; \
						double C3 = STORAGE_TYPE_CAST cov_latent_j[k]; \
						prior->Qfunc(thread_id, k, -1, values, prior->Qfunc_arg); \
						trace0 += values[0] * C2 * C3; \
						double trace_tmp = 0.0;	\
						double *valuesp1 = values + 1; \
						int *kk_a =  preopt->latent_graph->lnbs[k]; \
						for (int kkk = 0; kkk < preopt->latent_graph->lnnbs[k]; kkk++) { \
							int kk = kk_a[kkk]; \
							trace_tmp += valuesp1[kkk] * (C2 * STORAGE_TYPE_CAST cov_latent_j[kk] + C3 * STORAGE_TYPE_CAST cov_latent_i[kk]); \
						}			\
						trace1 += trace_tmp;	\
					}				\
					trace0 *= twoC1;		\
					trace1 *= twoC1;		\
					double val = (mell + 0.5 * (trace0 + trace1 + ldet)) * param_correction_i * param_correction_j; \
					if (ii == jj) {			\
						val = DMAX(0.0, val);	\
						gsl_matrix_set(hessian, (size_t) ii, (size_t) jj, val);	\
					} else {			\
						gsl_matrix_set(hessian, (size_t) ii, (size_t) jj, val);	\
						gsl_matrix_set(hessian, (size_t) jj, (size_t) ii, val);	\
					}				\
				}					\
			}						\
			if (enable_tref_a) {				\
				tref_a[4] += GMRFLib_timer();		\
			}						\
		}

		if (iter < hessian_update) {
			gsl_matrix_set_zero(hessian);
		}
		if (hessian_full && (iter < hessian_update)) {
			// as we in this case has a triagular double loop
			RUN_CODE_BLOCK_DYNAMIC(GMRFLib_MAX_THREADS(), 1, graph->n);
		} else {
			RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 1, graph->n);
		}
#undef CODE_BLOCK

		// GMRFLib_printf_gsl_matrix(stdout, hessian, "%.2g ");

		if (enable_tref_a) {
			tref_a[5] -= GMRFLib_timer();
		}

		grad_err = 0.0;
		for (int i = 0; i < vb_idx->n; i++) {
			grad_err += SQR(gsl_vector_get(gradient, (size_t) i));
		}
		grad_err = sqrt(grad_err / (double) vb_idx->n);

		if (iter < hessian_update) {
			GMRFLib_gsl_spd_inv(hessian, GSL_ROOT3_DBL_EPSILON);
		}

		double one = 1.0, zero = 0.0;
		gsl_blas_dgemv(CblasNoTrans, one, hessian, gradient, zero, delta);

		if (debug) {
			for (int ii = 0; ii < vb_idx->n; ii++) {
				fprintf(fp, "\t[%1d]Iter [%1d/%1d] delta[%1d] = %.12f\n", tn, iter, niter, ii, gsl_vector_get(delta, ii));
			}
		}
		for (int ii = 0; ii < vb_idx->n; ii++) {
			int i = vb_idx->idx[ii];
			theta[ii] -= gsl_vector_get(delta, ii);
			if (debug) {
				fprintf(fp, "\t[%1d]Iter [%1d/%1d] c_add[%1d] = %.12f\n", tn, iter, niter, ii, c_like[i] * FUN(theta[ii]));
			}
		}
		if (enable_tref_a) {
			tref_a[5] += GMRFLib_timer();
		}
	}

	if (c_corrected) {
		for (int ii = 0; ii < vb_idx->n; ii++) {
			int i = vb_idx->idx[ii];
			theta[ii] -= gsl_vector_get(delta, ii);
			c_corrected[i] = c_like[i] * FUN(theta[ii]);
			// printf("c_corrected[%1d] = %f\n", i, c_corrected[i]);
		}
	}

	if (enable_tref_a) {
		count_tref_a++;
		for (int i = 0; i < 6; i++) {
			printf("%s: accumulated time [%1d] = %.3f   mean time [%.8f]\n", __GMRFLib_FuncName, i, tref_a[i],
			       tref_a[i] / count_tref_a);
		}
	}

	GMRFLib_preopt_update(thread_id, preopt, like_b_save, like_c_save);
	if (problem) {
		GMRFLib_free_problem(ai_store->problem);
		ai_store->problem = problem;
		Memcpy(problem->mean_constr, mean_constr, graph->n * sizeof(double));
		for (int i = 0; i < graph->n; i++) {
			if (density[i][dens_count]) {
				GMRFLib_density_new_user_stdev(density[i][dens_count], sd_prev[i]);
			}
		}
	}

	if (cov_eta_latent_store) {
		for (int ii = 0; ii < vb_idx->n; ii++) {
			Free(cov_eta_latent_store[ii]);
		}
		Free(cov_eta_latent_store);
	}
	if (cov_latent_store) {
		for (int ii = 0; ii < vb_idx->n; ii++) {
			Free(cov_latent_store[ii]);
		}
		Free(cov_latent_store);
	}

#undef FUN
#undef FUN_DERIV

	gsl_matrix_free(hessian);
	gsl_vector_free(delta);
	gsl_vector_free(gradient);

	GMRFLib_free_tabulate_Qfunc(prior);
	GMRFLib_idx_free(vb_idx);
	Calloc_free();

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_ai_vb_fit_gaussian(int thread_id, double *ell, double *fitted_mean, double *fitted_prec, int idx, double d,
			       GMRFLib_logl_tp *loglFunc, void *loglFunc_arg, double *x_vec, double mean, double sd)
{
	/*
	 * fit a Gaussian to the log-likelihood using prior (mean,sd)
	 */

	if (ISZERO(d)) {
		return GMRFLib_SUCCESS;
	}

	static double *wp = NULL;
	static double *xp = NULL;
	static double *xp1 = NULL;
	static double *xp2 = NULL;
	static double *xp3 = NULL;
	static double *xp4 = NULL;
	static double *xp5 = NULL;

	const int np = 1 + 2 * GMRFLib_INT_GHQ_POINTS, nnp = np / 2L, nnp1 = nnp + 1;

	if (!wp) {
#pragma omp critical (Name_0ddd01862f572e8e2021d8c931021738790dccc7)
		{
			if (!wp) {
				double *wtmp = NULL;
				GMRFLib_ghq(&xp, &wtmp, np);   /* just give ptr to storage */
				int nn = GMRFLib_align((size_t) nnp1, sizeof(double));
				xp1 = Calloc(5 * nn, double);
				xp2 = xp1 + 1 * nn;
				xp3 = xp1 + 2 * nn;
				xp4 = xp1 + 3 * nn;
				xp5 = xp1 + 4 * nn;

				for (int i = 0; i < nnp1; i++) {
					double x = xp[nnp + i];
					double z2 = SQR(x);
					xp1[i] = x;	       // d mu
					xp2[i] = 0.5 * (z2 - 1.0);	// d var
					xp3[i] = z2 - 1.0;     // d mu mu
					xp4[i] = 0.25 * (3.0 - 6.0 * z2 + SQR(z2));	// d var var
					xp5[i] = 0.5 * x * (z2 - 3.0);	// d var mu
				}
				wp = wtmp;
			}
		}
	}


	int nn = GMRFLib_align((size_t) np, sizeof(double));
	double x_user[5 * nn];
	double *loglik = x_user + nn;
	double *wloglik = x_user + 2 * nn;
	double *wloglik_sym = x_user + 3 * nn;
	double *wloglik_asym = x_user + 4 * nn;
	// GMRFLib_spline_tp *spline = NULL;

	double fit_mean = mean, fit_log_var = log(SQR(sd)), prior_var_inv = 1.0 / SQR(sd), step = 0.0;
	int max_iter = 100;

	for (int iter = 0; iter < max_iter; iter++) {
		double s = exp(0.5 * fit_log_var), s2 = SQR(s);
		GMRFLib_daxpb(np, s, xp, fit_mean, x_user);
		loglFunc(thread_id, loglik, x_user, np, idx, x_vec, NULL, loglFunc_arg, NULL);

		// don't know if I should interpolate and use that one
		// spline = GMRFLib_spline_create_x(x_user, loglik, np, GMRFLib_INTPOL_TRANS_NONE, GMRFLib_INTPOL_CACHE_SIMPLE);
		// GMRFLib_spline_eval_x(np, x_user, spline, loglik);

		GMRFLib_mul(np, loglik, wp, wloglik);
		wloglik_sym[0] = wloglik[nnp];
		wloglik_asym[0] = wloglik[nnp];
		for (int i = 1; i < nnp1; i++) {
			double a = wloglik[nnp + i];
			double b = wloglik[nnp - i];
			wloglik_sym[i] = a + b;
			wloglik_asym[i] = a - b;
		}

		double G1, G2, H11, H12, H22, s_inv = 1.0 / s, s2_inv = SQR(s_inv);
		G1 = -d * s_inv * GMRFLib_ddot(nnp1, wloglik_asym, xp1);
		G2 = -d * s2_inv * GMRFLib_ddot(nnp1, wloglik_sym, xp2);
		H11 = -d * s2_inv * GMRFLib_ddot(nnp1, wloglik_sym, xp3);
		H22 = -d * SQR(s2_inv) * GMRFLib_ddot(nnp1, wloglik_sym, xp4);
		H12 = -d * s2_inv * s_inv * GMRFLib_ddot(nnp1, wloglik_asym, xp5);

		// convert G2, H22, H12 to gradients/hessians wrt log(var) instead of var
		// diff(f(exp(y)),y);
		// diff(f(exp(y)),y,y);
		G2 *= s2;
		H22 = H22 * SQR(s2) + G2;
		H12 *= s2;

		// add contribtion from KLD (now we are in (mu, var=exp(theta)) parameterisation)
		G1 += (fit_mean - mean) * prior_var_inv;
		G2 += 0.5 * (s2 * prior_var_inv - 1.0);
		H11 += prior_var_inv;
		H22 += 0.5 * s2 * prior_var_inv;
		// H12 += 0.0;

		double step_len = DMIN(1.0, (1.0 + iter) / 1.0);
		double idet = 1.0 / (H11 * H22 - SQR(H12));

		// spell out explictely the inverse of the 2x2 Hessian times gradient
		double d_fit_mean = -idet * (H22 * G1 - H12 * G2);
		double d_fit_log_var = -idet * (-H12 * G1 + H11 * G2);
		fit_mean += step_len * d_fit_mean;
		fit_log_var += step_len * d_fit_log_var;

		step = sqrt((SQR(d_fit_mean) + SQR(d_fit_log_var)) / 2.0);
		if (idx == 0)
			printf("idx=%1d iter %d diff (%.12g, %.12g) fit(mean=%.12g, sd=%.12g) step %.12g\n",
			       idx, iter, d_fit_mean, d_fit_log_var, fit_mean, exp(0.5 * fit_log_var), step);
		if (step < 1.0E-5)
			break;
	}

	*fitted_prec = 1.0 / exp(fit_log_var);
	*fitted_mean = fit_mean;
	*ell = GMRFLib_dsum(np, wloglik);

	return GMRFLib_SUCCESS;
}

int GMRFLib_ai_store_config_preopt(int thread_id, GMRFLib_ai_misc_output_tp *mo, int ntheta, double *theta, double log_posterior,
				   double log_posterior_orig, GMRFLib_problem_tp *problem, double *mean_corrected,
				   GMRFLib_preopt_tp *preopt, GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_arg, double *cpodens_moments,
				   double *gcpodens_moments, char **arg_str, double *ll_info, double *lpred_mean, double *lpred_variance,
				   double *c_corrected)
{
	if (!mo || !(mo->configs_preopt)) {
		return GMRFLib_SUCCESS;
	}
	int id = omp_get_thread_num();
	int lite = mo->config_lite;

	if (!(mo->configs_preopt[id])) {

		GMRFLib_graph_tp *g = NULL;
		mo->configs_preopt[id] = Calloc(1, GMRFLib_store_configs_preopt_tp);

		mo->configs_preopt[id]->lite = lite;
		mo->configs_preopt[id]->mpred = preopt->mpred;
		mo->configs_preopt[id]->npred = preopt->npred;
		mo->configs_preopt[id]->mnpred = preopt->mnpred;
		mo->configs_preopt[id]->Npred = preopt->Npred;
		mo->configs_preopt[id]->n = preopt->n;
		mo->configs_preopt[id]->ntheta = ntheta;

		int nelm;
		nelm = preopt->preopt_graph->n + (lite ? 0 : preopt->preopt_graph->nnz);
		mo->configs_preopt[id]->nz = (nelm - mo->configs_preopt[id]->n) / 2 + mo->configs_preopt[id]->n;

		nelm = preopt->latent_graph->n + (lite ? 0 : preopt->latent_graph->nnz);
		mo->configs_preopt[id]->prior_nz = (nelm - mo->configs_preopt[id]->n) / 2 + mo->configs_preopt[id]->n;

		mo->configs_preopt[id]->A = preopt->A;
		mo->configs_preopt[id]->pA = preopt->pA;
		GMRFLib_duplicate_constr(&(mo->configs_preopt[id]->constr), preopt->latent_constr, preopt->preopt_graph);

		int *i = NULL, *j = NULL, ii, jj, k, kk;

		i = mo->configs_preopt[id]->i = Calloc(mo->configs_preopt[id]->nz, int);
		j = mo->configs_preopt[id]->j = Calloc(mo->configs_preopt[id]->nz, int);
		g = preopt->preopt_graph;
		if (lite) {
			for (k = 0; k < g->n; k++) {
				i[k] = k;
				j[k] = k;
			}
		} else {
			for (ii = k = 0; ii < g->n; ii++) {
				i[k] = ii;
				j[k] = ii;
				k++;
				for (kk = 0; kk < g->lnnbs[ii]; kk++) {
					jj = g->lnbs[ii][kk];
					i[k] = ii;
					j[k] = jj;
					k++;
				}
			}
		}

		i = mo->configs_preopt[id]->iprior = Calloc(mo->configs_preopt[id]->prior_nz, int);
		j = mo->configs_preopt[id]->jprior = Calloc(mo->configs_preopt[id]->prior_nz, int);
		g = preopt->latent_graph;
		if (lite) {
			for (k = 0; k < g->n; k++) {
				i[k] = k;
				j[k] = k;
			}
		} else {
			for (ii = k = 0; ii < g->n; ii++) {
				i[k] = ii;
				j[k] = ii;
				k++;
				for (kk = 0; kk < g->lnnbs[ii]; kk++) {
					jj = g->lnbs[ii][kk];
					i[k] = ii;
					j[k] = jj;
					k++;
				}
			}
		}

		mo->configs_preopt[id]->nconfig = 0;
		mo->configs_preopt[id]->config = NULL;
	}

	int nconfig = mo->configs_preopt[id]->nconfig;
	mo->configs_preopt[id]->config = Realloc(mo->configs_preopt[id]->config, nconfig + 1, GMRFLib_store_config_preopt_tp *);
	mo->configs_preopt[id]->config[nconfig] = Calloc(1, GMRFLib_store_config_preopt_tp);

	int ii, jj, k, kk;
	double *Qinv = NULL, *Q = NULL, *Qprior = NULL, *mean = NULL, *imean = NULL;
	GMRFLib_graph_tp *g = NULL;

	// TODO This should be changed: The code below to store Q and Qprior, is a little sloppy, as we do not use the (i,j) arrays
	// mo->configs_preopt[id]->iprior etc, but we just store the matrices in a vector, in the same order as the indices were stored.

	Q = Calloc(mo->configs_preopt[id]->nz, double);
	g = preopt->preopt_graph;
	if (lite) {
		for (ii = k = 0; ii < g->n; ii++) {
			Q[k++] = Qfunc(thread_id, ii, ii, NULL, Qfunc_arg) + (c_corrected ? c_corrected[ii] : 0.0);
		}
	} else {
		for (ii = k = 0; ii < g->n; ii++) {
			Q[k++] = Qfunc(thread_id, ii, ii, NULL, Qfunc_arg) + (c_corrected ? c_corrected[ii] : 0.0);
			for (kk = 0; kk < g->lnnbs[ii]; kk++) {
				jj = g->lnbs[ii][kk];
				Q[k++] = Qfunc(thread_id, ii, jj, NULL, Qfunc_arg);
			}
		}
	}

	Qprior = Calloc(mo->configs_preopt[id]->prior_nz, double);
	assert(Qfunc == GMRFLib_preopt_Qfunc);
	g = preopt->latent_graph;
	if (lite) {
		for (ii = k = 0; ii < g->n; ii++) {
			Qprior[k++] = GMRFLib_preopt_Qfunc_prior(thread_id, ii, ii, NULL, Qfunc_arg);
		}
	} else {
		for (ii = k = 0; ii < g->n; ii++) {
			Qprior[k++] = GMRFLib_preopt_Qfunc_prior(thread_id, ii, ii, NULL, Qfunc_arg);
			for (kk = 0; kk < g->lnnbs[ii]; kk++) {
				jj = g->lnbs[ii][kk];
				Qprior[k++] = GMRFLib_preopt_Qfunc_prior(thread_id, ii, jj, NULL, Qfunc_arg);
			}
		}
	}

	mean = Calloc(g->n, double);
	imean = Calloc(g->n, double);
	Memcpy(mean, problem->mean_constr, g->n * sizeof(double));
	Memcpy(imean, mean_corrected, g->n * sizeof(double));

	Qinv = Calloc(mo->configs_preopt[id]->nz, double);
	for (k = 0; k < mo->configs_preopt[id]->nz; k++) {
		double *tmp = GMRFLib_Qinv_get(problem, mo->configs_preopt[id]->i[k], mo->configs_preopt[id]->j[k]);
		Qinv[k] = (tmp ? *tmp : NAN);
	}

	GMRFLib_store_config_preopt_tp *cfg = mo->configs_preopt[id]->config[mo->configs_preopt[id]->nconfig];

	cfg->Q = Q;
	cfg->Qinv = Qinv;
	cfg->Qprior = Qprior;
	cfg->mean = mean;
	cfg->improved_mean = imean;
	cfg->log_posterior = log_posterior;		       /* may include integration weights */
	cfg->log_posterior_orig = log_posterior_orig;	       /* do NOT include integration weights */
	cfg->cpodens_moments = cpodens_moments;
	cfg->gcpodens_moments = gcpodens_moments;
	cfg->ll_info = ll_info;
	cfg->lpred_mean = lpred_mean;
	cfg->lpred_variance = lpred_variance;
	if (ntheta) {
		cfg->theta = Calloc(ntheta, double);
		Memcpy(cfg->theta, theta, ntheta * sizeof(double));
	} else {
		cfg->theta = NULL;
	}
	cfg->arg_str = arg_str;
	mo->configs_preopt[id]->nconfig++;

	return GMRFLib_SUCCESS;
}

int GMRFLib_ai_compute_lincomb(GMRFLib_density_tp ***lindens, double **cross, int nlin, GMRFLib_lc_tp **Alin,
			       GMRFLib_ai_store_tp *ai_store, double *improved_mean, int lookup_tables)
{
	/*
	 * Compute the marginals for the linear combinations using just the Gaussians. The computations gets a bit messy since we will try to
	 * avoid dependency of n, in both a, in a^T x, and in remap(a). We only need the range of non-zero terms 'remap(a)' and the non-zero
	 * terms in 'a'.
	 */

	assert(ai_store);
	GMRFLib_problem_tp *problem = ai_store->problem;
	int *remap = problem->sub_sm_fact.remap;
	int n, nc = 0, one = 1;
	GMRFLib_density_tp **d = NULL;

	// yes, disable this with PARDISO as PARDISO do not have this feaure
	int disable_opt = (GMRFLib_smtp == GMRFLib_SMTP_PARDISO ? 1 : 0);

	typedef struct {
		double *v;
		int from_idx;
		int to_idx;
	} cross_tp;
	cross_tp *cross_store = NULL;

	if (GMRFLib_smtp == GMRFLib_SMTP_TAUCS || GMRFLib_smtp == GMRFLib_SMTP_BAND) {
		remap = problem->sub_sm_fact.remap;
	} else {
		remap = problem->sub_sm_fact.PARDISO_fact->pstore[GMRFLib_PSTORE_TNUM_REF]->perm;
	}
	assert(remap);

	int id = 0;
	GMRFLib_CACHE_SET_ID(id);

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

#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_inner)
	for (int i = 0; i < nlin; i++) {

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
		for (int j = 0; j < Alin[i]->n; j++) {
			a[Alin[i]->idx[j] - from_idx_a] = (double) Alin[i]->weight[j];
		}

		/*
		 * compute the first non-zero index (mapped) if not already there
		 */
		if (Alin[i]->tinfo[id].first_nonzero_mapped < 0) {
			int findx = n;

			for (int j = 0; j < Alin[i]->n; j++) {
				int k = remap[Alin[i]->idx[j]];
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
		for (int j = 0; j < Alin[i]->n; j++) {
			vv[remap[Alin[i]->idx[j]]] = (double) Alin[i]->weight[j];
		}
		GMRFLib_solve_l_sparse_matrix_special(vv, &(problem->sub_sm_fact), problem->sub_graph, from_idx, from_idx + len - 1, 1);
		v = Calloc(len, double);
		Memcpy(v, vv + from_idx, len * sizeof(double));
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
			for (int j = 0; j < nc; j++) {
				/*
				 * w = AA^T CONSTR_M 
				 */
				double *p = NULL, *pp = NULL, w, ww;

				w = ww = 0.0;
				p = &(problem->constr_m[j * n]);
				pp = &(problem->qi_at_m[j * n]);

				for (jj = 0; jj < Alin[i]->n; jj++) {
					int k = Alin[i]->idx[jj];
					weight = (double) Alin[i]->weight[jj];

					w += weight * p[k];
					ww += weight * pp[k];
				}
				var_corr += w * ww;
			}
		}

		mean = imean = 0.0;
		for (int j = 0; j < Alin[i]->n; j++) {
			int k = Alin[i]->idx[j];
			weight = (double) Alin[i]->weight[j];

			mean += weight * problem->mean_constr[k];
			imean += weight * improved_mean[k];
		}
		var = DMAX(DBL_EPSILON, var - var_corr);
		GMRFLib_density_create_normal(&d[i], (imean - mean) / sqrt(var), 1.0, mean, sqrt(var), lookup_tables);
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
		{
			int k, i;
			for (k = 0, i = 0; i < nlin; i++) {
				for (int j = i; j < nlin; j++) {
					arr[k].i = i;
					arr[k].j = j;
					k++;
				}
			}
			assert(k == klen);
		}

		/*
		 * this loop is quick in any case, so no need to make do it in parallel unless we have constraints ? 
		 */
		omp_set_num_threads(GMRFLib_openmp->max_threads_inner);
#pragma omp parallel for if(nc) num_threads(GMRFLib_openmp->max_threads_inner)
		for (int k = 0; k < klen; k++) {

			int i = arr[k].i;
			int j = arr[k].j;
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
		if (cross_store) {
			for (int i = 0; i < nlin; i++) {
				Free(cross_store[i].v);
			}
			Free(cross_store);
		}
	}

	*lindens = d;

	return GMRFLib_SUCCESS;
}

int GMRFLib_ai_correct_cpodens(double *logdens, double *x, int *n, GMRFLib_ai_param_tp *ai_par)
{
	/*
	 * remove 'local' maxima at the extremes and return !GMRFLib_SUCCESS if the density is monotone otherwise GMRFLib_SUCCESS. If maximum, or nearly, is on the
	 * extremes, then flag an error.
	 */
	int idx, i, j;
	char *code = Calloc(*n, char);
	double mode;

	for (i = j = 0; i < *n; i++) {
		if (!ISINF(logdens[i])) {
			logdens[j] = logdens[i];
			j++;
		}
	}
	*n = j;

	idx = *n;
	do {
		idx--;
	} while (idx > 0 && (logdens[idx] > logdens[idx - 1]));

	if (idx < *n)
		Memset(&code[idx], 1, *n - idx);

	idx = -1;
	do {
		idx++;
	} while (idx < *n && (logdens[idx] > logdens[idx + 1]));

	if (idx > 0)
		Memset(&code[0], 1, idx);

	for (idx = i = 0; i < *n; i++) {
		if (!code[i]) {
			logdens[idx] = logdens[i];
			x[idx] = x[i];
			idx++;
		}
	}
	*n = idx + 1;
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

double GMRFLib_ai_cpopit_integrate(int thread_id, double *cpo, double *pit, int idx, GMRFLib_density_tp *cpo_density, double d,
				   GMRFLib_logl_tp *loglFunc, void *loglFunc_arg, double *x_vec)
{
	/*
	 * cpo_density is the marginal for x_idx without y_idx, density: is the marginal for x_idx with y_idx.
	 */
	int retval, compute_cpo = 1, np = GMRFLib_INT_NUM_POINTS;
	double low, dx, dxi, *xp = NULL, *xpi = NULL, *dens = NULL, *prob = NULL, integral = 0.0, integral2 = 0.0, integral_one, *loglik = NULL;
	double fail = 0.0;

	static double *w = NULL;
	if (!w) {
#pragma omp critical (Name_743b7d82abb3dc313f542012c3a2640ccca29d15)
		if (!w) {
			double www[] = { 4.0, 2.0 };
			double *ww = Calloc(np, double);
			ww[0] = ww[np - 1] = 1.0;
			for (int i = 1, k = 0; i < np - 1; i++, k = (k + 1L) % 2L) {
				ww[i] = www[k];
			}
			w = ww;
		}
	}

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

	retval = loglFunc(thread_id, NULL, NULL, 0, idx, x_vec, NULL, loglFunc_arg, NULL);
	if (!(retval == GMRFLib_LOGL_COMPUTE_CDF)) {
		compute_cpo = 0;
	}

	GMRFLib_ASSERT_RETVAL(np > 3, GMRFLib_ESNH, 0.0);

	Calloc_init(5 * np, 5);
	xp = Calloc_get(np);
	xpi = Calloc_get(np);
	dens = Calloc_get(np);
	prob = Calloc_get(np);
	loglik = Calloc_get(np);

	dxi = (cpo_density->x_max - cpo_density->x_min) / (np - 1.0);
	low = GMRFLib_density_std2user(cpo_density->x_min, cpo_density);
	dx = (GMRFLib_density_std2user(cpo_density->x_max, cpo_density) - low) / (np - 1.0);

	xp[0] = low;
	xpi[0] = cpo_density->x_min;
#pragma omp simd
	for (int i = 1; i < np; i++) {
		xp[i] = xp[0] + i * dx;
	}
#pragma omp simd
	for (int i = 1; i < np; i++) {
		xpi[i] = xpi[0] + i * dxi;
	}
	GMRFLib_evaluate_ndensity(dens, xpi, np, cpo_density);

	if (compute_cpo) {
		loglFunc(thread_id, prob, xp, -np, idx, x_vec, NULL, loglFunc_arg, NULL);	/* no correction for 'd' here; should we? */
	} else {
		Memset(prob, 0, np * sizeof(double));
	}
	loglFunc(thread_id, loglik, xp, np, idx, x_vec, NULL, loglFunc_arg, NULL);
	GMRFLib_dscale(np, d, loglik);

	GMRFLib_mul(np, prob, dens, xp);
	GMRFLib_exp(np, loglik, loglik);
	GMRFLib_mul(np, loglik, dens, xpi);

	integral = GMRFLib_ddot(np, w, xp);
	integral2 = GMRFLib_ddot(np, w, xpi);
	integral_one = GMRFLib_ddot(np, w, dens);

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

	Calloc_free();
	return fail;
}

double GMRFLib_ai_po_integrate(int thread_id, double *po, double *po2, double *po3, int idx, GMRFLib_density_tp *po_density,
			       double d, GMRFLib_logl_tp *loglFunc, void *loglFunc_arg, double *x_vec)
{
	double fail = 0.0;
	double integral2 = 0.0, integral3 = 0.0, integral4 = 0.0;

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

	if (po_density->type == GMRFLib_DENSITY_TYPE_GAUSSIAN) {
		int np = GMRFLib_INT_GHQ_POINTS;
		double *xp = NULL, *wp = NULL;
		double mean = po_density->user_mean;
		double stdev = po_density->user_stdev;

		GMRFLib_ghq(&xp, &wp, np);

		Calloc_init(4 * np, 4);
		double *x = Calloc_get(np);
		double *ll = Calloc_get(np);
		double *mask = Calloc_get(np);
		double *ell = Calloc_get(np);

		GMRFLib_fill(np, 1.0, mask);
		GMRFLib_daxpb(np, stdev, xp, mean, x);
		loglFunc(thread_id, ll, x, np, idx, x_vec, NULL, loglFunc_arg, NULL);
		double dmax = GMRFLib_max_value(ll, np, NULL);
		double dmin = GMRFLib_min_value(ll, np, NULL);
		double limit = -0.5 * SQR(xp[0]);	       // prevent extreme values
		if (dmin - dmax < limit) {
			for (int i = 0; i < np; i++) {
				if (ll[i] - dmax < limit) {
					mask[i] = 0.0;
					ll[i] = 0.0;
				}
			}
		}

		integral3 = GMRFLib_ddot(np, ll, wp);
		GMRFLib_exp(np, ll, ell);
		GMRFLib_mul(np, ell, mask, ell);	       /* so that ell[i]=exp(ll[i])=0 if ll[i]=0 */
		integral2 = GMRFLib_ddot(np, ell, wp);
		GMRFLib_sqr(np, ll, ll);
		integral4 = GMRFLib_ddot(np, ll, wp);
		Calloc_free();
	} else {
		double low, dx, dxi, *xp = NULL, *xpi = NULL, *ldens = NULL, w[2] = { 4.0, 2.0 }, integral_one, *loglik = NULL;

		int np = GMRFLib_INT_NUM_POINTS;
		int npm = GMRFLib_INT_NUM_INTERPOL * np - (GMRFLib_INT_NUM_INTERPOL - 1);

		GMRFLib_ASSERT_RETVAL(np > 3, GMRFLib_ESNH, 0.0);
		Calloc_init(4 * np + 2 * npm, 6);
		xp = Calloc_get(np);
		xpi = Calloc_get(np);
		ldens = Calloc_get(np);
		loglik = Calloc_get(np);

		dxi = (po_density->x_max - po_density->x_min) / (np - 1.0);
		low = GMRFLib_density_std2user(po_density->x_min, po_density);
		dx = (GMRFLib_density_std2user(po_density->x_max, po_density) - low) / (np - 1.0);

		xp[0] = low;
		xpi[0] = po_density->x_min;
#pragma omp simd
		for (int i = 1; i < np; i++) {
			xp[i] = xp[0] + i * dx;
			xpi[i] = xpi[0] + i * dxi;
		}
		GMRFLib_evaluate_nlogdensity(ldens, xpi, np, po_density);
		loglFunc(thread_id, loglik, xp, np, idx, x_vec, NULL, loglFunc_arg, NULL);

		double *dens = Calloc_get(npm);
		double *llik = Calloc_get(npm);

		if (GMRFLib_INT_NUM_INTERPOL == 3) {
#pragma omp simd
			for (int i = 0; i < np - 1; i++) {
				llik[3 * i + 0] = loglik[i];
				llik[3 * i + 1] = (2.0 * loglik[i] + loglik[i + 1]) / 3.0;
				llik[3 * i + 2] = (loglik[i] + 2.0 * loglik[i + 1]) / 3.0;
			}
#pragma omp simd
			for (int i = 0; i < np - 1; i++) {
				dens[3 * i + 0] = exp(ldens[i]);
				dens[3 * i + 1] = exp((2.0 * ldens[i] + ldens[i + 1]) / 3.0);
				dens[3 * i + 2] = exp((ldens[i] + 2.0 * ldens[i + 1]) / 3.0);
			}
			llik[3 * (np - 2) + 3] = loglik[np - 1];
			dens[3 * (np - 2) + 3] = exp(ldens[np - 1]);
			assert(3 * (np - 2) + 3 == npm - 1);
		} else if (GMRFLib_INT_NUM_INTERPOL == 2) {
#pragma omp simd
			for (int i = 0; i < np - 1; i++) {
				llik[2 * i + 0] = loglik[i];
				llik[2 * i + 1] = (loglik[i] + loglik[i + 1]) / 2.0;
			}
#pragma omp simd
			for (int i = 0; i < np - 1; i++) {
				dens[2 * i + 0] = exp(ldens[i]);
				dens[2 * i + 1] = exp((ldens[i] + ldens[i + 1]) / 2.0);
			}
			llik[2 * (np - 2) + 2] = loglik[np - 1];
			dens[2 * (np - 2) + 2] = exp(ldens[np - 1]);
			assert(2 * (np - 2) + 2 == npm - 1);
		} else {
			assert(GMRFLib_INT_NUM_INTERPOL == 2 || GMRFLib_INT_NUM_INTERPOL == 3);
		}

		integral2 = exp(llik[0]) * dens[0] + exp(llik[npm - 1]) * dens[npm - 1];
		integral3 = llik[0] * dens[0] + llik[npm - 1] * dens[npm - 1];
		integral4 = SQR(llik[0]) * dens[0] + SQR(llik[npm - 1]) * dens[npm - 1];
		integral_one = dens[0] + dens[npm - 1];
		for (int i = 1, k = 0; i < npm - 1; i++, k = (k + 1) % 2) {
			integral2 += w[k] * exp(llik[i]) * dens[i];
			integral3 += w[k] * llik[i] * dens[i];
			integral4 += w[k] * SQR(llik[i]) * dens[i];
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
		Calloc_free();
	}


	if (po) {
		*po = exp(d * log(DMAX(DBL_EPSILON, integral2)));
	}
	if (po2) {
		*po2 = d * integral3;
	}
	if (po3) {
		*po3 = SQR(d) * integral4;
	}

	return fail;
}

double *GMRFLib_ai_dic_integrate(int thread_id, int idx, GMRFLib_density_tp *density, double d, GMRFLib_logl_tp *loglFunc, void *loglFunc_arg,
				 double *x_vec)
{
	/*
	 * compute the integral of -2*loglikelihood * density(x), wrt x. also return the saturated one
	 */
	double integral = 0.0;
	double integral_sat = 0.0;
	double sat_ll = inla_compute_saturated_loglik(thread_id, idx, loglFunc, x_vec, loglFunc_arg);

	if (density->type == GMRFLib_DENSITY_TYPE_GAUSSIAN) {
		int np = GMRFLib_INT_GHQ_POINTS;
		double *xp = NULL, *wp = NULL;
		double mean = density->user_mean;
		double stdev = density->user_stdev;

		GMRFLib_ghq(&xp, &wp, np);

		Calloc_init(3 * np, 3);
		double *x = Calloc_get(np);
		double *ll = Calloc_get(np);
		double *ll_sat = Calloc_get(np);

		GMRFLib_daxpb(np, stdev, xp, mean, x);
		loglFunc(thread_id, ll, x, np, idx, x_vec, NULL, loglFunc_arg, NULL);
		GMRFLib_daxpb(np, 1.0, ll, -sat_ll, ll_sat);

		double dmax = GMRFLib_max_value(ll, np, NULL);
		double dmin = GMRFLib_min_value(ll, np, NULL);
		double limit = -0.5 * SQR(xp[0]);	       // prevent extreme values
		if (dmin - dmax < limit) {
			for (int i = 0; i < np; i++) {
				if (ll[i] - dmax < limit) {
					ll[i] = ll_sat[i] = 0.0;
				}
			}
		}

		integral = -2.0 * d * GMRFLib_ddot(np, ll, wp);
		integral_sat = -2.0 * d * GMRFLib_ddot(np, ll_sat, wp);
		Calloc_free();
	} else {

		// THIS PART NEEDS TO BE REWRITTEN

		double low, dx, dxi, *xp = NULL, *xpi = NULL, *ldens = NULL, w[2] = { 4.0, 2.0 }, integral_one, *loglik = NULL;

		int np = GMRFLib_INT_NUM_POINTS;
		int npm = GMRFLib_INT_NUM_INTERPOL * np - (GMRFLib_INT_NUM_INTERPOL - 1);

		GMRFLib_ASSERT_RETVAL(np > 3, GMRFLib_ESNH, NULL);

		Calloc_init(4 * np + 3 * npm, 7);
		xp = Calloc_get(np);
		xpi = Calloc_get(np);
		ldens = Calloc_get(np);
		loglik = Calloc_get(np);

		dxi = (density->x_max - density->x_min) / (np - 1.0);
		low = GMRFLib_density_std2user(density->x_min, density);
		dx = (GMRFLib_density_std2user(density->x_max, density) - low) / (np - 1.0);

		xp[0] = low;
		xpi[0] = density->x_min;
		for (int i = 1; i < np; i++) {
			xp[i] = xp[0] + i * dx;
			xpi[i] = xpi[0] + i * dxi;
		}
		GMRFLib_evaluate_nlogdensity(ldens, xpi, np, density);
		loglFunc(thread_id, loglik, xp, np, idx, x_vec, NULL, loglFunc_arg, NULL);

		double *dens = Calloc_get(npm);
		double *llik = Calloc_get(npm);
		double *llik_sat = Calloc_get(npm);

		if (GMRFLib_INT_NUM_INTERPOL == 3) {
#pragma omp simd
			for (int i = 0; i < np - 1; i++) {
				llik[3 * i + 0] = loglik[i];
				llik[3 * i + 1] = (2.0 * loglik[i] + loglik[i + 1]) / 3.0;
				llik[3 * i + 2] = (loglik[i] + 2.0 * loglik[i + 1]) / 3.0;
			}
#pragma omp simd
			for (int i = 0; i < np - 1; i++) {
				llik_sat[3 * i + 0] = llik[3 * i + 0] - sat_ll;
				llik_sat[3 * i + 1] = llik[3 * i + 1] - sat_ll;
				llik_sat[3 * i + 2] = llik[3 * i + 2] - sat_ll;
			}
#pragma omp simd
			for (int i = 0; i < np - 1; i++) {
				dens[3 * i + 0] = exp(ldens[i]);
				dens[3 * i + 1] = exp((2.0 * ldens[i] + ldens[i + 1]) / 3.0);
				dens[3 * i + 2] = exp((ldens[i] + 2.0 * ldens[i + 1]) / 3.0);
			}
			llik[3 * (np - 2) + 3] = loglik[np - 1];
			llik_sat[3 * (np - 2) + 3] = loglik[np - 1] - sat_ll;

			dens[3 * (np - 2) + 3] = exp(ldens[np - 1]);
			assert(3 * (np - 2) + 3 == npm - 1);
		} else if (GMRFLib_INT_NUM_INTERPOL == 2) {
#pragma omp simd
			for (int i = 0; i < np - 1; i++) {
				llik[2 * i + 0] = loglik[i];
				llik[2 * i + 1] = (loglik[i] + loglik[i + 1]) / 2.0;
			}
#pragma omp simd
			for (int i = 0; i < np - 1; i++) {
				llik_sat[2 * i + 0] = llik[2 * i + 0] - sat_ll;
				llik_sat[2 * i + 1] = llik[2 * i + 1] - sat_ll;
			}
#pragma omp simd
			for (int i = 0; i < np - 1; i++) {
				dens[2 * i + 0] = exp(ldens[i]);
				dens[2 * i + 1] = exp((ldens[i] + ldens[i + 1]) / 2.0);
			}
			llik[2 * (np - 2) + 2] = loglik[np - 1];
			llik_sat[2 * (np - 2) + 2] = loglik[np - 1] - sat_ll;

			dens[2 * (np - 2) + 2] = exp(ldens[np - 1]);
			assert(2 * (np - 2) + 2 == npm - 1);
		} else {
			assert(GMRFLib_INT_NUM_INTERPOL == 2 || GMRFLib_INT_NUM_INTERPOL == 3);
		}

		// prevent extreme values
		double dmax = GMRFLib_max_value(llik, np, NULL);
		double limit = -0.5 * SQR(6.0);
		for (int i = 0; i < np; i++) {
			if (llik[i] - dmax < limit) {
				llik[i] = llik_sat[i] = dens[i] = 0.0;
			}
		}

		integral = llik[0] * dens[0] + llik[npm - 1] * dens[npm - 1];
		integral_sat = llik_sat[0] * dens[0] + llik_sat[npm - 1] * dens[npm - 1];
		integral_one = dens[0] + dens[npm - 1];
		for (int i = 1, k = 0; i < npm - 1; i++, k = (k + 1) % 2) {
			integral += w[k] * llik[i] * dens[i];
			integral_sat += w[k] * llik_sat[i] * dens[i];
			integral_one += w[k] * dens[i];
		}

		integral = -2.0 * d * (integral / integral_one);
		integral_sat = -2.0 * d * (integral_sat / integral_one);
		Calloc_free();
	}

	double *res = Calloc(2, double);
	res[0] = integral;
	res[1] = integral_sat;

	return res;
}

int GMRFLib_ai_cpo_free(GMRFLib_ai_cpo_tp *cpo)
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

int GMRFLib_ai_po_free(GMRFLib_ai_po_tp *po)
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

int GMRFLib_ai_add_Qinv_to_ai_store(GMRFLib_ai_store_tp *ai_store)
{
	if (!ai_store || !(ai_store->problem)) {
		return GMRFLib_SUCCESS;
	}

	if (!ai_store->problem->sub_inverse) {
		int n = ai_store->problem->n;

		taucs_ccs_matrix *L = ai_store->problem->sub_sm_fact.TAUCS_L;
		if (GMRFLib_taucs_sort_L && L) {
#define CODE_BLOCK							\
			for (int i = 0; i < n; i++) {			\
				int m = L->colptr[i + 1] - L->colptr[i]; \
				int j = L->colptr[i];			\
				my_sort2_id(L->rowind + j, (double *) L->values.d + j, m); \
			}

			RUN_CODE_BLOCK(4, 0, 0);
#undef CODE_BLOCK

			if (GMRFLib_opt_solve) {
				if (!(ai_store->problem->sub_sm_fact.TAUCS_LL)) {
					ai_store->problem->sub_sm_fact.TAUCS_LL = GMRFLib_ccs2crs(ai_store->problem->sub_sm_fact.TAUCS_L);
					taucs_crs_matrix *LL = ai_store->problem->sub_sm_fact.TAUCS_LL;
#define CODE_BLOCK							\
					for (int i = 0; i < n; i++) {	\
						int m = LL->rowptr[i + 1] - LL->rowptr[i]; \
						int j = LL->rowptr[i];	\
						my_sort2_id(LL->colind + j, (double *) LL->values.d + j, m); \
					}
					RUN_CODE_BLOCK(4, 0, 0);
#undef CODE_BLOCK
				}
			}
		}

		GMRFLib_Qinv(ai_store->problem);
		Free(ai_store->stdev);
		ai_store->stdev = Calloc(n, double);
		for (int i = 0; i < n; i++) {
			double *var = GMRFLib_Qinv_get(ai_store->problem, i, i);
			ai_store->stdev[i] = (var ? sqrt(*var) : 0.0);
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
			correction *= (GMRFLib_cdfnorm(x + dz / 2.0) - GMRFLib_cdfnorm(x - dz / 2.0)) / (dz * f * exp(-0.5 * SQR(x)));
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

int GMRFLib_ai_marginal_one_hyperparamter(GMRFLib_density_tp **density, int idx, int nhyper, int hyper_count, double *hyper_z,
					  double *hyper_ldens, double *theta_mode, gsl_vector *sqrt_eigen_values,
					  gsl_matrix *eigen_vectors, double *std_stdev_theta, double dz,
					  double *stdev_corr_pos, double *stdev_corr_neg,
					  GMRFLib_ai_interpolator_tp interpolator, GMRFLib_ai_param_tp *ai_par, double *covmat)
{
#define COV(i, j)  covmat[ (i) + (j)*nhyper ]
#define NEXTRA 19
	int i, j;
	double *points = NULL, *ldens_values = NULL, *theta_max = NULL, *theta_min = NULL, sd;
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
		double *dens = NULL, *theta_tmp = NULL, *npoints_j = NULL, ldens_max;

		npoints = 0;
		Calloc_init(5 * hyper_count, 5);
		points = Calloc_get(hyper_count);
		dens = Calloc_get(hyper_count);
		ldens_values = Calloc_get(hyper_count);
		theta_tmp = Calloc_get(hyper_count);
		npoints_j = Calloc_get(hyper_count);

		ldens_max = GMRFLib_max_value(hyper_ldens, hyper_count, NULL);
		for (i = 0; i < hyper_count; i++) {

			GMRFLib_ai_z2theta(theta_tmp, nhyper, theta_mode, &(hyper_z[i * nhyper]), sqrt_eigen_values, eigen_vectors);
			j = GMRFLib_which(theta_tmp[idx], points, npoints);
			if (j >= 0) {
				/*
				 * point we have already 
				 */
				dens[j] += exp(hyper_ldens[i] - ldens_max);
				npoints_j[j]++;
			} else {
				/*
				 * new point 
				 */
				points[npoints] = theta_tmp[idx];
				dens[npoints] = exp(hyper_ldens[i] - ldens_max);
				npoints_j[npoints] = 1.0;
				npoints++;
			}
		}

		if (nhyper > 1) {
			int new_npoints = 0;

			for (i = j = 0; i < npoints; i++) {
				if (npoints_j[i] > 0) {
					points[j] = points[i];
					dens[j] = dens[i];
					new_npoints++;
					j++;
				}
			}
			npoints = new_npoints;
		}

		sd = std_stdev_theta[idx];
		double ldens_min = log(GSL_DBL_EPSILON);
		for (i = 0; i < npoints; i++) {
			ldens_values[i] = log(dens[i]);
			ldens_values[i] = DMAX(ldens_values[i], ldens_min);
			points[i] = (points[i] - theta_mode[idx]) / sd;
		}
		GMRFLib_density_create(density, GMRFLib_DENSITY_TYPE_SCGAUSSIAN, npoints, points, ldens_values,
				       theta_mode[idx], std_stdev_theta[idx], GMRFLib_TRUE);
		Calloc_free();
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

		npoints = 51;
		double theta_fixed, *x = NULL, *xx = NULL, *xxx = NULL;

		GMRFLib_ghq_abscissas(&xx, npoints);
		xxx = Calloc(npoints + NEXTRA, double);
		Memcpy(xxx, xx, npoints * sizeof(double));
		Memcpy(xxx + npoints, extra_points, NEXTRA * sizeof(double));

		npoints += NEXTRA;
		QSORT_FUN((void *) xxx, (size_t) npoints, sizeof(double), GMRFLib_dcmp);
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

		int tmax = GMRFLib_MAX_THREADS();
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
#pragma omp critical (Name_8a7281e161252d30ef221c00a3554c5a81e762b4)
				{
					fprintf(stderr, "\n\tGMRFLib_ai_marginal_one_hyperparamter: warning:\n");
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
		Memcpy(a->theta, x, a->nhyper * sizeof(double));	/* ndim == nhyper */
	}

	GMRFLib_ai_theta2z(a->z, a->nhyper, a->theta_mode, a->theta, a->sqrt_eigen_values, a->eigen_vectors);
	switch (a->interpolator) {
	case GMRFLib_AI_INTERPOLATOR_CCD:
	case GMRFLib_AI_INTERPOLATOR_CCD_INTEGRATE:
	{
		val = GMRFLib_interpolator_ccd(a->nhyper, a->hyper_count, a->z, a->hyper_z, a->hyper_ldens, (void *) a);
	}
		break;

	case GMRFLib_AI_INTERPOLATOR_NEAREST:
	{
		val = GMRFLib_interpolator_nearest(a->nhyper, a->hyper_count, a->z, a->hyper_z, a->hyper_ldens, (void *) &(a->dz));
	}
		break;

	case GMRFLib_AI_INTERPOLATOR_LINEAR:
	{
		val = GMRFLib_interpolator_linear(a->nhyper, a->hyper_count, a->z, a->hyper_z, a->hyper_ldens, (void *) &(a->dz));
	}
		break;

	case GMRFLib_AI_INTERPOLATOR_QUADRATIC:
	{
		val = GMRFLib_interpolator_quadratic(a->nhyper, a->hyper_count, a->z, a->hyper_z, a->hyper_ldens, (void *) &(a->dz));
	}
		break;

	case GMRFLib_AI_INTERPOLATOR_WEIGHTED_DISTANCE:
	{
		val = GMRFLib_interpolator_wdistance(a->nhyper, a->hyper_count, a->z, a->hyper_z, a->hyper_ldens, (void *) &(a->dz));
	}
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


double GMRFLib_interpolator_ccd(int ndim, int UNUSED(nobs), double *x, double *UNUSED(xobs), double *UNUSED(yobs), void *arg)
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

GMRFLib_ai_store_tp *GMRFLib_duplicate_ai_store(GMRFLib_ai_store_tp *ai_store, int skeleton, int copy_ptr, int copy_pardiso_ptr)
{
	/*
	 * duplicate AI_STORE. 'skeleton' only duplicate 'required' features. 'copy_ptr' only copies pointers to some objects known to be 'read only'
	 */
#define DUPLICATE(name_, len_, tp_, skeleton_)				\
	if (1) {							\
		if (ai_store->name_ && (len_) && !(skeleton_)){		\
			new_ai_store->name_ = Calloc(len_, tp_);	\
			size_t len = (len_) * sizeof(tp_);		\
			assert(len < PTRDIFF_MAX);			\
			Memcpy(new_ai_store->name_, ai_store->name_, len); \
		} else {						\
			new_ai_store->name_ = NULL;			\
	 	}							\
	}

#define COPY(name_) new_ai_store->name_ = ai_store->name_

	GMRFLib_ENTER_ROUTINE;
	if (!ai_store) {
		GMRFLib_LEAVE_ROUTINE;
		return NULL;
	}
	GMRFLib_ai_store_tp *new_ai_store = Calloc(1, GMRFLib_ai_store_tp);
	int n = (ai_store->problem ? ai_store->problem->n : 0);
	int Npred = ai_store->Npred;

	new_ai_store->store = GMRFLib_duplicate_store(ai_store->store, skeleton, copy_ptr, copy_pardiso_ptr);
	new_ai_store->problem = GMRFLib_duplicate_problem(ai_store->problem, skeleton, copy_ptr, copy_pardiso_ptr);
	COPY(nidx);
	COPY(Npred);

	DUPLICATE(mode, n, double, 0);
	DUPLICATE(aa, Npred, double, skeleton);
	DUPLICATE(bb, Npred, double, skeleton);
	DUPLICATE(cc, Npred, double, skeleton);
	DUPLICATE(stdev, n, double, skeleton);
	DUPLICATE(correction_term, n, double, skeleton);
	DUPLICATE(derivative3, n, double, skeleton);
	DUPLICATE(derivative4, n, double, skeleton);
	DUPLICATE(correction_idx, n, int, skeleton);
	new_ai_store->d_idx = GMRFLib_idx_duplicate(ai_store->d_idx);

	char *tmp = Calloc(1, char);
	Free(tmp);

	GMRFLib_LEAVE_ROUTINE;
	return new_ai_store;

#undef DUPLICATE
#undef COPY
}

double GMRFLib_bfunc_eval(int thread_id, double *constant, GMRFLib_bfunc_tp *bfunc)
{
	// evaluate bfunc: b[idx] = sum_i Q[idx,i]*mean[i]. 'con' is the contribution to the constant m'(Q+c)m.

#define MAPIDX(_idx, _d) MOD(MOD(_idx, (_d)->n * (_d)->ngroup), (_d)->n)

	if (bfunc == NULL || bfunc->bdef == NULL || bfunc->idx < 0) {
		return 0.0;
	}

	double b = 0.0, mu0 = 0.0, mu = 0.0;
	int i, j, idx = bfunc->idx;
	GMRFLib_bfunc2_tp *d = bfunc->bdef;

	// fprintf(stderr, "idx %d mapidx %d n %d nr %d ng %d\n", idx, MAPIDX(idx, d), d->n, d->nreplicate, d->ngroup);
	mu0 = d->mfunc(thread_id, MAPIDX(idx, d), d->mfunc_arg);
	b = (mu0 ? (d->diagonal + d->Qfunc(thread_id, idx, idx, NULL, d->Qfunc_arg)) * mu0 : 0.0);
	for (i = 0; i < d->graph->nnbs[idx]; i++) {
		j = d->graph->nbs[idx][i];
		mu = d->mfunc(thread_id, MAPIDX(j, d), d->mfunc_arg);
		if (mu) {
			b += d->Qfunc(thread_id, idx, j, NULL, d->Qfunc_arg) * mu;
		}
	}

	*constant = (mu0 ? b * mu0 : 0.0);

#undef MAPIDX
	return b;
}

double GMRFLib_prior_mean_func_eval(int thread_id, GMRFLib_prior_mean_tp *prior_mean)
{
	// this code is extracted from bfunc...

#define MAPIDX(_idx, _d) MOD(MOD(_idx, (_d)->n * (_d)->ngroup), (_d)->n)
	GMRFLib_bfunc2_tp *d = prior_mean->bdef;
	if (d) {
		int idx = prior_mean->idx;
		return (d->mfunc(thread_id, MAPIDX(idx, d), d->mfunc_arg));
	} else {
		return (prior_mean->fixed_mean);
	}
#undef MAPIDX
	return 0.0;
}

int GMRFLib_prior_mean_get(int thread_id, double *pmean, int n, GMRFLib_prior_mean_tp **prior_mean)
{
	if (prior_mean) {
		for (int i = 0; i < n; i++) {
			pmean[i] = (prior_mean[i] ? GMRFLib_prior_mean_func_eval(thread_id, prior_mean[i]) : 0.0);
		}
	} else {
		Memset(pmean, 0, n * sizeof(double));
	}
	return GMRFLib_SUCCESS;
}
