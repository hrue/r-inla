
/* approx-inference.c
 * 
 * Copyright (C) 2006-2023 Havard Rue
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


static int pool_nhyper = -1;

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
	(*ai_par)->mode_known = 0;
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
	(*ai_par)->vb_nodes_mean = NULL;
	(*ai_par)->vb_nodes_variance = NULL;

	return GMRFLib_SUCCESS;
}

int GMRFLib_ai_param_duplicate(GMRFLib_ai_param_tp ** ai_par_new, GMRFLib_ai_param_tp * ai_par)
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
	fprintf(fp, "\t\tOptimise: try to be smart: %s\n", (ai_par->optimise_smart ? "Yes" : "No"));
	fprintf(fp, "\t\tOptimise: use directions: %s\n", (ai_par->optimise_use_directions ? "Yes" : "No"));
	fprintf(fp, "\t\tMode known: %s\n", (ai_par->mode_known ? "Yes" : "No"));
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
		fprintf(fp, "\t\tstrategy            = [%s]\n", (ai_par->vb_strategy == GMRFLib_AI_VB_MEAN ? "mean" :
								 (ai_par->vb_strategy ==
								  GMRFLib_AI_VB_VARIANCE ? "mean and variance" : "UNKNOWN")));
		fprintf(fp, "\t\tverbose             = [%s]\n", (ai_par->vb_verbose ? "Yes" : "No"));
		fprintf(fp, "\t\tf_enable_limit_mean = [%1d]\n", ai_par->vb_f_enable_limit_mean);
		fprintf(fp, "\t\tf_enable_limit_var  = [%1d]\n", ai_par->vb_f_enable_limit_variance);
		fprintf(fp, "\t\titer_max            = [%1d]\n", ai_par->vb_iter_max);
		fprintf(fp, "\t\temergency           = [%.2f]\n", ai_par->vb_emergency);
		fprintf(fp, "\t\thessian_update      = [%1d]\n", ai_par->vb_hessian_update);
		fprintf(fp, "\t\thessian_strategy    = [%s]\n", VB_HESSIAN_STRATEGY_NAME(ai_par->vb_hessian_strategy));
	} else {
		fprintf(fp, "\tVB-correction is [Disabled]\n");
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_ai_marginal_hyperparam(int thread_id,
				   double *logdens,
				   double *x, double *b, double *c, double *mean, double *d,
				   GMRFLib_logl_tp * loglFunc, void *loglFunc_arg,
				   GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg,
				   GMRFLib_constr_tp * constr, GMRFLib_ai_param_tp * ai_par, GMRFLib_ai_store_tp * ai_store,
				   GMRFLib_preopt_tp * preopt)
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
				nnr_step_factor_first_time_only = Calloc(GMRFLib_CACHE_LEN, int);
				for (int i = 0; i < GMRFLib_CACHE_LEN; i++) {
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
	optpar->abserr_step = ai_par->optpar_abserr_func;
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
			       (thread_id, &problem, ai_store->mode, b, c, mean, d, loglFunc, loglFunc_arg, graph, Qfunc,
				Qfunc_arg, constr, optpar, blockpar, ai_store->store, ai_store->aa, ai_store->bb, ai_store->cc,
				ai_par->gaussian_data, ai_par->cmin, ai_par->b_strategy, 0, preopt));
	} else {
		GMRFLib_EWRAP1(GMRFLib_init_GMRF_approximation_store__intern
			       (thread_id, &problem, x, b, c, mean, d, loglFunc, loglFunc_arg, graph, Qfunc, Qfunc_arg,
				constr, optpar, blockpar, ai_store->store, ai_store->aa, ai_store->bb, ai_store->cc,
				ai_par->gaussian_data, ai_par->cmin, ai_par->b_strategy, 0, preopt));
	}
	GMRFLib_ASSERT(problem, GMRFLib_EOPTNR);

	/*
	 * if store, then store the mode to use as the initial point at later calls 
	 */
	Free(ai_store->mode);
	ai_store->mode = Calloc(n, double);
	Memcpy(ai_store->mode, problem->mean_constr, n * sizeof(double));

	if (mean == NULL && 1) {
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

int GMRFLib_ai_log_posterior(int thread_id, double *logdens,
			     double *x, double *b, double *c, double *mean, double *d,
			     GMRFLib_logl_tp * loglFunc, void *loglFunc_arg,
			     GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg, GMRFLib_constr_tp * UNUSED(constr))
{
	/*
	 * compute the log posterior of configuration 'x' up to an additive constant and return the value in 'logdens'. 
	 */

	int n;
	double *xx = NULL, val, logl, result;

	GMRFLib_ENTER_ROUTINE;

	n = graph->n;
	xx = Calloc(n, double);				       /* xx = x - mean */

	if (mean) {
#pragma GCC ivdep
		for (int i = 0; i < n; i++) {
			xx[i] = x[i] - mean[i];
		}
	} else {
		Memcpy(xx, x, n * sizeof(double));
	}

	double tmp0 = 0.0, tmp1 = 0.0, tmp2 = 0.0, tmp3 = 0.0, tmp4 = 0.0;

	// c-term is now in here
	GMRFLib_xQx2(thread_id, &result, xx, graph, Qfunc, Qfunc_arg, c);
	tmp4 = -0.5 * result;

	/*
	 * add the linear term 
	 */
	if (b) {
		tmp1 = GMRFLib_ddot(n, b, x);
	}

	if (d) {
		/*
		 * do not include fixed points 
		 */
		if (1) {
			/*
			 * new code; better for omp 
			 */
			int *idxs = NULL, nidx;
			Calloc_init(n, 1);
			idxs = (int *) Calloc_get(n);
			assert(idxs);
			nidx = 0;

			// why isn't this precomputed as its the same all the time ?
			for (int ii = 0; ii < n; ii++) {
				if (d[ii]) {
					idxs[nidx++] = ii;
				}
			}

			tmp2 = 0.0;
#pragma omp parallel for reduction(+: tmp2) num_threads(GMRFLib_openmp->max_threads_inner)
			for (int iii = 0; iii < nidx; iii++) {
				int ii = idxs[iii];
				double ll = 0.0;
				loglFunc(thread_id, &ll, &x[ii], 1, ii, x, NULL, loglFunc_arg, NULL);
				tmp2 += d[ii] * ll;
			}
			Calloc_free();
		} else {
			/*
			 * old version 
			 */
			for (int ii = 0; ii < n; ii++) {
				if (d[ii]) {
					loglFunc(thread_id, &logl, &x[ii], 1, ii, x, NULL, loglFunc_arg, NULL);
					tmp2 += d[ii] * logl;
				}
			}
		}
	}

	val = tmp0 + tmp1 + tmp2 + tmp3 + tmp4;

	Free(xx);

	*logdens = val;

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_ai_log_posterior_restricted(int thread_id, double *logdens, double *x, double *x_mode, double *x_gradient, double delta,
					double *b, double *c, double *mean, double *d,
					GMRFLib_logl_tp * loglFunc, void *loglFunc_arg,
					GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg,
					GMRFLib_constr_tp * UNUSED(constr), GMRFLib_graph_tp * subgraph, GMRFLib_ai_store_tp * UNUSED(ai_store),
					int *node_map, double *ql)
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
	double xx, *f = NULL, *g = NULL, val, tmp, logl = 0.0, q_value;

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
			i = node_map[ii];
			xx = (mean ? x_mode[i] - mean[i] : x_mode[i]);

			if (c) {
				q_value = (Qfunc(thread_id, i, i, NULL, Qfunc_arg) + c[i]);
				f[ii] += x_gradient[i] * q_value;
				g[ii] += xx * q_value;
			} else {
				q_value = Qfunc(thread_id, i, i, NULL, Qfunc_arg);
				f[ii] += x_gradient[i] * q_value;
				g[ii] += xx * q_value;
			}
			if (mean) {
				for (jj = 0; jj < graph->nnbs[i]; jj++) {
					j = graph->nbs[i][jj];
					q_value = Qfunc(thread_id, i, j, NULL, Qfunc_arg);
					f[ii] += x_gradient[j] * q_value;
					g[ii] += (x_mode[j] - mean[j]) * q_value;
				}
			} else {
				for (jj = 0; jj < graph->nnbs[i]; jj++) {
					j = graph->nbs[i][jj];
					q_value = Qfunc(thread_id, i, j, NULL, Qfunc_arg);
					f[ii] += x_gradient[j] * q_value;
					g[ii] += x_mode[j] * q_value;
				}
			}
		}

		// ql[0] is linear, ql[1] is quadratic
		ql[0] = ql[1] = 0.0;
		for (ii = 0; ii < ns; ii++) {
			i = node_map[ii];
			ql[0] -= g[ii] * x_gradient[i];
			ql[1] += f[ii] * x_gradient[i];
		}
		if (b) {
			for (ii = 0; ii < ns; ii++) {
				i = node_map[ii];
				ql[0] += x_gradient[i] * b[i];
			}
		}
		Free(f);
		Free(g);
	} else {
		val = -0.5 * SQR(delta) * ql[1] + delta * ql[0];
		if (d) {
			tmp = 0.0;
			for (ii = 0; ii < ns; ii++) {
				i = node_map[ii];
				if (d[i]) {
					loglFunc(thread_id, &logl, &x[i], 1, i, x, NULL, loglFunc_arg, NULL);
					tmp += d[i] * logl;
				}
			}
			val += tmp;
		}
		*logdens = val;
	}

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_ai_marginal_hidden(int thread_id, GMRFLib_density_tp ** density, GMRFLib_density_tp ** cpo_density,
			       int lookup_tables,
			       int idx, double *x, double *b, double *c, double *mean, double *d,
			       GMRFLib_logl_tp * loglFunc, void *loglFunc_arg,
			       GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg,
			       GMRFLib_constr_tp * constr, GMRFLib_ai_param_tp * ai_par, GMRFLib_ai_store_tp * ai_store, GMRFLib_preopt_tp * preopt)
{
	/*
	 * compute the approximation to the marginal for the hidden field at index 'idx' and return the density in *density. if
	 * (cpo_density), then compute also the density when y_idx is removed.
	 */

	char *fix = NULL, *fixx = NULL;
	int i, j, k, nd = -1, n = -1, free_ai_par = 0, n_points, ii, free_ai_store = 0, i_idx, j_idx, one = 1, *node_map = NULL;
	double *x_points = NULL, x_sd, x_mean, *cond_mode = NULL, *fixed_mode = NULL, *log_density = NULL,
	    log_dens_cond = 0.0, deriv_log_dens_cond = 0.0, a, *derivative = NULL, *mean_and_variance = NULL, deldif =
	    GMRFLib_eps(1.0 / 6.0), inv_stdev, *cov = NULL, corr, corr_term, *covariances = NULL, alpha;

	GMRFLib_graph_tp *subgraph = NULL;
	GMRFLib_Qinv_tp *store_Qinv = NULL;
	GMRFLib_optimize_param_tp *optpar = NULL;
	GMRFLib_blockupdate_param_tp *blockpar = NULL;
	GMRFLib_store_tp *store = NULL;
	GMRFLib_ai_strategy_tp strategy;

	int Npred = (preopt ? preopt->Npred : graph->n);
	ai_store->Npred = Npred;

#define COMPUTE_CPO_DENSITY						\
	if (cpo_density) {						\
		if (d[idx]) {						\
			double *xp = NULL, *xp_tmp = NULL,		\
				*ld = NULL, *logcor = NULL, *x_user = NULL, _alpha=-1.0; \
			const int _debug =  0;				\
			int itry, flag, np, np_orig = GMRFLib_INT_GHQ_POINTS + 4, \
				_one = 1, _i, npx = 8, itmp, np_new = np_orig + 2*npx; \
			double cor_eps = GMRFLib_eps(0.75), cor_max, range;	\
									\
			Calloc_init(4*np_new, 4);			\
			for(itry = 0; itry < 2;	itry++)			\
			{						\
				np = np_orig;				\
				flag = 0;				\
				if (!ld) {				\
					ld = Calloc_get(np_new);	\
					logcor = Calloc_get(np_new);	\
					x_user = Calloc_get(np_new);	\
					xp = Calloc_get(np_new);	\
				}					\
				GMRFLib_ghq_abscissas(&xp_tmp, np);	\
				range = xp_tmp[np-1];			\
				Memcpy(xp+npx, xp_tmp, np*sizeof(double)); \
				for(itmp = 0; itmp < npx; itmp++) {	\
					xp[itmp] = xp[npx] - range * (npx - itmp)/(double)npx; \
					xp[np + npx + itmp] = xp[npx + np - 1]  + range * (itmp + 1.0)/(double)npx; \
				}					\
				np = np_new;				\
				if (_debug) {				\
					if (0) \
						for(itmp = 0; itmp < np; itmp++) \
							printf("xp[%1d] = %.3f\n", itmp, xp[itmp]); \
					GMRFLib_density_printf(stdout, *density); \
				}					\
				GMRFLib_evaluate_nlogdensity(ld, xp, np, *density); \
				GMRFLib_density_std2user_n(x_user, xp, np, *density); \
				loglFunc(thread_id, logcor, x_user, np, idx, fixed_mode, NULL, loglFunc_arg, NULL); \
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
							       (*density)->std_mean, (*density)->std_stdev, lookup_tables); \
					if (flag && cpo_density) GMRFLib_setbit(&((*cpo_density)->flags), DENSITY_FLAGS_FAILURE); \
				} else {				\
					*cpo_density = NULL;		\
				}					\
				if (*cpo_density || itry == 1)		\
					break;				\
			}						\
			Calloc_free();					\
		} else {						\
			if (ai_par->cpo_manual){			\
				GMRFLib_density_duplicate(cpo_density, *density); \
			} else {					\
				*cpo_density = NULL;			\
			}						\
		}							\
	}

	GMRFLib_ENTER_ROUTINE;

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
#pragma omp critical (Name_55d88ef833913b76c8f458812b76256a0492204c)
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
#pragma omp critical (Name_2e9be797bbf18c2c12a4e333407447077e7a0eae)
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
		ai_store->aa = Calloc(Npred, double);
		ai_store->bb = Calloc(Npred, double);
		ai_store->cc = Calloc(Npred, double);

		GMRFLib_EWRAP1(GMRFLib_init_GMRF_approximation_store__intern(thread_id,
									     &(ai_store->problem),
									     (ai_store->mode ? ai_store->mode : x),
									     b, c, mean, d, loglFunc, loglFunc_arg,
									     graph, Qfunc, Qfunc_arg, constr, optpar, blockpar,
									     ai_store->store, ai_store->aa, ai_store->bb,
									     ai_store->cc, ai_par->gaussian_data, ai_par->cmin,
									     ai_par->b_strategy, 0, preopt));
		GMRFLib_ai_add_Qinv_to_ai_store(ai_store);
		/*
		 * store the mode for later usage 
		 */
		Free(ai_store->mode);
		ai_store->mode = Calloc(n, double);
		Memcpy(ai_store->mode, ai_store->problem->mean_constr, n * sizeof(double));
	} else {
		GMRFLib_ai_add_Qinv_to_ai_store(ai_store);
	}

	/*
	 * ensure that the Qinv is present to be sure.
	 */
	GMRFLib_ai_add_Qinv_to_ai_store(ai_store);

	/*
	 * for internal use 
	 */
	assert(idx < ai_store->problem->sub_graph->n);
	assert(ai_store->stdev);
	x_mean = ai_store->problem->mean_constr[idx];
	x_sd = ai_store->stdev[idx];

	/*
	 * if just the Gaussian is needed, exit here 
	 */
	if (strategy == GMRFLib_AI_STRATEGY_GAUSSIAN) {
		/*
		 * build the density-object 
		 */
		GMRFLib_density_create_normal(density, 0.0, 1.0, x_mean, x_sd, lookup_tables);
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

	if (!ai_store) {
		FIXME("I assume ai_store is !NULL.");
		abort();
	}

	n_points = ai_par->n_points;
	n_points = 2 * (ai_par->n_points / 2) + 1;	       /* ensure its odd */
	GMRFLib_ghq_abscissas(&x_points, n_points);	       /* get the x-points */
	qsort(x_points, (size_t) n_points, sizeof(double), GMRFLib_dcmp_abs);	/* sort them using ABS() */
	x_points[0] = 0.0;

	log_density = Calloc(n_points, double);		       /* values of the log_density */
	cond_mode = Calloc(n, double);
	fixed_mode = Calloc(n, double);
	Memcpy(cond_mode, ai_store->problem->mean_constr, n * sizeof(double));
	Memcpy(fixed_mode, ai_store->problem->mean_constr, n * sizeof(double));

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
	if ((strategy == GMRFLib_AI_STRATEGY_MEANCORRECTED_GAUSSIAN || strategy == GMRFLib_AI_STRATEGY_MEANSKEWCORRECTED_GAUSSIAN)
	    && ai_par->linear_correction == GMRFLib_AI_LINEAR_CORRECTION_FAST) {
		// do nothing
	} else {
		if (!ISZERO(ai_par->cutoff)) {
			for (i = 0; i < n; i++) {
				a = x_sd * derivative[i] / ai_store->stdev[i];
				if (ABS(a) < ai_par->cutoff) {
					fix[i] = fixx[i] = 1;
					derivative[i] = 0.0;
				}
			}
		}
		derivative[idx] = 1.0;

		/*
		 * note that idx is included in subgraph 
		 */
		GMRFLib_graph_comp_subgraph(&subgraph, graph, fixx, &node_map);
	}

	fixx[idx] = 0;					       /* this is how 'fix' and 'fixx' differ */
	fix[idx] = 1;

	optpar->opt_type = GMRFLib_OPTTYPE_NR;		       /* force this option */
	store = Calloc(1, GMRFLib_store_tp);		       /* this can be used ;-) */

	if ((ai_par->linear_correction == GMRFLib_AI_LINEAR_CORRECTION_FAST) && !(ai_store->correction_term)) {
#pragma omp critical (Name_c92a8d89dec4a7d11c6665a243937d76391084de)
		{
			if (!(ai_store->correction_term)) {
				ai_store->correction_term = Calloc(n, double);	/* compute this */
				ai_store->derivative3 = Calloc(n, double);	/* and this */
				ai_store->derivative4 = Calloc(n, double);	/* and this */
				ai_store->correction_idx = Calloc(n, int);	/* and this one */

				/*
				 * the same code splitted 
				 */
				ai_store->nidx = 0;
				// printf("RECOMPUTE derivative3 for thread %d and idx %d\n", omp_get_thread_num(), idx);
				for (ii = 0; ii < ai_store->nd; ii++) {
					double dd;
					i = ai_store->d_idx[ii];
					ai_store->correction_idx[ai_store->nidx++] = i;
					GMRFLib_2order_approx(thread_id, NULL, NULL, NULL, &dd, d[i], fixed_mode[i] + deldif, i,
							      fixed_mode, loglFunc, loglFunc_arg, &(ai_par->step_len), &(ai_par->stencil), NULL);
					ai_store->derivative3[i] = dd;
					ai_store->correction_term[i] = -SQR(ai_store->stdev[i]) * ai_store->derivative3[i];
#if defined(INLA_RESEARCH1)
					if (1) {
						// have to redo this later if this gets serious
						double d3[2];
						double s = 1.0 / (2.0 * deldif);
						GMRFLib_2order_approx(thread_id, NULL, NULL, NULL, &(d3[0]), d[i], fixed_mode[i] - deldif, i,
								      fixed_mode, loglFunc, loglFunc_arg,
								      &(ai_par->step_len), &(ai_par->stencil), NULL);
						GMRFLib_2order_approx(thread_id, NULL, NULL, NULL, &(d3[1]), d[i], fixed_mode[i] + deldif, i,
								      fixed_mode, loglFunc, loglFunc_arg,
								      &(ai_par->step_len), &(ai_par->stencil), NULL);
						ai_store->derivative4[i] = (d3[1] - d3[0]) * s;
					}
#endif

				}
			}
		}
	}

	switch (ai_par->linear_correction) {
	case GMRFLib_AI_LINEAR_CORRECTION_OFF:
		break;					       /* including this; warnings in gcc */

	case GMRFLib_AI_LINEAR_CORRECTION_FAST:
	{
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
			assert(store_Qinv);
			i_idx = store_Qinv->mapping[idx];
			for (j = 0; j < ai_store->nidx; j++) {
				i = ai_store->correction_idx[j];
				j_idx = store_Qinv->mapping[i];
				cov = map_id_ptr(store_Qinv->Qinv[IMIN(i_idx, j_idx)], IMAX(i_idx, j_idx));
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
	}
		break;

	case GMRFLib_AI_LINEAR_CORRECTION_CENTRAL_DIFFERENCE:
	{

		FIXME("option GMRFLib_AI_LINEAR_CORRECTION_CENTRAL_DIFFERENCE is no longer supported.");
		assert(0 == 1);
	}
		break;

	default:
		break;
	}

	if (strategy != GMRFLib_AI_STRATEGY_MEANCORRECTED_GAUSSIAN && strategy != GMRFLib_AI_STRATEGY_MEANSKEWCORRECTED_GAUSSIAN) {
		/*
		 * this is the main loop, where we evaluate different values for x_i. 
		 */

		int debug_me = 0;
		double ql[2] = { 0.0, 0.0 };

		for (k = 0; k < n_points; k++) {
			/*
			 * find first the conditional mode: compute the initial value `cond_mode' 
			 */
			for (ii = 0; ii < subgraph->n; ii++) {
				i = node_map[ii];
				cond_mode[i] = fixed_mode[i] + x_points[k] * x_sd * derivative[i];
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
				FIXME("THIS OPTION IS NO LONGER SUPPORTED");
				assert(0 == 1);
			}

			/*
			 * this is the fast version that take into account that x = x_mode + delta * gradient for the
			 * quadratic term
			 * 
			 * we first initialise the routine computing the linear (ql[0]) and quadratic term (ql[1]), and then we can get the
			 * speedup for successive calls
			 */
			if (k == 0) {
				assert(ISZERO(x_points[k]));
				GMRFLib_ai_log_posterior_restricted(thread_id,
								    NULL,
								    cond_mode, fixed_mode, derivative,
								    0.0, b, c, mean, d, loglFunc,
								    loglFunc_arg, graph, Qfunc, Qfunc_arg, constr, subgraph, ai_store,
								    node_map, ql);
				// GMRFLib_ai_log_posterior(&logdens_ref, cond_mode, b, c, mean, d, loglFunc, loglFunc_arg, graph, Qfunc,
				// Qfunc_arg, NULL);
			}
			GMRFLib_ai_log_posterior_restricted(thread_id,
							    &log_density[k],
							    cond_mode, fixed_mode, derivative,
							    x_points[k] * x_sd, b, c, mean, d, loglFunc,
							    loglFunc_arg, graph, Qfunc, Qfunc_arg, constr, subgraph, ai_store, node_map, ql);
			// GMRFLib_ai_log_posterior(&logdens_new, cond_mode, b, c, mean, d, loglFunc, loglFunc_arg, graph, Qfunc, Qfunc_arg, NULL);

			log_density[k] -= log_dens_cond;
		}
		Free(node_map);
		GMRFLib_graph_free(subgraph);
	}

	GMRFLib_free_store(store);

	if (strategy == GMRFLib_AI_STRATEGY_MEANCORRECTED_GAUSSIAN) {
		GMRFLib_density_create_normal(density, -deriv_log_dens_cond, 1.0, x_mean, x_sd, lookup_tables);
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

#if defined(INLA_RESEARCH1)
		double fourth_order_derivative = 0.0;
		for (jjj = 0; jjj < ai_store->nidx; jjj++) {
			iii = ai_store->correction_idx[jjj];
			fourth_order_derivative += ai_store->derivative4[iii] * gsl_pow_4(derivative[iii]);
		}
		fourth_order_derivative *= gsl_pow_4(x_sd);
		fprintf(stdout, "RESEARCH1: idx= %1d mean %.8f sd %.8f d3 %.8f d4 %.8f\n", idx, x_mean, x_sd, third_order_derivative,
			fourth_order_derivative);
#endif

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

			const int debug = 0;
			if (!(ai_par->improved_simplified_laplace)) {
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
			} else {
				double mom[3] = { 0.0, 1.0, 0.0 };

				mom[0] = -deriv_log_dens_cond;
				mom[2] = GMRFLib_sn_d3_to_skew(third_order_derivative);
				if (ABS(mom[2]) > GMRFLib_SN_SKEWMAX) {
					mom[2] = 0.0;
					fail = 1;
				}
				if (debug) {
					printf("NEW mom: %f %f %f\n", mom[0], mom[1], mom[2]);
				}
				GMRFLib_sn_moments2par(&snp, &mom[0], &mom[1], &mom[2]);
			}
			if (debug && !fail) {
				double mm[3];
				GMRFLib_sn_par2moments(&mm[0], &mm[1], &mm[2], &snp);
				printf("NEW d3 %f\n", third_order_derivative);
				printf("NEW par: %f %f %f\n", snp.xi, snp.omega, snp.alpha);
				printf("NEW mm: %f %f %f\n\n", mm[0], mm[1], mm[2]);
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
		GMRFLib_density_create_sn(density, snp, x_mean, x_sd, lookup_tables);
	} else {
		GMRFLib_density_create(density, GMRFLib_DENSITY_TYPE_SCGAUSSIAN, n_points, x_points, log_density, x_mean, x_sd, lookup_tables);

		if (0) {
			printf("idx %d Gaussian mean %f sd %f internal mean %f sd %f user mean %f sd %f\n",
			       idx, x_mean, x_sd, (*density)->mean, (*density)->stdev, (*density)->user_mean, (*density)->user_stdev);
		}
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

int GMRFLib_ai_update_conditional_mean(int thread_id, GMRFLib_problem_tp * pproblem, double *x, double *mean,
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
		GMRFLib_EWRAP1(GMRFLib_init_problem(thread_id, &pproblem, x, bbb, ccc, mean, graph, Qfunc, Qfunc_args, c));
		if (covariances) {
			*covariances = NULL;
		}
	} else {
		/*
		 * this is extracted from the problem_setup routine compute the new conditional mean using a rank-1 update assuming the old
		 * problem is the same. 
		 */

		GMRFLib_problem_tp **problem = &pproblem;
		int i, j, sub_n, nc, k, kk, inc = 1, n;
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

		assert(qi_at_m_store);
		(*problem)->qi_at_m = Calloc(nc * sub_n, double);
		Memcpy((*problem)->qi_at_m, qi_at_m_store, (nc - 1) * sub_n * sizeof(double));
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
				(*covariances)[i] = (*problem)->qi_at_m[kk + i];
			}

			if ((*problem)->sub_constr->nc > 1) {
				for (i = 0; i < sub_n; i++) {
					if (idx == i) {
						break;
					}
				}
				idx_map = i;
				for (i = 0; i < sub_n; i++) {
					ii = i;
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
		if (GMRFLib_faster_constr) {
			dgemm_special(nc, sub_n, aqat_m, (*problem)->sub_constr->a_matrix, (*problem)->qi_at_m, (*problem)->sub_constr);
		} else {
			dgemm_("N", "N", &nc, &nc, &sub_n, &alpha, (*problem)->sub_constr->a_matrix, &nc, (*problem)->qi_at_m, &sub_n,
			       &beta, aqat_m, &nc, F_ONE, F_ONE);
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

		for (i = 0, k = 0; i < sub_n; i++) {
			for (j = 0; j < nc; j++) {
				tmp_vector[k++] = (*problem)->qi_at_m[i + j * sub_n];
			}
		}

		GMRFLib_EWRAP1(GMRFLib_solveAxb_posdef(tmp_vector, (*problem)->l_aqat_m, tmp_vector, nc, sub_n));

		for (i = 0, k = 0; i < sub_n; i++) {
			for (j = 0; j < nc; j++) {
				(*problem)->constr_m[i + j * sub_n] = tmp_vector[k++];
			}
		}

		Free(tmp_vector);
		nc = constr->nc;
		Memcpy((*problem)->sub_mean_constr, (*problem)->sub_mean, sub_n * sizeof(double));
		t_vector = Calloc(nc, double);

		GMRFLib_EWRAP1(GMRFLib_eval_constr(t_vector, NULL, (*problem)->sub_mean, (*problem)->sub_constr, (*problem)->sub_graph));

		alpha = -1.0;
		beta = 1.0;				       /* mean_constr = mean - cond_m*t_vector */
		dgemv_("N", &sub_n, &nc, &alpha, (*problem)->constr_m, &sub_n, t_vector, &inc, &beta, (*problem)->sub_mean_constr, &inc, F_ONE);
		Free(t_vector);

		for (i = 0; i < sub_n; i++) {
			(*problem)->mean_constr[i] = (*problem)->sub_mean_constr[i];
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

	// this one is FIXED by design and such that A[IDX(i,j,n)] = A_ij, i=0...n-1, j = 0..k-1 for n x k matrix A 
#define IDX(i, j, n) ((i) + (j)*(n))

	int k, n, nc, ncc, one = 1;
	double *c = NULL, *v = NULL, *w = NULL, *z = NULL, alpha = 0.0, beta = 0.0, b22 = 0.0, *constr_m_new = NULL, *t_vec =
	    NULL, *tmp_m = NULL, val;

	GMRFLib_ENTER_ROUTINE;

	n = problem->n;
	assert(n == problem->sub_graph->n);
	nc = (problem->sub_constr ? problem->sub_constr->nc : 0);
	ncc = nc + 1;

	/*
	 * setup workspace for small-mem's for the hole routine here. 
	 */
	Calloc_init(n + ncc + (nc ? n + nc + nc + ISQR(nc) : 0), 2 + (nc ? 4 : 0));
	c = Calloc_get(n);
	t_vec = Calloc_get(ncc);
	if (nc) {
		z = Calloc_get(n);
		v = Calloc_get(nc);
		w = Calloc_get(nc);
		tmp_m = Calloc_get(ISQR(nc));
	}

	c[idx] = 1.0;
	GMRFLib_solve_llt_sparse_matrix_special(c, &(problem->sub_sm_fact), problem->sub_graph, idx);

	if (covariances) {
		*covariances = Calloc(n, double);
		Memcpy(*covariances, c, n * sizeof(double));

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

		assert(w);
		assert(v);
		assert(z);

		if (!problem->inv_aqat_m) {
			double *m;
			m = Calloc(ISQR(nc), double);
			alpha = 1.0;
			beta = 0.0;
			if (GMRFLib_faster_constr) {
				dgemm_special(nc, n, m, problem->sub_constr->a_matrix, problem->qi_at_m, problem->sub_constr);
			} else {
				dgemm_("N", "N", &nc, &nc, &n, &alpha, problem->sub_constr->a_matrix, &nc, problem->qi_at_m,
				       &n, &beta, m, &nc, F_ONE, F_ONE);
			}
			GMRFLib_comp_posdef_inverse(m, nc);
#pragma omp critical (Name_c13803f68f931b7ce671391d062203d502409e15)
			{
				if (!problem->inv_aqat_m) {
					problem->inv_aqat_m = m;
				}
			}
		}

		/*
		 * v = A c.  w = inv(AQ^-1A^T) v.  z = (Q^-1 A^T) w. 
		 */
		alpha = 1.0;
		beta = 0.0;
		if (GMRFLib_faster_constr) {
			dgemv_special(v, c, problem->sub_constr);
		} else {
			dgemv_("N", &nc, &n, &alpha, problem->sub_constr->a_matrix, &nc, c, &one, &beta, v, &one, F_ONE);
		}
		dgemv_("N", &nc, &nc, &alpha, problem->inv_aqat_m, &nc, v, &one, &beta, w, &one, F_ONE);
		dgemv_("N", &n, &nc, &alpha, problem->qi_at_m, &n, w, &one, &beta, z, &one, F_ONE);

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
		dgemm_("N", "N", &n, &nc, &nc, &alpha, problem->constr_m, &n, tmp_m, &nc, &beta, constr_m_new, &n, F_ONE, F_ONE);

		for (k = 0; k < nc; k++) {
			val = -b22 * w[k];
			daxpy_(&n, &val, c, &one, &constr_m_new[IDX(0, k, n)], &one);
		}

		k = (ncc - 1) * n;
		double *c_tmp = constr_m_new + k;
#pragma GCC ivdep
		for (int i = 0; i < n; i++) {
			c_tmp[i] = b22 * (c[i] - z[i]);
		}

		// No longer needed: GMRFLib_eval_constr(t_vec, NULL, problem->sub_mean, problem->sub_constr, problem->sub_graph);
		assert(problem->sub_constr_value);
		Memcpy(t_vec, problem->sub_constr_value, nc * sizeof(double));
	} else {
		double cinv = 1.0 / c[idx];

		// for(i=0; i<n; i++) constr_m_new[k + i] = cinv*c[i];
		daxpy_(&n, &cinv, c, &one, &constr_m_new[nc * n], &one);
	}
	t_vec[ncc - 1] = problem->sub_mean[idx] - evalue;

	alpha = -1.0;
	beta = 1.0;					       /* mean_constr = mean - cond_m*t_vec */
	Memcpy(cond_mean, problem->sub_mean, n * sizeof(double));
	dgemv_("N", &n, &ncc, &alpha, constr_m_new, &n, t_vec, &one, &beta, cond_mean, &one, F_ONE);

	Calloc_free();
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
		Free(ai_store->derivative4);
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

	assert(z);
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

int GMRFLib_init_GMRF_approximation_store__intern(int thread_id,
						  GMRFLib_problem_tp ** problem, double *x, double *b, double *c, double *mean,
						  double *d, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg,
						  GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg,
						  GMRFLib_constr_tp * constr, GMRFLib_optimize_param_tp * optpar,
						  GMRFLib_blockupdate_param_tp * blockupdate_par, GMRFLib_store_tp * store,
						  double *aa, double *bb, double *cc, int gaussian_data, double cmin, int b_strategy, int nested,
						  GMRFLib_preopt_tp * preopt)
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

	int i, free_b = 0, free_c = 0, free_mean = 0, free_d = 0, free_blockpar = 0, free_aa = 0, free_bb = 0, free_cc =
	    0, n, *idxs = NULL, nidx = 0;
	int Npred = (GMRFLib_inla_mode == GMRFLib_MODE_TWOSTAGE_PART1 || GMRFLib_inla_mode == GMRFLib_MODE_COMPACT ? preopt->Npred : graph->n);
	double *mode = NULL;

#define FREE_ALL if (1) { if (free_b) Free(b); if (free_c) Free(c); if (free_d) Free(d); \
		if (free_mean) Free(mean); if (free_blockpar) Free(blockupdate_par); if (free_aa) Free(aa); if (free_bb) Free(bb); \
		if (free_cc) Free(cc); Free(idxs); }

	GMRFLib_ENTER_ROUTINE;

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

	/*
	 * the NEW implementation which do optimisation and GMRF_approximation in the same operation. this is tailored for INLA of'course, and only ment
	 * for that. 
	 */

	if (optpar && optpar->fp)
		fprintf(optpar->fp, "\n[%1d] Computing GMRF approximation\n------------------------------\n", omp_get_thread_num());
	nidx = 0;
	idxs = Calloc(Npred, int);
	for (i = 0; i < Npred; i++) {
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

	Memcpy(mode_initial, mode, n * sizeof(double));	       /* store the starting value */

	double *bcoof = NULL, *ccoof = NULL, *linear_predictor = NULL;
	int free_linear_predictor = 0;

	bcoof = Calloc(Npred, double);
	ccoof = Calloc(Npred, double);

	for (iter = 0; iter < itmax; iter++) {

		if (0)
			for (i = 0; i < n; i++) {
				printf("i mode %d %f \n", i, mode[i]);
			}

		Memset(aa, 0, Npred * sizeof(double));
		Memset(bb, 0, Npred * sizeof(double));
		Memset(cc, 0, Npred * sizeof(double));

		if (GMRFLib_inla_mode == GMRFLib_MODE_TWOSTAGE_PART1 || GMRFLib_inla_mode == GMRFLib_MODE_COMPACT) {
			if (!free_linear_predictor) {
				linear_predictor = Calloc(Npred, double);
				free_linear_predictor = 1;
			}
			GMRFLib_preopt_predictor(linear_predictor, mode, preopt);

			if (0)
				for (i = 0; i < preopt->Npred; i++) {
					printf("i predictor %d %f\n", i, linear_predictor[i]);
				}

		} else {
			linear_predictor = mode;
		}

		cc_is_negative = 0;
		Memset(bcoof, 0, Npred * sizeof(double));
		Memset(ccoof, 0, Npred * sizeof(double));

#define CODE_BLOCK							\
		for (int i_ = 0; i_ < nidx; i_++) {			\
			int idx = idxs[i_];				\
			GMRFLib_2order_approx(thread_id, &(aa[idx]), &(bcoof[idx]), &(ccoof[idx]), NULL, d[idx], \
					      linear_predictor[idx], idx, mode, loglFunc, loglFunc_arg, \
					      &(optpar->step_len), &(optpar->stencil), &cmin); \
		}

		RUN_CODE_BLOCK(GMRFLib_openmp->max_threads_inner, 0, 0);
#undef CODE_BLOCK

		for (i = 0; i < nidx; i++) {
			int idx;
			idx = idxs[i];
			cc_is_negative = (cc_is_negative || ccoof[idx] < 0.0);	/* this line IS OK! also for multithread.. */
			if (cc_positive) {
				if (ccoof[idx] == cmin) {
					// then cmin is in effect
					if (b_strategy == INLA_B_STRATEGY_SKIP) {
						bcoof[idx] = 0.0;
					} else if (b_strategy == INLA_B_STRATEGY_KEEP) {
						// do nothing
					} else {
						assert(0 == 1);
					}
				}
				bb[idx] += bcoof[idx];
				cc[idx] += ccoof[idx];
			} else {
				// 
				// this is not in use...
				// 
				if (ccoof[idx] > 0.0) {
					bb[idx] += bcoof[idx];
					cc[idx] += ccoof[idx];
				} else {
					bb[idx] += bcoof[idx];
					// bb[idx] += cc_factor*bcoof[idx]; /* what to use?? if any...*/
					cc[idx] += cc_factor * ccoof[idx];
				}
			}
		}

		double *bb_use = NULL, *cc_use = NULL;

		if (0)
			for (i = 0; i < Npred; i++) {
				printf("i bb cc %d %f %f\n", i, bb[i], cc[i]);
			}

		if (preopt) {
			GMRFLib_preopt_update(thread_id, preopt, bb, cc);
			bb_use = preopt->total_b[thread_id];
#pragma GCC ivdep
			for (i = 0; i < n; i++) {
				bb_use[i] += b[i];
			}
			cc_use = c;			       /* that what is there from before */
			if (0)
				for (i = 0; i < n; i++) {
					printf("i bb_preopt cc_preopt %d %f\n", i, preopt->total_b[thread_id][i]);
				}
		} else {
			assert(Npred == n);

#pragma GCC ivdep
			for (i = 0; i < n; i++) {
				bb[i] += b[i];
				cc[i] += c[i];
			}
			bb_use = bb;
			cc_use = cc;
			if (0)
				for (i = 0; i < n; i++) {
					printf("i bb_use cc_use %d %f %f\n", i, bb_use[i], cc_use[i]);
				}
		}

		if (!cc_positive) {
			cc_factor = DMIN(1.0, cc_factor * cc_factor_mult);
		}

		/*
		 * I thought this was quicker without store, as there is just reuse and no copy... but not.  I free lproblem below and set it to NULL, so
		 * it will always be lproblem = NULL 
		 */

		if (!lproblem) {
			if (preopt)
				assert((void *) preopt == (void *) Qfunc_arg);
			// GMRFLib_printf_Qfunc2(stdout, graph, Qfunc, Qfunc_arg);

			int ret;
			ret = GMRFLib_init_problem_store(thread_id, &lproblem, x, bb_use, cc_use, mean, graph, Qfunc, Qfunc_arg, constr, store);
			if (ret != GMRFLib_SUCCESS) {
				catch_error = 1;
			}
		} else {
			GMRFLib_init_problem_store(thread_id, &lproblem, x, bb_use, cc_use, mean, graph, Qfunc, Qfunc_arg, constr, store);
		}

		if (0)
			for (i = 0; i < n; i++) {
				printf("AAA iter i mean %d %d %f\n", iter, i, lproblem->mean_constr[i]);
			}

		if (catch_error) {
			lproblem = NULL;
			break;
		}

		int flag_cycle_behaviour = 0;
		double f;
		if (gaussian_data) {
			f = 1.0;
		} else {
			f = DMIN(1.0, (iter + 1.0) * optpar->nr_step_factor);
		}

		double err = 0.0;
#pragma GCC ivdep
		for (i = 0; i < n; i++) {
			err += SQR((lproblem)->mean_constr[i] - mode[i]);
		}
		err = sqrt(err / n);

#pragma GCC ivdep
		for (i = 0; i < n; i++) {
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
			fprintf(optpar->fp, "[%1d] iteration %d error %.4f\n", thread_id, iter, err);

		if (gaussian_data) {
			/*
			 * I need to update 'aa' as this is not evaluated in the mode! The sum of the a's are/might-be used later
			 */

			if (GMRFLib_inla_mode == GMRFLib_MODE_TWOSTAGE_PART1 || GMRFLib_inla_mode == GMRFLib_MODE_COMPACT) {
				GMRFLib_preopt_predictor(linear_predictor, mode, preopt);
			} else {
				linear_predictor = mode;
			}
			for (i = 0; i < nidx; i++) {
				int idx = idxs[i];
				GMRFLib_2order_approx(thread_id,
						      &(aa[idx]), NULL, NULL, NULL, d[idx], linear_predictor[idx], idx, mode, loglFunc,
						      loglFunc_arg, &(optpar->step_len), &(optpar->stencil), NULL);
			}
		}
		// about comparing with err_previous, then we're already in the good regime, and another iteration will not give anything new.
		// this assumes err_i = m * (err_{i-1})^2
		double m = err / SQR(err_previous);
		int almost_there = ((iter > 0) && (f >= 1.0) && (err > optpar->abserr_step) &&
				    (m <= 1.0) && (m * SQR(err) < 0.1 * optpar->abserr_step));

		if (gaussian_data || err < optpar->abserr_step || almost_there || flag_cycle_behaviour) {
			/*
			 * we're done!  unless we have negative elements on the diagonal...
			 */
			break;

			// NOT IN USE
			if (cc_is_negative && !cc_positive && (cc_factor < 1.0)) {
				/*
				 * do nothing 
				 */
			} else if (cc_is_negative && cc_positive) {
				cc_positive = 0;
			} else {
				break;
			}
		}

		if (gsl_isnan(err))
			break;

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
			for (i = 0; i < graph->n; i++)
				printf(" %.4f", (*problem)->sub_mean_constr[i]);
			printf("\n");
		}
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
						c_new[i] = lambda * Qfunc(thread_id, i, i, NULL, Qfunc_arg) + lambda * c[i] + lambda;
						if (x && (ISNAN(x[i]) || ISINF(x[i]))) {
							x[i] = mode[i];
						}
						if (x && (ISNAN(x[i]) || ISINF(x[i]))) {
							x[i] = 0.0;
						}
					}
					retval = GMRFLib_init_GMRF_approximation_store__intern(thread_id, problem, x, b, c_new, mean, d,
											       loglFunc, loglFunc_arg,
											       graph, Qfunc, Qfunc_arg, constr,
											       &new_optpar, blockupdate_par, store,
											       aa, bb, cc, gaussian_data, cmin, b_strategy, 2,
											       preopt);
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
	Memcpy(iz_local, iz, nhyper * sizeof(int));
	Memcpy(izz_local, izz, nhyper * sizeof(int));

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

int GMRFLib_ai_INLA(GMRFLib_density_tp *** density,
		    GMRFLib_density_tp *** density_transform, GMRFLib_transform_array_func_tp ** tfunc,
		    GMRFLib_density_tp *** density_hyper,
		    GMRFLib_ai_cpo_tp ** cpo, GMRFLib_ai_po_tp ** po, GMRFLib_ai_dic_tp * dic,
		    GMRFLib_ai_marginal_likelihood_tp * marginal_likelihood,
		    char *compute, double ***hyperparam, int nhyper,
		    GMRFLib_ai_log_extra_tp * log_extra, void *log_extra_arg,
		    double *x, double *b, double *c, double *mean,
		    GMRFLib_bfunc_tp ** bfunc, double *d,
		    GMRFLib_logl_tp * loglFunc, void *loglFunc_arg,
		    GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg,
		    GMRFLib_constr_tp * constr, GMRFLib_ai_param_tp * ai_par, GMRFLib_ai_store_tp * ai_store,
		    int nlin, GMRFLib_lc_tp ** Alin, GMRFLib_density_tp *** dlin, GMRFLib_ai_misc_output_tp * misc_output,
		    GMRFLib_preopt_tp * preopt, GMRFLib_preopt_res_tp * rpreopt)
{
	assert(omp_get_thread_num() == 0);

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

#define COMPUTE_LINDENS(_store, _lookup_tables)				\
	if (nlin) {							\
		double *_improved_mean = Calloc(graph->n, double);	\
		for(int _i = 0; _i<graph->n; _i++) {			\
			if (dens && dens[_i] && dens[_i][dens_count]){		\
				_improved_mean[_i] = dens[_i][dens_count]->user_mean; \
			} else {					\
				_improved_mean[_i] = ai_store->problem->mean_constr[_i]; \
			}						\
		}							\
		GMRFLib_ai_compute_lincomb(&(lin_dens[dens_count]), (lin_cross ? &(lin_cross[dens_count]) : NULL), nlin, Alin, _store, _improved_mean, _lookup_tables); \
		Free(_improved_mean);					\
	}

#define COMPUTE_LINDENS_LOCAL(_store, _lookup_tables)			\
	if (nlin) {							\
		double *_improved_mean = Calloc(graph->n, double);	\
		for(int _i = 0; _i<graph->n; _i++) {			\
			if (dens_local && dens_local[_i]){		\
				_improved_mean[_i] = dens_local[_i]->user_mean; \
			} else {					\
				_improved_mean[_i] = _store->problem->mean_constr[_i]; \
			}						\
		}							\
		GMRFLib_ai_compute_lincomb(&lin_dens_local, (lin_cross ? &lin_cross_local : NULL), nlin, Alin, _store, _improved_mean, _lookup_tables); \
		Free(_improved_mean);					\
	}

#define ADD_CONFIG(_store, _theta, _log_posterior, _log_posterior_orig)	\
	if (1) {							\
		double *_improved_mean = Calloc(graph->n, double);	\
		double *_skewness = Calloc(graph->n, double);		\
		for(int _i = 0; _i<graph->n; _i++) {			\
			_skewness[_i] = NAN;				\
			if (dens && dens[_i] && dens[_i][dens_count]){		\
				_improved_mean[_i] = dens[_i][dens_count]->user_mean; \
				_skewness[_i] = dens[_i][dens_count]->skewness;	\
			} else {					\
				_improved_mean[_i] = (_store)->problem->mean_constr[_i]; \
			}						\
		}							\
		assert(_store);						\
		GMRFLib_ai_store_config(thread_id, misc_output, nhyper, _theta, _log_posterior, _log_posterior_orig, _improved_mean, _skewness, (_store)->problem, Qfunc, Qfunc_arg, c, dens_count); \
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
		Memset(&(izs[old_dens_max]), 0, (num_) * sizeof(double *)); \
		for (kk_ = 0; kk_ < compute_n; kk_++) {			\
			ii_ = compute_idx[kk_];				\
			if (dens && dens[ii_]){					\
				dens[ii_] = Realloc(dens[ii_], dens_max, GMRFLib_density_tp *); \
				for(jj_ = old_dens_max; jj_ < dens_max; jj_++) \
					dens[ii_][jj_] = NULL;		\
			}						\
			if (dens_transform && dens_transform[ii_]){			\
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
				if (d[jj_]){				\
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
					deviance_theta[jj_] = Realloc(deviance_theta[jj_], dens_max, double *); \
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
			failure_theta[ii][dens_count] = GMRFLib_ai_cpopit_integrate(thread_id, &cpo_theta[ii][dens_count], \
										    &pit_theta[ii][dens_count], ii, cpodens, \
										    (ai_par->cpo_manual ? 1.0 : d[ii]), loglFunc, loglFunc_arg, xx_mode); \
			if (cpodens && GMRFLib_getbit(cpodens->flags, DENSITY_FLAGS_FAILURE)) { \
				failure_theta[ii][dens_count] = 1.0;	\
			}						\
		}							\
		if (dic) {						\
			deviance_theta[ii][dens_count] = GMRFLib_ai_dic_integrate(thread_id, ii, dens[ii][dens_count], \
										  (ai_par->cpo_manual ? 1.0 : d[ii]), loglFunc, loglFunc_arg, xx_mode); \
		}							\
	}

#define COMPUTE_PO							\
	if (d[ii]) {							\
		if (po) {						\
			GMRFLib_ai_po_integrate(thread_id, &po_theta[ii][dens_count], &po2_theta[ii][dens_count], &po3_theta[ii][dens_count], \
						ii, dens[ii][dens_count], d[ii], loglFunc, loglFunc_arg, xx_mode); \
		}							\
	}

#define COMPUTE_CPO_AND_DIC_LOCAL					\
	if (d[ii] || ai_par->cpo_manual) {				\
		if (cpo || ai_par->cpo_manual) {			\
			failure_theta_local[ii] +=			\
				GMRFLib_ai_cpopit_integrate(thread_id, &cpo_theta_local[ii], &pit_theta_local[ii],	\
							    ii, cpodens, (ai_par->cpo_manual ? 1.0 : d[ii]), loglFunc, loglFunc_arg, xx_mode); \
			if (cpodens && GMRFLib_getbit(cpodens->flags, DENSITY_FLAGS_FAILURE)) { \
				failure_theta_local[ii] += 1.0;		\
			}						\
		}							\
		if (dic) {						\
			deviance_theta_local[ii] =			\
				GMRFLib_ai_dic_integrate(thread_id, ii, dens_local[ii], (ai_par->cpo_manual ? 1.0 : d[ii]), loglFunc, loglFunc_arg, xx_mode); \
		}							\
	}

#define COMPUTE_PO_LOCAL						\
	if (d[ii]) {							\
		if (po) {						\
			GMRFLib_ai_po_integrate(thread_id, &po_theta_local[ii], &po2_theta_local[ii], &po3_theta_local[ii], \
						ii, dens_local[ii], d[ii], loglFunc, loglFunc_arg, xx_mode); \
		}							\
	}

/* 
   since all threads compute the same quantity, this is it well defined
*/
#define COMPUTE        COMPUTE_CPO_AND_DIC; COMPUTE_PO;
#define COMPUTE2       COMPUTE_CPO_AND_DIC; COMPUTE_PO;
#define COMPUTE_LOCAL  COMPUTE_CPO_AND_DIC_LOCAL; COMPUTE_PO_LOCAL;

	int i, j, k, *k_max = NULL, *k_min = NULL, *k_maxx = NULL, *k_minn = NULL, ierr, *iz = NULL, *izz = NULL, *len =
	    NULL, *iz_axes = NULL, len_length, free_ai_par = 0, config_count = 0, free_compute = 0, dens_count =
	    0, dens_max, hyper_len = 0, hyper_count = 0, *compute_idx = NULL, compute_n = 0, tmax, need_Qinv = 0;

	double *hessian = NULL, *theta = NULL, *theta_mode = NULL, *x_mode = NULL, log_dens_mode = 0, log_dens, *z = NULL, **izs =
	    NULL, *stdev_corr_pos = NULL, *stdev_corr_neg = NULL, f, w, w_origo, tref, tu, *weights = NULL, *adj_weights =
	    NULL, *hyper_z = NULL, *hyper_ldens = NULL, **userfunc_values = NULL, *inverse_hessian = NULL, *timer,
	    **cpo_theta = NULL, **po_theta = NULL, **po2_theta = NULL, **po3_theta = NULL, **pit_theta = NULL, ***deviance_theta =
	    NULL, **failure_theta = NULL;
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

	assert(GMRFLib_inla_mode != GMRFLib_MODE_COMPACT);
	if (GMRFLib_inla_mode == GMRFLib_MODE_CLASSIC) {
		assert(!preopt);
		assert(!rpreopt);
	}
	if (GMRFLib_inla_mode == GMRFLib_MODE_TWOSTAGE_PART1) {
		assert(preopt);
		assert(rpreopt);
	}
	if (GMRFLib_inla_mode == GMRFLib_MODE_TWOSTAGE_PART2) {
		assert(!preopt);
		assert(rpreopt);
	}

	if (!(mean == NULL)) {
		FIXME("\n\n\n\nGMRFLib_INLA() I think the vb assumes mean=NULL, please check.\n");
		abort();
	}

	tmax = GMRFLib_MAX_THREADS();
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

	nhyper = IMAX(0, nhyper);
	dens_max = 1;

	// do initial setup that we do not need with 'preopt'
	weights = Calloc(dens_max, double);
	izs = Calloc(dens_max, double *);

	if (!compute) {
		free_compute = 1;
		compute = Calloc(graph->n, char);
	}
	compute_idx = Calloc(graph->n, int);
	compute_n = 0;
	for (i = 0; i < graph->n; i++) {
		if (compute[i]) {
			compute_idx[compute_n++] = i;
		}
	}
	x_mode = Calloc(graph->n, double);
	map_strd_init_hint(&hash_table, dens_max);
	hash_table.alwaysdefault = 0;

	if (GMRFLib_inla_mode == GMRFLib_MODE_CLASSIC || GMRFLib_inla_mode == GMRFLib_MODE_TWOSTAGE_PART2) {
		dens = Calloc(graph->n, GMRFLib_density_tp **);
		dens_transform = Calloc(graph->n, GMRFLib_density_tp **);
		weights = Calloc(dens_max, double);
		izs = Calloc(dens_max, double *);

		if (density && cpo) {
			(*cpo) = Calloc(1, GMRFLib_ai_cpo_tp);
			(*cpo)->n = graph->n;
			(*cpo)->value = Calloc(graph->n, double *);
			(*cpo)->pit_value = Calloc(graph->n, double *);
			(*cpo)->failure = Calloc(graph->n, double *);
		} else {
			cpo = NULL;
		}
		if (density && po) {
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
		 * only one of the marginal are computed with cpo_manual is TRUE. The code depends on the this assumption I think. 
		 */
		if (ai_par->cpo_manual) {
			GMRFLib_ASSERT(compute_n > 0, GMRFLib_ESNH);
			for (i = 0; i < compute_n; i++) {
				GMRFLib_ASSERT(d[compute_idx[i]] == 0.0, GMRFLib_ESNH);
			}
			if (dic) {			       /* meaningless to compute the DIC in this case */
				dic = NULL;
			}
		}

		need_Qinv = (GMRFLib_inla_mode == GMRFLib_MODE_CLASSIC || GMRFLib_inla_mode == GMRFLib_MODE_TWOSTAGE_PART2)
		    && (compute_n || ai_par->compute_nparam_eff);

		for (i = 0; i < compute_n; i++) {
			j = compute_idx[i];
			dens[j] = Calloc(dens_max, GMRFLib_density_tp *);	/* storage for the marginals */
			if (tfunc && tfunc[j]) {
				dens_transform[j] = Calloc(dens_max, GMRFLib_density_tp *);
			}
		}

		if (cpo) {
			cpo_theta = Calloc(graph->n, double *);	/* cpo-value conditioned on theta */
			pit_theta = Calloc(graph->n, double *);	/* pit-value conditioned on theta */
			failure_theta = Calloc(graph->n, double *);	/* failure indicator on theta */
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
			po_theta = Calloc(graph->n, double *); /* po-value conditioned on theta */
			po2_theta = Calloc(graph->n, double *);	/* po-value conditioned on theta */
			po3_theta = Calloc(graph->n, double *);	/* po-value conditioned on theta */
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
			deviance_theta = Calloc(graph->n, double **);	/* mean of deviance conditioned on theta */
			for (i = 0; i < compute_n; i++) {
				j = compute_idx[i];
				if (d[j]) {
					deviance_theta[j] = Calloc(dens_max, double *);
				}
			}
		}
	}						       /* end of: if (!preopt) */

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
		 * the first step is to locate the mode of \pi(\theta | y). here we use the opt-optimiser routine.  NOTE that this
		 * '_setup' ensure that ai_store is changed for each call to _opt_f. this is a bit dirty programming, but there is no
		 * good way to get around it for the moment.
		 */
		GMRFLib_opt_setup(hyperparam, nhyper, log_extra, log_extra_arg, compute, x, b, c, mean, bfunc, d, loglFunc,
				  loglFunc_arg, graph, Qfunc, Qfunc_arg, constr, ai_par, ai_store, preopt);
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
			{
				int fd_save = ai_par->gradient_forward_finite_difference;
				if (ai_par->optimise_smart) {
					ai_par->gradient_forward_finite_difference = GMRFLib_TRUE;
					if (ai_par->fp_log) {
						fprintf(ai_par->fp_log, "Smart optimise part I: estimate gradient using forward differences\n");
					}
				}

				GMRFLib_gsl_optimize(ai_par);

				if (ai_par->optimise_smart) {
					ai_par->restart = IMAX(1, ai_par->restart);
					ai_par->gradient_forward_finite_difference = GMRFLib_FALSE;
					if (ai_par->fp_log) {
						fprintf(ai_par->fp_log, "Smart optimise part II: estimate gradient using central differences\n");
						fprintf(ai_par->fp_log, "Smart optimise part II: restart optimiser\n");
					}
				}

				if (ai_par->restart) {
					for (k = 0; k < IMAX(0, ai_par->restart); k++) {
						GMRFLib_gsl_optimize(ai_par);	/* restart */
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
			for (i = 0; i < nhyper; i++) {
				theta_mode[i] = hyperparam[i][0][0];
			}
			if (ai_par->fp_log) {
				fprintf(ai_par->fp_log, "Using known modal configuration = [");
				for (i = 0; i < nhyper; i++) {
					fprintf(ai_par->fp_log, " %6.3f", theta_mode[i]);
				}
				fprintf(ai_par->fp_log, " ]\n");
			}
			// this is not needed as we do that below. I am not quite sure if we need this in general, but...
			if (GMRFLib_inla_mode != GMRFLib_MODE_TWOSTAGE_PART2) {
				int thread_id = 0;
				assert(omp_get_thread_num() == 0);
				GMRFLib_opt_f(thread_id, theta_mode, &log_dens_mode, &ierr, NULL, NULL);
				log_dens_mode *= -1.0;
				if (ai_par->fp_log) {
					fprintf(ai_par->fp_log, "Compute mode: %10.3f\n", log_dens_mode);
				}
			}
		}

		SET_THETA_MODE;
		GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_HESSIAN, (void *) &nhyper, NULL);

		if (GMRFLib_inla_mode == GMRFLib_MODE_CLASSIC || GMRFLib_inla_mode == GMRFLib_MODE_TWOSTAGE_PART1) {
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

			// SMART MODE: we try to be smart. do a prerun using forward differences. if its ok, keep it.
			if (ai_par->optimise_smart) {
				ai_par->hessian_forward_finite_difference = GMRFLib_TRUE;
				smart_success = 1;
				if (ai_par->fp_log) {
					fprintf(ai_par->fp_log, "Smart optimise part III: estimate Hessian using forward differences\n");
				}
				while (GMRFLib_opt_estimate_hessian(hessian, theta_mode, &log_dens_mode, stupid_mode_iter) != GMRFLib_SUCCESS) {
					smart_success = 0;
					if (!stupid_mode_iter) {
						if (ai_par->fp_log)
							fprintf(ai_par->fp_log,
								"Mode not sufficient accurate; switch to a stupid local search strategy.\n");
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
					smart_success = 1;
				}
				ai_par->hessian_forward_finite_difference = fd_save;

				if (smart_success) {
					// check if the hessian is valid. if its ok, we accept, otherwise, we retry with central differences
					double *chol_tmp = NULL;
					int ecode = 99, ret_ecode;

					ret_ecode = GMRFLib_comp_chol_general(&chol_tmp, hessian, nhyper, NULL, ecode);
					Free(chol_tmp);
					if (ret_ecode == ecode) {
						// we failed, at least one eigenvalue is negative...
						smart_success = 0;
					}
				}

				/*
				 * do this again to get the ai_store set correctly.
				 */
				SET_THETA_MODE;
				if (x_mode) {
					Memcpy(x_mode, ai_store->mode, graph->n * sizeof(double));
				}

				if (stupid_mode_iter) {
					// FIXME("------------> do one function call");
					for (i = 0; i < nhyper; i++) {
						theta_mode[i] = hyperparam[i][0][0];
					}
					int thread_id = 0;
					assert(omp_get_thread_num() == 0);
					GMRFLib_opt_f(thread_id, theta_mode, &log_dens_mode, &ierr, NULL, NULL);
					log_dens_mode *= -1.0;
					SET_THETA_MODE;
					if (x_mode) {
						Memcpy(x_mode, ai_store->mode, graph->n * sizeof(double));
					}
				}

				if (ai_par->fp_log) {
					if (ai_par->optimise_smart) {
						if (smart_success) {
							fprintf(ai_par->fp_log, "Smart optimise part III: Hessian seems fine, keep it\n");
						} else {
							fprintf(ai_par->fp_log, "Smart optimise part III: detected trouble with the Hessian...\n");
							fprintf(ai_par->fp_log, "Smart optimise part III: try a restart before trying again.\n");
						}
					}
				}

				if (!smart_success) {
					// we'll try to compute the Hessian again, but before that, lets restart the optimizer
					GMRFLib_gsl_optimize(ai_par);	/* restart */
					GMRFLib_gsl_get_results(theta_mode, &log_dens_mode);
				}
			}

			stupid_mode_iter = 0;		       /* reset it */
			if (!(ai_par->optimise_smart) || !smart_success) {

				if (ai_par->optimise_smart) {
					ai_par->hessian_forward_finite_difference = GMRFLib_FALSE;
					if (ai_par->fp_log) {
						fprintf(ai_par->fp_log, "Smart optimise part IV: re-estimate Hessian using central differences\n");
					}
				}

				while (GMRFLib_opt_estimate_hessian(hessian, theta_mode, &log_dens_mode, stupid_mode_iter) != GMRFLib_SUCCESS) {
					if (!stupid_mode_iter) {
						if (ai_par->fp_log)
							fprintf(ai_par->fp_log,
								"Mode not sufficient accurate; switch to a stupid local search strategy.\n");
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

			ai_par->hessian_forward_finite_difference = fd_save;

			/*
			 * do this again to get the ai_store set correctly.
			 */
			SET_THETA_MODE;
			if (x_mode) {
				Memcpy(x_mode, ai_store->mode, graph->n * sizeof(double));
			}

			if (stupid_mode_iter) {
				// FIXME("------------> do one function call");
				for (i = 0; i < nhyper; i++) {
					theta_mode[i] = hyperparam[i][0][0];
				}
				int thread_id = 0;
				assert(omp_get_thread_num() == 0);
				GMRFLib_opt_f(thread_id, theta_mode, &log_dens_mode, &ierr, NULL, NULL);
				log_dens_mode *= -1.0;
				SET_THETA_MODE;
				if (x_mode) {
					Memcpy(x_mode, ai_store->mode, graph->n * sizeof(double));
				}
			}

			if (ai_par->fp_log) {
				for (i = 0; i < nhyper; i++) {
					for (j = 0; j < nhyper; j++) {
						fprintf(ai_par->fp_log, " %10.3f", hessian[i + j * nhyper]);
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
				GMRFLib_printf_gsl_matrix(ai_par->fp_log, eigen_vectors, "\t%7.3f");
				fprintf(ai_par->fp_log, "Eigenvalues of the Hessian\n");
				gsl_vector_fprintf(ai_par->fp_log, eigen_values, "\t%12.3f");
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
				min_pos_eigenvalue = 1.0;      /* if all are negative, zero included */
			}
			int a_change = 0, all_negative = 1;

			for (i = 0; i < nhyper; i++) {
				double eigv = gsl_vector_get(eigen_values, (unsigned int) i);

				all_negative = (all_negative && (eigv <= 0.0 || ISZERO(eigv)));
				if (eigv < 0.0) {
					fprintf(stderr, "\n");
					fprintf(stderr, "\t*** WARNING *** Eigenvalue %1d of the Hessian is %.6g < 0\n", i, eigv);
					fprintf(stderr, "\t*** WARNING *** This have consequence for the accurancy of the hyperpar\n");
					fprintf(stderr, "\t*** WARNING *** Continue with a diagonal Hessian.\n");
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
						"\n\t*** WARNING *** R-inla: All eigenvalues of the Hessian are negative. Move on with Hessian = Identity\n\n");
					Memset(hessian, 0, ISQR(nhyper) * sizeof(double));
					for (i = 0; i < nhyper; i++)
						hessian[i + i * nhyper] = 1.0;
				} else {
					// I have changed my mind. It is better to knock of all off-diagonal terms and just use the
					// diagonal, we do not control anything about negative eigenvalue(s).
					if (1) {
						// new. revert back to a diagonal hessian
						for (i = 0; i < nhyper; i++) {
							hessian[i + i * nhyper] = DMAX(DBL_EPSILON, hessian[i + i * nhyper]);
							for (j = i + 1; j < nhyper; j++) {
								hessian[i + j * nhyper] = hessian[j + i * nhyper] = 0.0;
							}
						}
						// need the new eigenvalues/vectors for futher calculations. its easy, we just compute
						// them again.
						for (i = 0; i < nhyper; i++) {
							for (j = 0; j < nhyper; j++) {
								gsl_matrix_set(H, (size_t) i, (size_t) j, hessian[i + nhyper * j]);
							}
						}
						work = gsl_eigen_symmv_alloc((size_t) nhyper);
						gsl_eigen_symmv(H, eigen_values, eigen_vectors, work);
						gsl_eigen_symmv_free(work);
					} else {
						// old 
						for (i = 0; i < nhyper; i++) {
							for (j = i; j < nhyper; j++) {
								double sum = 0.0;
								for (k = 0; k < nhyper; k++) {
									sum +=
									    gsl_matrix_get(eigen_vectors, i, k) * gsl_matrix_get(eigen_vectors, j,
																 k)
									    * gsl_vector_get(eigen_values, k);
								}
								hessian[i + j * nhyper] = hessian[j + i * nhyper] = sum;
							}
						}
					}
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
				for (i = 0; i < nhyper; i++) {
					misc_output->eigenvalues[i] = 1.0 / gsl_vector_get(eigen_values, i);	/* need the eigenvalues of the
														 * cov.mat not hessian */
				}
				GMRFLib_gsl_mat2plain(&(misc_output->eigenvectors), eigen_vectors);
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
							fprintf(ai_par->fp_log, " %7.3f", val);
						} else {
							fprintf(ai_par->fp_log, " %7s", "");
						}
					}
					fprintf(ai_par->fp_log, "\n");
				}
			}

			// we need to save what we need from this section....
			// 
			if (GMRFLib_inla_mode == GMRFLib_MODE_TWOSTAGE_PART1) {
				assert(rpreopt);
				rpreopt->hessian = hessian;
				rpreopt->inverse_hessian = inverse_hessian;
				rpreopt->H = H;
				rpreopt->eigen_vectors = eigen_vectors;
				rpreopt->eigen_values = eigen_values;
				rpreopt->sqrt_eigen_values = sqrt_eigen_values;
				rpreopt->cov_m = misc_output->cov_m;
			}
		} else {
			// ... and pick it up here
			// 
			assert(GMRFLib_inla_mode == GMRFLib_MODE_TWOSTAGE_PART2);
			assert(rpreopt);
			hessian = rpreopt->hessian;
			inverse_hessian = rpreopt->inverse_hessian;
			H = rpreopt->H;
			assert(misc_output);
			misc_output->nhyper = nhyper;
			eigen_vectors = rpreopt->eigen_vectors;
			eigen_values = rpreopt->eigen_values;
			GMRFLib_gsl_mat2plain(&(misc_output->eigenvectors), eigen_vectors);
			GMRFLib_gsl_vec2plain(&(misc_output->eigenvalues), eigen_values);
			sqrt_eigen_values = rpreopt->sqrt_eigen_values;
			misc_output->cov_m = rpreopt->cov_m;
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

		iz = Calloc(nhyper, int);
		Memset(iz, 0, nhyper * sizeof(int));

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

		if (GMRFLib_inla_mode == GMRFLib_MODE_CLASSIC || GMRFLib_inla_mode == GMRFLib_MODE_TWOSTAGE_PART1) {
			/*
			 * compute the corrected scalings/stdevs, if required. 
			 */
			if ((ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_CCD)
			    || (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_GRID && density_hyper &&
				(ai_par->interpolator == GMRFLib_AI_INTERPOLATOR_CCD
				 || ai_par->interpolator == GMRFLib_AI_INTERPOLATOR_CCD_INTEGRATE))
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
						int thread_id = omp_get_thread_num();
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
						Memset(zz, 0, nhyper * sizeof(double));

						if (GMRFLib_OPENMP_IN_PARALLEL_ONEPLUS_THREAD()) {
							if (!ais[thread_id]) {
								ais[thread_id] =
								    GMRFLib_duplicate_ai_store(ai_store, GMRFLib_TRUE, GMRFLib_TRUE, GMRFLib_FALSE);
							}
							s = ais[thread_id];
						} else {
							s = ai_store;	/* the common one */
						}

						if (opt == 0) {
							zz[kk] = 2.0;
							GMRFLib_ai_z2theta(ttheta, nhyper, theta_mode, zz, sqrt_eigen_values, eigen_vectors);
							GMRFLib_opt_f_intern(thread_id, ttheta, &llog_dens, &ierr, s, NULL, NULL);
							llog_dens *= -1.0;
							f0 = log_dens_mode - llog_dens;
							stdev_corr_pos[kk] = (f0 > 0.0 ? sqrt(2.0 / f0) : 1.0);
						} else {
							zz[kk] = -2.0;
							GMRFLib_ai_z2theta(ttheta, nhyper, theta_mode, zz, sqrt_eigen_values, eigen_vectors);
							GMRFLib_opt_f_intern(thread_id, ttheta, &llog_dens, &ierr, s, NULL, NULL);
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
						int thread_id = omp_get_thread_num();
						double f0, *zz = NULL, *ttheta = NULL, llog_dens;
						GMRFLib_ai_store_tp *s = NULL;
						zz = Calloc(nhyper, double);
						ttheta = Calloc(nhyper, double);
						Memset(zz, 0, nhyper * sizeof(double));

						if (GMRFLib_OPENMP_IN_PARALLEL_ONEPLUS_THREAD()) {
							if (!ais[thread_id]) {
								ais[thread_id] =
								    GMRFLib_duplicate_ai_store(ai_store, GMRFLib_TRUE, GMRFLib_TRUE, GMRFLib_FALSE);
							}
							s = ais[thread_id];
						} else {
							s = ai_store;	/* the common one */
						}

						zz[k] = 2.0;
						GMRFLib_ai_z2theta(ttheta, nhyper, theta_mode, zz, sqrt_eigen_values, eigen_vectors);
						GMRFLib_opt_f_intern(thread_id, ttheta, &llog_dens, &ierr, s, NULL, NULL);
						llog_dens *= -1.0;
						f0 = log_dens_mode - llog_dens;
						stdev_corr_pos[k] = (f0 > 0.0 ? sqrt(2.0 / f0) : 1.0);

						zz[k] = -2.0;
						GMRFLib_ai_z2theta(ttheta, nhyper, theta_mode, zz, sqrt_eigen_values, eigen_vectors);
						GMRFLib_opt_f_intern(thread_id, ttheta, &llog_dens, &ierr, s, NULL, NULL);
						llog_dens *= -1.0;
						f0 = log_dens_mode - llog_dens;
						stdev_corr_neg[k] = (f0 > 0.0 ? sqrt(2.0 / f0) : 1.0);

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
		} else {
			assert(GMRFLib_inla_mode == GMRFLib_MODE_TWOSTAGE_PART2);
		}

		if (GMRFLib_inla_mode == GMRFLib_MODE_TWOSTAGE_PART1) {
			rpreopt->stdev_corr_neg = stdev_corr_neg;
			rpreopt->stdev_corr_pos = stdev_corr_pos;
		} else if (GMRFLib_inla_mode == GMRFLib_MODE_TWOSTAGE_PART2) {
			stdev_corr_neg = rpreopt->stdev_corr_neg;
			stdev_corr_pos = rpreopt->stdev_corr_pos;
			if (misc_output && !(misc_output->stdev_corr_pos)) {
				misc_output->stdev_corr_pos = Calloc(nhyper, double);
				misc_output->stdev_corr_neg = Calloc(nhyper, double);
				Memcpy(misc_output->stdev_corr_pos, stdev_corr_pos, nhyper * sizeof(double));
				Memcpy(misc_output->stdev_corr_neg, stdev_corr_neg, nhyper * sizeof(double));
			}
		}

		for (k = 0; k < nhyper; k++) {
			if (ai_par->fp_log) {
				fprintf(ai_par->fp_log,
					"%s corrected stdev for theta[%1d]: negative %.3f  positive %.3f\n",
					(GMRFLib_inla_mode == GMRFLib_MODE_TWOSTAGE_PART2 ? "Retrive" : "Compute"),
					k, stdev_corr_neg[k], stdev_corr_pos[k]);
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
		if (x_mode && ai_store->mode) {
			Memcpy(x_mode, ai_store->mode, graph->n * sizeof(double));
		}
		GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_OPTIMIZE, (void *) &nhyper, NULL);
		GMRFLib_opt_f(0, theta_mode, &log_dens_mode, &ierr, NULL, NULL);
		log_dens_mode *= -1.0;
		misc_output->log_posterior_mode = log_dens_mode;

		SET_THETA_MODE;
		if (x_mode) {
			Memcpy(x_mode, ai_store->mode, graph->n * sizeof(double));
		}

		if (timer) {
			timer[1] = GMRFLib_cpu() - timer[1];
			timer[2] = GMRFLib_cpu();
		}
		GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_INTEGRATE_HYPERPAR, NULL, NULL);

		if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_EMPIRICAL_BAYES) {
			if (need_Qinv) {
				GMRFLib_ai_add_Qinv_to_ai_store(ai_store);	/* add Qinv if required */
			}
			GMRFLib_ai_store_tp **ai_store_id = Calloc(GMRFLib_MAX_THREADS(), GMRFLib_ai_store_tp *);

			int thread_id = 0;
			assert(omp_get_thread_num() == 0);
#pragma omp parallel for private(i) num_threads(GMRFLib_openmp->max_threads_outer)
			for (i = 0; i < compute_n; i++) {
				int ii = compute_idx[i];
				int id = omp_get_thread_num();
				GMRFLib_density_tp *cpodens = NULL;
				if (!ai_store_id[id]) {
					ai_store_id[id] = GMRFLib_duplicate_ai_store(ai_store, GMRFLib_FALSE, GMRFLib_TRUE, GMRFLib_TRUE);
				}
				GMRFLib_ai_marginal_hidden(thread_id, &dens[ii][dens_count], (cpo && (d[ii]
												      || ai_par->cpo_manual) ? &cpodens : NULL),
							   GMRFLib_TRUE,
							   ii, x, b, c, mean, d, loglFunc, loglFunc_arg, graph,
							   Qfunc, Qfunc_arg, constr, ai_par, ai_store_id[id], preopt);
				if (tfunc && tfunc[ii]) {
					GMRFLib_transform_density(&dens_transform[ii][dens_count], dens[ii][dens_count], tfunc[ii]);
				}
				double *xx_mode = ai_store_id[id]->mode;
				COMPUTE2;
				GMRFLib_free_density(cpodens);
			}
			for (i = 0; i < GMRFLib_MAX_THREADS(); i++) {
				if (!ai_store_id[i]) {
					GMRFLib_free_ai_store(ai_store_id[i]);
				}
			}
			Free(ai_store_id);

			if (ai_par->vb_enable) {
				GMRFLib_ai_vb_correct_mean(thread_id, dens, dens_count, NULL,
							   c, d, ai_par, ai_store, graph, Qfunc, Qfunc_arg, loglFunc, loglFunc_arg, preopt);
			}

			if (GMRFLib_ai_INLA_userfunc0) {
				userfunc_values[dens_count] = GMRFLib_ai_INLA_userfunc0(thread_id, ai_store->problem, theta, nhyper);
			}
			COMPUTE_LINDENS(ai_store, GMRFLib_TRUE);
			ADD_CONFIG(ai_store, theta_mode, log_dens_mode, log_dens_mode);

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
				GMRFLib_design_ccd(&design, nhyper);
			} else {
				design = ai_par->int_design;
			}

			f = DMAX(ai_par->f0, 1.0) * sqrt((double) nhyper);
			w = 1.0 / ((design->nexperiments - 1.0) * (1.0 + exp(-0.5 * SQR(f)) * (SQR(f) / nhyper - 1.0)));
			w_origo = 1.0 - (design->nexperiments - 1.0) * w;

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
				int thread_id = omp_get_thread_num();

				double *z_local, *theta_local, log_dens_orig;
				GMRFLib_ai_store_tp *ai_store_id = NULL;
				GMRFLib_tabulate_Qfunc_tp *tabQfunc = NULL;
				double *bnew = NULL;

				dens_count = k;
				hyper_count = k;

				if (GMRFLib_OPENMP_IN_PARALLEL_ONEPLUS_THREAD()) {
					if (!ais[thread_id]) {
						ais[thread_id] = GMRFLib_duplicate_ai_store(ai_store, GMRFLib_FALSE, GMRFLib_TRUE, GMRFLib_FALSE);
					}
					ai_store_id = ais[thread_id];
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
					Memcpy(theta_local, z_local, nhyper * sizeof(double));
				}
				GMRFLib_opt_f_intern(thread_id, theta_local, &log_dens, &ierr, ai_store_id, &tabQfunc, &bnew);
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
					if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_USER_EXPERT) {
						log_dens = log(design->int_weight[k]) + log_dens_mode;
					} else {
						log_dens += log(design->int_weight[k]);
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
				for (i = 0; i < compute_n; i++) {
					int ii = compute_idx[i];
					GMRFLib_density_tp *cpodens = NULL;

					GMRFLib_ai_marginal_hidden(thread_id, &dens[ii][dens_count],
								   (cpo && (d[ii] || ai_par->cpo_manual) ? &cpodens : NULL),
								   GMRFLib_FALSE,
								   ii, x, bnew, c, mean, d,
								   loglFunc, loglFunc_arg,
								   graph, tabQfunc->Qfunc, tabQfunc->Qfunc_arg,
								   constr, ai_par, ai_store_id, preopt);
					if (tfunc && tfunc[ii]) {
						GMRFLib_transform_density(&dens_transform[ii][dens_count], dens[ii][dens_count], tfunc[ii]);
					}
					double *xx_mode = ai_store_id->mode;
					COMPUTE;
					GMRFLib_free_density(cpodens);
				}
				if (ai_par->vb_enable) {
					GMRFLib_ai_vb_correct_mean(thread_id, dens, dens_count, NULL,
								   c, d, ai_par, ai_store_id, graph, tabQfunc->Qfunc,
								   tabQfunc->Qfunc_arg, loglFunc, loglFunc_arg, preopt);
				}

				if (GMRFLib_ai_INLA_userfunc0) {
					assert(ai_store_id);
					userfunc_values[dens_count] =
					    GMRFLib_ai_INLA_userfunc0(thread_id, ai_store_id->problem, theta_local, nhyper);
				}
				COMPUTE_LINDENS(ai_store_id, GMRFLib_FALSE);
				ADD_CONFIG(ai_store_id, theta_local, log_dens, log_dens_orig);
				tu = GMRFLib_cpu() - tref;
				if (ai_par->fp_log) {
#pragma omp critical (Name_ef0a51d4c8b775aad87f67bb17dd269e61620fab)
					{
						fprintf(ai_par->fp_log, "config %2d/%1d=[", config_count++, design->nexperiments);
						for (i = 0; i < nhyper; i++) {
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
			}

			/*
			 * set the values back 
			 */
			dens_count = design->nexperiments;
			hyper_count = design->nexperiments;

			/*
			 * END OF GMRFLib_AI_INT_STRATEGY_CCD / USER / USER_STD
			 */
		} else {
			/*
			 * integrate using GRID 
			 */
			/*
			 * new code which parallise over configurations 
			 */
			GMRFLib_ai_pool_tp *pool = NULL;
			GMRFLib_ai_pool_init(&pool, ai_par, nhyper);

			GMRFLib_ASSERT(dens_count == 0, GMRFLib_ESNH);
			GMRFLib_ASSERT(hyper_count == 0, GMRFLib_ESNH);

#pragma omp parallel for private(i, log_dens, tref, tu, ierr) num_threads(GMRFLib_openmp->max_threads_outer)
			for (size_t kk = 0; kk < pool->nconfig; kk++) {
				int thread_id = omp_get_thread_num();

				GMRFLib_ai_store_tp *ai_store_id = NULL;
				GMRFLib_density_tp **dens_local = NULL;
				GMRFLib_density_tp **dens_local_transform = NULL;
				double *z_local = NULL, *theta_local = NULL, *userfunc_values_local = NULL, weights_local, val;
				double *cpo_theta_local = NULL, *po_theta_local = NULL, *po2_theta_local =
				    NULL, *po3_theta_local = NULL, *pit_theta_local = NULL, *failure_theta_local =
				    NULL, **deviance_theta_local = NULL;
				int err, *iz_local = NULL;
				size_t idx = 0;
				GMRFLib_tabulate_Qfunc_tp *tabQfunc = NULL;
				double *bnew = NULL;

				iz_local = Calloc(nhyper, int);
				err = GMRFLib_ai_pool_get(pool, iz_local, &idx);
				/*
				 * if we get a new config, then go on, otherwise, do nothing 
				 */
				if (err == GMRFLib_SUCCESS) {
					tref = GMRFLib_cpu();
					if (GMRFLib_OPENMP_IN_PARALLEL_ONEPLUS_THREAD()) {
						if (!ais[thread_id]) {
							ais[thread_id] =
							    GMRFLib_duplicate_ai_store(ai_store, GMRFLib_FALSE, GMRFLib_TRUE, GMRFLib_FALSE);
						}
						ai_store_id = ais[thread_id];
					} else {
						ai_store_id = ai_store;	/* the common one */
					}

					z_local = Calloc(nhyper, double);
					theta_local = Calloc(nhyper, double);
					for (i = 0; i < nhyper; i++) {
						z_local[i] = iz_local[i] * ai_par->dz;
					}
					GMRFLib_ai_z2theta(theta_local, nhyper, theta_mode, z_local, sqrt_eigen_values, eigen_vectors);
					GMRFLib_opt_f_intern(thread_id, theta_local, &log_dens, &ierr, ai_store_id, &tabQfunc, &bnew);
					log_dens *= -1.0;

					val = log_dens - log_dens_mode;
					if ((ISINF(val) || ISNAN(val)) || -val > ai_par->diff_log_dens) {
						GMRFLib_ai_pool_set(pool, idx, val);
						if (ai_par->fp_log) {
#pragma omp critical (Name_241fb67d09a4cff6907d0608d097c348251195e9)
							{
								fprintf(ai_par->fp_log, "config %2d=[", config_count++);
								for (i = 0; i < nhyper; i++) {
									fprintf(ai_par->fp_log, " %6.3f", z_local[i]);
								}
								fprintf(ai_par->fp_log,
									" ] log(rel.dens)= %6.3f, reject, %.2fs\n", val, GMRFLib_cpu() - tref);
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
							deviance_theta_local = Calloc(graph->n, double *);
						}

						for (i = 0; i < compute_n; i++) {
							GMRFLib_density_tp *cpodens = NULL;
							int ii;
							double *xx_mode = NULL;

							ii = compute_idx[i];
							GMRFLib_ai_marginal_hidden(thread_id, &dens_local[ii], (cpo && (d[ii]
															|| ai_par->cpo_manual)
														? &cpodens : NULL),
										   GMRFLib_FALSE, ii,
										   x, bnew, c, mean, d, loglFunc,
										   loglFunc_arg, graph,
										   tabQfunc->Qfunc, tabQfunc->Qfunc_arg,
										   constr, ai_par, ai_store_id, preopt);
							if (tfunc && tfunc[ii]) {
								GMRFLib_transform_density(&dens_local_transform[ii], dens_local[ii], tfunc[ii]);
							}
							xx_mode = ai_store_id->mode;

							COMPUTE_LOCAL;
							GMRFLib_free_density(cpodens);
						}
						if (ai_par->vb_enable) {
							GMRFLib_ai_vb_correct_mean(thread_id, NULL, -1, dens_local, c, d,
										   ai_par, ai_store_id, graph, tabQfunc->Qfunc,
										   tabQfunc->Qfunc_arg, loglFunc, loglFunc_arg, preopt);
						}

						GMRFLib_density_tp **lin_dens_local = NULL;
						double *lin_cross_local = NULL;
						COMPUTE_LINDENS_LOCAL(ai_store_id, GMRFLib_FALSE);

						if (GMRFLib_ai_INLA_userfunc0) {
							userfunc_values_local =
							    GMRFLib_ai_INLA_userfunc0(thread_id, ai_store_id->problem, theta_local, nhyper);
						}
						tu = GMRFLib_cpu() - tref;

// not the best, but at least now its all about storing only...
#pragma omp critical (Name_36b1b7dfeb7a205ea072f283e7f5ed9408c3aca1)
						{
							int ii;
							if (ai_par->fp_log) {
								{
									fprintf(ai_par->fp_log, "config %2d=[", config_count++);
									for (i = 0; i < nhyper; i++) {
										fprintf(ai_par->fp_log, " %6.3f", z_local[i]);
									}
									fprintf(ai_par->fp_log,
										" ] log(rel.dens)= %6.3f, [%1d] accept, compute,",
										val, omp_get_thread_num());
									fprintf(ai_par->fp_log, " %.2fs\n", tu);
								}
							}
							CHECK_DENS_STORAGE;
							CHECK_HYPER_STORAGE;

							weights[dens_count] = weights_local;
							Memcpy(&hyper_z[hyper_count * nhyper], z_local, nhyper * sizeof(double));
							hyper_ldens[hyper_count] = log_dens - log_dens_mode;
							izs[dens_count] = Calloc(nhyper, double);
							Memcpy(izs[dens_count], z_local, nhyper * sizeof(double));

							for (i = 0; i < compute_n; i++) {
								ii = compute_idx[i];
								dens[ii][dens_count] = dens_local[ii];
								if (tfunc && tfunc[ii]) {
									dens_transform[ii][dens_count] = dens_local_transform[ii];
								}
							}
							if (nlin) {
								lin_dens[dens_count] = lin_dens_local;
								if (lin_cross) {
									lin_cross[dens_count] = lin_cross_local;
								}
							}
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

		int thread_id = 0;
		assert(omp_get_thread_num() == 0);
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

		/*
		 * In this case the contents of ai_store is NULL, we need to recompute the Gaussian approximation since
		 * the contents of ai_store is NULL in this case.
		 */
		double tmp_logdens;
		double *bnew = NULL, con = 0.0;
		GMRFLib_bnew(thread_id, &bnew, &con, graph->n, b, bfunc);
		GMRFLib_ai_marginal_hyperparam(thread_id, &tmp_logdens, x, bnew, c, mean, d,
					       loglFunc, loglFunc_arg, graph, Qfunc, Qfunc_arg, constr, ai_par, ai_store, preopt);
		log_dens_mode = tmp_logdens + con + log_extra(thread_id, NULL, nhyper, log_extra_arg);

		GMRFLib_ai_add_Qinv_to_ai_store(ai_store);
		Free(bnew);

		GMRFLib_ai_store_tp **ai_store_id = Calloc(GMRFLib_MAX_THREADS(), GMRFLib_ai_store_tp *);
		GMRFLib_bnew(thread_id, &bnew, &con, graph->n, b, bfunc);
#pragma omp parallel for private(i) num_threads(GMRFLib_openmp->max_threads_outer)
		for (i = 0; i < compute_n; i++) {
			int id = omp_get_thread_num();
			int ii = compute_idx[i];
			GMRFLib_density_tp *cpodens = NULL;

			if (!ai_store_id[id]) {
				ai_store_id[id] = GMRFLib_duplicate_ai_store(ai_store, GMRFLib_FALSE, GMRFLib_TRUE, GMRFLib_TRUE);
			}
			GMRFLib_ai_marginal_hidden(thread_id, &dens[ii][dens_count], (cpo && (d[ii]
											      || ai_par->cpo_manual) ? &cpodens : NULL),
						   GMRFLib_TRUE, ii, x, bnew, c, mean, d, loglFunc, loglFunc_arg, graph, Qfunc, Qfunc_arg, constr,
						   ai_par, ai_store_id[id], preopt);
			if (tfunc && tfunc[ii]) {
				GMRFLib_transform_density(&dens_transform[ii][dens_count], dens[ii][dens_count], tfunc[ii]);
			}
			double *xx_mode = ai_store_id[id]->mode;
			COMPUTE2;
			GMRFLib_free_density(cpodens);
		}

		int id_nz = 0;
		for (id_nz = 0; id_nz < GMRFLib_MAX_THREADS(); id_nz++) {
			if (ai_store_id[id_nz]) {
				break;
			}
		}
		if (id_nz < GMRFLib_MAX_THREADS()) {
			Memcpy(x_mode, ai_store_id[id_nz]->mode, graph->n * sizeof(double));
		} else {
			Memset(x_mode, 0, graph->n * sizeof(double));
		}

		for (i = 0; i < GMRFLib_MAX_THREADS(); i++) {
			if (!ai_store_id[i]) {
				GMRFLib_free_ai_store(ai_store_id[i]);
			}
		}
		Free(ai_store_id);
		Free(bnew);

		if (ai_par->vb_enable) {
			GMRFLib_ai_vb_correct_mean(thread_id, dens, dens_count, NULL, c, d, ai_par, ai_store, graph, Qfunc, Qfunc_arg,
						   loglFunc, loglFunc_arg, preopt);
		}

		if (GMRFLib_ai_INLA_userfunc0) {
			assert(userfunc_values);
			userfunc_values[dens_count] = GMRFLib_ai_INLA_userfunc0(thread_id, ai_store->problem, theta, nhyper);
		}
		COMPUTE_LINDENS(ai_store, GMRFLib_TRUE);
		ADD_CONFIG(ai_store, NULL, log_dens_mode, log_dens_mode);
		weights[dens_count] = 0.0;
		dens_count++;

		/*
		 * END OF nhyper == 0 
		 */
	}

	if (preopt) {
		// in this case, just save (x, theta) adding the predictors
		preopt->mode_theta = Calloc(nhyper, double);
		Memcpy(preopt->mode_theta, theta_mode, nhyper * sizeof(double));
		preopt->mode_x = Calloc(preopt->mnpred + preopt->n, double);
		// GMRFLib_opt_get_latent(&(preopt->mode_x[preopt->mnpred]));
		Memcpy(&(preopt->mode_x[preopt->mnpred]), x_mode, preopt->n * sizeof(double));
		GMRFLib_preopt_full_predictor(preopt->mode_x, &(preopt->mode_x[preopt->mnpred]), preopt);
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
		Memcpy(adj_weights, weights, dens_count * sizeof(double));
	}
	// need to adjust the weights in configs if we have adjusted the integration weights. this is a bit unfortunate, but the adjustment
	// needs too weights to be present before adjusting.
	if (misc_output->configs) {
		for (int dc = 0; dc < dens_count; dc++) {
			int found = 0;
			for (int id = 0; id < GMRFLib_MAX_THREADS(); id++) {
				if (misc_output->configs[id]) {
					for (i = 0; i < misc_output->configs[id]->nconfig; i++) {
						if (misc_output->configs[id]->config[i]->dens_count == dc) {
							misc_output->configs[id]->config[i]->log_posterior = log(adj_weights[dc]);
							found++;	/* for the check below */
						}
					}
				}
			}
			assert(found == 1);		       /* just a check... */
		}
	}

	if (ai_par->fp_log) {
		fprintf(ai_par->fp_log, "\nCombine the densities with relative weights:\n");
		for (j = 0; j < dens_count; j++) {
			fprintf(ai_par->fp_log, "config %2d/%2d=[", j, dens_count);
			for (k = 0; k < nhyper; k++) {
				fprintf(ai_par->fp_log, " %6.3f", izs[j][k]);
			}
			fprintf(ai_par->fp_log, " ] weight= %6.3f", weights[j]);
			if (ai_par->adjust_weights && (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_GRID)) {
				fprintf(ai_par->fp_log, " adjusted weight= %6.3f", adj_weights[j]);
			}
			fprintf(ai_par->fp_log, "\n");
		}
	}


	GMRFLib_normalize(dens_count, adj_weights);
	GMRFLib_idxval_tp *probs = GMRFLib_density_prune_weights(adj_weights, dens_count, GMRFLib_weight_prob);

	if (density) {
		/*
		 * need a separate strategy here, as this might take time if the number of points are large, and we essentially want to do this loop with max
		 * num_threads.
		 */
		GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_COMBINE, NULL, NULL);
#pragma omp parallel for private(j) num_threads(GMRFLib_openmp->max_threads_outer)
		for (j = 0; j < compute_n; j++) {
			int ii = compute_idx[j];
			GMRFLib_density_tp *dens_combine = NULL;
			GMRFLib_density_combine(&dens_combine, dens[ii], probs);
			if (density) {
				(*density)[ii] = dens_combine;
			}

			if (tfunc && tfunc[ii]) {
				GMRFLib_density_tp *dens_c = NULL;
				GMRFLib_density_combine((density_transform ? &dens_c : NULL), (density_transform ? dens_transform[ii] : NULL),
							probs);
				if (density_transform && *density_transform) {
					(*density_transform)[ii] = dens_c;
				}
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
			GMRFLib_density_combine(&dcombine, dtmp, probs);
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
			Memcpy(misc_output->cov_lin, ptmp, ISQR(nlin) * sizeof(double));

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
		fprintf(ai_par->fp_log, "\n");
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
		SET_THETA_MODE;
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
		SET_THETA_MODE;
		double mean_deviance = 0.0, mean_deviance_sat = 0.0, deviance_mean = 0.0, deviance_mean_sat = 0.0, *x_vec = NULL;

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
		    *deviance_e = Calloc(ndev, double), *deviance_e_sat = Calloc(ndev, double), *sign = Calloc(ndev, double);

		for (j = 0; j < ndev; j++) {
			e_deviance[j] = e_deviance_sat[j] = deviance_e[j] = deviance_e_sat[j] = sign[j] = NAN;
		}

		for (j = 0; j < compute_n; j++) {
			double md = 0.0, md_sat = 0.0, dm = 0.0, dm_sat = 0.0, logl_sat = 0.0;
			int ii = compute_idx[j];
			int thread_id = omp_get_thread_num();
			assert(thread_id == 0);

			if (d[ii]) {
				double evalue = 0.0, evalue_sat = 0.0, sum = 0.0, logl;
				for (int jj = 0; jj < dens_count; jj++) {
					evalue += deviance_theta[ii][jj][0] * adj_weights[jj];
					evalue_sat += deviance_theta[ii][jj][1] * adj_weights[jj];
					sum += adj_weights[jj];
				}
				md = evalue / sum;
				md_sat = evalue_sat / sum;

				if (!(density && (*density)[ii])) {
					fprintf(stderr, "\n\n\nFIXME FIXME!!!!!!!!\n\n\n");
					abort();
				}

				double x_tmp = (double) ((*density)[ii]->user_mean);
				loglFunc(thread_id, &logl, &x_tmp, 1, ii, x_vec, NULL, loglFunc_arg, NULL);
				logl_sat = inla_compute_saturated_loglik(thread_id, ii, loglFunc, x_vec, loglFunc_arg);
				dm = -2.0 * d[ii] * logl;
				dm_sat = -2.0 * d[ii] * (logl - logl_sat);
				e_deviance[ii] = md;
				e_deviance_sat[ii] = md_sat;
				deviance_e[ii] = dm;
				deviance_e_sat[ii] = dm_sat;

				// neither of these options are fail-safe. I cannot see how to do this fail-safe without really mapping to the
				// real data doing the comparison there. But this information is not available at this level
				double sig = 0.0;
				if (loglFunc(0, NULL, NULL, 0, ii, NULL, NULL, loglFunc_arg, NULL) == GMRFLib_LOGL_COMPUTE_CDF) {
					loglFunc(0, &sig, &((*density)[ii]->user_mean), -1, ii, NULL, NULL, loglFunc_arg, NULL);
					sig = (sig <= 0.5 ? -1.0 : 1.0);
				} else {
					double xx[2], ld[2] = { 0.0, 0.0 };
					xx[0] = (*density)[ii]->user_mean;
					xx[1] = xx[0] + 0.01 * (*density)[ii]->user_stdev;
					loglFunc(0, ld, xx, 2, ii, NULL, NULL, loglFunc_arg, NULL);
					sig = (ld[1] > ld[0] ? 1.0 : -1.0);
				}
				sign[ii] = sig;
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
			GMRFLib_density_create_normal(&(GMRFLib_ai_INLA_userfunc0_density[j]), 0.0, 1.0, mmean, ssd, GMRFLib_TRUE);

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
			for (k = 0; k < nhyper; k++) {
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

				izz = Calloc(nhyper, int);
				for (ntimes = 0; ntimes < 2; ntimes++) {
					int guard_count = 0;

					for (i = 0, len_length = 1; i < nhyper; i++) {
						k_maxx[i]++;   /* add one to construct the guard-zone */
						k_minn[i]--;   /* subtract one ... */
						len[i] = 1 + k_maxx[i] - k_minn[i];
						len_length *= len[i];
					}
					Memset(izz, 0, nhyper * sizeof(int));
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
							"\tLoop %1d: Added %1d guard-points with log-dens-diff equal to %.3f\n",
							ntimes, guard_count, hyper_ldens[hyper_count - 1]);
					}
				}
			}
			double *std_stdev_theta = Calloc(nhyper, double);

			for (k = 0; k < nhyper; k++) {
				std_stdev_theta[k] = sqrt(inverse_hessian[k + nhyper * k]);
			}

			/*
			 * write out the hole set 
			 */
			double *theta_tmp = NULL, log_jacobian = 0.0;
			theta_tmp = Calloc((int) IMAX(0, nhyper), double);

#define Amat(i_, j_) (rpreopt->int_design->A[ (i_) + (j_) * (rpreopt->int_design->nrow)])
			if (GMRFLib_inla_mode == GMRFLib_MODE_TWOSTAGE_PART1 && nhyper > 0) {
				rpreopt->int_design = Calloc(1, GMRFLib_matrix_tp);
				rpreopt->int_design->nrow = IMAX(1, hyper_count);
				rpreopt->int_design->ncol = nhyper + 1;
				rpreopt->int_design->A = Calloc(rpreopt->int_design->nrow * rpreopt->int_design->ncol, double);
				rpreopt->adj_weights = Calloc(dens_count, double);
				Memcpy(rpreopt->adj_weights, adj_weights, dens_count * sizeof(double));
			}

			if (eigen_values) {
				for (k = 0; k < nhyper; k++) {
					log_jacobian -= 0.5 * log(gsl_vector_get(eigen_values, (unsigned int) k));
				}
			}

			if (nhyper > 0 && hyper_count == 0) {
				if (rpreopt) {
					for (int kk = 0; kk < nhyper; kk++) {
						Amat(0, kk) = theta_mode[kk];
					}
					Amat(0, nhyper) = 1.0;
				}
			} else {
				for (k = 0; k < hyper_count; k++) {
					int kk;

					GMRFLib_ai_z2theta(theta_tmp, nhyper, theta_mode, &(hyper_z[k * nhyper]), sqrt_eigen_values, eigen_vectors);
					if (ai_par->fp_hyperparam) {
						for (kk = 0; kk < nhyper; kk++) {
							fprintf(ai_par->fp_hyperparam, " %.10g", theta_tmp[kk]);
						}
						fprintf(ai_par->fp_hyperparam, " %.10g %.10g\n", hyper_ldens[k] + log_dens_mode + log_jacobian,
							adj_weights[k]);
					}
					if (rpreopt) {
						for (kk = 0; kk < nhyper; kk++) {
							Amat(k, kk) = theta_tmp[kk];
						}
						Amat(k, nhyper) = adj_weights[k];
					}
				}
			}
#undef Amat
			if (ai_par->fp_hyperparam) {
				fflush(ai_par->fp_hyperparam);
			}
			Free(theta_tmp);

			if (ai_par->fp_log) {
				fprintf(ai_par->fp_log,
					"\tCompute the marginal for theta[%1d] to theta[%1d] using numerical integration...\n", 0, nhyper - 1);
			}

			for (k = 0; k < nhyper; k++) {
				GMRFLib_ai_marginal_one_hyperparamter(&((*density_hyper)[k]), k, nhyper, hyper_count, hyper_z,
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
	SET_THETA_MODE;
	if (x && x_mode) {
		Memcpy(x, x_mode, graph->n * sizeof(double));
	}

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

	/*
	 * userfunction1 
	 */
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

	if (GMRFLib_inla_mode == GMRFLib_MODE_CLASSIC) {
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
				Free(deviance_theta[j][0]);
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
		GMRFLib_opt_exit();
	}

	if (timer) {
		timer[3] = GMRFLib_cpu() - timer[3];
	}

	GMRFLib_LEAVE_ROUTINE;
#undef ADD_LINEAR_TERM
#undef ADD_LINEAR_TERM_LOCAL
#undef CHECK_DENS_STORAGE
#undef CHECK_DENS_STORAGE_FORCE
#undef CHECK_DENS_STORAGE_INTERN
#undef CHECK_HYPER_STORAGE
#undef COMPUTE
#undef COMPUTE2
#undef COMPUTE_CPO_AND_DIC
#undef COMPUTE_CPO_AND_DIC_LOCAL
#undef COMPUTE_LOCAL
#undef COMPUTE_NEFF
#undef COMPUTE_NEFF2
#undef COMPUTE_NEFF_LOCAL
#undef SET_THETA_MODE

	return GMRFLib_SUCCESS;
}

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
				 GMRFLib_preopt_tp * preopt)
{
#define SET_THETA_MODE							\
	if (theta_mode) {						\
		int i_, j_;						\
		for(j_=0; j_ < GMRFLib_MAX_THREADS(); j_++) {		\
			for(i_ = 0; i_ < nhyper; i_++){			\
				hyperparam[i_][j_][0] = theta_mode[i_]; \
			}						\
		}							\
	}

	int *k_max = NULL, *k_min = NULL, *k_maxx = NULL, *k_minn = NULL, ierr, *iz = NULL, *izz = NULL, *len =
	    NULL, *iz_axes = NULL, free_ai_par = 0, config_count = 0, dens_count = 0, dens_max, hyper_len = 0;

	double *hessian = NULL, *theta = NULL, *theta_mode = NULL, *x_mode = NULL, log_dens_mode = 0, log_dens, *z = NULL, **izs =
	    NULL, *stdev_corr_pos = NULL, *stdev_corr_neg = NULL, f, w, w_origo, tref, tu, *weights = NULL, *adj_weights =
	    NULL, *hyper_z = NULL, *hyper_ldens = NULL, **userfunc_values = NULL, *inverse_hessian = NULL, *timer,
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
	GMRFLib_idx_tp *d_idx = NULL;
	GMRFLib_gcpo_groups_tp *gcpo_groups = NULL;

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

	nhyper = IMAX(0, nhyper);
	ais = Calloc(GMRFLib_MAX_THREADS(), GMRFLib_ai_store_tp *);

	// need to determine dens_max
	GMRFLib_design_tp *tdesign = NULL;
	if (ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_CCD && nhyper > 0) {
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
			po_theta[j] = Calloc(dens_max, double);
			po2_theta[j] = Calloc(dens_max, double);
			po3_theta[j] = Calloc(dens_max, double);
		}
	}
	if (dic) {
		deviance_theta = Calloc(preopt->Npred, double **);	/* mean of deviance conditioned on theta */
		assert(d_idx);
		for (int jj = 0; jj < d_idx->n; jj++) {
			int j = d_idx->idx[jj];
			deviance_theta[j] = Calloc(dens_max, double *);
		}
	}

	if (timer) {
		timer[0] = GMRFLib_cpu() - timer[0];	       /* preparation */
		timer[1] = GMRFLib_cpu();
	}

	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_OPTIMIZE, (void *) &nhyper, NULL);
	GMRFLib_opt_setup(hyperparam, nhyper, log_extra, log_extra_arg, NULL, x, b, c, mean, bfunc, d, loglFunc,
			  loglFunc_arg, graph, Qfunc, Qfunc_arg, constr, ai_par, ai_store, preopt);

	if (nhyper > 0) {
		/*
		 * the first step is to locate the mode of \pi(\theta | y). here we use the opt-optimiser routine.  NOTE that this
		 * '_setup' ensure that ai_store is changed for each call to _opt_f. this is a bit dirty programming, but there is no
		 * good way to get around it for the moment.
		 */
		GMRFLib_opt_setup(hyperparam, nhyper, log_extra, log_extra_arg, NULL, x, b, c, mean, bfunc, d, loglFunc,
				  loglFunc_arg, graph, Qfunc, Qfunc_arg, constr, ai_par, ai_store, preopt);

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
			{
				int retval;
				int fd_save = ai_par->gradient_forward_finite_difference;
				if (ai_par->optimise_smart) {
					ai_par->gradient_forward_finite_difference = GMRFLib_TRUE;
					if (ai_par->fp_log) {
						fprintf(ai_par->fp_log, "Smart optimise part I: estimate gradient using forward differences\n");
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

		SET_THETA_MODE;
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

		if (ai_par->fixed_mode) {
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
						if (ai_par->fp_log)
							fprintf(ai_par->fp_log,
								"Mode not sufficient accurate; switch to a stupid local search strategy.\n");
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
		SET_THETA_MODE;
		if (x_mode) {
			Memcpy(x_mode, ai_store->mode, graph->n * sizeof(double));
		}

		if (stupid_mode_iter) {
			for (int i = 0; i < nhyper; i++) {
				theta_mode[i] = hyperparam[i][0][0];
			}
			int thread_id = 0;
			assert(omp_get_thread_num() == 0);
			GMRFLib_opt_f(thread_id, theta_mode, &log_dens_mode, &ierr, NULL, NULL);
			log_dens_mode *= -1.0;
			SET_THETA_MODE;
			if (x_mode) {
				Memcpy(x_mode, ai_store->mode, graph->n * sizeof(double));
			}
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
				fprintf(stderr, "\n");
				fprintf(stderr, "\t*** WARNING *** Eigenvalue %1d of the Hessian is %.6g < 0\n", i, eigv);
				fprintf(stderr, "\t*** WARNING *** This have consequence for the accurancy of the hyperpar\n");
				fprintf(stderr, "\t*** WARNING *** Continue with a diagonal Hessian.\n");
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

		if (ai_par->fixed_mode) {
			if (ai_par->fp_log) {
				fprintf(ai_par->fp_log, "fixed_mode=1, so artificially, scaling of sd is set to 1.0\n");
			}
			for (int i = 0; i < nhyper; i++) {
				stdev_corr_neg[i] = 1.0;
				stdev_corr_pos[i] = 1.0;
			}
		} else {
#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_outer)
			for (int k = 0; k < nhyper; k++) {
				int thread_id = omp_get_thread_num();

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

				zz[k] = 2.0;
				GMRFLib_ai_z2theta(ttheta, nhyper, theta_mode, zz, sqrt_eigen_values, eigen_vectors);
				GMRFLib_opt_f_intern(thread_id, ttheta, &llog_dens, &ierr, s, NULL, NULL);
				llog_dens *= -1.0;
				f0 = log_dens_mode - llog_dens;
				stdev_corr_pos[k] = (f0 > 0.0 ? sqrt(2.0 / f0) : 1.0);

				zz[k] = -2.0;
				GMRFLib_ai_z2theta(ttheta, nhyper, theta_mode, zz, sqrt_eigen_values, eigen_vectors);
				GMRFLib_opt_f_intern(thread_id, ttheta, &llog_dens, &ierr, s, NULL, NULL);
				llog_dens *= -1.0;
				f0 = log_dens_mode - llog_dens;
				stdev_corr_neg[k] = (f0 > 0.0 ? sqrt(2.0 / f0) : 1.0);

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
		SET_THETA_MODE;
		if (x_mode) {
			Memcpy(x_mode, ai_store->mode, graph->n * sizeof(double));
		}
	}

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

	if (x_mode && ai_store->mode) {
		Memcpy(x_mode, ai_store->mode, graph->n * sizeof(double));
	}

	if (timer) {
		timer[1] = GMRFLib_cpu() - timer[1];
		timer[2] = GMRFLib_cpu();
	}
	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_INTEGRATE_HYPERPAR, NULL, NULL);

	// if fixed_mode=1, then we need to use EB
	if (ai_par->fixed_mode) {
		if (ai_par->int_strategy != GMRFLib_AI_INT_STRATEGY_EMPIRICAL_BAYES) {
			ai_par->int_strategy = GMRFLib_AI_INT_STRATEGY_EMPIRICAL_BAYES;
			if (ai_par->fp_log) {
				fprintf(ai_par->fp_log, "int.strategy is set to EB, since fixed_mode=1\n");
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
	assert(dens_max == design->nexperiments);

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
					char *s =
					    GMRFLib_strdup("Using Predictor and/or Apredictor in lincomb in 'experimental-mode' is not supported.");
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
		gcpo_groups = GMRFLib_gcpo_build(thread_id, ai_store, preopt, gcpo_param);
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

		double *z_local, *theta_local, log_dens_orig;
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

		if (nhyper > 0 || GMRFLib_OPENMP_IN_PARALLEL_ONEPLUS_THREAD()) {
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

		tref = GMRFLib_cpu();
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
							  c, d, ai_par, ai_store_id, graph,
							  (tabQfunc ? tabQfunc->Qfunc : Qfunc), (tabQfunc ? tabQfunc->Qfunc_arg : Qfunc_arg),
							  loglFunc, loglFunc_arg, preopt);
		}

		if (ai_par->vb_enable && (ai_par->vb_strategy == GMRFLib_AI_VB_VARIANCE)) {
			GMRFLib_ai_vb_correct_variance_preopt(thread_id, dens, dens_count,
							      c, d, ai_par, ai_store_id, graph,
							      (tabQfunc ? tabQfunc->Qfunc : Qfunc), (tabQfunc ? tabQfunc->Qfunc_arg : Qfunc_arg),
							      loglFunc, loglFunc_arg, preopt);
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
				for (int ii = 0; ii < 3 * preopt->Npred; ii++) {
					gcpodens_moments[ii] = NAN;
				}
			}
			gcpo_theta[dens_count] = GMRFLib_gcpo(thread_id, ai_store_id, lpred_mean, lpred_mode, lpred_variance, preopt, gcpo_groups,
							      d, loglFunc, loglFunc_arg, ai_par, gcpo_param, gcpodens_moments);
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
			for (int ii = 0; ii < 3 * preopt->Npred; ii++) {
				cpodens_moments[ii] = NAN;
			}
		}

		if (cpo || dic || po) {
#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_inner)
			for (int ii = 0; ii < d_idx->n; ii++) {
				int i = d_idx->idx[ii];
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
			for (int j = 0; j < preopt->Npred; j++) {
				int jj = 3 * j;
				if (d[j]) {
					GMRFLib_2order_taylor(thread_id, NULL, &(ll_info[jj]), &(ll_info[jj + 1]), &(ll_info[jj + 2]), d[j],
							      lpred_mode[j], j, lpred_mode, loglFunc, loglFunc_arg, &ai_par->step_len,
							      &ai_par->stencil);
				} else {
					ll_info[jj] = ll_info[jj + 1] = ll_info[jj + 2] = NAN;
				}
			}
		}

		GMRFLib_ai_store_config_preopt(thread_id, misc_output, nhyper, theta_local, log_dens, log_dens_orig, ai_store_id->problem,
					       mean_corrected, preopt, Qfunc, Qfunc_arg, cpodens_moments, gcpodens_moments, arg_str, ll_info);

		tu = GMRFLib_cpu() - tref;
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
		Free(lpred_mean);
		Free(lpred_variance);
		Free(mean_corrected);
		// do not free cpodens_moments!!!
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
		timer[2] = GMRFLib_cpu() - timer[2];
		timer[3] = GMRFLib_cpu();
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
			GMRFLib_density_combine(&dens_combine, lpred[i], probs_combine);
			(*density)[i] = dens_combine;
		} else {
			i = ii - preopt->mnpred;
			GMRFLib_density_tp *dens_combine = NULL;
			GMRFLib_density_combine(&dens_combine, dens[i], probs);
			(*density)[ii] = dens_combine;	       /* yes, its 'ii' */
			if (tfunc && tfunc[i]) {
				GMRFLib_density_tp *dens_c = NULL;
				GMRFLib_density_combine(&dens_c, dens_transform[i], probs_combine);
				(*density_transform)[i] = dens_c;
			}
		}
	}
	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_DEFAULT, NULL, NULL);

	if (dlin && nlin) {
		GMRFLib_density_tp **dtmp, *dcombine;

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
			double *ptmp;
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

			double *ptmp_scale = Calloc(ISQR(nlin), double);
			for (int i = 0; i < nlin; i++) {
				ptmp_scale[i + i * nlin] = 1.0 / sqrt(ptmp[i + i * nlin]);
			}

			for (int i = 0; i < nlin; i++) {
				for (int j = i + 1; j < nlin; j++) {
					ptmp[i + j * nlin] = ptmp[i + j * nlin] * ptmp_scale[i + i * nlin] * ptmp_scale[j + j * nlin];
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

		for (int j = 0; j < preopt->Npred; j++) {
			int ii = j;
			if (cpo_theta[ii]) {
				for (int jjj = 0; jjj < probs->n; jjj++) {
					int jj = probs->idx[jjj];
					if (!ISNAN(cpo_theta[ii][jj]))	/* we ignore those that have failed */
						Z[ii] += probs->val[jjj] / cpo_theta[ii][jj];
				}
			}
		}

		for (int j = 0; j < preopt->Npred; j++) {
			double evalue, evalue2, evalue_one = 1.0;
			int ii = j;

			if (cpo_theta[ii]) {
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
		SET_THETA_MODE;
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
		SET_THETA_MODE;

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

		for (int j = 0; j < ndev; j++) {
			e_deviance[j] = e_deviance_sat[j] = deviance_e[j] = deviance_e_sat[j] = sign[j] = NAN;
		}

		for (int j = 0; j < d_idx->n; j++) {
			double md = 0.0, md_sat = 0.0, dm = 0.0, dm_sat = 0.0, logl_sat = 0.0;
			int ii = d_idx->idx[j];

			double evalue = 0.0, evalue_sat = 0.0;
			for (int jjj = 0; jjj < probs->n; jjj++) {
				int jj = probs->idx[jjj];
				evalue += deviance_theta[ii][jj][0] * probs->val[jjj];
				evalue_sat += deviance_theta[ii][jj][1] * probs->val[jjj];
			}
			md = evalue;
			md_sat = evalue_sat;

			if (!(density && (*density)[ii])) {
				fprintf(stderr, "\n\n\nFIXME FIXME!!!!!!!!\n\n\n");
				abort();
			}

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
	SET_THETA_MODE;
	if (x && x_mode) {
		Memcpy(x, x_mode, graph->n * sizeof(double));
	}

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
		for (int j = 0; j < dens_max; j++) {
			GMRFLib_free_density(dens[i][j]);
		}
		Free(dens[i]);
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
		timer[3] = GMRFLib_cpu() - timer[3];
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
#undef SET_THETA_MODE

	return GMRFLib_SUCCESS;
}

GMRFLib_gcpo_groups_tp *GMRFLib_gcpo_build(int thread_id, GMRFLib_ai_store_tp * ai_store, GMRFLib_preopt_tp * preopt,
					   GMRFLib_gcpo_param_tp * gcpo_param)
{
	GMRFLib_ENTER_ROUTINE;
#define A_idx(node_) (preopt->pAA_idxval ? preopt->pAA_idxval[node_] : preopt->A_idxval[node_])
#define EQUAL_COR(c1_, c2_) (ABS((c1_) - (c2_)) < gcpo_param->epsilon)

	int detailed_output = GMRFLib_DEBUG_IF();
	int Npred = preopt->Npred;
	int mnpred = preopt->mnpred;
	int N = IMAX(preopt->n, Npred);
	GMRFLib_idxval_tp **groups = NULL;

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

			double big = 1.0 / GMRFLib_eps(1.0);

			// default
#pragma GCC ivdep
			for (int i = 0; i < nn; i++) {
				mask[i] = 1.0;
			}
			Memset(diag, 0, nn * sizeof(double));

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
#pragma GCC ivdep
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
#pragma GCC ivdep
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
		groups = GMRFLib_idxval_ncreate_x(Npred, IABS(gcpo_param->num_level_sets));
		GMRFLib_ai_add_Qinv_to_ai_store(ai_store);
		GMRFLib_ai_add_Qinv_to_ai_store(build_ai_store);

		GMRFLib_preopt_predictor_moments(NULL, isd, preopt, build_ai_store->problem, NULL);
#pragma GCC ivdep
		for (int i = 0; i < Npred; i++) {
			isd[i] = 1.0 / sqrt(isd[i]);
		}

		GMRFLib_idx_tp *selection = NULL;
		int free_selection = 0;

		if (!(gcpo_param->selection)) {
			for (int i = 0; i < Npred; i++) {
				GMRFLib_idx_add(&selection, i);
			}
			free_selection = 1;
		} else {
			selection = gcpo_param->selection;
			assert(GMRFLib_imax_value(selection->idx, selection->n, NULL) < Npred);
			assert(GMRFLib_imin_value(selection->idx, selection->n, NULL) >= 0);
		}
		if (gcpo_param->verbose || detailed_output) {
			assert(selection);
			printf("%s[%1d]: Use selection of %1d indices and num.level.sets %1d\n", __GMRFLib_FuncName,
			       omp_get_thread_num(), selection->n, gcpo_param->num_level_sets);
		}
		assert(selection);

#define CODE_BLOCK							\
		for(int ii = 0;	ii < selection->n; ii++) {		\
			CODE_BLOCK_ALL_WORK_ZERO();			\
			int node = selection->idx[ii];			\
									\
			if (gcpo_param->num_level_sets == -1) {		\
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
			for (int k = 0; k < v->n; k++) {		\
				a[v->idx[k]] = v->val[k];		\
			}						\
			GMRFLib_Qsolve(Sa, a, build_ai_store->problem, -1); \
			cor[node] = 1.0;				\
									\
			for (int nnode = 0; nnode < Npred; nnode++) {	\
				if (nnode == node) continue;		\
				GMRFLib_idxval_tp *vv = A_idx(nnode);	\
				double sum = GMRFLib_dot_product(vv, Sa); \
				sum *= isd[node] * isd[nnode];		\
				cor[nnode] = TRUNCATE(sum, -1.0, 1.0);	\
				cor_abs[nnode] = ABS(cor[nnode]);	\
			}						\
									\
			int levels_ok = 0;				\
			int levels_magnify = 1;				\
			cor[node] = cor_abs[node] = 1.0;		\
			while (!levels_ok) {				\
				GMRFLib_idxval_free(groups[node]);	\
				groups[node] = NULL;			\
				int siz_g = IMIN(Npred, levels_magnify * (IABS(gcpo_param->num_level_sets) + 4L)); \
				levels_magnify *= 10;			\
				GMRFLib_DEBUG_i_v("node siz_g Npred num_level_sets levels_magnify", node, siz_g, Npred, gcpo_param->num_level_sets, levels_magnify); \
				gsl_sort_largest_index(largest, (size_t) siz_g, cor_abs, (size_t) 1, (size_t) Npred); \
				int nlevels = 1;			\
				int i_prev = 0;				\
				double cor_abs_prev = 1.0;		\
				GMRFLib_idxval_add(&(groups[node]), (int) largest[i_prev], cor_abs_prev); \
				for (int i = 1; i < siz_g && !levels_ok; i++) {	\
					int i_new = (int) largest[i];	\
					double cor_abs_new = cor_abs[i_new]; \
					if (!EQUAL_COR(cor_abs_new, cor_abs_prev)) { \
						nlevels++;		\
						i_prev = i;		\
						cor_abs_prev = cor_abs_new; \
						if (nlevels <= gcpo_param->num_level_sets) { \
							GMRFLib_DEBUG_id("add new level  i_new cor_abs_new", i_new, cor_abs_new); \
							GMRFLib_idxval_add(&(groups[node]), i_new , cor[i_new]); \
						} else {		\
							levels_ok = 1;	\
						}			\
					} else {			\
						cor_abs[i_new] = cor_abs_prev; \
						cor[i_new] = SIGN(cor[i_new]) * cor_abs_prev; \
						GMRFLib_idxval_add(&(groups[node]), i_new, cor[i_new]); \
						GMRFLib_DEBUG_id("add to old level  i_new cor_abs_prev", i_new, cor_abs_prev); \
					}				\
					if (gcpo_param->size_max > 0 && groups[node]->n >= gcpo_param->size_max) { \
						levels_ok = 1;		\
					}				\
				}					\
				if (siz_g == Npred) levels_ok = 1;	\
				if (levels_ok) {			\
					if (gcpo_param->verbose || detailed_output) { \
						printf("%s[%1d]: for node=%1d, number of levels is %d\n", __GMRFLib_FuncName, omp_get_thread_num(), node,  nlevels); \
						printf("%s[%1d]: either because there are no more levels, or size.max is reached.\n", __GMRFLib_FuncName, omp_get_thread_num()); } \
				}					\
				GMRFLib_DEBUG_i("found nlevels", nlevels); \
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
			if (0) P(node);					\
			if (0) GMRFLib_idxval_printf(stdout, groups[node], "after adding friends"); \
		}

		RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 5, N);
#undef CODE_BLOCK

		if (free_selection) {
			GMRFLib_idx_free(selection);
		}
		Free(isd);
		GMRFLib_free_ai_store(local_ai_store);
	} else {
		if (gcpo_param->verbose || detailed_output) {
			printf("%s[%1d]: Use user-defined groups\n", __GMRFLib_FuncName, omp_get_thread_num());
		}

		// if the number of given groups are to short compared to Npred, then pad with empty groups. if its longer then that is an error.
		if (gcpo_param->ngroups < Npred) {
			gcpo_param->groups = Realloc(gcpo_param->groups, Npred, GMRFLib_idxval_tp *);
			for (int i = gcpo_param->ngroups; i < Npred; i++) {
				GMRFLib_idxval_create_x(&(gcpo_param->groups[i]), IABS(gcpo_param->num_level_sets));
			}
		} else if (gcpo_param->ngroups > Npred) {
			assert(gcpo_param->ngroups > Npred);
		}
		groups = gcpo_param->groups;
	}


	// add first off-diagonals
	GMRFLib_idx2_tp **missing = GMRFLib_idx2_ncreate_x(Npred, IABS(gcpo_param->num_level_sets));
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
				char *msg;
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

#undef A_idx
#undef EQUAL_COR

	GMRFLib_LEAVE_ROUTINE;
	return ggroups;
}

GMRFLib_gcpo_elm_tp **GMRFLib_gcpo(int thread_id, GMRFLib_ai_store_tp * ai_store_id, double *lpred_mean, double *lpred_mode,
				   double *lpred_variance, GMRFLib_preopt_tp * preopt,
				   GMRFLib_gcpo_groups_tp * groups, double *d, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg,
				   GMRFLib_ai_param_tp * ai_par, GMRFLib_gcpo_param_tp * gcpo_param, double *gcpodens_moments)
{
#define A_idx(node_) (preopt->pAA_idxval ? preopt->pAA_idxval[node_] : preopt->A_idxval[node_])

	GMRFLib_ENTER_ROUTINE;

	int detailed_output = GMRFLib_DEBUG_IF();
	int Npred = preopt->Npred;
	int mnpred = preopt->mnpred;
	int nn = preopt->n;
	int N = IMAX(nn, Npred);
	int max_ng = -1;
	int corr_hypar = gcpo_param->correct_hyperpar;
	const int np = GMRFLib_INT_GHQ_POINTS;
	double zero = 0.0;
	double spd_eps = GMRFLib_eps(0.5);
	double diag_eps = GMRFLib_eps(0.236);		       /* gives 0.0002, so the stdev is scaled with 1.0001 */
	double diag_scale = 1.0 + diag_eps;

	if (gcpo_param->verbose || detailed_output) {
		printf("enter _gcpo with...\n");
		for (int i = 0; i < Npred; i++) {
			printf("\ti %d lpred_mean %f lpred_mode %f lpred_variance %f\n", i, lpred_mean[i], lpred_mode[i], lpred_variance[i]);
		}
	}

	Calloc_init(mnpred, 1);
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
					for (int kk = 0; kk < v->n; kk++) { \
						a[v->idx[kk]] = v->val[kk]; \
					}				\
					GMRFLib_Qsolve(Sa, a, ai_store_id->problem, -1); \
					need_Sa = 0;			\
				}					\
									\
				GMRFLib_idxval_tp *v = A_idx(nnode);	\
				double sum = GMRFLib_dot_product(v, Sa); \
				double f = sd[node] * sd[nnode];	\
				sum /= f;				\
				double cov = TRUNCATE(sum, -1.0, 1.0) * f; \
									\
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

	if (detailed_output) {
#pragma omp critical (Name_0139eb204165e8e82ee3aaaaff59eab1d5b3cc14)
		{
			for (int node = 0; node < Npred; node++) {
				if (gcpo[node]->cov_mat && gcpo[node]->cov_mat->size1 > 0) {
					printf("\ncov_mat for node=%d size=%d idx_node=%d\n", node, (int) gcpo[node]->cov_mat->size1,
					       (int) gcpo[node]->idx_node);
					GMRFLib_idxval_printf(stdout, gcpo[node]->idxs, GMRFLib_strdup("nodes in this group"));
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

	GMRFLib_idx_tp *node_idx2 = NULL;
	for (int node = 0; node < Npred; node++) {
		if (gcpo[node]->cov_mat) {
			GMRFLib_idx_add(&node_idx2, node);
		} else {
			gcpo[node]->value = NAN;
		}
	}

#define CODE_BLOCK							\
	for(int inode = 0; inode < node_idx2->n; inode++) {		\
		int node = node_idx2->idx[inode];			\
									\
		CODE_BLOCK_ALL_WORK_ZERO();				\
		int ng = gcpo[node]->idxs->n;				\
		int *idxs = gcpo[node]->idxs->idx;			\
		size_t idx_node = gcpo[node]->idx_node;			\
		double bb_idx_node = NAN, cc_idx_node = NAN;		\
		double *bb = CODE_BLOCK_WORK_PTR(0);			\
		double *cc = CODE_BLOCK_WORK_PTR(1);			\
		gsl_vector *mean_old = gsl_vector_calloc((size_t) ng);	\
		gsl_vector *mean = gsl_vector_calloc((size_t) ng);	\
		gsl_vector *b = gsl_vector_calloc((size_t) ng);		\
		gsl_matrix *S = GMRFLib_gsl_duplicate_matrix(gcpo[node]->cov_mat); \
		gsl_matrix *Q = S;					\
		if (detailed_output) {					\
			printf("node %d, idx_node %zu,cov mat\n", node, idx_node); \
			GMRFLib_printf_gsl_matrix(stdout, S, " %.8f ");	\
		}							\
									\
		gcpo[node]->marg_theta_correction = 0.0;		\
		for(int i = 0; i < ng; i++) {				\
			int nnode = idxs[i];				\
			gsl_vector_set(mean_old, (size_t) i, lpred_mode[nnode]); \
			double ll = 0.0, local_bb = 0.0, local_cc = 0.0; \
			if (d[nnode]) {					\
				if (corr_hypar) loglFunc(thread_id, &ll, &(lpred_mode[nnode]), 1, nnode, lpred_mode, NULL, loglFunc_arg, NULL); \
				GMRFLib_2order_approx(thread_id, NULL, &local_bb, &local_cc, NULL, d[nnode], lpred_mode[nnode], nnode, \
						      lpred_mode, loglFunc, loglFunc_arg, &ai_par->step_len, &ai_par->stencil, &zero); \
			}						\
			if (corr_hypar) {				\
				gcpo[node]->marg_theta_correction += ll; \
			}						\
			bb[i] += local_bb;				\
			cc[i] += local_cc;				\
			if (i == (int) idx_node) {			\
				bb_idx_node = local_bb;			\
				cc_idx_node = local_cc;			\
			}						\
		}							\
									\
		if (detailed_output) {					\
			printf("node %d, cov.mat and mean\n", node);	\
			GMRFLib_printf_gsl_matrix(stdout, Q, " %.8f ");	\
			GMRFLib_printf_gsl_vector(stdout, mean_old, " %.8f "); \
		}							\
									\
		if (1) { /* the new low-rank solution... */		\
			double low_rank_eps = GMRFLib_eps(0.3833);	\
			/* new low-rank approach */			\
			size_t n = (size_t) ng;				\
			gsl_matrix *Cov = S;				\
			gsl_matrix *B = GMRFLib_gsl_low_rank(Cov, low_rank_eps); \
			gsl_matrix *Bt = GMRFLib_gsl_transpose_matrix(B); \
			size_t m = B->size2;				\
			gsl_matrix *H = gsl_matrix_alloc(n, n);		\
			gsl_matrix_set_zero(H);				\
			for(size_t i = 0; i < n; i++) {			\
				gsl_matrix_set(H, i, i, cc[i]);		\
			}						\
									\
			gsl_vector *ztmp = gsl_vector_alloc(n);		\
			GMRFLib_gsl_mv(H, mean_old, ztmp);		\
			for(size_t i = 0; i < n; i++) {			\
				gsl_vector_set(ztmp, i, gsl_vector_get(ztmp, i) - bb[i]); \
			}						\
			gsl_vector *zb = gsl_vector_alloc(m);		\
			GMRFLib_gsl_mv(Bt, ztmp, zb);			\
									\
			gsl_matrix *BtHB = gsl_matrix_calloc(m, m);	\
			GMRFLib_gsl_mmm(Bt, H, B, BtHB);		\
			gsl_matrix *QQ = gsl_matrix_alloc(m, m);	\
			gsl_matrix_set_zero(QQ);			\
			if (corr_hypar) {				\
				gcpo[node]->marg_theta_correction -= GMRFLib_gsl_log_dnorm(NULL, NULL, QQ, NULL, 1); \
			}						\
			for(size_t i = 0; i < m; i++) {			\
				for(size_t j = 0; j < m; j++) {		\
					if (i == j) {			\
						gsl_matrix_set(QQ, i, j, DMAX(0, 1.0 - gsl_matrix_get(BtHB, i, j))); \
					} else {			\
						gsl_matrix_set(QQ, i, j, - gsl_matrix_get(QQ, i, j)); \
					}				\
				}					\
			}						\
			/* same ptr */					\
			gsl_matrix *SS = QQ;				\
			GMRFLib_gsl_ensure_spd_inverse(SS, low_rank_eps, NULL); \
									\
			gsl_vector *zmean = gsl_vector_alloc(m);	\
			GMRFLib_gsl_mv(SS, zb, zmean);			\
									\
			if (detailed_output) {				\
				printf("node %d, cov.mat and mean after correction (low-rank)\n", node); \
				GMRFLib_printf_gsl_matrix(stdout, SS, " %.8f "); \
				GMRFLib_printf_gsl_vector(stdout, zmean, " %.8f "); \
			}						\
									\
			gsl_vector *Bzmean = gsl_vector_alloc(n);	\
			GMRFLib_gsl_mv(B, zmean, Bzmean);		\
			for(size_t i = 0; i < n; i++) {			\
				gsl_vector_set(mean, i, gsl_vector_get(mean_old, i) + gsl_vector_get(Bzmean, i)); \
			}						\
			GMRFLib_gsl_mmm(B, SS, Bt, S);			\
									\
			if (corr_hypar) {				\
				gcpo[node]->marg_theta_correction += GMRFLib_gsl_log_dnorm(NULL, zmean, NULL, SS, 0); \
				/* we define the correction to be multiplicative */ \
				gcpo[node]->marg_theta_correction *= -1.0; \
			}						\
									\
			gsl_matrix_free(B);				\
			gsl_matrix_free(Bt);				\
			gsl_matrix_free(H);				\
			gsl_matrix_free(BtHB);				\
			gsl_matrix_free(QQ);				\
			gsl_vector_free(ztmp);				\
			gsl_vector_free(zb);				\
			gsl_vector_free(zmean);				\
			gsl_vector_free(Bzmean);			\
									\
		} else {						\
									\
			for(size_t i = 0; i < (size_t) ng; i++) {	\
				gsl_matrix_set(Q, i, i, gsl_matrix_get(Q, i, i) * diag_scale); \
			}						\
			GMRFLib_gsl_ensure_spd_inverse(Q, spd_eps, NULL); \
			GMRFLib_gsl_mv(Q, mean_old, b);			\
			if (corr_hypar) {				\
				gcpo[node]->marg_theta_correction -= GMRFLib_gsl_log_dnorm(NULL, NULL, Q, NULL, 0); \
			}						\
									\
			for(size_t i = 0; i < (size_t) ng; i++) {	\
				gsl_matrix_set(Q, i, i, DMAX(0.0, gsl_matrix_get(Q, i, i) - cc[i])); \
				gsl_vector_set(b, i, gsl_vector_get(b, i) - bb[i]); \
			}						\
									\
			if (detailed_output) {				\
				printf("node %d, prec mat and b after correction\n", node); \
				GMRFLib_printf_gsl_matrix(stdout, Q, " %.8f ");	\
				GMRFLib_printf_gsl_vector(stdout, b, " %.8f "); \
			}						\
									\
			GMRFLib_gsl_ensure_spd_inverse(Q, spd_eps, NULL); \
			GMRFLib_gsl_mv(S, b, mean);			\
			if (corr_hypar) {				\
				gcpo[node]->marg_theta_correction += GMRFLib_gsl_log_dnorm(mean_old, mean, NULL, S, 0); \
				/* we define the correction to be multiplicative */ \
				gcpo[node]->marg_theta_correction *= -1.0; \
			}						\
									\
			if (detailed_output) {				\
				printf("node %d, new cov mat and mean\n", node); \
				GMRFLib_printf_gsl_matrix(stdout, S, " %.8f ");	\
				GMRFLib_printf_gsl_vector(stdout, mean, " %.8f "); \
			}						\
									\
		}							\
		gsl_vector_set(mean, idx_node, gsl_vector_get(mean, idx_node) + (lpred_mean[node] - lpred_mode[node])); \
		gcpo[node]->lpred_mean = gsl_vector_get(mean, idx_node); \
		gcpo[node]->lpred_sd = sqrt(DMAX(DBL_EPSILON, gsl_matrix_get(S, idx_node, idx_node) / (1 + 0.0 * diag_scale))); \
		gcpo[node]->kld =  0.5 * (SQR(gcpo[node]->lpred_sd) / lpred_variance[node] - 1.0 + \
					  SQR(gcpo[node]->lpred_mean - lpred_mean[node]) / lpred_variance[node] + \
					  log(lpred_variance[node] / SQR(gcpo[node]->lpred_sd))); \
									\
		if (gcpodens_moments) {					\
			gcpodens_moments[node * 3 + 0] = gcpo[node]->lpred_mean; \
			gcpodens_moments[node * 3 + 1] = SQR(gcpo[node]->lpred_sd); \
			gcpodens_moments[node * 3 + 2] = gcpo[node]->marg_theta_correction; \
		}							\
									\
		if (d[node]) {						\
			/* do the integral by approximating phi(x)*exp(ll(x)) with a normal, and */ \
			/* then do the GHQ with respect to that normal as the kernel, with the */ \
			/* ``errors'' as the function */		\
			double *weights = NULL, *xx = NULL;		\
			GMRFLib_ghq(&xx, &weights, np);			\
									\
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
		} else {						\
			gcpo[node]->value = NAN;			\
		}							\
									\
		if (gcpo_param->verbose || detailed_output) {		\
			printf("%s[%1d]: node %d lpred_mean %f lpred_sd %f kld %f value %f\n", \
			       __GMRFLib_FuncName, omp_get_thread_num(), \
			       node, gcpo[node]->lpred_mean, gcpo[node]->lpred_sd, gcpo[node]->kld, gcpo[node]->value); \
		}							\
									\
		gsl_vector_free(mean);					\
		gsl_vector_free(mean_old);				\
		gsl_vector_free(b);					\
		gsl_matrix_free(Q);					\
	}

	RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 4, IMAX(np, max_ng));
#undef CODE_BLOCK

	GMRFLib_idx_free(node_idx2);
	Calloc_free();
#undef A_idx

	GMRFLib_LEAVE_ROUTINE;
	return gcpo;
}

int GMRFLib_compute_cpodens(int thread_id, GMRFLib_density_tp ** cpo_density, GMRFLib_density_tp * density,
			    int idx, double d, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, GMRFLib_ai_param_tp * ai_par)
{
	if (!d) {
		*cpo_density = NULL;
		return GMRFLib_SUCCESS;
	}

	const int debug = 0;
	int itry, flag, np, np_orig = GMRFLib_INT_GHQ_POINTS + 4, i, npx = 8, itmp, np_new = np_orig + 2 * npx, one = 1;
	double *xp = NULL, *xp_tmp = NULL, *ld = NULL, *logcor = NULL, *x_user = NULL, alpha = -1.0;
	double cor_eps = GMRFLib_eps(0.75), cor_max, range;

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
#pragma GCC ivdep
		for (i = 0; i < np; i++) {
			logcor[i] *= d;
		}
		if (debug && np) {
#pragma omp critical (Name_45542d32821a8fbfd2cec71e8219d7eeb4b423f2)
			{
				for (i = 0; i < np; i++)
					printf("CPO: %d BEFORE x_user %g xp %g ld %g logcor %g ld-logcor %g\n", idx,
					       x_user[i], xp[i], ld[i], logcor[i], ld[i] - logcor[i]);
			}
		}
		if (itry == 1 && cor_eps > 0.0) {
			flag = 1;
			cor_max = exp(log(cor_eps) + GMRFLib_max_value(logcor, np, NULL));
#pragma GCC ivdep
			for (i = 0; i < np; i++) {
				ld[i] += logcor[i] - 2.0 * GMRFLib_log_apbex(cor_max, logcor[i]);
			}
		} else {
			daxpy_(&np, &alpha, logcor, &one, ld, &one);	/* ld = ld + logcor */
		}
		GMRFLib_ai_correct_cpodens(ld, xp, &np, ai_par);
		if (debug && np) {
#pragma omp critical (Name_c6e59ebf504f17645e98f57731cc4de48bd2748a)
			{
				for (i = 0; i < np; i++)
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

int GMRFLib_ai_vb_prepare(int thread_id,
			  GMRFLib_vb_coofs_tp * coofs, int idx, GMRFLib_density_tp * density, double d, GMRFLib_logl_tp * loglFunc,
			  void *loglFunc_arg, double *x_vec)
{
	/*
	 * compute the Taylor-expansion of -loglikelihood * density(x), around the mean of x
	 */

	// Normal kernel: deriv: ... * (x-m)/s^2
	// dderiv: ... * ((x-m)^2 - s^2)/s^4
	// GMRFLib_density_type_tp type;
	/*
	 * params for the GMRFLib_DENSITY_TYPE_GAUSSIAN 
	 * params for the GMRFLib_DENSITY_TYPE_SKEWNORMAL 
	 * params for the GMRFLib_DENSITY_TYPE_SCGAUSSIAN 
	 */

	if (ISZERO(d)) {
		coofs->coofs[0] = coofs->coofs[1] = coofs->coofs[2] = 0.0;
		return GMRFLib_SUCCESS;
	}

	if (density->type == GMRFLib_DENSITY_TYPE_GAUSSIAN) {
		// life is simpler in this case

		const int np = GMRFLib_INT_GHQ_POINTS;
		double m = density->user_mean;
		double s = density->user_stdev;

		Calloc_init(3 * np, 3);
		double *xp = NULL;
		double *wp = NULL;
		double *x_user = Calloc_get(np);
		double *x_std = Calloc_get(np);
		double *loglik = Calloc_get(np);

		GMRFLib_ghq(&xp, &wp, np);		       /* just give ptr to storage */

#pragma GCC ivdep
		for (int i = 0; i < np; i++) {
			x_user[i] = m + s * xp[i];
		}
		GMRFLib_density_user2std_n(x_std, x_user, density, np);
		loglFunc(thread_id, loglik, x_user, np, idx, x_vec, NULL, loglFunc_arg, NULL);

		double A = 0.0, B = 0.0, C = 0.0, s_inv = 1.0 / s, s2_inv = 1.0 / SQR(s);
		for (int i = 0; i < np; i++) {
			double tmp = wp[i] * loglik[i];
			A -= tmp;
			B -= tmp * xp[i];
			C -= tmp * (SQR(xp[i]) - 1.0);
		}
		coofs->coofs[0] = d * A;
		coofs->coofs[1] = d * B * s_inv;
		coofs->coofs[2] = d * C * s2_inv;

		Calloc_free();
		return GMRFLib_SUCCESS;
	} else {
		int i, k, np = GMRFLib_INT_NUM_POINTS;
		double low, dx, dxi, *xp = NULL, *xpi = NULL, *dens = NULL, *loglik = NULL, w[2] =
		    { 4.0, 2.0 }, integral_one, integral, integral_p, integral_m;

		Calloc_init(4 * np, 4);
		xp = Calloc_get(np);
		xpi = Calloc_get(np);
		dens = Calloc_get(np);
		loglik = Calloc_get(np);

		dxi = (density->x_max - density->x_min) / (np - 1.0);
		low = GMRFLib_density_std2user(density->x_min, density);
		dx = (GMRFLib_density_std2user(density->x_max, density) - low) / (np - 1.0);

		xp[0] = low;
		xpi[0] = density->x_min;
#pragma GCC ivdep
		for (i = 1; i < np; i++) {
			xp[i] = xp[0] + i * dx;
			xpi[i] = xpi[0] + i * dxi;
		}

		GMRFLib_evaluate_ndensity(dens, xpi, np, density);
		loglFunc(thread_id, loglik, xp, np, idx, x_vec, NULL, loglFunc_arg, NULL);
#pragma GCC ivdep
		for (i = 0; i < np; i++) {
			loglik[i] *= d;
		}

		// those zero'ed out are included in the loop
		integral_one = dens[0] + dens[np - 1];
		integral = dens[0] * loglik[0] + dens[np - 1] * loglik[np - 1];
		integral_p = dens[1] * loglik[0] + 0 * dens[np - 1] * loglik[np - 2];
		integral_m = 0 * dens[0] * loglik[1] + dens[np - 2] * loglik[np - 1];
		for (i = 1, k = 0; i < np - 1; i++, k = (k + 1) % 2) {
			integral_one += w[k] * dens[i];
			integral += w[k] * loglik[i] * dens[i];
			integral_p += w[k] * loglik[i] * dens[i - 1];
			integral_m += w[k] * loglik[i] * dens[i + 1];
		}
		integral /= (-integral_one);		       /* -E() */
		integral_p /= (-integral_one);
		integral_m /= (-integral_one);

		// c[0] + c[1]*x + 1/2*c[2]*x^2 + c[3]*y + 1/2*c[4]*y2 + c[5]*x*y
		coofs->coofs[0] = integral;
		coofs->coofs[1] = (integral_p - integral_m) / (2.0 * dx);
		coofs->coofs[2] = (integral_p - 2.0 * integral + integral_m) / SQR(dx);

		Calloc_free();
		return GMRFLib_SUCCESS;
	}
}

int GMRFLib_ai_vb_prepare_mean(int thread_id,
			       GMRFLib_vb_coofs_tp * coofs, int idx, double d, GMRFLib_logl_tp * loglFunc,
			       void *loglFunc_arg, double *x_vec, double mean, double sd)
{
	/*
	 * compute the Taylor-expansion of integral of -loglikelihood * density(x), around the mean of x
	 */

	// Normal kernel: deriv: ... * (x-m)/s^2
	// dderiv: ... * ((x-m)^2 - s^2)/s^4
	// GMRFLib_density_type_tp type;

	if (ISZERO(d)) {
		coofs->coofs[0] = coofs->coofs[1] = coofs->coofs[2] = 0.0;
		return GMRFLib_SUCCESS;
	}

	static double *xp = NULL;
	static double *wp = NULL;
	static double *xp2 = NULL;

	if (!wp) {
#pragma omp critical (Name_00c5c0bab9ee4213c2351e3b2275ded2f8b87d22)
		{
			if (!wp) {
				double *wtmp = NULL;
				GMRFLib_ghq(&xp, &wtmp, GMRFLib_INT_GHQ_POINTS);	/* just give ptr to storage */
				xp2 = Calloc(GMRFLib_INT_GHQ_POINTS, double);
				for (int i = 0; i < GMRFLib_INT_GHQ_POINTS; i++) {
					xp2[i] = SQR(xp[i]) - 1.0;
				}
				wp = wtmp;
			}
		}
	}

	double x_user[2 * GMRFLib_INT_GHQ_POINTS];
	double *loglik = x_user + GMRFLib_INT_GHQ_POINTS;

#pragma GCC ivdep
	for (int i = 0; i < GMRFLib_INT_GHQ_POINTS; i++) {
		x_user[i] = mean + sd * xp[i];
	}
	loglFunc(thread_id, loglik, x_user, GMRFLib_INT_GHQ_POINTS, idx, x_vec, NULL, loglFunc_arg, NULL);

	double A, B, C, s_inv = 1.0 / sd, s2_inv = 1.0 / SQR(sd), tmp;

	// optimized version. since xp and wp are symmetric and xp[idx]=0
	int ni = GMRFLib_INT_GHQ_POINTS / 2L;
	tmp = wp[ni] * loglik[ni];
	A = tmp;
	B = 0.0;
	C = -tmp;
	for (int i = 0; i < ni; i++) {
		int ii = GMRFLib_INT_GHQ_POINTS - 1 - i;
		double tt = wp[i] * (loglik[i] + loglik[ii]);
		double tt2 = wp[i] * (loglik[i] - loglik[ii]);
		A += tt;
		C += tt * xp2[i];
		B += tt2 * xp[i];
	}
	coofs->coofs[0] = - d * A;
	coofs->coofs[1] = - d * B * s_inv;
	coofs->coofs[2] = - d * C * s2_inv;

	return GMRFLib_SUCCESS;
}

int GMRFLib_ai_vb_prepare_variance(int thread_id, GMRFLib_vb_coofs_tp * coofs, int idx, double d,
				   GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, double *x_vec, double mean, double sd)
{
	/*
	 * compute the Taylor-expansion of the integral of -loglikelihood * density(x), for the variance of x (assuming its N())
	 */

	if (ISZERO(d)) {
		coofs->coofs[0] = coofs->coofs[1] = coofs->coofs[2] = 0.0;
		return GMRFLib_SUCCESS;
	}

	static double *wp = NULL;
	static double *xp = NULL;
	static double *xp2 = NULL;
	static double *xp3 = NULL;

	if (!wp) {
#pragma omp critical (Name_0ddd01862f572e8e2021d8c931021738790dccc7)
		{
			if (!wp) {
				double *wtmp = NULL;
				GMRFLib_ghq(&xp, &wtmp, GMRFLib_INT_GHQ_POINTS);	/* just give ptr to storage */
				xp2 = Calloc(2 * GMRFLib_INT_GHQ_POINTS, double);
				xp3 = xp2 + GMRFLib_INT_GHQ_POINTS;
				for (int i = 0; i < GMRFLib_INT_GHQ_POINTS; i++) {
					double z2 = SQR(xp[i]);
					xp2[i] = z2 - 1.0;
					xp3[i] = 3.0 - 6.0 * z2 + SQR(z2);
				}
				wp = wtmp;
			}
		}
	}

	double x_user[2 * GMRFLib_INT_GHQ_POINTS];
	double *loglik = x_user + GMRFLib_INT_GHQ_POINTS;

#pragma GCC ivdep
	for (int i = 0; i < GMRFLib_INT_GHQ_POINTS; i++) {
		x_user[i] = mean + sd * xp[i];
	}
	loglFunc(thread_id, loglik, x_user, GMRFLib_INT_GHQ_POINTS, idx, x_vec, NULL, loglFunc_arg, NULL);

	double A, B, C, s2_inv = 1.0 / SQR(sd);

	// optimized version, as both xp and wp are symmetric and xp[idx]=0
	int ni = GMRFLib_INT_GHQ_POINTS / 2L;
	double tmp = wp[ni] * loglik[ni];
	A = tmp;
	B = -tmp;
	C = 3.0 * tmp;
	for (int i = 0; i < ni; i++) {
		int ii =  GMRFLib_INT_GHQ_POINTS - 1 - i;
		double tt = wp[i] * (loglik[i] + loglik[ii]);
		A += tt;
		B += tt * xp2[i];
		C += tt * xp3[i];
	}

	coofs->coofs[0] = - d * A;
	coofs->coofs[1] = - d * B * 0.5 * s2_inv;
	coofs->coofs[2] = - d * C * 0.25 * s2_inv;

	return GMRFLib_SUCCESS;
}

int GMRFLib_vb_fit_gaussian(int n, double *x, double *ld, double *mean, double *sd)
{
	// do a quick fit of a sequence of (x,ld) to get an approximate mean and sd.
	// assume 'x' is sorted

	int imax = -1;
	GMRFLib_max_value(ld, n, &imax);
	assert(imax >= 0);

	int istart = IMAX(0, imax - 2);
	int nn = IMIN(n - 1, imax + 2) - istart + 1;

	gsl_matrix *X = gsl_matrix_alloc(nn, 3);
	gsl_matrix *cov = gsl_matrix_alloc(3, 3);
	gsl_vector *c = gsl_vector_alloc(3);
	gsl_vector *y = gsl_vector_alloc(nn);

	for (int i = 0; i < nn; i++) {
		double xi = x[istart + i];
		gsl_matrix_set(X, i, 0, 1.0);
		gsl_matrix_set(X, i, 1, xi);
		gsl_matrix_set(X, i, 2, -0.5 * SQR(xi));
		gsl_vector_set(y, i, ld[istart + i]);
	}

	double chisq = 0.0;
	gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(nn, 3);
	gsl_multifit_linear(X, y, c, cov, &chisq, work);

	double c2 = gsl_vector_get(c, 2);
	if (c2 > 0.0) {
		*mean = gsl_vector_get(c, 1) / c2;
		*sd = sqrt(1.0 / c2);
	} else {
		*mean = NAN;
		*sd = NAN;
	}

	gsl_multifit_linear_free(work);
	gsl_matrix_free(X);
	gsl_vector_free(y);
	gsl_vector_free(c);
	gsl_matrix_free(cov);

	return GMRFLib_SUCCESS;
}
int GMRFLib_ai_vb_correct_mean(int thread_id, GMRFLib_density_tp *** density,	// need two types
			       int dens_count,
			       GMRFLib_density_tp ** dens_local,
			       double *c,
			       double *d,
			       GMRFLib_ai_param_tp * ai_par,
			       GMRFLib_ai_store_tp * ai_store,
			       GMRFLib_graph_tp * graph,
			       GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg,
			       GMRFLib_preopt_tp * UNUSED(preopt))
{
	if (GMRFLib_inla_mode == GMRFLib_MODE_TWOSTAGE_PART1 || GMRFLib_inla_mode == GMRFLib_MODE_COMPACT) {
		// nothing to do here
		return GMRFLib_SUCCESS;
	} else {
		return GMRFLib_ai_vb_correct_mean_std(thread_id, density, dens_count, dens_local,
						      c, d, ai_par, ai_store, graph, Qfunc, Qfunc_arg, loglFunc, loglFunc_arg);
	}
}

int GMRFLib_ai_vb_correct_mean_std(int thread_id, GMRFLib_density_tp *** density,	// need two types
				   int dens_count,
				   GMRFLib_density_tp ** dens_local,
				   double *c,
				   double *d,
				   GMRFLib_ai_param_tp * ai_par,
				   GMRFLib_ai_store_tp * ai_store,
				   GMRFLib_graph_tp * graph,
				   GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg)
{
	// this could be improved similarly to _preopt

	double one = 1.0, mone = -1.0, zero = 0.0;
	double max_correct = 1.0;

	if (!(ai_par->vb_enable && ai_par->vb_nodes_mean))
		return GMRFLib_SUCCESS;

	// need the idx's for the vb correction and the data locations
	GMRFLib_idx_tp *vb_idx = NULL, *d_idx = NULL;

	GMRFLib_idx_create(&vb_idx);
	assert(vb_idx);
	for (int i = 0; i < graph->n; i++) {
		if (ai_par->vb_nodes_mean[i]) {
			GMRFLib_idx_add(&vb_idx, i);
		}
	}
	for (int i = 0; i < graph->n; i++) {
		if (d[i]) {
			GMRFLib_idx_add(&d_idx, i);
		}
	}

	if (vb_idx->n > 0) {

		double *sd = Calloc(graph->n, double);
		double tref = GMRFLib_cpu();
		double *mode = Calloc(graph->n, double);
		GMRFLib_vb_coofs_tp **vb_coof = Calloc(graph->n, GMRFLib_vb_coofs_tp *);
		for (int i = 0; i < graph->n; i++) {
			vb_coof[i] = Calloc(1, GMRFLib_vb_coofs_tp);
		}

		GMRFLib_ai_add_Qinv_to_ai_store(ai_store);
		for (int i = 0; i < graph->n; i++) {
			double *var = GMRFLib_Qinv_get(ai_store->problem, i, i);
			sd[i] = (var ? sqrt(*var) : NAN);
		}

		Memcpy(mode, ai_store->mode, graph->n * sizeof(double));
		for (int i = 0; i < graph->n; i++) {
			if (density) {
				if (density[i][dens_count]) {
					mode[i] = density[i][dens_count]->user_mean;
				}
			} else {
				if (dens_local && dens_local[i]) {
					mode[i] = dens_local[i]->user_mean;
				}
			}
		}

#define CODE_BLOCK						\
		for (int ii = 0; ii < d_idx->n; ii++) {		\
			int i = d_idx->idx[ii];			\
			if (density) {					\
				GMRFLib_ai_vb_prepare(thread_id, vb_coof[i], i, density[i][dens_count], d[i], loglFunc, loglFunc_arg, mode); \
			} else {					\
				GMRFLib_ai_vb_prepare(thread_id, vb_coof[i], i, dens_local[i], d[i], loglFunc, loglFunc_arg, mode); \
			}						\
		}

		assert(d_idx);
		RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 0, 0);
#undef CODE_BLOCK

		double *c_diag = Calloc(graph->n, double);
		double *cmean = Calloc(graph->n, double);
		double *corr = Calloc(graph->n, double);
		double *tmp = Calloc(graph->n, double);
		gsl_matrix *QM = gsl_matrix_alloc(graph->n, vb_idx->n);
		gsl_matrix *M = gsl_matrix_alloc(graph->n, vb_idx->n);	// matrix with Cov()
		gsl_vector *B = gsl_vector_alloc(graph->n);

		gsl_matrix_set_zero(QM);
		gsl_matrix_set_zero(M);
		gsl_vector_set_zero(B);

		for (int jj = 0; jj < vb_idx->n; jj++) {
			int j = vb_idx->idx[jj];
			GMRFLib_ai_update_conditional_mean2(cmean, ai_store->problem, j, ai_store->problem->mean_constr[j] + 1.0, NULL);
			for (int i = 0; i < graph->n; i++) {
				// need correlation identical to 1 for i=j
				corr[i] = (i == j ? 1.0 : sd[j] * (cmean[i] - ai_store->problem->mean_constr[i]) / sd[i]);
				gsl_matrix_set(M, i, jj, corr[i] * sd[i] * sd[j]);	/* yes, it is 'jj' and not 'j' */
			}
		}

		GMRFLib_Qx2(thread_id, tmp, mode, graph, Qfunc, Qfunc_arg, c);	/* mode=mean, the point we expand around */
		for (int i = 0; i < graph->n; i++) {
			gsl_vector_set(B, i, tmp[i]);
		}
		for (int ii = 0; ii < d_idx->n; ii++) {
			int i = d_idx->idx[ii];
			gsl_vector_set(B, i, vb_coof[i]->coofs[1] + gsl_vector_get(B, i));
			c_diag[i] = vb_coof[i]->coofs[2];
		}

#pragma GCC ivdep
		for (int i = 0; i < graph->n; i++) {
			c_diag[i] += c[i] + 1E-6 * (d[i] ? ai_store->cc[i] : 1.0);
		}

#define CODE_BLOCK							\
		for (int j = 0; j < vb_idx->n; j++) {			\
			double *col = CODE_BLOCK_WORK_PTR(0);		\
			double *res = CODE_BLOCK_WORK_PTR(1);		\
			assert(col);					\
			assert(res);					\
			for (int i = 0; i < graph->n; i++) {		\
				col[i] = gsl_matrix_get(M, i, j);	\
			}						\
			GMRFLib_Qx2(thread_id, res, col, graph, Qfunc, Qfunc_arg, c_diag); \
			for (int i = 0; i < graph->n; i++) {		\
				gsl_matrix_set(QM, i, j, res[i]);	\
			}						\
		}

		RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 2, graph->n);
#undef CODE_BLOCK

		gsl_matrix *MM = gsl_matrix_alloc(vb_idx->n, vb_idx->n);
		gsl_permutation *perm = gsl_permutation_alloc(vb_idx->n);
		gsl_vector *MB = gsl_vector_alloc(vb_idx->n);
		gsl_vector *delta = gsl_vector_alloc(vb_idx->n);
		gsl_vector *delta_mu = gsl_vector_alloc(graph->n);

		gsl_matrix_set_zero(MM);
		gsl_vector_set_zero(MB);
		gsl_vector_set_zero(delta);
		gsl_vector_set_zero(delta_mu);

		gsl_blas_dgemm(CblasTrans, CblasNoTrans, one, M, QM, zero, MM);
		gsl_blas_dgemv(CblasTrans, mone, M, B, zero, MB);

		// need pivoting to solve the system
		GMRFLib_gsl_ensure_spd(MM, GMRFLib_eps(0.5), NULL);
		gsl_linalg_pcholesky_decomp(MM, perm);
		gsl_linalg_pcholesky_solve(MM, perm, MB, delta);
		gsl_blas_dgemv(CblasNoTrans, one, M, delta, zero, delta_mu);

		int num_trunc = 0;
		for (int i = 0; i < graph->n; i++) {
			if (ABS(gsl_vector_get(delta_mu, i) / sd[i]) > max_correct) {
				gsl_vector_set(delta_mu, i, max_correct * sd[i] * SIGN(gsl_vector_get(delta_mu, i)));
				num_trunc++;
			}
		}

		if (density) {
			for (int i = 0; i < graph->n; i++) {
				if (density[i][dens_count]) {
					GMRFLib_density_new_user_mean(density[i][dens_count],
								      density[i][dens_count]->user_mean + gsl_vector_get(delta_mu, i));
				}
			}
		} else {
			for (int i = 0; i < graph->n; i++) {
				if (dens_local && dens_local[i]) {
					GMRFLib_density_new_user_mean(dens_local[i], dens_local[i]->user_mean + gsl_vector_get(delta_mu, i));
				}
			}
		}

		if (ai_par->vb_verbose) {
			printf("\tVB correct with strategy [mean] in [%.3f]seconds\n", GMRFLib_cpu() - tref);
			printf("\t\tNumber of nodes corrected for [%1d]\n", (int) delta->size);
			printf("\t\tIter.max [%1d]\n", (int) ai_par->vb_iter_max);
			printf("\t\tNumber of corrections truncated [%1d] with max.correct[%.2f]\n", num_trunc, max_correct);
			for (int jj = 0; jj < vb_idx->n; jj++) {
				int j = vb_idx->idx[jj];
				printf("\t\tNode[%1d] delta[%.3f] correction[%.3f] correction/stdev[%.3f]\n",
				       j, gsl_vector_get(delta_mu, j), gsl_vector_get(delta_mu, j), gsl_vector_get(delta_mu, j) / sd[j]);
			}
			printf("\t\tImplied correction for [%1d] nodes\n", graph->n - vb_idx->n);
		}

		// clean-up
		for (int i = 0; i < graph->n; i++) {
			if (vb_coof[i]) {
				Free(vb_coof[i]);
			}
		}
		Free(vb_coof);

		Free(c_diag);
		Free(cmean);
		Free(corr);
		Free(mode);
		Free(sd);
		Free(tmp);

		gsl_matrix_free(M);
		gsl_matrix_free(MM);
		gsl_matrix_free(QM);
		gsl_permutation_free(perm);
		gsl_vector_free(B);
		gsl_vector_free(MB);
		gsl_vector_free(delta);
		gsl_vector_free(delta_mu);
	}

	GMRFLib_idx_free(d_idx);
	GMRFLib_idx_free(vb_idx);

	return GMRFLib_SUCCESS;
}

int GMRFLib_ai_vb_correct_mean_preopt(int thread_id,
				      GMRFLib_density_tp *** density,
				      int dens_count,
				      double *UNUSED(c),
				      double *d,
				      GMRFLib_ai_param_tp * ai_par,
				      GMRFLib_ai_store_tp * ai_store,
				      GMRFLib_graph_tp * graph,
				      GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg,
				      GMRFLib_preopt_tp * preopt)
{
	GMRFLib_ENTER_ROUTINE;
	Calloc_init(5 * graph->n + 2 * preopt->mnpred + 4 * preopt->Npred, 11);
	FILE *fp = (ai_par->fp_log ? ai_par->fp_log : stdout);
	int verbose = ai_par->vb_verbose && ai_par->fp_log;

#define SHOW_TIME(msg_)							\
	if (debug) {							\
		fprintf(fp, "[%1d] vb_preopt: %s %.3f\n", omp_get_thread_num(), msg_, GMRFLib_cpu()-_tref); \
		_tref = GMRFLib_cpu();					\
	}

	assert(GMRFLib_inla_mode == GMRFLib_MODE_COMPACT);

	// save time: only compute MM the first time, and keep MM and its factorisation fixed during the iterations. the motivation is that the
	// 2nd order properties will hardly change while the 1st order properties, ie the mean, will.
	int keep_MM = 1, ratio_ok = 0;

	int niter = ai_par->vb_iter_max;
	int debug = GMRFLib_DEBUG_IF();
	int emergency = 0;
	double one = 1.0, mone = -1.0, zero = 0.0;
	double max_correct = 5.001;
	double tref = GMRFLib_cpu();
	GMRFLib_tabulate_Qfunc_tp *tabQ = NULL;

	if (!(ai_par->vb_enable && ai_par->vb_nodes_mean)) {
		GMRFLib_LEAVE_ROUTINE;
		return GMRFLib_SUCCESS;
	}

	// need the idx's for the vb correction and the data locations
	GMRFLib_idx_tp *vb_idx = NULL, *d_idx = NULL;

	for (int i = 0; i < graph->n; i++) {
		if (ai_par->vb_nodes_mean[i]) {
			GMRFLib_idx_add(&vb_idx, i);
		}
	}
	if (!vb_idx) {
		GMRFLib_LEAVE_ROUTINE;
		return GMRFLib_SUCCESS;
	}

	for (int i = 0; i < preopt->Npred; i++) {
		if (d[i]) {
			GMRFLib_idx_add(&d_idx, i);
		}
	}
	assert(d_idx);

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
	gsl_matrix *QM = gsl_matrix_alloc(graph->n, vb_idx->n);
	gsl_matrix *M = gsl_matrix_alloc(graph->n, vb_idx->n); // matrix with Cov()
	gsl_vector *B = gsl_vector_alloc(graph->n);
	gsl_matrix *MM = gsl_matrix_alloc(vb_idx->n, vb_idx->n);
	gsl_permutation *perm = gsl_permutation_alloc(vb_idx->n);
	gsl_vector *MB = gsl_vector_alloc(vb_idx->n);
	gsl_vector *delta = gsl_vector_alloc(vb_idx->n);
	gsl_vector *delta_mu = gsl_vector_alloc(graph->n);

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

	GMRFLib_tabulate_Qfunc_core(thread_id, &tabQ, graph, Qfunc, Qfunc_arg, NULL, 1);
	GMRFLib_tabulate_Qfunc_core(thread_id, &prior, preopt->latent_graph, GMRFLib_preopt_Qfunc_prior, Qfunc_arg, NULL, 1);
	gsl_matrix_set_zero(M);
	gsl_matrix_set_zero(QM);

#define CODE_BLOCK							\
	for (int jj = 0; jj < vb_idx->n; jj++) {			\
		int j = vb_idx->idx[jj];				\
		double *b = CODE_BLOCK_WORK_PTR(0);			\
		double *cov = CODE_BLOCK_WORK_PTR(1);			\
		CODE_BLOCK_WORK_ZERO(0);				\
		b[j] = 1.0;						\
		GMRFLib_Qsolve(cov, b, ai_store->problem, j);		\
		for (int i = 0; i < graph->n; i++) {			\
			gsl_matrix_set(M, i, jj, cov[i]);		\
		}							\
	}

	RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 2, graph->n);
#undef CODE_BLOCK

	for (int iter = 0; iter < niter; iter++) {
		int update_MM = ((iter == 0) || !keep_MM);
		double err_dx = 0.0;
		double ratio = NAN;
		double time_grad = 0.0, time_hess = 0.0, time_ref_grad, time_ref_hess;

		if (ratio_ok) {
			// this override options, as the decision is that it is most efficient to update MM all the time
			update_MM = 1;
			keep_MM = 0;
		}

		gsl_vector_set_zero(B);
		gsl_vector_set_zero(MB);
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
		time_grad = 0.0;
		time_ref_grad = GMRFLib_cpu();
		GMRFLib_preopt_predictor_moments(pmean, NULL, preopt, ai_store->problem, x_mean);
		time_grad += GMRFLib_cpu() - time_ref_grad;
		time_ref_grad = GMRFLib_cpu();

#define CODE_BLOCK							\
		for (int ii = 0; ii < d_idx->n; ii++) {			\
			int i = d_idx->idx[ii];				\
			GMRFLib_vb_coofs_tp vb_coof = {.coofs = {NAN, NAN, NAN}}; \
			GMRFLib_ai_vb_prepare_mean(thread_id, &vb_coof, i, d[i], loglFunc, loglFunc_arg, x_mean, pmean[i], sqrt(pvar[i])); \
			if (debug) {					\
				fprintf(fp, "[%1d] i %d (mean,sd) = %.6f %.6f (A,B,C) = %.6f %.6f %.6f\n", omp_get_thread_num(), i, \
				       pmean[i], sqrt(pvar[i]), vb_coof.coofs[0], vb_coof.coofs[1], vb_coof.coofs[2]); \
			}						\
			BB[i] = vb_coof.coofs[1];			\
			CC[i] = DMAX(0.0, vb_coof.coofs[2]);		\
		}

		RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 0, 0);
#undef CODE_BLOCK

		GMRFLib_preopt_update(thread_id, preopt, BB, CC);
		GMRFLib_Qx(thread_id, tmp, x_mean, preopt->latent_graph, prior->Qfunc, prior->Qfunc_arg);
		for (int i = 0; i < graph->n; i++) {
			tmp[i] += preopt->total_b[thread_id][i];
			gsl_vector_set(B, i, tmp[i]);
		}
		time_grad += GMRFLib_cpu() - time_ref_grad;
		time_ref_grad = GMRFLib_cpu();

		if (update_MM) {
			time_ref_hess = GMRFLib_cpu();
			GMRFLib_QM(thread_id, QM, M, graph, tabQ->Qfunc, tabQ->Qfunc_arg);
			time_hess += GMRFLib_cpu() - time_ref_hess;
		}

		if (update_MM) {
			time_ref_hess = GMRFLib_cpu();
			gsl_blas_dgemm(CblasTrans, CblasNoTrans, one, M, QM, zero, MM);
			time_hess += GMRFLib_cpu() - time_ref_hess;
		}

		// no timing; this one is common
		gsl_blas_dgemv(CblasTrans, mone, M, B, zero, MB);

		if (debug) {
			FIXME("M");
			GMRFLib_printf_gsl_matrix(stdout, M, "%.6f ");
			FIXME("B");
			GMRFLib_printf_gsl_vector(stdout, B, "%.6f ");
		}

		// the system can be singular, like with intrinsic model components. its safe to invert the non-singular part only
		if (keep_MM) {
			// in this case, keep the inv of MM through the iterations
			if (update_MM) {
				time_ref_hess = GMRFLib_cpu();
				GMRFLib_gsl_spd_inv(MM, GMRFLib_eps(1.0 / 3.0));
				time_hess += GMRFLib_cpu() - time_ref_hess;
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
			GMRFLib_gsl_safe_spd_solve(MM, MB, delta, GMRFLib_eps(1.0 / 3.0));
		}

		if (iter == 0) {
			// assume we need twice the number of iterations, then this is the decision value
			ratio = time_hess / time_grad;
			ratio_ok = (ratio < 1.0);
		} else {
			ratio = NAN;
		}

		for (int i = 0; i < (int) delta->size; i++) {
			if (ISNAN(gsl_vector_get(delta, i))) {
				gsl_vector_set_zero(delta);
				break;
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
				dx[i] = max_correct * sd[i] * SIGN(dx[i]);
			}
		}

#pragma GCC ivdep
		for (int i = 0; i < graph->n; i++) {
			x_mean[i] += dx[i];
		}

		double max_correction = 0.0;
#pragma GCC ivdep
		for (int i = 0; i < graph->n; i++) {
			max_correction = DMAX(max_correction, ABS(x_mean[i] - x_mean_orig[i]) / sd[i]);
		}
		if (max_correction >= ai_par->vb_emergency) {
#pragma omp critical (Name_1169f76e685daed4d69fb5a745f9e95b4f5f633b)
			{
				fprintf(stderr, "\n\n\t*** max_correction = %.2f >= %.2f, so 'vb.correction' is aborted\n",
					max_correction, ai_par->vb_emergency);
				fprintf(stderr, "\t*** Please (re-)consider your model, priors, confounding, etc.\n");
				fprintf(stderr, "\t*** You can change the emergency value (current value=%.2f) by \n", ai_par->vb_emergency);
				fprintf(stderr, "\t*** \t'control.inla=list(control.vb=list(emergency=...))'\n\n");
				if (fp != stderr) {
					fprintf(fp, "\n\n\t*** max_correction = %.2f >= %.2f, so 'vb.correction' is aborted\n",
						max_correction, ai_par->vb_emergency);
					fprintf(fp, "\t*** Please (re-)consider your model, priors, confounding, etc.\n");
					fprintf(fp, "\t*** You can change the emergency value (current value=%.2f) by \n", ai_par->vb_emergency);
					fprintf(fp, "\t*** \t'control.inla=list(control.vb=list(emergency=...))'\n\n");
				}
			}
			emergency = 1;
			break;
		}

		// this is RMS standardized change between the iterations, otherwise, just run the max iterations.
		// test this here so we can adjust the verbose output below
		int do_break = 0;
		if (err_dx < 0.01 || iter == niter - 1) {
			do_break = 1;
		}

		if (verbose) {
#pragma omp critical (Name_d9343cf5e9cd69d222c869579102b5231d628874)
			{
				fprintf(fp, "\t[%1d]Iter [%1d/%1d] VB correct with strategy [MEAN] in total[%.3fsec] hess/grad[%.3f]\n",
					omp_get_thread_num(), iter, niter, GMRFLib_cpu() - tref, ratio);
				fprintf(fp, "\t\tNumber of nodes corrected for [%1d] max(dx/sd)[%.4f]\n", (int) delta->size, err_dx);
				if (do_break) {
					for (int jj = 0; jj < vb_idx->n; jj++) {
						int j = vb_idx->idx[jj];
						fprintf(fp, "\t\tNode[%1d] delta[%.3f] dx/sd[%.3f] |x-mode|/sd[%.3f]\n", j,
							gsl_vector_get(delta_mu, j), dx[j] / sd[j],
							(x_mean[j] - ai_store->problem->mean_constr[j]) / sd[j]);
					}
					fprintf(fp, "\t\tImplied correction for [%1d] nodes\n", preopt->mnpred + graph->n - vb_idx->n);
				}
			}
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

	GMRFLib_idx_free(d_idx);
	GMRFLib_idx_free(vb_idx);
	Calloc_free();

	GMRFLib_LEAVE_ROUTINE;
#undef SHOW_TIME
	return GMRFLib_SUCCESS;
}

int GMRFLib_ai_vb_correct_variance_preopt(int thread_id,
					  GMRFLib_density_tp *** density,
					  int dens_count,
					  double *c,
					  double *d,
					  GMRFLib_ai_param_tp * ai_par,
					  GMRFLib_ai_store_tp * ai_store,
					  GMRFLib_graph_tp * graph,
					  GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg,
					  GMRFLib_preopt_tp * preopt)
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

	double tref = GMRFLib_cpu();
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
	GMRFLib_idx_tp *vb_idx = NULL, *d_idx = NULL;

	for (int i = 0; i < graph->n; i++) {
		if (ai_par->vb_nodes_variance[i]) {
			GMRFLib_idx_add(&vb_idx, i);
		}
	}
	if (!vb_idx) {
		GMRFLib_LEAVE_ROUTINE;
		return GMRFLib_SUCCESS;
	}

	for (int i = 0; i < preopt->Npred; i++) {
		if (d[i]) {
			GMRFLib_idx_add(&d_idx, i);
		}
	}
	assert(d_idx);

	Calloc_init(7 * graph->n + 2 * preopt->mnpred + 4 * preopt->Npred + 1 * vb_idx->n, 13);
	double *x_mean = Calloc_get(graph->n);
	double *mean_constr = Calloc_get(graph->n);
	double *pmean = Calloc_get(preopt->mnpred);
	double *pvar = Calloc_get(preopt->mnpred);

#pragma GCC ivdep
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
	double diff_sigma_max_limit = 0.01;
	GMRFLib_problem_tp *problem = NULL;

	// main loop
	int iter = 0;
	for (iter = 0; iter < niter + 1; iter++) {	       /* yes, it should be '+1': we 'break' at the top */

		if (enable_tref_a) {
			tref_a[0] -= GMRFLib_cpu();
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
			tref_a[0] += GMRFLib_cpu();
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
						tn, iter, niter, GMRFLib_cpu() - tref);
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
			tref_a[1] -= GMRFLib_cpu();
		}

#define CODE_BLOCK							\
		for (int ii = 0; ii < d_idx->n; ii++) {			\
			int i = d_idx->idx[ii];				\
			GMRFLib_vb_coofs_tp vb_coof = {.coofs = {NAN, NAN, NAN}}; \
			GMRFLib_ai_vb_prepare_variance(thread_id, &vb_coof, i, d[i], loglFunc, loglFunc_arg, x_mean, pmean[i], sqrt(pvar[i])); \
			if (debug && 0) {				\
				fprintf(fp, "\t[%1d] i %d (mean,sd) = %.6f %.6f (A,B,C) = %.6f %.6f %.6f\n", tn, i, \
				       pmean[i], sqrt(pvar[i]), vb_coof.coofs[0], vb_coof.coofs[1], vb_coof.coofs[2]); \
			}						\
			BB[i] = vb_coof.coofs[1];			\
			CC[i] = DMAX(0.0, vb_coof.coofs[2]);		\
		}

		RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 0, 0);
#undef CODE_BLOCK

		if (enable_tref_a) {
			tref_a[1] += GMRFLib_cpu();
		}

		GMRFLib_idxval_tp **A = NULL;
		if (preopt->pA_idxval) {
			A = preopt->pAA_idxval;
		} else {
			A = preopt->A_idxval;
		}

		int UNUSED(integer_one) = 1;

#define COV_ETA_LATENT(value_, k_, cov_latent_) (value_) = GMRFLib_dot_product(A[k_], cov_latent_)

#define COMPUTE_COV_LATENT(cov_latent_, j_, b_)				\
		if (1) {						\
			int _j = j_;					\
			Memset(b_, 0, graph->n * sizeof(double));	\
			b_[_j] = 1.0;					\
			GMRFLib_Qsolve(cov_latent_, b_, problem, _j);	\
		}

		if (enable_tref_a) {
			tref_a[2] -= GMRFLib_cpu();
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
			tref_a[2] += GMRFLib_cpu();
		}

#define CODE_BLOCK							\
		for (int ii = 0; ii < vb_idx->n; ii++) {		\
			if (enable_tref_a) {				\
				tref_a[3] -= GMRFLib_cpu();		\
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
				tref_a[3] += GMRFLib_cpu();		\
			}						\
			if (debug) {					\
				fprintf(fp, "\t[%1d]Iter [%1d/%1d] gradient[%1d] = %f\n", tn, iter, niter, ii, gsl_vector_get(gradient, (size_t) ii)); \
			}						\
									\
			if (enable_tref_a) {				\
				tref_a[4] -= GMRFLib_cpu();		\
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
				tref_a[4] += GMRFLib_cpu();		\
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
			tref_a[5] -= GMRFLib_cpu();
		}

		grad_err = 0.0;
		for (int i = 0; i < vb_idx->n; i++) {
			grad_err += SQR(gsl_vector_get(gradient, (size_t) i));
		}
		grad_err = sqrt(grad_err / (double) vb_idx->n);

		if (iter < hessian_update) {
			GMRFLib_gsl_spd_inv(hessian, GMRFLib_eps(1.0 / 3.0));
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
			tref_a[5] += GMRFLib_cpu();
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
	GMRFLib_idx_free(d_idx);
	GMRFLib_idx_free(vb_idx);
	Calloc_free();

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_transform_density(GMRFLib_density_tp ** UNUSED(tdensity), GMRFLib_density_tp * UNUSED(density),
			      GMRFLib_transform_array_func_tp * UNUSED(func))
{
	fprintf(stderr, "\n\n\n");
	fprintf(stderr, "DISABLE THIS FEATURE FOR NOW, DO NOT KNOW HOW TO DO THIS WELL AT THE MOMENT.\n");
	fprintf(stderr, "SINCE THE SCALE OF THE Xs CAN BE SO DIFFERENT, WE WILL NEED A NEW APPROACH OF HOW\n");
	fprintf(stderr, "TO REPRESENT AND COMPUTE MIXTURES OF THESE DENSITIES.  SEE INFO ABOUT THIS AT: ISSUES\n\n.");
	GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);		       /* never enter this function */
	abort();

	return GMRFLib_SUCCESS;
}

int GMRFLib_ai_store_config(int thread_id, GMRFLib_ai_misc_output_tp * mo, int ntheta, double *theta, double log_posterior,
			    double log_posterior_orig, double *improved_mean, double *skewness, GMRFLib_problem_tp * gmrf_approx,
			    GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg, double *c, int dens_count)
{
	if (!mo || !(mo->configs)) {
		return GMRFLib_SUCCESS;
	}

	const int debug = 0;
	int id = omp_get_thread_num();

	if (!(mo->configs[id])) {
		mo->configs[id] = Calloc(1, GMRFLib_store_configs_tp);
		GMRFLib_graph_tp *g;
		GMRFLib_graph_duplicate(&g, gmrf_approx->sub_graph);
		if (debug) {
			printf("remapped graph\n");
			GMRFLib_printf_graph(stdout, g);
		}

		int nelm = gmrf_approx->sub_graph->n + gmrf_approx->sub_graph->nnz;
		mo->configs[id]->n = gmrf_approx->sub_graph->n;
		mo->configs[id]->nz = (nelm - mo->configs[id]->n) / 2 + mo->configs[id]->n;
		mo->configs[id]->ntheta = ntheta;

		if (debug) {
			printf("n nz ntheta %d %d %d\n", mo->configs[id]->n, mo->configs[id]->nz, mo->configs[id]->ntheta);
		}

		GMRFLib_constr_tp *cc;
		GMRFLib_duplicate_constr(&cc, gmrf_approx->sub_constr, gmrf_approx->sub_graph);	/* might or might not be mapped ???? */
		mo->configs[id]->constr = cc;

		if (debug) {
			printf("constraints\n");
			GMRFLib_printf_constr(stdout, cc, gmrf_approx->sub_graph);
		}

		int *i, *j, ii, jj, k, kk;
		i = Calloc(mo->configs[id]->nz, int);
		j = Calloc(mo->configs[id]->nz, int);

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
		mo->configs[id]->i = i;
		mo->configs[id]->j = j;
		mo->configs[id]->nconfig = 0;
		mo->configs[id]->config = NULL;
	}

	mo->configs[id]->config = Realloc(mo->configs[id]->config, mo->configs[id]->nconfig + 1, GMRFLib_store_config_tp *);
	mo->configs[id]->config[mo->configs[id]->nconfig] = Calloc(1, GMRFLib_store_config_tp);

	int ii, jj, k, kk, found = 0;
	double *Qinv = NULL, *Q = NULL, *Qprior = NULL, *mean = NULL, *imean = NULL, *skew = NULL;
	GMRFLib_graph_tp *g = gmrf_approx->sub_graph;

	Q = Calloc(mo->configs[id]->nz, double);
	if (gmrf_approx->tab->Qfunc == GMRFLib_tabulate_Qfunction_std) {
		GMRFLib_tabulate_Qfunc_arg_tp *aa;
		aa = (GMRFLib_tabulate_Qfunc_arg_tp *) gmrf_approx->tab->Qfunc_arg;
		if (aa->Q) {
			assert(mo->configs[id]->nz == aa->Q->s->na);
			Memcpy(Q, aa->Q->a, aa->Q->s->na * sizeof(double));
			found = 1;
		}
	}
	if (!found) {
		for (ii = k = 0; ii < g->n; ii++) {
			Q[k++] = gmrf_approx->tab->Qfunc(thread_id, ii, ii, NULL, gmrf_approx->tab->Qfunc_arg);
			for (kk = 0; kk < g->lnnbs[ii]; kk++) {
				jj = g->lnbs[ii][kk];
				Q[k++] = gmrf_approx->tab->Qfunc(thread_id, ii, jj, NULL, gmrf_approx->tab->Qfunc_arg);
			}
		}
	}

	if (g->n >= 0) {
		Qprior = Calloc(g->n, double);
		for (ii = 0; ii < g->n; ii++) {
			Qprior[ii] = Qfunc(thread_id, ii, ii, NULL, Qfunc_arg) + c[ii];
		}
	}

	mean = Calloc(g->n, double);
	imean = Calloc(g->n, double);
	skew = Calloc(g->n, double);
	Memcpy(mean, gmrf_approx->mean_constr, g->n * sizeof(double));
	Memcpy(imean, improved_mean, g->n * sizeof(double));
	Memcpy(skew, skewness, g->n * sizeof(double));

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
	mo->configs[id]->config[mo->configs[id]->nconfig]->Qprior = Qprior;
	mo->configs[id]->config[mo->configs[id]->nconfig]->mean = mean;
	mo->configs[id]->config[mo->configs[id]->nconfig]->improved_mean = imean;
	mo->configs[id]->config[mo->configs[id]->nconfig]->skewness = skew;
	mo->configs[id]->config[mo->configs[id]->nconfig]->log_posterior = log_posterior;	/* may include integration weights */
	mo->configs[id]->config[mo->configs[id]->nconfig]->log_posterior_orig = log_posterior_orig;	/* do NOT include integration weights */
	mo->configs[id]->config[mo->configs[id]->nconfig]->dens_count = dens_count;

	if (mo->configs[id]->ntheta) {
		mo->configs[id]->config[mo->configs[id]->nconfig]->theta = Calloc(mo->configs[id]->ntheta, double);
		Memcpy(mo->configs[id]->config[mo->configs[id]->nconfig]->theta, theta, mo->configs[id]->ntheta * sizeof(double));
	} else {
		mo->configs[id]->config[mo->configs[id]->nconfig]->theta = NULL;
	}
	mo->configs[id]->nconfig++;

	return GMRFLib_SUCCESS;
}

int GMRFLib_ai_store_config_preopt(int thread_id, GMRFLib_ai_misc_output_tp * mo, int ntheta, double *theta, double log_posterior,
				   double log_posterior_orig, GMRFLib_problem_tp * problem, double *mean_corrected,
				   GMRFLib_preopt_tp * preopt, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg, double *cpodens_moments,
				   double *gcpodens_moments, char **arg_str, double *ll_info)
{
	if (!mo || !(mo->configs_preopt)) {
		return GMRFLib_SUCCESS;
	}
	int id = omp_get_thread_num();

	if (!(mo->configs_preopt[id])) {

		GMRFLib_graph_tp *g;
		mo->configs_preopt[id] = Calloc(1, GMRFLib_store_configs_preopt_tp);

		mo->configs_preopt[id]->mnpred = preopt->mnpred;
		mo->configs_preopt[id]->Npred = preopt->Npred;
		mo->configs_preopt[id]->n = preopt->n;
		mo->configs_preopt[id]->ntheta = ntheta;

		int nelm = preopt->preopt_graph->n + preopt->preopt_graph->nnz;
		mo->configs_preopt[id]->nz = (nelm - mo->configs_preopt[id]->n) / 2 + mo->configs_preopt[id]->n;
		nelm = preopt->latent_graph->n + preopt->latent_graph->nnz;
		mo->configs_preopt[id]->prior_nz = (nelm - mo->configs_preopt[id]->n) / 2 + mo->configs_preopt[id]->n;
		mo->configs_preopt[id]->A = preopt->A;
		mo->configs_preopt[id]->pA = preopt->pA;
		GMRFLib_duplicate_constr(&(mo->configs_preopt[id]->constr), preopt->latent_constr, preopt->preopt_graph);

		int *i, *j, ii, jj, k, kk;
		i = Calloc(mo->configs_preopt[id]->nz, int);
		j = Calloc(mo->configs_preopt[id]->nz, int);
		g = preopt->preopt_graph;
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
		mo->configs_preopt[id]->i = i;
		mo->configs_preopt[id]->j = j;

		g = preopt->latent_graph;
		i = Calloc(mo->configs_preopt[id]->prior_nz, int);
		j = Calloc(mo->configs_preopt[id]->prior_nz, int);
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
		mo->configs_preopt[id]->iprior = i;
		mo->configs_preopt[id]->jprior = j;

		mo->configs_preopt[id]->nconfig = 0;
		mo->configs_preopt[id]->config = NULL;
	}

	int nconfig = mo->configs_preopt[id]->nconfig;
	mo->configs_preopt[id]->config = Realloc(mo->configs_preopt[id]->config, nconfig + 1, GMRFLib_store_config_preopt_tp *);
	mo->configs_preopt[id]->config[nconfig] = Calloc(1, GMRFLib_store_config_preopt_tp);

	int ii, jj, k, kk;
	double *Qinv = NULL, *Q = NULL, *Qprior = NULL, *mean = NULL, *imean = NULL;
	GMRFLib_graph_tp *g;

	Q = Calloc(mo->configs_preopt[id]->nz, double);
	g = preopt->preopt_graph;
	for (ii = k = 0; ii < g->n; ii++) {
		Q[k++] = Qfunc(thread_id, ii, ii, NULL, Qfunc_arg);
		for (kk = 0; kk < g->lnnbs[ii]; kk++) {
			jj = g->lnbs[ii][kk];
			Q[k++] = Qfunc(thread_id, ii, jj, NULL, Qfunc_arg);
		}
	}

	Qprior = Calloc(mo->configs_preopt[id]->prior_nz, double);
	assert(Qfunc == GMRFLib_preopt_Qfunc);
	g = preopt->latent_graph;
	for (ii = k = 0; ii < g->n; ii++) {
		Qprior[k++] = GMRFLib_preopt_Qfunc_prior(thread_id, ii, ii, NULL, Qfunc_arg);
		for (kk = 0; kk < g->lnnbs[ii]; kk++) {
			jj = g->lnbs[ii][kk];
			Qprior[k++] = GMRFLib_preopt_Qfunc_prior(thread_id, ii, jj, NULL, Qfunc_arg);
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

int GMRFLib_bnew(int thread_id, double **bnew, double *constant, int n, double *b, GMRFLib_bfunc_tp ** bfunc)
{
	/*
	 * bnew is a new alloced ptr for the new b. constant is the missing constant to be added due to b=b(theta).
	 */

	int i;
	double *bb = Calloc(n, double);
	double con = 0.0, con_add = 0.0;

	if (b) {
		Memcpy((void *) bb, (void *) b, n * sizeof(double));
	}

	if (bfunc) {
		for (i = 0; i < n; i++) {
			if (bfunc[i]) {
				bb[i] += GMRFLib_bfunc_eval(thread_id, &con_add, bfunc[i]);
				con += con_add;
			}
		}
	}

	*bnew = bb;
	*constant = -con / 2.0;

	return GMRFLib_SUCCESS;
}

int GMRFLib_ai_compute_lincomb(GMRFLib_density_tp *** lindens, double **cross, int nlin, GMRFLib_lc_tp ** Alin,
			       GMRFLib_ai_store_tp * ai_store, double *improved_mean, int lookup_tables)
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
	GMRFLib_density_tp **d;

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
				double *p, *pp, w, ww;

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

int GMRFLib_ai_correct_cpodens(double *logdens, double *x, int *n, GMRFLib_ai_param_tp * ai_par)
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

double GMRFLib_ai_cpopit_integrate(int thread_id, double *cpo, double *pit, int idx, GMRFLib_density_tp * cpo_density, double d,
				   GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, double *x_vec)
{
	/*
	 * cpo_density is the marginal for x_idx without y_idx, density: is the marginal for x_idx with y_idx.
	 */
	int retval, compute_cpo = 1, i, k, np = GMRFLib_INT_NUM_POINTS;
	double low, dx, dxi, *xp = NULL, *xpi = NULL, *dens = NULL, *prob = NULL,
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
	for (i = 1; i < np; i++) {
		xp[i] = xp[0] + i * dx;
		xpi[i] = xpi[0] + i * dxi;
	}
	GMRFLib_evaluate_ndensity(dens, xpi, np, cpo_density);

	if (compute_cpo) {
		loglFunc(thread_id, prob, xp, -np, idx, x_vec, NULL, loglFunc_arg, NULL);	/* no correction for 'd' here; should we? */
	} else {
		Memset(prob, 0, np * sizeof(double));
	}
	loglFunc(thread_id, loglik, xp, np, idx, x_vec, NULL, loglFunc_arg, NULL);
	for (i = 0; i < np; i++) {
		loglik[i] *= d;
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

	Calloc_free();
	return fail;
}

double GMRFLib_ai_po_integrate(int thread_id, double *po, double *po2, double *po3, int idx, GMRFLib_density_tp * po_density,
			       double d, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, double *x_vec)
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

		Calloc_init(3 * np, 3);
		double *x = Calloc_get(np);
		double *ll = Calloc_get(np);
		double *mask = Calloc_get(np);

		for (int i = 0; i < np; i++) {
			x[i] = mean + stdev * xp[i];
		}
		loglFunc(thread_id, ll, x, np, idx, x_vec, NULL, loglFunc_arg, NULL);
		double dmax = GMRFLib_max_value(ll, np, NULL);
		double limit = -0.5 * SQR(xp[0]);	       // prevent extreme values
		for (int i = 0; i < np; i++) {
			if (ll[i] - dmax < limit) {
				mask[i] = 0.0;
				ll[i] = 0.0;
			} else {
				mask[i] = 1.0;
			}
		}

		integral2 = 0.0;
		integral3 = 0.0;
		integral4 = 0.0;
		for (int i = 0; i < np; i++) {
			integral2 += mask[i] * exp(ll[i]) * wp[i];
			integral3 += ll[i] * wp[i];
			integral4 += SQR(ll[i]) * wp[i];
		}
		Calloc_free();
	} else {
		int i, k;
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
		for (i = 1; i < np; i++) {
			xp[i] = xp[0] + i * dx;
			xpi[i] = xpi[0] + i * dxi;
		}
		GMRFLib_evaluate_nlogdensity(ldens, xpi, np, po_density);
		loglFunc(thread_id, loglik, xp, np, idx, x_vec, NULL, loglFunc_arg, NULL);

		double *dens = Calloc_get(npm);
		double *llik = Calloc_get(npm);

		if (GMRFLib_INT_NUM_INTERPOL == 3) {
#pragma GCC ivdep
			for (i = 0; i < np - 1; i++) {
				llik[3 * i + 0] = loglik[i];
				llik[3 * i + 1] = (2.0 * loglik[i] + loglik[i + 1]) / 3.0;
				llik[3 * i + 2] = (loglik[i] + 2.0 * loglik[i + 1]) / 3.0;
			}
#pragma GCC ivdep
			for (i = 0; i < np - 1; i++) {
				dens[3 * i + 0] = exp(ldens[i]);
				dens[3 * i + 1] = exp((2.0 * ldens[i] + ldens[i + 1]) / 3.0);
				dens[3 * i + 2] = exp((ldens[i] + 2.0 * ldens[i + 1]) / 3.0);
			}
			llik[3 * (np - 2) + 3] = loglik[np - 1];
			dens[3 * (np - 2) + 3] = exp(ldens[np - 1]);
			assert(3 * (np - 2) + 3 == npm - 1);
		} else if (GMRFLib_INT_NUM_INTERPOL == 2) {
#pragma GCC ivdep
			for (i = 0; i < np - 1; i++) {
				llik[2 * i + 0] = loglik[i];
				llik[2 * i + 1] = (loglik[i] + loglik[i + 1]) / 2.0;
			}
#pragma GCC ivdep
			for (i = 0; i < np - 1; i++) {
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
		for (i = 1, k = 0; i < npm - 1; i++, k = (k + 1) % 2) {
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

double *GMRFLib_ai_dic_integrate(int thread_id, int idx, GMRFLib_density_tp * density, double d, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg,
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

		for (int i = 0; i < np; i++) {
			x[i] = mean + stdev * xp[i];
		}
		loglFunc(thread_id, ll, x, np, idx, x_vec, NULL, loglFunc_arg, NULL);
		for (int i = 0; i < np; i++) {
			ll_sat[i] = ll[i] - sat_ll;
		}

		double dmax = GMRFLib_max_value(ll, np, NULL);
		double limit = -0.5 * SQR(xp[0]);	       // prevent extreme values
		for (int i = 0; i < np; i++) {
			if (ll[i] - dmax < limit) {
				ll[i] = ll_sat[i] = 0.0;
			}
		}
		integral = integral_sat = 0.0;
		for (int i = 0; i < np; i++) {
			integral += ll[i] * wp[i];
			integral_sat += ll_sat[i] * wp[i];
		}
		integral = -2.0 * d * integral;
		integral_sat = -2.0 * d * integral_sat;
		Calloc_free();
	} else {
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
#pragma GCC ivdep
			for (int i = 0; i < np - 1; i++) {
				llik[3 * i + 0] = loglik[i];
				llik[3 * i + 1] = (2.0 * loglik[i] + loglik[i + 1]) / 3.0;
				llik[3 * i + 2] = (loglik[i] + 2.0 * loglik[i + 1]) / 3.0;

				llik_sat[3 * i + 0] = llik[3 * i + 0] - sat_ll;
				llik_sat[3 * i + 1] = llik[3 * i + 1] - sat_ll;
				llik_sat[3 * i + 2] = llik[3 * i + 2] - sat_ll;

			}
#pragma GCC ivdep
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
#pragma GCC ivdep
			for (int i = 0; i < np - 1; i++) {
				llik[2 * i + 0] = loglik[i];
				llik[2 * i + 1] = (loglik[i] + loglik[i + 1]) / 2.0;

				llik_sat[2 * i + 0] = llik[2 * i + 0] - sat_ll;
				llik_sat[2 * i + 1] = llik[2 * i + 1] - sat_ll;
			}
#pragma GCC ivdep
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
		int i, n;

		GMRFLib_Qinv(ai_store->problem);
		Free(ai_store->stdev);
		n = ai_store->problem->n;
		ai_store->stdev = Calloc(n, double);
		for (i = 0; i < n; i++) {
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

int GMRFLib_ai_marginal_one_hyperparamter(GMRFLib_density_tp ** density, int idx, int nhyper, int hyper_count, double *hyper_z,
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
		double ldens_min = log(GMRFLib_eps(1.0));
		for (i = 0; i < npoints; i++) {
			ldens_values[i] = log(dens[i]);
			ldens_values[i] = IMAX(ldens_values[i], ldens_min);
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

double GMRFLib_interpolator_nearest(int ndim, int nobs, double *x, double *xobs, double *yobs, void *UNUSED(arg))
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

double GMRFLib_interpolator_linear(int ndim, int nobs, double *x, double *xobs, double *yobs, void *UNUSED(arg))
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

GMRFLib_ai_store_tp *GMRFLib_duplicate_ai_store(GMRFLib_ai_store_tp * ai_store, int skeleton, int copy_ptr, int copy_pardiso_ptr)
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
	int nd = ai_store->nd;
	int Npred = ai_store->Npred;

	new_ai_store->store = GMRFLib_duplicate_store(ai_store->store, skeleton, copy_ptr, copy_pardiso_ptr);
	new_ai_store->problem = GMRFLib_duplicate_problem(ai_store->problem, skeleton, copy_ptr, copy_pardiso_ptr);
	COPY(nidx);
	COPY(nd);
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
	DUPLICATE(d_idx, nd, int, 0);

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
	ASSIGN(derivative4);
	ASSIGN(mode);
	ASSIGN(nc_orig);
	ASSIGN(nd);
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
	size_t i, j, k;
	const int debug = 0;
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
	p->configurations = Calloc((size_t) (p->nconfig * p->nhyper), int);
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
			p->configurations[k + j] = (int) izz[j];
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
	if (GMRFLib_MAX_THREADS() > 1) {
		qsort(p->configurations, p->nconfig, p->nhyper * sizeof(int), GMRFLib_pool_cmp);
	} else {
		/*
		 * alternative sorting: _pool_cmp1: seems like _pool_cmp runs faster (better wrt 'reject') 
		 */
		qsort(p->configurations, p->nconfig, p->nhyper * sizeof(int), GMRFLib_pool_cmp);
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
	const int *ia, *ib;
	int i, larger, dist_a, dist_b;

	GMRFLib_ASSERT(pool_nhyper > 0, GMRFLib_ESNH);

	ia = (const int *) a;
	ib = (const int *) b;
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
	const int *ia, *ib;
	int i, larger, eq;

	GMRFLib_ASSERT(pool_nhyper > 0, GMRFLib_ESNH);

	ia = (const int *) a;
	ib = (const int *) b;
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

int GMRFLib_ai_pool_get(GMRFLib_ai_pool_tp * pool, int *iz, size_t *idx)
{
	return GMRFLib_ai_pool_intern(pool, iz, idx, 0.0, GMRFLib_AI_POOL_GET);
}

int GMRFLib_ai_pool_set(GMRFLib_ai_pool_tp * pool, size_t idx, double logdens)
{
	return GMRFLib_ai_pool_intern(pool, NULL, &idx, logdens, GMRFLib_AI_POOL_SET);
}

int GMRFLib_ai_pool_intern(GMRFLib_ai_pool_tp * pool, int *iz, size_t *idx, double logdens, int action)
{
	const int debug = 0;
	int retval = 0;

	GMRFLib_ASSERT(idx, GMRFLib_EPARAMETER);
#pragma omp critical (Name_bbcc3bdca5c5de72f4e6fce04148a03434b2c017)
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

double GMRFLib_bfunc_eval(int thread_id, double *constant, GMRFLib_bfunc_tp * bfunc)
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
