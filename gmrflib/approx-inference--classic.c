
/* approx-inference--classic.c
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


static int pool_nhyper = -1;
int GMRFLib_ai_log_posterior(int thread_id, double *logdens,
			     double *x, double *b, double *c, double *mean, double *d,
			     GMRFLib_logl_tp *loglFunc, void *loglFunc_arg,
			     GMRFLib_graph_tp *graph, GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_arg, GMRFLib_constr_tp *UNUSED(constr))
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
			// THIS CODE IS VERY STUPID!!! FIX

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
					GMRFLib_logl_tp *loglFunc, void *loglFunc_arg,
					GMRFLib_graph_tp *graph, GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_arg,
					GMRFLib_constr_tp *UNUSED(constr), GMRFLib_graph_tp *subgraph, GMRFLib_ai_store_tp *UNUSED(ai_store),
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

int GMRFLib_ai_marginal_hidden(int thread_id, GMRFLib_density_tp **density, GMRFLib_density_tp **cpo_density,
			       int lookup_tables,
			       int idx, double *x, double *b, double *c, double *mean, double *d,
			       GMRFLib_logl_tp *loglFunc, void *loglFunc_arg,
			       GMRFLib_graph_tp *graph, GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_arg,
			       GMRFLib_constr_tp *constr, GMRFLib_ai_param_tp *ai_par, GMRFLib_ai_store_tp *ai_store, GMRFLib_preopt_tp *preopt)
{
	/*
	 * compute the approximation to the marginal for the hidden field at index 'idx' and return the density in *density. if
	 * (cpo_density), then compute also the density when y_idx is removed.
	 */

	char *fix = NULL, *fixx = NULL;
	int i, j, k, n = -1, free_ai_par = 0, n_points, ii, free_ai_store = 0, i_idx, j_idx, one = 1, *node_map = NULL;
	double *x_points = NULL, x_sd, x_mean, *cond_mode = NULL, *fixed_mode = NULL, *log_density = NULL,
	    log_dens_cond = 0.0, deriv_log_dens_cond = 0.0, a, *derivative = NULL, *mean_and_variance = NULL, deldif =
	    GSL_ROOT6_DBL_EPSILON, inv_stdev, *cov = NULL, corr, corr_term, *covariances = NULL, alpha;

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
			double cor_eps = (GSL_SQRT_DBL_EPSILON * GSL_ROOT4_DBL_EPSILON), cor_max, range;	\
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
			GMRFLib_idx_tp *dd = NULL;
			GMRFLib_idx_create(&dd);
#pragma omp critical (Name_55d88ef833913b76c8f458812b76256a0492204c)
			if (!ai_store->d_idx) {
				for (i = 0; i < n; i++) {
					if (d[i]) {
						GMRFLib_idx_add(&dd, i);
					}
				}
			}
			ai_store->d_idx = dd;
		}
	} else {
#pragma omp critical (Name_2f4281b0d3628f0c97f72ada8a31d4cad908cd4e)
		GMRFLib_idx_create(&(ai_store->d_idx));
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
									     b, c, mean, d, NULL, loglFunc, loglFunc_arg,
									     graph, Qfunc, Qfunc_arg, constr, optpar, blockpar,
									     ai_store->store, ai_store->aa, ai_store->bb,
									     ai_store->cc, ai_par->cmin, ai_par->b_strategy, 0, preopt, NULL));
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
	QSORT_FUN(x_points, (size_t) n_points, sizeof(double), GMRFLib_dcmp_abs);	/* sort them using ABS() */
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
				for (ii = 0; ii < ai_store->d_idx->n; ii++) {
					double dd;
					i = ai_store->d_idx->idx[ii];
					ai_store->correction_idx[ai_store->nidx++] = i;
					double aa, bb, cc;
					GMRFLib_2order_approx(thread_id, &aa, &bb, &cc, &dd, d[i], fixed_mode[i] + deldif, i,
							      fixed_mode, loglFunc, loglFunc_arg, &(ai_par->step_len), &(ai_par->stencil), NULL);
					ai_store->derivative3[i] = dd;
					ai_store->correction_term[i] = -SQR(ai_store->stdev[i]) * ai_store->derivative3[i];
#if defined(INLA_RESEARCH1)
					if (1) {
						// have to redo this later if this gets serious
						double d3[2];
						double s = 1.0 / (2.0 * deldif);
						GMRFLib_2order_approx(thread_id, &aa, &bb, &cc, &(d3[0]), d[i], fixed_mode[i] - deldif, i,
								      fixed_mode, loglFunc, loglFunc_arg,
								      &(ai_par->step_len), &(ai_par->stencil), NULL);
						GMRFLib_2order_approx(thread_id, &aa, &bb, &cc, &(d3[1]), d[i], fixed_mode[i] + deldif, i,
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
			third_order_derivative += ai_store->derivative3[iii] * POW3(derivative[iii]);
		}
		third_order_derivative *= POW3(x_sd);

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
				ld[k] = -0.5 * SQR(xx) - deriv_log_dens_cond * xx + (third_order_derivative / 6.0) * POW3(xx + deriv_log_dens_cond);
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

int GMRFLib_ai_update_conditional_mean(int thread_id, GMRFLib_problem_tp *pproblem, double *x, double *mean,
				       GMRFLib_graph_tp *graph, GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_args,
				       GMRFLib_constr_tp *constr, double *bbb, double *ccc, double **covariances, int idx)
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

int GMRFLib_ai_update_conditional_mean2(double *cond_mean, GMRFLib_problem_tp *problem, int idx, double evalue, double **covariances)
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

int GMRFLib_ai_skip_configurations(map_strd *hash_table, int k, int *iz, int *izz, int *len, int *k_max, int len_length, int nhyper)
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

int GMRFLib_ai_INLA(GMRFLib_density_tp ***density,
		    GMRFLib_density_tp ***density_transform, GMRFLib_transform_array_func_tp **tfunc,
		    GMRFLib_density_tp ***density_hyper,
		    GMRFLib_ai_cpo_tp **cpo, GMRFLib_ai_po_tp **po, GMRFLib_ai_dic_tp *dic,
		    GMRFLib_ai_marginal_likelihood_tp *marginal_likelihood,
		    char *compute, double ***hyperparam, int nhyper,
		    GMRFLib_ai_log_extra_tp *log_extra, void *log_extra_arg,
		    double *x, double *b, double *c, double *mean,
		    GMRFLib_bfunc_tp **bfunc, double *d,
		    GMRFLib_logl_tp *loglFunc, void *loglFunc_arg,
		    GMRFLib_graph_tp *graph, GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_arg,
		    GMRFLib_constr_tp *constr, GMRFLib_ai_param_tp *ai_par, GMRFLib_ai_store_tp *ai_store,
		    int nlin, GMRFLib_lc_tp **Alin, GMRFLib_density_tp ***dlin, GMRFLib_ai_misc_output_tp *misc_output,
		    GMRFLib_preopt_tp *preopt, GMRFLib_preopt_res_tp *rpreopt)
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

	GMRFLib_idx_tp *d_idx = NULL;
	for (i = 0; i < graph->n; i++) {
		if (d[i]) {
			GMRFLib_idx_add(&d_idx, i);
		}
	}

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
		timer[0] = GMRFLib_timer();
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
	assert(x_mode);
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
		timer[0] = GMRFLib_timer() - timer[0];
		timer[1] = GMRFLib_timer();
	}

	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_OPTIMIZE, (void *) &nhyper, NULL);
	if (nhyper > 0) {
		/*
		 * the first step is to locate the mode of \pi(\theta | y). here we use the opt-optimiser routine.  NOTE that this
		 * '_setup' ensure that ai_store is changed for each call to _opt_f. this is a bit dirty programming, but there is no
		 * good way to get around it for the moment.
		 */
		GMRFLib_opt_setup(hyperparam, nhyper, log_extra, log_extra_arg, compute, x, b, c, mean, bfunc, d, NULL, loglFunc,
				  loglFunc_arg, graph, Qfunc, Qfunc_arg, constr, ai_par, ai_store, preopt, d_idx);
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
			timer[1] = GMRFLib_timer() - timer[1];
			timer[2] = GMRFLib_timer();
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

			if (0 && ai_par->vb_enable) {
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

				tref = GMRFLib_timer();

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
				if (0 && ai_par->vb_enable) {
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
				tu = GMRFLib_timer() - tref;
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
					tref = GMRFLib_timer();
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
									" ] log(rel.dens)= %6.3f, reject, %.2fs\n", val, GMRFLib_timer() - tref);
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
						if (0 && ai_par->vb_enable) {
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
						tu = GMRFLib_timer() - tref;

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
			timer[2] = GMRFLib_timer();
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

		GMRFLib_ai_marginal_hyperparam(thread_id, &tmp_logdens, x, bnew, c, mean, d, NULL,
					       loglFunc, loglFunc_arg, graph, Qfunc, Qfunc_arg, constr, ai_par, ai_store, preopt, d_idx);
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

		if (0 && ai_par->vb_enable) {
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
		timer[2] = GMRFLib_timer() - timer[2];
		timer[3] = GMRFLib_timer();
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

	GMRFLib_idx_free(d_idx);
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
		timer[3] = GMRFLib_timer() - timer[3];
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
int GMRFLib_ai_vb_prepare(int thread_id,
			  GMRFLib_vb_coofs_tp *coofs, int idx, GMRFLib_density_tp *density, double d, GMRFLib_logl_tp *loglFunc,
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

		GMRFLib_daxpb(np, s, xp, m, x_user);
		GMRFLib_density_user2std_n(x_std, x_user, density, np);
		loglFunc(thread_id, loglik, x_user, np, idx, x_vec, NULL, loglFunc_arg, NULL);

		double A = 0.0, B = 0.0, C = 0.0, s_inv = 1.0 / s, s2_inv = 1.0 / SQR(s);
#pragma omp simd reduction(+: A, B, C)
		for (int i = 0; i < np; i++) {
			double tmp = wp[i] * loglik[i];
			A += tmp;
			B += tmp * xp[i];
			C += tmp * (SQR(xp[i]) - 1.0);
		}
		coofs->coofs[0] = -d * A;
		coofs->coofs[1] = -d * B * s_inv;
		coofs->coofs[2] = -d * C * s2_inv;

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

int GMRFLib_ai_store_config(int thread_id, GMRFLib_ai_misc_output_tp *mo, int ntheta, double *theta, double log_posterior,
			    double log_posterior_orig, double *improved_mean, double *skewness, GMRFLib_problem_tp *gmrf_approx,
			    GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_arg, double *c, int dens_count)
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

int GMRFLib_bnew(int thread_id, double **bnew, double *constant, int n, double *b, GMRFLib_bfunc_tp **bfunc)
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

int GMRFLib_ai_pool_init(GMRFLib_ai_pool_tp **pool, GMRFLib_ai_param_tp *ai_par, int nhyper)
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
	assert(iz);
	assert(izz);
	Memset(iz, 0, p->nhyper * sizeof(int));
	Memset(izz, 0, p->nhyper * sizeof(int));

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

int GMRFLib_ai_pool_free(GMRFLib_ai_pool_tp *pool)
{
	if (pool) {
		Free(pool->configurations);
		Free(pool->idx_mapping);
		Free(pool->out);
		Free(pool);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_ai_pool_get(GMRFLib_ai_pool_tp *pool, int *iz, size_t *idx)
{
	return GMRFLib_ai_pool_intern(pool, iz, idx, 0.0, GMRFLib_AI_POOL_GET);
}

int GMRFLib_ai_pool_set(GMRFLib_ai_pool_tp *pool, size_t idx, double logdens)
{
	return GMRFLib_ai_pool_intern(pool, NULL, &idx, logdens, GMRFLib_AI_POOL_SET);
}

int GMRFLib_ai_pool_intern(GMRFLib_ai_pool_tp *pool, int *iz, size_t *idx, double logdens, int action)
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

int GMRFLib_transform_density(GMRFLib_density_tp **UNUSED(tdensity), GMRFLib_density_tp *UNUSED(density),
			      GMRFLib_transform_array_func_tp *UNUSED(func))
{
	fprintf(stderr, "\n\n\n");
	fprintf(stderr, "DISABLE THIS FEATURE FOR NOW, DO NOT KNOW HOW TO DO THIS WELL AT THE MOMENT.\n");
	fprintf(stderr, "SINCE THE SCALE OF THE Xs CAN BE SO DIFFERENT, WE WILL NEED A NEW APPROACH OF HOW\n");
	fprintf(stderr, "TO REPRESENT AND COMPUTE MIXTURES OF THESE DENSITIES.  SEE INFO ABOUT THIS AT: ISSUES\n\n.");
	GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);		       /* never enter this function */
	abort();

	return GMRFLib_SUCCESS;
}

int GMRFLib_ai_vb_correct_mean_std(int thread_id, GMRFLib_density_tp ***density,	// need two types
				   int dens_count,
				   GMRFLib_density_tp **dens_local,
				   double *c,
				   double *d,
				   GMRFLib_ai_param_tp *ai_par,
				   GMRFLib_ai_store_tp *ai_store,
				   GMRFLib_graph_tp *graph, GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_arg, GMRFLib_logl_tp *loglFunc, void *loglFunc_arg)
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
		double tref = GMRFLib_timer();
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
		GMRFLib_gsl_ensure_spd(MM, GSL_SQRT_DBL_EPSILON, NULL);
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
			printf("\tVB correct with strategy [mean] in [%.3f]seconds\n", GMRFLib_timer() - tref);
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
