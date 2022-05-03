
/* blockupdate.c
 * 
 * Copyright (C) 2001-2021 Havard Rue
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


#ifndef GITCOMMIT
#define GITCOMMIT
#endif
static const char GitID[] = "file: " __FILE__ "  " GITCOMMIT;

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif

#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

int GMRFLib_default_blockupdate_param(GMRFLib_blockupdate_param_tp ** blockupdate_par)
{
	GMRFLib_ASSERT(blockupdate_par, GMRFLib_EINVARG);

	*blockupdate_par = Calloc(1, GMRFLib_blockupdate_param_tp);
	(*blockupdate_par)->modeoption = GMRFLib_MODEOPTION_MODE;
	(*blockupdate_par)->fp = NULL;
	(*blockupdate_par)->step_len = GMRFLib_eps(0.25);
	(*blockupdate_par)->stencil = 5;

	return GMRFLib_SUCCESS;
}

int GMRFLib_blockupdate(double *laccept,
			double *x_new, double *x_old,
			double *b_new, double *b_old,
			double *c_new, double *c_old,
			double *mean_new, double *mean_old,
			double *d_new, double *d_old,
			GMRFLib_logl_tp * loglFunc_new, void *loglFunc_arg_new,
			GMRFLib_logl_tp * loglFunc_old, void *loglFunc_arg_old,
			GMRFLib_graph_tp * graph,
			GMRFLib_Qfunc_tp * Qfunc_new, void *Qfunc_arg_new,
			GMRFLib_Qfunc_tp * Qfunc_old, void *Qfunc_arg_old,
			GMRFLib_Qfunc_tp * Qfunc_old2new, void *Qfunc_arg_old2new,
			GMRFLib_Qfunc_tp * Qfunc_new2old, void *Qfunc_arg_new2old,
			GMRFLib_constr_tp * constr_new, GMRFLib_constr_tp * constr_old,
			GMRFLib_optimize_param_tp * optpar, GMRFLib_blockupdate_param_tp * blockupdate_par)
{
	GMRFLib_ENTER_ROUTINE;
	GMRFLib_EWRAP1(GMRFLib_blockupdate_store
		       (laccept, x_new, x_old, b_new, b_old, c_new, c_old, mean_new, mean_old, d_new, d_old, loglFunc_new,
			loglFunc_arg_new, loglFunc_old, loglFunc_arg_old, graph, Qfunc_new, Qfunc_arg_new,
			Qfunc_old, Qfunc_arg_old, Qfunc_old2new, Qfunc_arg_old2new, Qfunc_new2old, Qfunc_arg_new2old,
			constr_new, constr_old, optpar, blockupdate_par, NULL));
	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_blockupdate_store(double *laccept,
			      double *x_new, double *x_old,
			      double *b_new, double *b_old,
			      double *c_new, double *c_old,
			      double *mean_new, double *mean_old,
			      double *d_new, double *d_old,
			      GMRFLib_logl_tp * loglFunc_new, void *loglFunc_arg_new,
			      GMRFLib_logl_tp * loglFunc_old, void *loglFunc_arg_old,
			      GMRFLib_graph_tp * graph,
			      GMRFLib_Qfunc_tp * Qfunc_new, void *Qfunc_arg_new,
			      GMRFLib_Qfunc_tp * Qfunc_old, void *Qfunc_arg_old,
			      GMRFLib_Qfunc_tp * Qfunc_old2new, void *Qfunc_arg_old2new,
			      GMRFLib_Qfunc_tp * Qfunc_new2old, void *Qfunc_arg_new2old,
			      GMRFLib_constr_tp * constr_new, GMRFLib_constr_tp * constr_old,
			      GMRFLib_optimize_param_tp * optpar, GMRFLib_blockupdate_param_tp * blockupdate_par, GMRFLib_store_tp * store)
{
	/*
	 * do a blockupdate, and return the proposed new state in 'x' and the corresponding log-acceptrate in 'laccept'
	 * 
	 * the density is
	 * 
	 * exp[ -1/2 (x-mean)'(Q+diag(c))(x-mean) + b'x + \sum d_i f(x_i) ]
	 * 
	 * where values are fixed if fixed_value[i] are true, and where a constraint, Ax=b can be spesified in 'constr' with
	 * optional noise
	 * 
	 * it's now ok that all values in fixed are 1, then none are updated, just the acceptrate is computed.
	 * 
	 */

	int n, i, free_block, id;
	double *mode = NULL, *bb = NULL, *cc = NULL, old2new, new2old, old, neww, *xx = NULL, *yy = NULL, logll;
	GMRFLib_problem_tp *problem = NULL;
	GMRFLib_blockupdate_param_tp *blockpar = NULL;

	GMRFLib_ASSERT(laccept, GMRFLib_EINVARG);
	GMRFLib_ASSERT(x_new, GMRFLib_EINVARG);
	GMRFLib_ASSERT(graph, GMRFLib_EINVARG);
	GMRFLib_ASSERT(Qfunc_new, GMRFLib_EINVARG);
	GMRFLib_ASSERT(Qfunc_old, GMRFLib_EINVARG);

	GMRFLib_ENTER_ROUTINE;

	id = GMRFLib_thread_id;

	if (constr_new) {
		GMRFLib_prepare_constr(constr_new, graph, 0);  /* no scaleing */
	}
	if (constr_old) {
		GMRFLib_prepare_constr(constr_old, graph, 0);  /* no scaleing */
	}

	/*
	 * use default choices if these are not specified. 
	 */
	if (!Qfunc_old2new) {
		Qfunc_old2new = Qfunc_new;
		if (!Qfunc_arg_old2new) {
			Qfunc_arg_old2new = Qfunc_arg_new;
		}
	}
	if (!Qfunc_new2old) {
		Qfunc_new2old = Qfunc_old;
		if (!Qfunc_arg_new2old) {
			Qfunc_arg_new2old = Qfunc_arg_old;
		}
	}

	if (blockupdate_par) {
		blockpar = blockupdate_par;
		free_block = 0;
	} else {
		GMRFLib_default_blockupdate_param(&blockpar);
		free_block = 1;
	}

	n = graph->n;
	mode = Calloc(n, double);
	xx = bb = Calloc(n, double);			       /* two names for the same storage */
	yy = cc = Calloc(n, double);			       /* two names for the same storage */

	/*
	 * fix storage accoring to reject or accept if use_more is ON 
	 */
	if (store && store->store_problems) {
		if (store->fixed_hyperparameters) {
			/*
			 * This is feature experimental for the moment. if the hyperparameters are the same, then keep both
			 * problems. by def, they are the same 
			 */
			if (store->decision == GMRFLib_STORE_ACCEPT) {
				if (store->old_logdens && store->new_logdens) {
					*(store->old_logdens) = *(store->new_logdens);
				}
			}
			if (store->decision == GMRFLib_STORE_REJECT) {
				/*
				 * nothing for the moment 
				 */
			}
		} else {
			if (store->decision == GMRFLib_STORE_ACCEPT) {
				if (store->problem_new2old)
					GMRFLib_free_problem(store->problem_new2old);
				store->problem_new2old = store->problem_old2new;
				store->problem_old2new = NULL;

				if (store->old_logdens && store->new_logdens) {
					*(store->old_logdens) = *(store->new_logdens);
				}
			}
			if (store->decision == GMRFLib_STORE_REJECT) {
				if (store->problem_old2new) {
					GMRFLib_free_problem(store->problem_old2new);
				}
				store->problem_old2new = NULL;
			}
		}
	}

	if (store && store->store_problems && store->fixed_hyperparameters && store->problem_old2new) {
		/*
		 * use copy 
		 */
		problem = store->problem_old2new;
	} else {
		/*
		 * first, find the point to expand around
		 * 
		 * note that i use Qfunc_old2new if we have constraints as the Qfunc_new can then be singular!
		 * 
		 * this step is always performed, not matter store->use_more 
		 */
		memcpy(mode, x_old, n * sizeof(double));
		if (blockpar->modeoption == GMRFLib_MODEOPTION_MODE && d_new) {
			GMRFLib_EWRAP1(GMRFLib_optimize_store(mode, b_new, c_new, mean_new, graph,
							      (constr_new ? Qfunc_old2new : Qfunc_new),
							      (constr_new ? Qfunc_arg_old2new : Qfunc_arg_new),
							      constr_new, d_new, loglFunc_new, loglFunc_arg_new, optpar, store));
		}

		/*
		 * compute the terms from loglFunc 
		 */
		if (d_new) {
#pragma omp parallel for private(i)
			for (i = 0; i < n; i++) {
				double cmin = 0.0;
				GMRFLib_thread_id = id;
				if (d_new[i]) {
					GMRFLib_2order_approx(NULL, &bb[i], &cc[i], NULL, d_new[i], mode[i], i,
							      mode, loglFunc_new, loglFunc_arg_new, &(blockpar->step_len), &(blockpar->stencil),
							      &cmin);
				}
			}
			GMRFLib_thread_id = id;
		}

		/*
		 * add the linear and quadratic term to the general model. note that we need to correct the b-term due to the mean. 
		 */
		if (b_new) {
			if (mean_new) {
				for (i = 0; i < n; i++) {
					bb[i] += b_new[i] - cc[i] * mean_new[i];
				}
			} else {
				for (i = 0; i < n; i++) {
					bb[i] += b_new[i];
				}
			}
		} else {
			if (mean_new) {
				for (i = 0; i < n; i++) {
					bb[i] -= cc[i] * mean_new[i];
				}
			}
		}

		if (c_new) {
			for (i = 0; i < n; i++) {
				cc[i] += c_new[i];
			}
		}

		GMRFLib_EWRAP1(GMRFLib_init_problem_store(&problem, x_old, bb, cc, mean_new, graph,
							  Qfunc_old2new, Qfunc_arg_old2new, constr_new, store));
	}

	if (problem) {
		GMRFLib_EWRAP1(GMRFLib_sample(problem));
		old2new = problem->sub_logdens;
		memcpy(x_new, problem->sample, n * sizeof(double));

		if (store && store->store_problems) {
			store->problem_old2new = problem;
		} else {
			GMRFLib_free_problem(problem);
		}
		problem = NULL;
	} else {
		/*
		 * nothing to do really, just make sure that the x_new equals x_old 
		 */
		old2new = 0.0;
		memcpy(x_new, x_old, n * sizeof(double));
	}

	/*
	 * now, go backwards
	 * 
	 * note: no need to optimize if there is no data.
	 * 
	 * note that i use Qfunc_new2old we have constraints as then the Qfunc_old can be singular! 
	 */

	if (store && store->store_problems && store->problem_new2old) {
		/*
		 * use copy 
		 */
		problem = store->problem_new2old;
	} else {
		memcpy(mode, x_new, n * sizeof(double));
		if (blockpar->modeoption == GMRFLib_MODEOPTION_MODE && d_old) {
			GMRFLib_EWRAP1(GMRFLib_optimize_store(mode, b_old, c_old, mean_old, graph,
							      (constr_old ? Qfunc_new2old : Qfunc_old),
							      (constr_old ? Qfunc_arg_new2old : Qfunc_arg_old),
							      constr_old, d_old, loglFunc_old, loglFunc_arg_old, optpar, store));
		}

		Memset(bb, 0, n * sizeof(double));
		Memset(cc, 0, n * sizeof(double));
		if (d_old) {
#pragma omp parallel for private(i)
			for (i = 0; i < n; i++) {
				double cmin = 0.0;
				GMRFLib_thread_id = id;
				if (d_old[i]) {
					GMRFLib_2order_approx(NULL, &bb[i], &cc[i], NULL, d_old[i], mode[i], i, mode,
							      loglFunc_old, loglFunc_arg_old, &(blockpar->step_len), &(blockpar->stencil), &cmin);
				}
			}
			GMRFLib_thread_id = id;
		}

		/*
		 * add the linear and quadratic term to the general model. note that we need to correct the b-term due to the
		 * mean. 
		 */
		if (b_old) {
			if (mean_old) {
				for (i = 0; i < n; i++) {
					bb[i] += b_old[i] - cc[i] * mean_old[i];
				}
			} else {
				for (i = 0; i < n; i++) {
					bb[i] += b_old[i];
				}
			}
		} else {
			if (mean_old) {
				for (i = 0; i < n; i++) {
					bb[i] -= cc[i] * mean_old[i];
				}
			}
		}
		if (c_old) {
			for (i = 0; i < n; i++) {
				cc[i] += c_old[i];
			}
		}

		GMRFLib_EWRAP1(GMRFLib_init_problem_store
			       (&problem, x_new, bb, cc, mean_old, graph, Qfunc_new2old, Qfunc_arg_new2old, constr_old, store));
	}

	if (problem) {
		memcpy(problem->sample, x_old, n * sizeof(double));
		GMRFLib_EWRAP1(GMRFLib_evaluate(problem));
		new2old = problem->sub_logdens;

		if (store && store->store_problems) {
			store->problem_new2old = problem;      /* just set it back [this is OK] */
		} else {
			GMRFLib_free_problem(problem);
		}
		problem = NULL;
	} else {
		/*
		 * nothing to do 
		 */
		new2old = 0.0;
	}

	bb = cc = NULL;					       /* not used anymore */

	/*
	 * compute the density at x and x_old.
	 * 
	 * FIXME: here i do not use subgraph, but that require big tests to see if the results would be the same. is it worth
	 * it really? 
	 */

	neww = 0.0;
	if (mean_new) {
		for (i = 0; i < n; i++) {
			xx[i] = x_new[i] - mean_new[i];
		}
	} else {
		memcpy(xx, x_new, n * sizeof(double));
	}
	GMRFLib_Qx(yy, xx, graph, Qfunc_new, Qfunc_arg_new);
	if (c_new) {
		for (i = 0; i < n; i++) {
			neww += yy[i] * xx[i] + c_new[i] * SQR(xx[i]);
		}
	} else {
		for (i = 0; i < n; i++) {
			neww += yy[i] * xx[i];
		}
	}
	neww *= -0.5;
	if (b_new) {
		for (i = 0; i < n; i++) {
			neww += b_new[i] * x_new[i];
		}
	}
	if (d_new) {
		double sum = 0.0;

#pragma omp parallel for private(i) reduction(+: sum)
		for (i = 0; i < n; i++) {
			GMRFLib_thread_id = id;
			if (d_new[i]) {
				loglFunc_new(&logll, &x_new[i], 1, i, x_new, NULL, loglFunc_arg_new);
				sum += d_new[i] * logll;
			}
		}
		GMRFLib_thread_id = id;
		neww += sum;
	}

	/*
	 * store the new value 
	 */
	if (store && store->store_problems) {
		if (!(store->new_logdens)) {
			store->new_logdens = Calloc(1, double);
		}
		*(store->new_logdens) = neww;
	}

	if (store && store->store_problems && store->old_logdens) {
		old = *(store->old_logdens);
	} else {
		old = 0.0;
		if (mean_old) {
			for (i = 0; i < n; i++) {
				xx[i] = x_old[i] - mean_old[i];
			}
		} else {
			memcpy(xx, x_old, n * sizeof(double));
		}
		GMRFLib_Qx(yy, xx, graph, Qfunc_old, Qfunc_arg_old);
		if (c_old) {
			for (i = 0; i < n; i++) {
				old += yy[i] * xx[i] + c_old[i] * SQR(xx[i]);
			}
		} else {
			for (i = 0; i < n; i++) {
				old += yy[i] * xx[i];
			}
		}
		old *= -0.5;
		if (b_old) {
			for (i = 0; i < n; i++) {
				old += b_old[i] * x_old[i];
			}
		}
		if (d_old) {
			double sum = 0.0;

#pragma omp parallel for private(i) reduction(+: sum)
			for (i = 0; i < n; i++) {
				GMRFLib_thread_id = id;
				if (d_old[i]) {
					loglFunc_old(&logll, &x_old[i], 1, i, x_old, NULL, loglFunc_arg_old);
					sum += d_old[i] * logll;
				}
			}
			GMRFLib_thread_id = id;
			old += sum;
		}
	}

	/*
	 * store the old value 
	 */
	if (store && store->store_problems) {
		if (!(store->old_logdens)) {
			store->old_logdens = Calloc(1, double);
		}
		*(store->old_logdens) = old;
	}

	/*
	 * finally.... 
	 */
	*laccept = neww + new2old - (old + old2new);

	if (blockpar->fp) {
		fprintf(blockpar->fp, "\n%s: laccept %f\n", __GMRFLib_FuncName, *laccept);
		fprintf(blockpar->fp, "\tnew_ldens %12.6f\n", neww);
		fprintf(blockpar->fp, "\told_ldens %12.6f\n", old);
		fprintf(blockpar->fp, "\tnew2old   %12.6f\n", new2old);
		fprintf(blockpar->fp, "\told2new   %12.6f\n", old2new);
	}

	Free(xx);
	Free(yy);
	Free(mode);
	if (free_block) {
		Free(blockpar);
	}

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_init_GMRF_approximation(GMRFLib_problem_tp ** problem, double *x, double *b, double *c, double *mean, double *d,
				    GMRFLib_logl_tp * loglFunc, void *loglFunc_arg,
				    GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg,
				    GMRFLib_constr_tp * constr, GMRFLib_optimize_param_tp * optpar, GMRFLib_blockupdate_param_tp * blockupdate_par)
{
	GMRFLib_ENTER_ROUTINE;
	GMRFLib_EWRAP1(GMRFLib_init_GMRF_approximation_store(problem, x, b, c, mean, d, loglFunc, loglFunc_arg,
							     graph, Qfunc, Qfunc_arg, constr, optpar, blockupdate_par, NULL));
	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_init_GMRF_approximation_store(GMRFLib_problem_tp ** problem, double *x, double *b, double *c, double *mean,
					  double *d, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg,
					  GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg,
					  GMRFLib_constr_tp * constr, GMRFLib_optimize_param_tp * optpar,
					  GMRFLib_blockupdate_param_tp * blockupdate_par, GMRFLib_store_tp * store)
{
	int i, j, free_x = 0, free_b = 0, free_c = 0, free_mean = 0, free_d = 0, free_blockpar = 0, n, id;
	double *bb = NULL, *cc = NULL, *mode = NULL;

#define FREE_ALL if (1) { if (free_x) Free(x); if (free_b) Free(b); if (free_c) Free(c); if (free_d) Free(d); \
	if (free_mean) Free(mean); if (free_blockpar) Free(blockupdate_par); Free(bb); Free(cc); Free(mode);}

	id = GMRFLib_thread_id;
	GMRFLib_ASSERT(problem, GMRFLib_EINVARG);
	GMRFLib_ASSERT(graph, GMRFLib_EINVARG);
	GMRFLib_ASSERT(Qfunc, GMRFLib_EINVARG);
	GMRFLib_ENTER_ROUTINE;

	n = graph->n;
	if (n == 0) {
		*problem = NULL;
		GMRFLib_LEAVE_ROUTINE;
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

	bb = Calloc(n, double);
	cc = Calloc(n, double);
	mode = Calloc(n, double);

	if (!blockupdate_par) {
		GMRFLib_default_blockupdate_param(&blockupdate_par);
		free_blockpar = 1;
	}

	memcpy(mode, x, n * sizeof(double));
	if (blockupdate_par->modeoption == GMRFLib_MODEOPTION_MODE && d)
		GMRFLib_EWRAP1(GMRFLib_optimize_store(mode, b, c, mean, graph, Qfunc, Qfunc_arg, constr, d, loglFunc, loglFunc_arg, optpar, store));

	if (!(store && store->sub_graph)) {
		/*
		 * compute the terms from loglFunc 
		 */
		if (d) {
#pragma omp parallel for private(i)
			for (i = 0; i < n; i++) {
				double cmin = 0.0;
				GMRFLib_thread_id = id;
				if (d[i]) {
					GMRFLib_2order_approx(NULL, &bb[i], &cc[i], NULL, d[i], mode[i], i, mode, loglFunc, loglFunc_arg,
							      &(blockupdate_par->step_len), &(blockupdate_par->stencil), &cmin);
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
	} else {
		/*
		 * do the same as above, but only for those 'i' which is not fixed. this is faster if we have fixed values and
		 * the sub_graph is available. 
		 */

		int ns = store->sub_graph->n, *mothergraph_idx = store->sub_graph->mothergraph_idx;

#pragma omp parallel for private(i, j)
		for (j = 0; j < ns; j++) {
			double cmin = 0.0;
			GMRFLib_thread_id = id;
			i = mothergraph_idx[j];
			if (d[i]) {
				GMRFLib_2order_approx(NULL, &bb[i], &cc[i], NULL, d[i], mode[i], i, mode, loglFunc, loglFunc_arg,
						      &(blockupdate_par->step_len), &(blockupdate_par->stencil), &cmin);
			}
		}
		GMRFLib_thread_id = id;

		if (b) {
			if (mean) {
				for (j = 0; j < ns; j++) {
					i = mothergraph_idx[j];
					bb[i] += b[i] - cc[i] * mean[i];
				}
			} else {
				for (j = 0; j < ns; j++) {
					i = mothergraph_idx[j];
					bb[i] += b[i];
				}
			}
		} else {
			if (mean) {
				for (j = 0; j < ns; j++) {
					i = mothergraph_idx[j];
					bb[i] -= cc[i] * mean[i];
				}
			}
		}

		if (c) {
			for (j = 0; j < ns; j++) {
				i = mothergraph_idx[j];
				cc[i] += c[i];
			}
		}
	}

	GMRFLib_EWRAP1(GMRFLib_init_problem_store(problem, x, bb, cc, mean, graph, Qfunc, Qfunc_arg, constr, store));

	FREE_ALL;
	GMRFLib_LEAVE_ROUTINE;

	return GMRFLib_SUCCESS;

#undef FREE_ALL
}

int GMRFLib_2order_taylor(double *a, double *b, double *c, double *dd, double d, double x0, int indx,
			  double *x_vec, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, double *step_len, int *stencil)
{
	/*
	 * compute a,b,c in the taylor expansion around x0 of d*loglFunc(x0,...)
	 * 
	 * a + b*(x-x0) + 0.5*c*(x-x0)^2 + 1/6*dd*(x-x0)^3
	 * 
	 */
	double f0 = 0.0, df = 0.0, ddf = 0.0, dddf = 0.0;

	if (ISZERO(d)) {
		f0 = df = ddf = 0.0;
	} else {
		GMRFLib_2order_approx_core(&f0, &df, &ddf, (dd ? &dddf : NULL), x0, indx, x_vec, loglFunc, loglFunc_arg, step_len, stencil);
	}

	if (a) {
		*a = d * f0;
	}
	if (b) {
		*b = d * df;
	}
	if (c) {
		*c = d * ddf;
	}
	if (dd) {
		*dd = d * dddf;
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_2order_approx(double *a, double *b, double *c, double *dd, double d, double x0, int indx,
			  double *x_vec, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, double *step_len, int *stencil, double *cmin)
{
	/*
	 * compute a,b,c in the taylor expansion around x0 of d*loglFunc(x0,...)
	 * 
	 * a + b*x- 0.5*c*x^2 + 1/6*dd*x^3
	 *
	 * where cmin is the minimum value of c.
	 */

	/*
	 * > A:=collect(expand(a + b * (x-x0) + 1/2 \                                      
	 * > * c * (x-x0)^2 + 1/6 * d * (x-x0)^3), [x,x^2, x^3]);
	 *          3                    /            2    \             2              3
	 *       d x   /      d x0\  2   |        d x0     |         c x0           d x0
	 * A := ---- + |c/2 - ----| x  + |-c x0 + ----- + b| x + a + ----- - b x0 - -----
	 *             \       2  /      \          2      /           2              6
	 *
	 * > coeff(A,x);                                                                   
	 *             2
	 *         d x0
	 * -c x0 + ----- + b
	 *           2
	 * 
	 * > coeff(A,x^2);
	 *       d x0
	 * c/2 - ----
	 *        2
	 * 
	 * > coeff(A,x^3);
	 * d/6
	 * 
	 */

	double f0 = 0.0, df = 0.0, ddf = 0.0, dddf = 0.0;

	if (ISZERO(d)) {
		f0 = df = ddf = dddf = 0.0;
	} else {
		GMRFLib_2order_approx_core(&f0, &df, &ddf, (dd ? &dddf : NULL), x0, indx, x_vec, loglFunc, loglFunc_arg, step_len, stencil);
	}

	// If there is a truncation, we have to do this here
	if (cmin) {
		ddf = DMIN(-(*cmin), ddf);
	}

	if (a) {
		*a = d * (f0 - df * x0 + 0.5 * ddf * SQR(x0) + 1.0 / 6.0 * dddf * gsl_pow_3(x0));
	}
	if (b) {
		*b = d * (df - ddf * x0 + dddf / 2.0 * SQR(x0));
	}
	if (c) {
		*c = -d * (ddf - dddf * x0);
	}
	if (dd) {
		*dd = d * dddf;
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_2order_approx_core(double *a, double *b, double *c, double *dd, double x0, int indx,
			       double *x_vec, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, double *step_len, int *stencil)
{

#define ERR if (dd) {							\
		fprintf(stderr, "2order_approx_core: 3rd derivative requested but there are to few points in the stencil\n"); \
		assert(dd == NULL);					\
		exit(1);						\
	}

	double step, df = 0.0, ddf = 0.0, dddf = 0.0, xx[9], f[9], f0 = 0.0;
	int code = loglFunc(f, &x0, 0, indx, x_vec, NULL, loglFunc_arg);

	if (step_len && *step_len < 0.0) {
		/*
		 * for internal use only! 
		 */
		step = -(*step_len);

		xx[0] = x0 - 2 * step;
		xx[1] = x0 - step;
		xx[2] = x0;
		xx[3] = x0 + step;
		xx[4] = x0 + 2 * step;

		loglFunc(&f[0], &xx[0], 1, indx, x_vec, NULL, loglFunc_arg);
		loglFunc(&f[1], &xx[1], 1, indx, x_vec, NULL, loglFunc_arg);
		loglFunc(&f[2], &xx[2], 1, indx, x_vec, NULL, loglFunc_arg);
		loglFunc(&f[3], &xx[3], 1, indx, x_vec, NULL, loglFunc_arg);
		loglFunc(&f[4], &xx[4], 1, indx, x_vec, NULL, loglFunc_arg);

		f0 = f[2];
		df = (1.0 / 12.0 * f[4] - 2.0 / 3.0 * f[3] + 0.0 * f[2] + 2.0 / 3.0 * f[1] - 1.0 / 12.0 * f[0]) / step;
		ddf = (-1.0 / 12.0 * f[4] - 4.0 / 3.0 * f[3] - 5.0 / 2.0 * f[2] + 4.0 / 3.0 * f[1] - 1.0 / 12.0 * f[0]) / SQR(step);
		dddf = (-1.0 / 2.0 * f[4] + 1.0 * f[3] + 0.0 * f[2] - 1.0 * f[1] + 1.0 / 2.0 * f[0]) / gsl_pow_3(step);

	} else if (code == GMRFLib_LOGL_COMPUTE_DERIVATIES || code == GMRFLib_LOGL_COMPUTE_DERIVATIES_AND_CDF) {
		/*
		 * this tells that exact calculations is carried out in loglFunc 
		 */

		/*
		 * this indicate that exact calculations can and are carried out in loglFunc 
		 */
		xx[0] = xx[1] = xx[2] = x0;
		loglFunc(f, xx, 3, indx, x_vec, NULL, loglFunc_arg);

		f0 = f[0];
		df = f[1];
		ddf = f[2];
		ERR;
	} else {
		int num_points = (stencil ? *stencil : 5);
		step = (step_len && *step_len > 0.0 ? *step_len : GMRFLib_eps(1.0 / 3.9134));
		switch (num_points) {
			/*
			 * see https://en.wikipedia.org/wiki/Finite_difference_coefficients
			 */
		case 3:
		{
			xx[0] = x0 - step;
			xx[1] = x0;
			xx[2] = x0 + step;

			loglFunc(f, xx, 3, indx, x_vec, NULL, loglFunc_arg);
			f0 = f[1];
			df = 0.5 * (f[2] - f[0]) / step;
			ddf = (f[2] - 2.0 * f[1] + f[0]) / SQR(step);
			ERR;
			break;
		}

		case 5:
		{
			double wf[] = { 1.0 / 12.0, -2.0 / 3.0, 0.0, 2.0 / 3.0, -1.0 / 12.0 };
			double wff[] = { -1.0 / 12.0, 4.0 / 3.0, -5.0 / 2.0, 4.0 / 3.0, -1.0 / 12.0 };
			double wfff[] = { -1.0 / 2.0, 1.0, 0.0, -1.0, 1.0 / 2.0 };


			xx[0] = x0 - 2.0 * step;
			xx[1] = x0 - step;
			xx[2] = x0;
			xx[3] = x0 + step;
			xx[4] = x0 + 2.0 * step;

			loglFunc(f, xx, 5, indx, x_vec, NULL, loglFunc_arg);
			f0 = f[2];
			df = (wf[0] * f[0] + wf[1] * f[1] + wf[2] * f[2] + wf[3] * f[3] + wf[4] * f[4]) / step;
			ddf = (wff[0] * f[0] + wff[1] * f[1] + wff[2] * f[2] + wff[3] * f[3] + wff[4] * f[4]) / SQR(step);
			if (dd) {
				dddf = (wfff[0] * f[0] + wfff[1] * f[1] + wfff[2] * f[2] + wfff[3] * f[3] + wfff[4] * f[4]) / gsl_pow_3(step);
			}
			break;
		}

		case 7:
		{
			double wf[] = { -1.0 / 60.0, 3.0 / 20.0, -3.0 / 4.0, 0.0, 3.0 / 4.0, -3.0 / 20.0, 1.0 / 60.0 };
			double wff[] = { 1.0 / 90.0, -3.0 / 20.0, 3.0 / 2.0, -49.0 / 18.0, 3.0 / 2.0, -3.0 / 20.0, 1.0 / 90.0 };
			double wfff[] = { 1.0 / 8.0, -1.0, 13.0 / 8.0, 0.0, -13.0 / 8.0, 1.0, -1.0 / 8.0 };

			xx[0] = x0 - 3.0 * step;
			xx[1] = x0 - 2.0 * step;
			xx[2] = x0 - step;
			xx[3] = x0;
			xx[4] = x0 + step;
			xx[5] = x0 + 2.0 * step;
			xx[6] = x0 + 3.0 * step;

			loglFunc(f, xx, 7, indx, x_vec, NULL, loglFunc_arg);
			f0 = f[3];
			df = (wf[0] * f[0] + wf[1] * f[1] + wf[2] * f[2] + wf[3] * f[3] + wf[4] * f[4] + wf[5] * f[5] + wf[6] * f[6]) / step;
			ddf = (wff[0] * f[0] + wff[1] * f[1] + wff[2] * f[2] + wff[3] * f[3] + wff[4] * f[4] + wff[5] * f[5] +
			       wff[6] * f[6]) / SQR(step);
			if (dd) {
				dddf = (wfff[0] * f[0] + wfff[1] * f[1] + wfff[2] * f[2] + wfff[3] * f[3] + wfff[4] * f[4] + wfff[5] * f[5] +
					wfff[6] * f[6]) / gsl_pow_3(step);
			}
			break;
		}

		case 9:
		{
			double wf[] = { 1.0 / 280.0, -4.0 / 105.0, 1.0 / 5.0, -4.0 / 5.0, 0.0, 4.0 / 5.0, -1.0 / 5.0, 4.0 / 105.0, -1.0 / 280.0 };
			double wff[] =
			    { -1.0 / 560.0, 8.0 / 315.0, -1.0 / 5.0, 8.0 / 5.0, -205.0 / 72.0, 8.0 / 5.0, -1.0 / 5.0, 8.0 / 315.0, -1.0 / 560.0 };

			double wfff[] =
			    { -7.0 / 240.0, 3.0 / 10.0, -169.0 / 120.0, 61.0 / 30.0, 0.0, -61.0 / 30.0, 169.0 / 120.0, -3.0 / 10.0, 7.0 / 240.0 };

			xx[0] = x0 - 4.0 * step;
			xx[1] = x0 - 3.0 * step;
			xx[2] = x0 - 2.0 * step;
			xx[3] = x0 - step;
			xx[4] = x0;
			xx[5] = x0 + step;
			xx[6] = x0 + 2.0 * step;
			xx[7] = x0 + 3.0 * step;
			xx[8] = x0 + 4.0 * step;

			loglFunc(f, xx, 9, indx, x_vec, NULL, loglFunc_arg);
			f0 = f[4];
			df = (wf[0] * f[0] + wf[1] * f[1] + wf[2] * f[2] + wf[3] * f[3] + wf[4] * f[4] + wf[5] * f[5] + wf[6] * f[6] +
			      wf[7] * f[7] + wf[8] * f[8]) / step;
			ddf = (wff[0] * f[0] + wff[1] * f[1] + wff[2] * f[2] + wff[3] * f[3] + wff[4] * f[4] + wff[5] * f[5] + wff[6] * f[6] +
			       wff[7] * f[7] + wff[8] * f[8]) / SQR(step);
			if (dd) {
				dddf = (wfff[0] * f[0] + wfff[1] * f[1] + wfff[2] * f[2] + wfff[3] * f[3] + wfff[4] * f[4] + wfff[5] * f[5] +
					wfff[6] * f[6] + wfff[7] * f[7] + wfff[8] * f[8]) / gsl_pow_3(step);
			}
			break;
		}

		default:
			GMRFLib_ASSERT(num_points == 3 || num_points == 5 || num_points == 7 || num_points == 9, GMRFLib_EINVARG);
			abort();
		}
	}
	*a = f0;
	*b = df;
	*c = ddf;
	if (dd) {
		*dd = dddf;
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_blockupdate_hidden(double *laccept,
			       double *x_new, double *x_old,
			       double *b_new, double *b_old,
			       double *c_new, double *c_old,
			       double *mean_new, double *mean_old,
			       double *d_new, double *d_old,
			       GMRFLib_logl_tp * loglFunc_new, void *loglFunc_arg_new,
			       GMRFLib_logl_tp * loglFunc_old, void *loglFunc_arg_old,
			       GMRFLib_graph_tp * graph,
			       GMRFLib_Qfunc_tp * Qfunc_new, void *Qfunc_arg_new,
			       GMRFLib_Qfunc_tp * Qfunc_old, void *Qfunc_arg_old,
			       GMRFLib_Qfunc_tp * Qfunc_old2new, void *Qfunc_arg_old2new,
			       GMRFLib_Qfunc_tp * Qfunc_new2old, void *Qfunc_arg_new2old, GMRFLib_optimize_param_tp * optpar,
			       GMRFLib_hidden_param_tp * hidden_par)
{
	GMRFLib_ENTER_ROUTINE;
	GMRFLib_EWRAP1(GMRFLib_blockupdate_hidden_store(laccept, x_new, x_old, b_new, b_old, c_new, c_old, mean_new, mean_old,
							d_new, d_old, loglFunc_new, loglFunc_arg_new, loglFunc_old,
							loglFunc_arg_old, graph, Qfunc_new, Qfunc_arg_new,
							Qfunc_old, Qfunc_arg_old, Qfunc_old2new, Qfunc_arg_old2new,
							Qfunc_new2old, Qfunc_arg_new2old, optpar, hidden_par, NULL));
	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_blockupdate_hidden_store(double *laccept,
				     double *x_new, double *x_old,
				     double *b_new, double *b_old,
				     double *c_new, double *c_old,
				     double *mean_new, double *mean_old,
				     double *d_new, double *d_old,
				     GMRFLib_logl_tp * loglFunc_new, void *loglFunc_arg_new,
				     GMRFLib_logl_tp * loglFunc_old, void *loglFunc_arg_old,
				     GMRFLib_graph_tp * graph,
				     GMRFLib_Qfunc_tp * Qfunc_new, void *Qfunc_arg_new,
				     GMRFLib_Qfunc_tp * Qfunc_old, void *Qfunc_arg_old,
				     GMRFLib_Qfunc_tp * Qfunc_old2new, void *Qfunc_arg_old2new,
				     GMRFLib_Qfunc_tp * Qfunc_new2old, void *Qfunc_arg_new2old,
				     GMRFLib_optimize_param_tp * optpar, GMRFLib_hidden_param_tp * hidden_par, GMRFLib_store_tp * store)
{
	/*
	 * do a blockupdate, and return the proposed new state in 'x' and the corresponding log-acceptrate in 'laccept'
	 * 
	 * the density is
	 * 
	 * exp[ -1/2 (x-mean)'(Q+diag(c))(x-mean) + b'x + \sum d_i f(x_i) ]
	 * 
	 * where values are fixed if fixed_value[i] are true, and where a constraint, Ax=b can be spesified in 'constr'
	 * 
	 */

	int n, i;
	double *mode = NULL, old2new, new2old, old, neww, *xx = NULL, *yy = NULL, logll;
	GMRFLib_hidden_problem_tp *hidden_problem = NULL;

	GMRFLib_ASSERT(laccept, GMRFLib_EINVARG);
	GMRFLib_ASSERT(x_new, GMRFLib_EINVARG);
	GMRFLib_ASSERT(graph, GMRFLib_EINVARG);
	GMRFLib_ASSERT(Qfunc_new, GMRFLib_EINVARG);
	GMRFLib_ASSERT(Qfunc_old, GMRFLib_EINVARG);

	GMRFLib_ENTER_ROUTINE;

	if (!Qfunc_old2new) {
		Qfunc_old2new = Qfunc_new;
	}
	if (!Qfunc_arg_old2new) {
		Qfunc_arg_old2new = Qfunc_arg_new;
	}
	if (!Qfunc_new2old) {
		Qfunc_new2old = Qfunc_old;
	}
	if (!Qfunc_arg_new2old) {
		Qfunc_arg_new2old = Qfunc_arg_old;
	}

	n = graph->n;
	mode = Calloc(n, double);
	xx = Calloc(n, double);				       /* two names for the same storage */
	yy = Calloc(n, double);				       /* two names for the same storage */

	GMRFLib_EWRAP1(GMRFLib_init_problem_hidden_store(&hidden_problem,
							 x_old, b_new, c_new, mean_new, graph, Qfunc_old2new, Qfunc_arg_old2new,
							 d_new, loglFunc_new, loglFunc_arg_new, optpar, hidden_par, store));
	GMRFLib_EWRAP1(GMRFLib_sample_hidden(hidden_problem));
	old2new = hidden_problem->sub_logdens;
	memcpy(x_new, hidden_problem->sample, n * sizeof(double));

	GMRFLib_free_hidden(hidden_problem);
	hidden_problem = NULL;

	GMRFLib_EWRAP1(GMRFLib_init_problem_hidden_store(&hidden_problem,
							 x_new, b_old, c_old, mean_old, graph, Qfunc_new2old, Qfunc_arg_new2old,
							 d_old, loglFunc_old, loglFunc_arg_old, optpar, hidden_par, store));
	memcpy(hidden_problem->sample, x_old, n * sizeof(double));
	GMRFLib_EWRAP1(GMRFLib_evaluate_hidden(hidden_problem));
	new2old = hidden_problem->sub_logdens;
	GMRFLib_free_hidden(hidden_problem);
	hidden_problem = NULL;

	/*
	 * compute the density at x and x_old.
	 * 
	 * FIXME: here i do not use subgraph, but that require big tests to see if the results would be the same. is it worth
	 * it really? 
	 */
	neww = 0.0;
	if (mean_new) {
		for (i = 0; i < n; i++) {
			xx[i] = x_new[i] - mean_new[i];
		}
	} else {
		memcpy(xx, x_new, n * sizeof(double));
	}
	GMRFLib_Qx(yy, xx, graph, Qfunc_new, Qfunc_arg_new);
	if (c_new) {
		for (i = 0; i < n; i++) {
			neww += yy[i] * xx[i] + c_new[i] * SQR(xx[i]);
		}
	} else {
		for (i = 0; i < n; i++) {
			neww += yy[i] * xx[i];
		}
	}
	neww *= -0.5;
	if (b_new) {
		for (i = 0; i < n; i++) {
			neww += b_new[i] * x_new[i];
		}
	}
	if (d_new) {
		for (i = 0; i < n; i++) {
			if (d_new[i]) {
				loglFunc_new(&logll, &x_new[i], 1, i, x_new, NULL, loglFunc_arg_new);
				neww += d_new[i] * logll;
			}
		}
	}

	old = 0.0;
	if (mean_old) {
		for (i = 0; i < n; i++) {
			xx[i] = x_old[i] - mean_old[i];
		}
	} else {
		memcpy(xx, x_old, n * sizeof(double));
	}
	GMRFLib_Qx(yy, xx, graph, Qfunc_old, Qfunc_arg_old);
	if (c_old) {
		for (i = 0; i < n; i++) {
			old += yy[i] * xx[i] + c_old[i] * SQR(xx[i]);
		}
	} else {
		for (i = 0; i < n; i++) {
			old += yy[i] * xx[i];
		}
	}
	old *= -0.5;
	if (b_old) {
		for (i = 0; i < n; i++) {
			old += b_old[i] * x_old[i];
		}
	}
	if (d_old) {
		for (i = 0; i < n; i++) {
			if (d_old[i]) {
				loglFunc_old(&logll, &x_old[i], 1, i, x_old, NULL, loglFunc_arg_old);
				old += d_old[i] * logll;
			}
		}
	}

	*laccept = neww + new2old - (old + old2new);

	if (0) {					       /* FIXME */
		fprintf(stdout, "\n%s: laccept %f\n", __GMRFLib_FuncName, *laccept);
		fprintf(stdout, "\tnew_ldens %12.6f\n", neww);
		fprintf(stdout, "\told_ldens %12.6f\n", old);
		fprintf(stdout, "\tnew2old   %12.6f\n", new2old);
		fprintf(stdout, "\told2new   %12.6f\n", old2new);
	}

	Free(xx);
	Free(yy);
	Free(mode);

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_blockupdate(double *laccept,
			double *x_new, double *x_old,
			double *b_new, double *b_old,
			double *c_new, double *c_old,
			double *mean_new, double *mean_old,
			double *d_new, double *d_old,
			GMRFLib_logl_tp * loglFunc_new, void *loglFunc_arg_new,
			GMRFLib_logl_tp * loglFunc_old, void *loglFunc_arg_old,
			GMRFLib_graph_tp * graph,
			GMRFLib_Qfunc_tp * Qfunc_new, void *Qfunc_arg_new,
			GMRFLib_Qfunc_tp * Qfunc_old, void *Qfunc_arg_old,
			GMRFLib_Qfunc_tp * Qfunc_old2new, void *Qfunc_arg_old2new,
			GMRFLib_Qfunc_tp * Qfunc_new2old, void *Qfunc_arg_new2old,
			GMRFLib_constr_tp * constr_new, GMRFLib_constr_tp * constr_old,
			GMRFLib_optimize_param_tp * optpar, GMRFLib_blockupdate_param_tp * blockupdate_par)
{
	GMRFLib_ENTER_ROUTINE;
	GMRFLib_EWRAP1(GMRFLib_blockupdate_store
		       (laccept, x_new, x_old, b_new, b_old, c_new, c_old, mean_new, mean_old, d_new, d_old, loglFunc_new,
			loglFunc_arg_new, loglFunc_old, loglFunc_arg_old, graph, Qfunc_new, Qfunc_arg_new,
			Qfunc_old, Qfunc_arg_old, Qfunc_old2new, Qfunc_arg_old2new, Qfunc_new2old, Qfunc_arg_new2old,
			constr_new, constr_old, optpar, blockupdate_par, NULL));
	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_blockupdate_store(double *laccept,
			      double *x_new, double *x_old,
			      double *b_new, double *b_old,
			      double *c_new, double *c_old,
			      double *mean_new, double *mean_old,
			      double *d_new, double *d_old,
			      GMRFLib_logl_tp * loglFunc_new, void *loglFunc_arg_new,
			      GMRFLib_logl_tp * loglFunc_old, void *loglFunc_arg_old,
			      GMRFLib_graph_tp * graph,
			      GMRFLib_Qfunc_tp * Qfunc_new, void *Qfunc_arg_new,
			      GMRFLib_Qfunc_tp * Qfunc_old, void *Qfunc_arg_old,
			      GMRFLib_Qfunc_tp * Qfunc_old2new, void *Qfunc_arg_old2new,
			      GMRFLib_Qfunc_tp * Qfunc_new2old, void *Qfunc_arg_new2old,
			      GMRFLib_constr_tp * constr_new, GMRFLib_constr_tp * constr_old,
			      GMRFLib_optimize_param_tp * optpar, GMRFLib_blockupdate_param_tp * blockupdate_par, GMRFLib_store_tp * store)
{
	/*
	 * do a blockupdate, and return the proposed new state in 'x' and the corresponding log-acceptrate in 'laccept'
	 * 
	 * the density is
	 * 
	 * exp[ -1/2 (x-mean)'(Q+diag(c))(x-mean) + b'x + \sum d_i f(x_i) ]
	 * 
	 * where values are fixed if fixed_value[i] are true, and where a constraint, Ax=b can be spesified in 'constr' with
	 * optional noise
	 * 
	 * it's now ok that all values in fixed are 1, then none are updated, just the acceptrate is computed.
	 * 
	 */

	int n, i, free_block, id;
	double *mode = NULL, *bb = NULL, *cc = NULL, old2new, new2old, old, neww, *xx = NULL, *yy = NULL, logll;
	GMRFLib_problem_tp *problem = NULL;
	GMRFLib_blockupdate_param_tp *blockpar = NULL;

	GMRFLib_ASSERT(laccept, GMRFLib_EINVARG);
	GMRFLib_ASSERT(x_new, GMRFLib_EINVARG);
	GMRFLib_ASSERT(graph, GMRFLib_EINVARG);
	GMRFLib_ASSERT(Qfunc_new, GMRFLib_EINVARG);
	GMRFLib_ASSERT(Qfunc_old, GMRFLib_EINVARG);

	GMRFLib_ENTER_ROUTINE;

	id = GMRFLib_thread_id;

	if (constr_new) {
		GMRFLib_prepare_constr(constr_new, graph, 0);  /* no scaleing */
	}
	if (constr_old) {
		GMRFLib_prepare_constr(constr_old, graph, 0);  /* no scaleing */
	}

	/*
	 * use default choices if these are not specified. 
	 */
	if (!Qfunc_old2new) {
		Qfunc_old2new = Qfunc_new;
		if (!Qfunc_arg_old2new) {
			Qfunc_arg_old2new = Qfunc_arg_new;
		}
	}
	if (!Qfunc_new2old) {
		Qfunc_new2old = Qfunc_old;
		if (!Qfunc_arg_new2old) {
			Qfunc_arg_new2old = Qfunc_arg_old;
		}
	}

	if (blockupdate_par) {
		blockpar = blockupdate_par;
		free_block = 0;
	} else {
		GMRFLib_default_blockupdate_param(&blockpar);
		free_block = 1;
	}

	n = graph->n;
	mode = Calloc(n, double);
	xx = bb = Calloc(n, double);			       /* two names for the same storage */
	yy = cc = Calloc(n, double);			       /* two names for the same storage */

	/*
	 * fix storage accoring to reject or accept if use_more is ON 
	 */
	if (store && store->store_problems) {
		if (store->fixed_hyperparameters) {
			/*
			 * This is feature experimental for the moment. if the hyperparameters are the same, then keep both
			 * problems. by def, they are the same 
			 */
			if (store->decision == GMRFLib_STORE_ACCEPT) {
				if (store->old_logdens && store->new_logdens) {
					*(store->old_logdens) = *(store->new_logdens);
				}
			}
			if (store->decision == GMRFLib_STORE_REJECT) {
				/*
				 * nothing for the moment 
				 */
			}
		} else {
			if (store->decision == GMRFLib_STORE_ACCEPT) {
				if (store->problem_new2old)
					GMRFLib_free_problem(store->problem_new2old);
				store->problem_new2old = store->problem_old2new;
				store->problem_old2new = NULL;

				if (store->old_logdens && store->new_logdens) {
					*(store->old_logdens) = *(store->new_logdens);
				}
			}
			if (store->decision == GMRFLib_STORE_REJECT) {
				if (store->problem_old2new) {
					GMRFLib_free_problem(store->problem_old2new);
				}
				store->problem_old2new = NULL;
			}
		}
	}

	if (store && store->store_problems && store->fixed_hyperparameters && store->problem_old2new) {
		/*
		 * use copy 
		 */
		problem = store->problem_old2new;
	} else {
		/*
		 * first, find the point to expand around
		 * 
		 * note that i use Qfunc_old2new if we have constraints as the Qfunc_new can then be singular!
		 * 
		 * this step is always performed, not matter store->use_more 
		 */
		Memcpy(mode, x_old, n * sizeof(double));
		if (blockpar->modeoption == GMRFLib_MODEOPTION_MODE && d_new) {
			GMRFLib_EWRAP1(GMRFLib_optimize_store(mode, b_new, c_new, mean_new, graph,
							      (constr_new ? Qfunc_old2new : Qfunc_new),
							      (constr_new ? Qfunc_arg_old2new : Qfunc_arg_new),
							      constr_new, d_new, loglFunc_new, loglFunc_arg_new, optpar, store));
		}

		/*
		 * compute the terms from loglFunc 
		 */
		if (d_new) {
#pragma omp parallel for private(i)
			for (i = 0; i < n; i++) {
				GMRFLib_thread_id = id;
				double cmin = 0.0;
				if (d_new[i]) {
					GMRFLib_2order_approx(NULL, &bb[i], &cc[i], NULL, d_new[i], mode[i], i,
							      mode, loglFunc_new, loglFunc_arg_new, &(blockpar->step_len), &(blockpar->stencil),
							      &cmin);
				}
			}
			GMRFLib_thread_id = id;
		}

		/*
		 * add the linear and quadratic term to the general model. note that we need to correct the b-term due to the mean. 
		 */
		if (b_new) {
			if (mean_new) {
				for (i = 0; i < n; i++) {
					bb[i] += b_new[i] - cc[i] * mean_new[i];
				}
			} else {
				for (i = 0; i < n; i++) {
					bb[i] += b_new[i];
				}
			}
		} else {
			if (mean_new) {
				for (i = 0; i < n; i++) {
					bb[i] -= cc[i] * mean_new[i];
				}
			}
		}

		if (c_new) {
			for (i = 0; i < n; i++) {
				cc[i] += c_new[i];
			}
		}

		GMRFLib_EWRAP1(GMRFLib_init_problem_store(&problem, x_old, bb, cc, mean_new, graph,
							  Qfunc_old2new, Qfunc_arg_old2new, constr_new, store));
	}

	if (problem) {
		GMRFLib_EWRAP1(GMRFLib_sample(problem));
		old2new = problem->sub_logdens;
		Memcpy(x_new, problem->sample, n * sizeof(double));

		if (store && store->store_problems) {
			store->problem_old2new = problem;
		} else {
			GMRFLib_free_problem(problem);
		}
		problem = NULL;
	} else {
		/*
		 * nothing to do really, just make sure that the x_new equals x_old 
		 */
		old2new = 0.0;
		Memcpy(x_new, x_old, n * sizeof(double));
	}

	/*
	 * now, go backwards
	 * 
	 * note: no need to optimize if there is no data.
	 * 
	 * note that i use Qfunc_new2old we have constraints as then the Qfunc_old can be singular! 
	 */

	if (store && store->store_problems && store->problem_new2old) {
		/*
		 * use copy 
		 */
		problem = store->problem_new2old;
	} else {
		Memcpy(mode, x_new, n * sizeof(double));
		if (blockpar->modeoption == GMRFLib_MODEOPTION_MODE && d_old) {
			GMRFLib_EWRAP1(GMRFLib_optimize_store(mode, b_old, c_old, mean_old, graph,
							      (constr_old ? Qfunc_new2old : Qfunc_old),
							      (constr_old ? Qfunc_arg_new2old : Qfunc_arg_old),
							      constr_old, d_old, loglFunc_old, loglFunc_arg_old, optpar, store));
		}

		Memset(bb, 0, n * sizeof(double));
		Memset(cc, 0, n * sizeof(double));
		if (d_old) {
#pragma omp parallel for private(i)
			for (i = 0; i < n; i++) {
				GMRFLib_thread_id = id;
				double cmin = 0.0;
				if (d_old[i]) {
					GMRFLib_2order_approx(NULL, &bb[i], &cc[i], NULL, d_old[i], mode[i], i, mode,
							      loglFunc_old, loglFunc_arg_old, &(blockpar->step_len), &(blockpar->stencil), &cmin);
				}
			}
			GMRFLib_thread_id = id;
		}

		/*
		 * add the linear and quadratic term to the general model. note that we need to correct the b-term due to the
		 * mean. 
		 */
		if (b_old) {
			if (mean_old) {
				for (i = 0; i < n; i++) {
					bb[i] += b_old[i] - cc[i] * mean_old[i];
				}
			} else {
				for (i = 0; i < n; i++) {
					bb[i] += b_old[i];
				}
			}
		} else {
			if (mean_old) {
				for (i = 0; i < n; i++) {
					bb[i] -= cc[i] * mean_old[i];
				}
			}
		}
		if (c_old) {
			for (i = 0; i < n; i++) {
				cc[i] += c_old[i];
			}
		}

		GMRFLib_EWRAP1(GMRFLib_init_problem_store
			       (&problem, x_new, bb, cc, mean_old, graph, Qfunc_new2old, Qfunc_arg_new2old, constr_old, store));
	}

	if (problem) {
		Memcpy(problem->sample, x_old, n * sizeof(double));
		GMRFLib_EWRAP1(GMRFLib_evaluate(problem));
		new2old = problem->sub_logdens;

		if (store && store->store_problems) {
			store->problem_new2old = problem;      /* just set it back [this is OK] */
		} else {
			GMRFLib_free_problem(problem);
		}
		problem = NULL;
	} else {
		/*
		 * nothing to do 
		 */
		new2old = 0.0;
	}

	bb = cc = NULL;					       /* not used anymore */

	/*
	 * compute the density at x and x_old.
	 * 
	 * FIXME: here i do not use subgraph, but that require big tests to see if the results would be the same. is it worth
	 * it really? 
	 */

	neww = 0.0;
	if (mean_new) {
		for (i = 0; i < n; i++) {
			xx[i] = x_new[i] - mean_new[i];
		}
	} else {
		Memcpy(xx, x_new, n * sizeof(double));
	}
	GMRFLib_Qx(yy, xx, graph, Qfunc_new, Qfunc_arg_new);
	if (c_new) {
		for (i = 0; i < n; i++) {
			neww += yy[i] * xx[i] + c_new[i] * SQR(xx[i]);
		}
	} else {
		for (i = 0; i < n; i++) {
			neww += yy[i] * xx[i];
		}
	}
	neww *= -0.5;
	if (b_new) {
		for (i = 0; i < n; i++) {
			neww += b_new[i] * x_new[i];
		}
	}
	if (d_new) {
		double sum = 0.0;

#pragma omp parallel for private(i) reduction(+: sum)
		for (i = 0; i < n; i++) {
			GMRFLib_thread_id = id;
			if (d_new[i]) {
				loglFunc_new(&logll, &x_new[i], 1, i, x_new, NULL, loglFunc_arg_new);
				sum += d_new[i] * logll;
			}
		}
		GMRFLib_thread_id = id;
		neww += sum;
	}

	/*
	 * store the new value 
	 */
	if (store && store->store_problems) {
		if (!(store->new_logdens)) {
			store->new_logdens = Calloc(1, double);
		}
		*(store->new_logdens) = neww;
	}

	if (store && store->store_problems && store->old_logdens) {
		old = *(store->old_logdens);
	} else {
		old = 0.0;
		if (mean_old) {
			for (i = 0; i < n; i++) {
				xx[i] = x_old[i] - mean_old[i];
			}
		} else {
			Memcpy(xx, x_old, n * sizeof(double));
		}
		GMRFLib_Qx(yy, xx, graph, Qfunc_old, Qfunc_arg_old);
		if (c_old) {
			for (i = 0; i < n; i++) {
				old += yy[i] * xx[i] + c_old[i] * SQR(xx[i]);
			}
		} else {
			for (i = 0; i < n; i++) {
				old += yy[i] * xx[i];
			}
		}
		old *= -0.5;
		if (b_old) {
			for (i = 0; i < n; i++) {
				old += b_old[i] * x_old[i];
			}
		}
		if (d_old) {
			double sum = 0.0;

#pragma omp parallel for private(i) reduction(+: sum)
			for (i = 0; i < n; i++) {
				GMRFLib_thread_id = id;
				if (d_old[i]) {
					loglFunc_old(&logll, &x_old[i], 1, i, x_old, NULL, loglFunc_arg_old);
					sum += d_old[i] * logll;
				}
			}
			GMRFLib_thread_id = id;
			old += sum;
		}
	}

	/*
	 * store the old value 
	 */
	if (store && store->store_problems) {
		if (!(store->old_logdens)) {
			store->old_logdens = Calloc(1, double);
		}
		*(store->old_logdens) = old;
	}

	/*
	 * finally.... 
	 */
	*laccept = neww + new2old - (old + old2new);

	if (blockpar->fp) {
		fprintf(blockpar->fp, "\n%s: laccept %f\n", __GMRFLib_FuncName, *laccept);
		fprintf(blockpar->fp, "\tnew_ldens %12.6f\n", neww);
		fprintf(blockpar->fp, "\told_ldens %12.6f\n", old);
		fprintf(blockpar->fp, "\tnew2old   %12.6f\n", new2old);
		fprintf(blockpar->fp, "\told2new   %12.6f\n", old2new);
	}

	Free(xx);
	Free(yy);
	Free(mode);
	if (free_block) {
		Free(blockpar);
	}

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_init_GMRF_approximation(GMRFLib_problem_tp ** problem, double *x, double *b, double *c, double *mean, double *d,
				    GMRFLib_logl_tp * loglFunc, void *loglFunc_arg,
				    GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg,
				    GMRFLib_constr_tp * constr, GMRFLib_optimize_param_tp * optpar, GMRFLib_blockupdate_param_tp * blockupdate_par)
{
	GMRFLib_ENTER_ROUTINE;
	GMRFLib_EWRAP1(GMRFLib_init_GMRF_approximation_store(problem, x, b, c, mean, d, loglFunc, loglFunc_arg,
							     graph, Qfunc, Qfunc_arg, constr, optpar, blockupdate_par, NULL));
	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_init_GMRF_approximation_store(GMRFLib_problem_tp ** problem, double *x, double *b, double *c, double *mean,
					  double *d, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg,
					  GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg,
					  GMRFLib_constr_tp * constr, GMRFLib_optimize_param_tp * optpar,
					  GMRFLib_blockupdate_param_tp * blockupdate_par, GMRFLib_store_tp * store)
{
	int i, j, free_x = 0, free_b = 0, free_c = 0, free_mean = 0, free_d = 0, free_blockpar = 0, n, id;
	double *bb = NULL, *cc = NULL, *mode = NULL;

#define FREE_ALL if (1) { if (free_x) Free(x); if (free_b) Free(b); if (free_c) Free(c); if (free_d) Free(d); \
	if (free_mean) Free(mean); if (free_blockpar) Free(blockupdate_par); Free(bb); Free(cc); Free(mode);}

	id = GMRFLib_thread_id;
	GMRFLib_ASSERT(problem, GMRFLib_EINVARG);
	GMRFLib_ASSERT(graph, GMRFLib_EINVARG);
	GMRFLib_ASSERT(Qfunc, GMRFLib_EINVARG);
	GMRFLib_ENTER_ROUTINE;

	n = graph->n;
	if (n == 0) {
		*problem = NULL;
		GMRFLib_LEAVE_ROUTINE;
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

	bb = Calloc(n, double);
	cc = Calloc(n, double);
	mode = Calloc(n, double);

	if (!blockupdate_par) {
		GMRFLib_default_blockupdate_param(&blockupdate_par);
		free_blockpar = 1;
	}

	Memcpy(mode, x, n * sizeof(double));
	if (blockupdate_par->modeoption == GMRFLib_MODEOPTION_MODE && d)
		GMRFLib_EWRAP1(GMRFLib_optimize_store(mode, b, c, mean, graph, Qfunc, Qfunc_arg, constr, d, loglFunc, loglFunc_arg, optpar, store));

	if (!(store && store->sub_graph)) {
		/*
		 * compute the terms from loglFunc 
		 */
		if (d) {
#pragma omp parallel for private(i)
			for (i = 0; i < n; i++) {
				GMRFLib_thread_id = id;
				double cmin = 0.0;
				if (d[i]) {
					GMRFLib_2order_approx(NULL, &bb[i], &cc[i], NULL, d[i], mode[i], i, mode, loglFunc, loglFunc_arg,
							      &(blockupdate_par->step_len), &(blockupdate_par->stencil), &cmin);
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
	} else {
		/*
		 * do the same as above, but only for those 'i' which is not fixed. this is faster if we have fixed values and
		 * the sub_graph is available. 
		 */

		int ns = store->sub_graph->n;

#pragma omp parallel for private(i, j)
		for (j = 0; j < ns; j++) {
			GMRFLib_thread_id = id;
			double cmin = 0.0;
			i = j;
			if (d[i]) {
				GMRFLib_2order_approx(NULL, &bb[i], &cc[i], NULL, d[i], mode[i], i, mode, loglFunc, loglFunc_arg,
						      &(blockupdate_par->step_len), &(blockupdate_par->stencil), &cmin);
			}
		}
		GMRFLib_thread_id = id;

		if (b) {
			if (mean) {
				for (j = 0; j < ns; j++) {
					i = j;
					bb[i] += b[i] - cc[i] * mean[i];
				}
			} else {
				for (j = 0; j < ns; j++) {
					i = j;
					bb[i] += b[i];
				}
			}
		} else {
			if (mean) {
				for (j = 0; j < ns; j++) {
					i = j;
					bb[i] -= cc[i] * mean[i];
				}
			}
		}

		if (c) {
			for (j = 0; j < ns; j++) {
				i = j;
				cc[i] += c[i];
			}
		}
	}

	GMRFLib_EWRAP1(GMRFLib_init_problem_store(problem, x, bb, cc, mean, graph, Qfunc, Qfunc_arg, constr, store));

	FREE_ALL;
	GMRFLib_LEAVE_ROUTINE;

	return GMRFLib_SUCCESS;

#undef FREE_ALL
}
