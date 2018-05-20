
/* optimize.c
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
  \file optimize.c
  \brief The optimising routines in GMRFLib.
*/

#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

/* Pre-hg-Id: $Id: optimize.c,v 1.65 2010/04/08 03:18:30 hrue Exp $ */

#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif

#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

/* 
   to flag various 'store' options for the use in the optimize-routines.

   i think these can be merged into the optimize routine, //fix later
*/
static int store_store_sub_graph = 0, store_use_sub_graph = 0, store_store_remap = 0, store_use_remap = 0, store_store_symb_fact =
    0, store_use_symb_fact = 0, store_smtp = 0;

#pragma omp threadprivate (store_store_sub_graph, store_use_sub_graph, store_store_remap, store_use_remap, store_store_symb_fact, store_use_symb_fact, store_smtp)

/*!
  \brief Creates a \c GMRFLib_optimize_param_tp -object holding the default settings.

  \param[out] optpar A pointer to a \c GMRFLib_optimize_param_tp pointer. 
  At output the \c GMRFLib_optimize_param_tp -object contains the default values.

  \par Description of elements in \c GMRFLib_optimize_param_tp -object:
  \n \em fp: A file for printing output from the optimizer. \n
  <b>Default value: \c NULL</b> \n\n
  \em opt_type: The type of optimizer to be used. Four methods are available.\n
  <b>Default value: #GMRFLib_OPTTYPE_SAFENR</b> \n\n
  \em nsearch_dir: Indicates the number of previously search gradient directions 
  on which the current search direction should be orthogonal (using the conjugate
  gradient (CG) method). \n
  <b>Default value: 1</b> \n\n
  \n \em nr_step_factor: Use reduced step-len in the Newton-Raphson iterations, where the step-length
  for iteration i, is MIN(1, (i+1)*nr_step_factor).\n
  <b>Default value: 1.0</b> \n\n
  \em restart_interval: If <em>restart_interval = r </em>,
  the CG search will be restarted every <em>r</em>'th iteration. \n
  <b>Default value: 10</b> \n\n
  \em max_iter:  The maximum number of iterations. \n
  <b>Default value: 200</b> \n\n
  \em fixed_iter: If > 0, then this fix the number of iterations, whatever 
  all other stopping-options. \n
  <b>Default value: 0</b> \n\n
  \em max_linesearch_iter: The maximum number of iterations within each search 
  direction (CG). \n
  <b>Default value: 200</b> \n\n
  \em step_len: Step length in the computation of a Taylor expansion or second 
  order approximation of the log-likelihood around a point 
  <em> \b x_0</em> (CG and NR). \n
  <b>Default value: 1.0e-4</b> \n\n
  \em abserr_func: The absolute error tolerance for the value of the function 
  to be optimized (CG and NR). \n
  <b>Default value: 0.5e-3</b> \n\n
  \em abserr_step: The absolute error tolerance, relative to the number of nodes, 
  for the size of one step of the optimizer (CG and NR). \n
  <b>Default value: 0.5e-3</b> \n
 */
int GMRFLib_default_optimize_param(GMRFLib_optimize_param_tp ** optpar)
{

	/*
	 * define default values for the optimizer. methods are
	 */

	*optpar = Calloc(1, GMRFLib_optimize_param_tp);

	(*optpar)->fp = NULL;
	// (*optpar)->fp = stdout;FIXME("set fp=stdout");
	(*optpar)->opt_type = GMRFLib_OPTTYPE_NR;
	(*optpar)->nr_step_factor = 1.0;
	(*optpar)->nsearch_dir = 1;
	(*optpar)->restart_interval = 10;
	(*optpar)->max_iter = 50;
	(*optpar)->fixed_iter = 0;
	(*optpar)->max_linesearch_iter = 50;
	(*optpar)->step_len = GMRFLib_eps(0.25);
	(*optpar)->stencil = 5;				       /* 3,5,7 */
	(*optpar)->abserr_func = 0.005;
	(*optpar)->abserr_step = 0.005;

	return GMRFLib_SUCCESS;
}

int GMRFLib_optimize_set_store_flags(GMRFLib_store_tp * store)
{
	/*
	 * set flags for 'store' 
	 */
	if (store) {
		/*
		 * create logicals 
		 */
		store_store_sub_graph = (store->sub_graph ? 0 : 1);
		store_use_sub_graph = !store_store_sub_graph;
		store_store_remap = (store->remap ? 0 : 1);
		store_use_remap = !store_store_remap;

		if (GMRFLib_valid_smtp((int) store->smtp) == GMRFLib_TRUE) {
			store_smtp = store->smtp;
		} else {
			store_smtp = GMRFLib_smtp;
		}
		if (store_smtp == GMRFLib_SMTP_TAUCS) {
			store_store_symb_fact = (store->TAUCS_symb_fact ? 0 : 1);
			store_use_symb_fact = !store_store_symb_fact;
		} else {
			store_store_symb_fact = 0;
			store_use_symb_fact = 0;
		}
	} else {
		store_store_sub_graph = 0;
		store_use_sub_graph = 0;
		store_store_remap = 0;
		store_use_remap = 0;
		store_store_symb_fact = 0;
		store_use_symb_fact = 0;
		store_smtp = GMRFLib_smtp;
	}
	return GMRFLib_SUCCESS;
}

/*!
  \brief The optimizer.

  The optimizer will use one out of four methods (given by \a optpar)
    - #GMRFLib_OPTTYPE_CG The conjugate gradient method
    - #GMRFLib_OPTTYPE_NR The Newton-Raphson method
    - #GMRFLib_OPTTYPE_SAFECG A two-step conjugate gradient method
    - #GMRFLib_OPTTYPE_SAFENR A two-step Newton-Raphson method

  The first two are the ordinary conjugate gradient (CG) and Newton-Raphson (NR) methods, while the
  last two are expected to be safer, doing the optimization in two steps: First finding the optimum
  using a diagonal precision matrix <em>\b Q</em>, an then adjusting for the correlations.\n\n
 */
int GMRFLib_optimize(double *mode, double *b, double *c, double *mean,
		     GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_args,
		     char *fixed_value, GMRFLib_constr_tp * constr, double *d, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg,
		     GMRFLib_optimize_param_tp * optpar)
{
	GMRFLib_ENTER_ROUTINE;
	GMRFLib_EWRAP1(GMRFLib_optimize_store
		       (mode, b, c, mean, graph, Qfunc, Qfunc_args, fixed_value, constr, d, loglFunc, loglFunc_arg, optpar, NULL));
	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}
int GMRFLib_optimize_store(double *mode, double *b, double *c, double *mean,
			   GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_args,
			   char *fixed_value, GMRFLib_constr_tp * constr,
			   double *d, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, GMRFLib_optimize_param_tp * optpar, GMRFLib_store_tp * store)
{
	/*
	 * locate the mode starting in 'mode' using a conjugate gradient algorithm where the search direction is Q-ortogonal to 
	 * the m previous ones.
	 * 
	 * exp(-0.5(x-mean)'Q(x-mean)+b'x + \sum_i d_i loglFunc(x_i,i,logl_arg)) && fixed-flags
	 * 
	 * and linear deterministic and stochastic constaints 
	 */
	int sub_n, i, j, node, nnode, free_optpar, id;
	double *cc = NULL, *initial_value = NULL;
	GMRFLib_store_tp *store_ptr;
	GMRFLib_optimize_problem_tp *opt_problem = NULL;

	GMRFLib_ENTER_ROUTINE;

	id = GMRFLib_thread_id;
	GMRFLib_optimize_set_store_flags(store);

	/*
	 * our first task, is to convert the opt_problem so we get rid of the fixed_values and then write the function in the
	 * canonical form
	 * 
	 * large parts here is adapted from `problem-setup.c'!!! 
	 */

	/*
	 * create new opt_problem 
	 */
	opt_problem = Calloc(1, GMRFLib_optimize_problem_tp);

	/*
	 * get the options 
	 */
	if (optpar) {
		opt_problem->optpar = optpar;
		free_optpar = 0;
	} else {
		GMRFLib_default_optimize_param(&(opt_problem->optpar));
		free_optpar = 1;
	}

	/*
	 * first, find the new graph. 
	 */
	if (store_use_sub_graph) {
		/*
		 * copy from store 
		 */
		GMRFLib_EWRAP1(GMRFLib_copy_graph(&(opt_problem->sub_graph), store->sub_graph));
	} else {
		/*
		 * compute it 
		 */
		GMRFLib_EWRAP1(GMRFLib_compute_subgraph(&(opt_problem->sub_graph), graph, fixed_value));

		/*
		 * store it in store if requested 
		 */
		if (store_store_sub_graph) {
			GMRFLib_EWRAP1(GMRFLib_copy_graph(&(store->sub_graph), opt_problem->sub_graph));
		}
	}

	sub_n = opt_problem->sub_graph->n;
	if (sub_n == 0) {				       /* fast return if there is nothing todo */
		GMRFLib_free_graph(opt_problem->sub_graph);
		Free(opt_problem);
		GMRFLib_LEAVE_ROUTINE;
		return GMRFLib_SUCCESS;
	}

	/*
	 * the mapping is there in the graph, make a new pointer in the opt_problem-definition 
	 */
	opt_problem->map = opt_problem->sub_graph->mothergraph_idx;

	/*
	 * setup space & misc 
	 */
	{
		{
			int ii;
			opt_problem->mode = Calloc(sub_n, double);

			for (ii = 0; ii < sub_n; ii++) {
				opt_problem->mode[ii] = mode[opt_problem->map[ii]];
			}
		}
		{
			opt_problem->b = Calloc(sub_n, double);

			if (b) {
				int ii;

				for (ii = 0; ii < sub_n; ii++) {
					opt_problem->b[ii] = b[opt_problem->map[ii]];
				}
			}
		}
		{
			opt_problem->d = Calloc(sub_n, double);

			if (d) {
				int ii;

				for (ii = 0; ii < sub_n; ii++) {
					opt_problem->d[ii] = d[opt_problem->map[ii]];
				}
			}
		}
		{
			cc = Calloc(sub_n, double);

			if (c) {
				int ii;

				for (ii = 0; ii < sub_n; ii++) {
					cc[ii] = c[opt_problem->map[ii]];
				}
			}
		}
		{
			opt_problem->x_vec = Calloc(graph->n, double);
			memcpy(opt_problem->x_vec, mode, graph->n * sizeof(double));
		}
	}

	/*
	 * FIXME: i might want to make a wrapper of this one, then REMEMBER TO CHANGE ->map[i] in all calls to
	 * GMRFLib_2order_approx to i, or what it should be.!!! 
	 */
	opt_problem->loglFunc = loglFunc;		       /* i might want to make a wrapper of this one! */
	opt_problem->loglFunc_arg = loglFunc_arg;

	/*
	 * make the arguments to the wrapper function 
	 */
	opt_problem->sub_Qfunc = GMRFLib_Qfunc_wrapper;
	opt_problem->sub_Qfunc_arg = Calloc(1, GMRFLib_Qfunc_arg_tp);
	opt_problem->sub_Qfunc_arg->map = opt_problem->map;    /* yes, this ptr is needed */
	opt_problem->sub_Qfunc_arg->diagonal_adds = cc;

	opt_problem->sub_Qfunc_arg->user_Qfunc = Qfunc;
	opt_problem->sub_Qfunc_arg->user_Qfunc_args = Qfunc_args;

	/*
	 * now compute the new 'effective' b, and then the mean. recall to add the 'c' term manually, since we're using the
	 * original Qfunc. 
	 */
	if (!fixed_value) {				       /* then sub_graph = graph and map=I */
		if (mean) {
			double *tmp = NULL;
			tmp = Calloc(graph->n, double);

			GMRFLib_Qx(tmp, mean, graph, Qfunc, Qfunc_args);
			for (i = 0; i < graph->n; i++) {
				opt_problem->b[i] += tmp[i] + cc[i] * mean[i];
			}
			Free(tmp);
		}
	} else {
		/*
		 * x=(x1,x2), then x1|x2 has b = Q11 \mu1 - Q12(x2-\mu2) 
		 */
#pragma omp parallel for private(i, j, node, nnode)
		for (i = 0; i < sub_n; i++) {		       /* loop over all sub_nodes */
			GMRFLib_thread_id = id;
			node = opt_problem->map[i];
			if (mean) {
				opt_problem->b[i] += ((*Qfunc) (node, node, Qfunc_args) + cc[i]) * mean[node];
			}

			for (j = 0; j < graph->nnbs[node]; j++) {	/* then over all neighbors */
				double qvalue;

				nnode = graph->nbs[node][j];
				qvalue = (*Qfunc) (node, nnode, Qfunc_args);

				if (fixed_value[nnode]) {
					/*
					 * nnode is fixed 
					 */
					if (mean) {
						opt_problem->b[i] -= qvalue * (mode[nnode] - mean[nnode]);
					} else {
						opt_problem->b[i] -= qvalue * mode[nnode];
					}
				} else {
					/*
					 * nnone is not fixed 
					 */
					if (mean) {
						opt_problem->b[i] += qvalue * mean[nnode];
					}
				}
			}
		}
		GMRFLib_thread_id = id;
	}

	/*
	 * save initial value and return that unchanged if fail to converge 
	 */
	initial_value = Calloc(sub_n, double);
	memcpy(initial_value, opt_problem->mode, sub_n * sizeof(double));

	GMRFLib_ASSERT(opt_problem->optpar->opt_type == GMRFLib_OPTTYPE_CG ||
		       opt_problem->optpar->opt_type == GMRFLib_OPTTYPE_NR ||
		       opt_problem->optpar->opt_type == GMRFLib_OPTTYPE_SAFECG
		       || opt_problem->optpar->opt_type == GMRFLib_OPTTYPE_SAFENR, GMRFLib_EPARAMETER);

	/*
	 * do stuff relating to the constraint 
	 */
	if (constr) {
		/*
		 * constraints only allowed when using the newton-raphson optimizer 
		 */
		double *b_add = NULL;

		GMRFLib_ASSERT(constr && (opt_problem->optpar->opt_type == GMRFLib_OPTTYPE_NR ||
					  opt_problem->optpar->opt_type == GMRFLib_OPTTYPE_SAFENR), GMRFLib_EPARAMETER);

		b_add = Calloc(sub_n, double);

		GMRFLib_EWRAP1(GMRFLib_recomp_constr(&(opt_problem->sub_constr), constr, mode, b_add, fixed_value, graph, opt_problem->sub_graph));

		for (i = 0; i < sub_n; i++) {
			opt_problem->b[i] += b_add[i];
		}
		Free(b_add);
	}

	if (opt_problem->optpar->opt_type == GMRFLib_OPTTYPE_SAFENR || opt_problem->optpar->opt_type == GMRFLib_OPTTYPE_SAFECG) {
		/*
		 * to be more safe regarding initial valus we proceed as follows: first we change the graph to a diagonal, so
		 * to optimize without any interactions using NR/CG. this optimum is then the startvalue for the NR/CG
		 * optimizer. 
		 */
		int nzero_idx = -1;
		GMRFLib_graph_tp *diag_graph = NULL, *tmp_graph;

		if (store && !(store->diag_store)) {
			/*
			 * use diag_store here 
			 */
			store->diag_store = Calloc(1, GMRFLib_store_tp);
		}
		store_ptr = (store ? store->diag_store : store);

		GMRFLib_EWRAP1(GMRFLib_compute_subgraph(&diag_graph, opt_problem->sub_graph, NULL));
		for (i = 0; i < diag_graph->n; i++)
			if (diag_graph->nnbs[i]) {
				diag_graph->nnbs[i] = 0;
				if (nzero_idx < 0) {
					nzero_idx = i;	       /* the first non-zero pointer */
				}
			}

		tmp_graph = opt_problem->sub_graph;
		opt_problem->sub_graph = diag_graph;

		if (opt_problem->optpar->opt_type == GMRFLib_OPTTYPE_SAFECG) {
			GMRFLib_EWRAP1(GMRFLib_optimize2(opt_problem, store_ptr));
		} else {
			GMRFLib_EWRAP1(GMRFLib_optimize3(opt_problem, store_ptr));
		}
		opt_problem->sub_graph = tmp_graph;

		if (nzero_idx >= 0) {
			diag_graph->nnbs[nzero_idx] = 1;       /* to make _free_graph work correct */
		}
		GMRFLib_free_graph(diag_graph);
	}

	if (store && fixed_value && !(store->sub_store)) {
		/*
		 * we cannot use the ordinary 'store' as optimize2/3 work on sub_graph. however, we can use the sub_store for
		 * this purpose. prepare it if not already done. 
		 */
		store->sub_store = Calloc(1, GMRFLib_store_tp);
	}

	/*
	 * which store to use, if any 
	 */
	if (store) {
		store_ptr = (fixed_value ? store->sub_store : store);
	} else {
		store_ptr = store;
	}

	switch (opt_problem->optpar->opt_type) {
	case GMRFLib_OPTTYPE_CG:
	case GMRFLib_OPTTYPE_SAFECG:
		GMRFLib_EWRAP1(GMRFLib_optimize2(opt_problem, store_ptr));
		if (0) {				       /* FIXME: if (fail) then do as follows */
			memcpy(opt_problem->mode, initial_value, sub_n * sizeof(double));
			Free(initial_value);
			GMRFLib_ERROR(GMRFLib_EOPTCG);
		}
		break;

	case GMRFLib_OPTTYPE_NR:
	case GMRFLib_OPTTYPE_SAFENR:
		GMRFLib_EWRAP1(GMRFLib_optimize3(opt_problem, store_ptr));
		if (0) {				       /* FIXME: if (fail) then do as follows */
			memcpy(opt_problem->mode, initial_value, sub_n * sizeof(double));
			Free(initial_value);
			GMRFLib_ERROR(GMRFLib_EOPTNR);
		}
		break;
	}
	/*
	 * copy back 
	 */
	for (i = 0; i < sub_n; i++) {
		mode[opt_problem->map[i]] = opt_problem->mode[i];
	}

	/*
	 * free what's malloced and return 
	 */
	Free(opt_problem->mode);
	Free(opt_problem->b);
	Free(opt_problem->d);
	Free(opt_problem->sub_Qfunc_arg);
	Free(opt_problem->x_vec);
	if (constr) {
		GMRFLib_free_constr(opt_problem->sub_constr);
	}
	if (free_optpar) {
		Free(opt_problem->optpar);
	}
	GMRFLib_free_graph(opt_problem->sub_graph);
	Free(opt_problem);
	Free(initial_value);
	Free(cc);

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_optimize2(GMRFLib_optimize_problem_tp * opt_problem, GMRFLib_store_tp * store)
{

	/*
	 * optimize, using GMRFLib_Qfunc_wrapper as Qfunc. NOTE: the loglFunc is NOT wrapped, so the indices are in the real
	 * world.
	 * 
	 * this is conjugate gradient optimizer 
	 */

	int sub_n, i, nsdir, iter, sdir_indx, n_correct, fail, id;
	double **sdir = NULL, *grad = NULL, *omode = NULL, err, *c_orig = NULL;

	id = GMRFLib_thread_id;
	nsdir = opt_problem->optpar->nsearch_dir + 1;	       /* also contain the current dir */
	sub_n = opt_problem->sub_graph->n;
	grad = Calloc(sub_n, double);
	omode = Calloc(sub_n, double);
	c_orig = Calloc(sub_n, double);

	memcpy(c_orig, opt_problem->sub_Qfunc_arg->diagonal_adds, sub_n * sizeof(double));

	sdir = Calloc(nsdir, double *);

	for (i = 0; i < nsdir; i++) {
		sdir[i] = Calloc(sub_n, double);
	}

	if (opt_problem->optpar->fp) {
		fprintf(opt_problem->optpar->fp, "\n%6s%22s%6s%22s\n", "Iter", "Value", "SubIt", "StepLength");
		fprintf(opt_problem->optpar->fp, " --------------------------------------------------------\n");
	}

	for (iter = 0, sdir_indx = 0, fail = 1;
	     iter < IMAX(opt_problem->optpar->max_iter, opt_problem->optpar->fixed_iter); sdir_indx = ((sdir_indx + 1) % nsdir), iter++) {
		if (opt_problem->optpar->fp) {
			fprintf(opt_problem->optpar->fp, "%6d", iter);
		}

		memcpy(opt_problem->sub_Qfunc_arg->diagonal_adds, c_orig, sub_n * sizeof(double));

		/*
		 * first compute the gradient, -Qx + b 
		 */
		GMRFLib_Qx(grad, opt_problem->mode, opt_problem->sub_graph, opt_problem->sub_Qfunc, (void *) (opt_problem->sub_Qfunc_arg));
		for (i = 0; i < sub_n; i++) {
			grad[i] = opt_problem->b[i] - grad[i];
		}

		/*
		 * then compute contribution from the loglFunc 
		 */
		for (i = 0; i < sub_n; i++) {
			opt_problem->x_vec[opt_problem->map[i]] = opt_problem->mode[i];
		}
#pragma omp parallel for private(i)
		for (i = 0; i < sub_n; i++) {
			GMRFLib_thread_id = id;
			if (opt_problem->d[i]) {
				double bcoof, ccoof;

				GMRFLib_2order_taylor(NULL, &bcoof, &ccoof, opt_problem->d[i],
						      opt_problem->mode[i], opt_problem->map[i], opt_problem->x_vec,
						      opt_problem->loglFunc, opt_problem->loglFunc_arg, &(opt_problem->optpar->step_len),
						      &(opt_problem->optpar->stencil));
				grad[i] += bcoof;
				opt_problem->sub_Qfunc_arg->diagonal_adds[i] += DMAX(0.0, -ccoof);
			}
		}
		GMRFLib_thread_id = id;

		/*
		 * subtract terms Q-ortogonal to the previous search-directions 
		 */
		if ((iter + 1) % opt_problem->optpar->restart_interval != 0) {
			n_correct = IMIN(opt_problem->optpar->nsearch_dir, iter);
			for (i = 0; i < n_correct; i++) {
				GMRFLib_Qadjust(grad, sdir[MOD(sdir_indx - i - 1, nsdir)],
						opt_problem->sub_graph, opt_problem->sub_Qfunc, (void *) (opt_problem->sub_Qfunc_arg));
			}
		} else {
			for (i = 0; i < opt_problem->optpar->nsearch_dir; i++) {
				memset(sdir[i], 0, sub_n * sizeof(double));
			}
		}
		memcpy(opt_problem->sub_Qfunc_arg->diagonal_adds, c_orig, sub_n * sizeof(double));	/* go back to original */
		memcpy(sdir[sdir_indx], grad, sub_n * sizeof(double));

		/*
		 * we got a search-direction, do a line-search in that direction and update 'mode' 
		 */
		memcpy(omode, opt_problem->mode, sub_n * sizeof(double));
		GMRFLib_EWRAP0(GMRFLib_linesearch(opt_problem, sdir[sdir_indx]));

		for (i = 0, err = 0.0; i < sub_n; i++) {
			err += SQR(omode[i] - opt_problem->mode[i]);
		}
		err = sqrt(DMAX(0.0, err) / sub_n);
		if (opt_problem->optpar->fp) {
			fprintf(opt_problem->optpar->fp, "%22.10e\n", err);
		}

		if (opt_problem->optpar->fixed_iter > 0) {
			if (opt_problem->optpar->fixed_iter == iter + 1) {
				fail = 0;
				break;
			}
		} else {
			if (err < opt_problem->optpar->abserr_func) {
				fail = 0;
				break;
			}
		}
	}

	memcpy(opt_problem->sub_Qfunc_arg->diagonal_adds, c_orig, sub_n * sizeof(double));

	for (i = 0; i < nsdir; i++) {
		Free(sdir[i]);
	}
	Free(sdir);
	Free(omode);
	Free(c_orig);
	Free(grad);

	return (fail ? GMRFLib_EOPTCG : GMRFLib_SUCCESS);
}
int GMRFLib_Qadjust(double *dir, double *odir, GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg)
{
	/*
	 * compute: dir := dir - (dir,odir)/(odir,odir), where (u,v) = u^TQv 
	 */

	int i;
	double a, b, c, *v = NULL;

	v = Calloc(graph->n, double);

	GMRFLib_Qx(v, odir, graph, Qfunc, Qfunc_arg);
	for (i = 0, a = b = 0.0; i < graph->n; i++) {
		a += v[i] * dir[i];
		b += v[i] * odir[i];
	}

	if (b > 0.0 && !ISZERO(a)) {
		c = a / b;				       /* FIX? */
		for (i = 0; i < graph->n; i++) {
			dir[i] -= c * odir[i];
		}
	}

	Free(v);
	return GMRFLib_SUCCESS;
}
double GMRFLib_linesearch_func(double length, double *dir, GMRFLib_optimize_problem_tp * opt_problem)
{
	int i, sub_n, id;
	double *v = NULL, *u = NULL, fval = 0.0;

	id = GMRFLib_thread_id;
	sub_n = opt_problem->sub_graph->n;
	v = Calloc(sub_n, double);
	u = Calloc(sub_n, double);

	if (length != 0.0) {
		for (i = 0; i < sub_n; i++) {
			u[i] = opt_problem->mode[i] + length * dir[i];
		}
	} else {
		memcpy(u, opt_problem->mode, sub_n * sizeof(double));
	}
	GMRFLib_Qx(v, u, opt_problem->sub_graph, opt_problem->sub_Qfunc, (void *) (opt_problem->sub_Qfunc_arg));

	for (i = 0; i < sub_n; i++) {
		opt_problem->x_vec[opt_problem->map[i]] = u[i];
	}
#if defined(_OPENMP)
	{
		double sum = 0.0;

#pragma omp parallel for private(i) reduction(+: sum)
		for (i = 0; i < sub_n; i++) {
			double logll;

			GMRFLib_thread_id = id;
			sum += (-0.5 * v[i] + opt_problem->b[i]) * u[i];
			if (opt_problem->d[i]) {
				(*(opt_problem->loglFunc)) (&logll, &u[i], 1, opt_problem->map[i], opt_problem->x_vec, NULL,
							    opt_problem->loglFunc_arg);
				sum += opt_problem->d[i] * logll;
			}
		}
		GMRFLib_thread_id = id;
		fval += sum;
	}
#else
	for (i = 0, fval = 0.0; i < sub_n; i++) {
		double logll;

		GMRFLib_thread_id = id;
		fval += (-0.5 * v[i] + opt_problem->b[i]) * u[i];
		if (opt_problem->d[i]) {
			(*(opt_problem->loglFunc)) (&logll, &u[i], 1, opt_problem->map[i], opt_problem->x_vec, NULL, opt_problem->loglFunc_arg);
			fval += opt_problem->d[i] * logll;
		}
	}
	GMRFLib_thread_id = id;
#endif

	if (STOCHASTIC_CONSTR(opt_problem->sub_constr)) {
		double sqr_val;

		GMRFLib_eval_constr(NULL, &sqr_val, u, opt_problem->sub_constr, opt_problem->sub_graph);
		fval += -0.5 * sqr_val;
	}

	Free(u);
	Free(v);

	return fval;
}
int GMRFLib_linesearch(GMRFLib_optimize_problem_tp * opt_problem, double *dir)
{
	int i, sub_n, iter = 0, fail, id;
	double len, deriv, dderiv, *u = NULL, err, *loglikgrad = NULL, *loglikggrad = NULL;
	double *len_history, len_eps = sqrt(FLT_EPSILON), periodic_flag = 0;

	id = GMRFLib_thread_id;
	sub_n = opt_problem->sub_graph->n;
	len_history = Calloc(opt_problem->optpar->max_linesearch_iter, double);
	loglikggrad = Calloc(sub_n, double);
	loglikgrad = Calloc(sub_n, double);
	u = Calloc(sub_n, double);

	for (iter = 0, fail = 1; iter < opt_problem->optpar->max_linesearch_iter; iter++) {
		for (i = 0; i < sub_n; i++) {
			opt_problem->x_vec[opt_problem->map[i]] = opt_problem->mode[i];
		}
#pragma omp parallel for private(i)
		for (i = 0; i < sub_n; i++) {
			GMRFLib_thread_id = id;
			if (opt_problem->d[i]) {
				GMRFLib_2order_taylor(NULL, &loglikgrad[i], &loglikggrad[i], 1.0,
						      opt_problem->mode[i], opt_problem->map[i], opt_problem->x_vec,
						      opt_problem->loglFunc, opt_problem->loglFunc_arg, &(opt_problem->optpar->step_len),
						      &(opt_problem->optpar->stencil));
			} else {
				loglikgrad[i] = loglikggrad[i] = 0.0;
			}
		}
		GMRFLib_thread_id = id;

		GMRFLib_Qx(u, dir, opt_problem->sub_graph, opt_problem->sub_Qfunc, (void *) opt_problem->sub_Qfunc_arg);

		deriv = dderiv = 0.0;
		for (i = 0; i < sub_n; i++) {
			deriv += -opt_problem->mode[i] * u[i] + opt_problem->b[i] * dir[i] + opt_problem->d[i] * dir[i] * loglikgrad[i];
			dderiv += -dir[i] * u[i] + opt_problem->d[i] * SQR(dir[i]) * loglikggrad[i];
		}

		if (!ISZERO(deriv) && !ISZERO(dderiv)) {
			double len_max = 1.0;

			len = len_history[iter] = DMIN(len_max, DMAX(-len_max, -deriv / dderiv));

			/*
			 * make a robustnes check for the behaviour len = a, -a, a, -a, a, -a, etc... 
			 */
			if (iter >= 4) {
				if (ABS(len_history[iter] - len_history[iter - 2]) < len_eps &&
				    ABS(len_history[iter - 1] - len_history[iter - 3]) < len_eps
				    && ABS(len_history[iter] + len_history[iter - 1]) < len_eps) {
					len = len_history[iter] / 2.;
					len_history[iter] = len;
					periodic_flag = 1;     /* !!!stop after this step!!! */
					if (0)
						printf("\nperiodic behaviour detected, set len = %f\n", len);
				}
			}

			for (i = 0, err = 0.0; i < sub_n; i++) {
				opt_problem->mode[i] += len * dir[i];
				err += SQR(dir[i]);
			}
		} else {
			len = err = 0.0;
		}
		err = sqrt(DMAX(err * SQR(len) / sub_n, 0.0));
		if (err < opt_problem->optpar->abserr_step) {
			fail = 0;
			break;
		}
		if (periodic_flag) {
			fail = 0;
			break;
		}
		if (0) {
			printf("\n%d %.12f %.12f\n", iter, len, GMRFLib_linesearch_func(0.0, dir, opt_problem));
		}
	}
	if (fail)
		GMRFLib_ERROR(GMRFLib_EOPTCGLINE);

	if (opt_problem->optpar->fp) {
		fprintf(opt_problem->optpar->fp, "%22.10e%6d", GMRFLib_linesearch_func(0.0, dir, opt_problem), iter + 1);
	}

	Free(loglikgrad);
	Free(loglikggrad);
	Free(u);
	Free(len_history);
	return GMRFLib_SUCCESS;
}

int GMRFLib_optimize3(GMRFLib_optimize_problem_tp * opt_problem, GMRFLib_store_tp * store)
{
	/*
	 * optimize, using GMRFLib_Qfunc_wrapper as Qfunc. NOTE: the loglFunc is NOT wrapped, so the indices are in the real
	 * world.
	 * 
	 * this is the newton-raphson optimizer 
	 */

	int sub_n, i, iter, fail = 1, id;
	int *idxs = NULL, nidx = 0;
	double *bb = NULL, err, *c_orig = NULL, step_factor, f;
	GMRFLib_problem_tp *problem = NULL;

	id = GMRFLib_thread_id;
	sub_n = opt_problem->sub_graph->n;
	bb = Calloc(sub_n, double);
	c_orig = Calloc(sub_n, double);

	if (opt_problem->optpar->nr_step_factor > 0) {
		step_factor = opt_problem->optpar->nr_step_factor;
	} else {
		/*
		 * default choice
		 */
		step_factor = 1.0;
	}

	memcpy(c_orig, opt_problem->sub_Qfunc_arg->diagonal_adds, sub_n * sizeof(double));

	if (opt_problem->optpar->fp) {
		fprintf(opt_problem->optpar->fp, "\n%6s%22s%22s\n", "Iter", "Value", "StepLength");
		fprintf(opt_problem->optpar->fp, " --------------------------------------------------\n");
	}

	for (iter = 0, fail = 1; iter < IMAX(opt_problem->optpar->max_iter, opt_problem->optpar->fixed_iter); iter++) {
		if (opt_problem->optpar->fp) {
			fprintf(opt_problem->optpar->fp, "%6d", iter);
		}
		memcpy(opt_problem->sub_Qfunc_arg->diagonal_adds, c_orig, sub_n * sizeof(double));

		/*
		 * compute contribution from the loglFunc 
		 */
		for (i = 0; i < sub_n; i++) {
			opt_problem->x_vec[opt_problem->map[i]] = opt_problem->mode[i];
		}

		if (!idxs) {
			/*
			 * be somewhat more clever in this loop to divide the work better 
			 */
			idxs = Calloc(sub_n, int);

			for (i = 0; i < sub_n; i++) {
				if (opt_problem->d[i]) {
					idxs[nidx++] = i;
				} else {
					bb[i] = opt_problem->b[i];
				}
			}
		}
#pragma omp parallel for private(i)
		for (i = 0; i < nidx; i++) {
			int idx;
			double bcoof, ccoof;

			GMRFLib_thread_id = id;
			idx = idxs[i];
			GMRFLib_2order_approx(NULL, &bcoof, &ccoof, opt_problem->d[idx],
					      opt_problem->mode[idx], opt_problem->map[idx], opt_problem->x_vec,
					      opt_problem->loglFunc, opt_problem->loglFunc_arg, &(opt_problem->optpar->step_len),
					      &(opt_problem->optpar->stencil));
			bb[idx] = opt_problem->b[idx] + bcoof;
			opt_problem->sub_Qfunc_arg->diagonal_adds[idx] += DMAX(0.0, ccoof);
		}
		GMRFLib_thread_id = id;


		GMRFLib_EWRAP0(GMRFLib_init_problem_store(&problem, opt_problem->mode, bb, NULL, NULL, opt_problem->sub_graph,
							  opt_problem->sub_Qfunc, (void *) opt_problem->sub_Qfunc_arg,
							  NULL, opt_problem->sub_constr,
							  (unsigned int) (iter ==
									  0 ? GMRFLib_NEW_PROBLEM : GMRFLib_KEEP_graph | GMRFLib_KEEP_constr),
							  store));
		f = DMIN(1.0, step_factor * (iter + 1));
		if (problem->sub_constr) {
			for (i = 0, err = 0.0; i < sub_n; i++) {
				err += SQR(problem->sub_mean_constr[i] - opt_problem->mode[i]);
				opt_problem->mode[i] = opt_problem->mode[i] + f * (problem->sub_mean_constr[i] - opt_problem->mode[i]);
			}
		} else {
			for (i = 0, err = 0.0; i < sub_n; i++) {
				err += SQR(problem->sub_mean[i] - opt_problem->mode[i]);
				opt_problem->mode[i] = opt_problem->mode[i] + f * (problem->sub_mean[i] - opt_problem->mode[i]);
			}
		}

		err = sqrt(DMAX(0.0, err) / sub_n);
		if (opt_problem->optpar->fp) {
			memcpy(opt_problem->sub_Qfunc_arg->diagonal_adds, c_orig, sub_n * sizeof(double));
			fprintf(opt_problem->optpar->fp, "%22.10e%22.10e\n", GMRFLib_linesearch_func(0.0, NULL, opt_problem), err);
		}

		if (opt_problem->optpar->fixed_iter > 0) {
			if (opt_problem->optpar->fixed_iter == iter + 1) {
				fail = 0;
				break;
			}
		} else {
			if (err < opt_problem->optpar->abserr_step) {
				fail = 0;
				break;
			}
		}
	}

	memcpy(opt_problem->sub_Qfunc_arg->diagonal_adds, c_orig, sub_n * sizeof(double));

	GMRFLib_free_problem(problem);
	Free(bb);
	Free(c_orig);
	Free(idxs);

	return (fail ? GMRFLib_EOPTNR : GMRFLib_SUCCESS);
}
