
/* problem-setup.c
 * 
 * Copyright (C) 2001-2018 Havard Rue
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
  \file problem-setup.c
  \brief Sampling etc from a GMRF
  
  Sampling, whether unconditionally or conditionally, from a GMRF on a general graph, is performed
  using one function, \ref GMRFLib_sample().  Similarly, \ref GMRFLib_evaluate() computes the
  log-density for a sample from a GMRF and \ref GMRFLib_Qinv() computes elements in the inverse of
  the precision matrix.  The functions operate on a data structure, \ref GMRFLib_problem_tp, holding
  all external and internal information needed by the sampling and evaluation algorithms. This data
  structure is initialised by \ref GMRFLib_init_problem().

*/
#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

/* Pre-hg-Id: $Id: problem-setup.c,v 1.166 2010/02/18 13:45:00 hrue Exp $ */

#if !defined(__FreeBSD__)
#include <malloc.h>
#endif

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

double GMRFLib_Qfunc_wrapper(int sub_node, int sub_nnode, void *arguments)
{
	int node, nnode;
	double val;
	GMRFLib_Qfunc_arg_tp *args = NULL;

	args = (GMRFLib_Qfunc_arg_tp *) arguments;

	node = args->map[sub_node];
	nnode = args->map[sub_nnode];

	if (node == nnode) {
		val = (*(args->user_Qfunc)) (node, nnode, args->user_Qfunc_args) + args->diagonal_adds[sub_node];
	} else {
		val = (*(args->user_Qfunc)) (node, nnode, args->user_Qfunc_args);
	}

	return val;
}

/*! \brief Initializes and specifies a \c GMRFLib_problem_tp -object holding all 
  information needed for sampling from a GMRF and evaluating the log-density 
  of such a sample.

  \param[in,out] problem At output, <em>(*problem)</em> is a pointer to a 
  \c GMRFLib_problem_tp -object, initialised and defined according to the problem 
  specification. The required value of <em>(*problem)</em> at input depends on the 
  argument \a keep. If \a keep = 0, <em>(*problem)</em> is allocated within the 
  routine, and thus expected to be \c NULL at input. If \a keep > 0, 
  <em>(*problem)</em> is assumed to be an already initialised \c GMRFLib_problem_tp 
  -pointer at input.
  See also \a keep and Remarks.
  \param[in] x  A length \em n array, where \em n is the number of nodes in 
  the graph, of initial values of the GMRF. If \em x = \c NULL then all elements
  are taken to be zero.  If \a fixed_value \f$ \neq \f$ \c NULL, the elements of 
  <em>\b x</em> corresponding to \a fixed_value=1 are the fixed values in a
  conditional simulation. The remaining elements of <em>\b x</em> are not used, 
  and can take arbitrary values. If \a fixed_value = \c NULL, all values can be
  arbitrary.
  \param[in] b If <tt>!NULL</tt>, a length \em n array holding the elements of 
  the vector <em>\b b</em> in the general expression of the density if a GMRF, 
  as given in <b>(GMRF-2)</b> in \ref description.
  \param[in] c If <tt>!NULL</tt>, a length \em n array of elements to add to 
  the diagonal of the precision matrix defined by the function \a Qfunc. The 
  argument \a c should hold the elements of the vector <em>\b c</em> in 
  <b>(GMRF-2)</b>.
  \param[in] mean If <tt>!NULL</tt>, a length \em n array holding the elements 
  of the vector \f$ \mbox{\boldmath $\mu$} \f$ in <b>(GMRF-2)</b>.  
  If <em>\b b</em> = 0, \a mean is the mean of the GMRF.
  \param[in] graph The graph on which the GMRF is defined.
  \param[in] Qfunc A pointer to a user-defined function defining the
  precision matrix <em>\b Q</em> of a GMRF <em>\b x</em>.
  \param[in] Qfunc_args The arguments to the function \a Qfunc defining the 
  precision matrix.
  \param[in] fixed_value If <tt>!NULL</tt>, the sampling is done conditionally 
  on fixed values. The elements of this array, of length \em n, should take 
  the value 0 or 1. The conditional mean and covariance matrix of 
  \f$ \{x_i:\mbox{\small\tt fixed\_value}[i]=0\} \f$ given 
  \f$ \{x_i:\mbox{\small\tt fixed\_value}[i]=1\} \f$ are computed, and by 
  sampling using \c GMRFLib_sample(), elements of the GMRF <em>\b x</em>
  corresponding to <tt>fixed_value = 1</tt>, are kept fixed at the 
  corresponding values specified in the argument \a x.
  \param[in] constr If <tt>!NULL</tt>, a pointer to a \c GMRFLib_constr_tp
  -variable holding information about a linear (deterministic or stochastic) 
  constraint. If \c NULL, the sampling is unconstrained.
  \param keep If 0, <em>(*problem)</em> is allocated within the routine, and thus 
  expected to be \c NULL} at input. If >0, the <em>(*problem)</em>-pointer is 
  supposed to be already initialised, pointing at a problem specification 
  object, of which parts of the problem specification is to be kept unchanged. 
  Which parts to keep fixed, are defined by the actual value of \a keep. 
  See Remarks below.

  \remarks This function is to be called prior to sampling and evaluation of 
  the log-density. If \a keep = 0 (or \c GMRFLib_NEW_PROBLEM()), the 
  \c GMRFLib_problem_tp -pointer <em>(*problem)</em> is dynamically allocated 
  within the function. If \a keep > 0, <em>(*problem)</em> is required to be
  non- \c NULL, and is not allocated within the routine. One or more of the 
  members <em>(*problem)</em> is to be kept fixed and not be re-allocated. 
  What parts of the problem specification to keep fixed during sampling, is 
  specified by the value of \a keep. It should take one out of four predefined 
  macro values, or combinations of these. The choices are 
  \n - \c GMRFLib_KEEP_mean : Keep the computed conditional mean value 
  \n - \c GMRFLib_KEEP_graph : Keep the specified graph  
  \n - \c GMRFLib_KEEP_chol : Keep the Cholesky factorization of the 
  <em>\b Q</em>-matrix.  
  \n - \c GMRFLib_KEEP_constr : Keep the constraints \n
  Two or more of the \a keep-values are combined by letting
  keep = \f$ \mbox{\small\tt val1}|\mbox{\small\tt val2}|\ldots \f$
  The allocated memory of a \c GMRFLib_problem_tp -pointer is freed using the 
  function \c GMRFLib_free_problem().  Letting \a fixed_value \f$ \neq \f$ 
  \c NULL, the problem is defined to be a conditional sampling problem, 
  conditioning on elements for which \a fixed_value = 1. With \a constraint
  \f$ \neq \f$ \c NULL, sampling will be done conditionally on a linear 
  constraint. If both \a fixed_value and \a constraint are \c NULL, 
  unconditional samples will be generated by \c GMRFLib_sample().  
  The arguments \a problem, \a x, \a graph and \a Qfunc are all required 
  to be <tt>!NULL</tt> at input.

  \par Example:
  See \ref ex_problem-setup,\n \ref ex_wa and\n \ref ex_blockupdate.

  \sa GMRFLib_Qfunc_tp, GMRFLib_sample, GMRFLib_evaluate, GMRFLib_free_problem.
*/
int GMRFLib_init_problem(GMRFLib_problem_tp ** problem,
			 double *x,
			 double *b,
			 double *c,
			 double *mean,
			 GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_args, char *fixed_value,
			 GMRFLib_constr_tp * constr, unsigned int keep)
{
	GMRFLib_ENTER_ROUTINE;
	GMRFLib_EWRAP1(GMRFLib_init_problem_store(problem, x, b, c, mean, graph, Qfunc, Qfunc_args, fixed_value, constr, keep, NULL));
	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}
int GMRFLib_init_problem_store(GMRFLib_problem_tp ** problem,
			       double *x,
			       double *b,
			       double *c,
			       double *mean,
			       GMRFLib_graph_tp * graph,
			       GMRFLib_Qfunc_tp * Qfunc,
			       void *Qfunc_args, char *fixed_value, GMRFLib_constr_tp * constr, unsigned int keep, GMRFLib_store_tp * store)
{
	double *bb = NULL;
	int i, j, sub_n, node, nnode, free_x = 0, id;
	GMRFLib_smtp_tp smtp;

	int store_store_sub_graph = 0, store_use_sub_graph = 0;
	int store_store_remap = 0, store_use_remap = 0;
	int store_store_symb_fact = 0, store_use_symb_fact = 0;

	GMRFLib_Qfunc_tp *sub_Qfunc;
	GMRFLib_Qfunc_arg_tp *sub_Qfunc_arg;

	GMRFLib_ENTER_ROUTINE;

	id = GMRFLib_thread_id;
	GMRFLib_ASSERT(problem, GMRFLib_EINVARG);
	GMRFLib_ASSERT(graph, GMRFLib_EINVARG);
	GMRFLib_ASSERT(Qfunc, GMRFLib_EINVARG);

	if (keep == GMRFLib_NEW_PROBLEM) {
		*problem = Calloc(1, GMRFLib_problem_tp);
	}

	/*
	 * whatever to be stored, the Qinv is no longer valid. 
	 */
	if (keep != GMRFLib_NEW_PROBLEM) {
		GMRFLib_free_Qinv(*problem);
	}

	/*
	 * define the events, 'create a store', 'use a store'. these are possible if and only if 'keep==NEW_PROBLEM'.
	 * 
	 * NOTE, 'keep' has priority over 'store'. 
	 */

	if ((keep == GMRFLib_NEW_PROBLEM) && store) {
		/*
		 * create logicals 
		 */
		store_store_sub_graph = (store->sub_graph ? 0 : 1);
		store_use_sub_graph = !store_store_sub_graph;
		store_store_remap = (store->remap ? 0 : 1);
		store_use_remap = !store_store_remap;

		/*
		 * if the smtp is defined in the store, we use that one. otherwise, we use the default one. 
		 */
		if (GMRFLib_valid_smtp((int) store->smtp) == GMRFLib_TRUE) {
			smtp = store->smtp;
		} else {
			smtp = GMRFLib_smtp;
		}
		if (smtp == GMRFLib_SMTP_TAUCS) {
			store_store_symb_fact = (store->TAUCS_symb_fact ? 0 : 1);
			store_use_symb_fact = !store_store_symb_fact;
		} else if (smtp == GMRFLib_SMTP_PARDISO) {
			store_store_symb_fact = 0;
			store_use_symb_fact = !store_store_symb_fact;
		} else {
			store_store_symb_fact = 0;
			store_use_symb_fact = 0;
		}
	} else {
		if (!(keep == GMRFLib_NEW_PROBLEM)) {
			/*
			 * we use an old problem 
			 */
			smtp = (*problem)->sub_sm_fact.smtp;
		} else {
			/*
			 * use the default one 
			 */
			smtp = GMRFLib_smtp;
		}
	}

	/*
	 * if x = NULL, make it zeros 
	 */
	if (!x) {
		x = Calloc(graph->n, double);

		free_x = 1;
	}

	/*
	 * create new problem or reuse the old one 
	 */
	if (keep) {
		GMRFLib_ASSERT(*problem, GMRFLib_EPTR);
	} else {
		*problem = Calloc(1, GMRFLib_problem_tp);
	}

	/*
	 * this is the sparse-matrix method 
	 */
	(*problem)->sub_sm_fact.smtp = smtp;

	/*
	 * compute the internal stuff if needed 
	 */
	if (constr) {
		GMRFLib_EWRAP1(GMRFLib_prepare_constr(constr, graph, 0));
	}

	/*
	 * first, define the new graph (ok if fixed_value = NULL) 
	 */
	if (!(keep & GMRFLib_KEEP_graph)) {
		// if (keep && (*problem)->sub_graph)
		GMRFLib_free_graph((*problem)->sub_graph);     /* free old sub_graph */

		if (store_use_sub_graph) {
			/*
			 * copy from store 
			 */
			GMRFLib_copy_graph(&((*problem)->sub_graph), store->sub_graph);
		} else {
			/*
			 * compute it 
			 */
			GMRFLib_compute_subgraph(&((*problem)->sub_graph), graph, fixed_value);

			/*
			 * store a copy, if requested 
			 */
			if (store_store_sub_graph) {
				GMRFLib_copy_graph(&(store->sub_graph), (*problem)->sub_graph);
			}
		}
	}
	sub_n = (*problem)->sub_graph->n;
	if (sub_n == 0) {				       /* fast return if there is nothing todo */
		GMRFLib_free_graph((*problem)->sub_graph);
		Free(*problem);
		if (free_x) {
			Free(x);
		}
		GMRFLib_LEAVE_ROUTINE;
		return GMRFLib_SUCCESS;
	}

	/*
	 * the mapping is there in the graph, make a new pointer in the problem-definition 
	 */
	(*problem)->map = (*problem)->sub_graph->mothergraph_idx;

	/*
	 * compute the reordering which is stored in sub_sm_fact 
	 */
	if (!(keep & (GMRFLib_KEEP_chol | GMRFLib_KEEP_graph))) {
		/*
		 * it's ok to try to free the reordering if this is a new problem, it this case remap=NULL. 
		 */
		GMRFLib_free_reordering(&((*problem)->sub_sm_fact));
		if (store_use_remap) {
			/*
			 * use the reordering in store 
			 */
			(*problem)->sub_sm_fact.remap = Calloc(sub_n, int);
			memcpy((*problem)->sub_sm_fact.remap, store->remap, sub_n * sizeof(int));
			if (smtp == GMRFLib_SMTP_BAND) {
				(*problem)->sub_sm_fact.bandwidth = store->bandwidth;
			}
		} else {
			/*
			 * compute it 
			 */
			GMRFLib_EWRAP1(GMRFLib_compute_reordering(&((*problem)->sub_sm_fact), (*problem)->sub_graph, NULL));

			/*
			 * store a copy, if requested 
			 */
			if (store_store_remap) {
				if ((*problem)->sub_sm_fact.remap != NULL) {
					store->remap = Calloc(sub_n, int);
					memcpy(store->remap, (*problem)->sub_sm_fact.remap, sub_n * sizeof(int));
					if (smtp == GMRFLib_SMTP_BAND) {
						store->bandwidth = (*problem)->sub_sm_fact.bandwidth;
					}
				} else {
					store->remap = NULL;
				}
			}
		}
	}

	/*
	 * setup space 
	 */
	if (!keep) {
		(*problem)->sample = Calloc(graph->n, double);
		(*problem)->mean = Calloc(graph->n, double);
		(*problem)->mean_constr = Calloc(graph->n, double);
		(*problem)->sub_sample = Calloc(sub_n, double);
		(*problem)->sub_mean = Calloc(sub_n, double);
		(*problem)->sub_mean_constr = Calloc(sub_n, double);

		(*problem)->n = graph->n;		       /* for internal use only */
	} else {
		GMRFLib_ASSERT((*problem)->sample, GMRFLib_EPTR);
		GMRFLib_ASSERT((*problem)->mean, GMRFLib_EPTR);
		GMRFLib_ASSERT((*problem)->mean_constr, GMRFLib_EPTR);
		GMRFLib_ASSERT((*problem)->sub_sample, GMRFLib_EPTR);
		GMRFLib_ASSERT((*problem)->sub_mean, GMRFLib_EPTR);
		GMRFLib_ASSERT((*problem)->sub_mean_constr, GMRFLib_EPTR);
	}
	memcpy((*problem)->sample, x, graph->n * sizeof(double));
	memcpy((*problem)->mean, x, graph->n * sizeof(double));
	memcpy((*problem)->mean_constr, x, graph->n * sizeof(double));

	/*
	 * tabulate the Qfunc on the sub_graph. first, make the Qfunc, then tabulate it. 
	 */
	if (keep & GMRFLib_KEEP_chol) {
		/*
		 * use stored tab 
		 */
		GMRFLib_ASSERT((*problem)->tab, GMRFLib_EPTR);
	} else {
		if (keep) {
			GMRFLib_free_tabulate_Qfunc((*problem)->tab);
		}

		sub_Qfunc = GMRFLib_Qfunc_wrapper;
		sub_Qfunc_arg = Calloc(1, GMRFLib_Qfunc_arg_tp);
		sub_Qfunc_arg->map = (*problem)->map;	       /* yes, this ptr is needed */
		sub_Qfunc_arg->diagonal_adds = Calloc(sub_n, double);

		if (c) {
			for (i = 0; i < sub_n; i++) {
				sub_Qfunc_arg->diagonal_adds[i] = c[(*problem)->map[i]];
			}
		}
		sub_Qfunc_arg->user_Qfunc = Qfunc;
		sub_Qfunc_arg->user_Qfunc_args = Qfunc_args;
		GMRFLib_EWRAP1(GMRFLib_tabulate_Qfunc
			       (&((*problem)->tab), (*problem)->sub_graph, sub_Qfunc, (void *) sub_Qfunc_arg, NULL, NULL, NULL));

		Free(sub_Qfunc_arg->diagonal_adds);
		Free(sub_Qfunc_arg);
	}

	if (!(keep & GMRFLib_KEEP_mean)) {
		/*
		 * now compute the new 'effective' b, and then the mean 
		 */

		/*
		 * i use bb as name, pointing to the same storage as sub_mean. this to avoid yet another arraw. but i assume
		 * that sub_mean[] is zero, which is the case for !keep, but NOT the case otherwise. hence, i have to set it to 
		 * zero unless it already is so. [this was a nasty bug...] 
		 */
		bb = (*problem)->sub_mean;
		if (keep) {
			memset(bb, 0, sub_n * sizeof(double));
		}

		if (b) {
			for (i = 0; i < sub_n; i++) {
				bb[i] = b[(*problem)->map[i]];
			}
		}

		if (!fixed_value) {			       /* then sub_graph = graph */
			if (mean) {
				double *tmp = NULL;

				tmp = Calloc(sub_n, double);

				GMRFLib_Qx(tmp, mean, (*problem)->sub_graph, (*problem)->tab->Qfunc, (*problem)->tab->Qfunc_arg);
				for (i = 0; i < sub_n; i++) {
					bb[i] += tmp[i];
				}
				Free(tmp);
			}
		} else {
			/*
			 * x=(x1,x2), then x1|x2 has b = Q11 \mu1 - Q12(x2-\mu2) 
			 */
#pragma omp parallel for private(i, j, node, nnode)
			for (i = 0; i < sub_n; i++) {	       /* loop over all sub_nodes */
				GMRFLib_thread_id = id;
				node = (*problem)->map[i];
				if (mean) {
					bb[i] += (*((*problem)->tab->Qfunc)) (i, i, (*problem)->tab->Qfunc_arg) * mean[node];
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
							bb[i] -= qvalue * (x[nnode] - mean[nnode]);
						} else {
							bb[i] -= qvalue * x[nnode];
						}
					} else {
						/*
						 * nnone is not fixed 
						 */
						if (mean) {
							bb[i] += qvalue * mean[nnode];
						}
					}
				}
			}
			GMRFLib_thread_id = id;
		}
	}

	/*
	 * now compute the Cholesky-factorisation
	 * 
	 * solve to obtain the mean (recall that bb=(*problem)->sub_mean) later! 
	 */
	if (!(keep & GMRFLib_KEEP_chol)) {
		/*
		 * 13/5/2005. the next is a small hack, to be fixed properly later, well, it is not decided (yet,) if this is
		 * really needed.
		 * 
		 * in the case where graph is fixed and we're using TAUCS, we can reuse the symbolic factorisation. the
		 * symbolic, is roughly, about 10% of the costs of the numeric factorisation for large problems. however, the
		 * 'free_fact_sparse_...' does free both, and there is currently no options to free only the numerical one. so
		 * therefore i do this here. 
		 */
		if ((keep & GMRFLib_KEEP_graph) && (smtp == GMRFLib_SMTP_TAUCS)) {
			void *hold = NULL;

			hold = (void *) ((*problem)->sub_sm_fact.TAUCS_symb_fact);
			(*problem)->sub_sm_fact.TAUCS_symb_fact = NULL;

			GMRFLib_free_fact_sparse_matrix(&((*problem)->sub_sm_fact));
			(*problem)->sub_sm_fact.TAUCS_symb_fact = (supernodal_factor_matrix *) hold;
		} else {
			GMRFLib_free_fact_sparse_matrix(&((*problem)->sub_sm_fact));
		}

		if (store_use_symb_fact && (smtp == GMRFLib_SMTP_TAUCS)) {
			(*problem)->sub_sm_fact.TAUCS_symb_fact = GMRFLib_sm_fact_duplicate_TAUCS(store->TAUCS_symb_fact);
		}

		if (store_use_symb_fact && (smtp == GMRFLib_SMTP_PARDISO)) {
			// FIXME1("ADDED NEW EXPERIMENTAL CODE");
			GMRFLib_pardiso_store_tp *s = Calloc(1, GMRFLib_pardiso_store_tp);
			s->graph = (*problem)->sub_graph;
			// use the internal cached storage
			GMRFLib_duplicate_pardiso_store(&((*problem)->sub_sm_fact.PARDISO_fact), s, GMRFLib_FALSE, GMRFLib_FALSE);
			Free(s);
		}

		if (GMRFLib_catch_error_for_inla) {
			/*
			 * special version for INLA 
			 */
			int ret;
			ret = GMRFLib_build_sparse_matrix(&((*problem)->sub_sm_fact), (*problem)->tab->Qfunc,
							  (char *) ((*problem)->tab->Qfunc_arg), (*problem)->sub_graph);
			if (ret != GMRFLib_SUCCESS) {
				return ret;
			}

			ret = GMRFLib_factorise_sparse_matrix(&((*problem)->sub_sm_fact), (*problem)->sub_graph);
			if (ret != GMRFLib_SUCCESS) {
				return ret;
			}

		} else {
			/*
			 * plain version 
			 */
			GMRFLib_EWRAP1(GMRFLib_build_sparse_matrix(&((*problem)->sub_sm_fact), (*problem)->tab->Qfunc,
								   (char *) ((*problem)->tab->Qfunc_arg), (*problem)->sub_graph));
			GMRFLib_EWRAP1(GMRFLib_factorise_sparse_matrix(&((*problem)->sub_sm_fact), (*problem)->sub_graph));
		}

		if (store_store_symb_fact && (smtp == GMRFLib_SMTP_TAUCS)) {
			store->TAUCS_symb_fact = GMRFLib_sm_fact_duplicate_TAUCS((*problem)->sub_sm_fact.TAUCS_symb_fact);
		}
	}

	/*
	 * the next step is to initialize the constraints and take that part into account. 
	 */
	if (!(keep & GMRFLib_KEEP_constr)) {
		int fail = 1;

		/*
		 * free the old stuff from an old problem, if it exists 
		 */
		if (keep) {
			GMRFLib_free_constr((*problem)->sub_constr);
			Free((*problem)->constr_m);
			Free((*problem)->l_aqat_m);

			/*
			 * this is very special. if this flag is on, then we should reuse this vector if its exists. 
			 */
			if (!(keep & GMRFLib_UPDATE_constr)) {
				Free((*problem)->qi_at_m);
			}
		}

		if (constr && constr->nc > 0) {
			int nc, k, kk, det_computed = 0;
			double *aat_m, alpha, beta, *b_add = NULL, *aqat_m;

			/*
			 * first make new constraints on the subgraph. ok to call this function if fixed_values is NULL, then
			 * we just get a copy of constr back. 
			 */
			b_add = Calloc((*problem)->sub_graph->n, double);

			GMRFLib_EWRAP1(GMRFLib_recomp_constr
				       (&((*problem)->sub_constr), constr, x, b_add, fixed_value, graph, (*problem)->sub_graph));
			/*
			 * if we should keep the mean, then do not add the correction-terms 
			 */
			if (!(keep & GMRFLib_KEEP_mean)) {
				for (k = 0; k < (*problem)->sub_graph->n; k++) {
					(*problem)->sub_mean[k] += b_add[k];
				}
			}
			Free(b_add);

			if ((*problem)->sub_constr && (*problem)->sub_constr->nc > 0) {
				/*
				 * go further only if the constraint is still there: it might go away!!! 
				 */

				double *tmp_vector;

				nc = (*problem)->sub_constr->nc;	/* shortname */

				if (!STOCHASTIC_CONSTR((*problem)->sub_constr)) {
					/*
					 * this is for deterministic constraints only. the stochastic version will most likely
					 * be ok. 
					 */
					do {
						/*
						 * check that the constaints are not singular, if so, remove them 
						 */
						int debug = 0;
						int *map = NULL, rank, ii, jj;
						double *a = NULL, *e = NULL, eps = GMRFLib_eps(1. / 2.);

						nc = (*problem)->sub_constr->nc;
						if (debug) {
							printf("enter check with nc = %d\n", nc);
						}

						/*
						 * compute |A*A'| 
						 */
						alpha = 1.0;
						beta = 0.0;
						aat_m = Calloc(nc * nc, double);

						dgemm_("N", "T", &nc, &nc, &sub_n, &alpha, (*problem)->sub_constr->a_matrix,
						       &nc, (*problem)->sub_constr->a_matrix, &nc, &beta, aat_m, &nc, 1, 1);
						GMRFLib_EWRAP1(GMRFLib_comp_chol_semidef
							       (NULL, &map, &rank, aat_m, nc, &((*problem)->logdet_aat), eps));
						GMRFLib_ASSERT(rank != 0, GMRFLib_EPARAMETER);

						if (rank == nc) {
							/*
							 * ok, the reduced constraints are fine 
							 */
							fail = 0;
							det_computed = 1;	/* flag that (*problem)->logdet_aat is computed */
						} else {
							/*
							 * oops, the reduced constraints are singular. this have to be fixed 
							 */

							fail = 1;

							if (debug) {
								printf("rank is estimated to be %d < nc=%d\n", rank, nc);
								for (ii = 0; ii < rank; ii++) {
									printf("map[%1d] = %1d\n", ii, map[ii]);
								}
							}

							/*
							 * now we need to store only the 'non-singular' columns 
							 */
							a = Calloc(sub_n * rank, double);
							e = Calloc(rank, double);

							for (jj = 0; jj < rank; jj++) {
								kk = map[jj];
								for (ii = 0; ii < sub_n; ii++) {
									a[ii * rank + jj] = (*problem)->sub_constr->a_matrix[ii * nc + kk];
								}
								e[jj] = (*problem)->sub_constr->e_vector[kk];
							}

							/*
							 * free the terms and then copy the reduced ones 
							 */
							Free((*problem)->sub_constr->a_matrix);
							(*problem)->sub_constr->a_matrix = a;
							Free((*problem)->sub_constr->e_vector);
							(*problem)->sub_constr->e_vector = e;
							(*problem)->sub_constr->nc = rank;

							if (debug) {
								GMRFLib_print_constr(stdout, (*problem)->sub_constr, (*problem)->sub_graph);
							}
						}
						Free(map);
						Free(aat_m);   /* free temp storage */
					}
					while (fail);
				}

				/*
				 * then compute the important constr_matrix 
				 */

				double *qi_at_m_store = NULL;  /* possible reuse old results */

				if ((keep & GMRFLib_UPDATE_constr)) {
					if (nc == 1) {
						Free((*problem)->qi_at_m);
					} else {
						if ((*problem)->qi_at_m) {
							qi_at_m_store = (*problem)->qi_at_m;
							(*problem)->qi_at_m = NULL;
						}
					}
				}

				(*problem)->qi_at_m = Calloc(nc * sub_n, double);
				if (qi_at_m_store == NULL) {
					/*
					 * compute it as usual 
					 */
					if (0) {
						for (k = 0; k < nc; k++) {
							kk = k * sub_n;
							for (i = 0; i < sub_n; i++) {
								(*problem)->qi_at_m[i + kk] = (*problem)->sub_constr->a_matrix[k + nc * i];
							}
							GMRFLib_solve_llt_sparse_matrix(&((*problem)->qi_at_m[kk]), 1,
											&((*problem)->sub_sm_fact), (*problem)->sub_graph);
						}
					} else {
						for (k = 0; k < nc; k++) {
							kk = k * sub_n;
							for (i = 0; i < sub_n; i++) {
								(*problem)->qi_at_m[i + kk] = (*problem)->sub_constr->a_matrix[k + nc * i];
							}
						}
						GMRFLib_solve_llt_sparse_matrix((*problem)->qi_at_m, nc,
										&((*problem)->sub_sm_fact), (*problem)->sub_graph);
					}
				} else {
					/*
					 * reuse 
					 */
					if (0) {
						memcpy((*problem)->qi_at_m, qi_at_m_store, (nc - 1) * sub_n * sizeof(double));
						for (k = nc - 1; k < nc; k++) {
							kk = k * sub_n;
							for (i = 0; i < sub_n; i++) {
								(*problem)->qi_at_m[i + kk] = (*problem)->sub_constr->a_matrix[k + nc * i];
							}
							GMRFLib_solve_llt_sparse_matrix(&((*problem)->qi_at_m[kk]), 1,
											&((*problem)->sub_sm_fact), (*problem)->sub_graph);
						}
					} else {
						memcpy((*problem)->qi_at_m, qi_at_m_store, (nc - 1) * sub_n * sizeof(double));
						for (k = nc - 1; k < nc; k++) {
							kk = k * sub_n;
							for (i = 0; i < sub_n; i++) {
								(*problem)->qi_at_m[i + kk] = (*problem)->sub_constr->a_matrix[k + nc * i];
							}
						}
						GMRFLib_solve_llt_sparse_matrix(&((*problem)->qi_at_m[(nc - 1) * sub_n]), 1,
										&((*problem)->sub_sm_fact), (*problem)->sub_graph);
					}
				}
				Free(qi_at_m_store);

				/*
				 * compute l_aqat_m = chol(AQ^{-1}A^T)^{-1}) = chol(A qi_at_m)^{-1}, size = nc x nc 
				 */
				aqat_m = Calloc(nc * nc, double);

				alpha = 1.0;
				beta = 0.0;
				dgemm_("N", "N", &nc, &nc, &sub_n, &alpha, (*problem)->sub_constr->a_matrix, &nc,
				       (*problem)->qi_at_m, &sub_n, &beta, aqat_m, &nc, 1, 1);

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

				if (0) {
					int iii, jjj;
					FIXME("print aqat_m");
					for (iii = 0; iii < nc; iii++) {
						for (jjj = 0; jjj < nc; jjj++) {
							printf(" %.8g", aqat_m[iii + jjj * nc]);
						}
						printf("\n");
					}
				}

				/*
				 * compute chol(aqat_m), recall that GMRFLib_comp_chol_general returns a new malloced L 
				 */
				GMRFLib_EWRAP1(GMRFLib_comp_chol_general
					       (&((*problem)->l_aqat_m), aqat_m, nc, &((*problem)->logdet_aqat), GMRFLib_ESINGCONSTR));
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

				if (!STOCHASTIC_CONSTR((*problem)->sub_constr) && !det_computed) {
					/*
					 * compute |A*A'| 
					 */
					alpha = 1.0;
					beta = 0.0;
					aat_m = Calloc(nc * nc, double);

					dgemm_("N", "T", &nc, &nc, &sub_n, &alpha, (*problem)->sub_constr->a_matrix,
					       &nc, (*problem)->sub_constr->a_matrix, &nc, &beta, aat_m, &nc, 1, 1);
					tmp_vector = NULL;
					GMRFLib_EWRAP1(GMRFLib_comp_chol_general
						       (&tmp_vector, aat_m, nc, &((*problem)->logdet_aat), GMRFLib_ESINGCONSTR));
					Free(aat_m);
					Free(tmp_vector);
				}
			}
		}
	}

	if (!(keep & GMRFLib_KEEP_mean)) {
		GMRFLib_EWRAP1(GMRFLib_solve_llt_sparse_matrix((*problem)->sub_mean, 1, &((*problem)->sub_sm_fact), (*problem)->sub_graph));
		if (!((*problem)->sub_mean_constr)) {
			(*problem)->sub_mean_constr = Calloc(sub_n, double);
		}
		memcpy((*problem)->sub_mean_constr, (*problem)->sub_mean, sub_n * sizeof(double));
	}

	if (((*problem)->sub_constr && (*problem)->sub_constr->nc > 0)) {
		if (((keep & GMRFLib_KEEP_constr)) && ((keep & GMRFLib_KEEP_mean))) {
			/*
			 * do nothing 
			 */
		} else {
			/*
			 * compute the mean after correcting for the constraint. this is the same as if the sample itself is
			 * sub_mean! i have copied parts of code from GMRFLib_sample into here...
			 */
			int nc = constr->nc, inc = 1;
			double alpha, beta, *t_vector;

			Free((*problem)->sub_constr_value);
			(*problem)->sub_constr_value = t_vector = Calloc(nc, double);

			GMRFLib_EWRAP1(GMRFLib_eval_constr(t_vector, NULL, (*problem)->sub_mean, (*problem)->sub_constr, (*problem)->sub_graph));

			/*
			 * sub_mean_constr is pr.default equal to sub_mean 
			 */
			alpha = -1.0;
			beta = 1.0;			       /* mean_constr = mean - cond_m*t_vector */
			/*
			 * workaround for a 'bug' in goto-blas and mkl(intel); which is triggered in one example. The failsafe version is just the fortran code for
			 * dgemv(), but renamed to it does not used the one in another blas library. 
			 */
			if (0) {
				dgemv_("N", &sub_n, &nc, &alpha, (*problem)->constr_m, &sub_n, t_vector, &inc, &beta,
				       (*problem)->sub_mean_constr, &inc, 1);
			} else {
				dgemv_failsafe_("N", &sub_n, &nc, &alpha, (*problem)->constr_m, &sub_n, t_vector, &inc, &beta,
						(*problem)->sub_mean_constr, &inc, 1);
			}
		}
	}

	/*
	 * make the ``mean variables'', which are easier accessible for the user. 
	 */
	for (i = 0; i < sub_n; i++) {
		j = (*problem)->map[i];
		(*problem)->mean[j] = (*problem)->sub_mean[i];
		(*problem)->mean_constr[j] = (*problem)->sub_mean_constr[i];
	}
	/*
	 * evaluate the log-likelihood and compute the constants therein 
	 */
	GMRFLib_EWRAP1(GMRFLib_evaluate__intern(*problem, 1));

	if (free_x) {
		Free(x);
	}

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

/*!
  \brief Samples one realization of the elements of a GMRF <em>\b x</em>.

  Whether the sampling is constrained or unconstrained, or should be done 
  conditionally on fixed values, is defined by the \c GMRFLib_problem_tp -argument.

  \param[in,out] problem At input \a problem should contain the problem 
  specification, as defined by a call to \c GMRFLib_init_problem(). At output, 
  a sample of the GMRF has been generated, and is stored in the \a sample
  member of the data structure \a problem. Also, the log-density of the sample 
  is computed, and stored in the \a sub_logdens member. If \a problem
  = \c NULL, nothing is done, and the routine returns 0.

  \remark To sample one realization from a GMRF <em>\b x</em>, first invoke
  \c GMRFLib_init_problem() generating the problem specification by initializing 
  a \c GMRFLib_problem_tp -object \a problem, and then call \c GMRFLib_sample() 
  using \a problem as an argument. Repeated samples are generated by 
  repeated calls to \c GMRFLib_sample(), using the same problem specification 
  object. There is no need to call \c GMRFLib_init_problem() more than once for 
  each sampling problem. The log-density of the sample is computed within the 
  sampling routine, but can also be computed explicitly for a 
  \c GMRFLib_problem_tp -object by calling the routine \c GMRFLib_evaluate(). \n
  Prior to sampling, the built-in random number generator has to be initialised. 
  This is done by calling the function (<tt>*GMRFLib_uniform_init</tt>), taking 
  one integer argument, specifying the seed of the generator. The user might 
  replace the default random generator by another, see globals.h
  for details. \n
  At output, the sampling results are stored in the \c GMRFLib_problem_tp
  -object. No functions are provided for extracting the results of the 
  sampling, so the results are retrieved by explicitly accessing the members 
  of \a problem after running the routine. The sampled values are stored in 
  <tt>problem->sample</tt> and the corresponding log-density in 
  <tt>problem->sub_logdens</tt>. The fixed values are copied into the 
  corresponding elements in \a sample.

  \par Example:
  See \ref ex_problem-setup and \n \ref ex_wa.
  \sa GMRFLib_init_problem, GMRFLib_evaluate.
*/
int GMRFLib_sample(GMRFLib_problem_tp * problem)
{
	int i, n;
	double sqrterm;

	if (!problem)
		return GMRFLib_SUCCESS;
	GMRFLib_ENTER_ROUTINE;

	n = problem->sub_graph->n;

	/*
	 * just solve L^Tx=z, then add mean and store in correct place 
	 */
	for (i = 0, sqrterm = 0.0; i < n; i++) {
		double z;

		z = GMRFLib_stdnormal();
		sqrterm += SQR(z);
		problem->sub_sample[i] = z;
	}

	GMRFLib_EWRAP1(GMRFLib_solve_lt_sparse_matrix(problem->sub_sample, 1, &(problem->sub_sm_fact), problem->sub_graph));
	for (i = 0; i < n; i++) {
		problem->sub_sample[i] += problem->sub_mean[i];
		problem->sample[problem->map[i]] = problem->sub_sample[i];	/* will be modified later if constraints */
	}

	/*
	 * if there is no constraints, then we can evaluate the log-density and return 
	 */
	if (!problem->sub_constr || (problem->sub_constr && problem->sub_constr == 0)) {
		problem->sub_logdens = -0.5 * n * log(2.0 * M_PI) + problem->log_normc - 0.5 * sqrterm;
		GMRFLib_LEAVE_ROUTINE;
		return GMRFLib_SUCCESS;
	} else {
		/*
		 * otherwise, we must account for the constraint, BUT, we need to call GMRFLib_evaluate() to compute the
		 * log-density.
		 * 
		 * FIXME: There is a slight overhead, since the quadratic form can be faster evaluated knowing the random
		 * bits... 
		 */

		double alpha, beta, *t_vector;
		int inc = 1, nc = problem->sub_constr->nc;

		Free(problem->sub_constr_value);
		problem->sub_constr_value = Calloc(nc, double);
		t_vector = Calloc(nc, double);		       /* t_vector = Ax-e */

		GMRFLib_EWRAP1(GMRFLib_eval_constr(problem->sub_constr_value, NULL, problem->sub_sample, problem->sub_constr, problem->sub_graph));
		memcpy(t_vector, problem->sub_constr_value, nc * sizeof(double));

		if (STOCHASTIC_CONSTR(problem->sub_constr)) {
			/*
			 * add (YES: do not set t_vector to zero!) noisy terms here 
			 */
			double *z = NULL;

			z = Calloc(nc, double);

			for (i = 0; i < nc; i++) {
				z[i] = GMRFLib_stdnormal();
			}

			if (problem->sub_constr->errcov_diagonal) {
				/*
				 * add L*z, where L=sqrt(diag(Q)) 
				 */
				for (i = 0; i < nc; i++) {
					t_vector[i] += problem->sub_constr->intern->chol[i] * z[i];
				}
			} else {
				/*
				 * add L*z 
				 */
				dtrmv_("L", "N", "N", &nc, problem->sub_constr->intern->chol, &nc, z, &inc, 1, 1, 1);
				for (i = 0; i < nc; i++) {
					t_vector[i] += z[i];
				}
			}
			Free(z);
		}
		alpha = -1.0;
		beta = 1.0;				       /* sample := sample - cond_m*t_vector */
		dgemv_("N", &n, &nc, &alpha, problem->constr_m, &n, t_vector, &inc, &beta, problem->sub_sample, &inc, 1);
		Free(t_vector);

		for (i = 0; i < n; i++) {
			problem->sample[problem->map[i]] = problem->sub_sample[i];	/* put into place */
		}
		GMRFLib_EWRAP1(GMRFLib_evaluate(problem));     /* to compute the log-density */
	}

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

/*! \brief Evaluates the log-density of a sample.

  Just compute the part from the sample in the log-likelihood, the constants 
  are already computed. \c GMRFLib_evaluate() is called from within the sampling 
  routine \c GMRFLib_sample(), see also \ref log_dens.

  \param[in,out] problem At input \a problem should contain the problem 
  specification and a sample, as defined by a call to \c GMRFLib_init_problem() 
  followed by one or more calls to \c GMRFLib_sample(). At output, the log-density
  of the current sample has been computed, and is stored in the member
  \a sub_logdens of \a problem.

  \remark To extract the value of the log-density, the user should access the 
  member \a sub_logdens of \a problem directly.

  \sa GMRFLib_init_problem, GMRFLib_sample
*/
int GMRFLib_evaluate(GMRFLib_problem_tp * problem)
{
	GMRFLib_ENTER_ROUTINE;
	GMRFLib_EWRAP1(GMRFLib_evaluate__intern(problem, 0));
	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}
int GMRFLib_evaluate__intern(GMRFLib_problem_tp * problem, int compute_const)
{
	/*
	 * evaluate the log-density in point 'sample' in the problem definition 
	 */

	int i, n;
	double sqrterm, *xx = NULL, *yy = NULL;

	if (!problem) {
		return GMRFLib_SUCCESS;
	}

	n = problem->sub_graph->n;
	xx = Calloc(n, double);
	yy = Calloc(n, double);

	/*
	 * user has altered the 'sample', put the correct subset into sub_sample and compute (x-\mu)^TQ(x-\mu)
	 */
	for (i = 0; i < n; i++) {
		problem->sub_sample[i] = problem->sample[problem->map[i]];
		xx[i] = problem->sub_sample[i] - problem->sub_mean[i];
	}
	GMRFLib_Qx(yy, xx, problem->sub_graph, problem->tab->Qfunc, (void *) problem->tab->Qfunc_arg);
	for (i = 0, sqrterm = 0.0; i < n; i++) {
		sqrterm += yy[i] * xx[i];
	}

	Free(xx);
	Free(yy);

	/*
	 * evaluate the normalization constant and add up 
	 */
	if (compute_const) {
		GMRFLib_EWRAP0(GMRFLib_log_determinant(&(problem->log_normc), &(problem->sub_sm_fact), problem->sub_graph));
		problem->log_normc /= 2.0;		       /* |Q|^1/2 */
	}

	problem->sub_logdens = -0.5 * n * log(2.0 * M_PI) + problem->log_normc - 0.5 * sqrterm;

	/*
	 * now correct for the constraint, if any 
	 */
	if (!problem->sub_constr || problem->sub_constr->nc == 0) {
		return GMRFLib_SUCCESS;
	} else if (STOCHASTIC_CONSTR(problem->sub_constr)) {
		/*
		 * stochastic constraints 
		 */
		int nc = problem->sub_constr->nc;
		double exp_corr2;

		if (compute_const) {
			/*
			 * t_vector = A mu -b tt_vector = work 
			 */
			double *t_vector = NULL, *tt_vector = NULL, exp_corr;

			Free(problem->sub_constr_value);
			problem->sub_constr_value = t_vector = Calloc(nc, double);
			tt_vector = Calloc(nc, double);

			GMRFLib_EWRAP0(GMRFLib_eval_constr(t_vector, NULL, problem->sub_mean, problem->sub_constr, problem->sub_graph));


			GMRFLib_EWRAP0(GMRFLib_solveAxb_posdef(tt_vector, problem->l_aqat_m, t_vector, nc, 1));
			for (i = 0, exp_corr = 0.0; i < nc; i++) {
				exp_corr += t_vector[i] * tt_vector[i];
			}
			problem->exp_corr = exp_corr;
			Free(tt_vector);
		}

		GMRFLib_EWRAP0(GMRFLib_eval_constr(NULL, &exp_corr2, problem->sub_sample, problem->sub_constr, problem->sub_graph));

		/*
		 * [x|Ax] = [x] [Ax|x] / [Ax] 
		 */
		problem->sub_logdens += -0.5 * *(problem->sub_constr->intern->logdet) - 0.5 * exp_corr2	/* [Ax|x] */
		    - (-0.5 * problem->logdet_aqat - 0.5 * problem->exp_corr);	/* [Ax] */
	} else {
		/*
		 * deterministic constraints 
		 */
		int nc = problem->sub_constr->nc;

		if (compute_const) {
			/*
			 * t_vector = A mu-b tt_vector = i_aqat_m*t_vector 
			 */
			double *t_vector = NULL, *tt_vector = NULL, exp_corr;

			Free(problem->sub_constr_value);
			problem->sub_constr_value = t_vector = Calloc(nc, double);
			tt_vector = Calloc(nc, double);

			GMRFLib_EWRAP0(GMRFLib_eval_constr(t_vector, NULL, problem->sub_mean, problem->sub_constr, problem->sub_graph));
			GMRFLib_EWRAP0(GMRFLib_solveAxb_posdef(tt_vector, problem->l_aqat_m, t_vector, nc, 1));
			for (i = 0, exp_corr = 0.0; i < nc; i++) {
				exp_corr += t_vector[i] * tt_vector[i];
			}
			problem->exp_corr = exp_corr;
			Free(tt_vector);
		}

		/*
		 * [x|Ax] = [x] [Ax|x] / [Ax] 
		 */
		problem->sub_logdens += (-0.5 * problem->logdet_aat)	/* [Ax|x] */
		    -(-0.5 * nc * log(2.0 * M_PI) - 0.5 * problem->logdet_aqat - 0.5 * problem->exp_corr);	/* [Ax] */
	}

	return GMRFLib_SUCCESS;
}

/*!
  \brief Free all malloced stuff in 'problem'
*/
int GMRFLib_free_problem(GMRFLib_problem_tp * problem)
{
	/*
	 * free all malloced stuff in 'problem' 
	 */
	if (!problem) {
		return GMRFLib_SUCCESS;
	}
	GMRFLib_free_Qinv(problem);			       /* do this here, so that `n' is still there... */

	Free(problem->sample);
	Free(problem->mean);
	Free(problem->mean_constr);
	Free(problem->sub_sample);
	Free(problem->sub_mean);
	Free(problem->sub_mean_constr);
	GMRFLib_free_fact_sparse_matrix(&(problem->sub_sm_fact));
	GMRFLib_free_reordering(&(problem->sub_sm_fact));
	Free(problem->sub_constr_value);
	Free(problem->constr_m);
	Free(problem->l_aqat_m);
	Free(problem->inv_aqat_m);
	Free(problem->qi_at_m);
	GMRFLib_free_graph(problem->sub_graph);
	GMRFLib_free_tabulate_Qfunc(problem->tab);

	GMRFLib_free_constr(problem->sub_constr);
	problem->sub_constr = NULL;

	Free(problem);
	return GMRFLib_SUCCESS;
}

int GMRFLib_free_Qinv(GMRFLib_problem_tp * problem)
{
	if (problem && problem->sub_inverse) {
		int i, n = problem->sub_graph->n;

		for (i = 0; i < n; i++) {
			map_id_free(problem->sub_inverse->Qinv[i]);
			Free(problem->sub_inverse->Qinv[i]);
		}
		Free(problem->sub_inverse->Qinv);

		map_ii_free(problem->sub_inverse->mapping);
		Free(problem->sub_inverse->mapping);
		Free(problem->sub_inverse);
	}
	return GMRFLib_SUCCESS;
}

/*!
  \brief Free storage \c store which holds temporary calculations which can be reused
  
  \param[in] store  The pointer to the (malloced) \c GMRFLib_store_tp object

  \sa GMRFLib_blockupdate, GMRFLib_init_problem, GMRFLib_optimize, GMRFLib_init_problem_hidden,
  GMRFLib_init_GMRF_approximation

  \note \a store itself also free'ed.
*/
int GMRFLib_free_store(GMRFLib_store_tp * store)
{
	/*
	 * free contents in store 
	 */
	if (!store) {
		return GMRFLib_SUCCESS;
	}

	Free(store->remap);
	if (store->copy_ptr) {
		/*
		 * do nothing 
		 */
	} else {
		GMRFLib_free_graph(store->sub_graph);
		if (store->TAUCS_symb_fact) {
			taucs_supernodal_factor_free(store->TAUCS_symb_fact);
		}
	}

	if (store->copy_pardiso_ptr) {
		/*
		 * do nothing 
		 */
	} else {
		if (store->PARDISO_fact) {
			GMRFLib_pardiso_free(&(store->PARDISO_fact));
			store->PARDISO_fact = NULL;
		}
	}

	store->sub_graph = NULL;
	store->TAUCS_symb_fact = NULL;

	/*
	 * free the diag and sub-store. its of the same type, therefore we can do this recursively. 
	 */
	GMRFLib_free_store(store->sub_store);
	GMRFLib_free_store(store->diag_store);

	if (store->store_problems) {
		GMRFLib_free_problem(store->problem_old2new);
		GMRFLib_free_problem(store->problem_new2old);
		Free(store->old_logdens);
		Free(store->new_logdens);
	}

	/*
	 * free the ptr itself 
	 */
	Free(store);

	return GMRFLib_SUCCESS;
}

/*!
  \brief Compute elements in the inverse of Q

  \param[in,out] problem Compute elements in the inverse of the Qmatrix defined in the \ref
  GMRFLib_problem_tp object \c problem, and store the result in \c problem. Hard and soft
  constraints are taken into account. Use repeatedly the function \ref GMRFLib_Qinv_get() to
  retrieve elements of the inverse of Q.

  \param[in] storage The user can decide how many of the terms in the inverse computed that are stored, and how they are
  computed. The options are
  
  - \c GMRFLib_QINV_ALL    Store all elements (default)
  - \c GMRFLib_QINV_NEIGB  Store elements (i,j) such that i=j or i~j
  - \c GMRFLib_QINV_DIAG   Store only the diagonal
  - \c GMRFLib_QINV_NO_CHECK Disable checking for ``complete'' \c L. Can be OR'ed with any of the
                           first three options and has only effect for the TAUCS-solver.
  - \c GMRFLib_QINV_CHECK_ONCE Check once only  for ``complete'' \c L. Can be OR'ed with any of the
                           first three options and has only effect for the TAUCS-solver.

  \note There are EXTRA CPU COST related to reduced storage, as all elements corresponding to \c
  GMRFLib_QINV_ALL must be computed in any case.

  \note For the TAUCS-solver only: There is reduced CPU cost using \c GMRFLib_QINV_NO_CHECK or \c GMRFLib_QINV_CHECK_ONCE, but a
  small error (usually) (less for \c GMRFLib_CHECK_ONCE) can be introduced in some of the covariances. The error is most notable
  for those ``far away'' but this is not always the case. The default behaviour is to check until a complete \c L is found, and
  this approach is without error.

  \note For the BAND-solver: The complete \c L is always used. 

  \note Example
  \verbatim
   GMRFLib_Qinv(problem, GMRFLib_QINV_ALL);
   for(i=0;i<graph->n;i++)
   {
       double *val = NULL;
       val = GMRFLib_Qinv_get(problem, i, i);
       printf("Marginal variance for x[%1d] = %f\n", i, *val);
   }
   \endverbatim

   \sa GMRFLib_Qinv_get()

*/
int GMRFLib_Qinv(GMRFLib_problem_tp * problem, int storage)
{
	GMRFLib_ENTER_ROUTINE;
	if (problem) {
		GMRFLib_EWRAP1(GMRFLib_compute_Qinv((void *) problem, storage));
	}
	GMRFLib_LEAVE_ROUTINE;

	return GMRFLib_SUCCESS;
}

/*!
  \brief Get entry (i,j) in the computed Qinv.

  \param[in] problem  The \c problem for which \c Qinv is stored.
  \param[in] i Index
  \param[in] j Index

  This function returns a pointer to the values of \c Qinv(i,j) and \c NULL if (i,j) is illegal, not
  computed or not stored. All entries corresponding to non-zero entries in \c L are computed.

  \note Se \ref GMRFLib_Qinv() for an example.
*/
double *GMRFLib_Qinv_get(GMRFLib_problem_tp * problem, int i, int j)
{
	int *ii = NULL, *jj = NULL;

	if (!problem || !problem->sub_inverse) {
		return NULL;
	}

	ii = map_ii_ptr(problem->sub_inverse->mapping, i);
	if (!ii) {
		return NULL;
	}
	jj = map_ii_ptr(problem->sub_inverse->mapping, j);
	if (!jj) {
		return NULL;
	}
	return map_id_ptr(problem->sub_inverse->Qinv[IMIN(*ii, *jj)], IMAX(*ii, *jj));
}

/*! \brief Create an empty \c GMRFLib_constr_tp -object.

  \param[in,out] constr A pointer to an \c GMRFLib_constr_tp object. At output, \a constr points to
  an empty \c GMRFLib_constr_tp -object. The internal pointer member \a intern points to an empty \c
  GMRFLib_constr__intern_tp -object, the other members are uninitialized.
  
  \remarks This functions allocates the memory needed by the \c GMRFLib_constr_tp -object and the \c
  GMRFLib_constr__intern_tp -object. To initialize the members of <em>(*constr)</em>, allocate \a
  a_matrix and \a e_vector directly, using \c calloc, define \em nc and call \c
  GMRFLib_prepare_constr().
  
  \sa GMRFLib_free_constr, GMRFLib_prepare_constr.
 */
int GMRFLib_make_empty_constr(GMRFLib_constr_tp ** constr)
{
	if (constr) {
		*constr = Calloc(1, GMRFLib_constr_tp);
		(*constr)->intern = Calloc(1, GMRFLib_constr__intern_tp);
	}
	return GMRFLib_SUCCESS;
}

/*! \brief Free the memory held by a \c GMRFLib_constr_tp -object. 

  \note To ensure safe memory handling, this function should only be applied 
  to pointers allocated by using the library routines \c GMRFLib_make_empty_constr() 
  and \c GMRFLib_prepare_constr(). Also, as long as \em nc > 0, the arrays 
  \a a_matrix and \a e_vector should both be <tt>!NULL</tt>, and should have
  been dynamically allocated (using \c calloc).

  \param[in,out] constr A pointer to a \c GMRFLib_constr_tp -object. At output, 
  the pointer \a constr and it's member pointers are deallocated.

  \sa GMRFLib_make_empty_constr
*/
int GMRFLib_free_constr(GMRFLib_constr_tp * constr)
{
	if (constr) {
		Free(constr->a_matrix);
		Free(constr->e_vector);
		Free(constr->errcov_diagonal);
		Free(constr->errcov_general);
		Free(constr->intern->chol);
		Free(constr->intern->Qpattern);
		Free(constr->intern->logdet);
		Free(constr->intern);
		Free(constr);
	}
	return GMRFLib_SUCCESS;
}

/*! \brief Prints the available information on a constraint held by an
  \c GMRFLib_constr_tp -object.

  \param fp The \c FILE* on which to print the constraints.
  \param[in] constr A pointer to a \c GMRFLib_constr_tp -object.
  \param[in] graph The graph on which the constraints are defined.

  \sa GMRFLib_make_empty_constr
*/
int GMRFLib_print_constr(FILE * fp, GMRFLib_constr_tp * constr, GMRFLib_graph_tp * graph)
{
	int i, j;
	FILE *fpp = NULL;

	if (!constr) {
		return GMRFLib_SUCCESS;
	}
	fpp = (fp ? fp : stdout);

	fprintf(fpp, "n_constr %d\n", constr->nc);
	for (i = 0; i < constr->nc; i++) {
		fprintf(fpp, "constraint %d, e= %f\na[%1d, ] = ", i, constr->e_vector[i], i);
		for (j = 0; j < graph->n; j++) {
			fprintf(fpp, " %9.6f", constr->a_matrix[i + j * constr->nc]);
		}
		fprintf(fpp, "\n");
	}
	if (STOCHASTIC_CONSTR(constr)) {
		if (constr->errcov_diagonal) {
			fprintf(fpp, "errcov_diagonal ");
			for (i = 0; i < constr->nc; i++) {
				fprintf(fpp, " %12.6f", constr->errcov_diagonal[i]);
			}
			fprintf(fpp, "\n");
			if (constr->intern && constr->intern->chol) {
				fprintf(fpp, "intern->chol ");
				for (i = 0; i < constr->nc; i++) {
					fprintf(fpp, " %12.6f", constr->intern->chol[i]);
				}
				fprintf(fpp, "\n");
			}
		} else {
			fprintf(fpp, "errcov_general\n");
			for (i = 0; i < constr->nc; i++) {
				fprintf(fpp, "%1d: ", i);
				for (j = 0; j <= i; j++) {
					fprintf(fpp, " %12.6f", constr->errcov_general[i + j * constr->nc]);
				}
				fprintf(fpp, "\n");
			}
			if (constr->intern && constr->intern->chol) {
				fprintf(fpp, "intern->chol\n");
				for (i = 0; i < constr->nc; i++) {
					fprintf(fpp, "%1d: ", i);
					for (j = 0; j <= i; j++) {
						fprintf(fpp, " %12.6f", constr->intern->chol[i + j * constr->nc]);
					}
					fprintf(fpp, "\n");
				}
			}
		}
	}

	return GMRFLib_SUCCESS;
}

/*! \brief Prepare a constraint by computing the information in the \em intern 
  object, if needed.

  \param[in,out] constr At input, a pointer to a \c GMRFLib_constr_tp -object, 
  created by \c GMRFLib_make_empty_constr().  At output, internal information needed
  by the sampling routines is computed and added to the object.
  \param[in] graph The graph on which the constraint is defined.
  \param scale_constr If TRUE, then scale the constraint so that
  \f$ \max_j A_{ij}=1 \f$ for each \em i, and correct the right hand 
  side/covaraince accordingly. This is only needed for numerical
  reasons in some cases. If FALSE, do not scale.

  \remarks The function computes the Cholesky factorization and the log of 
  the determinant of the covariance matrix 
  \f$ \mbox{\boldmath $\Sigma$}_{\epsilon} \f$.  The matrix <em>\b A</em> and 
  the vector <em>\b e</em> of the constraint <em>\b Ax = \b e</em> should be
  allocated and initialised explicitly. This can be done before or after 
  calling \c GMRFLib_prepare_constr().
  
  \sa GMRFLib_make_empty_constr, GMRFLib_free_constr
 */
int GMRFLib_prepare_constr(GMRFLib_constr_tp * constr, GMRFLib_graph_tp * graph, int scale_constr)
{
	/*
	 * prepare a constraint: compute the internal variables in constr->intern
	 * 
	 * NOTE: at the moment, the intern->Q is only used for checking if Q[i,j] is zero or not. 
	 */
	int i, nc, k, kk, n, debug = 0;
	double ldet, *scale = NULL;

	if (!constr) {
		return GMRFLib_SUCCESS;
	}

	nc = constr->nc;
	n = graph->n;

	if (scale_constr) {

		/*
		 * first, scale the constraints so that max(|A[i,]|)=1 
		 */
		int need_scaling = 0;

		scale = Calloc(nc, double);

		for (k = 0; k < nc; k++) {
			for (i = 0; i < n; i++) {
				scale[k] = DMAX(scale[k], ABS(constr->a_matrix[i * nc + k]));
			}
		}

		if (debug) {
			for (k = 0; k < nc; k++) {
				printf("scale A[%1d,] .... %f\n", k, scale[k]);
			}
		}

		for (k = 0, need_scaling = 0; k < nc; k++) {
			need_scaling = (need_scaling || (scale[k] != 1.0));
		}

		if (need_scaling) {
			if (debug) {
				printf("need scaling\n");
			}

			for (k = 0; k < nc; k++) {
				for (i = 0; i < n; i++) {
					constr->a_matrix[i * nc + k] /= scale[k];
				}
				constr->e_vector[k] /= scale[k];
			}

			if (STOCHASTIC_CONSTR(constr)) {
				if (constr->errcov_diagonal) {
					for (k = 0; k < nc; k++) {
						constr->errcov_diagonal[k] /= SQR(scale[k]);
					}
				}
				if (constr->errcov_general) {
					for (k = 0; k < nc; k++) {
						for (kk = 0; kk < nc; kk++) {
							constr->errcov_general[k + kk * nc] /= (scale[k] * scale[kk]);
						}
					}
				}

				if (constr->intern) {
					/*
					 * if we have scaled, then recompute these [can be avoided but....] 
					 */
					Free(constr->intern->logdet);
					Free(constr->intern->chol);
					Free(constr->intern->Qpattern);
				}
			}
		}
		Free(scale);
		if (debug) {
			printf("scaled constr\n");
			GMRFLib_print_constr(stdout, constr, graph);
		}
	}

	if (!STOCHASTIC_CONSTR(constr)) {
		return GMRFLib_SUCCESS;
	}
	if (constr->intern && constr->intern->chol && constr->intern->Qpattern && constr->intern->logdet) {
		return GMRFLib_SUCCESS;
	}

	if (!constr->intern) {
		constr->intern = Calloc(1, GMRFLib_constr__intern_tp);

	}

	if (constr->errcov_diagonal) {
		if (!constr->intern->chol) {
			constr->intern->chol = Calloc(nc, double);

			for (i = 0; i < nc; i++) {
				constr->intern->chol[i] = sqrt(constr->errcov_diagonal[i]);
			}
		}
		if (!constr->intern->Qpattern) {
			constr->intern->Qpattern = Calloc(nc, char);

			memset(constr->intern->Qpattern, 1, (unsigned) nc);	/* by definition they are.... */
		}
		if (!constr->intern->logdet) {
			constr->intern->logdet = Calloc(1, double);

			for (i = 0, ldet = 0.0; i < nc; i++) {
				ldet += log(constr->errcov_diagonal[i]);
			}
			*(constr->intern->logdet) = ldet;
		}
	} else {
		if (!constr->intern->chol) {
			GMRFLib_comp_chol_general(&(constr->intern->chol), constr->errcov_general, nc, NULL, GMRFLib_ESINGCONSTR);
		}
		if (!constr->intern->logdet) {
			constr->intern->logdet = Calloc(1, double);

			for (ldet = 0.0, i = 0; i < nc; i++) {
				ldet += log(constr->intern->chol[i + i * nc]);
			}
			*(constr->intern->logdet) = 2.0 * ldet;
		}
		if (!constr->intern->Qpattern) {
			double *Q = NULL, eps = GMRFLib_eps(1. / 2.);

			Q = Calloc(ISQR(nc), double);
			constr->intern->Qpattern = Calloc(ISQR(nc), char);

			memcpy(Q, constr->errcov_general, ISQR(nc) * sizeof(double));
			GMRFLib_comp_posdef_inverse(Q, nc);
			for (i = 0; i < ISQR(nc); i++) {
				constr->intern->Qpattern[i] = (ABS(Q[i]) < eps ? 0 : 1);
			}
			Free(Q);
		}
	}

	return GMRFLib_SUCCESS;
}

/*! \brief Evaluates the expressions 
  \f$ \mbox{\boldmath $Ax-e$} \f$ and 
  \f$ (\mbox{\boldmath $Ax-e$})^T\mbox{\boldmath $Q$}(\mbox{\boldmath $Ax-e$}) \f$ 
  for a given value of the GMRF <em>\b x</em>.

  \param[out] value A length <tt>constr->nc</tt> array holding the values of 
  <em>\b Ax - \b e</em>. If \c NULL, <em>\b Ax - \b e</em> is not evaluated.
  \param[out] sqr_value A length <tt>constr->nc</tt> array holding the value of 
  \f$ (\mbox{\boldmath $Ax-e$})^T\mbox{\boldmath $Q$}(\mbox{\boldmath $Ax-e$}) \f$
  If \c NULL,
  \f$ (\mbox{\boldmath $Ax-e$})^T\mbox{\boldmath $Q$}(\mbox{\boldmath $Ax-e$}) \f$
  is not evaluated.
  \param[in] x The instance of the GMRF, <em>\b x</em>, for which to compute 
  the constraints.
  \param[in] constr A pointer to a \c GMRFLib_constr_tp -object holding the 
  information on the type of linear constraint.
  \param[in] graph The graph on which the GMRF and the constraint are defined.
*/
int GMRFLib_eval_constr(double *value, double *sqr_value, double *x, GMRFLib_constr_tp * constr, GMRFLib_graph_tp * graph)
{
	/*
	 * eval a constraint `constr' at x-value `x':
	 * 
	 * if value, *value = Ax-e if sqr_value, *sqr_value = (Ax-e)'Q(Ax-e) 
	 */
	double *t_vector, alpha, beta;
	int inc = 1, nc;

	nc = constr->nc;
	t_vector = Calloc(nc, double);
	memcpy(t_vector, constr->e_vector, nc * sizeof(double));
	alpha = 1.0;
	beta = -1.0;
	dgemv_("N", &nc, &(graph->n), &alpha, constr->a_matrix, &nc, x, &inc, &beta, t_vector, &inc, 1);

	if (value) {
		memcpy(value, t_vector, nc * sizeof(double));
	}
	if (sqr_value) {
		*sqr_value = 0.0;
	}
	if (sqr_value && STOCHASTIC_CONSTR(constr)) {
		double *tt_vector, val;
		int i;

		tt_vector = Calloc(nc, double);

		if (constr->errcov_diagonal) {
			for (i = 0; i < nc; i++) {
				tt_vector[i] = t_vector[i] / constr->errcov_diagonal[i];
			}
		} else {
			if (constr->errcov_general) {
				GMRFLib_EWRAP0(GMRFLib_solveAxb_posdef(tt_vector, constr->intern->chol, t_vector, nc, 1));
			}
		}

		for (i = 0, val = 0.0; i < nc; i++) {
			val += t_vector[i] * tt_vector[i];
		}
		*sqr_value = val;

		Free(tt_vector);
	}
	Free(t_vector);

	return GMRFLib_SUCCESS;
}
int GMRFLib_duplicate_constr(GMRFLib_constr_tp ** new_constr, GMRFLib_constr_tp * constr, GMRFLib_graph_tp * graph)
{
	if (constr) {
		return GMRFLib_recomp_constr(new_constr, constr, NULL, NULL, NULL, graph, NULL);
	} else {
		if (new_constr) {
			*new_constr = NULL;
		}
		return GMRFLib_SUCCESS;
	}
}
int GMRFLib_recomp_constr(GMRFLib_constr_tp ** new_constr, GMRFLib_constr_tp * constr, double *x,
			  double *b_add, char *mask, GMRFLib_graph_tp * graph, GMRFLib_graph_tp * sub_graph)
{
	/*
	 * remap the constaints Ax=e, on graph to Ax=e on sub-graph
	 * 
	 * if stochastic constraints, then compute also ``b_add'', which is the additional terms to the 'b'-term: b^Tx.
	 * 
	 */

	GMRFLib_ENTER_ROUTINE;

	int i, ii, j, k, kk, n, ns, *in_use = NULL, *cmap = NULL, nc = 0;

	if (!constr) {
		if (new_constr) {
			*new_constr = NULL;
		}
		GMRFLib_LEAVE_ROUTINE;
		return GMRFLib_SUCCESS;
	}

	if (!mask) {
		/*
		 * new_constr is equal to the old one, make a copy 
		 */
		GMRFLib_make_empty_constr(new_constr);
		(*new_constr)->nc = constr->nc;

		(*new_constr)->a_matrix = Calloc(graph->n * constr->nc, double);
		memcpy((*new_constr)->a_matrix, constr->a_matrix, graph->n * constr->nc * sizeof(double));

		(*new_constr)->e_vector = Calloc(constr->nc, double);
		memcpy((*new_constr)->e_vector, constr->e_vector, constr->nc * sizeof(double));

		if (constr->errcov_diagonal) {
			(*new_constr)->errcov_diagonal = Calloc(constr->nc, double);
			memcpy((*new_constr)->errcov_diagonal, constr->errcov_diagonal, constr->nc * sizeof(double));
		} else {
			(*new_constr)->errcov_diagonal = NULL;
		}

		if (constr->errcov_general) {
			(*new_constr)->errcov_general = Calloc(ISQR(constr->nc), double);
			memcpy((*new_constr)->errcov_general, constr->errcov_general, ISQR(constr->nc) * sizeof(double));
		} else {
			(*new_constr)->errcov_general = NULL;
		}

		/*
		 * compute the internal stuff. scaling should not be needed.
		 */
		GMRFLib_prepare_constr(*new_constr, graph, 0);

		GMRFLib_LEAVE_ROUTINE;
		return GMRFLib_SUCCESS;
	}

	n = graph->n;
	ns = sub_graph->n;
	in_use = Calloc(constr->nc, int);
	cmap = Calloc(constr->nc, int);

	GMRFLib_make_empty_constr(new_constr);
	if (b_add) {
		memset(b_add, 0, sub_graph->n * sizeof(double));
	}


	/*
	 * find those constrs that are in use, part I 
	 */
	for (k = 0; k < constr->nc; k++) {
		for (i = 0; i < ns && in_use[k] == 0; i++) {
			if (constr->a_matrix[k + constr->nc * sub_graph->mothergraph_idx[i]]) {
				in_use[k] = 1;
			}
		}
	}

	/*
	 * find those constrs that are in use, part II.
	 * 
	 * if are in the errcov_general-case, then we need to do a detailed check: constraint k is out, iff Q_kj=0 for all
	 * other constraints j that are in (in_use[k]=1). (Q_kj can non-zero for a constraint j that is out (in_use[j]=0)). Q
	 * is the inverse of COV and its zero/non-zero pattern is stored in intern->Qpattern. 
	 */
	if (constr->errcov_general) {
		for (i = 0; i < constr->nc; i++) {
			if (!in_use[i]) {
				/*
				 * test if constraint i is really out 
				 */
				int really_out;

				for (j = 0, really_out = 1; j < constr->nc && really_out == 1; j++) {
					if (j != i) {
						if (in_use[j] == 1 && constr->intern->Qpattern[i + j * constr->nc]) {
							really_out = 0;
						}
					}
				}
				if (!really_out) {
					in_use[i] = 1;
				}
			}
		}
	}

	/*
	 * find those constrs that are in use, part III 
	 */
	for (k = nc = 0; k < constr->nc; k++) {
		if (in_use[k]) {
			cmap[nc++] = k;			       /* the mapping */
		}
	}
	(*new_constr)->nc = nc;

	if ((*new_constr)->nc == 0) {			       /* we do not have any constraints in sub_graph */
		Free(in_use);
		Free(cmap);
		Free(*new_constr);
		return GMRFLib_SUCCESS;
	}

	(*new_constr)->a_matrix = Calloc(nc * ns, double);
	(*new_constr)->e_vector = Calloc(nc, double);

	for (k = 0; k < nc; k++) {
		kk = cmap[k];
		(*new_constr)->e_vector[k] = constr->e_vector[kk];
		for (i = 0; i < n; i++) {
			if (mask[i]) {
				(*new_constr)->e_vector[k] -= constr->a_matrix[kk + constr->nc * i] * x[i];
			}
		}
		for (i = 0; i < ns; i++) {
			(*new_constr)->a_matrix[k + nc * i] = constr->a_matrix[kk + constr->nc * sub_graph->mothergraph_idx[i]];
		}
	}

	if (STOCHASTIC_CONSTR(constr)) {
		/*
		 * in this case we need to rebuild the covariance matrix, and compute 'b_add' which account also for setting
		 * 'b=0' in the constraint. 
		 */
		double *t_vec = NULL;

		t_vec = Calloc(nc, double);		       /* tmp-storage */

		/*
		 * compute new covariance-matrix for the new constraints 
		 */
		if (constr->errcov_diagonal) {
			(*new_constr)->errcov_diagonal = Calloc(nc, double);

			for (i = 0; i < nc; i++) {
				(*new_constr)->errcov_diagonal[i] = constr->errcov_diagonal[cmap[i]];
			}
		} else {				       /* the general case */
			(*new_constr)->errcov_general = Calloc(ISQR(nc), double);

			for (i = 0; i < nc; i++) {
				ii = cmap[i];
				for (j = 0; j < nc; j++) {
					(*new_constr)->errcov_general[i + j * nc] = constr->errcov_general[ii + cmap[j] * constr->nc];
				}
			}
		}

		/*
		 * compute the internal stuff & scale 
		 */
		GMRFLib_prepare_constr(*new_constr, sub_graph, 1);

		/*
		 * compute the contribution to b, ``b_add''.  t_vec = Q*(e-A_2 x_2). `b_add' account for `e_vector' in this
		 * case, so e_vector=0 
		 */
		if ((*new_constr)->errcov_diagonal) {
			for (i = 0; i < nc; i++) {
				t_vec[i] = (*new_constr)->e_vector[i] / (*new_constr)->errcov_diagonal[i];
			}
		} else {
			GMRFLib_solveAxb_posdef(t_vec, (*new_constr)->intern->chol, (*new_constr)->e_vector, nc, 1);
		}

		for (j = 0; j < nc; j++) {
			for (i = 0; i < ns; i++) {
				b_add[i] += (*new_constr)->a_matrix[i * nc + j] * t_vec[j];
			}
		}

		memset((*new_constr)->e_vector, 0, nc * sizeof(double));

		/*
		 * done ;-) 
		 */
		Free(t_vec);
	} else {
		/*
		 * scale it 
		 */
		GMRFLib_prepare_constr(*new_constr, sub_graph, 1);
	}

	Free(in_use);
	Free(cmap);
	GMRFLib_LEAVE_ROUTINE;

	return GMRFLib_SUCCESS;
}
int GMRFLib_info_problem(FILE * fp, GMRFLib_problem_tp * problem)
{
	if (problem) {
		return GMRFLib_fact_info_report(fp, &(problem->sub_sm_fact));
	}
	return GMRFLib_SUCCESS;
}
int GMRFLib_fact_info_report(FILE * fp, GMRFLib_sm_fact_tp * sm_fact)
{
	GMRFLib_fact_info_tp f;
	double possible_fillins;
	FILE *ffp = NULL;
	const char *sep = "-------------------------------------------------------------------------";

	if (!sm_fact)
		return GMRFLib_SUCCESS;

	ffp = (fp ? fp : stdout);
	f = sm_fact->finfo;

	fprintf(ffp, "\n\nGMRFLib report on factorisation\n%s\n", sep);

	fprintf(ffp, "Size ................................: %8d\n", f.n);
	if (f.n > 0.0) {
		fprintf(ffp, "Number of non-zeros in Q ............: %8d  (In percentage %.6f )\n", f.nnzero,
			(100.0 * f.nnzero / (SQR((double) f.n))));
		fprintf(ffp, "Number of non-zeros in L ............: %8d  (In percentage %.6f )\n",
			(f.nnzero - f.n) / 2 + f.n + f.nfillin,
			(100.0 * ((f.nnzero - f.n) / 2. + f.n + f.nfillin)) / (SQR((double) f.n) / 2. + (double) f.n / 2.));

		possible_fillins = SQR((double) f.n) / 2. + (double) f.n / 2. - f.nnzero;	/* terms in L - those in Q */

		if (possible_fillins > 0.0) {
			fprintf(ffp, "Number of fillins ...................: %8d  (In percentage %.6f )\n", f.nfillin,
				(100.0 * f.nfillin) / possible_fillins);
		}
		fprintf(ffp, "Terms in L / Minimum terms in L .....: %8.3f\n",
			((f.nnzero - f.n) / 2. + f.n + f.nfillin) / ((double) (f.nnzero - f.n) / 2. + f.n));
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_print_problem(FILE * fp, GMRFLib_problem_tp * problem)
{
	/*
	 * for verification purposes; just print the contents of 'problem' to 'fp' 
	 */

	int n, ns, nc;
	FILE *fpp = NULL;

	if (!problem) {
		return GMRFLib_SUCCESS;
	}

	n = problem->n;					       /* full graph */
	ns = problem->sub_graph->n;			       /* sub_graph */
	nc = (problem->sub_constr ? problem->sub_constr->nc : 0);
	fpp = (fp ? fp : stdout);

	fprintf(fpp, "\nContents of problem:\n");
	GMRFLib_print_darray(fpp, problem->sample, n, "sample");
	GMRFLib_print_darray(fpp, problem->mean, n, "mean");
	GMRFLib_print_darray(fpp, problem->mean_constr, n, "mean_constr");
	GMRFLib_print_iarray(fpp, problem->map, ns, "map");
	GMRFLib_print_darray(fpp, &(problem->sub_logdens), 1, "sub_logdens");
	GMRFLib_print_darray(fpp, problem->sub_sample, ns, "sub_sample");
	GMRFLib_print_darray(fpp, problem->sub_mean, ns, "sub_mean");
	GMRFLib_print_darray(fpp, problem->sub_mean_constr, ns, "sub_mean_constr");
	GMRFLib_print_darray(fpp, problem->constr_m, ns * nc, "constr_m");
	GMRFLib_print_darray(fpp, problem->l_aqat_m, nc * nc, "l_aqat_m");
	GMRFLib_print_darray(fpp, problem->qi_at_m, ns * nc, "qi_at_m");
	GMRFLib_print_darray(fpp, &(problem->logdet_aqat), 1, "logdet_aqat");
	GMRFLib_print_darray(fpp, &(problem->log_normc), 1, "log_normc");
	GMRFLib_print_darray(fpp, &(problem->exp_corr), 1, "exp_corr");

	GMRFLib_print_constr(fpp, problem->sub_constr, problem->sub_graph);

	fflush(fpp);
	return GMRFLib_SUCCESS;
}

GMRFLib_problem_tp *GMRFLib_duplicate_problem(GMRFLib_problem_tp * problem, int skeleton, int copy_ptr, int copy_pardiso_ptr)
{
	/*
	 * duplicate a problem 
	 */

#define DUPLICATE(name, len, tp, skeleton_)				\
	if (1) {							\
		if (problem->name && ((len)>0) && !skeleton_){		\
			np->name = Calloc((len), tp);			\
			memcpy(np->name, problem->name, (len)*sizeof(tp)); \
		} else {						\
			np->name = NULL;				\
		}							\
	}
#define COPY(name) np->name = problem->name

	if (!problem) {
		return NULL;
	}

	GMRFLib_problem_tp *np = Calloc(1, GMRFLib_problem_tp);

	int i;
	int n = problem->n;				       /* full graph */
	int ns = problem->sub_graph->n;			       /* sub_graph */
	int nc = (problem->sub_constr ? problem->sub_constr->nc : 0);

	DUPLICATE(sample, n, double, skeleton);
	DUPLICATE(mean, n, double, skeleton);
	DUPLICATE(mean_constr, n, double, skeleton);

	COPY(n);
	COPY(sub_logdens);
	DUPLICATE(sub_sample, ns, double, skeleton);
	DUPLICATE(sub_mean, ns, double, skeleton);
	DUPLICATE(sub_mean_constr, ns, double, skeleton);

	/*
	 * duplicate the sparse-matrix factorisation 
	 */
	DUPLICATE(sub_sm_fact.remap, ns, int, 0);
	DUPLICATE(sub_sm_fact.bchol, ns * (problem->sub_sm_fact.bandwidth + 1), double, 0);

	COPY(sub_sm_fact.bandwidth);
	COPY(sub_sm_fact.smtp);

	if (problem->sub_sm_fact.TAUCS_L && !skeleton) {
		np->sub_sm_fact.TAUCS_L = GMRFLib_L_duplicate_TAUCS(problem->sub_sm_fact.TAUCS_L, problem->sub_sm_fact.TAUCS_L->flags);
	} else {
		np->sub_sm_fact.TAUCS_L = NULL;
	}

	if (problem->sub_sm_fact.TAUCS_L_inv_diag && !skeleton) {
		DUPLICATE(sub_sm_fact.TAUCS_L_inv_diag, ns, double, skeleton);
	} else {
		np->sub_sm_fact.TAUCS_L_inv_diag = NULL;
	}
	np->sub_sm_fact.TAUCS_symb_fact = GMRFLib_sm_fact_duplicate_TAUCS(problem->sub_sm_fact.TAUCS_symb_fact);
	COPY(sub_sm_fact.finfo);

	if (problem->sub_sm_fact.PARDISO_fact) {
		GMRFLib_duplicate_pardiso_store(&(np->sub_sm_fact.PARDISO_fact), problem->sub_sm_fact.PARDISO_fact, copy_ptr, copy_pardiso_ptr);
	}

	/*
	 * then the constraint 
	 */
	if (skeleton) {
		np->sub_constr = NULL;
	} else {
		if (problem->sub_constr) {
			GMRFLib_make_empty_constr(&(np->sub_constr));
			COPY(sub_constr->nc);
			DUPLICATE(sub_constr->a_matrix, nc * ns, double, 0);
			DUPLICATE(sub_constr->e_vector, nc, double, 0);

			if (STOCHASTIC_CONSTR(problem->sub_constr)) {
				DUPLICATE(sub_constr->errcov_diagonal, nc, double, 0);
				DUPLICATE(sub_constr->errcov_general, ISQR(nc), double, 0);

				if (problem->sub_constr->intern) {
					np->sub_constr->intern = Calloc(1, GMRFLib_constr__intern_tp);
					DUPLICATE(sub_constr->intern->Qpattern, ISQR(nc), char, 0);
				}
			}
		} else {
			COPY(sub_constr);
		}
	}

	DUPLICATE(sub_constr_value, nc, double, skeleton);
	DUPLICATE(constr_m, ns * nc, double, skeleton);
	DUPLICATE(l_aqat_m, nc * nc, double, skeleton);
	DUPLICATE(inv_aqat_m, nc * nc, double, skeleton);
	DUPLICATE(qi_at_m, ns * nc, double, skeleton);

	COPY(logdet_aat);
	COPY(logdet_aqat);
	COPY(log_normc);
	COPY(exp_corr);
	GMRFLib_copy_graph(&(np->sub_graph), problem->sub_graph);
	np->map = np->sub_graph->mothergraph_idx;	       /* its kind of special */

	/*
	 * copy the tab 
	 */
	if (problem->tab && !skeleton) {
		GMRFLib_tabulate_Qfunc_tp *tab = Calloc(1, GMRFLib_tabulate_Qfunc_tp);
		GMRFLib_tabulate_Qfunc_arg_tp *Qfunc_arg = Calloc(1, GMRFLib_tabulate_Qfunc_arg_tp);
		GMRFLib_tabulate_Qfunc_arg_tp *tmp = (GMRFLib_tabulate_Qfunc_arg_tp *) (problem->tab->Qfunc_arg);

		tab->Qfunc = problem->tab->Qfunc;
		Qfunc_arg->n = tmp->n;
		Qfunc_arg->prec = tmp->prec;
		Qfunc_arg->log_prec = tmp->log_prec;
		if (tmp->log_prec_omp) {
			int tmax = GMRFLib_MAX_THREADS;
			Qfunc_arg->log_prec_omp = Calloc(tmax, double *);
			for (i = 0; i < tmax; i++) {
				Qfunc_arg->log_prec_omp[i] = tmp->log_prec_omp[i];
			}
		}
		Qfunc_arg->values = Calloc(ns, map_id *);
		for (i = 0; i < ns; i++) {
			Qfunc_arg->values[i] = GMRFLib_duplicate_map_id(tmp->values[i]);
		}
		tab->Qfunc_arg = (void *) Qfunc_arg;
		np->tab = tab;
	} else {
		np->tab = NULL;
	}

	/*
	 * copy the sub_inverse 
	 */
	if (problem->sub_inverse && !skeleton) {
		np->sub_inverse = Calloc(1, GMRFLib_Qinv_tp);
		map_id **Qinv = Calloc(n, map_id *);

		for (i = 0; i < n; i++) {
			Qinv[i] = GMRFLib_duplicate_map_id(problem->sub_inverse->Qinv[i]);
		}
		np->sub_inverse->Qinv = Qinv;
		np->sub_inverse->mapping = GMRFLib_duplicate_map_ii(problem->sub_inverse->mapping);
	} else {
		np->sub_inverse = NULL;
	}

#undef DUPLICATE
#undef COPY
	return np;
}

GMRFLib_sizeof_tp GMRFLib_sizeof_problem(GMRFLib_problem_tp * problem)
{
#define DUPLICATE(name, len, tp)  if (problem->name && ((len)>0)) siz += (len)*sizeof(tp);

	GMRFLib_sizeof_tp siz = 0;

	if (!problem) {
		return 0;
	}

	siz += sizeof(GMRFLib_problem_tp);

	int i;
	int n = problem->n;				       /* full graph */
	int ns = problem->sub_graph->n;			       /* sub_graph */
	int nc = (problem->sub_constr ? problem->sub_constr->nc : 0);

	DUPLICATE(sample, n, double);
	DUPLICATE(mean, n, double);
	DUPLICATE(mean_constr, n, double);

	DUPLICATE(sub_sample, ns, double);
	DUPLICATE(sub_mean, ns, double);
	DUPLICATE(sub_mean_constr, ns, double);

	/*
	 * duplicate the sparse-matrix factorisation 
	 */
	DUPLICATE(sub_sm_fact.remap, ns, int);
	DUPLICATE(sub_sm_fact.bchol, ns * (problem->sub_sm_fact.bandwidth + 1), double);

	if (problem->sub_sm_fact.TAUCS_L) {
		siz += GMRFLib_L_sizeof_TAUCS(problem->sub_sm_fact.TAUCS_L);
	}

	if (problem->sub_sm_fact.TAUCS_L_inv_diag) {
		DUPLICATE(sub_sm_fact.TAUCS_L_inv_diag, ns, double);
	}
	siz += GMRFLib_sm_fact_sizeof_TAUCS(problem->sub_sm_fact.TAUCS_symb_fact);

	/*
	 * then the constraint 
	 */
	if (problem->sub_constr) {
		siz += sizeof(GMRFLib_constr_tp);
		DUPLICATE(sub_constr->a_matrix, nc * ns, double);
		DUPLICATE(sub_constr->e_vector, nc, double);

		if (STOCHASTIC_CONSTR(problem->sub_constr)) {
			DUPLICATE(sub_constr->errcov_diagonal, nc, double);
			DUPLICATE(sub_constr->errcov_general, ISQR(nc), double);

			if (problem->sub_constr->intern) {
				DUPLICATE(sub_constr->intern->Qpattern, ISQR(nc), char);
			}
		}
	}

	DUPLICATE(sub_constr_value, nc, double);
	DUPLICATE(constr_m, ns * nc, double);
	DUPLICATE(l_aqat_m, nc * nc, double);
	DUPLICATE(inv_aqat_m, nc * nc, double);
	DUPLICATE(qi_at_m, ns * nc, double);

	siz += GMRFLib_sizeof_graph(problem->sub_graph);

	/*
	 * copy the tab 
	 */
	if (problem->tab) {
		GMRFLib_tabulate_Qfunc_arg_tp *tmp = (GMRFLib_tabulate_Qfunc_arg_tp *) (problem->tab->Qfunc_arg);

		siz += ns * sizeof(map_id *);
		for (i = 0; i < ns; i++) {
			siz += GMRFLib_sizeof_map_id(tmp->values[i]);
		}
	}

	/*
	 * copy the sub_inverse 
	 */
	if (problem->sub_inverse) {
		siz += sizeof(GMRFLib_Qinv_tp);
		for (i = 0; i < n; i++) {
			siz += GMRFLib_sizeof_map_id(problem->sub_inverse->Qinv[i]);
		}
		siz += GMRFLib_sizeof_map_ii(problem->sub_inverse->mapping);
	}
#undef DUPLICATE

	return siz;
}
GMRFLib_sizeof_tp GMRFLib_sizeof_store(GMRFLib_store_tp * store)
{
	/*
	 * return, approximately, the size of STORE 
	 */

	if (!store)
		return 0;

	GMRFLib_sizeof_tp siz = 0;
	int ns = store->sub_graph->n;

	siz += sizeof(GMRFLib_store_tp);
	siz += sizeof(double);
	siz += ns * sizeof(int);
	siz += GMRFLib_sizeof_graph(store->sub_graph);
	siz += GMRFLib_sm_fact_sizeof_TAUCS(store->TAUCS_symb_fact);
	siz += 5 * sizeof(double);
	siz += GMRFLib_sizeof_problem(store->problem_old2new);
	siz += GMRFLib_sizeof_problem(store->problem_new2old);

	if (store->diag_store) {
		siz += ns * sizeof(double);
	}
	if (store->sub_store) {
		siz += GMRFLib_sizeof_store(store->sub_store);
	}

	return siz;
}

GMRFLib_store_tp *GMRFLib_duplicate_store(GMRFLib_store_tp * store, int skeleton, int copy_ptr, int copy_pardiso_ptr)
{
	/*
	 * duplicate STORE 
	 */
#define DUPLICATE(name, len, tp, skeleton_)				\
	if (1) {							\
		if (store->name && ((len)>0) && !skeleton_){		\
			new_store->name = Calloc((len), tp);		\
			memcpy(new_store->name, store->name, (len)*sizeof(tp)); \
		} else {						\
			new_store->name = NULL;				\
		}							\
	}
#define COPY(name) new_store->name = store->name

	if (!store) {
		return NULL;
	}

	int ns = store->sub_graph->n;
	GMRFLib_store_tp *new_store = Calloc(1, GMRFLib_store_tp);

	COPY(bandwidth);
	DUPLICATE(remap, ns, int, 0);

	int id = omp_get_thread_num();
	GMRFLib_meminfo_thread_id = id;

	if (copy_ptr == GMRFLib_TRUE) {
		/*
		 * just copy ptr's; read only 
		 */
		new_store->sub_graph = store->sub_graph;
		new_store->TAUCS_symb_fact = store->TAUCS_symb_fact;
	} else {
		GMRFLib_copy_graph(&(new_store->sub_graph), store->sub_graph);
		new_store->TAUCS_symb_fact = GMRFLib_sm_fact_duplicate_TAUCS(store->TAUCS_symb_fact);
	}
	new_store->copy_ptr = copy_ptr;
	new_store->copy_pardiso_ptr = copy_pardiso_ptr;
	if (store->PARDISO_fact) {
		GMRFLib_duplicate_pardiso_store(&(new_store->PARDISO_fact), store->PARDISO_fact, copy_ptr, copy_pardiso_ptr);
	}

	GMRFLib_meminfo_thread_id *= -1;
	char *tmp = Calloc(1, char);
	Free(tmp);


	COPY(store_problems);
	COPY(fixed_hyperparameters);
	COPY(decision);
	DUPLICATE(old_logdens, 1, double, skeleton);
	DUPLICATE(new_logdens, 1, double, skeleton);

	if (!skeleton) {
		new_store->problem_old2new = GMRFLib_duplicate_problem(store->problem_old2new, skeleton, copy_ptr, copy_pardiso_ptr);
		new_store->problem_new2old = GMRFLib_duplicate_problem(store->problem_new2old, skeleton, copy_ptr, copy_pardiso_ptr);
	} else {
		new_store->problem_new2old = NULL;
		new_store->problem_old2new = NULL;
	}

	if (store->diag_store) {
		new_store->diag_store = GMRFLib_duplicate_store(store->diag_store, skeleton, copy_ptr, copy_pardiso_ptr);
	}
	if (store->sub_store) {
		new_store->sub_store = GMRFLib_duplicate_store(store->sub_store, skeleton, copy_ptr, copy_pardiso_ptr);
	}
#undef DUPLICATE
#undef COPY
	return new_store;
}
double GMRFLib_Qfunc_generic(int i, int j, void *arg)
{
	if (i != j) {
		return -1.0;
	} else {
		GMRFLib_graph_tp *g = (GMRFLib_graph_tp *) arg;
		return g->n;
	}
}

/*!
  \brief Set \c GMRFLib_reorder to the best reordering found for a given graph

  This functions factorise a precision matrix (symbolically) using several different reorderings techniques and chose the one with fewest fillins. If sizeof_L is
  non-NULL, then the sizeof_L in bytes, is returned.
*/
int GMRFLib_optimize_reorder(GMRFLib_graph_tp * graph, GMRFLib_sizeof_tp * nnz_opt, int *use_global, GMRFLib_global_node_tp * gn)
{
	if (!graph) {
		if (nnz_opt)
			*nnz_opt = 0;
		return GMRFLib_SUCCESS;
	}

	if (GMRFLib_smtp == GMRFLib_SMTP_BAND) {
		GMRFLib_reorder = GMRFLib_REORDER_DEFAULT;
		*nnz_opt = 0;
	} else if (GMRFLib_smtp == GMRFLib_SMTP_PARDISO) {
		GMRFLib_reorder = GMRFLib_REORDER_PARDISO;
		*nnz_opt = 0;
	} else {
		GMRFLib_sizeof_tp *nnzs = NULL, nnz_best;
		int k, debug = 0, n = -1, nk, r, i, ne = 0, use_global_nodes;
		GMRFLib_reorder_tp rs[] = { GMRFLib_REORDER_METIS, GMRFLib_REORDER_AMDC };
		taucs_ccs_matrix *Q = NULL;
		char *fixed = NULL;
		double *cputime = NULL;

		n = graph->n;

		/*
		 * build the Q-matrix; just symbolically, so I set Q_ij = 1. this matrix is common for all reordering
		 * schemes. 
		 */
		int ic, kk, j, nnz;

		for (i = 0, nnz = n; i < n; i++) {
			nnz += graph->nnbs[i];
		}
		Q = taucs_ccs_create(n, n, nnz, TAUCS_DOUBLE);
		GMRFLib_ASSERT(Q, GMRFLib_EMEMORY);
		Q->flags = (TAUCS_DOUBLE | TAUCS_SYMMETRIC | TAUCS_TRIANGULAR | TAUCS_LOWER | TAUCS_PATTERN);
		Q->colptr[0] = 0;

		for (i = 0, ic = 0; i < n; i++) {
			Q->rowind[ic] = i;
			Q->values.d[ic] = 1.0;
			ic++;
			ne = 1;

			for (kk = 0; kk < graph->nnbs[i]; kk++) {
				j = graph->nbs[i][kk];
				if (j > i) {
					break;
				}
				Q->rowind[ic] = j;
				Q->values.d[ic] = 1.0;
				ic++;
				ne++;
			}
			Q->colptr[i + 1] = Q->colptr[i] + ne;
		}

		/*
		 * do this for all reorderings, first with global node detection and then without global nodes detection 
		 */
		nk = 2 * (int) (sizeof(rs) / sizeof(int));     /* yes, twice... */
		nnzs = Calloc(nk, GMRFLib_sizeof_tp);
		cputime = Calloc(nk, double);

//#pragma omp parallel for private(k) schedule(dynamic)
		for (k = 0; k < nk; k++) {
			int *iperm = NULL, *perm = NULL, ii, kkk, use_global_nodes;
			supernodal_factor_matrix *TAUCS_symb_fact = NULL;
			taucs_ccs_matrix *L = NULL;

			GMRFLib_global_node_tp lgn;
			if (gn) {
				memcpy((void *) &lgn, (void *) gn, sizeof(GMRFLib_global_node_tp));
			} else {
				memcpy((void *) &lgn, (void *) &GMRFLib_global_node, sizeof(GMRFLib_global_node_tp));
			}

			/*
			 * first half is with, second half is without, the detection of global nodes 
			 */
			kkk = k % (nk / 2);
			use_global_nodes = (k < nk / 2);

			if (!use_global_nodes) {
				/*
				 * currently disable 
				 */
				lgn.factor = 2.0;
				lgn.degree = INT_MAX;
			}

			/*
			 * if we run with global nodes, it can be that the settings are so that they are not in effect. we check for this here, and if so, we do not
			 * need to try it, as it will be tried again in the second half where the global_nodes are disabled.
			 */
			if (!use_global_nodes || !(use_global_nodes && (lgn.factor > 1.0) && (lgn.degree > graph->n - 1))) {

				cputime[k] = GMRFLib_cpu();
				GMRFLib_compute_reordering_TAUCS(&iperm, graph, rs[kkk], &lgn);

				perm = Calloc(n, int);
				for (ii = 0; ii < n; ii++) {
					perm[iperm[ii]] = ii;
				}

				L = taucs_ccs_permute_symmetrically(Q, perm, iperm);	/* permute the matrix */
				TAUCS_symb_fact = (supernodal_factor_matrix *) taucs_ccs_factor_llt_symbolic(L);
				nnzs[k] = GMRFLib_sm_fact_nnz_TAUCS(TAUCS_symb_fact);
				Free(perm);
				Free(iperm);
				taucs_ccs_free(L);
				taucs_supernodal_factor_free(TAUCS_symb_fact);

				cputime[k] = GMRFLib_cpu() - cputime[k];

				if (debug) {
#pragma omp critical
					{
						printf("%s: reorder=[%s] \tnnz=%lu \tUseGlobalNodes=%1d cpu=%.4f\n",
						       __GMRFLib_FuncName, GMRFLib_reorder_name(rs[kkk]), nnzs[k], use_global_nodes, cputime[k]);
					}
				}

			} else {
				nnzs[k] = UINT_MAX;
			}
		}

		/*
		 * find the best one 
		 */
		r = 0;
		nnz_best = nnzs[r];
		for (k = 1; k < nk; k++) {
			if (nnzs[k] < nnz_best) {
				r = k;
				nnz_best = nnzs[r];
			}
		}

		/*
		 * find out which one this is 
		 */
		use_global_nodes = (r < (nk / 2));
		r = r % (nk / 2);

		if (!use_global_nodes) {
			GMRFLib_global_node_tp g = { 2.0, INT_MAX };
			GMRFLib_global_node = g;
			if (gn) {
				/*
				 * if so, copy the result over there as well
				 */
				memcpy((void *) gn, (void *) &g, sizeof(GMRFLib_global_node));
			}
		} else {
			if (gn) {
				memcpy((void *) &GMRFLib_global_node, (void *) gn, sizeof(GMRFLib_global_node));
			}
		}
		GMRFLib_reorder = rs[r];

		if (debug) {
			printf("%s: best reordering=[%s] UseGlobalNodes=%1d\n", __GMRFLib_FuncName,
			       GMRFLib_reorder_name(GMRFLib_reorder), use_global_nodes);
		}
		taucs_ccs_free(Q);
		if (nnz_opt)
			*nnz_opt = nnz_best;
		if (use_global)
			*use_global = use_global_nodes;

		Free(nnzs);
		Free(fixed);
	}
	return GMRFLib_SUCCESS;
}


/*



  Example for manual
 */

/*! \page ex_problem-setup Sampling from a GMRF

  In the example below, unconditional and conditional sampling from a
  GMRF <em>\b x</em> are illustrated. Th unconditional density of the GMRF is
  given by the general expression in <b>(GMRF-2)</b>, re-stated
  below for convenience:
  \f[\pi(\mbox{\boldmath $x$}) \propto \exp\left( -\frac{1}{2} 
  (\mbox{\boldmath $x$}-\mbox{\boldmath $\mu$})^T
  (\mbox{\boldmath $Q$} + \mbox{diag}(\mbox{\boldmath $c$}))
  (\mbox{\boldmath $x$}-\mbox{\boldmath $\mu$}) + 
  \mbox{\boldmath $b$}^T\mbox{\boldmath $x$}
  \right). \f]

  For all sampling methods illustrated in the example, we define the
  graph to be a lattice graph on a 6*6 lattice, using a 3*3 neighbourhood. 
  The <em>\b Q</em>-matrix is defined by
  \f[ \mbox{\boldmath $Q$} = \kappa \left( \begin{array}{cccccc}
  nnbs[1] & -1 & 0 & 0 & \cdots & 0 \\ 
  -1 & nnbs[2] & -1 & 0 & \cdots & 0 \\
  0 & -1 & nnbs[3] & -1 & \cdots & 0 \\
  \vdots & \vdots & \vdots & \vdots & \ddots & \vdots \\
  0 & 0 & 0 & 0 & \cdots &  nnbs[6]\\
  \end{array} \right), \f]

  where <em>nnbs[i]</em> is the number of neighbours of node \em i.
  This matrix is singular, but the singularity is resolved either by
  adding elements to the diagonal, or conditioning on fixed values.
  In the example, we use \f$ \kappa=1 \f$.

  Four sampling methods are implemented in the example. These are
  - <em>Unconditional sampling</em>. In this case, we generate samples
    from <b>(GMRF-2)</b>.  The remaining parameters are set to
    \f$ \mbox{\boldmath $b$}=\mbox{\boldmath $0$}, 
    \mbox{\boldmath $\mu$} = \mbox{\boldmath $0$} \mbox{ and }
    \mbox{\boldmath $c$} = \mbox{\boldmath $1$}_n \f$. The
    initial values of the GMRF (that are not used in this case) are
    set to <em>\b x = \b 0</em>.
  - <em>Conditioning on fixed values</em>. We sample from
    \f$ \mbox{\boldmath $x$}_{-\mathcal{F}}|
    \mbox{\boldmath $x$}_{\mathcal{F}} \f$, where \f$ \mathcal{F} \f$
    is the set of indices \em i for which \em x_i is fixed. In the example
    we let 
    \f$ \mathcal{F} = \{i:i\mbox{2} = 0\}, \mbox{ and } x[i] =
    \mbox{Unif}(0,1);\; i=0,\ldots,n-1 \f$.  
    The remaining parameters are assigned the values 
    \f$ \mbox{\boldmath $b$}=\mbox{\boldmath $0$}, 
    \mbox{\boldmath $\mu$} = \mbox{\boldmath $0$} \mbox{ and }
    \mbox{\boldmath $c$} = \mbox{\boldmath $0$} \f$.
  - <em>Conditioning on a deterministic linear constraint</em>,
    sampling from 
    \f$ \pi(\mbox{\boldmath $x$}|\mbox{\boldmath $Ax$} = 
    \mbox{\boldmath $e$}) \f$.  
    In the example, we use the constraints \f$ \sum_i x_i = 0 \mbox{ and }
    x_0 + 2 x_1 = 1 \f$, which correspond to
    \f[ \mbox{\boldmath $A$} = \left( \begin{array}{ccccc}
    1 & 1 & 1 & \cdots & 1 \\ 1 & 2 & 0 & \cdots & 0 \\ 
    \end{array} \right), \mbox{ and  } 
    \mbox{\boldmath $e$} = \left( \begin{array}{c} 0 \\ 1 \\ 
    \end{array} \right). \f]
    Also, we condition on fixed values, fixing the same elements of
    \f$ \mbox{\boldmath $x$} \f$ as for method 2. The values of the 
    GMRF are 
    \f$ \mbox{\boldmath $x$} = \mbox{Unif}(0,1);\; i=0,\ldots,n-1 \f$. 
    The remaining parameters are set to 
    \f$ \mbox{\boldmath $b$}=\mbox{Unif}(0,1), 
    \mbox{\boldmath $\mu$} = \mbox{\boldmath $0$} \mbox{ and }
    \mbox{\boldmath $c$} = \mbox{\boldmath $1$}_n \f$.
  - <em>Conditioning on a stochastic linear constraint</em>, sampling
    from 
    \f$ \pi(\mbox{\boldmath $x$}|\mbox{\boldmath $Ax$} = 
    \mbox{\boldmath $e$} + \mbox{\boldmath $\epsilon$}) \f$. 
    This is equivalent to sampling unconditionally from
    <b>(GMRF-11)</b> in \ref sampling. We use the constraint 
    \f$ \sum_i x_i = 0 \f$ and assume that 
    \f$ \epsilon \sim N(0, \sigma^2) \f$, where \f$ \sigma^2 = 1 \f$. 
    This corresponds to
    \f[ \mbox{\boldmath $A$} = \left( \begin{array}{cccc}
    1 & 1  & \cdots & 1 \\ \end{array} \right), \quad
    \mbox{\boldmath $e$} = e = 0 \mbox{  and  } 
    \mbox{\boldmath $\Sigma$}_{\epsilon} =  \sigma^2. \f]
    The other parameters are assign values equal to the values used
    for method 3, but we do not condition on fixed values in this
    case.

  \par Program code:

  \verbinclude example-doxygen-sample.txt

  \par Output: for <tt> method = 3 </tt>

  \verbinclude doxygen_problem_2.txt

*/
