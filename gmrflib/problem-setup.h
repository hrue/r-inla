
/* problem-setup.h
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
  \file problem-setup.h
  \brief Typedefs and defines for \ref problem-setup.c
*/
#ifndef __GMRFLib_PROBLEM_SETUP_H__
#define __GMRFLib_PROBLEM_SETUP_H__

#if !defined(__FreeBSD__)
#include <strings.h>
#include <malloc.h>
#endif
#include <stdlib.h>
#include <math.h>

#include "GMRFLib/hashP.h"

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

/*!
  \brief Define a \c new problem. No previous information is avilable or used
*/
#define GMRFLib_NEW_PROBLEM   0x0000

/*!
  \brief Keep the \c mean
*/
#define GMRFLib_KEEP_mean     0x0001

/*!
  \brief Keep the \c graph (or \c subgraph)
*/
#define GMRFLib_KEEP_graph    0x0002

/*!
  \brief Keep the Cholesky factorisation
*/
#define GMRFLib_KEEP_chol     0x0004

/*!
  \brief Keep the constraint
*/
#define GMRFLib_KEEP_constr   0x0008

/*
  Update the constraint with one extra. WARNING: THIS IS FOR INTERNAL USE ONLY AND NOT DOCUMENTED, as very very special
  conditions apply...
*/
#define GMRFLib_UPDATE_constr 0x0010

/*! 
  \struct GMRFLib_constr__intern_tp problem-setup.h

  \brief Mainly intended for internal usage, holding the Cholesky decomposition (\c chol}) and the
  log of the determinant (\c logdet}) of the covariance matrix \f$ \mbox{\boldmath
  $\Sigma$}_{\epsilon}\f$
*/

/*
 */
    typedef struct {

	/**
	 *  \brief = chol(errcov_*) 
	 */
	double *chol;

	/**
	 *  \brief = log(det(errcov_*)) 
	 */
	double *logdet;

	/**
	 *  \brief = the zero/non-zero pattern of inv(errcov_*) 
	 */
	char *Qpattern;
} GMRFLib_constr__intern_tp;

/*! 
  \struct GMRFLib_constr_tp problem-setup.h

  \brief Defining the constraint, of the form \f$ \mbox{\boldmath $Ax$} = \mbox{\boldmath
  $e$}+\mbox{\boldmath $\epsilon$} \f$ for conditional simulation on a linear deterministic or
  stochastic constraint.

  In the case of a deterministic constraint, the last three members of \c GMRFLib_constr_tp are not
  needed, and the members \a errcov_diagonal and \a errcov_general should be set to \c NULL. If one
  or both of these members are \c !NULL, the constraint is assumed to be stochastic, and the last
  three members specify the covariance matrix \f$ \mbox{\boldmath $\Sigma$}_{\epsilon} \f$ of \f$
  \mbox{\boldmath $\epsilon$} \f$ and it's decomposition. The matrix \f$ \mbox{\boldmath
  $\Sigma$}_{\epsilon} \f$ is required to be non-singular.

  To simplify allocation, initialization and deallocation of a \c GMRFLib_constr_tp -object, a set
  of functions operating on this data structure is implemented. To ensure safe memory handling, the
  user is recommended to use these function, which are documented in the next subsection.  A \c
  GMRFLib_constr_tp -object is created and initialized by the following steps:
  - Allocate memory by calling \c GMRFLib_make_empty_constr().
  - Set the number of constraints \em nc, and explicitly allocate the arrays 
  \a a_matrix and \a e_vector by using \c calloc, and fill in the values
  of their elements. The number of elements of the GMRF (the number
  of columns of <em>\b A</em>), will be defined within the corresponding
  \c GMRFLib_graph_tp -object, for which the constraints are defined.
  - Compute internal information, customizing the object for use in
  the sampling routines, by calling \c GMRFLib_prepare_constr(). If the constraint 
  is deterministic, nothing is to be done, and the function is empty.

  Two other operations on constraints are available:
  - To print the contents of the \c GMRFLib_constr_tp -object, call 
  \c GMRFLib_print_constr().
  - To evaluate the expressions 
  \f$ \mbox{\boldmath $Ax-e$} \f$ and \f$ (\mbox{\boldmath $Ax-e$})^T\mbox{\boldmath
  $Q$}(\mbox{\boldmath $Ax-e$}) \f$ for a given value of the GMRF <em>\b x</em>, call \c
  GMRFLib_eval_constr().
*/
typedef struct {

	/**
	 *  \brief Number of constaints, can be 0. 
	 */
	int nc;

	/**
	 *  \brief <em>\b A</em> matrix
	 * 
	 * A length <tt>(nc * n)</tt> array holding the elements of the <tt>(nc * n)</tt> matrix <em>\b A</em> in the linear
	 * constraint <em>\b Ax = \b e</em>. The matrix is assumed to be stored column by column, that is, element <tt>i + nc * 
	 * j</tt> of \a a_matrix should correspond to <em>\b A</em>[i,j]. \n\n 
	 */
	double *a_matrix;

	/**
	 *  \brief <em>\b e</em> vector
	 * 
	 * A length \c nc array holding the elements of the vector <em>\b e</em> of the constraint <em>\b Ax = \b e</em>. \n\n 
	 */
	double *e_vector;

	/**
	 *  \brief If non-NULL, a diagonal non-singular error cov-matrix
	 * 
	 * If <tt>!=NULL</tt>, the covariance matrix \f$ \mbox{\boldmath $\Sigma$}_{\epsilon} \f$ of the error \f$
	 * \mbox{\boldmath $\epsilon$} \f$ is assumed to be diagonal, and \a errcov_diagonal is a length \c nc -array
	 * containing the diagonal elements.\n\n 
	 */
	double *errcov_diagonal;

	/**
	 *  \brief If non-NULL, a general positive definite error cov-matrix
	 * 
	 * If <tt>!=NULL</tt>, \a errcov_general specifies a general positive definite covariance matrix \f$ \mbox{\boldmath
	 * $\Sigma$}_{\epsilon} \f$ of the error \f$ \mbox{\boldmath $\epsilon$} \f$. If both \a errcov_general and \a
	 * errcov_diagonal are <tt>!NULL</tt>, the matrix is assumed to be diagonal, and \a errcov_diagonal is used. \n\n 
	 */
	double *errcov_general;

	/**
	 *  \brief For internal use only 
	 */
	GMRFLib_constr__intern_tp *intern;
} GMRFLib_constr_tp;

typedef struct {
	void *user_Qfunc_args;
	GMRFLib_Qfunc_tp *user_Qfunc;
	int *map;					       /* Mapping to the real-world */
	double *diagonal_adds;
} GMRFLib_Qfunc_arg_tp;

/* 
   options for Qinv
*/

/*!
  \brief Store all values computed
*/
#define GMRFLib_QINV_ALL         0x0000

/*!
  \brief Store only marginal values and covariances for neigbours
*/
#define GMRFLib_QINV_NEIGB       0x0001

/*!
  \brief Store only marginal variances
*/
#define GMRFLib_QINV_DIAG        0x0002

/*!
  \brief Disable check for ``complete'' L.
*/
#define GMRFLib_QINV_NO_CHECK    0x0004

/*!
  \brief Check once only for a ``complete'' L
*/
#define GMRFLib_QINV_CHECK_ONCE  0x0008

typedef struct {

	/**
	 *  \brief The hash-table holding the inverse of \a Q.
	 * 
	 * The inverse of \a Q, \$S=(S_{ij})\$, are computed for those \a i and \a j where \$L_{ij}\$ is non-zero. The values
	 * of \$S_{ij}\$ are store in internal (reordered) coordinates.
	 * 
	 */
	map_id **Qinv;

	/**
	 *  \brief The mapping used to lookup values in \a Qinv 
	 */
	map_ii *mapping;
} GMRFLib_Qinv_tp;

/*! 
  \struct GMRFLib_problem_tp problem-setup.h
  \brief Specification of the sampling problem.

  The simulation routines operate on objects of type \c GMRFLib_problem_tp. 
  The members of the \em struct are all initialized and operated on within the
  library routines. Many of the members are internal, and only the members 
  - \a sample
  - \a mean
  - \a mean_constr
  - \a map
  - \a sub_logdens
  - \a sub_mean
  - \a sub_mean_constr
  - \a sub_graph 

  are needed by the user for extracting the results of the sampling routines. The \a sample is the
  sample itself (of length \a graph->n), and \a mean is the mean BEFORE any constraints and \a
  mean_constr is the mean AFTER (if any) constraints. If there is no constraints, then \a mean and
  \a mean_constr is identical. Note that both \a mean and \a mean_constr is of full length \a
  graph->n. The members \a sub_mean and \a sub_mean_constr are of length \a n_sub, where \a n_sub is
  the number of elements of the GMRF <em>\b x</em> that are updated during the sampling. Which
  elements of <em>\b x</em> that are updated, is given by the index vector \a map. More
  specifically, the elements <em> map[i]; i=0,..., n_sub-1</em> of <em>\b x</em> are updated, and
  these correspond to the elements <em>i=0,..., n_sub-1</em> of the mean vectors \a sub_mean and \a
  sub_mean_constr.
  
  \par Examples:
  See \ref ex_problem-setup
*/

typedef struct {

	/**
	 *  \brief The sample (full graph)
	 * 
	 * A length \em n array containing a sample of the GMRF. If some of the elements of <em>\b x</em> are fixed, these
	 * fixed values are copied into the corresponding positions in \c sample, while the rest are samples (or configuration
	 * to be evalutated) from the conditional distribution given the fixed values.
	 * 
	 * This array is just \ref GMRFLib_problem_tp::sub_sample copied into the correct locations, where the rest is fixed
	 * values.
	 * 
	 * \sa GMRFLib_problem_tp::sub_sample, GMRFLib_problem_tp::mean, GMRFLib_problem_tp::mean_constr 
	 */
	double *sample;

	/**
	 *  \brief The mean (full graph, BEFORE constraints)
	 * 
	 * A length \em n array containing a mean of the GMRF \c before correcting for constraints, if any.  If some of the
	 * elements of <em>\b x</em> are fixed, these fixed values are copied into the corresponding positions in \c mean,
	 * while the rest are the conditional mean from the conditional distribution given the fixed values.
	 * 
	 * This array is just \ref GMRFLib_problem_tp::sub_mean copied into the correct locations, where the rest is fixed
	 * values.
	 * 
	 * \sa GMRFLib_problem_tp::sub_mean, GMRFLib_problem_tp::mean_constr 
	 */
	double *mean;

	/**
	 *  \brief The mean (full graph, AFTER constraints)
	 * 
	 * A length \em n array containing a mean of the GMRF \c after correcting for constraints, if any. If there is no
	 * constraints, then \a mean_constr equals \a mean. If some of the elements of <em>\b x</em> are fixed, these fixed
	 * values are copied into the corresponding positions in \c mean_constr, while the rest are samples (or configuration
	 * to be evalutated) from the conditional distribution given the fixed values.
	 * 
	 * This array is just \ref GMRFLib_problem_tp::sub_mean_constr copied into the correct locations, where the rest is
	 * fixed values.
	 * 
	 * \sa GMRFLib_problem_tp::sub_mean_constr, GMRFLib_problem_tp::mean 
	 */
	double *mean_constr;

	/*
	 * for internal use only; the length of sample, mean and mean_constr. 
	 */
	int n;

	/**
	 *  \brief Mapping to the real-world
	 * 
	 * The mapping between the elements of \a sub_sample, \a sub_mean, \a sub_mean_constr and \a sub_graph, and the real
	 * world. The correspondence is defined by \f$ \mbox{\small\tt sub\_sample}[i] \leftrightarrow \mbox{\small\tt sample}
	 * [\mbox{\small\tt map}[i]] \f$, and similar with the other quantities. \n\n 
	 */
	int *map;

	/**
	 *  \brief The log-density of sub_sample, corresponding to the elements of <em>\b x</em> that are not fixed. 
	 */
	double sub_logdens;

	/**
	 *  \brief A new sample or the configuration which log-density have been evaluated.
	 * 
	 * This vector is mainly used internally, use instead GMRFLib_problem_tp::sample
	 * 
	 * \sa GMRFLib_problem_tp::map, GMRFLib_problem_tp::sample 
	 */
	double *sub_sample;

	/**
	 *  \brief The (conditional) mean before constraint
	 * 
	 * The computed conditional mean (given the values of the fixed elements). This is the mean values of an \em
	 * unconstrained sample, that is, before taking into account any constraints.
	 * 
	 * \sa GMRFLib_problem_tp::map, GMRFLib_problem_tp::sub_mean_constr 
	 */
	double *sub_mean;

	/**
	 *  \brief The mean (subgraph) after constraint
	 * 
	 * The mean of the variables corresponding to the sampled values in \a sample. This is the conditional mean of the
	 * sampled elements of <em>\b x</em>, given the values of the fixed elements and \em after taking any constraints into
	 * account.
	 * 
	 * \sa GMRFLib_problem_tp::sub_mean, GMRFLib_problem_tp::map 
	 */
	double *sub_mean_constr;

	/**
	 *  \brief The Choleskty triangle of Q for the subgraph 
	 */

	GMRFLib_sm_fact_tp sub_sm_fact;

	/**
	 *  \brief The constraint, if any 
	 */
	GMRFLib_constr_tp *sub_constr;

	/**
	 *  \brief The value of Ax - e
	 */
	double *sub_constr_value;

	/**
	 *  \brief The value of \f$ Q^{-1}A^T inv(AQ^{-1}A^T) \f$ 
	 */
	double *constr_m;

	/**
	 *  \brief The Cholesky triangle of \f$ AQ^{-1}A^T \f$ 
	 */
	double *l_aqat_m;

	/**
	 *  \brief  \f$ (AQ^{-1}A^T)^{-1} \f$.
	 *
	 *  This is only used for INLA() and added in there
	 */
	double *inv_aqat_m;

	/**
	 *  \brief The value of \f$ Q^{-1}A^T \f$ 
	 */
	double *qi_at_m;

	/**
	 *  \brief The value of \f$ \log (det(AA^T)) \f$ 
	 */
	double logdet_aat;

	/**
	 *  \brief The value of \f$ \log (det(AQA^T)) \f$ 
	 */
	double logdet_aqat;

	/**
	 *  \brief The value of \f$ \log (det(Q)^{1/2}) \f$ 
	 */
	double log_normc;

	/**
	 *  \brief The value of \f$x^TQx\f$ 
	 */
	double exp_corr;

	/**
	 *  \brief The subgraph
	 * 
	 * The subgraph defined on the nodes for which the elements of <em>\b x</em> are not fixed. This can be used to extract 
	 * \f$ \mbox{\small\tt sub\_graph} \rightarrow \mbox{\small\tt n} \f$, the number of elements of the graph, equal to
	 */
	GMRFLib_graph_tp *sub_graph;

	/**
	 *  \brief The tabulated Q-function on the sub_graph 
	 */
	GMRFLib_tabulate_Qfunc_tp *tab;

	/**
	 *  \brief The (structural) inverse of Q 
	 */
	GMRFLib_Qinv_tp *sub_inverse;
} GMRFLib_problem_tp;

/*!
  \brief To flag that the proposal was accepted
*/
#define GMRFLib_STORE_ACCEPT (1)

/*!
  \brief To flag that the proposal was rejected
*/
#define GMRFLib_STORE_REJECT (2)

/*!
  \struct GMRFLib_store_tp problem-setup.h
  \brief The structure to store intermediate calculations
*/
typedef struct GMRFLib_store_struct GMRFLib_store_tp;	       /* used recursively */
							       
struct GMRFLib_store_struct {
	GMRFLib_smtp_tp smtp;				       /* sparse matrix type */
	int bandwidth;					       /* for GMRFLib_smtp == GMRFLib_SMTP_BAND */
	int *remap;
	int copy_ptr;
	int copy_pardiso_ptr;
	GMRFLib_graph_tp *sub_graph;

	supernodal_factor_matrix *TAUCS_symb_fact;	       /* for GMRFLib_smtp == GMRFLib_SMTP_TAUCS */
	GMRFLib_pardiso_store_tp *PARDISO_fact;
	
	GMRFLib_store_tp *diag_store;			       /* store SAFE-optims in optimize */
	GMRFLib_store_tp *sub_store;			       /* store the same if fixed values in optimize */

	int store_problems;				       /* store problems or not? */
	int fixed_hyperparameters;			       /* if true, both problems are the same and fixed */
	int decision;					       /* what is the decision? accept or reject? */
	double *old_logdens;				       /* old log-density */
	double *new_logdens;				       /* new log-density */
	GMRFLib_problem_tp *problem_old2new;		       /* stored problem */
	GMRFLib_problem_tp *problem_new2old;		       /* stored problem */

};

#define STOCHASTIC_CONSTR(constr) ((constr) && ((constr)->errcov_diagonal || (constr)->errcov_general))

double *GMRFLib_Qinv_get(GMRFLib_problem_tp * problem, int i, int j);
double GMRFLib_Qfunc_wrapper(int sub_node, int sub_nnode, void *arguments);
int GMRFLib_Qinv(GMRFLib_problem_tp * problem, int storage);
int GMRFLib_eval_constr(double *value, double *sqr_value, double *x, GMRFLib_constr_tp * constr, GMRFLib_graph_tp * graph);
int GMRFLib_evaluate(GMRFLib_problem_tp * problem);
int GMRFLib_evaluate__intern(GMRFLib_problem_tp * problem, int compute_const);
int GMRFLib_fact_info_report(FILE * fp, GMRFLib_sm_fact_tp * sm_fact);
int GMRFLib_free_Qinv(GMRFLib_problem_tp * problem);
int GMRFLib_free_constr(GMRFLib_constr_tp * constr);
int GMRFLib_free_problem(GMRFLib_problem_tp * problem);
int GMRFLib_free_store(GMRFLib_store_tp * store);
int GMRFLib_info_problem(FILE * fp, GMRFLib_problem_tp * problem);
int GMRFLib_init_problem(GMRFLib_problem_tp ** problem, double *x, double *b, double *c, double *mean, GMRFLib_graph_tp * graph,
			 GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_args, char *fixed_value, GMRFLib_constr_tp * constraint, unsigned int keep);
int GMRFLib_init_problem_store(GMRFLib_problem_tp ** problem, double *x, double *b, double *c, double *mean,
			       GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_args, char *fixed_value,
			       GMRFLib_constr_tp * constr, unsigned int keep, GMRFLib_store_tp * store);
int GMRFLib_make_empty_constr(GMRFLib_constr_tp ** constr);
int GMRFLib_prepare_constr(GMRFLib_constr_tp * constr, GMRFLib_graph_tp * graph, int scale_constr);
int GMRFLib_print_constr(FILE * fp, GMRFLib_constr_tp * constr, GMRFLib_graph_tp * graph);
int GMRFLib_duplicate_constr(GMRFLib_constr_tp ** new_constr, GMRFLib_constr_tp * constr, GMRFLib_graph_tp * graph);
int GMRFLib_print_problem(FILE * fp, GMRFLib_problem_tp * problem);
int GMRFLib_recomp_constr(GMRFLib_constr_tp ** new_constr, GMRFLib_constr_tp * constr, double *x, double *b_add, char *mask,
			  GMRFLib_graph_tp * graph, GMRFLib_graph_tp * sub_graph);
int GMRFLib_sample(GMRFLib_problem_tp * problem);

GMRFLib_problem_tp *GMRFLib_duplicate_problem(GMRFLib_problem_tp * problem, int skeleton, int copy_ptr, int copy_pardiso_ptr);
GMRFLib_store_tp *GMRFLib_duplicate_store(GMRFLib_store_tp * store, int skeleton, int copy_ptr, int copy_pardiso_ptr);
double GMRFLib_Qfunc_generic(int i, int j, void *arg);
int GMRFLib_optimize_reorder(GMRFLib_graph_tp * graph, GMRFLib_sizeof_tp * nnz_opt, int *use_global, GMRFLib_global_node_tp *gn);
GMRFLib_sizeof_tp GMRFLib_sizeof_store(GMRFLib_store_tp * store);
GMRFLib_sizeof_tp GMRFLib_sizeof_problem(GMRFLib_problem_tp * problem);

__END_DECLS
#endif
