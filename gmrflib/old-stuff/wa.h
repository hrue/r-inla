#ifndef __GMRFLib_WA_H__
#define __GMRFLib_WA_H__

#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>
#include <stdio.h>

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
#include "GMRFLib/hashP.h"
    typedef struct {
	int node_nnode;					       /* flag if these are neigb or not */
	int nnode_node;					       /* flag if these are neigb or not */
	int *node_list;					       /* list of terms to include */
	int n_nodes;					       /* number of terms in the list */
} GMRFLib_node_list_tp;

typedef struct {
	GMRFLib_graph_tp *wagraph;
	GMRFLib_Qfunc_tp *waQfunc;
	void *waQfunc_arg;

	GMRFLib_graph_tp *waQgraph;			       /* graph before prune */
	spmatrix neigb_idx_hash;			       /* keep hash table, (i,j)->idx */
	int n_neigb_info;				       /* indexs */
	GMRFLib_node_list_tp **neigb_info;		       /* keep neigbors for each idx */
} GMRFLib_waQfunc_arg_tp;

/*!
  \struct GMRFLib_wa_problem_tp wa.h 
  \brief The resulting graph and Q-function (with arguments) for the equivalent formulation \f$x^TQx\f$

  The members of the data structure can be used the graph and the Q-function and its
  arguments for a general graph using the  \f$x^TQx\f$ formulation.
 */
typedef struct {

	/**
	 *  \brief The graph of \f$Q\f$ in the equivalent formulation \f$x^TQx\f$ 
	 */
	GMRFLib_graph_tp *graph;

	/**
	 *  \brief The <em>Qfunc</em>ion in the equivalent formulation \f$x^TQx\f$ 
	 */
	GMRFLib_Qfunc_tp *Qfunc;

	/**
	 *  \brief The arguments to GMRFLib_waQfunc_arg_tp::waQfunc in the equivalent formulation \f$x^TQx\f$
	 */
	void *Qfunc_arg;
} GMRFLib_wa_problem_tp;

/*!
  \struct GMRFLib_nwa_problem_tp wa.h 
  \brief  The resulting graph and Q-function (with arguments) for the equivalent formulation \f$x^TQx\f$

  The members of the data structure can be used the graph and the Q-function and its arguments for a
  general graph using the \f$x^TQx\f$ formulation.

 */
typedef struct {

	/**
	 *  \brief The graph of \f$Q\f$ in the equivalent formulation \f$x^TQx\f$ 
	 */
	GMRFLib_graph_tp *graph;

	/**
	 *  \brief The <em>Qfunc</em>ion in the equivalent formulation \f$x^TQx\f$ 
	 */
	GMRFLib_Qfunc_tp *Qfunc;

	/**
	 *  \brief The arguments to GMRFLib_waQfunc_arg_tp::waQfunc in the equivalent formulation \f$x^TQx\f$ 
	 */
	void *Qfunc_arg;

	/**
	 *  \brief Flag if this is a wa_problem or not. 
	 */
	int wa_problem;
} GMRFLib_nwa_problem_tp;

GMRFLib_Qfunc_tp GMRFLib_waQfunc, GMRFLib_nwaQfunc;
int GMRFLib_init_wa_problem(GMRFLib_wa_problem_tp ** wa_problem, GMRFLib_graph_tp * wagraph, GMRFLib_Qfunc_tp * wafunc, void *wafunc_arg);
int GMRFLib_free_wa_problem(GMRFLib_wa_problem_tp * wa_problem);
int GMRFLib_init_nwa_problem(GMRFLib_nwa_problem_tp ** nwa_problem,
			     int n_wa, GMRFLib_graph_tp ** wagraph, GMRFLib_Qfunc_tp ** wafunc, void **wafunc_arg,
			     int n_g, GMRFLib_graph_tp ** graph, GMRFLib_Qfunc_tp ** Qfunc, void **Qfunc_arg);
int GMRFLib_free_nwa_problem(GMRFLib_nwa_problem_tp * nwa_problem);

__END_DECLS
#endif
