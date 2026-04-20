
/*!
  \file tabulate-Qfunc.h
  \brief Typedefs for \ref tabulate-Qfunc.c
*/

#ifndef __GMRFLib_TABULATE_QFUNC_H__
#       define __GMRFLib_TABULATE_QFUNC_H__

#       include <stdlib.h>
#       include <zlib.h>

#       include "GMRFLib/hashP.h"
#       include "GMRFLib/GMRFLibP.h"
#       include "GMRFLib/GMRFLib.h"

#       undef __BEGIN_DECLS
#       undef __END_DECLS
#       ifdef __cplusplus
#              define __BEGIN_DECLS extern "C" {
#              define __END_DECLS }
#       else
#              define __BEGIN_DECLS			       /* empty */
#              define __END_DECLS			       /* empty */
#       endif

__BEGIN_DECLS

/* 
   dummy comment
 */
    typedef struct {
	int n;						       /* the size of the graph */
	map_id **values;				       /* hash-table for the values */
	double **log_prec_omp;				       /* log(prec) thread dependent */

	// new format
	GMRFLib_graph_tp *graph;
	GMRFLib_csr_tp *Q;
	map_ii **Q_idx;
} GMRFLib_tabulate_Qfunc_arg_tp;

/*!
   \struct GMRFLib_tabulate_Qfunc_tp

   \brief The structure retured by \c GMRFLib_tabulate_Qfunc()
*/

typedef struct {

	/**
	 *  \brief The Qfunction which returns the tabulated values 
	 */
	GMRFLib_Qfunc_tp *Qfunc;

	/**
	 *  \brief The arguments to GMRFLib_tabulate_Qfunc_tp::Qfunc 
	 */
	void *Qfunc_arg;
} GMRFLib_tabulate_Qfunc_tp;

double GMRFLib_tabulate_Qfunction(int thread_id, int node, int nnode, double *values, void *arg);
double GMRFLib_tabulate_Qfunction_std(int thread_id, int node, int nnode, double *values, void *arg);
int GMRFLib_free_tabulate_Qfunc(GMRFLib_tabulate_Qfunc_tp * tabulated_Qfunc);
int GMRFLib_tabulate_Qfunc(int thread_id, GMRFLib_tabulate_Qfunc_tp ** tabulated_Qfunc, GMRFLib_graph_tp * graph,
			   GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg, double **log_prec_omp);
int GMRFLib_tabulate_Qfunc_from_file(GMRFLib_tabulate_Qfunc_tp ** tabulate_Qfunc, GMRFLib_graph_tp ** graph, const char *filename, int dim,
				     double **log_prec_omp);
int GMRFLib_tabulate_Qfunc_from_list(GMRFLib_tabulate_Qfunc_tp ** tabulate_Qfunc, GMRFLib_graph_tp ** graph, int ntriples,
				     int *ilist, int *jlist, double *Qijlist, int dim, double **log_prec_omp);
int GMRFLib_tabulate_Qfunc_from_list2(GMRFLib_tabulate_Qfunc_tp ** tabulate_Qfunc, GMRFLib_graph_tp * graph,
				      int ntriples, int *ilist, int *jlist, double *Qijlist, int dim, double **log_prec_omp);
int GMRFLib_tabulate_Qfunc_core(int thread_id, GMRFLib_tabulate_Qfunc_tp ** tabulate_Qfunc, GMRFLib_graph_tp * graph,
				GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg, double **log_prec_omp, int force);
__END_DECLS
#endif
