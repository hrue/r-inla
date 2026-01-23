#ifndef __GMRFLib_MATERN_H__
#       define __GMRFLib_MATERN_H__

#       include <stdlib.h>

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

/*!
  \struct GMRFLib_matern2ddef_tp 
  \brief The typedef used to define Matern models

  This defines Matern models models on a regular lattice. This is used as an argument to \ref GMRFLib_matern2d(), to compute an entry in
  the presicion matrix.  \ref GMRFLib_make_matern2d_graph() can be used to create the appropirate graph.

  \sa \ref GMRFLib_matern2d(),  \ref GMRFLib_make_matern2d_graph()
*/

/*
 */
    typedef struct {

	/**
	 *  \brief Number of rows in the lattice
	 */
	int nrow;

	/**
	 *  \brief Number of columns in the lattice
	 */
	int ncol;

	/**
	 *  \brief The parameter controlling the degree of smoothing in the corresponding countinous indexed field. Must be 0, 1, 2 or 3.
	 */
	int nu;

	/**
	 *  \brief A flag for using cyclic boundary conditions or not.
	 * 
	 * Use cyclic boundary conditions if cyclic= #GMRFLib_TRUE, otherwise do not use cyclic boundary conditions. 
	 */
	int cyclic;

	double **log_prec_omp;
	double **log_range_omp;
} GMRFLib_matern2ddef_tp;

double GMRFLib_matern2d(int thread_id, int node, int nnode, double *values, void *def);
int GMRFLib_make_matern2d_graph(GMRFLib_graph_tp ** graph, GMRFLib_matern2ddef_tp * def);

__END_DECLS
#endif
