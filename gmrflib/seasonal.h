#ifndef __GMRFLib_SEASONAL_H__
#       define __GMRFLib_SEASONAL_H__

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
   \struct GMRFLib_seasonaldef_tp
   \brief The typedef used to define a seasonal model.

   \sa GMRFLib_seasonal(), GMRFLib_make_seasonal_graph()
 */
    typedef struct {

	/** 
	 * \brief The size of the graph
	 */
	int n;

	/** 
	 *\brief  The  length of the seasonal cycle
	 */
	int s;

	/**
	 * \brief A flag indicating wheather the graph is cyclical
	 *
	 * In a cyclic graph of lenght \c n  node \c n coincides with node 0.
	 * Use cyclic boundary conditions if cyclic=#GMRFLib_TRUE, otherwise do not use cyclic boundary conditions. 
	 */
	int cyclic;

	double **log_prec_omp;
	double *prec_scale;				       /* scaling of the precision for scale.model=TRUE */
} GMRFLib_seasonaldef_tp;

double GMRFLib_seasonal(int thread_id, int node, int nnode, double *values, void *seasonal_def);
int GMRFLib_make_seasonal_graph(GMRFLib_graph_tp ** graph, GMRFLib_seasonaldef_tp * def);
int GMRFLib_seasonal_scale(int thread_id, GMRFLib_seasonaldef_tp * def);

__END_DECLS
#endif
