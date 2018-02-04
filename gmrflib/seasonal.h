
/* seasonal.h
 * 
 * Copyright (C) 2007 Havard Rue
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
  \file seasonal.h
  \brief Typedefs for \c seasonal.c
*/

#ifndef __GMRFLib_SEASONAL_H__
#define __GMRFLib_SEASONAL_H__

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

	/**
	 *  \brief A (possible) pointer to the precision 
	 * 
	 * If prec!=\c NULL, then this pointer points to the precision. Alternatively, \c log_prec is tried.
	 */
	double *prec;

	/**
	 *  \brief A (possible) pointer to the log-precision 
	 * 
	 * If log_prec!=\c NULL, then this pointer points to the log-precision. Alternatively, \c log_prec_omp is
	 * tried.
	 */
	double *log_prec;

	/**
	 *  \brief A (possible) ppointer to the log-precision where each tread has its own value.
	 * 
	 * If \c log_prec_omp !=\c NULL, then \c log_prec_omp[ID] points to the log-precision, where ID is \c GMRFLib_thread_id.
	 * if \c log_prec_omp is \c NULL, then a unit precision is used.
	 */
	double **log_prec_omp;


	double *prec_scale;				       /* scaling of the precision for scale.model=TRUE */
} GMRFLib_seasonaldef_tp;

double GMRFLib_seasonal(int node, int nnode, void *seasonal_def);
int GMRFLib_make_seasonal_graph(GMRFLib_graph_tp ** graph, GMRFLib_seasonaldef_tp * def);
int GMRFLib_seasonal_scale(GMRFLib_seasonaldef_tp *def);

__END_DECLS
#endif
