
/* matern.h
 * 
 * Copyright (C) 2008 Havard Rue
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
  \file matern.h
  \brief Typedefs used to define Matern-models
*/

#ifndef __GMRFLib_MATERN_H__
#define __GMRFLib_MATERN_H__

#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>

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
	 * tried.	 * tried.
	 */
	double *log_prec;

	/**
	 *  \brief A (possible) ppointer to the log-precision where each tread has its own value.
	 * 
	 * If \c log_prec_omp !=\c NULL, then \c log_prec_omp[ID] points to the log-precision, where ID is \c GMRFLib_thread_id.
	 * if \c log_prec_omp is \c NULL, then a unit precision is used.
	 */
	double **log_prec_omp;

	/**
	 *  \brief A (possible) pointer to the range
	 * 
	 * If prec!=\c NULL, then this pointer points to the range. Alternatively, \c log_range is tried.
	 */
	double *range;

	/**
	 *  \brief A (possible) pointer to the log-range
	 * 
	 * If log_range!=\c NULL, then this pointer points to the log-range. Alternatively, \c log_range_omp is
	 * tried.	 
	 */
	double *log_range;

	/**
	 *  \brief A (possible) ppointer to the log-range where each tread has its own value.
	 * 
	 * If \c log_range_omp !=\c NULL, then \c log_range_omp[ID] points to the log-range, where ID is \c GMRFLib_thread_id.
	 * if \c log_range_omp is \c NULL, then \c log-range is 0.
	 */
	double **log_range_omp;
} GMRFLib_matern2ddef_tp;

double GMRFLib_matern2d(int node, int nnode, void *def);
int GMRFLib_make_matern2d_graph(GMRFLib_graph_tp ** graph, GMRFLib_matern2ddef_tp * def);


__END_DECLS
#endif
