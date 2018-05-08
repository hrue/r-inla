
/* rw.h
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
  \file rw.h
  \brief Typedefs used to define RW1, RW2, CRW1 and CRW2 models.
*/

#ifndef __GMRFLib_RW_H__
#define __GMRFLib_RW_H__

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
  \struct GMRFLib_rwdef_tp rw.h
  \brief The typedef used to define RW1 and RW2 models with regular locations.

  This defines the RW1 and RW2 models with regular locations, which is an argument to \ref GMRFLib_rw(), to compute an entry in
  the presicion matrix.  \ref GMRFLib_make_rw_graph() can be used to create the appropirate graph.

  \sa \ref GMRFLib_crwdef_tp, \ref GMRFLib_rw(),  \ref GMRFLib_make_rw_graph()
*/

/*
 */
    typedef struct {

	/**
	 *  \brief The size of the graph or the number of locations 
	 */
	int n;

	/**
	 *  \brief The order or the random walk. Must be 0, 1 or 2.
	 */
	int order;

	/**
	 *  \brief A flag for using cyclic boundary conditions or not.
	 * 
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

	/*
	 * scale the prec: only for order = 0. Not public
	 */
	double *scale0;

	double *prec_scale;

} GMRFLib_rwdef_tp;

/*!
  \brief Define the layout for the CRW2 model: Layout z[0], ..., z[n-1].
*/
#define GMRFLib_CRW_LAYOUT_SIMPLE 0

/*!
  \brief Define the layout for the CRW2 model: Layout z[0], z'[0], z[1], z'[1]....
*/
#define GMRFLib_CRW_LAYOUT_PAIRS  1

/*!
  \brief Define the layout for the CRW2 model: Layout z[0], ..., z[n-1], z'[0], ..., z'[n-1]
*/
#define GMRFLib_CRW_LAYOUT_BLOCK  2

/*!
  \struct GMRFLib_crwdef_tp rw.h

  \brief The typedef used to define CRW1 and CRW2 models with irregular locations.

  This defines the CRW1, CRW2 and approximate CRW2 models with irregular locations, which is
  an argument to to \ref GMRFLib_crw(), to compute an entry in the presicion matrix and \ref
  GMRFLib_make_crw_graph() which create the appropirate graph.

  \sa GMRFLib_rwdef_tp, GMRFLib_crw(), GMRFLib_make_crw_graph()
*/
typedef struct {

	/**
	*  \brief The size of the graph or the number of locations 
	 */
	int n;

	/**
	*  \brief The order or the random walk. Must be 0, 1 or 2.
	 */
	int order;

	/**
	 * \brief Use a dense representation to allow for 'si'. (EXPERT Option)
	 */
	int si;

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

	 * if \c log_prec_omp is \c NULL, then a unit precision is used.
	 */
	double **log_prec_omp;

	/**
	 *  \brief An array with the positions for x[0]...x[n-1]
	 * 
	 * \c position[i] is the position to x[i], for i=0...n-1.  <em> It is assumed, and not checked for, that the positions
	 * are increasing </em>, i.e. position[0] < position[1] < ... < position[n-1].  If position is \c NULL, then the
	 * positions are assumed to be 0, 1, ..., n-1. 
	 */
	double *position;

	/**
	 *  \brief Define the memory layout for the CRW2 model.
	 * 
	 * This argument defines the memory layout for the CRW2 model. If it is equal to #GMRFLib_CRW_LAYOUT_SIMPLE, then the
	 * approximative model is used where no augmentation with velocities are needed. If it is equal to
	 * #GMRFLib_CRW_LAYOUT_BLOCK or #GMRFLib_CRW_LAYOUT_PAIRS, then the exact augmented solution is used.
	 * 
	 * The exact solution of the CRW2 model augment the state-vector z with its derivatives z' (or velocities). Both
	 * vectors are stored in the general array x which is of length 2*n. The term ij in the presicion matrix, and also the
	 * graph itself, depends then on how z and z' are stored in x.
	 * 
	 * If \c layout = #GMRFLib_CRW_LAYOUT_PAIRS then <em>x = (z[0], z'[0], z[1], z'[1], ...., z[n-1], z'[n-1])</em>.
	 * 
	 * If \c layout = #GMRFLib_CRW_LAYOUT_BLOCK then <em>x = (z[0], z[1], ...., z[n-1], z'[0], z'[1], ..., z'[n-1])</em>.
	 * 
	 * If \c layout = #GMRFLib_CRW_LAYOUT_SIMPLE then <em>x = (z[0], z[1], ...., z[n-1])</em>.
	 * 
	 */
	int layout;

	/**
	*  \brief Work array for the \c GMRFLib_crw() function.
	 * 
	 * The array \c work is a working array for the \c GMRFLib_crw() function. It must be initialised to \c NULL, and then
	 * \c GMRFLib_rw() allocate the space needed and calculate the contents. This array can be free'ed using \c free().
	 * 
	 * The workarray depends ONLY of the positions, and does NOT depend of the precision, order or the layout.
	 * 
	 */
	double *work;

	/*
	 * scale the prec: only for order = 0. Not public
	 */
	double *scale0;

	/* 
	 * scale the prec: for any order. 
	 */
	double *prec_scale;
} GMRFLib_crwdef_tp;


typedef enum {
	GMRFLib_BVALUE_DEFAULT = 0,			       /* do not change this */
	GMRFLib_BVALUE_ZERO = 1
} GMRFLib_rw2d_bvalue_tp;


/*!
  \struct GMRFLib_rw2ddef_tp rw.h
  \brief The typedef used to define RW models in 2D on lattices.

  This defines the RW models on 2D on regular locations, which is an argument to \ref GMRFLib_rw2d(), to compute an entry in the
  presicion matrix.  \ref GMRFLib_make_rw2d_graph() can be used to create the appropirate graph. Currently only order=2 is
  supported, however, it does incorporate correct boundary conditions for the non-cyclic case.

  \sa \ref GMRFLib_rw2d(), \ref GMRFLib_make_lattice_graph()
*/
typedef struct {

	/**
	*  \brief The size of the lattice 
	 */
	int nrow, ncol;

	/**
	*  \brief The order or the random walk. Since only order=2 is currrently supported, this option is currently not in
	 * use. It might be so in the future. 
	 */
	int order;

	/**
	*  \brief A flag for using cyclic boundary conditions or not.
	 * 
	 * Use cyclic boundary conditions if cyclic=\c TRUE, otherwise do not use cyclic boundary conditions. 
	 */
	int cyclic;

	/**
	 * \brief Choice of boundary values: #GMRFLib_BVALUE_DEFAULT use the correct null-space, whereas #GMRFLib_BVALUE_ZERO condition on zero outside the lattice.
	 */
	int bvalue;

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

	double *prec_scale;
} GMRFLib_rw2ddef_tp;

double GMRFLib_rw(int node, int nnode, void *rwdef);
double GMRFLib_crw(int node, int nnode, void *crwdef);
double GMRFLib_rw2d(int node, int nnode, void *rw2ddef);
int GMRFLib_make_rw2d_graph(GMRFLib_graph_tp ** graph, GMRFLib_rw2ddef_tp * def);
int GMRFLib_make_rw_graph(GMRFLib_graph_tp ** graph, GMRFLib_rwdef_tp * def);
int GMRFLib_make_crw_graph(GMRFLib_graph_tp ** graph, GMRFLib_crwdef_tp * def);
int GMRFLib_crw_scale_OLD(void *def);
int GMRFLib_crw_scale(void *def);
int GMRFLib_rw_scale(void *def);
int GMRFLib_rw2d_scale(void *def);

__END_DECLS
#endif
