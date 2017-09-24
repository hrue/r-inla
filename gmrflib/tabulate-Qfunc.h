
/* tabulate-Qfunc.h
 * 
 * Copyright (C) 2004-2006 Havard Rue
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
  \file tabulate-Qfunc.h
  \brief Typedefs for \ref tabulate-Qfunc.c
*/

#ifndef __GMRFLib_TABULATE_QFUNC_H__
#define __GMRFLib_TABULATE_QFUNC_H__

#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>
#include <zlib.h>

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

/* 
   dummy comment
 */
    typedef struct {
	int n;						       /* the size of the graph */
	map_id **values;				       /* hash-table for the values */
	double *prec;					       /* precision */
	double *log_prec;				       /* log(prec) */
	double **log_prec_omp;				       /* log(prec) thread dependent */
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

double GMRFLib_tabulate_Qfunction(int node, int nnode, void *arg);
double GMRFLib_tabulate_Qfunction_std(int node, int nnode, void *arg);
int GMRFLib_free_tabulate_Qfunc(GMRFLib_tabulate_Qfunc_tp * tabulated_Qfunc);
int GMRFLib_tabulate_Qfunc(GMRFLib_tabulate_Qfunc_tp ** tabulated_Qfunc, GMRFLib_graph_tp * graph,
			   GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg, double *prec, double *log_prec, double **log_prec_omp);
int GMRFLib_tabulate_Qfunc_from_file_OLD(GMRFLib_tabulate_Qfunc_tp ** tabulate_Qfunc, GMRFLib_graph_tp * graph,
					 const char *filename, double *prec, double *log_prec, double **log_prec_omp);
int GMRFLib_tabulate_Qfunc_from_file(GMRFLib_tabulate_Qfunc_tp ** tabulate_Qfunc, GMRFLib_graph_tp ** graph, const char *filename, int dim,
				     double *prec, double *log_prec, double **log_prec_omp);
int GMRFLib_tabulate_Qfunc_from_list(GMRFLib_tabulate_Qfunc_tp ** tabulate_Qfunc, GMRFLib_graph_tp ** graph, int ntriples,
				     int *ilist, int *jlist, double *Qijlist, int dim, double *prec, double *log_prec, double **log_prec_omp);
__END_DECLS
#endif
