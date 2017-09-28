
/* geo.h
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

/* This code was initially written by Hanne T Wist, later modified by H. Rue */

/*!
  \file geo.h
  \brief Typedefs and defines for \ref geo.c
*/

#ifndef __GMRFLib_GEO_H__
#define __GMRFLib_GEO_H__

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
  \brief Use Matern CF
*/
#define GMRFLib_CORTP_MATERN  (1)

/*!
  \brief Use Exponential CF
*/
#define GMRFLib_CORTP_EXP     (2)

/*!
  \brief Use Gaussian CF
*/
#define GMRFLib_CORTP_GAUSS   (3)

/*!
  \struct GMRFLib_geo_problem_tp geo.h
  \brief Specification of the geo sampling problem.
*/
    typedef struct {

	/**
	 *  \brief The graph 
	 */
	GMRFLib_graph_tp *graph;

	/**
	 *  \brief The Q-function, computing Q(i,j). 
	 */
	GMRFLib_Qfunc_tp *Qfunc;

	/**
	 *  \brief The arguments to the Q-function. 
	 */
	void *Qfunc_arg;
} GMRFLib_geo_problem_tp;

/*!
  \struct GMRFLib_geoQfunc_arg_tp geo.h
*/
typedef struct {
	int name;
	int neigh;
	double param;
	double range;
	int nrow;
	int ncol;
	double *prec;
	int cyclic_flag;				       /* if != 0, the grid is made cyclic */
	double *coef;
	map_ii *hash_table;
} GMRFLib_geoQfunc_arg_tp;

typedef struct {
	map_stri coef2;
	map_stri coef3;
} GMRFLib_Global_geo_tp;

/*
  format for makeing keys for the table-lookup: name(as int) param range

  NOTE we're using only 2 digits of accurancy in the representation (i.e. %.2f) as we only have various ranges and various
  smoothing paramtersq in the Matern in that resolution.
*/
#define GMRFLib_GEO_FMT  "%1d|%.2f|%.2f"

const char *GMRFLib_geo_translate_cortp(int name);
const char *GMRFLib_geo_translate_neigh(int neigh);
double GMRFLib_geoQfunc(int node, int nnode, void *arg);
int GMRFLib_free_geo_problem(GMRFLib_geo_problem_tp * geo_problem);
int GMRFLib_get_geo_coefs(double **coef, int name, int neigh, double param, double range);
int GMRFLib_get_geo_coefs2(double **coef, int name, double param, double range);
int GMRFLib_get_geo_coefs3(double **coef, int name, double param, double range);
int GMRFLib_init_geo_problem(GMRFLib_geo_problem_tp ** geo_problem, int name, int neigh, double param,
			     double range, int nrow, int ncol, double *prec, int cyclic_flag);
int GMRFLib_is_geo_coefs(int name, int neigh, double param, double range);
int GMRFLib_print_geo_coefs(FILE * fp);
int GMRFLib_revise_geo_problem(GMRFLib_geo_problem_tp * geo_problem, int name, double param, double range, double *prec);

__END_DECLS
#endif
