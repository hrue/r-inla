
/* fgn.h
 * 
 * Copyright (C) 2016 Havard Rue
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
 */
#ifndef __INLA_FGN_H__
#define __INLA_FGN_H__
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
 *
 */
#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"
#include "inla.h"
#define FGN_KMAX (4L)					       /* maximum K in the tables */

//
    typedef struct {
	int N;						       /* total size = (k+1)*n */
	int n;						       /* size of each component */
	int k;						       /* number of components */
	double prec_eps;				       /* fixed tiny noise */
	double **log_prec;				       /* theta[0] */
	double **H_intern;				       /* theta[1] */
} inla_fgn_arg_tp;

typedef struct {
	int N;						       /* total size = k*n */
	int n;						       /* size of each component */
	int k;						       /* number of components */
	double **log_prec;				       /* theta[0] */
	double **H_intern;				       /* theta[1] */
} inla_fgn2_arg_tp;

double inla_fgn2_helper(int i, int j, int n, double phi);

int inla_make_fgn_graph(GMRFLib_graph_tp ** graph, inla_fgn_arg_tp * def);
int inla_make_fgn2_graph(GMRFLib_graph_tp ** graph, inla_fgn2_arg_tp * def);
double Qfunc_fgn(int i, int j, void *arg);
double Qfunc_fgn2(int i, int j, void *arg);

int inla_fgn_get(double *phi, double *w, double H_intern, int k);
double priorfunc_fgn_priorH(double *H_intern, double *param);

__END_DECLS
#endif
