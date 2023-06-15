
/* spde2.h
 * 
 * Copyright (C) 2011-2023 Havard Rue
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
#ifndef __INLA_SPDE2_H__
#define __INLA_SPDE2_H__
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
#include "GMRFLib/GMRFLib.h"

/* 
   
 */
    typedef enum {
	SPDE2_TRANSFORM_LOGIT = 1,			       /* cos(pi/(1+exp(-x)) */
	SPDE2_TRANSFORM_LOG,				       /* 2*exp(x)-1 */
	SPDE2_TRANSFORM_IDENTITY			       /* x */
} spde2_transform_tp;

typedef struct {
	int i;
	int need_transform;
	double *theta;
	double *vals;

	// need a buffer...
	double l1_cachline[4L];
} spde2_cache_tp;

typedef struct {
	int n;
	int ntheta;					       /* that is `p' in Finn's notes */
	int ntheta_used;				       /* how many are non-fixed. */
	int debug;
	int theta_first_idx;

	int *fixed;
	double *fixed_values;

	spde2_transform_tp transform;

	GMRFLib_matrix_tp **B;
	GMRFLib_matrix_tp **M;
	GMRFLib_matrix_tp *BLC;

	double ***theta;

	GMRFLib_Qfunc_tp *Qfunc;
	void *Qfunc_arg;
	GMRFLib_graph_tp *graph;

	GMRFLib_vmatrix_tp *vmatrix;
	spde2_cache_tp ***cache;
} inla_spde2_tp;

double inla_spde2_Qfunction_orig(int thread_id, int node, int nnode, double *values, void *arg);
double inla_spde2_Qfunction_cache(int thread_id, int ii, int jj, double *values, void *arg);
double inla_spde2_Qfunction_old(int thread_id, int ii, int jj, double *values, void *arg);
double inla_spde2_Qfunction_ref(int thread_id, int ii, int jj, double *values, void *arg);
double inla_spde2_Qfunction(int thread_id, int ii, int jj, double *values, void *arg);
double *inla_spde2_userfunc2(int number, double *theta, int nhyper, double *covmat, void *arg);
int inla_spde2_build_model(int thread_id, inla_spde2_tp ** smodel, const char *prefix, const char *transform);

__END_DECLS
#endif
