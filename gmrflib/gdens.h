
/* gdens.h
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
  \file gdens.h
  \brief Typedefs and defines for \ref gdens.c
*/

#ifndef __GMRFLib_GDENS_H__
#define __GMRFLib_GDENS_H__

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

__BEGIN_DECLS typedef struct {
	double x, xl, xm, xr, dx, b, c, lp, lfl, lfr, lfm, lfmax, lnormc;
} GMRFLib_gdens_Elm;

typedef struct {
	int n;
	int order;
	double xmin, xmax;
	GMRFLib_gdens_Elm *elm;
} GMRFLib_gdens_tp;

typedef double GMRFLib_gdens_Func_tp(double, void *);

double GMRFLib_erfi(double x);
double GMRFLib_lerf_diff(double x1, double x0);
double GMRFLib_lerfi_diff(double x1, double x0);
double *GMRFLib_cut_points(double fac, int n);
int GMRFLib_gdens_1LDens(GMRFLib_gdens_tp * ptr, double *x, double *ldens);
int GMRFLib_gdens_1Sample(GMRFLib_gdens_tp * ptr, double *x, double *ldens);
int GMRFLib_gdens_1Eval(GMRFLib_gdens_tp * ptr, double *x, double *ldens, int flag);
GMRFLib_gdens_tp *GMRFLib_gdens_1Init(double xmin, double xmax, int n, int n_min, GMRFLib_gdens_Func_tp * lfunc, void *lfunc_arg);
int GMRFLib_gdens_1Update(GMRFLib_gdens_tp * ptr);

int GMRFLib_gdens_2LDens(GMRFLib_gdens_tp * ptr, double *x, double *ldens);
int GMRFLib_gdens_2Sample(GMRFLib_gdens_tp * ptr, double *x, double *ldens);
int GMRFLib_gdens_2Eval(GMRFLib_gdens_tp * ptr, double *x, double *ldens, int flag);
GMRFLib_gdens_tp *GMRFLib_gdens_2Init(double xmin, double xmax, int n, int n_min, GMRFLib_gdens_Func_tp * lfunc, void *lfunc_arg);
GMRFLib_gdens_tp *GMRFLib_gdens_2InitNew(double fac, double mean, double stdev, int n, GMRFLib_gdens_Func_tp * lfunc, void *lfunc_arg);
int GMRFLib_gdens_2Update(GMRFLib_gdens_tp * ptr);

int GMRFLib_gdens_LDens(GMRFLib_gdens_tp * ptr, double *x, double *ldens);
int GMRFLib_gdens_Sample(GMRFLib_gdens_tp * ptr, double *x, double *ldens);
int GMRFLib_gdens_Eval(GMRFLib_gdens_tp * ptr, double *x, double *ldens, int flag);
GMRFLib_gdens_tp *GMRFLib_gdens_Init(double xmin, double xmax, int n, int n_min, int order, GMRFLib_gdens_Func_tp * lfunc, void *lfunc_arg);
GMRFLib_gdens_tp *GMRFLib_gdens_InitNew(double fac, double mean, double stdev, int n, GMRFLib_gdens_Func_tp * lfunc, void *lfunc_arg);

int GMRFLib_gdens_Free(GMRFLib_gdens_tp * ptr);
int GMRFLib_gdens_ElmCompare(const void *e, const void *ee);

__END_DECLS
#endif
