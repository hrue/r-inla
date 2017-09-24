
/* integrator.h
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

/* 
   This file is just for adapting the original code in integrator.c to GMRFLib.
 */

#ifndef __GMRFLib_INTEGRATOR_H__
#define __GMRFLib_INTEGRATOR_H__

#include <math.h>
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

/* 
   dummy comment
 */
typedef double (*integrand) (unsigned ndim, const double *x, void *);

/* Integrate the function f from xmin[dim] to xmax[dim], with at most maxEval function evaluations (0 for no limit), until the
   given absolute or relative error is achieved.  val returns the integral, and err returns the estimate for the absolute error
   in val.  The return value of the function is 0 on success and non-zero if there was an error. */

int adapt_integrate(integrand f, void *fdata, unsigned dim, const double *xmin, const double *xmax,
		    unsigned maxEval, double reqAbsError, double reqRelError, double *val, double *err);

__END_DECLS
#endif
