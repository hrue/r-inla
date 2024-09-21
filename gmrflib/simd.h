
/* simd.h
 * 
 * Copyright (C) 2024-2024 Havard Rue
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
  \file utils.h
  \brief Typedefs for \ref utils.c
*/

#ifndef __GMRFLib_SIMD_H__
#define __GMRFLib_SIMD_H__

#include <stdlib.h>
#include <stddef.h>
#include <math.h>

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
#include "GMRFLib/GMRFLibP.h"

#if defined(INLA_WITH_SIMD)
#define SLEEF_ENABLE_OMP_SIMD
#include <x86intrin.h>
#include <sleef.h>
#endif

void GMRFLib_exp(int, double *, double *);
void GMRFLib_exp_inc(int n, double *x, int inc, double *y);
void GMRFLib_log(int, double *, double *);
void GMRFLib_log1p(int, double *, double *);
void GMRFLib_sqr(int n, double *x, double *y);
void GMRFLib_sqrt(int n, double *x, double *y);
void GMRFLib_add(int n, double *x, double *y, double *z);
void GMRFLib_mul(int n, double *x, double *y, double *z);
void GMRFLib_daddto(int n, double *x, double *y);
void GMRFLib_cdaddto(int n, double *x, double cx, double *y);

__END_DECLS
#endif
