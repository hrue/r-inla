
/* sn-g.h
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

#ifndef __GMRFLib_SN_G_H__
#define __GMRFLib_SN_G_H__

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

__BEGIN_DECLS int GMRFLib_sn_g_get_order(void);

//
double *GMRFLib_sn_g_get_coof(double skew, double *cx);
double *GMRFLib_sn_ginv_get_coof(double skew, double *cx);
double GMRFLib_sn_g_eval(double x, double *cx);
double GMRFLib_sn_g_eval_deriv(double x, double *cx);
//

__END_DECLS
#endif
