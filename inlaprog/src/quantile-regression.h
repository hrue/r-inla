
/* quantile-regression.h
 * 
 * Copyright (C) 2016-2019 Havard Rue
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
#ifndef __INLA_QUANTILE_REGRESSION_H__
#define __INLA_QUANTILE_REGRESSION_H__
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

#include "GMRFLib/density.h"

struct inla_qgamma_cache_tp {
	double quantile;
	GMRFLib_spline_tp **s;
};

GMRFLib_spline_tp **inla_qcontpois_func(double alpha, int num);
double inla_pcontpois(double y, double lambda);
double inla_pcontpois_deriv(double y, double lambda);
double inla_qcontpois(double quantile, double alpha, double *initial_guess);
double inla_qcontpois_eta(double quantile, double alpha, double *initial_guess);
double inla_qgamma_cache(double shape, double quantile, int id);

__END_DECLS
#endif
