/* pc-bgev-link.c
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
 */

#include "inla.h"
#include "rmath.h"
#include "pc-bgev-link.h"

#define XI_MIN 1.0E-8
#define P_MIN 1.0E-8

double inla_log_pbgev(double y, double xi)
{
#define _log_pfrechet(y_) (-DMAX(0.0, pow(((y_) - qalpha) / (sbeta / (l2_xi - l3_xi)) + l1_xi, -1.0/xi)))
	xi = DMAX(XI_MIN, xi);
	double qalpha = 0, alpha = 0.5, beta = 0.5, sbeta = 1.0; 

	double l1_xi = pow(-log(alpha), -xi);
	double l2_xi = pow(-log(1.0 - beta / 2.0), -xi);
	double l3_xi = pow(-log(beta / 2.0), -xi);

	return (_log_pfrechet(y));
#undef _log_pfrechet
}

double inla_pbgev(double y, double xi)
{
	return (exp(inla_log_pbgev(y, xi)));
}

double inla_inv_pbgev(double p, double xi) 
{
#define _qfrechet(p_) ((pow(-log(p_), -xi) - l1_xi) * sbeta / (l2_xi - l3_xi) + qalpha)
	xi = DMAX(XI_MIN, xi);
	double qalpha = 0, alpha = 0.5, beta = 0.5, sbeta = 1.0; 
	double l1_xi = pow(-log(alpha), -xi);
	double l2_xi = pow(-log(1.0 - beta / 2.0), -xi);
	double l3_xi = pow(-log(beta / 2.0), -xi);

	return (_qfrechet(p));
#undef _qfrechet
}

double link_bgev(int thread_id, double arg, map_arg_tp typ, void *param, double *UNUSED(cov))
{
	Link_param_tp *par = (Link_param_tp *) param;
	double xi_intern = par->bgev_tail[thread_id][0];
	double intercept_intern = par->bgev_intercept[thread_id][0];
	double xi = map_interval(xi_intern, MAP_FORWARD, (void *) par->bgev_tail_interval);
	double intercept;

	// ...as the mapping using this reparameterisation
	if (!ISNAN(intercept_intern)) {
		intercept = inla_inv_pbgev(map_probability(intercept_intern, MAP_FORWARD, NULL), xi);
	} else {
		// this removes the intercept from the model
		intercept = 0.0;
	}

	double low = 0.0, high = 0.0; 

	if (typ == MAP_FORWARD || typ == MAP_DFORWARD) {
		low = inla_inv_pbgev(P_MIN, xi);
		high = inla_inv_pbgev(1.0 - P_MIN, xi);
		arg += intercept;
		arg = TRUNCATE(arg, low, high);
	}
	
	switch (typ) {
	case MAP_FORWARD:
		return (inla_pbgev(arg, xi));

	case MAP_BACKWARD:
		return (inla_inv_pbgev(arg, xi) - intercept);

	case MAP_DFORWARD:
		double h = 1.0e-5;
		return ((inla_pbgev(arg + h, xi) - inla_pbgev(arg - h, xi)) / (2.0 * h));

	case MAP_INCREASING:
		return 1.0;

	default:
		abort();
	}
	return 0.0;
}
