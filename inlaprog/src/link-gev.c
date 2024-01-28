
/* link-gev.c
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
#include "link-gev.h"

#define XI_MIN 1.0E-8
#define P_MIN 1.0E-12

// argument l_xi[3] = { l1_xi, l2_xi, l3_xi}, is defined as 
//
// const double qalpha = 0, alpha = 0.5, beta = 0.5, sbeta = 1.0;
// double l1_xi = pow(-log(alpha), -xi);
// double l2_xi = pow(-log(1.0 - beta / 2.0), -xi);
// double l3_xi = pow(-log(beta / 2.0), -xi);
// double l4_xi = l2_xi - l3_xi
// double l5_xi = 1/l4_xi

double inla_log_pgev(double y, double xi, double *l_xi)
{
#define _log_pfrechet(y_) (-DMAX(0.0, pow((y_) * l_xi[3] + l_xi[0], -1.0/xi)))
	return (_log_pfrechet(y));
#undef _log_pfrechet
}

double inla_log_pcgev(double y, double xi, double *l_xi)
{
#define _log_pfrechet(y_) (-DMAX(0.0, pow((y_) * l_xi[3] + l_xi[0], -1.0/xi)))
	return (log1p(-exp(_log_pfrechet(y))));
#undef _log_pfrechet
}

double inla_pgev(double y, double xi, double *l_xi)
{
	return (exp(inla_log_pgev(y, xi, l_xi)));
}

double inla_pcgev(double y, double xi, double *l_xi)
{
	return (exp(inla_log_pcgev(y, xi, l_xi)));
}

double inla_inv_pgev(double p, double xi, double *l_xi)
{
#define _qfrechet(p_) ((pow(-log(p_), -xi) - l_xi[0]) * l_xi[4])
	return (_qfrechet(p));
#undef _qfrechet
}

double inla_inv_pcgev(double p, double xi, double *l_xi)
{
#define _qfrechet(p_) ((pow(-log(p_), -xi) - l_xi[0]) * l_xi[4])
	return (_qfrechet(1.0 - p));
#undef _qfrechet
}


double link_gev(int thread_id, double arg, map_arg_tp typ, void *param, double *UNUSED(cov))
{
	if (typ == MAP_INCREASING) {
		return 1.0;
	} else {
		return link_gev_core(thread_id, arg, typ, param, 1);
	}
}

double link_cgev(int thread_id, double arg, map_arg_tp typ, void *param, double *UNUSED(cov))
{
	if (typ == MAP_INCREASING) {
		return -1.0;
	} else {
		return link_gev_core(thread_id, arg, typ, param, -1);
	}
}

double link_gev_core(int thread_id, double arg, map_arg_tp typ, void *param, int type)
{
	// type ==1 for gev, type == -1 for cgev

	Link_param_tp *par = (Link_param_tp *) param;
	double xi_intern = par->bgev_tail[thread_id][0];
	double intercept_intern = par->bgev_intercept[thread_id][0];
	double intercept = 0.0;
	double *l_xi = NULL;
	double xi = map_interval(xi_intern, MAP_FORWARD, (void *) par->bgev_tail_interval);
	xi = DMAX(XI_MIN, xi);

#define C_L_XI 0
#define C_XI 5
#define C_LOW 6
#define C_HIGH 7
#define C_INTERCEPT_INTERN 9
#define C_INTERCEPT 9
	typedef struct {
		double store[10];
	} cache_tp;
	static cache_tp **cache = NULL;

	if (!cache) {
#pragma omp critical (Name_05a52e84fdd11002deeb92bd0344f40bc9bac1a6)
		if (!cache) {
			cache = Calloc(GMRFLib_CACHE_LEN(), cache_tp *);
		}
	}

	int id = 0;
	GMRFLib_CACHE_SET_ID(id);

	double *c = NULL;
	if (!cache[id] || xi != cache[id]->store[C_XI]) {
		if (!cache[id]) {
			cache[id] = Calloc(1, cache_tp);
		}
		// double vec[3] = {-log(0.5), -log(1.0 - 0.5 / 2.0), -log(0.5 / 2.0)};
		double vec[3] = { 0.69314718055994528623, 0.28768207245178090137, 1.38629436111989057245 };
		c = cache[id]->store;
		l_xi = c + C_L_XI;
		GMRFLib_powx(3, vec, -xi, l_xi);
		l_xi[3] = l_xi[1] - l_xi[2];
		l_xi[4] = 1.0 / l_xi[3];

		c[C_XI] = xi;
		if (type > 0) {
			c[C_LOW] = inla_inv_pgev(P_MIN, xi, l_xi);
			c[C_HIGH] = inla_inv_pgev(1.0 - P_MIN, xi, l_xi);
		} else {
			c[C_LOW] = inla_inv_pcgev(1.0 - P_MIN, xi, l_xi);
			c[C_HIGH] = inla_inv_pcgev(P_MIN, xi, l_xi);
		}

		c[C_INTERCEPT_INTERN] = NAN;
		c[C_INTERCEPT] = NAN;
	} else {
		c = cache[id]->store;
		l_xi = c + C_L_XI;
	}

	// ...as the mapping using this reparameterisation
	if (!ISNAN(intercept_intern)) {
		if (c[C_INTERCEPT_INTERN] != intercept_intern) {
			if (type > 0) {
				intercept = inla_inv_pgev(map_probability(intercept_intern, MAP_FORWARD, NULL), xi, l_xi);
			} else {
				intercept = inla_inv_pcgev(map_probability(intercept_intern, MAP_FORWARD, NULL), xi, l_xi);
			}
			c[C_INTERCEPT_INTERN] = intercept_intern;
			c[C_INTERCEPT] = intercept;
		} else {
			intercept = c[C_INTERCEPT];
		}
	}

	if (typ == MAP_FORWARD || typ == MAP_DFORWARD) {
		arg += intercept;
		arg = TRUNCATE(arg, c[C_LOW], c[C_HIGH]);
	}

	double h = 1.0e-5;
	if (type > 0) {
		switch (typ) {
		case MAP_FORWARD:
			return (inla_pgev(arg, xi, l_xi));

		case MAP_BACKWARD:
			return (inla_inv_pgev(arg, xi, l_xi) - intercept);

		case MAP_DFORWARD:
			// same derivative, but we compute the numerical one in log scale, as with
			// f = exp(g) then f'= exp(g) * g' = f * g'
			return (inla_pgev(arg, xi, l_xi) * (inla_log_pgev(arg + h, xi, l_xi) - inla_log_pgev(arg - h, xi, l_xi)) / (2.0 * h));

		case MAP_INCREASING:
		default:
			abort();
		}
	} else {
		switch (typ) {
		case MAP_FORWARD:
			return (inla_pcgev(arg, xi, l_xi));

		case MAP_BACKWARD:
			return (inla_inv_pcgev(arg, xi, l_xi) - intercept);

		case MAP_DFORWARD:
			// same derivative, but we compute the numerical one in log scale, as with
			// f = exp(g) then f'= exp(g) * g' = f * g'
			return (inla_pcgev(arg, xi, l_xi) * (inla_log_pcgev(arg + h, xi, l_xi) - inla_log_pcgev(arg - h, xi, l_xi)) / (2.0 * h));

		case MAP_INCREASING:
		default:
			abort();
		}
	}

#undef C_L_XI
#undef C_XI
#undef C_LOW
#undef C_INTERCEPT_INTERN
#undef C_INTERCEPT

	return NAN;
}

#undef XI_MIN
#undef P_MIN
