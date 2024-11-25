
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

#define C_XI 5
#define C_LOW 6
#define C_HIGH 7
#define C_INTERCEPT_INTERN 8
#define C_INTERCEPT 9
#define C_LEN 10
#define P_MIN 1.0E-12
#define XI_MIN 1.0E-8
#define DERIV_H 1.0E-5
#define _log_pfrechet(y_) (-pow(DMAX(0.0, (y_) * l_xi[3] + l_xi[0]), -1.0/xi))
#define _qfrechet(p_) ((pow(-log(p_), -xi) - l_xi[0]) * l_xi[4])

// argument l_xi[3] = { l1_xi, l2_xi, l3_xi}, is defined as 
//
// const double qalpha = 0, alpha = 0.5, beta = 0.5, sbeta = 1.0;
// double l1_xi = pow(-log(alpha), -xi);
// double l2_xi = pow(-log(1.0 - beta / 2.0), -xi);
// double l3_xi = pow(-log(beta / 2.0), -xi);
// double l4_xi = l2_xi - l3_xi
// double l5_xi = 1/l4_xi

double link_gev_bound(double xi, double *l_xi)
{
	return (xi < 0.0 ? (-l_xi[0] / l_xi[3]) : INFINITY);
}

double inla_log_pgev(double y, double xi, double *l_xi)
{
	return (_log_pfrechet(y));
}

double inla_log_pcgev(double y, double xi, double *l_xi)
{
	return (log1p(-exp(_log_pfrechet(y))));
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
	return (_qfrechet(p));
}

double inla_inv_pcgev(double p, double xi, double *l_xi)
{
	return (_qfrechet(1.0 - p));
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
	double xi = map_interval(xi_intern, MAP_FORWARD, (void *) par->bgev_tail_interval);

	typedef struct {
		double store[C_LEN];
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

	xi = DSIGN(xi) * DMAX(XI_MIN, ABS(xi));
	double *c = NULL, *l_xi = NULL;			       // why keep both? they are the same...

	if (!cache[id] || xi != cache[id]->store[C_XI]) {
		if (!cache[id]) {
			cache[id] = Calloc(1, cache_tp);
		}
		// double vec[3] = {-log(0.5), -log(1.0 - 0.5 / 2.0), -log(0.5 / 2.0)};
		double vec[3] = { 0.69314718055994528623, 0.28768207245178090137, 1.38629436111989057245 };
		l_xi = c = cache[id]->store;
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
		l_xi = c = cache[id]->store;
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

	if (type > 0) {
		switch (typ) {
		case MAP_FORWARD:
			return (inla_pgev(arg, xi, l_xi));

		case MAP_BACKWARD:
			return (inla_inv_pgev(arg, xi, l_xi) - intercept);

		case MAP_DFORWARD:
			// same derivative, but we compute the numerical one in log scale, as with
			// f = exp(g) then f'= exp(g) * g' = f * g'
			if (arg + DERIV_H >= c[C_HIGH]) {
				return (inla_pgev(arg, xi, l_xi) *
					((inla_log_pgev(arg, xi, l_xi) - inla_log_pgev(arg - DERIV_H, xi, l_xi)) / DERIV_H));
			} else if (arg - DERIV_H <= c[C_LOW]) {
				return (inla_pgev(arg, xi, l_xi) *
					((inla_log_pgev(arg, xi + DERIV_H, l_xi) - inla_log_pgev(arg, xi, l_xi)) / DERIV_H));
			} else {
				return (inla_pgev(arg, xi, l_xi) *
					((inla_log_pgev(arg, xi + DERIV_H, l_xi) - inla_log_pgev(arg - DERIV_H, xi, l_xi)) / (2.0 * DERIV_H)));
			}

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
			if (arg + DERIV_H > c[C_HIGH]) {
				return (inla_pcgev(arg, xi, l_xi) *
					((inla_log_pcgev(arg, xi, l_xi) - inla_log_pcgev(arg - DERIV_H, xi, l_xi)) / DERIV_H));
			} else if (arg - DERIV_H < c[C_LOW]) {
				return (inla_pcgev(arg, xi, l_xi) *
					((inla_log_pcgev(arg, xi + DERIV_H, l_xi) - inla_log_pcgev(arg, xi, l_xi)) / DERIV_H));
			} else {
				return (inla_pcgev(arg, xi, l_xi) *
					((inla_log_pcgev(arg, xi + DERIV_H, l_xi) - inla_log_pcgev(arg - DERIV_H, xi, l_xi)) / (2.0 * DERIV_H)));
			}

		case MAP_INCREASING:
		default:
			abort();
		}
	}

	return NAN;
}

void link_gev_test(double xi, double intercept)
{
	Link_param_tp *param = Calloc(1, Link_param_tp);

	param->bgev_tail = Calloc(1, double *);
	param->bgev_tail[0] = Calloc(1, double);
	param->bgev_intercept = Calloc(1, double *);
	param->bgev_intercept[0] = Calloc(1, double);

	double intercept_intern = map_probability(intercept, MAP_BACKWARD, NULL);
	param->bgev_intercept[0][0] = intercept_intern;

	param->bgev_tail_interval = Calloc(2, double);
	param->bgev_tail_interval[0] = -0.5;
	param->bgev_tail_interval[1] = 0.5;

	double xi_intern = map_interval(xi, MAP_BACKWARD, (void *) param->bgev_tail_interval);
	param->bgev_tail[0][0] = xi_intern;

	printf("xi %.12g\n", xi);
	printf("xi_intern %.12g\n", xi_intern);
	printf("intercept %.12g\n", intercept);
	printf("intercept_intern %.12g\n", intercept_intern);

	for (double y = -3.0; y <= 3.0; y += 0.01) {
		printf("\ty = %.12g pgev = %.12g %.12g cgev = %.12g %.12g\n", y,
		       link_gev_core(0, y, MAP_FORWARD, (void *) param, 1),
		       link_gev_core(0, y, MAP_DFORWARD, (void *) param, 1),
		       link_gev_core(0, y, MAP_FORWARD, (void *) param, -1), link_gev_core(0, y, MAP_DFORWARD, (void *) param, -1));
	}
}

#undef C_XI
#undef C_LOW
#undef C_HIGH
#undef C_INTERCEPT_INTERN
#undef C_INTERCEPT
#undef C_LEN

#undef P_MIN
#undef XI_MIN
#undef DERIV_H
#undef _log_pfrechet
#undef _qfrechet
