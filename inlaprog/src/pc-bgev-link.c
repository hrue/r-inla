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

// 

double inla_pbgev_OLD(double y, double xi)
{
#define _log_pfrechet(y_) (-DMAX(0.0, pow(((y_) - qalpha) / (sbeta / (l2_xi - l3_xi)) + l1_xi, -1.0/xi)))
#define _qfrechet(p_) ((pow(-log(p_), -xi) - l1_xi) * sbeta / (l2_xi - l3_xi) + qalpha)
#define _log_pgumbel(y_) (-exp(-(((y_) - qtilde) / stilde)))
#define _pbeta(y_) MATHLIB_FUN(pbeta)(y_, c12, c12, 1, 0)

	if (ISZERO(xi)) {
		xi = 1.0e-8;
	}

	double qalpha = 0, alpha = 0.5, beta = 0.5, sbeta = 1.0, pa = 0.1, pb = 0.2, c12 = 5.0;
	double l1 = log(-log(alpha));
	double l2 = log(-log(1.0 - beta / 2.0));
	double l3 = log(-log(beta / 2.0));
	double l4 = log(-log(pa));
	double l5 = log(-log(pb));

	double l1_xi = pow(-log(alpha), -xi);
	double l2_xi = pow(-log(1.0 - beta / 2.0), -xi);
	double l3_xi = pow(-log(beta / 2.0), -xi);

	double a = _qfrechet(pa);
	double b = _qfrechet(pb);
	double px = _pbeta(y - a / (b - a));

	double qtilde = a - (b - a) * (l1 - l4) / (l4 - l5);
	double stilde = (b - a) * (l3 - l2) / (l4 - l5);
	double part1 = _log_pfrechet(y);
	double part2 = _log_pgumbel(y);
	double value = exp(px * part1 + (1 - px) * part2);

#undef _log_pfrechet
#undef _qfrechet
#undef _log_pgumbel
#undef _pbeta

	return (value);
}

double link_bgev_OLD(int thread_id, double arg, map_arg_tp typ, void *param, double *UNUSED(cov))
{
#define MAP(_x) log((_x)/(1.0 - (_x)))
#define iMAP(_x) (exp(_x)/(1.0+exp(_x)))
#define diMAP(_x) (exp(_x)/SQR(1.0+exp(_x)))

	static inla_link_bgev_table_tp **table = NULL;
	static char first = 1;

	int i, j, id = 0;
	const int debug = 1;
	double dx = 0.02, range = 10.0, p, pp; 
	double intercept;

	Link_param_tp *par = (Link_param_tp *) param;
	double xi_intern = par->bgev_tail[thread_id][0];
	double intercept_intern = par->bgev_intercept[thread_id][0];
	double xi = map_interval(xi_intern, MAP_FORWARD, (void *) par->bgev_tail_interval);

	if (first) {
#pragma omp critical (Name_f35fc78992e1ea9c433855032e75213c669a9284)
		if (first) {
			if (debug) {
				fprintf(stderr, "map_invbgev: build table\n");
			}
			table = Calloc(GMRFLib_CACHE_LEN(), inla_link_bgev_table_tp *);
			for (i = 0; i < GMRFLib_CACHE_LEN(); i++) {
				table[i] = Calloc(1, inla_link_bgev_table_tp);
				table[i]->xi = -INLA_REAL_BIG;
				table[i]->cdf = NULL;
				table[i]->icdf = NULL;
			}
			first = 0;
		}
	}

	GMRFLib_CACHE_SET_ID(id);

	if (!ISEQUAL(xi, table[id]->xi)) {

		if (debug) {
			printf("link_bgev: enter with arg= %g, xi= %g\n", arg, xi);
		}

		int len = (int) (2.0 * range / dx + 0.5) + 1, llen = 0;
		double *work, *x, *y, *yy, nc = 0.0, xx;

		if (debug) {
			fprintf(stderr, "map_invbgev: build new table for xi=%g id=%1d\n", xi, id);
		}

		work = Calloc(3 * len, double);
		x = work;
		y = work + len;
		yy = work + 2 * len;

		for (xx = -range, i = 0, llen = 0; xx <= range; xx += dx, i++) {
			x[i] = xx;
			y[i] = inla_pbgev(x[i], xi);
			llen++;
		}
		len = llen;

		for (i = j = 0; i < len; i++) {
			if (!ISNAN(y[i]) && !ISINF(y[i])) {
				x[j] = x[i];
				y[j] = y[i];
				j++;
			}
		}
		len = j;

		for (i = 0, nc = 0.0; i < len; i++) {
			nc += y[i];
		}

		nc = 1.0 / nc;
		for (i = 0; i < len; i++) {
			y[i] *= nc;
		}

		yy[0] = y[0];
		for (i = 1; i < len; i++) {
			yy[i] = yy[i - 1] + (y[i] + y[i - 1]) / 2.0;
		}
		for (i = 0; i < len; i++) {
			y[i] = MAP(yy[i]);
		}

		for (i = j = 0; i < len; i++) {
			if (!ISNAN(y[i]) && !ISINF(y[i])) {
				x[j] = x[i];
				y[j] = y[i];
				j++;
			}
		}
		len = j;
		// Remove values in 'y' that are to close (difference is to small)
		GMRFLib_unique_additive2(&len, y, x, (GSL_SQRT_DBL_EPSILON * GSL_ROOT4_DBL_EPSILON));

		table[id]->xi = xi;
		table[id]->xmin = x[0];
		table[id]->xmax = x[len - 1];
		table[id]->pmin = iMAP(y[0]);
		table[id]->pmax = iMAP(y[len - 1]);

		GMRFLib_spline_free(table[id]->cdf);	       /* ok if NULL */
		GMRFLib_spline_free(table[id]->icdf);	       /* ok if NULL */

		table[id]->cdf = GMRFLib_spline_create(x, y, len);
		table[id]->icdf = GMRFLib_spline_create(y, x, len);
		Free(work);
	}

	// ...as the mapping using this reparameterisation
	if (!ISNAN(intercept_intern)) {
		intercept = GMRFLib_spline_eval(intercept_intern, table[id]->icdf);
	} else {
		// this removes the intercept from the model
		intercept = 0.0;
	}

	switch (typ) {
	case MAP_FORWARD:
		/*
		 * extern = func(local) 
		 */
		arg = TRUNCATE(intercept + arg, table[id]->xmin, table[id]->xmax);
		p = GMRFLib_spline_eval(arg, table[id]->cdf);
		return iMAP(p);

	case MAP_BACKWARD:
		/*
		 * local = func(extern) 
		 */
		arg = TRUNCATE(arg, table[id]->pmin, table[id]->pmax);
		return GMRFLib_spline_eval(MAP(arg), table[id]->icdf) - intercept;

	case MAP_DFORWARD:
		/*
		 * d_extern / d_local 
		 */
		arg = TRUNCATE(intercept + arg, table[id]->xmin, table[id]->xmax);
		p = GMRFLib_spline_eval(arg, table[id]->cdf);
		pp = GMRFLib_spline_eval_deriv(arg, table[id]->cdf);
		return diMAP(p) * pp;

	case MAP_INCREASING:
		/*
		 * return 1.0 if montone increasing and 0.0 otherwise
		 */
		return 1.0;

	default:
		abort();
	}
	abort();

#undef MAP
#undef iMAP
#undef diMAP
	return 0.0;
}
