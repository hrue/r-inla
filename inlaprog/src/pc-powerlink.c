
/* pc-powerlink.c
 * 
 * Copyright (C) 2021 Havard Rue
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
#ifndef GITCOMMIT
#define GITCOMMIT
#endif

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

static const char GitID[] = "file: " __FILE__ "  " GITCOMMIT;

#include "inla.h"
#include "pc-powerlink.h"

double map_inv_powerlink_core(double arg, map_arg_tp typ, void *param, double *intercept_out)
{
	// if 'intercept' is !NULL, just return the contents. its a backdoor avoid duplicating code

#define MAP(_x) log((_x)/(1.0 - (_x)))
#define iMAP(_x) (exp(_x)/(1.0+exp(_x)))
#define diMAP(_x) (exp(_x)/SQR(1.0+exp(_x)))

	int id;
	GMRFLib_CACHE_SET_ID(id);
	
	static inla_powerlink_table_tp **table = NULL;
	static char first = 1;

	int i, debug = 1;
	double alpha;
	double **par, intercept_local, intercept_intern, power, power_intern;

	par = (double **) param;
	power_intern = *(par[0]);
	intercept_intern = *(par[1]);

	// power is the exponent
	// intercept is the quantile-level
	power = map_exp(power_intern, MAP_FORWARD, NULL);

	if (debug) {
		printf("map_inv_powerlink: enter with arg= %g, power= %g, intercept_intern= %g\n", arg, power, intercept_intern);
	}

	static double *lcdf_ref = NULL;
	static double *x_ref = NULL;
	static int x_len = 0;

	if (first) {
#pragma omp critical
		if (first) {
			if (debug) {
				fprintf(stderr, "map_inv_powerlink: build table\n");
			}
			table = Calloc(GMRFLib_CACHE_LEN, inla_powerlink_table_tp *);
			for (i = 0; i < GMRFLib_CACHE_LEN; i++) {
				table[i] = Calloc(1, inla_powerlink_table_tp);
				table[i]->power = INLA_REAL_BIG;
				table[i]->cdf = NULL;
				table[i]->icdf = NULL;
			}

			int len_add = 7;
			int len_extra = 2 * len_add + 1;
			x_len = 1000;
			lcdf_ref = Calloc(2 * (x_len + len_extra), double);
			x_ref = lcdf_ref + x_len + len_extra;

			// layout quantiles and some extra points
			for (i = 0; i < x_len; i++) {
				double p = (i + 0.5) / (double) x_len;
				x_ref[i] = gsl_cdf_ugaussian_Pinv(p);
			}
			for (i = -len_add; i <= len_add; i++) {
				x_ref[x_len++] = i;
			}
			qsort((void *) x_ref, (size_t) x_len, sizeof(double), GMRFLib_dcmp);
			for (i = 0; i < x_len; i++) {
				lcdf_ref[i] = inla_log_Phi(x_ref[i]);
			}
			first = 0;
		}
	}

	GMRFLib_CACHE_SET_ID(id);
	if (!ISEQUAL(power, table[id]->power)) {
		double *x, *cdf, yy;
		int len;

		if (debug) {
			fprintf(stderr, "map_invsn: build new table for alpha=%g id=%1d\n", alpha, id);
		}

		Calloc_init(2 * x_len);
		x = Calloc_get(x_len);
		cdf = Calloc_get(x_len);

		for (i = 0; i < x_len; i++) {
			x[i] = x_ref[i];
			cdf[i] = exp(power * lcdf_ref[i]);
		}
		len = x_len;

		/*
		  moments computed from the CDF, using:

		  E(x^k) = k * [ \int_{-inf}^0 x^{k-1} (0-F(x)) dx + \int_0^inf x^{k-1} * (1-F(x)) dx ]
		*/
		double mom[3] = { 0.0, 0.0, 0.0 }, w;
		for (i = 1; i < len - 1; i++) {
			// I remove the ()/2.0 for w, and rather correct at the end
			w = x[i + 1] - x[i - 1];
			yy = (x[i] < 0.0 ? -cdf[i] : 1.0 - cdf[i]);
			mom[1] += w * yy;
			mom[2] += w * x[i] * yy;
		}
		mom[1] /= 2.0;

		double sd = sqrt(mom[2] - SQR(mom[1]));

		P(mom[1]);
		P(sd);

		// Remove values in 'cdf' that are to close (difference is to small), as this will create issues later on in the interpolation
		GMRFLib_unique_additive2(&len, cdf, x, GMRFLib_eps(0.75));

		table[id]->power = power;
		table[id]->xmin = x[0];
		table[id]->xmax = x[len - 1];
		table[id]->pmin = cdf[0];
		table[id]->pmax = cdf[len-1];
		// transform before spline'ing
		for (i = 0; i < len; i++) {
			cdf[i] = MAP(cdf[i]);
		}
		GMRFLib_spline_free(table[id]->cdf);	       /* ok if NULL */
		GMRFLib_spline_free(table[id]->icdf);	       /* ok if NULL */

		table[id]->cdf = GMRFLib_spline_create(x, cdf, len);
		table[id]->icdf = GMRFLib_spline_create(cdf, x, len);

		Calloc_free();
	}

	if (!ISNAN(intercept_intern)) {
		intercept = GMRFLib_spline_eval(intercept_intern, table[id]->icdf);
	} else {
		intercept = 0.0;
	}

	if (debug) {
		printf("... intercept_intern= %g intercept_alpha=%g intercept= %g\n",
		       intercept_intern, iMAP(intercept_intern), intercept);
	}

	if (intercept_out) {
		*intercept_out = intercept;
		return 0.0;
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
