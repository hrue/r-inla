
/* pc-powerlink.c
 * 
 * Copyright (C) 2021-2023 Havard Rue
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

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

#include "inla.h"
#include "pc-powerlink.h"


// have to add the pc-prior for the log(power) later....

// we can use this code to make a generic function for any power-link

double map_inv_powerlink_core(double arg, map_arg_tp typ, void *param, double *intercept_out)
{
	// if 'intercept' is !NULL, just return the contents. its a backdoor avoid duplicating code

#define MAP(_x) log((_x)/(1.0 - (_x)))
#define iMAP(_x) (exp(_x)/(1.0+exp(_x)))
#define diMAP(_x) (exp(_x)/SQR(1.0+exp(_x)))

#define Probit_Pinv(_prob, _power) (-log1p(exp(-log(_prob)/(_power))-2.0))
#define Probit_P(_x, _power) exp(-(_power) * log1p(exp(-(_x))))

	int id = -1;
	GMRFLib_CACHE_SET_ID(id);

	static inla_powerlink_table_tp **table = NULL;
	static int first = 1, x_len = 256;

	int i, j;
	const int debug = 0;
	double **par, intercept_intern, power, power_intern, sd;
	double eps = GMRFLib_eps(0.5);

	par = (double **) param;
	power_intern = *(par[0]);
	intercept_intern = *(par[1]);
	power = map_exp(power_intern, MAP_FORWARD, NULL);

	if (debug) {
		printf("map_inv_powerlink: enter with arg= %g, power= %g, intercept_intern= %g\n", arg, power, intercept_intern);
	}

	if (first) {
#pragma omp critical (Name_9d79a06b70461ee05e66db9259a55b502b777c69)
		if (first) {
			if (1) {
				fprintf(stderr, "map_inv_powerlink: build table with power=%f\n", power);
			}
			table = Calloc(GMRFLib_CACHE_LEN, inla_powerlink_table_tp *);
			for (i = 0; i < GMRFLib_CACHE_LEN; i++) {
				table[i] = Calloc(1, inla_powerlink_table_tp);
				table[i]->power = INLA_REAL_BIG;
				table[i]->cdf = NULL;
				table[i]->icdf = NULL;
			}
			first = 0;
		}
	}

	GMRFLib_CACHE_SET_ID(id);
	if (!ISEQUAL(power, table[id]->power)) {
		double *x, *cdf, yy, p;
		int len;

		if (debug) {
			fprintf(stderr, "map_invsn: build new table for power=%g id=%1d\n", power, id);
		}

		double pp[] = {
			1.0E-3,
			1.0 - 1.0E-3,
			1.0E-4,
			1.0 - 1.0E-4,
			1.0E-5,
			1.0 - 1.0E-5,
			1.0E-6,
			1.0 - 1.0E-6,
			1.0E-7,
			1.0 - 1.0E-7,
			1.0E-8,
			1.0 - 1.0E-8,
			1.0E-9,
			1.0 - 1.0E-9,
			1.0E-10,
			1.0 - 1.0E-10
		};
		int x_len_extra = sizeof(pp) / sizeof(double);
		len = x_len + x_len_extra;

		Calloc_init(2 * len, 2);
		x = Calloc_get(len);
		cdf = Calloc_get(len);

		for (i = 0; i < x_len; i++) {
			p = (i + 0.5) / (double) x_len;
			x[i] = Probit_Pinv(p, power);
			cdf[i] = Probit_P(x[i], power);
		}
		for (j = 0; j < x_len_extra; j++) {
			p = pp[j];
			i = x_len + j;
			x[i] = Probit_Pinv(p, power);
			cdf[i] = Probit_P(x[i], power);
		}

		GMRFLib_qsorts((void *) x, (size_t) len, sizeof(double), (void *) cdf, sizeof(double), (void *) NULL, (size_t) 0, GMRFLib_dcmp);

		/*
		 * moments computed from the CDF, using:
		 *
		 * E(x^k) = k * [ \int_{-inf}^0 x^{k-1} (0-F(x)) dx + \int_0^inf x^{k-1} * (1-F(x)) dx ]
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
		sd = sqrt(mom[2] - SQR(mom[1]));

		for (i = 0; i < len; i++) {
			x[i] = (x[i] - mom[1]) / sd;
		}

		// Remove values in 'cdf' and 'x',that are to close (difference is to small), as this will create issues later on in the
		// interpolation
		GMRFLib_unique_additive2(&len, cdf, x, eps);
		GMRFLib_unique_additive2(&len, x, cdf, eps);

		table[id]->power = power;
		table[id]->xmin = x[0];
		table[id]->xmax = x[len - 1];
		table[id]->pmin = cdf[0];
		table[id]->pmax = cdf[len - 1];
		table[id]->mean = mom[1];
		table[id]->sd = sd;

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

	double intercept;

	if (!ISNAN(intercept_intern)) {
		intercept = GMRFLib_spline_eval(intercept_intern, table[id]->icdf);
	} else {
		intercept = 0.0;
	}

	if (debug) {
		static double intercept_save = 0.0;
		if (intercept != intercept_save) {
			P(intercept);
			P(intercept_intern);
			intercept_save = intercept;
		}
	}

	if (intercept_out) {
		*intercept_out = intercept;
		return 0.0;
	}

	double p, pp;

	switch (typ) {
	case MAP_FORWARD:
	{
		/*
		 * extern = func(local) 
		 */
		arg = TRUNCATE(intercept + arg, table[id]->xmin, table[id]->xmax);
		p = GMRFLib_spline_eval(arg, table[id]->cdf);
		return iMAP(p);
	}

	case MAP_BACKWARD:
	{
		/*
		 * local = func(extern) 
		 */
		arg = TRUNCATE(arg, table[id]->pmin, table[id]->pmax);
		return GMRFLib_spline_eval(MAP(arg), table[id]->icdf) - intercept;
	}

	case MAP_DFORWARD:
	{
		/*
		 * d_extern / d_local 
		 */
		arg = TRUNCATE(intercept + arg, table[id]->xmin, table[id]->xmax);
		p = GMRFLib_spline_eval(arg, table[id]->cdf);
		pp = GMRFLib_spline_eval_deriv(arg, table[id]->cdf);
		return diMAP(p) * pp;
	}

	case MAP_INCREASING:
	{
		/*
		 * return 1.0 if montone increasing and 0.0 otherwise
		 */
		return 1.0;
	}

	default:
		abort();
	}
	abort();

#undef MAP
#undef iMAP
#undef diMAP
	return 0.0;
}
