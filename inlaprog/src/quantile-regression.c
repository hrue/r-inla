
/* quantile-regression.c
 * 
 * Copyright (C) 2016-7 Havard Rue
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
#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = HGVERSION;

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

#include "quantile-regression.h"
#include "interpol.h"

double inla_pcontpois(double y, double lambda)
{
	// the cdf for the continous poisson
	return gsl_sf_gamma_inc_Q(y, lambda);
}
double inla_pcontpois_deriv(double y, double lambda)
{
	// the derivative of the cdf for the continous poisson, wrt lambda
	return (-exp((y - 1.0) * log(lambda) - lambda - gsl_sf_lngamma(y)));
}
double inla_qcontpois(double quantile, double alpha, double *initial_guess)
{
	double eta;
	if (initial_guess) {
		eta = log(*initial_guess);
	}
	return (exp(inla_qcontpois_eta(quantile, alpha, (initial_guess ? &eta : NULL))));
}
double inla_qcontpois_eta(double quantile, double alpha, double *initial_guess)
{
#define LOGIT(p) log((p)/(1.0-(p)))
#define DLOGIT(p) (1.0/((p)*(1.0 - (p))))
	int iter_max = 1000, verbose = 0, first_hit = 0;
	double eta_0, eta, max_step = 10, max_step_f = 0.8 * max_step, tol = GMRFLib_eps(0.5);
	double d, f, fd, lambda;

	eta_0 = (initial_guess ? *initial_guess : log(quantile));
	for (int i = 0; i < iter_max; i++) {
		lambda = exp(eta_0);
		f = LOGIT(inla_pcontpois(quantile, lambda)) - LOGIT(alpha);
		fd = DLOGIT((inla_pcontpois(quantile, lambda))) * inla_pcontpois_deriv(quantile, lambda) * lambda;
		d = -DMIN(max_step, DMAX(-max_step_f, f / fd));
		eta = eta_0 = eta_0 + d;
		if (verbose)
			printf("iter=%1d eta=%.6f f=%.8f fd=%.8f d=%.10f\n", i, eta, f, fd, d);
		if (ABS(d) < tol) {
			if (!first_hit)
				return (eta);
			first_hit++;
		}
	}
	GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
#undef LOGIT
#undef DLOGIT
	return eta;
}

GMRFLib_spline_tp **inla_qcontpois_func(double alpha, int num)
{
	/*
	 * return a spline function solving eta for a given log_quantile:
	 * F(exp(log_quantile); exp(eta)) = alpha
	 * alpha is kept fixed
	 */

	int n = 1024, verbose = 0;
	double lq_min = -5.0, lq_max = 10.0, lq_delta = (lq_max - lq_min) / n;
	double *lquantile, *eta;
	GMRFLib_spline_tp **spline;

	eta = Calloc(n, double);
	lquantile = Calloc(n, double);
	for (int i = 0; i < n; i++) {
		lquantile[i] = lq_min + i * lq_delta;
		eta[i] = inla_qcontpois_eta(exp(lquantile[i]), alpha, (i ? &eta[i - 1] : NULL));
		if (verbose)
			printf("i %d lquantile %g eta %g\n", i, lquantile[i], eta[i]);
	}
	spline = Calloc(num, GMRFLib_spline_tp *);
	for (int i = 0; i < num; i++) {
		spline[i] = inla_spline_create(lquantile, eta, n);
	}

	Free(eta);
	Free(lquantile);

	return (spline);
}
