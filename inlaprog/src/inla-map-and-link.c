
/* inla-map-and-links.c
 * 
 * Copyright (C) 2007-2024 Havard Rue
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *  * This program is distributed in the hope that it will be useful, but
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

double map_one(double UNUSED(arg), map_arg_tp typ, void *UNUSED(param))
{
	switch (typ) {
	case MAP_FORWARD:
		return 1.0;
	case MAP_BACKWARD:
		return 1.0;
	case MAP_DFORWARD:
		return 0.0;
	case MAP_INCREASING:
		return 1.0;
	default:
		abort();
	}
	return 0.0;
}

double map_identity(double arg, map_arg_tp typ, void *UNUSED(param))
{
	/*
	 * the idenity map-function
	 */
	switch (typ) {
	case MAP_FORWARD:
		/*
		 * extern = func(local) 
		 */
		return arg;
	case MAP_BACKWARD:
		/*
		 * local = func(extern) 
		 */
		return arg;
	case MAP_DFORWARD:
		/*
		 * d_extern / d_local 
		 */
		return 1.0;
	case MAP_INCREASING:
		/*
		 * return 1.0 if montone increasing and 0.0 otherwise 
		 */
		return 1.0;
	default:
		abort();
	}
	abort();
	return 0.0;
}

double map_inverse(double arg, map_arg_tp typ, void *UNUSED(param))
{
	/*
	 * the inverse map-function, assuming > 0
	 */
	switch (typ) {
	case MAP_FORWARD:
		/*
		 * extern = func(local) 
		 */
		return 1.0 / arg;
	case MAP_BACKWARD:
		/*
		 * local = func(extern) 
		 */
		return 1.0 / arg;
	case MAP_DFORWARD:
		/*
		 * d_extern / d_local 
		 */
		return -1.0 / SQR(arg);
	case MAP_INCREASING:
		/*
		 * return 1.0 if montone increasing and 0.0 otherwise. assuming positive...
		 */
		return 0.0;
	default:
		abort();
	}
	abort();
	return 0.0;
}

double map_identity_scale(double arg, map_arg_tp typ, void *param)
{
	/*
	 * the idenity map-function
	 */
	double scale = (param ? *((double *) param) : 1.0);

	switch (typ) {
	case MAP_FORWARD:
		/*
		 * extern = func(local) 
		 */
		return arg * scale;
	case MAP_BACKWARD:
		/*
		 * local = func(extern) 
		 */
		return arg / scale;
	case MAP_DFORWARD:
		/*
		 * d_extern / d_local 
		 */
		return scale;
	case MAP_INCREASING:
		/*
		 * return 1.0 if montone increasing and 0.0 otherwise 
		 */
		return (scale > 0 ? 1.0 : 0.0);
	default:
		abort();
	}
	abort();
	return 0.0;
}

double map_exp(double arg, map_arg_tp typ, void *UNUSED(param))
{
	/*
	 * the exp-map-function
	 */
	switch (typ) {
	case MAP_FORWARD:
		/*
		 * extern = func(local) 
		 */
		return exp(arg);
	case MAP_BACKWARD:
		/*
		 * local = func(extern) 
		 */
		return log(arg);
	case MAP_DFORWARD:
		/*
		 * d_extern / d_local 
		 */
		return exp(arg);
	case MAP_INCREASING:
		/*
		 * return 1.0 if montone increasing and 0.0 otherwise 
		 */
		return 1.0;
	default:
		abort();
	}
	abort();
	return 0.0;
}

double map_exp_scale(double arg, map_arg_tp typ, void *param)
{
	/*
	 * the exp-map-function with scaling
	 */
	double scale = *((double *) param);
	switch (typ) {
	case MAP_FORWARD:
		/*
		 * extern = func(local) 
		 */
		return exp(scale * arg);
	case MAP_BACKWARD:
		/*
		 * local = func(extern) 
		 */
		return log(arg) / scale;
	case MAP_DFORWARD:
		/*
		 * d_extern / d_local 
		 */
		return exp(scale * arg) * scale;
	case MAP_INCREASING:
		/*
		 * return 1.0 if montone increasing and 0.0 otherwise 
		 */
		return (scale > 0.0 ? 1.0 : 0.0);
	default:
		abort();
	}
	abort();
	return 0.0;
}

double map_negexp(double arg, map_arg_tp typ, void *UNUSED(param))
{
	/*
	 * the negexp-map-function
	 */
	switch (typ) {
	case MAP_FORWARD:
		/*
		 * extern = func(local) 
		 */
		return exp(-arg);
	case MAP_BACKWARD:
		/*
		 * local = func(extern) 
		 */
		return -log(arg);
	case MAP_DFORWARD:
		/*
		 * d_extern / d_local 
		 */
		return -exp(-arg);
	case MAP_INCREASING:
		/*
		 * return 1.0 if montone increasing and 0.0 otherwise 
		 */
		return 0.0;
	default:
		abort();
	}
	abort();
	return 0.0;
}

double map_exp_scale2(double arg, map_arg_tp typ, void *param)
{
	/*
	 * the exp-map-function with a scale
	 */
	double *p = (double *) param;
	double a = p[0], b = p[1];

	switch (typ) {
	case MAP_FORWARD:
		/*
		 * extern = func(local) 
		 */
		return a * exp(b * arg);
	case MAP_BACKWARD:
		/*
		 * local = func(extern) 
		 */
		return log(arg / a) / b;
	case MAP_DFORWARD:
		/*
		 * d_extern / d_local 
		 */
		return a * b * exp(b * arg);
	case MAP_INCREASING:
		/*
		 * return 1.0 if montone increasing and 0.0 otherwise 
		 */
		return (map_exp_scale2(1.0, MAP_FORWARD, param) > map_exp_scale2(0.0, MAP_FORWARD, param) ? 1 : -1);
	default:
		abort();
	}
	abort();
	return 0.0;
}

double map_invrobit(double arg, map_arg_tp typ, void *param)
{
	/*
	 * the inverse of the robit link: cdf of student t (scaled to have variance 1)
	 */
	double df_intern = *((double *) param);
	double df = map_dof(df_intern, MAP_FORWARD, NULL);
	double scale = sqrt(df / (df - 2.0));

	switch (typ) {
	case MAP_FORWARD:
		/*
		 * extern = func(local) 
		 */
		return MATHLIB_FUN(pt) (arg / scale, df, 1, 0);
	case MAP_BACKWARD:
		/*
		 * local = func(extern) 
		 */
		return MATHLIB_FUN(qt) (arg, df, 1, 0) * scale;
	case MAP_DFORWARD:
		/*
		 * d_extern / d_local 
		 */
		return MATHLIB_FUN(dt) (arg / scale, df, 0) / scale;
	case MAP_INCREASING:
		/*
		 * return 1.0 if montone increasing and 0.0 otherwise
		 */
		return 1.0;
	default:
		abort();
	}
	abort();
	return 0.0;
}
double inla_get_sn_param(inla_sn_arg_tp *output, double **param)
{
	// param = *(skew_intern, intercept_intern)
	return (map_invsn_core(0.0, MAP_FORWARD, (void *) param, output));
}

double map_invsn(double arg, map_arg_tp typ, void *param)
{
	return (map_invsn_core(arg, typ, param, NULL));
}

double map_invsn_core(double arg, map_arg_tp typ, void *param, inla_sn_arg_tp *output)
{
	// if 'output' is !NULL, just return the contents. its a backdoor avoid duplicating code

#define MAP(_x) log((_x)/(1.0 - (_x)))
#define iMAP(_x) (exp(_x)/(1.0+exp(_x)))
#define diMAP(_x) (exp(_x)/SQR(1.0+exp(_x)))

	/*
	 * the inverse link is the cdf of the sn, scaled to have E()=1 and Var()=1
	 */
	static inla_sn_table_tp **table = NULL;
	static char first = 1;

	int i, j, id = 0;
	const int debug = 0;
	double alpha, dx = 0.02, range = 10.0, p, pp, omega, delta, xi, skew, skew_intern, skew_max = GMRFLib_SN_SKEWMAX;
	double **par = NULL, intercept, intercept_intern, intercept_alpha;

	par = (double **) param;
	assert(par);
	assert(par[0]);
	assert(par[1]);
	skew_intern = *(par[0]);
	intercept_intern = *(par[1]);

	// parameters are SKEW and INTERCEPT
	skew = map_phi(skew_intern, MAP_FORWARD, (void *) &skew_max);
	if (!ISNAN(intercept_intern)) {
		intercept_alpha = map_probability(intercept_intern, MAP_FORWARD, NULL);
	} else {
		intercept_alpha = NAN;
	}

	if (debug) {
		printf("map_invsn: enter with arg= %g, skew= %g, intercept_alpha= %g\n", arg, skew, intercept_alpha);
	}

	if (first) {
#pragma omp critical (Name_375605f5a0a02853485931b7541aa2de0692bb42)
		if (first) {
			if (debug) {
				fprintf(stderr, "map_invsn: build table\n");
			}
			table = Calloc(GMRFLib_CACHE_LEN(), inla_sn_table_tp *);
			for (i = 0; i < GMRFLib_CACHE_LEN(); i++) {
				table[i] = Calloc(1, inla_sn_table_tp);
				table[i]->alpha = INLA_REAL_BIG;
				table[i]->cdf = NULL;
				table[i]->icdf = NULL;
			}
			first = 0;
		}
	}

	alpha = inla_pc_sn_skew2alpha(skew);
	delta = alpha / sqrt(1.0 + SQR(alpha));
	omega = 1.0 / sqrt(1.0 - 2.0 * SQR(delta) / M_PI);
	xi = -omega * delta * sqrt(2.0 / M_PI);
	GMRFLib_CACHE_SET_ID(id);
	if (debug) {
		printf("...this gives alpha= %g, delta= %g, omega= %g, xi= %g\n", alpha, delta, omega, xi);
	}

	if (!ISEQUAL(alpha, table[id]->alpha)) {
		int len = (int) (2.0 * range / dx + 0.5) + 1, llen = 0;
		double *work = NULL, *x = NULL, *y = NULL, *yy = NULL, nc = 0.0, xx;

		if (debug) {
			fprintf(stderr, "map_invsn: build new table for alpha=%g id=%1d\n", alpha, id);
		}

		work = Calloc(3 * len, double);
		x = work;
		y = work + len;
		yy = work + 2 * len;

		for (xx = -range, i = 0, llen = 0; xx <= range; xx += dx, i++) {
			x[i] = xx;
			if (alpha != 0.0) {
				double z = (xx - xi) / omega;
				y[i] = 2.0 / omega * MATHLIB_FUN(dnorm) (z, 0.0, 1.0, 0) * MATHLIB_FUN(pnorm) (alpha * z, 0.0, 1.0, 1, 0);
			} else {
				y[i] = MATHLIB_FUN(dnorm) (xx, 0.0, 1.0, 0);
			}
			llen++;
		}
		len = llen;

		if (debug) {
			// check that we have done it right...
			double mom[4] = { 0, 0, 0, 0 }, negative = 0;
			for (i = 0; i < len; i++) {
				mom[0] += y[i];
				mom[1] += y[i] * x[i];
				mom[2] += y[i] * SQR(x[i]);
				mom[3] += y[i] * POW3(x[i]);
				negative += y[i] * (x[i] < 0);
			}
			mom[1] /= mom[0];
			mom[2] /= mom[0];
			mom[3] /= mom[0];
			negative /= mom[0];
			printf("map_invsn: alpha= %.6g, negative= %.6g moments: %.6g %.6g %.6g\n", alpha, negative, mom[1], mom[2], mom[3]);
		}

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

		table[id]->alpha = alpha;
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

	if (debug) {
		printf("... intercept_alpha= %g intercept= %g\n", intercept_alpha, intercept);
	}

	if (output) {
		output->omega = omega;
		output->xi = xi;
		output->intercept = intercept;
		output->alpha = alpha;
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

double map_invprobit(double arg, map_arg_tp typ, void *UNUSED(param))
{
	/*
	 * the inverse probit function
	 */
	switch (typ) {
	case MAP_FORWARD:
		/*
		 * extern = func(local) 
		 */
		return GMRFLib_cdfnorm(arg);
	case MAP_BACKWARD:
		/*
		 * local = func(extern) 
		 */
		return GMRFLib_cdfnorm_inv(arg);
	case MAP_DFORWARD:
		/*
		 * d_extern / d_local 
		 */
		return gsl_ran_ugaussian_pdf(arg);

	case MAP_INCREASING:
		/*
		 * return 1.0 if montone increasing and 0.0 otherwise 
		 */
		return 1.0;
	default:
		abort();
	}
	abort();
	return 0.0;
}

double map_invloglog(double arg, map_arg_tp typ, void *UNUSED(param))
{
	/*
	 * the inverse loglog function
	 */
	switch (typ) {
	case MAP_FORWARD:
		/*
		 * extern = func(local) 
		 */
		return exp(-exp(-arg));
	case MAP_BACKWARD:
		/*
		 * local = func(extern) 
		 */
		return -log(-log(arg));
	case MAP_DFORWARD:
		/*
		 * d_extern / d_local 
		 */
		return exp(-arg - exp(-arg));
	case MAP_INCREASING:
		/*
		 * return 1.0 if montone increasing and 0.0 otherwise 
		 */
		return 1.0;
	default:
		abort();
	}
	abort();
	return 0.0;
}

double map_invcauchit(double arg, map_arg_tp typ, void *UNUSED(param))
{
	/*
	 * the inverse cauchit function
	 */
	switch (typ) {
	case MAP_FORWARD:
		/*
		 * extern = func(local) 
		 */
		return M_1_PI * atan(arg) + 0.5;
	case MAP_BACKWARD:
		/*
		 * local = func(extern) 
		 */
		return tan(M_PI * (arg - 0.5));
	case MAP_DFORWARD:
		/*
		 * d_extern / d_local 
		 */
		return M_1_PI / (SQR(arg) + 1.0);
	case MAP_INCREASING:
		/*
		 * return 1.0 if montone increasing and 0.0 otherwise 
		 */
		return 1.0;
	default:
		abort();
	}
	abort();
	return 0.0;
}

double link_log_invcloglog(double x)
{
	// log(invcloglog(x))
	return log1p(-exp(-exp(x)));
}

double link_log_1m_invcloglog(double x)
{
	// log(1 - invcloglog(x))
	return (-exp(x));
}

void link_log_invcloglog2(double x, double *r1, double *r2)
{
	// return the result of link_log_invcloglog(x) in r1, link_log_1m_invcloglog(x) in r2
	double v = -exp(x);
	*r2 = v;
	*r1 = log1p(-exp(v));
}

double map_invcloglog(double arg, map_arg_tp typ, void *UNUSED(param))
{
	/*
	 * the inverse cloglog function
	 */
	switch (typ) {
	case MAP_FORWARD:
		/*
		 * extern = func(local) 
		 */
		return ONE_mexp(-exp(arg));
	case MAP_BACKWARD:
		/*
		 * local = func(extern) 
		 */
		return log(-LOG_1mp(arg));
	case MAP_DFORWARD:
		/*
		 * d_extern / d_local 
		 */
		return exp(arg - exp(arg));
	case MAP_INCREASING:
		/*
		 * return 1.0 if montone increasing and 0.0 otherwise 
		 */
		return 1.0;
	default:
		abort();
	}
	abort();
	return 0.0;
}

double link_log_invccloglog(double x)
{
	// log(invccloglog(x))
	return (-exp(-x));
}

double link_log_1m_invccloglog(double x)
{
	// log(1 - invccloglog(x))
	return log1p(-exp(-exp(-x)));
}

void link_log_invccloglog2(double x, double *r1, double *r2)
{
	// return the result of link_log_invccloglog(x) in r1, link_log_1m_invccloglog(x) in r2
	double v = -exp(-x);
	*r1 = v;
	*r2 = log1p(-exp(v));
}

double map_invccloglog(double arg, map_arg_tp typ, void *UNUSED(param))
{
	/*
	 * the inverse complement cloglog function
	 */
	switch (typ) {
	case MAP_FORWARD:
		/*
		 * extern = func(local) 
		 */
		return exp(-exp(-arg));
	case MAP_BACKWARD:
		/*
		 * local = func(extern) 
		 */
		return -log(-log(arg));
	case MAP_DFORWARD:
		/*
		 * d_extern / d_local 
		 */
		return exp(-arg - exp(-arg));
	case MAP_INCREASING:
		/*
		 * return 1.0 if montone increasing and 0.0 otherwise 
		 */
		return 1.0;
	default:
		abort();
	}
	abort();
	return 0.0;
}

double map_beta(double x, map_arg_tp typ, void *param)
{
	/*
	 * the map for the beta parameter, which can have a lower and upper range as well. If range.low=range.high, then its
	 * interpreted as range.low = -INF and range.high = INF, ie the mapping is the identity. If range.high = INF and
	 * range.low != INF, then the mapping is range.low + exp(...).
	 */

	double *range = (double *) param;

	if (param == NULL || ISEQUAL(range[0], range[1])) {
		return map_identity(x, typ, param);
	}

	if (ISINF(range[1]) && !ISINF(range[0])) {
		switch (typ) {
		case MAP_FORWARD:
			/*
			 * extern = func(local) 
			 */
			return range[0] + exp(x);
		case MAP_BACKWARD:
			/*
			 * local = func(extern) 
			 */
			return log(x - range[0]);
		case MAP_DFORWARD:
			/*
			 * d_extern / d_local 
			 */
			return exp(x);
		case MAP_INCREASING:
			/*
			 * return 1.0 if montone increasing and 0.0 otherwise 
			 */
			return 1.0;
		default:
			abort();
		}
	} else if (ISINF(range[0]) && !ISINF(range[1])) {
		FIXME("the case: ISINF(range[0]) && !ISINF(range[1]), is not yet implemented.");
		exit(EXIT_FAILURE);
	} else {
		/*
		 * Then the mapping is
		 * 
		 * range[0] + exp(x)/(1 + exp(x)) * (range[1] - range[0]) 
		 */

		double d = range[1] - range[0], xx;

		switch (typ) {
		case MAP_FORWARD:
			return range[0] + d * exp(x) / (1.0 + exp(x));
		case MAP_BACKWARD:
			xx = (x - range[0]) / d;
			return log(xx / (1.0 - xx));
		case MAP_DFORWARD:
			xx = exp(x);
			return d * xx / SQR(1.0 + xx);
		case MAP_INCREASING:
			return 1.0;
		default:
			GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
		}
	}
	return 0.0;
}

double map_1exp(double arg, map_arg_tp typ, void *UNUSED(param))
{
	/*
	 * the 1/exp-map-function
	 */
	switch (typ) {
	case MAP_FORWARD:
		/*
		 * extern = func(local) 
		 */
		return exp(-arg);
	case MAP_BACKWARD:
		/*
		 * local = func(extern) 
		 */
		return -log(arg);
	case MAP_DFORWARD:
		/*
		 * d_extern / d_local 
		 */
		return -exp(-arg);
	case MAP_INCREASING:
		/*
		 * return 1.0 if montone increasing and 0.0 otherwise 
		 */
		return 0.0;
	default:
		abort();
	}
	abort();
	return 0.0;
}

double map_sqrt1exp(double arg, map_arg_tp typ, void *UNUSED(param))
{
	/*
	 * the sqrt(1/exp) map
	 */
	switch (typ) {
	case MAP_FORWARD:
		/*
		 * extern = func(local) 
		 */
		return exp(-0.5 * arg);
	case MAP_BACKWARD:
		/*
		 * local = func(extern) 
		 */
		return -2.0 * log(arg);
	case MAP_DFORWARD:
		/*
		 * d_extern / d_local 
		 */
		return -0.5 * exp(-0.5 * arg);
	case MAP_INCREASING:
		/*
		 * return 1.0 if montone increasing and 0.0 otherwise 
		 */
		return 0.0;
	default:
		abort();
	}
	abort();
	return 0.0;
}

double map_dof(double arg, map_arg_tp typ, void *UNUSED(param))
{
	/*
	 * the map-function for the degrees of freedom for the student-t 
	 */
	switch (typ) {
	case MAP_FORWARD:
		/*
		 * extern = func(local) 
		 */
		return 2.0 + exp(arg);
	case MAP_BACKWARD:
		/*
		 * local = func(extern) 
		 */
		return log(arg - 2.0);
	case MAP_DFORWARD:
		/*
		 * d_extern / d_local 
		 */
		return exp(arg);
	case MAP_INCREASING:
		/*
		 * return 1.0 if montone increasing and 0.0 otherwise 
		 */
		return 1.0;
	default:
		abort();
	}
	abort();
	return 0.0;
}

double map_phi(double arg, map_arg_tp typ, void *param)
{
	/*
	 * the map-function for the lag-1 correlation in the AR(1) model. The
	 * extra argument, if present, is the range (default = 1)
	 */
	double xx;
	double range = (param ? *((double *) param) : 1.0);

	switch (typ) {
	case MAP_FORWARD:
		/*
		 * extern = func(local) 
		 */
		xx = exp(arg);
		return (range * (2.0 * (xx / (1.0 + xx)) - 1.0));
	case MAP_BACKWARD:
		/*
		 * local = func(extern) 
		 */
		return log((1.0 + arg / range) / (1.0 - arg / range));
	case MAP_DFORWARD:
		/*
		 * d_extern / d_local 
		 */
		xx = exp(arg);
		return (range * 2.0 * xx / SQR(1.0 + xx));
	case MAP_INCREASING:
		/*
		 * return 1.0 if montone increasing and 0.0 otherwise 
		 */
		return (range >= 0 ? 1.0 : -1.0);
	default:
		abort();
	}
	abort();
	return 0.0;
}

double map_rho(double arg, map_arg_tp typ, void *UNUSED(param))
{
	/*
	 * the map-function for the lag-1 correlation in the AR(1) model. The
	 * extra argument, if present, is the range (default = 1)
	 */
	double xx;

	switch (typ) {
	case MAP_FORWARD:
		/*
		 * extern = func(local) 
		 */
		xx = exp(-arg);
		return (2.0 / (1.0 + xx) - 1.0);
	case MAP_BACKWARD:
		/*
		 * local = func(extern) 
		 */
		return log((1.0 + arg) / (1.0 - arg));
	case MAP_DFORWARD:
		/*
		 * d_extern / d_local 
		 */
		xx = exp(arg);
		return (2.0 * xx / SQR(1.0 + xx));
	case MAP_INCREASING:
		/*
		 * return 1.0 if montone increasing and 0.0 otherwise 
		 */
		return 1.0;
	default:
		abort();
	}
	abort();
	return 0.0;
}

double map_precision(double arg, map_arg_tp typ, void *param)
{
	/*
	 * the map-function for the precision variables 
	 */
	return map_exp(arg, typ, param);
}

double map_range(double arg, map_arg_tp typ, void *param)
{
	/*
	 * the map-function for the range
	 */
	return map_exp(arg, typ, param);
}

double map_alpha_weibull(double arg, map_arg_tp typ, void *UNUSED(param))
{
	/*
	 * the map-function for the range
	 */
	double scale = INLA_WEIBULL_ALPHA_SCALE;
	return map_exp_scale(arg, typ, (void *) &scale);
}

double map_alpha_gompertz(double arg, map_arg_tp typ, void *UNUSED(param))
{
	/*
	 * the map-function for the range
	 */
	double scale = INLA_GOMPERTZ_ALPHA_SCALE;
	return map_exp_scale(arg, typ, (void *) &scale);
}

double map_prec_qkumar(double arg, map_arg_tp typ, void *UNUSED(param))
{
	/*
	 * the map-function for the precision
	 */
	double scale = INLA_QKUMAR_PREC_SCALE;
	return map_exp_scale(arg, typ, (void *) &scale);
}

double map_invlogit(double x, map_arg_tp typ, void *UNUSED(param))
{
	/*
	 * extern = exp(local) / (1 + exp(local)) 
	 */
	switch (typ) {
	case MAP_FORWARD:
		return 1.0 / (1.0 + exp(-x));
	case MAP_BACKWARD:
		return log(x / (1.0 - x));
	case MAP_DFORWARD:
	{
		double xx = exp(x);
		return xx / SQR(1.0 + xx);
	}
	case MAP_INCREASING:
		return 1.0;
	default:
		GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
	}
	return 0.0;
}

double map_probability(double x, map_arg_tp typ, void *UNUSED(param))
{
	/*
	 * extern = exp(local) / (1 + exp(local)) 
	 */
	double xx;

	switch (typ) {
	case MAP_FORWARD:
		xx = exp(x);
		return xx / (1.0 + xx);
	case MAP_BACKWARD:
		return log(x / (1.0 - x));
	case MAP_DFORWARD:
		xx = exp(x);
		return xx / SQR(1.0 + xx);
	case MAP_INCREASING:
		return 1.0;
	default:
		GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
	}
	return 0.0;
}

double map_shape_svnig(double arg, map_arg_tp typ, void *param)
{
	return (map_one_plus_exp(arg, typ, param));
}

double map_one_plus_exp(double arg, map_arg_tp typ, void *UNUSED(param))
{
	/*
	 * the mapping for the shape-parameters in the stochvol-nig model. shape = 1 + exp(shape_intern)
	 */
	switch (typ) {
	case MAP_FORWARD:
		/*
		 * extern = func(local) 
		 */
		return 1.0 + exp(arg);
	case MAP_BACKWARD:
		/*
		 * local = func(extern) 
		 */
		return log(arg - 1.0);
	case MAP_DFORWARD:
		/*
		 * d_extern / d_local 
		 */
		return exp(arg);
	case MAP_INCREASING:
		/*
		 * return 1.0 if montone increasing and 0.0 otherwise 
		 */
		return 1.0;
	default:
		abort();
	}
	abort();
	return 0.0;
}

double map_H(double x, map_arg_tp typ, void *UNUSED(param))
{
	/*
	 * extern = 1/2  + 1/2 * exp(local) / (1 + exp(local)) 
	 */
	double xx;

	switch (typ) {
	case MAP_FORWARD:
		xx = exp(x);
		return 0.5 + 0.5 * xx / (1.0 + xx);
	case MAP_BACKWARD:
		return log((2.0 * x - 1.0) / (2.0 * (1.0 - x)));
	case MAP_DFORWARD:
		xx = exp(x);
		return 0.5 * xx / SQR(1.0 + xx);
	case MAP_INCREASING:
		return 1.0;
	default:
		GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
	}
	return 0.0;
}

double map_interval(double x, map_arg_tp typ, void *param)
{
	/*
	 * extern = A  + (B-A) * exp(C * local) / (1 + exp(C * local)) ,  B > A, C>0
	 */
	double *ABC = (double *) param;
	double A = ABC[0];
	double B = ABC[1];
	double C = 1.0;					       /* no longer in use */
	double ex;

	// printf("%g %g %g %g %g\n", x, A, B, C, (double) typ);

	switch (typ) {
	case MAP_FORWARD:
		ex = exp(C * x);
		return A + (B - A) * ex / (1.0 + ex);
	case MAP_BACKWARD:
		return log(-(A - x) / (B - x)) / C;
	case MAP_DFORWARD:
		ex = exp(C * x);
		return C * (B - A) * ex / SQR(1.0 + ex);
	case MAP_INCREASING:
		return 1.0;
	default:
		GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
	}
	return 0.0;
}

double map_group_rho(double x, map_arg_tp typ, void *param)
{
	/*
	 * extern = 
	 */
	assert(param != NULL);
	int ngroups = *((int *) param);
	double xx;

	switch (typ) {
	case MAP_FORWARD:
		xx = exp(x);
		return (xx - 1.0) / (xx + ngroups - 1.0);
	case MAP_BACKWARD:
		return log((1.0 + (ngroups - 1.0) * x) / (1.0 - x));
	case MAP_DFORWARD:
		xx = exp(x);
		return (xx * ngroups) / SQR(xx + ngroups - 1.0);
	case MAP_INCREASING:
		return 1.0;
	default:
		GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
	}
	return 0.0;
}

double map_invtan(double x, map_arg_tp typ, void *UNUSED(param))
{
	/*
	 * y = 2*atan(x), so that |y| <= Pi
	 */
	switch (typ) {
	case MAP_FORWARD:
		/*
		 * extern = func(local) 
		 */
		return 2.0 * atan(x);
	case MAP_BACKWARD:
		/*
		 * local = func(extern) 
		 */
		return tan(x / 2.0);
	case MAP_DFORWARD:
		return 2.0 / (1.0 + SQR(x));
	case MAP_INCREASING:
		return 1.0;
	default:
		GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
	}
	return 0.0;
}

double link_this_should_not_happen(int UNUSED(thread_id), double UNUSED(x), map_arg_tp UNUSED(typ), void *UNUSED(param), double *UNUSED(cov))
{
	/*
	 * the link-functions calls the inverse map-function 
	 */
	FIXMEstderr("This function is called because a wrong link function is used.");
	abort();
	return 0.0;
}

double link_probit(int UNUSED(thread_id), double x, map_arg_tp typ, void *param, double *UNUSED(cov))
{
	/*
	 * the link-functions calls the inverse map-function 
	 */
	return map_invprobit(x, typ, param);
}

double link_robit(int thread_id, double x, map_arg_tp typ, void *param, double *UNUSED(cov))
{
	/*
	 * the link-functions calls the inverse map-function 
	 */
	Link_param_tp *p = (Link_param_tp *) param;
	double dof_intern = p->dof_intern[thread_id][0];

	return map_invrobit(x, typ, (void *) &dof_intern);
}

double link_sn(int thread_id, double x, map_arg_tp typ, void *param, double *UNUSED(cov))
{
	/*
	 * the link-functions calls the inverse map-function 
	 */
	Link_param_tp *p = (Link_param_tp *) param;
	double skew = p->sn_skew[thread_id][0], intercept = p->sn_intercept[thread_id][0], *par[2];

	par[0] = &skew;
	par[1] = &intercept;

	return map_invsn(x, typ, (void *) par);
}

double link_power_logit(int thread_id, double x, map_arg_tp typ, void *param, double *UNUSED(cov))
{
	/*
	 * the link-functions calls the inverse map-function 
	 */
	Link_param_tp *p = (Link_param_tp *) param;
	double power = p->power_intern[thread_id][0];
	double intercept = p->intercept_intern[thread_id][0];
	double *par[2];

	par[0] = &power;
	par[1] = &intercept;

	return map_inv_powerlink_core(x, typ, (void *) par, NULL);
}

double link_tan(int UNUSED(thread_id), double x, map_arg_tp typ, void *param, double *UNUSED(cov))
{
	return map_invtan(x, typ, param);
}

double link_cloglog(int UNUSED(thread_id), double x, map_arg_tp typ, void *param, double *UNUSED(cov))
{
	/*
	 * the link-functions calls the inverse map-function 
	 */
	return map_invcloglog(x, typ, param);
}

double link_ccloglog(int UNUSED(thread_id), double x, map_arg_tp typ, void *param, double *UNUSED(cov))
{
	return map_invccloglog(x, typ, param);
}

double link_loglog(int UNUSED(thread_id), double x, map_arg_tp typ, void *param, double *UNUSED(cov))
{
	/*
	 * the link-functions calls the inverse map-function 
	 */
	return map_invloglog(x, typ, param);
}

double link_cauchit(int UNUSED(thread_id), double x, map_arg_tp typ, void *param, double *UNUSED(cov))
{
	/*
	 * the link-functions calls the inverse map-function 
	 */
	return map_invcauchit(x, typ, param);
}

double link_log(int UNUSED(thread_id), double x, map_arg_tp typ, void *param, double *UNUSED(cov))
{
	/*
	 * the link-functions calls the inverse map-function 
	 */
	return map_exp(x, typ, param);
}

double link_loga(int UNUSED(thread_id), double x, map_arg_tp typ, void *param, double *UNUSED(cov))
{
#define MAP(_x) log((_x)/(1.0 - (_x)))
#define iMAP(_x) (exp(_x)/(1.0+exp(_x)))
#define diMAP(_x) (exp(_x)/SQR(1.0+exp(_x)))

	Link_param_tp *par = (Link_param_tp *) param;
	double a = par->a;

	static inla_loga_table_tp **table = NULL;
	static char first = 1;

	int i, id = 0;
	const int debug = 0;
	double dx = 0.1, xx, range = 25.0, p, pp;

	if (first) {
#pragma omp critical (Name_a862df53446b09ea53b4bd233636df42f8acd1ff)
		if (first) {
			if (debug) {
				fprintf(stderr, "link_loga: init tables\n");
			}
			table = Calloc(GMRFLib_CACHE_LEN(), inla_loga_table_tp *);
			for (i = 0; i < GMRFLib_CACHE_LEN(); i++) {
				table[i] = Calloc(1, inla_loga_table_tp);
				table[i]->a = NAN;
				table[i]->cdf = NULL;
				table[i]->icdf = NULL;
			}
			first = 0;
		}
	}

	GMRFLib_CACHE_SET_ID(id);
	if (a != table[id]->a) {
		int len, llen;
		double *work = NULL, *x_ = NULL, *y = NULL, p_local;
		if (debug) {
			fprintf(stderr, "link_loga: build new table for a=%g [%1d]\n", a, id);
		}
		// count to find the length
		for (xx = -range, len = 0; xx <= 2 * range; xx += (ABS(xx) < 3.0 ? dx / 5.0 : dx), len++);
		work = x_ = Calloc(2 * len, double);
		y = work + len;

		for (xx = -range, i = 0, llen = 0; xx <= 2 * range; xx += (ABS(xx) < 3.0 ? dx / 5.0 : dx), i++, llen++) {
			p_local = iMAP(xx);
			if (p_local == 1.0) {
				// no point of going further if this happens. we need this as a small value of 'a' makes this
				// happen rather quickly
				break;
			}
			y[i] = xx;
			x_[i] = LOG_p(p_local) - a * LOG_1mp(p_local);	// log(p/(1-p)^a)
		}
		len = llen;

		table[id]->a = a;
		table[id]->p_intern_min = y[0];
		table[id]->p_intern_max = y[len - 1];
		table[id]->eta_min = x_[0];
		table[id]->eta_max = x_[len - 1];

		GMRFLib_spline_free(table[id]->cdf);	       /* ok if NULL */
		GMRFLib_spline_free(table[id]->icdf);	       /* ok if NULL */

		table[id]->cdf = GMRFLib_spline_create(x_, y, len);
		table[id]->icdf = GMRFLib_spline_create(y, x_, len);
		Free(work);
	}

	switch (typ) {
	case MAP_FORWARD:
		/*
		 * extern = func(local) 
		 */
		x = TRUNCATE(x, table[id]->eta_min, table[id]->eta_max);
		p = GMRFLib_spline_eval(x, table[id]->cdf);
		return (iMAP(p));

	case MAP_BACKWARD:
		/*
		 * local = func(extern) 
		 */
		x = TRUNCATE(MAP(x), table[id]->p_intern_min, table[id]->p_intern_max);
		return GMRFLib_spline_eval(x, table[id]->icdf);

	case MAP_DFORWARD:
		/*
		 * d_extern / d_local 
		 */
		x = TRUNCATE(x, table[id]->eta_min, table[id]->eta_max);
		p = GMRFLib_spline_eval(x, table[id]->cdf);
		pp = GMRFLib_spline_eval_deriv(x, table[id]->cdf);
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

	return 0.0;
}

double link_neglog(int UNUSED(thread_id), double x, map_arg_tp typ, void *param, double *UNUSED(cov))
{
	/*
	 * the link-functions calls the inverse map-function 
	 */
	return map_negexp(x, typ, param);
}

double link_logit(int UNUSED(thread_id), double x, map_arg_tp typ, void *param, double *UNUSED(cov))
{
	/*
	 * the link-functions calls the inverse map-function 
	 */
	return map_invlogit(x, typ, param);
}

double link_identity(int UNUSED(thread_id), double x, map_arg_tp typ, void *param, double *UNUSED(cov))
{
	/*
	 * the link-functions calls the inverse map-function 
	 */
	return map_identity(x, typ, param);
}

double link_inverse(int UNUSED(thread_id), double x, map_arg_tp typ, void *param, double *UNUSED(cov))
{
	/*
	 * the link-functions calls the inverse map-function 
	 */
	return map_inverse(x, typ, param);
}

double link_logoffset(int thread_id, double x, map_arg_tp typ, void *param, double *cov)
{
	/*
	 * the link-functions calls the inverse map-function 
	 */
	Link_param_tp *p = NULL;
	double beta, off, sign;

	if (!cov) {
		char *msg = NULL;
		GMRFLib_sprintf(&msg, "You need to pass the covariates to the link.model[logoffset] in the inla()-argument 'link.covariates'");
		inla_error_general(msg);
		exit(1);
	}
	if (cov[0] < 0.0) {
		char *msg = NULL;
		GMRFLib_sprintf(&msg, "The covariates to link.model[logoffset] must be all >= 0.0. Yours is [%g].", cov[0]);
		inla_error_general(msg);
		exit(1);
	}

	p = (Link_param_tp *) param;
	beta = exp(p->beta_intern[thread_id][0]);
	off = beta * cov[0];
	sign = (p->variant == 0 ? 1.0 : -1.0);

	switch (typ) {
	case MAP_FORWARD:
		return off + sign * exp(x);
	case MAP_BACKWARD:
		return log(sign * (x - off));
	case MAP_DFORWARD:
		return sign * exp(x);
	case MAP_INCREASING:
		return sign;
	default:
		abort();
	}
	return NAN;
}

double link_logitoffset(int thread_id, double x, map_arg_tp typ, void *param, double *UNUSED(cov))
{
	/*
	 * the link-functions calls the inverse map-function 
	 */
	Link_param_tp *p = NULL;
	double prob;

	p = (Link_param_tp *) param;
	prob = map_probability(p->prob_intern[thread_id][0], MAP_FORWARD, NULL);

	switch (typ) {
	case MAP_FORWARD:
		return prob + (1.0 - prob) * map_probability(x, MAP_FORWARD, NULL);
	case MAP_BACKWARD:
		return map_probability((x - prob) / (1.0 - prob), MAP_BACKWARD, NULL);
	case MAP_DFORWARD:
		return (1.0 - prob) * map_probability(x, MAP_DFORWARD, NULL);
	case MAP_INCREASING:
		return 1.0;
	default:
		abort();
	}
	return NAN;
}

double link_sslogit(int thread_id, double x, map_arg_tp typ, void *param, double *UNUSED(cov))
{
	Link_param_tp *p = NULL;
	double sens, spec, a, b, xx;

	p = (Link_param_tp *) param;
	sens = map_probability(p->sensitivity_intern[thread_id][0], MAP_FORWARD, NULL);
	spec = map_probability(p->specificity_intern[thread_id][0], MAP_FORWARD, NULL);

	// sens * pi + (1-spec)*(1-pi) = a * pi + b
	a = sens + spec - 1.0;
	b = 1.0 - spec;

	switch (typ) {
	case MAP_FORWARD:
		return a / (1.0 + exp(-x)) + b;

	case MAP_BACKWARD:
		assert((x - b) / (x - b - a) > 0.0);
		return log((x - b) / (x - b - a));

	case MAP_DFORWARD:
		xx = exp(-x);
		return a * xx / SQR(1.0 + xx);

	case MAP_INCREASING:
		return (a > 0 ? 1 : 0);
	default:
		abort();
	}
	return NAN;
}

double link_special2(int thread_id, double x, map_arg_tp typ, void *param, double *cov)
{
	/*
	 * the link-functions calls the inverse map-function 
	 */
	Link_param_tp *p = NULL;
	double beta, f;

	if (!cov) {
		char *msg = NULL;
		GMRFLib_sprintf(&msg, "You need to pass the covariate to the link.model[special2] in the inla()-argument 'link.covariates'");
		inla_error_general(msg);
		exit(1);
	}
	if (cov[0] <= 0.0 || cov[0] >= 1.0) {
		char *msg = NULL;
		GMRFLib_sprintf(&msg, "The covariate to link.model[special2] must be between 0 and 1. Your is [%g].", cov[0]);
		inla_error_general(msg);
		exit(1);
	}

	p = (Link_param_tp *) param;
	beta = p->beta[thread_id][0];
	f = (1.0 - cov[0] + cov[0] * exp(beta));

	switch (typ) {
	case MAP_FORWARD:
		return exp(x) * f;
	case MAP_BACKWARD:
		return log(x / f);
	case MAP_DFORWARD:
		return exp(x) * f;
	case MAP_INCREASING:
		return 1.0;
	default:
		abort();
	}
	return NAN;

}

double link_qpoisson(int thread_id, double x, map_arg_tp typ, void *param, double *cov)
{
	double shape, ret = 0.0;
	Link_param_tp *lparam = (Link_param_tp *) param;

	switch (typ) {
	case INVLINK:
	{
		shape = exp(x) + 1.0;
		if (shape > 400.0) {
			ret = SQR(sqrt(shape) + sqrt(0.25) * gsl_cdf_ugaussian_Qinv(lparam->quantile));
		} else {
			ret = gsl_cdf_gamma_Qinv(lparam->quantile, shape, 1.0);
		}
	}
		break;

	case LINK:
	{
		CODE_NEEDED;
	}
		break;

	case DINVLINK:
	{
		double dx = GSL_ROOT4_DBL_EPSILON;	       // about 0.0001 on my laptop
		double wf[] = { 1.0 / 12.0, -2.0 / 3.0, 0.0, 2.0 / 3.0, -1.0 / 12.0 };
		double wf_sum = 0.0;
		int i, nwf = sizeof(wf) / sizeof(double), nwf2 = (nwf - 1) / 2;	/* gives 5 and 2 */

		for (i = 0; i < nwf; i++) {
			wf_sum += wf[i] * link_qpoisson(thread_id, x + (i - nwf2) * dx, INVLINK, param, cov);
		}
		ret = wf_sum / dx;
	}
		break;

	case LINKINCREASING:
	{
		return 1.0;
	}
		break;

	default:
	{
		assert(0 == 1);
	}
		break;

	}

	return (ret);
}

double link_qweibull(int thread_id, double x, map_arg_tp typ, void *param, double *cov)
{
	Link_param_tp *lparam = (Link_param_tp *) param;
	double alpha = map_alpha_weibull(lparam->alpha_intern[thread_id][0], MAP_FORWARD, NULL);
	double ret = 0.0;

	switch (typ) {
	case INVLINK:
	{
		switch (lparam->variant) {
		case 0:
		{
			ret = -1.0 / pow(exp(x), alpha) * LOG_1mp(lparam->quantile);
		}
			break;
		case 1:
		{
			ret = 1.0 / exp(x) * pow(-LOG_1mp(lparam->quantile), 1.0 / alpha);
		}
			break;
		default:
			assert(0 == 1);
		}
	}
		break;

	case LINK:
	{
		CODE_NEEDED;
	}
		break;

	case DINVLINK:
	{
		double dx = GSL_ROOT4_DBL_EPSILON;	       // about 0.0001 on my laptop
		double wf[] = { 1.0 / 12.0, -2.0 / 3.0, 0.0, 2.0 / 3.0, -1.0 / 12.0 };
		double wf_sum = 0.0;
		int i, nwf = sizeof(wf) / sizeof(double), nwf2 = (nwf - 1) / 2;	/* gives 5 and 2 */

		for (i = 0; i < nwf; i++) {
			wf_sum += wf[i] * link_qweibull(thread_id, x + (i - nwf2) * dx, INVLINK, param, cov);
		}
		ret = wf_sum / dx;
	}
		break;

	case LINKINCREASING:
	{
		int ret_val = 0;
		static int do_check[2] = { 1, 1 };

		if (do_check[lparam->variant]) {
#pragma omp critical (Name_e77483d42a508f8162880242479cd58817992341)
			if (do_check[lparam->variant]) {
				if (ret_val !=
				    (link_qweibull(thread_id, x + 1.0, INVLINK, param, cov) >
				     link_qweibull(thread_id, x, INVLINK, param, cov) ? 1 : 0)) {
					FIXME("LINKINCREASING has error in link_qweibull");
					exit(EXIT_FAILURE);
				}
				do_check[lparam->variant] = 0;
			}
		}
		return (ret_val);
	}
		break;

	default:
	{
		assert(0 == 1);
	}
		break;

	}

	return (ret);
}

double link_qgamma(int thread_id, double x, map_arg_tp typ, void *param, double *cov)
{
	Link_param_tp *lparam = (Link_param_tp *) param;
	double s = (lparam->scale ? lparam->scale[lparam->idx] : 1.0);
	// double phi_param = map_exp(lparam->log_prec[thread_id][0], MAP_FORWARD, NULL);
	double phi_param = exp(lparam->log_prec[thread_id][0]);
	double shape = phi_param * s;
	double ret = 0.0;

	switch (typ) {
	case INVLINK:
	{
		// ret = exp(x) * shape / MATHLIB_FUN(qgamma) (lparam->quantile, shape/100, 1.0, 1, 0);
		ret = exp(x) * shape / inla_qgamma_cache(shape, lparam->quantile);
	}
		break;

	case LINK:
	{
		CODE_NEEDED;
	}
		break;

	case DINVLINK:
	{
		double dx = GSL_ROOT4_DBL_EPSILON;	       // about 0.0001 on my laptop
		double wf[] = { 1.0 / 12.0, -2.0 / 3.0, 0.0, 2.0 / 3.0, -1.0 / 12.0 };
		double wf_sum = 0.0;
		int i, nwf = sizeof(wf) / sizeof(double), nwf2 = (nwf - 1) / 2;	/* gives 5 and 2 */

		for (i = 0; i < nwf; i++) {
			wf_sum += wf[i] * link_qgamma(thread_id, x + (i - nwf2) * dx, INVLINK, param, cov);
		}
		ret = wf_sum / dx;
	}
		break;

	case LINKINCREASING:
	{
		ret = 1;
	}
		break;

	default:
	{
		assert(0 == 1);
	}
		break;

	}

	return (ret);
}

double link_qexppower(int thread_id, double x, map_arg_tp typ, void *param, double *UNUSED(cov))
{
	// do caching on two levels

	typedef struct {
		// level 1: if all these match
		double p;
		double w;
		double lprec;
		double lpower;
		double qval1;

		// level 2: if beta match only
		double beta;
		double lg_expr;
		double qval2;
	} lcache_t;

	static lcache_t **llcache = NULL;
	if (!llcache) {
#pragma omp critical (Name_2396789afcc20ddee4600d09ab8d0fe4a104e9f3)
		if (!llcache) {
			llcache = Calloc(GMRFLib_CACHE_LEN(), lcache_t *);
		}
	}

	int cidx = 0;
	GMRFLib_CACHE_SET_ID(cidx);
	if (!llcache[cidx]) {
#pragma omp critical (Name_c393bf22256042fb97a79700a66d05c333658625)
		if (!llcache[cidx]) {
			lcache_t *ptr = Calloc(1, lcache_t);
			ptr->qval1 = NAN;
			ptr->qval2 = NAN;
			llcache[cidx] = ptr;
		}
	}
	lcache_t *lc = llcache[cidx];

	Link_param_tp *lparam = (Link_param_tp *) param;
	int idx = lparam->idx;
	double ret = NAN;
	double p = lparam->quantile;
	double w = lparam->scale[idx];
	double lprec = lparam->log_prec[thread_id][0] + log(w);
	double lpower = lparam->log_power[thread_id][0];

	// level 1 caching
	if (lc->p == p && lc->w == w && lc->lprec == lprec && lc->lpower == lpower) {
		switch (typ) {
		case INVLINK:
			ret = x - lc->qval1;
			break;

		case LINK:
			CODE_NEEDED;
			break;

		case DINVLINK:
			ret = 1.0;
			break;

		case LINKINCREASING:
			ret = 1;
			break;

		default:
			assert(0 == 1);
			break;
		}
		return ret;
	}
	// level 2 caching
	lc->p = p;
	lc->w = w;
	lc->lprec = lprec;
	lc->lpower = lpower;

	double beta = 1.0 + exp(lpower);		       // map_one_plus_exp
	double shape = 1.0 / beta;
	double p2 = ABS(p - 0.5) * 2.0;

	if (lc->beta != beta) {
		lc->beta = beta;
		lc->lg_expr = exp(0.5 * (my_gsl_sf_lngamma(1.0 / beta) - my_gsl_sf_lngamma(3.0 / beta)));
		lc->qval2 = pow(MATHLIB_FUN(qgamma) (p2, shape, 1.0, 1, 0), 1.0 / beta);
	}

	double sign = DSIGN(p - 0.5);
	double sigma = exp(-0.5 * lprec);
	double alpha = sigma * lc->lg_expr;
	lc->qval1 = sign * alpha * lc->qval2;

	switch (typ) {
	case INVLINK:
		// double ialpha = 1.0 / alpha, lambda = pow(ialpha, beta), scale = 1.0/lambda;
		// ret = x - (sign * pow(scale * inla_qgamma_cache(shape, p2), 1.0/beta));
		ret = x - lc->qval1;
		break;

	case LINK:
		CODE_NEEDED;
		break;

	case DINVLINK:
		ret = 1.0;
		break;

	case LINKINCREASING:
		ret = 1;
		break;

	default:
		assert(0 == 1);
		break;
	}

	return (ret);
}

double link_qbinomial(int thread_id, double x, map_arg_tp typ, void *param, double *cov)
{
	// individual link
	Link_param_tp *lparam = (Link_param_tp *) param;
	double q, ret = 0.0;

	switch (typ) {
	case INVLINK:
	{
		q = 1.0 / (1.0 + exp(-x));
		ret = MATHLIB_FUN(qbeta) (lparam->quantile, q + 1.0, 1.0 - q, 0, 0);
	}
		break;

	case LINK:
	{
		CODE_NEEDED;
	}
		break;

	case DINVLINK:
	{
		double dx = GSL_ROOT4_DBL_EPSILON;	       // about 0.0001 on my laptop
		double wf[] = { 1.0 / 12.0, -2.0 / 3.0, 0.0, 2.0 / 3.0, -1.0 / 12.0 };
		double wf_sum = 0.0;
		int i, nwf = sizeof(wf) / sizeof(double), nwf2 = (nwf - 1) / 2;	/* gives 5 and 2 */

		for (i = 0; i < nwf; i++) {
			wf_sum += wf[i] * link_qbinomial(thread_id, x + (i - nwf2) * dx, INVLINK, param, cov);
		}
		ret = wf_sum / dx;
	}
		break;

	case LINKINCREASING:
	{
		return 1.0;
	}
		break;

	default:
	{
		assert(0 == 1);
	}
		break;

	}

	return (ret);
}

double link_pqbinomial(int thread_id, double x, map_arg_tp typ, void *param, double *cov)
{
	// population link
	Link_param_tp *lparam = (Link_param_tp *) param;
	double q, ret = 0.0;

	switch (typ) {
	case INVLINK:
	{
		q = 1.0 / (1.0 + exp(-x));
		ret = MATHLIB_FUN(qbeta) (lparam->quantile, lparam->Ntrial * q + 1.0, lparam->Ntrial * (1.0 - q), 0, 0);
	}
		break;

	case LINK:
	{
		CODE_NEEDED;
	}
		break;

	case DINVLINK:
	{
		double dx = GSL_ROOT4_DBL_EPSILON;	       // about 0.0001 on my laptop
		double wf[] = { 1.0 / 12.0, -2.0 / 3.0, 0.0, 2.0 / 3.0, -1.0 / 12.0 };
		double wf_sum = 0.0;
		int i, nwf = sizeof(wf) / sizeof(double), nwf2 = (nwf - 1) / 2;	/* gives 5 and 2 */

		for (i = 0; i < nwf; i++) {
			wf_sum += wf[i] * link_pqbinomial(thread_id, x + (i - nwf2) * dx, INVLINK, param, cov);
		}
		ret = wf_sum / dx;
	}
		break;

	case LINKINCREASING:
	{
		return 1.0;
	}
		break;

	default:
	{
		assert(0 == 1);
	}
		break;

	}

	return (ret);
}

double link_test1(int thread_id, double x, map_arg_tp typ, void *param, double *cov)
{
	/*
	 * the link-functions calls the inverse map-function 
	 */
	Link_param_tp *p = NULL;
	double beta;

	p = (Link_param_tp *) param;
	beta = p->beta[thread_id][0];

	return map_exp(x - beta * cov[0], typ, param);
}

double link_special1(int thread_id, double x, map_arg_tp typ, void *param, double *cov)
{
	Link_param_tp *p = NULL;
	int i;
	double covariate_contribution, h = 1.0E-4, prec;

	p = (Link_param_tp *) param;
	prec = map_precision(p->log_prec[thread_id][0], MAP_FORWARD, NULL);
	covariate_contribution = 0.0;
	for (i = 0; i < p->order; i++) {
		covariate_contribution += p->betas[i][thread_id][0] * cov[i];
	}
	assert(!ISNAN(covariate_contribution));

	switch (typ) {
	case MAP_FORWARD:
	{
		return gsl_cdf_lognormal_Pinv(inla_cdf_normal(x), covariate_contribution - 0.5 / prec, 1.0 / sqrt(prec));
	}
		break;

	case MAP_BACKWARD:
	{
		return GMRFLib_cdfnorm_inv(gsl_cdf_lognormal_P(x, covariate_contribution - 0.5 / prec, 1.0 / sqrt(prec)));
	}
		break;

	case MAP_DFORWARD:
	{
		return (gsl_cdf_lognormal_Pinv(inla_cdf_normal(x + h), covariate_contribution - 0.5 / prec, 1.0 / sqrt(prec)) -
			gsl_cdf_lognormal_Pinv(inla_cdf_normal(x - h), covariate_contribution - 0.5 / prec, 1.0 / sqrt(prec))) / (2.0 * h);
	}
		break;

	case MAP_INCREASING:
	{
		return 1.0;
	}
		break;

	default:
		abort();
	}

	abort();
	return 0.0;
}

double inla_boxcox_core(double y, double lambda)
{
	// simplify(subs(O=0,series((x^lam-1)/lam, lam=0, 10)));
	const double eps = 1.0E-4;
	double val = 0.0;

	if (ABS(lambda) < eps) {
		double ly = log(y);
		double a = ly * lambda;
		val = ly *
		    (1.0 + (a / 2.0) *
		     (1.0 + (a / 3.0) *
		      (1.0 + (a / 4.0) *
		       (1.0 + (a / 5.0) * (1.0 + (a / 6.0) * (1.0 + (a / 7.0) * (1.0 + (a / 8.0) * (1.0 + (a / 9.0) * (1.0 + (a / 10.0))))))))));
	} else {
		val = (pow(y, lambda) - 1.0) / lambda;
	}
	return val;
}
double inla_boxcox(double y, double mean, double lambda)
{
	return inla_boxcox_core(y, lambda) - (mean > 0.0 ? inla_boxcox_core(mean, lambda) : 0.0);
}
