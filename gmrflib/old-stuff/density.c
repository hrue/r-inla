	if (xx < s->xmin || xx > s->xmax) {
		if (extrapolate) {
			// maybe I should put this into the GMRFLib_spline_tp as a parameter...
			double deriv;
			if (xx > s->xmax) {
				deriv = GMRFLib_spline_eval_deriv_x(s->xmax, s);
				val = GMRFLib_spline_eval(s->xmax, s) + deriv * (xx - s->xmax);
			} else if (xx < s->xmin) {
				deriv = GMRFLib_spline_eval_deriv_x(s->xmin, s);
				val = GMRFLib_spline_eval(s->xmin, s) + deriv * (xx - s->xmin);
			} else {
				assert(0 == 1);
			}
		} else {
			val = NAN;
		}
	} else {













int GMRFLib_init_density_ORIG(GMRFLib_density_tp * density, int lookup_tables)
{
	/*
	 * initialize 'density': compute the mean, stdev and the norm_const (for the log spline fit) 
	 */
	int i, k, debug = 0, np = GMRFLib_INT_NUM_POINTS, npm = 2 * np;
	double result = 0.0, error, tmp, low = 0.0, high = 0.0, xval, ldens_max = -FLT_MAX, *xpm =
	    NULL, *ldm = NULL, *xp = NULL, integral, dx = 0.0, m1, m2, m3, x0, x1, d0, d1;

	if (!density) {
		return GMRFLib_SUCCESS;
	}

	Calloc_init(2 * npm + 4 * np);

	if (density->type == GMRFLib_DENSITY_TYPE_GAUSSIAN) {
		// density->mean = density->mean_gaussian;
		// density->stdev = density->stdev_gaussian;
		density->skewness = 0.0;
		density->x_min = -GMRFLib_DENSITY_INTEGRATION_LIMIT * density->stdev + density->mean;
		density->x_max = GMRFLib_DENSITY_INTEGRATION_LIMIT * density->stdev + density->mean;
	} else if (density->type == GMRFLib_DENSITY_TYPE_SKEWNORMAL) {
		/*
		 * for the skew-normal we know the moments 
		 */
		double mom[3] = { 0, 0, 0 };
		GMRFLib_sn_par2moments(&mom[0], &mom[1], &mom[2], density->sn_param);
		density->mean = mom[0];
		density->stdev = mom[1];
		density->skewness = mom[2];
		density->x_min = -GMRFLib_DENSITY_INTEGRATION_LIMIT * density->stdev + density->mean;
		density->x_max = GMRFLib_DENSITY_INTEGRATION_LIMIT * density->stdev + density->mean;
	} else {
		GMRFLib_ASSERT(density->type == GMRFLib_DENSITY_TYPE_SCGAUSSIAN, GMRFLib_ESNH);

		low = density->x_min;
		high = density->x_max;
		dx = (high - low) / (npm - 1.0);

		if (debug) {
			P(low);
			P(high);
			P(dx);
		}

		if (1) {
			double ldmax, log_integral;
			double w[2] = { 4.0, 2.0 };

			xpm = Calloc_get(npm);
			ldm = Calloc_get(npm);
			for (xval = low, i = 0; i < npm; xval += dx, i++) {
				xpm[i] = xval;
			}
			density->log_norm_const = 0.0;
			GMRFLib_evaluate_nlogdensity(ldm, xpm, npm, density);

			if (debug) {
				for (i = 0; i < npm; i++) {
					printf("INIT: i %d xpm %g ldm %g\n", i, xpm[i], ldm[i]);
				}
			}

			ldmax = GMRFLib_max_value(ldm, npm, NULL);
			GMRFLib_adjust_vector(ldm, npm);       /* so its well-behaved... */

			integral = exp(ldm[0]) + exp(ldm[npm - 1]);
			for (i = 1, k = 0; i < npm - 1; i++, k = (k + 1) % 2) {
				integral += exp(ldm[i]) * w[k];
			}
			integral *= dx / 3.0;

			log_integral = log(integral);
			for (i = 0; i < npm; i++) {
				ldm[i] -= log_integral;
			}
			density->log_norm_const = log_integral + ldmax;

			d0 = exp(ldm[0]);
			d1 = exp(ldm[npm - 1]);
			x0 = xpm[0];
			x1 = xpm[npm - 1];
			m1 = x0 * d0 + x1 * d1;
			m2 = SQR(x0) * d0 + SQR(x1) * d1;
			m3 = gsl_pow_3(x0) * d0 + gsl_pow_3(x1) * d1;
			for (i = 1, k = 0; i < npm - 1; i++, k = (k + 1) % 2) {
				double d, x, x2, x3;

				d = exp(ldm[i]) * w[k];
				x = xpm[i];
				x2 = x * x;
				x3 = x * x2;

				m1 += x * d;
				m2 += x2 * d;
				m3 += x3 * d;
			}
			m1 *= dx / 3.0;
			m2 *= dx / 3.0;
			m3 *= dx / 3.0;

			if (debug) {
				P(m1);
				P(m2);
				P(m3);
			}

			density->mean = m1;
			density->stdev = sqrt(DMAX(0.0, m2 - SQR(m1)));
			density->skewness = (m3 - 3.0 * m1 * SQR(density->stdev) - gsl_pow_3(m1)) / gsl_pow_3(density->stdev);
		} else {
			double ldens;
			double eps = GMRFLib_eps(1. / 2.);
			GMRFLib_density_properties_tp prop;
			gsl_function F;

			low = density->x_min;
			high = density->x_max;
			dx = (high - low) / (npm - 1.0);

			prop.density = density;
			F.function = GMRFLib_evaluate_density__intern;
			F.params = (void *) &prop;

			/*
			 * the __intern function *use* the norm_const, so first we need to set it temporary to 1.0 
			 */
			density->log_norm_const = 0.0;
			for (xval = low; xval < high; xval += (high - low) / 20.0) {
				GMRFLib_evaluate_logdensity(&ldens, xval, density);
				if (xval == low || ldens > ldens_max) {
					ldens_max = ldens;
				}
			}
			density->log_norm_const = ldens_max;
			GMRFLib_gsl_integration_wrapper(&F, low, high, eps, eps, &result, &error);
			density->log_norm_const = log(result) + ldens_max;

			/*
			 * ...and from now on, the density is normalised 
			 */
			F.function = GMRFLib_evaluate_density_power__intern;	/* use this function for f(x)*x^power */

			prop.power = 1;
			GMRFLib_gsl_integration_wrapper(&F, low, high, eps, eps, &result, &error);
			density->mean = result;

			prop.power = 2;
			GMRFLib_gsl_integration_wrapper(&F, low, high, eps, eps, &result, &error);
			tmp = result - SQR(density->mean);
			density->stdev = sqrt(DMAX(tmp, 0.0));

			/*
			 *  I do not care to add the skewness computations right now, as this code part is obsolete!
			 */
			fprintf(stderr, "\n\n\nTODO: added computation of skewness for a density, HERE!!!\n\n\n");
			abort();
		}
	}

	/*
	 * for convenience, here the mean and the stdev in the users scale 
	 */
	density->user_mean = density->std_stdev * density->mean + density->std_mean;
	density->user_stdev = density->std_stdev * density->stdev;

	if (density->type == GMRFLib_DENSITY_TYPE_GAUSSIAN) {
		density->user_mode = density->user_mean;
	} else {
		density->user_mode = NAN;		       /* yes, this is the value if its not computed */
	}

	/*
	 * new style speedup 
	 */
	if (lookup_tables) {
		/*
		 * build fast lookup tables for P() and Pinv() calculations using linear interpolation 
		 */
		double *p = NULL, *val = NULL, *dens;
		int some_nans = 0;

		if (!xpm) {
			/*
			 * in case we haven't done this before; same as above 
			 */
			low = density->x_min;
			high = density->x_max;
			dx = (high - low) / (npm - 1.0);

			xpm = Calloc_get(npm);
			ldm = Calloc_get(npm);
			for (xval = low, i = 0; i < npm; xval += dx, i++) {
				xpm[i] = xval;
			}
			GMRFLib_evaluate_nlogdensity(ldm, xpm, npm, density);
			GMRFLib_adjust_vector(ldm, npm);       /* so its well-behaved... */
		}

		/*
		 * find the mode fitting a quadratic around the best point, like a one-step Newton-Raphson, since we have so excellent initial value
		 */
		int imax;
		GMRFLib_max_value(ldm, npm - 1, &imax);
		density->user_mode = density->std_mean + density->std_stdev *
		    (xpm[imax] - (((ldm[imax + 1] - ldm[imax - 1]) / (2.0 * dx)) / ((ldm[imax + 1] - 2.0 * ldm[imax] + ldm[imax - 1]) / SQR(dx))));

		dens = Calloc_get(np);
		val = Calloc_get(np);
		p = Calloc_get(np);
		xp = Calloc_get(np);

		for (i = 0; i < np; i++) {
			xp[i] = xpm[2 * i];
			dens[i] = exp(ldm[2 * i]);
		}
		val[0] = 0.0;
		integral = val[0];
		for (i = 1; i < np; i++) {
			val[i] = (dens[i - 1] + dens[i] + 4.0 * exp(ldm[2 * i - 1])) * dx / 6.0;
			val[i] = DMAX(0.0, val[i]);
			integral += val[i];
		}

		p[0] = 0.0;
		for (i = 1; i < np - 1; i++) {
			p[i] = p[i - 1] + val[i] / integral;
			some_nans = some_nans || isnan(p[i]);
		}
		p[np - 1] = 1.0;

		if (!some_nans) {
			density->P = Calloc(1, GMRFLib_spline_tp);
			density->Pinv = Calloc(1, GMRFLib_spline_tp);
			GMRFLib_EWRAP0_GSL_PTR(density->P->accel = gsl_interp_accel_alloc());
			GMRFLib_EWRAP0_GSL_PTR(density->P->spline = gsl_spline_alloc(gsl_interp_linear, (unsigned int) np));
			GMRFLib_EWRAP0_GSL(gsl_spline_init(density->P->spline, xp, p, (unsigned int) np));


			GMRFLib_EWRAP0_GSL_PTR(density->Pinv->accel = gsl_interp_accel_alloc());
			GMRFLib_EWRAP0_GSL_PTR(density->Pinv->spline = gsl_spline_alloc(gsl_interp_linear, (unsigned int) np));
			GMRFLib_EWRAP0_GSL(gsl_spline_init(density->Pinv->spline, p, xp, (unsigned int) np));
		} else {
			density->P = NULL;
		}
	} else {
		density->P = NULL;
		density->Pinv = NULL;
	}

	Calloc_free();

	return GMRFLib_SUCCESS;
}
