
int inla_integrate_func_ORIG(double *d_mean, double *d_stdev, double *d_mode, GMRFLib_density_tp * density, map_func_tp * func,
			     void *func_arg, GMRFLib_transform_array_func_tp * tfunc)
{
	/*
	 * We need to integrate to get the transformed mean and variance. Use a simple Simpsons-rule.  The simple mapping we did before was not good enough,
	 * obviously... 
	 */
	if (!func && !tfunc) {
		*d_mean = density->user_mean;
		*d_stdev = density->user_stdev;
		*d_mode = density->user_mode;

		return GMRFLib_SUCCESS;
	}
#define _MAP_X(_x_user) (func ? func(_x_user, MAP_FORWARD, func_arg) : \
			(tfunc ? tfunc->func(_x_user, GMRFLib_TRANSFORM_FORWARD, tfunc->arg, tfunc->cov) : \
			 (_x_user)))
#define _MAP_STDEV(_stdev_user, _mean_user) (func ? (_stdev_user)*ABS(func(_mean_user, MAP_DFORWARD, func_arg)) : \
	(tfunc ? (_stdev_user)*ABS(tfunc->func(_mean_user, GMRFLib_TRANSFORM_DFORWARD, tfunc->arg, tfunc->cov)) : \
	 (_stdev_user)))

#define _TRANSFORMED_LOGDENS(_x, _logdens) (func ? ((_logdens) - log(ABS(func(_x, MAP_DFORWARD, func_arg)))) : \
					   (tfunc ? ((_logdens) - log(ABS(tfunc->func(_x, GMRFLib_TRANSFORM_DFORWARD, tfunc->arg, tfunc->cov)))) : \
				       (_logdens)))

	int i, k, debug = 0;
	int npm = 4 * GMRFLib_INT_NUM_POINTS;
	double dxx, d, xx_local, fval;
	double w[2] = { 4.0, 2.0 };
	double sum[3] = { 0.0, 0.0, 0.0 };
	double *xpm, *ldm, *work;

	double tlogdens, tlogdens_max = NAN;
	int idx_tlogdens_max = -1;

	work = Calloc(2 * npm, double);
	xpm = work;
	ldm = work + npm;

	dxx = (density->x_max - density->x_min) / (npm - 1.0);
	xpm[0] = density->x_min;
	for (i = 1; i < npm; i++) {
		xpm[i] = xpm[0] + i * dxx;
	}
	GMRFLib_evaluate_nlogdensity(ldm, xpm, npm, density);

	sum[0] = exp(ldm[0]) + exp(ldm[npm - 1]);
	xx_local = GMRFLib_density_std2user(xpm[0], density);
	sum[1] = exp(ldm[0]) * _MAP_X(xx_local);
	sum[2] = exp(ldm[0]) * SQR(_MAP_X(xx_local));
	xx_local = GMRFLib_density_std2user(xpm[npm - 1], density);
	sum[1] += exp(ldm[npm - 1]) * _MAP_X(xx_local);
	sum[2] += exp(ldm[npm - 1]) * SQR(_MAP_X(xx_local));

	for (i = 1, k = 0; i < npm - 1; i++, k = (k + 1) % 2) {
		d = exp(ldm[i]);
		xx_local = GMRFLib_density_std2user(xpm[i], density);
		fval = _MAP_X(xx_local);

		sum[0] += d * w[k];
		sum[1] += d * w[k] * fval;
		sum[2] += d * w[k] * SQR(fval);

		/*
		 * extract the mode
		 */
		tlogdens = _TRANSFORMED_LOGDENS(xx_local, ldm[i]);
		if (idx_tlogdens_max < 0 || tlogdens > tlogdens_max) {
			idx_tlogdens_max = i;
			tlogdens_max = tlogdens;
		}
	}
	sum[0] *= dxx / 3.0;				       /* this should be 1.0 */
	sum[1] *= dxx / 3.0;
	sum[2] *= dxx / 3.0;

	if (debug) {
		printf("NEW intergral %g\n", sum[0]);
		printf("NEW mean %g simple %g\n", sum[1] / sum[0], _MAP_X(density->user_mean));
		printf("NEW stdev %g simple %g\n", sqrt(sum[2] / sum[0] - SQR(sum[1] / sum[0])),
		       _MAP_STDEV(density->user_stdev, density->user_mean));
	}

	*d_mean = sum[1] / sum[0];
	*d_stdev = sqrt(DMAX(0.0, sum[2] / sum[0] - SQR(*d_mean)));

	/*
	 * adjust the  mode. use a 2.order polynomial
	 *
	 * > ans := solve({ a+b*x0+c*x0^2=d0, a+b*x1+c*x1^2=d1,a+b*x2+c*x2^2=d2}, [a,b,c]):
	 * > a := rhs(ans[1,1]):                                                           
	 * > b := rhs(ans[1,2]):
	 * > c := rhs(ans[1,3]):
	 * > simplify(-b/(2*c));
	 *         2        2        2        2        2        2
	 *    d0 x1  - d0 x2  - d1 x0  + d1 x2  + d2 x0  - d2 x1
	 *    ---------------------------------------------------
	 *      2 (d0 x1 - d0 x2 - d1 x0 + d1 x2 + d2 x0 - d2 x1)
	 */

	double xx[3], tld[3];

	idx_tlogdens_max = IMAX(1, IMIN(npm - 2, idx_tlogdens_max));	/* so that we can do interpolation */
	for (k = 0; k < 3; k++) {
		i = idx_tlogdens_max - 1 + k;
		xx[k] = GMRFLib_density_std2user(xpm[i], density);
		tld[k] = _TRANSFORMED_LOGDENS(xx[k], ldm[i]);
		xx[k] = _MAP_X(xx[k]);			       /* need it in the transformed scale */
	}
	*d_mode = (tld[0] * xx[1] * xx[1] - tld[0] * xx[2] * xx[2] - tld[1] * xx[0] * xx[0] +
		   tld[1] * xx[2] * xx[2] + tld[2] * xx[0] * xx[0] - tld[2] * xx[1] * xx[1]) / (tld[0] * xx[1] - tld[0] * xx[2] -
												tld[1] * xx[0] + tld[1] * xx[2] +
												tld[2] * xx[0] - xx[1] * tld[2]) / 0.2e1;

	Free(work);
#undef _MAP_X
#undef _MAP_STDEV
#undef _TRANSFORMED_LOGDENS
	return GMRFLib_SUCCESS;
}
































		/*
		 * THERE are up to STRATA_MAXTHETA of the probs, called prob 1... 10 
		 */
		ds->data_observations.probN_intern = Calloc(STRATA_MAXTHETA, double **);
		ds->data_nfixed = Calloc(STRATA_MAXTHETA, int);
		ds->data_nprior = Calloc(STRATA_MAXTHETA, Prior_tp);

		for (int count = 0; count < STRATA_MAXTHETA; count++) {
			char *ctmp = NULL;

			/*
			 * the zeroinflation prob-parameter
			 */
			GMRFLib_sprintf(&ctmp, "INITIAL%1d", count + 1);
			tmp = iniparser_getdouble(ini, inla_string_join(secname, ctmp), -1.0);

			GMRFLib_sprintf(&ctmp, "FIXED%1d", count + 1);
			ds->data_nfixed[count] = iniparser_getboolean(ini, inla_string_join(secname, ctmp), 0);
			if (!ds->data_nfixed[count] && mb->reuse_mode) {
				tmp = mb->theta_file[mb->theta_counter_file++];
			}
			HYPER_NEW(ds->data_observations.probN_intern[count], tmp);

			if (mb->verbose) {
				printf("\t\tinitialise prob%1d_intern[%g]\n", count + 1, ds->data_observations.probN_intern[count][0][0]);
				printf("\t\tfixed%1d=[%1d]\n", count + 1, ds->data_nfixed[count]);
			}
			inla_read_priorN(mb, ini, sec, &(ds->data_nprior[count]), "GAUSSIAN-std", count + 1, NULL);

			/*
			 * add theta 
			 */
			if (!ds->data_nfixed[count]) {
				mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
				mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
				mb->theta_hyperid[mb->ntheta] = ds->data_nprior[count].hyperid;
				mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
				mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
				mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);

				GMRFLib_sprintf(&ctmp, "intern zero-probability%1d for zero-inflated nbinomial_strata2", count + 1);
				mb->theta_tag[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				GMRFLib_sprintf(&ctmp, "zero-probability%1d for zero-inflated nbinomial_strata2", count + 1);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				GMRFLib_sprintf(&msg, "%s-parameter%1d", secname, count + 1);
				mb->theta_dir[mb->ntheta] = msg;

				mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
				mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
				mb->theta_from[mb->ntheta] = Strdup(ds->data_nprior[count].from_theta);
				mb->theta_to[mb->ntheta] = Strdup(ds->data_nprior[count].to_theta);

				mb->theta[mb->ntheta] = ds->data_observations.probN_intern[count];
				mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
				mb->theta_map[mb->ntheta] = map_probability;
				mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
				mb->theta_map_arg[mb->ntheta] = NULL;
				mb->ntheta++;
				ds->data_ntheta++;
			}
		}
		break;

		/*
		 * THERE are up to STRATA_MAXTHETA of the probs, called prob 1... 10 
		 */
		ds->data_observations.probN_intern = Calloc(STRATA_MAXTHETA, double **);
		ds->data_nfixed = Calloc(STRATA_MAXTHETA, int);
		ds->data_nprior = Calloc(STRATA_MAXTHETA, Prior_tp);

		for (int count = 0; count < STRATA_MAXTHETA; count++) {
			char *ctmp = NULL;

			/*
			 * the zeroinflation prob-parameter
			 */
			GMRFLib_sprintf(&ctmp, "INITIAL%1d", count + 1);
			tmp = iniparser_getdouble(ini, inla_string_join(secname, ctmp), -1.0);

			GMRFLib_sprintf(&ctmp, "FIXED%1d", count + 1);
			ds->data_nfixed[count] = iniparser_getboolean(ini, inla_string_join(secname, ctmp), 0);
			if (!ds->data_nfixed[count] && mb->reuse_mode) {
				tmp = mb->theta_file[mb->theta_counter_file++];
			}
			HYPER_NEW(ds->data_observations.probN_intern[count], tmp);

			if (mb->verbose) {
				printf("\t\tinitialise prob%1d_intern[%g]\n", count + 1, ds->data_observations.probN_intern[count][0][0]);
				printf("\t\tfixed%1d=[%1d]\n", count + 1, ds->data_nfixed[count]);
			}
			inla_read_priorN(mb, ini, sec, &(ds->data_nprior[count]), "GAUSSIAN-std", count + 1, NULL);

			/*
			 * add theta 
			 */
			if (!ds->data_nfixed[count]) {
				mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
				mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
				mb->theta_hyperid[mb->ntheta] = ds->data_nprior[count].hyperid;
				mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
				mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
				mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);

				GMRFLib_sprintf(&ctmp, "intern zero-probability%1d for zero-inflated nbinomial_strata2", count + 1);
				mb->theta_tag[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				GMRFLib_sprintf(&ctmp, "zero-probability%1d for zero-inflated nbinomial_strata2", count + 1);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				GMRFLib_sprintf(&msg, "%s-parameter%1d", secname, count + 1);
				mb->theta_dir[mb->ntheta] = msg;

				mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
				mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
				mb->theta_from[mb->ntheta] = Strdup(ds->data_nprior[count].from_theta);
				mb->theta_to[mb->ntheta] = Strdup(ds->data_nprior[count].to_theta);

				mb->theta[mb->ntheta] = ds->data_observations.probN_intern[count];
				mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
				mb->theta_map[mb->ntheta] = map_probability;
				mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
				mb->theta_map_arg[mb->ntheta] = NULL;
				mb->ntheta++;
				ds->data_ntheta++;
			}
		}
		break;
