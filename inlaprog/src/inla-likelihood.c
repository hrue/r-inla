
/* inla-likelihood.c
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


double inla_compute_saturated_loglik(int thread_id, int idx, GMRFLib_logl_tp *UNUSED(loglfunc), double *x_vec, void *arg)
{
	inla_tp *a = (inla_tp *) arg;
	// assert(loglfunc == loglikelihood_inla);
	return inla_compute_saturated_loglik_core(thread_id, idx, a->loglikelihood[idx], x_vec, a->loglikelihood_arg[idx]);
}

double inla_compute_saturated_loglik_core(int thread_id, int idx, GMRFLib_logl_tp *loglfunc, double *x_vec, void *arg)
{
	double prec_high = 1.0E4, prec_low = 1.0 / prec_high, eps = 1.0E-6;
	double log_prec_high = log(prec_high), log_prec_low = log(prec_low);
	double prec, x, xsol, xnew, f, deriv, dderiv, arr[3], arr_old[3], steplen = GSL_ROOT4_DBL_EPSILON, w;
	int niter, niter_min = 5, niter_max = 100, stencil = 5;
	const int debug = 0;

	(void) loglfunc(thread_id, NULL, NULL, 0, 0, NULL, NULL, NULL, NULL);
	x = xnew = xsol = 0.0;
	for (niter = 0; niter < niter_max; niter++) {
		w = DMIN(1.0, (double) niter / (double) niter_min);
		prec = exp(log_prec_high * (1.0 - w) + log_prec_low * w);

		Memcpy(arr_old, arr, sizeof(arr));
		GMRFLib_2order_taylor(thread_id, &arr[0], &arr[1], &arr[2], NULL, 1.0, x, idx, x_vec, loglfunc, arg, &steplen, &stencil);
		if (ISNAN(arr[0]) || ISINF(arr[0])) {
			Memcpy(arr, arr_old, sizeof(arr));
			break;
		}

		f = arr[0] - 0.5 * prec * SQR(x);
		deriv = arr[1] - prec * x;
		dderiv = DMIN(0.0, arr[2]) - prec;

		xnew = x - deriv / dderiv;
		if (debug) {
			printf("ITER%1d: idx %d x %.6g xnew %.6g f %.6g deriv %.6g dderiv %.6g prec %.6g\n", niter, idx, x, xnew, f, deriv, dderiv,
			       prec);
		}
		x = xnew;

		if (niter > 2 + niter_min && ABS(deriv / dderiv) < eps) {
			break;
		}
	}

	return (arr[0]);
}

int inla_read_data_likelihood(inla_tp *mb, dictionary *UNUSED(ini), int UNUSED(sec))
{
	/*
	 * read data from file 
	 */
#define _DIM_A  (4096L)

	double *x = NULL, *a[_DIM_A];
	int n = 0, na, i, j, ii, idiv = 0, k, ncol_data_all = -1;
	Data_section_tp *ds = &(mb->data_sections[mb->nds - 1]);

	double *attr = NULL;
	int n_attr = 0;
	inla_read_data_all(&attr, &n_attr, ds->attr_file.name, NULL);
	assert(n_attr >= 1);
	if (n_attr > 0) {
		n_attr--;
		attr++;
	} else {
		n_attr = 0;
		attr = NULL;
	}
	ds->data_observations.attr = attr;
	ds->data_observations.n_attr = n_attr;
	if (n_attr > 0 && mb->verbose) {
		printf("\t\tmdata.nattributes = %d\n", n_attr);
		for (i = 0; i < n_attr; i++) {
			printf("\t\tmdata.attribute[%1d] = %g\n", i, attr[i]);
		}
	}

	/*
	 * first read all entries in the file 
	 */
	inla_read_data_all(&x, &n, ds->data_file.name, &ncol_data_all);
	assert(ncol_data_all <= _DIM_A);

	if (mb->verbose) {
		printf("\t\tread n=[%1d] entries from file=[%s]\n", n, ds->data_file.name);
	}
	switch (ds->data_id) {
	case L_GAUSSIAN:
	case L_STDGAUSSIAN:
	case L_LOGNORMAL:
	case L_EXPPOWER:
	{
		idiv = 3;
		a[0] = ds->data_observations.weight_gaussian = Calloc(mb->predictor_ndata, double);
	}
		break;

	case L_BC_GAUSSIAN:
	{
		idiv = 4;
		a[0] = ds->data_observations.bc_mean = Calloc(mb->predictor_ndata, double);
		a[1] = ds->data_observations.bc_scale = Calloc(mb->predictor_ndata, double);
	}
		break;

	case L_SEM:
	{
		idiv = 2;
	}
		break;

	case L_SIMPLEX:
	{
		idiv = 3;
		a[0] = ds->data_observations.weight_simplex = Calloc(mb->predictor_ndata, double);
	}
		break;

	case L_IID_GAMMA:
	{
		idiv = 3;
		a[0] = ds->data_observations.iid_gamma_scale = Calloc(mb->predictor_ndata, double);
	}
		break;

	case L_LOGISTIC:
	{
		idiv = 3;
		a[0] = ds->data_observations.weight_logistic = Calloc(mb->predictor_ndata, double);
	}
		break;

	case L_SKEWNORMAL:
	{
		idiv = 3;
		a[0] = ds->data_observations.sn_scale = Calloc(mb->predictor_ndata, double);
	}
		break;

	case L_GEV:
	{
		idiv = 3;
		a[0] = ds->data_observations.weight_gev = Calloc(mb->predictor_ndata, double);
	}
		break;

	case L_T:
	{
		idiv = 3;
		a[0] = ds->data_observations.weight_t = Calloc(mb->predictor_ndata, double);
	}
		break;

	case L_TSTRATA:
	{
		idiv = 4;
		a[0] = ds->data_observations.weight_tstrata = Calloc(mb->predictor_ndata, double);
		a[1] = ds->data_observations.strata_tstrata = Calloc(mb->predictor_ndata, double);
	}
		break;

	case L_CENPOISSON:
	case L_CONTPOISSON:
	case L_GAMMACOUNT:
	case L_GPOISSON:
	case L_POISSON:
	case L_NPOISSON:
	case L_NZPOISSON:
	case L_QCONTPOISSON:
	case L_XPOISSON:
	case L_ZEROINFLATEDCENPOISSON0:
	case L_ZEROINFLATEDCENPOISSON1:
	case L_ZEROINFLATEDNBINOMIAL0:
	case L_ZEROINFLATEDNBINOMIAL1:
	case L_ZEROINFLATEDNBINOMIAL2:
	case L_ZEROINFLATEDPOISSON0:
	case L_ZEROINFLATEDPOISSON1:
	case L_ZEROINFLATEDPOISSON2:
	case L_POISSON_SPECIAL1:
	case L_BELL:
	{
		idiv = 3;
		a[0] = ds->data_observations.E = Calloc(mb->predictor_ndata, double);
	}
		break;

	case L_NBINOMIAL:
	{
		idiv = 4;
		a[0] = ds->data_observations.E = Calloc(mb->predictor_ndata, double);
		a[1] = ds->data_observations.S = Calloc(mb->predictor_ndata, double);
	}
		break;

	case L_CENPOISSON2:
	{
		idiv = 5;
		a[0] = ds->data_observations.E = Calloc(mb->predictor_ndata, double);
		a[1] = ds->data_observations.cen_low = Calloc(mb->predictor_ndata, double);
		a[2] = ds->data_observations.cen_high = Calloc(mb->predictor_ndata, double);
	}
		break;

	case L_CENNBINOMIAL2:
	{
		idiv = 6;
		a[0] = ds->data_observations.E = Calloc(mb->predictor_ndata, double);
		a[1] = ds->data_observations.S = Calloc(mb->predictor_ndata, double);
		a[2] = ds->data_observations.cen_low = Calloc(mb->predictor_ndata, double);
		a[3] = ds->data_observations.cen_high = Calloc(mb->predictor_ndata, double);
	}
		break;

	case L_CBINOMIAL:
	{
		idiv = 4;
		a[0] = ds->data_observations.cbinomial_k = Calloc(mb->predictor_ndata, double);
		a[1] = ds->data_observations.cbinomial_n = Calloc(mb->predictor_ndata, double);
	}
		break;

	case L_GAUSSIANJW:
	{
		idiv = 5;
		a[0] = ds->data_observations.gjw_n = Calloc(mb->predictor_ndata, double);
		a[1] = ds->data_observations.gjw_df = Calloc(mb->predictor_ndata, double);
		a[2] = ds->data_observations.gjw_var = Calloc(mb->predictor_ndata, double);
	}
		break;

	case L_GAMMA:
	case L_MGAMMA:
	{
		idiv = 3;
		a[0] = ds->data_observations.gamma_scale = Calloc(mb->predictor_ndata, double);
	}
		break;

	case L_GAMMAJW:
	{
		idiv = 2;
		a[0] = NULL;
	}
		break;

	case L_DGP:
	case L_EXPONENTIAL:
	case L_GP:
	case L_IID_LOGITBETA:
	case L_LOGGAMMA_FRAILTY:
	case L_LOGLOGISTIC:
	case L_LOGPERIODOGRAM:
	case L_POM:
	case L_QKUMAR:
	case L_QLOGLOGISTIC:
	case L_STOCHVOL:
	case L_STOCHVOL_LN:
	case L_STOCHVOL_SN:
	case L_STOCHVOL_NIG:
	case L_STOCHVOL_T:
	case L_WEIBULL:
	case L_GOMPERTZ:
	case L_EGP:
	{
		idiv = 2;
		a[0] = NULL;
	}
		break;

	case L_BETA:
	{
		idiv = 3;
		a[0] = ds->data_observations.beta_weight = Calloc(mb->predictor_ndata, double);
	}
		break;

	case L_BETABINOMIALNA:
	{
		idiv = 4;
		a[0] = ds->data_observations.nb = Calloc(mb->predictor_ndata, double);
		a[1] = ds->data_observations.betabinomialnb_scale = Calloc(mb->predictor_ndata, double);
	}
		break;

	case L_XBINOMIAL:
	{
		idiv = 4;
		a[0] = ds->data_observations.nb = Calloc(mb->predictor_ndata, double);
		a[1] = ds->data_observations.p_scale = Calloc(mb->predictor_ndata, double);
	}
		break;

	case L_BETABINOMIAL:
	case L_BINOMIAL:
	case L_NBINOMIAL2:
	case L_ZEROINFLATEDBETABINOMIAL0:
	case L_ZEROINFLATEDBETABINOMIAL1:
	case L_ZEROINFLATEDBETABINOMIAL2:
	case L_ZEROINFLATEDBINOMIAL0:
	case L_ZEROINFLATEDBINOMIAL1:
	case L_ZEROINFLATEDBINOMIAL2:
	case L_ZERO_N_INFLATEDBINOMIAL2:
	case L_ZERO_N_INFLATEDBINOMIAL3:
	{
		idiv = 3;
		a[0] = ds->data_observations.nb = Calloc(mb->predictor_ndata, double);
	}
		break;

	case L_ZEROINFLATEDNBINOMIAL1STRATA2:
	case L_ZEROINFLATEDNBINOMIAL1STRATA3:
	{
		idiv = 4;
		a[0] = ds->data_observations.E = Calloc(mb->predictor_ndata, double);
		a[1] = ds->data_observations.strata = Calloc(mb->predictor_ndata, double);
	}
		break;

	case L_EXPONENTIALSURV:
	case L_GAMMASURV:
	case L_MGAMMASURV:
	case L_GAMMAJWSURV:
	case L_LOGLOGISTICSURV:
	case L_LOGNORMALSURV:
	case L_QLOGLOGISTICSURV:
	case L_WEIBULLSURV:
	case L_FMRISURV:
	case L_GOMPERTZSURV:
	{
		ds->data_observations.weight_gaussian = NULL;
		ds->data_observations.fmri_scale = NULL;

		assert(ncol_data_all <= 6 + CURE_MAXTHETA && ncol_data_all >= 6);
		idiv = ncol_data_all;
		na = ds->data_observations.cure_ncov = ncol_data_all - 6;
		assert(na >= 0);

		a[0] = ds->data_observations.event = Calloc(mb->predictor_ndata, double);	/* the failure code */
		a[1] = ds->data_observations.truncation = Calloc(mb->predictor_ndata, double);
		a[2] = ds->data_observations.lower = Calloc(mb->predictor_ndata, double);
		a[3] = ds->data_observations.upper = Calloc(mb->predictor_ndata, double);

		if (na) {
			// we'll wrap this around at the end of this function
			ds->data_observations.cure_cov = Calloc(mb->predictor_ndata * na, double);
			for (i = 0; i < na; i++) {
				a[4 + i] = ds->data_observations.cure_cov + i * mb->predictor_ndata;
			}
		}
	}
		break;

	case L_CIRCULAR_NORMAL:
	{
		idiv = 3;
		a[0] = ds->data_observations.weight_circular_normal = Calloc(mb->predictor_ndata, double);
	}
		break;

	case L_WRAPPED_CAUCHY:
	{
		idiv = 3;
		a[0] = ds->data_observations.weight_wrapped_cauchy = Calloc(mb->predictor_ndata, double);
	}
		break;

	case L_TWEEDIE:
	{
		idiv = 3;
		a[0] = ds->data_observations.tweedie_w = Calloc(mb->predictor_ndata, double);
	}
		break;

	case L_FMRI:
	{
		idiv = 3;
		a[0] = ds->data_observations.fmri_scale = Calloc(mb->predictor_ndata, double);
	}
		break;


	case L_FL:
	{
		// the 'fl_c' matrix is transposed at a later stage
		int m = L_FL_NC;
		idiv = m + 2;
		ds->data_observations.fl_c = Calloc(m, double *);
		for (k = 0; k < m; k++) {
			a[k] = ds->data_observations.fl_c[k] = Calloc(mb->predictor_ndata, double);
		}

		// only this is supported for the moment
		assert((ds->data_id == L_FL) && (GMRFLib_inla_mode == GMRFLib_MODE_COMPACT));
	}
		break;

	case L_NMIX:
	case L_NMIXNB:
	{
		int dim_y;
		// this case is a bit special, as the real data 'y' is fake, and the list
		// of replicated data is in the 'a' below.
		assert(ncol_data_all >= 3L + NMIX_MMAX && ncol_data_all < _DIM_A);
		idiv = ncol_data_all;
		ds->data_observations.nmix_x = Calloc(NMIX_MMAX, double *);
		dim_y = ncol_data_all - NMIX_MMAX - 2L;
		ds->data_observations.nmix_y = Calloc(dim_y + 1, double *);	/* yes, its +1 */
		for (i = 0; i < NMIX_MMAX; i++) {
			a[i] = ds->data_observations.nmix_x[i] = Calloc(mb->predictor_ndata, double);
		}
		for (i = 0; i < dim_y; i++) {
			a[i + NMIX_MMAX] = ds->data_observations.nmix_y[i] = Calloc(mb->predictor_ndata, double);
		}
		// fill the fake column of NA's so we know when to stop
		ds->data_observations.nmix_y[dim_y] = Calloc(mb->predictor_ndata, double);
		for (i = 0; i < mb->predictor_ndata; i++) {
			ds->data_observations.nmix_y[dim_y][i] = NAN;
		}
	}
		break;

	case L_OCCUPANCY:
	{
		int ny = (int) attr[0];
		int nx = (int) attr[1];
		int m = nx / ny;
		assert(ny + nx + 2 == ncol_data_all);
		idiv = ncol_data_all;

		ds->data_observations.occ_nbeta = m;
		ds->data_observations.occ_ny_max = ny;
		for (i = 0; i < ny; i++) {
			a[i] = Calloc(mb->predictor_ndata, double);
		}
		for (i = 0; i < nx; i++) {
			a[i + ny] = Calloc(mb->predictor_ndata, double);
		}
	}
		break;

	case L_BGEV:
	{
		assert(ncol_data_all <= 3 + BGEV_MAXTHETA && ncol_data_all >= 3);
		idiv = ncol_data_all;
		na = ncol_data_all - 2;
		ds->data_observations.bgev_x = Calloc(na, double *);
		a[0] = ds->data_observations.bgev_scale = Calloc(mb->predictor_ndata, double);
		for (i = 1; i < na; i++) {
			a[i] = ds->data_observations.bgev_x[i - 1] = Calloc(mb->predictor_ndata, double);
		}
	}
		break;

	case L_AGAUSSIAN:
	{
		int four = 4;
		assert(ncol_data_all == four + 2);
		idiv = ncol_data_all;
		ds->data_observations.agaussian = Calloc(four, double *);
		for (i = 0; i < four; i++) {
			ds->data_observations.agaussian[i] = Calloc(mb->predictor_ndata, double);
			a[i] = ds->data_observations.agaussian[i];
		}
	}
		break;

	case L_RCPOISSON:
	{
		assert(ncol_data_all <= 5 + RCPOISSON_MAXTHETA && ncol_data_all >= 5);
		idiv = ncol_data_all;
		int nbeta = ncol_data_all - 5;
		na = 3 + nbeta;
		ds->data_observations.rcp_nbeta = nbeta;
		a[0] = ds->data_observations.rcp_E = Calloc(mb->predictor_ndata, double);
		a[1] = ds->data_observations.rcp_event = Calloc(mb->predictor_ndata, double);
		a[2] = ds->data_observations.rcp_offset = Calloc(mb->predictor_ndata, double);
		if (nbeta > 0) {
			ds->data_observations.rcp_x = Calloc(nbeta, double *);
			for (i = 0; i < nbeta; i++) {
				a[3 + i] = ds->data_observations.rcp_x[i] = Calloc(mb->predictor_ndata, double);
			}
		} else {
			ds->data_observations.rcp_x = NULL;
		}
	}
		break;

	case L_TPOISSON:
	{
		assert(ncol_data_all <= 5 + TPOISSON_MAXTHETA && ncol_data_all >= 5);
		idiv = ncol_data_all;
		int nbeta = ncol_data_all - 5;
		na = 3 + nbeta;
		ds->data_observations.tp_nbeta = nbeta;
		a[0] = ds->data_observations.tp_E = Calloc(mb->predictor_ndata, double);
		a[1] = ds->data_observations.tp_event = Calloc(mb->predictor_ndata, double);
		a[2] = ds->data_observations.tp_offset = Calloc(mb->predictor_ndata, double);
		if (nbeta > 0) {
			ds->data_observations.tp_x = Calloc(nbeta, double *);
			for (i = 0; i < nbeta; i++) {
				a[3 + i] = ds->data_observations.tp_x[i] = Calloc(mb->predictor_ndata, double);
			}
		} else {
			ds->data_observations.tp_x = NULL;
		}
	}
		break;

	case L_GGAUSSIAN:
	{
		assert(ncol_data_all <= 3 + GGAUSSIAN_MAXTHETA && ncol_data_all >= 3);
		idiv = ncol_data_all;
		na = ncol_data_all - 2;
		ds->data_observations.ggaussian_nbeta = na - 1;
		ds->data_observations.ggaussian_x = Calloc(na, double *);
		a[0] = ds->data_observations.ggaussian_scale = Calloc(mb->predictor_ndata, double);
		for (i = 1; i < na; i++) {
			a[i] = ds->data_observations.ggaussian_x[i - 1] = Calloc(mb->predictor_ndata, double);
		}
	}
		break;

	case L_GGAUSSIANS:
	{
		assert(ncol_data_all <= 3 + GGAUSSIAN_MAXTHETA && ncol_data_all >= 3);
		idiv = ncol_data_all;
		na = ncol_data_all - 2;
		ds->data_observations.ggaussian_nbeta = na - 1;
		ds->data_observations.ggaussian_x = Calloc(na, double *);
		a[0] = ds->data_observations.ggaussian_offset = Calloc(mb->predictor_ndata, double);
		for (i = 1; i < na; i++) {
			a[i] = ds->data_observations.ggaussian_x[i - 1] = Calloc(mb->predictor_ndata, double);
		}
	}
		break;

	case L_0POISSON:
	case L_0POISSONS:
	{
		assert(ncol_data_all <= 3 + POISSON0_MAXTHETA && ncol_data_all >= 3);
		idiv = ncol_data_all;
		na = ncol_data_all - 2;
		ds->data_observations.poisson0_nbeta = na - 1;
		ds->data_observations.poisson0_x = Calloc(na, double *);
		a[0] = ds->data_observations.poisson0_E = Calloc(mb->predictor_ndata, double);
		for (i = 1; i < na; i++) {
			a[i] = ds->data_observations.poisson0_x[i - 1] = Calloc(mb->predictor_ndata, double);
		}
	}
		break;

	case L_0BINOMIAL:
	case L_0BINOMIALS:
	{
		assert(ncol_data_all <= 3 + BINOMIAL0_MAXTHETA && ncol_data_all >= 3);
		idiv = ncol_data_all;
		na = ncol_data_all - 2;
		ds->data_observations.binomial0_nbeta = na - 1;
		ds->data_observations.binomial0_x = Calloc(na, double *);
		a[0] = ds->data_observations.binomial0_Ntrials = Calloc(mb->predictor_ndata, double);
		for (i = 1; i < na; i++) {
			a[i] = ds->data_observations.binomial0_x[i - 1] = Calloc(mb->predictor_ndata, double);
		}
	}
		break;

	case L_BINOMIALMIX:
	{
		assert(ncol_data_all == 16);
		idiv = ncol_data_all;
		na = ncol_data_all - 2;
		for (i = 0; i < na; i++) {
			a[i] = Calloc(mb->predictor_ndata, double);
		}
	}
		break;

	default:
		GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
	}

	na = idiv - 2;

	if (!inla_divisible(n, idiv)) {
		inla_error_file_numelm(__GMRFLib_FuncName, ds->data_file.name, n, idiv);
	}
	ds->data_observations.ndata = n / idiv;
	ds->data_observations.y = Calloc(mb->predictor_ndata, double);
	ds->data_observations.d = Calloc(mb->predictor_ndata, double);

	double *w = NULL;
	int nw = 0;
	inla_read_data_all(&w, &nw, ds->weight_file.name, NULL);
	if (nw) {
		// P(nw); P(ds->data_observations.ndata); P(mb->predictor_ndata);
		assert(nw == mb->predictor_ndata);
	}

	double *lp_scale = NULL;
	int n_lp_scale = 0;
	inla_read_data_all(&lp_scale, &n_lp_scale, ds->lp_scale_file.name, NULL);
	if (n_lp_scale) {
		assert(n_lp_scale == mb->predictor_ndata);
	}
	for (i = 0; i < n_lp_scale; i++) {
		lp_scale[i] = (int) (lp_scale[i] - 1.0);
	}
	mb->data_sections[0].lp_scale = lp_scale;

	for (i = j = 0; i < n; i += idiv, j++) {
		ii = (int) x[i];
		if (!LEGAL(ii, mb->predictor_ndata)) {
			inla_error_file_error(__GMRFLib_FuncName, ds->data_file.name, n, i, x[i]);
		}
		for (k = 0; k < na; k++) {
			a[k][ii] = x[i + k + 1];
		}
		ds->data_observations.y[ii] = x[i + idiv - 1];
		if (w) {
			ds->data_observations.d[ii] = w[ii];
		} else {
			ds->data_observations.d[ii] = 1.0;
		}
		if (mb->verbose && j < PREVIEW) {
			switch (na) {
			case 0:
			{
				printf("\t\t\t%1d/%1d  (idx,y,d) = (%1d, %g, %g)\n", j, ds->data_observations.ndata, ii,
				       ds->data_observations.y[ii], ds->data_observations.d[ii]);
			}
				break;

			case 1:
			{
				printf("\t\t\t%1d/%1d  (idx,a,y,d) = (%1d, %g, %g, %g)\n", j,
				       ds->data_observations.ndata, ii, a[0][ii], ds->data_observations.y[ii], ds->data_observations.d[ii]);
			}
				break;

			case 2:
			{
				printf("\t\t\t%1d/%1d (idx,a[0],a[1],y,d) = (%1d, %g, %g, %g, %g)\n", j,
				       ds->data_observations.ndata, ii, a[0][ii], a[1][ii], ds->data_observations.y[ii],
				       ds->data_observations.d[ii]);
			}
				break;

			case 3:
			{
				printf("\t\t\t%1d/%1d (idx,a[0],a[1],a[2],y,d) = (%1d, %g, %g, %g, %g, %g)\n", j,
				       ds->data_observations.ndata, ii, a[0][ii], a[1][ii], a[2][ii], ds->data_observations.y[ii],
				       ds->data_observations.d[ii]);
			}
				break;

			case 4:
			{
				printf("\t\t\t%1d/%1d (idx,a[0],a[1],a[2],a[3],y,d) = (%1d, %g, %g, %g, %g, %g, %g)\n", j,
				       ds->data_observations.ndata, ii, a[0][ii], a[1][ii], a[2][ii], a[3][ii],
				       ds->data_observations.y[ii], ds->data_observations.d[ii]);
			}
				break;

			case 5:
			{
				printf("\t\t\t%1d/%1d (idx,a[0],a[1],a[2],a[3],a[4],y,d) = (%1d, %g, %g, %g, %g, %g, %g, %g)\n", j,
				       ds->data_observations.ndata, ii, a[0][ii], a[1][ii], a[2][ii], a[3][ii], a[4][ii],
				       ds->data_observations.y[ii], ds->data_observations.d[ii]);
			}
				break;

			default:
				printf("\t\t\t%1d/%1d (idx,a[],y,d) = (%1d, ", j, ds->data_observations.ndata, ii);
				for (k = 0; k < na; k++) {
					printf("%6.4g, ", a[k][ii]);
				}
				printf("%g, %g)\n", ds->data_observations.y[ii], ds->data_observations.d[ii]);
				break;
			}
		}
	}

	if (ds->data_id == L_BINOMIALMIX) {
		// re-arrange data so they become more accessible
		ds->data_observations.binmix_dat = Calloc(mb->predictor_ndata, double *);
		for (i = 0; i < mb->predictor_ndata; i++) {
			ds->data_observations.binmix_dat[i] = Calloc(na, double);
			for (k = 0; k < na; k++) {
				ds->data_observations.binmix_dat[i][k] = a[k][i];
			}
			char *msg = NULL;
			if (!(ds->data_observations.y[i] >= 0 && ds->data_observations.y[i] <= ds->data_observations.binmix_dat[i][na - 1])) {
				GMRFLib_sprintf(&msg, "binomialmix Ntrials[%1d] = %g y[%1d] = %g is void\n", i,
						ds->data_observations.y[i], ds->data_observations.binmix_dat[i][na - 1]);
				inla_error_general(msg);
			}
			if ((double) ((int) ds->data_observations.y[i]) != ds->data_observations.y[i]) {
				GMRFLib_sprintf(&msg, "binomialmix y[%1d] = %g is not integer\n", i, ds->data_observations.y[i]);
				inla_error_general(msg);
			}
		}
		for (i = 0; i < na; i++) {
			Free(a[i]);
		}
	}
	// wrap it around so we can access all cure-covariates for one observation sequentially
	if (ds->data_observations.cure_cov) {
		int ncov = ds->data_observations.cure_ncov;
		n = mb->predictor_ndata;
		double *xx = Calloc(ncov * n, double);
		for (j = 0; j < ncov; j++) {
			for (i = 0; i < n; i++) {
				xx[i * ncov + j] = ds->data_observations.cure_cov[j * n + i];
			}
		}
		Free(ds->data_observations.cure_cov);
		ds->data_observations.cure_cov = xx;
	}
	// rearrange the data and covariates
	if (ds->data_id == L_OCCUPANCY) {
		int nb = ds->data_observations.occ_nbeta;
		int ny_max = ds->data_observations.occ_ny_max;
		int nd = mb->predictor_ndata;

		double **X = ds->data_observations.occ_x = Calloc(nd, double *);
		int **Y = ds->data_observations.occ_y = Calloc(nd, int *);
		int *ny = ds->data_observations.occ_ny = Calloc(nd, int);

		for (i = 0; i < nd; i++) {

			// count the number of observations
			int nyy = 0;
			for (j = 0; j < ny_max; j++) {
				double yy = a[j][i];
				if (!ISNAN(yy)) {
					nyy++;
				}
			}
			if (nyy == 0) {
				ny[i] = 0;
			} else {
				ny[i] = nyy;
				X[i] = Calloc(nb * nyy, double);
				Y[i] = Calloc(nyy, int);

				int jj = 0;
				for (j = 0; j < ny_max; j++) {
					double yy = a[j][i];
					if (!ISNAN(yy)) {
						Y[i][jj] = (int) yy;
						for (k = 0; k < nb; k++) {
							double xx = a[ny_max + j * nb + k][i];
							if (ISNAN(xx)) {
								xx = 0.0;
							}
							X[i][jj * nb + k] = xx;
						}
						jj++;
					}
				}
			}
		}
		for (j = 0; j < (1 + nb) * ny_max; j++) {
			Free(a[j]);
		}

		int *yzero = ds->data_observations.occ_yzero = Calloc(nd, int);
		for (i = 0; i < nd; i++) {
			k = 1;
			for (j = 0; j < ny[i]; j++) {
				k = (k && ISZERO(Y[i][j]));
			}
			yzero[i] = k;
		}

		// check
		for (i = 0; i < nd; i++) {
			if (ds->data_observations.d[i]) {
				for (int kk = 0; kk < ds->data_observations.occ_ny[i]; kk++) {
					if ((ds->data_observations.occ_y[i][kk] < 0) || (ds->data_observations.occ_y[i][kk] > 1)) {
						char *msg = NULL;
						GMRFLib_sprintf(&msg, "occupancy observation y[%1d,%1d] = %d is void\n", i, kk,
								ds->data_observations.occ_y[i][kk]);
						inla_error_general(msg);
					}
				}
			}
		}
	}

	Free(w);
	Free(x);
#undef _DIM_A

	return INLA_OK;
}

int loglikelihood_inla(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *x_vec, double *y_cdf, void *arg,
		       char **arg_str)
{
	inla_tp *a = (inla_tp *) arg;
	return a->loglikelihood[idx] (thread_id, logll, x, m, idx, x_vec, y_cdf, a->loglikelihood_arg[idx], arg_str);
}

double inla_dnchisq(double x, double df, double ncp)
{
	// code provided by L.Starke

	double ldens = 0.0;
	if (ISZERO(ncp)) {
		// besselI is not defined for x = 0
		ldens = MATHLIB_FUN(dnchisq) (x, df, ncp, 1);
	} else if (sqrt(ncp * x) < 8.0E4) {
		// the cutoff is due to besselI
		// see xlrg_IJ at https://github.com/atks/Rmath/blob/master/bessel.h
		// 
		// alternative form of pdf is used
		// https://en.wikipedia.org/wiki/Noncentral_chi-squared_distribution
		// 
		// note for C implementation with Rmath: 
		// besselI(..., ..., TRUE) -> bessel_i(..., ..., 2)
		double c1 = sqrt(x * ncp);
		ldens = -M_LN2 - (x + ncp) / 2.0 + (df / 4.0 - 0.5) * log(x / ncp) + c1 + log(MATHLIB_FUN(bessel_i) (c1, df / 2.0 - 1.0, 2.0));
	} else {
		// Approximation of Fraser (is perfect on right tail)
		double c1 = log(x / ncp);
		double c2 = sqrt(x) - sqrt(ncp);
		ldens = -0.5 * log(8.0 * M_PI * x) -
		    0.5 * SQR(c2 - (df - 1.0) / 4.0 * c1 / c2) + LOG_1mp((df - 1.0) / (2.0 * c2) * (1.0 / sqrt(x) - 0.5 * c1 / c2));
	}

	return (ldens);
}

int loglikelihood_gaussian(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			   void *arg, char **arg_str)
{
	/*
	 * y ~ Normal(x, stdev)
	 */
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}
	Data_section_tp *ds = (Data_section_tp *) arg;
	static double log_prec_limit = -log(INLA_REAL_SMALL);

	double y = ds->data_observations.y[idx];
	double w = ds->data_observations.weight_gaussian[idx];
	double lprec, prec;

	LINK_INIT;

	if (ds->data_observations.log_prec_gaussian_offset[thread_id][0] > log_prec_limit) {
		lprec = ds->data_observations.log_prec_gaussian[thread_id][0] + log(w);
		prec = exp(lprec);
	} else {
		double prec_offset = map_precision_forward(ds->data_observations.log_prec_gaussian_offset[thread_id][0], MAP_FORWARD, NULL);
		double prec_var = map_precision_forward(ds->data_observations.log_prec_gaussian[thread_id][0], MAP_FORWARD, NULL);
		double prec_tmp = 1.0 / (1.0 / prec_offset + 1.0 / prec_var);
		prec = prec_tmp * w;
		lprec = log(prec);
	}

	if (arg_str) {
		char *a = NULL;
		GMRFLib_sprintf(&a,
				"list(y = %.8g, family = \"gaussian\", scale = %.8g, link.model = \"%1s\", theta = %.8g, linear.predictor = %.8g)",
				y, w, ds->link_model, ds->data_observations.log_prec_gaussian[thread_id][0], x[0]);
		*arg_str = a;
		return GMRFLib_SUCCESS;
	}

	if (m > 0) {
		if (PREDICTOR_LINK_EQ(link_identity)) {
			if (PREDICTOR_SCALE == 1.0) {
				double a = -0.5 * prec;
				double b = LOG_NORMC_GAUSSIAN + 0.5 * lprec;
#pragma omp simd
				for (int i = 0; i < m; i++) {
					double res = y - (x[i] + off);
					logll[i] = b + a * SQR(res);
				}
			} else {
#pragma omp simd
				for (int i = 0; i < m; i++) {
					double ypred = PREDICTOR_INVERSE_IDENTITY_LINK(x[i] + off);
					logll[i] = LOG_NORMC_GAUSSIAN + 0.5 * (lprec - (SQR(ypred - y) * prec));
				}
			}
		} else {
#pragma omp simd
			for (int i = 0; i < m; i++) {
				double ypred = PREDICTOR_INVERSE_LINK(x[i] + off);
				logll[i] = LOG_NORMC_GAUSSIAN + 0.5 * (lprec - (SQR(ypred - y) * prec));
			}
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
		for (int i = 0; i < -m; i++) {
			double ypred = PREDICTOR_INVERSE_LINK(x[i] + off);
			logll[i] = inla_cdf_normal_fast((y - ypred) * sqrt(prec));
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_stdgaussian(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			      void *arg, char **arg_str)
{
	/*
	 * y ~ Normal(x, 1)
	 */
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}
	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx];
	double w = ds->data_observations.weight_gaussian[idx];
	double prec = w;
	double lprec = (w == 1.0 ? 0.0 : log(prec));

	LINK_INIT;
	if (arg_str) {
		char *a = NULL;
		GMRFLib_sprintf(&a,
				"list(y = %.8g, family = \"stdgaussian\", scale = %.8g, link.model = \"%1s\", theta = %.8g, linear.predictor = %.8g)",
				y, w, ds->link_model, ds->data_observations.log_prec_gaussian[thread_id][0], x[0]);
		*arg_str = a;
		return GMRFLib_SUCCESS;
	}

	if (m > 0) {
		if (PREDICTOR_LINK_EQ(link_identity)) {

			if (PREDICTOR_LINK_EQ(link_identity) && (PREDICTOR_SCALE == 1.0)) {
				double a = -0.5 * prec;
				double b = LOG_NORMC_GAUSSIAN + 0.5 * lprec;
#pragma omp simd
				for (int i = 0; i < m; i++) {
					double res = y - (x[i] + off);
					logll[i] = b + a * SQR(res);
				}
			} else {
#pragma omp simd
				for (int i = 0; i < m; i++) {
					double ypred = PREDICTOR_INVERSE_IDENTITY_LINK(x[i] + off);
					logll[i] = LOG_NORMC_GAUSSIAN + 0.5 * (lprec - (SQR(ypred - y) * prec));
				}
			}
		} else {
#pragma omp simd
			for (int i = 0; i < m; i++) {
				double ypred = PREDICTOR_INVERSE_LINK(x[i] + off);
				logll[i] = LOG_NORMC_GAUSSIAN + 0.5 * (lprec - (SQR(ypred - y) * prec));
			}
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
		for (int i = 0; i < -m; i++) {
			double ypred = PREDICTOR_INVERSE_LINK(x[i] + off);
			logll[i] = inla_cdf_normal_fast((y - ypred) * sqrt(prec));
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_exppower(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			   void *arg, char **UNUSED(arg_str))
{
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}
	Data_section_tp *ds = (Data_section_tp *) arg;

	typedef struct {
		double beta;
		double lgamma1;
		double lgamma3;
	} lcache_t;

	static lcache_t **llcache = NULL;
	if (!llcache) {
#pragma omp critical (Name_21bf489d4a1518aff583c47b3b79b92dc0329b87)
		if (!llcache) {
			llcache = Calloc(GMRFLib_CACHE_LEN(), lcache_t *);
		}
	}

	int cidx = 0;
	GMRFLib_CACHE_SET_ID(cidx);
	if (!llcache[cidx]) {
#pragma omp critical (Name_93af423a814c85a479569f0787bff31c76ef23bf)
		if (!llcache[cidx]) {
			llcache[cidx] = Calloc(1, lcache_t);
		}
	}

	double y = ds->data_observations.y[idx];
	double w = ds->data_observations.weight_gaussian[idx];
	double beta = map_one_plus_exp(ds->data_observations.log_power[thread_id][0], MAP_FORWARD, NULL);
	double lprec = ds->data_observations.log_prec_gaussian[thread_id][0] + log(w);
	double sigma = exp(-0.5 * lprec);

	lcache_t *lc = llcache[cidx];
	if (lc->beta != beta) {
		lc->beta = beta;
		lc->lgamma1 = my_gsl_sf_lngamma(1.0 / beta);
		lc->lgamma3 = my_gsl_sf_lngamma(3.0 / beta);
	}
	double alpha = sigma * exp(0.5 * (lc->lgamma1 - lc->lgamma3));
	double ialpha = 1.0 / alpha;

	LINK_INIT;
	if (m > 0) {
		double logZ = log(beta / (2.0 * alpha)) - lc->lgamma1;
		if (PREDICTOR_LINK_EQ(link_identity)) {
#pragma omp simd
			for (int i = 0; i < m; i++) {
				double ypred = PREDICTOR_INVERSE_IDENTITY_LINK(x[i] + off);
				double z = ABS(y - ypred) * ialpha;
				logll[i] = logZ - pow(z, beta);
			}
		} else {
			for (int i = 0; i < m; i++) {
				double ypred = PREDICTOR_INVERSE_LINK(x[i] + off);
				double z = ABS(y - ypred) * ialpha;
				logll[i] = logZ - pow(z, beta);
			}
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
		double shape = 1.0 / beta;
		double scale = 1.0 / pow(ialpha, beta);
		for (int i = 0; i < -m; i++) {
			double mu = PREDICTOR_INVERSE_LINK(x[i] + off);
			logll[i] = 0.5 * (1.0 + DSIGN(y - mu) * MATHLIB_FUN(pgamma) (pow(ABS(y - mu), beta), shape, scale, 1, 0));
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_sem(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
		      void *arg, char **UNUSED(arg_str))
{
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx];

	LINK_INIT;
	double prec = inla_eval_param_constraint(thread_id, ds);
	double lprec = log(prec);

	if (m > 0) {
		if (PREDICTOR_LINK_EQ(link_identity)) {
			double a = -0.5 * prec;
			double b = LOG_NORMC_GAUSSIAN + 0.5 * lprec;
#pragma omp simd
			for (int i = 0; i < m; i++) {
				double res = y - (x[i] + off);
				logll[i] = b + a * SQR(res);
			}
		} else {
#pragma omp simd
			for (int i = 0; i < m; i++) {
				double ypred = PREDICTOR_INVERSE_LINK(x[i] + off);
				logll[i] = LOG_NORMC_GAUSSIAN + 0.5 * (lprec - (SQR(ypred - y) * prec));
			}
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
		for (int i = 0; i < -m; i++) {
			double ypred = PREDICTOR_INVERSE_LINK(x[i] + off);
			logll[i] = inla_cdf_normal_fast((y - ypred) * sqrt(prec));
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_gaussianjw(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec),
			     double *UNUSED(y_cdf), void *arg, char **UNUSED(arg_str))
{
	if (m == 0) {
		return GMRFLib_SUCCESS;
	}
	Data_section_tp *ds = (Data_section_tp *) arg;

	double y = ds->data_observations.y[idx];
	double var_obs = ds->data_observations.gjw_var[idx];
	double log_n = log(ds->data_observations.gjw_n[idx]);
	double df = ds->data_observations.gjw_df[idx];
	double df2 = df / 2.0;
	double beta_0 = ds->data_observations.gjw_beta[0][thread_id][0];
	double beta_1 = ds->data_observations.gjw_beta[1][thread_id][0];
	double beta_2 = ds->data_observations.gjw_beta[2][thread_id][0];

	if (G_norm_const_compute[idx]) {
		G_norm_const[idx] = LOG_NORMC_GAUSSIAN - df2 * M_LN2 - gsl_sf_lngamma(df2);
		G_norm_const_compute[idx] = 0;
	}
	double normc = G_norm_const[idx];

	LINK_INIT;
	if (m > 0) {
#pragma omp simd
		for (int i = 0; i < m; i++) {
			double p = PREDICTOR_INVERSE_LINK(x[i] + off);
			double var = exp(beta_0 + beta_1 * log(p * (1.0 - p)) + beta_2 * log_n);
			double prec = 1.0 / var;
			logll[i] = normc + 0.5 * (log(prec) - (SQR(p - y) * prec));

			double chi_sqr = df * var_obs * prec;
			logll[i] += (df2 - 1.0) * log(chi_sqr) - chi_sqr / 2.0;
		}
	} else {
		GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_agaussian(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			    void *arg, char **UNUSED(arg_str))
{
	/*
	 * aggregated Gaussian
	 */
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}
	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double lprec, prec, ypred;
	double y, v, ldet_s, mm, nn;

	lprec = ds->data_observations.log_prec_gaussian[thread_id][0];
	prec = map_precision_forward(ds->data_observations.log_prec_gaussian[thread_id][0], MAP_FORWARD, NULL);
	y = ds->data_observations.y[idx];
	v = ds->data_observations.agaussian[0][idx];
	ldet_s = ds->data_observations.agaussian[1][idx];
	mm = ds->data_observations.agaussian[2][idx];
	nn = ds->data_observations.agaussian[3][idx];

	LINK_INIT;
	if (m > 0) {
		for (i = 0; i < m; i++) {
			ypred = PREDICTOR_INVERSE_LINK(x[i] + off);
			logll[i] = nn * LOG_NORMC_GAUSSIAN + nn / 2.0 * lprec + ldet_s - 0.5 * mm * prec * (SQR(y - ypred) + v);
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
		for (i = 0; i < -m; i++) {
			ypred = PREDICTOR_INVERSE_LINK(x[i] + off);
			logll[i] = inla_cdf_normal_fast((y - ypred) * mm * prec);
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_ggaussian(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			    void *arg, char **UNUSED(arg_str))
{
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	Data_section_tp *ds = (Data_section_tp *) arg;
	Data_tp *dtp = &(ds->data_observations);

	double y, w;
	y = dtp->y[idx];
	w = dtp->ggaussian_scale[idx];

	LINK_INIT;

	double lprec = 0.0;
	for (int i = 0; i < dtp->ggaussian_nbeta; i++) {
		lprec += dtp->ggaussian_beta[i][thread_id][0] * dtp->ggaussian_x[i][idx];
	}
	double prec = w * dtp->link_simple_invlinkfunc(thread_id, lprec, MAP_FORWARD, NULL, NULL);
	lprec = log(prec);

	if (m > 0) {
		double a = LOG_NORMC_GAUSSIAN + 0.5 * lprec;
		double b = -0.5 * prec;
		if (PREDICTOR_LINK_EQ(link_identity) && (PREDICTOR_SCALE == 1.0)) {
#pragma omp simd
			for (int i = 0; i < m; i++) {
				double err = (x[i] + off) - y;
				logll[i] = a + b * SQR(err);
			}
		} else {
#pragma omp simd
			for (int i = 0; i < m; i++) {
				double err = PREDICTOR_INVERSE_LINK(x[i] + off) - y;
				logll[i] = a + b * SQR(err);
			}
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
		double sprec = sqrt(prec);
		for (int i = 0; i < -m; i++) {
			double ypred = PREDICTOR_INVERSE_LINK(x[i] + off);
			logll[i] = inla_cdf_normal_fast((y - ypred) * sprec);
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_ggaussianS(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			     void *arg, char **UNUSED(arg_str))
{
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	Data_section_tp *ds = (Data_section_tp *) arg;
	Data_tp *dtp = &(ds->data_observations);

	double y = dtp->y[idx];

	LINK_INIT;

	double mean = dtp->ggaussian_offset[idx];
	for (int i = 0; i < dtp->ggaussian_nbeta; i++) {
		mean += dtp->ggaussian_beta[i][thread_id][0] * dtp->ggaussian_x[i][idx];
	}
	if (!PREDICTOR_SIMPLE_LINK_EQ(link_identity)) {
		mean = dtp->link_simple_invlinkfunc(thread_id, mean, MAP_FORWARD, NULL, NULL);
	}

	if (m > 0) {
		double err2 = SQR(y - mean);
		if (PREDICTOR_LINK_EQ(link_log)) {
#pragma omp simd
			for (int i = 0; i < m; i++) {
				double lprec = x[i] + off;
				logll[i] = LOG_NORMC_GAUSSIAN + 0.5 * lprec - 0.5 * exp(lprec) * err2;
			}
		} else {
#pragma omp simd
			for (int i = 0; i < m; i++) {
				double prec = PREDICTOR_INVERSE_LINK(x[i] + off);
				logll[i] = LOG_NORMC_GAUSSIAN + 0.5 * log(prec) - 0.5 * prec * err2;
			}
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
#pragma omp simd
		for (int i = 0; i < -m; i++) {
			double prec = PREDICTOR_INVERSE_LINK(x[i] + off);
			logll[i] = inla_cdf_normal_fast((y - mean) * sqrt(prec));
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_lognormal(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			    void *arg, char **UNUSED(arg_str))
{
	/*
	 * y ~ LogNormal. This is similar to the normal
	 */
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}
	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double ly, lprec, prec, w, lw, lypred;

	double *cache = NULL;
	if (G_norm_const_compute[idx]) {
		cache = Calloc(3, double);
		G_norm_const_v[idx] = (void *) cache;
		G_norm_const_compute[idx] = 0;

		// log(y)
		cache[0] = log(ds->data_observations.y[idx]);
		// w
		cache[1] = (ds->data_observations.weight_gaussian ? ds->data_observations.weight_gaussian[idx] : 1.0);
		// log(w)
		cache[2] = log(cache[1]);
	}
	cache = (double *) G_norm_const_v[idx];
	ly = cache[0];
	w = cache[1];
	lw = cache[2];

	lprec = ds->data_observations.log_prec_gaussian[thread_id][0] + lw;
	prec = map_precision_forward(ds->data_observations.log_prec_gaussian[thread_id][0], MAP_FORWARD, NULL) * w;

	LINK_INIT;
	if (m > 0) {
		for (i = 0; i < m; i++) {
			lypred = PREDICTOR_INVERSE_LINK(x[i] + off);
			logll[i] = LOG_NORMC_GAUSSIAN + 0.5 * (lprec - (SQR(lypred - ly) * prec)) - ly;
		}
	} else {
		if (y_cdf)
			ly = log(*y_cdf);
		for (i = 0; i < -m; i++) {
			lypred = PREDICTOR_INVERSE_LINK(x[i] + off);
			logll[i] = inla_cdf_normal_fast((ly - lypred) * sqrt(prec));
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_lognormalsurv(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *x_vec, double *y_cdf,
				void *arg, char **arg_str)
{
	return (m ==
		0 ? GMRFLib_SUCCESS : loglikelihood_generic_surv(thread_id, logll, x, m, idx, x_vec, y_cdf, arg, loglikelihood_lognormal, arg_str));
}

int loglikelihood_bcgaussian(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			     void *arg, char **UNUSED(arg_str))
{
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}
	Data_section_tp *ds = (Data_section_tp *) arg;
	double yo, y, lprec, prec, w, lambda, mean;

	FIXME("BCGAUSSIAN: THIS IS NOT YET DONE AND I DO NOT KNOW IF THIS WAY OF DOING IS CORRECT, EVEN THOUGH ITS WHAT HAS BEEN DONE EARLIER....");
	exit(1);

	lambda = ds->data_observations.bc_lambda[thread_id][0];
	mean = ds->data_observations.bc_mean[idx];
	yo = ds->data_observations.y[idx];
	y = inla_boxcox(yo, mean, lambda);
	w = ds->data_observations.bc_scale[idx];
	lprec = ds->data_observations.log_prec_gaussian[thread_id][0] + log(w) - 2.0 * (lambda - 1.0) * log(mean);
	prec = exp(lprec);

	LINK_INIT;

	if (m > 0) {
		double lcorr = LOG_NORMC_GAUSSIAN + (lambda - 1.0) * log(yo);
#pragma omp simd
		for (int i = 0; i < m; i++) {
			double ypred = PREDICTOR_INVERSE_LINK(x[i] + off);
			logll[i] = lcorr + 0.5 * (lprec - (SQR(ypred - y) * prec));
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
		for (int i = 0; i < -m; i++) {
			double ypred = PREDICTOR_INVERSE_LINK(x[i] + off);
			logll[i] = inla_cdf_normal_fast((y - ypred) * sqrt(prec));
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_fl(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *UNUSED(y_cdf),
		     void *arg, char **UNUSED(arg_str))
{
	// return c[0] + c[1] * x - 1/2 * c[2] * (c[3] - x)^2 - c[4] exp(c[5] + c[6] * x) - c[7] * log((exp(c[8]*x)-1.0)/(sign(c[8])x))

	if (m == 0) {
		return GMRFLib_SUCCESS;
	}

	Data_section_tp *ds = (Data_section_tp *) arg;
	LINK_INIT;
	if (m > 0) {
		double *c = ds->data_observations.fl_c[idx];

		if (0) {
			for (int i = 0; i < m; i++) {
				double eta = PREDICTOR_INVERSE_LINK(x[i] + off);
				logll[i] = c[0] + c[1] * eta - 0.5 * c[2] * SQR(c[3] - eta) - c[4] * exp(c[5] + c[6] * eta)
				    + c[7] * expm1(c[8] * eta) / PUSH_AWAY(eta);
				// + c[7] * (exp(c[8] * eta) - 1.0) / PUSH_AWAY(eta);
			}
		}

		size_t mm = GMRFLib_align_simple((size_t) m, sizeof(double));
		double eta[mm];
		for (int i = 0; i < m; i++) {
			eta[i] = PREDICTOR_INVERSE_LINK(x[i] + off);
			logll[i] = c[0] + c[1] * eta[i];
		}

		if (!ISZERO(c[2])) {
#pragma omp simd
			for (int i = 0; i < m; i++) {
				logll[i] += (-0.5 * c[2] * SQR(c[3] - eta[i]));
			}
		}

		if (!ISZERO(c[4])) {
#pragma omp simd
			for (int i = 0; i < m; i++) {
				logll[i] += (-c[4] * exp(c[5] + c[6] * eta[i]));
			}
		}

		if (!ISZERO(c[7])) {
			double sign = DSIGN(c[8]);
			for (int i = 0; i < m; i++) {
				logll[i] += (-c[7] * log(expm1(c[8] * eta[i]) / (sign * PUSH_AWAY(eta[i]))));
			}
		}
	} else {
		assert(0 == 1);
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_simplex(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec),
			  double *UNUSED(y_cdf), void *arg, char **UNUSED(arg_str))
{
	/*
	 * y ~ simplex
	 */
	if (m == 0) {
		return GMRFLib_SUCCESS;
	}
	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double y, yy, ypyp, lprec, prec, w, ypred;

	y = ds->data_observations.y[idx];
	w = ds->data_observations.weight_simplex[idx];
	lprec = ds->data_observations.log_prec_simplex[thread_id][0] + log(w);
	prec = map_precision_forward(ds->data_observations.log_prec_simplex[thread_id][0], MAP_FORWARD, NULL) * w;

	LINK_INIT;
	if (m > 0) {
		for (i = 0; i < m; i++) {
			ypred = PREDICTOR_INVERSE_LINK(x[i] + off);
			yy = y * (1.0 - y);
			ypyp = ypred * (1.0 - ypred);
			logll[i] = 0.5 * lprec - 0.5 * (log(2.0 * M_PI) + 3.0 * log(yy))
			    - prec * SQR(y - ypred) / (2.0 * yy * SQR(ypyp));
		}
	} else {
		GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);

	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_circular_normal(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec),
				  double *UNUSED(y_cdf), void *arg, char **UNUSED(arg_str))
{

	// this needs to be redone....
	assert(0 == 1);

	/*
	 * y ~ circular normal
	 */
	if (m == 0) {
		return GMRFLib_SUCCESS;
	}

	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double y, prec, w, ypred;

	LINK_INIT;
	y = ds->data_observations.y[idx];
	w = ds->data_observations.weight_circular_normal[idx];
	prec = map_precision_forward(ds->data_observations.log_prec_circular_normal[thread_id][0], MAP_FORWARD, NULL) * w;

	/*
	 * store the normalising constant as it involves bessel_I0: -log(2 Pi BesselI0(kappa)),
	 * which is ok as long as the scalings 'w' do not change to often.
	 */
	double log_norm_const = 0.0, log_norm_const_arg = INLA_REAL_BIG;
	if (!ISEQUAL(prec, log_norm_const_arg)) {
		log_norm_const_arg = prec;
		log_norm_const = -log(2.0 * M_PI * gsl_sf_bessel_I0(prec));
	}

	for (i = 0; i < m; i++) {
		ypred = PREDICTOR_INVERSE_LINK(x[i] + off);

		/*
		 * we need |y-ypred| <= Pi, but this might not be the case...  so we add a penalty if this condition is not met
		 */
		if (ABS(y - ypred) <= M_PI) {
			logll[i] = log_norm_const + prec * cos(y - ypred);
		} else {
			double penalty = 1.0e8 * prec;
			if (y - ypred > M_PI) {
				logll[i] = log_norm_const + prec * cos(M_PI) - penalty / 2.0 * SQR(y - ypred - M_PI);
			} else {
				logll[i] = log_norm_const + prec * cos(-M_PI) - penalty / 2.0 * SQR(y - ypred + M_PI);
			}
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_wrapped_cauchy(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec),
				 double *UNUSED(y_cdf), void *arg, char **UNUSED(arg_str))
{
	/*
	 * y ~ wrapped cauchy. DOES NOT WORK WELL OF'COURSE...
	 */
	if (m == 0) {
		return GMRFLib_SUCCESS;
	}

	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double y, rho, rho2, w, ypred;
	double mlog2pi = -1.8378770664093454836;	       /* -log(2*pi) */

	LINK_INIT;
	y = ds->data_observations.y[idx];
	w = ds->data_observations.weight_wrapped_cauchy[idx];
	rho = map_probability_forward(ds->data_observations.log_prec_wrapped_cauchy[thread_id][0], MAP_FORWARD, NULL) * w;
	rho2 = SQR(rho);

	for (i = 0; i < m; i++) {
		ypred = PREDICTOR_INVERSE_LINK(x[i] + off);

		/*
		 * we need |y-ypred| <= Pi, but this might not be the case...  so we add a penalty if this condition is not met
		 */
		if (ABS(y - ypred) <= M_PI) {
			logll[i] = mlog2pi + LOG_1mp(rho2) - log1p(rho2 - 2.0 * rho * cos(y - ypred));
		} else {
			double penalty = 1.0e6 * (2.0 * rho / SQR(1.0 - rho));	/* The -Hessian in the mode... */
			if (y - ypred > M_PI) {
				logll[i] = mlog2pi + LOG_1mp(rho2) - log1p(rho2 - 2.0 * rho * cos(M_PI))
				    - penalty / 2.0 * SQR(y - ypred - M_PI);
			} else {
				logll[i] = mlog2pi + LOG_1mp(rho2) - log1p(rho2 - 2.0 * rho * cos(-M_PI))
				    - penalty / 2.0 * SQR(y - ypred + M_PI);
			}
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_stochvol(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			   void *arg, char **UNUSED(arg_str))
{
	/*
	 * y ~ N(0, var = exp(x) + 1/tau) 
	 */
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}
	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx], var;
	double tau = map_precision_forward(ds->data_observations.log_offset_prec[thread_id][0], MAP_FORWARD, NULL);
	double var_offset;

	LINK_INIT;
	var_offset = ((ISINF(tau) || ISNAN(tau)) ? 0.0 : 1.0 / tau);
	if (m > 0) {
		for (int i = 0; i < m; i++) {
			var = PREDICTOR_INVERSE_LINK(x[i] + off) + var_offset;
			logll[i] = LOG_NORMC_GAUSSIAN - 0.5 * log(var) - 0.5 * SQR(y) / var;
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
		for (int i = 0; i < -m; i++) {
			var = PREDICTOR_INVERSE_LINK(x[i] + off) + var_offset;
			logll[i] = 1.0 - 2.0 * (1.0 - inla_cdf_normal_fast(ABS(y) / sqrt(var)));
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_stochvolln(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			     void *arg, char **UNUSED(arg_str))
{
	/*
	 * y ~ N(c - 1/2 * var, var)
	 */
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx];
	double c = ds->data_observations.stochvolln_c[thread_id][0];	// identity mapping

	LINK_INIT;
	if (m > 0) {
		for (int i = 0; i < m; i++) {
			double var = PREDICTOR_INVERSE_LINK(x[i] + off);
			double mean = c - 0.5 * var;
			logll[i] = LOG_NORMC_GAUSSIAN - 0.5 * log(var) - 0.5 * SQR((y - mean)) / var;
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
		for (int i = 0; i < -m; i++) {
			double var = PREDICTOR_INVERSE_LINK(x[i] + off);
			double mean = c - 0.5 * var;
			logll[i] = 1.0 - 2.0 * (1.0 - inla_cdf_normal_fast(ABS((y - mean) / sqrt(var))));
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_stochvol_t(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			     void *arg, char **UNUSED(arg_str))
{
	/*
	 * y / exp(x/2)  ~ Student-t_dof(0, ***var = 1***)
	 *
	 * Note that Student-t_dof has variance dof/(dof-2), so we need to scale it.
	 */
	int i;
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}
	Data_section_tp *ds = (Data_section_tp *) arg;
	double dof, y, sd, sd2, obs, var_u;

	dof = map_dof_forward(ds->data_observations.dof_intern_svt[thread_id][0], MAP_FORWARD, NULL);
	y = ds->data_observations.y[idx];
	sd2 = dof / (dof - 2.0);
	sd = sqrt(sd2);

	LINK_INIT;
	if (m > 0) {
		double lg1, lg2, f;
		lg1 = gsl_sf_lngamma(dof / 2.0);
		lg2 = gsl_sf_lngamma((dof + 1.0) / 2.0);
		for (i = 0; i < m; i++) {
			var_u = PREDICTOR_INVERSE_LINK(x[i] + off);
			f = sqrt(var_u) / sd;
			obs = y / f;
			logll[i] = lg2 - lg1 - 0.5 * log(M_PI * dof) - (dof + 1.0) / 2.0 * log1p(SQR(obs) / dof) - log(f);
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
		for (i = 0; i < -m; i++) {
			var_u = PREDICTOR_INVERSE_LINK(x[i] + off);
			logll[i] = 1.0 - 2.0 * gsl_cdf_tdist_Q(ABS(y) * sd / sqrt(var_u), dof);
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_stochvol_nig(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec),
			       double *UNUSED(y_cdf), void *arg, char **UNUSED(arg_str))
{
	/*
	 * y / exp(x/2)  ~ NIG with skew and shape parameter. beta = skew, psi = shape. Note that E=1 and Var=1.
	 *
	 *
	 * density: gamma
	 *          * exp[ psi^2 + beta*(gamma*x + beta) ]
	 *          * K_1[sqrt(beta^2+psi^2)*sqrt((gamma*x+beta)^2 + psi^2)]
	 *
	 */
	if (m == 0) {
		return GMRFLib_SUCCESS;
	}
	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double skew, skew2, shape, shape2, y, gam, gam2, tmp, obs, a, var_u;

	skew = ds->data_observations.skew_intern_svnig[thread_id][0];
	skew2 = SQR(skew);
	shape = map_shape_svnig(ds->data_observations.shape_intern_svnig[thread_id][0], MAP_FORWARD, NULL);
	shape2 = SQR(shape);
	gam2 = 1.0 + SQR(skew) / SQR(shape);
	gam = sqrt(gam2);
	y = ds->data_observations.y[idx];
	a = log(gam * shape / M_PI) + 0.5 * log(skew2 + shape2) + shape2;

	LINK_INIT;
	if (m > 0) {
		for (i = 0; i < m; i++) {
			var_u = PREDICTOR_INVERSE_LINK(x[i] + off);
			obs = y / sqrt(var_u);
			tmp = SQR(gam * obs + skew) + shape2;
			logll[i] = a - 0.5 * log(tmp) + skew * (gam * obs + skew)
			    + gsl_sf_bessel_lnKnu(1.0, sqrt((skew2 + shape2) * tmp)) - log(var_u) / 2.0;
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_iid_gamma(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec),
			    double *UNUSED(y_cdf), void *arg, char **UNUSED(arg_str))
{
	/*
	 * y ~ iid_gamma
	 */
	if (m == 0) {
		return GMRFLib_SUCCESS;
	}

	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double shape, rate, w, xx, penalty = 1.0 / INLA_REAL_SMALL, cons;

	LINK_INIT;
	w = ds->data_observations.iid_gamma_scale[idx];
	shape = map_exp_forward(ds->data_observations.iid_gamma_log_shape[thread_id][0], MAP_FORWARD, NULL);
	rate = map_exp_forward(ds->data_observations.iid_gamma_log_rate[thread_id][0], MAP_FORWARD, NULL) * w;
	cons = -shape * log(rate) - gsl_sf_lngamma(shape);

	if (m > 0) {
		for (i = 0; i < m; i++) {
			xx = PREDICTOR_INVERSE_LINK(x[i] + off);
			if (xx > INLA_REAL_SMALL) {
				logll[i] = cons + (shape - 1.0) * log(xx) - rate * xx;
			} else {
				/*
				 * this is the penalty, and should not happen in the end... 
				 */
				logll[i] =
				    cons + (shape - 1.0) * log(INLA_REAL_SMALL) - rate * INLA_REAL_SMALL - penalty * SQR(INLA_REAL_SMALL - xx);
			}
		}
	}
	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_iid_logitbeta(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec),
				double *UNUSED(y_cdf), void *arg, char **UNUSED(arg_str))
{
	/*
	 * y ~ iid_logitbeta
	 */
	if (m == 0) {
		return GMRFLib_SUCCESS;
	}

	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double a, b, cons;

	LINK_INIT;
	a = map_exp_forward(ds->data_observations.iid_logitbeta_log_a[thread_id][0], MAP_FORWARD, NULL);
	b = map_exp_forward(ds->data_observations.iid_logitbeta_log_b[thread_id][0], MAP_FORWARD, NULL);
	cons = gsl_sf_lngamma(a + b) - (gsl_sf_lngamma(a) + gsl_sf_lngamma(b));

	if (m > 0) {
		for (i = 0; i < m; i++) {
			double xx = PREDICTOR_INVERSE_LINK(x[i] + off);
			logll[i] = cons + (a - 1.0) * log(xx) + (b - 1.0) * LOG_1mp(xx) + PREDICTOR_INVERSE_LINK_LOGJACOBIAN(x[i] + off);
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_loggamma_frailty(int thread_id, double *__restrict logll, double *__restrict x, int m, int UNUSED(idx), double *UNUSED(x_vec),
				   double *UNUSED(y_cdf), void *arg, char **UNUSED(arg_str))
{
	/*
	 * Log-gamma frailty Gamma(a,a), a = exp(log_prec...)
	 */
	if (m == 0) {
		return GMRFLib_SUCCESS;
	}
	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double lprec, prec;
	double log_gamma;

	// LINK_INIT;
	lprec = ds->data_observations.log_prec_loggamma_frailty[thread_id][0];
	prec = map_precision_forward(lprec, MAP_FORWARD, NULL);
	log_gamma = gsl_sf_lngamma(prec);

	if (m > 0) {
		for (i = 0; i < m; i++) {
			logll[i] = -log_gamma + prec * (lprec + x[i] - exp(x[i]));
		}
	}
	// LINK_END;

	return GMRFLib_SUCCESS;
}

int loglikelihood_logistic(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			   void *arg, char **UNUSED(arg_str))
{
	/*
	 * y ~ Logisistc. scaled so that prec = 1 gives variance = 1
	 *
	 * A := Pi/sqrt(3)
	 *
	 * > F(x);             
	 *                                                     1
	 *                                          ------------------------
	 *                                          1 + exp(-tau A (x - mu))
	 *
	 *
	 * > solve(F(x) = p,x);
	 *                                                         -1 + p
	 *                                           tau A mu - ln(- ------)
	 *                                                             p
	 *                                           -----------------------
	 *                                                    tau A
	 * > diff(F(x),x);
	 *                                          tau A exp(-tau A (x - mu))
	 *                                         ---------------------------
	 *                                                                   2
	 *                                         (1 + exp(-tau A (x - mu)))
	 *
	 */
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}
	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double y, prec, w, A = M_PI / sqrt(3.0), precA, lprecA, eta;

	LINK_INIT;
	y = ds->data_observations.y[idx];
	w = ds->data_observations.weight_logistic[idx];
	prec = map_precision_forward(ds->data_observations.log_prec_logistic[thread_id][0], MAP_FORWARD, NULL) * w;
	precA = prec * A;
	lprecA = log(precA);

	if (m > 0) {
		for (i = 0; i < m; i++) {
			eta = PREDICTOR_INVERSE_LINK(x[i] + off);
			logll[i] = lprecA - precA * (y - eta) - 2.0 * log1p(exp(-precA * (y - eta)));
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
		for (i = 0; i < -m; i++) {
			eta = PREDICTOR_INVERSE_LINK(x[i] + off);
			logll[i] = 1.0 / (1.0 + exp(-precA * (y - eta)));
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_sn(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf, void *arg,
		     char **UNUSED(arg_str))
{
	/*
	 * y ~ Skew_Normal(x, stdev)
	 */
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}
	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double y, lprec, sprec, w, xarg, ypred, *param[2], nan = NAN;
	inla_sn_arg_tp sn_arg = { 0.0, 0.0, 0.0, 0.0 };


	LINK_INIT;
	y = ds->data_observations.y[idx];
	w = ds->data_observations.sn_scale[idx];
	lprec = ds->data_observations.sn_lprec[thread_id][0] + log(w);
	sprec = exp(lprec / 2.0);

	param[0] = ds->data_observations.sn_skew[thread_id];
	param[1] = &nan;
	inla_get_sn_param(&sn_arg, param);
	assert(sn_arg.intercept == 0);

	if (m > 0) {
		for (i = 0; i < m; i++) {
			ypred = PREDICTOR_INVERSE_LINK(x[i] + off);
			xarg = ((y - ypred - sn_arg.intercept) * sprec - sn_arg.xi) / sn_arg.omega;
			logll[i] =
			    LOG_NORMC_GAUSSIAN + M_LN2 + 0.5 * lprec - log(sn_arg.omega) - 0.5 * SQR(xarg) +
			    inla_logcdf_normal_fast(sn_arg.alpha * xarg);
		}
	} else {
		double yy = (y_cdf ? *y_cdf : y);
		param[1] = &nan;			       /* this will remove the internal intercept */
		for (i = 0; i < -m; i++) {
			ypred = PREDICTOR_INVERSE_LINK(x[i] + off);
			xarg = (yy - ypred - sn_arg.intercept) * sprec;
			logll[i] = map_invsn_core(xarg, MAP_FORWARD, param, NULL);
		}
	}
	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_stochvol_sn(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec),
			      double *UNUSED(y_cdf), void *arg, char **UNUSED(arg_str))
{
	/*
	 * y ~ Skew_Normal(0, var= exp(x)+offset, skew)
	 */
	if (m == 0) {
		return GMRFLib_SUCCESS;
	}
	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double y, sprec, xarg, *param[2], nan = NAN, var_offset, var, lomega;
	inla_sn_arg_tp sn_arg = {.xi = 0.0,.omega = 0.0,.intercept = 0.0,.alpha = 0.0 };

	LINK_INIT;
	y = ds->data_observations.y[idx];
	var_offset = 1.0 / map_precision_forward(ds->data_observations.log_offset_prec[thread_id][0], MAP_FORWARD, NULL);
	param[0] = ds->data_observations.sn_skew[thread_id];
	param[1] = &nan;
	inla_get_sn_param(&sn_arg, param);
	assert(sn_arg.intercept == 0.0);
	lomega = log(sn_arg.omega);

	if (m > 0) {
		for (i = 0; i < m; i++) {
			var = PREDICTOR_INVERSE_LINK(x[i] + off) + var_offset;
			sprec = sqrt(1.0 / var);
			xarg = (y * sprec - sn_arg.xi) / sn_arg.omega;
			logll[i] =
			    LOG_NORMC_GAUSSIAN + M_LN2 + log(sprec) - lomega - 0.5 * SQR(xarg) + inla_logcdf_normal_fast(sn_arg.alpha * xarg);
		}
	}
	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_gev(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
		      void *arg, char **UNUSED(arg_str))
{
	/*
	 * y ~ GEV
	 */
	int i;
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}
	Data_section_tp *ds = (Data_section_tp *) arg;
	double y, sprec, w, xi, xx, ypred;

	LINK_INIT;
	y = ds->data_observations.y[idx];
	w = ds->data_observations.weight_gev[idx];
	sprec = sqrt(map_precision_forward(ds->data_observations.log_prec_gev[thread_id][0], MAP_FORWARD, NULL) * w);
	/*
	 * map_identity_scale(theta, MAP_FORWARD, arg...);
	 */
	xi = ds->data_observations.gev_scale_xi * ds->data_observations.xi_gev[thread_id][0];

	if (m > 0) {
		if (ISZERO(xi)) {
			for (i = 0; i < m; i++) {
				ypred = PREDICTOR_INVERSE_LINK(x[i] + off);
				xx = sprec * (y - ypred);
				logll[i] = -xx - exp(-xx) + log(sprec);
			}
		} else {
			for (i = 0; i < m; i++) {
				ypred = PREDICTOR_INVERSE_LINK(x[i] + off);
				xx = 1.0 + xi * sprec * (y - ypred);
				if (xx > INLA_REAL_SMALL) {
					logll[i] = (-1.0 / xi - 1.0) * log(xx) - pow(xx, -1.0 / xi) + log(sprec);
				} else {
					logll[i] =
					    (-1.0 / xi - 1.0) * log(INLA_REAL_SMALL) - pow(INLA_REAL_SMALL,
											   -1.0 / xi) + log(sprec) -
					    1e6 * SQR(sprec * (INLA_REAL_SMALL - xx));
				}
			}
		}
	} else {
		double yy = (y_cdf ? *y_cdf : y);
		if (ISZERO(xi)) {
			for (i = 0; i < -m; i++) {
				ypred = PREDICTOR_INVERSE_LINK(x[i] + off);
				xx = sprec * (yy - ypred);
				logll[i] = exp(-exp(-xx));
			}
		} else {
			for (i = 0; i < -m; i++) {
				ypred = PREDICTOR_INVERSE_LINK(x[i] + off);
				xx = sprec * (yy - ypred);
				double a = 1.0 + xi * xx;
				if (a > 0.0) {
					logll[i] = exp(-pow(a, -1.0 / xi));
				} else {
					logll[i] = (xi > 0.0 ? 0.0 : 1.0);
				}
			}
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_bgev(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
		       void *arg, char **UNUSED(arg_str))
{
#define f3_BETA_STD(_x) (30.0 * SQR(_x) * SQR(1.0-(_x)))
#define F3_BETA_STD(_x) (POW3(_x) * (10.0 + (-15.0 + 6.0 * (_x)) * (_x)))
#define f3_BETA(_x, _a, _b) (f3_BETA_STD((((_x) - (_a)) / ((_b) - (_a)))) / ((_b) - (_a)))
#define F3_BETA(_x, _a, _b) F3_BETA_STD((((_x) - (_a)) / ((_b) - (_a))))
#define f4_BETA_STD(_x) (140.0 * POW3(_x) * POW3(1.0 - (_x)))
#define F4_BETA_STD(_x) (gsl_pow_4(_x) * (35.0 + (-84.0 + (70.0 - 20.0 * (_x)) * (_x)) * (_x)))
#define f4_BETA(_x, _a, _b) (f4_BETA_STD(((_x) - (_a)) / ((_b) - (_a))) / ((_b) - (_a)))
#define F4_BETA(_x, _a, _b) F4_BETA_STD((((_x) - (_a)) / ((_b) - (_a))))
#define f5_BETA_STD(_x) (630.0 * gsl_pow_4(_x) * gsl_pow_4(1.0 - (_x)))
#define F5_BETA_STD(_x) (gsl_pow_5(_x) * ((_x) * ((_x) * ((_x) * ((_x) * 70.0 - 315.0) + 540.0) - 420.0) + 126.0))
#define f5_BETA(_x, _a, _b) (f5_BETA_STD((((_x) - (_a)) / ((_b) - (_a)))) / ((_b) - (_a)))
#define F5_BETA(_x, _a, _b) F5_BETA_STD((((_x) - (_a)) / ((_b) - (_a))))
#define f6_BETA_STD(_x) (2772.0 * gsl_pow_5(_x) * gsl_pow_5(1.0 - (_x)))
#define F6_BETA_STD(_x) (gsl_pow_6(_x) * ((_x) * ((_x) * ((_x) * ((_x) * ((_x) * (-252.0) + 1386.0) - 3080.0) + 3465.0) - 1980.0) + 462.0))
#define f6_BETA(_x, _a, _b) (f6_BETA_STD((((_x) - (_a)) / ((_b) - (_a)))) / ((_b) - (_a)))
#define F6_BETA(_x, _a, _b) F6_BETA_STD((((_x) - (_a)) / ((_b) - (_a))))
#define f_BETA(_x, _a, _b) f5_BETA(_x, _a, _b)
#define F_BETA(_x, _a, _b) F5_BETA(_x, _a, _b)
#define p(_x, _a, _b) F_BETA(_x, _a, _b)
#define p_deriv(_x, _a, _b) f_BETA(_x, _a, _b)
#define T1_SCALED_M1(_x, _xi) pow( 1.0 + (_xi) * (_x), -1.0 / (_xi) - 1.0)
#define T1_SCALED(_x, _xi) pow( 1.0 + (_xi) * (_x), -1.0 / (_xi))
#define T1(_x, _loc, _scale, _xi) T1_SCALED((((_x) - (_loc)) / (_scale)), _xi)
#define T2_SCALED(_x) exp(-(_x))
#define T2(_x, _loc, _scale) T2_SCALED((((_x) - (_loc)) / (_scale)))
#define log_G(_x, _loc, _scale, _xi) (-T1(_x, _loc, _scale, _xi))
#define G(_x, _loc, _scale, _xi) exp(log_G(_x, _loc, _scale, _xi))
#define log_H(_x, _loc, _scale)	(-T2(_x, _loc, _scale))
#define H(_x, _loc, _scale) exp(log_H(_x, _loc, _scale))
#define h(_x, _loc, _scale) (H(_x, _loc, _scale) * T2(_x, _loc, _scale) / (_scale))
#define log_h(_x, _loc, _scale) (log_H(_x, _loc, _scale) + log(T2(_x, _loc, _scale)) - log(_scale))
#define g(_x, _loc, _scale, _xi) (G(_x, _loc, _scale, _xi) * T1_SCALED_M1(((_x) - (_loc))/(_scale), _xi) / (_scale))
#define log_g(_x, _loc, _scale, _xi) (log_G(_x, _loc, _scale, _xi) + log(T1_SCALED_M1(((_x) - (_loc))/(_scale), _xi)) - log(_scale))

	if (0) {
		double loc = 1.2, scale = 1.5, xi = 0.5, _x;
		FILE *fp = fopen("out.txt", "w");
		for (_x = 0.001 + loc - scale / xi; _x <= loc + 4.0 * scale / xi; _x += 0.01)
			fprintf(fp, "%g %g %g %g %g %g %g\n", _x, H(_x, loc, scale), h(_x, loc, scale), exp(log_h(_x, loc, scale)),
				G(_x, loc, scale, xi), g(_x, loc, scale, xi), exp(log_g(_x, loc, scale, xi)));
		fclose(fp);
		exit(0);
	}

	if (0) {
		double a = 1.2, b = 2.5, _x;
		FILE *fp = fopen("out.txt", "w");
		for (_x = a; _x <= b; _x += (b - a) / 1000)
			fprintf(fp, "%f %g %g %g %g %g %g %g %g\n", _x,
				f3_BETA(_x, a, b), F3_BETA(_x, a, b), f4_BETA(_x, a, b),
				F4_BETA(_x, a, b), f5_BETA(_x, a, b), F5_BETA(_x, a, b), f6_BETA(_x, a, b), F6_BETA(_x, a, b));
		fclose(fp);
		exit(1);
	}

	/*
	 * y ~ GEV
	 */
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}
	Data_section_tp *ds = (Data_section_tp *) arg;
	double location, spread, log_spread, log_xi, xi, sigma, sigmaH, d, mu, muH, sprec, ypred, xx, y, w, ld;
	double qlocation = ds->data_observations.bgev_qlocation;
	double qspread = ds->data_observations.bgev_qspread;
	double mix_a, qmix_a = ds->data_observations.bgev_qmix[0];
	double mix_b, qmix_b = ds->data_observations.bgev_qmix[1];

	double ab = ds->data_observations.bgev_beta_ab;
	static double count[3] = { 0.0, 0.0, 0.0 };

	if (ab != 5.0) {
		static char first = 1;
		if (first) {
			fprintf(stderr, "*** Warning ***: loglikelihood_bgev: argument beta_ab=[%.2f] is not used\n", ab);
			fprintf(stdout, "*** Warning ***: loglikelihood_bgev: argument beta_ab=[%.2f] is not used\n", ab);
		}
		first = 0;
	}

	LINK_INIT;
	y = ds->data_observations.y[idx];
	w = ds->data_observations.bgev_scale[idx];

	int ioff = 0;
	log_spread = ds->data_observations.bgev_log_spread[thread_id][0];
	for (int i = 0; i < ds->data_observations.bgev_nbetas[0]; i++) {
		log_spread += ds->data_observations.bgev_betas[i + ioff][thread_id][0] * ds->data_observations.bgev_x[i + ioff][idx];
	}
	spread = map_exp_forward(log_spread, MAP_FORWARD, NULL) / sqrt(w);

	ioff = ds->data_observations.bgev_nbetas[0];
	log_xi = ds->data_observations.bgev_intern_tail[thread_id][0];
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wbool-compare"
	if (ISINF(log_xi) < 0) {
#pragma GCC diagnostic pop
		xi = 0.0;
		assert(ds->data_observations.bgev_nbetas[1] == 0);
	} else {
		for (int i = 0; i < ds->data_observations.bgev_nbetas[1]; i++) {
			log_xi += ds->data_observations.bgev_betas[i + ioff][thread_id][0] * ds->data_observations.bgev_x[i + ioff][idx];
		}
		xi = map_interval(log_xi, MAP_FORWARD, (void *) (ds->data_observations.bgev_tail_interval));
	}

	if (m > 0) {
		if (ISZERO(xi)) {
			d = log(-log(qspread / 2.0)) - log(-LOG_1mp(qspread / 2.0));
			for (int i = 0; i < m; i++) {
				location = PREDICTOR_INVERSE_LINK(x[i] + off);
				sigma = spread / d;
				mu = location + sigma * log(-log(qlocation));
				sprec = 1.0 / sigma;
				ypred = mu;
				xx = sprec * (y - ypred);
				logll[i] = -xx - exp(-xx) + log(sprec);
			}
		} else {
			int left = 0;
			d = (pow(-LOG_1mp(qspread / 2.0), -xi) - pow(-log(qspread / 2.0), -xi)) / xi;
			for (int i = 0; i < m; i++) {
				sigma = spread / d;
				location = PREDICTOR_INVERSE_LINK(x[i] + off);
				mu = location - sigma * ((pow(-log(qlocation), -xi) - 1.0) / xi);
				mix_a = (pow(-log(qmix_a), -xi) - 1.0) / xi * sigma + mu;
				mix_b = (pow(-log(qmix_b), -xi) - 1.0) / xi * sigma + mu;
				sigmaH = (mix_b - mix_a) / log(log(qmix_a) / log(qmix_b));
				muH = mix_a + sigmaH * log(-log(qmix_a));

				if (y >= mix_b) {
					count[0]++;
					ld = log_g(y, mu, sigma, xi);
				} else if (y <= mix_a) {
					left = 1;
					count[1]++;
					ld = log_h(y, muH, sigmaH);
				} else {
					double value_p, value_p_deriv, value_g, value_log_G, value_h, value_log_H;
					count[2]++;

					value_p = p(y, mix_a, mix_b);
					value_p_deriv = p_deriv(y, mix_a, mix_b);
					value_g = g(y, mu, sigma, xi);
					value_log_G = log_G(y, mu, sigma, xi);
					value_h = h(y, muH, sigmaH);
					value_log_H = log_H(y, muH, sigmaH);

					ld = (value_p * value_log_G + (1.0 - value_p) * value_log_H) +
					    log(value_p_deriv * value_log_G + value_p * value_g / exp(value_log_G)
						- value_p_deriv * value_log_H + (1.0 - value_p) * value_h / exp(value_log_H));
				}
				if (ISNAN(ld) || ISINF(ld)) {
					printf("bgev: idx x ld y %d %g %g %g\n", idx, x[i], ld, y);
				}
				logll[i] = ld;
			}

			if (0 && left)
				printf("right %.3g left %.3g mix %.3g\n",
				       count[0] / (count[0] + count[1] + count[2]),
				       count[1] / (count[0] + count[1] + count[2]), count[2] / (count[0] + count[1] + count[2]));
		}
	} else {
		double yy = (y_cdf ? *y_cdf : y);
		if (ISZERO(xi)) {
			d = log(-log(qspread / 2.0)) - log(-LOG_1mp(qspread / 2.0));
			for (int i = 0; i < -m; i++) {
				location = PREDICTOR_INVERSE_LINK(x[i] + off);
				sigma = spread / d;
				mu = location + sigma * log(-log(qlocation));
				sprec = 1.0 / sigma;
				ypred = mu;
				xx = sprec * (yy - ypred);
				logll[i] = exp(-exp(-xx));
			}
		} else {
			d = (pow(-LOG_1mp(qspread / 2.0), -xi) - pow(-log(qspread / 2.0), -xi)) / xi;
			for (int i = 0; i < -m; i++) {
				sigma = spread / d;
				location = PREDICTOR_INVERSE_LINK(x[i] + off);
				mu = location - sigma * ((pow(-log(qlocation), -xi) - 1.0) / xi);
				mix_a = (pow(-log(qmix_a), -xi) - 1.0) / xi * sigma + mu;
				mix_b = (pow(-log(qmix_b), -xi) - 1.0) / xi * sigma + mu;
				sigmaH = (mix_b - mix_a) / log(log(qmix_a) / log(qmix_b));
				muH = mix_a + sigmaH * log(-log(qmix_a));

				if (yy >= mix_b) {
					ld = log_G(yy, mu, sigma, xi);
				} else if (yy <= mix_a) {
					ld = log_H(yy, muH, sigmaH);
				} else {
					double value_p = p(yy, mix_a, mix_b);
					ld = value_p * log_G(yy, mu, sigma, xi) + (1.0 - value_p) * log_H(yy, muH, sigmaH);
				}
				logll[i] = exp(ld);
			}
		}
	}

#undef f3_BETA_STD
#undef F3_BETA_STD
#undef f3_BETA
#undef F3_BETA
#undef f4_BETA_STD
#undef F4_BETA_STD
#undef f4_BETA
#undef F4_BETA
#undef f5_BETA_STD
#undef F5_BETA_STD
#undef f5_BETA
#undef F5_BETA
#undef f6_BETA_STD
#undef F6_BETA_STD
#undef f6_BETA
#undef F6_BETA
#undef f_BETA
#undef F_BETA
#undef p
#undef p_deriv
#undef T1_SCALED
#undef T1
#undef T2_SCALED
#undef T2
#undef logG
#undef G
#undef logH
#undef H
#undef h
#undef log_h
#undef g
#undef log_g

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_t(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf, void *arg,
		    char **UNUSED(arg_str))
{
	/*
	 * y -x ~ (Student_t with variance 1) times 1/sqrt(precision * weight)
	 */
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	Data_section_tp *ds = (Data_section_tp *) arg;
	double y, prec, w, dof, y_std, fac, ypred;

	LINK_INIT;
	dof = map_dof_forward(ds->data_observations.dof_intern_t[thread_id][0], MAP_FORWARD, NULL);
	y = ds->data_observations.y[idx];
	w = ds->data_observations.weight_t[idx];
	prec = exp(ds->data_observations.log_prec_t[thread_id][0]) * w;
	fac = sqrt((dof / (dof - 2.0)) * prec);

	double lg1 = gsl_sf_lngamma(dof / 2.0);
	double lg2 = gsl_sf_lngamma((dof + 1.0) / 2.0);
	if (m > 0) {
		double c1 = lg2 - lg1 - 0.5 * log(M_PI * dof) + log(fac);
		double c2 = (dof + 1.0) / 2.0;
		for (int i = 0; i < m; i++) {
			ypred = PREDICTOR_INVERSE_LINK(x[i] + off);
			y_std = (y - ypred) * fac;
			logll[i] = c1 - c2 * log1p(SQR(y_std) / dof);
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
		for (int i = 0; i < -m; i++) {
			ypred = PREDICTOR_INVERSE_LINK(x[i] + off);
			logll[i] = gsl_cdf_tdist_P((y - ypred) * fac, dof);
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_tstrata(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			  void *arg, char **UNUSED(arg_str))
{
	/*
	 * y -x ~ (Student_t with variance 1) times 1/sqrt(precision * weight)
	 */
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	int strata;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double y, prec, w, dof, y_std, fac, ypred;

	LINK_INIT;
	dof = map_dof_forward(ds->data_observations.dof_intern_tstrata[thread_id][0], MAP_FORWARD, NULL);
	y = ds->data_observations.y[idx];
	w = ds->data_observations.weight_tstrata[idx];
	strata = (int) (ds->data_observations.strata_tstrata[idx] + INLA_REAL_SMALL);
	prec = exp(ds->data_observations.log_prec_tstrata[strata][thread_id][0]) * w;
	fac = sqrt((dof / (dof - 2.0)) * prec);

	double lg1 = gsl_sf_lngamma(dof / 2.0);
	double lg2 = gsl_sf_lngamma((dof + 1.0) / 2.0);
	if (m > 0) {
		double c1 = lg2 - lg1 - 0.5 * log(M_PI * dof) + log(fac);
		double c2 = (dof + 1.0) / 2.0;
		for (int i = 0; i < m; i++) {
			ypred = PREDICTOR_INVERSE_LINK(x[i] + off);
			y_std = (y - ypred) * fac;
			logll[i] = c1 - c2 * log1p(SQR(y_std) / dof);
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
		for (int i = 0; i < -m; i++) {
			ypred = PREDICTOR_INVERSE_LINK(x[i] + off);
			logll[i] = gsl_cdf_tdist_P((y - ypred) * fac, dof);
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_gpoisson(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			   void *arg, char **UNUSED(arg_str))
{
	/*
	 * y ~ GPoisson(E*exp(x))
	 */
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	int i, yy;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx];
	double phi = map_exp_forward(ds->data_observations.gpoisson_overdispersion[thread_id][0], MAP_FORWARD, NULL);
	double p = map_identity_forward(ds->data_observations.gpoisson_p[thread_id][0], MAP_FORWARD, NULL);
	double E = ds->data_observations.E[idx];
	double a, b, lambda, mu;

	if (G_norm_const_compute[idx]) {
		G_norm_const[idx] = my_gsl_sf_lnfact(y);
		G_norm_const_compute[idx] = 0;
	}
	double log_y_fact = G_norm_const[idx];

	LINK_INIT;
	if (m > 0) {
		for (i = 0; i < m; i++) {
			lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
			mu = E * lambda;
			a = mu + phi * pow(mu, p - 1.0) * y;
			b = 1 + phi * pow(mu, p - 1.0);
			logll[i] = log(mu) + (y - 1.0) * log(a) - y * log(b) - log_y_fact - a / b;
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
		for (i = 0; i < -m; i++) {
			lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
			mu = E * lambda;
			b = 1 + phi * pow(mu, p - 1.0);
			logll[i] = 0.0;
			for (yy = 0; yy <= (int) y; yy++) {
				a = mu + phi * pow(mu, p - 1.0) * yy;
				logll[i] += exp(log(mu) + (yy - 1.0) * log(a) - yy * log(b) - my_gsl_sf_lnfact(yy) - a / b);
			}
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_poisson(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			  void *arg, char **arg_str)
{
#define _logE(E_) (E_ > 0.0 ? log(E_) : 0.0)

	/*
	 * y ~ Poisson(E*exp(x)), also accept E=0, giving the likelihood y * x.
	 */
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx], E = ds->data_observations.E[idx];
	double normc;

	if (G_norm_const_compute[idx]) {
		G_norm_const[idx] = y * _logE(E) - my_gsl_sf_lnfact((int) y);
		G_norm_const_compute[idx] = 0;
	}
	normc = G_norm_const[idx];

	LINK_INIT;
	if (arg_str) {
		char *a = NULL;
		GMRFLib_sprintf(&a,
				"list(y = %.8g, family = \"poisson\", E = %.8g, link.model = \"%1s\", linear.predictor = %.8g)",
				y, E, ds->link_model, x[0]);
		*arg_str = a;
		return GMRFLib_SUCCESS;
	}

	if (m > 0) {
		double ylEmn = normc;
		if (PREDICTOR_LINK_EQ(link_log)) {
			if ((PREDICTOR_SCALE == 1.0)) {
				const int mkl_lim = 4L;
				if (m >= mkl_lim) {
					size_t mm = GMRFLib_align_simple((size_t) m, sizeof(double));
					double xx[mm];
					double exp_x[mm];
#pragma omp simd
					for (int i = 0; i < m; i++) {
						xx[i] = x[i] + off;
					}
					GMRFLib_exp(m, xx, exp_x);
					if (y > 0.0) {
#pragma omp simd
						for (int i = 0; i < m; i++) {
							logll[i] = ylEmn + y * xx[i] - E * exp_x[i];
						}
					} else {
#pragma omp simd
						for (int i = 0; i < m; i++) {
							logll[i] = ylEmn - E * exp_x[i];
						}
					}
				} else {
					if (y > 0.0) {
#pragma omp simd
						for (int i = 0; i < m; i++) {
							double xx = x[i] + off;
							logll[i] = y * xx + ylEmn - E * exp(xx);
						}
					} else {
#pragma omp simd
						for (int i = 0; i < m; i++) {
							double xx = x[i] + off;
							logll[i] = ylEmn - E * exp(xx);
						}
					}
				}
			} else {
				if (y > 0.0) {
#pragma omp simd
					for (int i = 0; i < m; i++) {
						double log_lambda = PREDICTOR_INVERSE_IDENTITY_LINK(x[i] + off);
						logll[i] = y * log_lambda + ylEmn - E * exp(log_lambda);
					}
				} else {
#pragma omp simd
					for (int i = 0; i < m; i++) {
						double log_lambda = PREDICTOR_INVERSE_IDENTITY_LINK(x[i] + off);
						logll[i] = ylEmn - E * exp(log_lambda);
					}
				}
			}
		} else {
			// general case
#pragma omp simd
			for (int i = 0; i < m; i++) {
				double lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				logll[i] = y * log(lambda) + ylEmn - E * lambda;
			}
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
		for (int i = 0; i < -m; i++) {
			double lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
			if (ISZERO(E * lambda)) {
				if (ISZERO(y)) {
					logll[i] = 1.0;
				} else {
					assert(!ISZERO(y));
				}
			} else {
				double mean = E * lambda;
				if (mean > 10000.0) {
					logll[i] = GMRFLib_cdfnorm((y + 0.5 - mean) / sqrt(mean));
				} else {
					logll[i] = gsl_cdf_poisson_P((unsigned int) y, mean);
				}
			}
		}
	}

	LINK_END;
#undef _logE

	return GMRFLib_SUCCESS;
}

int loglikelihood_npoisson(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			   void *arg, char **UNUSED(arg_str))
{
	/*
	 * y ~ Poisson(E*exp(x)) using the Normal approximation
	 */

	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx];
	double E = ds->data_observations.E[idx];

	LINK_INIT;

	if (m > 0) {
		const double l6 = 1.791759469228054957313;     /* log(6.0) */
		for (int i = 0; i < m; i++) {
			double mean = E * PREDICTOR_INVERSE_LINK(x[i] + off);
			double prec = 1.0 / mean;
			double lprec = log(prec);
			double res = y - mean;
			double right = -0.5 * prec * SQR(res + 0.5);
			double mid = -0.5 * prec * SQR(res);
			double left = -0.5 * prec * SQR(res - 0.5);
			logll[i] = LOG_NORMC_GAUSSIAN + 0.5 * lprec - l6 + mid + log(4.0 + exp(left - mid) + exp(right - mid));
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
		for (int i = 0; i < -m; i++) {
			double mean = E * PREDICTOR_INVERSE_LINK(x[i] + off);
			double sd = sqrt(mean);
			logll[i] = GMRFLib_cdfnorm((y + 0.5 - mean) / sd);
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_nzpoisson(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			    void *arg, char **arg_str)
{
	/*
	 * y ~ nzPoisson(E*exp(x))
	 */
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx], E = ds->data_observations.E[idx];
	double normc;

	if (G_norm_const_compute[idx]) {
		G_norm_const[idx] = y * log(E) - my_gsl_sf_lnfact((int) y);
		G_norm_const_compute[idx] = 0;
	}
	normc = G_norm_const[idx];

	LINK_INIT;
	if (arg_str) {
		char *a = NULL;
		GMRFLib_sprintf(&a,
				"list(y = %.8g, family = \"nzpoisson\", E = %.8g, link.model = \"%1s\", linear.predictor = %.8g)",
				y, E, ds->link_model, x[0]);
		*arg_str = a;
		return GMRFLib_SUCCESS;
	}

	if (m > 0) {
		double ylEmn = normc;
		if (PREDICTOR_LINK_EQ(link_log) && (PREDICTOR_SCALE == 1.0)) {
#pragma omp simd
			for (int i = 0; i < m; i++) {
				double xx = x[i] + off;
				double lambda = exp(xx);
				double mu = E * lambda;
				double p0 = exp(-mu);
				logll[i] = y * xx + ylEmn - mu - LOG_1mp(p0);
			}
		} else {
			for (int i = 0; i < m; i++) {
				double lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				double mu = E * lambda;
				double p0 = exp(-mu);
				logll[i] = y * log(lambda) + ylEmn - mu - LOG_1mp(p0);
			}
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
		for (int i = 0; i < -m; i++) {
			double lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
			double mu = E * lambda;
			double p0 = exp(-mu);
			logll[i] = gsl_cdf_poisson_P((unsigned int) y, mu) / (1.0 - p0);
		}
	}

	LINK_END;

	return GMRFLib_SUCCESS;
}

int loglikelihood_rcpoisson(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec),
			    double *UNUSED(y_cdf), void *arg, char **UNUSED(arg_str))
{
	if (m == 0) {
		return GMRFLib_SUCCESS;
	}

	Data_section_tp *ds = (Data_section_tp *) arg;
	Data_tp *dtp = &(ds->data_observations);
	double y = dtp->y[idx];
	double E = dtp->rcp_E[idx];
	int event = (int) round(dtp->rcp_event[idx]);
	const int iszero = ISZERO(y);

	if (G_norm_const_compute[idx]) {
		G_norm_const[idx] = my_gsl_sf_lnfact((int) y);
		G_norm_const_compute[idx] = 0;
	}
	double normc = G_norm_const[idx];

	LINK_INIT;

	if (m > 0) {
		if (event == 1) {
			for (int i = 0; i < m; i++) {
				double lambda = E * PREDICTOR_INVERSE_LINK(x[i] + off);
				logll[i] = -normc + y * log(lambda) - lambda;
			}
		} else if (event == 0) {
			for (int i = 0; i < m; i++) {
				double lambda = E * PREDICTOR_INVERSE_LINK(x[i] + off);
				logll[i] = (iszero ? 0.0 : log(gsl_cdf_poisson_Q((unsigned int) (y - 1.0), lambda)));
			}
		} else {
			double lprob = dtp->rcp_offset[idx];
			for (int i = 0; i < dtp->rcp_nbeta; i++) {
				lprob += dtp->rcp_beta[i][thread_id][0] * dtp->rcp_x[i][idx];
			}
			double log_prob = -log1p(exp(-lprob));
			double log_1mprob = -lprob + log_prob;

			// truncate if we are in the extreme limits
			double lim = -7.0;
			if (log_prob < lim) {
				// event = 0
				if (iszero) {
					GMRFLib_fill(m, 0.0, logll);
				} else {
					for (int i = 0; i < m; i++) {
						double lambda = E * PREDICTOR_INVERSE_LINK(x[i] + off);
						logll[i] = log(gsl_cdf_poisson_Q((unsigned int) (y - 1.0), lambda));
					}
				}
			} else if (log_1mprob < lim) {
				// event = 1
				if (iszero) {
					for (int i = 0; i < m; i++) {
						double lambda = E * PREDICTOR_INVERSE_LINK(x[i] + off);
						logll[i] = -normc - lambda;
					}
				} else {
					for (int i = 0; i < m; i++) {
						double lambda = E * PREDICTOR_INVERSE_LINK(x[i] + off);
						logll[i] = -normc + y * log(lambda) - lambda;
					}
				}
			} else {
				// mixture
				if (iszero) {
					for (int i = 0; i < m; i++) {
						double lambda = E * PREDICTOR_INVERSE_LINK(x[i] + off);
						double logA = log_prob - normc - lambda;
						double logB = log_1mprob;
						logll[i] = GMRFLib_logsum(logA, logB);
					}
				} else {
					for (int i = 0; i < m; i++) {
						double lambda = E * PREDICTOR_INVERSE_LINK(x[i] + off);
						double logA = log_prob - normc + y * log(lambda) - lambda;
						double logB = log_1mprob + log(gsl_cdf_poisson_Q((unsigned int) (y - 1.0), lambda));
						logll[i] = GMRFLib_logsum(logA, logB);
					}
				}
			}
		}
	}

	LINK_END;

	return GMRFLib_SUCCESS;
}

int loglikelihood_tpoisson(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			   void *arg, char **UNUSED(arg_str))
{
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	Data_section_tp *ds = (Data_section_tp *) arg;
	Data_tp *dtp = &(ds->data_observations);
	double y = dtp->y[idx];
	double E = dtp->tp_E[idx];
	int event = (int) round(dtp->tp_event[idx]);

	if (G_norm_const_compute[idx]) {
		G_norm_const[idx] = my_gsl_sf_lnfact((int) y);
		G_norm_const_compute[idx] = 0;
	}
	double normc = G_norm_const[idx];

	LINK_INIT;
	if (m > 0) {
		double prob = 1.0;
		if (event != 1) {
			double lprob = dtp->tp_offset[idx];
			for (int i = 0; i < dtp->tp_nbeta; i++) {
				lprob += dtp->tp_beta[i][thread_id][0] * dtp->tp_x[i][idx];
			}
			prob = 1.0 / (1.0 + exp(-lprob));
		}
		for (int i = 0; i < m; i++) {
			double lambda = prob * E * PREDICTOR_INVERSE_LINK(x[i] + off);
			logll[i] = -normc + y * log(lambda) - lambda;
		}
	} else {
		double *yy = (y_cdf ? y_cdf : &y);
		double prob = 1.0;
		if (event != 1) {
			double lprob = dtp->tp_offset[idx];
			for (int i = 0; i < dtp->tp_nbeta; i++) {
				lprob += dtp->tp_beta[i][thread_id][0] * dtp->tp_x[i][idx];
			}
			prob = 1.0 / (1.0 + exp(-lprob));
		}
		for (int i = 0; i < -m; i++) {
			double lambda = prob * E * PREDICTOR_INVERSE_LINK(x[i] + off);
			logll[i] = gsl_cdf_poisson_P((unsigned int) *yy, lambda);
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_bell(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
		       void *arg, char **UNUSED(arg_str))
{
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx], E = ds->data_observations.E[idx];
	double normc;

	if (G_norm_const_compute[idx]) {
		G_norm_const[idx] = my_lbell((int) y) + 1.0;
		G_norm_const_compute[idx] = 0;
	}
	normc = G_norm_const[idx];

	LINK_INIT;
	if (m > 0) {
		size_t mm = GMRFLib_align_simple((size_t) m, sizeof(double));
		double mean[mm];
		double lambda[mm];
#pragma omp simd
		for (int i = 0; i < m; i++) {
			mean[i] = E * PREDICTOR_INVERSE_LINK(x[i] + off);
		}
		my_lambert_W0s(m, mean, lambda);
#pragma omp simd
		for (int i = 0; i < m; i++) {
			logll[i] = y * log(lambda[i]) - exp(lambda[i]) + normc;
		}
	} else {
		int yy = (int) (y_cdf ? *y_cdf : y);
#pragma omp simd
		for (int i = 0; i < -m; i++) {
			double mean = E * PREDICTOR_INVERSE_LINK(x[i] + off);
			double lambda = my_lambert_W0(mean);
			double t1 = exp(1.0 - exp(lambda));
			double cdf = 1.0;
			double p = 1.0;
			for (int iy = 1; iy <= yy; iy++) {
				p *= lambda;
				cdf += p * exp(my_lbell(iy));
			}
			cdf *= t1;
			logll[i] = cdf;
		}
	}

	LINK_END;

	return GMRFLib_SUCCESS;
}

int loglikelihood_0poisson(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			   void *arg, char **UNUSED(arg_str))
{
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx];
	double E = ds->data_observations.poisson0_E[idx];

	if (G_norm_const_compute[idx]) {
		G_norm_const[idx] = y * log(E) - my_gsl_sf_lnfact((int) y);
		G_norm_const_compute[idx] = 0;
	}
	double normc = G_norm_const[idx];

	LINK_INIT;

	double prob_intern = 0.0;
	for (int i = 0; i < ds->data_observations.poisson0_nbeta; i++) {
		if (0)
			printf("idx %d i %d beta %g x %g\n", idx, i, ds->data_observations.poisson0_beta[i][thread_id][0],
			       ds->data_observations.poisson0_x[i][idx]);
		prob_intern += ds->data_observations.poisson0_beta[i][thread_id][0] * ds->data_observations.poisson0_x[i][idx];
	}
	double prob = ds->data_observations.link_simple_invlinkfunc(thread_id, prob_intern, MAP_FORWARD, NULL, NULL);

	if (m > 0) {
		double ylEmn = normc;
		if (y > 0) {
			double l1mp = log(1.0 - prob);
#pragma omp simd
			for (int i = 0; i < m; i++) {
				double lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				logll[i] = l1mp + y * log(lambda) + ylEmn - E * lambda;
			}
		} else {
#pragma omp simd
			for (int i = 0; i < m; i++) {
				double lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				logll[i] = log(prob + (1.0 - prob) * exp(-lambda));
			}
		}
	} else {
		double *yy = (y_cdf ? y_cdf : &y);
		for (int i = 0; i < -m; i++) {
			double lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
			double mean = E * lambda;
			logll[i] = prob + (1.0 - prob) * gsl_cdf_poisson_P((unsigned int) *yy, mean);
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_0poissonS(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			    void *arg, char **UNUSED(arg_str))
{
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx];
	double E = ds->data_observations.poisson0_E[idx];

	if (G_norm_const_compute[idx]) {
		G_norm_const[idx] = y * log(E) - my_gsl_sf_lnfact((int) y);
		G_norm_const_compute[idx] = 0;
	}
	double normc = G_norm_const[idx];

	LINK_INIT;

	double log_lambda = 0.0;
	for (int i = 0; i < ds->data_observations.poisson0_nbeta; i++) {
		log_lambda += ds->data_observations.poisson0_beta[i][thread_id][0] * ds->data_observations.poisson0_x[i][idx];
	}
	double lambda = ds->data_observations.link_simple_invlinkfunc(thread_id, log_lambda, MAP_FORWARD, NULL, NULL);

	if (m > 0) {
		double lpois = y * log_lambda + normc - E * lambda;
		if (y > 0) {
#pragma omp simd
			for (int i = 0; i < m; i++) {
				double prob = PREDICTOR_INVERSE_LINK(x[i] + off);
				logll[i] = log(1.0 - prob) + lpois;
			}
		} else {
			double pois = exp(lpois);
#pragma omp simd
			for (int i = 0; i < m; i++) {
				double prob = PREDICTOR_INVERSE_LINK(x[i] + off);
				logll[i] = log(prob + (1.0 - prob) * pois);
			}
		}
	} else {
		double *yy = (y_cdf ? y_cdf : &y);
		double mean = E * lambda;
		double pois = gsl_cdf_poisson_P((unsigned int) *yy, mean);
		for (int i = 0; i < -m; i++) {
			double prob = PREDICTOR_INVERSE_LINK(x[i] + off);
			logll[i] = prob + (1.0 - prob) * pois;
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_occupancy(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec),
			    double *UNUSED(y_cdf), void *arg, char **UNUSED(arg_str))
{
	if (m == 0) {
		return GMRFLib_SUCCESS;
	}

	Data_section_tp *ds = (Data_section_tp *) arg;
	int ny = ds->data_observations.occ_ny[idx];

	if (ny == 0) {
		GMRFLib_fill(IABS(m), 0.0, logll);
		return GMRFLib_SUCCESS;
	}

	LINK_INIT;

	// code is hard-coded for this case... more work is needed to make it general
	assert(PREDICTOR_LINK_EQ(link_logit));

	int nb = ds->data_observations.occ_nbeta;
	int yzero = ds->data_observations.occ_yzero[idx];

	double *X = ds->data_observations.occ_x[idx];
	int *Y = ds->data_observations.occ_y[idx];

	double beta[nb];
	for (int i = 0; i < nb; i++) {
		beta[i] = ds->data_observations.occ_beta[i][thread_id][0];
	}

	if (m > 0) {
		const int mkl_lim = 4L;
		double logll0 = 0.0;

		if (ds->data_observations.link_simple_invlinkfunc == link_logit) {
			if (ny >= mkl_lim) {
				size_t mm = GMRFLib_align_simple((size_t) ny, sizeof(double));
				double w[mm], ww[mm];
				for (int i = 0; i < ny; i++) {
					double *xx = X + i * nb;
					double Xbeta = GMRFLib_ddot(nb, beta, xx);
					w[i] = (Y[i] ? -Xbeta : Xbeta);
				}
				GMRFLib_exp(ny, w, ww);
				GMRFLib_log1p(ny, ww, w);
				logll0 -= GMRFLib_dsum(ny, w);
			} else {
				for (int i = 0; i < ny; i++) {
					double *xx = X + i * nb;
					double Xbeta = GMRFLib_ddot(nb, beta, xx);
					logll0 += (Y[i] ? -log1p(exp(-Xbeta)) : -log1p(exp(Xbeta)));
				}
			}
		} else {
			for (int i = 0; i < ny; i++) {
				double *xx = X + i * nb;
				double Xbeta = GMRFLib_ddot(nb, beta, xx);
				double prob = ds->data_observations.link_simple_invlinkfunc(thread_id, Xbeta, MAP_FORWARD, NULL, NULL);
				logll0 += (Y[i] ? LOG_p(prob) : LOG_1mp(prob));
			}
		}

		double x_critical = -0.5 * logll0;
		double x0 = 0.900 * x_critical;
		double x1 = 0.999 * x_critical;
		int tail = (GMRFLib_max_value(x, m, NULL) + off > x0);

		if (!tail && PREDICTOR_SCALE == 1.0) {
			if (yzero) {
				double elogll0 = exp(logll0);

				if (m >= mkl_lim) {
					size_t mm = GMRFLib_align_simple((size_t) m, sizeof(double));
					double xx[mm], exx[mm];
#pragma omp simd
					for (int i = 0; i < m; i++) {
						xx[i] = x[i] + off;
					}
					GMRFLib_exp(m, xx, exx);
#pragma omp simd
					for (int i = 0; i < m; i++) {
						double a = 1.0 + 1.0 / exx[i];
						double b = (1.0 + exx[i]) * elogll0;
						xx[i] = 1.0 / a + 1.0 / b;
					}
					GMRFLib_log(m, xx, exx);
#pragma omp simd
					for (int i = 0; i < m; i++) {
						logll[i] = logll0 + exx[i];
					}
				} else {
#pragma omp simd
					for (int i = 0; i < m; i++) {
						double ex = exp(x[i] + off);
						double exd = 1.0 / ex;
						logll[i] = logll0 + log(1.0 / (1.0 + exd) + 1.0 / ((1.0 + ex) * elogll0));
					}
				}
			} else {
				if (m >= mkl_lim) {
					size_t mm = GMRFLib_align_simple((size_t) m, sizeof(double));
					double w[mm], ww[mm];
#pragma omp simd
					for (int i = 0; i < m; i++) {
						w[i] = -(x[i] + off);
					}
					GMRFLib_exp(m, w, ww);
					GMRFLib_log1p(m, ww, w);
#pragma omp simd
					for (int i = 0; i < m; i++) {
						logll[i] = logll0 - w[i];
					}
				} else {
#pragma omp simd
					for (int i = 0; i < m; i++) {
						double exd = exp(-(x[i] + off));
						logll[i] = logll0 - log1p(exd);
					}
				}
			}
		} else {
			if (yzero) {
				if (tail) {
					// we fix the tail in the low-likelihood area to avoid issues
					// this is documented in 'internal-doc/occupancy/description.pdf'
					double a = exp(logll0);
					for (int i = 0; i < m; i++) {
						double xoff = x[i] + off;
						double xx = PREDICTOR_INVERSE_IDENTITY_LINK(xoff);
						if (xx <= x0) {
							double phi = PREDICTOR_INVERSE_LINK(xoff);
							logll[i] = GMRFLib_logsum(logll0 + LOG_p(phi), LOG_1mp(phi));
						} else {
							double z = (xx - x0) / (x1 - x0);
							double xx0 = x0 + (x1 - x0) * z / (1.0 + z);
							double dx = xx - xx0;
							double exx0 = exp(-xx0);
							double t1 = a + exx0;
							double t2 = 1.0 + exx0;
							double c0 = log(t1 / t2);
							double c1 = (exx0 * (a - 1.0)) / (t1 * t2);
							double c2 = ((SQR(exx0) - a) * (a - 1.0) * exx0) / SQR(t1 * t2);
							logll[i] = c0 + dx * (c1 + 0.5 * c2 * dx);
						}
					}
				} else {
					for (int i = 0; i < m; i++) {
						double phi = PREDICTOR_INVERSE_LINK(x[i] + off);
						logll[i] = GMRFLib_logsum(logll0 + LOG_p(phi), LOG_1mp(phi));
					}
				}
			} else {
				for (int i = 0; i < m; i++) {
					double phi = PREDICTOR_INVERSE_LINK(x[i] + off);
					logll[i] = logll0 + LOG_p(phi);
				}
			}
		}
	} else {
		GMRFLib_fill(IABS(m), 0.0, logll);
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_0binomial(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			    void *arg, char **UNUSED(arg_str))
{
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx];
	double n = ds->data_observations.binomial0_Ntrials[idx];
	double ny = n - y;

	if (ISZERO(y) && ISZERO(n)) {
		double val = (m > 0 ? 0.0 : 1.0);
		for (int i = 0; i < m; i++) {
			logll[i] = val;
		}
		return GMRFLib_SUCCESS;
	}

	gsl_sf_result res = { 0, 0 };
	if (G_norm_const_compute[idx]) {
		gsl_sf_lnchoose_e((unsigned int) n, (unsigned int) y, &res);
		G_norm_const[idx] = res.val;
		G_norm_const_compute[idx] = 0;
	}
	res.val = G_norm_const[idx];
	LINK_INIT;

	double prob_intern = 0.0;
	for (int i = 0; i < ds->data_observations.binomial0_nbeta; i++) {
		prob_intern += ds->data_observations.binomial0_beta[i][thread_id][0] * ds->data_observations.binomial0_x[i][idx];
	}
	double prob = ds->data_observations.link_simple_invlinkfunc(thread_id, prob_intern, MAP_FORWARD, NULL, NULL);

	if (m > 0) {
		if (y > 0) {
			double tmp = log(1.0 - prob) + res.val;
			if (PREDICTOR_LINK_EQ(link_logit)) {
#pragma omp simd
				for (int i = 0; i < m; i++) {
					double eta = PREDICTOR_INVERSE_IDENTITY_LINK(x[i] + off);
					double ee = exp(eta);
					double log_1mp = -log1p(ee);
					double log_p = -log1p(1.0 / ee);
					logll[i] = tmp + y * log_p + ny * log_1mp;
				}
			} else {
#pragma omp simd
				for (int i = 0; i < m; i++) {
					double p = PREDICTOR_INVERSE_LINK(x[i] + off);
					logll[i] = tmp + y * LOG_p(p) + ny * LOG_1mp(p);
				}
			}
		} else {
			double ltmp = log((1.0 - prob) / prob) + res.val;
			double lprob = LOG_p(prob);
#pragma omp simd
			for (int i = 0; i < m; i++) {
				double p = PREDICTOR_INVERSE_LINK(x[i] + off);
				double tt = exp(ltmp + n * LOG_1mp(p));
				logll[i] = lprob + log1p(tt);
			}
		}
	} else {
		double *yy = (y_cdf ? y_cdf : &y);

		for (int i = 0; i < -m; i++) {
			double p = PREDICTOR_INVERSE_LINK(x[i] + off);
			logll[i] = prob + (1.0 - prob) * gsl_cdf_binomial_P((unsigned int) *yy, p, (unsigned int) n);
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_0binomialS(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			     void *arg, char **UNUSED(arg_str))
{
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx];
	double n = ds->data_observations.binomial0_Ntrials[idx];
	double ny = n - y;

	if (ISZERO(y) && ISZERO(n)) {
		double val = (m > 0 ? 0.0 : 1.0);
		for (int i = 0; i < m; i++) {
			logll[i] = val;
		}
		return GMRFLib_SUCCESS;
	}

	gsl_sf_result res = { 0, 0 };
	if (G_norm_const_compute[idx]) {
		gsl_sf_lnchoose_e((unsigned int) n, (unsigned int) y, &res);
		G_norm_const[idx] = res.val;
		G_norm_const_compute[idx] = 0;
	}
	res.val = G_norm_const[idx];
	LINK_INIT;

	double p_intern = 0.0;
	for (int i = 0; i < ds->data_observations.binomial0_nbeta; i++) {
		p_intern += ds->data_observations.binomial0_beta[i][thread_id][0] * ds->data_observations.binomial0_x[i][idx];
	}
	double p = ds->data_observations.link_simple_invlinkfunc(thread_id, p_intern, MAP_FORWARD, NULL, NULL);

	// assert(ds->data_observations.link_simple_invlinkfunc == link_logit);
	// assert(ds->predictor_invlinkfunc == link_cloglog);

	if (m > 0) {
		if (y > 0) {
			double tmp = res.val + y * LOG_p(p) + ny * LOG_1mp(p);
#pragma omp simd
			for (int i = 0; i < m; i++) {
				double prob = PREDICTOR_INVERSE_LINK(x[i] + off);
				logll[i] = tmp + LOG_1mp(prob);
			}
		} else {
			double a = exp(res.val + n * LOG_1mp(p));
#pragma omp simd
			for (int i = 0; i < m; i++) {
				double prob = PREDICTOR_INVERSE_LINK(x[i] + off);
				double tt = (1.0 - prob) / prob * a;
				logll[i] = LOG_p(prob) + log1p(tt);
			}
		}
	} else {
		double *yy = (y_cdf ? y_cdf : &y);
		double val = gsl_cdf_binomial_P((unsigned int) *yy, p, (unsigned int) n);
		for (int i = 0; i < -m; i++) {
			double prob = PREDICTOR_INVERSE_LINK(x[i] + off);
			logll[i] = prob + (1.0 - prob) * val;
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_binomialmix(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			      void *arg, char **UNUSED(arg_str))
{
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx];
	double *dat = ds->data_observations.binmix_dat[idx];
	double n = dat[BINOMIALMIX_NBETA + 4];		       /* z_1, ..., z_11, w_1, w_2 */
	double ny = n - y;

	if (ISZERO(y) && ISZERO(n)) {
		double val = (m > 0 ? 0.0 : 1.0);
		for (int i = 0; i < m; i++) {
			logll[i] = val;
		}
		return GMRFLib_SUCCESS;
	}

	if (G_norm_const_compute[idx]) {
		gsl_sf_result res = { 0, 0 };
		gsl_sf_lnchoose_e((unsigned int) n, (unsigned int) y, &res);
		G_norm_const[idx] = res.val;
		G_norm_const_compute[idx] = 0;
	}
	double normc = G_norm_const[idx];
	LINK_INIT;

	double beta[BINOMIALMIX_NBETA];
	for (int i = 0; i < BINOMIALMIX_NBETA; i++) {
		beta[i] = ds->data_observations.binmix_beta[i][thread_id][0];
	}
	double beta9 = beta[BINOMIALMIX_NBETA - 1];

	int nbeta2 = BINOMIALMIX_NBETA / 2L;
	double p1_intern = beta9 * dat[BINOMIALMIX_NBETA - 1] + GMRFLib_ddot(nbeta2, beta, dat);
	double p2_intern = beta9 * dat[BINOMIALMIX_NBETA + 0] + GMRFLib_ddot(nbeta2, beta + nbeta2, dat + nbeta2);
	double p1 = PREDICTOR_INVERSE_LINK_NO_SCALE(p1_intern);
	double p2 = PREDICTOR_INVERSE_LINK_NO_SCALE(p2_intern);
	double w1 = dat[BINOMIALMIX_NBETA + 2];
	double w2 = dat[BINOMIALMIX_NBETA + 3];
	double p12 = w1 * p1 + w2 * p2;
	double w3 = 1.0 - (w1 + w2);
	assert(w3 >= 0 && w3 <= 1.0);

	double offset = off + beta9 * dat[BINOMIALMIX_NBETA + 1];
	if (m > 0) {
		if (ISZERO(y)) {
#pragma omp simd
			for (int i = 0; i < m; i++) {
				double p = p12 + w3 * PREDICTOR_INVERSE_LINK(x[i] + offset);
				logll[i] = normc + ny * LOG_1mp(p);
			}
		} else if (ISZERO(ny)) {
#pragma omp simd
			for (int i = 0; i < m; i++) {
				double p = p12 + w3 * PREDICTOR_INVERSE_LINK(x[i] + offset);
				logll[i] = normc + y * LOG_p(p);
			}
		} else {
#pragma omp simd
			for (int i = 0; i < m; i++) {
				double p = p12 + w3 * PREDICTOR_INVERSE_LINK(x[i] + offset);
				logll[i] = normc + y * LOG_p(p) + ny * LOG_1mp(p);
			}
		}
	} else {
		double *yy = (y_cdf ? y_cdf : &y);
		for (int i = 0; i < -m; i++) {
			double p = p12 + w3 * PREDICTOR_INVERSE_LINK(x[i] + offset);
			logll[i] = gsl_cdf_binomial_P((unsigned int) *yy, p, (unsigned int) n);
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

double eval_log_contpoisson(double y, double lambda)
{
	// use that f.cont(x+1/2) approx f.poisson(x) for integer x's. I wasn't able to get the incomplete_gamma version
	// accurate enough to work. maybe the connection with the Gamma-distribution is a nice way to go?
#define _R 2
#define _L 3
#define _LEN (_R + _L + 2)
	int i, istart, iy, low, high, len;
	double work[2 * _LEN], *xx = NULL, *yy = NULL, lval;
	GMRFLib_spline_tp *spline = NULL;

	low = IMAX(0, (int) y - _L);
	high = (int) y + _R;
	len = high - low + 1;
	xx = work;
	yy = work + _LEN;

	if (low == 0) {
		xx[0] = 0.0;
		yy[0] = log(INLA_REAL_SMALL);
		istart = 1;
		len++;
	} else {
		istart = 0;
	}

	for (iy = low, i = istart; iy <= high; iy++, i++) {
		xx[i] = iy + 0.5;
		yy[i] = iy * log(lambda) - lambda - my_gsl_sf_lnfact(iy);
	}
	spline = GMRFLib_spline_create(xx, yy, len);
	lval = GMRFLib_spline_eval(y, spline);
	GMRFLib_spline_free(spline);
#undef _L
#undef _R
#undef _LEN

	return (lval);
}

int loglikelihood_contpoisson(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			      void *arg, char **UNUSED(arg_str))
{
	// this model is disabled
	assert(0 == 1);

	/*
	 * y ~ ContPoisson(E*exp(x))
	 */
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx], E = ds->data_observations.E[idx], lambda;

	LINK_INIT;
	if (m > 0) {
		for (i = 0; i < m; i++) {
			lambda = E * PREDICTOR_INVERSE_LINK(x[i] + off);
			logll[i] = eval_log_contpoisson(y + 1.0, lambda);
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);

		// slight inconsistency, as we use the 'exact' expression here, and an (good) approximation above.
		double normc = exp(gsl_sf_lngamma(y + 1.0));
		for (i = 0; i < -m; i++) {
			lambda = E * PREDICTOR_INVERSE_LINK(x[i] + off);
			logll[i] = gsl_sf_gamma_inc(y + 1.0, lambda) / normc;
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_qcontpoisson(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			       void *arg, char **UNUSED(arg_str))
{
	// this model is disabled
	assert(0 == 1);

	/*
	 * y ~ ContPoisson(E*exp(x)), also accept E=0, giving the likelihood y * x.  quantile version
	 */
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	int i, id = 0;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx], E = ds->data_observations.E[idx], lambda, q;

	GMRFLib_CACHE_SET_ID(id);
	LINK_INIT;
	if (m > 0) {
		for (i = 0; i < m; i++) {
			q = PREDICTOR_INVERSE_LINK(x[i] + off);
			lambda = E * exp(GMRFLib_spline_eval(log(q), ds->data_observations.qcontpoisson_func[id]));
			logll[i] = eval_log_contpoisson(y + 1.0, lambda);
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);

		// slight inconsistency, as we use the 'exact' expression here, and an (good) approximation above.
		double normc = exp(gsl_sf_lngamma(y + 1.0));
		for (i = 0; i < -m; i++) {
			q = PREDICTOR_INVERSE_LINK(x[i] + off);
			lambda = E * exp(GMRFLib_spline_eval(log(q), ds->data_observations.qcontpoisson_func[id]));
			logll[i] = gsl_sf_gamma_inc(y + 1.0, lambda) / normc;
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

double inla_poisson_interval(double mean, int ifrom, int ito)
{
	// Compute Prob(ifrom <= Y <= ito) for the poisson with given mean.
	// NOTE1: Both ends of the interval are included.
	// NOTE2: if ito < 0, then 'ito' is interpreted as INFINITY
	// NOTE3: if ifrom < 0, then 'ifrom' is interpreted as INFINITY

	double prob, prob_sum = 0.0;

	if (ifrom < 0) {
		// ifrom=INFINITE
	} else if (ito < 0) {
		if (ifrom == 0) {
			prob_sum = 1.0;
		} else {
			prob_sum = 1.0 - gsl_cdf_poisson_P((unsigned int) (ifrom - 1), mean);
		}
	} else {
		assert(ito >= ifrom);
		prob_sum = prob = inla_ipow(mean, ifrom) * exp(-mean - my_gsl_sf_lnfact(ifrom));
		for (int y = ifrom + 1; y <= ito; y++) {
			prob *= mean / (double) y;
			prob_sum += prob;
		}
	}
	return (prob_sum);
}

double inla_negative_binomial_interval(double size, double mu, int y_from, int y_to)
{
	// Compute Prob(y_from <= Y <= y_to) for nbinomial(size,prob).
	// NOTE1: Both ends of the interval are included.
	// NOTE2: if y_to < 0, then 'y_to' is interpreted as INFINITY
	// NOTE3: if y_from < 0, then 'y_from' is interpreted as INFINITY

	double prob = size / (size + mu);
	double p, p_sum = 0.0;

	if (y_from < 0) {
		p_sum = 0.0;
	} else if (y_to < 0) {
		if (y_from == 0) {
			p_sum = 1.0;
		} else {
			p_sum = 1.0 - gsl_cdf_negative_binomial_P((unsigned int) (y_from - 1), prob, size);
		}
	} else {
		assert(y_to >= y_from);
		p = p_sum = gsl_ran_negative_binomial_pdf((unsigned int) y_from, prob, size);
		double pp = 1.0 - prob;
		for (int y = y_from + 1; y <= y_to; y++) {
			double yy = (double) y;
			p *= (yy + size - 1.0) / yy * pp;
			p_sum += p;
		}
	}
	return (p_sum);
}

int loglikelihood_cenpoisson2(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			      void *arg, char **UNUSED(arg_str))
{
	/*
	 * y ~ Poisson(E*exp(x)), [cen_low,cen_high] is cencored. cen_high<0 means Inf
	 */
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double mu;
	double lambda;
	double y = ds->data_observations.y[idx];
	double E = ds->data_observations.E[idx];
	double cen_low = ds->data_observations.cen_low[idx];
	double cen_high = ds->data_observations.cen_high[idx];
	double normc = my_gsl_sf_lnfact((int) y);
	int int_low = (int) cen_low;
	int int_high = (int) cen_high;

	LINK_INIT;
	if (m > 0) {
		if (int_low < 0) {
			for (i = 0; i < m; i++) {
				lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				mu = E * lambda;
				logll[i] = y * log(mu) - mu - normc;
			}
		} else if (y >= int_low && (int_high < 0 || y <= int_high)) {
			for (i = 0; i < m; i++) {
				lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				mu = E * lambda;
				logll[i] = log(inla_poisson_interval(mu, int_low, int_high));
			}
		} else {
			for (i = 0; i < m; i++) {
				lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				mu = E * lambda;
				logll[i] = y * log(mu) - mu - normc;
			}
		}
	} else {
		double *yy = (y_cdf ? y_cdf : &y);
		int iy = (int) (*yy);

		for (i = 0; i < -m; i++) {
			lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
			mu = E * lambda;
			if (ISZERO(mu)) {
				if (ISZERO(iy)) {
					logll[i] = 1.0;
				} else {
					assert(!ISZERO(iy));
				}
			} else {
				if (int_low < 0) {
					logll[i] = gsl_cdf_poisson_P((unsigned int) iy, mu);
				} else if (iy < int_low || (int_high >= 0 && iy > int_high)) {
					logll[i] = gsl_cdf_poisson_P((unsigned int) iy, mu);
				} else {
					if (int_low > 0) {
						logll[i] = gsl_cdf_poisson_P((unsigned int) (int_low - 1), mu);
					} else {
						logll[i] = 0.0;
					}
					logll[i] += 0.5 * inla_poisson_interval(mu, int_low, int_high);
				}
			}
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_cenpoisson(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			     void *arg, char **UNUSED(arg_str))
{
	/*
	 * y ~ Poisson(E*exp(x)), also accept E=0, giving the likelihood y * x. values in CENINTERVAL is cencored
	 *
	 * interval[1] < 0 means infinity
	 */
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double *interval = ds->data_observations.cenpoisson_interval, mu;
	double y = ds->data_observations.y[idx], E = ds->data_observations.E[idx], normc = my_gsl_sf_lnfact((int) y), lambda;
	int int_low = (int) interval[0], int_high = (int) interval[1];

	// must use 'double' to store an INF
	if (ISINF(interval[1]) || interval[1] < 0) {
		int_high = -1;
	}

	LINK_INIT;
	if (m > 0) {
		for (i = 0; i < m; i++) {
			lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
			mu = E * lambda;
			if (y >= int_low && (int_high < 0 || y <= int_high)) {
				logll[i] = log(inla_poisson_interval(mu, int_low, int_high));
			} else {
				logll[i] = y * log(mu) - mu - normc;
			}
		}
	} else {
		double *yy = (y_cdf ? y_cdf : &y);
		int iy = (int) (*yy);

		for (i = 0; i < -m; i++) {
			lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
			mu = E * lambda;
			if (ISZERO(mu)) {
				if (ISZERO(iy)) {
					logll[i] = 1.0;
				} else {
					assert(!ISZERO(iy));
				}
			} else {
				if (iy < int_low || (int_high >= 0 && iy > int_high)) {
					logll[i] = gsl_cdf_poisson_P((unsigned int) iy, mu);
				} else {
					if (int_low > 0) {
						logll[i] = gsl_cdf_poisson_P((unsigned int) (int_low - 1), mu);
					} else {
						logll[i] = 0.0;
					}
					logll[i] += 0.5 * inla_poisson_interval(mu, int_low, int_high);
				}
			}
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_zeroinflated_cenpoisson0(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec),
					   double *y_cdf, void *arg, char **UNUSED(arg_str))
{
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double *interval = ds->data_observations.cenpoisson_interval;
	double mu, p0, fac, p = map_probability_forward(ds->data_observations.prob_intern[thread_id][0], MAP_FORWARD, NULL);
	double y = ds->data_observations.y[idx], E = ds->data_observations.E[idx], normc = my_gsl_sf_lnfact((int) y), lambda;

	LINK_INIT;
	/*
	 * '|y>0', and Prob(y=0) = exp(-mean) 
	 */

	if (m > 0) {
		if ((int) y == 0) {
			for (i = 0; i < m; i++) {
				logll[i] = LOG_p(p);
			}
		} else {
			for (i = 0; i < m; i++) {
				lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				mu = E * lambda;
				p0 = exp(-mu);
				fac = (1.0 - p) / (1.0 - p0);
				if (y >= interval[0] && y <= interval[1]) {
					logll[i] = log(fac) + log(gsl_cdf_poisson_P((unsigned int) interval[1], mu)
								  - gsl_cdf_poisson_P((unsigned int) (interval[0] - 1L), mu));
				} else {
					logll[i] = log(fac) + y * log(mu) - mu - normc;
				}
			}
		}

	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
		int iy = (int) y;

		if ((int) iy == 0) {
			for (i = 0; i < -m; i++) {
				logll[i] = p;
			}
		} else {
			for (i = 0; i < -m; i++) {
				lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				mu = E * lambda;
				p0 = exp(-mu);
				fac = (1.0 - p) / (1.0 - p0);
				if (iy == 0) {
					logll[i] = p;
				} else if (iy < interval[0] || iy > interval[1]) {
					// not censored
					logll[i] = p + fac * gsl_cdf_poisson_P((unsigned int) iy, E * lambda);
				} else {
					int ii;
					double sum = 0.0, prob, one = 0.0;
					for (ii = interval[0]; ii <= interval[1]; ii++) {
						prob = gsl_ran_poisson_pdf((unsigned int) ii, mu);
						sum += prob * (p + fac * gsl_cdf_poisson_P((unsigned int) ii, mu));
						one += prob;
					}
					logll[i] = sum / one;
				}
			}
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_zeroinflated_cenpoisson1(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec),
					   double *y_cdf, void *arg, char **UNUSED(arg_str))
{
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	Data_section_tp *ds = (Data_section_tp *) arg;
	double *interval = ds->data_observations.cenpoisson_interval;
	double mu, p = map_probability_forward(ds->data_observations.prob_intern[thread_id][0], MAP_FORWARD, NULL);
	double y = ds->data_observations.y[idx], E = ds->data_observations.E[idx], normc = my_gsl_sf_lnfact((int) y), lambda;

	LINK_INIT;
	/*
	 * '|y>0', and Prob(y=0) = exp(-mean) 
	 */

	if (m > 0) {
		if ((int) y == 0) {
			for (int i = 0; i < m; i++) {
				lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				mu = E * lambda;
				logll[i] = log(p + (1.0 - p) * gsl_ran_poisson_pdf((unsigned int) y, mu));
			}
		} else if (y >= interval[0] && y <= interval[1]) {
			for (int i = 0; i < m; i++) {
				lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				mu = E * lambda;
				logll[i] = log((1.0 - p) * (gsl_cdf_poisson_P((unsigned int) interval[1], mu)
							    - gsl_cdf_poisson_P((unsigned int) (interval[0] - 1L), mu)));
			}
		} else {
			for (int i = 0; i < m; i++) {
				lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				mu = E * lambda;
				logll[i] = LOG_1mp(p) + y * log(mu) - mu - normc;
			}
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
		int iy = (int) y;

		if (iy < interval[0] || iy > interval[1]) {
			for (int i = 0; i < -m; i++) {
				lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				mu = E * lambda;
				// not censored
				logll[i] = p + (1.0 - p) * gsl_cdf_poisson_P((unsigned int) iy, mu);
			}
		} else {
			for (int i = 0; i < -m; i++) {
				lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				mu = E * lambda;
				double sum = 0.0, prob, one = 0.0;
				for (int ii = interval[0]; ii <= interval[1]; ii++) {
					prob = gsl_ran_poisson_pdf((unsigned int) ii, mu);
					sum += prob * (p + (1.0 - p) * gsl_cdf_poisson_P((unsigned int) ii, mu));
					one += prob;
				}
				logll[i] = sum / one;
			}
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_pom(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *UNUSED(y_cdf),
		      void *arg, char **UNUSED(arg_str))
{
#define _F_CORE_LOGIT(_x) (1.0/(1.0 + exp(-(_x))))
#define _P_LOGIT(_class, _eta) ((_class) == 1 ? _F_CORE_LOGIT(alpha[(_class)] - (_eta)) : \
				((_class) == nclasses ? (1.0 - _F_CORE_LOGIT(alpha[(_class) -1] - (_eta))) : \
				 (_F_CORE_LOGIT(alpha[(_class)] - (_eta)) - _F_CORE_LOGIT(alpha[(_class) -1] - (_eta)))))
#define _F_CORE_PROBIT(_x) inla_cdf_normal_fast(_x)
#define _P_PROBIT(_class, _eta) ((_class) == 1 ? _F_CORE_PROBIT(alpha[(_class)] - (_eta)) : \
				 ((_class) == nclasses ? (1.0 - _F_CORE_PROBIT(alpha[(_class) -1] - (_eta))) : \
				  (_F_CORE_PROBIT(alpha[(_class)] - (_eta)) - _F_CORE_PROBIT(alpha[(_class) -1] - (_eta)))))

#define _F_CORE_PROBIT_FAST(_x) inla_cdf_normal_fast(_x)
#define _P_PROBIT_FAST(_class, _eta) ((_class) == 1 ? _F_CORE_PROBIT_FAST(alpha[(_class)] - (_eta)) : \
				      ((_class) == nclasses ? (1.0 - _F_CORE_PROBIT_FAST(alpha[(_class) -1] - (_eta))) : \
				       (_F_CORE_PROBIT_FAST(alpha[(_class)] - (_eta)) - _F_CORE_PROBIT_FAST(alpha[(_class) -1] - (_eta)))))

	/*
	 * y ~ POM(alpha_k + eta)
	 */
	if (m == 0) {
		return GMRFLib_SUCCESS;
	}

	double eta, theta, *alpha = NULL;
	Data_section_tp *ds = (Data_section_tp *) arg;
	int i, k, iy = (int) ds->data_observations.y[idx], nclasses = ds->data_observations.pom_nclasses;
	int use_logit = (ds->data_observations.pom_cdf == POM_CDF_LOGIT);
	int fast_probit = ds->data_observations.pom_fast_probit;

	static double **calpha = NULL;
	static int *nclass = NULL;

	if (!calpha) {
#pragma omp critical (Name_0e94df4562241e37d016a0edfddb0588df8765f1)
		{
			if (!calpha) {
				nclass = Calloc(GMRFLib_CACHE_LEN(), int);
				calpha = Calloc(GMRFLib_CACHE_LEN(), double *);
				for (i = 0; i < GMRFLib_CACHE_LEN(); i++) {
					nclass[i] = 4;
					calpha[i] = Calloc(nclass[i], double);
				}
			}
		}
	}

	int id = 0;
	GMRFLib_CACHE_SET_ID(id);

	if (nclasses > nclass[id]) {
		nclass[id] = nclasses;
		calpha[id] = Realloc(calpha[id], nclass[id], double);
	}
	alpha = calpha[id];

	for (i = 0; i < nclasses - 1; i++) {
		k = 1 + i;
		theta = map_identity_forward(ds->data_observations.pom_theta[i][thread_id][0], MAP_FORWARD, NULL);
		alpha[k] = (k == 1 ? theta : alpha[k - 1] + exp(theta));
	}

	LINK_INIT;
	if (m > 0) {
		if (use_logit) {
			for (i = 0; i < m; i++) {
				eta = PREDICTOR_INVERSE_LINK(x[i] + off);
				logll[i] = log(_P_LOGIT(iy, eta));
			}
		} else {
			if (fast_probit) {
				for (i = 0; i < m; i++) {
					eta = PREDICTOR_INVERSE_LINK(x[i] + off);
					logll[i] = log(_P_PROBIT_FAST(iy, eta));
				}
			} else {
				for (i = 0; i < m; i++) {
					eta = PREDICTOR_INVERSE_LINK(x[i] + off);
					logll[i] = log(_P_PROBIT(iy, eta));
				}
			}
		}
	} else {
		assert(0 == 1);
	}

	LINK_END;
#undef _P_LOGIT
#undef _P_PROBIT
#undef _P_PROBIT_FAST
#undef _F_CORE_LOGIT
#undef _F_CORE_PROBIT
#undef _F_CORE_PROBIT_FAST

	return GMRFLib_SUCCESS;
}

int loglikelihood_zeroinflated_poisson0(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec),
					double *y_cdf, void *arg, char **UNUSED(arg_str))
{
	/*
	 * zeroinflated Poission: y ~ p*1[y=0] + (1-p)*Poisson(E*exp(x) | y > 0)
	 */
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx], E = ds->data_observations.E[idx], normc = my_gsl_sf_lnfact((int) y),
	    p = map_probability_forward(ds->data_observations.prob_intern[thread_id][0], MAP_FORWARD, NULL), mu, lambda;

	LINK_INIT;
	if ((int) y == 0) {
		/*
		 * this is just the point-mass at zero 
		 */
		if (m > 0) {
			for (i = 0; i < m; i++) {
				logll[i] = LOG_p(p);
			}
		} else {
			for (i = 0; i < -m; i++) {
				logll[i] = p;
			}
		}
	} else {
		/*
		 * As for the Poisson but '|y>0', and Prob(y=0) = exp(-mean) 
		 */
		double p0;
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
		if (m > 0) {
			double c1 = LOG_1mp(p) - normc;
			for (i = 0; i < m; i++) {
				lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				mu = E * lambda;
				p0 = exp(-mu);
				logll[i] = c1 + y * log(mu) - mu - LOG_1mp(p0);
			}
		} else {
			for (i = 0; i < -m; i++) {
				lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				mu = E * lambda;
				p0 = exp(-mu);
				logll[i] = p + (1.0 - p) * (gsl_cdf_poisson_P((unsigned int) y, mu) - p0) / (1.0 - p0);
			}
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_zeroinflated_poisson1(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec),
					double *y_cdf, void *arg, char **UNUSED(arg_str))
{
	/*
	 * zeroinflated Poission: y ~ p*1[y=0] + (1-p)*Poisson(E*exp(x))
	 */
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx], E = ds->data_observations.E[idx], normc = my_gsl_sf_lnfact((int) y),
	    p = map_probability_forward(ds->data_observations.prob_intern[thread_id][0], MAP_FORWARD, NULL), mu, lambda, logA, logB;

	LINK_INIT;
	if ((int) y == 0) {
		if (m > 0) {
			for (i = 0; i < m; i++) {
				lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				mu = E * lambda;
				logA = LOG_p(p);
				logB = LOG_1mp(p) + y * log(mu) - mu - normc;
				// logll[i] = log(p + (1.0 - p) * gsl_ran_poisson_pdf((unsigned int) y, mu));
				logll[i] = GMRFLib_logsum(logA, logB);
			}
		} else {
			GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
			for (i = 0; i < -m; i++) {
				lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				mu = E * lambda;
				logll[i] = p + (1.0 - p) * gsl_cdf_poisson_P((unsigned int) y, mu);
			}
		}
	} else {
		if (m > 0) {
			for (i = 0; i < m; i++) {
				lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				mu = E * lambda;
				logll[i] = LOG_1mp(p) + y * log(mu) - mu - normc;
			}
		} else {
			GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
			for (i = 0; i < -m; i++) {
				lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				mu = E * lambda;
				logll[i] = p + (1.0 - p) * gsl_cdf_poisson_P((unsigned int) y, mu);
			}
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_zeroinflated_poisson2(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec),
					double *y_cdf, void *arg, char **UNUSED(arg_str))
{
	/*
	 * zeroinflated Poission: y ~ p*1[y=0] + (1-p)*Poisson(E*exp(x)), where p=p(x; alpha)
	 */

#define _PROB(xx, EE) (1.0-pow(EE*exp(xx)/(1.0+EE*exp(xx)), alpha))

	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx], E = ds->data_observations.E[idx], normc = my_gsl_sf_lnfact((int) y),
	    alpha = map_exp_forward(ds->data_observations.zeroinflated_alpha_intern[thread_id][0], MAP_FORWARD, NULL), mu, log_mu, p, lambda;

	LINK_INIT;
	if ((int) y == 0) {
		if (m > 0) {
			for (i = 0; i < m; i++) {
				p = _PROB(x[i] + off, E);
				if (gsl_isnan(p)) {
					// P(p);
					// P(x[i]+off);
					logll[i] = 0.0;
				} else {
					lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
					mu = E * lambda;
					log_mu = log(mu);

					// better expression I hope

					if (ISEQUAL(p, 1.0)) {
						logll[i] = 0.0;
					} else if (p < 1e-10) {
						logll[i] = 0 * log_mu - mu - normc;
					} else {
						logll[i] =
						    0 * log_mu - mu - normc + log(p / (gsl_ran_poisson_pdf((unsigned int) y, mu)) + (1.0 - p));
					}
					// logll[i] = log(p + (1.0 - p) * gsl_ran_poisson_pdf((unsigned int) y, mu));

					/*
					 * if all fails... 
					 */
					if (gsl_isnan(logll[i])) {
						P(p);
						P(logll[i]);
						P(x[i] + off);
						fprintf(stderr, "inla.c: Don't know what to do. Please report problem...");
						exit(EXIT_FAILURE);
					}
				}
			}
		} else {
			GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
			for (i = 0; i < -m; i++) {
				lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				mu = E * lambda;
				p = _PROB(x[i] + off, E);
				logll[i] = p + (1.0 - p) * gsl_cdf_poisson_P((unsigned int) y, mu);
			}
		}
	} else {
		if (m > 0) {
			for (i = 0; i < m; i++) {
				p = _PROB(x[i] + off, E);
				if (gsl_isnan(p)) {
					logll[i] = -INLA_REAL_BIG;
				} else {
					lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
					mu = E * lambda;
					log_mu = log(mu);
					logll[i] = LOG_1mp(p) + y * log_mu - mu - normc;
				}
			}
		} else {
			GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
			for (i = 0; i < -m; i++) {
				lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				mu = E * lambda;
				p = _PROB(x[i] + off, E);
				logll[i] = p + (1.0 - p) * gsl_cdf_poisson_P((unsigned int) y, mu);
			}
		}
	}

	for (i = 0; i < IABS(m); i++) {
		if (gsl_isinf(logll[i]))
			logll[i] = ((double) INLA_REAL_BIG) * gsl_isinf(logll[i]);
	}

	LINK_END;
#undef _PROB
	return GMRFLib_SUCCESS;
}

int loglikelihood_poisson_special1(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec),
				   double *y_cdf, void *arg, char **UNUSED(arg_str))
{
	/*
	 * poisson special 1 : y ~ p*1[y=1] + (1-p)*Poisson(E*exp(x) | y > 0)
	 */
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx], E = ds->data_observations.E[idx], normc = my_gsl_sf_lnfact((int) y),
	    p = map_probability_forward(ds->data_observations.prob_intern[thread_id][0], MAP_FORWARD, NULL), mu, p0, pp0;

	LINK_INIT;

	if (m > 0) {
		if (y == 1.0) {
			for (i = 0; i < m; i++) {
				mu = E * PREDICTOR_INVERSE_LINK(x[i] + off);
				p0 = exp(-mu);
				pp0 = 1.0 - p0;
				logll[i] = log(p + (1.0 - p) / pp0 * exp(y * log(mu) - mu - normc));
			}
		} else {
			for (i = 0; i < m; i++) {
				mu = E * PREDICTOR_INVERSE_LINK(x[i] + off);
				p0 = exp(-mu);
				pp0 = 1.0 - p0;
				logll[i] = log((1.0 - p) / pp0) + (y * log(mu) - mu - normc);
			}
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
		if (y < 1.0) {
			for (i = 0; i < -m; i++) {
				logll[i] = 0.0;
			}
		} else {
			for (i = 0; i < -m; i++) {
				mu = E * PREDICTOR_INVERSE_LINK(x[i] + off);
				p0 = exp(-mu);
				pp0 = 1.0 - p0;
				logll[i] = p + (1.0 - p) * (gsl_cdf_poisson_P((unsigned int) y, mu) - p0) / pp0;
			}
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

double exp_taylor(double x, double x0, int order)
{
	int i;
	double val = 1.0, xx = x - x0;

	assert(order > 0);
	for (i = order; i >= 1; i--) {
		val = 1.0 + val * xx / (double) i;
	}
	val *= exp(x0);
	return val;
}

double dexp_taylor(double x, double x0, int order)
{
	return exp_taylor(x, x0, order - 1);
}

double ddexp_taylor(double x, double x0, int order)
{
	return exp_taylor(x, x0, order - 2);
}

int loglikelihood_logperiodogram(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec),
				 double *UNUSED(y_cdf), void *arg, char **UNUSED(arg_str))
{
	/*
	 * y = x -log(2) + log\chi^2
	 *
	 * use the exact expression or a taylor-series for the exp-term of a given order ??
	 */
	int i;

	if (m == 0) {
		return 0;
	}
	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx], v, ypred;

	LINK_INIT;
	if (m > 0) {
		for (i = 0; i < m; i++) {
			ypred = PREDICTOR_INVERSE_LINK(x[i] + off);
			v = y - ypred + M_LN2;
			logll[i] = -M_LN2 + v - 0.5 * exp(v);
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_negative_binomial(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec),
				    double *y_cdf, void *arg, char **UNUSED(arg_str))
{
	/*
	 * y ~ NegativeBinomial(size, p) where E(y) = E*exp(x); same definition as in R and GSL, similar parameterisation as for the Poisson.
	 */

	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx];
	double E = ds->data_observations.E[idx];
	double S = ds->data_observations.S[idx];
	double size = (ds->variant == 0 ? 1.0 : (ds->variant == 1 ? E : S)) * exp(ds->data_observations.log_size[thread_id][0]);

	LINK_INIT;
	if (G_norm_const_compute[idx]) {
		double *cache = Calloc(2, double);
		G_norm_const_v[idx] = (void *) cache;
		cache[0] = my_gsl_sf_lnfact((int) y);
		cache[1] = y * log(E);
		G_norm_const_compute[idx] = 0;
	}
	double *cache = (double *) G_norm_const_v[idx];
	double normc = cache[0];
	double y_log_E = cache[1];

	// there is a tradeoff in the computations, either we can use lgamma() functions or a sum of log()'s.
	static int calibrate = 1;
	static int ylim = 8L;
	if (calibrate) {
#pragma omp critical (Name_e618c7278d96ebc883f4ddb21a27897f1dbbed07)
		if (calibrate) {
			const int ntimes = 16L;
			const int verbose = 0;
			const int yymax = 100;
			double s[ntimes];
			for (int yy = 4, dy = 1; yy < yymax; yy += dy) {
				double t[] = { 0, 0 }, tmp0 = 0.0, tmp1 = 0.0;

				for (int time = 0; time < ntimes; time++) {
					s[time] = exp(2.0 * (GMRFLib_uniform() - 0.5));
				}

				t[0] = -GMRFLib_timer();
				for (int time = 0; time < ntimes; time++) {
					tmp0 += gsl_sf_lngamma(yy + s[time]) - gsl_sf_lngamma(s[time]);
				}
				t[0] += GMRFLib_timer();

				t[1] -= GMRFLib_timer();
				for (int time = 0; time < ntimes; time++) {
					double ss = s[time];
#pragma omp simd reduction(+: tmp1)
					for (int y1 = 0; y1 < yy; y1++) {
						tmp1 += log(y1 + ss);
					}
				}
				t[1] += GMRFLib_timer();

				assert(ABS(((tmp0 - tmp1)) / (tmp0 + tmp1)) < FLT_EPSILON);
				if (verbose) {
					printf("Optimize nbinomial: yy %d sf=%.3f prod=%.3f\n", yy, t[0] / (t[0] + t[1]), t[1] / (t[0] + t[1]));
				}
				if (t[1] > t[0]) {
					ylim = yy - dy / 2L;
					if (verbose) {
						printf("Optimize nbinomial: chose ylim = %1d\n", ylim);
					}
					break;
				}
			}
			calibrate = 0;
		}
	}


	if (m > 0) {
		// the expression lgamma(y+s)-lgamm(s) reduces using Gamma(1+z)=z*Gamma(z)
		double lnorm = -normc;
		if (y >= ylim) {
			lnorm += gsl_sf_lngamma(y + size) - gsl_sf_lngamma(size);
		} else {
#pragma omp simd reduction(+: lnorm)
			for (int yy = 0; yy < (int) y; yy++) {
				lnorm += log(yy + size);
			}
		}

		if (PREDICTOR_LINK_EQ(link_log)) {

			if (0) {
				// old code
				for (int i = 0; i < m; i++) {
					double lambda = exp(PREDICTOR_INVERSE_IDENTITY_LINK(x[i] + off));
					double mu = E * lambda;
					double p = size / (size + mu);
					logll[i] = lnorm + size * LOG_p(p) + y * LOG_1mp(p);
				}
			}
			// optimised code
			double lsize = log(size);
			double t2 = lnorm + size * log(size) + y_log_E;
			double t3 = -(size + y);

			if (PREDICTOR_SCALE == 1.0) {
				double tt2 = t2 + t3 * lsize;
				if (0) {
					double b = E / size;
					if (y > 0) {
#pragma omp simd
						for (int i = 0; i < m; i++) {
							double xx = x[i] + off;
							logll[i] = tt2 + t3 * log1p(b * exp(xx)) + y * xx;
						}
					} else {
#pragma omp simd
						for (int i = 0; i < m; i++) {
							logll[i] = tt2 + t3 * log1p(b * exp(x[i] + off));
						}
					}
				} else {
					size_t mm = GMRFLib_align_simple((size_t) m, sizeof(double));
					double xx[mm], ex[mm], lx[mm];

					GMRFLib_cdaddto(m, x, off, xx);
					GMRFLib_exp(m, xx, ex);
					GMRFLib_dscale(m, E / size, ex);
					GMRFLib_log1p(m, ex, lx);

					if (y > 0) {
						// logll[i] = tt2 + t3 * lx[i] + y * x[i]);
						GMRFLib_daxpbypcz(m, t3, lx, y, xx, tt2, logll);
					} else {
						// logll[i] = tt2 + t3 * lx[i];
						GMRFLib_daxpb(m, t3, lx, tt2, logll);
					}
				}
			} else {
				double lEsize = log(E) - lsize;
#pragma omp simd
				for (int i = 0; i < m; i++) {
					double xx = PREDICTOR_INVERSE_IDENTITY_LINK(x[i] + off);
					double t1 = lsize + log1p(exp(lEsize + xx));
					logll[i] = t2 + t3 * t1 + y * xx;
				}
			}
		} else {
#pragma omp simd
			for (int i = 0; i < m; i++) {
				double lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				double mu = E * lambda;
				double p = size / (size + mu);
				logll[i] = lnorm + size * LOG_p(p) + y * LOG_1mp(p);
			}
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
#pragma omp simd
		for (int i = 0; i < -m; i++) {
			double lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
			double mu = E * lambda;
			double p = size / (size + mu);
			logll[i] = gsl_cdf_negative_binomial_P((unsigned int) y, p, size);
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_negative_binomial_cen2(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec),
					 double *y_cdf, void *arg, char **UNUSED(arg_str))
{
	/*
	 * y ~ NegativeBinomial(size, p) where E(y) = E*exp(x); same definition as in R and GSL, similar parameterisation as for the Poisson.
	 *
	 * this version allow for cencoring in the interval cen_low, cen_high, like 'cenpoisson2'
	 */

	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	Data_section_tp *ds = (Data_section_tp *) arg;
	double size = 0.0;
	double y = ds->data_observations.y[idx];
	double E = ds->data_observations.E[idx];
	double S = ds->data_observations.S[idx];
	double cen_low = ds->data_observations.cen_low[idx];
	double cen_high = ds->data_observations.cen_high[idx];
	int int_low = (int) cen_low;
	int int_high = (int) cen_high;

	switch (ds->variant) {
	case 0:
	{
		size = exp(ds->data_observations.log_size[thread_id][0]);
	}
		break;
	case 1:
	{
		size = E * exp(ds->data_observations.log_size[thread_id][0]);
	}
		break;
	case 2:
	{
		size = S * exp(ds->data_observations.log_size[thread_id][0]);
	}
		break;
	default:
		GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
	}

	LINK_INIT;
	if (m > 0) {
		double lnorm = gsl_sf_lngamma(y + size) - gsl_sf_lngamma(size) - gsl_sf_lngamma(y + 1.0);
		if ((y >= int_low && int_low >= 0) && (int_high < 0 || y <= int_high)) {
			for (int i = 0; i < m; i++) {
				double lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				double mu = E * lambda;
				logll[i] = log(inla_negative_binomial_interval(size, mu, int_low, int_high));
			}
		} else {
			for (int i = 0; i < m; i++) {
				double lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				double mu = E * lambda;
				double p = size / (size + mu);
				logll[i] = lnorm + size * LOG_p(p) + y * LOG_1mp(p);
			}
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
		for (int i = 0; i < -m; i++) {
			double lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
			double mu = E * lambda;
			double p = size / (size + mu);
			logll[i] = gsl_cdf_negative_binomial_P((unsigned int) y, p, size);
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_zeroinflated_negative_binomial0(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx,
						  double *UNUSED(x_vec), double *y_cdf, void *arg, char **UNUSED(arg_str))
{
	/*
	 * y ~ NegativeBinomial(size, p) where E(y) = E*exp(x); same definition as in R and GSL, similar parameterisation as for the Poisson.  This version is
	 * zeroinflated type 0.
	 */

	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double size = exp(ds->data_observations.log_size[thread_id][0]);
	double p_zeroinflated = map_probability_forward(ds->data_observations.prob_intern[thread_id][0], MAP_FORWARD, NULL);
	double y = ds->data_observations.y[idx];
	double E = ds->data_observations.E[idx];
	double lnorm, mu, p, prob_y_is_zero, lambda;
	double cutoff = 1.0e-4;				       /* switch to Poisson if mu/size < cutoff */

	LINK_INIT;
	if (m > 0) {
		if ((int) y == 0) {
			for (i = 0; i < m; i++) {
				logll[i] = log(p_zeroinflated);
			}
		} else {
			/*
			 * this is constant for the NegativeBinomial 
			 */
			lnorm = gsl_sf_lngamma(y + size) - gsl_sf_lngamma(size) - gsl_sf_lnfact((int) y);

			for (i = 0; i < m; i++) {
				lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				mu = E * lambda;
				if (mu / size > cutoff) {
					/*
					 * NegativeBinomial 
					 */
					p = size / (size + mu);
					prob_y_is_zero = gsl_ran_negative_binomial_pdf((unsigned int) 0, p, size);

					logll[i] = log((1.0 - p_zeroinflated) / (1.0 - prob_y_is_zero))
					    + lnorm + size * LOG_p(p) + y * LOG_1mp(p);
				} else {
					/*
					 * the Poission limit 
					 */
					prob_y_is_zero = gsl_ran_poisson_pdf((unsigned int) 0, mu);
					logll[i] = log((1.0 - p_zeroinflated) / (1.0 - prob_y_is_zero))
					    + y * log(mu) - mu - my_gsl_sf_lnfact((int) y);
				}
			}
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
		double p0;
		for (i = 0; i < -m; i++) {
			lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
			mu = E * lambda;
			if (mu / size > cutoff) {
				/*
				 * NegativeBinomial 
				 */
				p = size / (size + mu);
				p0 = gsl_cdf_negative_binomial_P((unsigned int) 0, p, size);
				logll[i] = p_zeroinflated + (1.0 - p_zeroinflated) *
				    (gsl_cdf_negative_binomial_P((unsigned int) y, p, size) - p0) / (1.0 - p0);
			} else {
				/*
				 * The Poission limit 
				 */
				p0 = gsl_cdf_poisson_P((unsigned int) 0, mu);
				logll[i] = p_zeroinflated + (1.0 - p_zeroinflated) * (gsl_cdf_poisson_P((unsigned int) y, mu) - p0) / (1.0 - p0);
			}
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_zeroinflated_negative_binomial1(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx,
						  double *UNUSED(x_vec), double *y_cdf, void *arg, char **UNUSED(arg_str))
{
	/*
	 * y ~ NegativeBinomial(size, p) where E(y) = E*exp(x); same definition as in R and GSL, similar parameterisation as for the Poisson.  This version is
	 * zeroinflated type 1.
	 */

	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double size = exp(ds->data_observations.log_size[thread_id][0]);
	double p_zeroinflated = map_probability_forward(ds->data_observations.prob_intern[thread_id][0], MAP_FORWARD, NULL);
	double y = ds->data_observations.y[idx];
	double E = ds->data_observations.E[idx];
	double lnorm, mu, p, lambda;
	double cutoff = 1.0e-4;				       /* switch to Poisson if mu/size < cutoff */

	LINK_INIT;
	if (m > 0) {
		/*
		 * this is constant for the NegativeBinomial 
		 */
		lnorm = gsl_sf_lngamma(y + size) - gsl_sf_lngamma(size) - gsl_sf_lngamma(y + 1.0);

		if ((int) y == 0) {
			for (i = 0; i < m; i++) {
				lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				mu = E * lambda;
				if (mu / size > cutoff) {
					/*
					 * NegativeBinomial 
					 */
					p = size / (size + mu);
					logll[i] =
					    log(p_zeroinflated + (1.0 - p_zeroinflated) * gsl_ran_negative_binomial_pdf((unsigned int) y, p, size));
				} else {
					/*
					 * the Poission limit 
					 */
					logll[i] = log(p_zeroinflated + (1.0 - p_zeroinflated) * gsl_ran_poisson_pdf((unsigned int) y, mu));
				}
			}
		} else {
			for (i = 0; i < m; i++) {
				lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				mu = E * lambda;
				if (mu / size > cutoff) {
					/*
					 * NegativeBinomial 
					 */
					p = size / (size + mu);
					logll[i] = LOG_1mp(p_zeroinflated) + lnorm + size * LOG_p(p) + y * LOG_1mp(p);
				} else {
					/*
					 * the Poission limit 
					 */
					logll[i] = LOG_1mp(p_zeroinflated) + y * log(mu) - mu - my_gsl_sf_lnfact((int) y);
				}
			}
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
		for (i = 0; i < -m; i++) {
			lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
			mu = E * lambda;
			if (mu / size > cutoff) {
				/*
				 * NegativeBinomial 
				 */
				p = size / (size + mu);
				logll[i] = p_zeroinflated + (1.0 - p_zeroinflated) * gsl_cdf_negative_binomial_P((unsigned int) y, p, size);
			} else {
				/*
				 * The Poission limit 
				 */
				logll[i] = p_zeroinflated + (1.0 - p_zeroinflated) * gsl_cdf_poisson_P((unsigned int) y, mu);
			}
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_zeroinflated_negative_binomial1_strata2(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx,
							  double *UNUSED(x_vec), double *y_cdf, void *arg, char **UNUSED(arg_str))
{
	/*
	 * y ~ NegativeBinomial(size, p) where E(y) = E*exp(x); same definition as in R and GSL, similar parameterisation as for the Poisson.  This version is
	 * zeroinflated type 1.
	 */

	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	Data_section_tp *ds = (Data_section_tp *) arg;
	int i, strata = (int) ds->data_observations.strata[idx];
	double size = exp(ds->data_observations.log_size[thread_id][0]);
	double y = ds->data_observations.y[idx];
	double E = ds->data_observations.E[idx];
	double lnorm, mu, p, lambda;
	double cutoff = 1.0e-4;				       /* switch to Poisson if mu/size < cutoff */
	double p_zeroinflated;

	p_zeroinflated = map_probability_forward(ds->data_observations.probN_intern[strata][thread_id][0], MAP_FORWARD, NULL);

	LINK_INIT;
	if (m > 0) {
		/*
		 * this is constant for the NegativeBinomial 
		 */
		lnorm = gsl_sf_lngamma(y + size) - gsl_sf_lngamma(size) - gsl_sf_lngamma(y + 1.0);

		if ((int) y == 0) {
			for (i = 0; i < m; i++) {
				lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				mu = E * lambda;
				if (mu / size > cutoff) {
					/*
					 * NegativeBinomial 
					 */
					p = size / (size + mu);
					logll[i] =
					    log(p_zeroinflated + (1.0 - p_zeroinflated) * gsl_ran_negative_binomial_pdf((unsigned int) y, p, size));
				} else {
					/*
					 * the Poission limit 
					 */
					logll[i] = log(p_zeroinflated + (1.0 - p_zeroinflated) * gsl_ran_poisson_pdf((unsigned int) y, mu));
				}
			}
		} else {
			for (i = 0; i < m; i++) {
				lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				mu = E * lambda;
				if (mu / size > cutoff) {
					/*
					 * NegativeBinomial 
					 */
					p = size / (size + mu);
					logll[i] = LOG_1mp(p_zeroinflated) + lnorm + size * LOG_p(p) + y * LOG_1mp(p);
				} else {
					/*
					 * the Poission limit 
					 */
					logll[i] = LOG_1mp(p_zeroinflated) + y * log(mu) - mu - my_gsl_sf_lnfact((int) y);
				}
			}
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
		for (i = 0; i < -m; i++) {
			lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
			mu = E * lambda;
			if (mu / size > cutoff) {
				/*
				 * NegativeBinomial 
				 */
				p = size / (size + mu);
				logll[i] = p_zeroinflated + (1.0 - p_zeroinflated) * gsl_cdf_negative_binomial_P((unsigned int) y, p, size);
			} else {
				/*
				 * The Poission limit 
				 */
				logll[i] = p_zeroinflated + (1.0 - p_zeroinflated) * gsl_cdf_poisson_P((unsigned int) y, mu);
			}
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_zeroinflated_negative_binomial1_strata3(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx,
							  double *UNUSED(x_vec), double *y_cdf, void *arg, char **UNUSED(arg_str))
{
	/*
	 * y ~ NegativeBinomial(size, p) where E(y) = E*exp(x); same definition as in R and GSL, similar parameterisation as for the Poisson.  This version is
	 * zeroinflated type 1.
	 */

	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	Data_section_tp *ds = (Data_section_tp *) arg;
	int i, strata = (int) ds->data_observations.strata[idx];
	double size;
	double p_zeroinflated;
	double y = ds->data_observations.y[idx];
	double E = ds->data_observations.E[idx];
	double lnorm, mu, p, lambda;
	double cutoff = 1.0e-4;				       /* switch to Poisson if mu/size < cutoff */

	p_zeroinflated = map_probability_forward(ds->data_observations.prob_intern[thread_id][0], MAP_FORWARD, NULL);
	size = exp(ds->data_observations.log_sizes[strata][thread_id][0]);

	LINK_INIT;
	if (m > 0) {
		/*
		 * this is constant for the NegativeBinomial 
		 */
		lnorm = gsl_sf_lngamma(y + size) - gsl_sf_lngamma(size) - gsl_sf_lngamma(y + 1.0);

		if ((int) y == 0) {
			for (i = 0; i < m; i++) {
				lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				mu = E * lambda;
				if (mu / size > cutoff) {
					/*
					 * NegativeBinomial 
					 */
					p = size / (size + mu);
					logll[i] =
					    log(p_zeroinflated + (1.0 - p_zeroinflated) * gsl_ran_negative_binomial_pdf((unsigned int) y, p, size));
				} else {
					/*
					 * the Poission limit 
					 */
					logll[i] = log(p_zeroinflated + (1.0 - p_zeroinflated) * gsl_ran_poisson_pdf((unsigned int) y, mu));
				}
			}
		} else {
			for (i = 0; i < m; i++) {
				lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				mu = E * lambda;
				if (mu / size > cutoff) {
					/*
					 * NegativeBinomial 
					 */
					p = size / (size + mu);
					logll[i] = LOG_1mp(p_zeroinflated) + lnorm + size * LOG_p(p) + y * LOG_1mp(p);
				} else {
					/*
					 * the Poission limit 
					 */
					logll[i] = LOG_1mp(p_zeroinflated) + y * log(mu) - mu - my_gsl_sf_lnfact((int) y);
				}
			}
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
		for (i = 0; i < -m; i++) {
			lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
			mu = E * lambda;
			if (mu / size > cutoff) {
				/*
				 * NegativeBinomial 
				 */
				p = size / (size + mu);
				logll[i] = p_zeroinflated + (1.0 - p_zeroinflated) * gsl_cdf_negative_binomial_P((unsigned int) y, p, size);
			} else {
				/*
				 * The Poission limit 
				 */
				logll[i] = p_zeroinflated + (1.0 - p_zeroinflated) * gsl_cdf_poisson_P((unsigned int) y, mu);
			}
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_zeroinflated_negative_binomial2(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx,
						  double *UNUSED(x_vec), double *y_cdf, void *arg, char **UNUSED(arg_str))
{
	/*
	 * y ~ NegativeBinomial(size, p) where E(y) = E*exp(x); same definition as in R and GSL, similar parameterisation as for the Poisson.  This version is
	 * zeroinflated type 3.
	 */

	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double size = exp(ds->data_observations.log_size[thread_id][0]);
	double alpha = map_exp_forward(ds->data_observations.zeroinflated_alpha_intern[thread_id][0], MAP_FORWARD, NULL);
	double p_zeroinflated = 0.0;
	double y = ds->data_observations.y[idx];
	double E = ds->data_observations.E[idx];
	double lnorm, mu, p, lambda;
	double cutoff = 1.0e-4;				       /* switch to Poisson if mu/size < cutoff */
	double normc = my_gsl_sf_lnfact(y);

	LINK_INIT;
	if (m > 0) {
		/*
		 * this is constant for the NegativeBinomial 
		 */
		lnorm = gsl_sf_lngamma(y + size) - gsl_sf_lngamma(size) - gsl_sf_lngamma(y + 1.0);

		if ((int) y == 0) {
			for (i = 0; i < m; i++) {
				lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				mu = E * lambda;
				p_zeroinflated = 1.0 - pow(mu / (1.0 + mu), alpha);

				if (gsl_isnan(p_zeroinflated)) {
					logll[i] = -INLA_REAL_BIG;
				} else {
					if (mu / size > cutoff) {
						/*
						 * NegativeBinomial 
						 */
						p = size / (size + mu);
						logll[i] =
						    log(p_zeroinflated +
							(1.0 - p_zeroinflated) * gsl_ran_negative_binomial_pdf((unsigned int) y, p, size));
					} else {
						/*
						 * the Poission limit 
						 */
						logll[i] = log(p_zeroinflated + (1.0 - p_zeroinflated) * gsl_ran_poisson_pdf((unsigned int) y, mu));
					}
				}
			}
		} else {
			for (i = 0; i < m; i++) {
				lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				mu = E * lambda;
				p_zeroinflated = 1.0 - pow(mu / (1.0 + mu), alpha);
				if (gsl_isnan(p_zeroinflated)) {
					logll[i] = -INLA_REAL_BIG;
				} else {
					if (mu / size > cutoff) {
						/*
						 * NegativeBinomial 
						 */
						p = size / (size + mu);
						logll[i] = LOG_1mp(p_zeroinflated) + lnorm + size * LOG_p(p) + y * LOG_1mp(p);
					} else {
						/*
						 * the Poission limit 
						 */
						logll[i] = LOG_1mp(p_zeroinflated) + y * log(mu) - mu - normc;
					}
				}
			}
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
		for (i = 0; i < -m; i++) {
			lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
			mu = E * lambda;
			p_zeroinflated = 1.0 - pow(mu / (1.0 + mu), alpha);
			if (mu / size > cutoff) {
				/*
				 * NegativeBinomial 
				 */
				p = size / (size + mu);
				logll[i] = p_zeroinflated + (1.0 - p_zeroinflated) * gsl_cdf_negative_binomial_P((unsigned int) y, p, size);
			} else {
				/*
				 * The Poission limit 
				 */
				logll[i] = p_zeroinflated + (1.0 - p_zeroinflated) * gsl_cdf_poisson_P((unsigned int) y, mu);
			}
		}
	}

	for (i = 0; i < IABS(m); i++) {
		if (gsl_isinf(logll[i]))
			logll[i] = ((double) INLA_REAL_BIG) * gsl_isinf(logll[i]);
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_binomial(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			   void *arg, char **UNUSED(arg_str))
{
	/*
	 * y ~ Binomial(n, p)
	 */
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	int status;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx];
	double n = ds->data_observations.nb[idx];
	double ny = n - y;
	/*
	 * this is a special case that should just return 0 or 1
	 */
	if (ISZERO(n) && ISZERO(y)) {
		if (m > 0) {
			for (int i = 0; i < m; i++) {
				logll[i] = 0.0;		       /* log(1) = 0 */
			}
		} else {
			for (int i = 0; i < -m; i++) {
				logll[i] = 1.0;
			}
		}
		return GMRFLib_SUCCESS;
	}

	/*
	 * this is the normal case...
	 */
	LINK_INIT;
	if (m > 0) {
		gsl_sf_result res = { 0, 0 };
		if (G_norm_const_compute[idx]) {
			if (ds->variant == 0) {
				// binomial
				status = gsl_sf_lnchoose_e((unsigned int) n, (unsigned int) y, &res);
			} else {
				// neg binomial
				status = gsl_sf_lnchoose_e((unsigned int) (n - 1.0), (unsigned int) (y - 1.0), &res);
			}
			assert(status == GSL_SUCCESS);
			G_norm_const[idx] = res.val;
			G_norm_const_compute[idx] = 0;
		}
		res.val = G_norm_const[idx];

		const int mkl_lim = 4L;
		int fast = (PREDICTOR_SCALE == 1.0);

		// special code for this case
		if (PREDICTOR_LINK_EQ(link_logit)) {
			if (ISZERO(y)) {
				if (m >= mkl_lim) {
					size_t mm = GMRFLib_align_simple((size_t) m, sizeof(double));
					double v_eta[mm], v_ee[mm], v_lee[mm];
					if (fast) {
#pragma omp simd
						for (int i = 0; i < m; i++) {
							v_eta[i] = x[i] + off;
						}
					} else {
#pragma omp simd
						for (int i = 0; i < m; i++) {
							v_eta[i] = PREDICTOR_INVERSE_IDENTITY_LINK(x[i] + off);
						}
					}
					GMRFLib_exp(m, v_eta, v_ee);
					GMRFLib_log1p(m, v_ee, v_lee);

#pragma omp simd
					for (int i = 0; i < m; i++) {
						logll[i] = res.val - ny * v_lee[i];
					}
				} else {
#pragma omp simd
					for (int i = 0; i < m; i++) {
						double eta = PREDICTOR_INVERSE_IDENTITY_LINK(x[i] + off);
						double ee = exp(eta);
						double log_1mp = -log1p(ee);
						logll[i] = res.val + ny * log_1mp;
					}
				}
			} else if (ISZERO(ny)) {
				if (m >= mkl_lim) {
					size_t mm = GMRFLib_align_simple((size_t) m, sizeof(double));
					double v_eta[mm], v_ee[mm], v_lee[mm];

					if (fast) {
#pragma omp simd
						for (int i = 0; i < m; i++) {
							v_eta[i] = -(x[i] + off);
						}
					} else {
#pragma omp simd
						for (int i = 0; i < m; i++) {
							v_eta[i] = -PREDICTOR_INVERSE_IDENTITY_LINK(x[i] + off);
						}
					}

					GMRFLib_exp(m, v_eta, v_ee);
					GMRFLib_log1p(m, v_ee, v_lee);
#pragma omp simd
					for (int i = 0; i < m; i++) {
						logll[i] = res.val - y * v_lee[i];
					}
				} else {
#pragma omp simd
					for (int i = 0; i < m; i++) {
						double eta = PREDICTOR_INVERSE_IDENTITY_LINK(x[i] + off);
						double ee = exp(-eta);
						double log_p = -log1p(ee);
						logll[i] = res.val + y * log_p;
					}
				}
			} else {
				if (m >= mkl_lim) {
					size_t mm = GMRFLib_align_simple((size_t) m, sizeof(double));
					double v_eta[mm], v_meta[mm], v_ee[mm], v_iee[mm], v_lee[mm], v_liee[mm];
#pragma omp simd
					for (int i = 0; i < m; i++) {
						v_eta[i] = PREDICTOR_INVERSE_IDENTITY_LINK(x[i] + off);
						v_meta[i] = -v_eta[i];
					}

					GMRFLib_exp(m, v_eta, v_ee);
					GMRFLib_log1p(m, v_ee, v_lee);
					GMRFLib_exp(m, v_meta, v_iee);
					GMRFLib_log1p(m, v_iee, v_liee);
#pragma omp simd
					for (int i = 0; i < m; i++) {
						logll[i] = res.val - ny * v_lee[i] - y * v_liee[i];
					}
				} else {
#pragma omp simd
					for (int i = 0; i < m; i++) {
						double eta = PREDICTOR_INVERSE_IDENTITY_LINK(x[i] + off);
						double ee = exp(eta);
						double log_1mp = -log1p(ee);
						double log_p = -log1p(1.0 / ee);
						logll[i] = res.val + y * log_p + ny * log_1mp;
					}
				}
			}
		} else if (PREDICTOR_LINK_EQ(link_cloglog)) {
			// this one is numerically unstable in the tails without special treatment
			if (ISZERO(y)) {
#pragma omp simd
				for (int i = 0; i < m; i++) {
					double eta = PREDICTOR_INVERSE_IDENTITY_LINK(x[i] + off);
					double log_1mp = link_log_1m_invcloglog(eta);
					logll[i] = res.val + ny * log_1mp;
				}
			} else if (ISZERO(ny)) {
#pragma omp simd
				for (int i = 0; i < m; i++) {
					double eta = PREDICTOR_INVERSE_IDENTITY_LINK(x[i] + off);
					double log_p = link_log_invcloglog(eta);
					logll[i] = res.val + y * log_p;
				}
			} else {
#pragma omp simd
				for (int i = 0; i < m; i++) {
					double eta = PREDICTOR_INVERSE_IDENTITY_LINK(x[i] + off);
					double log_p, log_1mp;
					link_log_invcloglog2(eta, &log_p, &log_1mp);
					logll[i] = res.val + y * log_p + ny * log_1mp;
				}
			}
		} else if (PREDICTOR_LINK_EQ(link_ccloglog)) {
			// this one is numerically unstable in the tails without special treatment
			if (ISZERO(y)) {
#pragma omp simd
				for (int i = 0; i < m; i++) {
					double eta = PREDICTOR_INVERSE_IDENTITY_LINK(x[i] + off);
					double log_1mp = link_log_1m_invccloglog(eta);
					logll[i] = res.val + ny * log_1mp;
				}
			} else if (ISZERO(ny)) {
#pragma omp simd
				for (int i = 0; i < m; i++) {
					double eta = PREDICTOR_INVERSE_IDENTITY_LINK(x[i] + off);
					double log_p = link_log_invccloglog(eta);
					logll[i] = res.val + y * log_p;
				}
			} else {
#pragma omp simd
				for (int i = 0; i < m; i++) {
					double eta = PREDICTOR_INVERSE_IDENTITY_LINK(x[i] + off);
					double log_p, log_1mp;
					link_log_invccloglog2(eta, &log_p, &log_1mp);
					logll[i] = res.val + y * log_p + ny * log_1mp;
				}
			}
		} else if (PREDICTOR_LINK_EQ(link_probit)) {
			// this one is numerically unstable in the tails without special treatment
			if (ISZERO(y)) {
#pragma omp simd
				for (int i = 0; i < m; i++) {
					double eta = PREDICTOR_INVERSE_IDENTITY_LINK(x[i] + off);
					double log_1mp = GMRFLib_log_cdfnorm(-eta);
					logll[i] = res.val + ny * log_1mp;
				}
			} else if (ISZERO(ny)) {
#pragma omp simd
				for (int i = 0; i < m; i++) {
					double eta = PREDICTOR_INVERSE_IDENTITY_LINK(x[i] + off);
					double log_p = GMRFLib_log_cdfnorm(eta);
					logll[i] = res.val + y * log_p;
				}
			} else {
#pragma omp simd
				for (int i = 0; i < m; i++) {
					double eta = PREDICTOR_INVERSE_IDENTITY_LINK(x[i] + off);
					double log_p = GMRFLib_log_cdfnorm(eta);
					double log_1mp = GMRFLib_log_cdfnorm(-eta);
					logll[i] = res.val + y * log_p + ny * log_1mp;
				}
			}
		} else {
#pragma omp simd
			for (int i = 0; i < m; i++) {
				// log(p) = LOG_p(p) and log(1-p) = LOG_1mp(p)
				double p = PREDICTOR_INVERSE_LINK(x[i] + off);
				logll[i] = res.val + y * LOG_p(p) + ny * LOG_1mp(p);
			}
		}
	} else {
		double *yy = (y_cdf ? y_cdf : &y);
		for (int i = 0; i < -m; i++) {
			double p = PREDICTOR_INVERSE_LINK((x[i] + off));
			p = DMIN(1.0, p);
			logll[i] = gsl_cdf_binomial_P((unsigned int) *yy, p, (unsigned int) n);
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_xbinomial(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			    void *arg, char **UNUSED(arg_str))
{
	/*
	 * y ~ xBinomial(n, p)
	 */
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	int status, i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx];
	double n = ds->data_observations.nb[idx], p;
	double p_scale = ds->data_observations.p_scale[idx];

	/*
	 * this is a special case that should just return 0 or 1
	 */
	if (ISZERO(y) && ISZERO(n)) {
		if (m > 0) {
			for (i = 0; i < m; i++) {
				logll[i] = 0.0;		       /* log(1) = 0 */
			}
		} else {
			for (i = 0; i < -m; i++) {
				logll[i] = 1.0;
			}
		}
		return GMRFLib_SUCCESS;
	}

	LINK_INIT;
	if (m > 0) {
		gsl_sf_result res = { 0, 0 };
		if (G_norm_const_compute[idx]) {
			if (ds->variant == 0) {
				// binomial
				status = gsl_sf_lnchoose_e((unsigned int) n, (unsigned int) y, &res);
			} else {
				// neg binomial
				status = gsl_sf_lnchoose_e((unsigned int) (n - 1.0), (unsigned int) (y - 1.0), &res);
			}
			assert(status == GSL_SUCCESS);
			G_norm_const[idx] = res.val;
			G_norm_const_compute[idx] = 0;
		}
		res.val = G_norm_const[idx];

		for (i = 0; i < m; i++) {
			p = p_scale * PREDICTOR_INVERSE_LINK(x[i] + off);
			p = DMIN(1.0 - FLT_EPSILON, p);
			logll[i] = res.val + y * LOG_p(p) + (n - y) * LOG_1mp(p);
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
		for (i = 0; i < -m; i++) {
			p = p_scale * PREDICTOR_INVERSE_LINK((x[i] + off));
			p = DMIN(1.0 - FLT_EPSILON, p);
			if (ds->variant == 0) {
				logll[i] = gsl_cdf_binomial_P((unsigned int) y, p, (unsigned int) n);
			} else {
				logll[i] = gsl_cdf_negative_binomial_P((unsigned int) (n - y), p, y);
			}
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_nbinomial2(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			     void *arg, char **UNUSED(arg_str))
{
	/*
	 * y ~ nBinomial2. y is the number of failures to get n successes with a success in the last trial
	 */
	int i;

	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}
	int status;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx];
	double n = ds->data_observations.nb[idx], p;

	/*
	 * this is the normal case...
	 */
	LINK_INIT;
	if (m > 0) {
		gsl_sf_result res = { 0, 0 };
		if (G_norm_const_compute[idx]) {
			status = gsl_sf_lnchoose_e((unsigned int) (y + n - 1.0), (unsigned int) (n - 1.0), &res);
			assert(status == GSL_SUCCESS);
			G_norm_const[idx] = res.val;
			G_norm_const_compute[idx] = 0;
		}
		res.val = G_norm_const[idx];

		for (i = 0; i < m; i++) {
			p = PREDICTOR_INVERSE_LINK(x[i] + off);
			p = TRUNCATE(p, 0.0, 1.0);
			logll[i] = res.val + y * LOG_1mp(p) + n * LOG_p(p);
		}
	} else {
		double *yy = (y_cdf ? y_cdf : &y);
		for (i = 0; i < -m; i++) {
			p = PREDICTOR_INVERSE_LINK((x[i] + off));
			p = TRUNCATE(p, 0.0, 1.0);
			logll[i] = gsl_cdf_negative_binomial_P((unsigned int) *yy, p, n);
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_nmix(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *UNUSED(y_cdf),
		       void *arg, char **UNUSED(arg_str))
{
	/*
	 * y ~ Binomial(n, p) * poisson(n, lambda), log(lambda) = X'beta
	 */
	if (m == 0) {
		return GMRFLib_SUCCESS;
	}

	int i, j, k;
	Data_section_tp *ds = (Data_section_tp *) arg;
	int n, nmax, ny;
	double *y = NULL, log_lambda, lambda, normc_poisson, fac, tt, tmp, p;

	assert(ds->data_observations.nmix_m > 0);
	for (i = 0, log_lambda = 0.0; i < ds->data_observations.nmix_m; i++) {
		log_lambda += ds->data_observations.nmix_beta[i][thread_id][0] * ds->data_observations.nmix_x[i][idx];
	}
	lambda = exp(log_lambda);

	LINK_INIT;

	static double **cy = NULL;
	static int *ncy = NULL;

	if (!cy) {
#pragma omp critical (Name_1d71960e99b4e67a36228891d1af914edbfd3dc3)
		{
			if (!cy) {
				ncy = Calloc(GMRFLib_CACHE_LEN(), int);
				cy = Calloc(GMRFLib_CACHE_LEN(), double *);
				for (i = 0; i < GMRFLib_CACHE_LEN(); i++) {
					ncy[i] = 8;
					cy[i] = Calloc(ncy[i], double);
				}
			}
		}
	}

	int id = 0;
	GMRFLib_CACHE_SET_ID(id);

	if (m > 0) {
		n = ds->data_observations.nmix_y[0][idx];
		assert(!gsl_isnan(ds->data_observations.nmix_y[0][idx]));
		ny = 1;
		for (i = 1; i > -1; i++) {
			if (gsl_isnan(ds->data_observations.nmix_y[i][idx]))
				break;
			ny++;
			n = IMAX(n, ds->data_observations.nmix_y[i][idx]);
		}
		normc_poisson = my_gsl_sf_lnfact(n);

		if (ny > ncy[id]) {
			ncy[id] = ny;
			cy[id] = Realloc(cy[id], ncy[id], double);
		}
		y = cy[id];

		for (i = 0; i < ny; i++) {
			y[i] = ds->data_observations.nmix_y[i][idx];
		}

		for (i = 0; i < m; i++) {
			gsl_sf_result res = { 0, 0 };
			p = PREDICTOR_INVERSE_LINK(x[i] + off);
			p = TRUNCATE(p, 0.0, 1.0);
			logll[i] = n * log_lambda - lambda - normc_poisson;
			for (j = 0; j < ny; j++) {
				gsl_sf_lnchoose_e((unsigned int) n, (unsigned int) y[j], &res);
				logll[i] += res.val + y[j] * LOG_p(p) + (n - y[j]) * LOG_1mp(p);
			}
			tt = lambda * inla_ipow(1.0 - p, ny);
			nmax = (int) DMAX(n + 50.0, DMIN(n + tt / 0.01, n + 500.0));	/* just to be sure */
			for (k = nmax, fac = 1.0; k > n; k--) {
				double kd = (double) k;
				for (j = 0, tmp = 1.0; j < ny; j++) {
					tmp *= kd / (kd - y[j]);
				}
				fac = 1.0 + fac * tt * tmp / kd;
			}
			logll[i] += log(fac);
		}
	} else {
		GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_nmixnb(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec),
			 double *UNUSED(y_cdf), void *arg, char **UNUSED(arg_str))
{
	/*
	 * y ~ Binomial(n, p) * NegBinom(n, mu=lambda, size=1/overdispersion), log(lambda) = X'beta
	 */
	if (m == 0) {
		return GMRFLib_SUCCESS;
	}

	int i, j, k;
	Data_section_tp *ds = (Data_section_tp *) arg;
	int n, nmax, ny;
	double *y = NULL, log_lambda, lambda, normc_nb, fac, tt, tmp, p, q, size;

	assert(ds->data_observations.nmix_m > 0);
	for (i = 0, log_lambda = 0.0; i < ds->data_observations.nmix_m; i++) {
		log_lambda += ds->data_observations.nmix_beta[i][thread_id][0] * ds->data_observations.nmix_x[i][idx];
	}
	lambda = exp(log_lambda);
	size = 1.0 / map_exp_forward(ds->data_observations.nmix_log_overdispersion[thread_id][0], MAP_FORWARD, NULL);

	LINK_INIT;

	static double **cy = NULL;
	static int *ncy = NULL;

	if (!cy) {
#pragma omp critical (Name_0b245bce3bb8f2007cc26fbbb141a5a0c7559165)
		{
			if (!cy) {
				ncy = Calloc(GMRFLib_CACHE_LEN(), int);
				cy = Calloc(GMRFLib_CACHE_LEN(), double *);
				for (i = 0; i < GMRFLib_CACHE_LEN(); i++) {
					ncy[i] = 8;
					cy[i] = Calloc(ncy[i], double);
				}
			}
		}
	}

	int id = 0;
	GMRFLib_CACHE_SET_ID(id);

	if (m > 0) {
		n = ds->data_observations.nmix_y[0][idx];
		assert(!gsl_isnan(ds->data_observations.nmix_y[0][idx]));
		ny = 1;
		for (i = 1; i > -1; i++) {
			if (gsl_isnan(ds->data_observations.nmix_y[i][idx]))
				break;
			ny++;
			n = IMAX(n, ds->data_observations.nmix_y[i][idx]);
		}
		normc_nb = gsl_sf_lngamma(n + size) - gsl_sf_lngamma(size) - my_gsl_sf_lnfact(n);

		if (ny > ncy[id]) {
			ncy[id] = ny;
			cy[id] = Realloc(cy[id], ncy[id], double);
		}
		y = cy[id];

		for (i = 0; i < ny; i++) {
			y[i] = ds->data_observations.nmix_y[i][idx];
		}

		for (i = 0; i < m; i++) {
			gsl_sf_result res = { 0, 0 };
			p = PREDICTOR_INVERSE_LINK(x[i] + off);
			p = TRUNCATE(p, 0.0, 1.0);
			q = size / (size + lambda);
			logll[i] = normc_nb + size * log(q) + n * LOG_1mp(q);
			for (j = 0; j < ny; j++) {
				gsl_sf_lnchoose_e((unsigned int) n, (unsigned int) y[j], &res);
				logll[i] += res.val + y[j] * LOG_p(p) + (n - y[j]) * LOG_1mp(p);
			}
			tt = lambda * sqrt(1.0 + lambda / size) * inla_ipow(1.0 - p, ny);
			nmax = (int) DMAX(n + 50.0, DMIN(n + tt / 0.01, n + 500.0));	/* just to be sure */
			tt = (1.0 - q) * inla_ipow(1.0 - p, ny);
			for (k = nmax, fac = 1.0; k > n; k--) {
				double kd = (double) k;
				for (j = 0, tmp = 1.0; j < ny; j++) {
					tmp *= kd / (kd - y[j]);
				}
				fac = 1.0 + fac * (kd + size - 1.0) * tt * tmp / kd;
			}
			logll[i] += log(fac);
		}
	} else {
		GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int inla_mix_int_quadrature_gaussian(int thread_id, double **x, double **w, int *n, void *arg)
{
	Data_section_tp *ds = (Data_section_tp *) arg;
	double prec = map_precision_forward(ds->data_observations.mix_log_prec_gaussian[thread_id][0], MAP_FORWARD, NULL);
	double sd = sqrt(1.0 / prec);
	double *xx = NULL, *ww = NULL, wmin;
	int i, j;
	GMRFLib_ghq(&xx, &ww, *n);

	wmin = MIX_INT_EPS * GMRFLib_max_value(ww, *n, NULL);
	*x = Calloc(*n, double);
	*w = Calloc(*n, double);
	for (i = j = 0; i < *n; i++) {
		if (ww[i] > wmin) {			       /* avoid to small weights */
			(*x)[j] = sd * xx[i];
			(*w)[j] = ww[i];
			j++;
		}
	}
	*n = j;

	return GMRFLib_SUCCESS;
}

int inla_mix_int_quadrature_loggamma(int UNUSED(thread_id), double **UNUSED(x), double **UNUSED(w), int *UNUSED(n), void *UNUSED(arg))
{
	char *msg = Strdup("This function is not yet implemented.");
	inla_error_general(msg);
	exit(1);
	return GMRFLib_SUCCESS;
}

int inla_mix_int_quadrature_mloggamma(int UNUSED(thread_id), double **UNUSED(x), double **UNUSED(w), int *UNUSED(n), void *UNUSED(arg))
{
	char *msg = Strdup("This function is not yet implemented.");
	inla_error_general(msg);
	exit(1);
	return GMRFLib_SUCCESS;
}

int inla_mix_int_simpson_gaussian(int thread_id, double **x, double **w, int *n, void *arg)
{
#define DENS(_x) exp(-0.5 * SQR(_x))
	Data_section_tp *ds = (Data_section_tp *) arg;
	double prec = map_precision_forward(ds->data_observations.mix_log_prec_gaussian[thread_id][0], MAP_FORWARD, NULL);
	double sd = sqrt(1.0 / prec);

	typedef struct {
		int n;					       /* is the requested length */
		int np;					       /* is the pruned length */
		double *__restrict x, *w;
	} lcache_t;

	static lcache_t **llcache = NULL;

	if (!llcache) {
#pragma omp critical (Name_f0fed30114239788701d492c4202a46b68cc060a)
		{
			if (!llcache) {
				llcache = Calloc(GMRFLib_CACHE_LEN(), lcache_t *);
			}
		}
	}
	int idx = 0;
	GMRFLib_CACHE_SET_ID(idx);
	if (!llcache) {
		llcache[idx] = Calloc(1, lcache_t);
	}

	lcache_t *lcache = llcache[idx];

	if (lcache->n != *n) {

		if (lcache->n > 0) {
			inla_error_general("Ask <help@r-inla.org> to rewrite inla_mix_int_simpson_gaussian()");
			exit(1);
		}

		Free(lcache->x);
		Free(lcache->w);
		lcache->n = *n;

		double *xx = Calloc(*n, double), *ww = Calloc(*n, double);
		double weight[2] = { 4.0, 2.0 }, limit = sqrt(-2.0 * log(MIX_INT_EPS)), dx = 2.0 * limit / (*n - 1.0), wmin;
		int i, j, np;

		xx[0] = -limit;
		xx[*n - 1] = limit;
		ww[0] = ww[*n - 1] = DENS(xx[0]);
		for (i = 1, j = 0; i < *n - 1; i++, j = (j + 1) % 2L) {
			xx[i] = xx[i - 1] + dx;
			ww[i] = weight[j] * DENS(xx[i]);
		}

		wmin = MIX_INT_EPS * GMRFLib_max_value(ww, *n, NULL);
		for (i = j = 0; i < *n; i++) {
			if (ww[i] > wmin) {
				xx[j] = xx[i];
				ww[j] = ww[i];
				j++;
			}
		}
		np = j;

		// make sure that integral(1) = 1 
		double corr = 0.0;
		for (i = 0; i < np; i++) {
			corr += ww[i];
		}
		corr = 1.0 / corr;
		for (i = 0; i < np; i++) {
			ww[i] *= corr;
		}

		lcache->x = xx;
		lcache->w = ww;
		lcache->np = np;
	}

	int i, np;

	np = lcache->np;
	*x = Calloc(np, double);
	*w = Calloc(np, double);
	for (i = 0; i < np; i++) {
		(*x)[i] = sd * lcache->x[i];
	}
	Memcpy(*w, lcache->w, np * sizeof(double));
	*n = np;

#undef DENS
	return GMRFLib_SUCCESS;
}

int inla_mix_int_simpson_loggamma(int thread_id, double **x, double **w, int *n, void *arg)
{
// Gamma(a,a) propto x^(a-1) * exp(-a*x). The density is normalized in any case, so we do not need to add the normalizing
// constant. z = log(x), so its the density for log(Gamma(a,a)). The '+1.0' is to stabilize it, so the modal value is around 0.
// shape = precision
#define DENS(_z, _a) exp( (_a)*((_z) -exp(_z) + 1.0) )

	Data_section_tp *ds = (Data_section_tp *) arg;
	double shape = map_precision_forward(ds->data_observations.mix_log_prec_loggamma[thread_id][0], MAP_FORWARD, NULL);

	typedef struct {
		int n, np;				       /* 'n' is the requested length and 'np' is the pruned length */
		double shape, *x, *w;
	} lcache_t;

	static lcache_t **llcache = NULL;

	if (!llcache) {
#pragma omp critical (Name_ba31fe1f20db2ec14750ea49f488e614a1920cb7)
		{
			if (!llcache) {
				llcache = Calloc(GMRFLib_CACHE_LEN(), lcache_t *);
			}
		}
	}
	int idx = 0;
	GMRFLib_CACHE_SET_ID(idx);

	if (!llcache[idx]) {
		llcache[idx] = Calloc(1, lcache_t);
	}
	lcache_t *lcache = llcache[idx];

	if (lcache->n != *n || lcache->shape != shape) {
		Free(lcache->x);
		double *work = Calloc(2 * *n, double);
		double *xx = NULL, *ww = NULL, weight[2] = { 4.0, 2.0 }, alpha = 0.001, low_limit, high_limit, dx, wmin;
		int i, j, np;

		low_limit = log(MATHLIB_FUN(qgamma) (alpha / 2.0, shape, 1.0 / shape, 1, 0));
		high_limit = log(MATHLIB_FUN(qgamma) (1.0 - alpha / 2.0, shape, 1.0 / shape, 1, 0));
		dx = (high_limit - low_limit) / (*n - 1.0);
		xx = work;				       /* use same storage */
		ww = work + *n;
		xx[0] = low_limit;
		ww[0] = DENS(low_limit, shape);
		xx[*n - 1] = high_limit;
		ww[*n - 1] = DENS(high_limit, shape);
		for (i = 1, j = 0; i < *n - 1; i++, j = (j + 1) % 2L) {
			xx[i] = xx[i - 1] + dx;
			ww[i] = weight[j] * DENS(xx[i], shape);
		}

		wmin = MIX_INT_EPS * GMRFLib_max_value(ww, *n, NULL);
		for (i = j = 0; i < *n; i++) {
			if (ww[i] > wmin) {
				xx[j] = xx[i];
				ww[j] = ww[i];
				j++;
			}
		}
		np = j;

		// make sure that integral(1) = 1 
		double corr = 0.0;
		for (i = 0; i < np; i++) {
			corr += ww[i];
		}
		corr = 1.0 / corr;
		for (i = 0; i < np; i++) {
			ww[i] *= corr;
		}

		lcache->n = *n;
		lcache->np = np;
		lcache->shape = shape;
		lcache->x = xx;
		lcache->w = ww;
	}

	*n = lcache->np;
	*x = Calloc(*n, double);
	*w = Calloc(*n, double);
	Memcpy(*x, lcache->x, *n * sizeof(double));
	Memcpy(*w, lcache->w, *n * sizeof(double));

#undef DENS
	return GMRFLib_SUCCESS;
}

int inla_mix_int_simpson_mloggamma(int thread_id, double **x, double **w, int *n, void *arg)
{
	inla_mix_int_simpson_loggamma(thread_id, x, w, n, arg);
	for (int i = 0; i < *n; i++) {
		/*
		 * just swap the sign
		 */
		(*x)[i] = -(*x)[i];
	}
	return GMRFLib_SUCCESS;
}

int loglikelihood_mix_loggamma(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *x_vec, double *y_cdf,
			       void *arg, char **arg_str)
{
	return (loglikelihood_mix_core
		(thread_id, logll, x, m, idx, x_vec, y_cdf, arg, inla_mix_int_quadrature_loggamma, inla_mix_int_simpson_loggamma, arg_str));
}

int loglikelihood_mix_mloggamma(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *x_vec, double *y_cdf,
				void *arg, char **arg_str)
{
	return (loglikelihood_mix_core
		(thread_id, logll, x, m, idx, x_vec, y_cdf, arg, inla_mix_int_quadrature_mloggamma, inla_mix_int_simpson_mloggamma, arg_str));
}

int loglikelihood_mix_gaussian(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *x_vec, double *y_cdf,
			       void *arg, char **arg_str)
{
	return (loglikelihood_mix_core
		(thread_id, logll, x, m, idx, x_vec, y_cdf, arg, inla_mix_int_quadrature_gaussian, inla_mix_int_simpson_gaussian, arg_str));
}

int loglikelihood_mix_core(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *x_vec, double *y_cdf, void *arg,
			   int (*func_quadrature)(int, double **, double **, int *, void *arg),
			   int(*func_simpson)(int, double **, double **, int *, void *arg), char **arg_str)
{
	Data_section_tp *ds =(Data_section_tp *) arg;
	if (m == 0) {
		if (arg) {
			return (ds->mix_loglikelihood(thread_id, NULL, NULL, 0, 0, NULL, NULL, arg, arg_str));
		} else {
			return (GMRFLib_LOGL_COMPUTE_CDF);
		}
	}

	int i, k, kk, mm, np;
	double *val = NULL, val_max, sum, *xx = NULL, *ll = NULL, *storage = NULL, *points = NULL, *weights = NULL;

	np = ds->mix_npoints;				       /* 'np' is the maximum number of points */
	switch (ds->mix_integrator) {
	case MIX_INT_QUADRATURE:
	{
		if (func_quadrature) {
			/*
			 * we can get back few points, as small weights are removed
			 */
			func_quadrature(thread_id, &points, &weights, &np, arg);
		} else {
			assert(0 == 1);
		}
	}
		break;
	case MIX_INT_DEFAULT:
	case MIX_INT_SIMPSON:
	{
		if (func_simpson) {
			/*
			 * we can get back few points, as small weights are removed
			 */
			func_simpson(thread_id, &points, &weights, &np, arg);
		} else {
			assert(0 == 1);
		}
	}
		break;
	default:
		assert(0 == 1);
	}

	mm = np * IABS(m);
	storage = Calloc(np + 2 * mm, double);		       /* use just one longer vector */
	val = storage;
	xx = storage + np;
	ll = storage + np + mm;

	if (m > 0) {
		for (i = 0, kk = 0; i < m; i++) {
			for (k = 0; k < np; k++) {
				xx[kk++] = x[i] + points[k];
			}
		}
		assert(kk == mm);
		ds->mix_loglikelihood(thread_id, ll, xx, mm, idx, x_vec, NULL, arg, arg_str);
		for (i = 0, kk = 0; i < m; i++) {
			for (k = 0; k < np; k++) {
				val[k] = log(weights[k]) + ll[kk++];
			}
			val_max = GMRFLib_max_value(val, np, NULL);
			for (k = 0, sum = 0.0; k < np; k++) {
				if (!ISNAN(val[k])) {
					sum += exp(val[k] - val_max);
				}
			}
			assert(sum > 0.0);
			logll[i] = log(sum) + val_max;
		}
		assert(kk == mm);
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
		for (i = 0, kk = 0; i < -m; i++) {
			for (k = 0; k < np; k++) {
				xx[kk++] = x[i] + points[k];
			}
		}
		assert(kk == mm);
		ds->mix_loglikelihood(thread_id, ll, xx, mm, idx, x_vec, NULL, arg, arg_str);
		for (i = 0, kk = 0; i < -m; i++) {
			for (k = 0, sum = 0.0; k < np; k++) {
				sum += weights[k] * ll[kk++];
			}
			logll[i] = sum;
		}
		assert(kk == mm);
	}

	Free(storage);
	Free(points);
	Free(weights);

	return GMRFLib_SUCCESS;
}

int loglikelihood_cbinomial(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			    void *arg, char **UNUSED(arg_str))
{
	/*
	 * z ~ CBinomial(k, n, p) == Binomial(k, 1-(1-p)^n)
	 */
	int i;

	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}
	int status;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx], k = ds->data_observations.cbinomial_k[idx], p;
	int n = ds->data_observations.cbinomial_n[idx];

	LINK_INIT;
	if (m > 0) {
		gsl_sf_result res = { 0, 0 };
		status = gsl_sf_lnchoose_e((unsigned int) k, (unsigned int) y, &res);	/* Yes, its 'k' */
		assert(status == GSL_SUCCESS);
		for (i = 0; i < m; i++) {
			p = PREDICTOR_INVERSE_LINK(x[i] + off);
			p = 1.0 - inla_ipow(1.0 - p, n);
			p = TRUNCATE(p, 0.0, 1.0);
			if (ISEQUAL(p, 1.0)) {
				/*
				 * this is ok if we get a 0*log(0) expression for the reminder 
				 */
				if (k == (int) y) {
					logll[i] = res.val + y * LOG_p(p);
				} else {
					logll[i] = -INLA_REAL_BIG;
				}
			} else if (ISZERO(p)) {
				/*
				 * this is ok if we get a 0*log(0) expression for the reminder 
				 */
				if ((int) y == 0) {
					logll[i] = res.val + (k - y) * LOG_1mp(p);
				} else {
					logll[i] = -INLA_REAL_BIG;
				}
			} else {
				logll[i] = res.val + y * LOG_p(p) + (k - y) * LOG_1mp(p);
			}
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
		for (i = 0; i < -m; i++) {
			p = PREDICTOR_INVERSE_LINK((x[i] + off));
			p = 1.0 - inla_ipow(1.0 - p, n);
			p = TRUNCATE(p, 0.0, 1.0);
			logll[i] = gsl_cdf_binomial_P((unsigned int) y, p, (unsigned int) k);
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_zeroinflated_binomial0(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec),
					 double *y_cdf, void *arg, char **UNUSED(arg_str))
{
	/*
	 * zeroinflated Binomial : y ~ p*1[y=0] + (1-p) Binomial(n, p | y > 0), where logit(p) = x. 
	 */
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx], n = ds->data_observations.nb[idx],
	    p = map_probability_forward(ds->data_observations.prob_intern[thread_id][0], MAP_FORWARD, NULL), prob = 0.0;

	LINK_INIT;
	if ((int) y == 0) {
		/*
		 * this is just the point-mass at zero 
		 */
		if (m > 0) {
			for (i = 0; i < m; i++) {
				logll[i] = LOG_p(p);
			}
		} else {
			for (i = 0; i < -m; i++) {
				logll[i] = p;
			}
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);

		gsl_sf_result res = { 0, 0 };
		gsl_sf_lnchoose_e((unsigned int) n, (unsigned int) y, &res);
		double p0;
		if (m > 0) {
			for (i = 0; i < m; i++) {
				prob = PREDICTOR_INVERSE_LINK(x[i] + off);
				p0 = gsl_cdf_binomial_P((unsigned int) 0, prob, (unsigned int) n);
				logll[i] = LOG_1mp(p) + res.val + y * LOG_p(prob) + (n - y) * LOG_1mp(prob) - LOG_1mp(p0);
			}
		} else {
			for (i = 0; i < -m; i++) {
				prob = PREDICTOR_INVERSE_LINK(x[i] + off);
				p0 = gsl_cdf_binomial_P((unsigned int) 0, prob, (unsigned int) n);
				logll[i] = p + (1.0 - p) * (gsl_cdf_binomial_P((unsigned int) y, prob, (unsigned int) n) - p0) / (1.0 - p0);
			}
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_zeroinflated_binomial1(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec),
					 double *y_cdf, void *arg, char **UNUSED(arg_str))
{
	/*
	 * zeroinflated Binomial : y ~ p*1[y=0] + (1-p)*Binomial(n, p), where logit(p) = x. 
	 */
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx], n = ds->data_observations.nb[idx],
	    p = map_probability_forward(ds->data_observations.prob_intern[thread_id][0], MAP_FORWARD, NULL), prob = 0.0, logA, logB;

	gsl_sf_result res = { 0, 0 };
	gsl_sf_lnchoose_e((unsigned int) n, (unsigned int) y, &res);

	LINK_INIT;
	if ((int) y == 0) {
		if (m > 0) {
			for (i = 0; i < m; i++) {
				prob = PREDICTOR_INVERSE_LINK(x[i] + off);
				logA = LOG_p(p);
				logB = LOG_1mp(p) + res.val + y * LOG_p(prob) + (n - y) * LOG_1mp(prob);
				// logll[i] = log(p + (1.0 - p) * gsl_ran_binomial_pdf((unsigned int) y, prob, (unsigned int) n));
				logll[i] = GMRFLib_logsum(logA, logB);
			}
		} else {
			GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
			for (i = 0; i < -m; i++) {
				prob = PREDICTOR_INVERSE_LINK(x[i] + off);
				logll[i] = p + (1.0 - p) * gsl_cdf_binomial_P((unsigned int) y, prob, (unsigned int) n);
			}
		}
	} else {
		if (m > 0) {
			for (i = 0; i < m; i++) {
				prob = PREDICTOR_INVERSE_LINK(x[i] + off);
				logll[i] = LOG_1mp(p) + res.val + y * LOG_p(prob) + (n - y) * LOG_1mp(prob);
			}
		} else {
			GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
			for (i = 0; i < -m; i++) {
				prob = PREDICTOR_INVERSE_LINK(x[i] + off);
				logll[i] = p + (1.0 - p) * gsl_cdf_binomial_P((unsigned int) y, prob, (unsigned int) n);
			}
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_zeroinflated_binomial2(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec),
					 double *y_cdf, void *arg, char **UNUSED(arg_str))
{
	/*
	 * zeroinflated Binomial : y ~ prob*1[y=0] + (1-prob)*Binomial(n, p), where logit(p) = x, and prob = 1-p^alpha.
	 */
#define _PROB(xx) (exp(xx)/(1.0+exp(xx)))
#define _PROBZERO(xx) (1.0-pow(_PROB(xx), alpha))

	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx], n = ds->data_observations.nb[idx], pzero, p,
	    alpha = map_exp_forward(ds->data_observations.zeroinflated_alpha_intern[thread_id][0], MAP_FORWARD, NULL), logA, logB;

	gsl_sf_result res = { 0, 0 };
	gsl_sf_lnchoose_e((unsigned int) n, (unsigned int) y, &res);

	LINK_INIT;
	if ((int) y == 0) {
		if (m > 0) {
			for (i = 0; i < m; i++) {
				pzero = _PROBZERO(x[i] + off);
				p = PREDICTOR_INVERSE_LINK(x[i] + off);
				if (gsl_isinf(pzero) || gsl_isinf(p)) {
					logll[i] = -INLA_REAL_BIG;
				} else {
					if (ISZERO(pzero)) {
						logll[i] = res.val + y * LOG_p(p) + (n - y) * LOG_1mp(p);
					} else {
						logA = log(pzero);
						logB = LOG_1mp(pzero) + res.val + y * LOG_p(p) + (n - y) * LOG_1mp(p);
						// logll[i] = log(pzero + (1.0 - pzero) * gsl_ran_binomial_pdf((unsigned int) y, p, 
						// 
						// 
						// 
						// 
						// (unsigned int) n));
						logll[i] = GMRFLib_logsum(logA, logB);
					}
				}
			}
		} else {
			GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
			for (i = 0; i < -m; i++) {
				pzero = _PROBZERO(x[i] + off);
				p = PREDICTOR_INVERSE_LINK(x[i] + off);
				if (gsl_isinf(pzero) || gsl_isinf(p)) {
					logll[i] = -INLA_REAL_BIG;
				} else {
					logll[i] = pzero + (1.0 - pzero) * gsl_cdf_binomial_P((unsigned int) y, p, (unsigned int) n);
				}
			}
		}
	} else {
		if (m > 0) {
			for (i = 0; i < m; i++) {
				pzero = _PROBZERO(x[i] + off);
				p = PREDICTOR_INVERSE_LINK(x[i] + off);
				if (gsl_isinf(pzero) || gsl_isinf(p)) {
					logll[i] = -INLA_REAL_BIG;
				} else {
					logll[i] = LOG_1mp(pzero) + res.val + y * LOG_p(p) + (n - y) * LOG_1mp(p);
				}
			}
		} else {
			GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
			for (i = 0; i < -m; i++) {
				pzero = _PROBZERO(x[i] + off);
				p = PREDICTOR_INVERSE_LINK(x[i] + off);
				if (gsl_isinf(pzero) || gsl_isinf(p)) {
					logll[i] = -INLA_REAL_BIG;
				} else {
					logll[i] = pzero + (1.0 - pzero) * gsl_cdf_binomial_P((unsigned int) y, p, (unsigned int) n);
				}
			}
		}
	}

	LINK_END;
#undef _PROB
#undef _PROBZERO
	return GMRFLib_SUCCESS;
}

#pragma omp declare simd
double GMRFLib_logsum(double lA, double lB)
{
	// evaluate log( exp(lA) + exp(lB) ) in a safe way 
	// return (lA > lB ? lA + log1p(exp(lB - lA)) : lB + log1p(exp(lA - lB)));

	return (fmax(lA, lB) + log1p(exp(-fabs(lB - lA))));
}

int loglikelihood_zero_n_inflated_binomial2(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec),
					    double *UNUSED(y_cdf), void *arg, char **UNUSED(arg_str))
{
	/*
	 * zeroNinflated Binomial : see doc from JS.
	 */
#define _P1(xx_) pow(1.0/(1.0+exp(-(xx_))), alpha1)
#define _P2(xx_) pow(1.0 - 1.0/(1.0+exp(-(xx_))), alpha2)

	if (m == 0) {
		return GMRFLib_SUCCESS;
	}

	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx], n = ds->data_observations.nb[idx],
	    alpha1 = map_exp_forward(ds->data_observations.zero_n_inflated_alpha1_intern[thread_id][0], MAP_FORWARD, NULL),
	    alpha2 = map_exp_forward(ds->data_observations.zero_n_inflated_alpha2_intern[thread_id][0], MAP_FORWARD, NULL), p, p1, p2, logA, logB;

	gsl_sf_result res = { 0, 0 };
	gsl_sf_lnchoose_e((unsigned int) n, (unsigned int) y, &res);

	LINK_INIT;
	if ((int) y == 0) {
		if (m > 0) {
			for (i = 0; i < m; i++) {
				p1 = _P1(x[i] + off);
				p2 = _P2(x[i] + off);
				p = PREDICTOR_INVERSE_LINK(x[i] + off);
				if (ISINF(p1) || ISINF(p2) || ISINF(p)) {
					logll[i] = -INLA_REAL_BIG;
				} else {
					if (ISZERO(1.0 - p1)) {
						logll[i] = LOG_p(p2) + res.val + y * LOG_p(p) + (n - y) * LOG_1mp(p);
					} else {
						logA = LOG_1mp(p1) + LOG_p(p2);
						logB = LOG_p(p1) + LOG_p(p2) + res.val + y * LOG_p(p) + (n - y) * LOG_1mp(p);
						// logll[i] = log((1.0 - p1) * p2 + p1 * p2 * gsl_ran_binomial_pdf((unsigned int)
						// y, p, (unsigned int) n));
						logll[i] = GMRFLib_logsum(logA, logB);
					}
				}
			}
		} else {
			GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
		}
	} else if ((int) y == (int) n) {
		if (m > 0) {
			for (i = 0; i < m; i++) {
				p1 = _P1(x[i] + off);
				p2 = _P2(x[i] + off);
				p = PREDICTOR_INVERSE_LINK(x[i] + off);
				if (ISINF(p1) || ISINF(p2) || ISINF(p)) {
					logll[i] = -INLA_REAL_BIG;
				} else {
					if (ISZERO(1.0 - p2)) {
						logll[i] = LOG_p(p1) + res.val + y * LOG_p(p) + (n - y) * LOG_1mp(p);
					} else {
						logA = LOG_1mp(p2) + LOG_p(p1);
						logB = LOG_p(p1) + LOG_p(p2) + res.val + y * LOG_p(p) + (n - y) * LOG_1mp(p);
						// logll[i] = log((1.0 - p2) * p1 + p1 * p2 * gsl_ran_binomial_pdf((unsigned int)
						// y, p, (unsigned int) n));
						logll[i] = GMRFLib_logsum(logA, logB);
					}
				}
			}
		} else {
			GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
		}
	} else {
		if (m > 0) {
			for (i = 0; i < m; i++) {
				p1 = _P1(x[i] + off);
				p2 = _P2(x[i] + off);
				p = PREDICTOR_INVERSE_LINK(x[i] + off);
				if (ISINF(p1) || ISINF(p2) || ISINF(p)) {
					logll[i] = -INLA_REAL_BIG;
				} else {
					logll[i] = LOG_p(p1) + LOG_p(p2) + res.val + y * LOG_p(p) + (n - y) * LOG_1mp(p);
				}
			}
		} else {
			GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
		}
	}

	LINK_END;
#undef _P1
#undef _P2
	return GMRFLib_SUCCESS;
}

int loglikelihood_zero_n_inflated_binomial3(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec),
					    double *UNUSED(y_cdf), void *arg, char **UNUSED(arg_str))
{
	/*
	 * zeroNinflated Binomial : see doc from JS.
	 */
#define _P0(_p) (pow(_p,        alpha0)/(1.0 + pow((_p), alpha0) + pow(1.0- (_p), alphaN)))
#define _PN(_p) (pow(1.0- (_p), alphaN)/(1.0 + pow((_p), alpha0) + pow(1.0- (_p), alphaN)))

	if (m == 0) {
		return GMRFLib_SUCCESS;
	}

	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx], n = ds->data_observations.nb[idx],
	    alpha0 = map_exp_forward(ds->data_observations.zero_n_inflated_alpha0_intern[thread_id][0], MAP_FORWARD, NULL),
	    alphaN = map_exp_forward(ds->data_observations.zero_n_inflated_alphaN_intern[thread_id][0], MAP_FORWARD, NULL), p, p0, pN, logA, logB;

	gsl_sf_result res = { 0, 0 };
	gsl_sf_lnchoose_e((unsigned int) n, (unsigned int) y, &res);

	LINK_INIT;
	if (m > 0) {
		for (i = 0; i < m; i++) {
			p = PREDICTOR_INVERSE_LINK(x[i] + off);
			p0 = _P0(p);
			pN = _PN(p);
			if (ISINF(p)) {
				logll[i] = -INLA_REAL_BIG;
			} else {
				logB = LOG_1mp((p0 + pN)) + res.val + y * LOG_p(p) + (n - y) * LOG_1mp(p);
				if ((int) y == 0) {
					logA = LOG_p(p0);
					logll[i] = GMRFLib_logsum(logA, logB);
				} else if ((int) y == (int) n) {
					logA = LOG_p(pN);
					logll[i] = GMRFLib_logsum(logA, logB);
				} else {
					logll[i] = logB;
				}
			}
		}
	} else {
		GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
	}

	LINK_END;
#undef _P0
#undef _PN
	return GMRFLib_SUCCESS;
}

int loglikelihood_gamma(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			void *arg, char **UNUSED(arg_str))
{
	/*
	 * Gamma
	 */

	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx];
	// 'scale' is not used for the survival version
	double s = (ds->data_observations.gamma_scale ? ds->data_observations.gamma_scale[idx] : 1.0);
	double phi_param = map_exp_forward(ds->data_observations.gamma_log_prec[thread_id][0], MAP_FORWARD, NULL);
	double phi = phi_param * s;
	double c = -gsl_sf_lngamma(phi) + (phi - 1.0) * log(y) + phi * log(phi);

	LINK_INIT;
	if (m > 0) {
		for (int i = 0; i < m; i++) {
			double mu = PREDICTOR_INVERSE_LINK(x[i] + off);
			logll[i] = c - phi * (log(mu) + y / mu);
		}
	} else {
		double yy = (y_cdf ? *y_cdf : y);
		for (int i = 0; i < -m; i++) {
			double mu = PREDICTOR_INVERSE_LINK(x[i] + off);
			double a = phi;
			double b = mu / phi;
			logll[i] = gsl_cdf_gamma_P(yy, a, b);
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_mgamma(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			 void *arg, char **UNUSED(arg_str))
{
	/*
	 * mGamma
	 */

	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx];
	double s = (ds->data_observations.gamma_scale ? ds->data_observations.gamma_scale[idx] : 1.0);
	double phi_param = map_exp_forward(ds->data_observations.gamma_log_prec[thread_id][0], MAP_FORWARD, NULL);
	double phi = phi_param * s;
	double delta = 0.5 * (sqrt(phi * (phi + 4.0)) + phi);
	double c = delta * log(y) - gsl_sf_lngamma(delta + 1.0);

	LINK_INIT;

	if (m > 0) {
		for (int i = 0; i < m; i++) {
			double mode = PREDICTOR_INVERSE_LINK(x[i] + off);
			double dmode = delta / mode;
			logll[i] = c + (1.0 + delta) * log(dmode) - y * dmode;
		}
	} else {
		double yy = (y_cdf ? *y_cdf : y);
		for (int i = 0; i < -m; i++) {
			double mode = PREDICTOR_INVERSE_LINK(x[i] + off);
			double dmode = delta / mode;
			logll[i] = gsl_cdf_gamma_P(yy, 1.0 + delta, 1.0 / dmode);
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_gammasurv(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *x_vec, double *y_cdf, void *arg,
			    char **arg_str)
{
	return (m ==
		0 ? GMRFLib_SUCCESS : loglikelihood_generic_surv(thread_id, logll, x, m, idx, x_vec, y_cdf, arg, loglikelihood_gamma, arg_str));
}

int loglikelihood_mgammasurv(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *x_vec, double *y_cdf, void *arg,
			     char **arg_str)
{
	return (m ==
		0 ? GMRFLib_SUCCESS : loglikelihood_generic_surv(thread_id, logll, x, m, idx, x_vec, y_cdf, arg, loglikelihood_mgamma, arg_str));
}

int loglikelihood_gammajw(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			  void *arg, char **UNUSED(arg_str))
{
	/*
	 * Gammajw
	 */

	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx];
	double ly = log(y);
	double mu;

	LINK_INIT;

	if (m > 0) {
		for (i = 0; i < m; i++) {
			mu = PREDICTOR_INVERSE_LINK(x[i] + off);
			logll[i] = -gsl_sf_lngamma(mu) + (mu - 1.0) * ly - y;
		}
	} else {
		double yy = (y_cdf ? *y_cdf : y);
		for (i = 0; i < -m; i++) {
			mu = PREDICTOR_INVERSE_LINK(x[i] + off);
			logll[i] = gsl_cdf_gamma_P(yy, mu, 1.0);
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_gammajwsurv(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *x_vec, double *y_cdf,
			      void *arg, char **arg_str)
{
	return (m ==
		0 ? GMRFLib_SUCCESS : loglikelihood_generic_surv(thread_id, logll, x, m, idx, x_vec, y_cdf, arg, loglikelihood_gammajw, arg_str));
}

int loglikelihood_gammacount(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			     void *arg, char **UNUSED(arg_str))
{
	/*
	 * Gammacount
	 */
#define _G(_alpha, _beta) gsl_cdf_gamma_P((_beta), (_alpha), 1.0)

	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx];
	double E = ds->data_observations.E[idx];
	double alpha = map_exp_forward(ds->data_observations.gammacount_log_alpha[thread_id][0], MAP_FORWARD, NULL);
	double beta, mu, p, logp;

	LINK_INIT;

	if (m > 0) {
		for (i = 0; i < m; i++) {
			mu = E * PREDICTOR_INVERSE_LINK(x[i] + off);
			beta = alpha * mu;
			p = _G(y * alpha, beta) - _G((y + 1.0) * alpha, beta);
			logp = LOG_p(p);
			// this can go in over/underflow...
			if (ISINF(logp) || ISNAN(logp)) {
				logll[i] = log(GSL_DBL_EPSILON) + PENALTY * SQR(x[i] + off);
			} else {
				logll[i] = logp;
			}
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
		for (i = 0; i < -m; i++) {
			mu = E * PREDICTOR_INVERSE_LINK(x[i] + off);
			beta = alpha * mu;
			logll[i] = _G((y + 1.0) * alpha, beta);
		}
	}

	LINK_END;
#undef _G
	return GMRFLib_SUCCESS;
}

int loglikelihood_qkumar(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			 void *arg, char **UNUSED(arg_str))
{
	/*
	 * qKumar-distr
	 */

	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx];
	double phi = map_prec_qkumar(ds->data_observations.qkumar_log_prec[thread_id][0], MAP_FORWARD, NULL);
	double q = ds->data_observations.quantile;
	double beta = LOG_1mp(q) / LOG_1mp(exp(-phi));
	double ibeta = 1.0 / beta;
	double tt1 = LOG_1mp(pow(1.0 - q, ibeta));

	LINK_INIT;
	if (m > 0) {
		double beta1 = beta - 1.0;
		double lbeta = log(beta);
		double ly = log(y);
		for (int i = 0; i < m; i++) {
			double kappa = PREDICTOR_INVERSE_LINK(x[i] + off);
			double alpha = tt1 / log(kappa);
			logll[i] = log(alpha) + lbeta + (alpha - 1.0) * ly + beta1 * LOG_1mp(pow(y, alpha));
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);
		for (int i = 0; i < -m; i++) {
			double kappa = PREDICTOR_INVERSE_LINK(x[i] + off);
			double alpha = tt1 / log(kappa);
			logll[i] = 1.0 - pow(1.0 - pow(y, alpha), beta);
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_gp(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf, void *arg,
		     char **UNUSED(arg_str))
{
	/*
	 * genPareto
	 */

	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx];
	double xi = map_interval(ds->data_observations.gp_intern_tail[thread_id][0], MAP_FORWARD,
				 (void *) (ds->data_observations.gp_tail_interval));
	double alpha = ds->data_observations.quantile;
	double q, sigma, fac;

	fac = xi / (pow(1.0 - alpha, -xi) - 1.0);

	LINK_INIT;
	if (m > 0) {
		for (i = 0; i < m; i++) {
			q = PREDICTOR_INVERSE_LINK(x[i] + off);
			sigma = q * fac;
			logll[i] = -log(sigma) - (1.0 / xi + 1.0) * log1p(xi * y / sigma);
		}
	} else {
		double yy = (y_cdf ? *y_cdf : y);
		for (i = 0; i < -m; i++) {
			q = PREDICTOR_INVERSE_LINK(x[i] + off);
			sigma = q * fac;
			logll[i] = 1.0 - pow(1.0 + xi * yy / sigma, -1.0 / xi);
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_dgp(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
		      void *arg, char **UNUSED(arg_str))
{
#define F(_y, _sigma, _xi) (1.0 - pow(1.0 + (_xi) * ((_y) + 1.0)/(_sigma), -1.0/(_xi)))
	/*
	 * discrete genPareto
	 */

	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx];
	double xi = map_interval(ds->data_observations.gp_intern_tail[thread_id][0], MAP_FORWARD,
				 (void *) (ds->data_observations.gp_tail_interval));
	double alpha = ds->data_observations.quantile;
	double q, sigma, fac;

	fac = xi / (pow(1.0 - alpha, -xi) - 1.0);

	LINK_INIT;
	if (m > 0) {
		for (i = 0; i < m; i++) {
			q = PREDICTOR_INVERSE_LINK(x[i] + off);
			sigma = q * fac;
			logll[i] = log(F(y, sigma, xi) - F(y - 1.01, sigma, xi));
		}
	} else {
		double yy = (y_cdf ? *y_cdf : y);
		for (i = 0; i < -m; i++) {
			q = PREDICTOR_INVERSE_LINK(x[i] + off);
			sigma = q * fac;
			logll[i] = F(yy, sigma, xi);
		}
	}
	LINK_END;

#undef F
	return GMRFLib_SUCCESS;
}

int loglikelihood_egp(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
		      void *arg, char **UNUSED(arg_str))
{
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx];
	double xi = map_interval(ds->data_observations.egp_intern_tail[thread_id][0], MAP_FORWARD,
				 (void *) (ds->data_observations.egp_tail_interval));
	double kappa = map_exp(ds->data_observations.egp_intern_shape[thread_id][0], MAP_FORWARD, NULL);
	double alpha = ds->data_observations.quantile;

	double eps = FLT_EPSILON;
	if (ABS(xi) < eps) {
		xi = DSIGN(xi) * eps;
	}

	LINK_INIT;
	assert(PREDICTOR_SCALE == 1.0);			       /* need to check... */

	double a = pow(1.0 - pow(alpha, 1.0 / kappa), -xi) - 1.0;
	double ia = 1.0 / a;
	double xii = -1.0 / xi;
	double lkappa = log(kappa);

	if (m > 0) {
		if (xi > 0.0) {
			for (int i = 0; i < m; i++) {
				double q = PREDICTOR_INVERSE_LINK(x[i] + off);
				double sigma = xi * q * ia;
				double yy = DMAX(DBL_EPSILON, 1.0 + xi * y / sigma);
				logll[i] = lkappa + (kappa - 1.0) * log1p(-pow(yy, xii)) - log(sigma) + (xii - 1.0) * log(yy);
			}
		} else {
			double f[] = { 0.90, 0.99 };
			double eta_c = log(-y * a);
			double eta_L, eta_H;
			if (eta_c < 0.0) {
				eta_H = f[0] * eta_c;
				eta_L = f[1] * eta_c;
			} else {
				eta_L = (1.0 / f[1]) * eta_c;
				eta_H = (1.0 / f[0]) * eta_c;
			}

			for (int i = 0; i < m; i++) {
				double xx = x[i] + off;
				if (xx > eta_H) {
					double q = PREDICTOR_INVERSE_LINK(xx);
					double sigma = xi * q * ia;
					double yy = 1.0 + xi * y / sigma;
					logll[i] = lkappa + (kappa - 1.0) * log1p(-pow(yy, xii)) - log(sigma) + (xii - 1.0) * log(yy);
				} else {
					double zz = (eta_H - xx) / (eta_H - eta_L);
					double eta = eta_H + (eta_L - eta_H) * zz / (1.0 + zz);
					double c0 = 0.0, c1 = 0.0, c2 = 0.0;
					double z = xx - eta;
					// needs 'alpha', 'kappa', 'xi', 'y' and 'eta' 
#include "egp-c012.h"
					logll[i] = c0 + z * (c1 + 0.5 * c2 * z);
				}
			}
		}
	} else {
		double yy = (y_cdf ? *y_cdf : y);
		for (int i = 0; i < -m; i++) {
			double q = PREDICTOR_INVERSE_LINK(x[i] + off);
			double sigma = xi * q * ia;
			logll[i] = pow(1.0 - pow(1.0 + xi * yy / sigma, xii), kappa);
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_beta(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
		       void *arg, char **UNUSED(arg_str))
{
	/*
	 * Beta : y ~ Beta(y; a, b) = BetaFunction(a,b)^{-1} y^{a-1} (1-y)^{b-1}. mu = a/(a+b), phi = a+b = exp(theta).
	 */

	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx];
	double w = ds->data_observations.beta_weight[idx];
	double phi = map_exp_forward(ds->data_observations.beta_precision_intern[thread_id][0], MAP_FORWARD, NULL) * w;
	double a, b, mu, lbeta;
	double censor_value = ds->data_observations.beta_censor_value;
	int no_censoring = (censor_value <= 0.0 || censor_value >= 0.5);

	LINK_INIT;
	if (m > 0) {
		for (i = 0; i < m; i++) {
			mu = PREDICTOR_INVERSE_LINK(x[i] + off);
			a = mu * phi;
			b = -mu * phi + phi;
			// If y is close to 0 then 'b' is tiny. Use the asymptotic expansion from `asympt(log(Beta(a,1/bb)), bb,
			// 1)'. If y is close to 1 then 'a'
			// is
			// tiny, do similarly
			if (DMIN(a, b) < INLA_REAL_SMALL) {
				lbeta = -log(DMIN(a, b));
			} else {
				lbeta = gsl_sf_lnbeta(a, b);
			}

			if (no_censoring) {
				// in most cases, we'll be here
				logll[i] = -lbeta + (a - 1.0) * log(y) + (b - 1.0) * LOG_1mp(y);
			} else {
				// if we have censoring, we have to be more careful
				if (y <= censor_value) {
					logll[i] = MATHLIB_FUN(pbeta) (censor_value, a, b, 1, 1);
				} else if (y < 1.0 - censor_value) {
					logll[i] = -lbeta + (a - 1.0) * log(y) + (b - 1.0) * LOG_1mp(y);
				} else {
					logll[i] = MATHLIB_FUN(pbeta) (1.0 - censor_value, a, b, 0, 1);
				}
			}
		}
	} else {
		double yy = (y_cdf ? *y_cdf : y);
		for (i = 0; i < -m; i++) {
			mu = PREDICTOR_INVERSE_LINK(x[i] + off);
			a = mu * phi;
			b = -mu * phi + phi;
			if (no_censoring) {
				logll[i] = gsl_cdf_beta_P(yy, a, b);
			} else {
				if (yy <= censor_value) {
					// use the expected prob instead
					logll[i] = MATHLIB_FUN(pbeta) (censor_value, a, b, 1, 0) / 2.0;
				} else if (yy < 1.0 - censor_value) {
					logll[i] = gsl_cdf_beta_P(yy, a, b);
				} else {
					// ... and also here
					logll[i] = 1.0 - MATHLIB_FUN(pbeta) (1.0 - censor_value, a, b, 0, 0) / 2.0;
				}
			}
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_betabinomial(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			       void *arg, char **UNUSED(arg_str))
{
	/*
	 * BetaBinomial : y ~ BetaBinomial(n, a, b), where logit(p) = a/(a+b), overdispertsion = 1/(a+b+1)
	 */
#define _LOGGAMMA_INT(xx) my_gsl_sf_lnfact(((xx) - 1))

	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	Data_section_tp *ds = (Data_section_tp *) arg;
	int y = (int) ds->data_observations.y[idx];
	int n = (int) ds->data_observations.nb[idx];

	double rho = map_probability_forward(ds->data_observations.betabinomial_overdispersion_intern[thread_id][0], MAP_FORWARD, NULL);
	double p, a, b;
	double normc;

	if (G_norm_const_compute[idx]) {
		G_norm_const[idx] = _LOGGAMMA_INT(n + 1) - _LOGGAMMA_INT(y + 1) - _LOGGAMMA_INT(n - y + 1);
		G_norm_const_compute[idx] = 0;
	}
	normc = G_norm_const[idx];

	LINK_INIT;
	if (m > 0) {
		// issues occur when x[i] is to large
		double p_upper = 0.999, xmax;
		double work[n];
		xmax = GMRFLib_max_value(x, m, NULL) + off;
		p = PREDICTOR_INVERSE_LINK(xmax);
		if (p < p_upper) {
			for (int i = 0; i < m; i++) {
				p = PREDICTOR_INVERSE_LINK(x[i] + off);
				a = p * (1.0 - rho) / rho;
				b = (p * rho - p - rho + 1.0) / rho;

				// joinig the beta-expressions, we can do better and be more robust
				// logll[i] = normc + gsl_sf_lnbeta(y + a, n - y + b) - gsl_sf_lnbeta(a, b);
				logll[i] = normc + my_betabinomial(y, n, a, b, work);
			}
		} else {
			// extrapolate linearly
			double xx[3], ll[3], h = 1.0E-4, diff, ddiff, dx;
			xx[1] = PREDICTOR_LINK(p_upper);
			xx[0] = xx[1] - h;
			xx[2] = xx[1] + h;
			for (int i = 0; i < 3; i++) {
				p = PREDICTOR_INVERSE_LINK(xx[i]);
				a = p * (1.0 - rho) / rho;
				b = (p * rho - p - rho + 1.0) / rho;
				ll[i] = normc + my_betabinomial(y, n, a, b, work);
			}
			diff = (ll[2] - ll[0]) / (2.0 * h);
			ddiff = (ll[2] - 2.0 * ll[1] + ll[0]) / SQR(h);
			diff = DMIN(0.0, diff);		       /* must have */
			ddiff = DMIN(0.0, ddiff);	       /* must have */
			for (int i = 0; i < m; i++) {
				dx = (x[i] + off) - xx[1];
				logll[i] = ll[1] + dx * diff + 0.5 * SQR(dx) * ddiff;
			}
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);

		for (int i = 0; i < -m; i++) {
			/*
			 * if 'n' is ``large'', then we need to do something approximate here... I guess we get the question back if so is the case. So we issue a
			 * warning for the ``extreme case'' only.
			 */
			static char give_warning = 1;
			if (n > 500 && give_warning) {
				give_warning = 0;
				printf("\n*** Warning ***  Version [%s]", GITCOMMIT);
				printf("\n*** Warning ***  The PIT calculations for the BetaBinomial can be time-consuming when Ntrials is large.");
				printf("\n*** Warning ***  Please contact <help@r-inla.org> if this becomes an issue.\n");
			}

			double normc2;

			p = PREDICTOR_INVERSE_LINK(x[i] + off);
			a = p * (1.0 - rho) / rho;
			b = (p * rho - p - rho + 1.0) / rho;
			normc2 = _LOGGAMMA_INT(n + 1) - my_gsl_sf_lnbeta(a, b);
			logll[i] = 0.0;

			if (y <= n / 2) {		       /* integer arithmetic is ok */
				for (int yy = y; yy >= 0; yy--) {
					logll[i] +=
					    exp(normc2 - _LOGGAMMA_INT(yy + 1) - _LOGGAMMA_INT(n - yy + 1) + my_gsl_sf_lnbeta(yy + a, n - yy + b));
				}
			} else {
				for (int yy = y + 1; yy <= n; yy++) {
					logll[i] +=
					    exp(normc2 - _LOGGAMMA_INT(yy + 1) - _LOGGAMMA_INT(n - yy + 1) + my_gsl_sf_lnbeta(yy + a, n - yy + b));
				}
				logll[i] = 1.0 - logll[i];
			}
		}
	}

	LINK_END;
#undef _LOGGAMMA_INT
	return GMRFLib_SUCCESS;
}

int loglikelihood_betabinomialna(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec),
				 double *y_cdf, void *arg, char **UNUSED(arg_str))
{
	/*
	 * BetaBinomialNA : y ~ BetaBinomial(n, a, b), where logit(p) = a/(a+b), overdispertsion = 1/(a+b+1), and use the normal
	 * approximation to it
	 */

	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx];
	double n = ds->data_observations.nb[idx];
	double s = ds->data_observations.betabinomialnb_scale[idx];
	double rho = map_probability_forward(ds->data_observations.betabinomial_overdispersion_intern[thread_id][0], MAP_FORWARD, NULL);
	double p, prec, lprec, ypred;

	LINK_INIT;
	if (m > 0) {
		for (i = 0; i < m; i++) {
			p = PREDICTOR_INVERSE_LINK(x[i] + off);
			ypred = n * p;
			prec = 1.0 / (n * p * (1.0 - p) * (1.0 + s * (n - 1.0) * rho));
			lprec = log(prec);
			logll[i] = LOG_NORMC_GAUSSIAN + 0.5 * (lprec - (SQR(ypred - y) * prec));
		}
	} else {
		double yy = (y_cdf ? *y_cdf : y);
		for (i = 0; i < -m; i++) {
			p = PREDICTOR_INVERSE_LINK(x[i] + off);
			ypred = n * p;
			prec = 1.0 / (n * p * (1.0 - p) * (1.0 + s * (n - 1.0) * rho));
			logll[i] = inla_cdf_normal_fast((yy - ypred) * sqrt(prec));
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_tweedie(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			  void *arg, char **UNUSED(arg_str))
{
	/*
	 * Tweedie
	 */

	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx];
	double w = (ds->data_observations.tweedie_w ? ds->data_observations.tweedie_w[idx] : 1.0);
	double interval[2] = { 1.0, 2.0 };
	double p = map_interval(ds->data_observations.tweedie_p_intern[thread_id][0], MAP_FORWARD, (void *) interval);
	double phi = map_exp_forward(ds->data_observations.tweedie_phi_intern[thread_id][0], MAP_FORWARD, NULL);

	phi /= w;
	LINK_INIT;

	static double **cmu = NULL;
	static int *ncmu = NULL;

	if (!cmu) {
#pragma omp critical (Name_f541b1464beaa9132d8c3f70fc8dc2de724ab8a5)
		{
			if (!cmu) {
				ncmu = Calloc(GMRFLib_CACHE_LEN(), int);
				cmu = Calloc(GMRFLib_CACHE_LEN(), double *);
				for (i = 0; i < GMRFLib_CACHE_LEN(); i++) {
					ncmu[i] = 8;
					cmu[i] = Calloc(ncmu[i], double);
				}
			}
		}
	}

	int id = 0;
	GMRFLib_CACHE_SET_ID(id);

	if (m > 0) {
		if (m > ncmu[id]) {
			ncmu[id] = m;
			cmu[id] = Realloc(cmu[id], ncmu[id], double);
		}
		double *mu = cmu[id];

		for (i = 0; i < m; i++) {
			mu[i] = PREDICTOR_INVERSE_LINK(x[i] + off);
		}
		dtweedie(m, y, mu, phi, p, logll);
	} else {
		double yy = (y_cdf ? *y_cdf : y);
		for (i = 0; i < -m; i++) {
			double mu = PREDICTOR_INVERSE_LINK(x[i] + off);
			logll[i] = ptweedie(yy, mu, phi, p);
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_zeroinflated_betabinomial0(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec),
					     double *y_cdf, void *arg, char **UNUSED(arg_str))
{
	/*
	 * zeroinflated BetaBinomial : y ~ BetaBinomial(n, a, b), where logit(p) = a/(a+b), overdispertsion = 1/(a+b+1)
	 */
#define _LOGGAMMA_INT(xx) my_gsl_sf_lnfact(((xx) - 1))

	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	int i, yzero = 0;
	Data_section_tp *ds = (Data_section_tp *) arg;
	int y = (int) ds->data_observations.y[idx];
	int n = (int) ds->data_observations.nb[idx];

	double rho = map_probability_forward(ds->data_observations.zeroinflated_rho_intern[thread_id][0], MAP_FORWARD, NULL);
	double pzero = map_probability_forward(ds->data_observations.prob_intern[thread_id][0], MAP_FORWARD, NULL);
	double p, a, b, prob_zero;
	double normc = _LOGGAMMA_INT(n + 1) - _LOGGAMMA_INT(y + 1) - _LOGGAMMA_INT(n - y + 1);
	double normc_zero = _LOGGAMMA_INT(n + 1) - _LOGGAMMA_INT(yzero + 1) - _LOGGAMMA_INT(n - yzero + 1);

	LINK_INIT;
	if (m > 0) {
		if (y == 0) {
			for (i = 0; i < m; i++) {
				logll[i] = LOG_p(pzero);
			}
		} else {
			for (i = 0; i < m; i++) {
				p = PREDICTOR_INVERSE_LINK(x[i] + off);
				a = p * (1.0 - rho) / rho;
				b = (p * rho - p - rho + 1.0) / rho;
				prob_zero = exp(normc_zero + gsl_sf_lnbeta(yzero + a, n - yzero + b) - gsl_sf_lnbeta(a, b));
				logll[i] = LOG_1mp(pzero) + normc + gsl_sf_lnbeta(y + a, n - y + b) - gsl_sf_lnbeta(a, b) - LOG_1mp(prob_zero);
			}
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);

		for (i = 0; i < -m; i++) {
			if (y == 0) {
				logll[i] = pzero;
			} else {
				int yy;
				double normc2;

				p = PREDICTOR_INVERSE_LINK(x[i] + off);
				a = p * (1.0 - rho) / rho;
				b = (p * rho - p - rho + 1.0) / rho;
				normc2 = _LOGGAMMA_INT(n + 1) - gsl_sf_lnbeta(a, b);
				prob_zero = normc + gsl_sf_lnbeta(yzero + a, n - yzero + b) - gsl_sf_lnbeta(a, b);
				logll[i] = 0.0;
				for (yy = y; yy > 0; yy--) {
					logll[i] +=
					    exp(normc2 - _LOGGAMMA_INT(yy + 1) - _LOGGAMMA_INT(n - yy + 1) + gsl_sf_lnbeta(yy + a, n - yy + b));
				}
				logll[i] = pzero + (1.0 - pzero) * logll[i] / (1.0 - prob_zero);
			}
		}
	}

	LINK_END;
#undef _LOGGAMMA_INT
	return GMRFLib_SUCCESS;
}

int loglikelihood_zeroinflated_betabinomial1(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec),
					     double *y_cdf, void *arg, char **UNUSED(arg_str))
{
	/*
	 * zeroinflated BetaBinomial : y ~ BetaBinomial(n, a, b), where logit(p) = a/(a+b), overdispertsion = 1/(a+b+1)
	 */
#define _LOGGAMMA_INT(xx) my_gsl_sf_lnfact(((xx) - 1))

	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	int y = (int) ds->data_observations.y[idx];
	int n = (int) ds->data_observations.nb[idx];

	double rho = map_probability_forward(ds->data_observations.zeroinflated_rho_intern[thread_id][0], MAP_FORWARD, NULL);
	double pzero = map_probability_forward(ds->data_observations.prob_intern[thread_id][0], MAP_FORWARD, NULL);
	double p, a, b, tmp;
	double normc = _LOGGAMMA_INT(n + 1) - _LOGGAMMA_INT(y + 1) - _LOGGAMMA_INT(n - y + 1);

	LINK_INIT;
	if (m > 0) {
		for (i = 0; i < m; i++) {
			p = PREDICTOR_INVERSE_LINK(x[i] + off);
			a = p * (1.0 - rho) / rho;
			b = (p * rho - p - rho + 1.0) / rho;
			tmp = LOG_1mp(pzero) + normc + gsl_sf_lnbeta(y + a, n - y + b) - gsl_sf_lnbeta(a, b);
			if (y == 0) {
				logll[i] = GMRFLib_log_apbex(pzero, tmp);
			} else {
				logll[i] = tmp;
			}
		}
	} else {
		GMRFLib_ASSERT(y_cdf == NULL, GMRFLib_ESNH);

		for (i = 0; i < -m; i++) {
			int yy;
			double normc2;

			p = PREDICTOR_INVERSE_LINK(x[i] + off);
			a = p * (1.0 - rho) / rho;
			b = (p * rho - p - rho + 1.0) / rho;
			normc2 = _LOGGAMMA_INT(n + 1) - gsl_sf_lnbeta(a, b);
			logll[i] = 0.0;
			for (yy = y; yy >= 0; yy--) {
				logll[i] += exp(normc2 - _LOGGAMMA_INT(yy + 1) - _LOGGAMMA_INT(n - yy + 1) + gsl_sf_lnbeta(yy + a, n - yy + b));
			}
			logll[i] = pzero + (1.0 - pzero) * logll[i];
		}
	}

	LINK_END;
#undef _LOGGAMMA_INT
	return GMRFLib_SUCCESS;
}

int loglikelihood_zeroinflated_betabinomial2(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec),
					     double *UNUSED(y_cdf), void *arg, char **UNUSED(arg_str))
{
	/*
	 * zeroinflated BetaBinomial : y ~ prob*1[y=0] + (1-prob)*BetaBinomial(n, p, delta), where logit(p) = x, and prob = 1-p^alpha.
	 */
#define _PROB(xx)         (exp(xx)/(1.0+exp(xx)))
#define _PROBZERO(xx)     (1.0-pow(_PROB(xx), alpha))
#define _LOGGAMMA(xx)     gsl_sf_lngamma(xx)
#define _LOGGAMMA_INT(xx) my_gsl_sf_lnfact(((xx) - 1))

	if (m == 0) {
		return 0;
	}

	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	int y = (int) ds->data_observations.y[idx];
	int n = (int) ds->data_observations.nb[idx];
	double pzero, p;
	double alpha = map_exp_forward(ds->data_observations.zeroinflated_alpha_intern[thread_id][0], MAP_FORWARD, NULL);
	double delta = map_exp_forward(ds->data_observations.zeroinflated_delta_intern[thread_id][0], MAP_FORWARD, NULL);
	double logA, logB;

	LINK_INIT;
	if ((int) y == 0) {
		for (i = 0; i < m; i++) {
			pzero = _PROBZERO(x[i] + off);
			p = PREDICTOR_INVERSE_LINK(x[i] + off);
			if (gsl_isinf(pzero) || gsl_isinf(p)) {
				logll[i] = -INLA_REAL_BIG;
			} else {

				logA = LOG_p(pzero);
				logB = LOG_1mp(pzero) + (_LOGGAMMA_INT(n + 1) - _LOGGAMMA_INT(y + 1) - _LOGGAMMA_INT(n - y + 1)
							 + _LOGGAMMA(delta * p + y) + _LOGGAMMA(n + delta * (1.0 - p) - y) - _LOGGAMMA(delta + n)
							 + _LOGGAMMA(delta) - _LOGGAMMA(delta * p) - _LOGGAMMA(delta * (1.0 - p)));
				logll[i] = GMRFLib_logsum(logA, logB);
			}
		}
	} else {
		for (i = 0; i < m; i++) {
			pzero = _PROBZERO(x[i] + off);
			p = PREDICTOR_INVERSE_LINK(x[i] + off);
			if (gsl_isinf(pzero) || gsl_isinf(p)) {
				logll[i] = -INLA_REAL_BIG;
			} else {
				logll[i] = LOG_1mp(pzero) + (_LOGGAMMA_INT(n + 1) - _LOGGAMMA_INT(y + 1) - _LOGGAMMA_INT(n - y + 1)
							     + _LOGGAMMA(delta * p + y) + _LOGGAMMA(n + delta * (1.0 - p) - y) -
							     _LOGGAMMA(delta + n)
							     + _LOGGAMMA(delta) - _LOGGAMMA(delta * p) - _LOGGAMMA(delta * (1.0 - p)));
			}
		}
	}

	LINK_END;
#undef _PROB
#undef _PROBZERO
#undef _LOGGAMMA
#undef _LOGGAMMA_INT

	return GMRFLib_SUCCESS;
}

int loglikelihood_exp(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
		      void *arg, char **arg_str)
{
	/*
	 * y ~ Exponential
	 */
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	Data_section_tp *ds = (Data_section_tp *) arg;
	int i;
	double y, lambda;

	y = ds->data_observations.y[idx];

	LINK_INIT;
	if (arg_str) {
		char *a = NULL;
		GMRFLib_sprintf(&a,
				"list(y = %.8g, family = \"exponential\", link.model = \"%1s\", linear.predictor = %.8g)", y, ds->link_model, x[0]);
		*arg_str = a;
		return GMRFLib_SUCCESS;
	}

	if (m > 0) {
		for (i = 0; i < m; i++) {
			lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
			logll[i] = log(lambda) - lambda * y;
		}
	} else {
		double yy = (y_cdf ? *y_cdf : y);
		for (i = 0; i < -m; i++) {
			lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
			// logll[i] = 1.0 - exp(-lambda * yy);
			logll[i] = ONE_mexp(-lambda * yy);
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_expsurv(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *x_vec, double *y_cdf, void *arg,
			  char **arg_str)
{
	return (m == 0 ? GMRFLib_SUCCESS : loglikelihood_generic_surv(thread_id, logll, x, m, idx, x_vec, y_cdf, arg, loglikelihood_exp, arg_str));
}

int loglikelihood_generic_surv(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *x_vec, double *y_cdf,
			       void *arg, GMRFLib_logl_tp *loglfun, char **arg_str)
{
#define SAFEGUARD1(value_) value_ = TRUNCATE(value_, eps, 1.0 - eps)
#define SAFEGUARD(value_)						\
	for(int i_ = 0; i_ < (m); i_++) {				\
		SAFEGUARD1(value_[i_]);					\
	}

	Data_section_tp *ds = (Data_section_tp *) arg;
	int ncov = ds->data_observations.cure_ncov;
	double pcure = 0.0, l_1mpcure = NAN;
	double eps = 0.001 * FLT_EPSILON;

	assert(y_cdf == NULL);				       // I do not think this should be used.
	if (ncov) {
		double *cov = ds->data_observations.cure_cov + idx * ncov;
		double sum = 0.0;
		for (int i = 0; i < ncov; i++) {
			sum += cov[i] * ds->data_observations.cure_beta[i][thread_id][0];
		}
		pcure = map_probability_forward(sum, MAP_FORWARD, NULL);
		SAFEGUARD1(pcure);
		l_1mpcure = LOG_p(1.0 - pcure);
	}

	double event = ds->data_observations.event[idx];
	double truncation = ds->data_observations.truncation[idx];
	double lower = ds->data_observations.lower[idx];
	double upper = ds->data_observations.upper[idx];
	int ievent = (int) event;
	int have_truncation = (truncation > 0.0);

	if (arg_str) {
		char *a = NULL, *b = NULL;
		double dummy;
		loglfun(thread_id, &dummy, x, 1, idx, x_vec, NULL, arg, &b);
		char *str_cov = Strdup(""), *str_beta = Strdup("");
		if (ncov) {
			double *cov = ds->data_observations.cure_cov + idx * ncov;
			str_cov = GMRFLib_vec2char(cov, ncov);
			double *beta = Calloc(ncov, double);
			for (int i = 0; i < ncov; i++) {
				beta[i] = ds->data_observations.cure_beta[i][thread_id][0];
			}
			str_beta = GMRFLib_vec2char(beta, ncov);
			Free(beta);
		}

		GMRFLib_sprintf(&a, "list(y.surv = list(time = %.8g, lower = %.8g, upper = %.8g, truncation = %.8g, event = %1d), \
family = \"inla.surv\", cure.prob = %.8g, cure.beta = c(%s), cure.covariates = c(%s), \
family.arg.str = %s)", ds->data_observations.y[idx], lower, upper, truncation, ievent, pcure, str_beta, str_cov, (b ? b : ""));
		*arg_str = a;
		return GMRFLib_SUCCESS;
	}

	if (m > 0) {
		Calloc_init(5 * m, 5);
		double *lf = Calloc_get(m);
		double *F_lower = Calloc_get(m);
		double *F_upper = Calloc_get(m);
		double *F_trunc = Calloc_get(m);
		double *FF_trunc = Calloc_get(m);

		GMRFLib_fill(m, 1.0, F_upper);
		if (have_truncation) {
			loglfun(thread_id, F_trunc, x, -m, idx, x_vec, &truncation, arg, arg_str);
			SAFEGUARD(F_trunc);
#pragma omp simd
			for (int i = 0; i < m; i++) {
				FF_trunc[i] = 1.0 / (1.0 - F_trunc[i]);
			}
		}

		switch (ievent) {
		case SURV_EVENT_FAILURE:
		{
			loglfun(thread_id, lf, x, m, idx, x_vec, NULL, arg, arg_str);
			if (have_truncation) {
#pragma omp simd
				for (int i = 0; i < m; i++) {
					lf[i] += LOG_p(FF_trunc[i]);
				}
			}

			if (pcure) {
#pragma omp simd
				for (int i = 0; i < m; i++) {
					logll[i] = l_1mpcure + lf[i];
				}
			} else {
				Memcpy(logll, lf, m * sizeof(double));
			}
		}
			break;

		case SURV_EVENT_RIGHT:
		{
			if (!ISZERO(lower)) {
				assert(lower >= truncation);
				loglfun(thread_id, F_lower, x, -m, idx, x_vec, &lower, arg, arg_str);
				SAFEGUARD(F_lower);
			}

			if (have_truncation) {
#pragma omp simd
				for (int i = 0; i < m; i++) {
					F_lower[i] = (F_lower[i] - F_trunc[i]) * FF_trunc[i];
				}
			}

			if (pcure) {
#pragma omp simd
				for (int i = 0; i < m; i++) {
					// logll[i] = log(pcure + (1.0 - pcure) * (1.0 - F_lower[i]));
					double A = LOG_p(pcure);
					double B = LOG_p((1.0 - pcure) * (1.0 - F_lower[i]));
					logll[i] = GMRFLib_logsum(A, B);
				}
			} else {
#pragma omp simd
				for (int i = 0; i < m; i++) {
					// logll[i] = log(1.0 - F_lower[i]);
					logll[i] = LOG_p(1.0 - F_lower[i]);
				}
			}
		}
			break;

		case SURV_EVENT_LEFT:
		{
			if (!ISINF(upper)) {
				assert(upper >= truncation);
				loglfun(thread_id, F_upper, x, -m, idx, x_vec, &upper, arg, arg_str);
				SAFEGUARD(F_upper);
			}

			if (have_truncation) {
#pragma omp simd
				for (int i = 0; i < m; i++) {
					F_upper[i] = (F_upper[i] - F_trunc[i]) * FF_trunc[i];
				}
			}

			if (pcure) {
#pragma omp simd
				for (int i = 0; i < m; i++) {
					logll[i] = LOG_p((1.0 - pcure) * F_upper[i]);
				}
			} else {
#pragma omp simd
				for (int i = 0; i < m; i++) {
					logll[i] = LOG_p(F_upper[i]);
				}
			}
		}
			break;

		case SURV_EVENT_INTERVAL:
		{
			if (!ISZERO(lower)) {
				assert(lower >= truncation);
				loglfun(thread_id, F_lower, x, -m, idx, x_vec, &lower, arg, arg_str);
				SAFEGUARD(F_lower);
			}
			if (!ISINF(upper)) {
				assert(upper >= truncation);
				loglfun(thread_id, F_upper, x, -m, idx, x_vec, &upper, arg, arg_str);
				SAFEGUARD(F_upper);
			}

			if (have_truncation) {
#pragma omp simd
				for (int i = 0; i < m; i++) {
					F_upper[i] = (F_upper[i] - F_trunc[i]) * FF_trunc[i];
					F_lower[i] = (F_lower[i] - F_trunc[i]) * FF_trunc[i];
				}
			}

			if (pcure) {
#pragma omp simd
				for (int i = 0; i < m; i++) {
					logll[i] = LOG_p((1.0 - pcure) * (F_upper[i] - F_lower[i]));
				}
			} else {
#pragma omp simd
				for (int i = 0; i < m; i++) {
					logll[i] = LOG_p(F_upper[i] - F_lower[i]);
				}
			}
		}
			break;

		case SURV_EVENT_ININTERVAL:
		{
			assert(lower >= truncation);
			if (!ISZERO(lower)) {
				assert(lower >= truncation);
				loglfun(thread_id, F_lower, x, -m, idx, x_vec, &lower, arg, arg_str);
				SAFEGUARD(F_lower);
			}
			if (!ISINF(upper)) {
				assert(upper >= truncation);
				loglfun(thread_id, F_upper, x, -m, idx, x_vec, &upper, arg, arg_str);
				SAFEGUARD(F_upper);
			}
			loglfun(thread_id, lf, x, m, idx, x_vec, NULL, arg, arg_str);

			if (have_truncation) {
#pragma omp simd
				for (int i = 0; i < m; i++) {
					F_upper[i] = (F_upper[i] - F_trunc[i]) * FF_trunc[i];
					F_lower[i] = (F_lower[i] - F_trunc[i]) * FF_trunc[i];
				}
#pragma omp simd
				for (int i = 0; i < m; i++) {
					lf[i] += LOG_p(FF_trunc[i]);
				}
			}

			if (pcure) {
#pragma omp simd
				for (int i = 0; i < m; i++) {
					logll[i] = l_1mpcure + lf[i] - LOG_p((F_upper[i] - F_lower[i]));
				}
			} else {
#pragma omp simd
				for (int i = 0; i < m; i++) {
					logll[i] = lf[i] - LOG_p((F_upper[i] - F_lower[i]));
				}
			}
		}
			break;

		default:
			GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
		}

		Calloc_free();
		return GMRFLib_SUCCESS;
	} else {
		GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
	}

#undef SAFEGUARD
#undef SAFEGUARD1
	return GMRFLib_SUCCESS;
}

int loglikelihood_weibull(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			  void *arg, char **arg_str)
{
	/*
	 * y ~ Weibull.
	 */
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	Data_section_tp *ds = (Data_section_tp *) arg;
	double y = ds->data_observations.y[idx];
	double ly = log(y);
	double alpha = map_alpha_weibull(ds->data_observations.alpha_intern[thread_id][0], MAP_FORWARD, NULL);
	double lalpha = log(alpha);

	LINK_INIT;
	if (arg_str) {
		char *a = NULL;
		GMRFLib_sprintf(&a,
				"list(y = %.8g, family = \"weibull\", theta = %.8g, alpha = %.8g, variant = %1d, link.model = \"%1s\", linear.predictor = %.8g)",
				y, ds->data_observations.alpha_intern[thread_id][0], alpha, ds->variant, ds->link_model, x[0]);
		*arg_str = a;
		return GMRFLib_SUCCESS;
	}

	switch (ds->variant) {
	case 0:
	{
		if (m > 0) {
			double ypow = pow(y, alpha);
			double t1 = lalpha + (alpha - 1.0) * ly;
			for (int i = 0; i < m; i++) {
				double lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				logll[i] = log(lambda) + t1 - lambda * ypow;
			}
		} else {
			double ypow = pow((y_cdf ? *y_cdf : y), alpha);
			for (int i = 0; i < -m; i++) {
				double lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				// logll[i] = 1.0 - exp(-lambda * ypow);
				logll[i] = ONE_mexp(-lambda * ypow);
			}
		}
	}
		break;
	case 1:
	{
		if (m > 0) {
			double t1 = lalpha - ly;
			for (int i = 0; i < m; i++) {
				double lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				double ypow = pow(lambda * y, alpha);
				logll[i] = log(ypow) + t1 - ypow;
			}
		} else {
			double yy = (y_cdf ? *y_cdf : y);
			for (int i = 0; i < -m; i++) {
				double lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				double ypow = pow(lambda * yy, alpha);
				// logll[i] = 1.0 - exp(-ypow);
				logll[i] = ONE_mexp(-ypow);
			}
		}
	}
		break;
	default:
		assert(0 == 1);
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_weibullsurv(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *x_vec, double *y_cdf,
			      void *arg, char **arg_str)
{
	return (m ==
		0 ? GMRFLib_SUCCESS : loglikelihood_generic_surv(thread_id, logll, x, m, idx, x_vec, y_cdf, arg, loglikelihood_weibull, arg_str));
}

int loglikelihood_gompertz(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			   void *arg, char **UNUSED(arg_str))
{
	/*
	 * y ~ gompertz
	 */
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}
	Data_section_tp *ds = (Data_section_tp *) arg;
	int i;
	double y, alpha, mu;

	y = ds->data_observations.y[idx];
	// yes, use the same mapping as weibull
	alpha = map_alpha_gompertz(ds->data_observations.alpha_intern[thread_id][0], MAP_FORWARD, NULL);

	LINK_INIT;
	if (m > 0) {
		for (i = 0; i < m; i++) {
			mu = PREDICTOR_INVERSE_LINK(x[i] + off);
			logll[i] = log(mu) + alpha * y - mu * (exp(alpha * y) - 1.0) / alpha;
			// if (i == 0)printf("idx %d x %f mu %f logll %f y %f alpha %f\n", idx, x[i], mu, logll[i], y, alpha);
		}
	} else {
		double yy = (y_cdf ? *y_cdf : y);
		for (i = 0; i < -m; i++) {
			mu = PREDICTOR_INVERSE_LINK(x[i] + off);
			// logll[i] = 1.0 - exp(-mu * (exp(alpha * yy) - 1.0) / alpha);
			logll[i] = ONE_mexp(-mu * (exp(alpha * yy) - 1.0) / alpha);
		}
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_gompertzsurv(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *x_vec, double *y_cdf,
			       void *arg, char **arg_str)
{
	return (m ==
		0 ? GMRFLib_SUCCESS : loglikelihood_generic_surv(thread_id, logll, x, m, idx, x_vec, y_cdf, arg, loglikelihood_gompertz, arg_str));
}

int loglikelihood_loglogistic(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			      void *arg, char **UNUSED(arg_str))
{
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	Data_section_tp *ds = (Data_section_tp *) arg;
	int i;
	double y, ly = 0, lambda, alpha, lalpha = 0;

	y = ds->data_observations.y[idx];
	alpha = map_exp_forward(ds->data_observations.alpha_intern[thread_id][0], MAP_FORWARD, NULL);
	LINK_INIT;

	if (m > 0) {
		ly = log(y);
		lalpha = log(alpha);
	}

	switch (ds->variant) {
	case 0:
	{
		if (m > 0) {
			for (i = 0; i < m; i++) {
				lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				logll[i] = log(lambda) + (-alpha - 1.0) * ly + lalpha - 2.0 * log1p(lambda * pow(y, -alpha));
			}
		} else {
			double yy = (y_cdf ? *y_cdf : y);
			for (i = 0; i < -m; i++) {
				lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				logll[i] = 1.0 / (1.0 + lambda * pow(yy, -alpha));
			}
		}
	}
		break;

	case 1:
	{
		if (m > 0) {
			double lam_y;
			for (i = 0; i < m; i++) {
				lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				lam_y = lambda * y;
				logll[i] = -alpha * log(lam_y) + lalpha - ly - 2.0 * log1p(pow(lam_y, -alpha));
			}
		} else {
			double yy = (y_cdf ? *y_cdf : y);
			for (i = 0; i < -m; i++) {
				lambda = PREDICTOR_INVERSE_LINK(x[i] + off);
				logll[i] = 1.0 / (1.0 + pow(lambda * yy, -alpha));
			}
		}
	}
		break;

	default:
		assert(0 == 1);
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_loglogisticsurv(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *x_vec, double *y_cdf,
				  void *arg, char **arg_str)
{
	return (m ==
		0 ? GMRFLib_SUCCESS : loglikelihood_generic_surv(thread_id, logll, x, m, idx, x_vec, y_cdf, arg, loglikelihood_loglogistic,
								 arg_str));
}

int loglikelihood_qloglogistic(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
			       void *arg, char **UNUSED(arg_str))
{
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double y, yq, ly = NAN, lambda, alpha, lalpha = NAN, q, qq;

	y = ds->data_observations.y[idx];
	alpha = map_exp_forward(ds->data_observations.alpha_intern[thread_id][0], MAP_FORWARD, NULL);
	q = ds->data_observations.quantile;
	qq = 1.0 / q - 1.0;
	LINK_INIT;

	if (m > 0) {
		ly = log(y);
		lalpha = log(alpha);
	}

	switch (ds->variant) {
	case 0:
	{
		if (m > 0) {
			for (i = 0; i < m; i++) {
				yq = PREDICTOR_INVERSE_LINK(x[i] + off);
				lambda = qq * pow(yq, alpha);
				logll[i] = log(lambda) + (-alpha - 1.0) * ly + lalpha - 2.0 * log1p(lambda * pow(y, -alpha));
			}
		} else {
			double yy = (y_cdf ? *y_cdf : y);
			for (i = 0; i < -m; i++) {
				yq = PREDICTOR_INVERSE_LINK(x[i] + off);
				lambda = qq * pow(yq, alpha);
				logll[i] = 1.0 / (1.0 + lambda * pow(yy, -alpha));
			}
		}
	}
		break;

	case 1:
	{
		if (m > 0) {
			double lam_y, qqinv = 1.0 / qq;
			for (i = 0; i < m; i++) {
				yq = PREDICTOR_INVERSE_LINK(x[i] + off);
				lambda = 1.0 / yq * pow(qqinv, 1.0 / alpha);
				lam_y = lambda * y;
				logll[i] = -alpha * log(lam_y) + lalpha - ly - 2.0 * log1p(pow(lam_y, -alpha));
			}
		} else {
			double yy = (y_cdf ? *y_cdf : y), qqinv = 1.0 / qq;
			for (i = 0; i < -m; i++) {
				yq = PREDICTOR_INVERSE_LINK(x[i] + off);
				lambda = 1.0 / yq * pow(qqinv, 1.0 / alpha);
				logll[i] = 1.0 / (1.0 + pow(lambda * yy, -alpha));
			}
		}
	}
		break;

	default:
		assert(0 == 1);
	}

	LINK_END;
	return GMRFLib_SUCCESS;
}

int loglikelihood_qloglogisticsurv(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *x_vec, double *y_cdf,
				   void *arg, char **arg_str)
{
	return (m ==
		0 ? GMRFLib_SUCCESS : loglikelihood_generic_surv(thread_id, logll, x, m, idx, x_vec, y_cdf, arg, loglikelihood_qloglogistic,
								 arg_str));
}

int loglikelihood_fmrisurv(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *x_vec, double *y_cdf, void *arg,
			   char **arg_str)
{

	return (m == 0 ? GMRFLib_SUCCESS : loglikelihood_generic_surv(thread_id, logll, x, m, idx, x_vec, y_cdf, arg, loglikelihood_fmri, arg_str));
}

int loglikelihood_fmri(int thread_id, double *__restrict logll, double *__restrict x, int m, int idx, double *UNUSED(x_vec), double *y_cdf,
		       void *arg, char **UNUSED(arg_str))
{
	/*
	 * y ~ fmri (noncentral-chi distribution).
	 */
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	int i;
	Data_section_tp *ds = (Data_section_tp *) arg;
	double eta, y, yy, prec, dof, ncp, scale, l2y, y2, yy2;

	y = ds->data_observations.y[idx];
	prec = map_exp_forward(ds->data_observations.fmri_lprec[thread_id][0], MAP_FORWARD, NULL);
	dof = map_identity_forward(ds->data_observations.fmri_ldof[thread_id][0], MAP_FORWARD, NULL);
	scale = (ds->data_observations.fmri_scale ? ds->data_observations.fmri_scale[idx] : 1.0);

	prec *= scale;
	LINK_INIT;

	if (m > 0) {
		y2 = prec * SQR(y);
		l2y = log(2.0 * prec * y);
		for (i = 0; i < m; i++) {
			eta = PREDICTOR_INVERSE_LINK(x[i] + off);
			ncp = prec * SQR(eta);
			// more robust implementation provided by L.Starke
			// logll[i] = l2y + MATHLIB_FUN(dnchisq) (y2, dof, ncp, 1);
			logll[i] = l2y + inla_dnchisq(y2, dof, ncp);
		}
	} else {
		yy = (y_cdf ? *y_cdf : y);
		yy2 = prec * SQR(yy);
		for (i = 0; i < -m; i++) {
			eta = PREDICTOR_INVERSE_LINK(x[i] + off);
			ncp = prec * SQR(eta);
			logll[i] = MATHLIB_FUN(pnchisq) (yy2, dof, ncp, 1, 0);
		}

	}

	LINK_END;
	return GMRFLib_SUCCESS;
}
