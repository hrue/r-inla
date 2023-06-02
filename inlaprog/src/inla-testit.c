
/* inla-testit.c
 * 
 * Copyright (C) 2007-2023 Havard Rue
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


int loglikelihood_testit(int UNUSED(thread_id), double *logll, double *x, int m, int UNUSED(idx), double *x_vec, double *UNUSED(y_cdf),
			 void *UNUSED(arg), char **UNUSED(arg_str))
{
	if (m == 0) {
		return GMRFLib_LOGL_COMPUTE_CDF;
	}

	int i;
	double a = 1.0, b = 2.0, c = 3.0, d = 4.0, x0;
	x0 = x_vec[0];

	if (m > 0) {
		for (i = 0; i < m; i++) {
			double xx = x[i] - x0;
			logll[i] = a + b * xx - c / 2.0 * SQR(xx) + d / 6.0 * gsl_pow_3(xx);
		}
	} else {
		abort();
	}
	return GMRFLib_SUCCESS;
}

int inla_testit_timer(void)
{
	GMRFLib_ENTER_ROUTINE;
	int ret = system("sleep 1");
	if (ret != 0)
		exit(1);
	GMRFLib_LEAVE_ROUTINE;
	return 0;
}

double testit_Qfunc(int UNUSED(thread_id), int i, int j, double *UNUSED(values), void *UNUSED(arg))
{
	return (i == j ? 100.0 : -1.0);
}

int testit(int argc, char **argv)
{
	int test_no = -1;
	char **args = NULL;
	int nargs = 0;

	if (argc > 0) {
		test_no = atoi(argv[0]);
		nargs = argc - 1;
		args = &(argv[1]);
	}
	printf("test_no = %1d  nargs = %1d\n", test_no, nargs);
	for (int i = 0; i < nargs; i++) {
		printf("\targs[%d] = %s\n", i, args[i]);
	}

	switch (test_no) {
	case -1:
	case 0:
	{
		double s, phi, mu, yyy, shape, rate, q, mmu, scale, alpha;
		s = 1.2;
		phi = 2.3;
		mu = 3.4;
		yyy = 4.5;
		shape = s * phi;
		rate = s * phi / mu;
		scale = 1.0 / rate;
		alpha = 0.5;
		alpha = MATHLIB_FUN(pgamma) (yyy, shape, scale, 1, 0);
		q = MATHLIB_FUN(qgamma) (alpha, shape, 1.0, 1, 0);
		mmu = yyy * s * phi / q;
		printf("alpha %f q %f mu %f mmu %f diff %f\n", alpha, q, mu, mmu, mu - mmu);

		P(MATHLIB_FUN(pgamma) (yyy, shape, scale, 1, 0));
		P(MATHLIB_FUN(qgamma) (alpha, shape, 1.0, 1, 0));
		P(gsl_cdf_gamma_P(yyy, shape, scale));
		P(gsl_cdf_gamma_Pinv(alpha, shape, 1.0));

		if (nargs) {
			int n = atoi(args[0]);
			double *x = Calloc(n, double);
			double *y = Calloc(n, double);
			double *yy = Calloc(n, double);

			for (int i = 0; i < n; i++) {
				x[i] = GMRFLib_uniform();
			}

			double tref[] = { 0, 0 };
			tref[0] -= GMRFLib_cpu();
			for (int i = 0; i < n; i++) {
				y[i] = MATHLIB_FUN(qgamma) (x[i], exp(x[i]), 1.0, 1, 0);
			}
			tref[0] += GMRFLib_cpu();

			tref[1] -= GMRFLib_cpu();
			for (int i = 0; i < n; i++) {
				yy[i] = gsl_cdf_gamma_Pinv(x[i], exp(x[i]), 1.0);
			}
			tref[1] += GMRFLib_cpu();

			P(y[0] - yy[0]);
			printf("MATHLIB %f GSL %f\n", tref[0] / (tref[0] + tref[1]), tref[1] / (tref[0] + tref[1]));

			Free(x);
			Free(y);
			Free(yy);
		}
	}
		break;

	case 1:
	{
		for (int i = 0; i < 10; i++) {
			P(GMRFLib_rng_uniform());
		}
	}
		break;

	case 2:
	{
		double par[] = { 0.8, 0.5 };
		double theta = 1.234;

		P(priorfunc_pc_cor1(&theta, par));
		exit(0);
	}
		break;

	case 3:
	{
		double x;
		double lambda = 1.2345;
		P(lambda);
		for (x = 0.1; x < 10; x += 0.1) {
			printf("F(%g) = %.12g\n", x, gsl_sf_gamma_inc_Q(x, lambda));
			printf("F(%g) = %.12g\n", x, gsl_sf_gamma_inc(x, lambda) / gsl_sf_gamma(x));
			printf("F(%g) = %.12g\n", x, gsl_sf_gamma_inc(x, lambda) / exp(gsl_sf_lngamma(x)));
		}
		exit(0);
	}
		break;

	case 4:
	{
		double lambda = 1.234;
		double x;
		for (x = -5; x < 5; x += 0.1) {
			printf("%f %f\n", x, priorfunc_pc_gammacount(&x, &lambda));
		}
		exit(0);
	}
		break;

	case 5:
	{
		// this force a race-condition
#define NN 10
		int x[NN];
#pragma omp parallel for
		for (int i = 0; i < NN; i++) {
			*(x + i) = i;
		}
		P((double) x[0]);
		P((double) x[NN - 1]);
#undef NN
	}
		break;

	case 6:
	{
		double y, lambda;
		for (lambda = 1.1; lambda < 5.1; lambda++) {
			for (y = 2.2; y < 8.3; y++) {
				printf("y %f lambda %f cdf cdf.deriv = %f %f\n",
				       y, lambda, inla_pcontpois(y, lambda), inla_pcontpois_deriv(y, lambda));
			}
		}

		double quantile, alpha;
		for (quantile = 1.1; quantile < 10; quantile++) {
			for (alpha = 0.1; alpha < 0.99; alpha += 0.1) {
				printf("quantile=%f alpha=%f eta=%f\n", quantile, alpha, inla_qcontpois_eta(quantile, alpha, NULL));
			}
		}

		exit(0);
	}
		break;

	case 7:
	{

		inla_fgn_arg_tp *arg = Calloc(1, inla_fgn_arg_tp);
		arg->n = 10;
		arg->k = 3;
		arg->N = arg->n * (arg->k + 1);
		arg->prec_eps = 100;

		GMRFLib_graph_tp *g;
		inla_make_fgn_graph(&g, arg);
		GMRFLib_printf_graph(stdout, g);
		exit(0);
	}
		break;

	case 8:
	{
#pragma omp critical (Name_5f0e63e7e6c4127b4e2e50c9002880440da56212)
		{
#define _MODEL "rgeneric.model"
			printf("test rgeneric\n");
			inla_R_library("INLA");
			inla_R_load("rgeneric.RData");
			inla_R_source("/home/hrue/p/inla/r-inla/rinla/R/rgeneric.R");

			double theta[] = { 1.0, 2.0 };
			int ntheta = sizeof(theta) / sizeof(double);
			int i;

			int n_out;
			double *x_out;

#define _PPP(cmd)							\
			printf("\ncmd [%s] n_out [%1d]\n", cmd, n_out);	\
			for(i=0; i<n_out; i++) {			\
				printf("x[ %1d ] = %g\n", i, x_out[i]);	\
			};						\
			Free(x_out)

			inla_R_rgeneric(&n_out, &x_out, "graph", _MODEL, &ntheta, theta);
			_PPP("graph");

			inla_R_rgeneric(&n_out, &x_out, "Q", _MODEL, &ntheta, theta);
			_PPP("Q");

			inla_R_rgeneric(&n_out, &x_out, "initial", _MODEL, &ntheta, theta);
			_PPP("initial");

			inla_R_rgeneric(&n_out, &x_out, "log.norm.const", _MODEL, &ntheta, theta);
			_PPP("log.norm.const");

			inla_R_rgeneric(&n_out, &x_out, "log.prior", _MODEL, &ntheta, theta);
			_PPP("log.prior");

			inla_R_rgeneric(&n_out, &x_out, "quit", _MODEL, &ntheta, theta);
			_PPP("quit");
		}

#undef _MODEL
#undef _PPP
	}
		break;

	case 9:
	{
		printf("test R, source %s\n", argv[0]);
		inla_R_source(argv[0]);

		double x[] = { 1.123, 2.234, 3.345 };
		int nx = sizeof(x) / sizeof(x[1]);

		double *xx = NULL;
		int nxx;
		int i;

		inla_R_assign("x", &nx, x);
		inla_R_get(&nxx, &xx, "x");
		for (i = 0; i < nxx; i++) {
			printf("xx[%1d] = %f\n", i, xx[i]);
		}

		exit(0);
	}
		break;

	case 10:
	{
		// checking the expression and the jacobian for this prior
		double x, xx, xxx, dx = 0.01, sum = 0.0, parameters[2], low = -4.001, high = 4.0;
		int i;

		parameters[0] = 2.123;			       /* lambda */
		parameters[1] = 1;			       /* p */
		sum = 0;
#pragma omp parallel for private(i, x, xx) reduction(+: sum)
		for (i = 0; i < (int) ((high - low) / dx + 1); i++) {
			x = low + dx * i;
			double x2[1];
			x2[0] = x;
			sum += exp(priorfunc_pc_ar(x2, parameters));
		}
		P(sum * pow(dx, 1.0));

		parameters[0] = 2.123;			       /* lambda */
		parameters[1] = 2;			       /* p */
		sum = 0;
#pragma omp parallel for private(i, x, xx) reduction(+: sum)
		for (i = 0; i < (int) ((high - low) / dx + 1); i++) {
			x = low + dx * i;
			for (xx = low; xx < high; xx += dx) {
				double x2[2];
				x2[0] = x;
				x2[1] = xx;
				sum += exp(priorfunc_pc_ar(x2, parameters));
			}
		}
		P(sum * pow(dx, 2.0));

		parameters[0] = 3.123;			       /* lambda */
		parameters[1] = 3;			       /* p */
		sum = 0;
#pragma omp parallel for private(i, x, xx, xxx) reduction(+: sum)
		for (i = 0; i < (int) ((high - low) / dx + 1); i++) {
			x = low + dx * i;
			for (xx = low; xx < high; xx += dx) {
				for (xxx = low; xxx < high; xxx += dx) {
					double x2[3];
					x2[0] = x;
					x2[1] = xx;
					x2[2] = xxx;
					sum += exp(priorfunc_pc_ar(x2, parameters));
				}
			}
		}
		P(sum * pow(dx, 3.0));
	}
		break;

	case 11:
	{
		// test the new R-interface

		printf("TESTIT!\n");
		inla_R_source("example-code.R");
		double x[] = { 1, 2, 3 };
		int nx = sizeof(x) / sizeof(x[1]);

		double *xx = NULL;
		int nxx;
#pragma omp parallel for
		for (int i = 0; i < 10; i++) {
			inla_R_funcall2(&nxx, &xx, "lprior2", "ThisIsTheTag", &nx, x);
		}
		inla_R_funcall2(&nxx, &xx, "lprior2", NULL, &nx, x);

		for (int i = 0; i < nxx; i++) {
			printf("lprior2[%1d] = %f\n", i, xx[i]);
		}

		exit(0);
	}
		break;

	case 12:
	{
		// testing spde3

		inla_spde3_tp *smodel = NULL;
		int i, j, jj, n, t, p;
		int thread_id = 0;

		inla_spde3_build_model(thread_id, &smodel, "./", "identity");
		n = smodel->graph->n;
		p = smodel->ntheta;

		GMRFLib_matrix_tp *theta = GMRFLib_read_fmesher_file("theta", 0, -1);

		for (i = 0; i < p; i++) {
			printf("theta[%1d] = %g\n", i, GMRFLib_matrix_get(i, 0, theta));
			for (t = 0; t < GMRFLib_MAX_THREADS(); t++) {
				smodel->theta[i][t][0] = GMRFLib_matrix_get(i, 0, theta);
			}
		}

		GMRFLib_matrix_tp *Q = Calloc(1, GMRFLib_matrix_tp);
		Q->nrow = Q->ncol = n;
		Q->elems = ISQR(n);
		Q->A = Calloc(ISQR(n), double);

		int ntimes = 10, itim;

		for (itim = 0; itim < ntimes; itim++) {
#pragma omp parallel for private(i, jj, j, thread_id)
			for (i = 0; i < n; i++) {
				thread_id = omp_get_thread_num();
				Q->A[i + i * n] = inla_spde3_Qfunction(thread_id, i, i, NULL, (void *) smodel);
				for (jj = 0; jj < smodel->graph->nnbs[i]; jj++) {
					j = smodel->graph->nbs[i][jj];
					Q->A[i + j * n] = Q->A[j + i * n] = inla_spde3_Qfunction(thread_id, i, j, NULL, (void *) smodel);
				}
			}
		}
		GMRFLib_write_fmesher_file(Q, "Q", 0, -1);
	}
		break;

	case 13:
	{
		int n = 100;
		if (nargs) {
			n = atoi(args[0]);
		}

		printf("Build matrix with dim = %1d\n", n);
		double *A = Calloc(SQR(n), double);
		for (int i = 0; i < n; i++) {
			for (int j = i + 1; j < n; j++) {
				int k = i + j * n;
				int kk = j + i * n;
				A[k] = A[kk] = GMRFLib_uniform();
			}
		}
		for (int i = 0; i < n; i++) {
			int k = i + i * n;
			A[k] = n + 1.0;
		}

		printf("Call ...ensure_spd ");
		double tref = GMRFLib_cpu();
		GMRFLib_ensure_spd(A, n, FLT_EPSILON, NULL);
		printf("%f seconds\n", GMRFLib_cpu() - tref);

		printf("Call ...chol ");
		tref = GMRFLib_cpu();
		double *chol = NULL;
		GMRFLib_comp_chol_general(&chol, A, n, NULL, 1);
		printf("%f seconds\n", GMRFLib_cpu() - tref);

		Free(A);
		Free(chol);
	}
		break;

	case 14:
	{
		double a_data[] = { 0.18, 0.41, 0.14, 0.51,
			0.60, 0.24, 0.30, 0.13,
			0.57, 0.99, 0.97, 0.19,
			0.96, 0.58, 0.66, 0.85
		};

		gsl_matrix_view m = gsl_matrix_view_array(a_data, 4, 4);
		gsl_matrix *A = GMRFLib_gsl_duplicate_matrix(&m.matrix);

		GMRFLib_printf_gsl_matrix(stdout, A, " %.12f");
		printf("\n");
		GMRFLib_gsl_ginv(A, GSL_SQRT_DBL_EPSILON, -1);
		GMRFLib_printf_gsl_matrix(stdout, A, " %.12f");

		exit(EXIT_SUCCESS);
	}
		break;

	case 15:
	{
		P(sizeof(void *));
		P(sizeof(double));
		P(sizeof(double *));
		P(sizeof(int *));
	}
		break;

	case 16:
	{
		ar_test1();
		exit(EXIT_SUCCESS);
	}
		break;

	case 17:
		break;

	case 18:
	{
		double x[] = { -2, -1, 0, 1, 2 };
		double ld[] = { -2.0, -0.5, 0.0, -0.5, -2.0 };
		double m = 0.0, sd = 0.0;
		GMRFLib_vb_fit_gaussian(5, x, ld, &m, &sd);
		P(m);
		P(sd);
	}
		break;

	case 19:
	{
		printf("physical= %1d logical= %1d\n", UTIL_countPhysicalCores(), UTIL_countLogicalCores());
	}
		break;

	case 20:
	{
		inla_file_contents_tp *fc;

		fc = inla_read_file_contents("aa.dat");
		inla_write_file_contents("bb.dat", fc);
		exit(EXIT_SUCCESS);
	}
		break;

	case 21:
	{
		GMRFLib_matrix_tp *M = NULL;

		int i, j, k, kk;

		printf("read file %s\n", argv[3]);
		M = GMRFLib_read_fmesher_file(argv[3], 0L, -1);

		if (1)
			if (M->i)
				for (k = 0; k < M->elems; k++)
					printf("k %d %d %d %g\n", k, M->i[k], M->j[k], M->values[k]);

		if (M->graph) {
			printf("n %d\n", M->graph->n);
			for (k = 0; k < M->graph->n; k++) {
				printf("%d nnbs %d:\n", k, M->graph->nnbs[k]);
				for (kk = 0; kk < M->graph->nnbs[k]; kk++)
					printf("\t\t%d\n", M->graph->nbs[k][kk]);
			}
		}

		for (i = 0; i < M->nrow; i++)
			for (j = 0; j < M->ncol; j++)
				printf("%d %d %g\n", i, j, GMRFLib_matrix_get(i, j, M));

		printf("\n\ntranspose...\n\n\n");
		GMRFLib_matrix_tp *N = GMRFLib_matrix_transpose(M);

		if (1)
			if (N->i)
				for (k = 0; k < N->elems; k++)
					printf("k %d %d %d %g\n", k, N->i[k], N->j[k], N->values[k]);

		if (1)
			for (i = 0; i < N->nrow; i++)
				for (j = 0; j < N->ncol; j++)
					printf("%d %d %g\n", i, j, GMRFLib_matrix_get(i, j, N));

		if (N->graph) {
			printf("n %d\n", N->graph->n);
			for (k = 0; k < N->graph->n; k++) {
				printf("%d nnbs %d:\n", k, N->graph->nnbs[k]);
				for (kk = 0; kk < N->graph->nnbs[k]; kk++)
					printf("\t\t%d\n", N->graph->nbs[k][kk]);
			}
		}
		GMRFLib_matrix_free(M);
		GMRFLib_matrix_free(N);
	}
		break;

	case 22:
	{
		double x;
		for (x = -100.0; x < 100.0; x = x + 0.1) {
			printf("x %.12g log(Phi(x)) %.12g %.12g\n", x, inla_log_Phi_fast(x), inla_log_Phi(x));
		}
		exit(0);
	}
		break;

	case 23:
	{
		// test ghq
#define FUN0(x) (1)
#define FUN1(x) (x)
#define FUN2(x) SQR(x)
#define FUN3(x) (SQR(x)*(x))
#define FUN4(x) SQR(SQR(x))

		double *xp, *wp, integral0 = 0, integral1 = 0, integral2 = 0, integral3 = 0, integral4 = 0.0;
		int np = GMRFLib_INT_GHQ_POINTS, i;
		if (nargs) {
			np = atoi(args[0]);
		}
		GMRFLib_ghq(&xp, &wp, np);
		for (i = 0; i < np; i++) {
			integral0 += wp[i] * FUN0(xp[i]);
			integral1 += wp[i] * FUN1(xp[i]);
			integral2 += wp[i] * FUN2(xp[i]);
			integral3 += wp[i] * FUN3(xp[i]);
			integral4 += wp[i] * FUN4(xp[i]);
			printf("x[%1d]= %.12f  w[%1d] = %.12f\n", i, xp[i], i, wp[i]);
		}
		printf("integral of x^0 = 1 ?  %.12f\n", integral0);
		printf("integral of x^1 = 0 ?  %.12f\n", integral1);
		printf("integral of x^2 = 1 ?  %.12f\n", integral2);
		printf("integral of x^3 = 0 ?  %.12f\n", integral3);
		printf("integral of x^4 = 3 ?  %.12f\n", integral4);

		double mean = GMRFLib_uniform();
		double stdev = exp(GMRFLib_uniform());
		GMRFLib_ghq_ms(&xp, &wp, np, mean, stdev);
		integral0 = 0;
		integral1 = 0;
		integral2 = 0;
		integral3 = 0;
		integral4 = 0;
		for (i = 0; i < np; i++) {
			double xx = (xp[i] - mean) / stdev;
			integral0 += wp[i] * FUN0(xx);
			integral1 += wp[i] * FUN1(xx);
			integral2 += wp[i] * FUN2(xx);
			integral3 += wp[i] * FUN3(xx);
			integral4 += wp[i] * FUN4(xx);
			printf("x[%1d]= %.12f  w[%1d] = %.12f\n", i, xp[i], i, wp[i]);
		}
		printf("integral of z^0 = %.12f ?  %.12f\n", 1.0, integral0);
		printf("integral of z^1 = %.12f ?  %.12f\n", 0.0, integral1);
		printf("integral of z^2 = %.12f ?  %.12f\n", 1.0, integral2);
		printf("integral of z^3 = %.12f ?  %.12f\n", 0.0, integral3);
		printf("integral of z^4 = %.12f ?  %.12f\n", 3.0, integral4);

		double xx = 2 * (GMRFLib_uniform() - 0.5);

		printf("compute CDF xx=%f true = %.12f\n", xx, inla_Phi(xx));

		integral0 = 0.0;
		for (i = 0; i < np; i++) {
			if (xp[i] < xx) {
				integral0 += wp[i];
			} else {
				integral0 += wp[i] * (1.0 - (xx - xp[i - 1]) / (xp[i] - xp[i - 1]));
				break;
			}
		}
		printf("estimate %.12f\n", integral0);

		exit(0);
#undef FUN2
#undef FUN4
	}
		break;

	case 24:
	{
		int n = 1024 * 1024 * 32;
		double *xx = Calloc(n, double);
		double *yy = Calloc(n, double);

		for (int i = 0; i < n; i++) {
			xx[i] = GMRFLib_uniform();
			yy[i] = GMRFLib_uniform();
		}

		int one = 1;
		double sum1 = 0.0, sum2 = 0.0;
		double tref1 = 0.0, tref2 = 0.0;
		for (int k = 0; k < 10; k++) {
			sum1 = sum2 = 0.0;
			tref1 -= GMRFLib_cpu();
#pragma GCC ivdep
			for (int i = 0; i < n; i++) {
				sum1 += xx[i] * yy[i];
			}
			tref1 += GMRFLib_cpu();
			tref2 -= GMRFLib_cpu();
			sum2 += ddot_(&n, xx, &one, yy, &one);
			tref2 += GMRFLib_cpu();
			if (k == 0)
				P(sum1 - sum2);
		}
		printf("loop %.3f ddot %.3f (%.3f, %.3f)\n", tref1, tref2, tref1 / (tref1 + tref2), tref2 / (tref1 + tref2));
		Free(xx);
		Free(yy);
	}
		break;

	case 25:
	{
		int n = 256 * 512;
		double *xx = Calloc(n, double);

		for (int i = 0; i < n; i++) {
			xx[i] = GMRFLib_uniform();
		}

		GMRFLib_idxval_tp *h = NULL;
		for (int i = 0, j = 0; i < n; i++) {
			j += 1 + (GMRFLib_uniform() < 0.9 ? 0.0 : (int) (GMRFLib_uniform() * 8));
			if (j >= n)
				break;
			GMRFLib_idxval_add(&h, j, xx[j]);
		}
		GMRFLib_idxval_sort(h);
		assert(h);
		P(h->g_n);
		P(h->n / h->g_n);
		double sum1 = 0.0, sum2 = 0.0;
		double tref1 = 0.0, tref2 = 0.0;
		for (int k = 0; k < 10000; k++) {
			sum1 = sum2 = 0.0;
			tref1 -= GMRFLib_cpu();
			sum1 = GMRFLib_dot_product_serial(h, xx);
			tref1 += GMRFLib_cpu();
			tref2 -= GMRFLib_cpu();
			sum2 = GMRFLib_dot_product_group(h, xx);
			tref2 += GMRFLib_cpu();
			if (ABS(sum1 - sum2) > 1e-8) {
				P(sum1);
				P(sum2);
				exit(88);
			}
		}
		printf("serial %.3f group %.3f (%.3f, %.3f)\n", tref1, tref2, tref1 / (tref1 + tref2), tref2 / (tref1 + tref2));
		Free(xx);
	}
		break;

	case 26:
	{
		int nrow = 10;
		GMRFLib_vmatrix_tp *m = NULL;
		GMRFLib_vmatrix_init(&m, nrow, NULL);

		for (int i = 0; i < nrow; i++) {
			for (int j = i; j < nrow; j++) {
				double *val = Calloc(1, double);
				*val = (double) j;
				GMRFLib_vmatrix_set(m, i, j, val);
			}
		}
		for (int i = 0; i < nrow; i++) {
			for (int j = i; j < nrow; j++) {
				double *val = GMRFLib_vmatrix_get(m, i, j);
				printf("i %d j %d val %g\n", i, j, *val);
			}
		}
		GMRFLib_vmatrix_free(m, 1);
	}
		break;

	case 27:
	{
		double eps = atof(args[0]);
		for (double val = 0.99; val < 1.0; val += eps / 1000.0) {
			int eq = ISEQUAL_x(val, 1.0, eps);
			if (eq) {
				printf("eps %.12f val %.12f %d\n", eps, val, ISEQUAL_x(val, 1.0, eps));
				break;
			}
		}
	}
		break;

	case 28:
	{
		double range = 1.9;
		double nu = 0.94;
		double kappa = sqrt(8.0 * nu) / range;
		double d, dd;
		double corf;

		for (d = 0.0; d < 3.0 * range; d += range / 10.0) {
			dd = kappa * d;
			corf = (dd <= 0.0 ? 1.0 : 1.0 / pow(2.0, nu - 1.0) / MATHLIB_FUN(gammafn) (nu) *
				pow(dd, nu) * MATHLIB_FUN(bessel_k) (dd, nu, 1.0));
			printf("dmatern nu %.3f range %.3f dist %.3f dd %.5f corf %.5f\n", nu, range, d, dd, corf);
		}
	}
		break;

	case 29:
	{
		GMRFLib_idx_tp *h = NULL;
		int i;

		for (i = 0; i < 10; i++)
			GMRFLib_idx_add(&h, i);
		GMRFLib_idx_prune(h);
		GMRFLib_idx_printf(stdout, h, "IDX-test");
		GMRFLib_idx_free(h);

		GMRFLib_idx2_tp *h2 = NULL;
		GMRFLib_idx2_create(&h2);
		for (i = 0; i < 10; i++)
			GMRFLib_idx2_add(&h2, i, -i);
		GMRFLib_idx2_prune(h2);
		GMRFLib_idx2_printf(stdout, h2, "IDX2-test");
		GMRFLib_idx2_free(h2);

		GMRFLib_val_tp *hh = NULL;
		for (i = 0; i < 10; i++)
			GMRFLib_val_add(&hh, (double) i);
		GMRFLib_val_prune(hh);
		GMRFLib_val_printf(stdout, hh, "VAL-test");
		GMRFLib_val_free(hh);

		GMRFLib_idxval_tp *h3 = NULL;
		for (i = 0; i < 10; i++)
			GMRFLib_idxval_add(&h3, i, (double) i);
		GMRFLib_idxval_prune(h3);
		GMRFLib_idxval_printf(stdout, h3, "VAL-test");
		GMRFLib_idxval_free(h3);
	}
		break;

	case 30:
	{
		double ta[] = {
			2.8457954755, -0.4965452301, -1.645446141, -1.128319792, 1.262602638,
			-0.4965452301, 9.5273534221, 6.998890429, 4.924730852, -2.979131713,
			-1.6454461413, 6.9988904291, 6.648388858, 3.208607495, -2.034296532,
			-1.1283197920, 4.9247308523, 3.208607495, 3.515410801, -2.453892405,
			1.2626026384, -2.9791317133, -2.034296532, -2.453892405, 1.833029623
		};

		gsl_matrix_view m = gsl_matrix_view_array(ta, 5, 5);
		gsl_matrix *A = GMRFLib_gsl_duplicate_matrix(&m.matrix);

		GMRFLib_printf_gsl_matrix(stdout, A, " %.12f");
		printf("\n");
		gsl_matrix *B = GMRFLib_gsl_low_rank(A, 1.0E-8);
		GMRFLib_printf_gsl_matrix(stdout, B, " %.12f");
		gsl_matrix_free(A);
		gsl_matrix_free(B);
	}
		break;

	case 31:
	{
		double xx, df;
		Link_param_tp *param = Calloc(1, Link_param_tp);
		int thread_id = 0;
		assert(omp_get_thread_num() == 0);

		for (df = 4.0; df <= 125; df += 60.0) {
			printf("df = %.4g\n", df);
			for (xx = 0.01; xx <= 0.99; xx += 0.1) {
				double cdf = MATHLIB_FUN(pt) (xx, df, 1, 0);
				double icdf = MATHLIB_FUN(qt) (cdf, df, 1, 0);
				double ldens = MATHLIB_FUN(dt) (xx, df, 1);
				printf("\txx= %.6f  cdf= %.6f icdf= %.6f ldens = %.6f\n", xx, cdf, icdf, ldens);
			}
		}
		for (df = 4.0; df <= 125; df += 60.0) {
			double scale = sqrt(df / (df - 2.0));
			printf("Normalized df = %.4g\n", df);
			for (xx = 0.01; xx <= 0.99; xx += 0.1) {
				double cdf = MATHLIB_FUN(pt) (xx / scale, df, 1, 0);
				double icdf = MATHLIB_FUN(qt) (cdf, df, 1, 0) * scale;
				double ldens = MATHLIB_FUN(dt) (xx / scale, df, 1) - log(scale);
				printf("\txx= %.6f  cdf= %.6f icdf= %.6f ldens = %.6f\n", xx, cdf, icdf, ldens);
			}
		}

		for (df = 4.0; df <= 125; df += 60.0) {
			printf("Link probit\n");
			HYPER_NEW(param->dof_intern, log(df - 2.0));

			for (xx = 0.01; xx <= 0.99; xx += 0.1) {
				double back, forw, dforw;
				back = link_probit(thread_id, xx, MAP_BACKWARD, (void *) param, NULL);
				forw = link_probit(thread_id, back, MAP_FORWARD, (void *) param, NULL);
				dforw = link_probit(thread_id, back, MAP_DFORWARD, (void *) param, NULL);
				printf("\txx= %.6f  back= %.6f forw= %.6f dforw = %.6f\n", xx, back, forw, dforw);
			}

			printf("Link df = %.4g\n", df);
			for (xx = 0.01; xx <= 0.99; xx += 0.1) {
				double back, forw, dforw;
				back = link_robit(thread_id, xx, MAP_BACKWARD, (void *) param, NULL);
				forw = link_robit(thread_id, back, MAP_FORWARD, (void *) param, NULL);
				dforw = link_robit(thread_id, back, MAP_DFORWARD, (void *) param, NULL);
				printf("\txx= %.6f  back= %.6f forw= %.6f dforw = %.6f\n", xx, back, forw, dforw);
			}
		}
	}
		break;

	case 32:
	{
		double xx, skew, intercept;
		skew = (!args[0] ? 0.25 : atof(args[0]));
		intercept = (!args[1] ? 0.5 : atof(args[1]));
		printf("skew = %g\n", skew);
		printf("intercept = %g\n", intercept);
		double range = 9.0, dx = 0.1;
		double *arg[2];
		arg[0] = &skew;
		arg[1] = &intercept;

		for (int i = 0; i < (int) (2.0 * range / dx); i++) {
			xx = -range + i * dx;
			double a, b, c, d, h = 1e-6;
			a = map_invsn(xx, MAP_FORWARD, (void *) arg);
			b = map_invsn(a, MAP_BACKWARD, (void *) arg);
			c = map_invsn(xx, MAP_DFORWARD, (void *) arg);
			d = (map_invsn(xx + h, MAP_FORWARD, (void *) arg) - map_invsn(xx - h, MAP_FORWARD, (void *) arg)) / (2.0 * h);
			printf("xx = %.8g forw=%.8g backw=%.8g dforw=%.8g fdiff=%.8g (derr=%.8g)\n", xx, a, b, c, d, c - d);

		}
	}
		break;

	case 33:
	{
		double xx, yy, h = 1.0e-4, range = 1.123;

		for (xx = 1.2; xx < 3.0; xx += 0.35) {
			yy = map_phi(xx, MAP_FORWARD, NULL);
			printf("xx %g yy %g  xx.inv %g deriv %g dderiv %g\n",
			       xx, yy,
			       map_phi(yy, MAP_BACKWARD, NULL),
			       map_phi(xx, MAP_DFORWARD, NULL),
			       (map_phi(xx + h, MAP_FORWARD, NULL) - map_phi(xx - h, MAP_FORWARD, NULL)) / 2.0 / h);
		}
		for (xx = 1.2; xx < 3.0; xx += 0.35) {
			yy = map_phi(xx, MAP_FORWARD, (void *) &range);
			printf("xx %g yy %g  xx.inv %g deriv %g dderiv %g\n",
			       xx, yy,
			       map_phi(yy, MAP_BACKWARD, (void *) &range),
			       map_phi(xx, MAP_DFORWARD, (void *) &range),
			       (map_phi(xx + h, MAP_FORWARD, (void *) &range) - map_phi(xx - h, MAP_FORWARD, (void *) &range)) / 2.0 / h);
		}
	}
		break;

	case 34:
	{
		double theta, lambda = 40;

		for (theta = -5; theta <= 5; theta += 0.01) {
			printf("theta %g logprior %g\n", theta, priorfunc_pc_sn(&theta, &lambda));
		}
	}
		break;

	case 35:
	{
		double y, lambda = 1.234;

		for (y = 0.0; y <= 3.0; y += 0.1) {
			printf("y=%g  eval_log_contpoisson= %g\n", y, eval_log_contpoisson(y + 1, lambda));
		}
	}
		break;

	case 36:
	{
		GMRFLib_design_tp *design = Calloc(1, GMRFLib_design_tp);
		int nf = atoi(args[0]);

		GMRFLib_design_ccd(&design, nf);
		GMRFLib_design_print(stdout, design);

		GMRFLib_design_eb(&design, nf);
		GMRFLib_design_print(stdout, design);
	}
		break;

	case 37:
	{
		double xx, aa;
		aa = (!args[0] ? 1.0 : atof(args[0]));
		printf("aa = %g\n", aa);
		double range = 23.0, dx = 0.1;
		Link_param_tp *arg = Calloc(1, Link_param_tp);
		arg->a = aa;
		// paralle for must have int as loop-index
		// #pragma omp parallel for private(xx)
		for (int i = 0; i < (int) (2.0 * range / dx); i++) {
			int thread_id = 0;
			xx = -range + i * dx;
			double a, b, c, d, h = 1e-6;
			a = link_loga(thread_id, xx, MAP_FORWARD, (void *) arg, NULL);
			b = link_loga(thread_id, a, MAP_BACKWARD, (void *) arg, NULL);
			c = link_loga(thread_id, xx, MAP_DFORWARD, (void *) arg, NULL);
			d = (link_loga(thread_id, xx + h, MAP_FORWARD, (void *) arg, NULL) -
			     link_loga(thread_id, xx - h, MAP_FORWARD, (void *) arg, NULL)) / (2.0 * h);
			printf("xx = %.8f forw=%.12f backw=%.8f dforw=%.12f fdiff=%.12f (derr=%.12f)\n", xx, a, b, c, d, c - d);
		}
	}
		break;

	case 38:
		break;

	case 39:
	{
		double a, b, c, dd, x0 = 2.0;
		int stencil;
		int thread_id = 0;

		stencil = 5;
		GMRFLib_2order_taylor(thread_id, &a, &b, &c, &dd, 1.0, x0, 0, &x0, loglikelihood_testit, NULL, NULL, &stencil);
		printf("taylor: stencil= %d a= %.10g b= %.10g c= %.10g dd= %.10g\n", stencil, a, b, c, dd);
		stencil = 7;
		GMRFLib_2order_taylor(thread_id, &a, &b, &c, &dd, 1.0, x0, 0, &x0, loglikelihood_testit, NULL, NULL, &stencil);
		printf("taylor: stencil= %d a= %.10g b= %.10g c= %.10g dd= %.10g\n", stencil, a, b, c, dd);
		stencil = 9;
		GMRFLib_2order_taylor(thread_id, &a, &b, &c, &dd, 1.0, x0, 0, &x0, loglikelihood_testit, NULL, NULL, &stencil);
		printf("taylor: stencil= %d a= %.10g b= %.10g c= %.10g dd= %.10g\n", stencil, a, b, c, dd);

		stencil = 5;
		GMRFLib_2order_approx(thread_id, &a, &b, &c, &dd, 1.0, x0, 0, &x0, loglikelihood_testit, NULL, NULL, &stencil, NULL);
		printf("approx: stencil= %d a= %.10g b= %.10g c= %.10g dd= %.10g\n", stencil, a, b, c, dd);
		stencil = 7;
		GMRFLib_2order_approx(thread_id, &a, &b, &c, &dd, 1.0, x0, 0, &x0, loglikelihood_testit, NULL, NULL, &stencil, NULL);
		printf("approx: stencil= %d a= %.10g b= %.10g c= %.10g dd= %.10g\n", stencil, a, b, c, dd);
		stencil = 9;
		GMRFLib_2order_approx(thread_id, &a, &b, &c, &dd, 1.0, x0, 0, &x0, loglikelihood_testit, NULL, NULL, &stencil, NULL);
		printf("approx: stencil= %d a= %.10g b= %.10g c= %.10g dd= %.10g\n", stencil, a, b, c, dd);
	}
		break;

	case 40:
	{
		printf("eps= %.12g\n", GSL_DBL_EPSILON);
	}
		break;

	case 41:
	{
		inla_sn_intercept(0.43, 0.123);
		inla_sn_intercept(0.823, -0.123);
	}
		break;

	case 42:
		break;

	case 43:
	{
		gsl_matrix *A;
		A = gsl_matrix_alloc(3, 3);
		gsl_matrix_set(A, 0, 0, 1.0);
		gsl_matrix_set(A, 1, 0, 2.0);
		gsl_matrix_set(A, 2, 0, 2.0);

		gsl_matrix_set(A, 0, 1, -1.0);
		gsl_matrix_set(A, 1, 1, 0.0);
		gsl_matrix_set(A, 2, 1, 2.0);

		gsl_matrix_set(A, 0, 2, 0.0);
		gsl_matrix_set(A, 1, 2, 0.0);
		gsl_matrix_set(A, 2, 2, 1.0);

		GMRFLib_printf_gsl_matrix(stdout, A, " %.4f");
		GMRFLib_gsl_mgs(A);
		printf("\n");
		GMRFLib_printf_gsl_matrix(stdout, A, " %.4f");
	}
		break;

	case 44:
	{
		double a, b, y;
		int i;

		for (i = 0; i < 10; i++) {

			a = 2.0 * GMRFLib_uniform();
			b = 2.0 * GMRFLib_uniform();
			y = GMRFLib_uniform();

			printf("a %f b %f y %f", a, b, y);
			printf("  pbeta %f ", MATHLIB_FUN(pbeta) (y, a, b, 1, 1));
			printf("  1-pbeta %f\n", MATHLIB_FUN(pbeta) (y, a, b, 0, 1));
		}
	}
		break;

	case 45:
	{
		int n = 10;
		GMRFLib_idxval_tp *h = NULL;
		for (int i = 0, j = 0; i < n; i++) {
			j = (i < n / 2 ? 0 : 1);
			GMRFLib_idxval_add(&h, j, (double) j);
		}
		GMRFLib_idxval_sort(h);
		GMRFLib_idxval_printf(stdout, h, "case 45");
	}
		break;

	case 46:
	{
		GMRFLib_crwdef_tp *rw = Calloc(1, GMRFLib_crwdef_tp);
		GMRFLib_graph_tp *g;
		int n = 10, i, j;
		int thread_id = 0;
		assert(omp_get_thread_num() == 0);

		rw->n = n;
		rw->order = 2;
		rw->layout = GMRFLib_CRW_LAYOUT_SIMPLE;
		rw->position = Calloc(n, double);
		rw->position[0] = 0;
		for (i = 1; i < n; i++) {
			rw->position[i] = rw->position[i - 1] + i;
		}

		GMRFLib_make_crw_graph(&g, rw);

		double *len = Calloc(n, double);
		double *null = Calloc(n, double);
		double *res = Calloc(n, double);

		null[0] = 1.0;
		null[n - 1] = 1.0;
		len[0] = (rw->position[1] - rw->position[0]) / 1.0;
		len[n - 1] = (rw->position[n - 1] - rw->position[n - 2]) / 1.0;
		for (i = 1; i < n - 1; i++) {
			null[i] = 1.0;
			len[i] = (rw->position[i + 1] - rw->position[i]);
		}
		GMRFLib_Qx(thread_id, res, null, g, GMRFLib_crw, (void *) rw);
		for (i = 0; i < n; i++)
			printf("constr i %1d  x %f null %f Qnull %f\n", i, rw->position[i], null[i], res[i]);
		if (rw->order > 1) {
			Memcpy(null, rw->position, n * sizeof(double));
			GMRFLib_Qx(thread_id, res, null, g, GMRFLib_crw, (void *) rw);
			for (i = 0; i < n; i++)
				printf("linear i %1d  x %f null %f Qnull %f\n", i, rw->position[i], null[i], res[i]);
		}

		if (0) {
			gsl_matrix *Q = gsl_matrix_calloc(n, n);
			for (i = 0; i < n; i++)
				for (j = 0; j < n; j++) {
					if (i == j || GMRFLib_graph_is_nb(i, j, g)) {
						gsl_matrix_set(Q, i, j, GMRFLib_crw(thread_id, i, j, rw->position, rw));
					} else {
						gsl_matrix_set(Q, i, j, 0.0);
					}
				}

			gsl_matrix *vec = gsl_matrix_calloc(n, n);
			gsl_vector *val = gsl_vector_calloc(n);
			gsl_eigen_symmv_workspace *work = gsl_eigen_symmv_alloc(n);
			gsl_eigen_symmv(Q, val, vec, work);

			GMRFLib_printf_gsl_matrix(stdout, Q, " %8.4f");
			printf("\n");
			GMRFLib_printf_gsl_matrix(stdout, vec, " %8.4f");
			printf("\n");
			GMRFLib_printf_gsl_vector(stdout, val, " %8.4f");
		}
	}
		break;

	case 47:
	{
		GMRFLib_idxval_tp *h = NULL;
		GMRFLib_idxval_add(&h, 23, 1.000000);
		GMRFLib_idxval_add(&h, 61, 1.000000);
		GMRFLib_idxval_add(&h, 67, 1.000000);
		GMRFLib_idxval_add(&h, 127, 1.000000);
		GMRFLib_idxval_add(&h, 128, 1.000000);
		GMRFLib_idxval_add(&h, 162, 1.000000);
		GMRFLib_idxval_add(&h, 177, 1.000000);
		if (1) {
			GMRFLib_idxval_add(&h, 189, 1.000000);
			GMRFLib_idxval_add(&h, 223, 1.000000);
			GMRFLib_idxval_add(&h, 235, 1.000000);
			GMRFLib_idxval_add(&h, 360, 1.000000);
			GMRFLib_idxval_add(&h, 423, 1.000000);
			GMRFLib_idxval_add(&h, 432, 1.000000);
			GMRFLib_idxval_add(&h, 433, 1.000000);
			GMRFLib_idxval_add(&h, 443, 1.000000);
			GMRFLib_idxval_add(&h, 455, 1.000000);
			GMRFLib_idxval_add(&h, 492, 1.000000);
			GMRFLib_idxval_add(&h, 581, 1.000000);
			GMRFLib_idxval_add(&h, 620, 1.000000);
			GMRFLib_idxval_add(&h, 646, 1.000000);
			GMRFLib_idxval_add(&h, 675, 1.000000);
			GMRFLib_idxval_add(&h, 694, 1.000000);
			GMRFLib_idxval_add(&h, 730, 1.000000);
			GMRFLib_idxval_add(&h, 734, 1.000000);
			GMRFLib_idxval_add(&h, 754, 1.000000);
			GMRFLib_idxval_add(&h, 769, 1.000000);
			GMRFLib_idxval_add(&h, 781, 1.000000);
			GMRFLib_idxval_add(&h, 792, 1.000000);
			GMRFLib_idxval_add(&h, 822, 1.000000);
			GMRFLib_idxval_add(&h, 826, 1.000000);
			GMRFLib_idxval_add(&h, 896, 1.000000);
			GMRFLib_idxval_add(&h, 947, 1.000000);
			GMRFLib_idxval_add(&h, 966, 1.000000);
			GMRFLib_idxval_add(&h, 998, 1.000000);
			GMRFLib_idxval_add(&h, 1013, 1.000000);
		}

		GMRFLib_idxval_sort(h);
		GMRFLib_idxval_printf(stdout, h, "test47");
	}
		break;

	case 48:
	{
		for (double x = 1.0;; x *= 10.0) {
			printf("x= %f log(gsl_sf_psi_1(x)= %f  -log(x)= %f diff= %f\n",
			       x, log(gsl_sf_psi_1(x)), -log(x), log(gsl_sf_psi_1(x)) + log(x));
		}
	}
		break;

	case 49:
	{
#define SPECIAL(x) ((x > 0 ?						\
		     -2.0 * log(x) - log(2.0) + 1.0/(3.0*(x)) - 1.0/(18.0*SQR(x)) : \
		     log(gsl_sf_psi_1(x) - 1.0/(x))))

		for (double x = 1.0;; x *= 10.0) {
			printf("x= %f log(gsl_sf_psi_1(x)-1/x)= %f  %f %f\n",
			       x, log(gsl_sf_psi_1(x) - 1 / x), SPECIAL(x), log(gsl_sf_psi_1(x) - 1 / x) - SPECIAL(x));
		}
	}
#undef SPECIAL
		break;

	case 50:
	{
		if (nargs != 3) {
			printf("X DF NCP\n");
			exit(1);
		}
		double x = atof(args[0]);
		double df = atof(args[1]);
		double ncp = atof(args[2]);
		printf("R --vanilla --quiet -e 'df=%.8f; x= %.8f; ncp= %.8f; dchisq(x,df,ncp,log=TRUE); pchisq(x,df,ncp)'\n", df, x, ncp);
		printf("log.density= %.8f\n", MATHLIB_FUN(dnchisq) (x, df, ncp, 1));
		printf("CDF(%.8f)= %.8f\n", x, MATHLIB_FUN(pnchisq) (x, df, ncp, 1, 0));
		printf("iCDF(CDF(%.8f))= %.8f\n", x, MATHLIB_FUN(qnchisq) (MATHLIB_FUN(pnchisq) (x, df, ncp, 1, 0), df, ncp, 1, 0));

	}
		break;

	case 51:
	{
		dtweedie_init_cache();

		if (1) {
			double phi = 1.0 + GMRFLib_uniform();
			double xi = 1.0 + GMRFLib_uniform();
			double mu = exp(GMRFLib_uniform());
			double y = exp(1 + GMRFLib_uniform());
			double ldens;
			printf("R --vanilla --quiet -e 'library(tweedie);phi=%f;xi=%f;mu=%f;y=%f;dtweedie(y,xi,mu,phi);ptweedie(y,xi,mu,phi)'",
			       phi, xi, mu, y);
			dtweedie(1, y, &mu, phi, xi, &ldens);
			P(exp(ldens));
			P(ptweedie(y, mu, phi, xi));

			break;
		}

		double mu = 7.986;
		double phi = 1.717755;
		double p = 1.476;
		double ldens;
		double y, dy = 0.001, sum = 0.0;
		int iy;

#pragma omp parallel for private(iy, y, ldens) reduction(+:sum)
		for (iy = 0; iy < 100000; iy++) {
			y = SQR(dy) + iy * dy;
			dtweedie(1, y, &mu, phi, p, &ldens);
			printf("LDENS %f %f\n", y, ldens);
			sum += dy * exp(ldens);
		}
		y = 0.0;
		dtweedie(1, y, &mu, phi, p, &ldens);
		P(sum);
		P(exp(ldens));
		P(sum + exp(ldens));

		double lmu;
		y = 1;
		p = 1.51;
		for (lmu = -10; lmu < 10; lmu += 0.01) {
			mu = exp(lmu);
			dtweedie(1, y, &mu, phi, p, &ldens);
			printf("LMU %f %f\n", lmu, ldens);
		}

	}
		break;

	case 52:
	{
		for (double x = 0.1; x < 20; x += 0.1) {
			printf("x %f lgamma %f lgamma.fast %f diff %f\n",
			       x, gsl_sf_lngamma(x), inla_lgamma_fast(x), gsl_sf_lngamma(x) - inla_lgamma_fast(x));
		}

	}
		break;

	case 53:
	{
		double mu = 17.986;
		double phi = 1.717755;
		double p = 1.476;
		double ldens;
		int iy;

		for (iy = 10; iy < 100; iy += 10) {
			double pphi = phi / iy;
			double y = (double) iy;
			dtweedie(1, y, &mu, pphi, p, &ldens);
			printf("LDENS %f %f %f\n", y, pphi, ldens);
		}

	}
		break;

	case 54:
	{

		double x0 = 0, v;
		for (v = 0.0; v <= 1.0; v += 0.1) {
			printf("x %f exp %.12f exp_taylor %.12f %.12f\n", v, exp(v), exp_taylor(v, x0, 6), exp_taylor(v, x0, 12));
		}

	}
		break;

	case 55:
	{
		double skew3 = GMRFLib_skew_to_skew3(0.3);
		GMRFLib_snq_tp *q;
		int n = 31;

		q = GMRFLib_snq(n, skew3);
		for (int i = 0; i < q->n; i++) {
			printf("i %d x %.8f w %.8f ww %.8f www %.8f\n", i, q->nodes[i], q->w[i], q->w_grad[i], q->w_hess[i]);
		}

		double fun = 0, fund = 0, fundd = 0, fval = 0;
		for (int i = 0; i < q->n; i++) {
			fval = sin(q->nodes[i]);
			fun += fval * q->w[i];
			fund += fval * q->w_grad[i];
			fundd += fval * q->w_hess[i];
		}
		printf("sin(x): value= %.8f deriv= %.8f dderiv= %.8f\n", fun, fund, fundd);

		GMRFLib_snq_free(q);
	}
		break;

	case 56:
	{
		for (int i = 0, j = 1; i < 10; i++, j = j + 2) {
			double lambda = exp(-1 + GMRFLib_uniform());
			double nnew = inla_poisson_interval(lambda, i, j);
			double gsl = (gsl_cdf_poisson_P((unsigned) j, lambda) - (i <= 0 ? 0.0 : gsl_cdf_poisson_P((unsigned) (i - 1), lambda)));
			printf("lambda %f from= %d to= %d: nnew %f gsl %f diff %.12f\n", lambda, i, j, nnew, gsl, nnew - gsl);
		}
		// j < 0 <==> j=INF
		for (int i = 0, j = -1; i < 10; i++) {
			double lambda = exp(-1 + GMRFLib_uniform());
			double nnew = inla_poisson_interval(lambda, i, j);
			double gsl = (gsl_cdf_poisson_P((unsigned) 1000, lambda) - (i <= 0 ? 0.0 : gsl_cdf_poisson_P((unsigned) (i - 1), lambda)));
			printf("lambda %f from= %d to= %d: nnew %f gsl %f diff %.12f\n", lambda, i, j, nnew, gsl, nnew - gsl);
		}
	}
		break;

	case 57:
	{
		// testing qpoisson

		double eta, shape;
		Link_param_tp *lparam = Calloc(1, Link_param_tp);
		lparam->quantile = 0.87;
		int thread_id = 0;
		assert(omp_get_thread_num() == 0);

		P(lparam->quantile);
		for (eta = 1.0; eta < 10.0; eta += 0.5) {
			shape = exp(eta) + 1;
			printf("eta %f shape %f qpoisson %f qgamma %f\n", eta, shape,
			       link_qpoisson(thread_id, eta, INVLINK, lparam, NULL), MATHLIB_FUN(qgamma) (lparam->quantile, shape, 1.0, 0, 0));
		}
	}
		break;

	case 58:
	{
		double df = 3.0;
		double x = SQR(50);
		double ncp, ncp_sqrt;

		for (ncp_sqrt = 0.0; ncp_sqrt < 2 * sqrt(x); ncp_sqrt += 2.0 * sqrt(x) / 1.0E4) {
			ncp = SQR(ncp_sqrt);
			printf("sqrt(ncp) inla_dnchisq %f %f\n", sqrt(ncp), inla_dnchisq(x, df, ncp));
		}
	}
		break;

	case 59:
	{
		double x = GMRFLib_uniform();

		printf("x= %.12f\n", x);
		printf("Phi= %.12f\n", GMRFLib_cdfnorm(x));
		printf("Phi_inv = %.12f\n", GMRFLib_cdfnorm_inv(x));
		printf("erf = %.12f\n", GMRFLib_erf(x));
		printf("erfinv = %.12f\n", GMRFLib_erf_inv(x));
		printf("erfc = %.12f\n", GMRFLib_erfc(x));
		printf("erfcinv = %.12f\n", GMRFLib_erfc_inv(x));

		printf("%s%.12f%s%s\n", "R --vanilla --quiet -e 'library(pracma);x=",
		       x, ";print(x); print(pnorm(x)); print(qnorm(x)); print(erf(x));", "print(erfinv(x)); print(erfc(x)); print(erfcinv(x))'\n");
	}
		break;

	case 60:
	{
		for (int i = 0; i < 10; i++) {
			double p = GMRFLib_uniform();
			double a = GMRFLib_uniform() - 0.5;
			printf("p= %.12f\n", p);
			printf("a= %.12f\n", a);
			printf("sn_inv= %.12f\n", GMRFLib_sn_Pinv(p, a));

			printf("%s%.12f%s%.12f%s\n", "R --vanilla --quiet -e 'library(sn);p=", p, "; a=", a, "; print(qsn(p,alpha=a))'\n");
		}
	}
		break;

	case 61:
	{
		GMRFLib_problem_tp *problem;
		GMRFLib_graph_tp *g;
		GMRFLib_graph_mk_linear(&g, 5, 5, 0);

		int thread_id = 0;
		assert(omp_get_thread_num() == 0);
		GMRFLib_init_problem(thread_id, &problem, NULL, NULL, NULL, NULL, g, testit_Qfunc, NULL, NULL);
		GMRFLib_evaluate(problem);
		GMRFLib_Qinv(problem);

		for (int i = 0; i < g->n; i++) {
			printf("Qinv[%1d]=  %f\n", i, *GMRFLib_Qinv_get(problem, i, i));
		}
	}
		break;

	case 62:
	{
		double mom[3] = { 1.123, 0.123, -0.213 };
		GMRFLib_sn_param_tp p;

		printf("mom %f %f %f\n", mom[0], mom[1], mom[2]);
		GMRFLib_sn_moments2par(&p, mom, mom + 1, mom + 2);
		GMRFLib_sn_par2moments(mom, mom + 1, mom + 2, &p);
		printf("sn  %f %f %f\n", p.xi, p.omega, p.alpha);
		printf("mom %f %f %f\n", mom[0], mom[1], mom[2]);
	}
		break;

	case 63:
	{
		int n = atoi(args[0]);
		P(n);
		double *x = Calloc(n, double);
		double tref;
		for (int i = 0; i < n; i++) {
			x[i] = GMRFLib_uniform();
		}

		tref = GMRFLib_cpu();
		FILE *fp = fopen("REMOVE_ME_1.dat", "wb");
		fwrite((void *) x, sizeof(double), (size_t) n, fp);
		fclose(fp);
		printf("Optimal %f\n", (GMRFLib_cpu() - tref));

		tref = GMRFLib_cpu();
		fp = fopen("REMOVE_ME_2.dat", "wb");
		for (int i = 0; i < n; i++)
			fwrite((void *) &x[i], sizeof(double), (size_t) 1, fp);
		fclose(fp);
		printf("One-by-one %f\n", (GMRFLib_cpu() - tref));

		tref = GMRFLib_cpu();
		{
			Dinit("REMOVE_ME_3.dat");
			for (int i = 0; i < n; i++) {
				D1W(x[i]);
			}
			Dclose();
			printf("D-cache %f\n", (GMRFLib_cpu() - tref));
		}

		tref = GMRFLib_cpu();
		{
			Dinit_s("REMOVE_ME_4.dat");
			for (int i = 0; i < n; i++) {
				D1W(x[i]);
			}
			Dclose();
			printf("D-cache short %f\n", (GMRFLib_cpu() - tref));
		}

		tref = GMRFLib_cpu();
		{
			fp = fopen("REMOVE_ME_5.dat", "wb");
			char *buff = (char *) Calloc(16777216L, double);
			setvbuf(stdout, buff, _IOFBF, 16777216L * sizeof(double));
			for (int i = 0; i < n; i++) {
				fwrite(x + i, sizeof(double), (size_t) 1, fp);
			}
			fclose(fp);
			printf("setvbuf %f\n", (GMRFLib_cpu() - tref));
		}
	}
		break;

	case 64:
	{
		int n = atoi(args[0]);
		gsl_bfgs4_test1((size_t) n);
	}
		break;

	case 65:
	{
		double x[] = { -0.200, -0.075, 0.000, 0.040, 0.160, 0.360, 0.640, 1.000 };
		double y[] = { 14960.675457, 14934.327851, 14927.976542, 14943.616530, 14945.530949, 15000.597367, 15140.861227, 15412.165108 };
		int n = sizeof(y) / sizeof(double);

		double xmin;
		double ymin;

		for (int order = 2; order <= 4; order += 2) {
			bfgs4_robust_minimize(&xmin, &ymin, n, x, y, 0, NULL, NULL, order);
			printf("xmin = %f ymin= %f when order = %d\n", xmin, ymin, order);
		}
	}
		break;

	case 66:
	{
		double power, power_intern;
		double intercept, intercept_intern;
		double **param;
		param = Calloc(2, double *);

		power = 1.5;
		power_intern = map_exp(power, MAP_BACKWARD, NULL);
		intercept = 0.75;
		intercept_intern = map_probability(intercept, MAP_BACKWARD, NULL);

		param[0] = &power_intern;
		param[1] = &intercept_intern;

		map_inv_powerlink_core(0.0, MAP_FORWARD, (void *) param, NULL);
	}
		break;

	case 67:
	{
		double xx, power, intercept;
		power = (!args[0] ? 0.01 : atof(args[0]));
		intercept = (!args[1] ? 0.01 : atof(args[1]));
		printf("power = %g\n", power);
		printf("intercept = %g\n", intercept);
		double range = 2.0, dx = 0.2;
		double *arg[2];
		arg[0] = &power;
		arg[1] = &intercept;

		for (int i = 0; i < (int) (2.0 * range / dx); i++) {
			xx = -range + i * dx;
			double a, b, c, d, h = 1e-6;
			a = map_inv_powerlink_core(xx, MAP_FORWARD, (void *) arg, NULL);
			b = map_inv_powerlink_core(a, MAP_BACKWARD, (void *) arg, NULL);
			c = map_inv_powerlink_core(xx, MAP_DFORWARD, (void *) arg, NULL);
			d = (map_inv_powerlink_core(xx + h, MAP_FORWARD, (void *) arg, NULL)
			     - map_inv_powerlink_core(xx - h, MAP_FORWARD, (void *) arg, NULL)) / (2.0 * h);
			printf("xx = %.8g forw=%.8g backw=%.8g dforw=%.8g fdiff=%.8g (derr=%.8g)\n", xx, a, b, c, d, c - d);

		}
	}
		break;

	case 68:
	{
		assert(nargs == 3);
		printf("Call 'double (*fun)(double)' function [%s] in [%s] with argument [%s]\n", args[0], args[1], args[2]);
		lt_dlhandle handle;
		typedef double fun_tp(double);
		fun_tp *fun = NULL;
		const char *error;

		lt_dlinit();
		handle = lt_dlopen(args[1]);
		if (!handle) {
			fprintf(stderr, "%s\n", lt_dlerror());
			exit(1);
		}
		lt_dlerror();

		fun = (fun_tp *) lt_dlsym(handle, args[0]);
		if ((error = lt_dlerror()) != NULL) {
			fprintf(stderr, "%s\n", error);
			exit(1);
		}
		lt_dlerror();

		double x = atof(args[2]);
		printf("fun(%g) = %g\n", x, fun(x));
		lt_dlclose(handle);
	}
		break;

	case 69:
	{
		int n = 10;

		printf("\n\n");
		FIXME("run with 4 threads");
#pragma omp parallel for num_threads(4)
		for (int i = 0; i < n; i++) {
			P(GMRFLib_OPENMP_IN_SERIAL());
			P(GMRFLib_OPENMP_IN_PARALLEL());
			P(GMRFLib_OPENMP_IN_PARALLEL_ONE_THREAD());
			P(GMRFLib_OPENMP_IN_PARALLEL_ONEPLUS_THREAD());
		}

		printf("\n\n");
		FIXME("run with 1 threads");
#pragma omp parallel for num_threads(1)
		for (int i = 0; i < n; i++) {
			P(GMRFLib_OPENMP_IN_SERIAL());
			P(GMRFLib_OPENMP_IN_PARALLEL());
			P(GMRFLib_OPENMP_IN_PARALLEL_ONE_THREAD());
			P(GMRFLib_OPENMP_IN_PARALLEL_ONEPLUS_THREAD());
		}

		printf("\n\n");
		FIXME("run serial");
		for (int i = 0; i < n; i++) {
			P(GMRFLib_OPENMP_IN_SERIAL());
			P(GMRFLib_OPENMP_IN_PARALLEL());
			P(GMRFLib_OPENMP_IN_PARALLEL_ONE_THREAD());
			P(GMRFLib_OPENMP_IN_PARALLEL_ONEPLUS_THREAD());
		}
	}
		break;

	case 70:
	{
		GMRFLib_design_tp *design = Calloc(1, GMRFLib_design_tp);
		int nf = atoi(args[0]);
		GMRFLib_design_grid(&design, nf);
		GMRFLib_design_print(stdout, design);
	}
		break;

	case 71:
	{
		const int n = 51;
		double x[n];
		double ld[n];

		for (int i = 0; i < n; i++) {
			x[i] = -5 + (double) i / (n - 1.0) * 10;
			ld[i] = (-0.5 * SQR(x[i]));
		}

		GMRFLib_density_tp *dens = NULL;
		GMRFLib_density_create(&dens, GMRFLib_DENSITY_TYPE_SCGAUSSIAN, n, x, ld, 0.0, 1.0, 1);
		GMRFLib_density_printf(stdout, dens);
	}
		break;

	case 72:
	{
		P(omp_get_num_threads());
		P(omp_get_max_threads());
		P(omp_get_thread_num());
		P(GMRFLib_OPENMP_IN_SERIAL());
	}
		break;

	case 73:
	{
		const int n = 51;
		double x[n];
		double ld[n];

		for (int i = 0; i < n; i++) {
			x[i] = (double) i / (n - 1.0) * 10;
			ld[i] = -0.5 * SQR(x[i]);
		}

		GMRFLib_density_tp *dens = NULL;
		GMRFLib_density_tp *dens_dup = NULL;
		GMRFLib_sn_param_tp sn_par = { 0., 1.0, 0.4 };

		GMRFLib_density_create_sn(&dens, sn_par, 1.1, 2.2, 1);
		FIXME("SN");
		GMRFLib_density_printf(stdout, dens);
		GMRFLib_density_duplicate(&dens_dup, dens);
		FIXME("SN dup");
		GMRFLib_density_printf(stdout, dens_dup);

		GMRFLib_free_density(dens);
		GMRFLib_free_density(dens_dup);

		GMRFLib_density_create(&dens, GMRFLib_DENSITY_TYPE_SCGAUSSIAN, n, x, ld, 0.0, 1.0, 1);
		FIXME("SC Gaussian");
		GMRFLib_density_duplicate(&dens_dup, dens);
		GMRFLib_density_printf(stdout, dens);
		FIXME("SC Gaussian dup");
		GMRFLib_density_printf(stdout, dens_dup);

		GMRFLib_free_density(dens);
		GMRFLib_free_density(dens_dup);
	}
		break;

	case 74:
	{
		double x[100];
		double *p, *pp;

		p = &x[0];
		pp = &x[10];
		P(OVERLAP(p, pp, 5));
		P(OVERLAP(p, pp, 9));
		P(OVERLAP(p, pp, 10));
		P(OVERLAP(p, pp, 11));
		P(OVERLAP(p, pp, 15));
	}
		break;

	case 75:
	{
		int n = 100, i;
		double *x = Calloc(n, double);
		double *pp = Calloc(n, double);
		double xx, dx;
		double xmax = 5;
		dx = 2.0 * xmax / (n - 1);

		for (i = 0; i < n; i++) {
			xx = -xmax + dx * i;
			x[i] = xx;
			pp[i] = inla_Phi(xx);
			// printf("%d %.20f %.20f\n", i, x[i], pp[i]); 
		}

		GMRFLib_spline_tp *P, *Pinv;

		P = GMRFLib_spline_create_x(x, pp, n, GMRFLib_INTPOL_TRANS_P, GMRFLib_INTPOL_CACHE_LEVEL12);
		Pinv = GMRFLib_spline_create_x(pp, x, n, GMRFLib_INTPOL_TRANS_Pinv, GMRFLib_INTPOL_CACHE_LEVEL12);

		for (xx = -(xmax + 2); xx <= (xmax + 2); xx += 0.5) {
			double p1 = inla_Phi(xx);
			double p2 = GMRFLib_spline_eval(xx, P);
			double xx2 = GMRFLib_spline_eval(p2, Pinv);
			printf("XX %.20f %.20f %.20f %.20f\n", xx, p1, p2, xx2);
		}
	}
		break;

	case 76:
	{
		for (int i = 0; i < 10; i++) {
			printf("%d %f %f\n", i, gsl_sf_lnfact((unsigned int) i), my_gsl_sf_lnfact(i));
		}
	}
		break;

	case 77:
	{
		int n = 111;
		double dx = 12.0 / (n - 1);
		Calloc_init(2 * n, 2);
		double *x = Calloc_get(n);
		double *y = Calloc_get(n);

		double z = 0.0;
		for (int i = 0; i < n; i++) {
			x[i] = -6.0 + i * dx;
			// this one is normalized
			y[i] = log(2.0) + log(1.0 / sqrt(2.0 * M_PI)) - 0.5 * SQR(x[i]) + inla_log_Phi(1.0 * x[i]);
			z += dx * exp(y[i]);
		}
		P(z);					       /* should be 1 */
		GMRFLib_density_tp *density;
		GMRFLib_density_create(&density, GMRFLib_DENSITY_TYPE_SCGAUSSIAN, n, x, y, 0.0, 1.0, 1);
		GMRFLib_density_printf(stdout, density);

		for (int i = 0; i < n; i++) {
			x[i] = -6.0 + i * dx;
			double yy;
			GMRFLib_evaluate_logdensity(&yy, x[i], density);
			printf("Evaluate x %f true %f scg %f\n", x[i], y[i], yy);
		}
	}
		break;

	case 78:
	{
		/*
		 * 
		 * library(mvtnorm) Q <- matrix(c(17, -1, -2, -1, 15, -3, -2, -3, 14), 3, 3) S <- solve(Q) x <- c(0.55, 1.55, 2.55) m <- c(2.05,
		 * 3.45, 1.25) zero <- rep(0, length(x)) print(dmvnorm(x = x, mean = m, sigma = S, log = TRUE)) print(dmvnorm(x = x, mean = zero,
		 * sigma = S, log = TRUE)) print(dmvnorm(x = zero, mean = m, sigma = S, log = TRUE)) print(dmvnorm(x = zero, mean = zero, sigma =
		 * S, log = TRUE))
		 * 
		 */

		gsl_matrix *Q = gsl_matrix_alloc(3, 3);
		gsl_matrix_set(Q, 0, 0, 17);
		gsl_matrix_set(Q, 1, 1, 15);
		gsl_matrix_set(Q, 2, 2, 14);
		gsl_matrix_set(Q, 0, 1, -1);
		gsl_matrix_set(Q, 1, 0, -1);
		gsl_matrix_set(Q, 0, 2, -2);
		gsl_matrix_set(Q, 2, 0, -2);
		gsl_matrix_set(Q, 1, 2, -3);
		gsl_matrix_set(Q, 2, 1, -3);

		gsl_matrix *S = GMRFLib_gsl_duplicate_matrix(Q);
		GMRFLib_gsl_spd_inverse(S);

		gsl_vector *x = gsl_vector_alloc(3);
		gsl_vector_set(x, 0, 0.55);
		gsl_vector_set(x, 1, 1.55);
		gsl_vector_set(x, 2, 2.55);

		gsl_vector *mean = gsl_vector_alloc(3);
		gsl_vector_set(mean, 0, 2.05);
		gsl_vector_set(mean, 1, 3.45);
		gsl_vector_set(mean, 2, 1.25);

		P(GMRFLib_gsl_log_dnorm(x, mean, NULL, S, 0));
		P(GMRFLib_gsl_log_dnorm(x, mean, Q, NULL, 0));

		P(GMRFLib_gsl_log_dnorm(x, NULL, NULL, S, 0));
		P(GMRFLib_gsl_log_dnorm(x, NULL, Q, NULL, 0));

		P(GMRFLib_gsl_log_dnorm(NULL, mean, NULL, S, 0));
		P(GMRFLib_gsl_log_dnorm(NULL, mean, Q, NULL, 0));

		P(GMRFLib_gsl_log_dnorm(NULL, NULL, NULL, S, 0));
		P(GMRFLib_gsl_log_dnorm(NULL, NULL, Q, NULL, 0));
	}
		break;

	case 79:
	{
		/*
		 * Q <- matrix(c(0.1, -1, -2, -1, 15, -3, -2, -3, 14), 3, 3) 
		 */

		gsl_matrix *Q = gsl_matrix_alloc(3, 3);
		gsl_matrix_set(Q, 0, 0, .1);
		gsl_matrix_set(Q, 1, 1, 15);
		gsl_matrix_set(Q, 2, 2, 14);
		gsl_matrix_set(Q, 0, 1, -1);
		gsl_matrix_set(Q, 1, 0, -1);
		gsl_matrix_set(Q, 0, 2, -2);
		gsl_matrix_set(Q, 2, 0, -2);
		gsl_matrix_set(Q, 1, 2, -3);
		gsl_matrix_set(Q, 2, 1, -3);

		gsl_matrix *S = GMRFLib_gsl_duplicate_matrix(Q);
		GMRFLib_gsl_ensure_spd_inverse(S, GSL_SQRT_DBL_EPSILON, NULL);
		GMRFLib_printf_gsl_matrix(stdout, Q, " %.8f");
		GMRFLib_printf_gsl_matrix(stdout, S, " %.8f");
		GMRFLib_gsl_ensure_spd_inverse(S, GSL_SQRT_DBL_EPSILON, NULL);
		GMRFLib_printf_gsl_matrix(stdout, S, " %.8f");
	}
		break;

	case 80:
	{
		/*
		 * fun = function(x) -x^2 + (x-1)^3 - 0.5 * (x+1)^4 plot(xx, fun(xx)) xx[which.max(fun(xx))] [1] 0.07507507508 
		 */

		double xx1[] = { 0, 1, 2 };
		double yy1[] = { -1.5, -9.0, -43.5 };
		double xm;
		xm = inla_interpolate_mode(xx1, yy1);
		P(xm);

		double xx2[] = { -1, 0, 0.5 };
		double yy2[] = { -9.00000, -1.50000, -2.90625 };

		xm = inla_interpolate_mode(xx2, yy2);
		P(xm);

		double xx3[] = { 0, 0.1, 0.12 };
		double yy3[] = { -1.50000000, -1.47105000, -1.48263168 };

		xm = inla_interpolate_mode(xx3, yy3);
		P(xm);

		double xx4[] = { 0.12, 0.0, 0.1 };
		double yy4[] = { -1.48263168, -1.50000000, -1.47105000 };

		xm = inla_interpolate_mode(xx4, yy4);
		P(xm);
	}
		break;

	case 81:
	{
		GMRFLib_idxval_tp *h = NULL;
		GMRFLib_idxval_add(&h, 22830, 1);
		GMRFLib_idxval_add(&h, 22832, 1);
		GMRFLib_idxval_add(&h, 22847, 1);
		GMRFLib_idxval_add(&h, 22850, 1);
		GMRFLib_idxval_add(&h, 22856, 1);
		GMRFLib_idxval_add(&h, 22861, 1);
		GMRFLib_idxval_add(&h, 22869, 1);
		GMRFLib_idxval_add(&h, 22877, 1);
		GMRFLib_idxval_add(&h, 22885, 1);
		GMRFLib_idxval_add(&h, 22892, 1);
		GMRFLib_idxval_add(&h, 22893, 1);
		GMRFLib_idxval_add(&h, 22904, 1);
		GMRFLib_idxval_add(&h, 22905, 1);
		GMRFLib_idxval_add(&h, 22918, 1);
		GMRFLib_idxval_add(&h, 22922, 1);
		GMRFLib_idxval_add(&h, 22933, 1);
		GMRFLib_idxval_add(&h, 22946, 1);
		GMRFLib_idxval_add(&h, 22949, 1);
		GMRFLib_idxval_add(&h, 22950, 1);
		GMRFLib_idxval_add(&h, 22965, 1);
		GMRFLib_idxval_add(&h, 22969, 1);
		GMRFLib_idxval_add(&h, 22980, 1);
		GMRFLib_idxval_add(&h, 22982, 1);
		GMRFLib_idxval_add(&h, 22983, 1);
		GMRFLib_idxval_add(&h, 22995, 1);
		GMRFLib_idxval_add(&h, 23009, 1);
		GMRFLib_idxval_add(&h, 23014, 1);
		GMRFLib_idxval_add(&h, 23015, 1);
		GMRFLib_idxval_add(&h, 23017, 1);
		GMRFLib_idxval_add(&h, 23032, 1);
		GMRFLib_idxval_add(&h, 23033, 1);
		GMRFLib_idxval_add(&h, 23045, 1);
		GMRFLib_idxval_add(&h, 23060, 1);
		GMRFLib_idxval_add(&h, 23070, 1);
		GMRFLib_idxval_add(&h, 23084, 1);
		GMRFLib_idxval_add(&h, 23093, 1);
		GMRFLib_idxval_add(&h, 23106, 1);
		GMRFLib_idxval_add(&h, 23107, 1);
		GMRFLib_idxval_add(&h, 23117, 1);
		GMRFLib_idxval_add(&h, 23124, 1);
		GMRFLib_idxval_add(&h, 23139, 1);
		GMRFLib_idxval_add(&h, 23143, 1);
		GMRFLib_idxval_add(&h, 23158, 1);
		GMRFLib_idxval_add(&h, 23173, 1);
		GMRFLib_idxval_add(&h, 23183, 1);
		GMRFLib_idxval_add(&h, 23197, 1);
		GMRFLib_idxval_add(&h, 23204, 1);
		GMRFLib_idxval_add(&h, 23214, 1);
		GMRFLib_idxval_add(&h, 23229, 1);
		GMRFLib_idxval_add(&h, 23232, 1);
		GMRFLib_idxval_add(&h, 23240, 1);
		GMRFLib_idxval_add(&h, 23252, 1);
		GMRFLib_idxval_add(&h, 23259, 1);
		GMRFLib_idxval_add(&h, 23262, 1);
		GMRFLib_idxval_add(&h, 23267, 1);
		GMRFLib_idxval_add(&h, 23271, 1);
		GMRFLib_idxval_add(&h, 23277, 1);
		GMRFLib_idxval_add(&h, 23287, 1);
		GMRFLib_idxval_add(&h, 23302, 1);
		GMRFLib_idxval_add(&h, 23309, 1);
		GMRFLib_idxval_add(&h, 23324, 1);
		GMRFLib_idxval_add(&h, 23339, 1);
		GMRFLib_idxval_add(&h, 23348, 1);
		GMRFLib_idxval_add(&h, 23351, 1);
		GMRFLib_idxval_add(&h, 23356, 1);
		GMRFLib_idxval_add(&h, 23370, 1);
		GMRFLib_idxval_add(&h, 23383, 1);
		GMRFLib_idxval_add(&h, 23391, 1);
		GMRFLib_idxval_add(&h, 23405, 1);
		GMRFLib_idxval_add(&h, 23409, 1);
		GMRFLib_idxval_add(&h, 23420, 1);
		GMRFLib_idxval_add(&h, 23434, 1);
		GMRFLib_idxval_add(&h, 23448, 1);
		GMRFLib_idxval_add(&h, 23462, 1);
		GMRFLib_idxval_add(&h, 23468, 1);
		GMRFLib_idxval_add(&h, 23472, 1);
		GMRFLib_idxval_add(&h, 23473, 1);
		GMRFLib_idxval_add(&h, 23474, 1);
		GMRFLib_idxval_add(&h, 23479, 1);
		GMRFLib_idxval_add(&h, 23492, 1);
		GMRFLib_idxval_add(&h, 23495, 1);
		GMRFLib_idxval_add(&h, 23503, 1);
		GMRFLib_idxval_add(&h, 23514, 1);
		GMRFLib_idxval_add(&h, 23528, 1);
		GMRFLib_idxval_add(&h, 23542, 1);
		GMRFLib_idxval_add(&h, 23552, 1);
		GMRFLib_idxval_add(&h, 23557, 1);
		GMRFLib_idxval_add(&h, 23558, 1);
		GMRFLib_idxval_add(&h, 23567, 1);
		GMRFLib_idxval_add(&h, 23572, 1);
		GMRFLib_idxval_add(&h, 23580, 1);
		GMRFLib_idxval_add(&h, 23582, 1);
		GMRFLib_idxval_add(&h, 23584, 1);
		GMRFLib_idxval_add(&h, 23598, 1);
		GMRFLib_idxval_add(&h, 23601, 1);
		GMRFLib_idxval_add(&h, 23602, 1);
		GMRFLib_idxval_add(&h, 23615, 1);
		GMRFLib_idxval_add(&h, 23617, 1);
		GMRFLib_idxval_add(&h, 23630, 1);
		GMRFLib_idxval_add(&h, 23643, 1);
		GMRFLib_idxval_add(&h, 23645, 1);
		GMRFLib_idxval_add(&h, 23658, 1);
		GMRFLib_idxval_add(&h, 23669, 1);
		GMRFLib_idxval_add(&h, 23670, 1);
		GMRFLib_idxval_add(&h, 23679, 1);
		GMRFLib_idxval_add(&h, 23688, 1);
		GMRFLib_idxval_add(&h, 23698, 1);
		GMRFLib_idxval_add(&h, 23710, 1);
		GMRFLib_idxval_add(&h, 23718, 1);
		GMRFLib_idxval_add(&h, 23728, 1);
		GMRFLib_idxval_add(&h, 23734, 1);
		GMRFLib_idxval_add(&h, 23741, 1);
		GMRFLib_idxval_add(&h, 23751, 1);
		GMRFLib_idxval_add(&h, 23754, 1);
		GMRFLib_idxval_add(&h, 23764, 1);
		GMRFLib_idxval_add(&h, 23776, 1);
		GMRFLib_idxval_add(&h, 23788, 1);
		GMRFLib_idxval_add(&h, 23792, 1);
		GMRFLib_idxval_add(&h, 23799, 1);
		GMRFLib_idxval_add(&h, 23801, 1);
		GMRFLib_idxval_add(&h, 23807, 1);
		GMRFLib_idxval_add(&h, 23808, 1);
		GMRFLib_idxval_add(&h, 23820, 1);
		GMRFLib_idxval_add(&h, 23823, 1);
		GMRFLib_idxval_add(&h, 23835, 1);
		GMRFLib_idxval_add(&h, 23843, 1);
		GMRFLib_idxval_add(&h, 23846, 1);
		GMRFLib_idxval_add(&h, 23858, 1);
		GMRFLib_idxval_add(&h, 23861, 1);
		GMRFLib_idxval_add(&h, 23873, 1);
		GMRFLib_idxval_add(&h, 23878, 1);
		GMRFLib_idxval_add(&h, 23881, 1);
		GMRFLib_idxval_add(&h, 23892, 1);
		GMRFLib_idxval_add(&h, 23901, 1);
		GMRFLib_idxval_add(&h, 23912, 1);
		GMRFLib_idxval_add(&h, 23923, 1);
		GMRFLib_idxval_add(&h, 23934, 1);
		GMRFLib_idxval_add(&h, 23945, 1);
		GMRFLib_idxval_add(&h, 23949, 1);
		GMRFLib_idxval_add(&h, 23956, 1);
		GMRFLib_idxval_add(&h, 23967, 1);
		GMRFLib_idxval_add(&h, 23978, 1);
		GMRFLib_idxval_add(&h, 23985, 1);
		GMRFLib_idxval_add(&h, 23988, 1);
		GMRFLib_idxval_add(&h, 23991, 1);
		GMRFLib_idxval_add(&h, 24002, 1);
		GMRFLib_idxval_add(&h, 24011, 1);
		GMRFLib_idxval_add(&h, 24021, 1);
		GMRFLib_idxval_add(&h, 24026, 1);
		GMRFLib_idxval_add(&h, 24029, 1);
		GMRFLib_idxval_add(&h, 24040, 1);
		GMRFLib_idxval_add(&h, 24051, 1);
		GMRFLib_idxval_add(&h, 24055, 1);
		GMRFLib_idxval_add(&h, 24066, 1);
		GMRFLib_idxval_add(&h, 24067, 1);
		GMRFLib_idxval_add(&h, 24078, 1);
		GMRFLib_idxval_add(&h, 24081, 1);
		GMRFLib_idxval_add(&h, 24092, 1);
		GMRFLib_idxval_add(&h, 24100, 1);
		GMRFLib_idxval_add(&h, 24105, 1);
		GMRFLib_idxval_add(&h, 24115, 1);
		GMRFLib_idxval_add(&h, 24125, 1);
		GMRFLib_idxval_add(&h, 24126, 1);
		GMRFLib_idxval_add(&h, 24132, 1);
		GMRFLib_idxval_add(&h, 24133, 1);
		GMRFLib_idxval_add(&h, 24137, 1);
		GMRFLib_idxval_add(&h, 24147, 1);
		GMRFLib_idxval_add(&h, 24152, 1);
		GMRFLib_idxval_add(&h, 24162, 1);
		GMRFLib_idxval_add(&h, 24166, 1);
		GMRFLib_idxval_add(&h, 24176, 1);
		GMRFLib_idxval_add(&h, 24186, 1);
		GMRFLib_idxval_add(&h, 24195, 1);
		GMRFLib_idxval_add(&h, 24205, 1);
		GMRFLib_idxval_add(&h, 24215, 1);
		GMRFLib_idxval_add(&h, 24225, 1);
		GMRFLib_idxval_add(&h, 24230, 1);
		GMRFLib_idxval_add(&h, 24240, 1);
		GMRFLib_idxval_add(&h, 24250, 1);
		GMRFLib_idxval_add(&h, 24260, 1);
		GMRFLib_idxval_add(&h, 24270, 1);
		GMRFLib_idxval_add(&h, 24280, 1);
		GMRFLib_idxval_add(&h, 24290, 1);
		GMRFLib_idxval_add(&h, 24297, 1);
		GMRFLib_idxval_add(&h, 24300, 1);
		GMRFLib_idxval_add(&h, 24309, 1);
		GMRFLib_idxval_add(&h, 24314, 1);
		GMRFLib_idxval_add(&h, 24317, 1);
		GMRFLib_idxval_add(&h, 24325, 1);
		GMRFLib_idxval_add(&h, 24334, 1);
		GMRFLib_idxval_add(&h, 24344, 1);
		GMRFLib_idxval_add(&h, 24345, 1);
		GMRFLib_idxval_add(&h, 24354, 1);
		GMRFLib_idxval_add(&h, 24357, 1);
		GMRFLib_idxval_add(&h, 24366, 1);
		GMRFLib_idxval_add(&h, 24375, 1);
		GMRFLib_idxval_add(&h, 24383, 1);
		GMRFLib_idxval_add(&h, 24392, 1);
		GMRFLib_idxval_add(&h, 24401, 1);
		GMRFLib_idxval_add(&h, 24410, 1);
		GMRFLib_idxval_add(&h, 24419, 1);
		GMRFLib_idxval_add(&h, 24428, 1);
		GMRFLib_idxval_add(&h, 24437, 1);
		GMRFLib_idxval_add(&h, 24446, 1);
		GMRFLib_idxval_add(&h, 24452, 1);
		GMRFLib_idxval_add(&h, 24458, 1);
		GMRFLib_idxval_add(&h, 24467, 1);
		GMRFLib_idxval_add(&h, 24476, 1);
		GMRFLib_idxval_add(&h, 24479, 1);
		GMRFLib_idxval_add(&h, 24488, 1);
		GMRFLib_idxval_add(&h, 24497, 1);
		GMRFLib_idxval_add(&h, 24506, 1);
		GMRFLib_idxval_add(&h, 24515, 1);
		GMRFLib_idxval_add(&h, 24522, 1);
		GMRFLib_idxval_add(&h, 24528, 1);
		GMRFLib_idxval_add(&h, 24532, 1);
		GMRFLib_idxval_add(&h, 24539, 1);
		GMRFLib_idxval_add(&h, 24542, 1);
		GMRFLib_idxval_add(&h, 24551, 1);
		GMRFLib_idxval_add(&h, 24560, 1);
		GMRFLib_idxval_add(&h, 24564, 1);
		GMRFLib_idxval_add(&h, 24572, 1);
		GMRFLib_idxval_add(&h, 24574, 1);
		GMRFLib_idxval_add(&h, 24575, 1);
		GMRFLib_idxval_add(&h, 24583, 1);
		GMRFLib_idxval_add(&h, 24591, 1);
		GMRFLib_idxval_add(&h, 24599, 1);
		GMRFLib_idxval_add(&h, 24602, 1);
		GMRFLib_idxval_add(&h, 24610, 1);
		GMRFLib_idxval_add(&h, 24611, 1);
		GMRFLib_idxval_add(&h, 24619, 1);
		GMRFLib_idxval_add(&h, 24623, 1);
		GMRFLib_idxval_add(&h, 24631, 1);
		GMRFLib_idxval_add(&h, 24639, 1);
		GMRFLib_idxval_add(&h, 24647, 1);
		GMRFLib_idxval_add(&h, 24654, 1);
		GMRFLib_idxval_add(&h, 24662, 1);
		GMRFLib_idxval_add(&h, 24668, 1);
		GMRFLib_idxval_add(&h, 24674, 1);
		GMRFLib_idxval_add(&h, 24676, 1);
		GMRFLib_idxval_add(&h, 24684, 1);
		GMRFLib_idxval_add(&h, 24687, 1);
		GMRFLib_idxval_add(&h, 24694, 1);
		GMRFLib_idxval_add(&h, 24697, 1);
		GMRFLib_idxval_add(&h, 24702, 1);
		GMRFLib_idxval_add(&h, 24710, 1);
		GMRFLib_idxval_add(&h, 24715, 1);
		GMRFLib_idxval_add(&h, 24718, 1);
		GMRFLib_idxval_add(&h, 24726, 1);
		GMRFLib_idxval_add(&h, 24734, 1);
		GMRFLib_idxval_add(&h, 24741, 1);
		GMRFLib_idxval_add(&h, 24748, 1);
		GMRFLib_idxval_add(&h, 24756, 1);
		GMRFLib_idxval_add(&h, 24763, 1);
		GMRFLib_idxval_add(&h, 24766, 1);
		GMRFLib_idxval_add(&h, 24773, 1);
		GMRFLib_idxval_add(&h, 24780, 1);
		GMRFLib_idxval_add(&h, 24787, 1);
		GMRFLib_idxval_add(&h, 24794, 1);
		GMRFLib_idxval_add(&h, 24801, 1);
		GMRFLib_idxval_add(&h, 24806, 1);
		GMRFLib_idxval_add(&h, 24813, 1);
		GMRFLib_idxval_add(&h, 24819, 1);
		GMRFLib_idxval_add(&h, 24823, 1);
		GMRFLib_idxval_add(&h, 24828, 1);
		GMRFLib_idxval_add(&h, 24832, 1);
		GMRFLib_idxval_add(&h, 24839, 1);
		GMRFLib_idxval_add(&h, 24840, 1);
		GMRFLib_idxval_add(&h, 24844, 1);
		GMRFLib_idxval_add(&h, 24850, 1);
		GMRFLib_idxval_add(&h, 24856, 1);
		GMRFLib_idxval_add(&h, 24863, 1);
		GMRFLib_idxval_add(&h, 24870, 1);
		GMRFLib_idxval_add(&h, 24875, 1);
		GMRFLib_idxval_add(&h, 24880, 1);
		GMRFLib_idxval_add(&h, 24886, 1);
		GMRFLib_idxval_add(&h, 24893, 1);
		GMRFLib_idxval_add(&h, 24900, 1);
		GMRFLib_idxval_add(&h, 24906, 1);
		GMRFLib_idxval_add(&h, 24913, 1);
		GMRFLib_idxval_add(&h, 24919, 1);
		GMRFLib_idxval_add(&h, 24920, 1);
		GMRFLib_idxval_add(&h, 24925, 1);
		GMRFLib_idxval_add(&h, 24932, 1);
		GMRFLib_idxval_add(&h, 24938, 1);
		GMRFLib_idxval_add(&h, 24944, 1);
		GMRFLib_idxval_add(&h, 24950, 1);
		GMRFLib_idxval_add(&h, 24955, 1);
		GMRFLib_idxval_add(&h, 24959, 1);
		GMRFLib_idxval_add(&h, 24962, 1);
		GMRFLib_idxval_add(&h, 24968, 1);
		GMRFLib_idxval_add(&h, 24971, 1);
		GMRFLib_idxval_add(&h, 24977, 1);
		GMRFLib_idxval_add(&h, 24983, 1);
		GMRFLib_idxval_add(&h, 24989, 1);
		GMRFLib_idxval_add(&h, 24992, 1);
		GMRFLib_idxval_add(&h, 24998, 1);
		GMRFLib_idxval_add(&h, 25000, 1);
		GMRFLib_idxval_add(&h, 25006, 1);
		GMRFLib_idxval_add(&h, 25012, 1);
		GMRFLib_idxval_add(&h, 25013, 1);
		GMRFLib_idxval_add(&h, 25019, 1);
		GMRFLib_idxval_add(&h, 25025, 1);
		GMRFLib_idxval_add(&h, 25030, 1);
		GMRFLib_idxval_add(&h, 25036, 0);
		GMRFLib_idxval_add(&h, 25037, 0);
		GMRFLib_idxval_add(&h, 25038, 0);
		GMRFLib_idxval_add(&h, 25039, 0);
		GMRFLib_idxval_add(&h, 25040, 0);
		GMRFLib_idxval_add(&h, 25042, 1);
		GMRFLib_idxval_add(&h, 25048, 1);
		GMRFLib_idxval_add(&h, 25054, 1);
		GMRFLib_idxval_add(&h, 25060, 1);
		GMRFLib_idxval_add(&h, 25065, 1);
		GMRFLib_idxval_add(&h, 25070, 1);
		GMRFLib_idxval_add(&h, 25075, 1);

		GMRFLib_idxval_prepare(&h, 1, 1);
	}
		break;

	case 82:
	{
		int n = atoi(args[0]);
		int ntimes = atoi(args[1]);
		double *x = Calloc(n, double);
		for (int i = 0; i < n; i++) {
			x[i] = GMRFLib_uniform();
			// x[i] = i+1;
		}

		P(n);
		P(ntimes);

		double tref[2] = { 0.0, 0.0 };
		double r = 0.0, rr = 0.0;

		for (int time = 0; time < ntimes; time++) {

			tref[0] -= GMRFLib_cpu();
			r += GMRFLib_dsum(n, x);
			tref[0] += GMRFLib_cpu();

			tref[1] -= GMRFLib_cpu();
#pragma omp simd reduction(+: rr)
			for (int i = 0; i < n; i++) {
				rr += x[i];
			}
			tref[1] += GMRFLib_cpu();
		}

		printf("dsum %.3f plain %.3f (r-rr=%.12f, %1d)\n", tref[0] / (tref[0] + tref[1]), tref[1] / (tref[0] + tref[1]), r - rr, r == rr);

		Free(x);
	}
		break;

	case 83:
	{
		FIXME("????????????????? FREE in idxval.c needs to disabled for this to run");
		int n = atoi(args[0]);
		int ntimes = atoi(args[1]);
		double *xx = Calloc(n, double);

		for (int i = 0; i < n; i++) {
			xx[i] = GMRFLib_uniform();
		}

		GMRFLib_idxval_tp *h = NULL;
		for (int i = 0, j = 0; i < n; i++) {
			j += 1 + (GMRFLib_uniform() < 0.9 ? 0.0 : 1 + (int) (GMRFLib_uniform() * 31));
			if (j >= n)
				break;
			GMRFLib_idxval_add(&h, j, xx[j]);
		}
		GMRFLib_idxval_prepare(&h, 1, 1);
		assert(h);
		P(n);
		P(h->g_n);
		P(h->n / h->g_n);

		double sum1 = 0.0, sum2 = 0.0;
		double tref1 = 0.0, tref2 = 0.0;
		for (int k = 0; k < ntimes; k++) {
			sum1 = sum2 = 0.0;
			tref1 -= GMRFLib_cpu();
			sum1 = GMRFLib_dot_product_serial(h, xx);
			tref1 += GMRFLib_cpu();

			tref2 -= GMRFLib_cpu();
			sum2 = GMRFLib_dot_product_group(h, xx);
			tref2 += GMRFLib_cpu();
			if (ABS(sum1 - sum2) > 1e-8) {
				P(sum1);
				P(sum2);
				exit(88);
			}
		}
		printf("serial %.3f group %.3f (%.3f, %.3f)\n", tref1, tref2, tref1 / (tref1 + tref2), tref2 / (tref1 + tref2));
		Free(xx);
	}
		break;

	case 84:
	{
		FIXME("FREE in idxval.c needs to disabled for this to run");
		int n = atoi(args[0]);
		int m = atoi(args[1]);
		double *xx = Calloc(n, double);

		for (int i = 0; i < n; i++) {
			xx[i] = GMRFLib_uniform();
		}

		GMRFLib_idxval_tp *h = NULL;
		for (int i = 0, j = 0; i < n; i++) {
			j += 1 + (GMRFLib_uniform() < 0.9 ? 0.0 : 1 + (int) (GMRFLib_uniform() * 31));
			if (j >= n)
				break;
			GMRFLib_idxval_add(&h, j, xx[j]);
		}
		GMRFLib_idxval_prepare(&h, 1, 1);
		P(n);
		P(m);
		P(h->g_n);
		P(h->n / h->g_n);

		double sum1 = 0.0, sum2 = 0.0;
		double tref1 = 0.0, tref2 = 0.0;
		for (int k = 0; k < m; k++) {
			sum1 = sum2 = 0.0;
			tref1 -= GMRFLib_cpu();
			sum1 = GMRFLib_ddot_idx(h->n, h->val, xx, h->idx);
			tref1 += GMRFLib_cpu();

			tref2 -= GMRFLib_cpu();
			for (int i = 0; i < h->n; i++) {
				sum2 += h->val[i] * xx[h->idx[i]];
			}
			tref2 += GMRFLib_cpu();
			if (ABS(sum1 - sum2) > 1e-8) {
				P(sum1);
				P(sum2);
				exit(88);
			}
		}
		printf("dot_idx %.3f serial %.3f (%.3f, %.3f)\n", tref1, tref2, tref1 / (tref1 + tref2), tref2 / (tref1 + tref2));
		Free(xx);
	}
		break;

	case 85:
	{
		int n = atoi(args[0]);
		int m = atoi(args[1]);

		double time = -GMRFLib_cpu();
		double sum = 0.0;
		int imin, imax;
		for (int k = 0; k < m; k++) {
			for (int i = 0; i < n; i++) {
				int ii = GMRFLib_uniform() * n;
				for (int j = 0; j < n; j++) {
					int jj = GMRFLib_uniform() * n;
					imin = IMIN(ii, jj);
					imax = IMAX(ii, jj);
					sum += imin - imax;
				}
			}
		}
		time += GMRFLib_cpu();
		printf("IMAX/MIN %.12f\n", time);

		sum = 0.0;
		time = -GMRFLib_cpu();
		for (int k = 0; k < m; k++) {
			for (int i = 0; i < n; i++) {
				int ii = GMRFLib_uniform() * n;
				for (int j = 0; j < n; j++) {
					int jj = GMRFLib_uniform() * n;
					if (ii <= jj) {
						imin = ii;
						imax = jj;
					} else {
						imin = jj;
						imax = ii;
					}
					sum += imin - imax;
				}
			}
		}
		time += GMRFLib_cpu();
		printf("if/ %.12f\n", time);
	}
		break;

	case 86:
	{
		double x = 0.0;
		double param[] = { 1, 0.001 };
		for (x = 0;; x++) {
			printf("x %f ldens %f\n", x, priorfunc_loggamma(&x, param));
		}
	}
		break;

	case 87:
	{
		int n = atoi(args[0]);
		int ntimes = atoi(args[1]);
		double *x = Calloc(n, double);
		for (int i = 0; i < n; i++) {
			x[i] = GMRFLib_uniform();
		}

		P(n);
		P(ntimes);

		double tref[2] = { 0.0, 0.0 };
		double r = 0.0, rr = 0.0;

		for (int time = 0; time < ntimes; time++) {

			tref[0] -= GMRFLib_cpu();
			r += GMRFLib_dsum(n, x);
			tref[0] += GMRFLib_cpu();

			tref[1] -= GMRFLib_cpu();
			rr += GMRFLib_dsum2(n, x);
			tref[1] += GMRFLib_cpu();
		}

		printf("dsum %.3f dsum2 %.3f\n", tref[0] / (tref[0] + tref[1]), tref[1] / (tref[0] + tref[1]));
		P((r - rr) / r);

		Free(x);
	}
		break;

	case 88:
	{
		int n = atoi(args[0]);
		int m = atoi(args[1]);
		double *xx = Calloc(n, double);

		for (int i = 0; i < n; i++) {
			xx[i] = GMRFLib_uniform();
		}

		P(n);
		P(m);

		double tref1 = 0.0, tref2 = 0.0;
		for (int k = 0; k < m; k++) {
			tref1 -= GMRFLib_cpu();
			double s = GMRFLib_uniform();
#pragma omp simd reduction(+: s)
			for (int i = 0; i < n; i++) {
				s += xx[i];
			}
			tref1 += GMRFLib_cpu();
			tref2 -= GMRFLib_cpu();
#pragma GCC ivdep
			for (int i = 0; i < n; i++) {
				s += xx[i];
			}
			tref2 += GMRFLib_cpu();
		}
		printf("simd %.3f opt %.3f (%.3f, %.3f)\n", tref1, tref2, tref1 / (tref1 + tref2), tref2 / (tref1 + tref2));
		Free(xx);
	}
		break;

	case 89:
	{
		int n = atoi(args[0]);
		int m = atoi(args[1]);
		P(n);
		P(m);
		GMRFLib_idxval_tp *h = NULL;
		double *xx = Calloc(n + 1, double);
		for (int i = 0; i < n + 1; i++) {
			xx[i] = GMRFLib_uniform();
		}

		for (int i = 0, j = 0; i < ISQR(n); i++) {
			j += 1 + (GMRFLib_uniform() < 0.8 ? 0 : 1 + (int) (GMRFLib_uniform() * 63));
			GMRFLib_idxval_add(&h, j, GMRFLib_uniform());
			if (h->n >= n)
				break;
		}
		GMRFLib_idxval_nsort_x(&h, 1, 1, 0, 0);
		P(h->n);
		double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, sum4 = 0.0;
		double tref1 = 0.0, tref2 = 0.0, tref3 = 0.0, tref4 = 0.0;
		for (int k = 0; k < m; k++) {

			sum1 = sum2 = sum3 = sum4 = 0.0;
			tref1 -= GMRFLib_cpu();
			sum1 = GMRFLib_ddot_idx(h->n, h->val, xx, h->idx);
			tref1 += GMRFLib_cpu();

			tref2 -= GMRFLib_cpu();
			sum2 = GMRFLib_ddot_idx_mkl_OLD(h->n, h->val, xx, h->idx);
			tref2 += GMRFLib_cpu();

			tref3 -= GMRFLib_cpu();
			sum3 = GMRFLib_ddot_idx_mkl_NEW(h->n, h->val, xx, h->idx);
			tref3 += GMRFLib_cpu();

			tref4 -= GMRFLib_cpu();
			sum4 = GMRFLib_ddot_idx_mkl(h->n, h->val, xx, h->idx);
			tref4 += GMRFLib_cpu();

			if (ABS(sum1 - sum2) > 1e-8 || ABS(sum1 - sum3) > 1e-8 || ABS(sum1 - sum4) > 1e-8) {
				P(sum1);
				P(sum2);
				P(sum3);
				P(sum4);
				exit(88);
			}
		}
		printf("dot_idx %.3f mkl_OLD %.3f mkl_NEW %.3f mkl %.3f (%.3f, %.3f, %.3f, %.3f)\n",
		       tref1, tref2, tref3, tref4,
		       tref1 / (tref1 + tref2 + tref3 + tref4),
		       tref2 / (tref1 + tref2 + tref3 + tref4), tref3 / (tref1 + tref2 + tref3 + tref4), tref4 / (tref1 + tref2 + tref3 + tref4));
		Free(xx);
	}
		break;

	case 90:
	{
		int n = atoi(args[0]);
		int m = atoi(args[1]);
		double *xx = Calloc(n, double);
		for (int i = 0; i < n; i++) {
			xx[i] = GMRFLib_uniform();
		}

		P(n);
		P(m);
		double sum1 = 0.0, sum2 = 0.0;
		double tref1 = 0.0, tref2 = 0.0;
		for (int k = 0; k < m; k++) {
			sum1 = sum2 = 0.0;
			tref1 -= GMRFLib_cpu();
			sum1 = GMRFLib_dsum(n, xx);
			tref1 += GMRFLib_cpu();

			tref2 -= GMRFLib_cpu();
			sum2 = 0.0;
			for (int i = 0; i < n; i++) {
				sum2 += xx[i];
			}
			tref2 += GMRFLib_cpu();
			if (ABS(sum1 - sum2) > 1e-8) {
				P(sum1);
				P(sum2);
				exit(88);
			}
		}
		printf("_dsum %.3f loop %.3f (%.3f, %.3f)\n", tref1, tref2, tref1 / (tref1 + tref2), tref2 / (tref1 + tref2));
		Free(xx);
	}
		break;

	case 91:
	{
		void *monobound_bsearch(const void *key, const void *array, size_t nmemb, size_t size);

		int n = atoi(args[0]);
		int m = atoi(args[1]);
		P(n);
		P(m);
		double start = 0, finish = 0;
		double start2 = 0, finish2 = 0;
		int *idx = Calloc(n, int);
		int key = 1;
		for (int i = 0; i < n; i++) {
			key += i;
			// printf("%d %d \n", i, key);
			idx[i] = key;
		}
		for (int k = 0; k < m; k++) {
			double sum = 0.0;
			start += omp_get_wtime();
			int low = 0;
			for (int i = 0; i < key + 1; i++) {
				int p = GMRFLib_iwhich_sorted_g(i, idx, n, &low);
				if (p >= 0)
					sum += idx[p];
			}
			finish += omp_get_wtime();

			double sum2 = 0.0;
			start2 += omp_get_wtime();
			int guess[2] = { 0, 0 };
			for (int i = 0; i < key + 1; i++) {
				int p = GMRFLib_iwhich_sorted_g2(i, idx, n, guess);
				if (p >= 0)
					sum2 += idx[p];
			}
			finish2 += omp_get_wtime();
			if (k == m - 1)
				printf("n.lookups= %1d  Time for iwhich_g= %.4g iwhich_g2= %.4g ratio g/g2= %.4f\n",
				       key, (finish - start) / (k + 1.0), (finish2 - start2) / (k + 1.0), (finish - start) / (finish2 - start2));
		}
	}
		break;

	case 92:
	{
		double a = 2.0;
		double n = 1.0;
		while (1) {
			n *= 10.0;
			a *= log(n);
			printf("a = %g n = %g log(Beta(a, a/n) = %g\n", a, n, my_gsl_sf_lnbeta(a, a / n));
		}
	}
		break;

	case 93:
	{
		int n = atoi(args[0]);
		int m = atoi(args[1]);
		double *a = Calloc(n, double);
		for (int i = 0; i < n; i++) {
			a[i] = GMRFLib_uniform();
		}

		P(n);
		P(m);
		double sum1 = 0.0, sum2 = 0.0;
		double start = 0, start2 = 0, finish = 0, finish2 = 0;

		for (int k = 0; k < m; k++) {
			start += omp_get_wtime();
			for (int i = 0; i < n; i++) {
				double aa = a[i];
				for (int j = 0; j < i; j++) {
					sum1 += log(aa + j);
				}
			}
			finish += omp_get_wtime();

			start2 += omp_get_wtime();
			for (int i = 0; i < n; i++) {
				sum2 += my_betabinomial_helper(i, a[i]);
			}
			finish2 += omp_get_wtime();
		}
		printf("plain= %.4g helper= %.4g ratio plain/helper= %.4f\n",
		       (finish - start) / m, (finish2 - start2) / m, (finish - start) / (finish2 - start2));
		P((sum1 - sum2) / (sum1 + sum2));
	}
		break;

	case 94:
	{
		int n = atoi(args[0]);
		int nt = atoi(args[1]);
		P(n);
		P(nt);
		double *x = Calloc(n, double);
		for (int i = 0; i < n; i++) {
			x[i] = GMRFLib_uniform();
		}

		double tref;
		tref = -GMRFLib_cpu();
#pragma omp parallel for num_threads(nt)
		for (int i = 0; i < nt; i++) {
			char *fnm;
			GMRFLib_sprintf(&fnm, "REMOVE_ME_%1d.dat", i);
			FILE *fp = fopen(fnm, "wb");
			fwrite((void *) x, sizeof(double), (size_t) n, fp);
			fclose(fp);
			Free(fnm);
		}
		tref += GMRFLib_cpu();
		P(tref);

		double tref2;
		tref2 = -GMRFLib_cpu();
		for (int i = 0; i < nt; i++) {
			char *fnm;
			GMRFLib_sprintf(&fnm, "REMOVE_ME_ALSO_%1d.dat", i);
			FILE *fp = fopen(fnm, "wb");
			fwrite((void *) x, sizeof(double), (size_t) n, fp);
			fclose(fp);
			Free(fnm);
		}
		tref2 += GMRFLib_cpu();
		P(tref2);
		P(tref2 / tref);
		Free(x);
	}
		break;

	case 95:
	{
		int n = atoi(args[0]);
		int m = atoi(args[1]);
		P(n);
		P(m);
		int *ix = Calloc(n, int);
		double *x = Calloc(n, double);

		double tref = 0.0;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				ix[j] = (int) ((2 * n) * GMRFLib_uniform());
				x[j] = GMRFLib_uniform();
			}

			tref -= GMRFLib_cpu();
			my_insertionSort_id(ix, x, n);
			tref += GMRFLib_cpu();

			int errors = 0;
			for (int j = 1; j < n; j++) {
				errors += (ix[j] < ix[j - 1]);
			}
			assert(errors == 0);
		}

		double tref2 = 0.0;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				ix[j] = (int) ((2 * n) * GMRFLib_uniform());
				x[j] = GMRFLib_uniform();
			}

			tref2 -= GMRFLib_cpu();
			gsl_sort2_id(ix, x, n);
			tref2 += GMRFLib_cpu();

			int errors = 0;
			for (int j = 1; j < n; j++) {
				errors += (ix[j] < ix[j - 1]);
			}
			assert(errors == 0);
		}

		printf("insert %g sort2 %g  insert/sort2 =  %g\n", tref, tref2, tref / tref2);
	}
		break;

	case 96:
	{
		int n = atoi(args[0]);
		int m = atoi(args[1]);
		P(n);
		P(m);
		double *x = Calloc(n, double);
		double *xx = Calloc(n, double);

		double tref = 0.0;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				xx[j] = GMRFLib_uniform();
				x[j] = GMRFLib_uniform();
			}

			tref -= GMRFLib_cpu();
			my_insertionSort_dd(xx, x, n);
			tref += GMRFLib_cpu();

			int errors = 0;
			for (int j = 1; j < n; j++) {
				errors += (xx[j] < xx[j - 1]);
			}
			assert(errors == 0);
		}

		double tref2 = 0.0;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				xx[j] = GMRFLib_uniform();
				x[j] = GMRFLib_uniform();
			}

			tref2 -= GMRFLib_cpu();
			gsl_sort2_dd(xx, x, n);
			tref2 += GMRFLib_cpu();

			int errors = 0;
			for (int j = 1; j < n; j++) {
				errors += (xx[j] < xx[j - 1]);
			}
			assert(errors == 0);
		}

		printf("insert %g sort2 %g  insert/sort2 =  %g\n", tref, tref2, tref / tref2);
	}
		break;

	case 97:
	{
		GMRFLib_idxval_tp *v = NULL;
		int idx[] = { 0, 1, 2, 4, 5, 6, 7, 8, 10, 11, 13, 15, 17, 18, 19, 21, 22, 24, 26, 27, 30 };
		// int idx[] = {0, 2, 3, 4, 5, 6, 7, 8, 10}; 

		for (int i = 0; i < (int) (sizeof(idx) / sizeof(int)); i++) {
			// GMRFLib_idxval_add(&v, idx[i], (double) 1.0);
			GMRFLib_idxval_add(&v, idx[i], 1.0 + (double) i + 0 * GMRFLib_uniform());
		}

		GMRFLib_idxval_prepare(&v, 1, 1);
	}
		break;

	case 98:
	{
		int n = atoi(args[0]);
		int ntimes = atoi(args[1]);
		double *xx = Calloc(n, double);

		// GMRFLib_rng_init(1);
		for (int i = 0; i < n; i++) {
			xx[i] = GMRFLib_uniform();
		}
		P(xx[0]);

		GMRFLib_idxval_tp *h = NULL;
		GMRFLib_idxval_tp *hh = NULL;
		for (int i = 0, j = 0; i < n; i++) {
			j += 1 + (GMRFLib_uniform() < 0.95 ? 0.0 : 1 + (int) (GMRFLib_uniform() * 31));
			if (j >= n)
				break;
			GMRFLib_idxval_add(&h, j, xx[j]);
			GMRFLib_idxval_add(&hh, j, xx[j]);
			// GMRFLib_idxval_add(&h, j, 1.0);
			// GMRFLib_idxval_add(&hh, j, 1.0);
		}
		if (!h || !hh)
			exit(0);

		double tref1 = 0.0, tref2 = 0.0;
		tref1 -= GMRFLib_cpu();
		GMRFLib_idxval_prepare(&h, 1, 1);
		tref1 += GMRFLib_cpu();
		tref2 -= GMRFLib_cpu();
		GMRFLib_idxval_prepare(&hh, 1, 1);
		tref2 += GMRFLib_cpu();

		P(tref1 / (tref1 + tref2));

		P(n);
		P(ntimes);

		double sum1 = 0.0, sum2 = 0.0;
		tref1 = 0.0;
		tref2 = 0.0;
		for (int k = 0; k < ntimes; k++) {
			sum1 = sum2 = 0.0;
			tref1 -= GMRFLib_cpu();
			sum1 = GMRFLib_dot_product(h, xx);
			tref1 += GMRFLib_cpu();

			tref2 -= GMRFLib_cpu();
			sum2 = GMRFLib_dot_product(hh, xx);
			tref2 += GMRFLib_cpu();
			if (ABS(sum1 - sum2) > 1e-8) {
				P(sum1);
				P(sum2);
				exit(88);
			}
		}
		GMRFLib_idxval_free(h);
		GMRFLib_idxval_free(hh);
		printf("new %.3f old %.3f (%.3f, %.3f)\n", tref1, tref2, tref1 / (tref1 + tref2), tref2 / (tref1 + tref2));
		Free(xx);
	}
		break;

	case 99:
	{
		int n = atoi(args[0]);
		int m = atoi(args[1]);
		int nm = n * m;

		double *x = Calloc(nm, double);
		double *xx = Calloc(nm, double);
		for (int i = 0; i < nm; i++) {
			x[i] = xx[i] = GMRFLib_uniform();
		}
		double tref1 = 0.0, tref2 = 0.0;
		for (int nt = 0; nt < 100; nt++) {
			tref1 -= GMRFLib_cpu();
			for (int i = 0; i < n; i++) {
				int offset = (n - i - 1) * m;
				Memcpy(x + offset, xx + offset, m * sizeof(double));
			}
			tref1 += GMRFLib_cpu();
			tref2 -= GMRFLib_cpu();
			for (int i = 0; i < n; i++) {
				int offset = (n - i - 1) * m;
				for (int j = 0; j < m; j++) {
					x[offset + j] = xx[offset + j];
				}
			}
			tref2 += GMRFLib_cpu();
		}
		printf("memcpy %.3f plain %.3f (%.3f, %.3f)\n", tref1, tref2, tref1 / (tref1 + tref2), tref2 / (tref1 + tref2));
	}
		break;

	case 100:
	{
		int n = atoi(args[0]);
		double *x = Calloc(n, double);
		for (int i = 0; i < n; i++) {
			x[i] = GMRFLib_uniform();
		}
		printf("%s\n", GMRFLib_vec2char(x, n));
	}
		break;

	case 101:
	{
		P(log(0.0));
		P(sqrt(-1.0));
	}
		break;

	case 102:
	{
		double x = GMRFLib_uniform(), y;
		P(x);
		y = map_invcloglog(x, MAP_FORWARD, NULL);
		P(y);
		y = map_invcloglog(y, MAP_BACKWARD, NULL);
		P(x - y);

		P(x);
		y = map_invccloglog(x, MAP_FORWARD, NULL);
		P(y);
		y = map_invccloglog(y, MAP_BACKWARD, NULL);
		P(x - y);

		double h = 1.0e-4;
		y = (map_invcloglog(x + h, MAP_FORWARD, NULL) - map_invcloglog(x - h, MAP_FORWARD, NULL)) / 2.0 / h;
		P(map_invcloglog(x, MAP_DFORWARD, NULL) - y);
		y = (map_invccloglog(x + h, MAP_FORWARD, NULL) - map_invccloglog(x - h, MAP_FORWARD, NULL)) / 2.0 / h;
		P(map_invccloglog(x, MAP_DFORWARD, NULL) - y);

	}
		break;

	case 103:
	{
		char a;
		unsigned char b;
		signed char c;

		a = b = c = 1;
		printf("char %d unsigned %d signed %d\n", (int) a, (int) b, (int) c);
		a = b = c = -1;
		printf("char %d unsigned %d signed %d\n", (int) a, (int) b, (int) c);
	}
		break;

	case 104:
	{
		P(sqrt(-1.0));
		P(log(-1.0));
		P(1.0 / 0.0);
		P(exp(exp(100.0)));
	}
		break;

	case 105:
	{
		P(ISZERO(0));
		P(ISZERO(1));
		P(ISZERO(0.0));
		P(ISZERO(1.0));
		P(ISZERO(0.0f));
		P(ISZERO(1.0f));
	}
		break;

	case 106:
	{
		printf("## check with\nlibrary(VGAM)\nfor(i in 0:200) print(c(i,  log(bell(i)) - lfactorial(i)))\n\n");
#pragma omp parallel for
		for (int i = -1; i <= 2057; i++) {
#pragma omp critical (Name_69969525cb4835f6178baf1c8599321a9419c0a8)
			{
				printf("%d %.20g\n", i, my_lbell(i));
			}
		}
	}
		break;

	case 107:
	{
		double tref[2] = { 0, 0 };
		int n = atoi(args[0]);
		double *y = Calloc(2 * n, double);
		double *res = y + n;

		double rel_err = 0.0;
		for (int i = 0; i < n; i++) {
			y[i] = exp(5.0 * GMRFLib_stdnormal());
			double ref = gsl_sf_lambert_W0(y[i]);
			rel_err += ABS((ref - my_lambert_W0(y[i])) / ref);
		}
		P(rel_err / n);

		tref[0] = -GMRFLib_cpu();
		double sum = 0.0;
		for (int i = 0; i < n; i++) {
			sum += gsl_sf_lambert_W0(y[i]);
		}
		tref[0] += GMRFLib_cpu();

		tref[1] = -GMRFLib_cpu();
		my_lambert_W0s(n, y, res);
		tref[1] += GMRFLib_cpu();
		printf("random arguments: GSL:  %.4f  Cache:  %.4f\n", tref[0] / (tref[0] + tref[1]), tref[1] / (tref[0] + tref[1]));

		qsort((void *) y, (size_t) n, sizeof(double), GMRFLib_dcmp);
		tref[0] = -GMRFLib_cpu();
		for (int i = 0; i < n; i++) {
			sum += gsl_sf_lambert_W0(y[i]);
		}
		tref[0] += GMRFLib_cpu();

		tref[1] = -GMRFLib_cpu();
		my_lambert_W0s(n, y, res);
		tref[1] += GMRFLib_cpu();
		printf("sorted arguments: GSL:  %.4f  Cache:  %.4f\n", tref[0] / (tref[0] + tref[1]), tref[1] / (tref[0] + tref[1]));
		P(sum);
	}
		break;

	case 108:
	{
		int n = atoi(args[0]);
		my_lbell(n);
		for (int y = 0; y <= n; y++) {
			printf("y %d my_lbell %.12g\n", y, my_lbell(y));
		}
	}
		break;

	case 109:
	{
		priorfunc_fgn_priorH_extract();
	}
		break;

	case 110:
	{
		double tref[3] = { 0, 0, 0 };
		int n = atoi(args[0]);
		double *y = Calloc(2 * n, double);

		double rel_err = 0.0;
		for (int i = 0; i < n; i++) {
			y[i] = exp(2 * GMRFLib_stdnormal());
			double ref = gsl_sf_lngamma(y[i]);
			rel_err += ABS((ref - lgamma(y[i])));
		}
		P(rel_err / n);

		tref[0] = -GMRFLib_cpu();
		double sum = 0.0;
		for (int i = 0; i < n; i++) {
			sum += gsl_sf_lngamma(y[i]);
		}
		tref[0] += GMRFLib_cpu();
		P(sum);

		sum = 0.0;
		tref[1] = -GMRFLib_cpu();
		for (int i = 0; i < n; i++) {
			sum += lgamma(y[i]);
		}
		tref[1] += GMRFLib_cpu();
		P(sum);

		printf("GSL:  %.4f  libm: %.4f NULL:  %.4f\n", tref[0] / (tref[0] + tref[1] + tref[2]),
		       tref[1] / (tref[0] + tref[1] + tref[2]), tref[2] / (tref[0] + tref[1] + tref[2]));

		P(sum);
	}
		break;


	case 111:
	{
		int n = atoi(args[0]);
		int m = atoi(args[1]);
		P(n);
		P(m);
		double *x = Calloc(n, double);
		double *xx = Calloc(n, double);
		for (int i = 0; i < n; i++) {
			x[i] = xx[i] = GMRFLib_uniform();
		}

		double tref[] = { 0, 0 };
		for (int i = 0; i < m; i++) {
			tref[0] -= GMRFLib_cpu();
			qsort(x, (size_t) n, sizeof(double), GMRFLib_dcmp);
			tref[0] += GMRFLib_cpu();
			Memcpy(x, xx, n * sizeof(double));

			tref[1] -= GMRFLib_cpu();
			qsort(x, (size_t) n, sizeof(double), GMRFLib_dcmp);
			tref[1] += GMRFLib_cpu();
			Memcpy(x, xx, n * sizeof(double));
		}
		printf("sorted arguments: qsort:  %.4f  mkl:  %.4f\n", tref[0] / (tref[0] + tref[1]), tref[1] / (tref[0] + tref[1]));
	}
		break;

	case 112:
	{
		int n = atoi(args[0]);
		double *x = Calloc(2 * n, double);
		double *y = x + n;

		for (int i = 0; i < n; i++) {
			x[i] = i;
			y[i] = sin(x[i] / n * 6.0 * M_PI);
		}

		GMRFLib_spline_tp *s = GMRFLib_spline_create(x, y, n);

		for (int i = 0; i < n; i++) {
			printf("X:  %g %g\n", x[i], y[i]);
		}
		printf("X:\n");
		for (double xx = -n / 4.0; xx < n + n / 4.0; xx += 0.01) {
			printf("X:  %g %g\n", xx, GMRFLib_spline_eval(xx, s));
		}
	}
		break;

	case 113:
	{
		int n = atoi(args[0]);
		double *x = Calloc(2 * n, double);
		double *y = x + n;

		for (int i = 0; i < n; i++) {
			x[i] = (i > 0 ? x[i - 1] + 2.0 * GMRFLib_uniform() : 0);
			y[i] = i * (GSL_IS_ODD(i) ? 1.0 : -1.0);
		}

		GMRFLib_spline_tp *s = GMRFLib_spline_create(x, y, n);

		for (int i = 0; i < n; i++) {
			printf("X:  %g %g\n", x[i], y[i]);
		}
		printf("X:\n");
		for (double xx = -5; xx < n + 1; xx += 0.01) {
			printf("X:  %.12g %.12g\n", xx, GMRFLib_spline_eval(xx, s));
		}
	}
		break;

	case 114:
	{
		int n = atoi(args[0]);
		int m = atoi(args[1]);
		P(n);
		P(m);
		double *x = Calloc(3 * n, double);
		double *y = x + n;
		double *yy = x + 2 * n;

		for (int i = 0; i < n; i++) {
			x[i] = GMRFLib_uniform();
		}

		double tref[] = { 0, 0 };
		for (int i = 0; i < m; i++) {
			double a = GMRFLib_uniform();
			double b = GMRFLib_uniform();

			tref[0] -= GMRFLib_cpu();
			GMRFLib_daxpb(n, a, x, b, y);
			tref[0] += GMRFLib_cpu();

			tref[1] -= GMRFLib_cpu();
#pragma omp simd
			for (int j = 0; j < n; j++) {
				yy[j] = a * x[j] + b;
			}
			tref[1] += GMRFLib_cpu();

			double err = 0.0;
			for (int j = 0; j < n; j++) {
				err = DMAX(err, ABS(y[j] - yy[j]));
			}
			assert(err < FLT_EPSILON);
		}
		printf("mod:  %.4f  plain:  %.4f\n", tref[0] / (tref[0] + tref[1]), tref[1] / (tref[0] + tref[1]));
	}
		break;

	case 115:
	{
		P(GMRFLib_CACHE_LEN());
		P(sizeof(gsl_interp_accel *));
	}
		break;

	case 116:
	{
		int n = atoi(args[0]);
		double *x = Calloc(2 * n, double);
		double *y = x + n;

		for (int i = 0; i < n; i++) {
			x[i] = (i > 0 ? x[i - 1] + 2.0 * GMRFLib_uniform() : 0);
			y[i] = i * (GSL_IS_ODD(i) ? 1.0 : -1.0);
		}

		GMRFLib_spline_tp *s = GMRFLib_spline_create(x, y, n);

		for (int i = 0; i < n; i++) {
			printf("X:  %g %g\n", x[i], y[i]);
		}
		printf("X:\n");
#pragma omp parallel for
		for (int i = 0; i < 100; i++) {
			int t = omp_get_thread_num();
#pragma omp parallel for
			for (int ii = 0; ii < 100; ii++) {
				int tt = omp_get_thread_num();
				double xx = -5 + ii * 0.1;
				printf("(%1d %1d):  %.12g %.12g\n", t, tt, xx, GMRFLib_spline_eval(xx, s));
			}
		}
	}
		break;

	case 117:
	{
		GMRFLib_sn_param_tp p;
		p.xi = GMRFLib_uniform();
		p.omega = 1 + exp(GMRFLib_uniform() - 0.5);
		p.alpha = 3.0 * (GMRFLib_uniform() - 0.5);

		int n = 100;
		double mom[3];
		GMRFLib_sn_par2moments(mom, mom + 1, mom + 2, &p);

		double xlow = mom[0] - 5 * mom[1];
		double xhigh = mom[0] + 5 * mom[1];
		double dx = mom[1] / 5.0;

		n = (xhigh - xlow) / dx + 2;

		double *x = Calloc(n, double);
		double *ld = Calloc(n, double);

		int k = 0;
		for (double xx = xlow; xx <= xhigh; xx += dx) {
			x[k] = xx;
			double xs = (xx - p.xi) / p.omega;
			ld[k] = -0.5 * SQR(xs) + inla_log_Phi(p.alpha * xs);
			k++;
		}
		n = k;

		double std_mean = GMRFLib_uniform();
		double std_sd = GMRFLib_uniform();

		printf("std mean sd %g %g\n", std_mean, std_sd);
		GMRFLib_density_tp *d1 = NULL, *d2 = NULL;
		GMRFLib_density_create_sn(&d1, p, std_mean, std_sd, 1);
		GMRFLib_density_create(&d2, GMRFLib_DENSITY_TYPE_SCGAUSSIAN, n, x, ld, std_mean, std_sd, 1);

		for (int i = 0; i < 10; i++) {
			double z = xlow + GMRFLib_uniform() * (xhigh - xlow);
			double pp = GMRFLib_uniform();
			printf("z %g pp %g\n", z, pp);

			double dd1, dd2;
			GMRFLib_evaluate_logdensity(&dd1, z, d1);
			GMRFLib_evaluate_logdensity(&dd2, z, d2);
			printf("\tld %g %g\n", dd1, dd2);

			double q1, q2;
			GMRFLib_density_Pinv(&q1, pp, d1);
			GMRFLib_density_Pinv(&q2, pp, d2);
			printf("\tq  %g %g\n", q1, q2);
		}
	}
		break;

	case 118:
	{
		typedef double fun_tp(double arg);

		typedef struct {
			fun_tp *fun[2];
			char *name;
		} cmp_tp;

		cmp_tp cmp[] = {
			{ { GMRFLib_erf, gsl_sf_erf}, "erf" },
			{ { GMRFLib_erfc, gsl_sf_erfc}, "erfc" },
			{ { GMRFLib_cdfnorm, gsl_cdf_ugaussian_P}, "cdfnorm" },
			{ { GMRFLib_cdfnorm_inv, gsl_cdf_ugaussian_Pinv}, "cdfnorm_inv" },
		};

		int n = atoi(args[0]);
		double *x = Calloc(3 * n, double);
		double *y[] = { x + n, x + 2 * n };

		for (int i = 0; i < n; i++) {
			x[i] = GMRFLib_uniform();
		}

		int K = sizeof(cmp) / sizeof(cmp_tp);
		for (int k = 0; k < K; k++) {
			double tref[] = { 0, 0 };
			for (int j = 0; j < 2; j++) {
				tref[j] -= GMRFLib_cpu();
				double *yy = y[j];
				for (int i = 0; i < n; i++) {
					yy[i] = cmp[k].fun[j] (x[i]);
				}
				tref[j] += GMRFLib_cpu();
			}
			double max_err = 0.0;
			for (int i = 0; i < n; i++) {
				max_err = DMAX(max_err, y[0][i] - y[1][i]);
			}
			printf("%s %.4f %.4f (max.err %g)\n", cmp[k].name, tref[0] / (tref[0] + tref[1]), tref[1] / (tref[0] + tref[1]), max_err);
		}
	}
		break;

	case 119:
	{
		double p = GMRFLib_uniform();
		P(p - GMRFLib_erf(GMRFLib_erf_inv(p)));
		P(p - GMRFLib_erfc(GMRFLib_erfc_inv(p)));
		P(p - GMRFLib_cdfnorm(GMRFLib_cdfnorm_inv(p)));

		for (double x = -20.0; x < 20.0; x += 0.5) {
			p = 1.0 / (1.0 + exp(-x));
			printf("%g gsl %g cdfnorm_inv %g\n", p, gsl_cdf_ugaussian_Pinv(p), GMRFLib_cdfnorm_inv(p));
		}
	}
		break;

	case 120:
	{
		int n = atoi(args[0]);
		int m = atoi(args[1]);
		P(n);
		P(m);
		double *x = Calloc(3 * n, double);
		double *y = x + n;
		double *yy = x + 2 * n;

		for (int i = 0; i < n; i++) {
			x[i] = GMRFLib_uniform();
			y[i] = GMRFLib_uniform();
			yy[i] = y[i];
		}

		double tref[] = { 0, 0 };
		for (int i = 0; i < m; i++) {
			tref[0] -= GMRFLib_cpu();
			GMRFLib_daddto(n, x, y);	       /* y += x */
			tref[0] += GMRFLib_cpu();

			tref[1] -= GMRFLib_cpu();
#pragma omp simd
			for (int j = 0; j < n; j++) {
				yy[j] += x[j];
			}
			tref[1] += GMRFLib_cpu();

			double err = 0.0;
			for (int j = 0; j < n; j++) {
				err = DMAX(err, ABS(y[j] - yy[j]));
			}
			assert(err < FLT_EPSILON);
		}
		printf("mod:  %.4f  plain:  %.4f\n", tref[0] / (tref[0] + tref[1]), tref[1] / (tref[0] + tref[1]));
	}
		break;

	case 121:
	{
		int n = atoi(args[0]);
		int m = atoi(args[1]);
		int inc = atoi(args[2]);

		P(n);
		P(m);
		P(inc);

		double *x = Calloc(3 * inc * n, double);
		double *y = x + inc * n;
		double *yy = x + 2 * inc * n;

		for (int i = 0; i < inc * n; i++) {
			x[i] = GMRFLib_uniform();
			y[i] = GMRFLib_uniform();
			yy[i] = y[i];
		}

		double tref[] = { 0, 0 };
		for (int i = 0; i < m; i++) {
			tref[0] -= GMRFLib_cpu();
			if (inc == 1) {
#pragma omp simd
				for (int j = 0; j < n; j++) {
					y[j] = x[j] - y[j];
				}
			} else {
#pragma omp simd
				for (int j = 0; j < n; j++) {
					y[j] = x[j * inc] - y[j];
				}
			}
			tref[0] += GMRFLib_cpu();

			tref[1] -= GMRFLib_cpu();
			if (inc == 1) {
				for (int j = 0; j < n; j++) {
					yy[j] = x[j] - yy[j];
				}
			} else {
				for (int j = 0; j < n; j++) {
					yy[j] = x[j * inc] - yy[j];
				}
			}
			tref[1] += GMRFLib_cpu();

			double err = 0.0;
			for (int j = 0; j < n; j++) {
				err = DMAX(err, ABS(y[j] - yy[j]));
			}
			assert(err < FLT_EPSILON);
		}
		printf("mod:  %.4f  plain:  %.4f\n", tref[0] / (tref[0] + tref[1]), tref[1] / (tref[0] + tref[1]));
	}
		break;

	case 122:
	{
		int n = atoi(args[0]);
		int m = atoi(args[1]);

		P(n);
		P(m);

		double *x = Calloc(4 * n, double);
		double *y = x + n;
		double *yy = x + 2 * n;
#if defined(INLA_LINK_WITH_MKL)
		double *z = x + 3 * n;
#endif
		for (int i = 0; i < n; i++) {
			x[i] = GMRFLib_uniform();
		}

		double tref[] = { 0, 0 };
		for (int i = 0; i < m; i++) {
			tref[0] -= GMRFLib_cpu();
#pragma omp simd
			for (int j = 0; j < n; j++) {
				y[j] = x[j] + exp(x[j]);
			}
			tref[0] += GMRFLib_cpu();

			tref[1] -= GMRFLib_cpu();
#if defined(INLA_LINK_WITH_MKL)
			vdExp(n, x, z);
			GMRFLib_daxpbyz(n, 1.0, x, 1.0, z, yy);
#else
#pragma omp simd
			for (int j = 0; j < n; j++) {
				y[j] = x[j] + GMRFLib_exp(x[j]);
			}
#endif
			tref[1] += GMRFLib_cpu();

			double err = 0.0;
			for (int j = 0; j < n; j++) {
				err = DMAX(err, ABS(y[j] - yy[j]));
			}
			assert(err < FLT_EPSILON);
		}
		printf("plain:  %.4f  MKL:  %.4f\n", tref[0] / (tref[0] + tref[1]), tref[1] / (tref[0] + tref[1]));
	}
		break;

	case 999:
	{
		GMRFLib_pardiso_check_install(0, 0);
	}
		break;

	default:
	{
		printf("\nNo such test: %d\n", test_no);
	}
		break;
	}
	exit(EXIT_SUCCESS);
}
