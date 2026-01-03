#include <float.h>
#include <stdio.h>
#include <math.h>
#include "GMRFLib/GMRFLib.h"

// 2 * phi((x-mu)*exp(lsinv)) * Phi(a*(x-mu)*exp(lsinv)) * exp(lsinv) * exp(Z), where Z allows for un-normalized density

#define NPARAM 4
#define IDX(i_, j_) ((i_) + (j_) * NPARAM)
#define SQR2(a_, b_) SQR(a_)

double ERF(double x)
{
	// robust version
	const double lim = 1.0 - 10.0 * DBL_EPSILON;
	return TRUNCATE(erf(x), -lim, lim);
}

void fitsn_ld(int n, double *x, double *param, double *ld)
{
	double mu = param[0];
	double lsinv = param[1];
	double a = param[2];
	double Z = param[3];

	for (int i = 0; i < n; i++) {
		double cg = -log(0.2e1) / 0.2e1 - log(M_PI) / 0.2e1 + log(exp(-SQR2(-x[i] + mu, 2) * exp((2 * lsinv)) / 0.2e1)) +
		    log1p(-ERF(sqrt(0.2e1) * a * (-x[i] + mu) * exp(lsinv) / 0.2e1)) + lsinv + Z;
		ld[i] = cg;
	}
}

void fitsn_gradhess(double x, double *param, double *grad, double *hess)
{
	// param = [mu, lsinv, a, Z]

	double mu = param[0];
	double lsinv = param[1];
	double a = param[2];
	double POSSIBLY_UNUSED(Z) = param[3];

	double t1, t10, t11, t12, t13, t14, t15, t17, t18, t19, t2, t20, t21, t22, t23, t24,
	    t25, t26, t27, t28, t29, t3, t30, t31, t32, t33, t34, t35, t36, t38, t39, t4,
	    t40, t41, t43, t44, t45, t46, t47, t49, t5, t50, t52, t54, t57, t58, t6, t61,
	    t62, t64, t68, t69, t7, t70, t76, t77, t8, t84, t88, t89, t9, t91, t92;


	t1 = exp(lsinv);
	t2 = sqrt(M_PI);
	t3 = t2 * t1;
	t4 = M_SQRT2;
	t6 = -x + mu;
	t10 = ERF(t1 * t6 * a * t4 / 0.2e1);
	t17 = a * a;
	t18 = t6 * t6;
	t21 = exp(0.2e1 * lsinv);
	t24 = exp(-t21 * t18 * t17 / 0.2e1);
	t33 = 0.1e1 / (t10 - 1) / t2 * (t10 * x * t3 - x * t3 - t10 * mu * t3 + mu * t3 + a * t4 * t24) * t1;
	grad[0] = t33;

	t2 = exp((2 * lsinv));
	t3 = sqrt(M_PI);
	t4 = t3 * t2;
	t5 = x * x;
	t6 = M_SQRT2;
	t9 = exp(lsinv);
	t13 = ERF(t9 * (-x + mu) * a * t6 / 0.2e1);
	t17 = x * mu;
	t23 = mu * mu;
	t27 = a * a;
	t28 = t2 * t27;
	t35 = exp(-t5 * t28 / 0.2e1 + t17 * t28 - t23 * t28 / 0.2e1 + lsinv);
	t36 = t6 * t35;
	t47 =
	    0.1e1 / (t13 - 1) / t3 * (-t13 * t5 * t4 + t5 * t4 + 0.2e1 * t13 * t17 * t4 - 0.2e1 * t17 * t4 - t4 * t13 * t23 + t23 * t4 -
				      a * x * t36 + a * mu * t36 + t13 * t3 - t3);
	grad[1] = t47;


	t1 = a * a;
	t3 = exp((2 * lsinv));
	t4 = t3 * t1;
	t5 = x * x;
	t10 = mu * mu;
	t14 = exp(-t5 * t4 / 0.2e1 + x * mu * t4 - t10 * t4 / 0.2e1 + lsinv);
	t15 = M_SQRT2;
	t17 = -x + mu;
	t18 = sqrt(M_PI);
	t22 = exp(lsinv);
	t26 = ERF(t22 * t17 * a * t15 / 0.2e1);
	t30 = 0.1e1 / (t26 - 1) / t18 * t17 * t15 * t14;
	grad[2] = t30;

	t1 = 1;
	grad[3] = t1;

	t1 = sqrt(M_PI);
	t2 = t1 * M_PI;
	t3 = M_SQRT2;
	t5 = -x + mu;
	t6 = exp(lsinv);
	t10 = ERF(t6 * t5 * a * t3 / 0.2e1);
	t11 = t10 * t10;
	t15 = a * a;
	t18 = exp(0.2e1 * lsinv);
	t19 = t18 * t15;
	t20 = x * x;
	t25 = mu * mu;
	t29 = exp(-t20 * t19 / 0.2e1 + x * mu * t19 - t25 * t19 / 0.2e1 + lsinv);
	t30 = t29 * t15 * a;
	t31 = t3 * t30;
	t35 = M_PI * t3;
	t43 = t5 * t5;
	t46 = exp(-t18 * t43 * t15);
	t54 = SQR2((t10 - 1), 2);
	t58 =
	    -0.1e1 / t54 / t2 * t18 * (t11 * t2 - 0.2e1 * t10 * t2 + t2 - t10 * x * M_PI * t31 + x * t35 * t30 + t10 * mu * M_PI * t31 -
				       mu * t35 * t30 + 0.2e1 * t1 * t15 * t46);
	hess[IDX(0, 0)] = t58;

	t1 = exp(lsinv);
	t2 = sqrt(M_PI);
	t3 = t2 * M_PI;
	t4 = t3 * t1;
	t5 = M_SQRT2;
	t7 = -x + mu;
	t11 = ERF(t1 * t7 * a * t5 / 0.2e1);
	t12 = t11 * t11;
	t24 = t11 * mu;
	t29 = a * a;
	t31 = 0.2e1 * lsinv;
	t32 = exp(t31);
	t33 = t32 * t29;
	t34 = x * x;
	t35 = t34 * t33;
	t38 = x * mu * t33;
	t39 = mu * mu;
	t40 = t39 * t33;
	t43 = exp(t31 - t35 / 0.2e1 + t38 - t40 / 0.2e1);
	t44 = t43 * t29 * a;
	t45 = t5 * t44;
	t49 = M_PI * t5;
	t52 = M_PI * x;
	t64 = t7 * t7;
	t68 = exp(-t32 * t64 * t29 / 0.2e1);
	t69 = t5 * t68;
	t70 = M_PI * a;
	t76 = exp(-t35 + 0.2e1 * t38 - t40 + lsinv);
	t77 = t29 * t76;
	t84 =
	    0.2e1 * t12 * x * t4 - 0.4e1 * t11 * x * t4 + 0.2e1 * x * t4 - 0.2e1 * t12 * mu * t4 + 0.4e1 * t4 * t24 - 0.2e1 * mu * t4 -
	    t11 * t34 * M_PI * t45 + t34 * t49 * t44 + 0.2e1 * t24 * t52 * t45 - 0.2e1 * mu * t52 * t45 - t11 * t39 * M_PI * t45 + t39 * t49 * t44 +
	    t11 * t70 * t69 - t70 * t69 + 0.2e1 * x * t2 * t77 - 0.2e1 * mu * t2 * t77;
	t88 = SQR2(t11 - 0.1e1, 2);
	t91 = 0.1e1 / t88 / t3 * t1 * t84;
	hess[IDX(0, 1)] = hess[IDX(1, 0)] = t91;

	t1 = a * a;
	t3 = exp((2 * lsinv));
	t4 = t3 * t1;
	t5 = M_SQRT2;
	t6 = t5 * t4;
	t7 = x * x;
	t11 = exp(lsinv);
	t15 = ERF(t11 * (-x + mu) * a * t5 / 0.2e1);
	t18 = M_PI * t5;
	t21 = M_PI * x;
	t29 = mu * mu;
	t43 = exp(-t7 * t4 / 0.2e1 + x * mu * t4 - t29 * t4 / 0.2e1 + lsinv);
	t44 = a * t43;
	t45 = sqrt(M_PI);
	t57 = SQR2((t15 - 1), 2);
	t61 =
	    -0.1e1 / t57 / t45 / M_PI * t43 * (t15 * t7 * M_PI * t6 - t7 * t18 * t4 - 0.2e1 * t15 * mu * t21 * t6 + 0.2e1 * mu * t21 * t6 +
					       t15 * t29 * M_PI * t6 - t29 * t18 * t4 - t18 * t15 + t18 - 0.2e1 * x * t45 * t44 +
					       0.2e1 * mu * t45 * t44);
	hess[IDX(0, 2)] = hess[IDX(2, 0)] = t61;

	t1 = 0.0;
	hess[IDX(0, 3)] = hess[IDX(3, 0)] = t1;

	t1 = -x + mu;
	t2 = exp(lsinv);
	t3 = sqrt(M_PI);
	t4 = M_PI * t3;
	t5 = t4 * t2;
	t6 = M_SQRT2;
	t11 = ERF(t2 * t1 * a * t6 / 0.2e1);
	t12 = t11 * t11;
	t24 = t11 * mu;
	t29 = a * a;
	t31 = 0.2e1 * lsinv;
	t32 = exp(t31);
	t33 = t32 * t29;
	t34 = x * x;
	t35 = t34 * t33;
	t38 = x * mu * t33;
	t39 = mu * mu;
	t40 = t39 * t33;
	t43 = exp(t31 - t35 / 0.2e1 + t38 - t40 / 0.2e1);
	t44 = t43 * t29 * a;
	t45 = t6 * t44;
	t49 = M_PI * t6;
	t52 = M_PI * x;
	t64 = t1 * t1;
	t68 = exp(-t32 * t64 * t29 / 0.2e1);
	t69 = t6 * t68;
	t70 = M_PI * a;
	t76 = exp(-t35 + 0.2e1 * t38 - t40 + lsinv);
	t77 = t29 * t76;
	t84 =
	    0.2e1 * t12 * x * t5 - 0.4e1 * t11 * x * t5 + 0.2e1 * x * t5 - 0.2e1 * t12 * mu * t5 + 0.4e1 * t24 * t5 - 0.2e1 * mu * t5 -
	    t11 * t34 * M_PI * t45 + t34 * t49 * t44 + 0.2e1 * t24 * t52 * t45 - 0.2e1 * mu * t52 * t45 - t11 * t39 * M_PI * t45 + t39 * t49 * t44 +
	    t11 * t70 * t69 - t70 * t69 + 0.2e1 * x * t3 * t77 - 0.2e1 * mu * t3 * t77;
	t89 = SQR2(t11 - 0.1e1, 2);
	t92 = 0.1e1 / t89 / t4 * t2 * t1 * t84;
	hess[IDX(1, 1)] = t92;

	t1 = -x + mu;
	t2 = a * a;
	t4 = exp((2 * lsinv));
	t5 = t4 * t2;
	t6 = M_SQRT2;
	t7 = t6 * t5;
	t8 = x * x;
	t11 = exp(lsinv);
	t15 = ERF(t11 * t1 * a * t6 / 0.2e1);
	t18 = M_PI * t6;
	t21 = M_PI * x;
	t29 = mu * mu;
	t43 = exp(-t8 * t5 / 0.2e1 + x * mu * t5 - t29 * t5 / 0.2e1 + lsinv);
	t44 = a * t43;
	t45 = sqrt(M_PI);
	t58 = SQR2((t15 - 1), 2);
	t62 =
	    -0.1e1 / t58 / t45 / M_PI * t43 * (t15 * t8 * M_PI * t7 - t8 * t18 * t5 - 0.2e1 * t15 * mu * t21 * t7 + 0.2e1 * mu * t21 * t7 +
					       t15 * t29 * M_PI * t7 - t29 * t18 * t5 - t18 * t15 + t18 - 0.2e1 * x * t45 * t44 +
					       0.2e1 * mu * t45 * t44) * t1;
	hess[IDX(1, 2)] = hess[IDX(2, 1)] = t62;

	t1 = 0.0;
	hess[IDX(1, 3)] = hess[IDX(3, 1)] = t1;

	t1 = exp(lsinv);
	t2 = t1 * a;
	t3 = M_SQRT2;
	t4 = t3 * t2;
	t7 = -x + mu;
	t11 = ERF(t1 * t7 * a * t3 / 0.2e1);
	t14 = M_PI * t3;
	t22 = a * a;
	t23 = t7 * t7;
	t25 = 0.2e1 * lsinv;
	t26 = exp(t25);
	t29 = exp(-t26 * t23 * t22 / 0.2e1);
	t30 = sqrt(M_PI);
	t35 = t26 * t22;
	t36 = x * x;
	t41 = mu * mu;
	t45 = exp(t25 - t36 * t35 / 0.2e1 + x * mu * t35 - t41 * t35 / 0.2e1);
	t50 = SQR2((t11 - 1), 2);
	t54 =
	    -0.1e1 / t50 / t30 / M_PI * t45 * t23 * (-t11 * x * M_PI * t4 + x * t14 * t2 + t11 * mu * M_PI * t4 - mu * t14 * t2 +
						     0.2e1 * t30 * t29);
	hess[IDX(2, 2)] = t54;

	t1 = 0.0;
	hess[IDX(2, 3)] = hess[IDX(3, 2)] = t1;

	t1 = 0.0;
	hess[IDX(3, 3)] = t1;
}

__attribute__((target_clones(INLA_CLONE_TARGETS "default")))
void fitsn_fit(int n, double *w, double *x, double *y, GMRFLib_sn_param_tp *sn)
{
	// param = arg min 1/2 * \sum_i=0^n-1 w_i * (ld(x_i, ...) - y_i)^2

	int verbose = 1;
	double *ld = Malloc(n, double);
	double *p = Malloc(n, double);
	double *param = Calloc(NPARAM, double);
	double m0 = 0.0, m1 = 0.0, m2 = 0.0, m3 = 0.0;
	double eps_f = 1.0e-6;

	// compute initial values
	Memcpy(ld, y, n * sizeof(double));
	GMRFLib_adjust_vector(ld, n);
	GMRFLib_exp(n, ld, p);

	for (int i = 0; i < n; i++) {
		m0 += p[i];
		m1 += p[i] * x[i];
		m2 += p[i] * SQR(x[i]);
	}
	double mean, stdev, skew;
	mean = m1 / m0;
	stdev = sqrt(DMAX(FLT_EPSILON, m2 / m0 - SQR(mean)));

	for (int i = 0; i < n; i++) {
		m3 += p[i] * POW3(x[i] - mean);
	}
	skew = TRUNCATE((m3 / m0) / POW3(stdev), -0.9, 0.9);
	GMRFLib_sn_moments2par(sn, &mean, &stdev, &skew);

	param[0] = sn->xi;
	param[1] = log(1.0 / sn->omega);
	param[2] = sn->alpha;

	double aa = 0.0;
	param[3] = 0.0;
	fitsn_ld(n, x, param, ld);
	for (int i = 0; i < n; i++) {
		aa += w[i] * (y[i] - ld[i]);
	}
	double wsi = 1.0 / GMRFLib_dsum(n, w);
	param[3] = aa * wsi;

	int iter_max = 100;
	double *grad = Malloc(NPARAM, double);
	double *Grad = Malloc(NPARAM, double);
	double *hess = Malloc(ISQR(NPARAM), double);
	double *Hess = Malloc(ISQR(NPARAM), double);
	double *delta = Malloc(NPARAM, double);
	double func_val = 0.0;
	double func_val_p = 0.0;

	for (int iter = 0; iter <= iter_max; iter++) {	       /* yes: <= */
		assert(iter < iter_max);
		func_val_p = func_val;
		func_val = 0.0;

		GMRFLib_dfill(NPARAM, 0.0, Grad);
		GMRFLib_dfill(ISQR(NPARAM), 0.0, Hess);
		fitsn_ld(n, x, param, ld);
		for (int i = 0; i < n; i++) {
			double res = w[i] * (y[i] - ld[i]);
			func_val += w[i] * SQR(y[i] - ld[i]);
			fitsn_gradhess(x[i], param, grad, hess);
			for (int k = 0; k < NPARAM; k++) {
				Grad[k] += res * (-grad[k]);
				for (int kk = k; kk < NPARAM; kk++) {
					Hess[IDX(k, kk)] += w[i] * grad[k] * grad[kk] - res * hess[IDX(k, kk)];
					Hess[IDX(kk, k)] = Hess[IDX(k, kk)];
				}
			}
		}

		double *chol = NULL;
		GMRFLib_ensure_spd(Hess, NPARAM, -1.0, NULL);
		GMRFLib_comp_chol_general(&chol, Hess, NPARAM, NULL, !GMRFLib_SUCCESS);
		GMRFLib_solveAxb_posdef(delta, chol, Grad, NPARAM, 1);
		Free(chol);

		for (int k = 0; k < NPARAM; k++) {
			param[k] -= delta[k];
		}
		func_val = sqrt(func_val * wsi);
		if (verbose) {
			printf("\tfitsn: iter %1d param = %.4g %.4g %.4g %.4g |diff.fval| %.6g\n",
			       iter, param[0], param[1], param[2], param[3], ABS(func_val - func_val_p));
		}
		if (iter > 0 && ABS(func_val - func_val_p) < eps_f) {
			break;
		}
	}

	sn->xi = param[0];
	sn->omega = exp(-param[1]);
	sn->alpha = param[2];

	if (verbose) {
		printf("\tFit found xi = %.4g omega %.4g alpha %.4g\n", sn->xi, sn->omega, sn->alpha);
	}

	Free(grad);
	Free(Grad);
	Free(hess);
	Free(Hess);
	Free(delta);
	Free(ld);
	Free(p);
}

#if defined(INLA_WITH_DEVEL)

void fitsn_test_grad(void)
{
	double param[NPARAM] = { 0.1, 1.2, 2.3, -1.0 };
	double x = 0.123;
	double h = 1e-4;
	double grad[NPARAM];
	double hess[ISQR(NPARAM)];
	double *pparam = Malloc(NPARAM, double);

	fitsn_gradhess(x, param, grad, hess);

	for (int k = 0; k < NPARAM; k++) {
		double ld[2];
		Memcpy(pparam, param, NPARAM * sizeof(double));

		pparam[k] = param[k] - h;
		fitsn_ld(1, &x, pparam, ld);
		pparam[k] = param[k] + h;
		fitsn_ld(1, &x, pparam, ld + 1);

		printf("Grad %d %f %f\n", k, (ld[1] - ld[0]) / (2.0 * h), grad[k]);
	}

	Free(pparam);
}

void fitsn_test_hess(void)
{
	double param[NPARAM] = { 0.1, 1.1, 2.2, -1.0 };
	double x = 0.123;
	double h = 1e-4;
	double grad[NPARAM];
	double hess[ISQR(NPARAM)];
	double *pparam = Malloc(NPARAM, double);

	fitsn_gradhess(x, param, grad, hess);

	for (int k = 0; k < NPARAM; k++) {
		double ld[3];
		Memcpy(pparam, param, NPARAM * sizeof(double));

		pparam[k] = param[k] - h;
		fitsn_ld(1, &x, pparam, ld);
		pparam[k] = param[k] + h;
		fitsn_ld(1, &x, pparam, ld + 2);
		pparam[k] = param[k];
		fitsn_ld(1, &x, pparam, ld + 1);

		printf("Hess[%d %d] %f %f\n", k, k, (ld[2] - 2 * ld[1] + ld[0]) / SQR(h), hess[IDX(k, k)]);
	}

	Free(pparam);
}

void fitsn_test(void)
{
	/*
	 * library(sn) x <- c(-2, -1, 0, 1, 2.2) y <- c(-2.0, -0.5, 0, -0.54, -3.3) param <- list(xi = 0.687, omega = 1.212, alpha = -1.313) Z <-
	 * 0.8858 plot(x, y) xx <- seq(-4, 4, by = 0.01) lines(xx, dsn(xx, log = TRUE, xi = param$xi, omega = param$omega, alpha = param$alpha) +
	 * Z) 
	 */

	double w[] = { 1, 1, 1, 1, 1 };
	double x[] = { -2, -1, 0, 1, 2.2 };
	double y[] = { -2.0, -0.5, 0, -0.54, -3.3 };

	GMRFLib_sn_param_tp sn;
	// fitsn_test_grad();
	// fitsn_test_hess();
	fitsn_fit(5, w, x, y, &sn);
}

#else

void fitsn_test_grad(void)
{
}
void fitsn_test_hess(void)
{
}
void fitsn_test(void)
{
}

#endif
