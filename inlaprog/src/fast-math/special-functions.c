#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <strings.h>
#include <stdio.h>

#include "GMRFLib/GMRFLib.h"
#include "inla.h"
#include "fast-math/special-functions.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((target_clones(INLA_CLONE_TARGETS "default")))
double inla_logcdf_normal(double x)
{
	// return the log of the cummulative distribution function for a standard normal.
	// This version is ok for all x 
	if (ABS(x) <= 7.0) {
		return (log(GMRFLib_cdfnorm(x)));
	} else {
		double t1, t4, t3, t8, t9, t13, t27, t28, t31, t47;

		if (x > 7.0) {
			t1 = 1.77245385090551602729816748334;
			t3 = M_SQRT2;
			t4 = t3 / t1;
			t8 = x * x;
			t9 = t8 * x;
			t13 = t8 * t8;
			t27 = exp(t8);
			t28 = sqrt(t27);
			t31 = 0.1e1 / M_PI;
			t47 =
			    0.1e1 / t28 * (-0.1e1 / x * t4 / 0.2e1 + 0.1e1 / t9 * t4 / 0.2e1 - 0.3e1 / 0.2e1 / t13 / x * t4 +
					   0.15e2 / 0.2e1 / t13 / t9 * t4)
			    + 0.1e1 / t27 * (-0.1e1 / t8 * t31 / 0.4e1 + 0.1e1 / t13 * t31 / 0.2e1 - 0.7e1 / 0.4e1 / t13 / t8 * t31);
			return t47;
		} else {
			// x < -7.0
			double xx = -x, cg1;
			cg1 =
			    -(pow(xx, 0.6e1) + log(0.2e1) * pow(xx, 0.4e1) + log(0.3141592653589793e1) * pow(xx, 0.4e1) +
			      0.2e1 * log(xx) * pow(xx, 0.4e1)
			      - 0.5e1 + 0.2e1 * xx * xx) * pow(xx, -0.4e1) / 0.2e1;
			return (cg1);
		}
	}
	abort();
	return 0;
}
#pragma GCC diagnostic pop

double inla_cdf_normal(double x)
{
	/*
	 * the un-log version of inla_logcdf_normal 
	 */
	if (ABS(x) < 7.0) {
		return GMRFLib_cdfnorm(x);
	} else {
		return exp(inla_logcdf_normal(x));
	}
}

double inla_cdf_normal_fast(double x)
{
	// a faster approximation, see misc/doc/doc/approximate-cdf-normal.pdf
	if (ABS(x) <= 7.0) {
		// see misc/doc/doc/approximate-cdf-normal.pdf
		// sqrt(M_PI / 8.0) = 0.6266570686577502....
		if (x > 0.0) {
			return (0.5 + 0.5 * sqrt(ONE_mexp(-0.6266570686577502 * SQR(x))));
		} else {
			return (1.0 - (0.5 + 0.5 * sqrt(ONE_mexp(-0.6266570686577502 * SQR(x)))));
		}
		abort();
		return (0.5 + 0.5 * sqrt(ONE_mexp(-0.6266570686577502 * SQR(x))));
	} else {
		return inla_cdf_normal(x);
	}
}

double inla_logitcdf_normal(double x)
{
	// return log(Phi(x)/(1-Phi(x)))
#define M_LN_SQRT_2PI       0.918938533204672741780329736406

	if (ABS(x) < 7.0) {
		double y = inla_cdf_normal(x);
		return (log(y / (1.0 - y)));
	} else {
		// > asympt(log(Phi(x)/(1-Phi(x))), x, 16); 
		// 2
		// x 1/2 1/2 1
		// ---- + ln(x) + ln(2 Pi ) + O(----)
		// 2 2
		// 
		double val = (SQR(x) / 2.0 + log(x) + M_LN_SQRT_2PI);
		return (x > 0.0 ? val : -val);
	}
#undef M_LN_SQRT_2PI
}

double inla_logcdf_normal_fast(double x)
{
	// a faster approximation, see misc/doc/doc/approximate-cdf-normal.pdf
	// sqrt(M_PI / 8.0) = 0.6266570686577502....
	// log(1.0/4.0) = -1.386294361119891...
	if (ABS(x) < 7.0) {
		return (log(inla_cdf_normal_fast(x)));
	} else {
		if (x > 7.0) {
			return (-0.25 * exp(-0.6266570686577502 * SQR(x)));
		} else {
			// return (log(1.0 / 4.0) - 0.6266570686577502 * SQR(x));
			return (-1.386294361119891 - 0.6266570686577502 * SQR(x));
		}
	}
}

forceinline double inla_ipow(double x, int k)
{
	// x^k
	return gsl_sf_pow_int(x, k);
}

forceinline double inla_lgamma(double x)
{
	return lgamma(x);
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((target_clones(INLA_CLONE_TARGETS "default")))
double inla_lgamma_fast(double x)
{
	if (unlikely(x <= 0.0)) {
		return lgamma(x);
	}
#define G 5
#define N 7
	static const double p[N] = {
		1.0000000001900148240,
		76.180091729471463483,
		-86.505320329416767652,
		24.014098240830910490,
		-1.2317395724501553875,
		0.0012086509738661785061,
		-5.3952393849531283785e-6
	};

	double tmp = x + G + 0.5;
	tmp -= (x + 0.5) * log(tmp);

	double ser = p[0];
	for (int i = 1; i < N; ++i) {
		ser += p[i] / (x + i);
	}
#undef G
#undef N
	return -tmp + log(2.5066282746310005 * ser / x);
}
#pragma GCC diagnostic pop


void inla_lgamma_m(size_t m, double *RESTRICT x, double *RESTRICT res)
{
	for (size_t i = 0; i < m; i++) {
		res[i] = lgamma(x[i]);
	}
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((target_clones(INLA_CLONE_TARGETS "default")))
void XXXinla_lgamma_fast_m(size_t m, double *RESTRICT x, double *RESTRICT res)
{
	// evaluate M calls to lgamma_fast together, assume all x[] > 0. this is what is used for the lbeta(a,b) function, for
	// which all arguments are positive: lbeta(a,b) := lgamma(a)+lgamma(b)-lgamma(a+b)
#define G 5
#define N 7
	static const double p[N] = {
		1.0000000001900148240,
		76.180091729471463483,
		-86.505320329416767652,
		24.014098240830910490,
		-1.2317395724501553875,
		0.0012086509738661785061,
		-5.3952393849531283785e-6
	};
	double tmp[m];
	for (size_t i = 0; i < m; i++) {
		double tt = x[i] + G + 0.5;
		tmp[i] = tt - (x[i] + 0.5) * log(tt);
	}

	double ser[m];
	// this is for N==7
#pragma omp simd
	for (size_t i = 0; i < m; i++) {
		double xi = x[i];
		double s = p[0];
		s += p[1] / (xi + 1.0);
		s += p[2] / (xi + 2.0);
		s += p[3] / (xi + 3.0);
		s += p[4] / (xi + 4.0);
		s += p[5] / (xi + 5.0);
		s += p[6] / (xi + 6.0);
		ser[i] = s;
	}
#pragma omp simd
	for (size_t i = 0; i < m; i++) {
		res[i] = -tmp[i] + log(2.5066282746310005 * ser[i] / x[i]);
	}
#undef G
#undef N
}
#pragma GCC diagnostic pop

#include <stddef.h>
#include <math.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((target_clones(INLA_CLONE_TARGETS "default")))
void inla_lgamma_fast_m(size_t m, double *RESTRICT x, double *RESTRICT res)
{
#define G 5
#define N 7
	static const double p[N] = {
		1.0000000001900148240,
		76.180091729471463483,
		-86.505320329416767652,
		24.014098240830910490,
		-1.2317395724501553875,
		0.0012086509738661785061,
		-5.3952393849531283785e-6
	};

	// Fusing into a single SIMD loop completely eliminates 'tmp' and 'ser' arrays in the old version of the function
#pragma omp simd
	for (size_t i = 0; i < m; i++) {
		double xi = x[i];
		double tt = xi + (G + 0.5);

		// Compute the Lanczos series approximation
		double s = p[0];
		s += p[1] / (xi + 1.0);
		s += p[2] / (xi + 2.0);
		s += p[3] / (xi + 3.0);
		s += p[4] / (xi + 4.0);
		s += p[5] / (xi + 5.0);
		s += p[6] / (xi + 6.0);

		// Inline computation of combined math terms
		double log_term = log((2.5066282746310005 * s) / xi);
		double tmp_term = tt - (xi + 0.5) * log(tt);

		res[i] = -tmp_term + log_term;
	}
#undef G
#undef N
}
#pragma GCC diagnostic pop


forceinline double inla_gamma(double x)
{
	return (exp(lgamma(x)));
}

forceinline double inla_gamma_fast(double x)
{
	return (exp(inla_lgamma_fast(x)));
}

forceinline double inla_beta(double a, double b)
{
	return exp(inla_lbeta(a, b));
}

double inla_lbeta(double a, double b)
{
	double x[3] = { a, b, a + b };
	double res[3] = { 0 };
	LGAMMAfn_m(3, x, res);
	return res[0] + res[1] - res[2];
}

void inla_lbeta_m(size_t m, double *RESTRICT a, double *RESTRICT b, double *RESTRICT llbeta)
{
	double x[3 * m];
	for (size_t i = 0, j = 0; i < 3 * m; i += 3, j++) {
		x[i] = a[j];
		x[i + 1] = b[j];
		x[i + 2] = a[j] + b[j];
	}
	double r[3 * m];
	LGAMMAfn_m(3 * m, x, r);
	for (size_t i = 0, j = 0; i < 3 * m; i += 3, j++) {
		llbeta[j] = r[i] + r[i + 1] - r[i + 2];
	}
}
