#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>

#include "GMRFLib/GMRFLib.h"
#include "inla.h"
#include "fast-math/special-functions.h"

FORCEINLINE double inla_lgamma(double x)
{
	return lgamma(x);
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
NOINLINE __attribute__((target_clones(INLA_CLONE_TARGETS "default")))
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
NOINLINE __attribute__((target_clones(INLA_CLONE_TARGETS "default")))
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

FORCEINLINE double inla_gamma(double x)
{
	return (exp(lgamma(x)));
}

FORCEINLINE double inla_gamma_fast(double x)
{
	return (exp(inla_lgamma_fast(x)));
}

FORCEINLINE double inla_beta(double a, double b)
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
	size_t m3 = 3 * m;
	double x[m3];
	for (size_t i = 0, j = 0; i < m3; i += 3, j++) {
		x[i] = a[j];
		x[i + 1] = b[j];
		x[i + 2] = a[j] + b[j];
	}

	double r[m3];
	LGAMMAfn_m(m3, x, r);

	for (size_t i = 0, j = 0; i < m3; i += 3, j++) {
		llbeta[j] = r[i] + r[i + 1] - r[i + 2];
	}
}
