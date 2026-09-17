#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#define LIM 32

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
NOINLINE __attribute__((target_clones(INLA_CLONE_TARGETS "default")))
void GMRFLib_exp(int n, double *RESTRICT x, double *RESTRICT y)
{
#if 0
	static double tref = 0;
	static double trefn = 0;
	tref -= GMRFLib_timer();
#endif

#if defined(INLA_WITH_MKL)
	vdExp(n, x, y);
#elif defined(INLA_WITH_FRAMEWORK_ACCELERATE)
	vvexp(y, x, &n);
#else
#       pragma omp simd
	for (int i = 0; i < n; i++) {
		y[i] = exp(x[i]);
	}
#endif

#if 0
	tref += GMRFLib_timer();
	trefn++;
	if ((int) trefn % 100000 == 0) {
		printf("_exp 1E-6 * %.6f\n", 1.0E6 * tref / trefn);
	}
#endif
}
#pragma GCC diagnostic pop

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
NOINLINE __attribute__((target_clones(INLA_CLONE_TARGETS "default")))
void GMRFLib_exp_inc(int n, double *RESTRICT x, int inc, double *RESTRICT y)
{
	// y = exp(x) with inc
#if defined(INLA_WITH_MKL)
	vdExpI(n, x, inc, y, inc);
#else
#       pragma omp simd
	for (int i = 0; i < n * inc; i += inc) {
		y[i] = exp(x[i]);
	}
#endif
}
#pragma GCC diagnostic pop

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
NOINLINE __attribute__((target_clones(INLA_CLONE_TARGETS "default")))
void GMRFLib_log(int n, double *RESTRICT x, double *RESTRICT y)
{
	// y = log(x)
#if defined(INLA_WITH_MKL)
	vdLn(n, x, y);
#elif defined(INLA_WITH_FRAMEWORK_ACCELERATE)
	vvlog(y, x, &n);
#else
#       pragma omp simd
	for (int i = 0; i < n; i++) {
		y[i] = log(x[i]);
	}
#endif
}
#pragma GCC diagnostic pop

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
NOINLINE __attribute__((target_clones(INLA_CLONE_TARGETS "default")))
void GMRFLib_log1p(int n, double *RESTRICT x, double *RESTRICT y)
{
	// y = log1p(x)
#if defined(INLA_WITH_MKL)
	vdLog1p(n, x, y);
#elif defined(INLA_WITH_FRAMEWORK_ACCELERATE)
	vvlog1p(y, x, &n);
#else
#       pragma omp simd
	for (int i = 0; i < n; i++) {
		y[i] = log1p(x[i]);
	}
#endif
}
#pragma GCC diagnostic pop

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((target_clones(INLA_CLONE_TARGETS "default")))
void GMRFLib_abs(int n, double *RESTRICT x, double *RESTRICT y)
{
#if defined(INLA_WITH_MKL)
	vdAbs(n, x, y);
#elif defined(INLA_WITH_FRAMEWORK_ACCELERATE)
	vvfabs(y, x, &n);
#else
#       pragma omp simd
	for (int i = 0; i < n; i++) {
		y[i] = fabs(x[i]);
	}
#endif
}
#pragma GCC diagnostic pop

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
NOINLINE __attribute__((target_clones(INLA_CLONE_TARGETS "default")))
void GMRFLib_sqr(int n, double *RESTRICT x, double *RESTRICT y)
{
	// y = x * x
#if defined(INLA_WITH_MKL)
	vdSqr(n, x, y);
#else
#       pragma omp simd
	for (int i = 0; i < n; i++) {
		y[i] = SQR(x[i]);
	}
#endif
}
#pragma GCC diagnostic pop

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
NOINLINE __attribute__((target_clones(INLA_CLONE_TARGETS "default")))
void GMRFLib_sqrt(int n, double *RESTRICT x, double *RESTRICT y)
{
	// y = sqrt(x)
#if defined(INLA_WITH_MKL)
	vdSqrt(n, x, y);
#elif defined(INLA_WITH_FRAMEWORK_ACCELERATE)
	vvsqrt(y, x, &n);
#else
#       pragma omp simd
	for (int i = 0; i < n; i++) {
		y[i] = sqrt(x[i]);
	}
#endif
}
#pragma GCC diagnostic pop

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
NOINLINE __attribute__((target_clones(INLA_CLONE_TARGETS "default")))
void GMRFLib_add(int n, double *RESTRICT x, double *RESTRICT y, double *RESTRICT z)
{
	// z = x + y
	if (n <= LIM) {
#pragma omp simd
		for (int i = 0; i < n; i++) {
			z[i] = x[i] + y[i];
		}
		return;
	}
#if defined(INLA_WITH_MKL)
	vdAdd(n, x, y, z);
#else
#       pragma omp simd
	for (int i = 0; i < n; i++) {
		z[i] = x[i] + y[i];
	}
#endif
}
#pragma GCC diagnostic pop

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
NOINLINE __attribute__((target_clones(INLA_CLONE_TARGETS "default")))
void GMRFLib_mul(int n, double *RESTRICT x, double *RESTRICT y, double *RESTRICT z)
{
	// z = x * y
	if (n <= LIM) {
#pragma omp simd
		for (int i = 0; i < n; i++) {
			z[i] = x[i] * y[i];
		}
		return;
	}
#if defined(INLA_WITH_MKL)
	vdMul(n, x, y, z);
#else
#       pragma omp simd
	for (int i = 0; i < n; i++) {
		z[i] = x[i] * y[i];
	}
#endif
}
#pragma GCC diagnostic pop

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
NOINLINE __attribute__((target_clones(INLA_CLONE_TARGETS "default")))
void GMRFLib_daddto(int n, double *RESTRICT x, double *RESTRICT y)
{
	// y = y + x
#pragma omp simd
	for (int i = 0; i < n; i++) {
		y[i] += x[i];
	}
}
#pragma GCC diagnostic pop

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
NOINLINE __attribute__((target_clones(INLA_CLONE_TARGETS "default")))
void GMRFLib_cdaddto(int n, double *RESTRICT x, double cx, double *RESTRICT y)
{
	// y = x + const.x
#pragma omp simd
	for (int i = 0; i < n; i++) {
		y[i] = x[i] + cx;
	}
}
#pragma GCC diagnostic pop

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
NOINLINE __attribute__((target_clones(INLA_CLONE_TARGETS "default")))
void GMRFLib_fm_po_1(int n, double *RESTRICT mask, double *RESTRICT ll, double *RESTRICT wp,
		     double *RESTRICT integral2, double *RESTRICT integral3, double *RESTRICT integral4)
{
	// integral2 = sum mask * exp(ll) * wp
	// integral3 = sum ll * wp
	// integral4 = sum ll^2 * wp

	double r2 = 0.0, r3 = 0.0, r4 = 0.0;
#pragma omp simd reduction(+: r2, r3, r4)
	for (int i = 0; i < n; i++) {
		r2 += exp(ll[i]) * mask[i] * wp[i];
		r3 += ll[i] * wp[i];
		r4 += ll[i] * ll[i] * wp[i];
	}

	*integral2 = r2;
	*integral3 = r3;
	*integral4 = r4;
}
#pragma GCC diagnostic pop
