// this function should be renamed to something else, as the name is not that representative anymore

#include <math.h>
#include <omp.h>
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"

void GMRFLib_exp(int n, double *x, double *y)
{
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
}

void GMRFLib_exp_inc(int n, double *x, int inc, double *y)
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

void GMRFLib_log(int n, double *x, double *y)
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

void GMRFLib_log1p(int n, double *x, double *y)
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

void GMRFLib_sqr(int n, double *x, double *y)
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

void GMRFLib_sqrt(int n, double *x, double *y)
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

void GMRFLib_add(int n, double *x, double *y, double *z)
{
	// z = x + y
#if defined(INLA_WITH_MKL)
	vdAdd(n, x, y, z);
#else
#       pragma omp simd
	for (int i = 0; i < n; i++) {
		z[i] = x[i] + y[i];
	}
#endif
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((target_clones(INLA_CLONE_TARGETS "default")))
void GMRFLib_mul(int n, double *x, double *y, double *z)
{
	// z = x * y
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

void GMRFLib_daddto(int n, double *x, double *y)
{
	// y = y + x
	int inc = 1;
	double one = 1.0;
	daxpy_(&n, &one, x, &inc, y, &inc);
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((target_clones(INLA_CLONE_TARGETS "default")))
void GMRFLib_cdaddto(int n, double *x, double cx, double *y)
{
	// y = x + const.x
#pragma omp simd
	for (int i = 0; i < n; i++) {
		y[i] = x[i] + cx;
	}
}
#pragma GCC diagnostic pop
