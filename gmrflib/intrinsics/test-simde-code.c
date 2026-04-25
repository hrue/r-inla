#include <assert.h>
#include <math.h>
#include <stdalign.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SIMDE_ENABLE_NATIVE_ALIASES
#define SIMDE_FAST_MATH
#if !defined(_OPENMP)
#       define SIMDE_ENABLE_OPENMP
#endif
#
#include <simde/simde-common.h>
#include <simde/x86/sse2.h>
#
#if defined(TEST_AVX2)
#       include <simde/x86/avx2.h>
#       include <simde/x86/fma.h>
#endif
#
#if defined(TEST_AVX512F)
#       include <simde/x86/avx512.h>
#       include <simde/x86/fma.h>
#endif

#define Calloc(n_, type_)  (type_ *)calloc((n_), sizeof(type_))
#define Malloc(n_, type_)  (type_ *)malloc((n_) * sizeof(type_))
#define Sqr(x_) ((x_)*(x_))

int test_ddot2(int n, double *x, double *y, double *z)
{
	double *a = Calloc(1, double);
	double *b = Calloc(1, double);

#if defined(TEST_SSE2)
#       include "simde/ddot2-sse2.h"
#elif defined(TEST_AVX2)
#       include "simde/ddot2-avx2.h"
#elif defined(TEST_AVX512F)
#       include "simde/ddot2-avx512f.h"
#else
#       error "NOT VALID CHOICE"
#endif

	double a_ref = 0.0, b_ref = 0.0;
	for (int i = 0; i < n; i++) {
		a_ref += x[i] * y[i];
		b_ref += x[i] * z[i];
	}

	assert(*a == a_ref);
	assert(*b == b_ref);
	printf("ddot2 pass\n");
	return 1;
}

double test_dsum_sse2(int n, double *x)
{
	double r0 = 0;
#if defined(TEST_SSE2)
#       include "simde/dsum-sse2.h"
#endif
	return r0;
}

double test_dsum_avx2(int n, double *x)
{
	double r0 = 0;
#if defined(TEST_AVX2)
#       include "simde/dsum-avx2.h"
#endif
	return r0;
}

double test_dsum_avx512f(int n, double *x)
{
	double r0 = 0;
#if defined(TEST_AVX512F)
#       include "simde/dsum-avx512f.h"
#endif
	return r0;
}

int test_dsum(int n, double *x)
{
	double sum = 0.0;
#if defined(TEST_SSE2)
	sum = test_dsum_sse2(n, x);
#elif defined(TEST_AVX2)
	sum = test_dsum_avx2(n, x);
#elif defined(TEST_AVX512F)
	sum = test_dsum_avx512f(n, x);
#else
#       error "NOT VALID CHOICE"
#endif

	double sum_ref = 0.0;
	for (int i = 0; i < n; i++) {
		sum_ref += x[i];
	}

	assert(sum == sum_ref);
	printf("dsum pass\n");
	return 1;
}

int test_dscale2(int n, double a, double *x, double *y)
{
	double *yy = Calloc(n, double);
#if defined(TEST_SSE2)
#       include "simde/dscale2-sse2.h"
#elif defined(TEST_AVX2)
#       include "simde/dscale2-avx2.h"
#elif defined(TEST_AVX512F)
#       include "simde/dscale2-avx512f.h"
#else
#       error "NOT VALID CHOICE"
#endif

	for (int i = 0; i < n; i++) {
		yy[i] = a * x[i];
	}

	int ok = 1;
	for (int i = 0; i < n; i++) {
		ok = ok && (y[i] == yy[i]);
	}
	assert(ok == 1 && "dscale2 does not pass");
	printf("dscale2 pass\n");
	return 0;
}

int test_is_sorted_double_sse2(int n, double *a)
{
#if defined(TEST_SSE2)
#       include "simde/is-sorted-double-sse2.h"
#endif
	return -1;
}

int test_is_sorted_double_avx2(int n, double *a)
{
#if defined(TEST_AVX2)
#       include "simde/is-sorted-double-avx2.h"
#endif
	return -1;
}

int test_is_sorted_double_avx512f(int n, double *a)
{
#if defined(TEST_AVX512F)
#       include "simde/is-sorted-double-avx512.h"
#endif
	return -1;
}

int test_is_sorted_double(int n, double *x)
{
	int value = 0.0;
#if defined(TEST_SSE2)
	value = test_is_sorted_double_sse2(n, x);
#elif defined(TEST_AVX2)
	value = test_is_sorted_double_avx2(n, x);
#elif defined(TEST_AVX512F)
	value = test_is_sorted_double_avx512f(n, x);
#else
#       error "NOT VALID CHOICE"
#endif

	int ref_value = 1;
	for (int i = 0; i < n - 1; i++) {
		ref_value = ref_value && (x[i + 1] >= x[i]);
	}

	if (value != ref_value) {
		printf("value %d ref_value %d\n", value, ref_value);
		for (int i = 0; i < n; i++)
			printf("x[%1d] =  %g\n", i, x[i]);
	}
	assert(value == ref_value);
	printf("is_sorted_double pass\n");
	return 1;
}

int test_is_sorted_int_sse2(int n, int *a)
{
#if defined(TEST_SSE2)
#       include "simde/is-sorted-int-sse2.h"
#endif
	return -1;
}

int test_is_sorted_int_avx2(int n, int *a)
{
#if defined(TEST_AVX2)
#       include "simde/is-sorted-int-avx2.h"
#endif
	return -1;
}

int test_is_sorted_int_avx512f(int n, int *a)
{
#if defined(TEST_AVX512F)
#       include "simde/is-sorted-int-avx512.h"
#endif
	return -1;
}

int test_is_sorted_int(int n, int *x)
{
	int value = 0.0;
#if defined(TEST_SSE2)
	value = test_is_sorted_int_sse2(n, x);
#elif defined(TEST_AVX2)
	value = test_is_sorted_int_avx2(n, x);
#elif defined(TEST_AVX512F)
	value = test_is_sorted_int_avx512f(n, x);
#else
#       error "NOT VALID CHOICE"
#endif

	int ref_value = 1;
	for (int i = 0; i < n - 1; i++) {
		ref_value = ref_value && (x[i + 1] >= x[i]);
	}

	assert(value == ref_value);
	printf("is_sorted_int pass\n");
	return 1;
}

int test_zero_small(int n, double eps, double *x)
{
	double *xx = Calloc(n, double);
	memcpy(xx, x, n * sizeof(double));


#if defined(TEST_SSE2)
#       include "simde/zero-small-sse2.h"
#elif defined(TEST_AVX2)
#       include "simde/zero-small-avx2.h"
#elif defined(TEST_AVX512F)
#       include "simde/zero-small-avx512f.h"
#else
#       error "NOT VALID CHOICE"
#endif

	for (int i = 0; i < n; i++) {
		if (fabs(xx[i]) <= eps)
			xx[i] = 0.0;
	}

	int ok = 1;
	for (int i = 0; i < n; i++) {
		ok = ok && x[i] == xx[i];
	}
	if (!ok) {
		printf("eps %g\n", eps);
		for (int i = 0; i < n; i++) {
			printf("x[%1d] = %g xx =  %g\n", i, x[i], xx[i]);
		}
	}
	assert(ok && "zero_small test");
	printf("test_zero_test pass\n");
	return 1;
}

int test_isum_sse2(int n, int *x)
{
#if defined(TEST_SSE2)
#       include "simde/isum-sse2.h"
#endif
	return -1;
}

int test_isum(int n, int *x)
{
#if defined(TEST_SSE2)
	int sum = 0;
	sum = test_isum_sse2(n, x);
	int ref_sum = 0;
	for (int i = 0; i < n; i++)
		ref_sum += x[i];
	if (sum != ref_sum) {
		printf("sum %d ref_sum %d\n", sum, ref_sum);
	}
	assert(sum == ref_sum && "isum failed");
	printf("test_isum pass\n");
#endif
	return 0;
}


int main(int argc, char **argv)
{
	assert(argc > 1 && "argc > 1");
	
	int n = atoi(argv[1]);
	int nn = (n < 16 ? 16 : n);
	
	double *x = Calloc(nn, double);
	double *xx = Calloc(nn, double);
	double *xxx = Calloc(nn, double);
	double *y = Calloc(nn, double);
	double *z = Calloc(nn, double);

	double *xs = Calloc(nn, double);
	double *xns = Calloc(nn, double);

	int *ixs = Calloc(nn, int);
	int *ixns = Calloc(nn, int);

	for (int i = 0; i < nn; i++) {
		x[i] = i;
		xx[i] = i;
		xxx[i] = -i;
		y[i] = 2 * nn - i;
		z[i] = 1 - i;
		xs[i] = i;
		xns[i] = (i == nn - 1 ? -1 : i);
		ixs[i] = i;
		ixns[i] = (i == nn - 1 ? -1 : i);
	}

	if (n >= 16) {
		test_ddot2(16, x, y, z);
	} else {
		printf("test_ddot2 require n >= 16\n");
	}

	test_dsum(n, x);
	test_dscale2(n, 3.0, x, y);

	printf("test sorted\n");
	test_is_sorted_double(n, xs);
	test_is_sorted_int(n, ixs);

	printf("test not sorted\n");
	test_is_sorted_double(n, xns);
	test_is_sorted_int(n, ixns);

	test_zero_small(n, n / 2.0, xx);
	test_zero_small(n, -n / 2.0, xxx);

	test_isum(n, ixs);
	test_isum(n, ixns);
}
