	__m128d sum0 = _mm_setzero_pd();
	__m128d sum1 = _mm_setzero_pd();
	__m128d sum2 = _mm_setzero_pd();
	__m128d sum3 = _mm_setzero_pd();
	int i = 0;
	int limit = n & ~7;
	for (; i < limit; i += 8) {
		__m128d data0 = _mm_load_pd(&x[i]);
		__m128d data1 = _mm_load_pd(&x[i + 2]);
		__m128d data2 = _mm_load_pd(&x[i + 4]);
		__m128d data3 = _mm_load_pd(&x[i + 6]);

		sum0 = _mm_add_pd(sum0, data0);
		sum1 = _mm_add_pd(sum1, data1);
		sum2 = _mm_add_pd(sum2, data2);
		sum3 = _mm_add_pd(sum3, data3);
	}
	sum0 = _mm_add_pd(sum0, sum1);
	sum2 = _mm_add_pd(sum2, sum3);
	sum0 = _mm_add_pd(sum0, sum2);
	__m128d sum_swapped = _mm_shuffle_pd(sum0, sum0, 1);
	__m128d sum_total = _mm_add_pd(sum0, sum_swapped);
	double r;
	_mm_store_sd(&r, sum_total);
#pragma omp simd reduction(+: r)
	for (int ii = i; ii < n; ii++) {
		r += x[ii];
	}
	return r + r0;
