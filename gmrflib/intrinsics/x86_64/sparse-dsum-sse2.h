	int i = 0;
	int limit = n & ~3;
	double r = 0.0;
	if (limit > 0) {
		__m128d sum0 = _mm_setzero_pd();
		__m128d sum1 = _mm_setzero_pd();
		for (; i < limit; i += 4) {
			__m128d val1 = _mm_load_sd(&a[idx[i]]);
			__m128d val2 = _mm_load_sd(&a[idx[i + 1]]);
			__m128d val3 = _mm_load_sd(&a[idx[i + 2]]);
			__m128d val4 = _mm_load_sd(&a[idx[i + 3]]);
			__m128d data0 = _mm_unpacklo_pd(val1, val2);
			__m128d data1 = _mm_unpacklo_pd(val3, val4);
			sum0 = _mm_add_pd(sum0, data0);
			sum1 = _mm_add_pd(sum1, data1);
		}
		sum0 = _mm_add_pd(sum0, sum1);
		__m128d sum_swapped = _mm_shuffle_pd(sum0, sum0, 1);
		__m128d sum_total = _mm_add_pd(sum0, sum_swapped);
		_mm_store_sd(&r, sum_total);
	}
#pragma omp simd reduction(+: r)
	for (int ii = i; ii < n; ii++) {
		r += a[idx[ii]];
	}
	return r;
