	int i = 0;
	int limit = n & ~3;
	double result = 0.0;
	if (limit > 0) {
		__m128d sum0 = _mm_setzero_pd();
		__m128d sum1 = _mm_setzero_pd();
		for (; i < limit; i += 4) {
			__m128d sparse1 = _mm_load_sd(&a[idx[i]]);
			__m128d sparse2 = _mm_load_sd(&a[idx[i+1]]);
			__m128d sparse3 = _mm_load_sd(&a[idx[i+2]]);
			__m128d sparse4 = _mm_load_sd(&a[idx[i+3]]);
			__m128d sparse_vals0 = _mm_unpacklo_pd(sparse1, sparse2);
			__m128d sparse_vals1 = _mm_unpacklo_pd(sparse3, sparse4);
			__m128d dense_vals0 = _mm_load_pd(&v[i]);
			__m128d dense_vals1 = _mm_load_pd(&v[i + 2]);
			__m128d prod0 = _mm_mul_pd(sparse_vals0, dense_vals0);
			__m128d prod1 = _mm_mul_pd(sparse_vals1, dense_vals1);
			sum0 = _mm_add_pd(sum0, prod0);
			sum1 = _mm_add_pd(sum1, prod1);
		}
		sum0 = _mm_add_pd(sum0, sum1);
		__m128d sum_swapped = _mm_shuffle_pd(sum0, sum0, 1);
		__m128d sum_total = _mm_add_pd(sum0, sum_swapped);
		_mm_store_sd(&result, sum_total);
	}
	for (; i < n; i++) {
		result += a[idx[i]] * v[i];
	}
	return result;
