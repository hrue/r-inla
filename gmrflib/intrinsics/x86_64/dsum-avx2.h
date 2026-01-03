	double r = 0.0;
	int i = 0;
	int m = n & ~15;
	if (m > 0) {
		__m256d sum0 = _mm256_setzero_pd();
		__m256d sum1 = _mm256_setzero_pd();
		__m256d sum2 = _mm256_setzero_pd();
		__m256d sum3 = _mm256_setzero_pd();
		for (; i < m; i += 16) {
			__m256d vec0 = _mm256_loadu_pd(&x[i]);
			__m256d vec1 = _mm256_loadu_pd(&x[i + 4]);
			__m256d vec2 = _mm256_loadu_pd(&x[i + 8]);
			__m256d vec3 = _mm256_loadu_pd(&x[i + 12]);
			sum0 = _mm256_add_pd(sum0, vec0);
			sum1 = _mm256_add_pd(sum1, vec1);
			sum2 = _mm256_add_pd(sum2, vec2);
			sum3 = _mm256_add_pd(sum3, vec3);
		}
		sum0 = _mm256_add_pd(sum0, sum1);
		sum2 = _mm256_add_pd(sum2, sum3);
		sum0 = _mm256_add_pd(sum0, sum2);
		__m128d sum_low = _mm256_castpd256_pd128(sum0);
		__m128d sum_high = _mm256_extractf128_pd(sum0, 1);
		__m128d sum_128 = _mm_add_pd(sum_low, sum_high);
		sum_128 = _mm_hadd_pd(sum_128, sum_128);
		r = _mm_cvtsd_f64(sum_128);
	}
	for (int ii = i; ii < n; ii++) {
		r += x[ii];
	}
	return r + r0;
