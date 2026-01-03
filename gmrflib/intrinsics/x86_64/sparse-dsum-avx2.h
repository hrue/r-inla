	int i = 0;
	int m = n & ~15;
	double r = 0.0;
	if (m > 0) {
		__m256d sum0 = _mm256_setzero_pd();
		__m256d sum1 = _mm256_setzero_pd();
		__m256d sum2 = _mm256_setzero_pd();
		__m256d sum3 = _mm256_setzero_pd();
		for (; i < m; i += 16) {
			__m128i idx0 = _mm_loadu_si128((__m128i *) & idx[i]);
			__m128i idx1 = _mm_loadu_si128((__m128i *) & idx[i + 4]);
			__m128i idx2 = _mm_loadu_si128((__m128i *) & idx[i + 8]);
			__m128i idx3 = _mm_loadu_si128((__m128i *) & idx[i + 12]);
			__m256d vec0 = _mm256_i32gather_pd(a, idx0, 8);
			__m256d vec1 = _mm256_i32gather_pd(a, idx1, 8);
			__m256d vec2 = _mm256_i32gather_pd(a, idx2, 8);
			__m256d vec3 = _mm256_i32gather_pd(a, idx3, 8);
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
		r += a[idx[ii]];
	}
	return r;
