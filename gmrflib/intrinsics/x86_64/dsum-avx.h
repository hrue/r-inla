	double r = 0.0;
	int remaining = n % 16;
	int limit = n - remaining;
	if (limit > 0) {
		__m256d sum0 = _mm256_setzero_pd();
		__m256d sum1 = _mm256_setzero_pd();
		__m256d sum2 = _mm256_setzero_pd();
		__m256d sum3 = _mm256_setzero_pd();
		for (int i = 0; i < limit; i += 16) {
			sum0 = _mm256_add_pd(sum0, _mm256_loadu_pd(x + i));
			sum1 = _mm256_add_pd(sum1, _mm256_loadu_pd(x + i + 4));
			sum2 = _mm256_add_pd(sum2, _mm256_loadu_pd(x + i + 8));
			sum3 = _mm256_add_pd(sum3, _mm256_loadu_pd(x + i + 12));
		}
		sum0 = _mm256_add_pd(sum0, sum1);
		sum0 = _mm256_add_pd(sum0, sum2);
		sum0 = _mm256_add_pd(sum0, sum3);
		double result[4];
		_mm256_storeu_pd(result, sum0);
		r = result[0] + result[1] + result[2] + result[3];
	}
	for (int i = limit; i < n; i++) {
		r += x[i];
	}
	return r + r0;
