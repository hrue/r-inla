	int i = 0;
	__m256d sum_vec = _mm256_setzero_pd();
	for (; i + 4 <= n; i += 4) {
		__m128i idx_vec = _mm_loadu_si128((__m128i*)&idx[i]);
		__m256d sparse_vec = _mm256_i32gather_pd(a, idx_vec, 8);
		__m256d dense_vec = _mm256_loadu_pd(&v[i]);
		__m256d prod_vec = _mm256_mul_pd(dense_vec, sparse_vec);
		sum_vec = _mm256_add_pd(sum_vec, prod_vec);
	}
	__m128d sum128 = _mm256_extractf128_pd(sum_vec, 0);
	__m128d sum256 = _mm256_extractf128_pd(sum_vec, 1);
	sum128 = _mm_add_pd(sum128, sum256);
	sum128 = _mm_hadd_pd(sum128, sum128);
	double result = _mm_cvtsd_f64(sum128);
	for (; i < n; ++i) {
		result += v[i] * a[idx[i]];
	}
	return result;
