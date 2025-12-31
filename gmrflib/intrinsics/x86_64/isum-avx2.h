	int i = 0;
	int m = n & ~31;
	int r = 0;
	if (m > 0) {
		__m256i sum0 = _mm256_setzero_si256();
		__m256i sum1 = _mm256_setzero_si256();
		__m256i sum2 = _mm256_setzero_si256();
		__m256i sum3 = _mm256_setzero_si256();
		for (; i < m; i += 32) {
			__m256i vec0 = _mm256_loadu_si256((__m256i *) & x[i]);
			__m256i vec1 = _mm256_loadu_si256((__m256i *) & x[i + 8]);
			__m256i vec2 = _mm256_loadu_si256((__m256i *) & x[i + 16]);
			__m256i vec3 = _mm256_loadu_si256((__m256i *) & x[i + 24]);
			sum0 = _mm256_add_epi32(sum0, vec0);
			sum1 = _mm256_add_epi32(sum1, vec1);
			sum2 = _mm256_add_epi32(sum2, vec2);
			sum3 = _mm256_add_epi32(sum3, vec3);
		}
		sum0 = _mm256_add_epi32(sum0, sum1);
		sum2 = _mm256_add_epi32(sum2, sum3);
		sum0 = _mm256_add_epi32(sum0, sum2);
		__m128i sum128 = _mm256_extracti128_si256(sum0, 0);
		__m128i sum128_high = _mm256_extracti128_si256(sum0, 1);
		sum128 = _mm_add_epi32(sum128, sum128_high);
		sum128 = _mm_hadd_epi32(sum128, sum128);
		sum128 = _mm_hadd_epi32(sum128, sum128);
		r = _mm_extract_epi32(sum128, 0);
	}
#pragma omp simd reduction(+: r)
	for (int ii = i; ii < n; ii++) {
		r += x[ii];
	}
	return r + r0;
