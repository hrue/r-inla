	int i = 0;
	for (; i <= n - 16; i += 16) {
		__m128i idx0_3 = _mm_loadu_si128((__m128i*)&ia[i]);
		__m128i idx4_7 = _mm_loadu_si128((__m128i*)&ia[i + 4]);
		__m128i idx8_11 = _mm_loadu_si128((__m128i*)&ia[i + 8]);
		__m128i idx12_15 = _mm_loadu_si128((__m128i*)&ia[i + 12]);
		__m256d y0_3 = _mm256_i32gather_pd(a, idx0_3, 8);
		__m256d y4_7 = _mm256_i32gather_pd(a, idx4_7, 8);
		__m256d y8_11 = _mm256_i32gather_pd(a, idx8_11, 8);
		__m256d y12_15 = _mm256_i32gather_pd(a, idx12_15, 8);
		_mm256_storeu_pd(&y[i], y0_3);
		_mm256_storeu_pd(&y[i + 4], y4_7);
		_mm256_storeu_pd(&y[i + 8], y8_11);
		_mm256_storeu_pd(&y[i + 12], y12_15);
	}
	for (; i < n; i++) {
		y[i] = a[ia[i]];
	}
