	int  i = 0;
	for (; i <= n - 4; i += 4) {
		__m128i vindex = _mm_loadu_si128((const __m128i*)(ia + i));
		__m256d gathered = _mm256_i32gather_pd(a, vindex, 8);
		_mm256_storeu_pd(y + i, gathered);
	}
	for (; i < n; i++) {
		y[i] = a[ia[i]];
	}
