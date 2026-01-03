	int i = 0;
	int limit = n & ~7;
	int r = 0;
	if (limit > 0) {
		__m128i sum0 = _mm_setzero_si128();
		__m128i sum1 = _mm_setzero_si128();
		for (; i < limit; i += 8) {
			__m128i data0 = _mm_load_si128((__m128i *) & x[i]);
			__m128i data1 = _mm_load_si128((__m128i *) & x[i + 4]);
			sum0 = _mm_add_epi32(sum0, data0);
			sum1 = _mm_add_epi32(sum1, data1);
		}
		sum0 = _mm_add_epi32(sum0, sum1);
		int sum_array[4];
		_mm_store_si128((__m128i *) sum_array, sum0);
		r = sum_array[0] + sum_array[1] + sum_array[2] + sum_array[3];
	}
	for (int ii = i; ii < n; ii++) {
		r += x[ii];
	}
	return r + r0;
