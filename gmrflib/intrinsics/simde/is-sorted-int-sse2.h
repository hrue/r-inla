	int i = 0;
	if (n >= 5) {
		for (; i <= n - 5; i += 4) {
			__m128i current = _mm_loadu_si128((__m128i*)&a[i]);
			__m128i next    = _mm_loadu_si128((__m128i*)&a[i + 1]);
			__m128i comparison = _mm_cmpgt_epi32(current, next);
			if (_mm_movemask_epi8(comparison) != 0) {
				return 0;
			}
		}
	}
	for (; i < n - 1; i++) {
		if (a[i] > a[i + 1])
			return 0;
	}
	return 1;
