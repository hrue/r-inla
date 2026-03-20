	int i = 0;
	if (n >= 17) {
		for (; i <= n - 17; i += 16) {
			__m512i current = _mm512_loadu_si512((__m512i*)&a[i]);
			__m512i next    = _mm512_loadu_si512((__m512i*)&a[i + 1]);
			__mmask16 is_unsorted = _mm512_cmpgt_epi32_mask(current, next);
			if (is_unsorted != 0) {
				return 0;
			}
		}
	}
	for (; i < n - 1; i++) {
		if (a[i] > a[i + 1])
			return 0;
	}
	return 1;
