{
	int i = 0;
	if (n >= 9) {
		for (; i <= n - 9; i += 8) {
			__m256i current = _mm256_loadu_si256((__m256i *) & a[i]);
			__m256i next = _mm256_loadu_si256((__m256i *) & a[i + 1]);
			__m256i comparison = _mm256_cmpgt_epi32(current, next);
			if (!_mm256_testz_si256(comparison, comparison)) {
				return 0;
			}
		}
	}
	for (; i < n - 1; i++) {
		if (a[i] > a[i + 1])
			return 0;
	}
	return 1;
}
