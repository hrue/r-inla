{
	int i = 0;
	if (n >= 5) {
		for (; i <= n - 5; i += 4) {
			__m256d current = _mm256_loadu_pd(&a[i]);
			__m256d next = _mm256_loadu_pd(&a[i + 1]);
			__m256d mask = _mm256_cmp_pd(current, next, _CMP_GT_OQ);
			if (!_mm256_testz_pd(mask, mask)) {
				return false;
			}
		}
	}
	for (; i < n - 1; i++) {
		if (a[i] > a[i + 1]) {
			return false;
		}
	}
	return true;
}
