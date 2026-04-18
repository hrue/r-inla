{
	int i = 0;
	if (n >= 5) {
		for (; i <= n - 5; i += 4) {
			simde__m256d current = simde_mm256_loadu_pd(&a[i]);
			simde__m256d next = simde_mm256_loadu_pd(&a[i + 1]);
			simde__m256d mask = simde_mm256_cmp_pd(current, next, _CMP_GT_OQ);
			if (!simde_mm256_testz_pd(mask, mask)) {
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
