{
	int i = 0;
	if (n >= 8) {
		for (; i <= n - 5; i += 4) {
			simde__m256d current;
			memcpy(&current, &a[i], sizeof(current));
			simde__m256d next;
			memcpy(&next, &a[i + 1], sizeof(next));
			simde__m256d mask = simde_mm256_cmp_pd(current, next, SIMDE_CMP_GT_OQ);
			if (!simde_mm256_testz_pd(mask, mask)) {
				return 0;
			}
		}
	}
	for (; i < n - 1; i++) {
		if (a[i] > a[i + 1])
			return false;
	}
	return true;
}
