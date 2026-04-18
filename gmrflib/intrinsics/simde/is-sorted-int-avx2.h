{
	int i = 0;
	if (n >= 9) {
		for (; i <= n - 9; i += 8) {
			simde__m256i current;
			memcpy(&current, &a[i], sizeof(current));
			simde__m256i next;
			memcpy(&next, &a[i + 1], sizeof(next));
			simde__m256i comparison = simde_mm256_cmpgt_epi32(current, next);
			if (!simde_mm256_testz_si256(comparison, comparison)) {
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
