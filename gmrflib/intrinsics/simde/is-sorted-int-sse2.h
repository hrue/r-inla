{
	int i = 0;
	if (n >= 6) {
		for (; i <= n - 5; i += 4) {
			simde__m128i current;
			memcpy(&current, &a[i], sizeof(current));
			simde__m128i next;
			memcpy(&next, &a[i + 1], sizeof(next));
			simde__m128i comparison = simde_mm_cmpgt_epi32(current, next);
			if (simde_mm_movemask_pd(simde_mm_castsi128_pd(comparison)) != 0) {
				return false;
			}
		}
	}
	for (; i < n - 1; i++) {
		if (a[i] > a[i + 1])
			return false;
	}

	for (; i < n - 1; i++) {
		if (a[i] > a[i + 1])
			return false;
	}
	return true;
}
