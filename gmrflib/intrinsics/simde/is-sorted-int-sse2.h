{
	int i = 0;
	if (n >= 5) {
		for (; i <= n - 5; i += 4) {
			simde__m128i current;
			memcpy(&current, &a[i], sizeof(current));
			simde__m128i next;
			memcpy(&next, &a[i + 1], sizeof(next));
			simde__m128i comparison = simde_mm_cmpgt_epi32(current, next);
			if (simde_mm_movemask_epi8(comparison) != 0) {
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
