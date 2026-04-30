{
	int i = 0;
	if (n >= 6) {
		for (; i <= n - 5; i += 4) {
			simde__m128i current, next;
			memcpy(&current, &a[i], sizeof(current));
			memcpy(&next,   &a[i + 1], sizeof(next));
			simde__m128i cmp = simde_mm_cmpgt_epi32(current, next);
			int mask = simde_mm_movemask_epi8(cmp);
			if (mask != 0)
				return 0;
		}
	}
	for (; i < n - 1; i++) {
		if (a[i] > a[i + 1])
			return 0;
	}
	return 1;
}
