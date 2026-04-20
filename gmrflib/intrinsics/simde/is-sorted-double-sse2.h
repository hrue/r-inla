{
	int i = 0;
	if (n >= 3) {
		for (; i <= n - 3; i += 2) {
			simde__m128d current;
			memcpy(&current, &a[i], sizeof(current));
			simde__m128d next;
			memcpy(&next, &a[i + 1], sizeof(next));
			simde__m128d comparison = simde_mm_cmpgt_pd(current, next);
			if (simde_mm_movemask_pd(comparison) != 0) {
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
