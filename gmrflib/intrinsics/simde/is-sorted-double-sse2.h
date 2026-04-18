{
	int i = 0;
	if (n >= 3) {
		for (; i <= n - 3; i += 2) {
			simde__m128d current = simde_mm_loadu_pd(&a[i]);
			simde__m128d next = simde_mm_loadu_pd(&a[i + 1]);
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
