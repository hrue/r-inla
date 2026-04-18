{
	int i = 0;
	if (n >= 9) {
		for (; i <= n - 9; i += 8) {
			simde__m512d current = simde_mm512_loadu_pd(&a[i]);
			simde__m512d next = simde_mm512_loadu_pd(&a[i + 1]);
			simde__mmask8 is_unsorted = simde_mm512_cmp_pd_mask(current, next, _CMP_GT_OQ);
			if (is_unsorted != 0) {
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
