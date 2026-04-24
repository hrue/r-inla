{
	int i = 0;
	if (n >= 16) {
		for (; i <= n - 9; i += 8) {
			simde__m512d current;
			memcpy(&current, &a[i], sizeof(current));
			simde__m512d next;
			memcpy(&next, &a[i + 1], sizeof(next));
			simde__mmask8 is_unsorted = simde_mm512_cmp_pd_mask(current, next, SIMDE_CMP_GT_OQ);
			if (is_unsorted != 0) {
				return false;
			}
		}
	}
	for (; i < n - 1; i++) {
		if (a[i] > a[i + 1])
			return false;
	}
	return true;
}
