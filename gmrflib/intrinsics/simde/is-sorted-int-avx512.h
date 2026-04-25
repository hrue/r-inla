{
	int i = 0;
	if (n >= 33) {
		for (; i <= n - 17; i += 16) {
			simde__m512i current;
			memcpy(&current, &a[i], sizeof(current));
			simde__m512i next;
			memcpy(&next, &a[i + 1], sizeof(next));
			simde__mmask16 is_unsorted = simde_mm512_cmpgt_epi32_mask(current, next);
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
