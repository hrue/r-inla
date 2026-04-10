{
	int i = 0;
	if (n >= 9) {
		for (; i <= n - 9; i += 8) {
			__m512d current = _mm512_loadu_pd(&a[i]);
			__m512d next = _mm512_loadu_pd(&a[i + 1]);
			__mmask8 is_unsorted = _mm512_cmp_pd_mask(current, next, _CMP_GT_OQ);
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
