{
	int i = 0;
	if (n >= 3) {
		for (; i <= n - 3; i += 2) {
			__m128d current = _mm_loadu_pd(&a[i]);
			__m128d next = _mm_loadu_pd(&a[i + 1]);
			__m128d comparison = _mm_cmpgt_pd(current, next);
			if (_mm_movemask_pd(comparison) != 0) {
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
