	int limit = n & ~3;
	for (int i = 0; i < limit; i += 4) {
		simde__m128i indices1 = simde_mm_loadu_si128((simde__m128i *) & ia[i]);
		int *idx1 = (int *) &indices1;
		simde__m128d values1 = simde_mm_set_pd(a[idx1[1]], a[idx1[0]]);
		simde_mm_storeu_pd(&y[i], values1);

		simde__m128i indices2 = simde_mm_loadu_si128((simde__m128i *) & ia[i + 2]);
		int *idx2 = (int *) &indices2;
		simde__m128d values2 = simde_mm_set_pd(a[idx2[1]], a[idx2[0]]);
		simde_mm_storeu_pd(&y[i + 2], values2);
	}
	for (int i = limit; i < n; i++) {
		y[i] = a[ia[i]];
	}
