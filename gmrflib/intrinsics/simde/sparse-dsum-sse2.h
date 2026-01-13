	double r = 0.0;
	int limit = n & ~3;
	if (limit > 0) {
		simde__m128d sum0 = simde_mm_setzero_pd();
		simde__m128d sum1 = simde_mm_setzero_pd();
		for (int i = 0; i < limit; i += 4) {
			simde__m128d val1 = simde_mm_load_sd(&a[idx[i]]);
			simde__m128d val2 = simde_mm_load_sd(&a[idx[i + 1]]);
			simde__m128d val3 = simde_mm_load_sd(&a[idx[i + 2]]);
			simde__m128d val4 = simde_mm_load_sd(&a[idx[i + 3]]);
			simde__m128d data0 = simde_mm_unpacklo_pd(val1, val2);
			simde__m128d data1 = simde_mm_unpacklo_pd(val3, val4);
			sum0 = simde_mm_add_pd(sum0, data0);
			sum1 = simde_mm_add_pd(sum1, data1);
		}
		sum0 = simde_mm_add_pd(sum0, sum1);
		simde__m128d sum_swapped = simde_mm_shuffle_pd(sum0, sum0, 1);
		simde__m128d sum_total = simde_mm_add_pd(sum0, sum_swapped);
		simde_mm_store_sd(&r, sum_total);
	}
	for (int i = limit; i < n; i++) {
		r += a[idx[i]];
	}
	return r;
