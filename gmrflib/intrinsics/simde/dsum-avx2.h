	double r = 0.0;
	int limit = n & ~3;  // Align to 4 doubles
	if (limit > 0) {
		simde__m256d sum0 = simde_mm256_setzero_pd();
		for (int i = 0; i < limit; i += 4) {
			simde__m256d data = simde_mm256_loadu_pd(&x[i]);
			sum0 = simde_mm256_add_pd(sum0, data);
		}
		simde__m128d sum_low = simde_mm256_castpd256_pd128(sum0);
		simde__m128d sum_high = simde_mm256_extractf128_pd(sum0, 1);
		sum_low = simde_mm_add_pd(sum_low, sum_high);
		simde__m128d sum_swapped = simde_mm_shuffle_pd(sum_low, sum_low, 1);
		simde__m128d sum_total = simde_mm_add_pd(sum_low, sum_swapped);
		simde_mm_store_sd(&r, sum_total);
	}
	for (int i = limit; i < n; i++) {
		r += x[i];
	}
	return r + r0;
