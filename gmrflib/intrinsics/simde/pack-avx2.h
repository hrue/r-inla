	int limit = n & ~7;
	for (int i = 0; i < limit; i += 8) {
		simde__m128i indices1 = simde_mm_loadu_si128((simde__m128i*)&ia[i]);
		simde__m256d values1 = simde_mm256_i32gather_pd(a, indices1, 8);
		simde_mm256_storeu_pd(&y[i], values1);
		simde__m128i indices2 = simde_mm_loadu_si128((simde__m128i*)&ia[i + 4]);
		simde__m256d values2 = simde_mm256_i32gather_pd(a, indices2, 8);
		simde_mm256_storeu_pd(&y[i + 4], values2);
	}
	for (int i = limit; i < n; i++) {
		y[i] = a[ia[i]];
	}
