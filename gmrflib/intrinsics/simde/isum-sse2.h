{
	int limit = n & ~7;
	int r = 0;
	if (limit > 0) {
		simde__m128i sum0 = simde_mm_setzero_si128();
		simde__m128i sum1 = simde_mm_setzero_si128();
		for (int i = 0; i < limit; i += 8) {
			simde__m128i data0 = simde_mm_loadu_si128((simde__m128i *) & x[i]);
			simde__m128i data1 = simde_mm_loadu_si128((simde__m128i *) & x[i + 4]);
			sum0 = simde_mm_add_epi32(sum0, data0);
			sum1 = simde_mm_add_epi32(sum1, data1);
		}
		sum0 = simde_mm_add_epi32(sum0, sum1);
		int sum_array[4];
		simde_mm_store_si128((simde__m128i *) sum_array, sum0);
		r = (sum_array[0] + sum_array[1]) + (sum_array[2] + sum_array[3]);
	}
	for (int i = limit; i < n; i++) {
		r += x[i];
	}
	return r;
}
