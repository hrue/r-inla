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
		r = (simde_mm_extract_epi32(sum0, 0) + simde_mm_extract_epi32(sum0, 1)) +
		    (simde_mm_extract_epi32(sum0, 2) + simde_mm_extract_epi32(sum0, 3));
	}
	for (int i = limit; i < n; i++) {
		r += x[i];
	}
	return r;
}
