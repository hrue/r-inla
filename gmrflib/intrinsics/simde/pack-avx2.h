{
	int i = 0;
	for (; i <= n - 16; i += 16) {
		// 1. Load all indices first
		simde__m128i idx0 = simde_mm_loadu_si128((const simde__m128i*)&ia[i]);
		simde__m128i idx1 = simde_mm_loadu_si128((const simde__m128i*)&ia[i + 4]);
		simde__m128i idx2 = simde_mm_loadu_si128((const simde__m128i*)&ia[i + 8]);
		simde__m128i idx3 = simde_mm_loadu_si128((const simde__m128i*)&ia[i + 12]);

		// 2. Fire ALL gathers simultaneously. 
		// This forces the CPU's memory units to track all 16 memory lookups in parallel.
		simde__m256d v0 = simde_mm256_i32gather_pd(a, idx0, 8);
		simde__m256d v1 = simde_mm256_i32gather_pd(a, idx1, 8);
		simde__m256d v2 = simde_mm256_i32gather_pd(a, idx2, 8);
		simde__m256d v3 = simde_mm256_i32gather_pd(a, idx3, 8);

		// 3. Independent Software Prefetching for the index stream
		simde_mm_prefetch((const char*)&ia[i + 32], _MM_HINT_T0);

		// 4. Store the values only AFTER all the gathers have had time to pipeline
		simde_mm256_storeu_pd(&y[i], v0);
		simde_mm256_storeu_pd(&y[i + 4], v1);
		simde_mm256_storeu_pd(&y[i + 8], v2);
		simde_mm256_storeu_pd(&y[i + 12], v3);
	}
	for (; i < n; i++)
		y[i] = a[ia[i]];
}
