{
	int i = 0;
	for (; i <= n - 16; i += 16) {
		simde__m128i idx0 = simde_mm_loadu_si128((const simde__m128i*)&ia[i]);
		simde__m128i idx1 = simde_mm_loadu_si128((const simde__m128i*)&ia[i + 4]);
		simde__m128i idx2 = simde_mm_loadu_si128((const simde__m128i*)&ia[i + 8]);
		simde__m128i idx3 = simde_mm_loadu_si128((const simde__m128i*)&ia[i + 12]);
		simde_mm_prefetch((const char*)&ia[i + 64], _MM_HINT_T0);
		simde__m256d v0 = simde_mm256_i32gather_pd(a, idx0, 8);
		simde__m256d v1 = simde_mm256_i32gather_pd(a, idx1, 8);
		simde__m256d v2 = simde_mm256_i32gather_pd(a, idx2, 8);
		simde__m256d v3 = simde_mm256_i32gather_pd(a, idx3, 8);
		simde_mm256_storeu_pd(&y[i], v0);
		simde_mm256_storeu_pd(&y[i + 4], v1);
		simde_mm256_storeu_pd(&y[i + 8], v2);
		simde_mm256_storeu_pd(&y[i + 12], v3);
	}
	for (; i < n; i++) y[i] = a[ia[i]];
}

