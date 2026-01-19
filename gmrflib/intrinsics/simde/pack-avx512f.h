	int limit = n & ~7;
	for (int i = 0; i < limit; i += 8) {
		simde__m256i indices = simde_mm256_loadu_si256((simde__m256i*)&ia[i]);
		simde__m512d values = simde_mm512_i32gather_pd(a, indices, 8);
		simde_mm512_storeu_pd(&y[i], values);
	}
	for (int i = limit; i < n; i++) {
		y[i] = a[ia[i]];
	}
