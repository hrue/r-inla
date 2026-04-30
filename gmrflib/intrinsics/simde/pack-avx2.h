{
	int i = 0;
	for (; i + 3 < n; i += 4) {
		simde__m128i v_ia = simde_mm_loadu_si128((const simde__m128i*)&ia[i]);
		simde__m256d v_a = simde_mm256_i32gather_pd(a, v_ia, 8);
		simde_mm256_storeu_pd(&y[i], v_a);
	}
	for (; i < n; i++) {
		y[i] = a[ia[i]];
	}
}
