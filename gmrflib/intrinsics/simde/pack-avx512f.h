{
	int i = 0;
	for (; i <= n - 16; i += 16) {
		simde__m256i v_idx0 = simde_mm256_loadu_si256((const simde__m256i*)&ia[i]);
		simde__m256i v_idx1 = simde_mm256_loadu_si256((const simde__m256i*)&ia[i + 8]);
		simde__m512d v_val0 = simde_mm512_i32gather_pd(v_idx0, a, 8);
		simde__m512d v_val1 = simde_mm512_i32gather_pd(v_idx1, a, 8);
		simde_mm512_storeu_pd(&y[i], v_val0);
		simde_mm512_storeu_pd(&y[i + 8], v_val1);
	}
	for (; i < n; i++) {
		y[i] = a[ia[i]];
	}
}
