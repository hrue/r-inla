{
	int i = 0;
	for (; i + 7 < n; i += 8) {
		simde__m256i v_ia = simde_mm256_loadu_si256((constsimde__m256i*)&ia[i]);
		simde__m512d v_a = simde_mm512_i32gather_pd(a, v_ia, 8);
		simde_mm512_storeu_pd(&y[i], v_a);
	}
	if (i < n) {
		simde__mmask8 mask = (simde__mmask8)((1 << (n - i)) - 1);
		simde__m256i v_ia = simde_mm256_maskz_loadu_epi32(0, mask, &ia[i]);
		simde__m512d v_a = simde_mm512_mask_i32gather_pd(simde_mm512_setzero_pd(), mask, v_ia, a, 8);
		simde_mm512_mask_storeu_pd(&y[i], mask, v_a);
	}
}
