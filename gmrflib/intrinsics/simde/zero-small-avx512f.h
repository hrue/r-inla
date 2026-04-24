{
	simde__m512d eps_vec = simde_mm512_set1_pd(eps);
	int i = 0;
	int simd_n = n & ~7;
	for (; i < simd_n; i += 8) {
		simde__m512d v = simde_mm512_loadu_pd(&x[i]);
		simde__m512d abs_v = simde_mm512_and_pd(v, simde_mm512_castps_pd(simde_mm512_set1_ps(0x7FFFFFFF)));
		simde__mmask8 mask = simde_mm512_cmp_pd_mask(abs_v, eps_vec, _CMP_LT_OS);
		simde_mm512_mask_storeu_pd(&x[i], mask, simde_mm512_setzero_pd());
	}
	for (; i < n; i++) {
		if (fabs(x[i]) < eps)
			x[i] = 0;
	}
}
