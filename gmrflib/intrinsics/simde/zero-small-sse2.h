{
	simde__m128d eps_vec = simde_mm_set1_pd(eps);
	int i = 0;
	int simd_n = n & ~1;
	for (; i < simd_n; i += 2) {
		simde__m128d v = simde_mm_loadu_pd(&x[i]);
		simde__m128d abs_v = simde_mm_and_pd(v, simde_mm_castsi128_pd(simde_mm_set1_epi64x(0x7FFFFFFFFFFFFFFFULL)));
		simde__m128d mask = simde_mm_cmplt_pd(abs_v, eps_vec);
		simde_mm_storeu_pd(&x[i], simde_mm_andnot_pd(mask, v));
	}
	for (; i < n; i++) {
		if (fabs(x[i]) < eps)
			x[i] = 0;
	}
}
