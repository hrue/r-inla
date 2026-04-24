{
	simde__m256d eps_vec = simde_mm256_set1_pd(eps);
	int i = 0;
	int simd_n = n & ~3;
	for (; i < simd_n; i += 4) {
		simde__m256d v = simde_mm256_loadu_pd(&x[i]);
		simde__m256d abs_v = simde_mm256_andnot_pd(simde_mm256_set1_pd(-0.0), v);
		simde__m256d mask = simde_mm256_cmp_pd(abs_v, eps_vec, _CMP_LE_OS);
		simde_mm256_storeu_pd(&x[i], simde_mm256_andnot_pd(mask, v));
	}
	for (; i < n; i++) {
		if (fabs(x[i]) <= eps)
			x[i] = 0;
	}
}
