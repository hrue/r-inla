{
	simde__m256d va = simde_mm256_set1_pd(a);
	int i = 0;
	for (; i <= n - 16; i += 16) {
		simde__m256d x0 = simde_mm256_loadu_pd(&x[i]);
		simde__m256d x1 = simde_mm256_loadu_pd(&x[i + 4]);
		simde__m256d x2 = simde_mm256_loadu_pd(&x[i + 8]);
		simde__m256d x3 = simde_mm256_loadu_pd(&x[i + 12]);
		simde_mm256_storeu_pd(&y[i], simde_mm256_mul_pd(x0, va));
		simde_mm256_storeu_pd(&y[i + 4], simde_mm256_mul_pd(x1, va));
		simde_mm256_storeu_pd(&y[i + 8], simde_mm256_mul_pd(x2, va));
		simde_mm256_storeu_pd(&y[i + 12], simde_mm256_mul_pd(x3, va));
	}
	for (; i <= n - 4; i += 4) {
		simde__m256d vx = simde_mm256_loadu_pd(&x[i]);
		simde_mm256_storeu_pd(&y[i], simde_mm256_mul_pd(vx, va));
	}
	for (; i < n; i++) {
		y[i] = x[i] * a;
	}
}
