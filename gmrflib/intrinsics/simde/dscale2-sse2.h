{
	__m128d va = simde_mm_set1_pd(a);
	int i = 0;
	for (; i <= n - 16; i += 16) {
		__m128d x0 = simde_mm_loadu_pd(&x[i]);
		__m128d x1 = simde_mm_loadu_pd(&x[i + 2]);
		__m128d x2 = simde_mm_loadu_pd(&x[i + 4]);
		__m128d x3 = simde_mm_loadu_pd(&x[i + 6]);
		__m128d x4 = simde_mm_loadu_pd(&x[i + 8]);
		__m128d x5 = simde_mm_loadu_pd(&x[i + 10]);
		__m128d x6 = simde_mm_loadu_pd(&x[i + 12]);
		__m128d x7 = simde_mm_loadu_pd(&x[i + 14]);
		simde_mm_storeu_pd(&y[i],      simde_mm_mul_pd(x0, va));
		simde_mm_storeu_pd(&y[i + 2],  simde_mm_mul_pd(x1, va));
		simde_mm_storeu_pd(&y[i + 4],  simde_mm_mul_pd(x2, va));
		simde_mm_storeu_pd(&y[i + 6],  simde_mm_mul_pd(x3, va));
		simde_mm_storeu_pd(&y[i + 8],  simde_mm_mul_pd(x4, va));
		simde_mm_storeu_pd(&y[i + 10], simde_mm_mul_pd(x5, va));
		simde_mm_storeu_pd(&y[i + 12], simde_mm_mul_pd(x6, va));
		simde_mm_storeu_pd(&y[i + 14], simde_mm_mul_pd(x7, va));
	}
	for (; i <= n - 2; i += 2) {
		simde_mm_storeu_pd(&y[i], simde_mm_mul_pd(simde_mm_loadu_pd(&x[i]), va));
	}
	for (; i < n; i++) {
		y[i] = x[i] * a;
	}
}
