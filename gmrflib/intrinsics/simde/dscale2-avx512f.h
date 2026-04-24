{
	simde__m512d va = simde_mm512_set1_pd(a);
	int i = 0;
	for (; i <= n - 32; i += 32) {
		simde__m512d x0 = simde_mm512_loadu_pd(&x[i]);
		simde__m512d x1 = simde_mm512_loadu_pd(&x[i + 8]);
		simde__m512d x2 = simde_mm512_loadu_pd(&x[i + 16]);
		simde__m512d x3 = simde_mm512_loadu_pd(&x[i + 24]);
		simde_mm512_storeu_pd(&y[i], simde_mm512_mul_pd(x0, va));
		simde_mm512_storeu_pd(&y[i + 8], simde_mm512_mul_pd(x1, va));
		simde_mm512_storeu_pd(&y[i + 16], simde_mm512_mul_pd(x2, va));
		simde_mm512_storeu_pd(&y[i + 24], simde_mm512_mul_pd(x3, va));
	}
	if (i < n) {
		int rem = n - i;
		if (rem >= 8) {
			simde_mm512_storeu_pd(&y[i], simde_mm512_mul_pd(simde_mm512_loadu_pd(&x[i]), va));
			i += 8;
		}
		if (rem > 0) {
			simde__mmask8 mask = (simde__mmask8) ((1ULL << rem) - 1);
			simde__m512d vx = simde_mm512_maskz_loadu_pd(mask, &x[i]);
			simde__m512d res = simde_mm512_mul_pd(vx, va);
			simde_mm512_mask_storeu_pd(&y[i], mask, res);
		}
	}
}
