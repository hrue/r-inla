{
	__m512d va = simde_mm512_set1_pd(a);
	int i = 0;
	for (; i <= n - 32; i += 32) {
		__m512d x0 = simde_mm512_loadu_pd(&x[i]);
		__m512d x1 = simde_mm512_loadu_pd(&x[i + 8]);
		__m512d x2 = simde_mm512_loadu_pd(&x[i + 16]);
		__m512d x3 = simde_mm512_loadu_pd(&x[i + 24]);
		simde_mm512_storeu_pd(&y[i], simde_mm512_mul_pd(x0, va));
		simde_mm512_storeu_pd(&y[i + 8], simde_mm512_mul_pd(x1, va));
		simde_mm512_storeu_pd(&y[i + 16], simde_mm512_mul_pd(x2, va));
		simde_mm512_storeu_pd(&y[i + 24], simde_mm512_mul_pd(x3, va));
	}
	if (i < n) {
		while (i < n) {
			if (n - i >= 8) {
				simde_mm512_storeu_pd(&y[i], simde_mm512_mul_pd(simde_mm512_loadu_pd(&x[i]), va));
				i += 8;
			} else {
				__mmask8 final_mask = (__mmask8) ((1U << (n - i)) - 1);
				__m512d last_x = simde_mm512_loadu_pd(&x[i]);
				__m512d last_res = simde_mm512_mul_pd(last_x, va);
				simde_mm512_mask_storeu_pd(&y[i], final_mask, last_res);
				i = n;
			}
		}
	}
}

#if 0
{
	simde__m512d scalar = simde_mm512_set1_pd(a);
	int limit = n & ~7;
	for (int i = 0; i < limit; i += 8) {
		simde__m512d yvec = simde_mm512_loadu_pd(&x[i]);
		simde__m512d xvec = simde_mm512_mul_pd(yvec, scalar);
		simde_mm512_storeu_pd(&y[i], xvec);
	}
#       pragma omp simd
	for (int i = limit; i < n; i++) {
		y[i] = x[i] * a;
	}
}
#endif
