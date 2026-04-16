{
	double alignas(16) total_sum = 0.0;
	int i = 0;
	if (likely(n >= 2)) {
		__m128d sum_reg = simde_mm_setzero_pd();
		for (; i <= n - 2; i += 2) {
			sum_reg = simde_mm_add_pd(sum_reg, simde_mm_loadu_pd(&x[i]));
		}
		double alignas(16) res[2];
		simde_mm_storeu_pd(res, sum_reg);
		total_sum = res[0] + res[1];
	}
	for (; i < n; i++) {
		total_sum += x[i];
	}
	return total_sum;
}
