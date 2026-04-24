{
	alignas(16) double total_sum = r0;
	int i = 0;
	if (likely(n >= 2)) {
		simde__m128d sum_reg = simde_mm_setzero_pd();
		for (; i <= n - 2; i += 2) {
			sum_reg = simde_mm_add_pd(sum_reg, simde_mm_loadu_pd(&x[i]));
		}
		alignas(16) double res[2];
		simde_mm_store_pd(res, sum_reg);
		total_sum += res[0] + res[1];
	}
	for (; i < n; i++) {
		total_sum += x[i];
	}
	return total_sum;
}
