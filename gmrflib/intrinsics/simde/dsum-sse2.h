{
	double alignas(16) total_sum = r0;
	int limit = n & ~7;
	if (limit > 0) {
		simde__m128d sum0 = simde_mm_setzero_pd();
		simde__m128d sum1 = simde_mm_setzero_pd();
		simde__m128d sum2 = simde_mm_setzero_pd();
		simde__m128d sum3 = simde_mm_setzero_pd();
		for (int i = 0; i < limit; i += 8) {
			simde__m128d data0 = simde_mm_loadu_pd(&x[i]);
			simde__m128d data1 = simde_mm_loadu_pd(&x[i + 2]);
			simde__m128d data2 = simde_mm_loadu_pd(&x[i + 4]);
			simde__m128d data3 = simde_mm_loadu_pd(&x[i + 6]);
			sum0 = simde_mm_add_pd(sum0, data0);
			sum1 = simde_mm_add_pd(sum1, data1);
			sum2 = simde_mm_add_pd(sum2, data2);
			sum3 = simde_mm_add_pd(sum3, data3);
		}
		sum0 = simde_mm_add_pd(sum0, sum1);
		sum2 = simde_mm_add_pd(sum2, sum3);
		sum0 = simde_mm_add_pd(sum0, sum2);
		simde__m128d sum_swapped = simde_mm_shuffle_pd(sum0, sum0, 1);
		simde__m128d sum_total = simde_mm_add_pd(sum0, sum_swapped);
		double alignas(16) tmp_sum;
		simde_mm_store_sd(&tmp_sum, sum_total);
		total_sum += tmp_sum;
	}
	for (int i = limit; i < n; i++) {
		total_sum += x[i];
	}
	total_sum += r0;
	return total_sum;
}
