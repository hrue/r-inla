{
	alignas(32) double total_sum = r0;
	alignas(32) double temp[4];
	int i = 0;
	if (n >= 16) {
		simde__m256d sum0 = simde_mm256_setzero_pd();
		simde__m256d sum1 = simde_mm256_setzero_pd();
		simde__m256d sum2 = simde_mm256_setzero_pd();
		simde__m256d sum3 = simde_mm256_setzero_pd();
		for (; i <= n - 16; i += 16) {
			sum0 = simde_mm256_add_pd(sum0, simde_mm256_loadu_pd(&x[i]));
			sum1 = simde_mm256_add_pd(sum1, simde_mm256_loadu_pd(&x[i + 4]));
			sum2 = simde_mm256_add_pd(sum2, simde_mm256_loadu_pd(&x[i + 8]));
			sum3 = simde_mm256_add_pd(sum3, simde_mm256_loadu_pd(&x[i + 12]));
		}
		simde__m256d sum01 = simde_mm256_add_pd(sum0, sum1);
		simde__m256d sum23 = simde_mm256_add_pd(sum2, sum3);
		simde__m256d final_sum_reg = simde_mm256_add_pd(sum01, sum23);
		simde_mm256_store_pd(temp, final_sum_reg);
		total_sum += temp[0] + temp[1] + temp[2] + temp[3];
	}
	if (n - i >= 4) {
		simde__m256d rem_sum = simde_mm256_setzero_pd();
		for (; i <= n - 4; i += 4) {
			rem_sum = simde_mm256_add_pd(rem_sum, simde_mm256_loadu_pd(&x[i]));
		}
		simde_mm256_store_pd(temp, rem_sum);
		total_sum += temp[0] + temp[1] + temp[2] + temp[3];
	}
	if (n - i >= 2) {
		simde__m128d sse_sum = simde_mm_loadu_pd(&x[i]);
		simde_mm_store_pd(temp, sse_sum);
		total_sum += temp[0] + temp[1];
		i += 2;
	}
	if (i < n) {
		total_sum += x[i];
	}
	return total_sum;
}
