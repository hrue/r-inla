{
	if (n < 4) {
		double sum = 0.0;
		for (int i = 0; i < n; i++) sum += x[i];
		return sum;
	}
	if (n < 32) {
		__m256d sum_vec = simde_mm256_setzero_pd();
		int i = 0;
		for (; i <= n - 4; i += 4) {
			sum_vec = simde_mm256_add_pd(sum_vec, simde_mm256_loadu_pd(&x[i]));
		}
		__m128d low = simde_mm256_castpd256_pd128(sum_vec);
		__m128d high = simde_mm256_extractf128_pd(sum_vec, 1);
		__m128d combined = simde_mm_add_pd(low, high);
		__m128d shuffled = simde_mm_unpackhi_pd(combined, combined);
		double total_sum = simde_mm_cvtsd_f64(simde_mm_add_pd(combined, shuffled));
		for (; i < n; i++)
			total_sum += x[i];
		return total_sum;
	}
	__m256d s0 = simde_mm256_setzero_pd();
	__m256d s1 = simde_mm256_setzero_pd();
	__m256d s2 = simde_mm256_setzero_pd();
	__m256d s3 = simde_mm256_setzero_pd();
	int i = 0;
	for (; i <= n - 16; i += 16) {
		s0 = simde_mm256_add_pd(s0, simde_mm256_loadu_pd(&x[i]));
		s1 = simde_mm256_add_pd(s1, simde_mm256_loadu_pd(&x[i + 4]));
		s2 = simde_mm256_add_pd(s2, simde_mm256_loadu_pd(&x[i + 8]));
		s3 = simde_mm256_add_pd(s3, simde_mm256_loadu_pd(&x[i + 12]));
	}
	__m256d final_v = simde_mm256_add_pd(simde_mm256_add_pd(s0, s1), 
					simde_mm256_add_pd(s2, s3));
	__m128d low = simde_mm256_castpd256_pd128(final_v);
	__m128d high = simde_mm256_extractf128_pd(final_v, 1);
	__m128d combined = simde_mm_add_pd(low, high);
	__m128d shuffled = simde_mm_unpackhi_pd(combined, combined);
	double total_sum = simde_mm_cvtsd_f64(simde_mm_add_pd(combined, shuffled));
	for (; i < n; i++)
		total_sum += x[i];
	return total_sum;
}
