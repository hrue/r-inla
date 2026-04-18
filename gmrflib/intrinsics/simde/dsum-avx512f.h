{
	double alignas(64) total_sum = r0;
	double alignas(64) res[2];
	int i = 0;
	if (n >= 32) {
		simde__m512d sum0 = simde_mm512_setzero_pd();
		simde__m512d sum1 = simde_mm512_setzero_pd();
		simde__m512d sum2 = simde_mm512_setzero_pd();
		simde__m512d sum3 = simde_mm512_setzero_pd();
		for (; i <= n - 32; i += 32) {
			sum0 = simde_mm512_add_pd(sum0, simde_mm512_loadu_pd(&x[i]));
			sum1 = simde_mm512_add_pd(sum1, simde_mm512_loadu_pd(&x[i + 8]));
			sum2 = simde_mm512_add_pd(sum2, simde_mm512_loadu_pd(&x[i + 16]));
			sum3 = simde_mm512_add_pd(sum3, simde_mm512_loadu_pd(&x[i + 24]));
		}
		simde__m512d combined = simde_mm512_add_pd(simde_mm512_add_pd(sum0, sum1), simde_mm512_add_pd(sum2, sum3));
		simde__m256d low256 = simde_mm512_extractf64x4_pd(combined, 0);
		simde__m256d high256 = simde_mm512_extractf64x4_pd(combined, 1);
		simde__m256d sum256 = simde_mm256_add_pd(low256, high256);
		simde__m128d low128 = simde_mm256_extractf128_pd(sum256, 0);
		simde__m128d high128 = simde_mm256_extractf128_pd(sum256, 1);
		simde__m128d sum128 = simde_mm_add_pd(low128, high128);
		simde_mm_store_pd(res, sum128);
		total_sum += res[0] + res[1];
	}
	if (n - i >= 8) {
		simde__m512d rem_sum = simde_mm512_loadu_pd(&x[i]);
		simde__m256d low256 = simde_mm512_extractf64x4_pd(rem_sum, 0);
		simde__m256d high256 = simde_mm512_extractf64x4_pd(rem_sum, 1);
		simde__m256d sum256 = simde_mm256_add_pd(low256, high256);
		simde__m128d low128 = simde_mm256_extractf128_pd(sum256, 0);
		simde__m128d high128 = simde_mm256_extractf128_pd(sum256, 1);
		simde__m128d sum128 = simde_mm_add_pd(low128, high128);
		simde_mm_store_pd(res, sum128);
		total_sum += res[0] + res[1];
		i += 8;
	}
	if (n - i >= 4) {
		simde__m256d rem_sum = simde_mm256_loadu_pd(&x[i]);
		simde__m128d low128 = simde_mm256_extractf128_pd(rem_sum, 0);
		simde__m128d high128 = simde_mm256_extractf128_pd(rem_sum, 1);
		simde__m128d sum128 = simde_mm_add_pd(low128, high128);
		simde_mm_store_pd(res, sum128);
		total_sum += res[0] + res[1];
		i += 4;
	}
	if (n - i >= 2) {
		simde__m128d rem_sum = simde_mm_loadu_pd(&x[i]);
		simde_mm_store_pd(res, rem_sum);
		total_sum += res[0] + res[1];
		i += 2;
	}
	for (; i < n; i++) {
		total_sum += x[i];
	}
	return total_sum;
}
