{
	if (n >= 8) {
		__m128d s0 = simde_mm_setzero_pd();
		__m128d s1 = simde_mm_setzero_pd();
		__m128d s2 = simde_mm_setzero_pd();
		__m128d s3 = simde_mm_setzero_pd();
		int i = 0;
		for (; i <= n - 8; i += 8) {
			s0 = simde_mm_add_pd(s0, simde_mm_loadu_pd(&x[i]));
			s1 = simde_mm_add_pd(s1, simde_mm_loadu_pd(&x[i + 2]));
			s2 = simde_mm_add_pd(s2, simde_mm_loadu_pd(&x[i + 4]));
			s3 = simde_mm_add_pd(s3, simde_mm_loadu_pd(&x[i + 6]));
		}
		__m128d final_sum = simde_mm_add_pd(simde_mm_add_pd(s0, s1), simde_mm_add_pd(s2, s3));
		__m128d shuffled = simde_mm_unpackhi_pd(final_sum, final_sum);
		__m128d res_vec = simde_mm_add_pd(final_sum, shuffled);
		double r = simde_mm_cvtsd_f64(res_vec);
		for (; i < n; i++) {
			r += x[i];
		}
		return r;
	} else {
		double r = 0.0;
		for (int i = 0; i < n; i++) {
			r += x[i];
		}
		return r;
	}
}

#if 0
{
	double r = 0.0;
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
		simde_mm_store_sd(&r, sum_total);
	}
	for (int i = limit; i < n; i++) {
		r += x[i];
	}
	return r;
}
#endif
