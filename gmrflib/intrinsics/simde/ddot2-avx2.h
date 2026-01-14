	double aa = 0.0, bb = 0.0;
	int limit = n & ~3;
	if (limit > 0) {
		simde__m256d sum_a = simde_mm256_setzero_pd();
		simde__m256d sum_b = simde_mm256_setzero_pd();
		for (int i = 0; i < limit; i += 4) {
			simde__m256d xvec = simde_mm256_loadu_pd(&x[i]);
			simde__m256d yvec = simde_mm256_loadu_pd(&y[i]);
			simde__m256d zvec = simde_mm256_loadu_pd(&z[i]);
			sum_a = simde_mm256_add_pd(sum_a, simde_mm256_mul_pd(xvec, yvec));
			sum_b = simde_mm256_add_pd(sum_b, simde_mm256_mul_pd(xvec, zvec));
		}
		sum_a = simde_mm256_hadd_pd(sum_a, sum_a);  // [a0+a1, a0+a1, a2+a3, a2+a3]
		simde__m128d low_a = simde_mm256_castpd256_pd128(sum_a);
		simde__m128d high_a = simde_mm256_extractf128_pd(sum_a, 1);
		low_a = simde_mm_add_pd(low_a, high_a);  // [total, total]
		aa += simde_mm_cvtsd_f64(low_a);
		sum_b = simde_mm256_hadd_pd(sum_b, sum_b);
		simde__m128d low_b = simde_mm256_castpd256_pd128(sum_b);
		simde__m128d high_b = simde_mm256_extractf128_pd(sum_b, 1);
		low_b = simde_mm_add_pd(low_b, high_b);
		bb += simde_mm_cvtsd_f64(low_b);
	}
	for (int i = limit; i < n; i++) {
		aa += x[i] * y[i];
		bb += x[i] * z[i];
	}
	*a = aa;
	*b = bb;
