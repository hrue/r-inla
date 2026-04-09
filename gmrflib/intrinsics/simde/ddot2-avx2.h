#if 1
	__m256d sum_xy0 = simde_mm256_setzero_pd();
	__m256d sum_xy1 = simde_mm256_setzero_pd();
	__m256d sum_xy2 = simde_mm256_setzero_pd();
	__m256d sum_xy3 = simde_mm256_setzero_pd();

	__m256d sum_xz0 = simde_mm256_setzero_pd();
	__m256d sum_xz1 = simde_mm256_setzero_pd();
	__m256d sum_xz2 = simde_mm256_setzero_pd();
	__m256d sum_xz3 = simde_mm256_setzero_pd();

	int i = 0;
	for (; i <= n - 16; i += 16) {
		__m256d vx0 = simde_mm256_loadu_pd(&x[i]);
		__m256d vx1 = simde_mm256_loadu_pd(&x[i + 4]);
		__m256d vx2 = simde_mm256_loadu_pd(&x[i + 8]);
		__m256d vx3 = simde_mm256_loadu_pd(&x[i + 12]);

		__m256d vy0 = simde_mm256_loadu_pd(&y[i]);
		__m256d vy1 = simde_mm256_loadu_pd(&y[i + 4]);
		__m256d vy2 = simde_mm256_loadu_pd(&y[i + 8]);
		__m256d vy3 = simde_mm256_loadu_pd(&y[i + 12]);

		__m256d vz0 = simde_mm256_loadu_pd(&z[i]);
		__m256d vz1 = simde_mm256_loadu_pd(&z[i + 4]);
		__m256d vz2 = simde_mm256_loadu_pd(&z[i + 8]);
		__m256d vz3 = simde_mm256_loadu_pd(&z[i + 12]);

		sum_xy0 = simde_mm256_add_pd(sum_xy0, simde_mm256_mul_pd(vx0, vy0));
		sum_xy1 = simde_mm256_add_pd(sum_xy1, simde_mm256_mul_pd(vx1, vy1));
		sum_xy2 = simde_mm256_add_pd(sum_xy2, simde_mm256_mul_pd(vx2, vy2));
		sum_xy3 = simde_mm256_add_pd(sum_xy3, simde_mm256_mul_pd(vx3, vy3));

		sum_xz0 = simde_mm256_add_pd(sum_xz0, simde_mm256_mul_pd(vx0, vz0));
		sum_xz1 = simde_mm256_add_pd(sum_xz1, simde_mm256_mul_pd(vx1, vz1));
		sum_xz2 = simde_mm256_add_pd(sum_xz2, simde_mm256_mul_pd(vx2, vz2));
		sum_xz3 = simde_mm256_add_pd(sum_xz3, simde_mm256_mul_pd(vx3, vz3));
	}

	__m256d f_xy = simde_mm256_add_pd(simde_mm256_add_pd(sum_xy0, sum_xy1), 
					  simde_mm256_add_pd(sum_xy2, sum_xy3));
    
	__m256d f_xz = simde_mm256_add_pd(simde_mm256_add_pd(sum_xz0, sum_xz1), 
					  simde_mm256_add_pd(sum_xz2, sum_xz3));

#define HORIZONTAL_ADD(r_, v_) {					\
		__m128d low = simde_mm256_castpd256_pd128(v_);		\
		__m128d high = simde_mm256_extractf128_pd(v_, 1);	\
		__m128d combined = simde_mm_add_pd(low, high);		\
		__m128d shuffled = simde_mm_unpackhi_pd(combined, combined); \
		r_ = simde_mm_cvtsd_f64(simde_mm_add_pd(combined, shuffled)); \
	}								

	double aa, bb;
	HORIZONTAL_ADD(aa, f_xy);
	HORIZONTAL_ADD(bb, f_xz);
#undef HORIZONTAL_ADD

	for (; i < n; i++) {
		aa += x[i] * y[i];
		bb += x[i] * z[i];
	}
	*a = aa;
	*b = bb;
#else 
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
#endif
