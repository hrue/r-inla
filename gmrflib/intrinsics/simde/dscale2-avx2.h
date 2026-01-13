	simde__m256d scalar = simde_mm256_set1_pd(a);
	int limit = n & ~7;
	for (int i = 0; i < limit; i += 8) {
		simde__m256d yvec1 = simde_mm256_loadu_pd(&x[i]);
		simde__m256d xvec1 = simde_mm256_mul_pd(yvec1, scalar);
		simde_mm256_storeu_pd(&y[i], xvec1);
		simde__m256d yvec2 = simde_mm256_loadu_pd(&x[i + 4]);
		simde__m256d xvec2 = simde_mm256_mul_pd(yvec2, scalar);
		simde_mm256_storeu_pd(&y[i + 4], xvec2);
	}
	for (int i = limit; i < n; i++) {
		y[i] = x[i] * a;
	}
