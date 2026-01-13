	simde__m256d scalar = simde_mm256_set1_pd(a);
	int limit = n & ~3;
	for (int i = 0; i < limit; i += 4) {
		simde__m256d yvec = simde_mm256_loadu_pd(&x[i]);
		simde__m256d xvec = simde_mm256_mul_pd(yvec, scalar);
		simde_mm256_storeu_pd(&y[i], xvec);
	}
	for (int i = limit; i < n; i++) {
		y[i] = x[i] * a;
	}
