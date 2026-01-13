	simde__m256d scalar = simde_mm256_set1_pd(a);
	int limit = n & ~3;
	for (int i = 0; i < limit; i += 4) {
		__m256d yvec = _mm256_loadu_pd(&x[i]);
		__m256d xvec = _mm256_mul_pd(yvec, scalar);
		_mm256_storeu_pd(&y[i], xvec);
	}
	for (int i = limit; i < n; i++) {
		y[i] = x[i] * a;
	}
