#if 1
	__m256d va = _mm256_set1_pd(a);
	int i = 0;
	for (; i <= n - 16; i += 16) {
		__m256d x0 = _mm256_loadu_pd(&x[i]);
		__m256d x1 = _mm256_loadu_pd(&x[i + 4]);
		__m256d x2 = _mm256_loadu_pd(&x[i + 8]);
		__m256d x3 = _mm256_loadu_pd(&x[i + 12]);
		_mm256_storeu_pd(&y[i],      _mm256_mul_pd(x0, va));
		_mm256_storeu_pd(&y[i + 4],  _mm256_mul_pd(x1, va));
		_mm256_storeu_pd(&y[i + 8],  _mm256_mul_pd(x2, va));
		_mm256_storeu_pd(&y[i + 12], _mm256_mul_pd(x3, va));
	}
	for (; i <= n - 4; i += 4) {
		__m256d vx = _mm256_loadu_pd(&x[i]);
		_mm256_storeu_pd(&y[i], _mm256_mul_pd(vx, va));
	}
	for (; i < n; i++) {
		y[i] = x[i] * a;
	}
#else
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
#endif
