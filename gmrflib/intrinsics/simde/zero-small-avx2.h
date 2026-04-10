{
	__m256d eps_vec = _mm256_set1_pd(eps);
	int i = 0;
	int simd_n = n & ~3;
	for (; i < simd_n; i += 4) {
		__m256d v = _mm256_loadu_pd(&x[i]);
		__m256d abs_v = _mm256_andnot_pd(_mm256_set1_pd(-0.0), v);
		__m256d mask = _mm256_cmp_pd(abs_v, eps_vec, _CMP_LT_OS);
		_mm256_storeu_pd(&x[i], _mm256_andnot_pd(mask, v));
	}
	for (; i < n; i++) {
		if (fabs(x[i]) < eps)
			x[i] = 0;
	}
}
