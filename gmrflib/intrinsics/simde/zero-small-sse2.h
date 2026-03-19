	__m128d eps_vec = _mm_set1_pd(eps);
	int i = 0;
	int simd_n = n & ~1;
	for (; i < simd_n; i += 2) {
		__m128d v = _mm_loadu_pd(&x[i]);
		__m128d abs_v = _mm_andnot_pd(_mm_set1_pd(-0.0), v);
		__m128d mask = _mm_cmplt_pd(abs_v, eps_vec);
		_mm_storeu_pd(&x[i], _mm_andnot_pd(mask, v));
	}
	for (; i < n; i++) {
		if (fabs(x[i]) < eps) x[i] = 0;
	}
