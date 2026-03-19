	__m512d eps_vec = _mm512_set1_pd(eps);
	int i = 0;
	int simd_n = n & ~7;
	for (; i < simd_n; i += 8) {
		__m512d v = _mm512_loadu_pd(&x[i]);
		__m512d abs_v = _mm512_andnot_pd(_mm512_set1_pd(-0.0), v);
		__mmask8 mask = _mm512_cmp_pd_mask(abs_v, eps_vec, _CMP_LT_OS);
		_mm512_mask_storeu_pd(&x[i], mask, _mm512_setzero_pd());
	}
	for (; i < n; i++) {
		if (fabs(x[i]) < eps) x[i] = 0;
	}
