#if 1
	__m512d s_xy0 = _mm512_setzero_pd();
	__m512d s_xy1 = _mm512_setzero_pd();
	__m512d s_xz0 = _mm512_setzero_pd();
	__m512d s_xz1 = _mm512_setzero_pd();

	int i = 0;
	for (; i <= n - 16; i += 16) {
		__m512d vx0 = _mm512_loadu_pd(&x[i]);
		__m512d vx1 = _mm512_loadu_pd(&x[i + 8]);
		s_xy0 = _mm512_fmadd_pd(vx0, _mm512_loadu_pd(&y[i]), s_xy0);
		s_xy1 = _mm512_fmadd_pd(vx1, _mm512_loadu_pd(&y[i + 8]), s_xy1);
		s_xz0 = _mm512_fmadd_pd(vx0, _mm512_loadu_pd(&z[i]), s_xz0);
		s_xz1 = _mm512_fmadd_pd(vx1, _mm512_loadu_pd(&z[i + 8]), s_xz1);
	}
	__m512d final_xy = _mm512_add_pd(s_xy0, s_xy1);
	__m512d final_xz = _mm512_add_pd(s_xz0, s_xz1);
	if (i < n) {
		int rem = n - i;
		if (rem >= 8) {
			__m512d vx = _mm512_loadu_pd(&x[i]);
			final_xy = _mm512_add_pd(final_xy, _mm512_mul_pd(vx, _mm512_loadu_pd(&y[i])));
			final_xz = _mm512_add_pd(final_xz, _mm512_mul_pd(vx, _mm512_loadu_pd(&z[i])));
			i += 8;
			rem -= 8;
		}
		if (rem > 0) {
			__mmask8 mask = (__mmask8)((1ULL << rem) - 1);
			__m512d vx = _mm512_maskz_loadu_pd(mask, &x[i]);
			final_xy = _mm512_add_pd(final_xy, _mm512_mul_pd(vx, _mm512_maskz_loadu_pd(mask, &y[i])));
			final_xz = _mm512_add_pd(final_xz, _mm512_mul_pd(vx, _mm512_maskz_loadu_pd(mask, &z[i])));
		}
	}

#define REDUCE_512(r_, v_) {						\
		__m256d low = _mm256_castpd512_pd256(v_);		\
		__m256d high = _mm256_extractf128_pd(_mm256_castpd512_pd256(_mm512_alignr_epi64(v_, v_, 32)), 1); \
		double temp[8];						\
		_mm512_storeu_pd(temp, v_);				\
		double s = 0;						\
		for(int k = 0; k < 8; k++) s += temp[k];		\
		r_ = s;							\
	}

	REDUCE_512(*a, final_xy);
	REDUCE_512(*b, final_xz);
#undef REDUCE_512

#else

	double aa = 0.0, bb = 0.0;
	int limit = n & ~7;
	if (limit > 0) {
		simde__m512d sum_a = simde_mm512_setzero_pd();
		simde__m512d sum_b = simde_mm512_setzero_pd();
		for (int i = 0; i < limit; i += 8) {
			simde__m512d xvec = simde_mm512_loadu_pd(&x[i]);
			simde__m512d yvec = simde_mm512_loadu_pd(&y[i]);
			simde__m512d zvec = simde_mm512_loadu_pd(&z[i]);
			sum_a = simde_mm512_add_pd(sum_a, simde_mm512_mul_pd(xvec, yvec));
			sum_b = simde_mm512_add_pd(sum_b, simde_mm512_mul_pd(xvec, zvec));
		}
		aa += _mm512_reduce_add_pd(sum_a);
		bb += _mm512_reduce_add_pd(sum_b);
	}
#pragma omp simd reduction(+: aa, bb)
	for (int i = limit; i < n; i++) {
		aa += x[i] * y[i];
		bb += x[i] * z[i];
	}
	*a = aa;
        *b = bb;
#endif
