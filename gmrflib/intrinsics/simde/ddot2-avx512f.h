// optimize this code for n=16, which is the case used

#if 1

	__m512d s_xy0 = simde_mm512_setzero_pd();
	__m512d s_xy1 = simde_mm512_setzero_pd();
	__m512d s_xz0 = simde_mm512_setzero_pd();
	__m512d s_xz1 = simde_mm512_setzero_pd();
	int i = 0;
	for (; i <= n - 16; i += 16) {
		__m512d vx0 = simde_mm512_loadu_pd(&x[i]);
		__m512d vx1 = simde_mm512_loadu_pd(&x[i + 8]);
		s_xy0 = simde_mm512_fmadd_pd(vx0, simde_mm512_loadu_pd(&y[i]), s_xy0);
		s_xy1 = simde_mm512_fmadd_pd(vx1, simde_mm512_loadu_pd(&y[i + 8]), s_xy1);
		s_xz0 = simde_mm512_fmadd_pd(vx0, simde_mm512_loadu_pd(&z[i]), s_xz0);
		s_xz1 = simde_mm512_fmadd_pd(vx1, simde_mm512_loadu_pd(&z[i + 8]), s_xz1);
	}
	__m512d final_xy = simde_mm512_add_pd(s_xy0, s_xy1);
	__m512d final_xz = simde_mm512_add_pd(s_xz0, s_xz1);
	if (i < n) {
		int rem = n - i;
		if (rem >= 8) {
			__m512d vx = simde_mm512_loadu_pd(&x[i]);
			final_xy = simde_mm512_add_pd(final_xy, simde_mm512_mul_pd(vx, simde_mm512_loadu_pd(&y[i])));
			final_xz = simde_mm512_add_pd(final_xz, simde_mm512_mul_pd(vx, simde_mm512_loadu_pd(&z[i])));
			i += 8;
			rem -= 8;
		}
		if (rem > 0) {
			__mmask8 mask = (__mmask8)((1ULL << rem) - 1);
			__m512d vx = simde_mm512_maskz_loadu_pd(mask, &x[i]);
			final_xy = simde_mm512_add_pd(final_xy, simde_mm512_mul_pd(vx, simde_mm512_maskz_loadu_pd(mask, &y[i])));
			final_xz = simde_mm512_add_pd(final_xz, simde_mm512_mul_pd(vx, simde_mm512_maskz_loadu_pd(mask, &z[i])));
		}
	}
	double buf_xy[8];
	double buf_xz[8];
	double sum_xy;
	double sum_xz;

	simde_mm512_storeu_pd(buf_xy, final_xy);
	sum_xy = buf_xy[0] + buf_xy[1] + buf_xy[2] + buf_xy[3];
	sum_xy += buf_xy[4] + buf_xy[5] + buf_xy[6] + buf_xy[7];

	simde_mm512_storeu_pd(buf_xz, final_xz);
	sum_xz = buf_xz[0] + buf_xz[1] + buf_xz[2] + buf_xz[3];
	sum_xz += buf_xz[4] + buf_xz[5] + buf_xz[6] + buf_xz[7];

	*a = sum_xy;
	*b = sum_xz;

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
