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
	for (int i = limit; i < n; i++) {
		aa += x[i] * y[i];
		bb += x[i] * z[i];
	}
	*a = aa;
        *b = bb;
