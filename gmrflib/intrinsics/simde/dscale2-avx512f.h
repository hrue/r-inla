	simde__m512d scalar = simde_mm512_set1_pd(a);
	int limit = n & ~7;
	for (int i = 0; i < limit; i += 8) {
		simde__m512d yvec = simde_mm512_loadu_pd(&x[i]);
		simde__m512d xvec = simde_mm512_mul_pd(yvec, scalar);
		simde_mm512_storeu_pd(&y[i], xvec);
	}
	for (int i = limit; i < n; i++) {
		y[i] = x[i] * a;
	}
