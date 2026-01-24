	double r = 0.0;
	int limit = n & ~7;  // Align to 8 doubles
	if (limit > 0) {
		simde__m512d sum0 = simde_mm512_setzero_pd();
		for (int i = 0; i < limit; i += 8) {
			simde__m512d data = simde_mm512_loadu_pd(&x[i]);
			sum0 = simde_mm512_add_pd(sum0, data);
		}
		r += _mm512_reduce_add_pd(sum0);
	}
	for (int i = limit; i < n; i++) {
		r += x[i];
	}
	return r;
