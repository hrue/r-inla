{
	if (n < 8) {
		double sum = 0.0;
		for (int i = 0; i < n; i++) sum += x[i];
		return sum;
	}
	if (n < 32) {
		__m512d total_v = simde_mm512_setzero_pd();
		int i = 0;
		while (i < n) {
			int rem = n - i;
			if (rem >= 8) {
				total_v = simde_mm512_add_pd(total_v, simde_mm512_loadu_pd(&x[i]));
				i += 8;
			} else {
				__mmask8 mask = (__mmask8)((1ULL << rem) - 1);
				__m512d tail = simde_mm512_maskz_loadu_pd(mask, &x[i]);
				total_v = simde_mm512_add_pd(total_v, tail);
				i = n;
			}
		}
		double temp[8];
		simde_mm512_storeu_pd(temp, total_v);
		double sum = 0;
		for (int k = 0; k < 8; k++) sum += temp[k];
		return sum;
	}

	__m512d s0 = simde_mm512_setzero_pd();
	__m512d s1 = simde_mm512_setzero_pd();
	__m512d s2 = simde_mm512_setzero_pd();
	__m512d s3 = simde_mm512_setzero_pd();
	int i = 0;
	for (; i <= n - 32; i += 32) {
		s0 = simde_mm512_add_pd(s0, simde_mm512_loadu_pd(&x[i]));
		s1 = simde_mm512_add_pd(s1, simde_mm512_loadu_pd(&x[i + 8]));
		s2 = simde_mm512_add_pd(s2, simde_mm512_loadu_pd(&x[i + 16]));
		s3 = simde_mm512_add_pd(s3, simde_mm512_loadu_pd(&x[i + 24]));
	}
	__m512d final_v = simde_mm512_add_pd(simde_mm512_add_pd(s0, s1), 
					     simde_mm512_add_pd(s2, s3));
	if (i < n) {
		int rem = n - i;
		while (rem > 0) {
			if (rem >= 8) {
				final_v = simde_mm512_add_pd(final_v, simde_mm512_loadu_pd(&x[i]));
				i += 8;
				rem -= 8;
			} else {
				__mmask8 mask = (__mmask8)((1ULL << rem) - 1);
				final_v = simde_mm512_add_pd(final_v, simde_mm512_maskz_loadu_pd(mask, &x[i]));
				rem = 0;
			}
		}
	}
	double temp[8];
	simde_mm512_storeu_pd(temp, final_v);
	double total_sum = 0;
	for (int k = 0; k < 8; k++) total_sum += temp[k];
	return total_sum;
}

