	simde__m128d scalar = simde_mm_set1_pd(a);
	int limit = n & ~3;
	for (int i = 0; i < limit; i += 4) {
		simde__m128d xvec1 = simde_mm_loadu_pd(&x[i]);
		simde__m128d yvec1 = simde_mm_mul_pd(xvec1, scalar);
		simde_mm_storeu_pd(&y[i], yvec1);
		simde__m128d xvec2 = simde_mm_loadu_pd(&x[i + 2]);
		simde__m128d yvec2 = simde_mm_mul_pd(xvec2, scalar);
		simde_mm_storeu_pd(&y[i + 2], yvec2);
	}
	for (int i = limit; i < n; i++) {
		y[i] = a * x[i];
	}
