// a = ddot(x,y) && b = ddot(x,z)

{
	double aa = 0.0, bb = 0.0;
	int limit = n & ~3;
	if (limit > 0) {
		simde__m128d sum_a = simde_mm_setzero_pd();
		simde__m128d sum_b = simde_mm_setzero_pd();
		for (int i = 0; i < limit; i += 4) {
			simde__m128d xvec1 = simde_mm_loadu_pd(&x[i]);
			simde__m128d yvec1 = simde_mm_loadu_pd(&y[i]);
			simde__m128d zvec1 = simde_mm_loadu_pd(&z[i]);
			sum_a = simde_mm_add_pd(sum_a, simde_mm_mul_pd(xvec1, yvec1));
			sum_b = simde_mm_add_pd(sum_b, simde_mm_mul_pd(xvec1, zvec1));
			simde__m128d xvec2 = simde_mm_loadu_pd(&x[i + 2]);
			simde__m128d yvec2 = simde_mm_loadu_pd(&y[i + 2]);
			simde__m128d zvec2 = simde_mm_loadu_pd(&z[i + 2]);
			sum_a = simde_mm_add_pd(sum_a, simde_mm_mul_pd(xvec2, yvec2));
			sum_b = simde_mm_add_pd(sum_b, simde_mm_mul_pd(xvec2, zvec2));
		}
		sum_a = simde_mm_add_pd(sum_a, simde_mm_shuffle_pd(sum_a, sum_a, 1));
		sum_b = simde_mm_add_pd(sum_b, simde_mm_shuffle_pd(sum_b, sum_b, 1));
		aa += simde_mm_cvtsd_f64(sum_a);
		bb += simde_mm_cvtsd_f64(sum_b);
	}
	for (int i = limit; i < n; ++i) {
		aa += x[i] * y[i];
		bb += x[i] * z[i];
	}
	*a = aa;
	*b = bb;
}
