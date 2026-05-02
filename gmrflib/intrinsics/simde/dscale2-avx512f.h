{
	simde__m512d va = simde_mm512_set1_pd(a);
	double* p_x = x;
	double* p_y = y;
	int unroll_block = 32;
	double* end_unrolled = x + (n / unroll_block) * unroll_block;
	while (p_x < end_unrolled) {
		simde__m512d x0 = simde_mm512_loadu_pd(p_x);
		simde__m512d x1 = simde_mm512_loadu_pd(p_x + 8);
		simde__m512d x2 = simde_mm512_loadu_pd(p_x + 16);
		simde__m512d x3 = simde_mm512_loadu_pd(p_x + 24);
		simde_mm512_storeu_pd(p_y,      simde_mm512_mul_pd(x0, va));
		simde_mm512_storeu_pd(p_y + 8,  simde_mm512_mul_pd(x1, va));
		simde_mm512_storeu_pd(p_y + 16, simde_mm512_mul_pd(x2, va));
		simde_mm512_storeu_pd(p_y + 24, simde_mm512_mul_pd(x3, va));
		p_x += 32;
		p_y += 32;
	}
	int rem = n - (p_x - x);
	if (rem > 0) {
		simde__mmask8 mask = (simde__mmask8) ((1ULL << rem) - 1);
		while (rem >= 8) {
			simde_mm512_storeu_pd(p_y, simde_mm512_mul_pd(simde_mm512_loadu_pd(p_x), va));
			p_x += 8;
			p_y += 8;
			rem -= 8;
		}
		if (rem > 0) {
			mask = (simde__mmask8) ((1ULL << rem) - 1);
			simde__m512d vx = simde_mm512_maskz_loadu_pd(mask, p_x);
			simde__m512d res = simde_mm512_mul_pd(vx, va);
			simde_mm512_mask_storeu_pd(p_y, mask, res);
		}
	}
}
