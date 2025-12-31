	double *end = x + n;
	int rem = n % 8;
	double *x_end = end - rem;
	float64x2_t sum_vec0 = vdupq_n_f64(0.0);
	float64x2_t sum_vec1 = vdupq_n_f64(0.0);
	float64x2_t sum_vec2 = vdupq_n_f64(0.0);
	float64x2_t sum_vec3 = vdupq_n_f64(0.0);
	for (double *p = x; p < x_end; p += 8) {
		float64x2_t vec0 = vld1q_f64(p);
		float64x2_t vec1 = vld1q_f64(p + 2);
		float64x2_t vec2 = vld1q_f64(p + 4);
		float64x2_t vec3 = vld1q_f64(p + 6);
		sum_vec0 = vaddq_f64(sum_vec0, vec0);
		sum_vec1 = vaddq_f64(sum_vec1, vec1);
		sum_vec2 = vaddq_f64(sum_vec2, vec2);
		sum_vec3 = vaddq_f64(sum_vec3, vec3);
	}
	sum_vec0 = vaddq_f64(sum_vec0, sum_vec1);
	sum_vec2 = vaddq_f64(sum_vec2, sum_vec3);
	sum_vec0 = vaddq_f64(sum_vec0, sum_vec2);
	double r = vgetq_lane_f64(sum_vec0, 0) + vgetq_lane_f64(sum_vec0, 1);
#pragma omp simd reduction(+: r)
	for (double *p = x_end; p < end; p++) {
		r += *p;
	}
	return r + r0;
