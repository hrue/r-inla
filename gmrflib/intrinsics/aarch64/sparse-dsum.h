	int i = 0;
	int rem = n % 4;
	int end = n - rem;
	double r = 0.0;
	if (end > 0) {
		float64x2_t sum_vec0 = vdupq_n_f64(0.0);
		float64x2_t sum_vec1 = vdupq_n_f64(0.0);
		for (; i < end; i += 4) {
			double d0[2] = { a[idx[i]], a[idx[i + 1]] };
			double d2[2] = { a[idx[i + 2]], a[idx[i + 3]] };
			float64x2_t vec0 = vld1q_f64(d0);
			float64x2_t vec1 = vld1q_f64(d2);
			sum_vec0 = vaddq_f64(sum_vec0, vec0);
			sum_vec1 = vaddq_f64(sum_vec1, vec1);
		}
		sum_vec0 = vaddq_f64(sum_vec0, sum_vec1);
		r = vgetq_lane_f64(sum_vec0, 0) + vgetq_lane_f64(sum_vec0, 1);
	}
#pragma omp simd reduction(+: r)
	for (int ii = i; ii < n; ii++) {
		r += a[idx[ii]];
	}
	return r;
