// a = ddot(x,y) && b = ddot(x,z)

{
	__m128d vx0 = simde_mm_loadu_pd(&x[0]);
	__m128d vx1 = simde_mm_loadu_pd(&x[2]);
	__m128d vx2 = simde_mm_loadu_pd(&x[4]);
	__m128d vx3 = simde_mm_loadu_pd(&x[6]);
	__m128d vx4 = simde_mm_loadu_pd(&x[8]);
	__m128d vx5 = simde_mm_loadu_pd(&x[10]);
	__m128d vx6 = simde_mm_loadu_pd(&x[12]);
	__m128d vx7 = simde_mm_loadu_pd(&x[14]);

	__m128d sum_xy = simde_mm_setzero_pd();
	sum_xy = simde_mm_add_pd(sum_xy, simde_mm_mul_pd(vx0, simde_mm_loadu_pd(&y[0])));
	sum_xy = simde_mm_add_pd(sum_xy, simde_mm_mul_pd(vx1, simde_mm_loadu_pd(&y[2])));
	sum_xy = simde_mm_add_pd(sum_xy, simde_mm_mul_pd(vx2, simde_mm_loadu_pd(&y[4])));
	sum_xy = simde_mm_add_pd(sum_xy, simde_mm_mul_pd(vx3, simde_mm_loadu_pd(&y[6])));
	sum_xy = simde_mm_add_pd(sum_xy, simde_mm_mul_pd(vx4, simde_mm_loadu_pd(&y[8])));
	sum_xy = simde_mm_add_pd(sum_xy, simde_mm_mul_pd(vx5, simde_mm_loadu_pd(&y[10])));
	sum_xy = simde_mm_add_pd(sum_xy, simde_mm_mul_pd(vx6, simde_mm_loadu_pd(&y[12])));
	sum_xy = simde_mm_add_pd(sum_xy, simde_mm_mul_pd(vx7, simde_mm_loadu_pd(&y[14])));

	__m128d sum_xz = simde_mm_setzero_pd();
	sum_xz = simde_mm_add_pd(sum_xz, simde_mm_mul_pd(vx0, simde_mm_loadu_pd(&z[0])));
	sum_xz = simde_mm_add_pd(sum_xz, simde_mm_mul_pd(vx1, simde_mm_loadu_pd(&z[2])));
	sum_xz = simde_mm_add_pd(sum_xz, simde_mm_mul_pd(vx2, simde_mm_loadu_pd(&z[4])));
	sum_xz = simde_mm_add_pd(sum_xz, simde_mm_mul_pd(vx3, simde_mm_loadu_pd(&z[6])));
	sum_xz = simde_mm_add_pd(sum_xz, simde_mm_mul_pd(vx4, simde_mm_loadu_pd(&z[8])));
	sum_xz = simde_mm_add_pd(sum_xz, simde_mm_mul_pd(vx5, simde_mm_loadu_pd(&z[10])));
	sum_xz = simde_mm_add_pd(sum_xz, simde_mm_mul_pd(vx6, simde_mm_loadu_pd(&z[12])));
	sum_xz = simde_mm_add_pd(sum_xz, simde_mm_mul_pd(vx7, simde_mm_loadu_pd(&z[14])));
    
#define REDUCE_SSE(r_, v_) {						\
		__m128d shuffled = simde_mm_unpackhi_pd(v_, v_);	\
		__m128d res = simde_mm_add_pd(v_, shuffled);		\
		r_ = simde_mm_cvtsd_f64(res);				\
	}

	REDUCE_SSE(*a, sum_xy);
	REDUCE_SSE(*b, sum_xz);
#undef REDUCE_SSE
}
