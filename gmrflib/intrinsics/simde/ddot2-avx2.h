// a = ddot(x,y) && b = ddot(x,z)
// THIS IS ONLY USED FOR n=16
{
	simde__m256d sum_xy = simde_mm256_setzero_pd();
	simde__m256d sum_xz = simde_mm256_setzero_pd();

	simde__m256d vx0 = simde_mm256_loadu_pd(&x[0]);
	simde__m256d vx1 = simde_mm256_loadu_pd(&x[4]);
	simde__m256d vx2 = simde_mm256_loadu_pd(&x[8]);
	simde__m256d vx3 = simde_mm256_loadu_pd(&x[12]);

	simde__m256d vy0 = simde_mm256_loadu_pd(&y[0]);
	simde__m256d vz0 = simde_mm256_loadu_pd(&z[0]);
	sum_xy = simde_mm256_fmadd_pd(vx0, vy0, sum_xy);
	sum_xz = simde_mm256_fmadd_pd(vx0, vz0, sum_xz);

	simde__m256d vy1 = simde_mm256_loadu_pd(&y[4]);
	simde__m256d vz1 = simde_mm256_loadu_pd(&z[4]);
	sum_xy = simde_mm256_fmadd_pd(vx1, vy1, sum_xy);
	sum_xz = simde_mm256_fmadd_pd(vx1, vz1, sum_xz);

	simde__m256d vy2 = simde_mm256_loadu_pd(&y[8]);
	simde__m256d vz2 = simde_mm256_loadu_pd(&z[8]);
	sum_xy = simde_mm256_fmadd_pd(vx2, vy2, sum_xy);
	sum_xz = simde_mm256_fmadd_pd(vx2, vz2, sum_xz);

	simde__m256d vy3 = simde_mm256_loadu_pd(&y[12]);
	simde__m256d vz3 = simde_mm256_loadu_pd(&z[12]);
	sum_xy = simde_mm256_fmadd_pd(vx3, vy3, sum_xy);
	sum_xz = simde_mm256_fmadd_pd(vx3, vz3, sum_xz);

#define REDUCE_VEC(r_, v_) {						\
		simde__m128d low = simde_mm256_castpd256_pd128(v_);		\
		simde__m128d high = simde_mm256_extractf128_pd(v_, 1);		\
		simde__m128d combined = simde_mm_add_pd(low, high);		\
		simde__m128d shuffled = simde_mm_unpackhi_pd(combined, combined);	\
		r_ = simde_mm_cvtsd_f64(simde_mm_add_pd(combined, shuffled));	\
	}

	REDUCE_VEC(*a, sum_xy);
	REDUCE_VEC(*b, sum_xz);
#undef REDUCE_VEC
}
