// a = ddot(x,y) && b = ddot(x,z)
// THIS IS ONLY USED FOR n=16
{
	simde__m512d vx0 = simde_mm512_loadu_pd(&x[0]);
	simde__m512d vx1 = simde_mm512_loadu_pd(&x[8]);
	simde__m512d sum_xy = simde_mm512_setzero_pd();
	sum_xy = simde_mm512_fmadd_pd(vx0, simde_mm512_loadu_pd(&y[0]), sum_xy);
	sum_xy = simde_mm512_fmadd_pd(vx1, simde_mm512_loadu_pd(&y[8]), sum_xy);
	simde__m512d sum_xz = simde_mm512_setzero_pd();
	sum_xz = simde_mm512_fmadd_pd(vx0, simde_mm512_loadu_pd(&z[0]), sum_xz);
	sum_xz = simde_mm512_fmadd_pd(vx1, simde_mm512_loadu_pd(&z[8]), sum_xz);
	alignas(64) double temp_xy[8];
	alignas(64) double temp_xz[8];
	simde_mm512_store_pd(temp_xy, sum_xy);
	simde_mm512_store_pd(temp_xz, sum_xz);
	double t_xy = 0, t_xz = 0;
	for (int i = 0; i < 8; i++) {
		t_xy += temp_xy[i];
		t_xz += temp_xz[i];
	}
	*a = t_xy;
	*b = t_xz;
}
