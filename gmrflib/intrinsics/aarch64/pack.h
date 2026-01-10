// looks like the SVE version is slower then NEON, so I disable it for the moment
#if 0 && defined(__ARM_FEATURE_SVE) && defined(__ARM_NEON_SVE_BRIDGE)
	size_t i = 0;
	const size_t vl = svcntd();
	svbool_t pg_all = svptrue_b64();
	while (i + (2 * vl) <= n) {
		svint64_t idx0 = svld1sw_s64(pg_all, ia + i);
		svint64_t idx1 = svld1sw_s64(pg_all, ia + i + vl);
		svfloat64_t val0 = svld1_gather_s64index_f64(pg_all, a, idx0);
		svfloat64_t val1 = svld1_gather_s64index_f64(pg_all, a, idx1);
		svst1_f64(pg_all, y + i,      val0);
		svst1_f64(pg_all, y + i + vl, val1);
		i += 2 * vl;
	}
	while (i < n) {
		svbool_t pg = svwhilelt_b64(i, n);
		svint64_t idx = svld1sw_s64(pg, ia + i);
		svfloat64_t val = svld1_gather_s64index_f64(pg, a, idx);
		svst1_f64(pg, y + i, val);
		i += vl;
	}
#else
	size_t i = 0;
	size_t remaining = n % 2;
	size_t limit = n - remaining;
	for (; i < limit; i += 2) {
		int idx0 = ia[i];
		int idx1 = ia[i + 1];
		float64x2_t result = vsetq_lane_f64(a[idx0], vdupq_n_f64(0.0), 0);
		result = vsetq_lane_f64(a[idx1], result, 1);
		vst1q_f64(y + i, result);
	}
	if (i < n) {
		y[i] = a[ia[i]];
	}
#endif
