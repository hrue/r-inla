// looks like the SVE version is slower then NEON, so I disable it for the moment
#if 0 && defined(__ARM_FEATURE_SVE) && defined(__ARM_NEON_SVE_BRIDGE)
	int i = 0;
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
	int i = 0;
	for (; i <= n - 4; i += 4) {
		int32x4_t idxs = vld1q_s32(ia + i);
		int i0 = vgetq_lane_s32(idxs, 0);
		int i1 = vgetq_lane_s32(idxs, 1);
		int i2 = vgetq_lane_s32(idxs, 2);
		int i3 = vgetq_lane_s32(idxs, 3);
		double d0 = a[i0];
		double d1 = a[i1];
		double d2 = a[i2];
		double d3 = a[i3];
		vst1q_f64(y + i, vcombine_f64(vld1_f64(&d0), vld1_f64(&d1)));
		vst1q_f64(y + i + 2, vcombine_f64(vld1_f64(&d2), vld1_f64(&d3)));
	}
	for (; i < n; i++) {
		y[i] = a[ia[i]];
	}
#endif
