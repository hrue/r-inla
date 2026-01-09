// disabled: not yet tested
#if 0 && defined(__ARM_FEATURE_SVE) && defined(__ARM_NEON_SVE_BRIDGE)
	size_t i = 0;
	while (i < n) {
		svbool_t pg = svwhilelt_b64(i, n);
		svuint32_t idx32 = svld1_u32(pg, (uint32_t*)(ia + i));
		svuint64_t idx = svunziplo_u64(svreinterpret_u64_u32(idx32), svdup_n_u64(0));
		svuint64_t offset = svmul_n_u64_z(pg, idx, 8);
		svfloat64_t gathered = svld1_gather_u64offset_f64(pg, a, offset);
		svst1_f64(pg, y + i, gathered);
		i += svcntd();
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
