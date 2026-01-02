	int i = 0;
	for (; i <= n - 8; i += 8) {
		__m128i indices0 = _mm_loadl_epi64((__m128i*)&ia[i]);
		__m128i indices1 = _mm_loadl_epi64((__m128i*)&ia[i+2]);
		__m128i indices2 = _mm_loadl_epi64((__m128i*)&ia[i+4]);
		__m128i indices3 = _mm_loadl_epi64((__m128i*)&ia[i+6]);
		long long idx64_0 = _mm_cvtsi128_si64(indices0);
		long long idx64_1 = _mm_cvtsi128_si64(indices1);
		long long idx64_2 = _mm_cvtsi128_si64(indices2);
		long long idx64_3 = _mm_cvtsi128_si64(indices3);
		int idx0 = (int)(idx64_0 & 0xFFFFFFFF);
		int idx1 = (int)((idx64_0 >> 32) & 0xFFFFFFFF);
		int idx2 = (int)(idx64_1 & 0xFFFFFFFF);
		int idx3 = (int)((idx64_1 >> 32) & 0xFFFFFFFF);
		int idx4 = (int)(idx64_2 & 0xFFFFFFFF);
		int idx5 = (int)((idx64_2 >> 32) & 0xFFFFFFFF);
		int idx6 = (int)(idx64_3 & 0xFFFFFFFF);
		int idx7 = (int)((idx64_3 >> 32) & 0xFFFFFFFF);
		_mm_storeu_pd(&y[i],   _mm_set_pd(a[idx1],  a[idx0]));
		_mm_storeu_pd(&y[i+2], _mm_set_pd(a[idx3],  a[idx2]));
		_mm_storeu_pd(&y[i+4], _mm_set_pd(a[idx5],  a[idx4]));
		_mm_storeu_pd(&y[i+6], _mm_set_pd(a[idx7],  a[idx6]));
	}
	for (; i <= n - 2; i += 2) {
		__m128i indices = _mm_loadl_epi64((__m128i*)&ia[i]);
		long long idx64 = _mm_cvtsi128_si64(indices);
		int idx0 = (int)(idx64 & 0xFFFFFFFF);
		int idx1 = (int)((idx64 >> 32) & 0xFFFFFFFF);
		__m128d result = _mm_set_pd(a[idx1], a[idx0]);
		_mm_storeu_pd(&y[i], result);
	}
	for (; i < n; i++) {
		y[i] = a[ia[i]];
	}
