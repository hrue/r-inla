	int i = 0;
	__m128d sum_vec = _mm_setzero_pd();
	for (; i + 2 <= n; i += 2) {
		__m128i idx_vec = _mm_loadu_si128((__m128i*)&idx[i]);
		int32_t idx0 = _mm_cvtsi128_si32(idx_vec);
		int32_t idx1 = _mm_cvtsi128_si32(_mm_srli_si128(idx_vec, 4));
		__m128d sparse_vec = _mm_set_pd(a[idx1], a[idx0]);
		__m128d dense_vec = _mm_loadu_pd(&v[i]);
		__m128d prod_vec = _mm_mul_pd(dense_vec, sparse_vec);
		sum_vec = _mm_add_pd(sum_vec, prod_vec);
	}
	__m128d shuffle = _mm_shuffle_pd(sum_vec, sum_vec, 1);
	sum_vec = _mm_add_pd(sum_vec, shuffle);
	double result = _mm_cvtsd_f64(sum_vec);
	for (; i < n; ++i) {
		result += v[i] * a[idx[i]];
	}
	return result;
