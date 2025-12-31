        int i = 0;
	__m512d sum_vec = _mm512_setzero_pd();
	for (; i + 8 <= n; i += 8) {
		__m256i idx_vec = _mm256_loadu_si256((__m256i*)&idx[i]);
		__m512d sparse_vec = _mm512_i32gather_pd(idx_vec, a, 8);
		__m512d dense_vec = _mm512_loadu_pd(&v[i]);
		__m512d prod_vec = _mm512_mul_pd(dense_vec, sparse_vec);
		sum_vec = _mm512_add_pd(sum_vec, prod_vec);
	}
	double result = _mm512_reduce_add_pd(sum_vec);
	for (; i < n; ++i) {
		result += v[i] * a[idx[i]];
	}
	return result;
