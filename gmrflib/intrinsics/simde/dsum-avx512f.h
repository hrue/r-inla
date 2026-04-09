{
	if (n >= 32) {
		__m512d s0 = _mm512_setzero_pd();
		__m512d s1 = _mm512_setzero_pd();
		__m512d s2 = _mm512_setzero_pd();
		__m512d s3 = _mm512_setzero_pd();
		int i = 0;
		for (; i <= n - 32; i += 32) {
			s0 = _mm512_add_pd(s0, _mm512_loadu_pd(&x[i]));
			s1 = _mm512_add_pd(s1, _mm512_loadu_pd(&x[i + 8]));
			s2 = _mm512_add_pd(s2, _mm512_loadu_pd(&x[i + 16]));
			s3 = _mm512_add_pd(s3, _mm512_loadu_pd(&x[i + 24]));
		}
		__m512d final_sum = _mm512_add_pd(_mm512_add_pd(s0, s1), 
						  _mm512_add_pd(s2, s3));
		__m256d low256 = _mm256_castpd512_pd256(final_sum);
		__m256d high256 = _mm256_extractf128_pd(final_sum, 1); 
		__m256d sum256 = _mm256_add_pd(low256, _mm256_extractf128_pd(_mm256_castpd512_pd256(_mm512_alignr_epi64(final_sum, final_sum, 32)), 1));
		double temp[8];
		_mm512_storeu_pd(temp, final_sum);
		double r = 0;
		for(int k=0; k<8; k++) {
			r += temp[k];
		}
		for (; i < n; i++) {
			r += x[i];
		}
		return r;
	} else {
#include "dsum-sse2.h"
	}

#if 0
	double r = 0.0;
	int limit = n & ~7;  // Align to 8 doubles
	if (limit > 0) {
		simde__m512d sum0 = simde_mm512_setzero_pd();
		for (int i = 0; i < limit; i += 8) {
			simde__m512d data = simde_mm512_loadu_pd(&x[i]);
			sum0 = simde_mm512_add_pd(sum0, data);
		}
		r += _mm512_reduce_add_pd(sum0);
	}
	for (int i = limit; i < n; i++) {
		r += x[i];
	}
	return r;
#endif
