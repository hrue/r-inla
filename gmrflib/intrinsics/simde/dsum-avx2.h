	if (n >= 16) {
		__m256d sum0 = _mm256_setzero_pd();
		__m256d sum1 = _mm256_setzero_pd();
		__m256d sum2 = _mm256_setzero_pd();
		__m256d sum3 = _mm256_setzero_pd();
		int i = 0;
		for (; i <= n - 16; i += 16) {
			sum0 = _mm256_add_pd(sum0, _mm256_loadu_pd(&x[i]));
			sum1 = _mm256_add_pd(sum1, _mm256_loadu_pd(&x[i + 4]));
			sum2 = _mm256_add_pd(sum2, _mm256_loadu_pd(&x[i + 8]));
			sum3 = _mm256_add_pd(sum3, _mm256_loadu_pd(&x[i + 12]));
		}
		__m256d final_sum = _mm256_add_pd(_mm256_add_pd(sum0, sum1), 
						  _mm256_add_pd(sum2, sum3));

		__m128d low = _mm256_castpd256_pd128(final_sum);
		__m128d high = _mm256_extractf128_pd(final_sum, 1);
		__m128d combined = _mm_add_pd(low, high);
    
		__m128d shuffled = _mm_unpackhi_pd(combined, combined);
		__m128d result_vec = _mm_add_pd(combined, shuffled);
		double r = _mm_cvtsd_f64(result_vec);
		for (; i < n; i++) {
			r += x[i];
		}
		return r;
	} else {
#include "dsum-sse2.h"
	}
	
#if 0
	double r = 0.0;
	int limit = n & ~3;  // Align to 4 doubles
	if (limit > 0) {
		simde__m256d sum0 = simde_mm256_setzero_pd();
		for (int i = 0; i < limit; i += 4) {
			simde__m256d data = simde_mm256_loadu_pd(&x[i]);
			sum0 = simde_mm256_add_pd(sum0, data);
		}
		simde__m128d sum_low = simde_mm256_castpd256_pd128(sum0);
		simde__m128d sum_high = simde_mm256_extractf128_pd(sum0, 1);
		sum_low = simde_mm_add_pd(sum_low, sum_high);
		simde__m128d sum_swapped = simde_mm_shuffle_pd(sum_low, sum_low, 1);
		simde__m128d sum_total = simde_mm_add_pd(sum_low, sum_swapped);
		simde_mm_store_sd(&r, sum_total);
	}
	for (int i = limit; i < n; i++) {
		r += x[i];
	}
	return r;
#endif
