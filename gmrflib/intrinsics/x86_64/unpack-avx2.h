	int i = 0;
	for (; i <= n - 8; i += 8) {
		if (iy[i+1] == iy[i] + 1 && iy[i+2] == iy[i] + 2 && iy[i+3] == iy[i] + 3) {
			__m256d vals0 = _mm256_loadu_pd(&a[i]);
			_mm256_storeu_pd(&y[iy[i]], vals0);
		} else {
			__m256d vals0 = _mm256_loadu_pd(&a[i]);
			__m128d low0 = _mm256_castpd256_pd128(vals0);
			__m128d high0 = _mm256_extractf128_pd(vals0, 1);
			double val0 = _mm_cvtsd_f64(low0);
			double val1 = _mm_cvtsd_f64(_mm_shuffle_pd(low0, low0, 1));
			double val2 = _mm_cvtsd_f64(high0);
			double val3 = _mm_cvtsd_f64(_mm_shuffle_pd(high0, high0, 1));
			y[iy[i]] = val0;
			y[iy[i+1]] = val1;
			y[iy[i+2]] = val2;
			y[iy[i+3]] = val3;
		}
		if (iy[i+5] == iy[i+4] + 1 && iy[i+6] == iy[i+4] + 2 && iy[i+7] == iy[i+4] + 3) {
			__m256d vals1 = _mm256_loadu_pd(&a[i+4]);
			_mm256_storeu_pd(&y[iy[i+4]], vals1);
		} else {
			__m256d vals1 = _mm256_loadu_pd(&a[i+4]);
			__m128d low1 = _mm256_castpd256_pd128(vals1);
			__m128d high1 = _mm256_extractf128_pd(vals1, 1);
			double val4 = _mm_cvtsd_f64(low1);
			double val5 = _mm_cvtsd_f64(_mm_shuffle_pd(low1, low1, 1));
			double val6 = _mm_cvtsd_f64(high1);
			double val7 = _mm_cvtsd_f64(_mm_shuffle_pd(high1, high1, 1));
			y[iy[i+4]] = val4;
			y[iy[i+5]] = val5;
			y[iy[i+6]] = val6;
			y[iy[i+7]] = val7;
		}
	}
	for (; i <= n - 4; i += 4) {
		if (iy[i+1] == iy[i] + 1 && iy[i+2] == iy[i] + 2 && iy[i+3] == iy[i] + 3) {
			__m256d vals = _mm256_loadu_pd(&a[i]);
			_mm256_storeu_pd(&y[iy[i]], vals);
		} else {
			__m256d vals = _mm256_loadu_pd(&a[i]);
			__m128d low = _mm256_castpd256_pd128(vals);
			__m128d high = _mm256_extractf128_pd(vals, 1);
			double val0 = _mm_cvtsd_f64(low);
			double val1 = _mm_cvtsd_f64(_mm_shuffle_pd(low, low, 1));
			double val2 = _mm_cvtsd_f64(high);
			double val3 = _mm_cvtsd_f64(_mm_shuffle_pd(high, high, 1));
			y[iy[i]] = val0;
			y[iy[i+1]] = val1;
			y[iy[i+2]] = val2;
			y[iy[i+3]] = val3;
		}
	}
	for (; i < n; i++) {
		y[iy[i]] = a[i];
	}
