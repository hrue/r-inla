	int i = 0;
	for (; i <= n - 8; i += 8) {
		if (iy[i+1] == iy[i] + 1 && iy[i+3] == iy[i+2] + 1 && 
		    iy[i+5] == iy[i+4] + 1 && iy[i+7] == iy[i+6] + 1) {
			__m128d vals0 = _mm_loadu_pd(&a[i]);
			__m128d vals1 = _mm_loadu_pd(&a[i+2]);
			__m128d vals2 = _mm_loadu_pd(&a[i+4]);
			__m128d vals3 = _mm_loadu_pd(&a[i+6]);
			_mm_storeu_pd(&y[iy[i]], vals0);
			_mm_storeu_pd(&y[iy[i+2]], vals1);
			_mm_storeu_pd(&y[iy[i+4]], vals2);
			_mm_storeu_pd(&y[iy[i+6]], vals3);
		} else {
			__m128d vals0 = _mm_loadu_pd(&a[i]);
			__m128d vals1 = _mm_loadu_pd(&a[i+2]);
			__m128d vals2 = _mm_loadu_pd(&a[i+4]);
			__m128d vals3 = _mm_loadu_pd(&a[i+6]);
			double val0 = _mm_cvtsd_f64(vals0);
			double val1 = _mm_cvtsd_f64(_mm_shuffle_pd(vals0, vals0, _MM_SHUFFLE2(0,1)));
			double val2 = _mm_cvtsd_f64(vals1);
			double val3 = _mm_cvtsd_f64(_mm_shuffle_pd(vals1, vals1, _MM_SHUFFLE2(0,1)));
			double val4 = _mm_cvtsd_f64(vals2);
			double val5 = _mm_cvtsd_f64(_mm_shuffle_pd(vals2, vals2, _MM_SHUFFLE2(0,1)));
			double val6 = _mm_cvtsd_f64(vals3);
			double val7 = _mm_cvtsd_f64(_mm_shuffle_pd(vals3, vals3, _MM_SHUFFLE2(0,1)));
			y[iy[i]] = val0;
			y[iy[i+1]] = val1;
			y[iy[i+2]] = val2;
			y[iy[i+3]] = val3;
			y[iy[i+4]] = val4;
			y[iy[i+5]] = val5;
			y[iy[i+6]] = val6;
			y[iy[i+7]] = val7;
		}
	}
	for (; i <= n - 2; i += 2) {
		if (iy[i+1] == iy[i] + 1) {
			__m128d vals = _mm_loadu_pd(&a[i]);
			_mm_storeu_pd(&y[iy[i]], vals);
		} else {
			__m128d vals = _mm_loadu_pd(&a[i]);
			double val0 = _mm_cvtsd_f64(vals);
			double val1 = _mm_cvtsd_f64(_mm_shuffle_pd(vals, vals, _MM_SHUFFLE2(0,1)));
			y[iy[i]] = val0;
			y[iy[i+1]] = val1;
		}
	}
	for (; i < n; i++) {
		y[iy[i]] = a[i];
	}
