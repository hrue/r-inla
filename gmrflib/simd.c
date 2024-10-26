
/* simd.c
 * 
 * Copyright (C) 2024-2024 Havard Rue
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * The author's contact information:
 *
 *        Haavard Rue
 *        CEMSE Division
 *        King Abdullah University of Science and Technology
 *        Thuwal 23955-6900, Saudi Arabia
 *        Email: haavard.rue@kaust.edu.sa
 *        Office: +966 (0)12 808 0640
 *
 */

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

void GMRFLib_exp(int n, double *x, double *y)
{
#if defined(INLA_WITH_SIMD)
	// y = exp(x)
	if (SIMD_ALIGNED(x) && SIMD_ALIGNED(y)) {
		const int len = 4L;
		div_t d = div(n, len);
		int m = d.quot * len;

		for (int i = 0; i < m; i += len) {
			double *xp = x + i;
			double *yp = y + i;
			__m256d vx = _mm256_loadu_pd(xp);
			__m256d vy = Sleef_expd4_u10(vx);
			_mm256_storeu_pd(yp, vy);
		}
		int k = n - m;
		if (k) {
			if (k <= 2) {
				double xx[2] = { 0, 0 };
				Memcpy(xx, x + m, k * sizeof(double));

				__m128d vx = _mm_loadu_pd(xx);
				__m128d vy = Sleef_expd2_u10(vx);

				if (k == 2) {
					_mm_storeu_pd(y + m, vy);
				} else {
					double yy[2];
					_mm_storeu_pd(yy, vy);
					Memcpy(y + m, yy, k * sizeof(double));
				}
			} else {
				double xx[] = { 0, 0, 0, 0 };
				Memcpy(xx, x + m, k * sizeof(double));

				__m256d vx = _mm256_loadu_pd(xx);
				__m256d vy = Sleef_expd4_u10(vx);

				if (k == 4L) {
					_mm256_storeu_pd(y + m, vy);
				} else {
					double yy[len];
					_mm256_storeu_pd(yy, vy);
					Memcpy(y + m, yy, k * sizeof(double));
				}
			}
		}
	} else {
		if (SIMD_ALIGNED(x + 1) && SIMD_ALIGNED(y + 1)) {
			y[0] = exp(x[0]);
			return (GMRFLib_exp(n - 1, x + 1, y + 1));
		} else {
			FIXME("miss alignment in _exp");
		}

#if defined(INLA_WITH_MKL)
		vdExp(n, x, y);
#elif defined(INLA_WITH_FRAMEWORK_ACCELERATE)
		vvexp(y, x, &n);
#else
		_Pragma("omp simd")
		    for (int i = 0; i < n; i++) {
			y[i] = exp(x[i]);
		}
#endif
	}
#else
#if defined(INLA_WITH_MKL)
	vdExp(n, x, y);
#elif defined(INLA_WITH_FRAMEWORK_ACCELERATE)
	vvexp(y, x, &n);
#else
	_Pragma("omp simd")
	    for (int i = 0; i < n; i++) {
		y[i] = exp(x[i]);
	}
#endif
#endif
}

void GMRFLib_exp_inc(int n, double *x, int inc, double *y)
{
	// y = exp(x) with inc
#if defined(INLA_WITH_MKL)
	vdExpI(n, x, inc, y, inc);
#else
	_Pragma("omp simd")
	    for (int i = 0; i < n * inc; i += inc) {
		y[i] = exp(x[i]);
	}
#endif
}

void GMRFLib_log(int n, double *x, double *y)
{
	// y = log(x)
#if defined(INLA_WITH_SIMD)
	if (SIMD_ALIGNED(x) && SIMD_ALIGNED(y)) {
		const int len = 4L;
		div_t d = div(n, len);
		int m = d.quot * len;

		for (int i = 0; i < m; i += len) {
			double *xp = x + i;
			double *yp = y + i;
			__m256d vx = _mm256_loadu_pd(xp);
			__m256d vy = Sleef_logd4_u10(vx);
			_mm256_storeu_pd(yp, vy);
		}
		int k = n - m;
		if (k) {
			if (k <= 2) {
				double xx[2] = { 1, 1 };
				Memcpy(xx, x + m, k * sizeof(double));
				__m128d vx = _mm_loadu_pd(xx);
				__m128d vy = Sleef_logd2_u10(vx);

				if (k == 2L) {
					_mm_storeu_pd(y + m, vy);
				} else {
					double yy[2];
					_mm_storeu_pd(yy, vy);
					Memcpy(y + m, yy, k * sizeof(double));
				}
			} else {
				double xx[] = { 1, 1, 1, 1 };
				Memcpy(xx, x + m, k * sizeof(double));
				__m256d vx = _mm256_loadu_pd(xx);
				__m256d vy = Sleef_logd4_u10(vx);

				if (k == len) {
					_mm256_storeu_pd(y + m, vy);
				} else {
					double yy[len];
					_mm256_storeu_pd(yy, vy);
					Memcpy(y + m, yy, k * sizeof(double));
				}
			}
		}
	} else {
		if (SIMD_ALIGNED(x + 1) && SIMD_ALIGNED(y + 1)) {
			y[0] = log(x[0]);
			return (GMRFLib_log(n - 1, x + 1, y + 1));
		} else {
			FIXME("miss alignment in _log");
		}

#if defined(INLA_WITH_MKL)
		vdLn(n, x, y);
#elif defined(INLA_WITH_FRAMEWORK_ACCELERATE)
		vvlog(y, x, &n);
#else
		_Pragma("omp simd")
		    for (int i = 0; i < n; i++) {
			y[i] = log(x[i]);
		}
#endif
	}
#else
#if defined(INLA_WITH_MKL)
	vdLn(n, x, y);
#elif defined(INLA_WITH_FRAMEWORK_ACCELERATE)
	vvlog(y, x, &n);
#else
	_Pragma("omp simd")
	    for (int i = 0; i < n; i++) {
		y[i] = log(x[i]);
	}
#endif
#endif
}

void GMRFLib_log1p(int n, double *x, double *y)
{
	// y = log1p(x)
#if defined(INLA_WITH_SIMD)
	if (SIMD_ALIGNED(x) && SIMD_ALIGNED(y)) {
		const int len = 4L;
		div_t d = div(n, len);
		int m = d.quot * len;

		for (int i = 0; i < m; i += len) {
			double *xp = x + i;
			double *yp = y + i;
			__m256d vx, vy;
			vx = _mm256_loadu_pd(xp);
			vy = Sleef_log1pd4_u10(vx);
			_mm256_storeu_pd(yp, vy);
		}
		int k = n - m;
		if (k) {
			if (k <= 2) {
				double xx[2] = { 0, 0 };
				Memcpy(xx, x + m, k * sizeof(double));

				__m128d vx = _mm_loadu_pd(xx);
				__m128d vy = Sleef_log1pd2_u10(vx);

				if (k == 2L) {
					_mm_storeu_pd(y + m, vy);
				} else {
					double yy[2];
					_mm_storeu_pd(yy, vy);
					Memcpy(y + m, yy, k * sizeof(double));
				}
			} else {
				double xx[] = { 0, 0, 0, 0 };
				Memcpy(xx, x + m, k * sizeof(double));

				__m256d vx = _mm256_loadu_pd(xx);
				__m256d vy = Sleef_log1pd4_u10(vx);

				if (k == len) {
					_mm256_storeu_pd(y + m, vy);
				} else {
					double yy[len];
					_mm256_storeu_pd(yy, vy);
					Memcpy(y + m, yy, k * sizeof(double));
				}
			}
		}
	} else {
		if (SIMD_ALIGNED(x + 1) && SIMD_ALIGNED(y + 1)) {
			y[0] = log1p(x[0]);
			return (GMRFLib_log1p(n - 1, x + 1, y + 1));
		} else {
			FIXME("miss alignment in _log1p");
		}

#if defined(INLA_WITH_MKL)
		vdLog1p(n, x, y);
#elif defined(INLA_WITH_FRAMEWORK_ACCELERATE)
		vvlog1p(y, x, &n);
#else
		_Pragma("omp simd")
		    for (int i = 0; i < n; i++) {
			y[i] = log1p(x[i]);
		}
#endif
	}
#else
#if defined(INLA_WITH_MKL)
	vdLog1p(n, x, y);
#elif defined(INLA_WITH_FRAMEWORK_ACCELERATE)
	vvlog1p(y, x, &n);
#else
	_Pragma("omp simd")
	    for (int i = 0; i < n; i++) {
		y[i] = log1p(x[i]);
	}
#endif
#endif
}

void GMRFLib_sqr(int n, double *x, double *y)
{
	// y = x * x
#if defined(INLA_WITH_MKL)
	if (n < 96) {
		_Pragma("omp simd")
		    for (int i = 0; i < n; i++) {
			y[i] = SQR(x[i]);
		}
	} else {
		vdSqr(n, x, y);
	}
#else
	_Pragma("omp simd")
	    for (int i = 0; i < n; i++) {
		y[i] = SQR(x[i]);
	}
#endif
}

void GMRFLib_sqrt(int n, double *x, double *y)
{
	// y = sqrt(x)
#if defined(INLA_WITH_MKL)
	vdSqrt(n, x, y);
#elif defined(INLA_WITH_FRAMEWORK_ACCELERATE)
	vvsqrt(y, x, &n);
#else
	_Pragma("omp simd")
	    for (int i = 0; i < n; i++) {
		y[i] = sqrt(x[i]);
	}
#endif
}

void GMRFLib_add(int n, double *x, double *y, double *z)
{
	// z = x + y
#if defined(INLA_WITH_MKL)
	vdAdd(n, x, y, z);
#else
	_Pragma("omp simd")
	    for (int i = 0; i < n; i++) {
		z[i] = x[i] + y[i];
	}
#endif
}

void GMRFLib_mul(int n, double *x, double *y, double *z)
{
	// z = x * y
#if defined(INLA_WITH_MKL)
	vdMul(n, x, y, z);
#else
	_Pragma("omp simd")
	    for (int i = 0; i < n; i++) {
		z[i] = x[i] * y[i];
	}
#endif
}

void GMRFLib_daddto(int n, double *x, double *y)
{
	// y = y + x
	int inc = 1;
	double one = 1.0;
	daxpy_(&n, &one, x, &inc, y, &inc);
}

void GMRFLib_cdaddto(int n, double *x, double cx, double *y)
{
	// y = x + const.x
#pragma omp simd
	for (int i = 0; i < n; i++) {
		y[i] = x[i] + cx;
	}
}
