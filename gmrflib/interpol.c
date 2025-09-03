#include <assert.h>
#include <float.h>
#include <math.h>
#include <omp.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

GMRFLib_spline_tp *GMRFLib_spline_create(double *x, double *y, int n)
{
	// return GMRFLib_spline_create_x(x, y, n, GMRFLib_INTPOL_TRANS_NONE, GMRFLib_INTPOL_CACHE_LEVEL12, 0);
	return GMRFLib_spline_create_x(x, y, n, GMRFLib_INTPOL_TRANS_NONE, GMRFLib_INTPOL_CACHE_NONE, 0);
}

GMRFLib_spline_tp *GMRFLib_spline_create_x(double *x, double *y, int n, GMRFLib_intpol_transform_tp trans, GMRFLib_intpol_cache_tp cache,
					   int skip_checks)
{
	/*
	 * Return a spline interpolant for {(x,y)}
	 * cache=0:cache only on level 1, if cache=1: cache on both levels, cache=2: serial cache, cache=3: none.
	 * if skip_checks!=0, then skip check for sorted 'x'
	 */

	if (n < 3) {
		return NULL;
	}

	int nn = n, special = 0;
	double *xx = NULL, *yy = NULL;
	double eps = (GSL_SQRT_DBL_EPSILON * GSL_ROOT4_DBL_EPSILON);
	GMRFLib_spline_tp *s = Calloc(1, GMRFLib_spline_tp);

	Malloc_init(2 * n, 2);
	xx = Malloc_get(n);
	yy = Malloc_get(n);
	Memcpy(xx, x, n * sizeof(double));
	Memcpy(yy, y, n * sizeof(double));

	if (trans == GMRFLib_INTPOL_TRANS_P) {
		special = 1;
		for (int i = 0; i < n; i++) {
			yy[i] = TRUNCATE(yy[i], eps, 1.0 - eps);
			yy[i] = GMRFLib_logit(yy[i]);
		}
	} else if (trans == GMRFLib_INTPOL_TRANS_Pinv) {
		special = 1;
		for (int i = 0; i < n; i++) {
			xx[i] = TRUNCATE(xx[i], eps, 1.0 - eps);
			xx[i] = GMRFLib_logit(xx[i]);
		}
	}

	if (!skip_checks) {
		// normally, 'xx' is sorted, but...
		if (!GMRFLib_is_sorted_dinc(n, xx)) {
			my_sort2_dd(xx, yy, n);
		}
	}
	GMRFLib_unique_additive2(&nn, xx, yy, GSL_SQRT_DBL_EPSILON);

	if (nn < 3) {
		Malloc_free();
		return NULL;
	}

	s->trans = trans;
	s->xmin = xx[0];
	s->xmax = xx[nn - 1];

	switch (cache) {
	case GMRFLib_INTPOL_CACHE_LEVEL12:
		s->cache = GMRFLib_INTPOL_CACHE_LEVEL12;
		s->cache_len = GMRFLib_CACHE_LEN();
		break;
	case GMRFLib_INTPOL_CACHE_LEVEL1:
		s->cache = GMRFLib_INTPOL_CACHE_LEVEL1;
		s->cache_len = GMRFLib_CACHE_LEN_LEVEL1_ONLY();
		break;
	case GMRFLib_INTPOL_CACHE_SIMPLE:
		s->cache = GMRFLib_INTPOL_CACHE_SIMPLE;
		s->cache_len = 1;
		break;
	case GMRFLib_INTPOL_CACHE_NONE:
		s->cache = GMRFLib_INTPOL_CACHE_NONE;
		s->cache_len = 0;
		break;
	default:
		assert(0 == 1);
	}

	if (s->cache_len > 0) {
		// not GMRFLib_INTPOL_CACHE_NONE
		s->accel = Calloc(s->cache_len, gsl_interp_accel *);
		s->accel[0] = gsl_interp_accel_alloc();	       /* rest will be created if needed */
	} else {
		s->accel = NULL;
	}

	if (special) {
		// we know its monotone inbetween the datapoints
		s->spline = gsl_spline_alloc((nn <= 2 ? gsl_interp_linear : gsl_interp_steffen), (unsigned int) nn);
	} else {
		s->spline = gsl_spline_alloc((nn <= 2 ? gsl_interp_linear : gsl_interp_cspline), (unsigned int) nn);
		// s->spline = gsl_spline_alloc((nn <= 2 ? gsl_interp_linear : gsl_interp_akima), (unsigned int) nn);
	}
	gsl_spline_init(s->spline, xx, yy, (unsigned int) nn);
	Malloc_free();

	return s;
}

GMRFLib_spline_tp *GMRFLib_spline_create_from_matrix(GMRFLib_matrix_tp *M)
{
	/*
	 * Return a spline interpolant between {(x,y)} where x is the first column in M and y is the second column of M. 
	 */

	GMRFLib_ASSERT_RETVAL(M->ncol == 2, GMRFLib_EINVARG, NULL);
	GMRFLib_ASSERT_RETVAL(M->nrow > 0, GMRFLib_EINVARG, NULL);
	GMRFLib_ASSERT_RETVAL(M->A, GMRFLib_EINVARG, NULL);

	double *x = M->A;
	double *y = M->A + M->nrow;			       /* column-based storage */

	return GMRFLib_spline_create(x, y, M->nrow);
}

double GMRFLib_spline_eval(double x, GMRFLib_spline_tp *s)
{
	/*
	 * Evaluate a spline 's' in point 'x' 
	 */
	double xx, xx_raw, val;
	static double eps = (GSL_SQRT_DBL_EPSILON * GSL_ROOT4_DBL_EPSILON);

	if (s->trans == GMRFLib_INTPOL_TRANS_Pinv) {
		xx_raw = TRUNCATE(x, eps, 1.0 - eps);
		xx_raw = GMRFLib_logit(xx_raw);
	} else {
		xx_raw = x;
	}
	xx = TRUNCATE(xx_raw, s->xmin, s->xmax);

	int tnum = 0;
	switch (s->cache) {
	case GMRFLib_INTPOL_CACHE_LEVEL12:
		GMRFLib_CACHE_SET_IDX(tnum);
		break;
	case GMRFLib_INTPOL_CACHE_LEVEL1:
		GMRFLib_CACHE_SET_IDX_LEVEL1_ONLY(tnum);
		break;
	case GMRFLib_INTPOL_CACHE_SIMPLE:
	case GMRFLib_INTPOL_CACHE_NONE:
		break;
	default:
		assert(0 == 1);
	}

	gsl_interp_accel *acc = NULL;
	if (s->accel) {
		if (tnum >= 0 && tnum < s->cache_len) {
			if (!(s->accel[tnum])) {
#pragma omp critical (Name_4ebacac2070ee6e249766cf77276653b9f3b684d)
				if (!(s->accel[tnum])) {
					gsl_interp_accel *ac = gsl_interp_accel_alloc();
					s->accel[tnum] = ac;
				}
			}
			acc = s->accel[tnum];
		}
	}

	val = gsl_spline_eval(s->spline, xx, acc);
	if (xx > s->xmin && xx < s->xmax) {
		// we are all fine
	} else {
		double h = (s->xmax - s->xmin) * 1e-4;
		double x_mid = (s->xmax + s->xmin) / 2.0;
		double h_mid = s->xmax - x_mid;
		double vval = 0.0, grad = 0.0, grad_mid = 0.0;
		double val_mid = gsl_spline_eval(s->spline, x_mid, acc);
		int increasing = 1;

		if (ISEQUAL(xx, s->xmax)) {
			vval = gsl_spline_eval(s->spline, xx - h, acc);
			grad = (val - vval) / h;
			if (s->trans == GMRFLib_INTPOL_TRANS_P || s->trans == GMRFLib_INTPOL_TRANS_Pinv) {
				increasing = (val > val_mid);
				grad_mid = (val - val_mid) / h_mid;
			}
		} else if (ISEQUAL(xx, s->xmin)) {
			vval = gsl_spline_eval(s->spline, xx + h, acc);
			grad = (vval - val) / h;
			if (s->trans == GMRFLib_INTPOL_TRANS_P || s->trans == GMRFLib_INTPOL_TRANS_Pinv) {
				increasing = (val_mid > val);
				grad_mid = (val_mid - val) / h_mid;
			}
		} else {
			assert(0 == 1);
		}

		// if grad has the wrong sign, then use the crude one
		if ((increasing && grad < 0.0) || (!increasing && grad > 0.0)) {
			grad = grad_mid;
		}
		val = val + (xx_raw - xx) * grad;
	}

	if (s->trans == GMRFLib_INTPOL_TRANS_P) {
		val = GMRFLib_inv_logit(val);
		val = TRUNCATE(val, eps, 1.0 - eps);
	}

	return val;
}

int GMRFLib_spline_eval_x(int n, double *x, GMRFLib_spline_tp *s, double *values)
{
	/*
	 * Evaluate a spline 's' in point 'x' n times, for increasing 'x'. simple case only.
	 */

	if (x[0] < s->xmin || x[n - 1] > s->xmax || s->trans != GMRFLib_INTPOL_TRANS_NONE) {
		// do this the slow-way
		for (int i = 0; i < n; i++) {
			values[i] = GMRFLib_spline_eval(x[i], s);
		}
		return GMRFLib_SUCCESS;
	}

	int tnum = 0;
	switch (s->cache) {
	case GMRFLib_INTPOL_CACHE_LEVEL12:
		GMRFLib_CACHE_SET_IDX(tnum);
		break;
	case GMRFLib_INTPOL_CACHE_LEVEL1:
		GMRFLib_CACHE_SET_IDX_LEVEL1_ONLY(tnum);
		break;
	case GMRFLib_INTPOL_CACHE_SIMPLE:
	case GMRFLib_INTPOL_CACHE_NONE:
		break;
	default:
		assert(0 == 1);
	}

	gsl_interp_accel *acc = NULL;
	if (s->accel) {
		if (tnum >= 0 && tnum < s->cache_len) {
			if (!(s->accel[tnum])) {
#pragma omp critical (Name_ab9a02f89e7e7b03314b34ac0715d9a6a335e0e2)
				if (!(s->accel[tnum])) {
					gsl_interp_accel *ac = gsl_interp_accel_alloc();
					s->accel[tnum] = ac;
				}
			}
			acc = s->accel[tnum];
		}
	} else {
		// use temporary cache since we evaluate 'n' at the time
		acc = gsl_interp_accel_alloc();
	}

	for (int i = 0; i < n; i++) {
		values[i] = gsl_spline_eval(s->spline, x[i], acc);
	}

	if (!(s->accel) && acc) {
		gsl_interp_accel_free(acc);
	}

	return GMRFLib_SUCCESS;
}

double GMRFLib_spline_eval_deriv(double x, GMRFLib_spline_tp *s)
{
	/*
	 * Evaluate the derivative of the spline 's' in point 'x' 
	 */

	// not yet implemented, I'm not sure I need this for P and Pinv
	assert(s->trans == GMRFLib_INTPOL_TRANS_NONE || s->trans == GMRFLib_INTPOL_TRANS_Pinv);
	double val = 0.0;

	int tnum = 0;
	switch (s->cache) {
	case GMRFLib_INTPOL_CACHE_LEVEL12:
		GMRFLib_CACHE_SET_IDX(tnum);
		break;
	case GMRFLib_INTPOL_CACHE_LEVEL1:
		GMRFLib_CACHE_SET_IDX_LEVEL1_ONLY(tnum);
		break;
	case GMRFLib_INTPOL_CACHE_SIMPLE:
	case GMRFLib_INTPOL_CACHE_NONE:
		break;
	default:
		assert(0 == 1);
	}

	gsl_interp_accel *acc = NULL;
	if (s->accel) {
		if (tnum >= 0 && tnum < s->cache_len) {
			if (!(s->accel[tnum])) {
#pragma omp critical (Name_bcc8a7f7a416bde91e4459c229fc294985c3674c)
				if (!(s->accel[tnum])) {
					gsl_interp_accel *ac = gsl_interp_accel_alloc();
					s->accel[tnum] = ac;
				}
			}
			acc = s->accel[tnum];
		}
	}

	if (s->trans == GMRFLib_INTPOL_TRANS_NONE) {
		val = gsl_spline_eval_deriv(s->spline, TRUNCATE(x, s->xmin, s->xmax), acc);
	} else if (s->trans == GMRFLib_INTPOL_TRANS_Pinv) {
		double xx = GMRFLib_logit(x);
		val = gsl_spline_eval_deriv(s->spline, TRUNCATE(xx, s->xmin, s->xmax), acc);
		double em = exp(-xx);
		val *= (em + 2.0 + 1.0 / em);
	} else {
		assert(0 == 1);
	}
	return val;
}

double GMRFLib_spline_eval_deriv2(double x, GMRFLib_spline_tp *s)
{
	/*
	 * Evaluate the 2.derivative of the spline 's' in point 'x' 
	 */

	// not yet implemented, I'm not sure I need this for P and Pinv
	assert(s->trans == GMRFLib_INTPOL_TRANS_NONE);
	double val;
	if (x < s->xmin || x > s->xmax) {
		val = NAN;
	} else {
		int tnum = 0;
		switch (s->cache) {
		case GMRFLib_INTPOL_CACHE_LEVEL12:
			GMRFLib_CACHE_SET_IDX(tnum);
			break;
		case GMRFLib_INTPOL_CACHE_LEVEL1:
			GMRFLib_CACHE_SET_IDX_LEVEL1_ONLY(tnum);
			break;
		case GMRFLib_INTPOL_CACHE_SIMPLE:
		case GMRFLib_INTPOL_CACHE_NONE:
			break;
		default:
			assert(0 == 1);
		}

		gsl_interp_accel *acc = NULL;
		if (s->accel) {
			if (tnum >= 0 && tnum < s->cache_len) {
				if (!(s->accel[tnum])) {
#pragma omp critical (Name_7db308fb16056e07320f9aa74e5445c74a6f298f)
					if (!(s->accel[tnum])) {
						gsl_interp_accel *ac = gsl_interp_accel_alloc();
						s->accel[tnum] = ac;
					}
				}
				acc = s->accel[tnum];
			}
		}
		val = gsl_spline_eval_deriv2(s->spline, x, acc);
	}

	return val;
}

double GMRFLib_spline_eval_deriv_x(double x, GMRFLib_spline_tp *s)
{
	// this expert version do not check for 's->trans'
	double val;
	if (x < s->xmin || x > s->xmax) {
		val = NAN;
	} else {
		int tnum = 0;
		switch (s->cache) {
		case GMRFLib_INTPOL_CACHE_LEVEL12:
			GMRFLib_CACHE_SET_IDX(tnum);
			break;
		case GMRFLib_INTPOL_CACHE_LEVEL1:
			GMRFLib_CACHE_SET_IDX_LEVEL1_ONLY(tnum);
			break;
		case GMRFLib_INTPOL_CACHE_SIMPLE:
		case GMRFLib_INTPOL_CACHE_NONE:
			break;
		default:
			assert(0 == 1);
		}

		gsl_interp_accel *acc = NULL;
		if (s->accel) {
			if (tnum >= 0 && tnum < s->cache_len) {
				if (!(s->accel[tnum])) {
#pragma omp critical (Name_8c1f6a9b1676b904f0235f5d9f3817343bf0b5d3)
					if (!(s->accel[tnum])) {
						gsl_interp_accel *ac = gsl_interp_accel_alloc();
						s->accel[tnum] = ac;
					}
				}
				acc = s->accel[tnum];
			}
		}
		val = gsl_spline_eval_deriv(s->spline, x, acc);
	}
	return val;
}

double GMRFLib_spline_eval_deriv2_x(double x, GMRFLib_spline_tp *s)
{
	// this expert version do not check for 's->trans'
	double val;
	if (x < s->xmin || x > s->xmax) {
		val = NAN;
	} else {
		int tnum = 0;
		switch (s->cache) {
		case GMRFLib_INTPOL_CACHE_LEVEL12:
			GMRFLib_CACHE_SET_IDX(tnum);
			break;
		case GMRFLib_INTPOL_CACHE_LEVEL1:
			GMRFLib_CACHE_SET_IDX_LEVEL1_ONLY(tnum);
			break;
		case GMRFLib_INTPOL_CACHE_SIMPLE:
		case GMRFLib_INTPOL_CACHE_NONE:
			break;
		default:
			assert(0 == 1);
		}

		gsl_interp_accel *acc = NULL;
		if (s->accel) {
			if (tnum >= 0 && tnum < s->cache_len) {
				if (!(s->accel[tnum])) {
#pragma omp critical (Name_6c1aed3d698e547929f98757e8a8e32e2adc4b68)
					if (!(s->accel[tnum])) {
						gsl_interp_accel *ac = gsl_interp_accel_alloc();
						s->accel[tnum] = ac;
					}
				}
				acc = s->accel[tnum];
			}
		}
		val = gsl_spline_eval_deriv2(s->spline, x, acc);
	}
	return val;
}

int GMRFLib_spline_free(GMRFLib_spline_tp *s)
{
	if (s) {
		gsl_spline_free(s->spline);
		if (s->accel) {
			for (int i = 0; i < s->cache_len; i++) {
				if (s->accel[i]) {
					gsl_interp_accel_free(s->accel[i]);
				}
			}
			Free(s->accel);
		}
		Free(s);
	}

	return GMRFLib_SUCCESS;
}
