
/* dot.c
 * 
 * Copyright (C) 2022-2024 Havard Rue
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
#include "GMRFLib/dot.h"

double GMRFLib_dot_product(GMRFLib_idxval_tp *__restrict ELM_, double *__restrict ARR_)
{
	if (ELM_->dot_product_func) {
#if !defined(INLA_WITH_MKL)
		if (GMRFLib_dot_product_gain >= 0.0) {
#pragma omp atomic
			GMRFLib_dot_product_gain += ELM_->cpu_gain;
		}
#endif
		return (ELM_->dot_product_func((GMRFLib_idxval_tp * __restrict) ELM_, (double *__restrict) ARR_));
	} else {
		return GMRFLib_dot_product_serial_mkl(ELM_, ARR_);
	}
}

double GMRFLib_dot_product_group(GMRFLib_idxval_tp *__restrict ELM_, double *__restrict ARR_)
{
	double value_ = 0.0;
	for (int g_ = 0; g_ < ELM_->g_n; g_++) {
		int len_ = ELM_->g_len[g_];
		int *__restrict ii_ = ELM_->g_idx[g_];
		double *__restrict vv_ = ELM_->g_val[g_];

		if (len_ > 0) {
			double *__restrict aa_ = &(ARR_[0]);
			value_ += (ELM_->g_1[g_] ? GMRFLib_dsum_idx(len_, aa_, ii_) : GMRFLib_ddot_idx(len_, vv_, aa_, ii_));
		} else if (len_ < 0) {
			int llen_ = -len_;
			double *__restrict aa_ = &(ARR_[ii_[0]]);
			value_ += (ELM_->g_1[g_] ? GMRFLib_dsum(llen_, aa_) : GMRFLib_ddot(llen_, vv_, aa_));
		}
	}
	return value_;
}

double GMRFLib_dot_product_group_mkl(GMRFLib_idxval_tp *__restrict ELM_, double *__restrict ARR_)
{
	double value_ = 0.0;
	for (int g_ = 0; g_ < ELM_->g_n; g_++) {
		int len_ = ELM_->g_len[g_];
		int *__restrict ii_ = ELM_->g_idx[g_];
		double *__restrict vv_ = ELM_->g_val[g_];

		if (len_ > 0) {
			double *__restrict aa_ = &(ARR_[0]);
			value_ += (ELM_->g_1[g_] ? GMRFLib_dsum_idx(len_, aa_, ii_) : GMRFLib_ddot_idx_mkl(len_, vv_, aa_, ii_));
		} else if (len_ < 0) {
			int llen_ = -len_;
			double *__restrict aa_ = &(ARR_[ii_[0]]);
			value_ += (ELM_->g_1[g_] ? GMRFLib_dsum(llen_, aa_) : GMRFLib_ddot(llen_, vv_, aa_));
		}
	}
	return value_;
}

double GMRFLib_dot_product_group_mkl_alt(GMRFLib_idxval_tp *__restrict ELM_, double *__restrict ARR_)
{
	double value_ = 0.0;
	for (int g_ = 0; g_ < ELM_->g_n; g_++) {
		int len_ = ELM_->g_len[g_];
		int *__restrict ii_ = ELM_->g_idx[g_];
		double *__restrict vv_ = ELM_->g_val[g_];

		if (len_ > 0) {
			double *__restrict aa_ = &(ARR_[0]);
			value_ += (ELM_->g_1[g_] ? GMRFLib_dsum_idx(len_, aa_, ii_) : GMRFLib_ddot_idx_mkl_alt(len_, vv_, aa_, ii_));
		} else if (len_ < 0) {
			int llen_ = -len_;
			double *__restrict aa_ = &(ARR_[ii_[0]]);
			value_ += (ELM_->g_1[g_] ? GMRFLib_dsum(llen_, aa_) : GMRFLib_ddot(llen_, vv_, aa_));
		}
	}
	return value_;
}

forceinline double GMRFLib_dot_product_serial(GMRFLib_idxval_tp *__restrict ELM_, double *__restrict ARR_)
{
	double *__restrict vv_ = ELM_->val;
	double *__restrict aa_ = ARR_;
	int *__restrict idx_ = ELM_->idx;

	return GMRFLib_ddot_idx(ELM_->n, vv_, aa_, idx_);
}

forceinline double GMRFLib_dot_product_serial_mkl(GMRFLib_idxval_tp *__restrict ELM_, double *__restrict ARR_)
{
	double *__restrict vv_ = ELM_->val;
	double *__restrict aa_ = ARR_;
	int *__restrict idx_ = ELM_->idx;

	return GMRFLib_ddot_idx_mkl(ELM_->n, vv_, aa_, idx_);
}

forceinline double GMRFLib_dot_product_serial_mkl_alt(GMRFLib_idxval_tp *__restrict ELM_, double *__restrict ARR_)
{
	double *__restrict vv_ = ELM_->val;
	double *__restrict aa_ = ARR_;
	int *__restrict idx_ = ELM_->idx;

	return GMRFLib_ddot_idx_mkl_alt(ELM_->n, vv_, aa_, idx_);
}

forceinline double GMRFLib_ddot(int n, double *x, double *y)
{
	int one = 1;
	return ddot_(&n, x, &one, y, &one);
}

int GMRFLib_isum(int n, int *ix)
{
	int s = 0;
#pragma omp simd reduction(+: s)
	for (int i = 0; i < n; i++) {
		s += ix[i];
	}
	return s;
}

double GMRFLib_dsum(int n, double *x)
{
	double s = 0.0;
#pragma omp simd reduction(+: s)
	for (int i = 0; i < n; i++) {
		s += x[i];
	}
	return s;
}

double GMRFLib_dsum_idx(int n, double *__restrict a, int *__restrict idx)
{
	const int roll = 8L;
	double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
	div_t d = div(n, roll);
	int m = d.quot * roll;

#pragma GCC ivdep
	for (int i = 0; i < m; i += roll) {
		int *iidx = idx + i;

		s0 += a[iidx[0]];
		s1 += a[iidx[1]];
		s2 += a[iidx[2]];
		s3 += a[iidx[3]];

		s0 += a[iidx[4]];
		s1 += a[iidx[5]];
		s2 += a[iidx[6]];
		s3 += a[iidx[7]];
	}

#pragma GCC ivdep
	for (int i = m; i < n; i++) {
		s0 += a[idx[i]];
	}

	return s0 + s1 + s2 + s3;
}

double GMRFLib_ddot_idx(int n, double *__restrict v, double *__restrict a, int *__restrict idx)
{
	const int roll = 8L;
	double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
	div_t d = div(n, roll);
	int m = d.quot * roll;

#pragma GCC ivdep
	for (int i = 0; i < m; i += roll) {
		double *vv = v + i;
		int *iidx = idx + i;

		s0 += vv[0] * a[iidx[0]];
		s1 += vv[1] * a[iidx[1]];
		s2 += vv[2] * a[iidx[2]];
		s3 += vv[3] * a[iidx[3]];

		s0 += vv[4] * a[iidx[4]];
		s1 += vv[5] * a[iidx[5]];
		s2 += vv[6] * a[iidx[6]];
		s3 += vv[7] * a[iidx[7]];
	}

#pragma GCC ivdep
	for (int i = m; i < n; i++) {
		s0 += v[i] * a[idx[i]];
	}

	return s0 + s1 + s2 + s3;
}

#if defined(INLA_WITH_MKL)

double GMRFLib_ddot_idx_mkl_alt(int n, double *__restrict v, double *__restrict a, int *__restrict idx)
{
	int iarr[4] = { 1, 0, n, idx[n - 1] + 1 };
	double darr[3] = { 1.0, 0.0, 0.0 };

	// we need to define this with length 6. the fifth argument is not used, so we use it for the
	// argument 'trans', the first argument in the call, trans='N'
	const char matdescra[6] = { 'G', '.', '.', 'C', 'N', '.' };

	mkl_dcsrmv(matdescra + 4, iarr, iarr + 3, darr, matdescra, v, idx, iarr + 1, iarr + 2, a, darr + 1, darr + 2);
	return darr[2];
}

forceinline double GMRFLib_ddot_idx_mkl(int n, double *__restrict v, double *__restrict a, int *__restrict idx)
{
	return cblas_ddoti(n, v, idx, a);
}

#else							       /* defined(INLA_WITH_MKL) */

forceinline double GMRFLib_ddot_idx_mkl_alt(int n, double *__restrict v, double *__restrict a, int *__restrict idx)
{
	return GMRFLib_ddot_idx(n, v, a, idx);
}

forceinline double GMRFLib_ddot_idx_mkl(int n, double *__restrict v, double *__restrict a, int *__restrict idx)
{
	return GMRFLib_ddot_idx(n, v, a, idx);
}

#endif							       /* if defined(INLA_WITH_MKL) */
