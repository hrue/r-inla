
/* ghq.c
 * 
 * Copyright (C) 2006-2024 Havard Rue
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

/*!
  \file ghq.c
  \brief Utilities for Gauss-Hermite quadrature.

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"
#include "GMRFLib/hashP.h"

int GMRFLib_ghq__intern(double *x, double *w, int n)
{
	int i, its, j, m;
	double p1, p2, p3, pp, z = 0, z1;
	double eps = 3.0e-14;				       /* shown to work fine for n < 199 */

	/*
	 * this function breaks down for n >= 199 
	 */
	GMRFLib_ASSERT(n < 199, GMRFLib_EINVARG);

	x--;						       /* indices from 1...n */
	if (w) {
		w--;					       /* in this routine */
	}

	m = (n + 1) / 2;
	for (i = 1; i <= m; i++) {
		switch (i) {
		case 1:
		{
			z = sqrt(2.0 * n + 1.0) - 1.85575 * pow(2.0 * n + 1.0, -0.16667);
		}
			break;
		case 2:
		{
			z -= 1.14 * pow((double) n, 0.426) / z;
		}
			break;
		case 3:
		{
			z = 1.86 * z - 0.86 * x[1];
		}
			break;
		case 4:
		{
			z = 1.91 * z - 0.91 * x[2];
		}
			break;
		default:
			z = 2.0 * z - x[i - 2];
		}
		its = 0;
		while (1) {
			its++;

			p1 = 0.7511255444649425;
			p2 = 0.0;
			for (j = 1; j <= n; j++) {
				p3 = p2;
				p2 = p1;
				p1 = z * sqrt(2.0 / j) * p2 - sqrt((j - 1.0) / j) * p3;
			}
			pp = sqrt(2.0 * n) * p2;
			z1 = z;
			z = z1 - p1 / pp;
			if (ABS(z - z1) <= eps)
				break;			       /* converged */
			GMRFLib_ASSERT(its <= 500, GMRFLib_ESNH);
		}

		x[i] = z;
		x[n + 1 - i] = -z;

		if (w) {
			w[i] = 2.0 / (pp * pp);
			w[n + 1 - i] = w[i];
		}
	}
	return 0;
}

int GMRFLib_ghq_abscissas(double **xp, int n)
{
	return GMRFLib_ghq(xp, NULL, n);
}

int GMRFLib_ghq_weights(double **wp, int n)
{
	return GMRFLib_ghq(NULL, wp, n);
}

int GMRFLib_ghq_ms(double **xp, double **wp, int n, double mean, double stdev)
{
	// the same for a given mean and stdev. Allocated new memory for xp and wp
	int i;
	double *xxp = NULL, *wwp = NULL;
	GMRFLib_ghq(&xxp, &wwp, n);
	assert(xxp);
	assert(wwp);

	*xp = Calloc(n, double);
	*wp = Calloc(n, double);
	assert(*xp);
	assert(*wp);
	for (i = 0; i < n; i++) {
		(*xp)[i] = xxp[i] * stdev + mean;
		(*wp)[i] = wwp[i];
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_ghq(double **xp, double **wp, int n)
{
	/*
	 * return a pointer to the abscissas and the weights for the Gauss-Hermite quadrature for given `n' with kernel
	 * exp(-x^2/2)/sqrt(2*pi)
	 * 
	 * both abscissas and the weights are optional
	 * 
	 * an internal storage is used to store previously computed abscissas and weights 
	 */

	static map_ivp **abscissas = NULL;		       /* keep previous computed elements here */
	static map_ivp **weights = NULL;		       /* keep previous computed elements here */

	if (!abscissas) {
#pragma omp critical (Name_57dd787c76d8e98b908fbe47a9af2b183bc7a84a)
		{
			if (!abscissas) {
				weights = Calloc(GMRFLib_CACHE_LEN(), map_ivp *);
				abscissas = Calloc(GMRFLib_CACHE_LEN(), map_ivp *);
			}
		}
	}
	int idx = 0;
	GMRFLib_CACHE_SET_ID(idx);

	if (!abscissas[idx]) {
		abscissas[idx] = Calloc(1, map_ivp);
		weights[idx] = Calloc(1, map_ivp);
		map_ivp_init(abscissas[idx]);
		map_ivp_init(weights[idx]);
	}

	int i;
	double *x, *w;
	void *ptr, *pptr;

	GMRFLib_ASSERT(n > 0, GMRFLib_EINVARG);
	if ((ptr = map_ivp_ptr(abscissas[idx], n))) {
		/*
		 * use previously computed. note that map_ivp_ptr returns a ptr to the stored ptr. 
		 */
		if (xp) {
			*xp = *((double **) ptr);
		}

		/*
		 * the weights should now be stored as well 
		 */
		if (wp) {
			pptr = map_ivp_ptr(weights[idx], n);
			GMRFLib_ASSERT(pptr, GMRFLib_ESNH);
			*wp = *((double **) pptr);
		}
	} else {
		/*
		 * compute new ones 
		 */

		x = Calloc(n, double);
		w = Calloc(n, double);

		GMRFLib_EWRAP0(GMRFLib_ghq__intern(x, w, n));

		/*
		 * the Gauss-Hermite is with kernel exp(-x^2), transform to kernel exp(-x^2/2)/sqrt(2*pi)
		 */
		for (i = 0; i < n; i++) {
			x[i] *= M_SQRT2;
			w[i] *= M_SQRT2 / sqrt(2.0 * M_PI);
		}

		/*
		 * reverse the order so its small to large, and sort the weights along 
		 */
		// GMRFLib_qsort2((void *) x, (size_t) n, sizeof(double), w, sizeof(double), NULL, 0, GMRFLib_dcmp);
		my_sort2_dd(x, w, n);
		map_ivp_set(abscissas[idx], n, (void *) x);
		map_ivp_set(weights[idx], n, (void *) w);

		if (xp) {
			*xp = x;			       /* return a ptr only */
		}
		if (wp) {
			*wp = w;			       /* return a ptr only */
		}
	}
	return GMRFLib_SUCCESS;
}

GMRFLib_snq_tp *GMRFLib_snq(int n, double skew3)
{
	// RATIO is the skew-normal density divided by the normal, each with mean zero and unit variance
	// skew3 is the skewness^(1/3)
#define RATIO_CORE(x_, z_) (2.0 / omega * exp(-0.5 * SQR(z_) + 0.5 * SQR(x_) + inla_log_Phi(alpha * (z_))))
#define RATIO(x_) RATIO_CORE(x_, (((x_)-xi)/omega))

	// Return a new allocted _snq_tp object with the required information. Note that 'n' in the object can be less than the
	// 'n' in the function call.

	// the weights wp[i+n] and wp[i+2*n] gives the weight to get the first and second derivative wrt the skew3

	// New memory for xp and wp is allocated

	int i, j, k;
	double *xxp = NULL, *wwp = NULL;

	double c1 = 1.2533141373155001208;		       // = sqrt(M_PI/2)
	double c2 = 0.63661977236758138243;		       // 2/M_PI
	double c3 = 0.568995659411924537;		       // pow((4 - M_PI)/2.0, 2.0/3.0)));
	double v1, delta, alpha, omega, xi;

	double *work = Calloc(4 * n, double);
	double *nodes = work;
	double *w = work + n;
	double *w_grad = work + 2 * n;
	double *w_hess = work + 3 * n;

	GMRFLib_ghq(&xxp, &wwp, n);
	assert(xxp);
	assert(wwp);
	Memcpy(nodes, xxp, n * sizeof(double));
	Memcpy(w, wwp, n * sizeof(double));

	// stencil for the first and second derivative. degree 5, 7 and 9
	// double wf[] = { 1.0 / 12.0, -2.0 / 3.0, 0.0, 2.0 / 3.0, -1.0 / 12.0 };
	// double wff[] = { -1.0 / 12.0, 4.0 / 3.0, -5.0 / 2.0, 4.0 / 3.0, -1.0 / 12.0 };
	// double wf[] = { -1.0 / 60.0, 3.0 / 20.0, -3.0 / 4.0, 0.0, 3.0 / 4.0, -3.0 / 20.0, 1.0 / 60.0 };
	// double wff[] = { 1.0 / 90.0, -3.0 / 20.0, 3.0 / 2.0, -49.0 / 18.0, 3.0 / 2.0, -3.0 / 20.0, 1.0 / 90.0 };
	// double wfff[] = { 1.0 / 8.0, -1.0, 13.0 / 8.0, 0.0, -13.0 / 8.0, 1.0, -1.0 / 8.0 };
	double wf[] = { 1.0 / 280.0, -4.0 / 105.0, 1.0 / 5.0, -4.0 / 5.0, 0.0, 4.0 / 5.0, -1.0 / 5.0, 4.0 / 105.0, -1.0 / 280.0 };
	double wff[] = { -1.0 / 560.0, 8.0 / 315.0, -1.0 / 5.0, 8.0 / 5.0, -205.0 / 72.0, 8.0 / 5.0, -1.0 / 5.0, 8.0 / 315.0, -1.0 / 560.0 };
	// double wfff[] = { -7.0 / 240.0, 3.0 / 10.0, -169.0 / 120.0, 61.0 / 30.0, 0.0, -61.0 / 30.0, 169.0 / 120.0, -3.0 / 10.0, 7.0 / 240.0 };

	double skew3s[sizeof(wf) / sizeof(double)], s;

	double ds = GSL_ROOT4_DBL_EPSILON;
	double **ww = NULL;
	int ns = sizeof(skew3s) / sizeof(double);
	int n2 = ns / 2;

	ww = Calloc(ns, double *);
	ww[0] = Calloc(ns * n, double);
	for (j = 1; j < ns; j++) {
		ww[j] = ww[0] + j * n;
	}

	for (j = 0, i = -n2; j < ns; j++, i++) {
		skew3s[j] = skew3 + i * ds;
	}

	for (j = 0; j < ns; j++) {
		s = GMRFLib_skew3_to_skew(skew3s[j]);
		v1 = pow(ABS(s), 2.0 / 3.0);
		delta = c1 * sqrt(v1 / (v1 + c3)) * SIGN(s);
		alpha = delta / sqrt(1.0 - SQR(delta));
		omega = sqrt(1.0 / (1.0 - c2 * SQR(delta)));
		xi = 0.0 - omega * delta / c1;

		for (i = 0; i < n; i++) {
			ww[j][i] = wwp[i] * RATIO(xxp[i]);
		}
	}

	for (i = 0; i < n; i++) {
		w[i] = ww[n2][i];
		w_grad[i] = 0.0;
		w_hess[i] = 0.0;
		for (j = 0; j < ns; j++) {
			w_grad[i] += wf[j] * ww[j][i];
			w_hess[i] += wff[j] * ww[j][i];
		}
		w_grad[i] /= ds;
		w_hess[i] /= SQR(ds);
	}

	double w_max = GMRFLib_max_value(w, n, NULL);
	double limit = GSL_SQRT_DBL_EPSILON;

	for (i = k = 0; i < n; i++) {
		if (w[i] / w_max > limit) {
			nodes[k] = nodes[i];
			w[k] = w[i];
			w_grad[k] = w_grad[i];
			w_hess[k] = w_hess[i];
			k++;
		}
	}

	GMRFLib_snq_tp *snq = Calloc(1, GMRFLib_snq_tp);
	snq->n = k;
	snq->skew3 = skew3;
	snq->nodes = Calloc(4 * snq->n, double);
	snq->w = snq->nodes + snq->n;
	snq->w_grad = snq->nodes + 2 * snq->n;
	snq->w_hess = snq->nodes + 3 * snq->n;

	Memcpy(snq->nodes, nodes, snq->n * sizeof(double));
	Memcpy(snq->w, w, snq->n * sizeof(double));
	Memcpy(snq->w_grad, w_grad, snq->n * sizeof(double));
	Memcpy(snq->w_hess, w_hess, snq->n * sizeof(double));

	double tmp = 0.0;
	for (i = 0; i < snq->n; i++) {
		tmp += snq->w[i];
	}
	tmp = 1.0 / tmp;
	for (i = 0; i < snq->n; i++) {
		snq->w[i] *= tmp;
		snq->w_grad[i] *= tmp;
		snq->w_hess[i] *= tmp;
	}

	Free(ww[0]);
	Free(ww);
	Free(work);

#undef RATIO_CORE
#undef RATIO
	return snq;
}

int GMRFLib_snq_free(GMRFLib_snq_tp *q)
{
	if (q) {
		Free(q->nodes);
		Free(q);
	}
	return GMRFLib_SUCCESS;
}
