
/* bessel.c
 * 
 * Copyright (C) 2012 Havard Rue
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
#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = HGVERSION;

#include <float.h>
#include <math.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"
#include "bessel.h"

/*
 * All this code is from R, slighty edited so its standalone
 */

double bessel_Knu(double alpha, double x)
{
	long nb, ncalc, ize;
	double *bk;

	ize = (long) 1;
	if (alpha < 0)
		alpha = -alpha;
	nb = 1 + (long) floor(alpha);			       /* nb-1 <= |alpha| < nb */
	alpha -= (nb - 1);
	bk = (double *) calloc(nb, sizeof(double));
	bessel_K(&x, &alpha, &nb, &ize, bk, &ncalc);
	x = bk[nb - 1];
	free(bk);
	return x;
}
double bessel_ftrunc(double x)
{
	if (x >= 0)
		return floor(x);
	else
		return ceil(x);
}
void bessel_K(double *x, double *alpha, long *nb, long *ize, double *bk, long *ncalc)
{
#define xmax_BESS_K	705.342				       /* maximal x for UNscaled answer */
#define sqxmin_BESS_K	1.49e-154
#define ML_POSINF DBL_MAX
#define M_SQRT_2dPI	0.797884560802865355879892119869       /* sqrt(2/pi) */

	/*---------------------------------------------------------------------
	 * Mathematical constants
	 *	A = LOG(2) - Euler's constant
	 *	D = SQRT(2/PI)
	 ---------------------------------------------------------------------*/
	static double a = .11593151565841244881;

	/*---------------------------------------------------------------------
	  P, Q - Approximation for LOG(GAMMA(1+ALPHA))/ALPHA + Euler's constant
	  Coefficients converted from hex to decimal and modified
	  by W. J. Cody, 2/26/82 */
	static double p[8] = { .805629875690432845, 20.4045500205365151,
		157.705605106676174, 536.671116469207504, 900.382759291288778,
		730.923886650660393, 229.299301509425145, .822467033424113231
	};
	static double q[7] = { 29.4601986247850434, 277.577868510221208,
		1206.70325591027438, 2762.91444159791519, 3443.74050506564618,
		2210.63190113378647, 572.267338359892221
	};
	/*
	 * R, S - Approximation for (1-ALPHA*PI/SIN(ALPHA*PI))/(2.D0*ALPHA) 
	 */
	static double r[5] = { -.48672575865218401848, 13.079485869097804016,
		-101.96490580880537526, 347.65409106507813131,
		3.495898124521934782e-4
	};
	static double s[4] = { -25.579105509976461286, 212.57260432226544008,
		-610.69018684944109624, 422.69668805777760407
	};
	/*
	 * T - Approximation for SINH(Y)/Y 
	 */
	static double t[6] = { 1.6125990452916363814e-10,
		2.5051878502858255354e-8, 2.7557319615147964774e-6,
		1.9841269840928373686e-4, .0083333333333334751799,
		.16666666666666666446
	};

	/*---------------------------------------------------------------------*/
	static double estm[6] = { 52.0583, 5.7607, 2.7782, 14.4303, 185.3004, 9.3715 };
	static double estf[7] = { 41.8341, 7.1075, 6.4306, 42.511, 1.35633, 84.5096, 20. };

	/*
	 * Local variables 
	 */
	long iend, i, j, k, m, ii, mplus1;
	double x2by4, twox, c, blpha, ratio, wminf;
	double d1, d2, d3, f0, f1, f2, p0, q0, t1, t2, twonu;
	double dm, ex, bk1, bk2, nu;

	ii = 0;						       /* -Wall */

	ex = *x;
	nu = *alpha;
	*ncalc = IMIN(*nb, 0) - 2;
	if (*nb > 0 && (0. <= nu && nu < 1.) && (1 <= *ize && *ize <= 2)) {
		if (ex <= 0 || (*ize == 1 && ex > xmax_BESS_K)) {
			if (ex <= 0) {
				if (ex < 0)
					abort();
				for (i = 0; i < *nb; i++)
					bk[i] = ML_POSINF;
			} else				       /* would only have underflow */
				for (i = 0; i < *nb; i++)
					bk[i] = 0.;
			*ncalc = *nb;
			return;
		}
		k = 0;
		if (nu < sqxmin_BESS_K) {
			nu = 0.;
		} else if (nu > .5) {
			k = 1;
			nu -= 1.;
		}
		twonu = nu + nu;
		iend = *nb + k - 1;
		c = nu * nu;
		d3 = -c;
		if (ex <= 1.) {
			/*
			 * ------------------------------------------------------------ Calculation of P0 = GAMMA(1+ALPHA) * (2/X)**ALPHA Q0 = GAMMA(1-ALPHA) *
			 * (X/2)**ALPHA ------------------------------------------------------------ 
			 */
			d1 = 0.;
			d2 = p[0];
			t1 = 1.;
			t2 = q[0];
			for (i = 2; i <= 7; i += 2) {
				d1 = c * d1 + p[i - 1];
				d2 = c * d2 + p[i];
				t1 = c * t1 + q[i - 1];
				t2 = c * t2 + q[i];
			}
			d1 = nu * d1;
			t1 = nu * t1;
			f1 = log(ex);
			f0 = a + nu * (p[7] - nu * (d1 + d2) / (t1 + t2)) - f1;
			q0 = exp(-nu * (a - nu * (p[7] + nu * (d1 - d2) / (t1 - t2)) - f1));
			f1 = nu * f0;
			p0 = exp(f1);
			/*
			 * ----------------------------------------------------------- Calculation of F0 =
			 * ----------------------------------------------------------- 
			 */
			d1 = r[4];
			t1 = 1.;
			for (i = 0; i < 4; ++i) {
				d1 = c * d1 + r[i];
				t1 = c * t1 + s[i];
			}
			/*
			 * d2 := sinh(f1)/ nu = sinh(f1)/(f1/f0) = f0 * sinh(f1)/f1 
			 */
			if (fabs(f1) <= .5) {
				f1 *= f1;
				d2 = 0.;
				for (i = 0; i < 6; ++i) {
					d2 = f1 * d2 + t[i];
				}
				d2 = f0 + f0 * f1 * d2;
			} else {
				d2 = sinh(f1) / nu;
			}
			f0 = d2 - nu * d1 / (t1 * p0);
			if (ex <= 1e-10) {
				/*
				 * --------------------------------------------------------- X <= 1.0E-10 Calculation of K(ALPHA,X) and X*K(ALPHA+1,X)/K(ALPHA,X)
				 * --------------------------------------------------------- 
				 */
				bk[0] = f0 + ex * f0;
				if (*ize == 1) {
					bk[0] -= ex * bk[0];
				}
				ratio = p0 / f0;
				c = ex * DBL_MAX;
				if (k != 0) {
					/*
					 * --------------------------------------------------- Calculation of K(ALPHA,X) and X*K(ALPHA+1,X)/K(ALPHA,X), ALPHA >=
					 * 1/2 --------------------------------------------------- 
					 */
					*ncalc = -1;
					if (bk[0] >= c / ratio) {
						return;
					}
					bk[0] = ratio * bk[0] / ex;
					twonu += 2.;
					ratio = twonu;
				}
				*ncalc = 1;
				if (*nb == 1)
					return;

				/*
				 * ----------------------------------------------------- Calculate K(ALPHA+L,X)/K(ALPHA+L-1,X), L = 1, 2, ... , NB-1
				 * ----------------------------------------------------- 
				 */
				*ncalc = -1;
				for (i = 1; i < *nb; ++i) {
					if (ratio >= c)
						return;

					bk[i] = ratio / ex;
					twonu += 2.;
					ratio = twonu;
				}
				*ncalc = 1;
				goto L420;
			} else {
				/*
				 * ------------------------------------------------------ 10^-10 < X <= 1.0 ------------------------------------------------------ 
				 */
				c = 1.;
				x2by4 = ex * ex / 4.;
				p0 = .5 * p0;
				q0 = .5 * q0;
				d1 = -1.;
				d2 = 0.;
				bk1 = 0.;
				bk2 = 0.;
				f1 = f0;
				f2 = p0;
				do {
					d1 += 2.;
					d2 += 1.;
					d3 = d1 + d3;
					c = x2by4 * c / d2;
					f0 = (d2 * f0 + p0 + q0) / d3;
					p0 /= d2 - nu;
					q0 /= d2 + nu;
					t1 = c * f0;
					t2 = c * (p0 - d2 * f0);
					bk1 += t1;
					bk2 += t2;
				} while (fabs(t1 / (f1 + bk1)) > DBL_EPSILON || fabs(t2 / (f2 + bk2)) > DBL_EPSILON);
				bk1 = f1 + bk1;
				bk2 = 2. * (f2 + bk2) / ex;
				if (*ize == 2) {
					d1 = exp(ex);
					bk1 *= d1;
					bk2 *= d1;
				}
				wminf = estf[0] * ex + estf[1];
			}
		} else if (DBL_EPSILON * ex > 1.) {
			/*
			 * ------------------------------------------------- X > 1./EPS ------------------------------------------------- 
			 */
			*ncalc = *nb;
			bk1 = 1. / (M_SQRT_2dPI * sqrt(ex));
			for (i = 0; i < *nb; ++i)
				bk[i] = bk1;
			return;

		} else {
			/*
			 * ------------------------------------------------------- X > 1.0 ------------------------------------------------------- 
			 */
			twox = ex + ex;
			blpha = 0.;
			ratio = 0.;
			if (ex <= 4.) {
				/*
				 * ---------------------------------------------------------- Calculation of K(ALPHA+1,X)/K(ALPHA,X), 1.0 <= X <= 4.0
				 * ----------------------------------------------------------
				 */
				d2 = bessel_ftrunc(estm[0] / ex + estm[1]);
				m = (long) d2;
				d1 = d2 + d2;
				d2 -= .5;
				d2 *= d2;
				for (i = 2; i <= m; ++i) {
					d1 -= 2.;
					d2 -= d1;
					ratio = (d3 + d2) / (twox + d1 - ratio);
				}
				/*
				 * ----------------------------------------------------------- Calculation of I(|ALPHA|,X) and I(|ALPHA|+1,X) by backward
				 * recurrence and K(ALPHA,X) from the wronskian -----------------------------------------------------------
				 */
				d2 = bessel_ftrunc(estm[2] * ex + estm[3]);
				m = (long) d2;
				c = fabs(nu);
				d3 = c + c;
				d1 = d3 - 1.;
				f1 = DBL_MIN;
				f0 = (2. * (c + d2) / ex + .5 * ex / (c + d2 + 1.)) * DBL_MIN;
				for (i = 3; i <= m; ++i) {
					d2 -= 1.;
					f2 = (d3 + d2 + d2) * f0;
					blpha = (1. + d1 / d2) * (f2 + blpha);
					f2 = f2 / ex + f1;
					f1 = f0;
					f0 = f2;
				}
				f1 = (d3 + 2.) * f0 / ex + f1;
				d1 = 0.;
				t1 = 1.;
				for (i = 1; i <= 7; ++i) {
					d1 = c * d1 + p[i - 1];
					t1 = c * t1 + q[i - 1];
				}
				p0 = exp(c * (a + c * (p[7] - c * d1 / t1) - log(ex))) / ex;
				f2 = (c + .5 - ratio) * f1 / ex;
				bk1 = p0 + (d3 * f0 - f2 + f0 + blpha) / (f2 + f1 + f0) * p0;
				if (*ize == 1) {
					bk1 *= exp(-ex);
				}
				wminf = estf[2] * ex + estf[3];
			} else {
				/*
				 * --------------------------------------------------------- Calculation of K(ALPHA,X) and K(ALPHA+1,X)/K(ALPHA,X), by backward
				 * recurrence, for X > 4.0 ----------------------------------------------------------
				 */
				dm = bessel_ftrunc(estm[4] / ex + estm[5]);
				m = (long) dm;
				d2 = dm - .5;
				d2 *= d2;
				d1 = dm + dm;
				for (i = 2; i <= m; ++i) {
					dm -= 1.;
					d1 -= 2.;
					d2 -= d1;
					ratio = (d3 + d2) / (twox + d1 - ratio);
					blpha = (ratio + ratio * blpha) / dm;
				}
				bk1 = 1. / ((M_SQRT_2dPI + M_SQRT_2dPI * blpha) * sqrt(ex));
				if (*ize == 1)
					bk1 *= exp(-ex);
				wminf = estf[4] * (ex - fabs(ex - estf[6])) + estf[5];
			}
			/*
			 * --------------------------------------------------------- Calculation of K(ALPHA+1,X) from K(ALPHA,X) and K(ALPHA+1,X)/K(ALPHA,X)
			 * --------------------------------------------------------- 
			 */
			bk2 = bk1 + bk1 * (nu + .5 - ratio) / ex;
		}

		/*--------------------------------------------------------------------
		  Calculation of 'NCALC', K(ALPHA+I,X),	I  =  0, 1, ... , NCALC-1,
		  &	  K(ALPHA+I,X)/K(ALPHA+I-1,X),	I = NCALC, NCALC+1, ... , NB-1
		  -------------------------------------------------------------------*/
		*ncalc = *nb;
		bk[0] = bk1;
		if (iend == 0)
			return;

		j = 1 - k;
		if (j >= 0)
			bk[j] = bk2;

		if (iend == 1)
			return;

		m = IMIN((long) (wminf - nu), iend);
		for (i = 2; i <= m; ++i) {
			t1 = bk1;
			bk1 = bk2;
			twonu += 2.;
			if (ex < 1.) {
				if (bk1 >= DBL_MAX / twonu * ex)
					break;
			} else {
				if (bk1 / ex >= DBL_MAX / twonu)
					break;
			}
			bk2 = twonu / ex * bk1 + t1;
			ii = i;
			++j;
			if (j >= 0) {
				bk[j] = bk2;
			}
		}

		m = ii;
		if (m == iend) {
			return;
		}
		ratio = bk2 / bk1;
		mplus1 = m + 1;
		*ncalc = -1;
		for (i = mplus1; i <= iend; ++i) {
			twonu += 2.;
			ratio = twonu / ex + 1. / ratio;
			++j;
			if (j >= 1) {
				bk[j] = ratio;
			} else {
				if (bk2 >= DBL_MAX / ratio)
					return;

				bk2 *= ratio;
			}
		}
		*ncalc = IMAX(1, mplus1 - k);
		if (*ncalc == 1)
			bk[0] = bk2;
		if (*nb == 1)
			return;

	      L420:
		for (i = *ncalc; i < *nb; ++i) {	       /* i == *ncalc */
#ifndef IEEE_754
			if (bk[i - 1] >= DBL_MAX / bk[i])
				return;
#endif
			bk[i] *= bk[i - 1];
			(*ncalc)++;
		}
	}
}
