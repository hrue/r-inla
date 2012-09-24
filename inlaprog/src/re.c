
/* re.c
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
 *       H{\aa}vard Rue
 *       Department of Mathematical Sciences
 *       The Norwegian University of Science and Technology
 *       N-7491 Trondheim, Norway
 *       Voice: +47-7359-3533    URL  : http://www.math.ntnu.no/~hrue  
 *       Fax  : +47-7359-3524    Email: havard.rue@math.ntnu.no
 *
 */
#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = HGVERSION;

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

#include "re.h"

#define BESSEL_KNU(_alpha, _x) bessel_Knu(_alpha, _x)
//#define BESSEL_KNU(_alpha, _x) gsl_sf_bessel_Knu(_alpha, _x)

#define Pq(_arg) (exp(0.25)/sqrt(8.0*M_PI) * (BESSEL_KNU(((_arg)+1.0)/2.0, 0.25) + BESSEL_KNU(ABS(((_arg)-1.0)/2.0), 0.25)))
#define M1(_epsilon, _delta) (sinh((_epsilon)/(_delta))*Pq(1.0/(_delta)))
#define M2(_epsilon, _delta) (0.5*(cosh(2.0*(_epsilon)/(_delta))*Pq(2.0/(_delta)) -1.0))
#define M3(_epsilon, _delta) (0.25*(sinh(3.0*(_epsilon)/(_delta))*Pq(3.0/(_delta)) - \
				    3.0*sinh((_epsilon)/(_delta))*Pq(1.0/(_delta))))
#define M4(_epsilon, _delta) (0.125*(cosh(4.0*(_epsilon)/(_delta))*Pq(4.0/(_delta)) - \
				     4.0*cosh(2.0*(_epsilon)/(_delta))*Pq(2.0/(_delta)) + 3.0))

#define XMATCH(x0,x1) (fabs((x0)-(x1)) == 0)
#define YMATCH(y0,y1) (fabs((y0)-(y1)) == 0)

#define REV(x, n) \
	if (1){						\
		double *tmp = Calloc(n, double);	\
		memcpy(tmp, x, (n)*sizeof(double));	\
		int ii;					\
		for(ii=0; ii < (n); ii++){		\
			x[ii] = tmp[(n)-1-ii];		\
		}					\
		Free(tmp);				\
	}

static re_sas_prior_tp *sas_prior_table = NULL;

double bessel_Knu(double alpha, double x)
{
	/*
	 * This is from R 
	 */
	long nb, ncalc, ize;
	double *bk;

	ize = (long) 1;
	if (alpha < 0)
		alpha = -alpha;
	nb = 1 + (long) floor(alpha);			       /* nb-1 <= |alpha| < nb */
	alpha -= (nb - 1);
	bk = (double *) calloc(nb, sizeof(double));
	K_bessel(&x, &alpha, &nb, &ize, bk, &ncalc);
	x = bk[nb - 1];
	free(bk);
	return x;
}
double ftrunc(double x)
{
	/*
	 * This is from R 
	 */
	if (x >= 0)
		return floor(x);
	else
		return ceil(x);
}
void K_bessel(double *x, double *alpha, long *nb, long *ize, double *bk, long *ncalc)
{
	/*
	 * This is from R 
	 */

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
				d2 = ftrunc(estm[0] / ex + estm[1]);
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
				d2 = ftrunc(estm[2] * ex + estm[3]);
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
				dm = ftrunc(estm[4] / ex + estm[5]);
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


int re_valid_skew_kurt(double *dist, double skew, double kurt)
{
#define KURT_LIMIT(s) (2.15 + SQR((s)/0.8))
	int retval;

	retval = (kurt > KURT_LIMIT(skew) ? GMRFLib_TRUE : GMRFLib_FALSE);
	if (dist) {
		if (retval == GMRFLib_FALSE) {
			*dist = ABS(kurt - KURT_LIMIT(skew));
		} else {
			*dist = 0.0;
		}
	}
#undef KURT_LIMIT
	return retval;
}
int re_shash_skew_kurt(double *skew, double *kurt, double epsilon, double delta)
{
	double m1 = M1(epsilon, delta);
	double m2 = M2(epsilon, delta);
	double m3 = M3(epsilon, delta);
	double m4 = M4(epsilon, delta);

	if (skew) {
		*skew = (m3 - 3.0 * m1 * m2 + 3.0 * gsl_pow_3(m1) - gsl_pow_3(m1)) / gsl_pow_3(sqrt(m2 - SQR(m1)));
	}
	if (kurt) {
		*kurt = (m4 - 4.0 * m1 * m3 + 6.0 * SQR(m1) * m2 - 4.0 * gsl_pow_4(m1) + gsl_pow_4(m1)) / gsl_pow_4(sqrt(m2 - SQR(m1)));
	}

	return GMRFLib_SUCCESS;
}


int re_shash_fit_parameters(re_shash_param_tp * param, double *mean, double *prec, double *skew, double *kurt)
{
	/*
	 * for given mean, prec, skew and kurt, fit the parameters in the shash-model. Either none or both of 'skew' and 'kurt' can be NULL. 
	 */

	int use_lookup = 1;

	static double **pin = NULL;
	static double **pout = NULL;
	static int *perr = NULL;

	assert((!skew && !kurt) || (skew && kurt));	       /* either non or both */

	if (use_lookup) {
		if (!pin && !pout && !perr) {
#pragma omp critical
			{
				FIXME("ENTER WAIT GO");
				if (!pin && !pout && !perr) {
					pin = Calloc(ISQR(GMRFLib_MAX_THREADS), double *);
					pout = Calloc(ISQR(GMRFLib_MAX_THREADS), double *);
					perr = Calloc(ISQR(GMRFLib_MAX_THREADS), int);

					int i;
					for (i = 0; i < ISQR(GMRFLib_MAX_THREADS); i++) {
						pin[i] = Calloc(2, double);
						pout[i] = Calloc(2, double);
					}
				}
			}
		}
	}

	if (!pin && !pout && !perr) {
		abort();
	}

	double npin[2], pout_tmp[2];
	int id, err = 0, debug = 0;

	id = omp_get_thread_num() + GMRFLib_thread_id * GMRFLib_MAX_THREADS;
	npin[0] = (skew ? *skew : 0.0);
	npin[1] = (kurt ? *kurt : 3.0);

	if (use_lookup && memcmp(npin, pin[id], sizeof(npin)) == 0) {
		/*
		 * we have already computed these values
		 */
		memcpy(pout_tmp, pout[id], sizeof(pout_tmp));
	} else if (!skew && !kurt) {
		pout_tmp[0] = 0.0;			       /* epsilon */
		pout_tmp[1] = log(1.0);			       /* log delta */
	} else {
		gsl_vector *x = gsl_vector_alloc(2);
		gsl_vector_set(x, 0, 0);		       /* initial values */
		gsl_vector_set(x, 1, log(1.0));

		double target[2];
		target[0] = *skew;
		target[1] = *kurt;

		const gsl_multifit_fdfsolver_type *T;
		gsl_multifit_fdfsolver *s;
		int n = 2;
		int p = 2;
		int status;

		gsl_multifit_function_fdf f;
		f.f = re_shash_f;
		f.df = re_shash_df;
		f.fdf = re_shash_fdf;
		f.n = n;
		f.p = p;
		f.params = (void *) target;


		T = gsl_multifit_fdfsolver_lmsder;
		s = gsl_multifit_fdfsolver_alloc(T, n, p);
		gsl_multifit_fdfsolver_set(s, &f, (const gsl_vector *) x);

		int iter = 0, iter_max = 1000;
		double eps = 1e-6;
		do {
			iter++;
			status = gsl_multifit_fdfsolver_iterate(s);
			if (status) {
				break;
			}
			if (debug) {
				printf("status = %s\n", gsl_strerror(status));
				printf("iter %1d: x %g %g |f| = %g\n", iter, gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1), gsl_blas_dnrm2(s->f));
			}
			status = gsl_multifit_test_delta(s->dx, s->x, eps, eps);
		}
		while (status == GSL_CONTINUE && iter < iter_max);

		err = (gsl_blas_dnrm2(s->f) > 1e-5);
		if (err) {
			if (re_valid_skew_kurt(NULL, target[0], target[1]) == GMRFLib_TRUE) {
				fprintf(stderr, "SHASH fail to fit target skew=%g kurt=%g err=%g, but VALID FAIL!\n", target[0], target[1], gsl_blas_dnrm2(s->f));
			}
		}

		pout_tmp[0] = gsl_vector_get(s->x, 0);
		pout_tmp[1] = gsl_vector_get(s->x, 1);
		gsl_vector_free(x);
		gsl_multifit_fdfsolver_free(s);
	}

	double m1, m2, mu, stdev, ss;

	pout_tmp[1] = exp(pout_tmp[1]);			       /* transform it */
	m1 = M1(pout_tmp[0], pout_tmp[1]);
	m2 = M2(pout_tmp[0], pout_tmp[1]);
	ss = sqrt(m2 - SQR(m1));
	mu = *mean - m1 / sqrt(*prec) / ss;
	stdev = 1.0 / sqrt(*prec) / ss;

	param->mu = mu;
	param->stdev = stdev;
	param->epsilon = pout_tmp[0];
	param->delta = pout_tmp[1];

	if (use_lookup) {
		perr[id] = err;
		memcpy(pin[id], npin, sizeof(npin));
		memcpy(pout[id], pout_tmp, sizeof(pout_tmp));
	}

	return (err ? !GMRFLib_SUCCESS : GMRFLib_SUCCESS);
}
int re_shash_f(const gsl_vector * x, void *data, gsl_vector * f)
{
	double epsilon, delta, skew, sskew, kurt, kkurt, *ddata = (double *) data;

	epsilon = gsl_vector_get(x, 0);
	delta = exp(gsl_vector_get(x, 1));
	sskew = ddata[0];
	kkurt = ddata[1];

	re_shash_skew_kurt(&skew, &kurt, epsilon, delta);
	gsl_vector_set(f, 0, skew - sskew);
	gsl_vector_set(f, 1, kurt - kkurt);

	return GSL_SUCCESS;
}
int re_shash_df(const gsl_vector * x, void *data, gsl_matrix * J)
{
	double h = GMRFLib_eps(1. / 3.), epsilon, log_delta, skew_h, skew_mh, kurt_h, kurt_mh;

	epsilon = gsl_vector_get(x, 0);
	log_delta = gsl_vector_get(x, 1);

	re_shash_skew_kurt(&skew_h, &kurt_h, epsilon + h, exp(log_delta));
	re_shash_skew_kurt(&skew_mh, &kurt_mh, epsilon - h, exp(log_delta));
	gsl_matrix_set(J, 0, 0, (skew_h - skew_mh) / 2.0 / h);
	gsl_matrix_set(J, 1, 0, (kurt_h - kurt_mh) / 2.0 / h);

	re_shash_skew_kurt(&skew_h, &kurt_h, epsilon, exp(log_delta + h));
	re_shash_skew_kurt(&skew_mh, &kurt_mh, epsilon, exp(log_delta - h));
	gsl_matrix_set(J, 0, 1, (skew_h - skew_mh) / 2.0 / h);
	gsl_matrix_set(J, 1, 1, (kurt_h - kurt_mh) / 2.0 / h);

	return GSL_SUCCESS;
}
int re_shash_fdf(const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J)
{
	re_shash_f(x, data, f);
	re_shash_df(x, data, J);

	return GSL_SUCCESS;
}

/* 
 *  The countourLine-code is taken from R: plot3d.c, with modifications.
 */
int ctr_intersect(double z0, double z1, double zc, double *f)
{
	if ((z0 - zc) * (z1 - zc) < 0.0) {
		*f = (zc - z0) / (z1 - z0);
		return 1;
	}
	return 0;
}
CLP ctr_newseg(double x0, double y0, double x1, double y1, CLP prev)
{
	CLP seg = (CLP) calloc(1, sizeof(CL));
	seg->x0 = x0;
	seg->y0 = y0;
	seg->x1 = x1;
	seg->y1 = y1;
	seg->next = prev;
	return seg;
}

void ctr_swapseg(CLP seg)
{
	double x, y;
	x = seg->x0;
	y = seg->y0;
	seg->x0 = seg->x1;
	seg->y0 = seg->y1;
	seg->x1 = x;
	seg->y1 = y;
}

int ctr_segdir(double xend, double yend, double *x, double *y, int *i, int *j, int nx, int ny)
{
	if (YMATCH(yend, y[*j])) {
		if (*j == 0)
			return 0;
		*j = *j - 1;
		return 3;
	}
	if (XMATCH(xend, x[*i])) {
		if (*i == 0)
			return 0;
		*i = *i - 1;
		return 4;
	}
	if (YMATCH(yend, y[*j + 1])) {
		if (*j >= ny - 1)
			return 0;
		*j = *j + 1;
		return 1;
	}
	if (XMATCH(xend, x[*i + 1])) {
		if (*i >= nx - 1)
			return 0;
		*i = *i + 1;
		return 2;
	}
	return 0;
}

CLP ctr_segupdate(double xend, double yend, int dir, int tail, CLP seglist, CLP * seg)
{
	if (seglist == NULL) {
		*seg = NULL;
		return NULL;
	}
	switch (dir) {
	case 1:
	case 3:
		if (YMATCH(yend, seglist->y0)) {
			if (!tail)
				ctr_swapseg(seglist);
			*seg = seglist;
			return seglist->next;
		}
		if (YMATCH(yend, seglist->y1)) {
			if (tail)
				ctr_swapseg(seglist);
			*seg = seglist;
			return seglist->next;
		}
		break;
	case 2:
	case 4:
		if (XMATCH(xend, seglist->x0)) {
			if (!tail)
				ctr_swapseg(seglist);
			*seg = seglist;
			return seglist->next;
		}
		if (XMATCH(xend, seglist->x1)) {
			if (tail)
				ctr_swapseg(seglist);
			*seg = seglist;
			return seglist->next;
		}
		break;
	}
	seglist->next = ctr_segupdate(xend, yend, dir, tail, seglist->next, seg);
	return seglist;
}

CLP *contourLines1(double *x, int nx, double *y, int ny, double *z, double zc)
{
	double f, xl, xh, yl, yh, zll, zhl, zlh, zhh, xx[4], yy[4];
	double atom = GMRFLib_eps(0.5);
	int i, j, k, l, m, nacode;
	CLP seglist;
	CLP *segmentDB;
	segmentDB = (CLP *) calloc(nx * ny, sizeof(CLP));
	for (i = 0; i < nx; i++)
		for (j = 0; j < ny; j++)
			segmentDB[i + j * nx] = NULL;
	for (i = 0; i < nx - 1; i++) {
		xl = x[i];
		xh = x[i + 1];
		for (j = 0; j < ny - 1; j++) {
			yl = y[j];
			yh = y[j + 1];
			k = i + j * nx;
			zll = z[k];
			zhl = z[k + 1];
			zlh = z[k + nx];
			zhh = z[k + nx + 1];

			/*
			 * If the value at a corner is exactly equal to a contour level, change that value by a tiny amount 
			 */
			if (zll == zc)
				zll += atom;
			if (zhl == zc)
				zhl += atom;
			if (zlh == zc)
				zlh += atom;
			if (zhh == zc)
				zhh += atom;
			/*
			 * Check for intersections with sides 
			 */
			nacode = 0;
			if (isfinite(zll))
				nacode += 1;
			if (isfinite(zhl))
				nacode += 2;
			if (isfinite(zlh))
				nacode += 4;
			if (isfinite(zhh))
				nacode += 8;

			k = 0;
			switch (nacode) {
			case 15:
				if (ctr_intersect(zll, zhl, zc, &f)) {
					xx[k] = xl + f * (xh - xl);
					yy[k] = yl;
					k++;
				}
				if (ctr_intersect(zll, zlh, zc, &f)) {
					yy[k] = yl + f * (yh - yl);
					xx[k] = xl;
					k++;
				}
				if (ctr_intersect(zhl, zhh, zc, &f)) {
					yy[k] = yl + f * (yh - yl);
					xx[k] = xh;
					k++;
				}
				if (ctr_intersect(zlh, zhh, zc, &f)) {
					xx[k] = xl + f * (xh - xl);
					yy[k] = yh;
					k++;
				}
				break;
			case 14:
				if (ctr_intersect(zhl, zhh, zc, &f)) {
					yy[k] = yl + f * (yh - yl);
					xx[k] = xh;
					k++;
				}
				if (ctr_intersect(zlh, zhh, zc, &f)) {
					xx[k] = xl + f * (xh - xl);
					yy[k] = yh;
					k++;
				}
				if (ctr_intersect(zlh, zhl, zc, &f)) {
					xx[k] = xl + f * (xh - xl);
					yy[k] = yh + f * (yl - yh);
					k++;
				}
				break;
			case 13:
				if (ctr_intersect(zll, zlh, zc, &f)) {
					yy[k] = yl + f * (yh - yl);
					xx[k] = xl;
					k++;
				}
				if (ctr_intersect(zlh, zhh, zc, &f)) {
					xx[k] = xl + f * (xh - xl);
					yy[k] = yh;
					k++;
				}
				if (ctr_intersect(zll, zhh, zc, &f)) {
					xx[k] = xl + f * (xh - xl);
					yy[k] = yl + f * (yh - yl);
					k++;
				}
				break;
			case 11:
				if (ctr_intersect(zhl, zhh, zc, &f)) {
					yy[k] = yl + f * (yh - yl);
					xx[k] = xh;
					k++;
				}
				if (ctr_intersect(zll, zhl, zc, &f)) {
					xx[k] = xl + f * (xh - xl);
					yy[k] = yl;
					k++;
				}
				if (ctr_intersect(zll, zhh, zc, &f)) {
					xx[k] = xl + f * (xh - xl);
					yy[k] = yl + f * (yh - yl);
					k++;
				}
				break;
			case 7:
				if (ctr_intersect(zll, zlh, zc, &f)) {
					yy[k] = yl + f * (yh - yl);
					xx[k] = xl;
					k++;
				}
				if (ctr_intersect(zll, zhl, zc, &f)) {
					xx[k] = xl + f * (xh - xl);
					yy[k] = yl;
					k++;
				}
				if (ctr_intersect(zlh, zhl, zc, &f)) {
					xx[k] = xl + f * (xh - xl);
					yy[k] = yh + f * (yl - yh);
					k++;
				}
				break;
			}

			/*
			 * We now have k(=2,4) endpoints 
			 * Decide which to join 
			 */
			seglist = NULL;

			if (k > 0) {
				if (k == 2) {
					seglist = ctr_newseg(xx[0], yy[0], xx[1], yy[1], seglist);
				} else if (k == 4) {
					for (k = 3; k >= 1; k--) {
						m = k;
						xl = xx[k];
						for (l = 0; l < k; l++) {
							if (xx[l] > xl) {
								xl = xx[l];
								m = l;
							}
						}
						if (m != k) {
							xl = xx[k];
							yl = yy[k];
							xx[k] = xx[m];
							yy[k] = yy[m];
							xx[m] = xl;
							yy[m] = yl;
						}
					}
					seglist = ctr_newseg(xx[0], yy[0], xx[1], yy[1], seglist);
					seglist = ctr_newseg(xx[2], yy[2], xx[3], yy[3], seglist);

				} else {
					int inla_error_general(const char *msg);
					inla_error_general("k != 2 or 4");
				}
			}
			segmentDB[i + j * nx] = seglist;
		}
	}
	return segmentDB;
}
inla_contour_tp *contourLines2(double *x, int nx, double *y, int ny, double *z, double zc, CLP * segmentDB)
{
	assert(segmentDB);

	double xend, yend;
	int i, ii, j, jj, ns, dir;
	CLP seglist, seg, s, start, end;

	/*
	 * Begin following contours: Grab a segment, Follow its tail, Follow its head, Save the contour 
	 */
	inla_contour_tp *c = Calloc(1, inla_contour_tp);
	c->level = zc;

	for (i = 0; i < nx - 1; i++) {
		for (j = 0; j < ny - 1; j++) {
			while ((seglist = segmentDB[i + j * nx])) {
				ii = i;
				jj = j;
				start = end = seglist;
				segmentDB[i + j * nx] = seglist->next;
				xend = seglist->x1;
				yend = seglist->y1;
				while ((dir = ctr_segdir(xend, yend, x, y, &ii, &jj, nx, ny))) {
					segmentDB[ii + jj * nx]
					    = ctr_segupdate(xend, yend, dir, GMRFLib_TRUE,	/* = tail */
							    segmentDB[ii + jj * nx], &seg);
					if (!seg)
						break;
					end->next = seg;
					end = seg;
					xend = end->x1;
					yend = end->y1;
				}
				end->next = NULL;	       /* <<< new for 1.2.3 */
				ii = i;
				jj = j;
				xend = seglist->x0;
				yend = seglist->y0;
				while ((dir = ctr_segdir(xend, yend, x, y, &ii, &jj, nx, ny))) {
					segmentDB[ii + jj * nx]
					    = ctr_segupdate(xend, yend, dir, GMRFLib_FALSE,	/* ie. head */
							    segmentDB[ii + jj * nx], &seg);
					if (!seg)
						break;
					seg->next = start;
					start = seg;
					xend = start->x0;
					yend = start->y0;
				}

				/*
				 * ns := #{segments of polyline} -- need to allocate 
				 */
				s = start;
				ns = 0;
				while (s && ns < INT_MAX) {
					ns++;
					s = s->next;
				}
				/*
				 * "write" the contour locations into the list of contours
				 */
				double *x, *y;

				x = Calloc(ns + 1, double);
				y = Calloc(ns + 1, double);
				s = start;
				x[0] = s->x0;
				y[0] = s->y0;
				ns = 1;
				while (s->next) {
					CLP s_free = s;	       /* need to save it, to free it */
					s = s->next;
					Free(s_free);
					x[ns] = s->x0;
					y[ns] = s->y0;
					ns++;
				}
				x[ns] = s->x1;
				y[ns] = s->y1;
				Free(s);

				if (0) {
					int ii;
					for (ii = 0; ii < ns + 1; ii++) {
						printf("(x,y) = ( %g %g )\n", x[ii], y[ii]);
					}
				}

				int ii;
				double len = 0.0;
				for (ii = 1; ii < ns + 1; ii++) {
					len += sqrt(SQR(x[ii] - x[ii - 1]) + SQR(y[ii] - y[ii - 1]));
				}

				c->nc++;
				c->ns = Realloc(c->ns, c->nc, int);
				c->cyclic = Realloc(c->cyclic, c->nc, int);
				c->length = Realloc(c->length, c->nc, double);
				c->x = Realloc(c->x, c->nc, double *);
				c->y = Realloc(c->y, c->nc, double *);
				c->ns[c->nc - 1] = ns;
				c->cyclic[c->nc - 1] = (x[0] == x[ns]) && (y[0] == y[ns]);
				c->length[c->nc - 1] = len;
				c->x[c->nc - 1] = x;
				c->y[c->nc - 1] = y;
			}
		}
	}
	return c;
}
inla_contour_tp *contourLines(double *x, int nx, double *y, int ny, double *z, double zc)
{
	CLP *seg;
	inla_contour_tp *c;
	seg = contourLines1(x, nx, y, ny, z, zc);
	c = contourLines2(x, nx, y, ny, z, zc, seg);
	Free(seg);

	return c;
}
int inla_print_contourLines(FILE * fp, inla_contour_tp * c)
{
	int i, j;

	if (fp == NULL) {
		fp = stdout;
	}
	fprintf(fp, "Number of contours %1d for level %g\n", c->nc, c->level);
	for (i = 0; i < c->nc; i++) {
		fprintf(fp, "\tContour[%1d]\n", i);
		fprintf(fp, "\tNumber of segments = %1d\n", c->ns[i]);
		fprintf(fp, "\tLength = %.12g\n", c->length[i]);
		fprintf(fp, "\tCyclic = %s\n", (c->cyclic[i] ? "Yes" : "No"));
		fprintf(fp, "\tPoints:\n");
		for (j = 0; j < c->ns[i] + 1; j++) {
			fprintf(fp, "\t\t(x,y) = %g %g\n", c->x[i][j], c->y[i][j]);
		}
	}

	return GMRFLib_SUCCESS;
}
int inla_free_contourLines(inla_contour_tp * c)
{
	if (!c)
		return GMRFLib_SUCCESS;

	int i;
	for (i = 0; i < c->nc; i++) {
		Free(c->x[i]);
		Free(c->y[i]);
	}
	Free(c->x);
	Free(c->y);
	Free(c->ns);
	Free(c->length);
	Free(c->cyclic);
	Free(c);

	return GMRFLib_SUCCESS;
}


// 

int re_join_contourLines(inla_contour_tp * c)
{
	/*
	 * modify the contourLines in 'c' so its one object 
	 */

	if (c->nc <= 1) {
		return GMRFLib_SUCCESS;
	}

	int debug = 0;

	if (debug) {
		printf("ENTER WITH\n");
		inla_print_contourLines(stdout, c);
	}

	int i, j, k;
	double *dists = Calloc(ISQR(c->nc), double);
	int *idir = Calloc(ISQR(c->nc), int);
	int *jdir = Calloc(ISQR(c->nc), int);

	if (debug) {
		P(c->nc);
	}

	for (i = 0; i < c->nc; i++) {
		for (j = 0; j < c->nc; j++) {
			k = i + j * c->nc;
			if (i < j) {
				double d, dd;

				// i ... j
				d = sqrt(SQR(c->x[i][c->ns[i]] - c->x[j][0]) + SQR(c->y[i][c->ns[i]] - c->y[j][0]));
				idir[k] = 1;
				jdir[k] = 1;

				// i ... -j
				dd = sqrt(SQR(c->x[i][c->ns[i]] - c->x[j][c->ns[j]]) + SQR(c->y[i][c->ns[i]] - c->y[j][c->ns[j]]));
				if (dd < d) {
					idir[k] = 1;
					jdir[k] = -1;
					d = dd;
				}
				// -i .. j
				dd = sqrt(SQR(c->x[i][0] - c->x[j][0]) + SQR(c->y[i][0] - c->y[j][0]));
				if (dd < d) {
					idir[k] = -1;
					jdir[k] = 1;
					d = dd;
				}
				// -i .. -j
				dd = sqrt(SQR(c->x[i][0] - c->x[j][c->ns[j]]) + SQR(c->y[i][0] - c->y[j][c->ns[j]]));
				if (dd < d) {
					idir[k] = -1;
					jdir[k] = -1;
					d = dd;
				}

				dists[k] = d;
				if (debug) {
					printf("compute dists[ %d %d ] = %d %d %g\n", i, j, idir[k], jdir[k], dists[k]);
				}
			}
		}
	}

	int imin, jmin, imindir, jmindir;
	double dmin = DBL_MAX;

	for (i = 0; i < c->nc; i++) {
		for (j = 0; j < c->nc; j++) {
			k = i + j * c->nc;
			if (i < j) {
				if (dists[k] < dmin) {
					dmin = dists[k];
					imin = i;
					jmin = j;
					imindir = idir[k];
					jmindir = jdir[k];
				}
			}
		}
	}
	if (debug) {
		printf("found min imin=%d jmin=%d imindir=%d jmindir=%d dist=%g\n", imin, jmin, imindir, jmindir, dmin);
	}

	c->x[imin] = Realloc(c->x[imin], c->ns[imin] + c->ns[jmin] + 2, double);
	c->y[imin] = Realloc(c->y[imin], c->ns[imin] + c->ns[jmin] + 2, double);

	if (debug) {
		P(c->ns[imin]);
		P(c->ns[jmin]);
	}
	if (imindir < 0) {
		REV(c->x[imin], c->ns[imin] + 1);
		REV(c->y[imin], c->ns[imin] + 1);
	}
	if (jmindir < 0) {
		REV(c->x[jmin], c->ns[jmin] + 1);
		REV(c->y[jmin], c->ns[jmin] + 1);
	}

	if (debug) {
		printf("IMIN\n");
		for (i = 0; i < c->ns[imin] + 1; i++)
			printf("imin %d (%g %g)\n", i, c->x[imin][i], c->y[imin][i]);

		printf("JMIN\n");
		for (i = 0; i < c->ns[jmin] + 1; i++)
			printf("jmin %d (%g %g)\n", i, c->x[jmin][i], c->y[jmin][i]);
	}

	memcpy(&(c->x[imin][c->ns[imin] + 1]), &(c->x[jmin][0]), (c->ns[jmin]+1)*sizeof(double));
	memcpy(&(c->y[imin][c->ns[imin] + 1]), &(c->y[jmin][0]), (c->ns[jmin]+1)*sizeof(double));
	c->ns[imin] = c->ns[imin] + c->ns[jmin] + 1;

	if (debug) {
		printf("JOINED IMIN\n");
		for (i = 0; i < c->ns[imin] + 1; i++)
			printf("imin %d (%g %g)\n", i, c->x[imin][i], c->y[imin][i]);
	}

	/*
	 *remove jmin
	 */
	Free(c->x[jmin]);
	Free(c->y[jmin]);
	for (j = jmin; j < c->nc - 1; j++) {
		c->x[j] = c->x[j + 1];
		c->y[j] = c->y[j + 1];
		c->ns[j] = c->ns[j + 1];
		c->cyclic[j] = c->cyclic[j + 1];
		c->length[j] = c->length[j + 1];
	}

	c->length[imin] = 0.0;
	for (i = 0; i < c->ns[imin]; i++) {
		c->length[imin] += sqrt(SQR(c->x[imin][i] - c->x[imin][i + 1]) + SQR(c->y[imin][i] - c->y[imin][i + 1]));
	}
	c->cyclic[imin] = (c->x[imin][0] == c->x[imin][c->ns[imin]]) && (c->y[imin][0] == c->y[imin][c->ns[imin]]);
	c->nc--;

	Free(dists);
	Free(idir);
	Free(jdir);

	return re_join_contourLines(c);
}


double re_point_on_countour(inla_contour_tp * c, double skew, double kurt)
{
#define WRAPIT(_i) (((_i) + c->ns[ic]+1) % (c->ns[ic]+1))
#define INBETWEEN(_x, _x0, _x1) (( (_x) >= (_x0) && (_x) <= (_x1)) || ((_x) >= (_x1) && (_x) <= (_x0)) )

	/*
	 * locate the corresponding 'point' on the contour; counted as the length around the contour from a well defined starting point. 
	 */

	int i, j, start_idx, ic, direction;
	double len;

	GMRFLib_max_value(c->length, c->nc, &ic);
	if (c->cyclic[ic]) {
		GMRFLib_min_value(c->y[ic], c->ns[ic] + 1, &start_idx);
	} else {
		start_idx = 0;
		if (!(c->x[ic][0] < c->x[ic][c->ns[ic] + 1])) {
			/*
			 * simply reverse the contour
			 */
			REV(c->x[ic], c->ns[ic] + 1);
			REV(c->y[ic], c->ns[ic] + 1);
			if (0) {
				FILE *fp = fopen("contour.dat", "w");
				inla_print_contourLines(fp, c);
				fprintf(stderr, "\n\n *** trouble; file written....\n\n");
				assert(0 == 1);
			}
		}
	}

	// P(start_idx);
	direction = (c->x[ic][WRAPIT(start_idx + 1)] < c->x[ic][start_idx] ? 1 : -1);

	len = 0.0;
	for (i = start_idx, j = 0; j < c->ns[ic] + 1; i = WRAPIT(i + direction), j++) {
		if (INBETWEEN(skew, c->x[ic][i], c->x[ic][WRAPIT(i + direction)]) && INBETWEEN(kurt, c->y[ic][i], c->y[ic][WRAPIT(i + direction)])) {
			len += sqrt(SQR(c->x[ic][i] - skew) + SQR(c->y[ic][i] - kurt));
			break;
		} else {
			len += sqrt(SQR(c->x[ic][i] - c->x[ic][WRAPIT(i + direction)]) + SQR(c->y[ic][i] - c->y[ic][WRAPIT(i + direction)]));
		}
	}
#undef WRAPIT
#undef INBETWEEN
	return len;
}
double re_sas_evaluate_log_prior(double skew, double kurt)
{
	if (GMRFLib_FALSE) {
		return -0.5 * 1 * SQR(skew) - 0.5 * 1 * SQR(kurt - 3);
	} else {
		double level;
		inla_contour_tp *c;
		double ldens_uniform, ldens_dist, length = 0, point;

		level = re_find_in_sas_prior_table(skew, kurt);	/* level for the contour */
		if (ISNAN(level)) {
			return 0.0;
		}
		c = contourLines(sas_prior_table->x, sas_prior_table->nx, sas_prior_table->y, sas_prior_table->ny, sas_prior_table->z, level);
		assert(c->nc);
		if (c->nc > 1) {
#pragma omp critical
			{
				static int count = 0;
				char *filename;
				FILE *fp;

				fprintf(stderr, "log_prior: nc=%1d, skew=%g kurt=%g\n", c->nc, skew, kurt);
				GMRFLib_sprintf(&filename, "c/contour-%1d.dat", count);
				fp = fopen(filename, "w");
				inla_print_contourLines(fp, c);
				fclose(fp);

				re_join_contourLines(c);
				GMRFLib_sprintf(&filename, "c/contour-NEW-%1d.dat", count);
				fp = fopen(filename, "w");
				inla_print_contourLines(fp, c);
				fclose(fp);

				Free(filename);
				count++;
			}
		}
		length = GMRFLib_max_value(c->length, c->nc, NULL);
		ldens_uniform = log(1.0 / length);
		point = re_point_on_countour(c, skew, kurt);
		inla_free_contourLines(c);

		double lambda = 20;
		ldens_dist = log(lambda) - lambda * level;

#define NEW(dskew, dkurt)						\
		if (1)							\
		{							\
			level = re_find_in_sas_prior_table(skew + dskew, kurt+dkurt); \
			if (ISNAN(level)) return 0.0;			\
			c = contourLines(sas_prior_table->x, sas_prior_table->nx, \
					 sas_prior_table->y, sas_prior_table->ny, \
					 sas_prior_table->z, level);	\
			assert(c->nc);					\
			int ic_max;					\
			GMRFLib_max_value(c->length, c->nc, &ic_max);	\
			length = c->length[ic_max];			\
			point = re_point_on_countour(c, skew + dskew, kurt+dkurt); \
			inla_free_contourLines(c);			\
		}

		double d_level_d_skew, d_level_d_kurt, d_point_d_skew, d_point_d_kurt;
		double lev[2], poi[2];
		double dskew, dkurt;

		if (skew <= 0.0) {
			dskew = -0.01;
			dkurt = 0.01;
		} else {
			dskew = 0.01;
			dkurt = 0.01;
		}

		lev[0] = level;
		poi[0] = point;

		NEW(dskew, 0);
		lev[1] = level;
		poi[1] = point;

		d_level_d_skew = (lev[1] - lev[0]) / 2.0 / dskew;
		d_point_d_skew = (poi[1] - poi[0]) / 2.0 / dskew;

		NEW(0, dkurt);
		lev[1] = level;
		poi[1] = point;

		d_level_d_kurt = (lev[1] - lev[0]) / 2.0 / dkurt;
		d_point_d_kurt = (poi[1] - poi[0]) / 2.0 / dkurt;

		double Jacobian = ABS(d_level_d_skew * d_point_d_kurt - d_point_d_skew * d_level_d_kurt);
#undef NEW
		return ldens_uniform + ldens_dist + log(Jacobian);
		// FIXME1("ONLY JACOBIAN");
		// return log(Jacobian);
	}
}
double re_find_in_sas_prior_table(double skew, double kurt)
{
	/*
	 * return the level based on 'skew' and 'kurt' 
	 */

	int ix = re_find_in_table_general(skew, sas_prior_table->x, sas_prior_table->nx);
	int iy = re_find_in_table_general(kurt, sas_prior_table->y, sas_prior_table->ny);
	// printf("skew %g kurt %g gives ix %d iy %d\n", skew, kurt, ix, iy);
	assert(LEGAL(ix, sas_prior_table->nx));
	assert(LEGAL(iy, sas_prior_table->ny));

	/*
	 * so the notation is like the one in https://en.wikipedia.org/wiki/Bilinear_interpolation 
	 */
#define LEVEL(_i, _j) sas_prior_table->z[(_i) + sas_prior_table->nx * (_j)]
	double x1, x2, y1, y2, x, y, fQ11, fQ12, fQ21, fQ22, level;

	x1 = sas_prior_table->x[ix];
	x2 = sas_prior_table->x[ix + 1];
	y1 = sas_prior_table->y[iy];
	y2 = sas_prior_table->y[iy + 1];
	fQ11 = LEVEL(ix, iy);
	fQ21 = LEVEL(ix + 1, iy);
	fQ12 = LEVEL(ix, iy + 1);
	fQ22 = LEVEL(ix + 1, iy + 1);
	x = skew;
	y = kurt;

	level = 1.0 / (DMAX(DBL_EPSILON, x2 - x1) * DMAX(DBL_EPSILON, y2 - y1))
	    * (fQ11 * (x2 - x) * (y2 - y) + fQ21 * (x - x1) * (y2 - y) + fQ12 * (x2 - x) * (y - y1) + fQ22 * (x - x1) * (y - y1));
#undef LEVEL
	return level;
}

int re_find_in_table_general(double value, double *x, int nx)
{
	// x[SOL] < value but the largest one
	int debug = 0;

	if (value < x[0]) {
		return -1;
	}
	if (value > x[nx - 1]) {
		return nx;
	}

	int low = 0, high = nx - 1, new_guess;

	while (GMRFLib_TRUE) {
		if (debug) {
			printf("locate %g low %d high %d\n", value, low, high);
		}
		if (high - low == 1) {
			if (debug) {
				printf("locate %g found %d\n", value, low);
			}
			return (low);
		}
		new_guess = (int) ((low + high) / 2);
		if (x[new_guess] <= value) {
			low = new_guess;
		} else {
			high = new_guess;
		}
	}
}

int re_make_sas_prior_table(void)
{
	if (sas_prior_table) {
		return GMRFLib_SUCCESS;
	}
//#pragma omp critical
	{
		if (!sas_prior_table) {
			printf("MAKE TABLE...\n");

			re_sas_prior_tp *s = Calloc(1, re_sas_prior_tp);
#define NX 500
#define NY 500
#define NZ (NX*NY)
#define SKEW 3
#define KURT_LOW 2.1
#define KURT_HIGH 10
			s->nx = NX;
			s->ny = NY;
			s->nz = NZ;

			s->x = Calloc(s->nx, double);
			s->y = Calloc(s->ny, double);
			s->z = Calloc(s->nz, double);

			int i, j;
			double dd;

			dd = 2 * SKEW / (s->nx - 1.0);
			for (i = 0; i < s->nx; i++) {
				s->x[i] = -SKEW + i * dd;
			}

			dd = (KURT_HIGH - KURT_LOW) / (s->ny - 1.0);
			for (j = 0; j < s->ny; j++) {
				s->y[j] = KURT_LOW + j * dd;
			}

#pragma omp parallel for private(j, i)
			for (j = 0; j < s->ny; j++) {
				for (i = 0; i < s->nx; i++) {
					int k = i + j * s->nx;
					if (re_valid_skew_kurt(NULL, s->x[i], s->y[j])) {
						s->z[k] = re_intrinsic_discrepancy_distance_map(re_intrinsic_discrepancy_distance(s->x[i], s->y[j]));
						// printf("valid skew %g kurt %g dist %g\n", s->x[i], s->y[j], s->z[k]);
						assert(!ISNAN(s->z[k]));
					} else {
						// printf("not valid skew %g kurt %g\n", s->x[i], s->y[j]);
						s->z[k] = NAN;
					}
				}
			}
			printf("MAKE TABLE... DONE\n");

			FILE *fp = fopen("table.dat", "w");
			fprintf(fp, "%g\n", (double) s->nx);
			fprintf(fp, "%g\n", (double) s->ny);
			fprintf(fp, "%g\n", (double) s->nz);
			for (i = 0; i < s->nx; i++)
				fprintf(fp, "%.12g\n", s->x[i]);
			for (i = 0; i < s->ny; i++)
				fprintf(fp, "%.12g\n", s->y[i]);
			for (i = 0; i < s->nz; i++)
				fprintf(fp, "%.12g\n", s->z[i]);
			fclose(fp);

			sas_prior_table = s;
		}
	}
#undef NX
#undef NY
#undef NZ
#undef SKEW
#undef KURT_LOW
#undef KURT_HIGH
	return GMRFLib_SUCCESS;
}
int re_dsas_intern(double *logdens, double *x, int n, double mu, double sigma, double delta, double epsilon)
{
	int i;
	double z, c, r, tmp;

	for (i = 0; i < n; i++) {
		z = (x[i] - mu) / sigma;
		tmp = delta * asinh(z) - epsilon;
		c = cosh(tmp);
		r = sinh(tmp);
		logdens[i] = -log(sigma) + log(delta) - log(2 * M_PI) / 2.0 - log(1.0 + SQR(z)) / 2.0 + log(c) - SQR(r) / 2.0;
		if (ISNAN(logdens[i])) {
			P(log(sigma));
			P(delta);
			P(log(delta));
			P(log(1.0 + SQR(z)));
			P(log(c));
			P(SQR(r));
			abort();
		}
	}

	return GMRFLib_SUCCESS;
}

int re_dsas(double *logdens, double *x, int n, double skew, double kurt)
{
	/*
	 * mean zero and precision one
	 */
	re_shash_param_tp param;
	double zero = 0.0, one = 1.0;

	re_shash_fit_parameters(&param, &zero, &one, &skew, &kurt);
	re_dsas_intern(logdens, x, n, param.mu, param.stdev, param.delta, param.epsilon);

	return GMRFLib_SUCCESS;
}
int re_dnorm(double *logdens, double *x, int n)
{

	int i;
	for (i = 0; i < n; i++) {
		logdens[i] = -0.91893853320467274178032973640560 - 0.5 * SQR(x[i]);
	}

	return GMRFLib_SUCCESS;
}
double re_intrinsic_discrepancy_distance_map(double distance)
{
	/*
	 * map distance to the eqivalent distance between N(0,1) and N(mu,1)
	 */
	return (sqrt(2.0 * distance));
}
double re_intrinsic_discrepancy_distance(double skew, double kurt)
{
#define N (2*10000+3)
#define LIM (10.0)

	double w[2] = { 4.0, 2.0 }, dx = (2.0 * LIM) / (N - 1.0), ldens_sas[N], integral[2] = {
	0.0, 0.0};
	int i, k;

	double x[N], ldens_norm[N];			       /* these are constants */
	for (i = 0; i < N; i++) {
		x[i] = -LIM + i * dx;
	}
	re_dnorm(ldens_norm, x, N);

#define F0(_idx) exp(ldens_norm[_idx]) * (ldens_norm[_idx] - ldens_sas[_idx])
#define F1(_idx) exp(ldens_sas[_idx]) * (ldens_sas[_idx] - ldens_norm[_idx])

	re_dsas(ldens_sas, x, N, skew, kurt);
	integral[0] = F0(0) + F0(N - 1);
	if (ISNAN(integral[0])) {
		P(F0(0));
		P(F0(N - 1));
		abort();
	}
	integral[1] = F1(0) + F1(N - 1);

	for (i = 1, k = 0; i < N - 1; i++, k = (k + 1) % 2) {
		integral[0] += F0(i) * w[k];
		integral[1] += F1(i) * w[k];
	}
	integral[0] *= dx / 3.0;
	integral[1] *= dx / 3.0;

#undef N
#undef LIM
#undef F0
#undef F1
	return (DMIN(integral[0], integral[1]));
}

int re_init()
{
	re_shash_param_tp param;
	double skew = 0, kurt = 3.0, mean = 0, s = 1;
	re_shash_fit_parameters(&param, &mean, &s, &skew, &kurt);

	re_make_sas_prior_table();

	return GMRFLib_SUCCESS;
}
