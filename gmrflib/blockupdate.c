#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>

#include "GMRFLib/GMRFLib.h"

FORCEINLINE double GMRFLib_prod_diff(double a, double b, double c, double d)
{
	// return a*b-c*d , see https://pharr.org/matt/blog/2019/11/03/difference-of-floats 
	double cd = c * d;
	return fma(a, b, -cd) + fma(-c, d, cd);
}

int GMRFLib_default_blockupdate_param(GMRFLib_blockupdate_param_tp **blockupdate_par)
{
	GMRFLib_ASSERT(blockupdate_par, GMRFLib_EINVARG);

	*blockupdate_par = Calloc(1, GMRFLib_blockupdate_param_tp);
	(*blockupdate_par)->modeoption = GMRFLib_MODEOPTION_MODE;
	(*blockupdate_par)->fp = NULL;
	(*blockupdate_par)->step_len = GSL_ROOT4_DBL_EPSILON;
	(*blockupdate_par)->stencil = 5;

	return GMRFLib_SUCCESS;
}


int GMRFLib_2order_taylor(int thread_id, int *lcache_idx, double *a, double *b, double *c, double *dd, double d, double x0, int idx,
			  double *x_vec, GMRFLib_logl_tp *loglFunc, void *loglFunc_arg, double *step_len, int *stencil)
{
	/*
	 * compute a,b,c in the taylor expansion around x0 of d*loglFunc(x0,...)
	 * 
	 * a + b*(x-x0) + 0.5*c*(x-x0)^2 + 1/6*dd*(x-x0)^3
	 * 
	 */
	double f0 = 0.0, df = 0.0, ddf = 0.0, dddf = 0.0;

	if (ISZERO(d)) {
		f0 = df = ddf = 0.0;
	} else {
		GMRFLib_2order_approx_core(thread_id, lcache_idx, &f0, &df, &ddf, (dd ? &dddf : NULL), x0, idx, x_vec, loglFunc, loglFunc_arg,
					   step_len, stencil);
	}

	*a = f0;
	*b = df;
	*c = ddf;
	if (dd)
		*dd = dddf;

	if (d != 1.0) {
		*a *= d;
		*b *= d;
		*c *= d;
		if (dd) {
			*dd *= d;
		}
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_2order_approx(int thread_id, int *lcache_idx, double *a, double *b, double *c, double *dd, double d, double x0, int idx,
			  double *x_vec, GMRFLib_logl_tp *loglFunc, void *loglFunc_arg, double *step_len, int *stencil, double *cmin)
{
	/*
	 * compute a,b,c in the taylor expansion around x0 of d*loglFunc(x0,...)
	 * 
	 * a + b*x - 0.5*c*x^2 + 1/6*dd*x^3
	 *
	 * where cmin is the minimum value of c.
	 */

	/*
	 * > A:=collect(expand(a + b * (x-x0) + 1/2 \                                      
	 * > * c * (x-x0)^2 + 1/6 * d * (x-x0)^3), [x,x^2, x^3]);
	 *          3                    /            2    \             2              3
	 *       d x   /      d x0\  2   |        d x0     |         c x0           d x0
	 * A := ---- + |c/2 - ----| x  + |-c x0 + ----- + b| x + a + ----- - b x0 - -----
	 *             \       2  /      \          2      /           2              6
	 *
	 * > coeff(A,x);                                                                   
	 *             2
	 *         d x0
	 * -c x0 + ----- + b
	 *           2
	 * 
	 * > coeff(A,x^2);
	 *       d x0
	 * c/2 - ----
	 *        2
	 * 
	 * > coeff(A,x^3);
	 * d/6
	 * 
	 */

#define INVALID(x_) (ISNAN(x_) || ISINF(x_))

	double f0 = 0.0, df = 0.0, ddf = 0.0, dddf = 0.0;
	int rescue = 0;
	static int give_warning_c = 0;
	static int give_warning_idx = -1;		       /* this is set at first call */

	if (give_warning_idx < 0) {
#pragma omp critical (Name_53f42442f89dc6478eaee39aa0766bbff846950c)
		if (give_warning_idx < 0) {
			give_warning_idx = idx;
		}
	}

	if (idx == give_warning_idx && give_warning_c > 1) {
		fprintf(stderr, "[%1d] *** WARNING *** GMRFLib_2order_approx: reset counter for %1d NAN/INF values in logl\n",
			omp_get_thread_num(), give_warning_c);
#pragma omp critical (Name_61a18063454b0e56bccffa14dda9ace39df612f8)
		give_warning_c = 0;
	}

	GMRFLib_2order_approx_core(thread_id, lcache_idx, &f0, &df, &ddf, (dd ? &dddf : NULL), x0, idx, x_vec, loglFunc, loglFunc_arg, step_len,
				   stencil);
	if (INVALID(ddf)) {
		if (give_warning_c == 0) {
			fprintf(stderr, " *** WARNING *** GMRFLib_2order_approx: rescue NAN/INF values in logl for idx=%1d\n", idx);
		}
#pragma omp atomic
		give_warning_c++;

		f0 = df = 0.0;
		ddf = -1.0;				       /* we try with this */
		if (dd) {
			dddf = 0.0;
		}
		rescue = 1;
	} else {
		if (cmin) {
			ddf = DMIN(-(*cmin), ddf);
		}
	}

	if (rescue) {
		*a = 0.0;
		*b = 0.0;
		*c = -d * ddf;
		if (dd) {
			*dd = 0.0;
		}
	} else {
		*a = f0 + x0 * (-df + 0.5 * x0 * (ddf + 0.3333333333333333333 * dddf * x0));
		*b = df + x0 * (-ddf + 0.5 * dddf * x0);
		*c = -ddf + dddf * x0;
		if (dd) {
			*dd = dddf;
		}

		if (d != 1.0) {
			*a *= d;
			*b *= d;
			*c *= d;
			if (dd) {
				*dd *= d;
			}
		}
	}

#undef INVALID
	return GMRFLib_SUCCESS;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((target_clones(INLA_CLONE_TARGETS "default")))
int GMRFLib_2order_approx_core(int thread_id, int *lcache_idx, double *a, double *b, double *c, double *dd, double x0, int idx,
			       double *x_vec, GMRFLib_logl_tp *loglFunc, void *loglFunc_arg, double *step_len, int *stencil)
{
	// default step-size is determined using test=151
	double step, df = 0.0, ddf = 0.0, dddf = 0.0, xx[9], f[9], f0 = 0.0, x00;
	int stenc = (stencil ? *stencil : 5);

	if (step_len && *step_len < 0.0) {
		/*
		 * for internal use only! 
		 */
		step = -(*step_len);

		xx[0] = x0 - 2 * step;
		xx[1] = x0 - step;
		xx[2] = x0;
		xx[3] = x0 + step;
		xx[4] = x0 + 2 * step;

		loglFunc(thread_id, lcache_idx, f, xx, 5, idx, x_vec, NULL, loglFunc_arg);

		f0 = f[2];
		df = (1.0 / 12.0 * f[4] - 2.0 / 3.0 * f[3] + 0.0 * f[2] + 2.0 / 3.0 * f[1] - 1.0 / 12.0 * f[0]) / step;
		ddf = (-1.0 / 12.0 * f[4] - 4.0 / 3.0 * f[3] - 5.0 / 2.0 * f[2] + 4.0 / 3.0 * f[1] - 1.0 / 12.0 * f[0]) / SQR(step);
		dddf = (-1.0 / 2.0 * f[4] + 1.0 * f[3] + 0.0 * f[2] - 1.0 * f[1] + 1.0 / 2.0 * f[0]) / POW3(step);
	} else {

		// this is the plain code
		// df=GMRFLib_ddot(n, wf, f);
		// ddf=GMRFLib_ddot(n, wff, f);

		switch (stenc) {
		case 3:
		{
			// special implementation: ONLY used for initial values
			step = 1.0e-4;
			const int n = 3;
			xx[0] = x0 - step;
			xx[1] = x0;
			xx[2] = x0 + step;

			loglFunc(thread_id, lcache_idx, f, xx, n, idx, x_vec, NULL, loglFunc_arg);

			f0 = f[1];
			df = 0.5 * (-f[0] + f[2]);
			ddf = f[0] - 2.0 * f[1] + f[2];
		}
			break;

		case 5:
		{
			if (unlikely(!step_len || ISZERO(*step_len))) {
				double ref = GSL_DBL_EPSILON / 2.220446049e-16;
				step = ref * 5.0E-4;
			} else {
				step = *step_len;
			}

			const int n = 5, nn = 2, wlength = 8;
			static const double wf[24] = {
				1.0 / 12.0, 
				- 2.0 / 3.0, 
				0,
				2.0 / 3.0, 
				-1.0 / 12.0, 
				0,
				0,
				0,

				- 1.0 / 12.0, 
				4.0 / 3.0, 
				-2.5,
				4.0 / 3.0, 
				- 1.0 / 12.0, 
				0,
				0,
				0,

				-0.5,
				1,
				0,
				-1,
				0.5,
				0,
				0,
				0
			};

			x00 = x0 - nn * step;
			for (int i = 0; i < n; i++) {
				xx[i] = x00 + i * step;
			}

			loglFunc(thread_id, lcache_idx, f, xx, n, idx, x_vec, NULL, loglFunc_arg);
			f0 = f[nn];

			double *wff = (double *) wf + wlength;
			double *wfff = (double *) wf + 2 * wlength;
			double *f_ref = f + nn;
			double *wf_ref = (double *) wf + nn;
			double *wff_ref = wff + nn;
#if 1
			if (!dd) {
				ddf = f_ref[0] * wff_ref [0];
				for(int i = 1; i <= nn ; i++) {
					df += wf_ref[i] * (f_ref[i] - f_ref[-i]);
					ddf += wff_ref[i] * (f_ref[i] + f_ref[-i]);
				}
			} else {
				double *wfff_ref = wfff + nn;
				ddf = f_ref[0] * wff_ref [0];
				for(int i = 1; i <= nn ; i++) {
					double dif = f_ref[i] - f_ref[-i];
					df += wf_ref[i] * dif;
					ddf += wff_ref[i] * (f_ref[i] + f_ref[-i]);
					dddf += wfff_ref[i] * dif;
				}
			}

#else
			df = GMRFLib_prod_diff(wf_ref[1], f_ref[1] - f_ref[-1], -wf_ref[2], f_ref[2] - f_ref[-2]);
			ddf = GMRFLib_prod_diff(wff_ref[1], f_ref[-1] + f_ref[1], -wff_ref[2], f_ref[-2] + f_ref[2]);
			ddf = fma(wff_ref[0], f_ref[0], ddf);
			if (dd) {
				double *wfff_ref = wfff + nn;
				dddf = GMRFLib_prod_diff(wfff_ref[1], f_ref[1] - f_ref[-1], -wfff_ref[2], f_ref[2] - f_ref[-2]);
			}
#endif
		}
			break;

		case 7:
		{
			if (!step_len || ISZERO(*step_len)) {
				double ref = GSL_DBL_EPSILON / 2.220446049e-16;
				step = ref * 100.0E-4;
			} else {
				step = *step_len;
			}

			const int n = 7, nn = 3, wlength = 8;
			static const double wf[24] = {
				-0.01666666666666667,
				0.15,
				-0.75,
				0,
				0.75,
				-0.15,
				0.01666666666666667,
				0,

				0.01111111111111111,
				-0.15,
				1.5,
				-2.722222222222222,
				1.5,
				-0.15,
				0.01111111111111111,
				0,

				0.125,
				-1,
				1.625,
				0,
				-1.625,
				1,
				-0.125,
				0
			};

			x00 = x0 - nn * step;
			for (int i = 0; i < n; i++) {
				xx[i] = x00 + i * step;
			}

			loglFunc(thread_id, lcache_idx, f, xx, n, idx, x_vec, NULL, loglFunc_arg);
			f0 = f[nn];

			double *wff = (double *) wf + wlength;
			double *wfff = (double *) wf + 2 * wlength;
			double *f_ref = f + nn;
			double *wf_ref = (double *) wf + nn;
			double *wff_ref = wff + nn;

			// we do not need to initialized df and dddf, as wf_ref[0]=0 and wfff_ref[0]=0
#if 1
			ddf = f_ref[0] * wff_ref [0];
			if (!dd) {
				for(int i = 1; i <= nn ; i++) {
					df += wf_ref[i] * (f_ref[i] - f_ref[-i]);
					ddf += wff_ref[i] * (f_ref[i] + f_ref[-i]);
				}
			} else {
				double *wfff_ref = wfff + nn;
				for(int i = 1; i <= nn ; i++) {
					double dif = f_ref[i] - f_ref[-i];
					df += wf_ref[i] * dif;
					ddf += wff_ref[i] * (f_ref[i] + f_ref[-i]);
					dddf += wfff_ref[i] * dif;
				}
			}
#else
			df = GMRFLib_prod_diff(wf_ref[1], f_ref[1] - f_ref[-1], -wf_ref[2], f_ref[2] - f_ref[-2]);
			df = fma(wf_ref[3], f_ref[3] - f_ref[-3], df);
			ddf = GMRFLib_prod_diff(wff_ref[0], f_ref[0], -wff_ref[1], f_ref[1] + f_ref[-1]) +
				GMRFLib_prod_diff(wff_ref[2], f_ref[2] + f_ref[-2], -wff_ref[3], f_ref[3] + f_ref[-3]);
			if (dd) {
				double *wfff_ref = wfff + nn;
				dddf = GMRFLib_prod_diff(wfff_ref[1], f_ref[1] - f_ref[-1], -wfff_ref[2], f_ref[2] - f_ref[-2]);
				dddf = fma(wfff_ref[3], f_ref[3] - f_ref[-3], dddf);
			}
#endif
		}
			break;

		default:
			assert(0 == 1);
		}
	}

	double istep = 1.0 / step;
	df *= istep;
	ddf *= SQR(istep);
	*a = f0;
	*b = df;
	*c = ddf;

	if (dd) {
		*dd = dddf * POW3(istep);
	}

	return GMRFLib_SUCCESS;
}
#pragma GCC diagnostic pop
