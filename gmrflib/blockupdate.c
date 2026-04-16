#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>

#include "GMRFLib/GMRFLib.h"

#pragma omp declare simd
static double GMRFLib_prod_diff(double a, double b, double c, double d)
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
		fprintf(stderr, " *** WARNING *** GMRFLib_2order_approx: reset counter for %1d NAN/INF values in logl\n", give_warning_c);
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
	// default step-size is determined using test=151. stencil=9 does not bring much...

	double step, df = 0.0, ddf = 0.0, dddf = 0.0, xx[9], f[9], f0 = 0.0, x00;
	int stenc = (stencil ? *stencil : 5);
	int numa = GMRFLib_numa_get_node();

	typedef struct {
		double **wf;
	} wf_tp;

	static wf_tp **lwork = NULL;
	if (!lwork) {
#pragma omp critical (Name_009f5f31299b4b554b667873ad6c4c874bfc9a77)
		if (!lwork) {
			wf_tp **tmp = Calloc(GMRFLib_numa_nodes(), wf_tp *);
			lwork = tmp;
		}
	}

	if (!lwork[numa]) {
#pragma omp critical (Name_b53c77704653d4b6a42cc3c6c8221441fac46a73)
		if (!lwork[numa]) {
			wf_tp *w = Calloc(1, wf_tp);
			w->wf = Calloc(10, double *);	       /* Must initialize to 0 */
			lwork[numa] = w;
		}
	}

	wf_tp *w = lwork[numa];

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
			int n = 3;
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
				static double ref = GSL_DBL_EPSILON / 2.220446049e-16;
				step = ref * 5.0e-4;
			} else {
				step = *step_len;
			}

			int n = 5, nn = 2, wlength = 8;

			if (unlikely(!(w->wf[stenc]))) {
#pragma omp critical (Name_4eb4719ffe22f0af964510f0aec612baccccbb0d)
				if (!(w->wf[stenc])) {
					int len_offset = GMRFLib_memory_alignment / sizeof(double);
					int len = 3 * wlength + len_offset;
					double *ww = Malloc(len, double);
					GMRFLib_dfill(len, 0.0, ww);
					GMRFLib_ENSURE_NUMA_PTR(ww, len, double);

					if (1 || GMRFLib_memory_alignment_enabled) {
						// ensure ww is aligned. we might change the ptr so we cannot free
						int ok = 0;
						for(int k = 0; k < len; k++) {
							if (GMRFLib_is_aligned(ww + k)) {
								ww += k;
								ok = 1;
								break;
							}
						}
						if (!ok) FIXME("Memory alignment failed");
					}

					ww[0] = 1.0 / 12.0;
					ww[1] = -2.0 / 3.0;
					ww[3] = 2.0 / 3.0;
					ww[4] = - 1.0 / 12.0;
					ww[8] = - 1.0 / 12.0;
					ww[9] = 4.0 / 3.0; 
					ww[10] = -2.5;
					ww[11] = 4.0 / 3.0; 
					ww[12] = - 1.0 / 12; 
					ww[16] = -0.5;
					ww[17] = 1.0;
					ww[19] = -1.0;
					ww[20] = 0.5;
					w->wf[stenc] = ww;
				}
			}

			double *wf = w->wf[stenc];
			double *wff = wf + wlength;
			double *wfff = wf + 2 * wlength;

			x00 = x0 - nn * step;
#pragma omp simd
			for (int i = 0; i < n; i++) {
				xx[i] = x00 + i * step;
			}

			loglFunc(thread_id, lcache_idx, f, xx, n, idx, x_vec, NULL, loglFunc_arg);
			f0 = f[nn];

			int iref = n / 2L;
			double *f_ref = f + iref;
			double *wf_ref = wf + iref;
			double *wff_ref = wff + iref;

			df = GMRFLib_prod_diff(wf_ref[1], f_ref[1] - f_ref[-1], -wf_ref[2], f_ref[2] - f_ref[-2]);
			ddf = GMRFLib_prod_diff(wff_ref[1], f_ref[-1] + f_ref[1], -wff_ref[2], f_ref[-2] + f_ref[2]);
			ddf = fma(wff_ref[0], f_ref[0], ddf);

			if (dd) {
				double *wfff_ref = wfff + iref;
				dddf = GMRFLib_prod_diff(wfff_ref[1], f_ref[1] - f_ref[-1], -wfff_ref[2], f_ref[2] - f_ref[-2]);
			}
		}
			break;

		case 7:
		{
			if (!step_len || ISZERO(*step_len)) {
				static double ref = GSL_DBL_EPSILON / 2.220446049e-16;
				step = ref * 100.0e-4;
			} else {
				step = *step_len;
			}

			int n = 7, nn = 3, wlength = 8;

			if (!(w->wf[stenc])) {
#pragma omp critical (Name_0eed179363c2b9a7edfda8a212fc6f63e8ec9741)
				if (!(w->wf[stenc])) {
					int len_offset = GMRFLib_memory_alignment / sizeof(double);
					int len = 3 * wlength + len_offset;
					double *ww = Malloc(len, double);
					GMRFLib_dfill(len, 0.0, ww);
					GMRFLib_ENSURE_NUMA_PTR(ww, len, double);

					if (1 || GMRFLib_memory_alignment_enabled) {
						// ensure ww is aligned. we might change the ptr so we cannot free
						int ok = 0;
						for(int k = 0; k < len_offset; k++) {
							if (GMRFLib_is_aligned(ww + k)) {
								ww += k;
								ok = 1;
								break;
							}
						}
						if (!ok) FIXME("Memory alignment failed");
					}

					ww[0] = - 1.0 / 60.0;
					ww[1] = 0.15;
					ww[2] = -0.75;
					ww[4] = 0.75;
					ww[5] = -0.15;
					ww[6] = 1.0 / 60.0; 
					ww[8] = 1.0 / 90.0;
					ww[9] = -0.15;
					ww[10] = 1.5;
					ww[11] = - 49.0 / 18.0;
					ww[12] = 1.5;
					ww[13] = -0.15;
					ww[14] = 1.0 / 90.0;
					ww[16] = 0.125;
					ww[17] = -1.0;
					ww[18] = 1.625;
					ww[20] = -1.625;
					ww[21] = 1.0;
					ww[22] = -0.125;
					w->wf[stenc] = ww;
				}
			}

			double *wf = w->wf[stenc];
			double *wff = wf + wlength;
			double *wfff = wf + 2 * wlength;

			x00 = x0 - nn * step;
#pragma omp simd
			for (int i = 0; i < n; i++) {
				xx[i] = x00 + i * step;
			}

			loglFunc(thread_id, lcache_idx, f, xx, n, idx, x_vec, NULL, loglFunc_arg);
			f0 = f[nn];

			int iref = n / 2L;
			double *f_ref = f + iref;
			double *wf_ref = wf + iref;
			double *wff_ref = wff + iref;

			df = GMRFLib_prod_diff(wf_ref[1], f_ref[1] - f_ref[-1], -wf_ref[2], f_ref[2] - f_ref[-2]);
			df = fma(wf_ref[3], f_ref[3] - f_ref[-3], df);

			ddf = GMRFLib_prod_diff(wff_ref[0], f_ref[0], -wff_ref[1], f_ref[1] + f_ref[-1]) +
			    GMRFLib_prod_diff(wff_ref[2], f_ref[2] + f_ref[-2], -wff_ref[3], f_ref[3] + f_ref[-3]);

			if (dd) {
				double *wfff_ref = wfff + iref;
				dddf = GMRFLib_prod_diff(wfff_ref[1], f_ref[1] - f_ref[-1], -wfff_ref[2], f_ref[2] - f_ref[-2]);
				dddf = fma(wfff_ref[3], f_ref[3] - f_ref[-3], dddf);
			}
		}
			break;

		default:
			assert(0 == 1);
		}
	}

	double istep = 1.0 / step;
	df *= istep;
	ddf *= SQR(istep);
	if (dd) {
		dddf *= POW3(istep);
	}

	*a = f0;
	*b = df;
	*c = ddf;
	if (dd) {
		*dd = dddf;
	}

	return GMRFLib_SUCCESS;
}
#pragma GCC diagnostic pop
