
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

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

#include "re.h"
#include "bessel.h"

/* 
   the version from R seems can handle more extreme cases without throwing an error..
 */
#define BESSEL_KNU(_alpha, _x) bessel_Knu(_alpha, _x)
//#define BESSEL_KNU(_alpha, _x) gsl_sf_bessel_Knu(_alpha, _x)


#define Pq(_arg) (exp(0.25)/sqrt(8.0*M_PI) * (BESSEL_KNU(((_arg)+1.0)/2.0, 0.25) + BESSEL_KNU(ABS(((_arg)-1.0)/2.0), 0.25)))
#define M1(_epsilon, _delta) (sinh((_epsilon)/(_delta))*Pq(1.0/(_delta)))
#define M2(_epsilon, _delta) (0.5*(cosh(2.0*(_epsilon)/(_delta))*Pq(2.0/(_delta)) -1.0))
#define M3(_epsilon, _delta) (0.25*(sinh(3.0*(_epsilon)/(_delta))*Pq(3.0/(_delta)) - \
				    3.0*sinh((_epsilon)/(_delta))*Pq(1.0/(_delta))))
#define M4(_epsilon, _delta) (0.125*(cosh(4.0*(_epsilon)/(_delta))*Pq(4.0/(_delta)) - \
				     4.0*cosh(2.0*(_epsilon)/(_delta))*Pq(2.0/(_delta)) + 3.0))

#define XMATCH(x0,x1) (fabs((x0)-(x1)) == 0.0)
#define YMATCH(y0,y1) (fabs((y0)-(y1)) == 0.0)

#define KURT_LIMIT(s) (2.15 + SQR((s)/0.8))

#define SKEW_MIN (sas_prior_table->skew[0])
#define SKEW_MAX (sas_prior_table->skew[sas_prior_table->nx-1])
#define KURT_MIN (sas_prior_table->kurt[0])
#define KURT_MAX (sas_prior_table->kurt[sas_prior_table->ny-1])
#define SCALEWITH(x) DMAX(1.0/1.5, DMIN(1.5, (x)))

#define REV(x, n) \
	if (1){						\
		double *_tmp = Calloc(n, double);	\
		memcpy(_tmp, (x), (n)*sizeof(double));	\
		int ii;					\
		for(ii=0; ii < (n); ii++){		\
			(x)[ii] = _tmp[(n)-1-ii];	\
		}					\
		Free(_tmp);				\
	}

static re_sas_prior_tp *sas_prior_table = NULL;
static double sas_prior_table_log_integral = 0.0;
static double contour_eps = 100.0 * FLT_EPSILON;
static double cyclic_eps = 0.01;

int re_valid_skew_kurt(double *dist, double skew, double kurt)
{
	int retval;

	retval = (kurt > KURT_LIMIT(skew) ? GMRFLib_TRUE : GMRFLib_FALSE);
	if (dist) {
		if (retval == GMRFLib_FALSE) {
			*dist = ABS(kurt - KURT_LIMIT(skew));
		} else {
			*dist = 0.0;
		}
	}
	return retval;
}
double re_valid_skew(double kurt)
{
	/*
	 * return valid skew 
	 */
	if (kurt > 2.15) {
		return (0.08 * sqrt(-215.0 + 100.0 * kurt));
	} else {
		return (NAN);
	}
}
double re_valid_kurt(double skew)
{
	/*
	 * return minimum valid kurt 
	 */
	return (KURT_LIMIT(skew));
}

int re_sas_skew_kurt(double *skew, double *kurt, double epsilon, double delta)
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

int re_sas_fit_parameters(re_sas_param_tp * param, double *mean, double *prec, double *skew, double *kurt)
{
	/*
	 * for given mean, prec, skew and kurt, fit the parameters in the sas-model. Either none or both of 'skew' and 'kurt' can be NULL. 
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
		if (debug) {
			printf("Have already %.12f %.12f --> %.12f %.12f\n", pin[id][0], pin[id][1], pout[id][0], pout[id][1]);
		}
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
		f.f = re_sas_f;
		f.df = re_sas_df;
		f.fdf = re_sas_fdf;
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
				printf("iter %1d: x %g %g |f| = %g\n", iter, gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1),
				       gsl_blas_dnrm2(s->f));
			}
			status = gsl_multifit_test_delta(s->dx, s->x, eps, eps);
		}
		while (status == GSL_CONTINUE && iter < iter_max);

		err = (gsl_blas_dnrm2(s->f) > 1e-5);
		if (err) {
			if (re_valid_skew_kurt(NULL, target[0], target[1]) == GMRFLib_TRUE) {
				fprintf(stderr, "SHASH fail to fit target skew=%g kurt=%g err=%g, but VALID FAIL!\n", target[0],
					target[1], gsl_blas_dnrm2(s->f));
			}
		}

		pout_tmp[0] = gsl_vector_get(s->x, 0);
		pout_tmp[1] = gsl_vector_get(s->x, 1);
		gsl_vector_free(x);
		gsl_multifit_fdfsolver_free(s);
	}

	double m1, m2, mu, stdev, ss, tmp1;

	tmp1 = exp(pout_tmp[1]);			       /* transform it */
	m1 = M1(pout_tmp[0], tmp1);
	m2 = M2(pout_tmp[0], tmp1);
	ss = sqrt(m2 - SQR(m1));
	mu = *mean - m1 / sqrt(*prec) / ss;
	stdev = 1.0 / sqrt(*prec) / ss;

	param->mu = mu;
	param->stdev = stdev;
	param->epsilon = pout_tmp[0];
	param->delta = tmp1;

	if (use_lookup) {
		perr[id] = err;
		memcpy(pin[id], npin, sizeof(npin));
		memcpy(pout[id], pout_tmp, sizeof(pout_tmp));
	}

	return (err ? !GMRFLib_SUCCESS : GMRFLib_SUCCESS);
}

int re_sas_f(const gsl_vector * x, void *data, gsl_vector * f)
{
	double epsilon, delta, skew, sskew, kurt, kkurt, *ddata = (double *) data;

	epsilon = gsl_vector_get(x, 0);
	delta = exp(gsl_vector_get(x, 1));
	sskew = ddata[0];
	kkurt = ddata[1];

	re_sas_skew_kurt(&skew, &kurt, epsilon, delta);
	gsl_vector_set(f, 0, skew - sskew);
	gsl_vector_set(f, 1, kurt - kkurt);

	return GSL_SUCCESS;
}

int re_sas_df(const gsl_vector * x, void *data, gsl_matrix * J)
{
	double h = GMRFLib_eps(1. / 3.), epsilon, log_delta, skew_h, skew_mh, kurt_h, kurt_mh;

	epsilon = gsl_vector_get(x, 0);
	log_delta = gsl_vector_get(x, 1);

	re_sas_skew_kurt(&skew_h, &kurt_h, epsilon + h, exp(log_delta));
	re_sas_skew_kurt(&skew_mh, &kurt_mh, epsilon - h, exp(log_delta));
	gsl_matrix_set(J, 0, 0, (skew_h - skew_mh) / 2.0 / h);
	gsl_matrix_set(J, 1, 0, (kurt_h - kurt_mh) / 2.0 / h);

	re_sas_skew_kurt(&skew_h, &kurt_h, epsilon, exp(log_delta + h));
	re_sas_skew_kurt(&skew_mh, &kurt_mh, epsilon, exp(log_delta - h));
	gsl_matrix_set(J, 0, 1, (skew_h - skew_mh) / 2.0 / h);
	gsl_matrix_set(J, 1, 1, (kurt_h - kurt_mh) / 2.0 / h);

	return GSL_SUCCESS;
}

int re_sas_fdf(const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J)
{
	re_sas_f(x, data, f);
	re_sas_df(x, data, J);

	return GSL_SUCCESS;
}

/* 
 *  The contourLine-code is taken from R: plot3d.c, with modifications.
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

re_contour_tp *contourLines2(double *x, int nx, double *y, int ny, double zc, CLP * segmentDB)
{
	assert(segmentDB);

	double xend, yend;
	int i, ii, j, jj, ns, dir;
	CLP seglist, seg, s, start, end;

	/*
	 * Begin following contours: Grab a segment, Follow its tail, Follow its head, Save the contour 
	 */
	re_contour_tp *c = Calloc(1, re_contour_tp);
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
				c->cyclic[c->nc - 1] = (ABS(x[0] - x[ns]) < cyclic_eps && ABS(y[0] - y[ns]) < cyclic_eps);
				c->length[c->nc - 1] = len;
				c->x[c->nc - 1] = x;
				c->y[c->nc - 1] = y;
			}
		}
	}
	return c;
}

re_contour_tp *contourLines(double *x, int nx, double *y, int ny, double *z, double zc)
{
	CLP *seg;
	re_contour_tp *c;
	seg = contourLines1(x, nx, y, ny, z, zc);
	c = contourLines2(x, nx, y, ny, zc, seg);
	Free(seg);

	return c;
}

int re_print_contourLines(FILE * fp, re_contour_tp * c)
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

int re_free_contourLines(re_contour_tp * c)
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

int re_join_contourLines(re_contour_tp * c)
{
	/*
	 * modify the contourLines in 'c' so its one object 
	 */

	int debug = 0;
	int i, j, k;

	if (!c || (c && c->nc == 0)) {
		return GMRFLib_SUCCESS;
	}

	if (c->nc == 1) {
		/*
		 * force them to be cyclic if close
		 */
		if ((c->x[0][0] != c->x[0][c->ns[0]] ||
		     c->y[0][0] != c->y[0][c->ns[0]]) && (ABS(c->x[0][0] - c->x[0][c->ns[0]]) < cyclic_eps
							  && ABS(c->y[0][0] - c->y[0][c->ns[0]]) < cyclic_eps)) {
			c->x[0] = Realloc(c->x[0], c->ns[0] + 2, double);
			c->y[0] = Realloc(c->y[0], c->ns[0] + 2, double);
			c->x[0][c->ns[0] + 1] = c->x[0][0];
			c->y[0][c->ns[0] + 1] = c->y[0][0];
			c->ns[0]++;

			/*
			 * update the length 
			 */
			c->length[0] = 0.0;
			for (i = 0; i < c->ns[0]; i++) {
				c->length[0] += sqrt(SQR(c->x[0][i] - c->x[0][i + 1]) + SQR(c->y[0][i] - c->y[0][i + 1]));
			}
		}

		return GMRFLib_SUCCESS;
	}

	if (c->nc <= 1) {
		return GMRFLib_SUCCESS;
	}

	if (debug) {
		printf("ENTER WITH\n");
		re_print_contourLines(stdout, c);
	}

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

	int imin = 0, jmin = 0, imindir = 0, jmindir = 0;
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

	memcpy(&(c->x[imin][c->ns[imin] + 1]), &(c->x[jmin][0]), (c->ns[jmin] + 1) * sizeof(double));
	memcpy(&(c->y[imin][c->ns[imin] + 1]), &(c->y[jmin][0]), (c->ns[jmin] + 1) * sizeof(double));
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
	c->cyclic[imin] = (ABS(c->x[imin][0] - c->x[imin][c->ns[imin]]) < cyclic_eps && ABS(c->y[imin][0] - c->y[imin][c->ns[imin]]) < cyclic_eps);
	c->nc--;

	Free(dists);
	Free(idir);
	Free(jdir);

	return re_join_contourLines(c);
}

int re_contour_get_direction(double *x, double *y, int nseg)
{
	int i;
	double area = 0.0;
	for (i = 0; i < nseg; i++) {
		area += (x[i + 1] - x[i]) * (y[i + 1] + y[i]) / 2.0;
	}
	return (area >= 0.0 ? 1 : -1);
}
double re_point_on_contour(re_contour_tp * c, double skew, double kurt)
{
#define WRAPIT(_i) (((_i) + c->ns[ic]+1)%(c->ns[ic]+1))
#define INBETWEEN(_x, _x0, _x1) (((_x) + contour_eps >= (_x0) && (_x) - contour_eps <= (_x1)) || ((_x) + contour_eps >= (_x1) && (_x) - contour_eps <= (_x0)))

	/*
	 * locate the corresponding 'point' on the contour; counted as the length around the contour from a well defined starting point. 
	 */

	int i, j, start_idx, ic, direction;
	double len;

	if (c->nc == 0)
		return 0.0;

	assert(c->nc == 1);
	GMRFLib_max_value(c->length, c->nc, &ic);
	assert(ic == 0);

	direction = re_contour_get_direction(c->x[ic], c->y[ic], c->ns[ic]);
	if (c->cyclic[ic]) {
		GMRFLib_min_value(c->y[ic], c->ns[ic] + 1, &start_idx);
	} else {
		start_idx = (direction == 1 ? 0 : c->ns[ic]);
	}

	if (0) {
		if (direction < 0) {
			direction = 1;
			if (!(c->x[ic][0] < c->x[ic][c->ns[ic] + 1])) {
				REV(c->x[ic], c->ns[ic] + 1);
				REV(c->y[ic], c->ns[ic] + 1);
			}
		}
	}

	len = 0.0;
	for (i = start_idx, j = 0; j < c->ns[ic] + 1; i = WRAPIT(i + direction), j++) {
		if (INBETWEEN(skew, c->x[ic][i], c->x[ic][WRAPIT(i + direction)])
		    && INBETWEEN(kurt, c->y[ic][i], c->y[ic][WRAPIT(i + direction)])) {
			len += sqrt(SQR(c->x[ic][i] - skew) + SQR(c->y[ic][i] - kurt));
			break;
		} else {
			len += sqrt(SQR(c->x[ic][i] - c->x[ic][WRAPIT(i + direction)]) + SQR(c->y[ic][i] - c->y[ic][WRAPIT(i + direction)]));
		}
	}
	if (!(len <= 1.1 * c->length[0])) {

		printf("SOMETHING WRONG\n");
		for (i = 0; i < c->ns[ic] + 1; i++) {
			printf("i,x,y %d %.12g %.12g BETWEEN %d\n", i, c->x[ic][i], c->y[ic][i],
			       (INBETWEEN(skew, c->x[ic][i], c->x[ic][WRAPIT(i + direction)])
				&& INBETWEEN(kurt, c->y[ic][i], c->y[ic][WRAPIT(i + direction)])));
		}
		P(start_idx);
		P(direction);
		P(len);
		P(c->length[ic]);
		P(skew);
		P(kurt);
		abort();
	}
#undef WRAPIT
#undef INBETWEEN
	return DMAX(0.0, DMIN(len, c->length[ic]));
}

double re_sas_table_log_integral(double *param)
{
	int i, j;
	double skew, kurt, sum = 0.0, val, *pri = NULL;

	for (i = 1; i < sas_prior_table->nx - 1; i++) {
		for (j = 1; j < sas_prior_table->ny - 1; j++) {
			skew = sas_prior_table->skew[i];
			kurt = sas_prior_table->kurt[i];
			pri = re_sas_evaluate_log_prior(skew, kurt, param);
			val = pri[0];
			Free(pri);
			if (!ISNAN(val)) {
				sum += exp(val) *
				    ((sas_prior_table->skew[i + 1] - sas_prior_table->skew[i - 1]) / 2.0) *
				    ((sas_prior_table->kurt[j + 1] - sas_prior_table->kurt[j - 1]) / 2.0);
			}
		}
	}
	return sum;
}


int re_find_in_sas_prior_table(double *output, double skew, double kurt)
{
	/*
	 * return the level based on 'skew' and 'kurt' 
	 */

	int ix, iy;

	ix = re_find_in_table_general(skew, sas_prior_table->skew, sas_prior_table->nx);
	iy = re_find_in_table_general(kurt, sas_prior_table->kurt, sas_prior_table->ny);

	/*
	 * so the notation is like the one in https://en.wikipedia.org/wiki/Bilinear_interpolation 
	 */
#define TABLE_GENERIC(_i, _j, what) sas_prior_table->what[(_i) + sas_prior_table->nx * (_j)]

	double x1, x2, y1, y2, x, y, fQ11, fQ12, fQ21, fQ22;
	x1 = sas_prior_table->skew[ix];
	x2 = sas_prior_table->skew[ix + 1];
	y1 = sas_prior_table->kurt[iy];
	y2 = sas_prior_table->kurt[iy + 1];
	x = skew;
	y = kurt;

	if (!(x >= x1 && x <= x2) || !(y >= y1 && y <= y2)) {
		output[0] = NAN;
		output[1] = NAN;
		output[2] = NAN;

		return GMRFLib_SUCCESS;
	}
#define INTERPOLATE(idx, what) if (1) {					\
		fQ11 = TABLE_GENERIC(ix, iy, what);			\
		fQ21 = TABLE_GENERIC(ix + 1, iy, what);			\
		fQ12 = TABLE_GENERIC(ix, iy + 1, what);			\
		fQ22 = TABLE_GENERIC(ix + 1, iy + 1, what);		\
		output[idx] = 1.0 / (DMAX(DBL_EPSILON, x2 - x1) * DMAX(DBL_EPSILON, y2 - y1)) *\
			(fQ11 * (x2 - x) * (y2 - y) + fQ21 * (x - x1) * (y2 - y) + fQ12 * (x2 - x) * (y - y1) + fQ22 * (x - x1) * (y - y1)); \
	}

	INTERPOLATE(0, level);
	INTERPOLATE(1, length);
	INTERPOLATE(2, point);
	if (sas_prior_table->logjac) {
		INTERPOLATE(3, logjac);
	} else {
		output[3] = NAN;
	}

	output[0] = DMAX(0.0, output[0]);
	output[1] = DMAX(0.0, output[1]);
	output[2] = DMAX(0.0, DMIN(output[1], output[2]));

#undef INTERPOLATE
#undef TABLE_GENERIC
	return GMRFLib_SUCCESS;
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

int re_sas_table_add_logjac(int debug)
{
	int i, j, off = 2;
	double *logjac;

	assert(sas_prior_table);
	logjac = Calloc(sas_prior_table->nz, double);
	for (i = 0; i < sas_prior_table->nz; i++) {
		logjac[i] = NAN;
	}

#pragma omp parallel for private(i, j)
	for (i = off; i < sas_prior_table->nx - off; i++) {
		if (debug) {
			printf(".");
			fflush(stdout);
		}
		for (j = off; j < sas_prior_table->ny - off; j++) {
			double output[4], ddefault = 0.01, target_level = 0.01, Jacobian, lev[3], poi[3], len[3],
			    dkurt, dskew, dlevel_dskew, dpoint_dskew, dlevel_dkurt, dpoint_dkurt, skew, kurt;
			int iter, iter_max = 100, adaptive = 1;

			skew = sas_prior_table->skew[i];
			kurt = sas_prior_table->kurt[j];

			re_find_in_sas_prior_table(output, skew, kurt);
			lev[0] = output[0];
			len[0] = output[1];
			poi[0] = output[2];

			if (!ISNAN(lev[0]) && !ISNAN(len[0]) && !ISNAN(poi[0])) {
				if (skew <= 0.0) {
					dskew = -ddefault;
					dkurt = 2 * ddefault;
				} else {
					dskew = ddefault;
					dkurt = 2 * ddefault;
				}

				if (skew + dskew > SKEW_MIN && skew + dskew < SKEW_MAX && kurt + dkurt > KURT_MIN && kurt + dkurt < KURT_MAX) {

					do {
						re_find_in_sas_prior_table(output, skew + dskew, kurt);
						if (ISNAN(output[0])) {
							dskew *= 0.9;
						}
					} while (ISNAN(output[0]));

					if (adaptive) {
						for (iter = 0; iter < iter_max; iter++) {
							double f = (iter < iter_max / 2 ? 0.9 : (iter < (iter_max * 3) / 4 ? 0.99 : 0.999));
							re_find_in_sas_prior_table(output, skew + dskew, kurt);
							if (ISNAN(output[0])) {
								dskew /= f;
								re_find_in_sas_prior_table(output, skew + dskew, kurt);
								break;
							}

							if (ABS(output[0] - lev[0]) > target_level) {
								dskew *= f;
							} else {
								dskew /= f;
							}
						}
					}
					// P(ABS(output[0]-lev[0])/target_level);

					lev[1] = output[0];
					len[1] = output[1];
					poi[1] = output[2];

					do {
						re_find_in_sas_prior_table(output, skew, kurt + dkurt);
						if (ISNAN(output[0])) {
							dkurt *= 0.9;
						}
					} while (ISNAN(output[0]));

					if (adaptive) {
						for (iter = 0; iter < iter_max; iter++) {
							double f = (iter < iter_max / 2 ? 0.9 : (iter < iter_max * 3 / 4 ? 0.99 : 0.999));
							re_find_in_sas_prior_table(output, skew, kurt + dkurt);
							if (ISNAN(output[0])) {
								dkurt /= f;
								re_find_in_sas_prior_table(output, skew, kurt + dkurt);
								break;
							}
							if (ABS(output[0] - lev[0]) > target_level) {
								dkurt *= f;
							} else {
								dkurt /= f;
							}
						}
					}
					// P(ABS(output[0]-lev[0])/target_level);

					lev[2] = output[0];
					len[2] = output[1];
					poi[2] = output[2];

					dlevel_dskew = (lev[1] - lev[0]) / dskew;
					dpoint_dskew = (poi[1] / len[1] - poi[0] / len[0]) / dskew;
					dlevel_dkurt = (lev[2] - lev[0]) / dkurt;
					dpoint_dkurt = (poi[2] / len[2] - poi[0] / len[0]) / dkurt;

					Jacobian = ABS(dlevel_dskew * dpoint_dkurt - dpoint_dskew * dlevel_dkurt);
					logjac[i + j * sas_prior_table->nx] = log(Jacobian);

					if (0 && ABS(skew) < 0.1) {
						printf("skew %g kurt %g\n", skew, kurt);
						printf("\t       lev %g len %g point %g\n", lev[0], len[0], poi[0]);
						printf("\tskew:  lev %g len %g point %g\n", lev[1], len[1], poi[1]);
						printf("\tkurt:  lev %g len %g point %g\n", lev[2], len[2], poi[2]);
						printf("\tdskew: lev %g point %g\n", dlevel_dskew, dpoint_dskew);
						printf("\tdkurt: lev %g point %g\n", dlevel_dkurt, dpoint_dkurt);
						printf("\tlog(a x b) = log(Jac) %g %g = %g\n", dlevel_dskew * dpoint_dkurt,
						       dpoint_dskew * dlevel_dkurt, log(Jacobian));
					}
				}
			}
		}
	}

	if (debug) {
		printf("\n");
		fflush(stdout);
	}

	sas_prior_table->logjac = logjac;

	FILE *fp;
	if (debug) {
		fprintf(stderr, "Add logjac tot the prior-table %s\n", SAS_PRIOR_TABLE);
	}
	fp = fopen(SAS_PRIOR_TABLE, "wb");
	fwrite(&(sas_prior_table->nx), sizeof(int), (size_t) 1, fp);
	fwrite(&(sas_prior_table->ny), sizeof(int), (size_t) 1, fp);
	fwrite(&(sas_prior_table->nz), sizeof(int), (size_t) 1, fp);
	fwrite(sas_prior_table->skew, sizeof(double), (size_t) sas_prior_table->nx, fp);
	fwrite(sas_prior_table->kurt, sizeof(double), (size_t) sas_prior_table->ny, fp);
	fwrite(sas_prior_table->level, sizeof(double), (size_t) sas_prior_table->nz, fp);
	fwrite(sas_prior_table->length, sizeof(double), (size_t) sas_prior_table->nz, fp);
	fwrite(sas_prior_table->point, sizeof(double), (size_t) sas_prior_table->nz, fp);
	fwrite(sas_prior_table->logjac, sizeof(double), (size_t) sas_prior_table->nz, fp);
	fclose(fp);

	return GMRFLib_SUCCESS;
}


int re_sas_table_init(double *param)
{
	/*
	 * only read it
	 */
	re_sas_prior_table_core(1, 0, 1, param);
	assert(sas_prior_table->logjac);
	assert(sas_prior_table_log_integral != 0.0);

	return GMRFLib_SUCCESS;
}

int re_sas_table_check(double *param)
{
	/*
	 * this will create the table if its not there, and add the logjac if its not there already
	 */
	re_sas_prior_table_core(1, 1, 1, param);

	int n = 500, i, j, ii;
	double *skew = Calloc(n, double);
	double *kurt = Calloc(n, double);
	double *ld = Calloc(ISQR(n), double);
	double *ldd = Calloc(ISQR(n), double);
	double *lddd = Calloc(ISQR(n), double);
	double *ldddd = Calloc(ISQR(n), double);
	double *lddddd = Calloc(ISQR(n), double);
	double s, k;

	for (s = -1.5, i = 0; i < n; s += 2 * 1.5 / (n - 1.0), i++) {
		skew[i] = s;
		if (i == 0 || i == n - 1)
			printf("skew[%d]  = %g\n", i, skew[i]);
	}
	for (k = 2.1, i = 0; i < n; k += (7 - 2.0) / (n - 1), i++) {
		kurt[i] = k;
		if (i == 0 || i == n - 1)
			printf("kurt[%d]  = %g\n", i, kurt[i]);
	}

#pragma omp parallel for private(j, i, ii)
	for (j = 0; j < n; j++) {
		for (i = 0; i < n; i++) {
			ii = i + j * n;
			// printf("skew kurt %g %g\n", skew[i], kurt[j]);
			if (re_valid_skew_kurt(NULL, skew[i], kurt[j])) {
				double *pri = re_sas_evaluate_log_prior(skew[i], kurt[j], param);
				ld[ii] = pri[0];
				ldd[ii] = pri[1];
				lddd[ii] = pri[2];
				ldddd[ii] = pri[3];
				lddddd[ii] = pri[4];
				Free(pri);
			} else {
				ld[ii] = NAN;
				ldd[ii] = NAN;
				lddd[ii] = NAN;
				ldddd[ii] = NAN;
				lddddd[ii] = NAN;
			}
		}
	}

	FILE *fp = fopen("testing1.dat", "w");
	fprintf(fp, "%d\n%d\n%d\n", n, n, ISQR(n));
	for (i = 0; i < n; i++)
		fprintf(fp, "%g\n", skew[i]);
	for (i = 0; i < n; i++)
		fprintf(fp, "%g\n", kurt[i]);
	for (j = ii = 0; j < n; j++)
		for (i = 0; i < n; i++)
			fprintf(fp, "%g\n", ld[ii++]);

	fclose(fp);
	printf("wrote file 1\n");

	fp = fopen("testing2.dat", "w");
	fprintf(fp, "%d\n%d\n%d\n", n, n, ISQR(n));
	for (i = 0; i < n; i++)
		fprintf(fp, "%g\n", skew[i]);
	for (i = 0; i < n; i++)
		fprintf(fp, "%g\n", kurt[i]);
	for (j = ii = 0; j < n; j++)
		for (i = 0; i < n; i++)
			fprintf(fp, "%g\n", ldd[ii++]);

	fclose(fp);
	printf("wrote file 2\n");

	fp = fopen("testing3.dat", "w");
	fprintf(fp, "%d\n%d\n%d\n", n, n, ISQR(n));
	for (i = 0; i < n; i++)
		fprintf(fp, "%g\n", skew[i]);
	for (i = 0; i < n; i++)
		fprintf(fp, "%g\n", kurt[i]);
	for (j = ii = 0; j < n; j++)
		for (i = 0; i < n; i++)
			fprintf(fp, "%g\n", lddd[ii++]);

	fclose(fp);
	printf("wrote file 3\n");

	fp = fopen("testing4.dat", "w");
	fprintf(fp, "%d\n%d\n%d\n", n, n, ISQR(n));
	for (i = 0; i < n; i++)
		fprintf(fp, "%g\n", skew[i]);
	for (i = 0; i < n; i++)
		fprintf(fp, "%g\n", kurt[i]);
	for (j = ii = 0; j < n; j++)
		for (i = 0; i < n; i++)
			fprintf(fp, "%g\n", ldddd[ii++]);

	fclose(fp);
	printf("wrote file 4\n");

	fp = fopen("testing5.dat", "w");
	fprintf(fp, "%d\n%d\n%d\n", n, n, ISQR(n));
	for (i = 0; i < n; i++)
		fprintf(fp, "%g\n", skew[i]);
	for (i = 0; i < n; i++)
		fprintf(fp, "%g\n", kurt[i]);
	for (j = ii = 0; j < n; j++)
		for (i = 0; i < n; i++)
			fprintf(fp, "%g\n", lddddd[ii++]);
	fclose(fp);
	printf("wrote file 5\n");

	return GMRFLib_SUCCESS;
}
int re_sas_table_create(double *param)
{
	re_sas_prior_table_core(0, 1, 1, param);
	return GMRFLib_SUCCESS;
}
int re_sas_prior_table_core(int read_only, int add_logjac, int debug, double *param)
{
	if (sas_prior_table) {
		return GMRFLib_SUCCESS;
	}

	FILE *fp = NULL;
	fp = fopen(SAS_PRIOR_TABLE, "rb");

	if (fp) {
		size_t ret;
		re_sas_prior_tp *s = Calloc(1, re_sas_prior_tp);

		if (debug) {
			fprintf(stderr, "read table from %s\n", SAS_PRIOR_TABLE);
		}

		ret = fread(&(s->nx), sizeof(int), (size_t) 1, fp);
		assert(ret == (size_t) 1);
		ret = fread(&(s->ny), sizeof(int), (size_t) 1, fp);
		assert(ret == (size_t) 1);
		ret = fread(&(s->nz), sizeof(int), (size_t) 1, fp);
		assert(ret == (size_t) 1);

		if (debug) {
			fprintf(stderr, "\tnx = %1d, ny = %1d, nz = %1d\n", s->nx, s->ny, s->nz);
		}
		s->skew = Calloc(s->nx, double);
		s->kurt = Calloc(s->ny, double);
		s->level = Calloc(s->nz, double);
		s->length = Calloc(s->nz, double);
		s->point = Calloc(s->nz, double);
		s->logjac = Calloc(s->nz, double);

		ret = fread(s->skew, sizeof(double), (size_t) s->nx, fp);
		assert(ret == (size_t) s->nx);
		ret = fread(s->kurt, sizeof(double), (size_t) s->ny, fp);
		assert(ret == (size_t) s->ny);
		ret = fread(s->level, sizeof(double), (size_t) s->nz, fp);
		assert(ret == (size_t) s->nz);
		ret = fread(s->length, sizeof(double), (size_t) s->nz, fp);
		assert(ret == (size_t) s->nz);
		ret = fread(s->point, sizeof(double), (size_t) s->nz, fp);
		assert(ret == (size_t) s->nz);
		ret = fread(s->logjac, sizeof(double), (size_t) s->nz, fp);
		if (ret == (size_t) 0) {
			if (debug) {
				fprintf(stderr, "\tlogjac is not available.\n");
			}
			Free(s->logjac);
		} else {
			if (debug) {
				fprintf(stderr, "\tlogjac is available.\n");
			}
			assert(ret == (size_t) s->nz);
		}
		fclose(fp);
		sas_prior_table = s;

		if (add_logjac && !(s->logjac)) {
			if (debug) {
				fprintf(stderr, "\tadd logjac\n");
			}
			re_sas_table_add_logjac(debug);
		}

		sas_prior_table_log_integral = re_sas_table_log_integral(param);
		if (debug) {
			fprintf(stderr, "sas_prior_table: log(integral) = %g for param = %g\n", sas_prior_table_log_integral, param[0]);
		}

		return GMRFLib_SUCCESS;
	}

	if (read_only) {
		return !GMRFLib_SUCCESS;
	}
//#pragma omp critical
	{
		if (!sas_prior_table) {
			printf("MAKE TABLE...\n");

			re_sas_prior_tp *s = Calloc(1, re_sas_prior_tp);
			s->nx = 2000;
			s->ny = 2000;
			s->nz = s->nx * s->ny;

			s->skew = Calloc(s->nx, double);
			s->kurt = Calloc(s->ny, double);
			s->level = Calloc(s->nz, double);
			s->length = Calloc(s->nz, double);
			s->point = Calloc(s->nz, double);

			int i, j;
			double dd, skew_limit = 4.0, kurt_low = 2, kurt_high = 20;

			dd = 2 * skew_limit / (s->nx - 1.0);
			for (i = 0; i < s->nx; i++) {
				s->skew[i] = -skew_limit + i * dd;
			}

			dd = (kurt_high - kurt_low) / (s->ny - 1.0);
			for (j = 0; j < s->ny; j++) {
				s->kurt[j] = kurt_low + j * dd;
			}

#pragma omp parallel for private(j, i) schedule(dynamic)
			for (j = 0; j < s->ny; j++) {
				for (i = 0; i < s->nx; i++) {
					int k = i + j * s->nx;
					if (re_valid_skew_kurt(NULL, s->skew[i], s->kurt[j])) {
						s->level[k] =
						    re_intrinsic_discrepancy_distance_map(re_intrinsic_discrepancy_distance
											  (s->skew[i], s->kurt[j]));
						assert(!ISNAN(s->level[k]));
					} else {
						s->level[k] = NAN;
					}
				}
			}
			printf("MAKE TABLE... Part 1 done\n");


#pragma omp parallel for private(j, i) schedule(dynamic)
			for (j = 0; j < s->ny; j++) {
				printf("...outer loop j = %d/%d\n", j, s->ny);
				for (i = 0; i < s->nx; i++) {
					int k = i + j * s->nx;
					if (re_valid_skew_kurt(NULL, s->skew[i], s->kurt[j])) {
						double level;

						level = s->level[k];
						if (ISNAN(level)) {
							s->length[k] = NAN;
						} else {
							re_contour_tp *c = NULL;
							double length, point;

							c = contourLines(s->skew, s->nx, s->kurt, s->ny, s->level, level);
							re_join_contourLines(c);
							if (c->nc == 0) {
								s->length[k] = FLT_EPSILON;
								s->point[k] = FLT_EPSILON / 2.0;
							} else {
								if (c->nc != 1) {
									re_print_contourLines(stdout, c);
								}
								assert(c->nc == 1);
								length = c->length[0];
								point = re_point_on_contour(c, s->skew[i], s->kurt[j]);
								re_free_contourLines(c);
								s->length[k] = length;
								s->point[k] = point;
							}
						}
					} else {
						s->length[k] = s->point[k] = NAN;
					}
				}
			}
			printf("MAKE TABLE... Part 2 done\n");

			FILE *fp;
			sas_prior_table = s;
			fprintf(stderr, "Write prior-table to %s\n", SAS_PRIOR_TABLE);
			fp = fopen(SAS_PRIOR_TABLE, "wb");
			assert(fp);
			fwrite(&(s->nx), sizeof(int), (size_t) 1, fp);
			fwrite(&(s->ny), sizeof(int), (size_t) 1, fp);
			fwrite(&(s->nz), sizeof(int), (size_t) 1, fp);
			fwrite(s->skew, sizeof(double), (size_t) s->nx, fp);
			fwrite(s->kurt, sizeof(double), (size_t) s->ny, fp);
			fwrite(s->level, sizeof(double), (size_t) s->nz, fp);
			fwrite(s->length, sizeof(double), (size_t) s->nz, fp);
			fwrite(s->point, sizeof(double), (size_t) s->nz, fp);
			fclose(fp);

			return re_sas_prior_table_core(read_only, add_logjac, debug, param);
		}
	}

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
	re_sas_param_tp param;
	double zero = 0.0, one = 1.0;

	re_sas_fit_parameters(&param, &zero, &one, &skew, &kurt);
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

int re_init(double *p)
{
	re_sas_param_tp param;
	double skew = 0, kurt = 3.0, mean = 0, s = 1;

	re_sas_fit_parameters(&param, &mean, &s, &skew, &kurt);
	re_sas_table_init(p);

	return GMRFLib_SUCCESS;
}

double re_sas_log_prior(double *val, double *param)
{
	double *p, rval;
	static int first = 1;

	if (first) {
#pragma omp critical
		{
			if (first) {
				fprintf(stderr, "INIT SAS PRIOR WITH LAMBDA = %g\n", param[0]);
				re_init(param);
				first = 0;
			}
		}
	}

	p = re_sas_evaluate_log_prior(val[0], val[1], param);
	rval = p[0];
	Free(p);

	return rval;
}

double *re_sas_evaluate_log_prior(double skew, double kurt, double *param)
{
	double output[4], level, length, point, ldens_uniform, ldens_dist, *pri, logjac, lambda = param[0];

	re_find_in_sas_prior_table(output, skew, kurt);
	level = output[0];
	length = output[1];
	point = output[2];
	logjac = output[3];

	if (ISNAN(level) || ISNAN(length) || ISNAN(point) || !re_valid_skew_kurt(NULL, skew, kurt)) {
		pri = Calloc(5, double);
		pri[0] = NAN;
		pri[1] = NAN;
		pri[2] = NAN;
		pri[3] = NAN;
		pri[4] = NAN;

		return pri;
	}


	ldens_uniform = log(1 / length);
	ldens_dist = log(lambda) - lambda * level;	       /* the prior for the distance */

	pri = Calloc(5, double);
	pri[0] = ldens_uniform + ldens_dist + logjac - sas_prior_table_log_integral;
	pri[1] = level;
	pri[2] = length;
	pri[3] = point;
	pri[4] = logjac;

	return pri;
}
