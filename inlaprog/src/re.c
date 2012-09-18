
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

#define Pq(_arg) (exp(0.25)/sqrt(8.0*M_PI) * (gsl_sf_bessel_Knu(((_arg)+1.0)/2.0, 0.25) + gsl_sf_bessel_Knu(ABS(((_arg)-1.0)/2.0), 0.25)))
#define M1(_epsilon, _delta) (sinh((_epsilon)/(_delta))*Pq(1.0/(_delta)))
#define M2(_epsilon, _delta) (0.5*(cosh(2.0*(_epsilon)/(_delta))*Pq(2.0/(_delta)) -1.0))
#define M3(_epsilon, _delta) (0.25*(sinh(3.0*(_epsilon)/(_delta))*Pq(3.0/(_delta)) - \
				    3.0*sinh((_epsilon)/(_delta))*Pq(1.0/(_delta))))
#define M4(_epsilon, _delta) (0.125*(cosh(4.0*(_epsilon)/(_delta))*Pq(4.0/(_delta)) - \
				     4.0*cosh(2.0*(_epsilon)/(_delta))*Pq(2.0/(_delta)) + 3.0))

#define XMATCH(x0,x1) (fabs((x0)-(x1)) == 0)
#define YMATCH(y0,y1) (fabs((y0)-(y1)) == 0)

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
				if (!pin && !pout && !perr) {
					pin = Calloc(GMRFLib_MAX_THREADS, double *);
					pout = Calloc(GMRFLib_MAX_THREADS, double *);
					perr = Calloc(GMRFLib_MAX_THREADS, int);

					int i;
					for (i = 0; i < GMRFLib_MAX_THREADS; i++) {
						pin[i] = Calloc(2, double);
						pout[i] = Calloc(2, double);
					}
				}
			}
		}
	}

	double npin[2], pout_tmp[2];
	int thread_id, err = 0, debug = 0;

	thread_id = omp_get_thread_num();
	npin[0] = (skew ? *skew : 0.0);
	npin[1] = (kurt ? *kurt : 3.0);

	if (use_lookup && memcmp(npin, pin[thread_id], sizeof(npin)) == 0) {
		/*
		 * we have already computed these values
		 */
		memcpy(pout_tmp, pout[thread_id], sizeof(pout_tmp));
	} else if (!skew && !kurt) {
		pout_tmp[0] = 0.0;			       /* epsilon */
		pout_tmp[1] = 1.0;			       /* delta */
	} else {
		double epsilon = 0.0;
		double delta = 1.0;

		gsl_vector *x = gsl_vector_alloc(2);
		gsl_vector_set(x, 0, 0.0);		       /* initial values */
		gsl_vector_set(x, 1, 1.0);

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
		double eps = 1e-4;
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

		err = (gsl_blas_dnrm2(s->f) > 1e-6);
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
		perr[thread_id] = err;
		memcpy(pin[thread_id], npin, sizeof(npin));
		memcpy(pout[thread_id], pout_tmp, sizeof(pout_tmp));
	}

	return (err ? !GMRFLib_SUCCESS : GMRFLib_SUCCESS);
}
int re_shash_f(const gsl_vector * x, void *data, gsl_vector * f)
{
	double epsilon, delta, skew, sskew, kurt, kkurt, *ddata = (double *) data;

	epsilon = gsl_vector_get(x, 0);
	delta = gsl_vector_get(x, 1);
	sskew = ddata[0];
	kkurt = ddata[1];

	re_shash_skew_kurt(&skew, &kurt, epsilon, delta);
	gsl_vector_set(f, 0, skew - sskew);
	gsl_vector_set(f, 1, kurt - kkurt);

	return GSL_SUCCESS;
}
int re_shash_df(const gsl_vector * x, void *data, gsl_matrix * J)
{
	double h = GMRFLib_eps(1. / 3.), epsilon, delta, skew_h, skew_mh, kurt_h, kurt_mh;

	epsilon = gsl_vector_get(x, 0);
	delta = gsl_vector_get(x, 1);

	re_shash_skew_kurt(&skew_h, &kurt_h, epsilon + h, delta);
	re_shash_skew_kurt(&skew_mh, &kurt_mh, epsilon - h, delta);
	gsl_matrix_set(J, 0, 0, (skew_h - skew_mh) / 2.0 / h);
	gsl_matrix_set(J, 1, 0, (kurt_h - kurt_mh) / 2.0 / h);

	re_shash_skew_kurt(&skew_h, &kurt_h, epsilon, delta + h);
	re_shash_skew_kurt(&skew_mh, &kurt_mh, epsilon, delta - h);
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
   The countourLine-code is taken from R: plot3d.c
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
			 */
			/*
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
inla_countour_tp *contourLines2(double *x, int nx, double *y, int ny, double *z, double zc, CLP * segmentDB)
{
	/*
	 * this function leaks memory, as the second contourLines2()-call, revise the 'segmentDB' object without free objects discared.
	 */

	if (segmentDB == NULL) {
		CLP *seg;
		inla_countour_tp *c;
		seg = contourLines1(x, nx, y, ny, z, zc);
		c = contourLines2(x, nx, y, ny, z, zc, seg);
		Free(seg);
		return c;
	} else {
		double xend, yend;
		int i, ii, j, jj, ns, dir;
		CLP seglist, seg, s, start, end;
		/*
		 * Begin following contours. 
		 */
		/*
		 * 1. Grab a segment 
		 */
		/*
		 * 2. Follow its tail 
		 */
		/*
		 * 3. Follow its head 
		 */
		/*
		 * 4. Save the contour 
		 */

		inla_countour_tp *c = Calloc(1, inla_countour_tp);
		c->level = zc;

		for (i = 0; i < nx - 1; i++)
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
					end->next = NULL;      /* <<< new for 1.2.3 */
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
						s = s->next;
						x[ns] = s->x0;
						y[ns] = s->y0;
						ns++;
					}
					x[ns] = s->x1;
					y[ns] = s->y1;

					if (0) {
						int ii;
						for (ii = 0; ii < ns + 1; ii++) {
							printf("(x,y) = ( %g %g )\n", x[ii], y[ii]);
						}
					}

					int ii;
					double len = 0.0;
					for (ii = 1; ii < ns + 1; ii++) {
						len += sqrt(SQR(x[ii] - x[ii-1]) + SQR(y[ii] - y[ii-1]));
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
		return c;
	}
}
inla_countour_tp *contourLines(double *x, int nx, double *y, int ny, double *z, double zc)
{
	return (contourLines2(x, nx, y, ny, z, zc, NULL));
}
int inla_print_contourLines(FILE *fp, inla_countour_tp *c)
{
	int i, j;

	if (fp == NULL){
		fp = stdout;
	}
	fprintf(fp, "Number of contours %1d for level %g\n", c->nc, c->level);
	for(i=0; i<c->nc; i++){
		fprintf(fp, "\tContour %1d\n", i);
		fprintf(fp, "\tNumber of segments %1d\n", c->ns[i]);
		fprintf(fp, "\tLength %g\n", c->length[i]);
		fprintf(fp, "\tCyclic = %1d\n", c->x[i][0] == c->x[i][c->ns[i]] && c->y[i][0] == c->y[i][c->ns[i]]);
		fprintf(fp, "\tPoints:\n");
		for(j=0; j<c->ns[i]+1; j++){
			fprintf(fp, "\t\t(x,y) = %g %g\n", c->x[i][j], c->y[i][j]);
		}
	}

	return GMRFLib_SUCCESS;
}
