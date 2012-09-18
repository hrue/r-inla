
/* re.h
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
#ifndef __INLA_RE_H__
#define __INLA_RE_H__
#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS					       /* empty */
#define __END_DECLS					       /* empty */
#endif
__BEGIN_DECLS

/* 
 *
 */
    typedef struct {
	double mu;
	double stdev;
	double delta;
	double epsilon;
} re_shash_param_tp;

typedef struct CL {
	struct CL *next;
	double x0;
	double y0;
	double x1;
	double y1;
} CL, *CLP;

typedef struct {
	int nc;
	int *ns;					       /* number of segments */
	int *cyclic;
	double *length;
	double level;
	double **x;
	double **y;
} inla_countour_tp;

int re_valid_skew_kurt(double *dist, double skew, double kurt);
int re_shash_skew_kurt(double *skew, double *kurt, double epsilon, double delta);
int re_shash_fit_parameters(re_shash_param_tp * param, double *mean, double *prec, double *skew, double *kurt);

int re_shash_f(const gsl_vector * x, void *data, gsl_vector * f);
int re_shash_df(const gsl_vector * x, void *data, gsl_matrix * J);
int re_shash_fdf(const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J);

int ctr_intersect(double z0, double z1, double zc, double *f);
CLP ctr_newseg(double x0, double y0, double x1, double y1, CLP prev);
void ctr_swapseg(CLP seg);
int ctr_segdir(double xend, double yend, double *x, double *y, int *i, int *j, int nx, int ny);
CLP ctr_segupdate(double xend, double yend, int dir, int tail, CLP seglist, CLP * seg);
CLP *contourLines1(double *x, int nx, double *y, int ny, double *z, double zc);
inla_countour_tp *contourLines2(double *x, int nx, double *y, int ny, double *z, double zc, CLP * segmentDB);
inla_countour_tp *contourLines(double *x, int nx, double *y, int ny, double *z, double zc);
int inla_print_contourLines(FILE *fp, inla_countour_tp *c);

__END_DECLS
#endif
