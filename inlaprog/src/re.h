
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
 *        Haavard Rue
 *        CEMSE Division
 *        King Abdullah University of Science and Technology
 *        Thuwal 23955-6900, Saudi Arabia
 *        Email: haavard.rue@kaust.edu.sa
 *        Office: +966 (0)12 808 0640
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
#define SAS_PRIOR_TABLE "sas-prior-table.dat"
    typedef struct {
	double mu;
	double stdev;
	double delta;
	double epsilon;
} re_sas_param_tp;

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
} re_contour_tp;

typedef struct {
	int nx;
	int ny;
	int nz;

	double *skew;					       /* x */
	double *kurt;					       /* y */
	double *level;					       /* various z's */
	double *length;
	double *point;
	double *logjac;
} re_sas_prior_tp;





CLP *contourLines1(double *x, int nx, double *y, int ny, double *z, double zc);
CLP ctr_newseg(double x0, double y0, double x1, double y1, CLP prev);
CLP ctr_segupdate(double xend, double yend, int dir, int tail, CLP seglist, CLP * seg);
double *re_sas_evaluate_log_prior(double skew, double kurt, double *param);
double bessel_Knu(double alpha, double x);
double re_intrinsic_discrepancy_distance(double skew, double kurt);
double re_intrinsic_discrepancy_distance_map(double distance);
double re_point_on_contour(re_contour_tp * c, double skew, double kurt);
double re_sas_log_prior(double *val, double *param);
double re_valid_kurt(double skew);
double re_valid_skew(double kurt);
int ctr_intersect(double z0, double z1, double zc, double *f);
int ctr_segdir(double xend, double yend, double *x, double *y, int *i, int *j, int nx, int ny);
int re_contour_get_direction(double *x, double *y, int nseg);
int re_dnorm(double *logdens, double *x, int n);
int re_dsas(double *logdens, double *x, int n, double skew, double kurt);
int re_dsas_intern(double *logdens, double *x, int n, double mu, double sigma, double delta, double epsilon);
int re_find_in_sas_prior_table(double *result, double skew, double kurt);
int re_find_in_table_general(double value, double *x, int nx);
int re_free_contourLines(re_contour_tp * c);
int re_init(double *param);
int re_join_contourLines(re_contour_tp * c);
int re_print_contourLines(FILE * fp, re_contour_tp * c);
int re_prior_check(double *param);
int re_read_sas_prior_table(double *param);
int re_sas_df(const gsl_vector * x, void *data, gsl_matrix * J);
int re_sas_f(const gsl_vector * x, void *data, gsl_vector * f);
int re_sas_fdf(const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J);
int re_sas_fit_parameters(re_sas_param_tp * param, double *mean, double *prec, double *skew, double *kurt);
int re_sas_prior_table_core(int read_only, int add_logjac, int debug, double *param);
int re_sas_skew_kurt(double *skew, double *kurt, double epsilon, double delta);
int re_sas_table_add_logjac(int debug);
int re_sas_table_check(double *param);
int re_sas_table_create(double *param);
int re_sas_table_init(double *param);
int re_valid_skew_kurt(double *dist, double skew, double kurt);
re_contour_tp *contourLines(double *x, int nx, double *y, int ny, double *z, double zc);
re_contour_tp *contourLines2(double *x, int nx, double *y, int ny, double zc, CLP * segmentDB);
void K_bessel(double *x, double *alpha, long *nb, long *ize, double *bk, long *ncalc);
void ctr_swapseg(CLP seg);

__END_DECLS
#endif
