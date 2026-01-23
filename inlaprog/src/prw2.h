#ifndef __INLA_PRW2_H__
#       define __INLA_PRW2_H__

#       undef __BEGIN_DECLS
#       undef __END_DECLS
#       ifdef __cplusplus
#              define __BEGIN_DECLS extern "C" {
#              define __END_DECLS }
#       else
#              define __BEGIN_DECLS			       /* empty */
#              define __END_DECLS			       /* empty */
#       endif

#       include <stdio.h>
#       include <stdlib.h>
#       include "GMRFLib/GMRFLib.h"

__BEGIN_DECLS

/*
 * same storage as lapack routines
 */
// general band
#       define BM_IDX(m_, i_, j_) (((i_) - (j_) + (m_)->ubw) + (j_) * (m_)->ldim)
#       define BM_PTR(m_, i_, j_) ((m_)->band + BM_IDX(m_, i_, j_))
#       define BM_VAL(m_, i_, j_) (*BM_PTR(m_, i_, j_))
// symmetric ones, lower storage, can call with (i,j) and (j,i)
#       define sBM_IDX(m_, i_, j_) BM_IDX(m_, IMAX(i_, j_), IMIN(i_, j_))
#       define sBM_PTR(m_, i_, j_) ((m_)->band + sBM_IDX(m_, i_, j_))
#       define sBM_VAL(m_, i_, j_) (*sBM_PTR(m_, i_, j_))

/*
 *
 */
    typedef struct {
	int n;
	int lbw;
	int ubw;
	int ldim;
	int memlen;
	double *band;
} inla_bm_tp;


typedef struct {
	double range;
	inla_bm_tp *Q;
} inla_prw2_cache_tp;

typedef struct {
	int n;
	double *loc;
	double *h;
	double h_size;
	inla_bm_tp *C;
	inla_bm_tp *C_tilde;
	inla_bm_tp *G;
	inla_bm_tp *M;

	double **log_prec_omp;
	double **log_range_omp;

	GMRFLib_graph_tp *graph;
	inla_prw2_cache_tp **cache;
} inla_prw2_arg_tp;


double inla_Qfunc_prw2(int thread_id, int i, int j, double *values, void *arg);
double priorfunc_prw2_pcprior_range(double *x, double *parameters);
double priorfunc_prw2_pcprior_range_calibrate(double r0, double alpha, double h_size);
double priorfunc_prw2_pcprior_range_calibrate_helper(double lambda, double r0, double alpha, double h_size);
inla_bm_tp *inla_bm_alloc(int n, int lbw, int ubw);
inla_bm_tp *inla_bm_chol(inla_bm_tp * A, inla_bm_tp * chol);
inla_bm_tp *inla_bm_duplicate(inla_bm_tp * A);
inla_bm_tp *inla_bm_mm(inla_bm_tp * A, inla_bm_tp * B, inla_bm_tp * AB);
inla_bm_tp *inla_bm_mmm(inla_bm_tp * A, inla_bm_tp * B, inla_bm_tp * C, inla_bm_tp * ABC);
inla_bm_tp *inla_bm_partial_inv(inla_bm_tp * chol, inla_bm_tp * inv);
inla_bm_tp *inla_bm_sum(double a, inla_bm_tp * A, double b, inla_bm_tp * B, double c, inla_bm_tp * C, inla_bm_tp * ABC);
inla_bm_tp *inla_bm_sym(inla_bm_tp * A, inla_bm_tp * Asym);
inla_bm_tp *inla_bm_trans(inla_bm_tp * A, inla_bm_tp * At);
inla_bm_tp *inla_prw2_build_Q(int thread_id, inla_prw2_arg_tp * arg);
inla_prw2_arg_tp *inla_prw2_create(int n, double *loc);
void inla_bm_fprintf(FILE * fp, inla_bm_tp * A, char *msg);
void inla_bm_free(inla_bm_tp * A);
void inla_bm_nsolve(inla_bm_tp * chol, double *b, double *x, int nrhs);
void inla_bm_scale(double a, inla_bm_tp * A);
void inla_bm_solve(inla_bm_tp * chol, double *b, double *x);
void inla_bm_test();
void inla_prw2_pcprior_cdf_range(double *range, int n, double lambda, double h_size, double *cdf);
void inla_prw2_pcprior_dist(double *rho, int n, double *d);
void inla_prw2_pcprior_range2rho(double *range, int n, double h_size, double *rho);
void inla_prw2_test(void);


__END_DECLS
#endif
