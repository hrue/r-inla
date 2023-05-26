
/* pre-opt.h
 * 
 * Copyright (C) 2021-2023 Havard Rue
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
 *
 */

#ifndef __GMRFLib_PREOPT_H__
#define __GMRFLib_PREOPT_H__

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
    typedef enum {
	GMRFLib_PREOPT_TP_F = 1,
	GMRFLib_PREOPT_TP_BETA,
	GMRFLib_PREOPT_TP___VOID = -1
} GMRFLib_preopt_type_types_tp;

typedef struct {
	GMRFLib_preopt_type_types_tp tp;
	int idx;
	int tp_idx;
} GMRFLib_preopt_type_tp;

typedef struct {
	GMRFLib_matrix_tp *int_design;
	double *adj_weights;
	double *hessian;
	double *inverse_hessian;
	double *stdev_corr_neg;
	double *stdev_corr_pos;
	double *cov_m;
	gsl_matrix *H;
	gsl_matrix *eigen_vectors;
	gsl_vector *eigen_values;
	gsl_vector *sqrt_eigen_values;
} GMRFLib_preopt_res_tp;

typedef struct {

	/*
	 * eta* = pA %*% eta length(eta*) = mpred eta = A %*% x, length(eta) = npred
	 * 
	 * B = mpred x npred might be given, this is pA... A = npred x n this matrix is constructed online
	 * 
	 * mnpred = mpred + npred total length for transfering mode_x etc
	 * 
	 * length of data is Npred: Npred = mpred if mpred > 0,  Npred = npred if mpred = 0
	 * 
	 * the AtA matrix is the matrix that goes into the likelihood, and its either buildt from t(A)A or t(pAA)pAA
	 * 
	 */
	int mpred;
	int npred;
	int mnpred;
	int Npred;
	int n;
	int nf;
	int nbeta;

	char **preopt_graph_latent_is_nb;
	char **preopt_graph_like_is_nb;

	GMRFLib_graph_tp *preopt_graph;
	GMRFLib_Qfunc_tp *preopt_Qfunc;
	void *preopt_Qfunc_arg;

	GMRFLib_Qfunc_tp *gcpo_Qfunc;

	GMRFLib_graph_tp *latent_graph;
	GMRFLib_Qfunc_tp *latent_Qfunc;
	void *latent_Qfunc_arg;
	GMRFLib_constr_tp *latent_constr;

	GMRFLib_graph_tp *like_graph;
	GMRFLib_Qfunc_tp *like_Qfunc;
	GMRFLib_Qfunc_tp *like_Qfunc_k;
	void *like_Qfunc_arg;

	GMRFLib_bfunc_tp **bfunc;
	double **like_c;
	double **like_b;
	double **total_b;
	double *total_const;

	int *idx_map_beta;
	int *idx_map_f;

	double **covariate;
	double *prior_precision;

	GMRFLib_Qfunc_tp ***ff_Qfunc;			       /* interaction */
	GMRFLib_Qfunc_tp **f_Qfunc;
	GMRFLib_graph_tp **f_graph;
	double *f_diag;
	GMRFLib_preopt_type_tp *what_type;
	void ***ff_Qfunc_arg;
	void **f_Qfunc_arg;

	GMRFLib_idxval_tp **pA_idxval;
	GMRFLib_idxval_tp **pAA_idxval;
	GMRFLib_idxval_tp **pAAt_idxval;
	GMRFLib_idxval_tp **A_idxval;
	GMRFLib_idxval_tp **At_idxval;
	GMRFLib_idxval_tp ***AtA_idxval;		       /* this is the total (pA%*%A)^T%*%(pA%*%A) */

	GMRFLib_matrix_tp *A;				       /* the model matrix to construct Predictor */
	GMRFLib_matrix_tp *pA;				       /* the matrix to construct APredictor from Predictor */

	double *mode_theta;
	double *mode_x;

	// this is for gcpo with strategy = prior. DO NOT FREE!
	double *gcpo_mask;
	double *gcpo_diag;

} GMRFLib_preopt_tp;

GMRFLib_preopt_type_tp GMRFLib_preopt_what_type(int node, GMRFLib_preopt_tp * a);
double *GMRFLib_preopt_measure_time(int thread_id, GMRFLib_preopt_tp * preopt, double *res, double *test_vector);
double *GMRFLib_preopt_measure_time2(GMRFLib_preopt_tp * preopt);
double GMRFLib_preopt_Qfunc(int thread_id, int node, int nnode, double *UNUSED(values), void *arg);
double GMRFLib_preopt_Qfunc_prior(int thread_id, int node, int nnode, double *UNUSED(values), void *arg);
double GMRFLib_preopt_gcpo_Qfunc(int thread_id, int node, int nnode, double *UNUSED(values), void *arg);
double GMRFLib_preopt_latent_Qfunc(int thread_id, int node, int nnode, double *values, void *arg);
double GMRFLib_preopt_like_Qfunc(int thread_id, int node, int nnode, double *values, void *arg);
double GMRFLib_preopt_like_Qfunc_k(int thread_id, int node, int k, double *UNUSED(values), void *arg);
int GMRFLib_preopt_free(GMRFLib_preopt_tp * preopt);
int GMRFLib_preopt_bnew(int thread_id, double *b, GMRFLib_preopt_tp * preopt);
int GMRFLib_preopt_bnew_like(double *bnew, double *blike, GMRFLib_preopt_tp * arg);
int GMRFLib_preopt_init(GMRFLib_preopt_tp ** preopt, int n, int nf, int **c, double **w,
			GMRFLib_graph_tp ** f_graph, GMRFLib_Qfunc_tp ** f_Qfunc,
			void **f_Qfunc_arg, char *f_sumzero, GMRFLib_constr_tp ** f_constr,
			double *f_diag,
			GMRFLib_Qfunc_tp *** ff_Qfunc, void ***ff_Qfunc_arg,
			int nbeta, double **covariate, double *prior_precision,
			GMRFLib_bfunc_tp ** bfunc, GMRFLib_ai_param_tp * UNUSED(ai_par), char *predictor_pA_fnm,
			GMRFLib_matrix_tp ** global_constr);
int GMRFLib_preopt_predictor(double *predictor, double *latent, GMRFLib_preopt_tp * preopt);
int GMRFLib_preopt_full_predictor(double *predictor, double *latent, GMRFLib_preopt_tp * preopt);
int GMRFLib_preopt_predictor_core(double *predictor, double *latent, GMRFLib_preopt_tp * preopt, int likelihood_only);
int GMRFLib_preopt_predictor_moments(double *mean, double *variance, GMRFLib_preopt_tp * preopt,
				     GMRFLib_problem_tp * problem, double *optional_mean);
int GMRFLib_preopt_test(GMRFLib_preopt_tp * preopt);
int GMRFLib_preopt_update(int thread_id, GMRFLib_preopt_tp * preopt, double *like_b, double *like_c);
int GMRFLib_preopt_test1(int n, int m);

__END_DECLS
#endif
