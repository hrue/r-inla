
/* inla.h
 * 
 * Copyright (C) 2007-2010 Havard Rue
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
 * RCSId: $Id: inla.h,v 1.238 2010/04/03 12:33:29 hrue Exp $
 *
 */
#ifndef __INLA_H__
#define __INLA_H__
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
#include "iniparser.h"
#include "dictionary.h"
#include "strlib.h"
#include "ar.h"
#define LOG_NORMC_GAUSSIAN (-0.91893853320467274178032973640560)	/* -1/2 * log(2*pi) */
#define INLA_FAIL  1
#define INLA_OK    0
#define FIFO_GET "inla-mcmc-fifo-get"
#define FIFO_PUT "inla-mcmc-fifo-put"
#define FIFO_GET_DATA "inla-mcmc-fifo-get-data"
#define FIFO_PUT_DATA "inla-mcmc-fifo-put-data"

/* 
 *
 */
    typedef enum {
	/*
	 * Failure time
	 */
	SURV_EVENT_FAILURE = 1,				       /* never change!!! */

	/*
	 * Right cencored
	 */
	SURV_EVENT_RIGHT = 0,				       /* never change!!! */

	/*
	 * Left cencored
	 */
	SURV_EVENT_LEFT = 2,

	/*
	 * Interval cencored
	 */
	SURV_EVENT_INTERVAL = 3
} inla_surv_event_tp;

typedef enum {
	INLA_MODE_DEFAULT = 0,
	INLA_MODE_MCMC,
	INLA_MODE_HYPER,
	INLA_MODE_QINV,
	INLA_MODE_QSOLVE,
	INLA_MODE_QREORDERING,
	INLA_MODE_QSAMPLE,
	INLA_MODE_FINN,
	INLA_MODE_GRAPH,
	INLA_MODE_SASPRIOR,
	INLA_MODE_TESTIT = 99
} inla_mode_tp;

typedef enum {
	/*
	 * map from internal representation to the user representation and link-functions. EX: precision = exp(log_prec), or the
	 * inverse-link-function
	 */
	MAP_FORWARD = 1,
	INVLINK = 1,
	/*
	 * the inverse of the _FORWARD map, EX: log_prec = log(precision), or the link-function
	 */
	MAP_BACKWARD = 2,
	LINK = 2,
	/*
	 * the derivative of the forward_map, EX: exp(log_prec). or the derivative of the inverse link-function.
	 */
	MAP_DFORWARD = 3,
	DINVLINK = 3,
	/*
	 * return 1 is monotone increasing and 0 otherwise
	 */
	MAP_INCREASING = 4
} map_arg_tp;

typedef double map_func_tp(double arg, map_arg_tp typ, void *param);

typedef struct {
	const char *name;
	map_func_tp *func;
} map_table_tp;


typedef struct {
	size_t len;
	char *contents;
} inla_file_contents_tp;

typedef struct {
	double *d;					       /* the d-array */
	int ndata;					       /* length of data (from file) */
	double *y;					       /* general responce */

	/*
	 * y ~ Poisson(E*exp(x)) 
	 */
	double *E;

	/*
	 * y ~ Binomial(nb, p(x))
	 */
	double *nb;

	/*
	 * y ~ BinomialRE()
	 */
	double **log_prec_binomialre;

	/*
	 * y ~ CBinomial(k n, p(x))
	 */
	double *cbinomial_k;
	double *cbinomial_n;

	/*
	 * y ~ Binomial test 
	 */
	double **link_psi;

	/*
	 * y ~ Neg.Binomial(n, p(x)), n=size is a hyperparameter (overdispersion)
	 */
	double **log_size;

	/*
	 * y ~ Normal(x, 1/(weight*prec)), also used for the log-normal
	 */
	double **log_prec_gaussian;
	double *weight_gaussian;			       /* weights for the gaussian: Variance = 1/(weight*prec) */

	/*
	 * y ~ Logistic, Variance = 1/ (weight * prec)
	 */
	double **log_prec_logistic;
	double *weight_logistic;

	/*
	 * y ~ T_dof(x, 1/(weight*prec)), where T_dof has Variance=1
	 */
	double **log_prec_t;
	double **dof_intern_t;
	double *weight_t;				       /* weights for the t: Variance = 1/(weight*prec) */

	/*
	 * y ~ Tstrata 
	 */
	double ***log_prec_tstrata;			       /* array of log_prec_t */
	double **dof_intern_tstrata;			       /* common dof */
	double *weight_tstrata;				       /* weight_t */
	double *strata_tstrata;				       /* strata_t (type int) */


	/*
	 * stochvol_t: y ~ Student-t_dof (x) with Variance=1. dof_intern_svt = log(dof - 2).
	 */
	double **dof_intern_svt;
	/*
	 * stochvol_nig: y ~ normal.inverse.gaussian(x) with Variance=1.
	 *
	 *  shew is the beta-parameter, which determine the skewness
	 *  shape is psi = 1+exp(shape), which determine the shape. 
	 */
	double **skew_intern_svnig;
	double **shape_intern_svnig;
	/*
	 * survival: event is 1 is y is the failure time or 0 if its right sensored 
	 */
	double **alpha_intern;				       /* For the Weibull, PS and loglogistic */
	double *truncation;				       /* for survival */
	double *event;					       /* one of inla_surv_event_tp */
	double *lower;					       /* for survival */
	double *upper;					       /* for survival */
	// dobule *time == y
	double **p_intern;				       /* For the L_WEIBULL_CURE */

	/*
	 * zero-inflated Poission/Binomial/inomial/BetaBinomial version 0/1/2...
	 */
	double **prob_intern;
	double **zeroinflated_alpha_intern;		       /* alpha = exp(alpha_intern) */
	double **zeroinflated_delta_intern;		       /* delta = exp(delta_intern) */


	/*
	 * the zero-n-inflated binomial 2 
	 */
	double **zero_n_inflated_alpha1_intern;
	double **zero_n_inflated_alpha2_intern;

	/*
	 * the overdispersion parameter for the betabinomial, \rho = 1/(a+b+1).
	 */
	double **betabinomial_overdispersion_intern;

	/*
	 * the precision parameter for the beta, \phi = exp(theta)
	 */
	double **beta_precision_intern;

	/*
	 * for the (asymmetric) laplace
	 */
	double **log_tau_laplace;
	double *weight_laplace;
	double alpha_laplace;
	double epsilon_laplace;
	double gamma_laplace;

	/*
	 * Skew-Normal
	 */
	double **shape_skew_normal;
	double shape_max_skew_normal;			       /* maximum value for |shape| allowed */
	double **log_prec_skew_normal;
	double *weight_skew_normal;			       /* weights for the skew_normal: Variance = 1/(weight*prec) [for a=0] */

	/*
	 * GEV 
	 */
	double *weight_gev;				       /* weights for the skew_normal: Variance propto 1/(weight*prec) */
	double **log_prec_gev;				       /* log prec for gev */
	double **xi_gev;				       /* the shape-parameter */
	double gev_scale_xi;				       /* scaling of the shape-parameter */

	/*
	 * Log gamma frailty
	 */
	double **log_prec_loggamma_frailty;

	/*
	 * iid gamma 
	 */
	double *iid_gamma_weight;
	double **iid_gamma_log_shape;
	double **iid_gamma_log_rate;

	/*
	 * iid logitbeta 
	 */
	double **iid_logitbeta_log_a;
	double **iid_logitbeta_log_b;

	/*
	 * Sinh-asinh
	 */
	double *sas_weight;
	double **sas_log_prec;
	double **sas_skew;
	double **sas_kurt;

	/*
	 * y ~ Circular Normal, with precision parameter: weight*prec
	 */
	double **log_prec_circular_normal;
	double *weight_circular_normal;

	/*
	 * y ~ Circular Cauchy, with precision parameter: weight*prec, related to the rho-parameter; see doc. 0<weight<=1
	 */
	double **log_prec_wrapped_cauchy;
	double *weight_wrapped_cauchy;

	/*
	 * generalised poisson
	 */
	double **gpoisson_overdispersion;
	double **gpoisson_p;

	/*
	 * test-implementations  
	 */
	double **test_binomial_1_s, **test_binomial_1_e;

	/*
	 * Gamma 
	 */
	double **gamma_log_prec;
	double *gamma_weight;				       /* the scalings 's' */

	/*
	 * MIX ~ Normal(x, 1/prec)
	 */
	double **mix_log_prec_gaussian;
} Data_tp;

/* 
   this is needed so we can identify each component in the model
 */
typedef enum {
	INVALID_COMPONENT = 0,
	L_GAUSSIAN,					       /* likelihood-models */
	L_LOGISTIC,
	L_SKEWNORMAL,
	L_GEV,
	L_T,
	L_TSTRATA,
	L_POISSON,
	L_GPOISSON,
	L_BINOMIAL,
	L_BINOMIALTEST,					       /* test-version of the binomial! */
	L_CBINOMIAL,					       /* clumped binomial */
	L_ZEROINFLATEDBINOMIAL0,
	L_ZEROINFLATEDBINOMIAL1,
	L_ZEROINFLATEDBINOMIAL2,
	L_ZEROINFLATEDBETABINOMIAL2,
	L_GAMMA,
	L_BETA,
	L_BETABINOMIAL,
	L_NBINOMIAL,
	L_ZEROINFLATEDNBINOMIAL0,
	L_ZEROINFLATEDNBINOMIAL1,
	L_ZEROINFLATEDNBINOMIAL2,
	L_STOCHVOL,
	L_STOCHVOL_T,
	L_STOCHVOL_NIG,
	L_LOGPERIODOGRAM,
	L_EXPONENTIAL,
	L_WEIBULL,
	L_LOGLOGISTIC,
	L_LOGNORMAL,
	L_ZEROINFLATEDPOISSON0,
	L_ZEROINFLATEDPOISSON1,
	L_ZEROINFLATEDPOISSON2,
	L_ZERO_N_INFLATEDBINOMIAL2,
	L_WEIBULL_CURE,					       /* Patrick and Silvia's model */
	L_LAPLACE,
	L_LOGGAMMA_FRAILTY,
	L_SAS,
	L_IID_GAMMA,
	L_IID_LOGITBETA,
	L_CIRCULAR_NORMAL,
	L_WRAPPED_CAUCHY,
	L_TEST_BINOMIAL_1,
	L_GAUSSIAN_WINDOW,
	F_RW2D,						       /* f-models */
	F_BESAG,
	F_BESAG2,					       /* the [a*x, x/a] model */
	F_BESAGPROPER,
	F_SEASONAL,
	F_IID,
	F_2DIID,
	F_IID1D,
	F_IID2D,
	F_IID3D,
	F_IID4D,
	F_IID5D,
	F_RW1,
	F_RW2,
	F_CRW2,
	F_AR1,
	F_AR,
	F_OU,
	F_Z,
	F_ZADD,
	F_BYM,
	F_GENERIC0,
	F_GENERIC1,
	F_GENERIC2,
	F_MATERN2D,
	F_SPDE,
	F_SPDE2,
	F_COPY,
	F_MEC,
	F_MEB,
	F_R_GENERIC,
	P_LOGGAMMA,					       /* priors */
	P_GAUSSIAN,
	P_MVGAUSSIAN,
	P_MVNORM,
	P_MINUSLOGSQRTRUNCGAUSSIAN,
	P_FLAT,
	P_WISHART1D,
	P_WISHART2D,
	P_WISHART3D,
	P_WISHART4D,
	P_WISHART5D,
	P_LOGFLAT,
	P_LOGIFLAT,
	P_NONE,
	P_BETACORRELATION,
	P_SASPRIOR,
	P_LOGITBETA,
	P_EXPRESSION,
	P_TABLE,
	P_JEFFREYS_T_DF,
	G_EXCHANGEABLE,					       /* group models */
	G_AR1,
	G_RW1,
	G_RW2,
	G_AR, 
	G_BESAG, 
	MIX_GAUSSIAN, 
} inla_component_tp;


/* 
   priors are defined using this template. return log(pi(precision, parameters....))
 */
typedef double inla_priorfunc_tp(double *param, double *parameters);

typedef struct {
	inla_component_tp id;				       /* prior Id */
	char *name;					       /* name of prior */
	double *parameters;				       /* the parameters */
	char *to_theta;					       /* R-code */
	char *from_theta;				       /* R-code */
	inla_priorfunc_tp *priorfunc;			       /* Either a priorfunction, or */
	char *expression;				       /* an alternative expression */
} Prior_tp;

/* 
   This is the macro to evaluate the prior. One and only one of `priorfunc' and `expression' is non-NULL, so we use that one
 */
#define PRIOR_EVAL(p_, arg_) ((p_).priorfunc ?				\
			      (p_).priorfunc(arg_, (p_).parameters)  :	\
			      inla_eval((p_).expression, arg_))

typedef struct {
	GMRFLib_tabulate_Qfunc_tp *tab;
	double **beta;
	double **log_prec;

	int n;
	double *eigenvalues;
	double max_eigenvalue;
	double min_eigenvalue;
} Generic1_tp;

typedef struct {
	GMRFLib_tabulate_Qfunc_tp *tab;
	double **log_prec;
	double **h2_intern;
	int n;						       /* size of graph */
	int N;						       /* total size: N=2*n */
} Generic2_tp;

typedef struct {
	char *name;
	char *type;
} File_tp;

typedef struct {
	int cpo;					       /* output CPO */
	int dic;					       /* output DIC */
	int summary;					       /* output marginal summaries (mean, stdev, etc) */
	int return_marginals;				       /* output detailed marginal density (even though they are computed) */
	int hyperparameters;				       /* compute also the marginal for the hyperparameters */
	int kld;					       /* output the (symmetric) kld between marginals */
	int mlik;					       /* compute the marginal likelihood? */
	int q;						       /* output image of the Q-matrix */
	int graph;					       /* output the graph */
	int config;					       /* output the configurations */
	int nquantiles;					       /* compute cdfs and/or quantiles; max 10 */
	int ncdf;
	double *quantiles;
	double *cdf;
} Output_tp;

typedef struct inla_tp_struct inla_tp;			       /* need it like this as they point to each other */

typedef struct {
	char *data_likelihood;
	GMRFLib_uchar variant;

	char *link;
	map_func_tp *predictor_invlinkfunc;

	inla_component_tp data_id;
	File_tp data_file;
	File_tp weight_file;
	Prior_tp data_prior;
	Prior_tp data_prior0;
	Prior_tp data_prior1;
	Prior_tp data_prior2;
	Prior_tp *data_nprior;
	Data_tp data_observations;
	int data_fixed;
	int data_fixed0;
	int data_fixed1;
	int data_fixed2;
	int *data_nfixed;
	int data_ntheta;
	GMRFLib_logl_tp *loglikelihood;
	double *offset;
	inla_tp *mb;					       /* to get the off_.... */

	/* 
	 * the re-extention
	 */
	int mix_use;
	inla_component_tp mix_id;
	GMRFLib_logl_tp *mix_loglikelihood;
	Prior_tp mix_prior;
	int mix_fixed;
	int mix_ntheta;
} Data_section_tp;


typedef struct {
	int n;
	GMRFLib_Qfunc_tp *Qfunc;
	void *Qfunc_arg;
	double **rho_intern;
} inla_rwc2_def_tp;

typedef struct
{
	int n;
	int p;

	double **log_prec;
	double ***pacf_intern;
	
	/* 
	 * these are the stored values, all for prec = 1
	 */
	double **hold_pacf_intern;			       /* [i][id][0] */
	double **hold_Q;				       /* dim = 2*p + 1 */
	double **hold_Qmarg;				       /* dim = p */
}
	ar_def_tp;

typedef struct {
	int N;
	int ngroup;
	int type;
	int cyclic;
	GMRFLib_graph_tp *graph;
	GMRFLib_Qfunc_tp *Qfunc;
	void *Qfunc_arg;
	double **group_rho_intern;
	double **group_prec_intern;
	GMRFLib_rwdef_tp *rwdef;
	ar_def_tp *ardef;
} inla_group_def_tp;

typedef struct {
	double precision;
	double **beta;

	GMRFLib_Qfunc_tp *Qfunc;
	void *Qfunc_arg;

	map_func_tp *map_beta;
	void *map_beta_arg;
} inla_copy_arg_tp;


struct inla_tp_struct {
	/*
	 * General stuff 
	 */
	int verbose;
	int strategy;

	/*
	 * parameters for global_nodes
	 */
	GMRFLib_global_node_tp gn;

	/*
	 * reuse the mode-stuff 
	 */
	unsigned char *sha1_hash;
	unsigned char *sha1_hash_file;
	int ntheta_file;
	int theta_counter_file;
	int reuse_mode;
	int reuse_mode_but_restart;
	int fixed_mode;					       /* if TRUE, then treat all thetas as known and fixed, otherwise, do as usual... */
	double *theta_file;
	double *x_file;
	int nx_file;

	/*
	 * Expert options 
	 */
	int expert_cpo_manual;
	int *expert_cpo_idx;				       /* list of indices */
	int expert_n_cpo_idx;				       /* the length */
	double expert_diagonal_emergencey;

	/*
	 * section Problem
	 */
	char *name;					       /* Id of the problem */
	char *dir;					       /* Where to store the files */
	Output_tp *output;				       /* Global output-options */

	/*
	 * type Predictor
	 */
	int predictor_n;				       /* dimension of \eta */
	int predictor_m;				       /* dimension of \tilde{\eta}, the extended predictor */
	int predictor_ndata;				       /* dimension of the data, ie _n if _m is 0, otherwise _m */

	char *predictor_tag;				       /* the tag */
	char *predictor_dir;				       /* the directory */
	double **predictor_log_prec;
	Prior_tp predictor_prior;
	int *predictor_cross_sumzero;
	int predictor_compute;
	int predictor_fixed;
	int predictor_user_scale;
	map_func_tp **predictor_invlinkfunc;		       /* these are rebuilt */
	Output_tp *predictor_output;
	double *offset;					       /* the offset y ~ f(eta + offset) */
	double *link_fitted_values;			       /* the index for the link function for missing observations */

	char *predictor_Aext_fnm;			       /* extension: filename for the Amatrix */
	double predictor_Aext_precision;		       /* extension: precision for the Amatrix */

	/*
	 * type Data 
	 */
	int ds;						       /* current data-section (when reading) */
	int nds;					       /* number of data-sections in total */
	int data_ntheta_all;
	int gaussian_data;				       /* corresponds to ai_par->gaussian_data */
	double *d;
	void **loglikelihood_arg;
	Data_section_tp *data_sections;
	GMRFLib_logl_tp **loglikelihood;

	/*
	 * type Ffield 
	 */
	int nf;
	char **f_tag;
	char **f_dir;
	char **f_modelname;
	int **f_c;
	int *f_n;					       /* basic size of f-model */
	int *f_N;					       /* size of f-model after internal augmentation */
	int *f_Ntotal;					       /* Grand total of f-model, including augmentation and replication */
	int *f_nrow;					       /* if 2d */
	int *f_ncol;					       /* if 2d */
	int *f_nrep;					       /* number of replicates */
	int *f_ngroup;
	int *f_group_model;
	int *f_group_cyclic;
	int *f_group_order;
	GMRFLib_graph_tp **f_group_graph;
	int *f_order;
	double **f_locations;
	double **f_weights;
	double **f_scale;
	GMRFLib_graph_tp **f_graph;
	GMRFLib_graph_tp **f_graph_orig;
	GMRFLib_Qfunc_tp **f_Qfunc;
	GMRFLib_Qfunc_tp **f_Qfunc_orig;
	void **f_Qfunc_arg;
	void **f_Qfunc_arg_orig;
	GMRFLib_bfunc2_tp **f_bfunc2;
	Prior_tp **f_prior;
	char *f_sumzero;
	GMRFLib_constr_tp **f_constr;
	GMRFLib_constr_tp **f_constr_orig;
	double *f_diag;
	double *f_rankdef;
	double ****f_theta;
	map_func_tp ***f_theta_map;
	void ***f_theta_map_arg;
	int *f_si;
	int *f_compute;
	int **f_fixed;
	double **f_initial;
	inla_component_tp *f_id;
	int *f_ntheta;
	int *f_cyclic;
	int *f_nu;
	int *f_Torder;
	char **f_Tmodel;
	int *f_Korder;
	char **f_Kmodel;
	void **f_model;
	char **f_of;
	char **f_same_as;
	double *f_precision;
	Output_tp **f_output;
	inla_file_contents_tp **f_id_names;

	GMRFLib_Qfunc_tp ***ff_Qfunc;			       /* interactions */
	void ***ff_Qfunc_arg;

	/*
	 * type Linear
	 */
	int nlinear;
	char **linear_tag;
	char **linear_dir;
	double **linear_covariate;
	double *linear_precision;
	double *linear_mean;
	int *linear_compute;
	Output_tp **linear_output;

	/*
	 * linear combinations 
	 */
	int nlc;					       /* number of linear combinations */
	int lc_derived_only;				       /* use only the derived lincombs ? */
	int lc_derived_correlation_matrix;		       /* compute correlations ? */
	char **lc_tag;					       /* the tags */
	double *lc_prec;				       /* the `high' precision */
	char **lc_dir;
	double *lc_order;
	Output_tp **lc_output;
	GMRFLib_lc_tp **lc_lc;
	double *lc_derived_c;				       /* optional: correlation for the lincombs (derived) */

	/*
	 * The final model 
	 */
	GMRFLib_hgmrfm_tp *hgmrfm;
	int ntheta;
	double ***theta;
	char **theta_tag;
	char **theta_tag_userscale;
	char **theta_dir;
	char **theta_from;
	char **theta_to;
	map_func_tp **theta_map;
	void **theta_map_arg;
	int *off_compute;
	char **off_modelname;
	int *off_id;

	/*
	 * INLA 
	 */
	GMRFLib_ai_param_tp *ai_par;
	/*
	 * results 
	 */
	GMRFLib_density_tp **density;
	GMRFLib_density_tp **gdensity;
	GMRFLib_density_tp **density_hyper;
	GMRFLib_density_tp **density_lin;
	GMRFLib_ai_cpo_tp *cpo;
	GMRFLib_ai_dic_tp *dic;
	GMRFLib_ai_marginal_likelihood_tp mlik;
	GMRFLib_ai_neffp_tp neffp;

	/*
	 * index-table 
	 */
	int idx_ntot;					       /* total length of x before lincomb */
	int idx_tot;					       /* total number of terms (ffield/linear) */
	char **idx_tag;					       /* the tags */
	int *idx_start;					       /* the starting index */
	int *idx_n;					       /* the length */
	map_stri idx_hash;				       /* the hash-table */

	/*
	 * misc output from INLA() 
	 */
	GMRFLib_ai_misc_output_tp *misc_output;
};

typedef struct {
	GMRFLib_graph_tp *graph;
	int si;						       /* expert option... */
	double **log_prec;
} inla_besag_Qfunc_arg_tp;

typedef struct {
	GMRFLib_graph_tp *graph;
	double **log_prec;
	double **log_a;					       /* the parameter a */
	double precision;				       /* the copy precision */
} inla_besag2_Qfunc_arg_tp;

typedef struct {
	int n;
	int N;
	double **log_prec_iid;
	inla_besag_Qfunc_arg_tp *besag_arg;
} inla_bym_Qfunc_arg_tp;

typedef struct {
	GMRFLib_graph_tp *graph;
	double **log_prec;				       /* precision */
	double **log_diag;				       /* the value added to the diagonal. diag = 0 is the besag model */
} inla_besag_proper_Qfunc_arg_tp;

typedef struct {
	/*
	 * the AR(1) model: X_t = phi * X_t-1 + Z_t. The arguments are Var(Z_t) = 1/exp(log_precision), and phi_intern =
	 * logit((phi+1)/2). 
	 */
	int n;
	int cyclic;
	double **log_prec;				       /* theta[0] */
	double **phi_intern;				       /* theta[1] */
} inla_ar1_arg_tp;

typedef struct {
	/*
	 * the OU model: dX_t = -phi * X_t-1 + sigma * dW_t. The arguments are phi > 0 and 2*phi/sigma^2 = precision
	 */
	int n;
	double *locations;
	double **log_prec;				       /* theta[0] */
	double **phi_intern;				       /* theta[1] */
} inla_ou_arg_tp;

typedef struct {
	/*
	 * 2D iid random effects. The coding is (x0,y0,x1,y1,...,xn-1,yn-1), so the total length is N=2*n. For the 2DIIDWISHART the coding is (x,y).
	 */
	int n;
	double **log_prec0;				       /* Precision for x */
	double **log_prec1;				       /* Precision for y */
	double **rho_intern;				       /* Corr(x,y) */
} inla_2diid_arg_tp;

typedef struct {
	/*
	 * 3D iid random effects. The coding is (x,y,z), so the total length is N=3*n.
	 */
	int n;
	double **log_prec0;				       /* Precision for x */
	double **log_prec1;				       /* Precision for y */
	double **log_prec2;				       /* Precision for z */
	double **rho_intern01;				       /* Corr(x,y) */
	double **rho_intern02;				       /* Corr(x,z) */
	double **rho_intern12;				       /* Corr(y,z) */
} inla_3diid_arg_tp;

typedef struct {
	/*
	 * iid2d random effects. The coding is (x0,x1,..., y0,y1...)
	 */
	int n;
	int N;
	double **log_prec0;				       /* Precision for x */
	double **log_prec1;				       /* Precision for y */
	double **rho_intern;				       /* Corr(x,y) */
} inla_iid2d_arg_tp;

typedef struct {
	/*
	 * iid3d random effects. The coding is (x0,x1,..., y0,y1,....,z0,z1,...), so the total length is N=3*n.
	 */
	int n;
	int N;
	double **log_prec0;				       /* Precision for x */
	double **log_prec1;				       /* Precision for y */
	double **log_prec2;				       /* Precision for z */
	double **rho_intern01;				       /* Corr(x,y) */
	double **rho_intern02;				       /* Corr(x,z) */
	double **rho_intern12;				       /* Corr(y,z) */
} inla_iid3d_arg_tp;

typedef struct {
	double *vec;
	gsl_matrix *Q;
} inla_wishart_hold_tp;


typedef struct {
	/*
	 * iid random effects with a Wishart prior The coding is (x0,x1,..., y0,y1,....,z0,z1,...), so the total length is N=dim*n.
	 */
	int n;
	int N;
	int dim;
	double ***log_prec;				       /* 0, 1, 2,.... */
	double ***rho_intern;				       /* 01,02,... 12, 13, ... , 23, 24... */

	inla_wishart_hold_tp **hold;
} inla_iid_wishart_arg_tp;

typedef struct {
	double **log_prec;
	int n;						       /* number of elements in F_Z and F_ZADD */
} inla_z_arg_tp;

typedef struct {
	int *array;
	int n;
} inla_iarray_tp;

typedef struct {
	GMRFLib_Qfunc_tp *Qfunc;
	void *Qfunc_arg;
	int n;						       /* size of the original problem */
} inla_replicate_tp;


/* 
   classic me model
 */
typedef struct {
	double **beta;
	double **log_prec_obs;
	double **mean_x;
	double **log_prec_x;
	double *x_obs;
	double *scale;
} inla_mec_tp;

/* 
   berkson me model
 */
typedef struct {
	double **beta;
	double **log_prec;
	double *x;
	double *scale;
} inla_meb_tp;

typedef struct {
	int Id;
	int R2c;
	int c2R;
	int ntheta;
	char *filename_R2c;
	char *filename_c2R;
	double ***theta;
	double **param;
	GMRFLib_tabulate_Qfunc_tp **Q;
} inla_rgeneric_tp;


#define R_GENERIC_Q "Q"
#define R_GENERIC_GRAPH "graph"
#define R_GENERIC_INITIAL "initial"
#define R_GENERIC_LOG_NORM_CONST "log.norm.const"
#define R_GENERIC_LOG_PRIOR "log.prior"
#define R_GENERIC_QUIT "quit"

#define INLA_LITTLE_ENDIAN 1
#define INLA_BIG_ENDIAN    2

/* 
   binary write macros
 */
#define DW(a) {double da = (a); fwrite(&da, sizeof(double), (size_t)1, fp); }
#define D2W(a, b) {DW(a); DW(b);}
#define D3W(a, b, c) {D2W(a, b); DW(c);}
#define D4W(a, b, c, d) {D2W(a, b); D2W(c, d);}
#define IW(a) {double d = (double)(a); DW(d);}		       /* OOPS: write ints as double! */
#define I2W(a, b) {IW(a); IW(b);}
#define I3W(a, b, c) {IW(a); IW(b); IW(c);}
#define I4W(a, b, c, d) {I2W(a, b); I2W(c, d);}
#define I5W(a, b, c, d, e) {I3W(a, b, c); I2W(d, e);}
#define IDW(a, b)  {IW(a); DW(b);}
#define ID2W(a, b, c)  {IW(a); D2W(b, c);}


/* 
   The chose-a-link function
 */
#define CHOSE_LINK(link)						\
	(strcasecmp(link, "identity") == 0 ? link_identity :		\
	 (strcasecmp(link, "log") == 0 ? link_log :			\
	  (strcasecmp(link, "probit") == 0 ? link_probit :		\
	   (strcasecmp(link, "cloglog") == 0 ? link_cloglog :		\
	    (strcasecmp(link, "logit") == 0 ? link_logit :		\
	     (strcasecmp(link, "tan") == 0 ? link_tan :			\
	      (strcasecmp(link, "h") == 0 ? link_h :			\
	       link_this_should_not_happen)))))))

/* 
   functions
 */

GMRFLib_constr_tp *inla_make_constraint(int n, int sumzero, GMRFLib_constr_tp * constr);
GMRFLib_constr_tp *inla_make_constraint2(int n, int replicate, int sumzero, GMRFLib_constr_tp * constr);
GMRFLib_constr_tp *inla_read_constraint(const char *filename, int n);
char *inla_fnmfix(char *name);
char *inla_make_tag(const char *string, int ds);
const char *inla_string_join(const char *a, const char *b);
double Qfunc_rgeneric(int i, int j, void *arg);
double Qfunc_2diid(int i, int j, void *arg);
double Qfunc_2diid_wishart(int i, int j, void *arg);
double Qfunc_ar1(int i, int j, void *arg);
double Qfunc_besag(int i, int j, void *arg);
double Qfunc_besag2(int i, int j, void *arg);
double Qfunc_bym(int i, int j, void *arg);
double Qfunc_copy_part00(int i, int j, void *arg);
double Qfunc_copy_part01(int i, int j, void *arg);
double Qfunc_copy_part11(int i, int j, void *arg);
double Qfunc_generic1(int i, int j, void *arg);
double Qfunc_generic2(int i, int j, void *arg);
double Qfunc_group(int i, int j, void *arg);
double Qfunc_ou(int i, int j, void *arg);
double Qfunc_replicate(int i, int j, void *arg);
double Qfunc_z(int i, int j, void *arg);
double ddexp_taylor(double x, double x0, int order);
double dexp_taylor(double x, double x0, int order);
double exp_taylor(double x, double x0, int order);
double inla_ar1_cyclic_logdet(int N_orig, double phi);
double extra(double *theta, int ntheta, void *argument);
double inla_compute_initial_value(int idx, GMRFLib_logl_tp * logl, double *x_vec, void *arg);
double laplace_likelihood_normalising_constant(double alpha, double gamma, double tau);
double link_cloglog(double x, map_arg_tp typ, void *param);
double link_identity(double x, map_arg_tp typ, void *param);
double link_log(double x, map_arg_tp typ, void *param);
double link_logit(double x, map_arg_tp typ, void *param);
double link_probit(double x, map_arg_tp typ, void *param);
double link_tan(double x, map_arg_tp typ, void *param);
double link_h(double x, map_arg_tp typ, void *param);
double link_this_should_not_happen(double x, map_arg_tp typ, void *param);
double log_apbex(double a, double b);
double map_1exp(double arg, map_arg_tp typ, void *param);
double map_alpha_loglogistic(double arg, map_arg_tp typ, void *param);
double map_alpha_weibull(double arg, map_arg_tp typ, void *param);
double map_alpha_weibull_cure(double arg, map_arg_tp typ, void *param);
double map_beta(double arg, map_arg_tp typ, void *param);
double map_binomialtest_psi(double x, map_arg_tp typ, void *param);
double map_dof(double arg, map_arg_tp typ, void *param);
double map_dof5(double arg, map_arg_tp typ, void *param);
double map_exp(double arg, map_arg_tp typ, void *param);
double map_group_rho(double x, map_arg_tp typ, void *param);
double map_identity(double arg, map_arg_tp typ, void *param);
double map_identity_scale(double arg, map_arg_tp typ, void *param);
double map_invcloglog(double arg, map_arg_tp typ, void *param);
double map_invlogit(double x, map_arg_tp typ, void *param);
double map_invprobit(double arg, map_arg_tp typ, void *param);
double map_invtan(double arg, map_arg_tp typ, void *param);
double map_p_weibull_cure(double arg, map_arg_tp typ, void *param);
double map_phi(double arg, map_arg_tp typ, void *param);
double map_precision(double arg, map_arg_tp typ, void *param);
double map_probability(double x, map_arg_tp typ, void *param);
double map_range(double arg, map_arg_tp typ, void *param);
double map_rho(double arg, map_arg_tp typ, void *param);
double map_shape_svnig(double arg, map_arg_tp typ, void *param);
double map_sqrt1exp(double arg, map_arg_tp typ, void *param);
double map_tau_laplace(double arg, map_arg_tp typ, void *param);
double priorfunc_beta(double *x, double *parameters);
double priorfunc_betacorrelation(double *x, double *parameters);
double priorfunc_logitbeta(double *x, double *parameters);
double priorfunc_bymjoint(double *logprec_besag, double *p_besag, double *logprec_iid, double *p_iid);
double priorfunc_flat(double *x, double *parameters);
double priorfunc_gamma(double *precision, double *parameters);
double priorfunc_gaussian(double *x, double *parameters);
double priorfunc_jeffreys_df_student_t(double *x, double *parameters);
double priorfunc_logflat(double *x, double *parameters);
double priorfunc_loggamma(double *x, double *parameters);
double priorfunc_logiflat(double *x, double *parameters);
double priorfunc_minuslogsqrtruncnormal(double *x, double *parameters);
double priorfunc_normal(double *x, double *parameters);
double priorfunc_wishart1d(double *x, double *parameters);
double priorfunc_wishart2d(double *x, double *parameters);
double priorfunc_wishart3d(double *x, double *parameters);
inla_iarray_tp *find_all_f(inla_tp * mb, inla_component_tp id);
inla_tp *inla_build(const char *dict_filename, int verbose, int make_dir);
int count_f(inla_tp * mb, inla_component_tp id);
int find_f(inla_tp * mb, inla_component_tp id);
int find_tag(inla_tp * mb, const char *name);
int fixup_zadd(inla_tp * mb);
int inla_INLA(inla_tp * mb);
int inla_MCMC(inla_tp * mb_old, inla_tp * mb_new);
int inla_add_copyof(inla_tp * mb);
int inla_computed(GMRFLib_density_tp ** d, int n);
int inla_divisible(int n, int by);
int inla_endian(void);
int inla_error_field_is_void(const char *funcname, const char *secname, const char *field, const char *value);
int inla_error_file_error(const char *funcname, const char *filename, int n, int element_number, double val);
int inla_error_file_error2(const char *funcname, const char *filename, int n, int element_number, double val, int element_number2, double val2);
int inla_error_file_numelm(const char *funcname, const char *filename, int n, int div);
int inla_error_file_totnumelm(const char *funcname, const char *filename, int n, int total);
int inla_error_general(const char *msg);
int inla_error_general2(const char *msg, const char *msg2);
int inla_error_missing_required_field(const char *funcname, const char *secname, const char *field);
int inla_error_open_file(const char *msg);
int inla_iid3d_adjust(double *rho);
int inla_initial_setup(inla_tp * mb);
int inla_integrate_func(double *d_mean, double *d_stdev, GMRFLib_density_tp * density, map_func_tp * func, void *func_arg);
int inla_is_NAs(int nx, const char *string);
int inla_layout_x(double **x, int *n, double xmin, double xmax, double mean);
int inla_make_2diid_graph(GMRFLib_graph_tp ** graph, inla_2diid_arg_tp * arg);
int inla_make_2diid_wishart_graph(GMRFLib_graph_tp ** graph, inla_2diid_arg_tp * arg);
int inla_make_3diid_graph(GMRFLib_graph_tp ** graph, inla_3diid_arg_tp * arg);
int inla_make_3diid_wishart_graph(GMRFLib_graph_tp ** graph, inla_3diid_arg_tp * arg);
int inla_make_ar1_graph(GMRFLib_graph_tp ** graph, inla_ar1_arg_tp * arg);
int inla_make_besag2_graph(GMRFLib_graph_tp ** graph_out, GMRFLib_graph_tp * graph);
int inla_make_bym_graph(GMRFLib_graph_tp ** new_graph, GMRFLib_graph_tp * graph);
int inla_make_group_graph(GMRFLib_graph_tp ** new_graph, GMRFLib_graph_tp * graph, int ngroup, int type, int cyclic, int order, GMRFLib_graph_tp *group_graph);
int inla_make_iid2d_graph(GMRFLib_graph_tp ** graph, inla_iid2d_arg_tp * arg);
int inla_make_iid3d_graph(GMRFLib_graph_tp ** graph, inla_iid3d_arg_tp * arg);
int inla_make_iid_wishart_graph(GMRFLib_graph_tp ** graph, inla_iid_wishart_arg_tp * arg);
int inla_make_ou_graph(GMRFLib_graph_tp ** graph, inla_ou_arg_tp * arg);
int inla_mkdir(const char *dirname);
int inla_ncpu(void);
int inla_output(inla_tp * mb);
int inla_output_Q(inla_tp * mb, const char *dir, GMRFLib_graph_tp * graph);
int inla_output_graph(inla_tp * mb, const char *dir, GMRFLib_graph_tp * graph);
int inla_output_id_names(const char *dir, const char *sdir, inla_file_contents_tp * fc);
int inla_output_matrix(const char *dir, const char *sdir, const char *filename, int n, double *matrix, int *order);
int inla_output_names(const char *dir, const char *sdir, int n, const char **names, const char *suffix);
int inla_output_detail(const char *dir, GMRFLib_density_tp ** density, GMRFLib_density_tp ** gdensity, double *locations, int n, int nrep, Output_tp * output,
		       const char *sdir, map_func_tp * func, void *func_arg, map_func_tp ** ffunc, const char *tag, const char *modelname, int verbose);
int inla_output_hgid(const char *dir);
int inla_output_detail_cpo(const char *dir, GMRFLib_ai_cpo_tp * cpo, int predictor_n, int verbose);
int inla_output_detail_dic(const char *dir, GMRFLib_ai_dic_tp * dic, int verbose);
int inla_output_detail_mlik(const char *dir, GMRFLib_ai_marginal_likelihood_tp * mlik, int verbose);
int inla_output_detail_neffp(const char *dir, GMRFLib_ai_neffp_tp * neffp, int verbose);
int inla_output_detail_theta(const char *dir, double ***theta, int n_theta);
int inla_output_detail_theta_sha1(unsigned char *sha1_hash, double ***theta, int n_theta);
int inla_output_detail_x(const char *dir, double *x, int n_x);
int inla_output_misc(const char *dir, GMRFLib_ai_misc_output_tp * mo, int ntheta, char **theta_tag, char **from_theta, char **to_theta, double *lc_order,
		     int verbose, inla_tp * mb);
int inla_output_size(const char *dir, const char *sdir, int n, int N, int Ntotal, int ngroup, int nrep);
int inla_output_ok(const char *dir);
int inla_parse_INLA(inla_tp * mb, dictionary * ini, int sec, int make_dir);
int inla_parse_data(inla_tp * mb, dictionary * ini, int sec);
int inla_parse_expert(inla_tp * mb, dictionary * ini, int sec);
int inla_parse_ffield(inla_tp * mb, dictionary * ini, int sec);
int inla_parse_lincomb(inla_tp * mb, dictionary * ini, int sec);
int inla_parse_linear(inla_tp * mb, dictionary * ini, int sec);
int inla_parse_mode(inla_tp * mb, dictionary * ini, int sec);
int inla_parse_output(inla_tp * mb, dictionary * ini, int sec, Output_tp ** out);
int inla_parse_predictor(inla_tp * mb, dictionary * ini, int sec);
int inla_parse_problem(inla_tp * mb, dictionary * ini, int sec, int mkdir);
int inla_print_sha1(FILE * fp, unsigned char *md);
int inla_qinv(const char *filename, const char *outfilename);
int inla_finn(const char *filename);
int inla_qreordering(const char *filename);
int inla_qsample(const char *filename, const char *outfile, const char *nsamples, const char *rngfile,
		 const char *samplefile, const char *bfile, const char *mufile, const char *constr_file);
int inla_qsolve(const char *Qfilename, const char *Afilename, const char *Bfilename, const char *method);
int inla_read_data_all(double **x, int *n, const char *filename);
int inla_read_data_general(double **xx, int **ix, int *nndata, const char *filename, int n, int column, int n_columns, int verbose, double default_value);
int inla_read_data_likelihood(inla_tp * mb, dictionary * ini, int sec);
int inla_read_fileinfo(inla_tp * mb, dictionary * ini, int sec, File_tp * file, const char *FILENAME);
int inla_read_prior(inla_tp * mb, dictionary * ini, int sec, Prior_tp * prior, const char *default_prior);
int inla_read_prior0(inla_tp * mb, dictionary * ini, int sec, Prior_tp * prior, const char *default_prior);
int inla_read_prior1(inla_tp * mb, dictionary * ini, int sec, Prior_tp * prior, const char *default_prior);
int inla_read_prior2(inla_tp * mb, dictionary * ini, int sec, Prior_tp * prior, const char *default_prior);
int inla_read_prior3(inla_tp * mb, dictionary * ini, int sec, Prior_tp * prior, const char *default_prior);
int inla_read_prior4(inla_tp * mb, dictionary * ini, int sec, Prior_tp * prior, const char *default_prior);
int inla_read_prior5(inla_tp * mb, dictionary * ini, int sec, Prior_tp * prior, const char *default_prior);
int inla_read_prior6(inla_tp * mb, dictionary * ini, int sec, Prior_tp * prior, const char *default_prior);
int inla_read_prior_generic(inla_tp * mb, dictionary * ini, int sec, Prior_tp * prior, const char *prior_tag, const char *param_tag,
			    const char *from_theta, const char *to_theta, const char *default_prior);
int inla_read_prior_mix(inla_tp * mb, dictionary * ini, int sec, Prior_tp * prior, const char *default_prior);
int inla_read_prior_group(inla_tp * mb, dictionary * ini, int sec, Prior_tp * prior, const char *default_prior);
int inla_read_prior_group0(inla_tp * mb, dictionary * ini, int sec, Prior_tp * prior, const char *default_prior);
int inla_read_prior_group1(inla_tp * mb, dictionary * ini, int sec, Prior_tp * prior, const char *default_prior);
int inla_read_theta_sha1(unsigned char **sha1_hash, double **theta, int *ntheta);
int inla_read_weightsinfo(inla_tp * mb, dictionary * ini, int sec, File_tp * file);
int inla_replicate_graph(GMRFLib_graph_tp ** g, int replicate);
int inla_setup_ai_par_default(inla_tp * mb);
int inla_sread(void *x, int nx, const char *str, int code);
int inla_sread_colon_ints(int *i, int *j, const char *str);
int inla_sread_doubles(double *x, int nx, const char *str);
int inla_sread_doubles_q(double **x, int *nx, const char *str);
int inla_sread_ints(int *x, int nx, const char *str);
int inla_sread_ints_q(int **x, int *nx, const char *str);
int inla_sread_q(void **x, int *nx, const char *str, int code);
int inla_tolower(char *string);
int inla_trim_family(char *family);
int inla_wishart3d_adjust(double *rho);
int loglikelihood_beta(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_betabinomial(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_binomial(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_binomialtest(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_cbinomial(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_exp(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_gaussian(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_gaussian_window(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_inla(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_laplace(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_logperiodogram(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_negative_binomial(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_poisson(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_gpoisson(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_stochvol(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_stochvol_nig(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_stochvol_t(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_tstrata(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_t(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_weibull(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_weibull_cure(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_logistic(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_loglogistic(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_zeroinflated_betabinomial2(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_zeroinflated_binomial0(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_zeroinflated_binomial1(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_zeroinflated_binomial2(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_zeroinflated_negative_binomial0(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_zeroinflated_negative_binomial1(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_zeroinflated_negative_binomial2(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_zeroinflated_poisson0(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_zeroinflated_poisson1(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_zeroinflated_poisson2(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_zeroinflated_poisson2_OLD(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_lognormal(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_iid_gamma(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_iid_logitbeta(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_zero_n_inflated_binomial2(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_circular_normal(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_wrapped_cauchy(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_test_binomial_1(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int my_setenv(char *str);
int testit(int argc, char **argv);
map_table_tp *mapfunc_find(const char *name);
unsigned char *inla_fp_sha1(FILE * fp);
unsigned char *inla_inifile_sha1(const char *filename);
void inla_signal(int sig);
int inla_read_graph(const char *filename);
int inla_iid_wishart_nparam(int dim);
int inla_iid_wishart_adjust(int dim, double *rho);
double priorfunc_wishart(int dim, double *x, double *parameters);
int loglikelihood_loggamma_frailty(double *logll, double *x, int m, int idx, double *x_vec, void *arg);

double inla_log_Phi(double x);
double inla_Phi(double x);
int loglikelihood_skew_normal(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_gev(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_sas(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
double eval_logsum_safe(double lA, double lB);

double priorfunc_mvnorm(double *x, double *parameters);
double priorfunc_sasprior(double *x, double *parameters);
double priorfunc_wishart1d(double *x, double *parameters);
double priorfunc_wishart2d(double *x, double *parameters);
double priorfunc_wishart3d(double *x, double *parameters);
double priorfunc_wishart4d(double *x, double *parameters);
double priorfunc_wishart5d(double *x, double *parameters);
double priorfunc_wishart_generic(int idim, double *x, double *parameters);
double Qfunc_iid_wishart(int node, int nnode, void *arg);
double Qfunc_besagproper(int i, int j, void *arg);
double Qfunc_iid2d(int i, int j, void *arg);
double Qfunc_mec(int i, int j, void *arg);
double mfunc_mec(int i, void *arg);

inla_file_contents_tp *inla_read_file_contents(const char *filename);
int inla_write_file_contents(const char *filename, inla_file_contents_tp * fc);

double iid_mfunc(int idx, void *arg);

int ar_pacf2phi(int p, double *pacf, double *phi);
int ar_phi2pacf(int p, double *phi, double *pacf);
int ar_test1();
int ar_marginal_distribution(int p, double *pacf, double *prec, double *Q);
double Qfunc_ar(int i, int j, void *arg);
double ar_map_pacf(double arg, map_arg_tp typ, void *param);


/* 
***
*/


typedef struct {
	int binary;					       /* use binary output-files */
	int fast_mode;					       /* avoid detailed calculations but use ok approximations */
	inla_mode_tp mode;				       /* which mode to run in */
	double log_prec_initial;			       /* inititial value for log-precisions */
	double mcmc_scale;				       /* scaling */
	int mcmc_thinning;				       /* thinning */
	int mcmc_niter;					       /* number of iterations: 0 is infinite */
	int reorder;					       /* reorder strategy: -1 for optimize */
	int mcmc_fifo;					       /* use fifo to communicate in mcmc mode */
	int mcmc_fifo_pass_data;			       /* use fifo to communicate in mcmc mode, pass also all data */
} G_tp;


#define HYPER_NEW2(name_, initial_, n_)  {				\
		int i_, j_;						\
		name_ = Calloc(n_, double **);				\
		for(j_=0; j_ < n_; j_++){				\
			name_[j_] = Calloc(GMRFLib_MAX_THREADS, double *); \
			for(i_ = 0; i_ < GMRFLib_MAX_THREADS; i_++) { \
				name_[j_][i_] = Calloc(1, double);	\
				name_[j_][i_][0] = initial_;		\
			}						\
		}							\
	}								\

#define HYPER_NEW(name_, initial_)  {					\
		int i_;							\
		name_ = Calloc(GMRFLib_MAX_THREADS, double *);	\
		for(i_ = 0; i_ < GMRFLib_MAX_THREADS; i_++) {	\
			name_[i_] = Calloc(1, double);			\
			name_[i_][0] = initial_;			\
		}							\
	}

#define HYPER_INIT(name_, initial_) {					\
		int i_;							\
		for(i_ = 0; i_ < GMRFLib_MAX_THREADS; i_++) {	\
			name_[i_][0] = initial_;			\
		}							\
	}

#define WISHART_DIM(idx) (mb->f_id[idx] == F_IID1D ? 1 :		\
			  (mb->f_id[idx] == F_IID2D ? 2 :		\
			   (mb->f_id[idx] == F_IID3D ? 3 :		\
			    (mb->f_id[idx] == F_IID4D ? 4 :		\
			     (mb->f_id[idx] == F_IID5D ? 5 : -1)))))


__END_DECLS
#endif
