
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
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS					       /* empty */
# define __END_DECLS					       /* empty */
#endif
__BEGIN_DECLS
#include "iniparser.h"
#include "dictionary.h"
#include "strlib.h"
#define INLA_FAIL  1
#define INLA_OK    0
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
	double *d;					       /* the d-array */
	int ndata;					       /* length of data (from file) */
	double *y;					       /* general responce */

	/*
	 * y ~ Poisson(E*exp(x)) 
	 */
	double *E;

	/*
	 * y ~ Poisson(theta*E1 + theta*E2 + E*exp(x)) 
	 */
	double *E1;
	double *E2;
	double **log_theta_E1;
	double **log_theta_E2;

	/*
	 * y ~ Binomial(nb, p(x))
	 */
	double *nb;

	/*
	 * y ~ Neg.Binomial(n, p(x)), n=size is a hyperparameter (overdispersion)
	 */
	double **log_size;

	/*
	 * y ~ Normal(x, 1/(weight*prec)) 
	 */
	double **log_prec_gaussian;
	double *weight_gaussian;			       /* weights for the gaussian: Variance = 1/(weight*prec) */
	/*
	 * y ~ T_dof(x, 1/(weight*prec)), where T_dof has Variance=1
	 */
	double **log_prec_t;
	double *weight_t;				       /* weights for the t: Variance = 1/(weight*prec) */
	double **dof_intern_t;
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
	double **alpha_intern;				       /* For the Weibull and PS */
	double *truncation;				       /* for survival */
	double *event;					       /* one of inla_surv_event_tp */
	double *lower;					       /* for survival */
	double *upper;					       /* for survival */
	// dobule *time == y
	double **p_intern;				       /* For the L_WEIBULL_CURE */

	/*
	 * zero-inflated Poission/Binomial/nbinomial version 0 and 1: prob
	 */
	double **prob_intern;
	double **zeroinflated_alpha_intern;		       /* alpha = exp(alpha_intern) */

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

} Data_tp;

/* 
   this is needed so we can identify each component in the model
 */
typedef enum {
	INVALID_COMPONENT = 0,
	L_GAUSSIAN,					       /* likelihood-models */
	L_SKEWNORMAL,
	L_T,
	L_POISSON,
	L_POISSONEXT,
	L_BINOMIAL,
	L_ZEROINFLATEDBINOMIAL0,
	L_ZEROINFLATEDBINOMIAL1,
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
	L_ZEROINFLATEDPOISSON0,
	L_ZEROINFLATEDPOISSON1,
	L_ZEROINFLATEDPOISSON2,
	L_WEIBULL_CURE,					       /* Patrick and Silvia's model */
	L_LAPLACE,
	F_RW2D,						       /* f-models */
	F_BESAG,
	F_BESAG2,					       /* the [a*x, x/a] model */
	F_BESAGMOD,
	F_SEASONAL,
	F_IID,
	F_2DIID,
	F_IID1D,
	F_IID2D,
	F_IID3D,
	F_RW1,
	F_RW2,
	F_CRW2,
	F_AR1,
	F_Z,
	F_ZADD,
	F_BYM,
	F_GENERIC0,
	F_GENERIC1,
	F_GENERIC2,
	F_MATERN2D,
	F_SPHERE,
	F_SPDE,
	F_COPY,
	P_LOGGAMMA,					       /* priors */
	P_GAUSSIAN,
	P_MINUSLOGSQRTRUNCGAUSSIAN,
	P_FLAT,
	P_BYMJOINT,
	P_WISHART1D,
	P_WISHART2D,
	P_WISHART3D,
	P_LOGFLAT, 
	O_POSITIVE /* offset-functions */ ,
	G_EXCHANGEABLE,					       /* group models */
	G_AR1
} inla_component_tp;


/* 
   priors are defined using this template. return log(pi(precision, parameters....))
 */
typedef double inla_priorfunc_tp(double *param, double *parameters);
typedef double inla_priorfunc2_tp(double *param0, double *parameters0, double *param1, double *parameters1);

typedef struct {
	inla_component_tp id;				       /* prior Id */
	char *name;					       /* name of prior */
	double *parameters;				       /* the parameters */
	inla_priorfunc_tp *priorfunc;			       /* priorfunction */
	inla_priorfunc2_tp *priorfunc2;			       /* priorfunction2 */
} Prior_tp;


typedef struct {
	GMRFLib_tabulate_Qfunc_tp *tab;
	double **beta;
	double **log_prec;

	int n;
	double *eigenvalues;
	double max_eigenvalue;
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
	int density;					       /* output detailed marginal density */
	int hyperparameters;				       /* compute also the marginal for the hyperparameters */
	int kld;					       /* output the (symmetric) kld between marginals */
	int mlik;					       /* compute the marginal likelihood? */
	int q;						       /* output image of the Q-matrix */
	int nquantiles;					       /* compute cdfs and/or quantiles; max 10 */
	int ncdf;
	double *quantiles;
	double *cdf;
} Output_tp;

typedef struct inla_tp_struct inla_tp;			       /* need it like this as they point to each other */

typedef struct {
	char *data_likelihood;
	inla_component_tp data_id;
	File_tp data_file;
	Prior_tp data_prior;
	Prior_tp data_prior0;
	Prior_tp data_prior1;
	Data_tp data_observations;
	int data_fixed;
	int data_fixed0;
	int data_fixed1;
	int data_ntheta;
	GMRFLib_logl_tp *loglikelihood;
	double *offset;
	map_func_tp *predictor_linkfunc;
	inla_tp *mb;					       /* to get the off_.... */
} Data_section_tp;


typedef struct {
	double **log_positive;
	double *weights;
	int n;
} inla_off_func_arg_tp;

typedef double inla_off_func_tp(int idx, void *arg);

typedef struct {
	int n;
	GMRFLib_Qfunc_tp *Qfunc;
	void *Qfunc_arg;
	double **rho_intern;
} inla_rwc2_def_tp;

typedef struct {
	int N;
	int ngroup;
	int type;
	GMRFLib_Qfunc_tp *Qfunc;
	void *Qfunc_arg;
	double **group_rho_intern;
} inla_group_def_tp;

typedef struct {
	double precision;
	double **beta;

	GMRFLib_Qfunc_tp *Qfunc;
	void *Qfunc_arg;
} inla_copy_arg_tp;

struct inla_tp_struct {
	/*
	 * General stuff 
	 */
	int verbose;

	/*
	 * reuse the mode-stuff 
	 */
	unsigned char *sha1_hash;
	unsigned char *sha1_hash_file;
	int ntheta_file;
	int theta_counter_file;
	int reuse_mode;
	int reuse_mode_but_restart;
	double *theta_file;
	double *x_file;
	int nx_file;

	/*
	 * Expert options 
	 */
	int expert_cpo_manual;
	int expert_cpo_idx;
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
	char *predictor_tag;				       /* the tag */
	char *predictor_dir;				       /* the directory */
	double **predictor_log_prec;
	Prior_tp predictor_prior;
	int *predictor_cross_sumzero;
	int predictor_compute;
	int predictor_fixed;
	int predictor_user_scale;
	map_func_tp **predictor_linkfunc;		       /* these are rebuilt */
	map_table_tp *predictor_usermap;
	Output_tp *predictor_output;
	double *offset;					       /* the offset y ~ f(eta + offset) */

	int predictor_n_ext;				       /* dimension of \eta_ext */
	char *predictor_Aext_fnm;			       /* extension: filename for the Amatrix  */
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
	double **f_locations;
	double **f_weights;
	GMRFLib_graph_tp **f_graph;
	GMRFLib_graph_tp **f_graph_orig;
	GMRFLib_Qfunc_tp **f_Qfunc;
	GMRFLib_Qfunc_tp **f_Qfunc_orig;
	void **f_Qfunc_arg;
	void **f_Qfunc_arg_orig;
	Prior_tp **f_prior;
	char *f_sumzero;
	GMRFLib_constr_tp **f_constr;
	GMRFLib_constr_tp **f_constr_orig;
	double *f_diag;
	double *f_rankdef;
	double ****f_theta;
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
	double *f_precision;
	Output_tp **f_output;
	map_table_tp **f_usermap;

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
	map_table_tp **linear_usermap;

	/*
	 * linear combinations 
	 */
	int nlc;					       /* number of linear combinations */
	int lc_derived_only;				       /* use only the derived lincombs ? */
	char **lc_tag;					       /* the tags */
	double **lc_w;					       /* their weights */
	double *lc_prec;				       /* the `high' precision */
	char **lc_dir;
	Output_tp **lc_output;
	map_table_tp **lc_usermap;
	double *Alc;

	/*
	 * offset functions
	 */
	Output_tp **off_output;
	Prior_tp **off_prior;
	char **off_dir;
	char **off_tag;
	double ****off_theta;
	double **off_weights;
	inla_off_func_tp **off_func;
	int **off_fixed;
	int *off_ntheta;
	int noff;
	map_table_tp **off_usermap;
	void **off_func_arg;

	/*
	 * The final model 
	 */
	GMRFLib_hgmrfm_tp *hgmrfm;
	int ntheta;
	double ***theta;
	char **theta_tag;
	char **theta_tag_userscale;
	char **theta_dir;
	map_func_tp **theta_map;
	void **theta_map_arg;
	map_table_tp **theta_usermap;
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


#define INLA_LITTLE_ENDIAN 1
#define INLA_BIG_ENDIAN    2

/* 
   binary write macros
 */
#define DW(a) {double da = (a); size_t ret;  ret = fwrite(&da, sizeof(double), (size_t)1, fp); }
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
   functions
 */

GMRFLib_constr_tp *inla_make_constraint(int n, int sumzero, GMRFLib_constr_tp * constr);
GMRFLib_constr_tp *inla_make_constraint2(int n, int replicate, int sumzero, GMRFLib_constr_tp * constr);
GMRFLib_constr_tp *inla_read_constraint(const char *filename, int n);
char *inla_fnmfix(char *name);
char *inla_make_tag(const char *string, int ds);
const char *inla_string_join(const char *a, const char *b);
double Qfunc_2diid(int i, int j, void *arg);
double Qfunc_ar1(int i, int j, void *arg);
double Qfunc_besag(int i, int j, void *arg);
double Qfunc_besag2(int i, int j, void *arg);
double Qfunc_besagmod(int i, int j, void *arg);
double Qfunc_bym(int i, int j, void *arg);
double Qfunc_copy_part00(int i, int j, void *arg);
double Qfunc_copy_part01(int i, int j, void *arg);
double Qfunc_copy_part11(int i, int j, void *arg);
double Qfunc_generic1(int i, int j, void *arg);
double Qfunc_generic2(int i, int j, void *arg);
double Qfunc_group(int i, int j, void *arg);
double Qfunc_iid2d(int i, int j, void *arg);
double Qfunc_iid3d(int i, int j, void *arg);
double Qfunc_replicate(int i, int j, void *arg);
double Qfunc_z(int i, int j, void *arg);
double ddexp_taylor(double x, double x0, int order);
double dexp_taylor(double x, double x0, int order);
double exp_taylor(double x, double x0, int order);
double extra(double *theta, int ntheta, void *argument);
double inla_all_offset(int idx, void *arg);
double inla_compute_initial_value(int idx, GMRFLib_logl_tp * logl, double *x_vec, void *arg);
double laplace_likelihood_normalising_constant(double alpha, double gamma, double tau);
double link_log(double x, map_arg_tp typ, void *param);
double link_logit(double x, map_arg_tp typ, void *param);
double log_apbex(double a, double b);
double map_1exp(double arg, map_arg_tp typ, void *param);
double map_alpha_weibull(double arg, map_arg_tp typ, void *param);
double map_alpha_weibull_cure(double arg, map_arg_tp typ, void *param);
double map_dof(double arg, map_arg_tp typ, void *param);
double map_exp(double arg, map_arg_tp typ, void *param);
double map_group_rho(double x, map_arg_tp typ, void *param);
double map_identity(double arg, map_arg_tp typ, void *param);
double map_invlogit(double x, map_arg_tp typ, void *param);
double map_p_weibull_cure(double arg, map_arg_tp typ, void *param);
double map_phi(double arg, map_arg_tp typ, void *param);
double map_precision(double arg, map_arg_tp typ, void *param);
double map_probability(double x, map_arg_tp typ, void *param);
double map_range(double arg, map_arg_tp typ, void *param);
double map_rho(double arg, map_arg_tp typ, void *param);
double map_shape_svnig(double arg, map_arg_tp typ, void *param);
double map_sqrt1exp(double arg, map_arg_tp typ, void *param);
double map_tau_laplace(double arg, map_arg_tp typ, void *param);
double offset_positive_func(int idx, void *arg);
double priorfunc_beta(double *x, double *parameters);
double priorfunc_bymjoint(double *logprec_besag, double *p_besag, double *logprec_iid, double *p_iid);
double priorfunc_flat(double *x, double *parameters);
double priorfunc_logflat(double *x, double *parameters);
double priorfunc_gamma(double *precision, double *parameters);
double priorfunc_gaussian(double *x, double *parameters);
double priorfunc_loggamma(double *x, double *parameters);
double priorfunc_minuslogsqrtruncnormal(double *x, double *parameters);
double priorfunc_normal(double *x, double *parameters);
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
int inla_make_bym_graph(GMRFLib_graph_tp ** new_graph, GMRFLib_graph_tp * graph);
int inla_make_besag2_graph(GMRFLib_graph_tp ** graph_out, GMRFLib_graph_tp * graph);
int inla_make_group_graph(GMRFLib_graph_tp ** new_graph, GMRFLib_graph_tp * graph, int ngroup, int type);
int inla_make_iid2d_graph(GMRFLib_graph_tp ** graph, inla_iid2d_arg_tp * arg);
int inla_make_iid3d_graph(GMRFLib_graph_tp ** graph, inla_iid3d_arg_tp * arg);
int inla_mkdir(const char *dirname);
int inla_ncpu(void);
int inla_output(inla_tp * mb);
int inla_output_Q(inla_tp * mb, const char *dir, GMRFLib_graph_tp * graph);
int inla_output_detail(const char *dir, GMRFLib_density_tp ** density, GMRFLib_density_tp ** gdensity, double *locations, int n, int nrep, Output_tp * output,
		       const char *sdir, map_func_tp * func, void *func_arg, map_func_tp ** ffunc, const char *tag, const char *modelname, int verbose);
int inla_output_detail_cpo(const char *dir, GMRFLib_ai_cpo_tp * cpo, int predictor_n, int verbose);
int inla_output_detail_dic(const char *dir, GMRFLib_ai_dic_tp * dic, int verbose);
int inla_output_detail_mlik(const char *dir, GMRFLib_ai_marginal_likelihood_tp * mlik, int verbose);
int inla_output_detail_neffp(const char *dir, GMRFLib_ai_neffp_tp * neffp, int verbose);
int inla_output_detail_theta(const char *dir, double ***theta, int n_theta);
int inla_output_detail_theta_sha1(unsigned char *sha1_hash, double ***theta, int n_theta);
int inla_output_detail_x(const char *dir, double *x, int n_x);
int inla_output_misc(const char *dir, GMRFLib_ai_misc_output_tp * mo, int verbose);
int inla_output_size(const char *dir, const char *sdir, int n, int N, int Ntotal, int ngroup, int nrep);
int inla_parse_INLA(inla_tp * mb, dictionary * ini, int sec, int make_dir);
int inla_parse_data(inla_tp * mb, dictionary * ini, int sec);
int inla_parse_expert(inla_tp * mb, dictionary * ini, int sec);
int inla_parse_ffield(inla_tp * mb, dictionary * ini, int sec);
int inla_parse_lincomb(inla_tp * mb, dictionary * ini, int sec);
int inla_parse_linear(inla_tp * mb, dictionary * ini, int sec);
int inla_parse_mode(inla_tp * mb, dictionary * ini, int sec);
int inla_parse_offset(inla_tp * mb, dictionary * ini, int sec);
int inla_parse_output(inla_tp * mb, dictionary * ini, int sec, Output_tp ** out);
int inla_parse_predictor(inla_tp * mb, dictionary * ini, int sec);
int inla_parse_problem(inla_tp * mb, dictionary * ini, int sec, int mkdir);
int inla_print_sha1(FILE * fp, unsigned char *md);
int inla_qinv(const char *filename);
int inla_read_data_all(double **x, int *n, const char *filename);
int inla_read_data_general(double **xx, int **ix, int *nndata, const char *filename, int n, int column, int n_columns, int verbose, double default_value);
int inla_read_data_likelihood(inla_tp * mb, dictionary * ini, int sec);
int inla_read_fileinfo(inla_tp * mb, dictionary * ini, int sec, File_tp * file);
int inla_read_prior(inla_tp * mb, dictionary * ini, int sec, Prior_tp * prior, const char *default_prior);
int inla_read_prior0(inla_tp * mb, dictionary * ini, int sec, Prior_tp * prior, const char *default_prior);
int inla_read_prior1(inla_tp * mb, dictionary * ini, int sec, Prior_tp * prior, const char *default_prior);
int inla_read_prior2(inla_tp * mb, dictionary * ini, int sec, Prior_tp * prior, const char *default_prior);
int inla_read_prior3(inla_tp * mb, dictionary * ini, int sec, Prior_tp * prior, const char *default_prior);
int inla_read_prior4(inla_tp * mb, dictionary * ini, int sec, Prior_tp * prior, const char *default_prior);
int inla_read_prior_generic(inla_tp * mb, dictionary * ini, int sec, Prior_tp * prior, const char *prior_tag, const char *param_tag, const char *default_prior);
int inla_read_prior_group(inla_tp * mb, dictionary * ini, int sec, Prior_tp * prior, const char *default_prior);
int inla_read_theta_sha1(unsigned char **sha1_hash, double **theta, int *ntheta);
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
int loglikelihood_binomial(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_exp(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_gaussian(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_inla(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_laplace(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_logperiodogram(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_negative_binomial(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_poisson(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_poisson_ext(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_stochvol(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_stochvol_nig(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_stochvol_t(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_t(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_weibull(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_weibull_cure(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_zeroinflated_binomial0(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_zeroinflated_binomial1(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_zeroinflated_negative_binomial0(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_zeroinflated_negative_binomial1(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_zeroinflated_negative_binomial2(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_zeroinflated_poisson0(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_zeroinflated_poisson1(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_zeroinflated_poisson2(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int loglikelihood_zeroinflated_poisson2_OLD(double *logll, double *x, int m, int idx, double *x_vec, void *arg);
int my_setenv(char *str);
int testit(void);
map_table_tp *mapfunc_find(const char *name);
unsigned char *inla_fp_sha1(FILE * fp);
unsigned char *inla_inifile_sha1(const char *filename);
void inla_signal(int sig);

double inla_log_Phi(double x);
int loglikelihood_skew_normal(double *logll, double *x, int m, int idx, double *x_vec, void *arg);


/* 
***
*/

typedef double inla_all_offset_func_tp(int idx, void *arg);

typedef struct {
	int tmax;					       /* maximum number of threads */
	int binary;					       /* use binary output-files */
	int fast_mode;					       /* avoid detailed calculations but use ok approximations */
	int hyper_mode;					       /* enable accurate computations of the hyperparameters */
	int mcmc_mode;					       /* enable MCMC mode */
	int qinv_mode;					       /* (expert) just compute marginal variance... */
	double log_prec_initial;			       /* inititial value for log-precisions */
	double dof_max;					       /* max dof for (additive) student-t */
	double mcmc_scale;				       /* scaling */
	int mcmc_thinning;				       /* thinning */
	int mcmc_niter;					       /* number of iterations: 0 is infinite */
	int reorder;					       /* reorder strategy: -1 for optimize */
	inla_all_offset_func_tp *all_offset;
} G_tp;


#define HYPER_NEW2(name_, initial_, n_)  {				\
		int i_, j_;						\
		name_ = Calloc(n_, double **);				\
		for(j_=0; j_ < n_; j_++){				\
			name_[j_] = Calloc(G.tmax, double *);		\
			for(i_ = 0; i_ < G.tmax; i_++) {		\
				name_[j_][i_] = Calloc(1, double);	\
				name_[j_][i_][0] = initial_;		\
			}						\
		}							\
	}								\

#define HYPER_NEW(name_, initial_)  {			\
		int i_;					\
		name_ = Calloc(G.tmax, double *);	\
		for(i_ = 0; i_ < G.tmax; i_++) {	\
			name_[i_] = Calloc(1, double);	\
			name_[i_][0] = initial_;	\
		}					\
	}

#define HYPER_INIT(name_, initial_) { \
		int i_;					\
		for(i_ = 0; i_ < G.tmax; i_++) {	\
			name_[i_][0] = initial_;	\
		}					\
	}


__END_DECLS
#endif
