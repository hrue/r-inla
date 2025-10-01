#ifndef __INLA_SPDE2_H__
#define __INLA_SPDE2_H__
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
#include "GMRFLib/GMRFLib.h"

/* 
   
 */
    typedef enum {
	SPDE2_TRANSFORM_LOGIT = 1,			       /* cos(pi/(1+exp(-x)) */
	SPDE2_TRANSFORM_LOG,				       /* 2*exp(x)-1 */
	SPDE2_TRANSFORM_IDENTITY			       /* x */
} spde2_transform_tp;

typedef struct {
	int i;
	int need_transform;
	double *theta;
	double *vals;

	// need a buffer...
	double l1_cachline[4L];
} spde2_cache_tp;

typedef struct {
	double *V;
	double *v;
} spde2_vV_tp;

typedef struct {
	int n;
	int ntheta;					       /* that is `p' in Finn's notes */
	int ntheta_used;				       /* how many are non-fixed. */
	int debug;
	int theta_first_idx;

	int *fixed;
	double *fixed_values;

	spde2_transform_tp transform;

	GMRFLib_matrix_tp **B;
	GMRFLib_matrix_tp **M;
	GMRFLib_matrix_tp *BLC;

	double ***theta;

	GMRFLib_Qfunc_tp *Qfunc;
	void *Qfunc_arg;
	GMRFLib_graph_tp *graph;

	GMRFLib_vmatrix_tp *Vmatrix;
	double **row_V;
	double **row_v;
} inla_spde2_tp;


double inla_spde2_Qfunction(int thread_id, int ii, int jj, double *values, void *arg);
double inla_spde2_Qfunction_ij_opt(int thread_id, int ii, int jj, double *UNUSED(values), void *arg);
void apply_exponentials(double *__restrict dij, int nb);
void apply_single_transform(int transform, double *d2);
void apply_transform_vectorized(int transform, double *__restrict dij, int nb);
void build_theta_vector(double *__restrict theta, int nc, double ***model_theta, int thread_id);
void compute_d_values_opt(double *__restrict d, double *__restrict vals, double *__restrict theta, int nc, int nc2, int use_ddot_lim);
void compute_diagonal_values(double *__restrict dij, double *__restrict v, double *__restrict values, int nb);
void perform_matrix_vector_mult(double *__restrict V, double *__restrict theta, double *__restrict dij, int nc, int n);

double inla_spde2_Qfunction___ORIG(int thread_id, int ii, int jj, double *values, void *arg);
double inla_spde2_Qfunction_ij(int thread_id, int ii, int jj, double *values, void *arg);
double *inla_spde2_userfunc2(int number, double *theta, int nhyper, double *covmat, void *arg);
int inla_spde2_build_model(int thread_id, inla_spde2_tp ** smodel, const char *prefix, const char *transform);

__END_DECLS
#endif
