#ifndef __INLA_SPDE3_H__
#       define __INLA_SPDE3_H__
#       undef __BEGIN_DECLS
#       undef __END_DECLS
#       ifdef __cplusplus
#              define __BEGIN_DECLS extern "C" {
#              define __END_DECLS }
#       else
#              define __BEGIN_DECLS			       /* empty */
#              define __END_DECLS			       /* empty */
#       endif

__BEGIN_DECLS
#       include "GMRFLib/GMRFLib.h"

/* 
   
 */
    typedef enum {
	SPDE3_TRANSFORM_IDENTITY = 1,			       /* x */
	SPDE3_TRANSFORM_LOG,				       /* exp(x) */
	SPDE3_TRANSFORM_SHIFTEDLOG,			       /* 2*exp(x)-1 */
	SPDE3_TRANSFORM_LOGIT,				       /* (1-exp(x))/(1+exp(x)) */
	SPDE3_TRANSFORM_OLDLOGIT			       /* cos(pi/(1+exp(-x))) */
} spde3_transform_tp;

typedef struct {
	double *theta;
	double *d3;
} inla_spde3_d3store_tp;

typedef struct {
	int n;
	int n3;						       /* for B[3] and M[3] */
	int ntheta;					       /* that is `p' in Finn's notes */
	int debug;
	int theta_first_idx;

	spde3_transform_tp transform;

	GMRFLib_matrix_tp **B;
	GMRFLib_matrix_tp **M;
	GMRFLib_matrix_tp *BLC;
	GMRFLib_matrix_tp *M3transpose;			       /* the transpose of M3 */

	inla_spde3_d3store_tp **store;

	double ***theta;

	GMRFLib_Qfunc_tp *Qfunc;
	void *Qfunc_arg;
	GMRFLib_graph_tp *graph;
} inla_spde3_tp;

double inla_spde3_Qfunction(int thread_id, int node, int nnode, double *values, void *arg);
int inla_spde3_build_model(int thread_id, inla_spde3_tp ** smodel, const char *prefix, const char *transform);
double *inla_spde3_userfunc3(int number, double *theta, int nhyper, double *covmat, void *arg);

__END_DECLS
#endif
