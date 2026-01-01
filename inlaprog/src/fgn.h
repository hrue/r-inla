#ifndef __INLA_FGN_H__
#       define __INLA_FGN_H__
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

/* 
 *
 */
#       include "GMRFLib/GMRFLib.h"
#       include "GMRFLib/GMRFLibP.h"
#       include "inla.h"
#       define FGN_KMAX (4L)				       /* maximum K in the tables */
//
    typedef struct {
	int N;						       /* total size = (k+1)*n */
	int n;						       /* size of each component */
	int k;						       /* number of components */
	double prec_eps;				       /* fixed tiny noise */
	double **log_prec;				       /* theta[0] */
	double **H_intern;				       /* theta[1] */
} inla_fgn_arg_tp;

typedef struct {
	int N;						       /* total size = k*n */
	int n;						       /* size of each component */
	int k;						       /* number of components */
	double **log_prec;				       /* theta[0] */
	double **H_intern;				       /* theta[1] */
} inla_fgn2_arg_tp;

double inla_fgn2_helper(int i, int j, int n, double phi);

int inla_make_fgn_graph(GMRFLib_graph_tp ** graph, inla_fgn_arg_tp * def);
int inla_make_fgn2_graph(GMRFLib_graph_tp ** graph, inla_fgn2_arg_tp * def);
double Qfunc_fgn(int thread_id, int i, int j, double *values, void *arg);
double Qfunc_fgn2(int thread_id, int i, int j, double *values, void *arg);

int inla_fgn_get(double *phi, double *w, double H_intern, int k);
double priorfunc_fgn_priorH(double *H_intern, double *param);
void priorfunc_fgn_priorH_extract(void);

__END_DECLS
#endif
