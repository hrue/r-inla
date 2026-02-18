#include <omp.h>
#include <stdbool.h>
#include <stdlib.h>
#include <strings.h>
#include <time.h>

#ifndef __GMRFLib_SMTP_STILES_H__
#       define __GMRFLib_SMTP_STILES_H__


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
//...
#       include "no-stiles.h"
// ...
    typedef struct {
	int verbose;
	int block_size;
	int param_len;
	int *param;
} GMRFLib_stiles_ctl_tp;

#       if 0
// this one is defined in problem-setup.h as we need it first time there
typedef struct {
	int in_group;
	int within_group;
	int nrhs;
} GMRFLib_stiles_idx_tp;
#       endif

typedef struct {
	int ng;
	int ng2;
	int ngt;
	int n_in_group;
	int *n_within_group;
	int *n_cores_group;
	int nt_outer;
	int nt_inner;
	int nt_special;
	int nt_max_threads;
	int rescale_on;
	int **perm;
	int **iperm;
	int *n;
	int *nnz;
	int nrhss;
	int *rhss;
	bool **Qinv_done;
	bool **bind_done;
	bool **chol_done;
	GMRFLib_ptr_tp *graphs;
	void *obj;
	double wtime;
} GMRFLib_stiles_store_tp;


typedef struct {
	GMRFLib_ptr_tp *graphs;
	GMRFLib_idx_tp *nrhss;
} GMRFLib_stiles_setup_tp;

GMRFLib_stiles_ctl_tp * GMRFLib_stiles_get_ctl(void);
GMRFLib_stiles_setup_tp *GMRFLib_stiles_get_setup(void *mb);
double GMRFLib_stiles_logdet(GMRFLib_stiles_idx_tp * stiles_idx);
int *GMRFLib_stiles_get_iperm(GMRFLib_stiles_idx_tp * stiles_idx);
int *GMRFLib_stiles_get_perm(GMRFLib_stiles_idx_tp * stiles_idx);
int GMRFLib_stiles_Qinv_INLA(GMRFLib_problem_tp * problem);
int GMRFLib_stiles_build(GMRFLib_stiles_idx_tp * stiles_idx, int thread_id, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg);
int GMRFLib_stiles_chol(GMRFLib_stiles_idx_tp * stiles_idx);
int GMRFLib_stiles_get_block_size(void);
int GMRFLib_stiles_get_tile_size(void);
int GMRFLib_stiles_get_verbose();
int GMRFLib_stiles_is_rescale(void);
int GMRFLib_stiles_rescale_group(void);
int GMRFLib_stiles_set_ctl(int verbose, int block_size, int len, int *param);
int GMRFLib_stiles_set_idx(GMRFLib_stiles_idx_tp * stiles_idx, int nrhs);
int GMRFLib_stiles_setup(GMRFLib_stiles_setup_tp * setup);
int GMRFLib_stiles_solve_L(GMRFLib_stiles_idx_tp * stiles_idx, double *rhs);
int GMRFLib_stiles_solve_LLT(GMRFLib_stiles_idx_tp * stiles_idx, double *rhs);
int GMRFLib_stiles_solve_LT(GMRFLib_stiles_idx_tp * stiles_idx, double *rhs);
void *GMRFLib_stiles_get_store_ptr(void);
void GMRFLib_stiles_Qinv(GMRFLib_stiles_idx_tp * stiles_idx);
void GMRFLib_stiles_bind(GMRFLib_stiles_idx_tp * stiles_idx);
void GMRFLib_stiles_free_setup(GMRFLib_stiles_setup_tp * setup);
void GMRFLib_stiles_print(FILE * fp);
void GMRFLib_stiles_print_ctl_param(FILE *fp, char *suf);
void GMRFLib_stiles_print_idx(GMRFLib_stiles_idx_tp * stiles_idx, FILE * fp);
void GMRFLib_stiles_quit(void);
void GMRFLib_stiles_rescale_end(void);
void GMRFLib_stiles_rescale_start(void);
void GMRFLib_stiles_unbind(GMRFLib_stiles_idx_tp * stiles_idx);
void GMRFLib_stiles_unbind_all(void);
void GMRFLib_stiles_unbind_group(int in_group);



// this function is not defined in 'stiles.h'
int get_auto_tile_size(void);

__END_DECLS
#endif
