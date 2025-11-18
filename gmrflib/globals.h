#ifndef __GMRFLib_GLOBALS_H__
#define __GMRFLib_GLOBALS_H__

#include <stdlib.h>

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
#if !defined(GMRFLib_FALSE)
#define GMRFLib_FALSE (0)
#endif
#if !defined(GMRFLib_TRUE)
#define GMRFLib_TRUE  (1)
#endif
typedef double GMRFLib_uniform_tp(void);
typedef int GMRFLib_uniform_init_tp(unsigned long int seed);
typedef void *GMRFLib_uniform_getstate_tp(size_t *siz);
typedef int GMRFLib_uniform_setstate_tp(void *state);

typedef double *GMRFLib_ai_INLA_userfunc0_tp(int thread_id, GMRFLib_problem_tp * problem, double *theta, int nhyper);
typedef double *GMRFLib_ai_INLA_userfunc1_tp(int thread_id, double *theta, int nhyper, double *covmat);
typedef double *GMRFLib_ai_INLA_userfunc2_tp(int number, double *theta, int nhyper, double *covmat, void *arg);
typedef double *GMRFLib_ai_INLA_userfunc3_tp(int number, double *theta, int nhyper, double *covmat, void *arg);

#ifndef __GMRFLib_DONT_DEFINE_GLOBALS

extern char GMRFLib_path[];

extern GMRFLib_smtp_tp GMRFLib_smtp;
extern GMRFLib_reorder_tp GMRFLib_reorder;

extern gsl_rng *GMRFLib_rng_ptr;
#pragma omp threadprivate(GMRFLib_rng_ptr)

extern GMRFLib_uniform_tp *GMRFLib_uniform;
extern GMRFLib_uniform_init_tp *GMRFLib_uniform_init;
extern GMRFLib_uniform_getstate_tp *GMRFLib_uniform_getstate;
extern GMRFLib_uniform_setstate_tp *GMRFLib_uniform_setstate;

extern GMRFLib_ai_INLA_userfunc0_tp *GMRFLib_ai_INLA_userfunc0;	/* points to the function */
extern int GMRFLib_ai_INLA_userfunc0_dim;		       /* dimension of func() */
extern GMRFLib_density_tp **GMRFLib_ai_INLA_userfunc0_density; /* return the marginal densities here */

extern GMRFLib_ai_INLA_userfunc1_tp *GMRFLib_ai_INLA_userfunc1;	/* points to the function */
extern int GMRFLib_ai_INLA_userfunc1_dim;		       /* dimension of func() */
extern GMRFLib_density_tp **GMRFLib_ai_INLA_userfunc1_density; /* return the marginal densities here */

extern GMRFLib_ai_INLA_userfunc2_tp **GMRFLib_ai_INLA_userfunc2;
extern void **GMRFLib_ai_INLA_userfunc2_args;
extern GMRFLib_density_tp ***GMRFLib_ai_INLA_userfunc2_density;
extern int GMRFLib_ai_INLA_userfunc2_n;
extern int *GMRFLib_ai_INLA_userfunc2_len;
extern char **GMRFLib_ai_INLA_userfunc2_tag;

extern GMRFLib_ai_INLA_userfunc3_tp **GMRFLib_ai_INLA_userfunc3;
extern void **GMRFLib_ai_INLA_userfunc3_args;
extern GMRFLib_density_tp ***GMRFLib_ai_INLA_userfunc3_density;
extern int GMRFLib_ai_INLA_userfunc3_n;
extern int *GMRFLib_ai_INLA_userfunc3_len;
extern char **GMRFLib_ai_INLA_userfunc3_tag;

extern int GMRFLib_bitmap_max_dimension;
extern int GMRFLib_bitmap_swap;
extern GMRFLib_openmp_tp *GMRFLib_openmp;
extern GMRFLib_global_node_tp GMRFLib_global_node;
extern GMRFLib_density_storage_strategy_tp GMRFLib_density_storage_strategy;
extern int GMRFLib_pardiso_ok;

extern int GMRFLib_faster_constr;
extern double GMRFLib_aqat_m_diag_add;
extern int GMRFLib_inla_mode;
extern int GMRFLib_Qx_strategy;				       // 0 = serial, 1 = parallel
extern int GMRFLib_preopt_predictor_strategy;		       // 0 = !data_rich, 1 = data_rich
extern double GMRFLib_weight_prob;
extern double GMRFLib_weight_prob_one;
extern double **GMRFLib_dot_product_optim_report;
extern int GMRFLib_internal_opt;
extern int GMRFLib_save_memory;

extern int GMRFLib_sort2_id_cut_off;
extern int GMRFLib_sort2_dd_cut_off;

extern int GMRFLib_write_state;
extern int GMRFLib_gaussian_data;

extern int GMRFLib_testit_mode;
extern int GMRFLib_testit_debug;
extern int GMRFLib_opt_solve;
extern int GMRFLib_opt_num_threads;

extern int GMRFLib_intern_flag;
extern int GMRFLib_cachelinesize;

extern int GMRFLib_model_idx;
extern int GMRFLib_model_n;

extern int GMRFLib_force_stiles;

extern char *GMRFLib_tmpdir;

extern double GMRFLib_overall_cpu[8];

#endif
__END_DECLS
#endif
