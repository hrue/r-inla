#ifndef __GMRFLib_HGMRFM_H__
#define __GMRFLib_HGMRFM_H__

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
	GMRFLib_HGMRFM_TP_ETA = 1,
	GMRFLib_HGMRFM_TP_F,
	GMRFLib_HGMRFM_TP_BETA,
	GMRFLib_HGMRFM_TP_LC,
	GMRFLib_HGMRFM_TP___VOID = -1
} GMRFLib_hgmrfm_type_types_tp;

typedef struct {
	GMRFLib_hgmrfm_type_types_tp tp;
	int idx;
	int tp_idx;
} GMRFLib_hgmrfm_type_tp;

typedef struct {
	int n;						       /* length of the linear predictor */
	int n_ext;					       /* length of the extended part of the linear predictor */
	int N;						       /* grand total */

	int *idx_map_f;
	int *idx_map_beta;
	int *idx_map_lc;

	int nf;
	GMRFLib_Qfunc_tp **f_Qfunc;
	void **f_Qfunc_arg;
	GMRFLib_graph_tp **f_graph;

	GMRFLib_Qfunc_tp ***ff_Qfunc;			       /* interaction */
	void ***ff_Qfunc_arg;

	int nbeta;
	double **covariate;
	double *prior_precision;

	GMRFLib_graph_tp *eta_graph;
	GMRFLib_tabulate_Qfunc_tp *eta_Q;

	GMRFLib_graph_tp *eta_ext_graph;
	GMRFLib_tabulate_Qfunc_tp *eta_ext_Q;

	GMRFLib_graph_tp *lc_graph;
	GMRFLib_tabulate_Qfunc_tp *lc_Q;

	GMRFLib_hgmrfm_type_tp *what_type;
} GMRFLib_hgmrfm_arg_tp;

typedef struct {

	/**
	 * The graph for the hgmrf-model
	 */
	GMRFLib_graph_tp *graph;

	/**
	 * The Qfunction for the hgmrf-model
	 */
	GMRFLib_Qfunc_tp *Qfunc;

	/**
	 * The arguments to GMRFLib_hgmrfm_tp::Qfunc
	 */
	void *Qfunc_arg;

	/**
	 * Linear constraints for the hgmrf-model (if any).
	 */
	GMRFLib_constr_tp *constr;
} GMRFLib_hgmrfm_tp;

int GMRFLib_init_hgmrfm(GMRFLib_hgmrfm_tp ** hgmrfm, int n, int n_ext,
			int *eta_sumzero, double **logprec_unstruct_omp,
			const char *Aext_fnm, double Aext_precision,
			int nf, int **c, double **w,
			GMRFLib_graph_tp ** f_graph, GMRFLib_Qfunc_tp ** f_Qfunc,
			void **f_Qfunc_arg, char *f_sumzero, GMRFLib_constr_tp ** f_constr,
			GMRFLib_Qfunc_tp *** ff_Qfunc, void ***ff_Qfunc_arg,
			int nbeta, double **covariate, double *prior_precision, GMRFLib_ai_param_tp * ai_par);
GMRFLib_hgmrfm_type_tp GMRFLib_hgmrfm_what_type(int node, GMRFLib_hgmrfm_arg_tp * a);
double GMRFLib_hgmrfm_Qfunc(int thread_id, int node, int nnode, double *values, void *arg);
int GMRFLib_free_hgmrfm(GMRFLib_hgmrfm_tp * hgmrfm);

__END_DECLS
#endif
