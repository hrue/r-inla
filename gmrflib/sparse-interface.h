
/*!
  \file sparse-interface.h
  \brief Typedefs and defines for the interface to the sparse-matrix libraries.
*/

#ifndef __GMRFLib_SPARSE_INTERFACE_H__
#define __GMRFLib_SPARSE_INTERFACE_H__

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

/* 
   
 */
// this is from smtp-stiles.h
    typedef struct {
	int in_group;
	int within_group;
	int sidx;
} GMRFLib_stiles_idx_tp;


typedef struct {
	double factor;
	int degree;
} GMRFLib_global_node_tp;

typedef enum {
	GMRFLib_SMTP_INVALID = -1,
	GMRFLib_SMTP_BAND = 1,
	GMRFLib_SMTP_TAUCS = 2,
	GMRFLib_SMTP_PARDISO = 3,
	GMRFLib_SMTP_STILES = 4,
	GMRFLib_SMTP_DEFAULT = 5
} GMRFLib_smtp_tp;

#define GMRFLib_SMTP_NAME(smtp)			     \
	((smtp) == GMRFLib_SMTP_BAND ? "band" :    \
	 ((smtp) == GMRFLib_SMTP_TAUCS ? "taucs" :	  \
	  ((smtp) == GMRFLib_SMTP_PARDISO ? "pardiso" :		\
	   ((smtp) == GMRFLib_SMTP_STILES ? "sTiles" :		\
		    ((smtp) == GMRFLib_SMTP_DEFAULT ? "default" : "THIS SHOULD NOT HAPPEN")))))

typedef enum {

	/**
	 * \brief The default reordering choice, which depends on the sparse-solver used.
	 * 
	 * Currently, it minimise the bandwidth for the band-solver, and nested dissection using the TAUCS-solver 
	 */
	GMRFLib_REORDER_AUTO = -1,
	GMRFLib_REORDER_DEFAULT = 0,
	GMRFLib_REORDER_IDENTITY,
	GMRFLib_REORDER_BAND,
	GMRFLib_REORDER_METIS,
	GMRFLib_REORDER_GENMMD,
	GMRFLib_REORDER_AMD,
	GMRFLib_REORDER_MD,
	GMRFLib_REORDER_MMD,
	GMRFLib_REORDER_AMDBAR,
	GMRFLib_REORDER_AMDC,
	GMRFLib_REORDER_AMDBARC,
	GMRFLib_REORDER_REVERSE_IDENTITY,
	GMRFLib_REORDER_PARDISO,
	GMRFLib_REORDER_STILES
} GMRFLib_reorder_tp;

/*! 
  \struct GMRFLib_fact_info_tp problem-setup.h
  \brief Description of Qmatrix
 */
typedef struct {

	/**
	 *  \brief Size of Q 
	 */
	int n;

	/**
	 *  \brief Number of non-zero elements 
	 */
	int nnzero;

	/**
	 *  \brief Fillin in L 
	 */
	int nfillin;
} GMRFLib_fact_info_tp;

typedef struct {
	int n;
	int nnz;
	int *len;
	int *rowind;
	int *rowind_sorted;
	GMRFLib_idx2_tp *sort2;
} GMRFLib_taucs_cache_tp;


typedef struct {

	/**
	 * \brief Sparse-matrix type.					       
	 */
	GMRFLib_smtp_tp smtp;

	/**
	 *  \brief Permutation 
	 */
	int *remap;

	/**
	 *  \brief The Cholesky factorisation (smtp == BAND)
	 */
	double *bchol;

	/**
	 *  \brief The bandwidth of the Cholesky factorisation (smtp == BAND)
	 */
	int bandwidth;

	/**
	 *  \brief The Cholesky factorisation (smtp == TAUCS)
	 */
	taucs_ccs_matrix *TAUCS_L;
	taucs_crs_matrix *TAUCS_LL;

	/**
	 *  \brief The symbolic factorisation (smtp == TAUCS)
	 */
	supernodal_factor_matrix *TAUCS_symb_fact;

	GMRFLib_taucs_cache_tp *TAUCS_cache;

	 /**
	 *  \brief Info about the factorization 
	 */
	GMRFLib_fact_info_tp finfo;

	 /**
	 *  \brief The factorisation of PARDISO
	 */
	GMRFLib_pardiso_store_tp *PARDISO_fact;

} GMRFLib_sm_fact_tp;

/* 
   
 */

const char *GMRFLib_reorder_name(GMRFLib_reorder_tp r);
int GMRFLib_bitmap_factorisation(const char *filename_body, GMRFLib_sm_fact_tp * sm_fact, GMRFLib_graph_tp * graph);
int GMRFLib_build_sparse_matrix(int thread_id, GMRFLib_sm_fact_tp * sm_fact, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg,
				GMRFLib_graph_tp * graph, GMRFLib_problem_tp * problem);
int GMRFLib_comp_cond_meansd(double *cmean, double *csd, int indx, double *x, int remapped, GMRFLib_sm_fact_tp * sm_fact, GMRFLib_graph_tp * graph);
int GMRFLib_compute_Qinv(void *problem);
int GMRFLib_compute_reordering(GMRFLib_sm_fact_tp * sm_fact, GMRFLib_graph_tp * graph, GMRFLib_global_node_tp * gn);
int GMRFLib_factorise_sparse_matrix(GMRFLib_sm_fact_tp * sm_fact, GMRFLib_graph_tp * graph, GMRFLib_problem_tp * problem);
int GMRFLib_free_fact_sparse_matrix(GMRFLib_sm_fact_tp * sm_fact);
int GMRFLib_free_reordering(GMRFLib_sm_fact_tp * sm_fact);
int GMRFLib_log_determinant(double *logdet, GMRFLib_sm_fact_tp * sm_fact, GMRFLib_graph_tp * graph, GMRFLib_problem_tp * problem);
int GMRFLib_reorder_id(const char *name);
int GMRFLib_solve_l_sparse_matrix(double *rhs, int nrhs, GMRFLib_sm_fact_tp * sm_fact, GMRFLib_graph_tp * graph, GMRFLib_problem_tp * problem);
int GMRFLib_solve_l_sparse_matrix_special(double *rhs, GMRFLib_sm_fact_tp * sm_fact, GMRFLib_graph_tp * graph, int findx,
					  int toindx, int remapped, GMRFLib_problem_tp * problem);
int GMRFLib_solve_llt_sparse_matrix(double *rhs, int nrhs, GMRFLib_sm_fact_tp * sm_fact, GMRFLib_graph_tp * graph,
				    GMRFLib_problem_tp * problem, GMRFLib_stiles_idx_tp * stiles_idx);
int GMRFLib_solve_llt_sparse_matrix_special(double *rhs, GMRFLib_sm_fact_tp * sm_fact, GMRFLib_graph_tp * graph, int idx,
					    GMRFLib_problem_tp * problem);
int GMRFLib_solve_lt_sparse_matrix(double *rhs, int nrhs, GMRFLib_sm_fact_tp * sm_fact, GMRFLib_graph_tp * graph, GMRFLib_problem_tp * problem);
int GMRFLib_solve_lt_sparse_matrix_special(double *rhs, GMRFLib_sm_fact_tp * sm_fact, GMRFLib_graph_tp * graph, int findx,
					   int toindx, int remapped, GMRFLib_problem_tp * problem);
int GMRFLib_valid_smtp(int smtp);

__END_DECLS
#endif
