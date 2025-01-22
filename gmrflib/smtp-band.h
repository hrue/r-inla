/*!
  \file smtp-band.h
  \brief GMRFLib interface to the band-solver in the LAPACK library
*/

#ifndef __GMRFLib_SMTP_BAND_H__
#define __GMRFLib_SMTP_BAND_H__

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
void gpskca_(int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *);
void dtbsvspecial_(const char *, const char *, const char *, int *, int *, double *, int *, double *, int *, int *, int *, fortran_charlen_t,
		   fortran_charlen_t, fortran_charlen_t);
int cmsd_(double *, double *, int *, int *, int *, double *, int *, double *);

int GMRFLib_compute_reordering_BAND(int **remap, GMRFLib_graph_tp * graph);
int GMRFLib_build_sparse_matrix_BAND(int thread_id, double **bandmatrix, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg, GMRFLib_graph_tp * graph,
				     int *remap, int bandwidth);
int GMRFLib_factorise_sparse_matrix_BAND(double *band, GMRFLib_fact_info_tp * finfo, GMRFLib_graph_tp * graph, int bandwidth);
int GMRFLib_free_fact_sparse_matrix_BAND(double *bchol);
int GMRFLib_solve_lt_sparse_matrix_BAND(double *rhs, double *bchol, GMRFLib_graph_tp * graph, int *remap, int bandwidth);
int GMRFLib_solve_llt_sparse_matrix_BAND(double *rhs, double *bchol, GMRFLib_graph_tp * graph, int *remap, int bandwidth, double *work);
int GMRFLib_solve_llt_sparse_matrix_special_BAND(double *rhs, double *bchol, GMRFLib_graph_tp * graph, int *remap, int bandwidth, int idx);
int GMRFLib_solve_lt_sparse_matrix_special_BAND(double *rhs, double *bchol, GMRFLib_graph_tp * graph, int *remap, int bandwidth,
						int findx, int toindx, int remapped);
int GMRFLib_solve_l_sparse_matrix_BAND(double *rhs, double *bchol, GMRFLib_graph_tp * graph, int *remap, int bandwidth);
int GMRFLib_solve_l_sparse_matrix_special_BAND(double *rhs, double *bchol, GMRFLib_graph_tp * graph, int *remap, int bandwidth,
					       int findx, int toindx, int remapped);
int GMRFLib_log_determinant_BAND(double *logdet, double *bchol, GMRFLib_graph_tp * graph, int bandwidth);
int GMRFLib_comp_cond_meansd_BAND(double *cmean, double *csd, int indx, double *x, int remapped, double *bchol,
				  GMRFLib_graph_tp * graph, int *remap, int bandwidth);
int GMRFLib_compute_Qinv_BAND(GMRFLib_problem_tp * problem);
int GMRFLib_bitmap_factorisation_BAND__iternal(const char *filename, double *band, GMRFLib_graph_tp * graph, int *remap, int bandwidth);
int GMRFLib_bitmap_factorisation_BAND(const char *filename_body, double *band, GMRFLib_graph_tp * graph, int *remap, int bandwidth);
int GMRFLib_bitmap_factorisation_BAND__intern(const char *filename, double *band, GMRFLib_graph_tp * graph, int *remap, int bandwidth);

__END_DECLS
#endif
