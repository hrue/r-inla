
/*!
  \file taucs_privateP.h
  \brief Private typedefs and defines for the TAUCS library version 2.0

  This file is included to make GMRFLib compile without refering to the TAUCS library,  
*/

#ifndef __GMRFLib_TAUCS_PRIVATEP_H__
#define __GMRFLib_TAUCS_PRIVATEP_H__

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
//
extern double taucs_dtl(zero_const);
extern double taucs_dtl(one_const);

double taucs_dtl(add_fn) (double a, double b);
double taucs_dtl(sub_fn) (double a, double b);
double taucs_dtl(mul_fn) (double a, double b);
double taucs_dtl(div_fn) (double a, double b);
double taucs_dtl(neg_fn) (double a);
double taucs_dtl(sqrt_fn) (double a);
double taucs_dtl(conj_fn) (double a);
double taucs_dtl(abs_fn) (double a);

taucs_ccs_matrix *taucs_dtl(ccs_create) (int m, int n, int nnz);
taucs_ccs_matrix *taucs_ccs_create(int m, int n, int nnz, int flags);
void taucs_dtl(ccs_free) (taucs_ccs_matrix * matrix);
void taucs_ccs_free(taucs_ccs_matrix * matrix);

void taucs_dtl(ccs_split) (taucs_ccs_matrix * A, taucs_ccs_matrix ** L, taucs_ccs_matrix ** R, int p);
void taucs_ccs_split(taucs_ccs_matrix * A, taucs_ccs_matrix ** L, taucs_ccs_matrix ** R, int p);

taucs_ccs_matrix *taucs_dtl(ccs_permute_symmetrically) (taucs_ccs_matrix * A, int *perm, int *invperm);
taucs_ccs_matrix *taucs_ccs_permute_symmetrically(taucs_ccs_matrix * A, int *perm, int *invperm);

void taucs_dtl(ccs_times_vec) (taucs_ccs_matrix * m, double *X, double *B);
void taucs_ccs_times_vec(taucs_ccs_matrix * m, void *X, void *B);


taucs_ccs_matrix *taucs_dtl(ccs_augment_nonpositive_offdiagonals) (taucs_ccs_matrix * A);
taucs_ccs_matrix *taucs_ccs_augment_nonpositive_offdiagonals(taucs_ccs_matrix * A);

void taucs_ccs_order(taucs_ccs_matrix * matrix, int **perm, int **invperm, char *which);

taucs_ccs_matrix *taucs_dtl(ccs_factor_llt) (taucs_ccs_matrix * A, double droptol, int modified);
taucs_ccs_matrix *taucs_ccs_factor_llt(taucs_ccs_matrix * A, double droptol, int modified);
taucs_ccs_matrix *taucs_ccs_factor_llt_partial(taucs_ccs_matrix * A, int p);
taucs_ccs_matrix *taucs_dtl(ccs_factor_llt_partial) (taucs_ccs_matrix * A, int p);
taucs_ccs_matrix *taucs_dtl(ccs_factor_ldlt) (taucs_ccs_matrix * A);
taucs_ccs_matrix *taucs_ccs_factor_ldlt(taucs_ccs_matrix * A);
taucs_ccs_matrix *taucs_dtl(ccs_factor_xxt) (taucs_ccs_matrix * A);

int taucs_ccs_solve_llt(void *L, void *x, void *b);
int taucs_dtl(ccs_solve_llt) (void *L, double *x, double *b);
int taucs_ccs_solve_ldlt(void *L, void *x, void *b);
int taucs_dtl(ccs_solve_ldlt) (void *L, double *x, double *b);

taucs_ccs_matrix *taucs_ccs_generate_mesh2d(int n, char *which);
taucs_ccs_matrix *taucs_ccs_generate_mesh2d_negative(int n);
taucs_ccs_matrix *taucs_ccs_generate_mesh3d(int X, int Y, int Z);
taucs_ccs_matrix *taucs_ccs_generate_dense(int m, int n, int flags);
taucs_ccs_matrix *taucs_ccs_generate_rrn(int X, int Y, int Z, double drop_probability, double rmin);
taucs_ccs_matrix *taucs_ccs_generate_discontinuous(int X, int Y, int Z, double jump);
double *taucs_vec_generate_continuous(int X, int Y, int Z, char *which);

int taucs_conjugate_gradients(taucs_ccs_matrix * A,
			      int (*precond_fn)(void *, void *x, void *b), void *precond_args, double *X, double *B, int itermax,
			      double convergetol);

int taucs_minres(taucs_ccs_matrix * A, int (*precond_fn)(void *, void *x, void *b), void *precond_args, double *X, double *B, int itermax,
		 double convergetol);

int taucs_sg_preconditioner_solve(void *P, double *z, double *r);

void *taucs_sg_preconditioner_create(taucs_ccs_matrix * A, int **perm, int **invperm, char *ordering, char *specification);
void taucs_sg_preconditioner_free(void *P);

taucs_ccs_matrix *taucs_amwb_preconditioner_create(taucs_ccs_matrix * symccs_mtxA, int rnd, double subgraphs);

void *taucs_recursive_amwb_preconditioner_create(taucs_ccs_matrix * A,
						 double c, double epsilon, int nsmall, int maxlevels, int innerits, double convratio, int **perm,
						 int **invperm);

int taucs_recursive_amwb_preconditioner_solve(void *P, void *Z, void *R);

int taucs_dtl(ccs_etree) (taucs_ccs_matrix * A, int *parent, int *l_colcount, int *l_rowcount, int *l_nnz);

int taucs_dtl(ccs_symbolic_elimination) (taucs_ccs_matrix * A, void *L, int do_order);

void *taucs_dtl(ccs_factor_llt_symbolic) (taucs_ccs_matrix * A);
int taucs_dtl(ccs_factor_llt_numeric) (taucs_ccs_matrix * A, void *L);

void *taucs_dtl(ccs_factor_llt_mf) (taucs_ccs_matrix * A);
void *taucs_dtl(ccs_factor_llt_ll) (taucs_ccs_matrix * A);
int taucs_dtl(supernodal_solve_llt) (void *vL, void *x, void *b);
void taucs_dtl(supernodal_factor_free) (void *L);
void taucs_dtl(supernodal_factor_free_numeric) (void *L);
taucs_ccs_matrix *taucs_dtl(supernodal_factor_to_ccs) (void *L);

int taucs_ccs_etree(taucs_ccs_matrix * A, int *parent, int *l_colcount, int *l_rowcount, int *l_nnz);

int taucs_ccs_symbolic_elimination(taucs_ccs_matrix * A, void *L, int do_order);

void *taucs_ccs_factor_llt_symbolic(taucs_ccs_matrix * A);
int taucs_ccs_factor_llt_numeric(taucs_ccs_matrix * A, void *L);

void *taucs_ccs_factor_llt_mf(taucs_ccs_matrix * A);
void *taucs_ccs_factor_llt_ll(taucs_ccs_matrix * A);
int taucs_supernodal_solve_llt(void *vL, void *x, void *b);
void taucs_supernodal_factor_free(void *L);
void taucs_supernodal_factor_free_numeric(void *L);
taucs_ccs_matrix *taucs_supernodal_factor_to_ccs(void *L);

void taucs_logfile(char *file_prefix);
int taucs_printf(char *fmt, ...);
double taucs_wtime();
double taucs_ctime();

__END_DECLS
#endif
