#include <metis.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#if !defined(INLA_WITH_PARDISO)

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
int METIS51PARDISO_NodeND(int *nvtxs, int *xadj, int *adjncy, int *vwgt, int *options, int *perm, int *iperm);
void pardiso(void *a, int *b, int *c, int *d, int *e, int *f, double *g, int *h, int *i, int *j, int *k, int *l, int *m, double *n, double *o, int *p, double *q);
void pardiso_chkmatrix(int *a, int *s, double *d, int *f, int *g, int *h);
void pardiso_chkvec(int *a, int *s, double *d, int *f);
void pardiso_copy_symbolic_factor_single(void *a, void *b, int *c, int *d, double *e, double *f, int *g, int *h, int *i);
void pardiso_delete_symbolic_factor_single(void *a, int *b, int *c);
void pardiso_get_factor_csc(void **a, double *s, int *d, int *f, double *g, int *h, int *j, int *k, int l);
void pardiso_get_inverse_factor_csc(void **a, double *s, int *d, int *f, int *g, int h);
void pardiso_printstats(int *a, int *s, double *d, int *f, int *g, int *h, double *j, int *k);
void pardiso_residual(int *mtype, int *n, double *a, int *ia, int *ja, double *b, double *x, double *y, double *norm_b, double *norm_res);
void pardisoinit(void *a, int *b, int *c, int *d, double *e, int *f);
__END_DECLS


/* 
   this creates an empty version of the pardiso functions with an appropriate wrapper to the relevant original metis function. Since
   the (pardiso modified) metis library is bundled into the libpardiso.so, I have to do the same.
 */

// do not change: also GMRFLib/smtp-pardiso.c uses this code
#define NOLIB_ECODE (270465)

#define NO_PARDISO_LIB						\
	{							\
		fprintf(stderr, "\n\n\t*** PARDISO library is not available. Exit.\n\n");	\
		exit(1);						\
	}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"

void pardisoinit(void *a, int *b, int *c, int *d, double *e, int *f)
{
	*f = NOLIB_ECODE;
	return;
}

void pardiso(void *a, int *b, int *c, int *d, int *e, int *f, double *g,
	     int *h, int *i, int *j, int *k, int *l, int *m, double *n, double *o, int *p, double *q) NO_PARDISO_LIB;
void pardiso_chkmatrix(int *a, int *s, double *d, int *f, int *g, int *h) NO_PARDISO_LIB;
void pardiso_chkvec(int *a, int *s, double *d, int *f) NO_PARDISO_LIB;
void pardiso_printstats(int *a, int *s, double *d, int *f, int *g, int *h, double *j, int *k) NO_PARDISO_LIB;
void pardiso_get_factor_csc(void **a, double *s, int *d, int *f, double *g, int *h, int *j, int *k, int l) NO_PARDISO_LIB;
void pardiso_get_inverse_factor_csc(void **a, double *s, int *d, int *f, int *g, int h) NO_PARDISO_LIB;
void pardiso_residual(int *mtype, int *n, double *a, int *ia, int *ja, double *b, double *x, double *y, double *norm_b,
		      double *norm_res) NO_PARDISO_LIB;
void pardiso_copy_symbolic_factor_single(void *a, void *b, int *c, int *d, double *e, double *f, int *g, int *h, int *i) NO_PARDISO_LIB;
void pardiso_delete_symbolic_factor_single(void *a, int *b, int *c) NO_PARDISO_LIB;

int METIS51PARDISO_NodeND(int *nvtxs, int *xadj, int *adjncy, int *vwgt, int *options, int *perm, int *iperm)
{
	int METIS_NodeND(int *, int *, int *, int *, int *, int *, int *);
	return METIS_NodeND(nvtxs, xadj, adjncy, vwgt, options, perm, iperm);
}

#pragma GCC diagnostic pop
#endif							       /* if !defined(INLA_WITH_PARDISO) */
