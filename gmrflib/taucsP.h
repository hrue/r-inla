#define TAUCS_CORE_DOUBLE				       /* we use the double version only */

#ifndef __GMRFLib_TAUCSP_H__
#define __GMRFLib_TAUCSP_H__

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

#define TAUCS_DOUBLE    2048
#define TAUCS_LOWER      1
#define TAUCS_UPPER      2
#define TAUCS_TRIANGULAR 4
#define TAUCS_SYMMETRIC  8
#define TAUCS_HERMITIAN  16
#define TAUCS_PATTERN    32

#define taucs_add(x,y) ((x)+(y))
#define taucs_sub(x,y) ((x)-(y))
#define taucs_mul(x,y) ((x)*(y))
#define taucs_div(x,y) ((x)/(y))
#define taucs_neg(x)   (-(x))
#define taucs_conj(x)  (x)
#define taucs_abs(x)   (fabs(x))
#define taucs_sqrt(x)  (sqrt(x))

#define taucs_im(x)   0.0
#define taucs_re(x)   (x)
#define taucs_minusone -1.0
#define taucs_one     1.0
#define taucs_zero    0.0

extern double taucs_dzero_const;
extern double taucs_done_const;
extern double taucs_dminusone_const;

#define taucs_zero_const     taucs_dzero_const
#define taucs_one_const      taucs_done_const
#define taucs_minusone_const taucs_dminusone_const

#define taucs_zero_real_const     taucs_dzero_const
#define taucs_one_real_const      taucs_done_const
#define taucs_minusone_real_const taucs_dminusone_const

#define taucs_gemm  dgemm_
#define taucs_potrf dpotrf_
#define taucs_herk  dsyrk_
#define taucs_trsm  dtrsm_

typedef struct {
	int n;
	int m;
	int flags;
	int *colptr;
	int *rowind;
	double *values;

} taucs_ccs_matrix;

typedef struct {
	int n;	
	int m;	
	int flags;
	int *rowptr;
	int *colind;
	double *values;
} taucs_crs_matrix;

#define taucs_dtl(X) taucs_d##X
#include "taucs_privateP.h"

#define TAUCS_CORE
#define TAUCS_CORE_REAL
#define TAUCS_CORE_DATATYPE TAUCS_DOUBLE

#define taucs_iszero(x) (((__typeof (x)) (x)) == 0)

#include <stdlib.h>

typedef struct {
	int flags;
	char uplo;					       /* 'u' for upper, 'l' for lower, ' ' don't know; prefer lower. */
	int n;						       /* size of matrix */
	int n_sn;					       /* number of supernodes */
	int *parent;					       /* supernodal elimination tree */
	int *first_child;
	int *next_child;
	int *sn_size;					       /* size of supernodes (diagonal block) */
	int *sn_up_size;				       /* size of subdiagonal update blocks */
	int **sn_struct;				       /* row structure of supernodes */
	int *sn_blocks_ld;				       /* lda of supernode blocks */
	double **sn_blocks;			       /* supernode blocks */
	int *up_blocks_ld;				       /* lda of update blocks */
	double **up_blocks;			       /* update blocks */
} supernodal_factor_matrix;

__END_DECLS
#endif
