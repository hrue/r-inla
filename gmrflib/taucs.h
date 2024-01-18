
/*!
  \file taucs.h
  \brief Typedefs and defines for the TAUCS library version 2.2

  This file is included to make GMRFLib compile without refering to the TAUCS library,  
*/

#define TAUCS_CORE_DOUBLE				       /* we use the double version only */

#ifndef __GMRFLib_TAUCS_H__
#define __GMRFLib_TAUCS_H__

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

/*********************************************************/

/* TAUCS                                                 */

/* Author: Sivan Toledo                                  */

/*********************************************************/
#define TAUCS_INT       1024
#define TAUCS_DOUBLE    2048
#define TAUCS_SINGLE    4096
#define TAUCS_DCOMPLEX  8192
#define TAUCS_SCOMPLEX 16384
#define TAUCS_LOWER      1
#define TAUCS_UPPER      2
#define TAUCS_TRIANGULAR 4
#define TAUCS_SYMMETRIC  8
#define TAUCS_HERMITIAN  16
#define TAUCS_PATTERN    32
typedef double taucs_double;
typedef float taucs_single;

/*
  hrue added the __STRICT__ANSI__
*/
#if defined(__GNUC__) && !defined(TAUCS_GENERIC_COMPLEX) && !defined(__STRICT_ANSI__)

typedef __complex__ double taucs_dcomplex;
typedef __complex__ float taucs_scomplex;

#define taucs_add(x,y) ((x)+(y))
#define taucs_sub(x,y) ((x)-(y))
#define taucs_mul(x,y) ((x)*(y))
#define taucs_div(x,y) ((x)/(y))
#define taucs_neg(x)   (-(x))

#if defined(TAUCS_CORE_DOUBLE)

#define taucs_conj(x)  (x)
#define taucs_im(x)    0.0
#define taucs_re(x)    (x)
#define taucs_minusone -1.0
#define taucs_one      1.0
#define taucs_zero     0.0
#define taucs_abs(x)   (fabs(x))
#define taucs_sqrt(x)  (sqrt(x))

#elif defined(TAUCS_CORE_SINGLE)

#define taucs_conj(x)  (x)
#define taucs_im(x)    0.0f
#define taucs_re(x)    (x)
#define taucs_minusone -1.0f
#define taucs_one      1.0f
#define taucs_zero     0.0f
#define taucs_abs(x)   ((taucs_single) fabs(x))
#define taucs_sqrt(x)  ((taucs_single) sqrt(x))

#elif defined(TAUCS_CORE_DCOMPLEX)

#define taucs_conj(x)  (~(x))
#define taucs_im(x)    (__imag__ (x))
#define taucs_re(x)    (__real__ (x))
#define taucs_minusone -1.0+0.0i
#define taucs_one      1.0+0.0i
#define taucs_zero     0.0+0.0i
#define taucs_abs(x)   taucs_zabs_fn(x)
#define taucs_sqrt(x)  taucs_zsqrt_fn(x)

#elif defined(TAUCS_CORE_SCOMPLEX)

#define taucs_conj(x)  (~(x))
#define taucs_im(x)    (__imag__ (x))
#define taucs_re(x)    (__real__ (x))
#define taucs_minusone -1.0f+0.0fi
#define taucs_one      1.0f+0.0fi
#define taucs_zero     0.0f+0.0fi
#define taucs_abs(x)   taucs_cabs_fn(x)
#define taucs_sqrt(x)  taucs_csqrt_fn(x)

#endif

#else							       /* __GNUC__ */

typedef struct {
	double r, i;
} taucs_dcomplex;
typedef struct {
	float r, i;
} taucs_scomplex;

#if defined(TAUCS_CORE_DOUBLE)

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

#elif defined(TAUCS_CORE_SINGLE)

#define taucs_add(x,y) ((x)+(y))
#define taucs_sub(x,y) ((x)-(y))
#define taucs_mul(x,y) ((x)*(y))
#define taucs_div(x,y) ((x)/(y))
#define taucs_neg(x)   (-(x))
#define taucs_conj(x)  (x)
#define taucs_abs(x)   ((taucs_single) fabs(x))
#define taucs_sqrt(x)  ((taucs_single) sqrt(x))

#define taucs_im(x)   0.0f
#define taucs_re(x)   (x)
#define taucs_minusone -1.0f
#define taucs_one     1.0f
#define taucs_zero    0.0f

#elif defined(TAUCS_CORE_DCOMPLEX)

#define taucs_add(x,y) taucs_zadd_fn(x,y)
#define taucs_sub(x,y) taucs_zsub_fn(x,y)
#define taucs_mul(x,y) taucs_zmul_fn(x,y)
#define taucs_div(x,y) taucs_zdiv_fn(x,y)
#define taucs_neg(x)   taucs_zneg_fn(x)
#define taucs_conj(x)  taucs_zconj_fn(x)
#define taucs_abs(x)   taucs_zabs_fn(x)
#define taucs_sqrt(x)  taucs_zsqrt_fn(x)

#define taucs_im(x)    ((x).i)
#define taucs_re(x)    ((x).r)
#define taucs_minusone taucs_zminusone_const
#define taucs_one      taucs_zone_const
#define taucs_zero     taucs_zzero_const

#elif defined(TAUCS_CORE_SCOMPLEX)

#define taucs_add(x,y) taucs_cadd_fn(x,y)
#define taucs_sub(x,y) taucs_csub_fn(x,y)
#define taucs_mul(x,y) taucs_cmul_fn(x,y)
#define taucs_div(x,y) taucs_cdiv_fn(x,y)
#define taucs_neg(x)   taucs_cneg_fn(x)
#define taucs_conj(x)  taucs_cconj_fn(x)
#define taucs_abs(x)   taucs_cabs_fn(x)
#define taucs_sqrt(x)  taucs_csqrt_fn(x)

#define taucs_im(x)    ((x).i)
#define taucs_re(x)    ((x).r)
#define taucs_minusone taucs_cminusone_const
#define taucs_one      taucs_cone_const
#define taucs_zero     taucs_czero_const

#endif

#endif

extern taucs_double taucs_dzero_const;
extern taucs_double taucs_done_const;
extern taucs_double taucs_dminusone_const;

extern taucs_single taucs_szero_const;
extern taucs_single taucs_sone_const;
extern taucs_single taucs_sminusone_const;

extern taucs_dcomplex taucs_zzero_const;
extern taucs_dcomplex taucs_zone_const;
extern taucs_dcomplex taucs_zminusone_const;

extern taucs_scomplex taucs_czero_const;
extern taucs_scomplex taucs_cone_const;
extern taucs_scomplex taucs_cminusone_const;

#define taucs_isnan(x) (isnan((double)(taucs_re(x))) || isnan((double)(taucs_im(x))))
#define taucs_isinf(x) (isinf((double)(taucs_re(x))) || isinf((double)(taucs_im(x))))

#ifdef TAUCS_CORE_SINGLE
#define taucs_zero_const     taucs_szero_const
#define taucs_one_const      taucs_sone_const
#define taucs_minusone_const taucs_sminusone_const

#define taucs_zero_real_const     taucs_szero_const
#define taucs_one_real_const      taucs_sone_const
#define taucs_minusone_real_const taucs_sminusone_const

#define taucs_gemm  sgemm_
#define taucs_potrf spotrf_
#define taucs_herk  ssyrk_
#define taucs_trsm  strsm_
#endif

#ifdef TAUCS_CORE_DOUBLE
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
#endif

#ifdef TAUCS_CORE_SCOMPLEX
#define taucs_zero_const     taucs_czero_const
#define taucs_one_const      taucs_cone_const
#define taucs_minusone_const taucs_cminusone_const

#define taucs_zero_real_const     taucs_szero_const
#define taucs_one_real_const      taucs_sone_const
#define taucs_minusone_real_const taucs_sminusone_const

#define taucs_gemm  cgemm_
#define taucs_potrf cpotrf_
#define taucs_herk  cherk_
#define taucs_trsm  ctrsm_
#endif

#ifdef TAUCS_CORE_DCOMPLEX
#define taucs_zero_const     taucs_zzero_const
#define taucs_one_const      taucs_zone_const
#define taucs_minusone_const taucs_zminusone_const

#define taucs_zero_real_const     taucs_dzero_const
#define taucs_one_real_const      taucs_done_const
#define taucs_minusone_real_const taucs_dminusone_const

#define taucs_gemm  zgemm_
#define taucs_potrf zpotrf_
#define taucs_herk  zherk_
#define taucs_trsm  ztrsm_
#endif

/*********************************************************/

/*                                                       */

/*********************************************************/

typedef struct {
	int n;						       /* columns */
	int m;						       /* rows; don't use if symmetric */
	int flags;
	int *colptr;					       /* pointers to where columns begin in rowind and values. */
	/*
	 * 0-based. Length is (n+1). 
	 */
	int *rowind;					       /* row indices */

	union {
		void *v;
		taucs_double *d;
		taucs_single *s;
		taucs_dcomplex *z;
		taucs_scomplex *c;
	} values;

} taucs_ccs_matrix;

typedef struct {
	int n;						       /* columns */
	int m;						       /* rows; don't use if symmetric */
	int flags;
	int *rowptr;					       /* pointers to where columns begin in rowind and values. */
	/*
	 * 0-based. Length is (n+1). 
	 */
	int *colind;					       /* row indices */

	union {
		void *v;
		taucs_double *d;
		taucs_single *s;
		taucs_dcomplex *z;
		taucs_scomplex *c;
	} values;

} taucs_crs_matrix;

typedef struct {
	int type;
	int nmatrices;
	void *type_specific;

	/*
	 * the following may change! do not rely on them. 
	 */
	double nreads, nwrites, bytes_read, bytes_written, read_time, write_time;
} taucs_io_handle;

/* generate all the prototypes */

#define taucs_datatype taucs_double
#define taucs_dtl(X) taucs_d##X
#include "taucs_private.h"
#undef taucs_datatype
#undef taucs_dtl

#define taucs_datatype taucs_single
#define taucs_dtl(X) taucs_s##X
#include "taucs_private.h"
#undef taucs_datatype
#undef taucs_dtl

#define taucs_datatype taucs_dcomplex
#define taucs_dtl(X) taucs_z##X
#include "taucs_private.h"
#undef taucs_datatype
#undef taucs_dtl

#define taucs_datatype taucs_scomplex
#define taucs_dtl(X) taucs_c##X
#include "taucs_private.h"
#undef taucs_datatype
#undef taucs_dtl

/*********************************************************/

/*                                                       */

/*********************************************************/

/* now define the data type for the file that we compile now */

#ifdef TAUCS_CORE_DOUBLE
#define TAUCS_CORE
#define TAUCS_CORE_REAL
#define TAUCS_CORE_DATATYPE TAUCS_DOUBLE
typedef taucs_double taucs_datatype;

#define taucs_dtl(X) taucs_d##X
#define taucs_values values.d
#define taucs_iszero(x) ((x) == 0.0)
#endif

#ifdef  TAUCS_CORE_SINGLE
#define TAUCS_CORE
#define TAUCS_CORE_REAL
#define TAUCS_CORE_DATATYPE TAUCS_SINGLE
typedef taucs_single taucs_datatype;

#define taucs_dtl(X) taucs_s##X
#define taucs_values values.s
#define taucs_iszero(x) ((x) == 0.0f)
#endif

#ifdef  TAUCS_CORE_DCOMPLEX
#define TAUCS_CORE
#define TAUCS_CORE_COMPLEX
#define TAUCS_CORE_DATATYPE TAUCS_DCOMPLEX
typedef taucs_dcomplex taucs_datatype;

#define taucs_dtl(X) taucs_z##X
#define taucs_values values.z
#define taucs_iszero(x) (taucs_re(x) == 0.0 && taucs_im(x) == 0.0)
#endif

#ifdef  TAUCS_CORE_SCOMPLEX
#define TAUCS_CORE
#define TAUCS_CORE_COMPLEX
#define TAUCS_CORE_DATATYPE TAUCS_SCOMPLEX
typedef taucs_scomplex taucs_datatype;

#define taucs_dtl(X) taucs_c##X
#define taucs_values values.c
#define taucs_iszero(x) (taucs_re(x) == 0.0f && taucs_im(x) == 0.0f)
#endif

#ifdef TAUCS_CORE_REAL
#endif

/*********************************************************/

/*                                                       */

/*********************************************************/

#ifdef OSTYPE_irix
#include <ieeefp.h>

#define isinf(x) (!finite((x)) && !isnan((x)))
#endif

/*********************************************************/

/*                                                       */

/*********************************************************/

#ifndef max
#define max(x,y) ( ((x) > (y)) ? (x) : (y) )
#endif

#ifndef min
#define min(x,y) ( ((x) < (y)) ? (x) : (y) )
#endif

/*********************************************************/

/*                                                       */

/*********************************************************/

/* 
   routines for testing memory allocation.
   Mostly useful for testing programs
   that hunt for memory leaks.
*/

double taucs_allocation_amount();
int taucs_allocation_count();
int taucs_allocation_attempts();
void taucs_allocation_assert_clean();
void taucs_allocation_mark_clean();
void taucs_allocation_induce_failure(int i);

/* 
   these are meant to allow allocation 
   and more importantly, deallocation,
   within the testing programs.
*/

#include <stdlib.h>

void *taucs_MALLOC(size_t size);
void *taucs_CALLOC(size_t nmemb, size_t size);
void *taucs_REALLOC(void *ptr, size_t size);
void taucs_FREE(void *ptr);

#if defined(TAUCS_CORE) && defined(TAUCS_MEMORY_TEST_yes)

#include <stdlib.h>

void *taucs_internal_calloc(size_t nmemb, size_t size, char *file, int line);
void *taucs_internal_malloc(size_t size, char *file, int line);
void *taucs_internal_realloc(void *ptr, size_t size, char *file, int line);
void taucs_internal_free(void *ptr, char *file, int line);

#define realloc(x,y) taucs_internal_realloc(x,y,__FILE__,__LINE__)
#define malloc(x)    taucs_internal_malloc(x,__FILE__,__LINE__)
#define calloc(x,y)  taucs_internal_calloc(x,y,__FILE__,__LINE__)
#define free(x)      taucs_internal_free(x,__FILE__,__LINE__)

#endif

/*********************************************************/

/*                                                       */

/*********************************************************/

/* 
   this is needed, as its only define in TAUCS c-source.
*/
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
	taucs_datatype **sn_blocks;			       /* supernode blocks */

	int *up_blocks_ld;				       /* lda of update blocks */
	taucs_datatype **up_blocks;			       /* update blocks */
} supernodal_factor_matrix;

__END_DECLS
#endif
