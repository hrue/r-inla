#if !defined(TAUCS_TAUCS_H)
#       define TAUCS_TAUCS_H

#       include <assert.h>
#       include <stddef.h>
#       include <time.h>
#       include <omp.h>
#       if __has_include(<malloc.h>)
#              include <malloc.h>
#       endif

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
#       ifdef __GNUC__
#              define UNUSED(x) x __attribute__((__unused__))
#       else
#              define UNUSED(x) x
#       endif
#       define TAUCS_CONFIG_DREAL
#       define TAUCS_CONFIG_BASE
#       define TAUCS_CONFIG_METIS
#       define TAUCS_CONFIG_AMD
#       define TAUCS_CONFIG_COLAMD
#       define TAUCS_CONFIG_GENMMD
#       define TAUCS_CONFIG_ORDERING
#       define TAUCS_CONFIG_LLT
#       if __GNUC__ > 7
typedef size_t fortran_charlen_t;
#       else
typedef int fortran_charlen_t;
#       endif
#       if !defined(F_ONE)
#              define F_ONE (fortran_charlen_t)1
#       endif

#       ifdef __GNUC__
#              define POSSIBLY_UNUSED_FUNCTION(x) __attribute__((__unused__)) x
#       else
#              define POSSIBLY_UNUSED_FUNCTION(x) x
#       endif

#       pragma omp declare simd
static int POSSIBLY_UNUSED_FUNCTION(IMAX) (int a, int b) {
	return ((a) > (b) ? (a) : (b));
}

#       pragma omp declare simd
static int POSSIBLY_UNUSED_FUNCTION(IMIN) (int a, int b) {
	return ((a) < (b) ? (a) : (b));
}

#       define DMAX(a_, b_) fmax(a_, b_)
#       define DMIN(a_, b_) fmin(a_, b_)

static double UNUSED(taucs_zero_real_const) = 0.0;
static double UNUSED(taucs_one_real_const) = 1.0;
static double UNUSED(taucs_minusone_real_const) = -1.0;

#       define TAUCS_DOUBLE_IN_BUILD
#       define TAUCS_BLAS_UNDERSCORE 1
#       if defined(TAUCS_BLAS_UNDERSCORE)
#              define taucs_blas_name(x) (x##_)
#       elif defined(TAUCS_BLAS_NOUNDERSCORE)
#              define taucs_blas_name(x) (x)
#       else
#              error "taucs_blas_[no]underscore_test: linking with the BLAS failed both attempts"
#       endif

#       ifdef _WIN32
#       endif

#       define TAUCS_SUCCESS                       0
#       define TAUCS_ERROR                        -1
#       define TAUCS_ERROR_NOMEM                  -2
#       define TAUCS_ERROR_BADARGS                -3
#       define TAUCS_ERROR_INDEFINITE             -4
#       define TAUCS_ERROR_MAXDEPTH               -5

#       define TAUCS_DOUBLE    2048

#       define TAUCS_LOWER      1
#       define TAUCS_UPPER      2
#       define TAUCS_TRIANGULAR 4
#       define TAUCS_SYMMETRIC  8
#       define TAUCS_HERMITIAN  16
#       define TAUCS_PATTERN    32

#       define TAUCS_METHOD_LLT  1
#       define TAUCS_METHOD_LDLT 2
#       define TAUCS_METHOD_PLU  3

#       define TAUCS_VARIANT_SNMF 1
#       define TAUCS_VARIANT_SNLL 2

#       define taucs_dtl(X) taucs_d##X

#       if defined(TAUCS_CORE_DOUBLE) || defined(TAUCS_CORE_GENERAL)

#              define taucs_conj(x)  (x)
#              define taucs_im(x)    0.0
#              define taucs_re(x)    (x)
#              define taucs_minusone -1.0
#              define taucs_one      1.0
#              define taucs_zero     0.0
#              define taucs_abs(x)   (fabs(x))
#              define taucs_sqrt(x)  (sqrt(x))

#              define taucs_zero_const     taucs_zero_real_const
#              define taucs_one_const      taucs_one_real_const
#              define taucs_minusone_const taucs_minusone_real_const

#              define taucs_gemm  taucs_blas_name(dgemm)
#              define taucs_potrf taucs_blas_name(dpotrf)
#              define taucs_herk  taucs_blas_name(dsyrk)
#              define taucs_trsm  taucs_blas_name(dtrsm)
#       endif

typedef struct {
	int n;						       /* columns */
	int m;						       /* rows; don't use if symmetric */
	int flags;
	int *colptr;					       /* pointers to where columns begin in rowind and values. */
	int *rowind;					       /* row indices */
	double *values;
} taucs_ccs_matrix;

#       include "taucs_private.h"

#       define taucs_isnan(x) (isnan(x) != 0)
#       define taucs_isinf(x) (isinf(x))

#       ifdef TAUCS_CORE_DOUBLE
#              define TAUCS_CORE
#              define taucs_dtl(X) taucs_d##X
#              define taucs_iszero(x) ((x) == 0.0)
#       endif

#       ifdef TAUCS_CORE_GENERAL
#              define TAUCS_CORE
#       endif

#       define taucs_add(x,y) ((x)+(y))
#       define taucs_sub(x,y) ((x)-(y))
#       define taucs_mul(x,y) ((x)*(y))
#       define taucs_div(x,y) ((x)/(y))
#       define taucs_neg(x)   (-(x))

#       define taucs_dadd(x,y) ((x)+(y))
#       define taucs_dsub(x,y) ((x)-(y))
#       define taucs_dmul(x,y) ((x)*(y))
#       define taucs_ddiv(x,y) ((x)/(y))
#       define taucs_dneg(x)   (-(x))
#       define taucs_dconj(x)  (x)
#       define taucs_dimag(x)    0.0
#       define taucs_dreal(x)    (x)
#       define taucs_dminusone -1.0
#       define taucs_done      1.0
#       define taucs_dzero     0.0
#       define taucs_dabs(x)   (fabs(x))
#       define taucs_dsqrt(x)  (sqrt(x))

double taucs_get_nan(void);

double taucs_allocation_amount(void);
int taucs_allocation_count(void);
int taucs_allocation_attempts(void);
void taucs_allocation_assert_clean(void);
void taucs_allocation_mark_clean(void);
void taucs_allocation_induce_failure(int i);

#       include <stdlib.h>

void *taucs_malloc(size_t size);
void *taucs_calloc(size_t nmemb, size_t size);
void *taucs_realloc(void *ptr, size_t size);
void taucs_free(void *ptr);

#       if defined(TAUCS_CORE)

void *taucs_calloc_stub(size_t nmemb, size_t size);
void *taucs_malloc_stub(size_t size);
void *taucs_realloc_stub(void *ptr, size_t size);
void taucs_free_stub(void *ptr);

#              define realloc(x,y) taucs_must_not_call_realloc_directly(x,y)
#              define malloc(x)    taucs_must_not_call_malloc_directly(x)
#              define calloc(x,y)  taucs_must_not_call_calloc_directly(x,y)
#              define free(x)      taucs_must_not_call_free_directly(x)

#              define taucs_realloc(x,y) taucs_realloc_stub(x,y)
#              define taucs_malloc(x)    taucs_malloc_stub(x)
#              define taucs_calloc(x,y)  taucs_calloc_stub(x,y)
#              define taucs_free(x)      taucs_free_stub(x)

#       endif

#       ifndef max
#              define max(x,y) ( ((x) > (y)) ? (x) : (y) )
#       endif

#       ifndef min
#              define min(x,y) ( ((x) < (y)) ? (x) : (y) )
#       endif

extern int dreadhb_(char *, int *, int *, int *, int *, int *, double *);
extern int amdexa_(int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *);
extern int amdtru_(int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *);
extern int amdbar_(int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *);
extern int amdbarnew_(int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *);
extern int amdnew_(int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *);
extern int genmmd_(int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *);

#       if (defined(OSTYPE_irix) || defined(OSTYPE_solaris))

#              include <math.h>
#              include <ieeefp.h>
#              define isinf(x) (!finite((x)) && !isnan((x)))

#       elif defined(_WIN32)

#              include <float.h>
//#define isnan(x)  (_isnan(x))
//#define isinf(x)  (!(_finite(x)) && !(_isnan(x)))
#              define finite(x) (_finite(x))

#       endif

#       include <math.h>
#       include <float.h>

extern int taucs_potrf(char *, int *, double *, int *, int *, fortran_charlen_t);
extern int taucs_trsm(char *, char *, char *, char *,
		      int *, int *, double *, double *, int *, double *, int *,
		      fortran_charlen_t, fortran_charlen_t, fortran_charlen_t, fortran_charlen_t);
extern int taucs_gemm(char *, char *, int *, int *, int *,
		      double *, double *, int *, double *, int *, double *, double *, int *, fortran_charlen_t, fortran_charlen_t);
extern int taucs_herk(char *, char *, int *, int *, double *, double *, int *, double *, double *, int *, fortran_charlen_t, fortran_charlen_t);

double taucs_blas_name(dnrm2) (int *, double *, int *);

__END_DECLS
#endif
