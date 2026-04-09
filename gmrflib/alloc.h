#ifndef __GMRFLib_ALLOC_H__
#       define __GMRFLib_ALLOC_H__
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
#       include <stddef.h>
#       include <string.h>
#       include <stdlib.h>
#       include <stdint.h>
extern unsigned int GMRFLib_memory_alignment;

static inline int GMRFLib_mem_align_test(void *p)
{
	return (((uintptr_t) (p)) % GMRFLib_memory_alignment == 0 ? 1 : 0);
}

static inline int GMRFLib_is_aligned(void *p)
{
	return GMRFLib_mem_align_test(p);
}

static inline int GMRFLib_is_aligned2(void *p, void *pp)
{
	return (GMRFLib_is_aligned(p) && GMRFLib_is_aligned(pp));
}

static inline int GMRFLib_is_aligned3(void *p, void *pp, void *ppp)
{
	return (GMRFLib_is_aligned2(p, pp) && GMRFLib_is_aligned(ppp));
}

size_t GMRFLib_align_len(size_t n, size_t size);
void *calloc_intern(size_t nmemb, size_t size);
void *malloc_intern(size_t size);
void *realloc_intern(void *ptr, size_t size);

int GMRFLib_mem_align_test(void *p);
int GMRFLib_is_aligned(void *p);
int GMRFLib_is_aligned2(void *p, void *pp);
int GMRFLib_is_aligned3(void *p, void *pp, void *ppp);

#       define GMRFLib_ALLOC_SAFE_SIZE(n_, type_) ((size_t)(n_)*sizeof(type_) <= PTRDIFF_MAX ? (size_t)(n_) : (size_t)1)

#       if 0
#              define Calloc(n, type)         (type *)GMRFLib_calloc(GMRFLib_ALLOC_SAFE_SIZE(n, type), sizeof(type), __FILE__, __GMRFLib_FuncName, __LINE__)
#              define Malloc(n, type)         (type *)GMRFLib_malloc(GMRFLib_ALLOC_SAFE_SIZE((n) * sizeof(type), char), __FILE__, __GMRFLib_FuncName, __LINE__)
#              define Realloc(ptr, n, type)   (type *)GMRFLib_realloc((void *)ptr, GMRFLib_ALLOC_SAFE_SIZE((n)*sizeof(type), char), __FILE__, __GMRFLib_FuncName, __LINE__)
#              define Free(ptr)               if (ptr) {GMRFLib_free((void *)(ptr), __FILE__, __GMRFLib_FuncName, __LINE__); ptr=NULL;}
#       else
#              undef  GMRFLib_TRACE_MEMORY
#              define Calloc(n, type)         (type *)calloc_intern(GMRFLib_ALLOC_SAFE_SIZE(n, type), sizeof(type))
#              define Malloc(n, type)         (type *)malloc_intern(GMRFLib_ALLOC_SAFE_SIZE((n) * sizeof(type), type))
#              define Realloc(ptr, n, type)   (type *)realloc_intern((void *)ptr, GMRFLib_ALLOC_SAFE_SIZE((n) * sizeof(type), char))
#              define Free(ptr)               if (ptr) {free((void *)(ptr)); ptr=NULL;}
#       endif

#       define Memcpy(dest, src, n)    memcpy((void *) (dest), (void *) (src), GMRFLib_ALLOC_SAFE_SIZE(n, char))
#       define Memset(dest, value, n)  memset((void *) (dest), (int) (value), (size_t) (n))

__END_DECLS
#endif
