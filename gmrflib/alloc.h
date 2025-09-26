#ifndef __GMRFLib_ALLOC_H__
#define __GMRFLib_ALLOC_H__

#include <stddef.h>
#include <stdlib.h>

#if defined(__linux) && defined(__AVX2__)
#define GMRFLib_MEM_ALIGN 32u
#else
#define GMRFLib_MEM_ALIGN 16u
#endif

void * malloc_intern(size_t size);
void * calloc_intern(size_t nmemb, size_t size);

size_t GMRFLib_align_len(size_t n, size_t size);
int GMRFLib_is_aligned(void *ptr);
int GMRFLib_is_aligned2(void *ptr, void *pptr);
int GMRFLib_is_aligned3(void *ptr, void *pptr, void *ppptr);

#endif
