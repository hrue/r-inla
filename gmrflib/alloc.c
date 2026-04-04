#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "GMRFLib/GMRFLib.h"

void *malloc_intern(size_t size)
{
	void *p = NULL;
	if (GMRFLib_memory_alignment) {
		p = aligned_alloc(GMRFLib_MEM_ALIGN, size);
	} else {
		p = malloc(size);
	}
	assert(p);
	return p;
}

void *calloc_intern(size_t nmemb, size_t size)
{
	void *p = NULL;
	if (GMRFLib_memory_alignment) {
		size_t n = nmemb * size;
		p = aligned_alloc(GMRFLib_MEM_ALIGN, n);
		assert(p);
		Memset(p, 0, n);
	} else {
		p = calloc(nmemb, size);
	}
	assert(p);
	return p;
}

void *realloc_intern(void *ptr, size_t size)
{
	void *p = NULL;
	if (!ptr) {
		p = malloc_intern(size);
	} else {
		if (GMRFLib_memory_alignment) {
#       if defined(INLA_WITH_MIMALLOC)
			void *mi_realloc_aligned(void *, size_t, size_t);
			p = mi_realloc_aligned(ptr, size, GMRFLib_MEM_ALIGN);
#       elif defined(INLA_WITH_JEMALLOC)
			// void *rallocx(void *, size_t, int); 
			p = rallocx(ptr, size, MALLOCX_ALIGN(GMRFLib_MEM_ALIGN));
#       else
			p = realloc(ptr, size);
#       endif
		} else {
			p = realloc(ptr, size);
		}
	}

	assert(p);
	return p;
}

size_t GMRFLib_align_len(size_t n, size_t size)
{
	// return 'N >= n' so that the endpoint is aligned at
	// GMRFLib_MEM_ALIGN bytes boundary

	int mm = (size_t) GMRFLib_MEM_ALIGN / size;
	div_t d = div(n, mm);

	return n + (d.rem == 0 ? 0 : mm - d.rem);
}
