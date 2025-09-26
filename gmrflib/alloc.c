#include <stdint.h>
#include <assert.h>
#include <stddef.h>
#include <string.h>
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"

void * malloc_intern(size_t size) 
{
	void *p = aligned_alloc(GMRFLib_MEM_ALIGN, size); assert(p);
	return p;
}

void * calloc_intern(size_t nmemb, size_t size) 
{
	void *p = aligned_alloc(GMRFLib_MEM_ALIGN, nmemb * size); assert(p);
	memset(p, 0, nmemb * size);
	return p;
}

size_t GMRFLib_align_len(size_t n, size_t size)
{
	// return 'N >= n' so that the endpoint is aligned at GMRFLib_MEM_ALIGN bytes boundary
	int mm = (size_t) GMRFLib_MEM_ALIGN / size;
	div_t d = div(n, mm);
	return n + (d.rem == 0 ? 0 : mm - d.rem);
}

forceinline int GMRFLib_is_aligned(void *ptr)
{
	return ((uintptr_t)ptr % GMRFLib_MEM_ALIGN == 0);
}

forceinline int GMRFLib_is_aligned2(void *ptr, void *pptr)
{
	return ((uintptr_t)ptr % GMRFLib_MEM_ALIGN == 0 && (uintptr_t)pptr % GMRFLib_MEM_ALIGN == 0);
}

forceinline int GMRFLib_is_aligned3(void *ptr, void *pptr, void *ppptr)
{
	return ((uintptr_t)ptr % GMRFLib_MEM_ALIGN == 0 && (uintptr_t)pptr % GMRFLib_MEM_ALIGN == 0 && (uintptr_t) ppptr % GMRFLib_MEM_ALIGN == 0);
}
