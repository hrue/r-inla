#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "GMRFLib/GMRFLib.h"

// set to 1 to enable and 0 to disable
#define USE_ALIGNMENT 1

void *malloc_intern(size_t size)
{
	void *p = NULL;

#if USE_ALIGNMENT && GMRFLib_MEM_ALIGN == 16u
	// quick alternative as most allocators align by default at 16
	p = malloc(size);
	assert(p);
	if (GMRFLib_is_aligned(p)) {
		return p;
	} else {
		free(p);
	}
#endif

#if USE_ALIGNMENT
	size_t newsize = sMAX(size, GMRFLib_MEM_ALIGN);
	int rem = (newsize % GMRFLib_MEM_ALIGN);
	newsize += (rem > 0 ? GMRFLib_MEM_ALIGN - rem : 0);
	p = aligned_alloc(GMRFLib_MEM_ALIGN, newsize);
	assert(p);
#else
	p = malloc(size);
	assert(p);
#endif
	return p;
}

void *calloc_intern(size_t nmemb, size_t size)
{
	void *p = NULL;

#if USE_ALIGNMENT && GMRFLib_MEM_ALIGN == 16u
	p = calloc(nmemb, size);
	assert(p);
	if (GMRFLib_is_aligned(p)) {
		return p;
	} else {
		free(p);
	}
#endif

#if USE_ALIGNMENT
	size_t newsize = sMAX(nmemb * size, GMRFLib_MEM_ALIGN);
	int rem = (newsize % GMRFLib_MEM_ALIGN);
	newsize += (rem > 0 ? GMRFLib_MEM_ALIGN - rem : 0);
	p = aligned_alloc(GMRFLib_MEM_ALIGN, newsize);
	assert(p);
	memset(p, 0, newsize);
#else
	p = calloc(nmemb, size);
	assert(p);
#endif
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
