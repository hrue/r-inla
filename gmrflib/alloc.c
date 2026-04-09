#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "GMRFLib/GMRFLib.h"

#if defined(_WIN32)
#       define aligned_alloc(a_, b_) malloc(b_)
bool not_windows = 0;
#else
bool not_windows = 1;
#endif

unsigned int GMRFLib_memory_alignment = 64;

void *malloc_intern(size_t size)
{
	void *p = NULL;
	if (not_windows) {
		p = aligned_alloc(GMRFLib_memory_alignment, size);
	} else {
		p = malloc(size);
	}
	assert(p);
	return p;
}

void *calloc_intern(size_t nmemb, size_t size)
{
	void *p = NULL;
	if (not_windows) {
		size_t n = nmemb * size;
		p = aligned_alloc(GMRFLib_memory_alignment, n);
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
		if (not_windows) {
			p = realloc(ptr, size);
			// this is a workaround for not having aligned realloc.
			// good thing is that 'size' is valid so we can do Memcpy!
			if (!GMRFLib_is_aligned(p)) {
				void *pp = malloc_intern(size);
				assert(pp);
				Memcpy(pp, p, size);
				Free(p);
				p = pp;
			}
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
	// GMRFLib_memory_alignment bytes boundary

	int mm = (size_t) GMRFLib_memory_alignment / size;
	div_t d = div(n, mm);

	return n + (d.rem == 0 ? 0 : mm - d.rem);
}
