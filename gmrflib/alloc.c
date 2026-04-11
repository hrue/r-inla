#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "GMRFLib/GMRFLib.h"

#if defined(_WIN32)
#       define aligned_alloc(a_, b_) malloc(b_)
int GMRFLib_memory_alignment_enabled = 0;
#else
int GMRFLib_memory_alignment_enabled = 1;
#endif

// will be over-rided by an adaptive default in inla-parse.c
unsigned int GMRFLib_memory_alignment = 16;

void *malloc_intern(size_t size)
{
	void *p = malloc(size);
	assert(p);
	return p;
}

void *amalloc_intern(size_t size)
{
	void *p = NULL;
	if (GMRFLib_memory_alignment_enabled) {
		size_t new_size = GMRFLib_align_len(size, 1);
		p = aligned_alloc(GMRFLib_memory_alignment, new_size);
	} else {
		p = malloc(size);
	}
	assert(p);
	return p;
}

void *calloc_intern(size_t nmemb, size_t size)
{
	void *p = calloc(nmemb, size);
	assert(p);
	return p;
}

void *acalloc_intern(size_t nmemb, size_t size)
{
	void *p = NULL;
	if (GMRFLib_memory_alignment_enabled) {
		size_t n = nmemb * size;
		size_t nn = GMRFLib_align_len(n, 1);
		p = aligned_alloc(GMRFLib_memory_alignment, nn);
		assert(p);
		Memset(p, 0, nn);
	} else {
		p = calloc(nmemb, size);
	}
	assert(p);
	return p;
}

void *realloc_intern(void *ptr, size_t size)
{
	void *p = NULL;
	if (ptr) {
		p = realloc(ptr, size);
	} else {
		p = malloc(size);
	}
	assert(p);
	return p;
}

void *arealloc_intern(void *ptr, size_t size)
{
	void *p = NULL;
	if (!ptr) {
		p = malloc_intern(size);
	} else {
		if (GMRFLib_memory_alignment_enabled) {
			p = realloc(ptr, size);
			// this is a workaround for not having aligned realloc.
			// good thing is that 'size' is valid so we can do Memcpy!
			if (!GMRFLib_is_aligned(p)) {
				size_t new_size = GMRFLib_align_len(size, 1);
				void *pp = malloc_intern(new_size);
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
	// GMRFLib_memory_alignment bytes boundary.

	// assume that _memory_alignment divisible by size (oops: no check)
	int mm = (size_t) GMRFLib_memory_alignment / size;
	div_t d = div(n, mm);

	return n + (d.rem == 0 ? 0 : mm - d.rem);
}
