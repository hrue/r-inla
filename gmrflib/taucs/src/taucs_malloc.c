#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#if __has_include(<malloc.h>)
#       include <malloc.h>
#endif
#include "taucs.h"

#undef malloc
#undef calloc
#undef realloc
#undef free

int taucs_printf(const char *UNUSED(fmt), ...)
{
	return 0;
}

void *taucs_malloc_stub(size_t size)
{
	if (size) {
		assert(size < PTRDIFF_MAX);
		void *p = malloc(size);
		if (!p) {
			fprintf(stderr, "\n\ntaucs_malloc_stub fail to malloc %zu bytes\n", size);
			abort();
		}
		return p;
	} else {
		return NULL;
	}
	// return (size ? calloc(1, size) : NULL);
}

void *taucs_calloc_stub(size_t nmemb, size_t size)
{
	if (nmemb) {
		assert(nmemb * size < PTRDIFF_MAX);
		void *p = calloc(nmemb, size);
		if (!p) {
			fprintf(stderr, "\n\ntaucs_calloc_stub fail to calloc %zu elements of size %zu bytes\n", nmemb, size);
			abort();
		}
		return p;
	} else {
		return NULL;
	}
	// return (nmemb ? calloc(nmemb, size) : NULL);
}

void *taucs_realloc_stub(void *ptr, size_t size)
{
	assert(size < PTRDIFF_MAX);
	void *p = realloc(ptr, size);
	if (!p) {
		fprintf(stderr, "\n\ntaucs_realloc_stub fail to realloc %zu bytes\n", size);
		abort();
	}
	return p;
	// return realloc(ptr, size);
}

void taucs_free_stub(void *ptr)
{
	if (ptr) {
		free(ptr);
		ptr = NULL;
	}
}
