#include <stddef.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#if !defined(__cplusplus)
#pragma GCC diagnostic ignored "-Wimplicit-function-declaration"
#endif

#include "GMRFLib/fsort.h"
#include "GMRFLib/fsort/fluxsort.h"

void quadfluxsort(void *array, size_t nmemb, size_t size, int (*cmp)(const void *, const void *))
{
	if(nmemb <= 64L) {
		quadsort(array, nmemb, size, cmp);
	} else {
		fluxsort(array, nmemb, size, cmp);
	}
}

#pragma GCC diagnostic pop
