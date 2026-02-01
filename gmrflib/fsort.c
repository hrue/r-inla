#include <stddef.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#if !defined(__cplusplus)
#       pragma GCC diagnostic ignored "-Wimplicit-function-declaration"
#endif

#include "GMRFLib/fsort.h"
#include "GMRFLib/fsort/fluxsort.h"

void quadfluxsort(void *array, size_t nmemb, size_t size, int (*cmp)(const void *, const void *))
{
#if 0
	void *a = (void *) malloc(nmemb * size);
	void *aa = (void *) malloc(nmemb * size);

	Memcpy(a, array, nmemb * size);
	Memcpy(aa, array, nmemb * size);

	static double win = 0.0;
	double tq = 0.0, tf = 0.0;

	tq = -GMRFLib_timer();
	quadsort(a, nmemb, size, cmp);
	tq += GMRFLib_timer();

	tf = -GMRFLib_timer();
	fluxsort(aa, nmemb, size, cmp);
	tf += GMRFLib_timer();

	double lwin = ABS(tf - tq);
	win += lwin;
	printf("SORT n %1d winner %s win %g total.win %g\n", (int) nmemb, (tf < tq ? "f" : "q"), lwin, win);

	Free(a);
	Free(aa);
#endif

	// quadsort(array, nmemb, size, cmp);
	fluxsort(array, nmemb, size, cmp);
}

#pragma GCC diagnostic pop
