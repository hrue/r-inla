#ifndef __GMRFLib_FSORT_H__
#       define __GMRFLib_FSORT_H__

#       include <stddef.h>
#       include <stdlib.h>

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

/*
 */
void fluxsort(void *array, size_t nmemb, size_t size, int (*cmp)(const void *, const void *));
void quadsort(void *array, size_t nmemb, size_t size, int (*cmp)(const void *, const void *));
void quadfluxsort(void *array, size_t nmemb, size_t size, int (*cmp)(const void *, const void *));

//#define QSORT_FUN quadsort
//#define QSORT_FUN fluxsort
#       define QSORT_FUN quadfluxsort
//#define QSORT_FUN qsort

__END_DECLS
#endif
