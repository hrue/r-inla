#ifndef __GMRFLib_REMAP_CACHE_H__
#define __GMRFLib_REMAP_CACHE_H__

#include <math.h>
#include <strings.h>
#include <stdlib.h>
#include <stdio.h>

#include "GMRFLib/sha.h"

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS					       /* empty */
#define __END_DECLS					       /* empty */
#endif

__BEGIN_DECLS

//

typedef struct 
{
	int n;
	int nrhs;
	int *remap;					       /* length n * nrhs */
	unsigned char *sha;
}
	GMRFLib_remap_tp;

int *GMRFLib_remap_get(int *remap, int n, int nrhs);
int GMRFLib_remap_init_store(void);
unsigned char *GMRFLib_remap_sha(int *remap, int n, int nrhs);
void icopy(int *n, int *x, int *ix, int *y, int *iy);


__END_DECLS
#endif
