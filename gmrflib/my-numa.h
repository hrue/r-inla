#ifndef __GMRFLib_MY_NUMA_H__
#define __GMRFLib_MY_NUMA_H__

#include <stdlib.h>
#include <stddef.h>
#include <math.h>

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
#include "GMRFLib/GMRFLibP.h"
void GMRFLib_numa_get(int *cpu, int *numa);
int GMRFLib_numa_nodes(void);
int GMRFLib_numa(void);

__END_DECLS
#endif
