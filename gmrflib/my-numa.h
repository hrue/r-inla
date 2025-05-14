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

__BEGIN_DECLS int GMRFLib_numa_have(void);
int GMRFLib_numa_get_node(void);
int GMRFLib_numa_node_of_ptr(void *ptr);
int GMRFLib_numa_nodes(void);
size_t GMRFLib_get_L3_cache(void);
size_t GMRFLib_numa_get_L3_cache(int nnode);
void *GMRFLib_numa_alloc_onnode(size_t size, int node);
void GMRFLib_numa_free(void *start, size_t size);
void GMRFLib_numa_get(int *cpu, int *numa);
void GMRFLib_numa_init(void);
void GMRFLib_numa_set_ctl(int enable);

__END_DECLS
#endif
