#ifndef __GMRFLib_HA_IDX_INTERFACE_H__
#define __GMRFLib_HA_IDX_INTERFACE_H__

#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif

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

int ha_idx_q(void *ha, int key);
void *ha_idx_init(void);
void *ha_idx_init_hint(int);
void ha_idx_free(void *);
void ha_idx_set(void *ha, int key);
void ha_idx_sets(void *ha, int n, int *keys);
void ha_idx_stats(void *, int print, double *mb);

__END_DECLS
#endif
