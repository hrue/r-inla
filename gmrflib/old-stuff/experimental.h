#ifndef __GMRFLib_EXPERIMENTAL_H__
#define __GMRFLib_EXPERIMENTAL_H__

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

/* ... */
int GMRFLib_graph_read_binary_EXPERIMENTAL(GMRFLib_graph_tp ** graph, const char *filename);
int GMRFLib_graph_write_b_EXPERIMENTAL(const char *filename, GMRFLib_graph_tp * graph);

__END_DECLS
#endif
