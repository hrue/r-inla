#ifndef __GMRFLib_BITMAP_H__
#       define __GMRFLib_BITMAP_H__

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

__BEGIN_DECLS int GMRFLib_bitmap_graph__intern(GMRFLib_graph_tp * graph, const char *filename, int *mapping);
int GMRFLib_bitmap_image(const char *filename, GMRFLib_uchar * image, int nx, int ny);
int GMRFLib_bitmap_graph(const char *filename_body, int *remap, GMRFLib_graph_tp * graph);
int GMRFLib_bitmap_problem(const char *filename_body, GMRFLib_problem_tp * problem);

__END_DECLS
#endif
