#ifndef __GMRFLib_GRAPH_EDIT_H__
#       define __GMRFLib_GRAPH_EDIT_H__

#       include <stdlib.h>
#       include <stddef.h>
#       include <math.h>

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

/*!
   \brief The various values for a tag
*/

/*
 */
    typedef struct {
	int n;
	int n_alloc;
	map_ii *Q;
} GMRFLib_ged_tp;

int GMRFLib_ged_add(GMRFLib_ged_tp * ged, int node, int nnode);
int GMRFLib_ged_append_graph(GMRFLib_ged_tp * ged, GMRFLib_graph_tp * graph);
int GMRFLib_ged_build(GMRFLib_graph_tp ** graph, GMRFLib_ged_tp * ged);
int GMRFLib_ged_free(GMRFLib_ged_tp * ged);
int GMRFLib_ged_init(GMRFLib_ged_tp ** ged, GMRFLib_graph_tp * graph);
int GMRFLib_ged_init2(GMRFLib_ged_tp ** ged, int n);
int GMRFLib_ged_insert_graph(GMRFLib_ged_tp * ged, GMRFLib_graph_tp * graph, int at_node);
int GMRFLib_ged_insert_graph2(GMRFLib_ged_tp * ged, GMRFLib_graph_tp * graph, int at_i_node, int at_j_node);

__END_DECLS
#endif
