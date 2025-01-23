#ifndef __GMRFLib_DESIGN_H__
#define __GMRFLib_DESIGN_H__

#include <stdlib.h>
#include <stddef.h>
#include <math.h>

#include "GMRFLib/fmesher-io.h"

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
    typedef struct {
	double **experiment;
	double *int_weight;
	int nexperiments;
	int nfactors;
	int std_scale;					       /* if true, then the weights are on a standardized scale. this is the normal case */
} GMRFLib_design_tp;

int GMRFLib_design_ccd(GMRFLib_design_tp ** design, int nfactors);
int GMRFLib_design_eb(GMRFLib_design_tp ** design, int nhyper);
int GMRFLib_design_free(GMRFLib_design_tp * design);
int GMRFLib_design_grid(GMRFLib_design_tp ** design, int nhyper);
int GMRFLib_design_print(FILE * fp, GMRFLib_design_tp * design);
int GMRFLib_design_prune(GMRFLib_design_tp * design, double prob);
int GMRFLib_design_read(GMRFLib_design_tp ** design, GMRFLib_matrix_tp * D, int std_scale);


__END_DECLS
#endif
