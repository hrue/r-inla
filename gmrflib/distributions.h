#ifndef __GMRFLib_DISTRIBUTIONS_H__
#       define __GMRFLib_DISTRIBUTIONS_H__

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
**
 */
double GMRFLib_stdnormal(void);
double GMRFLib_scale_proposal(double F);
double GMRFLib_Wishart_logdens(gsl_matrix * Q, double r, gsl_matrix * R);

__END_DECLS
#endif
