#ifndef __INLA_TWEEDIE_H__
#       define __INLA_TWEEDIE_H__
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
 *
 */
#       include "GMRFLib/GMRFLib.h"
#       include "GMRFLib/GMRFLibP.h"
#       include "inla.h"
#       include "my.h"
// ***
void dtweedie_init_cache(void);
void dtweedie(int n, double y, double *mu, double phi, double p, double *ldens);
double ptweedie(double y, double mu, double phi, double p);

__END_DECLS
#endif
