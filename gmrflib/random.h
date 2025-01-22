#ifndef __GMRFLib_RANDOM_H__
#define __GMRFLib_RANDOM_H__

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
#include "GMRFLib/GMRFLib.h"
#include <stdlib.h>
#define GMRFLib_rng (GMRFLib_rng_ptr ? GMRFLib_rng_ptr :  GMRFLib_rng_init_default())
    gsl_rng * GMRFLib_rng_init_default(void);
double GMRFLib_rng_uniform(void);
int GMRFLib_rng_init(unsigned long int seed);
int GMRFLib_rng_set_default_seed(void);
int GMRFLib_rng_setstate(void *saved_state);
void *GMRFLib_rng_getstate(size_t *siz);

__END_DECLS
#endif
