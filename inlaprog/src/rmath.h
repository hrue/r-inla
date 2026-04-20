#ifndef __INLA_RMATH_H__
#       define __INLA_RMATH_H__

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
#       define MATHLIB_STANDALONE
#       define MATHLIB_FUN(_fun) _fun
#       if defined(ISNAN)
#              undef ISNAN
#       endif
#       include <Rmath.h>
#       if defined(ISNAN)
#              undef ISNAN
// same definition in GMRFLibP.h
#              define ISNAN(x) (isnan(x) != 0)
#       endif
    __END_DECLS
#endif
