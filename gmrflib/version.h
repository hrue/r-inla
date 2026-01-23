
/*!
  \file version.h
  \brief Typedefs for \ref version.c
*/

#ifndef __GMRFLib_VERSION_H__
#       define __GMRFLib_VERSION_H__

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

/*
 */
int GMRFLib_version(FILE * fp);

__END_DECLS
#endif
