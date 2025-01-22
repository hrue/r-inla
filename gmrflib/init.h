#ifndef __GMRFLib_INIT_H__
#define __GMRFLib_INIT_H__

#include <stdlib.h>

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
#if !defined(GMRFLib_FALSE)

/*!
  \brief Define GMRFLib_FALSE as (0)
*/
#define GMRFLib_FALSE (0)
#endif
#if !defined(GMRFLib_TRUE)

/*!
  \brief Define GMRFLib_TRUE as (1)
*/
#define GMRFLib_TRUE  (1)
#endif
    __END_DECLS
#endif
