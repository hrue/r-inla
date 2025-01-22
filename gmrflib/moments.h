#ifndef __GMRFLib_MOMENTS_H__
#define __GMRFLib_MOMENTS_H__

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
// ...
#include "GMRFLib/GMRFLib.h"
double GMRFLib_noncentral_moment(int i, int j, double mi, double mj, double Sii, double Sij, double Sjj);
double GMRFLib_central_moment(int i, int j, double Sii, double Sij, double Sjj);

__END_DECLS
#endif
