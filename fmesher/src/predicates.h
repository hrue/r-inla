#include <cstddef>

#ifndef _FMESH_PREDICATES_
#define _FMESH_PREDICATES_ 1

namespace fmesh {
namespace predicates {

/* #define SINGLE */
#ifdef SINGLE
typedef float REAL;
#else  /* not SINGLE */
typedef double REAL;
#endif /* not SINGLE */
typedef const REAL CREAL;

REAL orient2dfast(CREAL *pa, CREAL *pb, CREAL *pc);
REAL orient2d(CREAL *pa, CREAL *pb, CREAL *pc);

REAL orient3dfast(CREAL *pa, CREAL *pb, CREAL *pc, CREAL *pd);
REAL orient3d(CREAL *pa, CREAL *pb, CREAL *pc, CREAL *pd);

REAL incirclefast(CREAL *pa, CREAL *pb, CREAL *pc, CREAL *pd);
REAL incircle(CREAL *pa, CREAL *pb, CREAL *pc, CREAL *pd);

REAL inspherefast(CREAL *pa, CREAL *pb, CREAL *pc, CREAL *pd, CREAL *pe);
REAL insphere(CREAL *pa, CREAL *pb, CREAL *pc, CREAL *pd, CREAL *pe);

} // namespace predicates
} // namespace fmesh

#endif
