#include <cstddef>

#ifndef _FMESH_PREDICATES_
#define _FMESH_PREDICATES_ 1

namespace fmesh {
namespace predicates {

/* #define SINGLE */
#ifdef SINGLE
  typedef float REAL;
#else /* not SINGLE */
  typedef double REAL;
#endif /* not SINGLE */

  REAL orient2dfast(REAL *pa, REAL *pb, REAL *pc);
  REAL orient2d(REAL *pa, REAL *pb, REAL *pc);

  REAL orient3dfast(REAL *pa, REAL *pb, REAL *pc, REAL *pd);
  REAL orient3d(REAL *pa, REAL *pb, REAL *pc, REAL *pd);
  
  REAL incirclefast(REAL *pa, REAL *pb, REAL *pc, REAL *pd);
  REAL incircle(REAL *pa, REAL *pb, REAL *pc, REAL *pd);

  REAL inspherefast(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL *pe);
  REAL insphere(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL *pe);

} /* namespace fmesh */
} /* namespace predicates */

#endif
