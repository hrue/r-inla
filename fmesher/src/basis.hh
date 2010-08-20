#ifndef _FMESH_BASIS_
#define _FMESH_BASIS_ 1

#include <cstddef>
#include <iostream>
#include <iomanip>
#include <vector>
#include <set>
#include <map>
#include <list>
#include <string>

#include "vector.hh"

#ifndef NOT_IMPLEMENTED
#define NOT_IMPLEMENTED (std::cout					\
			 << __FILE__ << "(" << __LINE__ << ")\t"	\
			 << "NOT IMPLEMENTED: "				\
			 << __PRETTY_FUNCTION__ << std::endl);
#endif

namespace fmesh {

  Matrix<double> spherical_harmonics(Matrix3<double>& S,
				     int max_order,
				     bool rotationally_symmetric);

} /* namespace fmesh */

#endif
