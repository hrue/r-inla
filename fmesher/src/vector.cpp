#include <cstddef>
#include <cstring>
#include <set>
#include <map>
#include <sstream>
#include <cmath>

#include "vector.h"

#define WHEREAMI __FILE__ << "(" << __LINE__ << ")\t"

#ifdef DEBUG
#define VECTOR_LOG(msg) std::cout << WHEREAMI << msg;
#else
#define VECTOR_LOG(msg)
#endif


using std::cout;
using std::endl;

namespace fmesh {

  template<>
  std::ostream& operator<<(std::ostream& output,
			   const SparseMatrixTriplet<int>& MT);
  template<>
  std::ostream& operator<<(std::ostream& output,
			   const SparseMatrixTriplet<double>& MT);

  template<>
  std::istream& operator>>(std::istream& input,
			   SparseMatrixTriplet<int>& MT);
  template<>
  std::istream& operator>>(std::istream& input,
			   SparseMatrixTriplet<double>& MT);

} /* namespace fmesh */
