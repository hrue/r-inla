#include <cstddef>
#include <cstring>
#include <set>
#include <map>
#include <sstream>
#include <cmath>

#include "locator.hh"

#define WHEREAMI __FILE__ << "(" << __LINE__ << ")\t"

#define LOG_(msg) std::cout << WHEREAMI << msg;
#ifdef DEBUG
#define LOG(msg) LOG_(msg)
#else
#define LOG(msg)
#endif


using std::cout;
using std::endl;

namespace fmesh {


  std::ostream& TriangleLocator::print(std::ostream& output)
  {
    return bbox_locator_.print(output);
  }
 
  std::ostream& operator<<(std::ostream& output, TriangleLocator& locator)
  {
    return locator.print(output);
  }


} /* namespace fmesh */
