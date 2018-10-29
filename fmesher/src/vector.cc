#include <cstddef>
#include <cstring>
#include <set>
#include <map>
#include <sstream>
#include <cmath>

#include "vector.hh"

#define WHEREAMI __FILE__ << "(" << __LINE__ << ")\t"

#ifdef DEBUG
#define VECTOR_LOG(msg) std::cout << WHEREAMI << msg;
#else
#define VECTOR_LOG(msg)
#endif


using std::cout;
using std::endl;

namespace fmesh {


  double Vec::length(const Point& s0)
  {
    return s0.length();
  }

  /*!
    Calculate an arbitrary perpendicular vector.
    
    Michael M. Stark, Efficient Construction of Perpendicular
    Vectors without Branching, Journal of graphics, gpu, and game
    tools, Vol. 14, No. 1: 55-62, 2009
  */
#define ABS(X) std::fabs(X)
#define SIGNBIT(X) ((unsigned int)(std::signbit(X) != 0))
  void arbitrary_perpendicular(Vector3<double>& n,
			       const Vector3<double>& v)
  {
    double s_0 = v.s[0];
    double s_1 = v.s[1];
    double s_2 = v.s[2];
    const unsigned int uyx = SIGNBIT(ABS(s_0) - ABS(s_1));
    const unsigned int uzx = SIGNBIT(ABS(s_0) - ABS(s_2));
    const unsigned int uzy = SIGNBIT(ABS(s_1) - ABS(s_2));
    const unsigned int xm = uyx & uzx;
    const unsigned int ym = (1^xm) & uzy;
    const unsigned int zm = 1^(xm & ym);
    n.s[0] = zm*s_1 - ym*s_2;
    n.s[1] = xm*s_2 - zm*s_0;
    n.s[2] = ym*s_0 - xm*s_1;
  }



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
