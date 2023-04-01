#ifndef _FMESH_RCPP_FMESHER_
#define _FMESH_RCPP_FMESHER_ 1
#ifdef FMESHER_WITH_R

#include "RcppCommon.h"
#include "RcppEigen.h"

namespace fmesh {
template <class T> class Matrix;
template <class T> class Matrix1;
template <class T> class Matrix3;
template <class T> class SparseMatrix;
class MatrixC;
}

/* forward declarations */
namespace Rcpp {
/* support for wrap */

#define FM_DEFINE_WRAP(__thetype__)                     \
  template<> inline SEXP wrap(const fmesh::__thetype__& obj);

FM_DEFINE_WRAP(Matrix<double>);
FM_DEFINE_WRAP(Matrix<int>);
FM_DEFINE_WRAP(Matrix1<double>);
FM_DEFINE_WRAP(Matrix1<int>);
FM_DEFINE_WRAP(Matrix3<double>);
FM_DEFINE_WRAP(Matrix3<int>);
FM_DEFINE_WRAP(SparseMatrix<double>);
FM_DEFINE_WRAP(MatrixC);

// TODO:
/* support for as */
// template<typename T> class Exporter< fmesh::Matrix<T> >;

}

#endif
#endif
