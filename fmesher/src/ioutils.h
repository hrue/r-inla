#ifndef _FMESH_IOUTILS_
#define _FMESH_IOUTILS_ 1

#include <cstddef>
#include <iostream>
#include <iomanip>
#include <vector>
#include <set>
#include <map>
#include <list>
#include <string>

#ifndef NOT_IMPLEMENTED
#define NOT_IMPLEMENTED (std::cout					\
			 << __FILE__ << "(" << __LINE__ << ")\t"	\
			 << "NOT IMPLEMENTED: "				\
			 << __PRETTY_FUNCTION__ << std::endl);
#endif

namespace fmesh {

  struct IOHeader {
    enum Datatype {Datatype_dense=0,
		   Datatype_sparse=1,
		   Datatype_map=2};
    enum Valuetype {Valuetype_int=0,
		    Valuetype_double=1};
    enum Matrixtype {Matrixtype_general=0,
		     Matrixtype_symmetric=1,
		     Matrixtype_diagonal=2};
    int elems;
    int rows;
    int cols;
    int datatype;
    int valuetype;
    int matrixtype;
  };

  std::ostream& IOHeaderOStream(std::ostream& output,
				const IOHeader& IO_header,
				bool ascii);
  std::istream& IOHeaderIStream(std::istream& input,
				IOHeader& IO_header,
				bool ascii);
  std::ostream& MatrixIntOStream(std::ostream& output,
				 const Matrix<int>& M,
				 bool ascii);
  std::istream& MatrixIntIStream(std::istream& input,
				 Matrix<int>& M,
				 bool ascii);
  std::ostream& MatrixDoubleOStream(std::ostream& output,
				    const Matrix<double>& M,
				    bool ascii);
  std::istream& MatrixDoubleIStream(std::istream& input,
				    Matrix<double>& M,
				    bool ascii);
  std::ostream& SparseMatrixDoubleOStream(std::ostream& output,
					  const SparseMatrix<double>& M,
					  bool ascii);
  std::istream& SparseMatrixDoubleIStream(std::istream& input,
					  SparseMatrix<double>& M,
					  bool ascii);


} /* namespace fmesh */

#endif
