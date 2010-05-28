#include <cstddef>
#include <cstring>
#include <set>
#include <map>
#include <sstream>
#include <cmath>

#include "vector.h"
#include "ioutils.h"

#define WHEREAMI __FILE__ << "(" << __LINE__ << ")\t"

#ifdef DEBUG
#define VECTOR_LOG(msg) std::cout << WHEREAMI << msg;
#else
#define VECTOR_LOG(msg)
#endif


using std::cout;
using std::endl;

namespace fmesh {

  struct SparseDatablock {
    int r;
    int c;
    double value;
  };


  std::ostream& IOHeaderOStream(std::ostream& output,
				const IOHeader& IO_header,
				bool ascii)
  {
    if (ascii) {
      output << IO_header.elems << " "
	     << IO_header.rows << " "
	     << IO_header.cols << " "
	     << IO_header.datatype << " "
	     << IO_header.valuetype << " "
	     << IO_header.matrixtype;
    } else {
      int sz = sizeof(IOHeader);
      output.write((char*)&sz, sizeof(sz));
      output.write((char*)&IO_header, sizeof(IOHeader));
    }
    return output;
  }

  std::istream& IOHeaderIStream(std::istream& input,
				IOHeader& IO_header,
				bool ascii)
  {
    if (ascii) {
    } else {
      int sz;
      input.read((char*)&sz, sizeof(sz));
      if (sz < sizeof(IOHeader)) {
	IO_header.elems = 0;
	IO_header.rows = 0;
	IO_header.cols = 0;
	IO_header.datatype = IOHeader::Datatype_dense;
	IO_header.valuetype = IOHeader::Valuetype_int;
	IO_header.matrixtype = IOHeader::Matrixtype_general;
	input.read((char*)&IO_header, sz);
      } else {
	input.read((char*)&IO_header, sizeof(IOHeader));
	if (sz > sizeof(IOHeader)) {
	  char* buf = new char[sizeof(IOHeader)-sz];
	  input.read(buf, sizeof(IOHeader)-sz);
	  delete[] buf;
	}
      }
    }
    return input;
  }

  std::istream& MatrixIntIStream(std::istream& input,
				 const Matrix<int>& M,
				 bool ascii);
  std::ostream& MatrixIntOStream(std::ostream& output,
				 const Matrix<int>& M,
				 bool ascii);
  std::ostream& MatrixDoubleOStream(std::ostream& output,
				    const Matrix<double>& M,
				    bool ascii);
  std::istream& MatrixDoubleIStream(std::istream& input,
				    const Matrix<double>& M,
				    bool ascii);
  std::ostream& SparseMatrixDoubleOStream(std::ostream& output,
					  const SparseMatrix<double>& M,
					  bool ascii);
  std::istream& SparseMatrixDoubleIStream(std::istream& input,
					  const SparseMatrix<double>& M,
					  bool ascii);


} /* namespace fmesh */
