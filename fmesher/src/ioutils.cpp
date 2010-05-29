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

   // enum Datatype {Datatype_dense=0,
   //                 Datatype_sparse=1,
   //                 Datatype_map=2};
   //  /*! int/double */
   //  enum Valuetype {Valuetype_int=0,
   //                  Valuetype_double=1};
   //  /*! general/symmentric/diagonal */
   //  enum Matrixtype {Matrixtype_general=0,
   //                   Matrixtype_symmetric=1,
   //                   Matrixtype_diagonal=2};
   //  /*! rowmajor/colmajor */
   //  enum Storagetype {Storagetype_rowmajor=0,
   //                    Storagetype_colmajor=1};

  IOHelper::IOHelper(const Matrix<int>& M,
		     IOHeader::Matrixtype matrixt) : ascii_(ASCII_DEFAULT)
  {
    h_.DefaultDense(M,matrixt);
    h_.valuetype = IOHeader::Valuetype_int;
  }
  IOHelper::IOHelper(const Matrix<double>& M,
		     IOHeader::Matrixtype matrixt) : ascii_(ASCII_DEFAULT)
  {
    h_.DefaultDense(M,matrixt);
    h_.valuetype = IOHeader::Valuetype_double;
  }

  IOHelper::IOHelper(const SparseMatrix<int>& M,
		     IOHeader::Matrixtype matrixt) : ascii_(ASCII_DEFAULT)
  {
    h_.DefaultSparse(M,matrixt);
    h_.valuetype = IOHeader::Valuetype_int;
  }
  IOHelper::IOHelper(const SparseMatrix<double>& M,
		     IOHeader::Matrixtype matrixt) : ascii_(ASCII_DEFAULT)
  {
    h_.DefaultSparse(M,matrixt);
    h_.valuetype = IOHeader::Valuetype_double;
  }

  IOHelper::IOHelper(const Matrix1int& m,
		     const Matrix<int>& M) : ascii_(ASCII_DEFAULT)
  {
    h_.DefaultMap(m,M);
    h_.valuetype = IOHeader::Valuetype_int;
  }
  IOHelper::IOHelper(const Matrix1int& m,
		     const Matrix<double>& M) : ascii_(ASCII_DEFAULT)
  {
    h_.DefaultMap(m,M);
    h_.valuetype = IOHeader::Valuetype_double;
  }


  IOHelper& IOHelper::O(std::ostream& output)
  {
    if (ascii_) {
      int header_length = sizeof(IOHeader)/sizeof(int);
      output << header_length << " ";
      const int* ioheader_p = (const int*)&h_;
      for (int i=0; i<header_length; i++) {
	output << ioheader_p[i];
	if (i+1<header_length)
	  output << " ";
      }
      output << endl;
    } else {
      int header_size = sizeof(IOHeader);
      output.write((char*)&header_size, sizeof(header_size));
      output.write((const char*)&h_, header_size);
    }
    return *this;
  }

  IOHelper& IOHelper::I(std::istream& input)
  {
    if (ascii_) {
      int header_length = sizeof(IOHeader)/sizeof(int);
      int file_header_length;
      input >> file_header_length;
      if (file_header_length < header_length) {
	h_.DefaultDense(Matrix<int>(0),IOHeader::Matrixtype_general);
	int* ioheader_p = (int*)&h_;
	for (int i=0; i<file_header_length; i++)
	  input >> ioheader_p[i];
      } else {
	int* ioheader_p = (int*)&h_;
	for (int i=0; i<header_length; i++)
	  input >> ioheader_p[i];
	if (file_header_length > header_length) {
	  int* buf = new int[header_length-file_header_length];
	  for (int i=0; i<header_length-file_header_length; i++)
	    input >> buf[i];
	  delete[] buf;
	}
      }
    } else {
      int header_size = sizeof(IOHeader);
      int file_header_size;
      input.read((char*)&file_header_size, sizeof(file_header_size));
      if (file_header_size < header_size) {
	h_.DefaultDense(Matrix<int>(0),IOHeader::Matrixtype_general);
	input.read((char*)&h_, file_header_size);
      } else {
	input.read((char*)&h_, header_size);
	if (file_header_size > header_size) {
	  char* buf = new char[header_size-file_header_size];
	  input.read(buf, header_size-file_header_size);
	  delete[] buf;
	}
      }
    }
    return *this;
  }


} /* namespace fmesh */
