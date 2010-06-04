#include <cstddef>
#include <cstring>
#include <set>
#include <map>
#include <sstream>
#include <cmath>

#include "vector.hh"
#include "ioutils.hh"

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
   //                 Datatype_sparse=1};
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

  IOHeader& IOHeader::def(const int& ref) {
    def();
    valuetype = IOValuetype_int;
    return *this;
  }

  IOHeader& IOHeader::def(const double& ref) {
    def();
    valuetype = IOValuetype_double;
    return *this;
  }


  IOHeader& IOHeader::def() {
    version = IOHEADER_VERSION;
    elems = 0;
    rows = 0;
    cols = 0;
    datatype = -1;
    valuetype = -1;
    matrixtype = -1;
    storagetype = IOStoragetype_rowmajor;
    return *this;
  }



  std::ostream& operator<<(std::ostream& output, const IOHeader& h)
  {
    int header_length = sizeof(IOHeader)/sizeof(int);
    output << header_length << " ";
    const int* ioheader_p = (const int*)&h;
    for (int i=0; i<header_length; i++) {
      output << ioheader_p[i];
      if (i+1<header_length)
	output << " ";
    }
    return output;
  }


  std::istream& operator>>(std::istream& input, IOHeader& h)
  {
    int header_length = sizeof(IOHeader)/sizeof(int);
    int file_header_length;
    input >> file_header_length;
    if (file_header_length < header_length) {
      h.dense(Matrix<int>(0),IOMatrixtype_general);
      int* ioheader_p = (int*)&h;
      for (int i=0; i<file_header_length; i++)
	input >> ioheader_p[i];
    } else {
      int* ioheader_p = (int*)&h;
      for (int i=0; i<header_length; i++)
	input >> ioheader_p[i];
      if (file_header_length > header_length) {
	int* buf = new int[header_length-file_header_length];
	for (int i=0; i<header_length-file_header_length; i++)
	  input >> buf[i];
	delete[] buf;
      }
    }
    return input;
  }


} /* namespace fmesh */
