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

#include "vector.h"

#ifndef NOT_IMPLEMENTED
#define NOT_IMPLEMENTED (std::cout					\
			 << __FILE__ << "(" << __LINE__ << ")\t"	\
			 << "NOT IMPLEMENTED: "				\
			 << __PRETTY_FUNCTION__ << std::endl);
#endif

#define IOHEADER_VERSION 0
#define BINARY_DEFAULT false

namespace fmesh {

  /*! dense/sparse/map */
  enum IODatatype {IODatatype_dense=0,
		   IODatatype_sparse=1};
  /*! int/double */
  enum IOValuetype {IOValuetype_int=0,
		    IOValuetype_double=1};
  /*! general/symmentric/diagonal */
  enum IOMatrixtype {IOMatrixtype_general=0,
		     IOMatrixtype_symmetric=1,
		     IOMatrixtype_diagonal=2};
  /*! rowmajor/colmajor */
  enum IOStoragetype {IOStoragetype_rowmajor=0,
		      IOStoragetype_colmajor=1};

  /*! Header for input and output file formats. */
  class IOHeader {
  public:
    int version; /*!< Format version */
    int elems; /*!< The number of data units
		 
		 For dense matrices, the total number of elements.
		 For sparse matrices, the number of elements contained in
		 the file.
		*/
    int rows; /*!< The number of data rows. */
    int cols; /*!< The number of data columns. */
    int datatype; /*!< The IODatatype. */
    int valuetype; /*!< The IOValuetype. */
    int matrixtype; /*!< The IOMatrixtype. */
    int storagetype; /*!< The IOStoragetype. */

    /* Default values: */
    template <class T>
    IOHeader& DefaultDense(const Matrix<T>& M,
			   IOMatrixtype matrixt = IOMatrixtype_general);
    template <class T>
    IOHeader& DefaultSparse(const SparseMatrix<T>& M,
			    IOMatrixtype matrixt = IOMatrixtype_general);
    
    /* Constructor, that sets the valuetype matching T: */
    template <class T>
    IOHeader(const T& ref);
  };

  template <>
  IOHeader::IOHeader(const int& ref);
  template <>
  IOHeader::IOHeader(const double& ref);

  std::ostream& operator<<(std::ostream& output, const IOHeader& h);
  std::istream& operator>>(std::istream& output, IOHeader& h);

  /*! Helper for input and output. */
  class IOHelper {
  public:
    IOHeader h_;
    bool binary_;
  public:
    /* Constructors: */
    IOHelper() : h_(0), binary_(BINARY_DEFAULT) {};
    IOHelper(const IOHeader& h) : h_(h), binary_(BINARY_DEFAULT) {};
    template <class T>
    IOHelper(const Matrix<T>& M,
	     IOMatrixtype matrixt = IOMatrixtype_general,
	     bool set_binary = BINARY_DEFAULT)
      : h_(T()), binary_(set_binary)
    { h_.DefaultDense(M,matrixt); };
    template <class T>
    IOHelper(const SparseMatrix<T>& M,
	     IOMatrixtype matrixt = IOMatrixtype_general,
	     bool set_binary = BINARY_DEFAULT)
      : h_(T()), binary_(set_binary)
    { h_.DefaultSparse(M,matrixt); };

    bool binaryformat() const {
      return binary_;
    };
    IOStoragetype storage() const {
      return (IOStoragetype)h_.storagetype;
    };
    IOHelper& ascii(bool set_ascii = true) {
      return binary(!set_ascii);
    };
    IOHelper& binary(bool set_binary = true) {
      binary_ = set_binary;
      return *this;
    };
    IOHelper& storage(IOStoragetype set_storage) {
      h_.storagetype = set_storage;
      return *this;
    };
    IOHelper& rowmajor(bool set_rowmajor = true) {
      h_.storagetype = (set_rowmajor
			? IOStoragetype_rowmajor
			: IOStoragetype_rowmajor);
      return *this;
    };
    IOHelper& colmajor(bool set_colmajor = true) {
      h_.storagetype = (set_colmajor
			? IOStoragetype_colmajor
			: IOStoragetype_colmajor);
      return *this;
    };

    /* Output/Input: */
    IOHelper& O(std::ostream& output);
    IOHelper& I(std::istream& input);
    IOHelper& H(const IOHeader& h);
    template <class T>
    IOHelper& O(std::ostream& output,
		const Matrix<T>& M);
    template <class T>
    IOHelper& I(std::istream& input,
		Matrix<T>& M);
    template <class T>
    IOHelper& O(std::ostream& output,
		const SparseMatrix<T>& M);
    template <class T>
    IOHelper& I(std::istream& input,
		SparseMatrix<T>& M);

  };



} /* namespace fmesh */

#include "ioutils.tpp"

#endif
