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
#define ASCII_DEFAULT true

namespace fmesh {

  /*! Header for input and output file formats. */
  class IOHeader {
  public:
    /*! dense/sparse/map */
    enum Datatype {Datatype_dense=0,
		   Datatype_sparse=1,
		   Datatype_map=2};
    /*! int/double */
    enum Valuetype {Valuetype_int=0,
		    Valuetype_double=1};
    /*! general/symmentric/diagonal */
    enum Matrixtype {Matrixtype_general=0,
		     Matrixtype_symmetric=1,
		     Matrixtype_diagonal=2};
    /*! rowmajor/colmajor */
    enum Storagetype {Storagetype_rowmajor=0,
		      Storagetype_colmajor=1};
    int version; /*!< Format version */
    int elems; /*!< The number of data units
		 
		 For dense matrices, the total number of elements.
		 For sparse matrices, the number of elements contained in
		 the file.
		 For maps, the number of mapped vectors.
		*/
    int rows; /*!< The number of data rows. */
    int cols; /*!< The number of data columns. */
    int datatype; /*!< The Datatype. */
    int valuetype; /*!< The Valuetype. */
    int matrixtype; /*!< The Matrixtype. */
    int storagetype; /*!< The Storagetype. */

    /* Default values: */
    template <class T>
    IOHeader& DefaultDense(const Matrix<T>& M,
			   Matrixtype matrixt = Matrixtype_general);
    template <class T>
    IOHeader& DefaultSparse(const SparseMatrix<T>& M,
			    Matrixtype matrixt = Matrixtype_general);
    template <class T>
    IOHeader& DefaultMap(const Matrix1int& m,
			 const Matrix<T>& M);
  };

  /*! Helper for input and output. */
  class IOHelper {
    IOHeader h_;
    bool ascii_;
  public:
    /* Constructors: */
    IOHelper() : ascii_(ASCII_DEFAULT) {};
    IOHelper(const Matrix<int>& M,
	     IOHeader::Matrixtype matrixt = IOHeader::Matrixtype_general);
    IOHelper(const Matrix<double>& M,
	     IOHeader::Matrixtype matrixt = IOHeader::Matrixtype_general);
    IOHelper(const SparseMatrix<int>& M,
	     IOHeader::Matrixtype matrixt = IOHeader::Matrixtype_general);
    IOHelper(const SparseMatrix<double>& M,
	     IOHeader::Matrixtype matrixt = IOHeader::Matrixtype_general);
    IOHelper(const Matrix1int& m, const Matrix<int>& M);
    IOHelper(const Matrix1int& m, const Matrix<double>& M);

    IOHelper& ascii() {
      ascii_ = !ascii_;
      return *this;
    };
    IOHelper& storage() {
      if (h_.storagetype == IOHeader::Storagetype_rowmajor)
	h_.storagetype = IOHeader::Storagetype_colmajor;
      else
	h_.storagetype = IOHeader::Storagetype_rowmajor;
      return *this;
    };

    /* Output: */
    IOHelper& O(std::ostream& output);
    IOHelper& I(std::istream& input);
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
    template <class T>
    IOHelper& O(std::ostream& output,
		const Matrix1int& m,
		const Matrix<T>& M);
    template <class T>
    IOHelper& I(std::istream& input,
		Matrix1int& m,
		Matrix<T>& M);

  };



} /* namespace fmesh */

#include "ioutils.tpp"

#endif
