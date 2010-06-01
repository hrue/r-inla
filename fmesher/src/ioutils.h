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

    /* Sets defaults, and the valuetype matching T: */
    template <class T>
    IOHeader& def(const T& ref);
    IOHeader& def();
    /* Default values: */
    template <class T>
    IOHeader& dense(const Matrix<T>& M,
		    IOMatrixtype matrixt = IOMatrixtype_general);
    template <class T>
    IOHeader& sparse(const SparseMatrix<T>& M,
		     IOMatrixtype matrixt = IOMatrixtype_general);
    
    /* Constructor, that sets the valuetype matching T: */
    template <class T>
    IOHeader(const T& ref) { def(ref); };
  };

  /*
  template <>
  IOHeader::def(const int& ref);
  template <>
  IOHeader::def(const double& ref);
  */

  std::ostream& operator<<(std::ostream& output, const IOHeader& h);
  std::istream& operator>>(std::istream& output, IOHeader& h);

  /*! Base helper for input and output. */
  template <class T>
  class IOHelper {
  public:
    IOHeader h_;
    bool binary_;
  public:
    /* Constructors: */
    IOHelper() : h_(T()), binary_(BINARY_DEFAULT) {};
    IOHelper(const IOHeader& h)
      : h_(h), binary_(BINARY_DEFAULT) {};

    bool binaryformat() const {
      return binary_;
    };
    IOMatrixtype matrixtype() const {
      return (IOMatrixtype)h_.matrixtype;
    };
    IOStoragetype storagetype() const {
      return (IOStoragetype)h_.storagetype;
    };

    IOHelper<T>& ascii(bool set_ascii = true) {
      return binary(!set_ascii);
    };
    IOHelper<T>& binary(bool set_binary = true) {
      binary_ = set_binary;
      return *this;
    };
    IOHelper<T>& storagetype(IOStoragetype set_storage) {
      h_.storagetype = set_storage;
      return *this;
    };
    IOHelper<T>& rowmajor(bool set_rowmajor = true) {
      h_.storagetype = (set_rowmajor
			? IOStoragetype_rowmajor
			: IOStoragetype_rowmajor);
      return *this;
    };
    IOHelper<T>& colmajor(bool set_colmajor = true) {
      h_.storagetype = (set_colmajor
			? IOStoragetype_colmajor
			: IOStoragetype_colmajor);
      return *this;
    };

    /* Output/Input: */
    IOHelper<T>& OH(std::ostream& output);
    IOHelper<T>& IH(std::istream& input);
    IOHelper<T>& IH(const IOHeader& h);

  };

  /*! Helper for Matrix input and output. */
  template <class T>
  class IOHelperM : public IOHelper<T> {
  public:
    const Matrix<T> *cM_;
    Matrix<T> *M_;
  public:
    /* Constructors: */
    IOHelperM() : IOHelper<T>() {};
    IOHelperM(const IOHeader& h) : IOHelper<T>(h) {};
    IOHelperM<T>& cD(const Matrix<T>* M) {
      cM_ = M;
      M_ = NULL;
      IOHelper<T>::h_.dense(*M);
      return *this;
    };
    IOHelperM<T>& D(Matrix<T>* M)
    {
      cM_ = M;
      M_ = M;
      IOHelper<T>::h_.dense(*M);
      return *this;
    };

    IOHelperM<T>& matrixtype(IOMatrixtype matrixt) {
      IOHelper<T>::h_.dense(*cM_,matrixt);
      return *this;
    };

    /* Output/Input: */
    IOHelperM<T>& OD(std::ostream& output);
    IOHelperM<T>& ID(std::istream& input);

    /* Overloaded from IOHelper: */
    IOHelperM<T>& ascii(bool set_ascii = true) {
      IOHelper<T>::ascii(set_ascii); return *this; };
    IOHelperM<T>& binary(bool set_binary = true) {
      IOHelper<T>::binary(set_binary); return *this; };
    IOHelperM<T>& general() { return matrixtype(IOMatrixtype_general); };
    IOHelperM<T>& symmetric() { return matrixtype(IOMatrixtype_symmetric); };
    IOHelperM<T>& diagonal() { return matrixtype(IOMatrixtype_diagonal); };
    IOHelperM<T>& storagetype(IOStoragetype set_storage) {
      IOHelper<T>::storagetype(set_storage); return *this; };
    IOHelperM<T>& rowmajor(bool set_rowmajor = true) {
      IOHelper<T>::rowmajor(set_rowmajor); return *this; };
    IOHelperM<T>& colmajor(bool set_colmajor = true) {
      IOHelper<T>::colmajor(set_colmajor); return *this; };
    IOHelperM<T>& OH(std::ostream& output) {
      IOHelper<T>::OH(output); return *this; };
    IOHelperM<T>& IH(std::istream& input) {
      IOHelper<T>::IH(input); return *this; };
    IOHelperM<T>& IH(const IOHeader& h) {
      IOHelper<T>::IH(h); return *this; };
  };

  /*! Helper for SparseMatrix input and output. */
  template <class T>
  class IOHelperSM : public IOHelper<T> {
  public:
    const SparseMatrix<T> *cM_;
    SparseMatrix<T> *M_;
  public:
    /* Constructors: */
    IOHelperSM() : IOHelper<T>() {};
    IOHelperSM(const IOHeader& h) : IOHelper<T>(h) {};
    IOHelperSM<T>& cD(const SparseMatrix<T>* M) {
      cM_ = M;
      M_ = NULL;
      IOHelper<T>::h_.sparse(*M);
      IOHelper<T>::colmajor();
      return *this;
    };
    IOHelperSM<T>& D(SparseMatrix<T>* M)
    {
      cM_ = M;
      M_ = M;
      IOHelper<T>::h_.sparse(*M);
      IOHelper<T>::colmajor();
      return *this;
    };

    IOHelperSM<T>& matrixtype(IOMatrixtype matrixt) {
      IOHelper<T>::h_.sparse(*cM_,matrixt);
      return *this;
    };

    /* Output/Input: */
    IOHelperSM<T>& OD(std::ostream& output);
    IOHelperSM<T>& ID(std::istream& input);

    /* Overloaded from IOHelper: */
    IOHelperSM<T>& ascii(bool set_ascii = true) {
      IOHelper<T>::ascii(set_ascii); return *this; };
    IOHelperSM<T>& binary(bool set_binary = true) {
      IOHelper<T>::binary(set_binary); return *this; };
    IOHelperSM<T>& general() { return matrixtype(IOMatrixtype_general); };
    IOHelperSM<T>& symmetric() { return matrixtype(IOMatrixtype_symmetric); };
    IOHelperSM<T>& diagonal() { return matrixtype(IOMatrixtype_diagonal); };
    IOHelperSM<T>& storagetype(IOStoragetype set_storage) {
      IOHelper<T>::storagetype(set_storage); return *this; };
    IOHelperSM<T>& rowmajor(bool set_rowmajor = true) {
      IOHelper<T>::rowmajor(set_rowmajor); return *this; };
    IOHelperSM<T>& colmajor(bool set_colmajor = true) {
      IOHelper<T>::colmajor(set_colmajor); return *this; };
    IOHelperSM<T>& OH(std::ostream& output) {
      IOHelper<T>::OH(output); return *this; };
    IOHelperSM<T>& IH(std::istream& input) {
      IOHelper<T>::IH(input); return *this; };
    IOHelperSM<T>& IH(const IOHeader& h) {
      IOHelper<T>::IH(h); return *this; };
  };




} /* namespace fmesh */

#include "ioutils.tpp"

#endif
