#ifndef _FMESH_IOUTILS_T_
#define _FMESH_IOUTILS_T_ 1

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


  template <class T>
  IOHeader& IOHeader::DefaultDense(const Matrix<T>& M,
				   Matrixtype matrixt)
  {
    version = IOHEADER_VERSION;
    switch (matrixt) {
    case IOHeader::Matrixtype_general:
      elems = M.rows()*M.cols();
      rows = M.rows();
      cols = M.cols();
      break;
    case IOHeader::Matrixtype_symmetric:
      if (M.rows() <= M.cols()) {
	elems = (M.rows()*(M.rows()+1))/2;
	rows = M.rows();
	cols = M.rows();
      } else {
	elems = (M.cols()*(M.cols()+1))/2;
	rows = M.cols();
	cols = M.cols();
      }
      break;
    case IOHeader::Matrixtype_diagonal:
      if (M.rows() <= M.cols()) {
	elems = M.rows();
	rows = M.rows();
	cols = M.rows();
      } else {
	elems = M.cols();
	rows = M.cols();
	cols = M.cols();
      }
      break;
    }
    datatype = Datatype_dense;
    valuetype = Valuetype_double; /* Must be overridden! */
    matrixtype = matrixt;
    storagetype = Storagetype_rowmajor;
    return *this;
  }
  template <class T>
  IOHeader& IOHeader::DefaultSparse(const SparseMatrix<T>& M,
				    Matrixtype matrixt)
  {
    version = IOHEADER_VERSION;
    switch (matrixt) {
    case IOHeader::Matrixtype_general:
      elems = M.nnz();
      rows = M.rows();
      cols = M.cols();
      break;
    case IOHeader::Matrixtype_symmetric:
      if (M.rows() <= M.cols()) {
	elems = 0;
	rows = M.rows();
	cols = M.rows();
      } else {
	elems = 0;
	rows = M.cols();
	cols = M.cols();
      }
      for (int i=0; i<rows; i++)
	for (int j=i; j<cols; j++)
	  if (M.non_zero(i,j)) elems++;
      break;
    case IOHeader::Matrixtype_diagonal:
      if (M.rows() <= M.cols()) {
	elems = 0;
	rows = M.rows();
	cols = M.rows();
      } else {
	elems = 0;
	rows = M.cols();
	cols = M.cols();
      }
      for (int i=0; i<rows; i++)
	if (M.non_zero(i,i)) elems++;
      break;
    }
    datatype = Datatype_sparse;
    valuetype = Valuetype_double; /* Must be overridden! */
    matrixtype = matrixt;
    storagetype = Storagetype_colmajor;
    return *this;
  }
  template <class T>
  IOHeader& IOHeader::DefaultMap(const Matrix1int& m,
				 const Matrix<T>& M)
  {
    version = IOHEADER_VERSION;
    elems = M.rows();
    rows = M.rows();
    cols = M.cols();
    datatype = Datatype_map;
    valuetype = Valuetype_double; /* Must be overridden! */
    matrixtype = Matrixtype_general;
    storagetype = Storagetype_colmajor;
    return *this;
  }




  struct SparseDatablock {
    int r;
    int c;
    double value;
  };

  template <class T>
  IOHelper& IOHelper::O(std::ostream& output,
			    const Matrix<T>& M)
  {
    if (ascii_) {
      output << std::setprecision(15) << std::scientific;
      switch (h_.matrixtype) {
      case IOHeader::Matrixtype_general:
	if (h_.storagetype == IOHeader::Storagetype_rowmajor) {
	  for (int i=0; i<M.rows(); i++) {
	    const T* Mrow = M[i];
	    for (int j=0; j+1<M.cols(); j++) {
	      output << Mrow[j] << " ";
	    }
	    output << Mrow[M.cols()-1] << std::endl;
	  }
	} else {
	  for (int j=0; j<M.cols(); j++) {
	    for (int i=0; i+1<M.rows(); i++) {
	      output << M[i][j] << " ";
	    }
	    output << M[M.rows()-1][j] << std::endl;
	  }
	}
	break;
      case IOHeader::Matrixtype_symmetric:
	if (h_.storagetype == IOHeader::Storagetype_rowmajor) {
	  for (int i=0; i<h_.rows; i++) {
	    const T* Mrow = M[i];
	    for (int j=i; j+1<h_.cols; j++) {
	      output << Mrow[j] << " ";
	    }
	    output << Mrow[h_.cols-1] << std::endl;
	  }
	} else {
	  for (int j=0; j<h_.cols; j++) {
	    for (int i=0; i<j; i++) {
	      output << M[i][j] << " ";
	    }
	    output << M[j][j] << std::endl;
	  }
	}
	break;
      case IOHeader::Matrixtype_diagonal:
	if (h_.storagetype == IOHeader::Storagetype_rowmajor) {
	  for (int i=0; i<h_.rows; i++) {
	    output << M[i][i] << std::endl;
	  }
	} else {
	  for (int i=0; i+1<h_.rows; i++) {
	    output << M[i][i] << " ";
	  }
	  output << M[h_.rows-1][h_.rows-1] << std::endl;
	}
	break;
      }
    } else {
      switch (h_.matrixtype) {
      case IOHeader::Matrixtype_general:
	if (h_.storagetype == IOHeader::Storagetype_rowmajor) {
	  output.write((char*)M.raw(), sizeof(T)*M.rows()*M.cols());
	} else {
	  for (int j=0; j<M.cols(); j++) {
	    for (int i=0; i<M.rows(); i++) {
	      output.write((char*)&M[i][j], sizeof(T));
	    }
	  }
	}
	break;
      case IOHeader::Matrixtype_symmetric:
	if (h_.storagetype == IOHeader::Storagetype_rowmajor) {
	  for (int i=0; i<h_.rows; i++) {
	    const T* Mrow = M[i];
	    output.write((char*)&Mrow[i], sizeof(T)*(h_.cols-i));
	  }
	} else {
	  for (int j=0; j<h_.cols; j++) {
	    for (int i=0; i<j; i++) {
	      output.write((char*)&M[i][j], sizeof(T));
	    }
	  }
	}
	break;
      case IOHeader::Matrixtype_diagonal:
	for (int i=0; i<h_.rows; i++) {
	  output.write((char*)&M[i][i], sizeof(T));
	}
	break;
      }
    }
    return *this;
  }

  template <class T>
  IOHelper& IOHelper::I(std::istream& input,
			Matrix<T>& M)
  {
    NOT_IMPLEMENTED;
    return *this;
  }
  
  template <class T>
  IOHelper& IOHelper::O(std::ostream& output,
			const SparseMatrix<T>& M)
  {
    NOT_IMPLEMENTED;
    return *this;
  }
  
  template <class T>
  IOHelper& IOHelper::I(std::istream& input,
			SparseMatrix<T>& M)
  {
    NOT_IMPLEMENTED;
    return *this;
  }
  
  template <class T>
  IOHelper& IOHelper::O(std::ostream& output,
			const Matrix1int& m,
			const Matrix<T>& M)
  {
    NOT_IMPLEMENTED;
    return *this;
  }
  
  template <class T>
  IOHelper& IOHelper::I(std::istream& input,
			Matrix1int& m,
			Matrix<T>& M)
  {
    NOT_IMPLEMENTED;
    return *this;
  }


  
} /* namespace fmesh */

#endif
