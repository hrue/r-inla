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
  IOHeader::IOHeader(const T& ref) { valuetype = -1; };

  template <class T>
  IOHeader& IOHeader::DefaultDense(const Matrix<T>& M,
				   IOMatrixtype matrixt)
  {
    version = IOHEADER_VERSION;
    switch (matrixt) {
    case IOMatrixtype_general:
      elems = M.rows()*M.cols();
      rows = M.rows();
      cols = M.cols();
      break;
    case IOMatrixtype_symmetric:
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
    case IOMatrixtype_diagonal:
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
    datatype = IODatatype_dense;
    /* valuetype is set by the constructor */
    matrixtype = matrixt;
    storagetype = IOStoragetype_rowmajor;
    return *this;
  }
  template <class T>
  IOHeader& IOHeader::DefaultSparse(const SparseMatrix<T>& M,
				    IOMatrixtype matrixt)
  {
    version = IOHEADER_VERSION;
    switch (matrixt) {
    case IOMatrixtype_general:
      elems = M.nnz();
      rows = M.rows();
      cols = M.cols();
      break;
    case IOMatrixtype_symmetric:
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
    case IOMatrixtype_diagonal:
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
    datatype = IODatatype_sparse;
    /* valuetype is set by the constructor */
    matrixtype = matrixt;
    storagetype = IOStoragetype_colmajor;
    return *this;
  }



  template <class T>
  IOHelper& IOHelper::O(std::ostream& output,
			const Matrix<T>& M)
  {
    if (!((h_.rows>0) && (h_.cols>0))) {
      return *this;
    }
    if (binary_) {
      switch (h_.matrixtype) {
      case IOMatrixtype_general:
	if ((h_.storagetype == IOStoragetype_rowmajor) ||
	    (M.cols()==1)) {
	  output.write((char*)M.raw(), sizeof(T)*h_.rows*h_.cols);
	} else {
	  for (int j=0; j<h_.cols; j++) {
	    for (int i=0; i<h_.rows; i++) {
	      output.write((char*)&M[i][j], sizeof(T));
	    }
	  }
	}
	break;
      case IOMatrixtype_symmetric:
	if (h_.storagetype == IOStoragetype_rowmajor) {
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
      case IOMatrixtype_diagonal:
	for (int i=0; i<h_.rows; i++) {
	  output.write((char*)&M[i][i], sizeof(T));
	}
	break;
      }
    } else { /* Text format. */
      output << std::setprecision(15) << std::scientific;
      switch (h_.matrixtype) {
      case IOMatrixtype_general:
	if (h_.storagetype == IOStoragetype_rowmajor) {
	  for (int i=0; i<h_.rows; i++) {
	    const T* Mrow = M[i];
	    for (int j=0; j+1<h_.cols; j++) {
	      output << Mrow[j] << " ";
	    }
	    output << Mrow[h_.cols-1] << std::endl;
	  }
	} else {
	  for (int j=0; j<h_.cols; j++) {
	    for (int i=0; i+1<h_.rows; i++) {
	      output << M[i][j] << " ";
	    }
	    output << M[h_.rows-1][j] << std::endl;
	  }
	}
	break;
      case IOMatrixtype_symmetric:
	if (h_.storagetype == IOStoragetype_rowmajor) {
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
      case IOMatrixtype_diagonal:
	if (h_.storagetype == IOStoragetype_rowmajor) {
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
    }
    return *this;
  }
  
  template <class T>
  IOHelper& IOHelper::I(std::istream& input,
			Matrix<T>& M)
  {
    M.clear();
    M.cols(h_.cols);
    M.capacity(h_.rows);
    M(h_.rows-1,h_.cols-1,T()); /* Initialize last element. */
    if (binary_) {
      switch (h_.matrixtype) {
      case IOMatrixtype_general:
	if ((h_.storagetype == IOStoragetype_rowmajor) ||
	    (h_.cols==1)) {
	  input.read((char*)M.raw(), sizeof(T)*h_.rows*h_.cols);
	} else {
	  for (int j=0; j<h_.cols; j++) {
	    for (int i=0; i<h_.rows; i++) {
	      input.read((char*)&M(i)[j], sizeof(T));
	    }
	  }
	}
	break;
      case IOMatrixtype_symmetric:
	if (h_.storagetype == IOStoragetype_rowmajor) {
	  for (int i=0; i<h_.rows; i++) {
	    T* Mrow = M(i);
	    input.read((char*)&Mrow[i], sizeof(T)*(h_.cols-i));
	  }
	} else {
	  for (int j=0; j<h_.cols; j++) {
	    for (int i=0; i<j; i++) {
	      input.read((char*)&M(i)[j], sizeof(T));
	    }
	  }
	}
	break;
      case IOMatrixtype_diagonal:
	for (int i=0; i<h_.rows; i++) {
	  input.read((char*)&M(i)[i], sizeof(T));
	}
	break;
      }
    } else { /* Text format. */
      switch (h_.matrixtype) {
      case IOMatrixtype_general:
	if (h_.storagetype == IOStoragetype_rowmajor) {
	  for (int i=0; i<h_.rows; i++) {
	    T* Mrow = M(i);
	    for (int j=0; j<h_.cols; j++) {
	      input >> Mrow[j];
	    }
	  }
	} else {
	  for (int j=0; j<h_.cols; j++) {
	    for (int i=0; i<h_.rows; i++) {
	      input >> M(i)[j];
	    }
	  }
	}
	break;
      case IOMatrixtype_symmetric:
	if (h_.storagetype == IOStoragetype_rowmajor) {
	  for (int i=0; i<h_.rows; i++) {
	    T* Mrow = M(i);
	    for (int j=i; j<h_.cols; j++) {
	      input >> Mrow[j];
	    }
	  }
	} else {
	  for (int j=0; j<h_.cols; j++) {
	    for (int i=0; i<j+1; i++) {
	      input >> M(i)[j];
	    }
	  }
	}
	break;
      case IOMatrixtype_diagonal:
	if (h_.storagetype == IOStoragetype_rowmajor) {
	  for (int i=0; i<h_.rows; i++) {
	    input >> M(i)[i];
	  }
	} else {
	  for (int i=0; i<h_.rows; i++) {
	    input >> M(i)[i];
	  }
	}
	break;
      }
    }
    return *this;
  }
  
  template <class T>
  IOHelper& IOHelper::O(std::ostream& output,
			const SparseMatrix<T>& M)
  {
    if (h_.storagetype == IOStoragetype_rowmajor) {
      if (h_.matrixtype == IOMatrixtype_diagonal) {
	Matrix1< SparseMatrixDuplet<T> > MT;
	M.tolist(MT);
	IOHelper(MT).binary(binary_).rowmajor().O(output,MT);
      } else {
	Matrix1< SparseMatrixTriplet<T> > MT;
	M.tolist(MT,(h_.matrixtype == IOMatrixtype_symmetric));
	IOHelper(MT).binary(binary_).rowmajor().O(output,MT);
      }
    } else {
      if (h_.matrixtype == IOMatrixtype_diagonal) {
	Matrix1int Mr;
	Matrix1<T> Mv;
	M.tolist(Mr,Mv);
	IOHelper(Mr).binary(binary_).colmajor().O(output,Mr);
	IOHelper(Mv).binary(binary_).colmajor().O(output,Mv);
      } else {
	Matrix1int Mr;
	Matrix1int Mc;
	Matrix1<T> Mv;
	M.tolist(Mr,Mc,Mv,(h_.matrixtype == IOMatrixtype_symmetric));
	IOHelper(Mr).binary(binary_).colmajor().O(output,Mr);
	IOHelper(Mc).binary(binary_).colmajor().O(output,Mc);
	IOHelper(Mv).binary(binary_).colmajor().O(output,Mv);
      }
    }
    return *this;
  }
  
  template <class T>
  IOHelper& IOHelper::I(std::istream& input,
			SparseMatrix<T>& M)
  {
    NOT_IMPLEMENTED;
    M.clear();
    if (h_.storagetype == IOStoragetype_rowmajor) {
      if (h_.matrixtype == IOMatrixtype_diagonal) {
	Matrix1< SparseMatrixDuplet<T> > MT;
	IOHelper(MT).binary(binary_).rowmajor().I(input,MT);
	M.fromlist(MT);
      } else {
	Matrix1< SparseMatrixTriplet<T> > MT;
	IOHelper(MT).binary(binary_).rowmajor().I(input,MT);
	M.fromlist(MT,(h_.matrixtype == IOMatrixtype_symmetric));
      }
    } else {
      if (h_.matrixtype == IOMatrixtype_diagonal) {
	Matrix1int Mr;
	Matrix1<T> Mv;
	IOHelper(Mr).binary(binary_).colmajor().I(input,Mr);
	IOHelper(Mv).binary(binary_).colmajor().I(input,Mv);
	M.fromlist(Mr,Mv);
      } else {
	Matrix1int Mr;
	Matrix1int Mc;
	Matrix1<T> Mv;
	IOHelper(Mr).binary(binary_).colmajor().I(input,Mr);
	IOHelper(Mc).binary(binary_).colmajor().I(input,Mc);
	IOHelper(Mv).binary(binary_).colmajor().I(input,Mv);
	M.fromlist(Mr,Mc,Mv,(h_.matrixtype == IOMatrixtype_symmetric));
      }
    }
    return *this;
  }
  


  
} /* namespace fmesh */

#endif
