#ifndef _FMESH_VECTOR_T_
#define _FMESH_VECTOR_T_ 1

#include "vector.hh"

namespace fmesh {



  template <class T>
  Matrix<T>::Matrix(size_t set_rows, size_t set_cols, const T* vals)
      : data_(NULL), rows_(0), cols_(0), cap_(0) {
      cols(set_cols);
      capacity(set_rows);
      rows_ = set_rows;
      if (vals) {
	std::memcpy(data_,vals,sizeof(T)*rows_*cols_);
      }
    };
  template <class T>
  
  const Matrix<T>& Matrix<T>::operator=(const Matrix<T>& from) {
    clear();
    cols(from.cols_);
    capacity(from.cap_);
    rows_ = from.rows_;
    if (data_) {
      std::memcpy(data_,from.data_,sizeof(T)*rows_*cols_);
      }
    return *this;
  };

  template <class T>
  bool Matrix<T>::capacity(size_t cap) {
    if (cap <= cap_) {
      return true;
    }
    size_t old_cap = cap_;
    if ((cap_==0) && 
	(cap < capacity_step_size_))
      cap_ = cap;
    while (cap > cap_) {
      if (cap_ < capacity_step_size_)
	cap_ = capacity_step_size_;
      else {
	if (cap_ < capacity_doubling_limit_)
	  cap_ *= 2;
	else
	  cap_ += capacity_step_size_;
      }
    }
    
    T* data_new_ = new T[cap_*cols_];
    
    if (data_) { /* Copy existing data: */
      std::memcpy(data_new_,data_,
		  sizeof(T)*old_cap*cols_);
      delete[] data_;
    }
    
    data_ = data_new_;
    zeros(old_cap);
    return true;
  };

  template <class T>
  bool Matrix<T>::append(const Matrix<T>& toappend) {
    if (cols_ != toappend.cols_) return false;
    if (!capacity(rows_+toappend.rows_)) return false;
    std::memcpy(data_+rows_*cols_,toappend.data_,
		sizeof(T)*toappend.rows_*cols_);
    rows_ += toappend.rows_;
    return true;
  };
  
  template <class T>
  Matrix<T>& Matrix<T>::rows(size_t set_rows) {
    if (set_rows>rows_)
      capacity(set_rows);
    else if (set_rows<rows_)
      truncate(set_rows);
    return *this;
  };
  template <class T>
  Matrix<T>& Matrix<T>::cols(size_t set_cols) {
    /* When set, cannot alter number of columns,
       unless we only need to expand a single row. */
    if ((cols_ > 0) &&
	(!((rows_<=1) && cols_<=set_cols))) {
      return *this;
    }
    if ((cols_>0) && (rows_>0)) { /* We already have some data, and
				     need to carefully make sure the
				     capacity is enough */
      /* Pre-data-size: cap_*cols_
	 Post-data_size needed cap_*set_cols
	 capacity is set as r*cols_
	 Requirement: r*cols_ >= cap_*set_cols
	 r = (cap_*set_cols)/cols_+1
      */
      capacity((cap_*set_cols)/cols_+1);
      cols_ = set_cols;
      cap_ = rows_; /* This makes sure we dont' overestimate the true cap_. */
    } else
      cols_ = set_cols;
    return *this;
  };





  template <class T>
  SparseMatrix<T> diag(const Matrix<T>& M1) {
    SparseMatrix<T> SM(M1.rows(),M1.rows());
    for (int i=0; i<M1.rows(); i++) {
      SM(i,i) = M1[i][0];
    }
    return SM;
  };

  template <class T>
  Matrix<T> diag(const SparseMatrix<T>& M1) {
    Matrix<T> M(M1.rows(),1);
    for (int i=0; ((i<M1.rows()) && (i<M1.cols())); i++) {
      M(i,0) = M1[i][i];
    }
    return M;
  };



  template <class T>
  SparseMatrix<T> operator*(const SparseMatrix<T>& M1,
			    const SparseMatrix<T>& M2)
  {
    SparseMatrix<T> M;
    int M1rows = M1.rows();
    int M2rows = M2.rows();
    M.cols(M2.cols()).rows(M1rows);
    for (int i=0; i<M1rows; i++) {
      SparseMatrixRow<T>& Mi = M(i);
      const SparseMatrixRow<T>& M1i = M1[i];
      if (M1i.size() > 0) {
	for (typename SparseMatrixRow<T>::ColCIter M1k = M1i.begin();
	     (M1k != M1i.end()) && (M1k->first < M2rows);
	     M1k++) {
	  int k = M1k->first;
	  const T& M1ik = M1i[k];
	  const SparseMatrixRow<T>& M2k = M2[k];
	  for (typename SparseMatrixRow<T>::ColCIter M2j = M2k.begin();
	       (M2j != M2k.end());
	       M2j++) {
	    Mi(M2j->first) += (M1ik * M2j->second);
	  }
	}
      }
    }
    return M;
  }

  template <class T>
  SparseMatrix<T> operator-(const SparseMatrix<T>& M1,
			    const SparseMatrix<T>& M2)
  {
    SparseMatrix<T> M(M1);
    for (int r=0; (r<M1.rows()) && (r<M2.rows()); r++) {
      SparseMatrixRow<T>& Mr = M(r);
      const SparseMatrixRow<T>& M2r = M2[r];
      for (typename SparseMatrixRow<T>::ColCIter c = M2r.begin();
	   (c != M2r.end()) && (c->first < M1.cols());
	   c++) {
	Mr(c->first) -= c->second;
      }
    }
    return M;
  }

  template <class T>
  SparseMatrix<T> inverse(const SparseMatrix<T>& M1,
			  bool diagonal)
  {
    SparseMatrix<T> M;
    M.cols(M1.cols()).rows(M1.rows());
    if (!diagonal) {
      /* NOT IMPLEMENTED */
      return M;
    }
    for (int r=0; (r<M1.rows()) && (r<M1.cols()); r++) {
      const T& val = M1[r][r];
      if (!(val==T()))
	M(r,r) = 1/val;
    }
    return M;
  }





  template <class T>
  const T fmesh::Matrix<T>::zero_ = T();
  template <class T>
  const T fmesh::SparseMatrixRow<T>::zero_ = T();
  template <class T>
  const T fmesh::SparseMatrix<T>::zero_ = T();


} /* namespace fmesh */

#endif
