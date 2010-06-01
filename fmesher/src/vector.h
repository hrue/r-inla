#ifndef _FMESH_VECTOR_
#define _FMESH_VECTOR_ 1

#include <cstddef>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <vector>
#include <set>
#include <map>
#include <list>
#include <string>
#include <cmath>

#ifndef NOT_IMPLEMENTED
#define NOT_IMPLEMENTED (std::cout					\
			 << __FILE__ << "(" << __LINE__ << ")\t"	\
			 << "NOT IMPLEMENTED: "				\
			 << __PRETTY_FUNCTION__ << std::endl);
#endif

namespace fmesh {



  template <class T> class Matrix;
  template <class T> class Vector3;
  template <class T> class Matrix3;
  template <class T> class SparseMatrix;

  template<class T>
  std::ostream& operator<<(std::ostream& output,
			   const Matrix<T>& M);
  template<class T>
  std::ostream& operator<<(std::ostream& output,
			   const SparseMatrix<T>& M);


  template <class T>
  class Matrix {
    friend
    std::ostream& operator<< <> (std::ostream& output,
				 const Matrix<T>& M);
    
  protected:
    static const size_t capacity_step_size_ = 1024;
    static const size_t capacity_doubling_limit_ = 8192;
    T* data_;
    size_t rows_;
    size_t cols_;
    size_t cap_;
    const T zero_;
  public:
    Matrix() : data_(NULL), rows_(0), cols_(0), cap_(0), zero_() {};
    Matrix(size_t set_cols) : data_(NULL), rows_(0), cols_(0),
			      cap_(0), zero_() {
      cols(set_cols);
    };
    Matrix(size_t set_rows, size_t set_cols, const T* vals = NULL)
      : data_(NULL), rows_(0), cols_(0), cap_(0), zero_() {
      cols(set_cols);
      capacity(set_rows);
      rows_ = set_rows;
      if (vals) {
	std::memcpy(data_,vals,sizeof(T)*rows_*cols_);
      }
    };
    ~Matrix() {
      if (data_)
	delete[] data_;
    };
    void clear(void) {
      if (data_) {
	delete[] data_;
	data_ = NULL;
      }
      cap_ = 0;
      rows_ = 0;
      cols_ = 0;
    };
    void zeros(const size_t from_row = 0) {
      for (size_t i=from_row*cols_; i<cap_*cols_; i++)
	data_[i] = zero_;
    }
    void truncate(const size_t new_rows) {
      if (new_rows < rows_)
	rows_ = new_rows;
    }

    size_t capacity() const { return cap_; };
    bool capacity(size_t cap) {
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

    bool append(const Matrix<T>& toappend) {
      if (cols_ != toappend.cols_) return false;
      if (!capacity(rows_+toappend.rows_)) return false;
      std::memcpy(data_+rows_*cols_,toappend.data_,
		  sizeof(T)*toappend.rows_*cols_);
      rows_ += toappend.rows_;
      return true;
    };

    int rows(void) const { return rows_; };
    int cols(void) const { return cols_; };
    int cols(size_t set_cols) {
      if (cols_ > 0) {
	/* Cannot alter number of columns */
	return cols_;
      }
      cols_ = set_cols;
      return cols_;
    };

    const T (* operator[](const int r) const) {
      if (r >= (int)rows_) {
	return NULL;
      }
      return &data_[r*cols_];
    };

    T* operator()(const int r) {
      if (r >= (int)rows_) {
	capacity(r+1);
	rows_ = r+1;
      }
      return &data_[r*cols_]; 
    };

    T& operator()(const int r, const int c) {
      return operator()(r)[c]; 
      /*
      if (r >= (int)rows_) {
	capacity(r+1);
	rows_ = r+1;
      }
      return data_[r*cols_+c]; 
      */
    };

    const T& operator()(const int r, const int c, const T& val) {
      operator()(r,c) = val; 
      return data_[r*cols_+c]; 
    };

    const T (* raw(void) const) { return data_; }
    T* raw(void) { return data_; }

  };


  template <class T>
  class Vector3 {
  public:
    typedef Vector3<T> selfT;
    typedef T Raw[3];
  private:
    T s[3];
  public:
    Vector3() {
      s[0] = T();
      s[1] = T();
      s[2] = T();
    };
    Vector3(const T& val0,
	    const T& val1,
	    const T& val2) {
      s[0] = val0;
      s[1] = val1;
      s[2] = val2;
    };
    Vector3(const T val[3]) {
      s[0] = val[0];
      s[1] = val[1];
      s[2] = val[2];
    };
    Vector3(const selfT& vec) {
      s[0] = vec.s[0];
      s[1] = vec.s[1];
      s[2] = vec.s[2];
    };

    selfT& operator=(const selfT vec) {
      return copy(vec);
    };

    const T& operator[](const int i) const {
      return s[i];
    };

    T& operator[](const int i) {
      return s[i];
    };

    const Raw& raw(void) const { return s; }
    


    selfT& copy(const selfT& s0)
    {
      if (this != &s0) {
	s[0] = s0.s[0];
	s[1] = s0.s[1];
	s[2] = s0.s[2];
      }
      return *this;
    };
    selfT& rescale(T s1)
    {
      s[0] *= s1;
      s[1] *= s1;
      s[2] *= s1;
      return *this;
    };
    selfT& scale(const selfT& s0, T s1)
    {
      s[0] = s0.s[0]*s1;
      s[1] = s0.s[1]*s1;
      s[2] = s0.s[2]*s1;
      return *this;
    };
    selfT& diff(const selfT& s0, const selfT& s1)
    {
      s[0] = s0.s[0]-s1.s[0];
      s[1] = s0.s[1]-s1.s[1];
      s[2] = s0.s[2]-s1.s[2];
      return *this;
    };
    selfT& sum(const selfT& s0, const selfT& s1)
    {
      s[0] = s0.s[0]+s1.s[0];
      s[1] = s0.s[1]+s1.s[1];
      s[2] = s0.s[2]+s1.s[2];
      return *this;
    };
    selfT& accum(const selfT& s0, T s1)
    {
      s[0] += s0.s[0]*s1;
      s[1] += s0.s[1]*s1;
      s[2] += s0.s[2]*s1;
      return *this;
    };
    T scalar(const selfT& s1) const
    {
      return (s[0]*s1.s[0]+s[1]*s1.s[1]+s[2]*s1.s[2]);
    };
    double length() const;
    selfT& cross(const selfT& s0,const selfT& s1)
    {
      if ((this == &s0) || (this == &s0)) {
	T s_0 = s0.s[1]*s1.s[2]-s0.s[2]*s1.s[1];
	T s_1 = s0.s[2]*s1.s[0]-s0.s[0]*s1.s[2];
	T s_2 = s0.s[0]*s1.s[1]-s0.s[1]*s1.s[0];
	s[0] = s_0;
	s[1] = s_1;
	s[2] = s_2;
      } else {
	s[0] = s0.s[1]*s1.s[2]-s0.s[2]*s1.s[1];
	s[1] = s0.s[2]*s1.s[0]-s0.s[0]*s1.s[2];
	s[2] = s0.s[0]*s1.s[1]-s0.s[1]*s1.s[0];
      }
    };
    T cross2(const selfT& s1) const
    {
      return (s[0]*s1.s[1]-s[1]*s1.s[0]);
    };
    /*!
      "Volume product" = scalar(cross(s0,s1),s2)
     */
    T volume(const selfT& s1, const selfT& s2) const
    {
      return ((s1.s[1]*s2.s[2]-s1.s[2]*s2.s[1])*s[0]+
	      (s1.s[2]*s2.s[0]-s1.s[0]*s2.s[2])*s[1]+
	      (s1.s[0]*s2.s[1]-s1.s[1]*s2.s[0])*s[2]);
    };
    double angle(const selfT& s0, const selfT& s1) const
    {
      Vector3<T> s0xs1;
      s0xs1.cross(s0,s1);
      return std::atan2((double)length(s0xs1),(double)scalar(s0,s1));
    };

    friend
    void arbitrary_perpendicular(Vector3<double>& n,
				 const Vector3<double>& v);
  };



  template <class T>
  class Matrix1 : public Matrix<T> {
  public:
    typedef T ValueRaw;
    Matrix1() : Matrix<T>(1) {};
    Matrix1(size_t set_rows, const ValueRaw* vals = NULL)
      : Matrix<T>(set_rows,1,(T*)vals) {};
    void clear(void) {
      Matrix<T>::clear();
      Matrix<T>::cols(1);
    };

    const T& operator[](const int r) const {
      if (r >= (int)Matrix<T>::rows_) {
	/* ERROR */
      }
      return Matrix<T>::data_[r];
    };

    T& operator()(const int r) {
      return Matrix<T>::operator()(r,0);
    };
    
  };


  template <class T>
  class Matrix3 : public Matrix<T> {
  public:
    typedef T ValueRaw[3];
    typedef Vector3<T> ValueRow;
    Matrix3() : Matrix<T>(3) {};
    Matrix3(const Matrix<T>& M) : Matrix<T>(3) {
      for (int r=0; r < M.rows(); r++) {
	for (int c=0; (c < 3) && (c < M.cols()); c++) {
	  operator()(r,c,M[r][c]);
	}
      }
    };
    Matrix3(size_t set_rows, const ValueRow* vals)
      : Matrix<T>(set_rows,3,(T*)vals) {};
    Matrix3(size_t set_rows, const ValueRaw* vals = NULL)
      : Matrix<T>(set_rows,3,(T*)vals) {};
    void clear(void) {
      Matrix<T>::clear();
      Matrix<T>::cols(3);
    };

    const ValueRow& operator[](const int r) const {
      return ((ValueRow*)Matrix<T>::data_)[r];
    };
    
    ValueRow& operator()(const int r) {
      return *(ValueRow*)(&Matrix<T>::operator()(r,0));
    };
    
    T& operator()(const int r, const int c) {
      return Matrix<T>::operator()(r,c);
    };
    
    const T& operator()(const int r, const int c, const T& val) {
      return Matrix<T>::operator()(r,c,val);
    };

  };


  typedef Matrix3<double>::ValueRaw PointRaw;
  typedef Matrix3<int>::ValueRaw Int3Raw;
  typedef Matrix3<double>::ValueRow Point;
  typedef Matrix3<int>::ValueRow Int3;
  typedef Matrix1<int> Matrix1int;
  typedef Matrix1<double> Matrix1double;
  typedef Matrix3<int> Matrix3int;
  typedef Matrix3<double> Matrix3double;




  template <class T>
  class SparseMatrixDuplet {
  public:
    int r;
    T value;
    SparseMatrixDuplet()
      : r(0), value(T()) {};
    SparseMatrixDuplet(int set_r, const T& set_value)
      : r(set_r), value(set_value) {};
  };
  template <class T>
  class SparseMatrixTriplet {
  public:
    int r;
    int c;
    T value;
    SparseMatrixTriplet()
      : r(0), c(0), value(T()) {};
    SparseMatrixTriplet(int set_r, int set_c, const T& set_value)
      : r(set_r), c(set_c), value(set_value) {};
  };



  template <class T>
  class SparseMatrix {
    friend
    std::ostream& operator<< <> (std::ostream& output,
				 const SparseMatrix<T>& M);
    
  public:
    typedef typename std::map<int, T> RowType;
    typedef typename std::map<int, RowType> DataType;
    typedef typename RowType::const_iterator ColConstIter;
    typedef typename DataType::const_iterator RowConstIter;
    typedef typename RowType::iterator ColIter;
    typedef typename DataType::iterator RowIter;

  private:
    DataType data_;
    const T zero_;
  public:
    SparseMatrix() : data_(), zero_() {};
    void clear() { data_.clear(); };

    int rows(void) const {
      if (data_.size() == 0)
	return 0;
      else {
	return data_.rbegin()->first+1;
      }
    };

    int cols(void) const {
      int cols_ = -1;
      int cols_row_;
      for (RowConstIter row = data_.begin();
	   row != data_.end();
	   row++) {
	if (row->second.size() > 0) {
	  cols_row_ = row->second.rbegin()->first;
	  if (cols_row_ > cols_)
	    cols_ = cols_row_;
	}
      }
      return cols_+1;
    };

    int nnz(void) const {
      int nnz_ = 0;
      for (RowConstIter row = data_.begin();
	   row != data_.end();
	   row++) {
	nnz_ += row->second.size();
       }
      return nnz_;
    };

    bool non_zero(const int r,  const int c) const {
      RowConstIter row;
      if ((row = data_.find(r)) != data_.end()) {
	return (row->second.find(c) != row->second.end());
      } else {
	return false; 
      }
    };

    const T& operator()(const int r,  const int c) const {
      RowConstIter row;
      if ((row = data_.find(r)) != data_.end()) {
	ColConstIter col;
	if ((col = row->second.find(c)) != row->second.end()) {
	  return col->second;
	} else {
	  return zero_; 
	}
      } else {
	return zero_; 
      }
    };

    const RowType& operator[](const int r) const {
      return data_[r];
    };

    RowType& operator()(const int r) {
      return data_[r];
    };

    RowConstIter begin() const {
      return data_.begin();
    };

    RowConstIter end() const {
      return data_.end();
    };



    T& operator()(const int r,  const int c) {
      return data_[r][c];
    };

    const T& operator()(const int r,  const int c, const T& val) {
      if (val == zero_) {
	RowIter row;
	if ((row = data_.find(r)) != data_.end()) {
	  ColIter col;
	  if ((col = row->second.find(c)) != row->second.end()) {
	    row->second.erase(col);
	    if (row->second.size() == 0) {
	      data_.erase(row);
	    }
	  }
	}
	return zero_;
      } else {
	data_[r][c] = val;
	return data_[r][c];
      }
    };

    /*! To list, assuming diagonal. */
    void tolist(Matrix1< SparseMatrixDuplet<T> >& MT) const {
      int elem = 0;
      for (RowConstIter row = data_.begin();
	  row != data_.end();
	  row++) {
	ColConstIter col;
	if ((col = row->second.find(row->first)) != row->second.end()) {
	  MT(elem) = SparseMatrixDuplet<T>(row->first,col->second);
	  elem++;
	}
      }
    };
    /*! To list, general or symmetric. */
    void tolist(Matrix1< SparseMatrixTriplet<T> >& MT,
		bool assume_symmetric = false) const {
      int elem = 0;
      for (RowConstIter row = data_.begin();
	  row != data_.end();
	  row++) {
	for (ColConstIter col = row->second.begin();
	    col != row->second.end();
	    col++) {
	  if ((!assume_symmetric) ||
	      (row->first <= col->first)) {
	    MT(elem) = SparseMatrixTriplet<T>(row->first,col->first,
					      col->second);
	    elem++;
	  }
	}
      }
    };
    /*! To list, assuming diagonal. */
    void tolist(Matrix1< int >& Tr,
		Matrix1< T >& Tv) const {
      int elem = 0;
      for (RowConstIter row = data_.begin();
	  row != data_.end();
	  row++) {
	ColConstIter col;
	if ((col = row->second.find(row->first)) != row->second.end()) {
	  Tr(elem) = row->first;
	  Tv(elem) = col->second;
	  elem++;
	}
      }
    };
    /*! To list, general or symmetric. */
    void tolist(Matrix1< int >& Tr,
		Matrix1< int >& Tc,
		Matrix1< T >& Tv,
		bool assume_symmetric = false) const {
      int elem = 0;
      for (RowConstIter row = data_.begin();
	  row != data_.end();
	  row++) {
	for (ColConstIter col = row->second.begin();
	    col != row->second.end();
	    col++) {
	  if ((!assume_symmetric) ||
	      (row->first <= col->first)) {
	    Tr(elem) = row->first;
	    Tc(elem) = col->first;
	    Tv(elem) = col->second;
	    elem++;
	  }
	}
      }
    };

    /*! From list, assuming diagonal. */
    void fromlist(const Matrix1< SparseMatrixDuplet<T> >& MT) {
      for (int i=0; i<MT.rows(); i++)
	operator()(MT[i].r,MT[i].r,MT[i].value);
    };
    /*! From list, general or symmetric. */
    void fromlist(const Matrix1< SparseMatrixTriplet<T> >& MT,
		bool assume_symmetric = false) {
      for (int i=0; i<MT.rows(); i++)
	operator()(MT[i].r,MT[i].c,MT[i].value);
    };
    /*! From list, assuming diagonal. */
    void fromlist(const Matrix1< int >& Tr,
		  const Matrix1< T >& Tv) {
      for (int i=0; i<Tr.rows(); i++)
	operator()(Tr[i],Tr[i],Tv[i]);
    };
    /*! From list, general or symmetric. */
    void fromlist(const Matrix1< int >& Tr,
		  const Matrix1< int >& Tc,
		  const Matrix1< T >& Tv,
		  bool assume_symmetric = false) {
      if (assume_symmetric) {
	for (int i=0; i<Tr.rows(); i++) {
	  operator()(Tr[i],Tc[i],Tv[i]);
	  operator()(Tc[i],Tr[i],Tv[i]);
	}
      } else {
	for (int i=0; i<Tr.rows(); i++) {
	  operator()(Tr[i],Tc[i],Tv[i]);
	}
      }
    };

  };


  template <class T>
  double Vector3<T>::length() const
  {
    return 0.0;
  };


  struct Vec {  
    static void copy(Point& s, const Point& s0)
    {
      s[0] = s0[0];
      s[1] = s0[1];
      s[2] = s0[2];
    };
    static void rescale(Point& s, double s1)
    {
      s[0] *= s1;
      s[1] *= s1;
      s[2] *= s1;
    };
    static void scale(Point& s, const Point& s0, double s1)
    {
      s[0] = s0[0]*s1;
      s[1] = s0[1]*s1;
      s[2] = s0[2]*s1;
    };
    static void diff(Point& s,const Point& s0, const Point& s1)
    {
      s[0] = s0[0]-s1[0];
      s[1] = s0[1]-s1[1];
      s[2] = s0[2]-s1[2];
    };
    static void sum(Point& s,const Point& s0, const Point& s1)
    {
      s[0] = s0[0]+s1[0];
      s[1] = s0[1]+s1[1];
      s[2] = s0[2]+s1[2];
    };
    static void accum(Point& s, const Point& s0, double s1 = 1.0)
    {
      s[0] += s0[0]*s1;
      s[1] += s0[1]*s1;
      s[2] += s0[2]*s1;
    };
    static double scalar(const Point& s0, const Point& s1)
    {
      return (s0[0]*s1[0]+s0[1]*s1[1]+s0[2]*s1[2]);
    };
    static double length(const Point& s0)
    {
      return (std::sqrt(s0[0]*s0[0]+s0[1]*s0[1]+s0[2]*s0[2]));
    };
    static void cross(Point& s, const Point& s0, const Point& s1)
    {
      s[0] = s0[1]*s1[2]-s0[2]*s1[1];
      s[1] = s0[2]*s1[0]-s0[0]*s1[2];
      s[2] = s0[0]*s1[1]-s0[1]*s1[0];
    };
    static double cross2(const Point& s0, const Point& s1)
    {
      return (s0[0]*s1[1]-s0[1]*s1[0]);
    };
    /*!
      "Volume product" = scalar(cross(s0,s1),s2)
     */
    static double volume(const Point& s0, const Point& s1, const Point& s2)
    {
      return ((s0[1]*s1[2]-s0[2]*s1[1])*s2[0]+
	      (s0[2]*s1[0]-s0[0]*s1[2])*s2[1]+
	      (s0[0]*s1[1]-s0[1]*s1[0])*s2[2]);
    };
    static double angle(const Point& s0, const Point& s1)
    {
      Point s0xs1;
      cross(s0xs1,s0,s1);
      return std::atan2(length(s0xs1),scalar(s0,s1));
    };
    /*!
      Calculate an arbitrary perpendicular vector.

      Michael M. Stark, Efficient Construction of Perpendicular
      Vectors without Branching, Journal of graphics, gpu, and game
      tools, Vol. 14, No. 1: 55-62, 2009
    */
#define ABS(X) std::fabs(X)
#define SIGNBIT(X) ((unsigned int)(std::signbit(X) != 0))
    static void arbitrary_perpendicular(Point& n, const Point& v)
    {
      const unsigned int uyx = SIGNBIT(ABS(v[0]) - ABS(v[1]));
      const unsigned int uzx = SIGNBIT(ABS(v[0]) - ABS(v[2]));
      const unsigned int uzy = SIGNBIT(ABS(v[1]) - ABS(v[2]));
      const unsigned int xm = uyx & uzx;
      const unsigned int ym = (1^xm) & uzy;
      const unsigned int zm = 1^(xm & ym);
      std::cout << uyx << ' ' << uzx << ' ' << uzy << std::endl;
      std::cout << xm << ' ' << ym << ' ' << zm << std::endl;
      n[0] =  zm*v[1] - ym*v[2];
      n[1] =  xm*v[2] - zm*v[0];
      n[2] =  ym*v[0] - xm*v[1];
    };
  };



  template<class T>
  std::ostream& operator<<(std::ostream& output,
			   const Matrix<T>& M)
  {
    output << M.rows_ << " "
	   << M.cols_ << std::endl;
    for (int r=0; r<(int)M.rows_; r++) {
      for (int c=0; c<(int)M.cols_; c++) {
	output << M.data_[r*M.cols_+c] << " ";
      }
      output << std::endl;
    }
    return output;
  }


  template<class T>
  std::ostream& operator<<(std::ostream& output,
			   const SparseMatrixDuplet<T>& MT)
  {
    output << MT.r << " "
	   << MT.value;
    return output;
  }

  template<class T>
  std::istream& operator>>(std::istream& input,
			   SparseMatrixDuplet<T>& MT)
  {
    input >> MT.r
	  >> MT.value;
    return input;
  }

  template<class T>
  std::ostream& operator<<(std::ostream& output,
			   const SparseMatrixTriplet<T>& MT)
  {
    output << MT.r << " "
	   << MT.c << " "
	   << MT.value;
    return output;
  }

  template<class T>
  std::istream& operator>>(std::istream& input,
			   SparseMatrixTriplet<T>& MT)
  {
    input >> MT.r
	  >> MT.c
	  >> MT.value;
    return input;
  }


  template<class T>
  std::ostream& operator<<(std::ostream& output,
			   const SparseMatrix<T>& M)
  {
    output << M.rows() << " "
	   << M.cols() << " "
	   << M.nnz() << std::endl;
    for (typename SparseMatrix<T>::RowConstIter row
	   = M.data_.begin();
	 row != M.data_.end();
	 row++) {
      for (typename SparseMatrix<T>::ColConstIter col
	     = row->second.begin();
	   col != row->second.end();
	   col++) {
	output << row->first << " "
	       << col->first << " "
	       << col->second << std::endl;
      }
    }
    return output;
  }


} /* namespace fmesh */

#endif
