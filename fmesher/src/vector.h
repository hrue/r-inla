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



  template <class ValueType> class Matrix;
  template <class ValueType> class Vector3;
  template <class ValueType> class Matrix3;
  template <class ValueType> class SparseMatrix;

  template<class ValueType>
  std::ostream& operator<<(std::ostream& output,
			   const Matrix<ValueType>& M);
  template<class ValueType>
  std::ostream& operator<<(std::ostream& output,
			   const SparseMatrix<ValueType>& M);


  template <class ValueType>
  class Matrix {
    friend
    std::ostream& operator<< <> (std::ostream& output,
				 const Matrix<ValueType>& M);
    
  protected:
    static const size_t capacity_step_size_ = 1024;
    static const size_t capacity_doubling_limit_ = 8192;
    ValueType* data_;
    size_t rows_;
    size_t cols_;
    size_t cap_;
  public:
    Matrix() : data_(NULL), rows_(0), cols_(0), cap_(0) {};
    Matrix(size_t set_cols) : data_(NULL), rows_(0), cols_(0), cap_(0) {
      cols(set_cols);
    };
    Matrix(size_t set_rows, size_t set_cols, const ValueType* vals = NULL)
      : data_(NULL), rows_(0), cols_(0), cap_(0) {
      cols(set_cols);
      capacity(set_rows);
      rows_ = set_rows;
      if (vals) {
	std::memcpy(data_,vals,sizeof(ValueType)*rows_*cols_);
      }
    };
    ~Matrix() {
      if (data_)
	delete[] data_;
      if (cap_>0)
	std::cout << "Capability = " << cap_
		  << ", rows = " << rows_
		  << ", cols = " << cols_
		  << std::endl;
    };
    void clear(void) {
      if (cap_>0)
	std::cout << "Capability = " << cap_
		  << ", rows = " << rows_
		  << ", cols = " << cols_
		  << std::endl;
      if (data_) {
	delete[] data_;
	data_ = NULL;
      }
      cap_ = 0;
      rows_ = 0;
      cols_ = 0;
    };
    void truncate(const size_t new_rows) {
      if (new_rows < rows_)
	rows_ = new_rows;
    }

    size_t capacity() const { return cap_; };
    bool capacity(size_t cap) {
      if (cap <= cap_) {
	return true;
      }
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

      ValueType* data_new_ = new ValueType[cap_*cols_];

      if ((data_) && (rows_>0)) { /* Copy existing data: */
	std::memcpy(data_new_,data_,
		    sizeof(ValueType)*rows_*cols_);
	delete[] data_;
      }

      data_ = data_new_;
      return true;
    };

    bool append(const Matrix<ValueType>& toappend) {
      if (cols_ != toappend.cols_) return false;
      if (!capacity(rows_+toappend.rows_)) return false;
      std::memcpy(data_+rows_*cols_,toappend.data_,
		  sizeof(ValueType)*toappend.rows_*cols_);
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

    const ValueType (* operator[](const int r) const) {
      if (r >= (int)rows_) {
	return NULL;
      }
      return &data_[r*cols_];
    };

    ValueType& operator()(const int r, const int c) {
      if (r >= (int)rows_) {
	capacity(r+1);
	rows_ = r+1;
      }
      return data_[r*cols_+c]; 
    };

  };


  template <class ValueType>
  class Vector3 {
  public:
    typedef ValueType Raw[3];
  private:
    ValueType data_[3];
  public:
    Vector3() {};
    Vector3(const ValueType& val0,
	    const ValueType& val1,
	    const ValueType& val2) {
      data_[0] = val0;
      data_[1] = val1;
      data_[2] = val2;
    };
    Vector3(const ValueType val[3]) {
      data_[0] = val[0];
      data_[1] = val[1];
      data_[2] = val[2];
    };
    Vector3(const Vector3<ValueType>& vec) {
      data_[0] = vec.data_[0];
      data_[1] = vec.data_[1];
      data_[2] = vec.data_[2];
    };

    Vector3<ValueType>& operator=(const Vector3<ValueType> vec) {
      if (this != &vec) {
	data_[0] = vec.data_[0];
	data_[1] = vec.data_[1];
	data_[2] = vec.data_[2];
      }
      return *this;
    };

    const ValueType& operator[](const int i) const {
      return data_[i];
    };

    ValueType& operator[](const int i) {
      return data_[i];
    };

    const Raw& raw(void) const { return data_; }
    
  };

  template <class ValueType>
  class Matrix1 : public Matrix<ValueType> {
  public:
    typedef ValueType ValueRaw;
    Matrix1() : Matrix<ValueType>(1) {};
    Matrix1(size_t set_rows, const ValueRaw* vals = NULL)
      : Matrix<ValueType>(set_rows,1,(ValueType*)vals) {};
    void clear(void) {
      Matrix<ValueType>::clear();
      Matrix<ValueType>::cols(1);
    };

    const ValueType& operator[](const int r) const {
      if (r >= (int)Matrix<ValueType>::rows_) {
	/* ERROR */
      }
      return Matrix<ValueType>::data_[r];
    };

    ValueType& operator()(const int r) {
      return Matrix<ValueType>::operator()(r,0);
    };
    
  };


  template <class ValueType>
  class Matrix3 : public Matrix<ValueType> {
  public:
    typedef ValueType ValueRaw[3];
    typedef Vector3<ValueType> ValueRow;
    Matrix3() : Matrix<ValueType>(3) {};
    Matrix3(size_t set_rows, const ValueRow* vals)
      : Matrix<ValueType>(set_rows,3,(ValueType*)vals) {};
    Matrix3(size_t set_rows, const ValueRaw* vals = NULL)
      : Matrix<ValueType>(set_rows,3,(ValueType*)vals) {};
    void clear(void) {
      Matrix<ValueType>::clear();
      Matrix<ValueType>::cols(3);
    };

    const ValueRow& operator[](const int r) const {
      return ((ValueRow*)Matrix<ValueType>::data_)[r];
    };
    
    ValueRow& operator()(const int r) {
      return *(ValueRow*)(&Matrix<ValueType>::operator()(r,0));
    };
    
    ValueType& operator()(const int r, const int c) {
      return Matrix<ValueType>::operator()(r,c);
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







  template <class ValueType>
  class SparseMatrix {
    friend
    std::ostream& operator<< <> (std::ostream& output,
				 const SparseMatrix<ValueType>& M);
    
    typedef typename std::map<int, ValueType> RowType;
    typedef typename std::map<int, RowType> DataType;
    typedef typename RowType::const_iterator ColConstIter;
    typedef typename DataType::const_iterator RowConstIter;
    typedef typename RowType::iterator ColIter;
    typedef typename DataType::iterator RowIter;
  private:
    static const ValueType zero_ = ValueType();
    DataType data_;
  public:
    SparseMatrix() : data_() {};

    int rows(void) const {
      if (data_.size() == 0)
	return 0;
      else {
	return data_.rbegin()->first;
      }
    };

    int cols(void) const {
      int cols_ = 0;
      int cols_row_;
       for (RowConstIter row = data_.begin();
	    row != data_.end();
	    row++) {
	 cols_row_ = row->second.rbegin()->first;
	 if (cols_row_ > cols_)
	   cols_ = cols_row_;
       }
      return cols_;
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

    const ValueType& operator()(const int r,  const int c) const {
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

    ValueType& operator()(const int r,  const int c) {
      return data_[r][c];
    };

    void assign(const int r,  const int c, const ValueType& val) {
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
      } else {
	data_[r][c] = val;
      }
    };

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



  template<class ValueType>
  std::ostream& operator<<(std::ostream& output,
			   const Matrix<ValueType>& M)
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


  template<class ValueType>
  std::ostream& operator<<(std::ostream& output,
			   const SparseMatrix<ValueType>& M)
  {
    output << M.rows() << " "
	   << M.cols() << " "
	   << M.nnz() << std::endl;
    for (typename SparseMatrix<ValueType>::RowConstIter row
	   = M.data_.begin();
	 row != M.data_.end();
	 row++) {
      for (typename SparseMatrix<ValueType>::ColConstIter col
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
