#ifndef _FMESH_VECTOR_
#define _FMESH_VECTOR_ 1

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

  typedef double Point[3];

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


  template <class ValueType> class SparseMatrix;

  template<class ValueType>
  std::ostream& operator<<(std::ostream& output,
			   const SparseMatrix<ValueType>& M);


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
    const ValueType zero_value_;
    DataType data_;
  public:
    SparseMatrix() : zero_value_(), data_() {};

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
	  return zero_value_; 
	}
      } else {
	return zero_value_; 
      }
    };

    void assign(const int r,  const int c, const ValueType& val) {
      if (val == zero_value_) {
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
