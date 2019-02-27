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

#ifndef WHEREAMI
#define WHEREAMI __FILE__ << "(" << __LINE__ << ")\t"
#endif

#ifndef LOG_
#define LOG_(msg) std::cout << WHEREAMI << msg;
#endif
#ifndef LOG
#ifdef DEBUG
#define LOG(msg) LOG_(msg)
#else
#define LOG(msg)
#endif
#endif

#ifndef NOT_IMPLEMENTED
#define NOT_IMPLEMENTED (std::cout					\
			 << __FILE__ << "(" << __LINE__ << ")\t"	\
			 << "NOT IMPLEMENTED: "				\
			 << __PRETTY_FUNCTION__ << std::endl);
#endif

namespace fmesh {

  /* Note: This definition really belongs in ioutils.hh, but is placed
     here to avoid chicken-and-egg definition problem. */
  /*! general/symmentric/diagonal */
  enum IOMatrixtype {IOMatrixtype_general=0,
                     IOMatrixtype_symmetric=1,
                     IOMatrixtype_diagonal=2};

  template <class T> class IOHelperM;
  template <class T> class IOHelperSM;
  template <class T> class Matrix;
  template <class T> class Vector3;
  template <class T> class Matrix3;
  template <class T> class SparseMatrixRow;
  template <class T> class SparseMatrix;

  template <class T>
  Matrix<T> operator*(const SparseMatrix<T>& M1,
		      double M2);

  template <class T>
  SparseMatrix<T> operator*(const SparseMatrix<T>& M1,
			    const SparseMatrix<T>& M2);
  template <class T>
  SparseMatrix<T> operator-(const SparseMatrix<T>& M1,
			    const SparseMatrix<T>& M2);
  template <class T>
  SparseMatrix<T> inverse(const SparseMatrix<T>& M1,
			  bool diagonal = false);
  template <class T>
  SparseMatrix<T> diag(const Matrix<T>& M1);
  template <class T>
  Matrix<T> diag(const SparseMatrix<T>& M1);



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
    static const T zero_;
    T* data_;
    size_t rows_;
    size_t cols_;
    size_t cap_;
  public:
    Matrix() : data_(NULL), rows_(0), cols_(0), cap_(0) {};
    Matrix(size_t set_cols) : data_(NULL), rows_(0), cols_(0),
			      cap_(0) {
      cols(set_cols);
    };
    Matrix(size_t set_rows, size_t set_cols, const T* vals = NULL);
    const Matrix<T>& operator=(const Matrix<T>& from);
    ~Matrix() {
      if (data_)
	delete[] data_;
    };
    Matrix<T>& clear(void) {
      if (data_) {
	delete[] data_;
	data_ = NULL;
      }
      cap_ = 0;
      rows_ = 0;
      cols_ = 0;
      return *this;
    };
    void zeros(const size_t from_row = 0,
	       const size_t num_rows = 0) {
      size_t num_rows_ = ((num_rows != 0)
			  ? ((num_rows>cap_) ? cap_ : num_rows)
			  : rows_);
      for (size_t i=from_row*cols_; i<num_rows_*cols_; i++)
	data_[i] = zero_;
    }

    size_t capacity() const { return cap_; };
    bool capacity(size_t cap);

    bool append(const Matrix<T>& toappend);

    Matrix<T>& rows(size_t set_rows);
    size_t rows(void) const { return rows_; };
    size_t cols(void) const { return cols_; };
    Matrix<T>& cols(size_t set_cols);

    const T (* operator[](const size_t r) const) {
      if (r >= rows_) {
	return NULL;
      }
      return &data_[r*cols_];
    };

    T* operator()(const size_t r) {
      if (r >= rows_) {
	rows(r+1);
      }
      return &data_[r*cols_]; 
    };

    T& operator()(const size_t r, const size_t c) {
      if (c >= cols_)
	cols(c+1);
      return operator()(r)[c]; 
    };

    const T& operator()(const size_t r, const size_t c, const T& val) {
      return (operator()(r,c) = val); 
    };

    const T (* raw(void) const) { return data_; }
    T* raw(void) { return data_; }


    /*! \brief Store the matrix in a file. */
    bool save(std::string filename,
	      IOMatrixtype matrixt = IOMatrixtype_general,
	      bool binary = true) const;
    /*! \brief Read a matrix from a file. */
    bool load(std::string filename, bool binary = true);
    /*! \brief Store the matrix in a file in old headerless ascii format. */
    bool save_ascii_2009(std::string filename,
	      IOMatrixtype matrixt = IOMatrixtype_general) const;
    /*! \brief Read a matrix from a file in old headerless ascii format. */
    bool load_ascii_2009(std::string filename);
    /*! \brief Read a matrix from a stream in old headerless ascii format. */
    void load_ascii_2009(std::istream& input);

    friend
    SparseMatrix<T> diag<T>(const Matrix<T>& M1);
    friend
    Matrix<T> diag<T>(const SparseMatrix<T>& M1);

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

    const T& operator[](const size_t i) const {
      return s[i];
    };

    T& operator[](const size_t i) {
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
      return *this;
    };
    T cross2(const selfT& s1) const
    {
      return (s[0]*s1.s[1]-s[1]*s1.s[0]);
    };
    /*!
      "Volume product" = scalar(s0,cross(s1,s2))
     */
    T volume(const selfT& s1, const selfT& s2) const
    {
      return ((s1.s[1]*s2.s[2]-s1.s[2]*s2.s[1])*s[0]+
	      (s1.s[2]*s2.s[0]-s1.s[0]*s2.s[2])*s[1]+
	      (s1.s[0]*s2.s[1]-s1.s[1]*s2.s[0])*s[2]);
    };
    double angle(const selfT& s1) const
    {
      Vector3<T> s0xs1;
      s0xs1.cross(*this,s1);
      return std::atan2((double)s0xs1.length(),(double)scalar(s1));
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
    Matrix1<T>& clear(void) {
      Matrix<T>::clear();
      Matrix<T>::cols(1);
      return *this;
    };

    const T& operator[](const size_t r) const {
      if (r >= Matrix<T>::rows_) {
	/* ERROR */
	LOG_("Error: Index out of bounds.")
      }
      return Matrix<T>::data_[r];
    };

    T& operator()(const size_t r) {
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
      for (size_t r=0; r < M.rows(); r++) {
	for (size_t c=0; (c < 3) && (c < M.cols()); c++) {
	  operator()(r,c,M[r][c]);
	}
      }
    };
    Matrix3(size_t set_rows, const ValueRow* vals)
      : Matrix<T>(set_rows,3,(T*)vals) {};
    Matrix3(size_t set_rows, const ValueRaw* vals = NULL)
      : Matrix<T>(set_rows,3,(T*)vals) {};
    Matrix3<T>& clear(void) {
      Matrix<T>::clear();
      Matrix<T>::cols(3);
      return *this;
    };

    const ValueRow& operator[](const size_t r) const {
      return ((ValueRow*)Matrix<T>::data_)[r];
    };
    
    ValueRow& operator()(const size_t r) {
      return *(ValueRow*)(&Matrix<T>::operator()(r,0));
    };
    
    T& operator()(const size_t r, const size_t c) {
      return Matrix<T>::operator()(r,c);
    };
    
    const T& operator()(const size_t r, const size_t c, const T& val) {
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
  class SparseMatrixRow {
    friend class SparseMatrix<T>;
  public:
    typedef typename std::map<int, T> DataType;
    typedef typename DataType::const_iterator ColCIter;
    typedef typename DataType::const_reverse_iterator ColCRIter;
    typedef typename DataType::iterator ColIter;
    typedef typename DataType::reverse_iterator ColRIter;
  protected:
    static const T zero_;
    SparseMatrix<T>* M_;
    DataType data_;
  public:
    SparseMatrixRow()
      : M_(NULL), data_() { };
    SparseMatrixRow(const SparseMatrixRow<T>& from)
      : M_(from.M_), data_(from.data_) { };
    SparseMatrixRow(SparseMatrix<T>* M)
      : M_(M), data_() { };

    size_t size() const { return data_.size(); };

    void cols(size_t set_cols) {
      ColRIter col;
      if (data_.size()>0) {
	for (col = data_.rbegin();
	     (col != data_.rend()) && (!(col->first < (int)set_cols));
	     col = data_.rbegin()) {
	  data_.erase((++col).base());
	}
      }
    };

    size_t nnz(int r, int matrixt = 0) const {
      size_t nnz_ = 0;
      if (matrixt == 2) {
	ColCIter col;
	if ((col = data_.find(r)) != data_.end()) {
	  nnz_ = 1;
	}
      } else if (matrixt == 1) {
	for (ColCIter c = data_.begin();
	     c != data_.end();
	     c++)
	  if (r <= c->first) nnz_++;
      } else {
	nnz_ = data_.size();
      }
      return nnz_;
    };

    int tolist(int offset,
	       int row,
	       Matrix1< SparseMatrixTriplet<T> >& MT,
	       int matrixt = 0) const {
      int elem = 0;
      if (matrixt == 2) {
	ColCIter col;
	if ((col = data_.find(row)) != data_.end()) {
	  MT(offset+elem) =
	    SparseMatrixTriplet<T>(row, row, col->second);
	  elem++;
	}
      } else {
	for (ColCIter col = data_.begin();
	     col != data_.end();
	     col++) {
	  if ((matrixt==0) ||
	      (row <= col->first)) {
	    MT(offset+elem) =
	      SparseMatrixTriplet<T>(row, col->first, col->second);
	    elem++;
	  }
	}
      }
      return elem;
    };
    /*! To list, general, symmetric, or diagonal. */
    int tolist(int offset,
	       int row,
	       Matrix1< int >& Tr,
	       Matrix1< int >& Tc,
	       Matrix1< T >& Tv,
	       int matrixt = 0) const {
      int elem = 0;
      if (matrixt==2) {
	ColCIter col;
	if ((col = data_.find(row)) != data_.end()) {
	  Tr(offset+elem) = row;
	  Tc(offset+elem) = row;
	  Tv(offset+elem) = col->second;
	  elem++;
	}
      } else {
	for (ColCIter col = data_.begin();
	     col != data_.end();
	     col++) {
	  if ((matrixt==0) ||
	      (row <= col->first)) {
	    Tr(offset+elem) = row;
	    Tc(offset+elem) = col->first;
	    Tv(offset+elem) = col->second;
	    elem++;
	  }
	}
      }
      return elem;
    };

    ColIter begin() { return data_.begin(); };
    ColRIter rbegin() { return data_.rbegin(); };
    ColIter end() { return data_.end(); };
    ColRIter rend() { return data_.rend(); };
    ColIter find(int c) { return data_.find(c); };

  public:

    ColCIter begin() const { return data_.begin(); };
    ColCRIter rbegin() const { return data_.rbegin(); };
    ColCIter end() const { return data_.end(); };
    ColCRIter rend() const { return data_.rend(); };
    ColCIter find(int c) const { return data_.find(c); };

    const T& operator[](const size_t c) const {
      if (!(c < M_->cols())) {
	/* Range error. */
	LOG_("Error: Column index out of bounds.")
	return zero_;
      }
      ColCIter col;
      if ((col = data_.find(c)) != data_.end())
	return col->second;
      return zero_;
    };
    
    T& operator()(const size_t c) {
      if (!(c < M_->cols()))
	M_->cols(c+1);
      return data_[c];
    };

    void erase(ColIter& col) {
      data_.erase(col);
    };
    
    void erase(size_t c) {
      ColIter col;
      if ((col = data_.find(c)) != data_.end())
	data_.erase(col);
    };
    
  };


  template <class T>
  class SparseMatrix {
    friend
    std::ostream& operator<< <> (std::ostream& output,
				 const SparseMatrix<T>& M);
    
  public:
    typedef typename fmesh::SparseMatrixRow<T> RowType;
    typedef typename std::vector<RowType> DataType;
    typedef typename RowType::ColCIter ColCIter;
    typedef typename RowType::ColCRIter ColCRIter;
    typedef typename RowType::ColIter ColIter;
    typedef typename RowType::ColRIter ColRIter;

  private:
    static const T zero_;
    size_t cols_;
    DataType data_;
  public:
    SparseMatrix(size_t set_rows = 0, size_t set_cols = 0)
      : cols_(set_cols), data_() {
      rows(set_rows);
    };
    SparseMatrix(const SparseMatrix<T>& from)
      : cols_(from.cols_), data_(from.data_) {
      for (size_t r=0; r<rows(); r++) {
	data_[r].M_ = this;
      }
      //      std::cout << "SM copy" << std::endl;
    };
    const SparseMatrix<T>& operator=(const SparseMatrix<T>& from) {
      cols_ = from.cols_;
      data_ = from.data_;
      for (size_t r=0; r<rows(); r++) {
	data_[r].M_ = this;
      }
      //      std::cout << "SM assignment" << std::endl;
      return *this;
    };
    SparseMatrix<T>& clear() {
      data_.clear();
      return *this;
    };

    SparseMatrix<T>& rows(size_t set_rows) {
      data_.resize(set_rows, RowType(this));
      return *this;
    };

    SparseMatrix<T>& cols(size_t set_cols) {
      if (!(cols_<set_cols)) {
	for (size_t row=0; row<rows(); row++) {
	  data_[row].cols(set_cols);
	}
      }
      cols_ = set_cols;
      return *this;
    };

    size_t rows(void) const {
      return data_.size();
    };

    size_t cols(void) const {
      return cols_;
    };

    size_t nnz(int matrixt = 0) const {
      size_t nnz_ = 0;
      for (size_t row=0; row<rows(); row++)
	nnz_ += data_[row].nnz(row,matrixt);
      return nnz_;
    };

    bool non_zero(const size_t r,  const size_t c) const {
      if (r < rows()) {
	return (data_[r].find(c) != data_[r].end());
      } else {
	return false; 
      }
    };

    const T& operator()(const size_t r,  const size_t c) const {
      if (r<rows()) {
	return data_[r][c];
      } else {
	return zero_; 
      }
    };

    const RowType& operator[](const size_t r) const {
      if (!(r<rows())) {
	/* Range error. */
	LOG_("Error: Row index out ouf bounds.")
      }
      return data_[r];
    };

    RowType& operator()(const size_t r) {
      if (!(r<rows()))
	rows(r+1); /* Expand. */
      return data_[r];
    };

    T& operator()(const size_t r,  const size_t c) {
      return operator()(r)(c);
    };

    const T& operator()(const size_t r,  const size_t c, const T& val) {
      if (val == zero_) {
	if (r < rows()) {
	  data_[r].erase(r);
	}
	return zero_;
      } else {
	return (operator()(r)(c) = val);
      }
    };


    /* Linear algebra */
    friend
    SparseMatrix<T> operator*<T>(const SparseMatrix<T>& M1,
				 const SparseMatrix<T>& M2);
    friend
    SparseMatrix<T> operator-<T>(const SparseMatrix<T>& M1,
				 const SparseMatrix<T>& M2);
    friend
    SparseMatrix<T> inverse<T>(const SparseMatrix<T>& M1,
			       bool diagonal);
    friend
    SparseMatrix<T> diag<T>(const Matrix<T>& M1);
    friend
    Matrix<T> diag<T>(const SparseMatrix<T>& M1);


    /*! To list, general, symmetric, or diagonal. */
    int tolist(Matrix1< SparseMatrixTriplet<T> >& MT,
	       int matrixt = 0) const {
      int elem = 0;
      for (size_t row=0; row<rows(); row++) {
	elem += data_[row].tolist(elem,row,MT,matrixt);
      }
      return elem;
    };
    /*! To list, general, symmetric, or diagonal. */
    int tolist(Matrix1< int >& Tr,
	       Matrix1< int >& Tc,
	       Matrix1< T >& Tv,
	       int matrixt = 0) const {
      int elem = 0;
      for (size_t row=0; row < rows(); row++)
	elem += data_[row].tolist(elem,row,Tr,Tc,Tv,matrixt);
      return elem;
    };

    /*! From list, general, symmetric, or diagonal. */
    void fromlist(const Matrix1< SparseMatrixTriplet<T> >& MT,
		  int matrixt = 0) {
      if (matrixt==1) {
	for (size_t i=0; i<MT.rows(); i++) {
	  operator()(MT[i].r,MT[i].c,MT[i].value);
	  operator()(MT[i].c,MT[i].r,MT[i].value);
	}
      } else if (matrixt==2) {
	for (size_t i=0; i<MT.rows(); i++)
	  operator()(MT[i].r,MT[i].r,MT[i].value);
      } else {
	for (size_t i=0; i<MT.rows(); i++)
	  operator()(MT[i].r,MT[i].c,MT[i].value);
      }
    };
    /*! From list, general or symmetric. */
    void fromlist(const Matrix1< int >& Tr,
		  const Matrix1< int >& Tc,
		  const Matrix1< T >& Tv,
		  int matrixt = 0) {
      if (matrixt==1) {
	for (size_t i=0; i<Tr.rows(); i++) {
	  operator()(Tr[i],Tc[i],Tv[i]);
	  operator()(Tc[i],Tr[i],Tv[i]);
	}
      } else if (matrixt==1) {
	for (size_t i=0; i<Tr.rows(); i++) {
	  operator()(Tr[i],Tr[i],Tv[i]);
	}
      } else {
	for (size_t i=0; i<Tr.rows(); i++) {
	  operator()(Tr[i],Tc[i],Tv[i]);
	}
      }
    };

    /*! \brief Store the matrix in a file. */
    bool save(std::string filename,
	      IOMatrixtype matrixt = IOMatrixtype_general,
	      bool binary = true) const;
    /*! \brief Read a matrix from a file. */
    bool load(std::string filename, bool binary = true);
    /*! \brief Store the matrix in a file in old headerless ascii format. */
    bool save_ascii_2009(std::string filename,
			 IOMatrixtype matrixt = IOMatrixtype_general) const;
    
  };


  struct Vec {  
    static void copy(Point& s, const Point& s0)
    { s.copy(s0); };
    static void rescale(Point& s, double s1)
    { s.rescale(s1); };
    static void scale(Point& s, const Point& s0, double s1)
    { s.scale(s0,s1); };
    static void diff(Point& s,const Point& s0, const Point& s1)
    { s.diff(s0,s1); };
    static void sum(Point& s,const Point& s0, const Point& s1)
    { s.sum(s0,s1); };
    static void accum(Point& s, const Point& s0, double s1 = 1.0)
    { s.accum(s0,s1); };
    static double scalar(const Point& s0, const Point& s1)
    { return s0.scalar(s1); };
    static double length(const Point& s0);
    static void cross(Point& s, const Point& s0, const Point& s1)
    { s.cross(s0,s1); };
    static double cross2(const Point& s0, const Point& s1)
    { return s0.cross2(s1); };
    static double volume(const Point& s0, const Point& s1, const Point& s2)
    { return s0.volume(s1,s2); };
    static double angle(const Point& s0, const Point& s1)
    { return s0.angle(s1); };
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

  template <class T>
  double Vector3<T>::length() const
  {
    return (std::sqrt(s[0]*s[0]+s[1]*s[1]+s[2]*s[2]));
  }




  template<class T>
  std::ostream& operator<<(std::ostream& output,
			   const Matrix<T>& M)
  {
    output << M.rows_ << " "
	   << M.cols_ << std::endl;
    for (size_t r=0; r<M.rows(); r++) {
      for (size_t c=0; c<M.cols(); c++) {
	output << M.data_[r*M.cols()+c] << " ";
      }
      output << std::endl;
    }
    return output;
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
    for (size_t row=0; row<M.rows(); row++)
      for (typename SparseMatrix<T>::ColCIter col
	     = M[row].begin();
	   col != M[row].end();
	   col++) {
	output << row << " "
	       << col->first << " "
	       << col->second << std::endl;
      }
    return output;
  }









} /* namespace fmesh */

#include "vector.tcc"

#endif
