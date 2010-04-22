#ifndef _FMESH_MESH_
#define _FMESH_MESH_ 1

#include <cstddef>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <set>
#include <map>

#define Mesh_V_capacity_step_size 128

namespace fmesh {

  typedef double Point[3];

  class Dart;
  class M3intO;
  class M3doubleO;
  
  class Mesh {
    friend class Dart;
    friend std::ostream& operator<<(std::ostream& output, const Mesh& M);
  private:
    size_t Vcap_;
    size_t Tcap_;
    size_t nV_;
    size_t nT_;
    bool use_TTi_;
    int (*TV_)[3];  /* TV[t]  : {v1,v2,v3} */
    int (*TT_)[3];  /* TT[t]  : {t1,t2,t3} */
    int (*TTi_)[3]; /* TTi[t] : {vi1,vi2,vi3},
		       t == TT[ TT[t][i] ][ TTi[t][i] ] */
    double (*S_)[3];
    
    Mesh& check_capacity(int nVc, int nTc);

  private:
    Mesh& rebuildTT();
    Mesh& rebuildTTi();
    
  public:
    Mesh(void) : Vcap_(0), Tcap_(0), nV_(0), nT_(0), use_TTi_(false),
      TV_(NULL), TT_(NULL), TTi_(NULL), S_(NULL) {};
    Mesh(size_t Vcapacity, bool use_TTi);
    Mesh(const Mesh& M) : Vcap_(0), Tcap_(0), nV_(0), nT_(0), use_TTi_(false),
      TV_(NULL), TT_(NULL), TTi_(NULL), S_(NULL) {
      *this = M;
    };
    Mesh& operator=(const Mesh& M) {
      clear();
      useTTi(M.use_TTi_);
      S_set(M.S_,M.nV_);
      TV_set(M.TV_,M.nT_);
    };
    ~Mesh();
    Mesh& clear();

    bool useTTi() const { return use_TTi_; }
    Mesh& useTTi(bool use_TTi);

    const size_t nV() { return nV_; };
    const size_t nT() { return nT_; };
    const int (*TV() const)[3] { return TV_; };
    const int (*TT() const)[3] { return TT_; };
    const int (*TTi() const)[3] { return TTi_; };
    const double (*S() const)[3] { return S_; };
    M3intO TVO() const;
    M3intO TTO() const;
    M3intO TTiO() const;
    M3doubleO SO() const;
    
    Mesh& S_set(double (*S)[3], int nV);
    Mesh& TV_set(int (*TV)[3], int nT); 
    
    void swapEdge(Dart& d);
    void splitEdge(Dart& d, const Point& p);
    void splitTriangle(Dart& d);
  };
  
  
  class Dart {
    friend std::ostream& operator<<(std::ostream& output, const Dart& d);
  private:
    const Mesh *M_;
    int vi_;
    int edir_;
    size_t t_;
    
  public:
    Dart(const Mesh& M, int t=0)
      : M_(&M), vi_(0), edir_(1), t_(t) {};
    Dart(const Dart& d) : M_(d.M_), vi_(d.vi_),
			  edir_(d.edir_), t_(d.t_) {};
    int vi() const { return vi_; };
    int edir() const { return edir_; };
    int t() const { return t_; };
    Dart& operator=(const Dart& d) {
      M_ = d.M_ ;
      vi_ = d.vi_;
      edir_ = d.edir_;
      t_ = d.t_;
    };
    bool operator==(const Dart& d) const {
      return ((d.t_ == t_) &&
	      (d.vi_ == vi_) &&
	      (d.edir_ == edir_));
    };
    bool operator<(const Dart& d) const {
      return ((d.t_ < t_) ||
	      ((d.t_ == t_) &&
	       ((d.edir_ < edir_) ||
		((d.edir_ == edir_) &&
		 (d.vi_ < vi_)))));
    };
    bool operator!=(const Dart& d) const {
      return !(d == *this);
    };

    /* Graph traversal algebra. */
    Dart& alpha0(void);
    Dart& alpha1(void);
    Dart& alpha2(void);
    Dart& orbit0(void);
    Dart& orbit1(void);
    Dart& orbit2(void);
    Dart& orbit0cw(void);
    Dart& orbit1cw(void);
    Dart& orbit2cw(void);

  };

  class MCtri {
  public:
    Dart d_;
    double value_;
    MCtri(Dart d, double value) : d_(d), value_(value) {};
    MCtri(const MCtri &T) : d_(T.d_), value_(T.value_) {};
    bool operator<(const MCtri& tb) {
      return ((value_ < tb.value_) ||
	      ((value_ == tb.value_) &&
	       (d_ < tb.d_)));
    }
  };

  typedef std::multiset<MCtri> MCtriSet;

  class MeshConstructor {
  private:
    Mesh *M_;
    MCtriSet Tset;

  public:
  };


  class M3intO {
    friend std::ostream& operator<<(std::ostream& output, const M3intO& MO);
  private:
    size_t n_;
    const int (*M_)[3];
  public:
    M3intO(const int (*M)[3],size_t n) : M_(M), n_(n) {};
  };

  class M3doubleO {
    friend std::ostream& operator<<(std::ostream& output, const M3doubleO& MO);
  private:
    size_t n_;
    const double (*M_)[3];
  public:
    M3doubleO(const double (*M)[3],size_t n) : M_(M), n_(n) {};
  };

} /* namespace fmesh */

#endif
