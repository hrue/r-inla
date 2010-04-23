#ifndef _FMESH_MESH_
#define _FMESH_MESH_ 1

#include <cstddef>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <vector>
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
    
  private:
    Mesh& rebuildTT();
    Mesh& rebuildTTi();
    
  public:
    Mesh(void) : Vcap_(0), Tcap_(0), nV_(0), nT_(0), use_TTi_(false),
      TV_(NULL), TT_(NULL), TTi_(NULL), S_(NULL) {};
    Mesh(size_t Vcapacity, bool use_TTi=false);
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

    /*!
      \brief Check the storage capacity, and increase if necessary
    */
    Mesh& check_capacity(int nVc, int nTc);

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
    
    Mesh& S_set(const double (*S)[3], int nV);
    Mesh& TV_set(const int (*TV)[3], int nT); 
    Mesh& S_append(const double (*S)[3], int nV);
    Mesh& TV_append(const int (*TV)[3], int nT); 
    
    Dart swapEdge(const Dart& d);
    Dart splitEdge(const Dart& d, int v);
    Dart splitTriangle(const Dart& d, int v);
  };
  
  /*! \breif Darts */
  class Dart {
    friend std::ostream& operator<<(std::ostream& output, const Dart& d);
  private:
    const Mesh *M_;
    int vi_;
    int edir_;
    size_t t_;
    
  public:
    /*! Test 1 */
    Dart(void)
      : M_(NULL), vi_(0), edir_(1), t_(0) {};
    Dart(const Mesh& M, int t=0, int edir=1, size_t vi=0)
      : M_(&M), vi_(vi), edir_(edir), t_(t) {};/*!< Test 2 */
    Dart(const Dart& d) : M_(d.M_), vi_(d.vi_),
			  edir_(d.edir_), t_(d.t_) {};
    /*!< Test 3 */
    int vi() const { return vi_; };
    int edir() const { return edir_; };
    int t() const { return t_; };
    Dart& operator=(const Dart& d) {
      M_ = d.M_ ;
      vi_ = d.vi_;
      edir_ = d.edir_;
      t_ = d.t_;
    };
    bool isnull() const { return (!M_) || (edir_==0); };
    bool operator==(const Dart& d) const {
      return ((d.t_ == t_) &&
	      (d.vi_ == vi_) &&
	      (d.edir_ == edir_));
    };
    bool operator<(const Dart& d) const {
      /* TODO: Add debug check for M_==d.M_ */
      return ((d.t_ < t_) ||
	      ((d.t_ == t_) &&
	       ((d.edir_ < edir_) ||
		((d.edir_ == edir_) &&
		 (d.vi_ < vi_)))));
    };
    bool operator!=(const Dart& d) const {
      return !(d == *this);
    };

    bool onBoundary() const {
      return (M_->TT_[t_][(vi_+2)%3] < 0);
    }

    /* Graph traversal algebra. */
    Dart& alpha0(void);
    Dart& alpha1(void);
    Dart& alpha2(void);
    Dart& orbit0(void);
    Dart& orbit1(void);
    Dart& orbit2(void);
    Dart& orbit0rev(void);
    Dart& orbit1rev(void);
    Dart& orbit2rev(void);

  };

  class MCdv {
  public:
    Dart d_;
    double value_;
    MCdv(Dart d, double value) : d_(d), value_(value) {};
    MCdv(const MCdv &T) : d_(T.d_), value_(T.value_) {};
    bool operator<(const MCdv& tb) {
      return ((value_ < tb.value_) ||
	      ((value_ == tb.value_) &&
	       (d_ < tb.d_)));
    }
  };

  typedef std::set<MCdv> MCdvSet;
  typedef std::set<Dart,double> MCdvMap;

  /*!
    \brief Class for constructing Delaunay triangulations
  */
  class MeshConstructor {
  private:
    enum State {State_noT, State_DT, State_CDT, State_refinedCDT};
    Mesh *M_;
    MCdvMap segm_; /*!< Segment darts, mapped to metadata */
    MCdvSet segm_encr_; /*!< Set of encroached segment darts */
    MCdvMap skinny_; /*!< Skinny triangles, mapped to metadata */
    MCdvSet skinny_sorted_;
    MCdvMap big_;
    MCdvSet big_sorted_;
    double skinny_limit_;
    double big_limit_;
    State state_;

  public:
    MeshConstructor() : M_(NULL), state_(State_noT) {};

    Dart locatePoint(const Dart& d, const Point s) {return Dart();};
    /*!
      \brief Append vertices

      Return index of the first of the added points.
    */
    int addVertices(const Dart& d, const double (*S)[3], int nV)
    {
      M_->S_append(S,nV);
      return M_->nV()-nV;
    };
    Dart swapEdge(const Dart& d);
    Dart splitEdge(const Dart& d, int v);
    Dart splitTriangle(const Dart& d, int v);

    /*!
      \brief Build Delaunay triangulation (DT)

      If the input mesh contains triangles, it is assumed to be a DT
      including at least the boundary vertices of the convex hull.

      If PruneExterior is to be used later, any exterior points must
      be at the end of the vertex list.
     */
    void DT() {};
    /*!
      \brief Add segments to constraint list, preparing for CDT
    */
    void addSegments() {};
    /*!
      \brief Build constrained Delaunay triangulation (CDT)
    */
    void CDT() {};
    /*!
      \brief Remove exterior triangles from a CDT
    */
    void PruneExterior() {};
    /*!
      \brief Refine a CDT
    */
    void RefineCDT() {};
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
