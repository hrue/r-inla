#ifndef _FMESH_MESH_
#define _FMESH_MESH_ 1

#include <cstddef>
#include <cstddef>
//#include <cstring>
#include <iostream>
#include <iomanip>
#include <vector>
#include <set>
#include <map>
#include <list>
#include <string>

#include "xtmpl.h"

#define Mesh_V_capacity_step_size 128
#define MESH_EPSILON 1e-18

#define NOT_IMPLEMENTED (std::cout << "Not implemented: \""	\
			 << __PRETTY_FUNCTION__ << std::endl);

namespace fmesh {

  typedef double Point[3];
  typedef std::list<int> vertexListT;
  typedef std::list<int> triangleListT;
  typedef std::pair<int,int> constrT;
  typedef std::list<constrT> constrListT;

  class Xtmpl;
  class Mesh;
  class Dart;
  class MOAint;
  class MOAint3;
  class MOAdouble3;
  class MeshC;
  
  class Mesh {
    friend class Dart;
    friend std::ostream& operator<<(std::ostream& output, const Mesh& M);
  public:
    enum Mtype {Mtype_manifold=0,
		Mtype_plane,
		Mtype_sphere};
  private:
    Mtype type_;
    size_t Vcap_;
    size_t Tcap_;
    size_t nV_;
    size_t nT_;
    bool use_VT_;
    bool use_TTi_;
    int (*TV_)[3];  /* TV[t]  : {v1,v2,v3} */
    int (*TT_)[3];  /* TT[t]  : {t1,t2,t3} */
    int (*VT_);     /* VT[v]  : t,
		       v == TV[t][vi]  for some vi=0,1,2 */
    int (*TTi_)[3]; /* TTi[t] : {vi1,vi2,vi3},
		       t == TT[ TT[t][i] ][ TTi[t][i] ] */
    double (*S_)[3];
    Xtmpl (*X11_);
    
  private:
    Mesh& rebuildTT();

    Mesh& updateVT(const int v, const int t);
    /*!< Change VT[v] only if not linked to a triangle */
    Mesh& setVT(const int v, const int t);
    /* Overwerite current VT[v] info */
    Mesh& updateVTtri(const int t);
    Mesh& setVTtri(const int t);
    Mesh& updateVTtri_private(const int t0);
    Mesh& setVTv_private(const int t0);

    Mesh& rebuildVT();
    Mesh& rebuildTTi();

    void redrawX11(std::string str);
    
  public:
    Mesh(void) : type_(Mtype_manifold), Vcap_(0), Tcap_(0),
      nV_(0), nT_(0), use_VT_(false), use_TTi_(true),
		 TV_(NULL), TT_(NULL), TTi_(NULL), S_(NULL), X11_(NULL) {};
    Mesh(Mtype manifold_type, size_t Vcapacity, bool use_VT=true, bool use_TTi=false);
    Mesh(const Mesh& M) : type_(Mtype_manifold), Vcap_(0), Tcap_(0),
      nV_(0), nT_(0), use_VT_(true), use_TTi_(false),
      TV_(NULL), TT_(NULL), TTi_(NULL), S_(NULL), X11_(NULL) {
      *this = M;
    };
    Mesh& operator=(const Mesh& M);
    ~Mesh();
    Mesh& clear();

    /*!
      \brief Check the storage capacity, and increase if necessary
    */
    Mesh& check_capacity(size_t nVc, size_t nTc);

    bool useVT() const { return use_VT_; };
    Mesh& useVT(bool use_VT);
    bool useTTi() const { return use_TTi_; };
    Mesh& useTTi(bool use_TTi);

    bool useX11() const { return (X11_!=NULL); };
    Mesh& useX11(bool use_X11,
		 int sx = 500, int sy = 500,
		 double minx = -0.05,
		 double maxx = 1.05,
		 double miny = -0.05,
		 double maxy = 1.05,
		 std::string name = "fmesher::Mesh");

    Mtype type() const { return type_; };
    size_t nV() const { return nV_; };
    size_t nT() const { return nT_; };
    const int (*TV() const)[3] { return TV_; };
    const int (*TT() const)[3] { return TT_; };
    const int (*VT() const) { return VT_; };
    const int (*TTi() const)[3] { return TTi_; };
    const double (*S() const)[3] { return S_; };
    Xtmpl *X11() { return X11_; };
    MOAint3 TVO() const;
    MOAint3 TTO() const;
    MOAint VTO() const;
    MOAint3 TTiO() const;
    MOAdouble3 SO() const;
    
    Mesh& S_set(const double (*S)[3], int nV);
    Mesh& TV_set(const int (*TV)[3], int nT); 
    Mesh& S_append(const double (*S)[3], int nV);
    Mesh& TV_append(const int (*TV)[3], int nT); 

    Dart locatePoint(const Dart& d0, const Point s, double* delta_min) const;
    Dart locateVertex(const Dart& d0, const int v) const;
    
    Dart swapEdge(const Dart& d);
    Dart splitEdge(const Dart& d, int v);
    Dart splitTriangle(const Dart& d, int v);

    /* Traits: */
    double edgeLength(const Dart& d) const;
    double triangleArea(int t) const;
    double triangleCircumcircleRadius(int t) const;
    double triangleShortestEdge(int t) const;
    double triangleLongestEdge(int t) const;
    double edgeEncroached(const Dart& d, const double s[3]) const;
    
    double encroachedQuality(const Dart& d) const;
    double skinnyQuality(int t) const;
    double bigQuality(int t) const;

    /*!
      Compute dart half-space test for a point.
      positive if s is to the left of the edge defined by d.
     */
    double inLeftHalfspace(const Dart& d, const double s[3]) const;
    double inCircumcircle(const Dart& d, const double s[3]) const;
    bool circumcircleOK(const Dart& d) const;
  };



  class MOAint {
    friend std::ostream& operator<<(std::ostream& output, const MOAint& MO);
  private:
    size_t n_;
    const int (*M_);
  public:
    MOAint(const int (*M),size_t n) : n_(n), M_(M) {};
  };

  class MOAint3 {
    friend std::ostream& operator<<(std::ostream& output, const MOAint3& MO);
  private:
    size_t n_;
    const int (*M_)[3];
  public:
    MOAint3(const int (*M)[3],size_t n) : n_(n), M_(M) {};
  };

  class MOAdouble3 {
    friend std::ostream& operator<<(std::ostream& output, const MOAdouble3& MO);
  private:
    size_t n_;
    const double (*M_)[3];
  public:
   MOAdouble3(const double (*M)[3],size_t n) : n_(n), M_(M) {};
  };



  
  /*! \breif Darts */
  class Dart {
    friend std::ostream& operator<<(std::ostream& output, const Dart& d);
  private:
    const Mesh *M_;
    size_t vi_;
    int edir_;
    int t_;
    
  public:
    Dart(void)
      : M_(NULL), vi_(0), edir_(1), t_(0) {};
    Dart(const Mesh& M, int t=0, int edir=1, size_t vi=0)
      : M_(&M), vi_(vi), edir_(edir), t_(t) {};/*!< Test 2 */
    Dart(const Dart& d) : M_(d.M_), vi_(d.vi_),
			  edir_(d.edir_), t_(d.t_) {};
    Dart& operator=(const Dart& d) {
      M_ = d.M_ ;
      vi_ = d.vi_;
      edir_ = d.edir_;
      t_ = d.t_;
      return *this;
    };

    const Mesh* M() const { return M_; };
    int vi() const { return vi_; };
    int edir() const { return edir_; };
    int t() const { return t_; };

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
      return (M_->TT_[t_][(vi_+(3-edir_))%3] < 0);
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














  class MCQdv {
  public:
    Dart d_;
    double value_;
    MCQdv(Dart d, double value) : d_(d), value_(value) {};
    MCQdv(const MCQdv &T) : d_(T.d_), value_(T.value_) {};
    bool operator<(const MCQdv& tb) const {
      return ((value_ < tb.value_) ||
	      ((value_ == tb.value_) &&
	       (d_ < tb.d_)));
    }
  };

  class MCQ {
  protected:
    typedef std::map<Dart,double> map_type;
    typedef map_type::value_type map_key_type;
    typedef std::set<MCQdv> set_type;
    map_type darts_; /*!< Darts, mapped to quality */
    set_type darts_quality_;
    /*!< Set of "bad quality" (calcQ>0.0) segment darts */
    bool only_quality_;
    /*!< If true, only store darts that are of bad quality */
  public:
    MCQ(bool only_quality) : darts_(), darts_quality_(),
			     only_quality_(only_quality) {};
    virtual double calcQ(const Dart& d) const = 0;
    void clear() {
      darts_.clear();
      darts_quality_.clear();
    };
    void insert(const Dart& d); /*!< Insert dart if not existing. */
    void erase(const Dart& d); /*!< Remove dart if existing. */
    int count() const { return darts_.size(); };
    int countQ() const { return darts_quality_.size(); };
    bool empty() const { return darts_.empty(); };
    bool emptyQ() const { return darts_quality_.empty(); };
    bool found(const Dart& d) const;
    bool foundQ(const Dart& d) const;
    const double quality(const Dart& d) const;
    Dart quality() const;

    friend std::ostream& operator<<(std::ostream& output, const MCQ& Q);
  };

  class MCQtri : public MCQ {
  protected:
    double quality_limit_;
    /*!< Larger values are included in the quality set */
    virtual double calcQtri(const Dart& d) const = 0;
  public:
    MCQtri(bool only_quality, double quality_limit);
    void setQ(double quality_limit);
    void insert(const Dart& d) {
      MCQ::insert(Dart(*d.M(),d.t()));
    };
    void erase(const Dart& d) {
      MCQ::erase(Dart(*d.M(),d.t()));
    };
    bool found(const Dart& d) const {
      return MCQ::found(Dart(*d.M(),d.t()));
    };
    bool foundQ(const Dart& d) const {
      return MCQ::foundQ(Dart(*d.M(),d.t()));
    };
    const double quality(const Dart& d) const {
      return MCQ::quality(Dart(*d.M(),d.t()));
    };
    Dart quality() const {
      return MCQ::quality();
    };

    virtual double calcQ(const Dart& d) const;
  };


  class MCQskinny : public MCQtri {
  private:
  public:
    MCQskinny() : MCQtri(true,1.42) {};
    virtual double calcQtri(const Dart& d) const;
  };

  class MCQbig : public MCQtri {
  private:
  public:
    MCQbig() : MCQtri(true,1.0) {};
    virtual double calcQtri(const Dart& d) const;
  };

  class MCQsegm : public MCQ {
  private:
    double encroached_limit_;
    /*!< Larger values are included in the quality set */
  public:
    MCQsegm() : MCQ(false), encroached_limit_(10*MESH_EPSILON) {};
    double calcQ(const Dart& d) const;
    bool segm(const Dart& d) const; /*! true if d or d.orbit1() is found */
  };


  /*!
    \brief Class for constructing Delaunay triangulations
  */
  class MeshC {
  public:
    enum State {State_noT=0, /*!< No triangulation present */
		State_CHT, /*!< Convex hull triangulation */
		State_DT, /*!< Delaunay triangulation */
		State_CDT, /*!< Constrained DT,
			     segment data structures active. */
		State_RCDT /*!< Refined CDT, triangle quality
                                    data structures active. */
    }; /*!< The current triangulation and data structure state. */
  private:
    Mesh *M_;
    /* CDT Constraint and segment data structures: */
    constrListT constr_boundary_; /*! Boundary edge
				    constraints not yet
				    added as segments. */
    constrListT constr_interior_; /*! Interior edge
				    constraints not yet
				    added as segments. */
    MCQsegm boundary_; /*!< Boundary segment */
    MCQsegm interior_; /*!< Interior segment */
    /* RCDT triangle quality data structures: */
    MCQskinny skinny_; /*!< Skinny triangles */
    MCQbig big_;
    /* State variables: */
    State state_;
    bool is_pruned_;

    bool recSwapDelaunay(const Dart& d0);
    Dart splitTriangleDelaunay(const Dart& td, int v);
    Dart splitEdgeDelaunay(const Dart& ed, int v);
    bool insertNode(int v, const Dart& ed);

    bool isSegment(const Dart& d) const;

    /*!
      \brief Make a DT from a CHT, calling LOP.
    */
    bool prepareDT();
    /*!
      \brief Initialise the CDT data structures, and add the
      boundaries as constraint segments.
    */
    bool prepareCDT();
    /*!
      \brief Build a CDT from constraint edge lists. Called by prepareCDT,
      CDTBoundary and CDTInterior.
    */
    bool buildCDT();
    /*!
      \brief Initialise the RCDT data structures.
    */
    bool prepareRCDT(double skinny_limit, double big_limit);
    /*!
      \brief Build a RCDT.
    */
    bool buildRCDT();

  public:
    MeshC() : M_(NULL), state_(State_noT) {};
    MeshC(Mesh* M, bool with_conv_hull)
      : M_(M), state_(State_noT), is_pruned_(false) {
      if (with_conv_hull)
	state_ = State_CHT;
    };

    /*!
      \brief Append vertices

      Return index of the first of the added points.
    */
    int addVertices(const double (*S)[3], int nV)
    {
      M_->S_append(S,nV);
      return M_->nV()-nV;
    };
    Dart swapEdge(const Dart& d);
    Dart splitEdge(const Dart& d, int v);
    Dart splitTriangle(const Dart& d, int v);

    /*!
      \brief Local Optimisation Procedure (LOP)

      Perform LOP to make the input triangulation Delaunay.
     */
    bool LOP(const triangleListT& t_set);
    /*!
      \brief Build Delaunay triangulation (DT)

      The vertices must be covered by the input triangulation, which
      must be convex.  LOP is called to make sure it is Delaunay, if
      not already known by the MeshC to be Delaunay.
     */
    bool DT(const vertexListT& v_set);
    /*!
      \brief Build boundary edge constrained Delaunay triangulation (CDT)
      
      The boundary edge constraints define what regions should be
      removed by a later call to PruneExterior.
    */
    bool CDTBoundary(const constrListT& constr);
    /*!
      \brief Build interior edge constrained Delaunay triangulation (CDT)
    */
    bool CDTInterior(const constrListT& constr);
    /*!
      \brief Alias to CDTInterior
    */
    bool CDT(const constrListT& constr) {
      return CDTInterior(constr);
    };
    /*!
      \brief Remove exterior triangles from a CDT

      Exterior points at the end of the vertex list are removed.

      TODO: allow optional vertex reordering.
    */
    bool PruneExterior();
    /*!
      \brief Refine a CDT
    */
    bool RCDT(double skinny_limit, double big_limit);
  };


} /* namespace fmesh */

#endif
