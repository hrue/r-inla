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
#include <list>

#define Mesh_V_capacity_step_size 128
#define MESH_EPSILON 1e-18

namespace fmesh {

  typedef double Point[3];

  class Dart;
  class M3intO;
  class M3doubleO;
  
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
    Mesh(void) : type_(Mtype_manifold), Vcap_(0), Tcap_(0),
      nV_(0), nT_(0), use_TTi_(false),
      TV_(NULL), TT_(NULL), TTi_(NULL), S_(NULL) {};
    Mesh(Mtype manifold_type, size_t Vcapacity, bool use_TTi=false);
    Mesh(const Mesh& M) : type_(Mtype_manifold), Vcap_(0), Tcap_(0),
      nV_(0), nT_(0), use_TTi_(false),
      TV_(NULL), TT_(NULL), TTi_(NULL), S_(NULL) {
      *this = M;
    };
    Mesh& operator=(const Mesh& M) {
      clear();
      type_ = M.type_;
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

    Mtype type() const { return type_; };
    size_t nV() const { return nV_; };
    size_t nT() const { return nT_; };
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

    Dart locatePoint(const Dart& d0, const Point s, double& delta_min) const;
    
    Dart swapEdge(const Dart& d);
    Dart splitEdge(const Dart& d, int v);
    Dart splitTriangle(const Dart& d, int v);

    double encroachedQuality(const Dart& d) const;
    double skinnyQuality(int t) const;
    double bigQuality(int t) const;
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

  class MCdv {
  public:
    Dart d_;
    double value_;
    MCdv(Dart d, double value) : d_(d), value_(value) {};
    MCdv(const MCdv &T) : d_(T.d_), value_(T.value_) {};
    bool operator<(const MCdv& tb) const {
      return ((value_ < tb.value_) ||
	      ((value_ == tb.value_) &&
	       (d_ < tb.d_)));
    }
  };


  class DartQualitySet {
  private:
    typedef std::map<Dart,double>::value_type map_key_type;
    std::map<Dart,double> darts_; /*!< Darts, mapped to quality */
    std::set<MCdv> darts_quality_; /*!< Set of "bad quality" segment darts */
    double quality_limit_; /*!< Large quality values are included in the set */
  public:
    DartQualitySet() : quality_limit_(0.0) {};
    DartQualitySet(double quality_limit) : quality_limit_(quality_limit) {};
    void clear();
    void insert(const Dart& d, double quality);
    void remove(const Dart& d, double quality);
    void remove(const Dart& d);
    Dart quality(const Dart& d);
    bool empty() const {return darts_.empty();};
    bool empty_quality() const {return darts_quality_.empty();};
    bool found(const Dart& d) const;
    bool found_quality(const Dart& d) const;
    Dart get_quality() const;
  };

  class SegmSet {
  private:
    const Mesh *M_;
    DartQualitySet segm_;
  public:
    SegmSet() : M_(NULL) {};
    SegmSet(const Mesh& M) : M_(&M), segm_(10*MESH_EPSILON) {};
    void clear();
    void insert(const Dart& d) {
      double quality = M_->encroachedQuality(d);
      segm_.insert(d,quality);
    };
    void remove(const Dart& d) {
      segm_.remove(d);
    };
    bool empty() const {return segm_.empty();};
    bool empty_encroached() const {return segm_.empty_quality();};
    bool found(const Dart& d) const {return segm_.found(d);};
    bool found_encroached(const Dart& d) const {return segm_.found_quality(d);};
    Dart get_encroached() const {return segm_.get_quality();};
  };


  /*!
    \brief Class for constructing Delaunay triangulations
  */
  class MeshConstructor {
  public:
    enum State {State_noT=0, /*!< No triangulation present */
		State_CHT, /*!< Convex hull triangulation */
		State_DT, /*!< Delaunay triangulation */
		State_CDT, /*!< Constrained DT,
			     segment data structures active. */
		State_RCDT /*!< Refined CDT, triangle quality
                                    data structures active. */
    }; /*!< The current triangulation and data structure state. */
    typedef std::list<int> vertex_input_type;
    typedef std::list<int> triangle_input_type;
    typedef std::pair<int,int> constraint_type;
    typedef std::list<constraint_type> constraint_input_type;
    typedef std::list<constraint_type> constraint_list_type;
  private:
    Mesh *M_;
    /* CDT Constraint and segment data structures: */
    constraint_list_type constr_boundary_; /*! Boundary edge
                                              constraints not yet
                                              added as segments. */
    constraint_list_type constr_interior_; /*! Interior edge
                                              constraints not yet
                                              added as segments. */
    SegmSet boundary_; /*!< Boundary segment */
    SegmSet interior_; /*!< Interior segment */
    /* RCDT triangle quality data structures: */
    DartQualitySet skinny_; /*!< Skinny triangles */
    DartQualitySet big_;
    double skinny_limit_;
    double big_limit_;
    /* State variables: */
    State state_;
    bool is_pruned_;

    bool recSwapDelaunay(const Dart& d0);
    Dart splitTriangleDelaunay(const Dart& td, int v);
    Dart splitEdgeDelaunay(const Dart& ed, int v);
    bool insertNode(int v, const Dart& ed);

    bool isSegmentDart(const Dart& d) const;

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
    MeshConstructor() : M_(NULL), state_(State_noT) {};
    MeshConstructor(Mesh* M, bool with_conv_hull)
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
    bool LOP(const triangle_input_type& t_set);
    /*!
      \brief Build Delaunay triangulation (DT)

      The vertices must be covered by the input triangulation, which
      must be convex.  LOP is called to make sure it is Delaunay, if
      not already known by the MeshConstructor to be Delaunay.
     */
    bool DT(const vertex_input_type& v_set);
    /*!
      \brief Build boundary edge constrained Delaunay triangulation (CDT)
      
      The boundary edge constraints define what regions should be
      removed by a later call to PruneExterior.
    */
    bool CDTBoundary(const constraint_input_type& constr);
    /*!
      \brief Build interior edge constrained Delaunay triangulation (CDT)
    */
    bool CDTInterior(const constraint_input_type& constr);
    /*!
      \brief Alias to CDTInterior
    */
    bool CDT(const constraint_input_type& constr) {
      CDTInterior(constr);
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



  struct Traits {
    /*!
      Compute dart half-space test for a point.
      positive if s is to the left of the edge defined by d.
     */
    static double inLeftHalfspace(const Dart& d, const double s[3]);
    static double inCircumcircle(const Dart& d, const double s[3]);
    static bool circumcircleTest(const Dart& d);
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
