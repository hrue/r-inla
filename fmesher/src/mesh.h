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

namespace fmesh {

  typedef double Point[3];

  class Xtmpl;
  class Dart;
  class Mesh;
  class MeshConstructor;
  class MintO;
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
    bool use_VT_;
    bool use_TTi_;
    bool use_X11_;
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

    
  public:
    Mesh(void) : type_(Mtype_manifold), Vcap_(0), Tcap_(0),
      nV_(0), nT_(0), use_VT_(false), use_TTi_(true), use_X11_(false),
		 TV_(NULL), TT_(NULL), TTi_(NULL), S_(NULL), X11_(NULL) {};
    Mesh(Mtype manifold_type, size_t Vcapacity, bool use_VT=true, bool use_TTi=false);
    Mesh(const Mesh& M) : type_(Mtype_manifold), Vcap_(0), Tcap_(0),
      nV_(0), nT_(0), use_VT_(true), use_TTi_(false), use_X11_(false),
      TV_(NULL), TT_(NULL), TTi_(NULL), S_(NULL), X11_(NULL) {
      *this = M;
    };
    Mesh& operator=(const Mesh& M) {
      clear();
      type_ = M.type_;
      useVT(M.use_VT_);
      useTTi(M.use_TTi_);
      useX11(M.use_X11_);
      S_set(M.S_,M.nV_);
      TV_set(M.TV_,M.nT_);
      return *this;
    };
    ~Mesh();
    Mesh& clear();

    /*!
      \brief Check the storage capacity, and increase if necessary
    */
    Mesh& check_capacity(size_t nVc, size_t nTc);

    bool useVT() const { return use_VT_; }
    Mesh& useVT(bool use_VT);
    bool useTTi() const { return use_TTi_; }
    Mesh& useTTi(bool use_TTi);

    bool useX11() const { return use_X11_; }
    Mesh& useX11(bool use_X11);
    void redrawX11();

    Mtype type() const { return type_; };
    size_t nV() const { return nV_; };
    size_t nT() const { return nT_; };
    const int (*TV() const)[3] { return TV_; };
    const int (*TT() const)[3] { return TT_; };
    const int (*VT() const) { return VT_; };
    const int (*TTi() const)[3] { return TTi_; };
    const double (*S() const)[3] { return S_; };
    Xtmpl *X11() { return X11_; };
    M3intO TVO() const;
    M3intO TTO() const;
    MintO VTO() const;
    M3intO TTiO() const;
    M3doubleO SO() const;
    
    Mesh& S_set(const double (*S)[3], int nV);
    Mesh& TV_set(const int (*TV)[3], int nT); 
    Mesh& S_append(const double (*S)[3], int nV);
    Mesh& TV_append(const int (*TV)[3], int nT); 

    Dart locatePoint(const Dart& d0, const Point s, double* delta_min) const;
    Dart locateVertex(const Dart& d0, const int v) const;
    
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
    typedef std::set<MCdv> set_type;
    typedef std::map<Dart,double> map_type;
    typedef map_type::value_type map_key_type;
    map_type darts_; /*!< Darts, mapped to quality */
    set_type darts_quality_; /*!< Set of "bad quality" segment darts */
    double quality_limit_; /*!< Large quality values are included in the set */
  public:
    DartQualitySet() : quality_limit_(0.0) {};
    DartQualitySet(double quality_limit) : quality_limit_(quality_limit) {};
    void clear() {
      darts_.clear();
      darts_quality_.clear();
    };
    void insert(const Dart& d, double quality);
    void erase(const Dart& d);
    bool empty() const {return darts_.empty();};
    bool empty_quality() const {return darts_quality_.empty();};
    bool found(const Dart& d) const;
    bool found_quality(const Dart& d) const;
    const double quality(const Dart& d) const;
    Dart quality_dart() const;
  };

  class SegmSet {
  private:
    const Mesh *M_;
    DartQualitySet segm_;
  public:
    SegmSet() : M_(NULL) {};
    SegmSet(const Mesh& M) : M_(&M), segm_(10*MESH_EPSILON) {};
    void clear() { segm_.clear(); };
    void insert(const Dart& d) {
      double quality = M_->encroachedQuality(d);
      segm_.insert(d,quality);
    };
    void erase(const Dart& d) { segm_.erase(d); };
    bool empty() const {return segm_.empty();};
    bool empty_encroached() const {return segm_.empty_quality();};
    bool found(const Dart& d) const {return segm_.found(d);};
    bool found_encroached(const Dart& d) const {return segm_.found_quality(d);};
    Dart get_encroached() const {return segm_.quality_dart();};
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



  struct Traits {
    /*!
      Compute dart half-space test for a point.
      positive if s is to the left of the edge defined by d.
     */
    static double inLeftHalfspace(const Dart& d, const double s[3]);
    static double inCircumcircle(const Dart& d, const double s[3]);
    static bool circumcircleOK(const Dart& d);
  };




  class MintO {
    friend std::ostream& operator<<(std::ostream& output, const MintO& MO);
  private:
    size_t n_;
    const int (*M_);
  public:
    MintO(const int (*M),size_t n) : n_(n), M_(M) {};
  };

  class M3intO {
    friend std::ostream& operator<<(std::ostream& output, const M3intO& MO);
  private:
    size_t n_;
    const int (*M_)[3];
  public:
    M3intO(const int (*M)[3],size_t n) : n_(n), M_(M) {};
  };

  class M3doubleO {
    friend std::ostream& operator<<(std::ostream& output, const M3doubleO& MO);
  private:
    size_t n_;
    const double (*M_)[3];
  public:
   M3doubleO(const double (*M)[3],size_t n) : n_(n), M_(M) {};
  };









  class Xtmpl {
  private:
    int window_;
    char* name_char_;
    int sx_, sy_;
    double minx_, maxx_, miny_, maxy_;
  public:
    Xtmpl() : window_(-1), name_char_(NULL) {};
    void reopen(int sx, int sy) {
      if (!(window_<0))
	close();
      window_ = 0;
      sx_ = sx;
      sy_ = sy;
      xtmpl_window = window_;
      xtmpl_open(sx_,sy_,name_char_);
    };
    void open(std::string name,
	      int sx, int sy) {
      if (!(window_<0))
	close();
      window_ = 0;
      sx_ = sx;
      sy_ = sy;
      if (name_char_) delete[] name_char_;
      name_char_ = new char[name.length()+1];
      name.copy(name_char_,name.length(),0);
      name_char_[name.length()] = '\0';
      xtmpl_window = window_;
      xtmpl_open(sx,sy,name_char_);
      setAxis(-0.05,1.05,-0.05,1.05);
    };
    void close() {
      if (window_<0)
	return;
      xtmpl_window = window_;
      xtmpl_close();
    };
    ~Xtmpl() {
      close();
      if (name_char_) delete[] name_char_;
    };

    void clear() {
      xtmpl_window = window_;
      xtmpl_clear();
    }

    void setSize(int sx, int sy) {
      reopen(sx,sy);
    };
    void setAxis(double minx, double maxx,
		 double miny, double maxy) {
      clear();
      minx_ = minx;
      maxx_ = maxx;
      miny_ = miny;
      maxy_ = maxy;
    };

    void lineFG(const double* s0, const double* s1)
    {
      xtmpl_window = window_;
      xtmpl_draw_line((int)(sx_*(s0[0]-minx_)/(maxx_-minx_)),
		      (int)(sy_*(s0[1]-miny_)/(maxy_-miny_)),
		      (int)(sx_*(s1[0]-minx_)/(maxx_-minx_)),
		      (int)(sy_*(s1[1]-miny_)/(maxy_-miny_)));
    };
    void text(const double* s0, std::string str)
    {
      char* str_ = new char[str.length()+1];
      str.copy(str_,str.length(),0);
      str_[str.length()] = '\0';
      xtmpl_window = window_;
      xtmpl_text((int)(sx_*(s0[0]-minx_)/(maxx_-minx_)),
		 (int)(sy_*(s0[1]-miny_)/(maxy_-miny_)),
		 str_,str.length());
      delete[] str_;
    };

  };






} /* namespace fmesh */

#endif
