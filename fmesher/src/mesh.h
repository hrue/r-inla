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
#define MESH_EPSILON 1e-10

#ifndef NOT_IMPLEMENTED
#define NOT_IMPLEMENTED (std::cout << "Not implemented: "	\
			 << __FILE__ << "(" << __LINE__ << ") "	\
			 << __PRETTY_FUNCTION__ << std::endl);
#endif

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


  struct Vec {  
    static void copy(double* s, const Point s0)
    {
      s[0] = s0[0];
      s[1] = s0[1];
      s[2] = s0[2];
    };
    static void rescale(double* s, double s1)
    {
      s[0] *= s1;
      s[1] *= s1;
      s[2] *= s1;
    };
    static void scale(double* s, const Point s0, double s1)
    {
      s[0] = s0[0]*s1;
      s[1] = s0[1]*s1;
      s[2] = s0[2]*s1;
    };
    static void diff(double* s,const Point s0, const Point s1)
    {
      s[0] = s0[0]-s1[0];
      s[1] = s0[1]-s1[1];
      s[2] = s0[2]-s1[2];
    };
    static void sum(double* s,const Point s0, const Point s1)
    {
      s[0] = s0[0]+s1[0];
      s[1] = s0[1]+s1[1];
      s[2] = s0[2]+s1[2];
    };
    static void accum(double* s, const Point s0, double s1 = 1.0)
    {
      s[0] += s0[0]*s1;
      s[1] += s0[1]*s1;
      s[2] += s0[2]*s1;
    };
    static double scalar(const Point s0, const Point s1)
    {
      return (s0[0]*s1[0]+s0[1]*s1[1]+s0[2]*s1[2]);
    };
    static double length(const Point s0)
    {
      return (std::sqrt(s0[0]*s0[0]+s0[1]*s0[1]+s0[2]*s0[2]));
    };
    static void cross(double* s, const Point s0, const Point s1)
    {
      s[0] = s0[1]*s1[2]-s0[2]*s1[1];
      s[1] = s0[2]*s1[0]-s0[0]*s1[2];
      s[2] = s0[0]*s1[1]-s0[1]*s1[0];
    };
    static double cross2(const Point s0, const Point s1)
    {
      return (s0[0]*s1[1]-s0[1]*s1[0]);
    };
  };


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
    Mesh& useX11(bool use_X11, bool draw_text,
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
    void triangleCircumcenter(int t, double* c) const;
    double triangleCircumcircleRadius(int t) const;
    double triangleShortestEdge(int t) const;
    double triangleLongestEdge(int t) const;
    double edgeEncroached(const Dart& d, const double s[3]) const;
    
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
    int v() const { return M_->TV_[t_][vi_]; };

    bool isnull() const { return (!M_); };
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


} /* namespace fmesh */

#endif
