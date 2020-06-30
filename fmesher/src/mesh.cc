#include <cstddef>
#include <cstring>
#include <set>
#include <map>
#include <sstream>
#include <cmath>
#include <ctime>
#include <cerrno>

#include "predicates.hh"
#include "x11utils.hh"
#include "meshc.hh"
#include "mesh.hh"

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


using std::cout;
using std::endl;

namespace fmesh {



  


  /*
    T-E+V=2
    
    closed 2-manifold triangulation:
    E = T*3/2
    T = 2*V-4

    simply connected 2-manifold triangulation:
    T <= 2*V-5

  */

  Mesh::Mesh(Mtype manifold_type,
	     size_t V_capacity,
	     bool use_VT,
	     bool use_TTi) : type_(manifold_type),
			     use_VT_(use_VT), use_TTi_(use_TTi),
			     TV_(), TT_(), VT_(), TTi_(), S_()
#ifndef FMESHER_NO_X
			   , X11_(NULL), X11_v_big_limit_(0)
#endif
  {
    if (V_capacity > 0) {
      TV_.capacity(V_capacity*2);
      TT_.capacity(V_capacity*2);
      if (use_VT_)
	VT_.capacity(V_capacity);
      if (use_TTi_)
	TTi_.capacity(V_capacity*2);
      S_.capacity(V_capacity);
    }
  }

  Mesh::~Mesh()
  {
    clear();
  }

  Mesh& Mesh::clear()
  {
    empty();
#ifndef FMESHER_NO_X
    if (X11_) { delete X11_; X11_ = NULL; }
#endif
    return *this;
  }

  Mesh& Mesh::empty()
  {
    TV_.clear();
    TT_.clear();
    VT_.clear();
    TTi_.clear();
    S_.clear();
    return *this;
  }

  Mesh& Mesh::operator=(const Mesh& M)
  {
    clear();
    type_ = M.type_;
    useVT(M.use_VT_);
    useTTi(M.use_TTi_);
#ifndef FMESHER_NO_X
    if (M.X11_) {
      X11_ = new Xtmpl(*M.X11_);
    } else {
      X11_ = NULL;
    }
    X11_v_big_limit_ = M.X11_v_big_limit_;
#endif
    S_set(M.S_);
    TV_set(M.TV_);
    return *this;
  }

  Mesh& Mesh::check_capacity(size_t nVc, size_t nTc)
  {
    if (nVc > S_.capacity()) {
      if (use_VT_) VT_.capacity(nVc);
      S_.capacity(nVc);
    }

    if (nTc > TV_.capacity()) {
      TV_.capacity(nTc);
      TT_.capacity(nTc);
      if (use_TTi_) TTi_.capacity(nTc);
    }

    return *this;
  }



  Mesh& Mesh::rebuildTT()
  {
    typedef IntPair E_Type;
    typedef std::map<E_Type,int> ET_Type;
    int t, vi;
    E_Type E0,E1;
    ET_Type::const_iterator Ei;
    ET_Type ET;
    /* Pass 1: */
    for (t=0; t<(int)nT(); t++) {
      const Int3& TVt = TV_[t];
      for (vi=0; vi<3; vi++) {
	E0 = IntPair(TVt[(vi+1)%3],TVt[(vi+2)%3]);
	E1 = IntPair(TVt[(vi+2)%3],TVt[(vi+1)%3]);
	Ei = ET.find(E1);
	if (Ei != ET.end()) { /* Found neighbour */
	  TT_(t)[vi] = Ei->second;
	} else { /* Either on boundary, or not added yet. */
	  TT_(t)[vi] = -1;
	}
	ET.insert(ET_Type::value_type(E0,t));
      }
    }

    /* Pass 2: */
    for (t=0; t<(int)nT(); t++) {
      const Int3& TVt = TV_[t];
      for (vi=0; vi<3; vi++) {
	if (TT_[t][vi]>=0) continue;
	E1 = IntPair(TVt[(vi+2)%3],TVt[(vi+1)%3]);
	Ei = ET.find(E1);
	if (Ei != ET.end()) { /* Found neighbour */
	  TT_(t)[vi] = Ei->second;
	}
      }
    }

    return *this;
  }



  Mesh& Mesh::update_VT(const int v, const int t)
  {
    if ((use_VT_) && (v<(int)nV()) & (t<(int)nT()) && (VT_[v]<0))
      VT_(v) = t;
    return *this;
  }

  Mesh& Mesh::set_VT(const int v, const int t)
  {
    if ((use_VT_) && (v<(int)nV()) & (t<(int)nT()))
      VT_(v) = t;
    return *this;
  }

  Mesh& Mesh::update_VT_triangle(const int t)
  {
    if ((use_VT_) && (t<(int)nT()) && (t>=0)) {
      const Int3& TVt = TV_[t];
      for (int vi=0; vi<3; vi++) {
	int v = TVt[vi];
	if (VT_[v]<0)
	  VT_(v) = t;
      }
    }
    return *this;
  }

  Mesh& Mesh::set_VT_triangle(const int t)
  {
    if ((use_VT_) && (t<(int)nT()) && (t>=0)) {
      const Int3& TVt = TV_[t];
      for (int vi=0; vi<3; vi++)
	VT_(TVt[vi]) = t;
    }
    return *this;
  }

  Mesh& Mesh::update_VT_triangles(const int t_start)
  {
    if (use_VT_) {
      for (int t=t_start; t<(int)nT(); t++)
	update_VT_triangle(t);
    }
    return *this;
  }

  Mesh& Mesh::reset_VT(const int v_start)
  {
    if (use_VT_) {
      for (int v=v_start; v<(int)nV(); v++)
	set_VT(v,-1);
    }
    return *this;
  }

  Mesh& Mesh::rebuild_VT()
  {
    if ((!use_VT_) || (!S_.capacity())) {
      VT_.clear();
    } else {
      VT_.rows(0);
      VT_.capacity(S_.capacity());
      reset_VT(0);
      update_VT_triangles(0);
    }
    return *this;
  }

  Mesh& Mesh::rebuildTTi()
  {
    int t, vi, v, t2, vi2;
    if (!use_TTi_) {
      TTi_.clear();
      return *this;
    }
    TTi_.rows(nT());
    if (!TV_.capacity())
      return *this;
    TTi_.capacity(TV_.capacity());
    for (t=0; t<(int)nT(); t++) {
      for (vi=0; vi<3; vi++) {
	v = TV_[t][vi];
	t2 = TT_[t][(vi+2)%3];
	if (t2>=0) {
	  for (vi2 = 0; (vi2<3) && (TV_[t2][vi2] != v); vi2++) { }
	  if (vi2<3) {
	    TTi_(t)[(vi+2)%3] = (vi2+1)%3;
	  } else {
	    /* Error! This should never happen! */
	    LOG("ERROR\n");
	  }
	} else {
	  TTi_(t)[(vi+2)%3] = -1;
	}
      }
    }
    return *this;
  }


  Mesh& Mesh::useVT(bool use_VT)
  {
    if (use_VT_ != use_VT) {
      use_VT_ = use_VT;
      rebuild_VT();
    }
    return *this;
  }

  Mesh& Mesh::useTTi(bool use_TTi)
  {
    if (use_TTi_ != use_TTi) {
      use_TTi_ = use_TTi;
      rebuildTTi();
    }
    return *this;
  }


  SparseMatrix<int> Mesh::VV() const
  {
    SparseMatrix<int> VV;
    for (int t=0; t<(int)nT(); t++) {
      VV(TV_[t][0],TV_[t][1]) = 1;
      VV(TV_[t][0],TV_[t][2]) = 1;
      VV(TV_[t][1],TV_[t][0]) = 1;
      VV(TV_[t][1],TV_[t][2]) = 1;
      VV(TV_[t][2],TV_[t][0]) = 1;
      VV(TV_[t][2],TV_[t][1]) = 1;
    };
    return VV;
  }
    

#ifndef FMESHER_NO_X
  void Mesh::setX11delay(double set_delay)
  {
    if (X11_) {
      X11_->delay(set_delay);
    }
  }

  Mesh& Mesh::useX11(bool use_X11,
		     bool draw_text,
		     int sx, int sy,
		     double minx,
		     double maxx,
		     double miny,
		     double maxy,
		     std::string name)
  {
    if (use_X11) {
      if (!X11_) { /* Init. */
	if (type_ == Mtype_sphere)
	  X11_ = new Xtmpl(draw_text,sx*2,sy,minx,2*maxx-minx,miny,maxy,name);
	else
	  X11_ = new Xtmpl(draw_text,sx,sy,minx,maxx,miny,maxy,name);
	redrawX11("");
      } else {
	if (type_ == Mtype_sphere) {
	  X11_->reopen(sx*2,sy,draw_text);
	  X11_->setAxis(minx,2*maxx-minx,miny,maxy);
	} else {
	  X11_->reopen(sx,sy,draw_text);
	  X11_->setAxis(minx,maxx,miny,maxy);
	}
	redrawX11("");
      }
    } else { /* Destroy. */
      if (X11_) { /* Destroy. */
	delete X11_;
	X11_ = NULL;
      }
    }
    return *this;
  }
#endif



  Mesh& Mesh::S_set(const Matrix3double& S)
  {
    S_.rows(0); /* Avoid possible unnecessary copy. */
    S_append(S);
    return *this;
  }
  
  Mesh& Mesh::TV_set(const Matrix3int& TV)
  {
    TV_.rows(0); /* Avoid possible unnecessary copy. */
    TV_append(TV);
    return *this;
  }

  Mesh& Mesh::S_append(const Point& s)
  {
    S_(nV()) = s;
    if (use_VT_)
      reset_VT(nV()-1);
    return *this;
  }

  Mesh& Mesh::S_append(const Matrix3double& S)
  {
    S_.append(S);
    if (use_VT_)
      reset_VT(nV()-S.rows());
    return *this;
  }

#ifndef FMESHER_NO_X
  void Mesh::drawX11point(int v, bool fg)
  {
    if (!X11_) return;

    int szbig = (nV()>50 ? (nV()>500 ? 1 : 3) : 5);
    int szsmall = 1;

    const Point& s = S_[v];
    bool otherside = (s[2]<0);

    if (type_==Mtype_sphere) {
      double offset(otherside ? X11_->width()/2.0 : 0.0);
      int sz = ((v<X11_v_big_limit_) ? szbig : szsmall);
      X11_->dot_on_sphere(fg,s,sz,offset);
    } else {
      int sz = ((v<X11_v_big_limit_) ? szbig : szsmall);
      X11_->dot(fg,s,sz);
    }
  }
#endif

#ifndef FMESHER_NO_X
  void Mesh::drawX11triangle(int t, bool fg)
  {
    if (!X11_) return;

    int szbig = (nV()>50 ? (nV()>500 ? 1 : 3) : 5);
    int szsmall = 1;

    const Int3& v = TV_[t];
    Point s[3];
    Point s0;
    bool arc_upper[3] = {false, false, false};
    bool arc_lower[3] = {false, false, false};
    Point arc_split_s_upper[3];
    Point arc_split_s[3];
    Point arc_split_s_lower[3];

    s0[0] = 0.0;
    s0[1] = 0.0;
    s0[2] = 0.0;
    for (int vi=0;vi<3;vi++) {
      Vec::copy(s[vi],S_[v[vi]]);
      Vec::accum(s0,s[vi],1.0/3.0);
    }
    if (type_==Mtype_sphere) {
      double l = std::sqrt(Vec::scalar(s0,s0));
      Vec::rescale(s0,l);

      if (s[0][2]>=0.) { arc_upper[1] = true; arc_upper[2] = true; }
      if (s[0][2]< 0.) { arc_lower[1] = true; arc_lower[2] = true; }
      if (s[1][2]>=0.) { arc_upper[0] = true; arc_upper[2] = true; }
      if (s[1][2]< 0.) { arc_lower[0] = true; arc_lower[2] = true; }
      if (s[2][2]>=0.) { arc_upper[0] = true; arc_upper[1] = true; }
      if (s[2][2]< 0.) { arc_lower[0] = true; arc_lower[1] = true; }

      for (int i=0; i<3; i++) {
	int i1 = (i+1)%3;
	int i2 = (i+2)%3;
	if (arc_upper[i] && arc_lower[i]) {
	  arc_split_s_lower[i] = ((s[i1][2] <  0.0) ? s[i1] : s[i2]);
	  arc_split_s_upper[i] = ((s[i1][2] >= 0.0) ? s[i1] : s[i2]);
	  /*  s[i1]z * (1-?) + s[i2]z * ? = 0
	      ? = - s[i1]z / (s[i2]z - s[i1]z)
	   */
	  Vec::copy(arc_split_s[i],s[i1]);
	  double beta = - s[i1][2] / (s[i2][2] - s[i1][2]);
	  Vec::rescale(arc_split_s[i],1-beta);
	  Vec::accum(arc_split_s[i],s[i2],beta);
	  Vec::rescale(arc_split_s[i],1.0/arc_split_s[i].length());
	} else if (arc_upper[i]) {
	  arc_split_s_upper[i] = s[i1];
	  arc_split_s[i] = s[i2];
	} else {
	  arc_split_s_lower[i] = s[i1];
	  arc_split_s[i] = s[i2];
	}
      }

      /*
      Point r0;
      Point r1;
      Point n;
      Vec::diff(r0,s[1],s[0]);
      Vec::diff(r1,s[2],s[0]);
      Vec::cross(n,r0,r1);
      otherside = (n[2]<0);

      Vec::copy(r0,s[0]);
      Vec::accum(r0,s[1]);
      Vec::accum(r0,s[2]);
      Vec::rescale(r0,1/3.);

      LOG_(arc_split << "\t" << r0[2] << "\t" << n[2] << endl);
      */
    }
    /* Draw triangle slightly closer to center. */
    if (type_==Mtype_sphere) {
      // for (int vi=0;vi<3;vi++) {
      // 	Vec::diff(s[vi],s[vi],s0);
      // 	Vec::rescale(s[vi],0.975);
      // 	Vec::accum(s[vi],s0);
      // 	Vec::rescale(s[vi],1./Vec::length(s[vi]));
      // }
      double offset = X11_->width()/2.;
      int sz = ((v[0]<X11_v_big_limit_) ? szbig : szsmall);
      X11_->dot_on_sphere(fg,s[0],sz,((s[0][2]<0.0) ? offset : 0.0));
      sz = ((v[1]<X11_v_big_limit_) ? szbig : szsmall);
      X11_->dot_on_sphere(fg,s[1],sz,((s[1][2]<0.0) ? offset : 0.0));
      sz = ((v[2]<X11_v_big_limit_) ? szbig : szsmall);
      X11_->dot_on_sphere(fg,s[2],sz,((s[2][2]<0.0) ? offset : 0.0));
      for (int i=0; i<3; i++) {
	if (arc_upper[i])
	  X11_->arc(fg,arc_split_s[i],arc_split_s_upper[i],0.0);
	if (arc_lower[i])
	  X11_->arc(fg,arc_split_s[i],arc_split_s_lower[i],offset);
      }
    } else {
      //      for (int vi=0;vi<3;vi++) {
      //       	Vec::diff(s[vi],s[vi],s0);
      //       	Vec::rescale(s[vi],0.975);
      //       	Vec::accum(s[vi],s0);
      //      }
      int sz = ((v[0]<X11_v_big_limit_) ? szbig : szsmall);
      X11_->dot(fg,s[0],sz);
      sz = ((v[1]<X11_v_big_limit_) ? szbig : szsmall);
      X11_->dot(fg,s[1],sz);
      sz = ((v[2]<X11_v_big_limit_) ? szbig : szsmall);
      X11_->dot(fg,s[2],sz);
      X11_->line(fg,s[0],s[1]);
      X11_->line(fg,s[1],s[2]);
      X11_->line(fg,s[2],s[0]);
    }
    /* Draw vertex indices even closer to center. */
    for (int vi=0;vi<3;vi++) {
      Vec::diff(s[vi],s[vi],s0);
      Vec::rescale(s[vi],0.9);
      Vec::accum(s[vi],s0);
    }
    for (int vi=0;vi<3;vi++) {
      std::ostringstream ss;
      //      ss << "(" << TV_[t][vi] << "," << TT_[t][vi] << ")";
      ss << "" << TV_[t][vi] << "";
      X11_->text(fg,s[vi],ss.str());
    }
    /* Draw triangle indices at center. */
    {
      std::ostringstream ss;
      ss << "t" << t << "";
      X11_->text(fg,s0,ss.str());
    }
  }
#endif

#ifndef FMESHER_NO_X
  void Mesh::redrawX11(std::string str)
  {
    if (!X11_) return;
    if (verbose_ > 0)
      std::cout << str << std::endl;
    
    X11_->clear();
    for (int v=0;v<(int)nV();v++)
      drawX11point(v,true);
    for (int t=0;t<(int)nT();t++)
      drawX11triangle(t,true);

    X11_->delay();
  }
#endif /* FMESHER_NO_X */
  
  Mesh& Mesh::TV_append(const Matrix3int& TV)
  {
    TV_.append(TV);
    if (use_VT_)
      update_VT_triangles(nT()-TV.rows());
    rebuildTT();
    rebuildTTi();
#ifndef FMESHER_NO_X
    redrawX11(std::string("TV appended"));
#endif
    return *this;
  }





  void Mesh::triangleBoundingBox(const Point& s0,
				 const Point& s1,
				 const Point& s2,
				 Point& mini, Point& maxi) const
  {
    for (int d=0; d<3; d++) {
      mini[d] = (s0[d] < s1[d] ?
		 (s0[d] < s2[d] ? s0[d] : s2[d] ) :
		 (s2[d] < s1[d] ? s2[d] : s1[d] ) );
      maxi[d] = (s0[d] > s1[d] ?
		 (s0[d] > s2[d] ? s0[d] : s2[d] ) :
		 (s2[d] > s1[d] ? s2[d] : s1[d] ) );
    }
  }



  /*!
   \brief Calculate the length of an edge. 

   For planes and triangular manifolds, the edge length is
   \f$L=\|s_1-s_0\|\f$.

   On the sphere, the "obvious" arccos-formula for the length of a
   geodesic may be numerically unstable for short and long edges.  Use
   the arctan-formula instead, that should handle all cases
   \f$L\in[0,\pi]\f$.
   \f{align*}{
   L &= \operatorname{acos}(s_0 \cdot s_1) \\
   \sin(L/2) &= \|s_1-s_0\|/2 \\
   \cos(L/2) &= \|s_0+s_1\|/2 \\
   L &= 2 \cdot \operatorname{atan2}(\|s_1-s_0\|,\|s_0+s_1\|)
   \f}
   */
  double Mesh::edgeLength(const Point& s0, const Point& s1) const
  {
    Point e;
    Vec::diff(e,s1,s0);
    double len = Vec::length(e);

    if (type_==Mesh::Mtype_sphere) {
      Point ssum;
      Vec::sum(ssum,s1,s0);
      len = 2.0*std::atan2(len,Vec::length(ssum));
    }

    return len;
  }

  /*!
    
    \see Mesh::edgeLength(const Dart& d)
  */
  double Mesh::edgeLength(const Dart& d) const
  {
    int t(d.t());
    if ((t<0) || (t>=(int)nT())) return 0.0;

    return edgeLength(S_[d.v()],S_[d.vo()]);
  }

  /*!
    \brief Calculate a triangle area.

    Notation:\n
    \f$s_k\f$ are the cartesian coordinates of vertex \f$k\f$.\n
    \f$e_0=s_2-s_1\f$ is the CCW edge vector opposite node 0,
    and similarly for the other nodes.\n
    \f$n_0=e_2\times(-e_1)=e_1\times e_2\f$ is an unnormalised
    outwards triangle normal vector.

    For planes, calculate the signed area, from the z-component of
    \f$(n_0+n_1+n_2)/6\f$.\n
    For manifolds, calculate the absolute area, from
    \f$A=\|n_0+n_1+n_3\|/6\f$.\n
    For spheres, calculate the CCW interior geodesic triangle area,
    with range \f$[0,4\pi)\f$.

    The spherical geodesic triangle area is given by the spherical
    excess formula \f$A = \theta_0+\theta_1+\theta_2-\pi\f$, where
    \f$\theta_i\f$ is the interior angle at node i.
   \f{align*}{
   d_1 &= e_1-s_0(s_0\cdot e_1) \\
   d_2 &= e_2-s_0(s_0\cdot e_2) \\
   \|d_1\| \|d_2\| \sin \theta_0 &= s_0\cdot(d_1\times d_2) \\
   \|d_1\| \|d_2\| \cos \theta_0 &= - (d_1 \cdot d_2) \\
   \theta_0 &= \operatorname{atan2}(s_0\cdot(d_1\times d_2),- (d_1 \cdot d_2))
               \mod [0,2\pi)
   \f}
   Rewrite in terms of the original edge vectors:
   \f{align*}{
   s_0\cdot(d_1\times d_2) &= s_0\cdot(e_1\times e_2
                             - (e_1\times s_0)(s_0\cdot e_2)
                             - (e_2\times s_0)(s_0\cdot e_1))
                       = s_0\cdot(e_1\times e_2) \\
   (d_1\cdot d_2) &= (e_1\cdot e_2) - (s_0\cdot e_1)(s_0\cdot e_2) \\
   \theta_0 &= \operatorname{atan2}(s_0\cdot(e_1\times e_2),
                      (s_0\cdot e_1)(s_0\cdot e_2) - (e_1 \cdot e_2))
	       \mod [0,2\pi)
   \f}
   The two remaining angles, \f$\theta_1\f$ and \f$\theta_2\f$, are
   obtained by permuting indices.

   If \f$s_0\cdot(e_1\times e_2)\f$ is negative, the triangle covers
   more than a hemisphere, and is non-convex (all interior angles are
   \f$>\pi\f$), with area \f$>2\pi\f$.  Since the total area of the
   sphere is \f$4\pi\f$, a geodesic triangulation can have at most one
   triangle of this type.

   If the vertices are co-planar with the origin, some special cases
   need to be analysed.  If the flat triangle spanned by the vertices
   contains the origin, the geodesic triangle covers a complete
   hemisphere with area \f$2\pi\f$, which is also calculated by the
   formula.  If the origin lies on the boundary of the flat triangle,
   two of the vertices are antipodes, say \f$s_1=-s_0\f$, and the
   formula breaks down, with both arguments to \p atan2 equal to zero.
   Such triangles are not uniquely defined, since the geodesic between
   \f$s_0\f$ and \f$-s_0\f$ is not unique.  For the other co-planar
   cases, the geodesic triangle is degenerate, and the exact result of
   the formula is \f$\pi+\pi-\pi-\pi=0\f$, reflecting the different
   signs of the second parameter to \p atan2.  In practice, the
   vertices may not be numerically co-planar, and the calculated area
   may become \f$4\pi\f$ instead.


  */
  /*
    Don't use the following, since they rely on the edge lengths being known.

    Heron's formula:
    a,b,c edge lengths
    s = (a+b+c)/2
    Area = sqrt(s(s-a)(s-b)(s-c)) 
    
    Numerically stable version from
    http://www.eecs.berkeley.edu/~wkahan/Triangle.pdf
    a >= b >= c
    Area = sqrt( (a+(b+c)) (c-(a-b)) (c+(a-b)) (a+(b-c)) )/4

    l'Huilier's Theorem for spherical triangle areas:
    a,b,c edge lengths
    s = (a+b+c)/2
    tan(E / 4) = sqrt(tan(s / 2) 
                      tan((s - a) / 2)
		      tan((s - b) / 2)
		      tan((s - c) / 2))
    Area = E  (E = spherical excess)
  */
  double Mesh::triangleArea(const Point& s0, const Point& s1, const Point& s2) const
  {
    Point e0, e1, e2;
    Vec::diff(e0,s2,s1);
    Vec::diff(e1,s0,s2);
    Vec::diff(e2,s1,s0);

    double area;

    switch (type_) {
    case Mesh::Mtype_manifold:
      {
	/* Calculate the upwards unscaled normal(s), calculate length. */
	Point n0, n1, n2;
	Vec::cross(n0,e1,e2);
	Vec::cross(n1,e2,e0);
	Vec::cross(n2,e0,e1);
	Vec::accum(n0,n1);
	Vec::accum(n0,n2);
	area = Vec::length(n0)/6.0;
      }
      break;
    case Mesh::Mtype_plane:
      {
	/* Calculate the upwards unscaled normal(s), extract z-component. */
	area = (Vec::cross2(e1,e2)+
		Vec::cross2(e2,e0)+
		Vec::cross2(e0,e1))/6.0;
      }
      break;
    case Mesh::Mtype_sphere:
      {
	if (false) {
	/* Calculate the spherical excess. */
	Point n0, n1, n2;
	Vec::cross(n0,e1,e2);
	Vec::cross(n1,e2,e0);
	Vec::cross(n2,e0,e1);
	/* Note: std::fmod calculates signed modulus, so use other trick. */
	int nneg = 0;
	double theta0 = std::atan2(Vec::scalar(s0,n0),
				   Vec::scalar(s0,e1)
				   *Vec::scalar(s0,e2)
				   -Vec::scalar(e1,e2));
	if (theta0<0) nneg++;
	double theta1 = std::atan2(Vec::scalar(s1,n1),
				   Vec::scalar(s1,e2)
				   *Vec::scalar(s1,e0)
				   -Vec::scalar(e2,e0));
	if (theta1<0) nneg++;
	double theta2 = std::atan2(Vec::scalar(s2,n2),
				   Vec::scalar(s2,e0)
				   *Vec::scalar(s2,e1)
				   -Vec::scalar(e0,e1));
	if (theta2<0) nneg++;
	area = theta0+theta1+theta2+static_cast<double>(2*nneg-1)*M_PI;
	}

	/*
	  New formula.
	*/
	double costh =
	  1.+Vec::scalar(s0,s1)+Vec::scalar(s1,s2)+Vec::scalar(s2,s0);
	double sinth = Vec::volume(s0,s1,s2);
	area = 2.*std::atan2(sinth,costh);
	if (area<0)
	  area += 4.*M_PI;
	
	//	cout << WHEREAMI << "Areas: (" << area << "," << area2
	//	     << ") a2/a1 = " << area2/area << endl;

	/*
	double costh =
	  1.+Vec::scalar(s0,s1)+Vec::scalar(s1,s2)+Vec::scalar(s2,s0);
	double sinth = 2.*Vec::volume(s0,s1,s2)/Vec::length(e2);
	double area2 = 2.*std::atan2(sinth,costh);
	
	cout << WHEREAMI << "Areas: (" << area << "," << area2
	     << ") a2/a1 = " << area2/area << endl;
	*/

	/*
	  // L'Huilier code, don't use:
	  Dart dh(*this,t);
	  double a(edgeLength(dh));
	  dh.orbit2();
	  double b(edgeLength(dh));
	  dh.orbit2();
	  double c(edgeLength(dh));
	  double s((a+b+c)/2.0);
	  double tanE4(std::sqrt(std::tan(s/2.0)*
	                         std::tan((s-a)/2.0)*
				 std::tan((s-b)/2.0)*
				 std::tan((s-c)/2.0)));
          area = 4.0*std::atan(tanE4);
	*/
      }
      break;
    default:
      /* ERROR: This should never be reached. */
      LOG("ERROR: unhandled mesh type.");
      area = 0.0;
    }

    return area;
  }

  void Mesh::triangleBoundingBox(int t, Point& mini, Point& maxi) const
  {
    if ((t<0) || (t>=(int)nT())) {
      return;
    }

    Dart dh(Dart(*this,t));
    int v0 = dh.v();
    dh.orbit2();
    int v1 = dh.v();
    dh.orbit2();
    int v2 = dh.v();
    const Point& s0 = S_[v0];
    const Point& s1 = S_[v1];
    const Point& s2 = S_[v2];
    Mesh::triangleBoundingBox(s0,s1,s2,mini,maxi);

    if (type_ == Mesh::Mtype_sphere) {
      /* Need to take curvature into account. */
      /* Calculate a tangent point inside the triangle. */
      Point e0, e1, e2;
      Point n0;
      Vec::sum(n0,s0,s1);
      Vec::accum(n0,s2,1.0);
      Vec::rescale(n0,1.0/Vec::length(n0));
      
      /* Radially project points onto tangent plane. */
      /* s0_ = x0*s0,  (x0*s0-n0)*n0 = 0,  x0 = 1/(s0*n0) */
      Vec::scale(e0, s0, 1.0/Vec::scalar(s0,n0));
      Vec::scale(e1, s1, 1.0/Vec::scalar(s1,n0));
      Vec::scale(e2, s2, 1.0/Vec::scalar(s2,n0));

      /* Calculate bounding box in tangent plane. */
      Point mini2, maxi2;
      Mesh::triangleBoundingBox(e0,e1,e2,mini2,maxi2);

      /* Calculate joint bounding box. */
      for (int d=0; d<3; d++) {
	mini[d] = ( mini[d] < mini2[d] ? mini[d] : mini2[d] );
	maxi[d] = ( maxi[d] > maxi2[d] ? maxi[d] : maxi2[d] );
      }
    }
  }

  double Mesh::triangleArea(int t) const
  {
    if ((t<0) || (t>=(int)nT())) return 0.0;

    Dart dh(Dart(*this,t));
    int v0 = dh.v();
    dh.orbit2();
    int v1 = dh.v();
    dh.orbit2();
    int v2 = dh.v();
    return Mesh::triangleArea(S_[v0],S_[v1],S_[v2]);
  }

  /*!
    Calculate triangle circumcenter

    For planes, we use a linear combination of the vertices, with
    weights obtained from
    http://en.wikipedia.org/wiki/Circumscribed_circle#Barycentric_coordinates_from_cross-_and_dot-products   \n
    Rewriting in our notation, the weights and circumcenter are given by
    \f{align*}{
    a_0  &= - \frac{\|e_0\|^2 (e_1\cdot e_2)}{2\|e_1\times e_2\|^2}
    = \frac{\|e_0\|^2}{4 A \tan(\theta_0)} \\
    c &= a_0s_0+a_1s_1+a_2s_2
    \f}
    where formulas for \f$a_1\f$ and \f$a_2\f$ are given by index
    permutation.
    
    On the sphere, the normalised flat triangle normal is the circumcenter.
    
    \see Mesh::triangleArea
    \see Mesh::triangleCircumcircleRadius
  */
  void Mesh::triangleCircumcenter(int t, Point& c) const
  {
    if ((t<0) || (t>=(int)nT())) {
      c[0] = 0.0;
      c[1] = 0.0;
      c[2] = 0.0;
      return;
    }

    int v0 = TV_[t][0];
    int v1 = TV_[t][1];
    int v2 = TV_[t][2];
    const Point& s0 = S_[v0];
    const Point& s1 = S_[v1];
    const Point& s2 = S_[v2];
    Point e0, e1, e2;
    Vec::diff(e0,s2,s1);
    Vec::diff(e1,s0,s2);
    Vec::diff(e2,s1,s0);

    switch (type_) {
    case Mesh::Mtype_manifold:
      /* TODO: Implement? Need more manifold theory for that! */
      NOT_IMPLEMENTED;
      {
	/* Calculate centroid instead of circumcenter. Not useful for RCDT. */
	Vec::copy(c,s0);
	Vec::rescale(c,1.0/3.0);
	Vec::accum(c,s1,1.0/3.0);
	Vec::accum(c,s2,1.0/3.0);
      }
      break;
    case Mesh::Mtype_plane:
      {
	Point n0, n1, n2;
	Vec::cross(n0,e1,e2);
	Vec::cross(n1,e2,e0);
	Vec::cross(n2,e0,e1);
	Vec::accum(n0,n1);
	Vec::accum(n0,n2);
	double scale(-4.5/Vec::scalar(n0,n0));
	Vec::scale(c,s0,scale*Vec::scalar(e0,e0)*Vec::scalar(e1,e2));
	Vec::accum(c,s1,scale*Vec::scalar(e1,e1)*Vec::scalar(e2,e0));
	Vec::accum(c,s2,scale*Vec::scalar(e2,e2)*Vec::scalar(e0,e1));
      }
      break;
    case Mesh::Mtype_sphere:
      {
	/* The triangle normal is equal to the circumcenter. */
	Point tmp;
	Vec::cross(c,e1,e2);
	Vec::cross(tmp,e2,e0);
	Vec::accum(c,tmp);
	Vec::cross(tmp,e0,e1);
	Vec::accum(c,tmp);
	Vec::rescale(c,1.0/Vec::length(c));
      }
      break;
    }

    return;
  }

  /*!
    \brief Calculate the radius of the triangle circumcircle

    We use the formula given at
    http://en.wikipedia.org/wiki/Circumscribed_circle#Barycentric_coordinates_from_cross-_and_dot-products \n
    Rewriting in our notation, the radius of the circumcircle
    is given by
    \f{align*}{
    r  &= \frac{3\|e_0\| \|e_1\| \|e_2\|}{2\|n_0+n_1+n_2\|}
   \f}

    \see Mesh::triangleArea
    \see Mesh::triangleCircumcenter
   */
  double Mesh::triangleCircumcircleRadius(const Point& s0,
					  const Point& s1,
					  const Point& s2) const
  {
    Point e0, e1, e2;
    Vec::diff(e0,s2,s1);
    Vec::diff(e1,s0,s2);
    Vec::diff(e2,s1,s0);

    Point n0, n1, n2;
    Vec::cross(n0,e1,e2);
    Vec::cross(n1,e2,e0);
    Vec::cross(n2,e0,e1);
    Vec::accum(n0,n1);
    Vec::accum(n0,n2);

    double radius = ((3.0
		      *Vec::length(e0)*Vec::length(e1)
		      *Vec::length(e2))
		     /(2.0*Vec::length(n0)));

    if (type_==Mesh::Mtype_sphere) {
      radius = std::asin(radius);
    }

    return radius;
  }

  /*!
    \brief Calculate the radius of the triangle circumcircle
   */
  double Mesh::triangleCircumcircleRadius(int t) const
  {
    if ((t<0) || (t>=(int)nT())) return -1.0;

    int v0 = TV_[t][0];
    int v1 = TV_[t][1];
    int v2 = TV_[t][2];
    const Point& s0 = S_[v0];
    const Point& s1 = S_[v1];
    const Point& s2 = S_[v2];

    return triangleCircumcircleRadius(s0,s1,s2);
  }

  /*!
    \brief Calculate intersection between two edges.
   */
  double Mesh::edgeIntersection(const Point& s00, const Point& s01,
				const Point& s10, const Point& s11,
				Point& c) const
  {
    Point e0, e1, e00, e01;
    double beta(0.5);

    Vec::diff(e0,s01,s00);
    Vec::diff(e1,s11,s10);

    if (type_ == Mesh::Mtype_sphere) {
      e00.cross(s00,s01);
    } else {
      Point n(0.0,0.0,1.0);
      e00.cross(n,e0);
    }

    e01.diff(s00,s10);
    c = s10;
    beta = Vec::scalar(e00,e01)/Vec::scalar(e00,e1);
    c.accum(e1, beta);

    if (type_ == Mesh::Mtype_sphere) {
      Vec::rescale(c,1./Vec::length(c));
      beta = s10.angle(c)/s10.angle(s11);
    }

    return beta;
  }
  
  bool Mesh::triangleEdgeLengths(int t, Point& len) const
  {
    if ((t<0) || (t>=(int)nT())) return 0.0;

    Dart dh(*this,t);
    len[2] = edgeLength(dh);
    dh.orbit2();
    len[0] = edgeLength(dh);
    dh.orbit2();
    len[1] = edgeLength(dh);

    return true;
  }

  int Mesh::triangleEdgeLengthsArgMin(int t, Point& len) const
  {
    if (!Mesh::triangleEdgeLengths(t,len)) return -1;

    return (len[0]<len[1] ?
	    (len[0]<len[2] ? 0 : 2) :
	    (len[1]<len[2] ? 1 : 2));
  }

  int Mesh::triangleEdgeLengthsArgMax(int t, Point& len) const
  {
    if (!Mesh::triangleEdgeLengths(t,len)) return -1;

    return (len[0]>len[1] ?
	    (len[0]>len[2] ? 0 : 2) :
	    (len[1]>len[2] ? 1 : 2));
  }

  double Mesh::triangleShortestEdge(int t) const
  {
    Point len;
    if (!Mesh::triangleEdgeLengths(t,len)) return -1;

    return (len[0]<len[1] ?
	    (len[0]<len[2] ? len[0] : len[2]) :
	    (len[1]<len[2] ? len[1] : len[2]));
  }

  double Mesh::triangleLongestEdge(int t) const
  {
    Point len;
    if (!Mesh::triangleEdgeLengths(t,len)) return -1;

    return (len[0]>len[1] ?
	    (len[0]>len[2] ? len[0] : len[2]) :
	    (len[1]>len[2] ? len[1] : len[2]));
  }


  double Mesh::edgeEncroached(const Dart& d, const Point& s) const
  /* > --> encroached */
  {
    int t(d.t());
    if ((t<0) || (t>=(int)nT())) return -1.0;

    const Point& s0 = S_[d.v()];
    const Point& s1 = S_[d.vo()];

    /* Edge is encorached if the distance to the midpoint is smaller
       than half the straight edge length. */
    Point e;
    Vec::diff(e,s1,s0);
    Point mid;
    Vec::copy(mid,s0);
    Vec::accum(mid,s1);
    Vec::rescale(mid,0.5);
    Point smid;
    Vec::diff(smid,s,mid);
    return (Vec::length(e)/2.0-Vec::length(smid));
  }





  double Mesh::inLeftHalfspace(const Point& s0,
			       const Point& s1,
			       const Point& s) const
  {
    switch (type_) {
    case Mesh::Mtype_manifold:
      //	return predicates::orient3d(M_->S[]);
      NOT_IMPLEMENTED;
      break;
    case Mesh::Mtype_plane:
      return predicates::orient2d(s0.raw(),s1.raw(),s.raw());
      break;
    case Mesh::Mtype_sphere:
      Point zero(0.,0.,0.);
      return -predicates::orient3d(s0.raw(),s1.raw(),zero.raw(),s.raw());
      break;
    }
    /* This should never be reached. */
    return 0.0;
  }

  double Dart::inLeftHalfspace(const Point& s) const
  {
    if (isnull()) return 0.0; /* TODO: should show a warning somewhere... */
    Dart dh(*this);
    int v0(dh.v());
    dh.orbit2();
    int v1(dh.v());
    return M_->inLeftHalfspace(M_->S_[v0],M_->S_[v1],s);
  }

  double Dart::inCircumcircle(const Point& s) const
  {
    if (isnull()) return 0.0; /* TODO: should show a warning somewhere... */
    Dart dh(*this);
    int v0(dh.v());
    dh.orbit2();
    int v1(dh.v());
    dh.orbit2();
    int v2(dh.v());
    switch (M_->type_) {
    case Mesh::Mtype_manifold:
      //	return predicates::orient3d(M_->S[]);
      break;
    case Mesh::Mtype_plane:
      return predicates::incircle(M_->S_[v0].raw(),
				  M_->S_[v1].raw(),
				  M_->S_[v2].raw(),
				  s.raw());
      break;
    case Mesh::Mtype_sphere:
      return -predicates::orient3d(M_->S_[v0].raw(),
				   M_->S_[v1].raw(),
				   M_->S_[v2].raw(),
				   s.raw());
      break;
    }
    /* This should never be reached. */
    return 0.0;
  }

  bool Dart::circumcircleOK(void) const
  {
    Dart dh(*this);
    if (isnull()) return true; /* TODO: should show a warning somewhere... */
    if (onBoundary()) return true; /* Locally optimal, OK. */
    dh.orbit0rev().orbit2();
    int v(dh.v());
    //    LOG("circumcircleOK? " << *this << endl);
    //    LOG("  result0 = "
    //	      << std::scientific << inCircumcircle(M_->S_[v]) << endl);
    if (inCircumcircle(M_->S_[v]) <= MESH_EPSILON) return true;
    /* For symmetric robusness, check with the reverse dart as well: */
    dh = *this;
    dh.orbit2rev();
    v = dh.v();
    dh.orbit2();
    dh.orbit1();
    //    LOG("  result1 = "
    //	      << std::scientific << dh.inCircumcircle(M_->S_[v]) << endl);
    return (dh.inCircumcircle(M_->S_[v]) <= MESH_EPSILON);
  }


  /*! \brief Swap an edge

     \verbatim
       2         2
      /0\       /|\
     0d--1 --> 00|11
      \1/       \d/
       3         3
     \endverbatim
     Dart 0-1 --> 3-2
    
  */
  Dart Mesh::swapEdge(const Dart& d)
  {
    Dart dh(d);
    int vi;
    int v_list[4];
    int t0, t1;
    int tt_list[4];
    int tti_list[4] = {-1,-1,-1};
    if (d.edir()<0) dh.alpha1(); /* Correct dart orientation */

    /* Step 1: Store geometry information. */
    t0 = dh.t();
    vi = dh.vi();
    v_list[0] = TV_[t0][vi];
    tt_list[0] = TT_[t0][vi];
    if (use_TTi_) tti_list[0] = TTi_[t0][vi];
    dh.orbit2();
    vi = dh.vi();
    v_list[1] = TV_[t0][vi];
    tt_list[1] = TT_[t0][vi];
    if (use_TTi_) tti_list[1] = TTi_[t0][vi];
    dh.orbit2();
    v_list[2] = TV_[t0][dh.vi()];
    dh.orbit2rev().orbit0();
    t1 = dh.t();
    if (t0 == t1) { dh = d; return dh; } /* ERROR: Boundary edge */
    vi = dh.vi();
    tt_list[2] = TT_[t1][vi];
    if (use_TTi_) tti_list[2] = TTi_[t1][vi];
    dh.orbit2();
    vi = dh.vi();
    tt_list[3] = TT_[t1][vi];
    if (use_TTi_) tti_list[3] = TTi_[t1][vi];
    dh.orbit2();
    v_list[3] = TV_[t1][dh.vi()];

#ifndef FMESHER_NO_X
    if (X11_) {
      drawX11triangle(t0,false);
      drawX11triangle(t1,false);
    }
#endif

    /* Step 2: Overwrite with new triangles. */
    TV_(t0)[0] = v_list[0];
    TV_(t0)[1] = v_list[3];
    TV_(t0)[2] = v_list[2];
    TT_(t0)[0] = t1;
    TT_(t0)[1] = tt_list[1];
    TT_(t0)[2] = tt_list[2];
    if (use_TTi_) {
      TTi_(t0)[0] = 0;
      TTi_(t0)[1] = tti_list[1];
      TTi_(t0)[2] = tti_list[2];
    }
    TV_(t1)[0] = v_list[1];
    TV_(t1)[1] = v_list[2];
    TV_(t1)[2] = v_list[3];
    TT_(t1)[0] = t0;
    TT_(t1)[1] = tt_list[3];
    TT_(t1)[2] = tt_list[0];
    if (use_TTi_) {
      TTi_(t1)[0] = 0;
      TTi_(t1)[1] = tti_list[3];
      TTi_(t1)[2] = tti_list[0];
    }

    /* Step 3: Relink neighbouring triangles. */
    if (use_TTi_) {
      if (TT_[t0][1]>=0) TTi_(TT_[t0][1])[TTi_[t0][1]] = 1;
      if (TT_[t0][2]>=0) TTi_(TT_[t0][2])[TTi_[t0][2]] = 2;
      if (TT_[t1][1]>=0) TTi_(TT_[t1][1])[TTi_[t1][1]] = 1;
      if (TT_[t1][2]>=0) TTi_(TT_[t1][2])[TTi_[t1][2]] = 2;
      if (TT_[t0][1]>=0) TT_(TT_[t0][1])[TTi_[t0][1]] = t0;
      if (TT_[t0][2]>=0) TT_(TT_[t0][2])[TTi_[t0][2]] = t0;
      if (TT_[t1][1]>=0) TT_(TT_[t1][1])[TTi_[t1][1]] = t1;
      if (TT_[t1][2]>=0) TT_(TT_[t1][2])[TTi_[t1][2]] = t1;
    } else {
      if (TT_[t0][1]>=0) {
	dh = Dart(*this,t0,1,2).orbit0rev();
	dh.orbit2();
	TT_(dh.t())[dh.vi()] = t0;
      }
      if (TT_[t0][2]>=0) {
	dh = Dart(*this,t0,1,0).orbit0rev();
	dh.orbit2();
	TT_(dh.t())[dh.vi()] = t0;
      }
      if (TT_[t1][1]>=0) {
	dh = Dart(*this,t1,1,2).orbit0rev();
	dh.orbit2();
	TT_(dh.t())[dh.vi()] = t1;
      }
      if (TT_[t1][2]>=0) {
	dh = Dart(*this,t1,1,0).orbit0rev();
	dh.orbit2();
	TT_(dh.t())[dh.vi()] = t1;
      }
    }

    /* Link vertices to triangles */
    if (use_VT_) {
      set_VT_triangle(t1);
      set_VT_triangle(t0);
    }

    /* Debug code: */
    /* 
    LOG("TT is \n" << TTO());
    rebuildTT();
    LOG("TT should be \n" << TTO());
    if (use_TTi_) {
      LOG("TTi is \n" << TTiO());
      rebuildTTi();
      LOG("TTi should be \n" << TTiO());
    }
    */

    LOG("Edge swapped" << endl);
#ifndef FMESHER_NO_X
    if (X11_) {
      drawX11triangle(t0,true);
      drawX11triangle(t1,true);
      X11_->delay();
    }
#endif
    
    return Dart(*this,t0,1,1);
  }
  
  /*!
     \verbatim
     2           2
    /|\         /|\
   / | \       /1d2\
  1 0|1 3 --> 1--v--3
   \ | /       \0|3/
    \d/         \|/
     0           0
     \endverbatim
   
     Dart 0-2 --> v-2
  */
  Dart Mesh::splitEdge(const Dart& d, int v)
  {
    Dart dh(d);
    int t, vi;
    int v0, v1, v2, v3;
    int t0, t1, t2, t3;
    int tt_list[4];
    int tti_list[4] = {-1,-1,-1};
    if (d.edir()<0) dh.alpha0(); /* Correct dart orientation */

    /* Step 1: Store geometry information. */
    /* Go through t0: */
    t0 = dh.t();
    vi = dh.vi();
    v0 = TV_[t0][vi];
    tt_list[1] = TT_[t0][vi];
    if (use_TTi_) tti_list[1] = TTi_[t0][vi];
    dh.orbit2();
    vi = dh.vi();
    v2 = TV_[t0][vi];
    tt_list[0] = TT_[t0][vi];
    if (use_TTi_) tti_list[0] = TTi_[t0][vi];
    dh.orbit2();
    vi = dh.vi();
    v1 = TV_[t0][vi];
    dh.orbit2();

#ifndef FMESHER_NO_X
    if (X11_)
      drawX11triangle(t0,false);
#endif

    bool on_boundary = dh.onBoundary();
    if (!on_boundary) {
      /* Go through t1: */
      dh.orbit1();
      t1 = dh.t();
      vi = dh.vi();
      tt_list[3] = TT_[t1][vi];
      if (use_TTi_) tti_list[3] = TTi_[t1][vi];
      dh.orbit2();
      vi = dh.vi();
      tt_list[2] = TT_[t1][vi];
      if (use_TTi_) tti_list[2] = TTi_[t1][vi];
      dh.orbit2();
      vi = dh.vi();
      v3 = TV_[t1][vi];

#ifndef FMESHER_NO_X
      if (X11_)
	drawX11triangle(t1,false);
#endif
    } else {
      v3 = -1;
      tt_list[2] = -1;
      tt_list[3] = -1;
      if (use_TTi_) {
	tti_list[2] = -1;
	tti_list[3] = -1;
      }
    }

    /* Step 2: Overwrite two/one triangles, create four/two new. */
    /* t0 = t0; */
    if (!on_boundary) {
      /* t1 = t1; */
      t2 = nT();
      t3 = nT()+1;
      check_capacity(0,nT()+2);
    } else {
      t1 = nT();
      check_capacity(0,nT()+1);
      t2 = -1;
      t3 = -1;
    }
    /* t0 */
    t = t0;
    TV_(t)[0] = v;
    TV_(t)[1] = v1;
    TV_(t)[2] = v0;
    TT_(t)[0] = tt_list[0];
    TT_(t)[1] = t3;
    TT_(t)[2] = t1;
    if (use_TTi_) {
      TTi_(t)[0] = tti_list[0];
      TTi_(t)[1] = 2;
      TTi_(t)[2] = 1;
    }
    /* t1 */
    t = t1;
    TV_(t)[0] = v;
    TV_(t)[1] = v2;
    TV_(t)[2] = v1;
    TT_(t)[0] = tt_list[1];
    TT_(t)[1] = t0;
    TT_(t)[2] = t2;
    if (use_TTi_) {
      TTi_(t)[0] = tti_list[1];
      TTi_(t)[1] = 2;
      TTi_(t)[2] = 1;
    }
    if (!on_boundary) {
      /* t2 */
      t = t2;
      TV_(t)[0] = v;
      TV_(t)[1] = v3;
      TV_(t)[2] = v2;
      TT_(t)[0] = tt_list[2];
      TT_(t)[1] = t1;
      TT_(t)[2] = t3;
      if (use_TTi_) {
	TTi_(t)[0] = tti_list[2];
	TTi_(t)[1] = 2;
	TTi_(t)[2] = 1;
      }
      /* t3 */
      t = t3;
      TV_(t)[0] = v;
      TV_(t)[1] = v0;
      TV_(t)[2] = v3;
      TT_(t)[0] = tt_list[3];
      TT_(t)[1] = t2;
      TT_(t)[2] = t0;
      if (use_TTi_) {
	TTi_(t)[0] = tti_list[3];
	TTi_(t)[1] = 2;
	TTi_(t)[2] = 1;
      }
    }

    /* Step 3: Relink neighbouring triangles. */
    if (use_TTi_) {
      if (TT_[t0][0]>=0) TTi_(TT_[t0][0])[TTi_[t0][0]] = 0;
      if (TT_[t1][0]>=0) TTi_(TT_[t1][0])[TTi_[t1][0]] = 0;
      if (TT_[t0][0]>=0) TT_(TT_[t0][0])[TTi_[t0][0]] = t0;
      if (TT_[t1][0]>=0) TT_(TT_[t1][0])[TTi_[t1][0]] = t1;
      if (!on_boundary) {
	if (TT_[t2][0]>=0) TTi_(TT_[t2][0])[TTi_[t2][0]] = 0;
	if (TT_[t3][0]>=0) TTi_(TT_[t3][0])[TTi_[t3][0]] = 0;
	if (TT_[t2][0]>=0) TT_(TT_[t2][0])[TTi_[t2][0]] = t2;
	if (TT_[t3][0]>=0) TT_(TT_[t3][0])[TTi_[t3][0]] = t3;
      }
    } else {
      if (TT_[t0][0]>=0) {
	dh = Dart(*this,t0,1,1).orbit0rev();
	dh.orbit2();
	TT_(dh.t())[dh.vi()] = t0;
      }
      if (TT_[t1][0]>=0) {
	dh = Dart(*this,t1,1,1).orbit0rev();
	dh.orbit2();
	TT_(dh.t())[dh.vi()] = t1;
      }
      if (!on_boundary) {
	if (TT_[t2][0]>=0) {
	  dh = Dart(*this,t2,1,1).orbit0rev();
	  dh.orbit2();
	  TT_(dh.t())[dh.vi()] = t2;
	}
	if (TT_[t3][0]>=0) {
	  dh = Dart(*this,t3,1,1).orbit0rev();
	  dh.orbit2();
	  TT_(dh.t())[dh.vi()] = t3;
	}
      }
    }

    /* Link vertices to triangles */
    if (use_VT_) {
      if (!on_boundary) {
	set_VT_triangle(t3);
	set_VT_triangle(t2);
      }
      set_VT_triangle(t1);
      set_VT_triangle(t0);
    }

    /* Debug code: */
    /*
    LOG("TT is \n" << TTO());
    rebuildTT();
    LOG("TT should be \n" << TTO())
    if (use_TTi_) {
      LOG("TTi is \n" << TTiO());
      rebuildTTi();
      LOG("TTi should be \n" << TTiO());
    }
    */

    LOG("Edge split" << endl);
#ifndef FMESHER_NO_X
    if (X11_) {
      if (!on_boundary) {
	drawX11triangle(t3,true);
	drawX11triangle(t2,true);
      }
      drawX11triangle(t1,true);
      drawX11triangle(t0,true);
      X11_->delay();
    }
#endif
    
    return Dart(*this,t1,1,0);
  }

  /*!
     \verbatim
      2          2
     |  \       |\ \
     |   \      |2\1\
     |    1 --> | v--1
     |   /      | d0/
     | d/       |/ /
      0          0
     \endverbatim
   
     Dart 0-1 --> v-1
  */
  Dart Mesh::splitTriangle(const Dart& d, int v)
  {
    Dart dh(d);
    int t, vi;
    int v0, v1, v2;
    int t0, t1, t2;
    int tt_list[3];
    int tti_list[3] = {-1,-1,-1};
    if (d.edir()<0) dh.alpha1(); /* Correct dart orientation */

    /* Step 1: Store geometry information. */
    t = dh.t();
    vi = dh.vi();
    v0 = TV_[t][vi];
    tt_list[1] = TT_[t][vi];
    if (use_TTi_) tti_list[1] = TTi_[t][vi];
    dh.orbit2();
    t = dh.t();
    vi = dh.vi();
    v1 = TV_[t][vi];
    tt_list[2] = TT_[t][vi];
    if (use_TTi_) tti_list[2] = TTi_[t][vi];
    dh.orbit2();
    t = dh.t();
    vi = dh.vi();
    v2 = TV_[t][vi];
    tt_list[0] = TT_[t][vi];
    if (use_TTi_) tti_list[0] = TTi_[t][vi];
    dh.orbit2();

#ifndef FMESHER_NO_X
    if (X11_)
      drawX11triangle(t,false);
#endif

    /* Step 2: Overwrite one triangle, create two new. */
    t0 = t;
    t1 = nT();
    t2 = nT()+1;
    check_capacity(0,nT()+2);
    
    LOG("Capacity (V,T) = ("
	     << S_.capacity() << "," << TV_.capacity() << "), T-indices = ("
	     << t0 << "," << t1 << "," << t2 << ")" << endl);

    TV_(t0)[0] = v;
    TV_(t0)[1] = v0;
    TV_(t0)[2] = v1;
    TT_(t0)[0] = tt_list[0];
    TT_(t0)[1] = t1;
    TT_(t0)[2] = t2;
    if (use_TTi_) {
      TTi_(t0)[0] = tti_list[0];
      TTi_(t0)[1] = 2;
      TTi_(t0)[2] = 1;
    }
    TV_(t1)[0] = v;
    TV_(t1)[1] = v1;
    TV_(t1)[2] = v2;
    TT_(t1)[0] = tt_list[1];
    TT_(t1)[1] = t2;
    TT_(t1)[2] = t0;
    if (use_TTi_) {
      TTi_(t1)[0] = tti_list[1];
      TTi_(t1)[1] = 2;
      TTi_(t1)[2] = 1;
    }
    TV_(t2)[0] = v;
    TV_(t2)[1] = v2;
    TV_(t2)[2] = v0;
    TT_(t2)[0] = tt_list[2];
    TT_(t2)[1] = t0;
    TT_(t2)[2] = t1;
    if (use_TTi_) {
      TTi_(t2)[0] = tti_list[2];
      TTi_(t2)[1] = 2;
      TTi_(t2)[2] = 1;
    }

    /* Step 3: Relink neighbouring triangles. */
    if (use_TTi_) {
      if (TT_[t0][0]>=0) TTi_(TT_[t0][0])[TTi_[t0][0]] = 0;
      if (TT_[t1][0]>=0) TTi_(TT_[t1][0])[TTi_[t1][0]] = 0;
      if (TT_[t2][0]>=0) TTi_(TT_[t2][0])[TTi_[t2][0]] = 0;
      if (TT_[t0][0]>=0) TT_(TT_[t0][0])[TTi_[t0][0]] = t0;
      if (TT_[t1][0]>=0) TT_(TT_[t1][0])[TTi_[t1][0]] = t1;
      if (TT_[t2][0]>=0) TT_(TT_[t2][0])[TTi_[t2][0]] = t2;
    } else {
      if (TT_[t0][0]>=0) {
	dh = Dart(*this,t0,1,1).orbit0rev();
	dh.orbit2();
	TT_(dh.t())[dh.vi()] = t0;
      }
      if (TT_[t1][0]>=0) {
	dh = Dart(*this,t1,1,1).orbit0rev();
	dh.orbit2();
	TT_(dh.t())[dh.vi()] = t1;
      }
      if (TT_[t2][0]>=0) {
	dh = Dart(*this,t2,1,1).orbit0rev();
	dh.orbit2();
	TT_(dh.t())[dh.vi()] = t2;
      }
    }

    /* Link vertices to triangles */
    if (use_VT_) {
      set_VT_triangle(t2);
      set_VT_triangle(t1);
      set_VT_triangle(t0);
    }

    /* Debug code: */
    /*
    LOG("TT is \n" << TTO());
    rebuildTT();
    LOG("TT should be \n" << TTO());
    if (use_TTi_) {
      LOG("TTi is \n" << TTiO());
      rebuildTTi();
      LOG("TTi should be \n" << TTiO());
    }
    */
    
    LOG("Triangle split" << endl);
#ifndef FMESHER_NO_X
    if (X11_) {
      drawX11triangle(t0,true);
      drawX11triangle(t1,true);
      drawX11triangle(t2,true);
      X11_->delay();
    }
#endif

    return Dart(*this,t0,1,0);
  }





  Mesh& Mesh::unlinkEdge(const Dart& d)
  {
    Dart dh(d);
    if (!d.onBoundary()) {
      dh.orbit0rev().orbit2();
      TT_(dh.t())[dh.vi()] = -1;
      if (use_TTi_)
	TTi_(dh.t())[dh.vi()] = -1;
      dh = d;
    }
    dh.orbit2rev();
    TT_(dh.t())[dh.vi()] = -1;
    if (use_TTi_)
      TTi_(dh.t())[dh.vi()] = -1;
    
    return *this;
  }
  
  /*!
    Unlink a triangle
   */
  Mesh& Mesh::unlinkTriangle(const int t)
  {
    Dart dh(*this,t);
    unlinkEdge(dh);
    dh.orbit2();
    unlinkEdge(dh);
    dh.orbit2();
    unlinkEdge(dh);
    return *this;
  }

  Mesh& Mesh::relocateTriangle(const int t_source, const int t_target)
  {
    if (t_target == t_source)
      return *this;
    if (t_target>t_source)
      check_capacity(0,t_target+1);
    TV_(t_target)[0] = TV_[t_source][0];
    TV_(t_target)[1] = TV_[t_source][1];
    TV_(t_target)[2] = TV_[t_source][2];
    TT_(t_target)[0] = TT_[t_source][0];
    TT_(t_target)[1] = TT_[t_source][1];
    TT_(t_target)[2] = TT_[t_source][2];
    if (use_VT_) {
      if (VT_[TV_[t_target][0]] == t_source)
	VT_(TV_[t_target][0]) = t_target;
      if (VT_[TV_[t_target][1]] == t_source)
	VT_(TV_[t_target][1]) = t_target;
      if (VT_[TV_[t_target][2]] == t_source)
	VT_(TV_[t_target][2]) = t_target;
    }
    if (use_TTi_) {
      TTi_(t_target)[0] = TTi_[t_source][0];
      TTi_(t_target)[1] = TTi_[t_source][1];
      TTi_(t_target)[2] = TTi_[t_source][2];
    }
    /* Relink neighbouring TT:s. TTi is not affected by the relocation. */
    Dart dh(*this,t_target,1,0);
    if (!dh.onBoundary()) {
      dh.orbit0rev().orbit2();
      TT_(dh.t())[dh.vi()] = t_target;
    }
    dh = Dart(*this,t_target,1,1);
    if (!dh.onBoundary()) {
      dh.orbit0rev().orbit2();
      TT_(dh.t())[dh.vi()] = t_target;
    }
    dh = Dart(*this,t_target,1,2);
    if (!dh.onBoundary()) {
      dh.orbit0rev().orbit2();
      TT_(dh.t())[dh.vi()] = t_target;
    }
    
    return *this;
  }

  /*!
    Remove a triangle.

    The current implementation is slow when useVT is true.  For better
    performance when removing many triangles, set to false while
    removing.
   */
  int Mesh::removeTriangle(const int t)
  {
    if ((t<0) || (t>=(int)nT()))
      return -1;

#ifndef FMESHER_NO_X
    if (X11_) {
      drawX11triangle(t,false);
      if (!(TT_[t][0]<0)) drawX11triangle(TT_[t][0],true);
      if (!(TT_[t][1]<0)) drawX11triangle(TT_[t][1],true);
      if (!(TT_[t][2]<0)) drawX11triangle(TT_[t][2],true);
      X11_->delay();
    }
#endif

    unlinkTriangle(t);
    relocateTriangle(nT()-1,t);
    TV_.rows(nT()-1);
    /* Note: nT() was changed by the alteration of TV_. */
    TT_.rows(nT());
    if (use_TTi_)
      TTi_.rows(nT());
    if (use_VT_)
      rebuild_VT();
    return nT();
  }












  /*!
    Calculate barycentric coordinates.

  */
  void Mesh::barycentric(const Dart& d, const Point& s, Point& bary) const
  {
    Dart dh(d);
    int v0(dh.v());
    dh.orbit2();
    int v1(dh.v());
    dh.orbit2();
    int v2(dh.v());
    bary[0] = triangleArea(S_[v1],S_[v2],s);
    bary[1] = triangleArea(S_[v2],S_[v0],s);
    bary[2] = triangleArea(S_[v0],S_[v1],s);

    switch (type_) {
    case Mesh::Mtype_manifold:
      break;
    case Mesh::Mtype_plane:
      break;
    case Mesh::Mtype_sphere:
      {
	double a(triangleArea(d.t()));
	if (a <= 2.0*M_PI) { // Regular triangle
	  if (bary[0] > 2.0*M_PI) bary[0] = bary[0]-4.0*M_PI;
	  if (bary[1] > 2.0*M_PI) bary[1] = bary[1]-4.0*M_PI;
	  if (bary[2] > 2.0*M_PI) bary[2] = bary[2]-4.0*M_PI;
	} else { // Inverted/big triangle
	  if (bary[0] > a) bary[0] = bary[0]-4.0*M_PI;
	  if (bary[1] > a) bary[1] = bary[1]-4.0*M_PI;
	  if (bary[2] > a) bary[2] = bary[2]-4.0*M_PI;
	}

	/* "Official Barycentric coordinates, normalised: */
	if (false) {
	  const Point& s0 = S_[v0];
	  const Point& s1 = S_[v1];
	  const Point& s2 = S_[v2];
	  Point vol;
	  vol[0] = Vec::volume(s1,s2,s);
	  vol[1] = Vec::volume(s2,s0,s);
	  vol[2] = Vec::volume(s0,s1,s);
	  Vec::rescale(bary,1.0/(bary[0]+bary[1]+bary[2]));
	  cout << WHEREAMI << "Barycentric:\t" << bary << endl;
	  Vec::rescale(vol,1.0/Vec::volume(s0,s1,s2));
	  cout << WHEREAMI << "Unnormalised:\t" << vol << endl;
	  Vec::rescale(vol,1.0/(vol[0]+vol[1]+vol[2]));
	  cout << WHEREAMI << "Normalised:\t" << vol << endl;
	}
      }
      break;
    }

    Vec::rescale(bary,1.0/(bary[0]+bary[1]+bary[2]));
  }


  /*!
    Find the edge opposite a vertex that a straight path will pass through.

    If the point is not found, a null Dart is returned.

    \verbatim
    findPathDirection(d0,s)
     1. d0 determines the starting vertex, v0 = d0.v
     2. d = d0
     4. d.orbit0rev
     5. while (d != d0) // Stop at edge or after one full orbit.
     6.   d.orbit0rev
     7. d.orbit2
     8. onleft0 = inLeftHalfSpace(S(v0),s,S(d.v))
     9. d.orbit2
    10. onleft1 = inLeftHalfSpace(S(v0),s,S(d.v))
    11. while !(!onleft0 & onleft1) & !d.onBoundary
    12.   d.orbit0rev
    13.   if d == d0 // We've gone a full orbit without success
    14.     return null
    15.   onleft0 = onleft1
    16.   d.orbit2
    17.   onleft1 = inLeftHalfSpace(S(v0),s,S(d.v))
    18. if (!onleft0 & onleft1)
    19.   return d
    20. return null // Not found (hit boundary).
    \endverbatim
   */
  Dart Mesh::find_path_direction(const Dart& d0,
				 const Point& s,
				 const int v) const
  {
    Dart d(d0);
    if (d.isnull())
      return Dart();

    int v0(d.v());
    if (d.v() == v) // Have we found a preexisting vertex?
      return d;
    d.orbit0rev();
    while ((d != d0) && (!d.onBoundary()))
      d.orbit0rev();
    Dart d00(d);

    LOG("Finding direction to point or v starting from d00, S:" << endl
	     << "\t\t" << s << endl
	     << "\t\t" << v << endl
	     << "\t\t" << d00 << endl
	     << "\t\t" << S_[d00.v()] << endl);

    d.orbit2();
    LOG("Finding direction to point or v starting from dart:" << endl
	     << "\t\t" << s << endl
	     << "\t\t" << v << endl
	     << "\t\t" << d << endl
	     << "\t\t" << S_[d.v()] << endl
	     << "\tiLHS\t" << inLeftHalfspace(S_[v0],s,S_[d.v()]) << endl);
    if (d.v() == v) // Have we found a preexisting vertex?
      return d;
    bool onleft0(inLeftHalfspace(S_[v0],s,S_[d.v()]) >= 0.0);
    bool onleft2(d.inLeftHalfspace(s) >= -MESH_EPSILON);
    LOG(d << endl);
    d.orbit2();
    LOG("Finding direction to point or v starting from dart:" << endl
	     << "\t\t" << s << endl
	     << "\t\t" << v << endl
	     << "\t\t" << d << endl
	     << "\t\t" << S_[d.v()] << endl
	     << "\tiLHS\t" << inLeftHalfspace(S_[v0],s,S_[d.v()]) << endl);
    if (d.v() == v) // Have we found a preexisting vertex?
      return d;
    bool onleft1(inLeftHalfspace(S_[v0],s,S_[d.v()]) >= 0.0);
    LOG("Locating direction "
	     << onleft0 << onleft1 << endl);
    while (!(!onleft0 && onleft1) && (!d.onBoundary())) {
      d.orbit0rev();
      LOG(d << endl);
      if (d.v()==d00.vo()) {
	if (onleft2) {
	  LOG("Went full circle. Point found." << endl);
	  return d;
	} else {
	  LOG("Went full circle. Point not found." << endl);
	  return Dart();
	}
      }
      onleft0 = onleft1;
      onleft2 = (onleft2 && (d.inLeftHalfspace(s) >= -MESH_EPSILON));
      d.orbit2();
      onleft1 = (inLeftHalfspace(S_[v0],s,S_[d.v()]) >= 0.0);
      if (d.v() == v) // Have we found a preexisting vertex?
	return d;
      LOG("Locating direction "
	       << onleft0 << onleft1 << endl);
    }
    if (!onleft0 && onleft1) {
      d.orbit2rev();
      return d;
    }
    return Dart();
  }



  /*!
    Find the edge that a straight path will pass through.

    If the path endpoint is inside the original triangle, a null Dart
    is returned.
  */
  Dart Mesh::find_path_direction(const Point& s0,
				 const Point& s1,
				 const Dart& d0) const
  {
    if (d0.isnull())
      return Dart();
    Dart d(*this, d0.t(), 1, 0);
    
    /* Check if we're starting on a vertex, and call alternative method */
    /* onleft[i] = is triangle vertex i to the left of the line? */
    /* inside[i] = is s1 inside triangle edge i? */
    bool onleft[3];
    for (int i = 0; i < 3; i++) {
      if ( edgeLength(S_[d.v()], s0) < MESH_EPSILON ) {
	d = find_path_direction(d, s1, -1);
	/* Check that the line actually crosses the dart. */
	if (inLeftHalfspace(S_[d.v()], S_[d.vo()], s1) < 0.0) {
	  LOG("Checkpoint 4" << endl);
	  return d;
	} else {
	  return Dart();
	}
      }
      onleft[i] = (inLeftHalfspace(s0, s1, S_[d.v()]) >= 0.0);
      LOG("D=" << d << endl);
      LOG("onleft[" << i << "] = " << onleft[i] << endl);
      d.orbit2();
    }

    for (int i0 = 0; i0 < 3; i0++) {
      int i1 = (i0+1) % 3;
      if (inLeftHalfspace(S_[d.v()], S_[d.vo()], s1) < 0.0) {
	if (!onleft[i1]) {
	  LOG("Checkpoint 1" << endl);
	  d.orbit2();
	  return d;
	}
	if (onleft[i0]) {
	  LOG("Checkpoint 2" << endl);
	  d.orbit2rev();
	  return d;
	}
	LOG("Checkpoint 3" << endl);
	return d;
      }
      d.orbit2();
    }
    return Dart();
  }



  /*!
    Trace the geodesic path from a vertex to a point or another vertex.

    Return a pair of darts identifying the starting vertex, together
    with the point containing triangle, or a dart originating at the
    endpoint vertex.  If the point/vertex is not found, a null Dart is
    returned.  Priority is given to finding the vertex; if found, the
    point is disregarded.  The trace only includes darts strictly
    intersected, i.e. not the initial and final darts.

    Alg 9.1 is non-robust, and does not take the shortest route.  This
    algorithm finds and follows the straight-line intersected edges
    instead, making it suitable for use in CDT segment insertion
    algorithms.

    \verbatim
    tracePath:
     1. d0 determines the starting vertex, v0 = d0.v
     2. d = findPathDirection(d0,s) // d is now opposite v0
     3. if inLeftHalfspace(d,s) return d
     4. while !d.onBoundary
     5.   store d in path-trace
     6.   d1 = d.orbit1.orbit2rev
     7.   found1 = inLeftHalfspace(d1,s)
     8.   leavethrough2 = inLeftHalfspace(S(v0),s,S(d.v))
     9.   d2 = d1.orbit2rev
    10.   found2 = inLeftHalfspace(d2,s)
    11.   if found1 & found2 return d2
    12.   if leavethrough2, d=d2, else d=d1
    13. store d in path-trace
    14. return null
    \endverbatim
   */
  DartPair Mesh::trace_path(const Dart& d0,
			    const Point& s1,
			    const int v1,
			    DartList* trace) const
  {
    Dart dh;
    bool found, other;
    int trace_index = 0;
    if (d0.isnull())
      dh = Dart(*this,0);
    else
      dh = Dart(*this,d0.t(),1,d0.vi());
    int v0(dh.v());
    LOG("Locating point " << s1
	     << " v0=" << v0
	     << " v1=" << v1
	     << endl);

    if (v1>=(int)nV()) { /* Vertex index out of range */
      return DartPair(dh,Dart()); 
    }

    Dart d(find_path_direction(dh,s1,v1));
    LOG("Path-direction " << d << endl);
    LOG("Starting triangle " << d.t() << " ("
	     << TV_[d.t()][0] << ","
	     << TV_[d.t()][1] << ","
	     << TV_[d.t()][2] << ")"
	     << endl);
    if (d.isnull()) {
      LOG("No direction found, so is in starting triangle" << endl);
      return DartPair(dh, dh);
    }
    Dart dstart = d;
    while (dstart.v() != d0.v())
      dstart.orbit2rev();
    LOG("Starting dart " << dstart << endl);
    if ((d.v() == v1) || (d.inLeftHalfspace(s1) >= -MESH_EPSILON)) {
      LOG("Found " << d << endl);
      return DartPair(dstart,d);
    }
    while (!d.onBoundary()) {
      if (trace) {
	trace->push_back(d);
	trace_index++;
      }
      d.orbit1().orbit2rev();
      LOG("In triangle " << d << endl);
      if (d.v() == v1) {
	LOG("Found vertex at " << d << endl);
	return DartPair(dstart,d);
      }
      found = (d.inLeftHalfspace(s1) >= -MESH_EPSILON);
      other = (inLeftHalfspace(S_[v0],s1,S_[d.v()]) > 0.0);
      d.orbit2rev();
      if (found && (d.inLeftHalfspace(s1) >= -MESH_EPSILON))
	return DartPair(dstart,d);
      else
	found = false;
      if (!other)
	d.orbit2();
      LOG("Go to next triangle, from " << d << endl);
    }
    LOG("Endpoint not found "
	     << dstart << " " << d << endl);
    return DartPair(dstart,Dart());
  }


  /*!
    Trace the geodesic path between two points.

    Return a pair of darts identifying the starting and end triangles.
    If the end point is not reached due to reaching a mesh boundary
    the second dart returned will be null, but the path trace up to an
    including the mesh boundary will be populated.  The trace only
    includes darts strictly intersected.

    \verbatim
    tracePath:
     2. d = findPathDirection(s0,s1,d0)
     3. if d is null (i.e. inLeftHalfspace(d,s1)), return (d0,d0)
     4. while !d.onBoundary
     5.   store d in path-trace
     6.   d1 = d.orbit1.orbit2rev
     7.   found1 = inLeftHalfspace(d1,s1)
     8.   leavethrough2 = inLeftHalfspace(s0,s1,S(d.v))
     9.   d2 = d1.orbit2rev
    10.   found2 = inLeftHalfspace(d2,s1)
    11.   if found1 & found2 return (d0,d2)
    12.   if leavethrough2, d=d2, else d=d1
    13. store d in path-trace
    14. return (d0,null)
    \endverbatim
   */
  DartPair Mesh::trace_path(const Point& s0,
			    const Point& s1,
			    const Dart& d0,
			    DartList* trace) const
  {
    Dart dh;
    bool found, other;
    int trace_index = 0;
    if (d0.isnull()) {
      return DartPair(Dart(), Dart());
    }
    dh = Dart(*this, d0.t(), 1, d0.vi());
    LOG("Locating point s1=" << s1
	<< "from s0=" << s0
	<< endl);
    Dart dstart = dh;
    LOG("Starting dart " << dstart << endl);

    Dart d(find_path_direction(s0, s1, dstart));
    LOG("Path-direction " << d << endl);
    LOG("Starting triangle " << d.t() << " ("
	     << TV_[d.t()][0] << ","
	     << TV_[d.t()][1] << ","
	     << TV_[d.t()][2] << ")"
	     << endl);
    if (d.isnull()) {
      LOG("No direction found, so is in starting triangle" << endl);
      return DartPair(dstart, dstart);
    }
    if (d.inLeftHalfspace(s1) >= -MESH_EPSILON) {
      LOG("Found " << d << endl);
      return DartPair(dstart, dstart);
    }
    while (!d.onBoundary()) {
      if (trace) {
	trace->push_back(d);
	trace_index++;
      }
      d.orbit1().orbit2rev();
      LOG("In triangle " << d << endl);
      found = (d.inLeftHalfspace(s1) >= -MESH_EPSILON);
      other = (inLeftHalfspace(s0, s1, S_[d.v()]) > 0.0);
      d.orbit2rev();
      if (found && (d.inLeftHalfspace(s1) >= -MESH_EPSILON))
	return DartPair(dstart, d);
      else
	found = false;
      if (!other)
	d.orbit2();
      LOG("Go to next triangle, from " << d << endl);
    }
    LOG("Endpoint not found "
	     << dstart << " " << d << endl);
    if (trace) {
      trace->push_back(d);
      trace_index++;
    }
    return DartPair(dstart, Dart());
  }


  /*!
    Locate a point in the graph.

    Return a dart identifying the containing triangle.

    If the point is not found, a null Dart is returned.
   */
  Dart Mesh::locate_point(const Dart& d0, const Point& s, const int v) const
  {
    Dart dh;
    if (d0.isnull())
      dh = Dart(*this,0); /* Another option here would be to pick a
			     random triangle instead of always triangle 0.
			     For now, we leave such options to the caller.
			  */
    else
      dh = Dart(*this,d0.t(),1,d0.vi());
    return trace_path(dh,s,v).second;
  }


  /*!
    Locate an existing vertex in the graph.

    Return a dart originating at the vertex.

    If the vertex is not found, a null Dart is returned.
   */
  Dart Mesh::locate_vertex(const Dart& d0,
			  const int v) const
  {
    if ((v<0) || (v>=(int)nV())) {
      return Dart(); /* Vertex index out of range */
    }

    if (use_VT_) {
      int t = VT_[v];
      if (t<0) /* Vertex not connected to any triangles. */
	return Dart();
      if (TV_[t][0] == v)
	return Dart(*this,t,1,0);
      if (TV_[t][1] == v)
	return Dart(*this,t,1,1);
      if (TV_[t][2] == v)
	return Dart(*this,t,1,2);
      LOG("ERROR: Inconsistent data structures!" << endl);
      return Dart(); /* ERROR: Inconsistent data structures! */
    }

    Dart dh;
    if (d0.isnull())
      dh = Dart(*this,0);
    else
      dh = Dart(*this,d0.t(),1,d0.vi());
    dh = trace_path(dh,S_[v],v).second;
    if (dh.v() != v) /* Point may be found, but not the actual vertex. */
      return Dart();
    return dh;
  }










  Mesh& Mesh::quad_tesselate(const Mesh& M)
  {
    NOT_IMPLEMENTED;
    // TODO: Implement;
    clear();
    std::cout << WHEREAMI << "M.nV = " << M.nV() << std::endl;
    return *this;
  }


  Matrix3double* make_globe_points(int subsegments)
  {
    int nT = 20*subsegments*subsegments;
    int nV = 2+nT/2;
    Matrix3double* S = new Matrix3double(nV);
    Matrix3double& S_ = *S;

    int offset = 0;
    S_(offset) = Point(0.,0.,1.);
    offset += 1;

    for (int i=1; i <= subsegments; i++) {
      /* #points in this ring: 5*i */
      double colatitude = i*M_PI/(subsegments*3.);
      for (int j=0; j < 5*i ; j++) {
	double longitude = j/(5.*i)*2.*M_PI;
	S_(offset+j) = Point(std::cos(longitude)*std::sin(colatitude),
			     std::sin(longitude)*std::sin(colatitude),
			     std::cos(colatitude));
      }
      offset += 5*i;
    }

    for (int i=1; i < subsegments; i++) {
      /* #points in this ring: 5*subsegments */
      double colatitude = (subsegments+i)*M_PI/(subsegments*3.);
      for (int j=0; j < 5*subsegments ; j++) {
	double longitude = (0.5*(i%2)+j)/(5.*subsegments)*2.*M_PI;
	S_(offset+j) = Point(std::cos(longitude)*std::sin(colatitude),
			     std::sin(longitude)*std::sin(colatitude),
			     std::cos(colatitude));
      }
      offset += 5*subsegments;
    }

    for (int i=subsegments; i>0; i--) {
      /* #points in this ring: 5*i */
      double colatitude = M_PI-i*M_PI/(subsegments*3.);
      for (int j=0; j < 5*i ; j++) {
	double longitude = (0.5*(i%2)+j)/(5.*i)*2.*M_PI;
	S_(offset+j) = Point(std::cos(longitude)*std::sin(colatitude),
			     std::sin(longitude)*std::sin(colatitude),
			     std::cos(colatitude));
      }
      offset += 5*i;
    }
    S_(offset) = Point(0.,0.,-1.);

    return S;
  }

  Mesh& Mesh::make_globe(int subsegments)
  {
    TV_set(Matrix3int());
    int nV0 = (*this).nV();
    type(Mtype_sphere);
    int nT = 20*subsegments*subsegments;
    int nV = 2+nT/2;
    check_capacity(nV0+nV,nT);

    Matrix3double* S = make_globe_points(subsegments);
    S_append(*S);
    delete S;

    MeshC MC(this);
    vertexListT vertices;
    for (int v=0;v<nV;v++)
      vertices.push_back(nV0+v);
    MC.DT(vertices);

    /*
    if (use_VT_)
      update_VT_triangles(0);
    rebuildTT();
    rebuildTTi();
    */

    return *this;
  }










  MOAint3 Mesh::TVO() const { return MOAint3(TV_,nT()); }
  MOAint3 Mesh::TTO() const { return MOAint3(TT_,nT()); }
  MOAint Mesh::VTO() const { return MOAint(VT_,nV()); }
  MOAint3 Mesh::TTiO() const { return MOAint3(TTi_,nT()); }
  MOAdouble3 Mesh::SO() const { return MOAdouble3(S_,nV()); }







  void Mesh::calcQblocks(SparseMatrix<double>& C0,
			 SparseMatrix<double>& C1,
			 SparseMatrix<double>& G1,
			 SparseMatrix<double>& B1,
			 Matrix<double>& Tareas) const
  {
    C0.clear().rows(nV()).cols(nV());
    C1.clear().rows(nV()).cols(nV());
    G1.clear().rows(nV()).cols(nV());
    B1.clear().rows(nV()).cols(nV());
    Tareas.clear().cols(1).rows(nT());
    Point e[3];
    for (int t = 0; t < (int)nT(); t++) {
      const Int3Raw& tv = TV_[t].raw();
      const Point& s0 = S_[tv[0]];
      const Point& s1 = S_[tv[1]];
      const Point& s2 = S_[tv[2]];
      e[0].diff(s2,s1);
      e[1].diff(s0,s2);
      e[2].diff(s1,s0);

      PointRaw eij[3];
      for (int i=0; i<3; i++) {
	eij[i][i] = Vec::scalar(e[i],e[i]);
	for (int j=i+1; j<3; j++) {
	  eij[i][j] = Vec::scalar(e[i],e[j]);
	  eij[j][i] = eij[i][j];
	}
      }
      
      bool b[3];
      b[0] = (TT_[t][0] < 0 ? true : false);
      b[1] = (TT_[t][1] < 0 ? true : false);
      b[2] = (TT_[t][2] < 0 ? true : false);

      double a = triangleArea(t);
      Tareas(t,0) = a;
      
      /* "Flat area" better approximation for use in G-calculation. */
      double fa = Point().cross(e[0],e[1]).length()/2.0;

      double vij;
      for (int i=0; i<3; i++) {
	C0(tv[i],tv[i]) += a/3.;
	C1(tv[i],tv[i]) += a/6.;
	G1(tv[i],tv[i]) += eij[i][i]/(4.*fa);
	for (int j=i+1; j<3; j++) {
	  C1(tv[i],tv[j]) += a/12.;
	  C1(tv[j],tv[i]) += a/12.;
	  vij = eij[i][j]/(4.*fa);
	  G1(tv[i],tv[j]) += vij;
	  G1(tv[j],tv[i]) += vij;
	}
      }

      if (b[0] || b[1] || b[2]) {
	vij = -1./(4.*fa);
	for (int i=0; i<3; i++) {
	  for (int j=0; j<3; j++) {
	    for (int k=0; k<3; k++) {
	      if (b[k] && (i != k)) {
		B1(tv[i],tv[j]) += eij[k][j]*vij;
	      }
	    }
	  }
	}
      }
    }
  }
  


  void crossmultiply(const Point* ax,
		     Point* H,
		     int n)
  {
    for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
	H[i][j] = 0.0;
	for (int k=0; k<n; k++) {
	  H[i][j] += ax[k][i]*ax[k][j];
	}
      }
    }
  }

  void adjugate(const Point* H,
		Point* aH)
  {
    aH[0][0] = H[1][1]*H[2][2]-H[1][2]*H[2][1];
    aH[0][1] = H[1][2]*H[2][0]-H[1][0]*H[2][2];
    aH[0][2] = H[1][0]*H[2][1]-H[1][1]*H[2][0];
    aH[1][1] = H[0][0]*H[2][2]-H[0][2]*H[2][0];
    aH[1][2] = H[0][1]*H[2][0]-H[0][0]*H[2][1];
    aH[2][2] = H[0][0]*H[1][1]-H[0][1]*H[1][0];
    aH[1][0] = aH[0][1];
    aH[2][0] = aH[0][2];
    aH[2][1] = aH[1][2];
  }


  void Mesh::calcQblocksAni(SparseMatrix<double>& G1,
			    const Matrix<double>& gamma,
			    const Matrix<double>& vec) const
  {
    G1.clear().rows(nV()).cols(nV());
    Matrix<double> Tareas;
    Tareas.clear().cols(1).rows(nT());
    Matrix3double vec_(vec);
    Point e[3];
    double t_gamma;
    Point t_vec;
    for (int t = 0; t < (int)nT(); t++) {
      const Int3Raw& tv = TV_[t].raw();
      const Point& s0 = S_[tv[0]];
      const Point& s1 = S_[tv[1]];
      const Point& s2 = S_[tv[2]];
      e[0].diff(s2,s1);
      e[1].diff(s0,s2);
      e[2].diff(s1,s0);
      t_gamma = (gamma[tv[0]][0]+gamma[tv[1]][0]+gamma[tv[2]][0])/3.0;
      t_vec.sum(vec_(tv[0]),vec_(tv[1])).accum(vec_(tv[2]),1.0);
      t_vec.rescale(1.0/3.0);
      
      Point H[3];
      Point aH[3];
      switch (type()) {
      case Mesh::Mtype_plane:
	H[0] = Point(t_vec[0]*t_vec[0]+t_gamma*t_gamma,
		     t_vec[1]*t_vec[0],
		     0.0);
	H[1] = Point(t_vec[0]*t_vec[1],
		     t_vec[1]*t_vec[1]+t_gamma*t_gamma,
		     0.0);
	H[2] = Point(0.0, 0.0, 0.0);

	aH[0] = Point(t_vec[1]*t_vec[1]+t_gamma*t_gamma,
		      -t_vec[1]*t_vec[0],
		      0.0);
	aH[1] = Point(-t_vec[0]*t_vec[1],
		      t_vec[0]*t_vec[0]+t_gamma*t_gamma,
		      0.0);
	aH[2] = Point(0.0, 0.0, 0.0);
	break;
      case Mesh::Mtype_sphere:
      case Mesh::Mtype_manifold:
	Point ax[4];
	Point n;
	Vec::cross(ax[0],e[1],e[2]);
	Vec::cross(ax[1],e[2],e[0]);
	Vec::cross(ax[2],e[0],e[1]);
	Vec::accum(ax[0],ax[1]);
	Vec::accum(ax[0],ax[2]);
	n.scale(ax[0],1.0/ax[0].length());
	ax[1] = e[0];
	ax[1].rescale(1.0/ax[1].length());
	Vec::cross(ax[2],n,ax[1]);
	ax[2].rescale(1.0/ax[2].length());
	ax[1].rescale(t_gamma);
	ax[2].rescale(t_gamma);
	ax[0] = n;

	t_vec.accum(n,-Vec::scalar(n,t_vec));
	ax[3] = t_vec;
	crossmultiply(ax,H,4);

	LOG("H[0]\t" << H[0] << endl);
	LOG("H[1]\t" << H[1] << endl);
	LOG("H[2]\t" << H[2] << endl);

	adjugate(H,aH);

	LOG("aH[0]\t" << aH[0] << endl);
	LOG("aH[1]\t" << aH[1] << endl);
	LOG("aH[2]\t" << aH[2] << endl);

	break;
      }

      PointRaw eij[3];
      for (int i=0; i<3; i++) {
	Point eiH(Vec::scalar(aH[0],e[i]),
		  Vec::scalar(aH[1],e[i]),
		  Vec::scalar(aH[2],e[i]));
	eij[i][i] = Vec::scalar(eiH,e[i]);
	for (int j=i+1; j<3; j++) {
	  eij[i][j] = Vec::scalar(eiH,e[j]);
	  eij[j][i] = eij[i][j];
	}
      }

      /* "Flat area" better approximation for use in G-calculation. */
      double fa = Point().cross(e[0],e[1]).length()/2.0;

      double vij;
      for (int i=0; i<3; i++) {
	G1(tv[i],tv[i]) += eij[i][i]/(4.*fa);
	for (int j=i+1; j<3; j++) {
	  vij = eij[i][j]/(4.*fa);
	  G1(tv[i],tv[j]) += vij;
	  G1(tv[j],tv[i]) += vij;
	}
      }

    }
  }
  


  void Mesh::calcGradientMatrices(SparseMatrix<double>** D) const
  {
    SparseMatrix<double> D_[3];
    D_[0].clear().rows(nV()).cols(nV());
    D_[1].clear().rows(nV()).cols(nV());
    D_[2].clear().rows(nV()).cols(nV());
    Matrix<double> weights(nV(),1);
    Point e[3];
    for (int t = 0; t < (int)nT(); t++) {
      const Int3Raw& tv = TV_[t].raw();
      const Point& s0 = S_[tv[0]];
      const Point& s1 = S_[tv[1]];
      const Point& s2 = S_[tv[2]];
      e[0].diff(s2,s1);
      e[1].diff(s0,s2);
      e[2].diff(s1,s0);

      PointRaw eij[3];
      for (int i=0; i<3; i++) {
	eij[i][i] = Vec::scalar(e[i],e[i]);
	for (int j=i+1; j<3; j++) {
	  eij[i][j] = Vec::scalar(e[i],e[j]);
	  eij[j][i] = eij[i][j];
	}
      }

      
      /*
	g0 = e1-e0*(e0*e1)/(e0*e0)
	|g0| = 2*|T|/|e0|
	g0 = g0/|g0|^2 = g0*(e0*e0)/(2*|T|)^2
      */
      /* "Flat area" . */
      double fa = Point().cross(e[0],e[1]).length()/2.0;

      Point gr[3];
      gr[0] = e[1];
      gr[1] = e[2];
      gr[2] = e[0];
      gr[0].accum(e[0],-eij[0][1]/eij[0][0]);
      gr[1].accum(e[1],-eij[1][2]/eij[1][1]);
      gr[2].accum(e[2],-eij[2][0]/eij[2][2]);
      gr[0].rescale(eij[0][0]/(4.0*fa*fa));
      gr[1].rescale(eij[1][1]/(4.0*fa*fa));
      gr[2].rescale(eij[2][2]/(4.0*fa*fa));

      for (size_t i=0; i<3; i++) {
	weights(tv[i],0) += fa;
	for (size_t j=0; j<3; j++) {
	  D_[0](tv[i],tv[j]) += gr[j][0]*fa;
	  D_[1](tv[i],tv[j]) += gr[j][1]*fa;
	  D_[2](tv[i],tv[j]) += gr[j][2]*fa;
	}
      }

    }

    for (size_t i=0; i<nV(); i++) {
      weights(i,0) = 1.0/weights(i,0);
    }
    SparseMatrix<double> w = diag(weights);
    D[0] = new SparseMatrix<double>(w*D_[0]);
    D[1] = new SparseMatrix<double>(w*D_[1]);
    D[2] = new SparseMatrix<double>(w*D_[2]);
 
  }
  


  bool Mesh::save(std::string filename_s,
		  std::string filename_tv,
		  bool binary) const
  {
    return (S_.save(filename_s,IOMatrixtype_general,binary) &&
	    TV_.save(filename_tv,IOMatrixtype_general,binary));
  }

  bool Mesh::load(std::string filename_s,
		  std::string filename_tv,
		  bool binary)
  {
    Matrix3<double> s;
    Matrix3<int> tv;

    if (!s.load(filename_s,binary) || !tv.load(filename_tv,binary))
      return false;

    S_set(s);
    TV_set(tv);

    return true;
  }

  bool Mesh::save_ascii_2009(std::string filename_s,
			     std::string filename_tv) const
  {
    return (S_.save_ascii_2009(filename_s) &&
	    TV_.save_ascii_2009(filename_tv));
  }

  bool Mesh::load_ascii_2009(std::string filename_s,
			     std::string filename_tv)
  {
    Matrix3<double> s;
    Matrix3<int> tv;

    if (!s.load_ascii_2009(filename_s) || !tv.load_ascii_2009(filename_tv))
      return false;

    S_set(s);
    TV_set(tv);

    return true;
  }








  Dart& Dart::alpha0()
  {
    vi_ = (vi_ + (3+edir_) ) % 3;
    edir_ = -edir_;
    return *this;
  }

  Dart& Dart::alpha1()
  {
    edir_ = -edir_;
    return *this;
  }

  Dart& Dart::alpha2()
  {
    if (!M_->use_TTi_) {
      int vi;
      int v = M_->TV_[t_][vi_];
      int t = M_->TT_[t_][(vi_+(3-edir_))%3];
      if (t<0) return *this;
      for (vi = 0; (vi<3) && (M_->TV_[t][vi] != v); vi++) { }
      if (vi>=3) return *this; /* Error! This should never happen! */
      vi_ = vi;
      edir_ = -edir_;
      t_ = t;
    } else {
      int vi = (vi_+(3-edir_))%3;
      int t = M_->TT_[t_][vi];
      if (t<0) return *this;
      vi_ = (M_->TTi_[t_][vi]+(3-edir_))%3;
      edir_ = -edir_;
      t_ = t;
    }
    return *this;
  }

  Dart& Dart::orbit0()
  {
    int t = t_;
    alpha1();
    alpha2();
    if (t == t_) alpha1(); /* Undo; boundary. */
    return *this;
  }

  Dart& Dart::orbit1()
  {
    int t = t_;
    alpha2();
    if (t != t_) alpha0(); /* Do only if not at boundary. */
    return *this;
  }

  Dart& Dart::orbit2()
  {
    /* "alpha0(); alpha1();" would be less efficient. */
    vi_ = (vi_+(3+edir_))%3;
    return *this;
  }

  Dart& Dart::orbit0rev()
  {
    int t = t_;
    alpha2();
    if (t != t_) alpha1(); /* Do only if not at boundary. */
    return *this;
  }

  Dart& Dart::orbit1rev() /* Equivalent to orbit1() */
  {
    orbit1();
    return *this;
  }

  Dart& Dart::orbit2rev()
  {
    /* "alpha1(); alpha0();" would be less efficient. */
    vi_ = (vi_+(3-edir_))%3;
    return *this;
  }


  std::ostream& operator<<(std::ostream& output, const Mesh& M)
  {
    output << "Mesh type:\t";
    switch (M.type()) {
    case Mesh::Mtype_manifold:
      output << "Manifold (Rd)";
      break;
    case Mesh::Mtype_plane:
      output << "Plane (R2)";
      break;
    case Mesh::Mtype_sphere:
      output << "Sphere (S2)";
      break;
    }
    output << endl;
    output << "Vertices:\t" << M.nV() << endl;
    output << "Triangles:\t" << M.nT() << endl;
    output << "Options:\t"
	   << (M.useVT() ? "VT " : "")
	   << (M.useTTi() ? "TTi " : "")
#ifndef FMESHER_NO_X
	   << (M.useX11() ? "X11 " : "")
#endif
	   << endl;
    return output;
  }

  std::ostream& operator<<(std::ostream& output, const MOAint& MO)
  {
    for (int i = 0; i < (int)MO.n_; i++) {
      output << ' ' << std::right << std::setw(4)
	     << MO.M_[i];
    }
    output << endl;
    return output;
  }

  std::ostream& operator<<(std::ostream& output, const MOAint3& MO)
  {
    for (int j = 0; j<3; j++) {
      for (int i = 0; i < (int)MO.n_; i++) {
	output << ' ' << std::right << std::setw(4)
	       << MO.M_[i][j];
      }
      output << endl;
    }
    return output;
  }

  std::ostream& operator<<(std::ostream& output, const MOAdouble3& MO)
  {
    for (int i = 0; i < (int)MO.n_; i++) {
      for (int j = 0; j<3; j++)
	output << ' ' << std::right << std::setw(10) << std::scientific
	       << MO.M_[i][j];
      output << endl;
    }
    return output;
  }

  std::ostream& operator<<(std::ostream& output, const Point& MO)
  {
    output << '(';
    for (int j = 0; j<3; j++) {
      output << std::right << std::setw(10) << std::scientific << MO[j];
      if (j<2)
	output << ',';
    }
    output << ')';
    return output;
  }

  std::ostream& operator<<(std::ostream& output, const Dart& d)
  {
    output << "D=("
	   <<std::right << std::setw(1) << d.t_
	   << std::right << std::setw(2) << d.edir_
	   << std::right << std::setw(2) << d.vi_
	   << ")";
    if ((!d.isnull()) && (d.t_<(int)d.M()->nT())) {
      output << " EV=("
	     << d.M()->TV(d.t_)[d.vi_]
	     << ","
	     << d.M()->TV(d.t_)[(d.vi_+(3+d.edir_))%3]
	     << ")";
      output << " TV=("
	     << d.M()->TV(d.t_)[0]
	     << ","
	     << d.M()->TV(d.t_)[1]
	     << ","
	     << d.M()->TV(d.t_)[2]
	     << ")";
      output << " TT=("
	     << d.M()->TT(d.t_)[0]
	     << ","
	     << d.M()->TT(d.t_)[1]
	     << ","
	     << d.M()->TT(d.t_)[2]
	     << ")";
    }
      
    return output;
  }


} /* namespace fmesh */
