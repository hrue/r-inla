#include <cstddef>
#include <cstring>
#include <set>
#include <map>
#include <sstream>

#include "predicates.h"

#include "mesh.h"

#define NOT_IMPLEMENTED std::cout << "Not implemented: \"" \
  << __PRETTY_FUNCTION__ << std::endl;

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
			     Vcap_(V_capacity), Tcap_(2*V_capacity),
			     nV_(0), nT_(0),
			     use_VT_(use_VT), use_TTi_(use_TTi)
  {
    if (Vcap_ > 0) {
      TV_ = new int[Vcap_][3];
      TT_ = new int[Tcap_][3];
      if (use_VT_)
	VT_ = new int[Vcap_];
      else
	VT_ = NULL;
      if (use_TTi_)
	TTi_ = new int[Tcap_][3];
      else
	TTi_ = NULL;
      S_ = new double[Vcap_][3];
    } else {
      TV_ = NULL;
      TT_ = NULL;
      VT_ = NULL;
      TTi_ = NULL;
      S_ = NULL;
    }
  };

  Mesh::~Mesh()
  {
    clear();
  }

  Mesh& Mesh::clear()
  {
    Vcap_ = 0;
    Tcap_ = 0;
    nV_ = 0;
    nT_ = 0;
    use_VT_ = false;
    use_TTi_ = false;
    if (TV_) { delete[] TV_; TV_ = NULL; }
    if (TT_) { delete[] TT_; TT_ = NULL; }
    if (TTi_) { delete[] TTi_; TTi_ = NULL; }
    if (S_) { delete[] S_; S_ = NULL; }
    return *this;
  }

  Mesh& Mesh::check_capacity(size_t nVc, size_t nTc)
  {
    if ((nVc <= Vcap_) && (nTc <= Tcap_))
      return *this;
    while ((nVc > Vcap_) || (nTc > Tcap_)) {
      Vcap_ = Vcap_+Mesh_V_capacity_step_size;
      Tcap_ = 2*Vcap_;
    }

    int (*TV)[3] = new int[Vcap_][3];
    int (*TT)[3] = new int[Tcap_][3];
    int (*VT) = NULL;
    if (use_VT_) VT = new int[Vcap_];
    int (*TTi)[3] = NULL;
    if (use_TTi_) TTi = new int[Tcap_][3];
    double (*S)[3] = new double[Vcap_][3];

    if (TV_) {
      if (TV_) memcpy(TV,TV_,sizeof(int)*nV_*3);
      if (TT_) memcpy(TT,TT_,sizeof(int)*nT_*3);
      if (VT_) memcpy(VT,VT_,sizeof(int)*nV_);
      if (TTi_) memcpy(TTi,TTi_,sizeof(int)*nT_*3);
      if (S_) memcpy(S,S_,sizeof(double)*nV_*3);
      if (TV_) delete[] TV_;
      if (TT_) delete[] TT_;
      if (VT_) delete[] VT_;
      if (TTi_) delete[] TTi_;
      if (S_) delete[] S_;
    }

    TV_ = TV;
    TT_ = TT;
    VT_ = VT;
    TTi_ = TTi;
    S_ = S;
    return *this;
  };



  Mesh& Mesh::rebuildTT()
  {
    typedef std::pair<int,int> E_Type;
    typedef std::map<E_Type,int> ET_Type;
    int t, vi;
    int* TVt;
    E_Type E0,E1;
    ET_Type::const_iterator Ei;
    ET_Type ET;
    /* Pass 1: */
    for (t=0; t<(int)nT_; t++) {
      TVt = TV_[t];
      for (vi=0; vi<3; vi++) {
	E0 = std::pair<int,int>(TVt[(vi+1)%3],TVt[(vi+2)%3]);
	E1 = std::pair<int,int>(TVt[(vi+2)%3],TVt[(vi+1)%3]);
	Ei = ET.find(E1);
	if (Ei != ET.end()) { /* Found neighbour */
	  TT_[t][vi] = Ei->second;
	} else { /* Either on boundary, or not added yet. */
	  TT_[t][vi] = -1;
	}
	ET.insert(ET_Type::value_type(E0,t));
      }
    }
    /*
    std::cout << TTO() << std::endl;
    for (Ei=ET.begin();Ei!=ET.end();Ei++) {
      std::cout << Ei->first.first << ' '
		<< Ei->first.second << ' '
		<< Ei->second << std::endl;
    }
    */

    /* Pass 2: */
    for (t=0; t<(int)nT_; t++) {
      TVt = TV_[t];
      for (vi=0; vi<3; vi++) {
	if (TT_[t][vi]>=0) continue;
	E1 = std::pair<int,int>(TVt[(vi+2)%3],TVt[(vi+1)%3]);
	Ei = ET.find(E1);
	if (Ei != ET.end()) { /* Found neighbour */
	  TT_[t][vi] = Ei->second;
	}
      }
    }

    return *this;
  }



  Mesh& Mesh::updateVT(const int v, const int t)
  {
    if ((!use_VT_) || (v>=(int)nV_) || (t>=(int)nT_) || (VT_[v]<0))
      return *this;
    VT_[v] = t;
    return *this;
  }

  Mesh& Mesh::setVT(const int v, const int t)
  {
    if ((!use_VT_) || (v>=(int)nV_) || (t>=(int)nT_))
      return *this;
    VT_[v] = t;
    return *this;
  }

  Mesh& Mesh::updateVTtri(const int t)
  {
    int vi;
    if ((!use_VT_) || (t>=(int)nT_) || (t<0))
      return *this;
    for (vi=0; vi<3; vi++)
      updateVT(TV_[t][vi],t);
    return *this;
  }

  Mesh& Mesh::setVTtri(const int t)
  {
    int vi;
    if ((!use_VT_) || (t>=(int)nT_) || (t<0))
      return *this;
    for (vi=0; vi<3; vi++)
      setVT(TV_[t][vi],t);
    return *this;
  }

  Mesh& Mesh::updateVTtri_private(const int t0)
  {
    if (!use_VT_) return *this;
    int t, vi;
    for (t=t0; t<(int)nT_; t++)
      for (vi=0; vi<3; vi++)
	updateVT(TV_[t][vi],t);
    return *this;
  }

  Mesh& Mesh::setVTv_private(const int v0)
  {
    if (!use_VT_) return *this;
    int v;
    for (v=v0; v<(int)nV_; v++)
      setVT(v,-1);
    return *this;
  }

  Mesh& Mesh::rebuildVT()
  {
    if (!use_VT_) {
      if (VT_) {
	delete[] VT_;
	VT_ = NULL;
      }
      return *this;
    }
    if (!Vcap_)
      return *this;
    if (!VT_)
      VT_ = new int[Vcap_];
    setVTv_private(0);
    updateVTtri_private(0);
    return *this;
  }

  Mesh& Mesh::rebuildTTi()
  {
    int t, vi, v, t2, vi2;
    if (!use_TTi_) {
      if (TTi_) {
	delete[] TTi_;
	TTi_ = NULL;
      }
      return *this;
    }
    if (!Tcap_)
      return *this;
    if (!TTi_)
      TTi_ = new int[Tcap_][3];
    for (t=0; t<(int)nT_; t++) {
      for (vi=0; vi<3; vi++) {
	v = TV_[t][vi];
	t2 = TT_[t][(vi+2)%3];
	if (t2>=0) {
	  for (vi2 = 0; (vi2<3) && (TV_[t2][vi2] != v); vi2++) { }
	  if (vi2<3) {
	    TTi_[t][(vi+2)%3] = (vi2+1)%3;
	  } else {
	    /* Error! This should never happen! */
	    std::cout << "ERROR\n";
	  }
	} else {
	  TTi_[t][(vi+2)%3] = -1;
	}
      }
    }
    return *this;
  }


  Mesh& Mesh::useVT(bool use_VT)
  {
    if (use_VT_ != use_VT) {
      if ((!use_VT_) && (VT_)) {
	/* This shouldn't happen. */
	delete[] VT_;
	VT_ = NULL;
      }
      use_VT_ = use_VT;
      rebuildVT();
    }
    return *this;
  }

  Mesh& Mesh::useTTi(bool use_TTi)
  {
    if (use_TTi_ != use_TTi) {
      if ((!use_TTi_) && (TTi_)) {
	/* This shouldn't happen. */
	delete[] TTi_;
	TTi_ = NULL;
      }
      use_TTi_ = use_TTi;
      rebuildTTi();
    }
    return *this;
  }

  Mesh& Mesh::useX11(bool use_X11)
  {
    if (use_X11_ != use_X11) {
      if (use_X11) { /* Init. */
	X11_ = new Xtmpl;
	X11_->open("fmesher::Mesh",500,500);
	use_X11_ = true;
      } else { /* Destroy. */
	X11_->close();
	use_X11_ = false;
      }
    }
    return *this;
  }



  Mesh& Mesh::S_set(const double (*S)[3], int nV)
  {
    nV_ = 0; /* Avoid possible unnecessary copy. */
    S_append(S,nV);
    return *this;
  }
  
  Mesh& Mesh::TV_set(const int (*TV)[3], int nT)
  {
    nT_ = 0; /* Avoid possible unnecessary copy. */
    TV_append(TV,nT);
    return *this;
  }

  Mesh& Mesh::S_append(const double (*S)[3], int nV)
  {
    check_capacity(nV_+nV,0);
    memcpy(S_+nV_,S,sizeof(double)*nV*3);
    nV_ += nV;
    if (use_VT_)
      setVTv_private(nV_-nV);
    return *this;
  }

  void Mesh::redrawX11()
  {
    if (!use_X11_) return;

    int v0, v1, v;
    double s[3][3];
    double s0[3];

    X11_->clear();
    for (int t=0;t<(int)nT_;t++) {
      s0[0] = 0.0;
      s0[1] = 0.0;
      s0[2] = 0.0;
      for (int vi=0;vi<3;vi++) {
	v = TV_[t][vi];
	for (int dim=0;dim<3;dim++) {
	  s[vi][dim] = S_[v][dim];
	  s0[dim] += s[vi][dim]/3;
	}
      }
      if (type_==Mtype_sphere) {
	double r0[3];
	double r1[3];
	double n[3];
	for (int dim=0;dim<3;dim++) {
	  r0[dim] = s[1][dim]-s[0][dim];
	  r1[dim] = s[2][dim]-s[0][dim];
	}
	n[0] = r0[1]*r1[2]-r0[2]*r1[1];
	n[1] = r0[2]*r1[0]-r0[0]*r1[2];
	n[2] = r0[0]*r1[1]-r0[1]*r1[0];
	if (n[2]<0) continue;
      }
      /* Draw triangle slightly closer to center. */
      for (int vi=0;vi<3;vi++)
	for (int dim=0;dim<3;dim++)
	  s[vi][dim] = (s[vi][dim]-s0[dim])*0.975+s0[dim];
      X11_->lineFG(s[0],s[1]);
      X11_->lineFG(s[1],s[2]);
      X11_->lineFG(s[2],s[0]);
      /* Draw vertex indices even closer to center. */
      for (int vi=0;vi<3;vi++)
	for (int dim=0;dim<3;dim++)
	  s[vi][dim] = (s[vi][dim]-s0[dim])*0.8+s0[dim];
      for (int vi=0;vi<3;vi++) {
	std::ostringstream ss;
	ss << "(" << TV_[t][vi] << "," << TT_[t][vi] << ")";
	X11_->text(s[vi],ss.str());
      }
      /* Draw triangle indices at center. */
      {
	std::ostringstream ss;
	ss << "(" << t << ")";
	X11_->text(s0,ss.str());
      }
    }
  }
  
  Mesh& Mesh::TV_append(const int (*TV)[3], int nT)
  {
    check_capacity(0,nT_+nT);
    memcpy(TV_+nT_,TV,sizeof(int)*nT*3);
    nT_ += nT;
    if (use_VT_)
      updateVTtri_private(nT_-nT);
    rebuildTT();
    rebuildTTi();
    if (use_X11_)
      redrawX11();
    xtmpl_press_ret("TV appended");
    return *this;
  }


  double Mesh::encroachedQuality(const Dart& d) const
  {
    /* TODO: Implement. */
    NOT_IMPLEMENTED;

    return -1.0; /* <=0 --> not encroached */
  }

  double Mesh::skinnyQuality(int t) const
  {
    /* TODO: Implement. */
    NOT_IMPLEMENTED;

    return 0.0;
  }

  double Mesh::bigQuality(int t) const
  {
    /* TODO: Implement. */
    NOT_IMPLEMENTED;

    return 0.0;
  }


  double Traits::inLeftHalfspace(const Dart& d, const double s[3])
  {
    Dart dhelper = d;
    const Mesh *M = d.M();
    int v0, v1;
    if (d.isnull()) return 0.0; /* TODO: should show a warning somewhere... */
    const int* tp = M->TV()[dhelper.t()];
    v0 = tp[dhelper.vi()];
    dhelper.orbit2();
    v1 = tp[dhelper.vi()];
    switch (M->type()) {
    case Mesh::Mtype_manifold:
      //	return predicates::orient3d(M_->S[]);
      break;
    case Mesh::Mtype_plane:
      return predicates::orient2d(M->S()[v0],M->S()[v1],s);
      break;
    case Mesh::Mtype_sphere:
      Point zero = {0.,0.,0.};
      return -predicates::orient3d(M->S()[v0],M->S()[v1],zero,s);
      break;
    }
    /* This should never be reached. */
    return 0.0;
  }

  double Traits::inCircumcircle(const Dart& d, const double s[3])
  {
    Dart dhelper = d;
    const Mesh *M = d.M();
    int v0, v1, v2;
    if (d.isnull()) return 0.0; /* TODO: should show a warning somewhere... */
    const int* tp = M->TV()[dhelper.t()];
    v0 = tp[dhelper.vi()];
    dhelper.orbit2();
    v1 = tp[dhelper.vi()];
    dhelper.orbit2();
    v2 = tp[dhelper.vi()];
    switch (M->type()) {
    case Mesh::Mtype_manifold:
      //	return predicates::orient3d(M_->S[]);
      break;
    case Mesh::Mtype_plane:
      return predicates::incircle(M->S()[v0],M->S()[v1],M->S()[v2],s);
      break;
    case Mesh::Mtype_sphere:
      return -predicates::orient3d(M->S()[v0],M->S()[v1],M->S()[v2],s);
      break;
    }
    /* This should never be reached. */
    return 0.0;
  }

  bool Traits::circumcircleOK(const Dart& d)
  {
    Dart dhelper = d;
    const Mesh *M = d.M();
    int v;
    double result;
    if (d.isnull()) return true; /* TODO: should show a warning somewhere... */
    if (d.onBoundary()) return true; /* Locally optimal, OK. */
    dhelper.orbit0rev().alpha0();
    v = M->TV()[dhelper.t()][dhelper.vi()];
    result = Traits::inCircumcircle(d,M->S()[v]);
    std::cout << "Dart=" << d
	      << " Node=" << v
	      << std::scientific << " result=" << result
	      << std::endl;
    if  (result > MESH_EPSILON)
      return false;
    /* For robusness, check with the reverse dart as well: */
    dhelper = d;
    dhelper.orbit2rev();
    v = M->TV()[dhelper.t()][dhelper.vi()];
    dhelper.orbit2();
    dhelper.orbit1();
    result = Traits::inCircumcircle(dhelper,M->S()[v]);
    std::cout << "Dart=" << dhelper
	      << " Node=" << v
	      << std::scientific << " result=" << result
	      << std::endl;
    if  (result > MESH_EPSILON)
      return false;
    return true;
  }


  /*! \brief Swap an edge

     \verbatim
       2         2
      /0\       /|\
     0---1 --> 00|11
      \1/       \|/
       3         3
     \endverbatim
     Dart 0-1 --> 3-2
    
  */
  Dart Mesh::swapEdge(const Dart& d)
  {
    Dart dhelper = d;
    int vi;
    int v_list[4];
    int t0, t1;
    int tt_list[4];
    int tti_list[4];
    if (d.edir()<0) dhelper.alpha1(); /* Correct dart orientation */

    /* Step 1: Store geometry information. */
    t0 = dhelper.t();
    vi = dhelper.vi();
    v_list[0] = TV_[t0][vi];
    tt_list[0] = TT_[t0][vi];
    if (use_TTi_) tti_list[0] = TTi_[t0][vi];
    dhelper.orbit2();
    vi = dhelper.vi();
    v_list[1] = TV_[t0][vi];
    tt_list[1] = TT_[t0][vi];
    if (use_TTi_) tti_list[1] = TTi_[t0][vi];
    dhelper.orbit2();
    v_list[2] = TV_[t0][dhelper.vi()];
    dhelper.orbit2rev().orbit0();
    t1 = dhelper.t();
    if (t0 == t1) { dhelper = d; return dhelper; } /* ERROR: Boundary edge */
    vi = dhelper.vi();
    tt_list[2] = TT_[t1][vi];
    if (use_TTi_) tti_list[2] = TTi_[t1][vi];
    dhelper.orbit2();
    vi = dhelper.vi();
    tt_list[3] = TT_[t1][vi];
    if (use_TTi_) tti_list[3] = TTi_[t1][vi];
    dhelper.orbit2();
    v_list[3] = TV_[t1][dhelper.vi()];

    /* Step 2: Overwrite with new triangles. */
    TV_[t0][0] = v_list[0];
    TV_[t0][1] = v_list[3];
    TV_[t0][2] = v_list[2];
    TT_[t0][0] = t1;
    TT_[t0][1] = tt_list[1];
    TT_[t0][2] = tt_list[2];
    if (use_TTi_) {
      TTi_[t0][0] = 0;
      TTi_[t0][1] = tti_list[1];
      TTi_[t0][2] = tti_list[2];
    }
    TV_[t1][0] = v_list[1];
    TV_[t1][1] = v_list[2];
    TV_[t1][2] = v_list[3];
    TT_[t1][0] = t0;
    TT_[t1][1] = tt_list[3];
    TT_[t1][2] = tt_list[0];
    if (use_TTi_) {
      TTi_[t1][0] = 0;
      TTi_[t1][1] = tti_list[3];
      TTi_[t1][2] = tti_list[0];
    }

    /* Step 3: Relink neighbouring triangles. */
    if (use_TTi_) {
      if (TT_[t0][1]>=0) TT_[TT_[t0][1]][TTi_[t0][1]] = t0;
      if (TT_[t0][2]>=0) TT_[TT_[t0][2]][TTi_[t0][2]] = t0;
      if (TT_[t1][1]>=0) TT_[TT_[t1][1]][TTi_[t1][1]] = t1;
      if (TT_[t1][2]>=0) TT_[TT_[t1][2]][TTi_[t1][2]] = t1;
    } else {
      if (TT_[t0][1]>=0) {
	dhelper = Dart(*this,t0,1,2).orbit0rev();
	dhelper.orbit2();
	TT_[dhelper.t()][dhelper.vi()] = t0;
      }
      if (TT_[t0][2]>=0) {
	dhelper = Dart(*this,t0,1,0).orbit0rev();
	dhelper.orbit2();
	TT_[dhelper.t()][dhelper.vi()] = t0;
      }
      if (TT_[t1][1]>=0) {
	dhelper = Dart(*this,t1,1,2).orbit0rev();
	dhelper.orbit2();
	TT_[dhelper.t()][dhelper.vi()] = t1;
      }
      if (TT_[t1][2]>=0) {
	dhelper = Dart(*this,t1,1,0).orbit0rev();
	dhelper.orbit2();
	TT_[dhelper.t()][dhelper.vi()] = t1;
      }
    }

    /* Link vertices to triangles */
    if (use_VT_) {
      setVTtri(t1);
      setVTtri(t0);
    }

    std::cout << "TT is \n" << TTO();
    rebuildTT();
    std::cout << "TT should be \n" << TTO();
    if (use_TTi_) {
      std::cout << "TTi is \n" << TTiO();
      rebuildTTi();
      std::cout << "TTi should be \n" << TTiO();
    }
    
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
    Dart dhelper = d;
    int t, vi;
    int v0, v1, v2, v3;
    int t0, t1, t2, t3;
    int tt_list[4];
    int tti_list[4];
    if (d.edir()<0) dhelper.alpha0(); /* Correct dart orientation */

    /* Step 1: Store geometry information. */
    /* Go through t0: */
    t0 = dhelper.t();
    vi = dhelper.vi();
    v0 = TV_[t0][vi];
    tt_list[1] = TT_[t0][vi];
    if (use_TTi_) tti_list[1] = TTi_[t0][vi];
    dhelper.orbit2();
    vi = dhelper.vi();
    v2 = TV_[t0][vi];
    tt_list[0] = TT_[t0][vi];
    if (use_TTi_) tti_list[0] = TTi_[t0][vi];
    dhelper.orbit2();
    vi = dhelper.vi();
    v1 = TV_[t0][vi];
    dhelper.orbit2();

    bool on_boundary = dhelper.onBoundary();
    if (!on_boundary) {
      /* Go through t1: */
      dhelper.orbit1();
      t1 = dhelper.t();
      vi = dhelper.vi();
      tt_list[3] = TT_[t1][vi];
      if (use_TTi_) tti_list[3] = TTi_[t1][vi];
      dhelper.orbit2();
      vi = dhelper.vi();
      tt_list[2] = TT_[t1][vi];
      if (use_TTi_) tti_list[2] = TTi_[t1][vi];
      dhelper.orbit2();
      vi = dhelper.vi();
      v3 = TV_[t1][vi];
    } else {
      v3 = -1;
      tt_list[2] = -1;
      tt_list[3] = -1;
      if (use_TTi_) {
	tti_list[2] = -1;
	tti_list[3] = -1;
      }
    }

    /* Step 2: Overwrite one/two triangles, create two/four new. */
    /* t0 = t0; */
    if (on_boundary) {
      t1 = nT_;
      check_capacity(0,nT_+1);
      t2 = -1;
      t3 = -1;
    } else {
      /* t1 = t1; */
      t2 = nT_;
      t3 = nT_+1;
      check_capacity(0,nT_+2);
    }
    /* t0 */
    t = t0;
    TV_[t][0] = v;
    TV_[t][1] = v1;
    TV_[t][2] = v0;
    TT_[t][0] = tt_list[0];
    TT_[t][1] = t3;
    TT_[t][2] = t1;
    if (use_TTi_) {
      TTi_[t][0] = tti_list[0];
      TTi_[t][1] = 2;
      TTi_[t][2] = 1;
    }
    /* t1 */
    t = t1;
    TV_[t][0] = v;
    TV_[t][1] = v2;
    TV_[t][2] = v1;
    TT_[t][0] = tt_list[1];
    TT_[t][1] = t0;
    TT_[t][2] = t2;
    if (use_TTi_) {
      TTi_[t][0] = tti_list[1];
      TTi_[t][1] = 2;
      TTi_[t][2] = 1;
    }
    if (!on_boundary) {
      /* t2 */
      t = t2;
      TV_[t][0] = v;
      TV_[t][1] = v3;
      TV_[t][2] = v2;
      TT_[t][0] = tt_list[2];
      TT_[t][1] = t1;
      TT_[t][2] = t3;
      if (use_TTi_) {
	TTi_[t][0] = tti_list[2];
	TTi_[t][1] = 2;
	TTi_[t][2] = 1;
      }
      /* t3 */
      t = t3;
      TV_[t][0] = v;
      TV_[t][1] = v0;
      TV_[t][2] = v3;
      TT_[t][0] = tt_list[3];
      TT_[t][1] = t2;
      TT_[t][2] = t0;
      if (use_TTi_) {
	TTi_[t][0] = tti_list[3];
	TTi_[t][1] = 2;
	TTi_[t][2] = 1;
      }
    }

    /* Step 3: Relink neighbouring triangles. */
    if (use_TTi_) {
      if (TT_[t0][0]>=0) TT_[TT_[t0][0]][TTi_[t0][0]] = t0;
      if (TT_[t1][0]>=0) TT_[TT_[t1][0]][TTi_[t1][0]] = t1;
      if (!on_boundary) {
	if (TT_[t2][0]>=0) TT_[TT_[t2][0]][TTi_[t2][0]] = t2;
	if (TT_[t3][0]>=0) TT_[TT_[t3][0]][TTi_[t3][0]] = t3;
      }
    } else {
      if (TT_[t0][0]>=0) {
	dhelper = Dart(*this,t0,1,1).orbit0rev();
	dhelper.orbit2();
	TT_[dhelper.t()][dhelper.vi()] = t0;
      }
      if (TT_[t1][0]>=0) {
	dhelper = Dart(*this,t1,1,1).orbit0rev();
	dhelper.orbit2();
	TT_[dhelper.t()][dhelper.vi()] = t1;
      }
      if (!on_boundary) {
	if (TT_[t2][0]>=0) {
	  dhelper = Dart(*this,t2,1,1).orbit0rev();
	  dhelper.orbit2();
	  TT_[dhelper.t()][dhelper.vi()] = t2;
	}
	if (TT_[t3][0]>=0) {
	  dhelper = Dart(*this,t3,1,1).orbit0rev();
	  dhelper.orbit2();
	  TT_[dhelper.t()][dhelper.vi()] = t3;
	}
      }
    }

    /* Step 4: Update triangle count. */
    if (on_boundary)
      nT_ = nT_+1;
    else
      nT_ = nT_+2;
  
    /* Link vertices to triangles */
    if (use_VT_) {
      if (!on_boundary) {
	setVTtri(t3);
	setVTtri(t2);
      }
      setVTtri(t1);
      setVTtri(t0);
    }

    std::cout << "TT is \n" << TTO();
    rebuildTT();
    std::cout << "TT should be \n" << TTO();
    if (use_TTi_) {
      std::cout << "TTi is \n" << TTiO();
      rebuildTTi();
      std::cout << "TTi should be \n" << TTiO();
    }
    
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
    Dart dhelper = d;
    int t, vi, i;
    //    int v0, v1, v2;
    int t0, t1, t2;
    int v_list[3];
    int tt_list[3];
    int tti_list[3];
    if (d.edir()<0) dhelper.alpha1(); /* Correct dart orientation */

    /* Step 1: Store geometry information. */
    t = dhelper.t();
    for (i=0;i<3;i++) {
      vi = dhelper.vi();
      v_list[i] = TV_[t][vi];
      tt_list[i] = TT_[t][vi];
      if (use_TTi_) tti_list[i] = TTi_[t][vi];
      dhelper.orbit2();
    }

    /* Step 2: Overwrite one triangles, create two new. */
    t0 = t;
    t1 = nT_;
    t2 = nT_+1;
    check_capacity(0,nT_+2);
    TV_[t0][0] = v;
    TV_[t0][1] = v_list[1];
    TV_[t0][2] = v_list[2];
    TT_[t0][0] = tt_list[0];
    TT_[t0][1] = t1;
    TT_[t0][2] = t2;
    if (use_TTi_) {
      TTi_[t0][0] = tti_list[0];
      TTi_[t0][1] = 1;
      TTi_[t0][2] = 2;
    }
    TV_[t1][0] = v;
    TV_[t1][1] = v_list[2];
    TV_[t1][2] = v_list[0];
    TT_[t1][0] = tt_list[1];
    TT_[t1][1] = t2;
    TT_[t1][2] = t0;
    if (use_TTi_) {
      TTi_[t1][0] = tti_list[1];
      TTi_[t1][1] = 2;
      TTi_[t1][2] = 0;
    }
    TV_[t2][0] = v;
    TV_[t2][1] = v_list[0];
    TV_[t2][2] = v_list[1];
    TT_[t2][0] = tt_list[2];
    TT_[t2][1] = t0;
    TT_[t2][2] = t1;
    if (use_TTi_) {
      TTi_[t2][0] = tti_list[2];
      TTi_[t2][1] = 0;
      TTi_[t2][2] = 1;
    }

    /* Step 3: Relink neighbouring triangles. */
    if (use_TTi_) {
      if (TT_[t0][0]>=0) TT_[TT_[t0][0]][TTi_[t0][0]] = t0;
      if (TT_[t1][0]>=0) TT_[TT_[t1][0]][TTi_[t1][0]] = t1;
      if (TT_[t2][0]>=0) TT_[TT_[t2][0]][TTi_[t2][0]] = t2;
    } else {
      if (TT_[t0][0]>=0) {
	dhelper = Dart(*this,t0,1,1).orbit0rev();
	dhelper.orbit2();
	TT_[dhelper.t()][dhelper.vi()] = t0;
      }
      if (TT_[t1][0]>=0) {
	dhelper = Dart(*this,t1,1,1).orbit0rev();
	dhelper.orbit2();
	TT_[dhelper.t()][dhelper.vi()] = t1;
      }
      if (TT_[t2][0]>=0) {
	dhelper = Dart(*this,t2,1,1).orbit0rev();
	dhelper.orbit2();
	TT_[dhelper.t()][dhelper.vi()] = t2;
      }
    }

    /* Step 4: Update triangle count. */
    nT_ = nT_+2;

    /* Link vertices to triangles */
    if (use_VT_) {
      setVTtri(t2);
      setVTtri(t1);
      setVTtri(t0);
    }

    std::cout << "TT is \n" << TTO();
    rebuildTT();
    std::cout << "TT should be \n" << TTO();
    if (use_TTi_) {
      std::cout << "TTi is \n" << TTiO();
      rebuildTTi();
      std::cout << "TTi should be \n" << TTiO();
    }
    
    return Dart(*this,t0,1,0);
  }





  std::ostream& operator<<(std::ostream& output, const Mesh& M)
  {
    //    output << "S =\n" << M.SO();
    output << "TV =\n" << M.TVO();
    output << "TT =\n" << M.TTO();
    if (M.useVT())
      output << "VT =\n" << M.VTO();
    if (M.useTTi())
      output << "TTi =\n" << M.TTiO();
    return output;
  }


  M3intO Mesh::TVO() const { return M3intO(TV_,nT_); };
  M3intO Mesh::TTO() const { return M3intO(TT_,nT_); };
  MintO Mesh::VTO() const { return MintO(VT_,nV_); };
  M3intO Mesh::TTiO() const { return M3intO(TTi_,nT_); };
  M3doubleO Mesh::SO() const { return M3doubleO(S_,nV_); };

  std::ostream& operator<<(std::ostream& output, const MintO& MO)
  {
    if (!MO.M_) return output;
    for (int i = 0; i < (int)MO.n_; i++) {
      output << ' ' << std::right << std::setw(4)
	     << MO.M_[i];
    }
    std::cout << std::endl;
    return output;
  }

  std::ostream& operator<<(std::ostream& output, const M3intO& MO)
  {
    if (!MO.M_) return output;
    for (int j = 0; j<3; j++) {
      for (int i = 0; i < (int)MO.n_; i++) {
	output << ' ' << std::right << std::setw(4)
	       << MO.M_[i][j];
      }
      std::cout << std::endl;
    }
    return output;
  }

  std::ostream& operator<<(std::ostream& output, const M3doubleO& MO)
  {
    if (!MO.M_) return output;
    for (int i = 0; i < (int)MO.n_; i++) {
      for (int j = 0; j<3; j++)
	output << ' ' << std::right << std::setw(10) << std::scientific
	       << MO.M_[i][j];
      std::cout << std::endl;
    }
    return output;
  }






  std::ostream& operator<<(std::ostream& output, const Dart& d)
  {
    output << std::right << std::setw(1) << d.t_
	   << std::right << std::setw(3) << d.edir_
	   << std::right << std::setw(2) << d.vi_;
    if ((!d.isnull()) && (d.t_<(int)d.M()->nV())) {
      output << " ("
	     << d.M()->TV()[d.t_][d.vi_]
	     << ","
	     << d.M()->TV()[d.t_][(d.vi_+(3+d.edir_))%3]
	     << ")";
    }
      
    return output;
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




  /*!
    Alg 9.1

    If the point is located within the triangulation domain,
    delta_min and the returned Dart correspond to the triangle edge
    with smallest distance, as measured by inLeftHalfspace.

    If the point is not found, a null Dart is returned.
   */
  Dart Mesh::locatePoint(const Dart& d0,
			 const Point s,
			 double* delta_min) const
  {
    Dart dart;
    if (d0.isnull())
      dart = Dart(*this,0);
    else
      dart = Dart(*this,d0.t(),1,d0.vi());
    Dart dart_start = dart;
    double delta;
    Dart dart_min = Dart();
    while (1) {
      std::cout << dart_start << ' '
		<< dart << ' '
		<< Traits::inLeftHalfspace(dart,s)
		<< std::endl;
      delta = Traits::inLeftHalfspace(dart,s);
      if (dart_min.isnull() || (delta<*delta_min)) {
	dart_min = dart;
	*delta_min = delta;
      }
      if (delta >= -MESH_EPSILON) {
	dart.orbit2();
	if (dart==dart_start)
	  return dart_min;
      } else {
	if (dart.onBoundary())
	  return Dart();
	dart.alpha2();
	dart_start = dart;
	dart_start.alpha0();
	dart.alpha1();
	dart_min = dart_start;
	*delta_min = -delta;
      }
    }

    return Dart();
  }

  /*!
    Alg 9.1 modified to locate a pre-existing vertex.

    If the vertex is not found, a null Dart is returned.
   */
  Dart Mesh::locateVertex(const Dart& d0,
			  const int v) const
  {
    if (use_VT_) {
      int t;
      t = VT_[v];
      if (t<0) /* Vertex not connected to any triangles. */
	return Dart();
      if (TV_[t][0] == v)
	return Dart(*this,t,1,0);
      if (TV_[t][1] == v)
	return Dart(*this,t,1,1);
      if (TV_[t][2] == v)
	return Dart(*this,t,1,2);
      return Dart(); /* ERROR: Inconsistent data structures! */
    }

    int i;
    Dart dart;
    if (d0.isnull())
      dart = Dart(*this,0);
    else
      dart = Dart(*this,d0.t(),1,d0.vi());
    Dart dart_start = dart;
    double delta;
    Dart dart_min = Dart();
    double* s = &(S_[v][0]);
    double delta_min = 0.0;
    while (1) {
      std::cout << dart_start << ' '
		<< dart << ' '
		<< Traits::inLeftHalfspace(dart,s)
		<< std::endl;
      for (i=0;i<3;i++) {
	if (TV_[dart.t()][dart.vi()] == v)
	  return dart;
	dart.orbit2();
      }

      delta = Traits::inLeftHalfspace(dart,s);
      if (dart_min.isnull() || (delta<delta_min)) {
	dart_min = dart;
	delta_min = delta;
      }
      if (delta >= -MESH_EPSILON) {
	dart.orbit2();
	if (dart==dart_start) {
	  for (i=0;i<3;i++) {
	    if (TV_[dart.t()][dart.vi()] == v)
	      return dart;
	    dart.orbit2();
	  }
	  return Dart(); /* ERROR: Point located, but not the vertex itself. */
	}
      } else {
	if (dart.onBoundary())
	  return Dart();
	dart.alpha2();
	dart_start = dart;
	dart_start.alpha0();
	dart.alpha1();
	dart_min = dart_start;
	delta_min = -delta;
      }
    }

    return Dart();
  }



  /*! Alg 4.3 */
  bool MeshConstructor::recSwapDelaunay(const Dart& d0)
  {
    Dart d1, d2;

    std::cout << "Trying to swap " << d0 << std::endl;

    if (d0.isnull() or d0.onBoundary())
      return true; /* OK. Not allowed to swap. */
    if (isSegmentDart(d0))
      return true ; /* OK. Not allowed to swap. */
    if (Traits::circumcircleOK(d0))
      return true; /* OK. Need not swap. */

    std::cout << "Swap " << d0 << std::endl;
    xtmpl_press_ret("swap edge");

    /* Get opposing darts. */
    d1 = d0;
    d1.alpha1();
    if (d1.onBoundary()) d1 = Dart(); else d1.alpha2();
    d2 = d0;
    d2.orbit2rev().alpha1(); 
    if (d2.onBoundary()) d2 = Dart(); else d2.alpha2();
    
    std::cout << "TVpre  = " << std::endl << M_->TVO();
    swapEdge(d0);
    std::cout << "TVpost = " << std::endl << M_->TVO();
    std::cout << "TTpost = " << std::endl << M_->TTO();
    M_->redrawX11();
    xtmpl_press_ret("edge swapped, next recSwapDelaunay");

    if (!d1.isnull()) recSwapDelaunay(d1);
    M_->redrawX11();
    xtmpl_press_ret("After d1-recSwapDelaunay");
    if (!d2.isnull()) recSwapDelaunay(d2);
    M_->redrawX11();
    xtmpl_press_ret("After d2-recSwapDelaunay");
    return true;
  }


  /*! Alg 9.3 */
  Dart MeshConstructor::splitTriangleDelaunay(const Dart& td, int v)
  {
    Dart d, d0, d1, d2;

    if (td.isnull()) { return Dart(); }; /* ERROR */
    /* Get opposing darts. */
    d = td;
    if (d.onBoundary()) d0 = Dart(); else {d0 = d; d0.orbit1();} 
    d.orbit2();
    if (d.onBoundary()) d1 = Dart(); else {d1 = d; d1.orbit1();} 
    d.orbit2();
    if (d.onBoundary()) d2 = Dart(); else {d2 = d; d2.orbit1();} 

    std::cout << "TV = " << std::endl << M_->TVO();
    std::cout << "Split triangle with vertex " << v << std::endl;
    d = splitTriangle(td,v);
    
    M_->redrawX11();
    xtmpl_press_ret("triangle split, next recSwapDelaunay");

    if (!d0.isnull()) recSwapDelaunay(d0);
    M_->redrawX11();
    xtmpl_press_ret("After d0-recSwapDelaunay");
    if (!d1.isnull()) recSwapDelaunay(d1);
    M_->redrawX11();
    xtmpl_press_ret("After d1-recSwapDelaunay");
    if (!d2.isnull()) recSwapDelaunay(d2);
    M_->redrawX11();
    xtmpl_press_ret("After d2-recSwapDelaunay");

    std::cout << "TV = " << std::endl << M_->TVO();
    
    return d;
  }

  /*! Modified Alg 9.3 */
  Dart MeshConstructor::splitEdgeDelaunay(const Dart& ed, int v)
  {
    Dart d, d0, d1, d2, d3;

    if (ed.isnull()) { return Dart(); }; /* ERROR */
    /* Get opposing darts. */
    d = ed;
    d.orbit2();
    if (d.onBoundary()) d0 = Dart(); else {d0 = d; d0.orbit1();} 
    d.orbit2();
    if (d.onBoundary()) d1 = Dart(); else {d1 = d; d1.orbit1();} 
    d = ed;
    if (d.onBoundary()) {
      d2 = Dart();
      d3 = Dart();
    } else {
      d.orbit0rev();
      if (d.onBoundary()) d2 = Dart(); else {d2 = d; d2.orbit1();} 
      d.orbit2();
      if (d.onBoundary()) d3 = Dart(); else {d3 = d; d3.orbit1();} 
    }

    std::cout << "TV = " << std::endl << M_->TVO();
    std::cout << "Split edge with vertex " << v << std::endl;
    d = splitEdge(ed,v);
    
    if (!d0.isnull()) recSwapDelaunay(d0);
    if (!d1.isnull()) recSwapDelaunay(d1);
    if (!d2.isnull()) recSwapDelaunay(d2);
    if (!d3.isnull()) recSwapDelaunay(d3);

    std::cout << "TV = " << std::endl << M_->TVO();
    
    return d;
  }

  /*! Alg 9.3 */
  bool MeshConstructor::insertNode(int v, const Dart& ed)
  {
    Dart td;
    double delta;

    std::cout << "Locating node " << v << std::endl;
    td = M_->locatePoint(ed,M_->S()[v],&delta);
    if (td.isnull()) { return false; }; /* ERROR, not found! */
    std::cout << "Closest dart " << td
	      << ' ' << delta << std::endl;

    if (delta>10*MESH_EPSILON) { /* Split triangle */
      splitTriangleDelaunay(td,v);
    } else { /* Split edge */
      splitEdgeDelaunay(td,v);
    }

    return true;
  }

  bool MeshConstructor::DT(const vertex_input_type& v_set)
  {
    if (is_pruned_) 
      return false; /* ERROR, cannot safely insert nodes into a pruned
		       triangulation. Call insertNode directly if known to
		       be visible/reachable from a given edge.  */

    if (state_ < State_CHT)
      return false; /* TODO: Add convex enclosure? */

    if (state_ < State_DT)
      if (!prepareDT()) /* Make sure we have a DT. */
	return false;

    int v;
    vertex_input_type::const_iterator v_iter;
    Dart td, d, d0, d1, d2;

    for (v_iter = v_set.begin(); v_iter != v_set.end(); v_iter++) {
      v = *v_iter;
      insertNode(v,Dart(*M_,0)); /* TODO: More clever starting edge? */
      std::cout << M_->VTO();
      
      M_->redrawX11();
      xtmpl_press_ret("next node insert");
    }

    state_ = State_DT;
    return true;
  }


  bool MeshConstructor::prepareDT()
  {
    if (state_<State_DT) {
      /* We need to build a DT first. */
      triangle_input_type t_set;
      for (int t=0;t<(int)M_->nT();t++)
	t_set.push_back(t);
      if (LOP(t_set))
	state_ = State_DT;
    }
    return (state_>=State_DT) && (!is_pruned_);
  }


  bool MeshConstructor::prepareCDT()
  {
    if (!prepareDT()) return false; /* Make sure we have a DT. */
    if (state_>=State_CDT)
      return true; /* Nothing to do. Data structures already active. */

    const int* tt;
    int vi;
    Dart d;
    for (int t=0;t<(int)M_->nT();t++) {
      tt = M_->TT()[t];
      for (vi=0;vi<3;vi++)
	if (tt[vi]<0) {
	  d = Dart(*M_,t,1,(vi+1)%3);
	  boundary_.insert(d);
	}
    }

    state_ = State_CDT;
    return true;
  }

  bool MeshConstructor::prepareRCDT(double skinny_limit, double big_limit)
  {
    if (!prepareCDT()) return false; /* Make sure we have a CDT. */

    skinny_limit_ = skinny_limit;
    big_limit_ = big_limit;

    skinny_ = DartQualitySet(skinny_limit_);
    big_ = DartQualitySet(big_limit_);

    double quality;
    for (int t=0;t<(int)M_->nT();t++) {
      quality = M_->skinnyQuality(t);
      if (quality>skinny_limit_)
	skinny_.insert(Dart(*M_,t),quality);
      quality = M_->bigQuality(t);
      if (quality>big_limit_)
	big_.insert(Dart(*M_,t),quality);
    }

    state_ = State_RCDT;
    return true;
  }



  bool MeshConstructor::CDTBoundary(const constraint_input_type& constr)
  {
    if (!prepareCDT()) return false;

    constr_boundary_ = constraint_list_type(constr.begin(),constr.end());
    
    return buildCDT();
  };

  bool MeshConstructor::CDTInterior(const constraint_input_type& constr)
  {
    if (!prepareCDT()) return false;

    constr_interior_ = constraint_list_type(constr.begin(),constr.end());
    
    return buildCDT();
  };






  bool MeshConstructor::LOP(const triangle_input_type& t_set)
  {
    /* TODO: Implement. */
    NOT_IMPLEMENTED;

    return true;
  }


  bool MeshConstructor::buildCDT()
  {
    if (!prepareCDT()) return false;

    /* TODO: Implement. */
    NOT_IMPLEMENTED;

    for (constraint_list_type::iterator ci = constr_boundary_.begin();
	 ci != constr_boundary_.end(); ci++) {
      if (true)
	ci = constr_boundary_.erase(ci);
    }
    for (constraint_list_type::iterator ci = constr_interior_.begin();
	 ci != constr_interior_.end(); ci++) {
      if (true)
	ci = constr_interior_.erase(ci);
    }

    return (constr_boundary_.empty() && constr_interior_.empty());
  };

  bool MeshConstructor::buildRCDT()
  {
    if (state_<State_RCDT)
      return false; /* ERROR: RCDT not initialised. */

    /* TODO: Implement. */
    NOT_IMPLEMENTED;

    return true;
  };

  bool MeshConstructor::RCDT(double skinny_limit, double big_limit)
  {
    if (!prepareRCDT(skinny_limit,big_limit)) return false;
    return buildRCDT();
  };


  bool MeshConstructor::PruneExterior()
  {
    if (state_ < State_CDT) {
      /* Since there are no constraints at this state, no exterior
	 needs to be pruned, but add the boundary to the constraint
	 set anyway, to get to state State_CDT. */
      prepareCDT();
      is_pruned_ = true;
      return true;
    }
    is_pruned_ = true;
    
    /* TODO: Implement. */
    NOT_IMPLEMENTED;

    return true;
  };





  Dart MeshConstructor::swapEdge(const Dart& d)
  {
    if (state_ < State_CDT) {
      return M_->swapEdge(d);
    }

    /* TODO: implement. */
    NOT_IMPLEMENTED;
    return Dart(d);
  }

  Dart MeshConstructor::splitEdge(const Dart& d, int v)
  {
    if (state_ < State_CDT) {
      return M_->splitEdge(d,v);
    }

    /* TODO: implement. */
    NOT_IMPLEMENTED;

    return Dart(d);
  }

  Dart MeshConstructor::splitTriangle(const Dart& d, int v)
  {
    if (state_ < State_CDT) {
      return M_->splitTriangle(d,v);
    }

    /* TODO: implement. */
    NOT_IMPLEMENTED;
    return Dart(d);
  }




  bool MeshConstructor::isSegmentDart(const Dart& d) const
  {
    if (state_<State_CDT) /* No segments */
      return false;

    return (boundary_.found(d) || interior_.found(d));
  }







  bool DartQualitySet::found(const Dart& d) const
  {
    return (darts_.find(d) != darts_.end());
  }

  bool DartQualitySet::found_quality(const Dart& d) const
  {
    map_type::const_iterator i = darts_.find(d);
    if (i == darts_.end())
      return false;
    return (darts_quality_.find(MCdv(i->first,i->second)) !=
	    darts_quality_.end());
  }

  const double DartQualitySet::quality(const Dart& d) const
  {
    if (empty())
      return 0.0;
    return darts_.find(d)->second;
  }

  Dart DartQualitySet::quality_dart() const
  {
    if (empty_quality())
      return Dart();
    return darts_quality_.begin()->d_;
  }

  void DartQualitySet::insert(const Dart& d, double quality)
  {
    darts_.insert(map_key_type(d,quality));
    if (quality>=quality_limit_)
      darts_quality_.insert(MCdv(d,quality));
  }

  void DartQualitySet::erase(const Dart& d)
  {
    double quality;
    map_type::iterator i = darts_.find(d);
    if (i != darts_.end()) {
      quality = i->second;
      darts_.erase(i);
      set_type::iterator j = darts_quality_.find(MCdv(d,quality));
      if (j != darts_quality_.end())
	darts_quality_.erase(j);
    }
  }


} /* namespace fmesh */
