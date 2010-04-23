#include <cstddef>
#include <cstring>
#include <set>
#include <map>

#include "mesh.h"

namespace fmesh {

  


  /*
    T-E+V=2
    
    closed 2-manifold triangulation:
    E = T*3/2
    T = 2*V-4

    simply connected 2-manifold triangulation:
    T <= 2*V-5

  */

  Mesh::Mesh(size_t V_capacity,
	     bool use_TTi) : Vcap_(V_capacity),
			     Tcap_(2*V_capacity),
			     nV_(0), nT_(0), use_TTi_(use_TTi)
  {
    if (Vcap_ > 0) {
      TV_ = new int[Vcap_][3];
      TT_ = new int[Tcap_][3];
      if (use_TTi_)
	TTi_ = new int[Tcap_][3];
      else
	TTi_ = NULL;
      S_ = new double[Vcap_][3];
    } else {
      TV_ = NULL;
      TT_ = NULL;
      TTi_ = NULL;
      S_ = NULL;
    }
  };

  Mesh::~Mesh()
  {
    if (TV_) delete[] TV_;
    if (TT_) delete[] TT_;
    if (TTi_) delete[] TTi_;
    if (S_) delete[] S_;
  }
  Mesh& Mesh::clear()
  {
    Vcap_ = 0;
    Tcap_ = 0;
    nV_ = 0;
    nT_ = 0;
    use_TTi_ = false;
    if (TV_) { delete[] TV_; TV_ = NULL; }
    if (TT_) { delete[] TT_; TT_ = NULL; }
    if (TTi_) { delete[] TTi_; TTi_ = NULL; }
    if (S_) { delete[] S_; S_ = NULL; }
    return *this;
  }

  Mesh& Mesh::check_capacity(int nVc, int nTc)
  {
    if ((nVc <= Vcap_) && (nTc <= Tcap_))
      return *this;
    while ((nVc > Vcap_) || (nTc > Tcap_)) {
      Vcap_ = Vcap_+Mesh_V_capacity_step_size;
      Tcap_ = 2*Vcap_;
    }

    int (*TV)[3] = new int[Vcap_][3];
    int (*TT)[3] = new int[Tcap_][3];
    int (*TTi)[3] = NULL;
    if (use_TTi_) TTi = new int[Tcap_][3];
    double (*S)[3] = new double[Vcap_][3];

    if (TV_) {
      if (TV_) memcpy(TV,TV_,sizeof(int)*nV_*3);
      if (TT_) memcpy(TT,TT_,sizeof(int)*nT_*3);
      if (TTi_) memcpy(TTi,TTi_,sizeof(int)*nT_*3);
      if (S_) memcpy(S,S_,sizeof(double)*nV_*3);
      if (TV_) delete[] TV_;
      if (TT_) delete[] TT_;
      if (TTi_) delete[] TTi_;
      if (S_) delete[] S_;
    }

    TV_ = TV;
    TT_ = TT;
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
    for (t=0; t<nT_; t++) {
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
    for (t=0; t<nT_; t++) {
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
    for (t=0; t<nT_; t++) {
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
    return *this;
  }
  
  Mesh& Mesh::TV_append(const int (*TV)[3], int nT)
  {
    check_capacity(0,nT_+nT);
    memcpy(TV_+nT_,TV,sizeof(int)*nT*3);
    nT_ += nT;
    rebuildTT();
    rebuildTTi();
    return *this;
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
    int t, vi;
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

    return Dart(*this,t0,1,1);
  }
  
  Dart Mesh::splitEdge(const Dart& d, int v)
  {
  }

  Dart Mesh::splitTriangle(const Dart& d, int v)
  /**
   *   2          2
   *  |  \       |\ \
   *  |   \      | \ \
   *  |    1 --> | vd-1
   *  |   /      | / /
   *  | d/       |/ /
   *   0          0
   *
   *  Dart 0-1 --> 3-1
   * 
   */
  {
    Dart dhelper = d;
    int t, vi, i;
    int v_list[3];
    int t0, t1, t2;
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

    return Dart(*this,t0,1,0);
  }





  std::ostream& operator<<(std::ostream& output, const Mesh& M)
  {
    output << "S =\n" << M.SO();
    output << "TV =\n" << M.TVO();
    output << "TT =\n" << M.TTO();
    if (M.useTTi())
      output << "TTi =\n" << M.TTiO();
    return output;
  }


  M3intO Mesh::TVO() const { return M3intO(TV_,nT_); };
  M3intO Mesh::TTO() const { return M3intO(TT_,nT_); };
  M3intO Mesh::TTiO() const { return M3intO(TTi_,nT_); };
  M3doubleO Mesh::SO() const { return M3doubleO(S_,nV_); };

  std::ostream& operator<<(std::ostream& output, const M3intO& MO)
  {
    if (!MO.M_) return output;
    for (int i = 0; i < MO.n_; i++) {
      for (int j = 0; j<3; j++)
	output << ' ' << std::right << std::setw(4)
	       << MO.M_[i][j];
      std::cout << '\n';
    }
    return output;
  }
  std::ostream& operator<<(std::ostream& output, const M3doubleO& MO)
  {
    if (!MO.M_) return output;
    for (int i = 0; i < MO.n_; i++) {
      for (int j = 0; j<3; j++)
	output << ' ' << std::right << std::setw(10) << std::scientific
	       << MO.M_[i][j];
      std::cout << '\n';
    }
    return output;
  }






  std::ostream& operator<<(std::ostream& output, const Dart& d)
  {
    output << std::right << std::setw(5) << d.t_
	   << std::right << std::setw(3) << d.edir_
	   << std::right << std::setw(2) << d.vi_;
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




  Dart Mesh::locatePoint(const Dart& d, const Point s) const
  {
    return Dart();
  }


  bool circumcircleTest(const Dart& d)
  {
    /* TODO: implement circle-test using predicates::incircle */
    return false;
  }

  /*! Alg 9.4 */
  void MeshConstructor::recSwapDelaunay(const Dart& d0)
  {
    Dart d1, d2;

    if (circumcircleTest(d0))
      return;

    /* Get opposing darts. */
    d1 = d0;
    d1.alpha1();
    if (!d1.onBoundary()) d1.alpha2();
    d2 = d0;
    d2.orbit2rev().alpha1(); 
    if (d2.onBoundary()) d2.alpha2();
    
    swapEdge(d0);

    if (!d1.onBoundary()) recSwapDelaunay(d1);
    if (!d2.onBoundary()) recSwapDelaunay(d2);
  }


  /*! Alg 9.3 */
  void MeshConstructor::insertNode(int v)
  {
    Dart td, d, d0, d1, d2;

    td = M_->locatePoint(Dart(*M_,0),M_->S()[v]);
    if (td.isnull()) { return; }; /* ERROR, not found! */
    
    /* Get opposing darts. */
    d = td;
    if (d.onBoundary()) d0 = Dart(); else {d0 = d; d0.orbit1();} 
    d.orbit2(); 
    if (d.onBoundary()) d1 = Dart(); else {d1 = d; d1.orbit1();} 
    d.orbit2(); 
    if (d.onBoundary()) d2 = Dart(); else {d2 = d; d1.orbit1();} 
    
    td = splitTriangle(td,v);
    
    recSwapDelaunay(d0);
    recSwapDelaunay(d1);
    recSwapDelaunay(d2);
  }

  void MeshConstructor::DT(const std::vector<int> v_set)
  {
    if (state_ > State_DT)
      return;

    int v;
    std::vector<int>::const_iterator v_iter;
    Dart td, d, d0, d1, d2;

    for (v_iter = v_set.begin(); v_iter != v_set.end(); v_iter++) {
      v = *v_iter;
      insertNode(v);
    }
  }


  Dart MeshConstructor::swapEdge(const Dart& d)
  {
    if (state_ < State_CDT_prepared) {
      return M_->swapEdge(d);
    }

    /* TODO: implement. */
  }

  Dart MeshConstructor::splitEdge(const Dart& d, int v)
  {
    if (state_ < State_CDT_prepared) {
      return M_->splitEdge(d,v);
    }

    /* TODO: implement. */
  }

  Dart MeshConstructor::splitTriangle(const Dart& d, int v)
  {
    if (state_ < State_CDT_prepared) {
      return M_->splitTriangle(d,v);
    }

    /* TODO: implement. */
  }





  /*

int point_in_dart_lefthalfplane(const trimesh_t* M,
				point_t* p,
				const dart_t* dart)
{
  return ( ((*p)[0]-(*p)[2])*
	   (trimesh_S((*M),dart->f,dart->vi,1)-
	    trimesh_S((*M),dart->f,dart->vi,2)) -
	   ((*p)[1]-(*p)[2])*
	   (trimesh_S((*M),dart->f,dart->vi,0)-
	    trimesh_S((*M),dart->f,dart->vi,2)) ) >= 0;
}


int locate_triangle(const trimesh_t* M,
		    point_t* p,
		    const dart_t* dart_init)
//
//  Algorithm from Hjelle & Daehlen, p. 209
//
{
  int f;
  dart_t dart_start;
  dart_t dart;
  dart_copy(dart_init,&dart);
  dart.edir = 1;
  dart_copy(&dart,&dart_start);
  while (1) {
    if (point_in_dart_lefthalfplane(M,p,&dart)) {
      dart_alpha1(M,dart_alpha0(M,&dart));
      if (dart_equal(M,&dart,&dart_start)) {
	return dart.f;
      }
    } else {
      dart_alpha2(M,dart_copy(&dart,&dart_start));
      if (dart_equal(M,&dart,&dart_start)) {
	return -1;
      }
      dart_alpha1(M,dart_copy(&dart_start,&dart));
      dart_alpha0(M,&dart_start);
    }
  }
}

  */

} /* namespace fmesh */
