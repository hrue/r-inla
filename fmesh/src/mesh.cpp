#include <cstddef>
#include <cstring>
#include <set>
#include <map>

#include "mesh.h"

#define Mesh_V_capacity_step_size 128

namespace fmesh {

  

  /*  
  bool isBoundaryEdge(const Dart& d)
  {
    Dart d2 = d;
    return (d2.alpha2() == d);
  }
  */



  /*
    T-E+V=2
    
    closed mesh in R3:
    E = T*3/2
    T = 2*V-4

    simple mesh in R2:
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



  Mesh& Mesh::S_set(double (*S)[3], int nV)
  {
    nV_ = 0; /* Avoid possible unnecessary copy. */
    check_capacity(nV,0);
    nV_ = nV;
    memcpy(S_,S,sizeof(double)*nV_*3);
    return *this;
  }
  
  Mesh& Mesh::TV_set(int (*TV)[3], int nT)
  {
    nT_ = 0; /* Avoid possible unnecessary copy. */
    check_capacity(0,nT);
    nT_ = nT;
    memcpy(TV_,TV,sizeof(int)*nT_*3);
    rebuildTT();
    rebuildTTi();
    return *this;
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
    /* "alpha1(); alpha0();" would be less efficient. */
    vi_ = (vi_+(3+edir_))%3;
    return *this;
  }

  Dart& Dart::orbit0cw()
  {
    int t = t_;
    alpha2();
    if (t != t_) alpha1(); /* Do only if not at boundary. */
    return *this;
  }

  Dart& Dart::orbit1cw() /* Equivalent to orbit1() */
  {
    orbit1();
    return *this;
  }

  Dart& Dart::orbit2cw()
  {
    /* "alpha0(); alpha1();" would be less efficient. */
    vi_ = (vi_+(3-edir_))%3;
    return *this;
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
