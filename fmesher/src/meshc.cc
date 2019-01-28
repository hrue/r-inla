#include <cstddef>
#include <cstring>
#include <set>
#include <map>
#include <sstream>
#include <cmath>
#include <vector>

#include "predicates.hh"
#include "meshc.hh"

#define WHEREAMI __FILE__ << "(" << __LINE__ << ")\t"

#define MESHC_LOG_(msg) cout << WHEREAMI << msg;
#ifdef DEBUG
#define MESHC_LOG(msg) MESHC_LOG_(msg)
#else
#define MESHC_LOG(msg)
#endif

using std::cout;
using std::endl;

namespace fmesh {

  bool MCQ::found(const Dart& d) const
  {
    return (darts_.find(d) != darts_.end());
  }

  bool MCQ::foundQ(const Dart& d) const
  {
    map_type::const_iterator i = darts_.find(d);
    if (i == darts_.end())
      return false;
    return (darts_quality_.find(MCQdv(i->first,i->second)) !=
	    darts_quality_.end());
  }

  double MCQ::quality(const Dart& d) const
  {
    if (empty())
      return 0.0;
    return darts_.find(d)->second;
  }

  Dart MCQ::quality() const
  {
    if (emptyQ())
      return Dart();
    return darts_quality_.rbegin()->d_;
  }

  void MCQ::insert(const Dart& d)
  {
    double quality_ = calcQ(d);
    if (quality_>0.0) {
      darts_.insert(map_key_type(d,quality_));
      darts_quality_.insert(MCQdv(d,quality_));
    } else if (!only_quality_)
      darts_.insert(map_key_type(d,quality_));
  }

  void MCQ::erase(const Dart& d)
  {
    double quality_;
    map_type::iterator i = darts_.find(d);
    if (i != darts_.end()) {
      quality_ = i->second;
      darts_.erase(i);
      set_type::iterator j = darts_quality_.find(MCQdv(d,quality_));
      if (j != darts_quality_.end())
	darts_quality_.erase(j);
    }
  }







  MCQtri::MCQtri(MeshC* MC, bool only_quality,
		 double quality_limit,
		 const double* quality_limits,
		 size_t nQL)
    : MCQ(MC, only_quality),
      quality_limit_(quality_limit),
      quality_limits_(NULL),
      quality_limits_cap_(0)
  {
    setQ(quality_limit, quality_limits, nQL);
  }

  void MCQtri::setQ(double quality_limit,
		    const double* quality_limits,
		    size_t nQL)
  {
    quality_limit_ = quality_limit;
    if (quality_limits) {
      if (quality_limits_cap_ < MC_->M_->Vcap()) {
	if (quality_limits_) {
	  delete[] quality_limits_;
	}
	quality_limits_cap_ = MC_->M_->Vcap();
	quality_limits_ = new double[quality_limits_cap_];
      }
      if (nQL>=MC_->M_->nV())
	memcpy(quality_limits_,quality_limits,sizeof(double)*MC_->M_->nV());
      else {
	memcpy(quality_limits_,quality_limits,sizeof(double)*nQL);
	for (int v=nQL; v<(int)MC_->M_->nV(); v++)
	  quality_limits_[v] = quality_limit_;
      }
    } else {
      if (quality_limits_) {
	delete[] quality_limits_;
	quality_limits_ = NULL;
      }
    }
  }

  void MCQtri::setQv(int v, double quality_limit)
  {
    if (quality_limits_cap_ < MC_->M_->Vcap()) {
      size_t old_quality_limits_cap_ = quality_limits_cap_;
      quality_limits_cap_ = MC_->M_->Vcap();
      double* ql = new double[quality_limits_cap_];
      if (quality_limits_) {
	memcpy(ql,quality_limits_,sizeof(double)*old_quality_limits_cap_);
	delete[] quality_limits_;
      }
      quality_limits_ = ql;
    }
    quality_limits_[v] = quality_limit;
  }

  double MCQtri::getQ(int t) const
  {
    if (quality_limits_ == NULL)
      return quality_limit_;
    else {
      double lim = quality_limits_[MC_->M_->TV(t)[0]];
      if (lim > quality_limits_[MC_->M_->TV(t)[1]])
	lim = quality_limits_[MC_->M_->TV(t)[1]];
      if (lim > quality_limits_[MC_->M_->TV(t)[2]])
	lim = quality_limits_[MC_->M_->TV(t)[2]];
      return lim;
    }
  }

  double MCQtri::calcQ(const Dart& d) const {
    return (calcQtri(d) - getQ(d.t()));
  }

  double MCQskinny::calcQtri(const Dart& d) const
  {
    return MC_->skinnyQuality(d.t());
  }
  
  double MCQbig::calcQtri(const Dart& d) const
  {
    return MC_->bigQuality(d.t());
  }
  
  double MCQsegm::calcQ(const Dart& d) const
  {
    double quality_ = MC_->encroachedQuality(d);
    Dart dh(d);
    dh.orbit1();
    if (d.t() != dh.t()) {
      double quality1_ = MC_->encroachedQuality(dh);
      if (quality1_>quality_)
	quality_ = quality1_;
    }
    return (quality_-encroached_limit_);
  }

  bool MCQsegm::segm(const Dart& d) const
  {
    if (found(d)) return true;
    Dart dh(d);
    dh.orbit1();
    return ((dh.t() != d.t()) && found(dh));
  }

  void MCQsegm::update(const Dart& d)
  {
    if (found(d)){
      insert(d,erase(d));
    }
    Dart dh(d);
    dh.orbit1();
    if ((dh.t() != d.t()) && found(dh)) {
      insert(dh,erase(dh));
    }
  }

  MCQsegm::meta_type MCQsegm::meta(const Dart& d) const
  {
    if (MCQ::empty())
      return meta_type();
    return meta_.find(d)->second;
  }
  
  void MCQsegm::clear()
  {
    meta_.clear();
    MCQ::clear();
  }
  
  void MCQsegm::insert(const Dart& d, const meta_type& meta)
  {
    MCQ::insert(d);
    meta_.insert(meta_map_key_type(d,meta));
  }

  MCQsegm::meta_type MCQsegm::erase(const Dart& d)
  {
    meta_type meta = meta_type();
    meta_map_type::iterator i = meta_.find(d);
    if (i != meta_.end()) {
      meta = i->second;
      meta_.erase(i);
    }
    MCQ::erase(d);
    return meta;
  }



  double MCQswapable::calcQ(const Dart& d) const
  {
    return (d.isSwapable() ? 1.0 : -1.0);
  }

  double MCQswapableD::calcQ(const Dart& d) const
  {
    return (d.isSwapableD() ? 1.0 : -1.0);
  }

  bool MCQswapable::found(const Dart& d) const
  {
    if (MCQ::found(d)) return true;
    Dart dh(d);
    dh.orbit1();
    if (dh.t() != d.t())
      if (MCQ::found(dh)) return true;
    return false;
  }

  bool MCQswapable::foundQ(const Dart& d) const
  {
    if (MCQ::foundQ(d)) return true;
    Dart dh(d);
    dh.orbit1();
    if (dh.t() != d.t())
      if (MCQ::foundQ(dh)) return true;
    return false;
  }

  double MCQswapable::quality(const Dart& d) const
  {
    if (MCQ::foundQ(d)) return MCQ::quality(d);
    Dart dh(d);
    dh.orbit1();
    if (dh.t() != d.t())
      return MCQ::quality(d);
    return 0.0;
  }

  void MCQswapable::insert(const Dart& d)
  {
    if (found(d)) return;  /* Don't add duplicates. */
      MCQ::insert(d);
  }

  void MCQswapable::erase(const Dart& d)
  {
    MCQ::erase(d);
    Dart dh(d);
    dh.orbit1();
    if (dh.t() != d.t()) {
      MCQ::erase(dh);
    }
  }

  bool MCQswapable::swapable(const Dart& d) const
  {
    if (foundQ(d)) return true;
    Dart dh(d);
    dh.orbit1();
    /* Boundary edges should never be in the set, but check
       for safety, and then check if found. */
    return ((dh.t() != d.t()) && foundQ(dh));
  }



  double MeshC::encroachedQuality(const Dart& d) const
  /* >0 --> encroached */
  {
    int t(d.t());
    if ((t<0) || (t>=(int)M_->nT())) return -1.0;

    Dart dh(d);
    dh.orbit2rev();
    
    double encr = M_->edgeEncroached(d,M_->S(dh.v()));

    dh.orbit2rev();
    MESHC_LOG("encroachedQ("
	      << d.v() << "," << dh.v()
	      << ") = " << encr << endl);

    return encr;
  }

  double MeshC::skinnyQuality(int t) const
  {
    if ((t<0) || (t>=(int)M_->nT())) return 0.0;

    double skinny;
    Point len;
    int argmin = M_->triangleEdgeLengthsArgMin(t,len);
    bool ok = true;

    if (state_>=State_CDT) {
      Dart dh(*M_,t,1,(argmin+2)%3);
      ok = ((!boundary_.segm(dh)) && (!interior_.segm(dh)));
      if (!ok) {
	dh.orbit2();
	ok = ((!boundary_.segm(dh)) && (!interior_.segm(dh)));
      }
    } else
      ok = !((M_->TT(t)[(argmin+1)%3] < 0) &&
	     (M_->TT(t)[(argmin+2)%3] < 0));

    if (ok)
      skinny = (M_->triangleCircumcircleRadius(t) / 
		len[argmin]);
    else
      skinny = 0.0;

    //    MESHC_LOG("skinnyQ(" << t << ") = " << skinny << endl);

    return skinny;
  }

  double MeshC::bigQuality(int t) const
  {
    return M_->triangleLongestEdge(t);
    //return M_->triangleArea(t);
  }




  typedef std::pair<int,Dart> intDartPairT;
  typedef std::multimap<int,Dart> intDartMapT;

  /* Remove dart from sets. */
  Dart erase_dart_from_set(intDartMapT::iterator i,
			   intDartMapT& map_v0_d)
  {
    std::pair<intDartMapT::iterator,intDartMapT::iterator> candidates;
    Dart d = i->second;
    map_v0_d.erase(i);
    /* Old code for handling the reverse mapping as well */
    /*
    candidates = map_v1_d.equal_range(d.vo());
    for (intDartMapT::iterator ic=candidates.first;
	 ic!=candidates.second;
	 ++ic) {
      if (ic->second == d) {
	map_v1_d.erase(ic);
	break;
      }
    }
    */
    return d;
  }


  /* Find "next" dart in v0-set. */
  intDartMapT::iterator find_next_dart_in_set(Dart d,
					      intDartMapT& map_v0_d)
  {
    intDartMapT::iterator candidate;
    candidate = map_v0_d.find(d.vo());
    if (candidate != map_v0_d.end())
      return candidate;
    return map_v0_d.end();
  }


  int extract_segments(const MCQsegm& seg,
		       Matrix<int>* segm,
		       Matrix<int>* segmgrp)
  {
    if (segm==NULL) {
      return seg.count();
    }

    intDartMapT map_v0_d;
    for (MCQsegm::const_iterator ci = seg.begin();
	 ci != seg.end();
	 ci++) {
      map_v0_d.insert(intDartPairT(ci->first.v(),ci->first));
    }

    int segm_initial = segm->rows();
    for (intDartMapT::iterator i = map_v0_d.begin();
	 i != map_v0_d.end();
	 i = map_v0_d.begin()) {
      Dart d;
      for (;
	   i != map_v0_d.end();
	   i = find_next_dart_in_set(d,map_v0_d)) {
	d = erase_dart_from_set(i,map_v0_d);
	int segm_i = segm->rows();
	(*segm)(segm_i,0) = d.v();
	(*segm)(segm_i,1) = d.vo();
	if (segmgrp) {
	  (*segmgrp)(segm_i,0) = seg.meta(d);
	}
      }
    }

    return segm->rows()-segm_initial;
  }

  int MeshC::segments(bool boundary,
		      Matrix<int>* segm,
		      Matrix<int>* segmgrp) const
  {
    if (boundary)
      return extract_segments(boundary_,segm,segmgrp);
    else
      return extract_segments(interior_,segm,segmgrp);
  }






  /*! Alg 4.3 */
  bool MeshC::recSwapDelaunay(const Dart& d0)
  {
    Dart d1, d2;

    MESHC_LOG("Trying to swap " << d0 << endl);

    if (d0.isnull() or d0.onBoundary()) {
      MESHC_LOG("Not allowed to swap, boundary" << endl);
      return true; /* OK. Not allowed to swap. */
    }
    if (isSegment(d0)) {
      MESHC_LOG("Not allowed to swap, segment" << endl);
      return true ; /* OK. Not allowed to swap. */
    }
    if (d0.circumcircleOK()) {
      MESHC_LOG("No need to swap, circumcircle OK" << endl);
      return true; /* OK. Need not swap. */
    }

    MESHC_LOG("Swap " << d0 << endl);

    /* Get opposing darts. */
    d1 = d0;
    d1.alpha1();
    if (d1.onBoundary()) d1 = Dart(); else d1.alpha2();
    d2 = d0;
    d2.orbit2rev().alpha1(); 
    if (d2.onBoundary()) d2 = Dart(); else d2.alpha2();
    
    //    MESHC_LOG("TVpre  = " << endl << M_->TVO());
    swapEdge(d0);
    //    MESHC_LOG("TVpost = " << endl << M_->TVO());
    //    MESHC_LOG("TTpost = " << endl << M_->TTO());

    if (!d1.isnull()) recSwapDelaunay(d1);
    if (!d2.isnull()) recSwapDelaunay(d2);
    return true;
  }


  /*! Alg 9.3 */
  Dart MeshC::splitTriangleDelaunay(const Dart& td, int v)
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

    //    MESHC_LOG("TV = " << endl << M_->TVO());
    MESHC_LOG("Split triangle " << td << " with vertex " << v << endl);
    d = splitTriangle(td,v);
    
    if (!d0.isnull()) recSwapDelaunay(d0);
    if (!d1.isnull()) recSwapDelaunay(d1);
    if (!d2.isnull()) recSwapDelaunay(d2);

    //    MESHC_LOG("TV = " << endl << M_->TVO());
    
    return d;
  }

  /*! Modified Alg 9.3 */
  Dart MeshC::splitEdgeDelaunay(const Dart& ed, int v)
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

    //    MESHC_LOG("TV = " << endl << M_->TVO());
    MESHC_LOG("Split edge " << ed << " with vertex " << v << endl);
    d = splitEdge(ed,v);
    
    if (!d0.isnull()) recSwapDelaunay(d0);
    if (!d1.isnull()) recSwapDelaunay(d1);
    if (!d2.isnull()) recSwapDelaunay(d2);
    if (!d3.isnull()) recSwapDelaunay(d3);

    //    MESHC_LOG("TV = " << endl << M_->TVO());
    
    return d;
  }

  Dart MeshC::bisectEdgeDelaunay(const Dart& d)
  {
    Point s;
    Dart dh(d);
    int v0(d.v());
    int v1(d.vo());
    int v2;
    bool segments(false);
    bool boundary(d.onBoundary());

    if (boundary ||
	((state_>=State_CDT) &&
	 (boundary_.segm(d) || interior_.segm(d)))) {
      dh.orbit2();
      v2 = dh.vo();
      MESHC_LOG("Bisect, enchroached 0: "
		<< M_->edgeEncroached(d,M_->S(v2)) << endl);
      if (M_->edgeEncroached(d,M_->S(v2))>0.0) {
	if (dh.onBoundary() ||
	    ((state_>=State_CDT) &&
	     (boundary_.segm(dh) || interior_.segm(dh)))) {
	  segments = true;
	  v0 = d.vo();
	  v1 = d.v();
	} else {
	  dh.orbit2();
	  segments = (dh.onBoundary() ||
		      ((state_>=State_CDT) &&
		       (boundary_.segm(dh) || interior_.segm(dh))));
	}
      }
      if ((!segments) && (!boundary)) {
	dh = d;
	dh.orbit0rev();
	v2 = dh.vo();
	MESHC_LOG("Bisect, enchroached 1: "
		  << M_->edgeEncroached(dh,M_->S(v2)) << endl);
	if (M_->edgeEncroached(d,M_->S(v2))>0.0) {
	  dh.orbit2();
	  if (dh.onBoundary() ||
	      ((state_>=State_CDT) &&
	       (boundary_.segm(dh) || interior_.segm(dh)))) {
	    segments = true;
	    v0 = d.vo();
	    v1 = d.v();
	  } else {
	    dh.orbit2rev();
	    segments =  (dh.onBoundary() ||
			 ((state_>=State_CDT) &&
			  (boundary_.segm(dh) || interior_.segm(dh))));
	  }
	}
      }
    }

    double beta = 0.5;
    if (segments) {
      MESHC_LOG("Adjacent segments, bisect at offcenter point.");
      const Point& s0(M_->S(v0));
      const Point& s1(M_->S(v1));
      const Point& s2(M_->S(v2));
      double l01 = M_->edgeLength(s0,s1);
      double l02 = M_->edgeLength(s0,s2);
      beta = l02/l01;
      if (beta>1./3.)
	while (beta>2./3.) beta *= 0.5;
      else if (beta<2./3.)
	while (beta<1./3.) beta *= 2.0;

      switch (M_->type()) {
      case Mesh::Mtype_manifold:
	/* Fall through to Mtype_plane behaviour; we have no
	   manifold-specific algorithm. */
      case Mesh::Mtype_plane:
	Vec::scale(s,s0,1.-beta);
	Vec::accum(s,s1,beta);
	break;
      case Mesh::Mtype_sphere:
	Vec::scale(s,s0,std::sin((1.-beta)*l01)/l01);
	Vec::accum(s,s1,std::sin(beta*l01)/l01);
	Vec::rescale(s,1./Vec::length(s));
	break;
      }
    } else {
      Vec::sum(s,M_->S(v0),M_->S(v1));
      Vec::rescale(s,0.5);
      
      switch (M_->type()) {
      case Mesh::Mtype_manifold:
	/* Fall through to Mtype_plane behaviour; we have no
	   manifold-specific algorithm. */
      case Mesh::Mtype_plane:
	/* Nothing to do! */
	break;
      case Mesh::Mtype_sphere:
	Vec::rescale(s,1./Vec::length(s));
	break;
      }
    }

    int v = addVertex(s);
    if ((state_>=State_RCDT) && big_.usingQv())
      big_.setQv(v,std::exp(std::log(big_.getQv(v0))*(1.-beta)+
			    std::log(big_.getQv(v1))*beta));
    dh = splitEdgeDelaunay(d,v);
    return dh;
  }


  /*!
    \f{align*}{
    r &= CR(s_0,s_1,c)
    r/\|s1-s0\| &> \beta
    \f}
  */
  void MeshC::calcSteinerPoint(const Dart& d, Point& c)
  {
    M_->triangleCircumcenter(d.t(),c);
    if ((M_->type() != Mesh::Mtype_sphere) &&
	(options_ & Option_offcenter_steiner)) {
      const double beta = ((state_>=State_RCDT) ?
			   skinny_.getQ(d.t()) :
			   std::sqrt(2.));
      Point len;
      const int argmin = M_->triangleEdgeLengthsArgMin(d.t(),len);
      const Point& s0 = M_->S(M_->TV(d.t())[(argmin+1)%3]);
      const Point& s1 = M_->S(M_->TV(d.t())[(argmin+2)%3]);
      const double radius = M_->triangleCircumcircleRadius(s0,s1,c);
      if (radius/len[argmin] <= beta) {
	MESHC_LOG("Steiner point (offcenter=circumcenter): " << c << endl);
      } else {
	MESHC_LOG("Steiner point could have been " << c << endl);
	MESHC_LOG("r/l = " << radius/len[argmin] << endl);
	Point s01;
	Vec::scale(s01,s0,0.5);
	Vec::accum(s01,s1,0.5);
	Vec::accum(c,s01,-1.0);
	MESHC_LOG("Distance before move: " << Vec::length(c) << endl);
	Vec::rescale(c,(len[argmin]*(beta+std::sqrt(beta*beta-0.25))
			/Vec::length(c)));
	MESHC_LOG("Distance after move: " << Vec::length(c) << endl);
	Vec::accum(c,s01);
	MESHC_LOG("Steiner point (offcenter): " << c << endl);
      }
    } else {
      MESHC_LOG("Steiner point (circumcenter): " << c << endl);
    }
  }



  /*! Alg 9.3 */
  Dart MeshC::insertNode(int v, const Dart& ed)
  {
    Dart td;

    MESHC_LOG("Locating node " << v
	      << " " << M_->S(v) << endl);

    if (M_->useVT()) {
      if (M_->VT(v) != -1) {/* Node already inserted! */
	MESHC_LOG("Node " << v << " already inserted." << endl);
	td = Dart(*M_,M_->VT(v));
	for (int k=0; k<3; k++) {
	  if (td.v() == v)
	    return td;
	  td.orbit2();
	}
      }
    }

    Dart ed0 = ed;
    MESHC_LOG("Trying, starting from dart " << ed0
	      << " " << M_->S(ed0.v())
	      << endl);
    if (M_->type() == Mesh::Mtype_sphere) {
      MESHC_LOG("PI minus distance to target " <<
		M_PI - M_->edgeLength(M_->S(v),M_->S(ed0.v())) << endl);
      if (M_PI - M_->edgeLength(M_->S(v),M_->S(ed0.v())) < 1e-6) {
	ed0.orbit2();
	MESHC_LOG("Trying, starting from dart " << ed0
		  << " " << M_->S(ed0.v())
		  << endl);
      }
    }

    td = M_->locate_point(ed0,M_->S(v),v);
    MESHC_LOG("Done looking." << endl);
    if (td.isnull()) { /* ERROR, not found! */
      MESHC_LOG("Error, node not found");
      return Dart();
    };
    if (td.v() == v) { /* Node already inserted! */
      MESHC_LOG("Node already inserted" << endl);
      return td;
    }
    td = Dart(*M_,td.t());
    Point bary;
    M_->barycentric(td,M_->S(v),bary);
    if ((bary[0]<-1000*MESH_EPSILON) ||
	(bary[0]<-1000*MESH_EPSILON) ||
	(bary[0]<-1000*MESH_EPSILON)) {
      MESHC_LOG("Triangle dart " << td
		<< "\n\t S[t]=("
		<< M_->S(M_->TV(td.t())[0]) << ",\n\t       "
		<< M_->S(M_->TV(td.t())[1]) << ",\n\t       "
		<< M_->S(M_->TV(td.t())[2]) << ")"
		<< "\n\t bary=" << bary
		<< "\n\t S[v]=" << M_->S(v)
		<< endl);
      MESHC_LOG("ERROR: locate_point returned triangle with bad barycentric coordinates for point.");
      return Dart();
    }
    size_t pattern(size_t(bary[0]>MESH_EPSILON)*1+
		   size_t(bary[1]>MESH_EPSILON)*2+
		   size_t(bary[2]>MESH_EPSILON)*4);
    MESHC_LOG("Triangle dart " << td
	      << "\n\t S[t]=("
	      << M_->S(M_->TV(td.t())[0]) << ",\n\t       "
	      << M_->S(M_->TV(td.t())[1]) << ",\n\t       "
	      << M_->S(M_->TV(td.t())[2]) << ")"
	      << "\n\t bary=" << bary
	      << "\n\t pattern=" << pattern
	      << "\n\t S[v]=" << M_->S(v)
	      << endl);
    switch (pattern) {
    case 6: // -++ Split e0
      td.orbit2();
      break;
    case 5: // +-+ Split e1
      td.orbit2rev();
      break;
    case 3: // ++- Split e2
      break;
    }
    if ((state_>=State_RCDT) && big_.usingQv())
      big_.setQv(v,std::exp(std::log(big_.getQv(M_->TV(td.t())[0]))*bary[0]+
		 std::log(big_.getQv(M_->TV(td.t())[1]))*bary[1]+
		 std::log(big_.getQv(M_->TV(td.t())[2]))*bary[2]));
    switch (pattern) {
    case 7: // +++
      return splitTriangleDelaunay(td,v);
      break;
    case 6: // -++ Split e0
    case 5: // +-+ Split e1
    case 3: // ++- Split e2
      MESHC_LOG("Edge dart " << td << endl);
      return splitEdgeDelaunay(td,v);
      break;
    case 1: // +-- Close to node 0, not allowed
    case 2: // -+- Close to node 1, not allowed
    case 4: // --+ Close to node 2, not allowed
      MESHC_LOG("Triangle dart " << td
		 << "\n\t S[t]=("
		 << M_->S(M_->TV(td.t())[0]) << ",\n\t       "
		 << M_->S(M_->TV(td.t())[1]) << ",\n\t       "
		 << M_->S(M_->TV(td.t())[2]) << ")"
		 << "\n\t bary=" << bary
		 << "\n\t pattern=" << pattern
		 << "\n\t S[v]=" << M_->S(v)
		 << endl);
      MESHC_LOG("ERROR: Attempt to add a duplicate point in triangle " << td << endl);
      break;
    case 0: // --- Close to all nodes, should not happen!
      MESHC_LOG("Close to all nodes, this should not happen, in triangle " << td << endl);
      break;
    }

    return Dart();
  }



  /*!  Calculate a convex covering of the convex hull of the points,
    as a convex geodesic polygon.  If the covering is larger than a
    hemisphere, the whole spere is used as cover.  The algorithm may
    overestimate the size of the convex hull.
    
    Let \f$(n,d)\f$ denote a plane with normal vector \f$n\f$ at
    (signed) distance \f$d\f$ from the origin.
    
    CET algorithm:
    \verbatim
    1. Find an enclosing circle defined by a plane (n,d).
    2. If d<=0
    3.   Find an enclosing geodesic polygon with  n  as reference point.
    4. else (enclosure may be more than a hemisphere)
    5.   Cover the whole sphere.
    \endverbatim

    Find an enclosing circle:
    \verbatim
    1. Set (n_0,d_0)=(s_0,1)
    2. For each point s_k, k=1,...,V-1
    3.   If n_{k-1}*s_k >= d_{k-1}
    4.     Set (n_k,d_k) = (n_{k-1},d_{k-1})
    5.   else
    6.     Calculate the plane (n_k,d_k) that covers both
           s_k and the circle defined by (n_{k-1},d_{k-1})
    \endverbatim

    Find plane \f$(n_1,d_1)\f$ that covers \f$s_1\f$ and the circle
    defined by \f$(n_0,d_0)\f$, given that \f$n_0\cdot s_1 < d_0\f$:
    \f{align*}{
    c &= n_0\cdot s_1 \\
    1-c^2 &= \|n_0\times s_1\|^2 \\
    s'_1 &= n_0d_0 -(s_1-n_0 c) b, \quad 0\leq b\leq 1 \\
    1 &= d_0^2 + (1-c^2) b^2 \\
    b &= \sqrt{1-d_0^2}/\|n_0\times s_1\| \\
    n'_1 &= (s_1-s'_1)\times (n_0\times s_1) \\
    n_1 &= n'_1/\|n'_1\| \\
    d_1 &= \min( n_1 \cdot s_1, n_1 \cdot s'_1 )
    \f}
    Since \f$n_0\cdot s_1 < d_0\f$, only the case \f$s_1=-n_0\f$ is
    degenerate, in which case \f$s'_1\f$ is set to
    \f$n_0d_0+n'_0\sqrt{1-d_0^2}\f$, where \f$n'_0\f$ is an arbitrary
    perpendicular point to \f$n_0\f$.

    \see MeshC::CET
  */
  bool MeshC::CETsphere(int sides, double margin)
  {
    if (state_ != State_noT) {
      MESHC_LOG("Cannot add convex enclosure to existing triangulation.");
      return false;
    }

    if (M_->type() != Mesh::Mtype_sphere) {
      MESHC_LOG("Mesh type mismatch: "
		<< M_->type() << " should be "
		<< "Mesh::Mtype_sphere" << endl);
      return false;
    }

    int nV = (int)M_->nV();

    if (nV<1) {
      MESHC_LOG("Need at least one vertex to calculate enclosure.");
      return false;
    }

    if (sides<3)
      sides = 3;

    int i;

    /* Calculate a covering circle. */
    MESHC_LOG("Calculate a covering circle.");

    Point n0, n0s1, s1prime, sh;
    Vec::copy(n0,M_->S(0));
    double d0 = 1.0;
    double nc,ns,b;
    
    for (i=1;i<nV;i++) {
      nc = Vec::scalar(n0,M_->S(i));
      Vec::cross(n0s1,n0,M_->S(i));
      ns = Vec::length(n0s1);
      if (nc < d0) {
	if (ns > 0.0) {
	  b = std::sqrt(1.0-d0*d0)/ns;
	  Vec::scale(s1prime,n0,d0+nc*b);
	  Vec::accum(s1prime,M_->S(i),-b);
	} else {
	  Point n0prime;
	  Vec::arbitrary_perpendicular(n0prime,n0);
	  Vec::scale(s1prime,n0,d0);
	  Vec::accum(s1prime,n0prime,std::sqrt(1.0-d0*d0));
	}
	Vec::diff(sh,M_->S(i),s1prime);
	Vec::cross(n0,sh,n0s1);
	Vec::rescale(n0,1.0/Vec::length(n0));
	d0 = Vec::scalar(n0,M_->S(i));
	/* For robustness, check s1prime as well: */
	b = Vec::scalar(n0,s1prime);
	if (b<d0) {
	  MESHC_LOG("min(" << d0 << ", " << b << ")" << endl);
	  d0 = b;
	}

	MESHC_LOG("n0: " << n0 << ", " << d0 << endl);
      }
    }

    /* Calculate margin */    
    double diameter = 2.0*std::acos(d0);
    if (margin<0.0) {
      margin = -diameter*margin;
    }
    
    MESHC_LOG("diameter = " << diameter << endl);
    MESHC_LOG("margin = " << margin << endl);

    if (diameter+2*margin+margin>=M_PI) {
      /* The whole sphere needs to be covered. */
      MESHC_LOG("Cover the whole sphere." << endl);

      if (nV<3) { /* Not enough points for even one triangle,
		     needs special treatment. */
	NOT_IMPLEMENTED;
	MESHC_LOG("nV=" << nV << endl);
	return false;
      }
      
      /*
	1. Pick a point v0.
	2. Find the point v1 that minimises (s0s1+1/3)^2
	3. Find the point v2 that minimises
	     (s0s2+1/3)^2+(s1s2+1/3)^2
	4. Reorder to a CCW triangle (s0,s1,s2)
	5. Find any point v3 in the triangle (-s0,-s1,-s2)
	6. If no point found in 5., add a vertex at -(s0+s1+s2)/\|s0+s1+s2\|
	6. Create triangles (v0,v1,v2), (v2,v1,v3), (v0,v2,v3), (v1,v0,v3)
       */

      int v0 = 0;
      Point const * s0 = &(M_->S(v0));
      int v1 = -1;
      int v2 = -1;
      int v3 = -1;
      Point const * s1 = NULL;
      Point const * s2 = NULL;

      MESHC_LOG("First point," << 
		" v0=" << v0 <<
		" s0=" << *s0 << endl);

      /* Find suitable v1: */
      double loss = 16.0/9.0+1.0;
      for (int v=0; v<nV; v++) {
	if (v==v0) continue;
	double loss_ = (Vec::scalar(*s0,M_->S(v))+1.0/3.0);
	loss_ *= loss_;
	if (loss_ < loss) {
	  loss = loss_;
	  v1 = v;
	}
      }
      s1 = &(M_->S(v1));

      MESHC_LOG("Second point," << 
		" v1=" << v1 <<
		" s1=" << *s1 << endl);

      /* Find suitable v2: */
      loss = 16.0/9.0+16.0/9.0+1.0;
      for (int v=0; v<nV; v++) {
	if ((v==v0) || (v==v1)) continue;
	double loss0_ = (Vec::scalar(*s0,M_->S(v))+1.0/3.0);
	double loss1_ = (Vec::scalar(*s1,M_->S(v))+1.0/3.0);
	double loss_ = loss0_*loss0_+loss1_*loss1_;
	if (loss_ < loss) {
	  loss = loss_;
	  v2 = v;
	}
      }
      s2 = &(M_->S(v2));

      MESHC_LOG("Third point," << 
		" v2=" << v2 <<
		" s2=" << *s2 << endl);

      /* Make sure we have a CCW triangle: */
      if (Vec::volume(*s0,*s1,*s2) < 0.0) {
	int v = v2;
	v2 = v1;
	v1 = v;
	s1 = &(M_->S(v1));
	s2 = &(M_->S(v2));

	MESHC_LOG("Swapped second and third point." << endl);
      }

      /* Calculate the inward normals of the triangle. */
      Point s0xs1; Vec::cross(s0xs1,*s0,*s1);
      Point s1xs2; Vec::cross(s1xs2,*s1,*s2);
      Point s2xs0; Vec::cross(s2xs0,*s2,*s0);

      /* Find a point in the opposing triangle: */
      Point outside;
      outside[0] = 0.0;
      outside[1] = 0.0;
      outside[2] = 0.0;
      for (int v=0; v<nV; v++) {
	if ((v==v0) || (v==v1) || (v==v2)) continue;
	Point outside_;
	outside_[0] = Vec::scalar(s1xs2,M_->S(v));
	outside_[1] = Vec::scalar(s2xs0,M_->S(v));
	outside_[2] = Vec::scalar(s0xs1,M_->S(v));
	if ((outside_[0] < outside[0]) &&
	    (outside_[1] < outside[1]) &&
	    (outside_[2] < outside[2])) {
	  Vec::copy(outside,outside_);
	  v3 = v;
	}
      }
      if (v3<0) {
	Point s3;
	Vec::sum(s3,*s0,*s1);
	Vec::accum(s3,*s2);
	Vec::rescale(s3,-1.0/Vec::length(s3));
	v3 = addVertex(s3);
	MESHC_LOG("Needed to add an extra vertex." << endl);
      }

      MESHC_LOG("Fourth point," << 
		" v3=" << v3 <<
		" s3=" << M_->S(v3) << endl);

      /* Create triangles: */
      Matrix3int TV(4);
      TV(0) = Int3(v0,v1,v2);
      TV(1) = Int3(v3,v2,v1);
      TV(2) = Int3(v3,v0,v2);
      TV(3) = Int3(v3,v1,v0);
      M_->TV_append(TV);
    } else {
      /* Calculate tight enclosure. */
      MESHC_LOG("Calculate tight enclosure." << endl);

      /* Construct interior boundary normals. */
      /* This initialises the enclosure. */
      Point n1, n2;
      //      Vec::arbitrary_perpendicular(n1,n0);
      //     Vec::cross(n2,n0,n1);
      if (n0[2]>0.9) {
	n1[0] = 0.;
	n1[1] = 1.;
	n1[2] = 0.;
	Vec::cross(n2,n0,n1);
	Vec::rescale(n2,1.0/Vec::length(n2));
	Vec::cross(n1,n2,n0);
	Vec::rescale(n1,1.0/Vec::length(n1));
      } else if (n0[2]<-0.9) {
	n1[0] = 0.;
	n1[1] = 1.;
	n1[2] = 0.;
	Vec::cross(n2,n0,n1);
	Vec::rescale(n2,1.0/Vec::length(n2));
	Vec::cross(n1,n2,n0);
	Vec::rescale(n1,1.0/Vec::length(n1));
      } else {
	n2[0] = 0.;
	n2[1] = 0.;
	n2[2] = 1.;
	Vec::cross(n1,n2,n0);
	Vec::rescale(n1,1.0/Vec::length(n1));
	Vec::cross(n2,n0,n1);
	Vec::rescale(n2,1.0/Vec::length(n2));
      }

      Matrix3double n(sides); /* Normal vectors. */
      double th;
      for (i=0;i<sides;i++) {
	th = 2.0*M_PI*double(i)/double(sides);
	Vec::scale(n(i),n1,-std::sin(th));
	Vec::accum(n(i),n2,std::cos(th));
      }
      
      double dist;
      for (int v=0;v<nV;v++) {
	for (i=0;i<sides;i++) {
	  dist = Vec::scalar(n[i],M_->S(v));
	  if (dist < 0.0) { /* Update enclosure. */
	    Vec::cross(sh,n0,n[i]);
	    Vec::cross(n(i),sh,M_->S(v));
	    Vec::rescale(n(i),1.0/Vec::length(n[i]));
	  }
	}
      }

      MESHC_LOG("Add margin." << endl);
      {
	double th;
	double margini;
	for (i=0;i<sides;i++) {
	  nc = Vec::scalar(n0,n[i]);
	  Vec::cross(sh,n0,n[i]);
	  ns = Vec::length(sh);
	  Vec::rescale(sh,1.0/ns);
	  th = std::atan2(ns,nc);
	  margini = margin*ns/d0;
	  if (th-margini > 0.0) {
	    nc = std::cos(th-margini);
	    ns = std::sin(th-margini);
	  } else {
	    cout << WHEREAMI << "Oops! Perhaps some NA input?" << endl;
	    nc = 1.0;
	    ns = 0.0;
	  }
	  Vec::cross(n(i),sh,n0);
	  Vec::rescale(n(i),ns);
	  Vec::accum(n(i),n0,nc);
	  Vec::rescale(n(i),1.0/Vec::length(n[i]));
	}
      }

      /* Calculate intersections. */
      MESHC_LOG("Calculate enclosure boundary." << endl);
      Matrix3double S(sides);
      {
	Point nip, nipp;
	double nip_nj, nipp_nj;
	double bi;
	int j;
	for (i=0;i<sides;i++) {
	  j = (i+1)%sides;
	  Vec::cross(nip,n[i],n0);
	  Vec::rescale(nip,1.0/Vec::length(nip));
	  Vec::cross(nipp,nip,n[i]);
	  Vec::rescale(nipp,1.0/Vec::length(nipp));
	  nip_nj = Vec::scalar(nip,n[j]);
	  nipp_nj = Vec::scalar(nipp,n[j]);
	  bi = std::sqrt(nip_nj*nip_nj+nipp_nj*nipp_nj);
	  Vec::scale(S(j),nip,nipp_nj/bi);
	  Vec::accum(S(j),nipp,-nip_nj/bi);
	}
      }

      /* Construct enclosure triangles. */
      MESHC_LOG("Construct enclosure triangles." << endl);
      Matrix3int TV(sides-2);
      for (i=0;i<sides-2;i++) {
	TV(i) = Int3(nV+(0),
		     nV+(i+1),
		     nV+((i+2)%sides));
      }

      M_->S_append(S);
      M_->TV_append(TV);

    }
    
    MESHC_LOG("CET finished" << endl << *this);
#ifndef FMESHER_NO_X
    M_->redrawX11("CET finished");
#endif
    
    state_ = State_CET;
    return true;
  }


  /*!
    Calculate a convex covering of the convex hull of the points, as a
    convex polygon.

    Let \f$(n,d)\f$ denote a line with normal vector \f$n\f$ at
    (signed) distance \f$d\f$ from the origin.
    
    Intersection \f$s_{01}\f$ between two planes \f$(n_0,d_0)\f$ and
    \f$(n_1,d_1)\f$:
    \f{align*}{
    s_{01} &= n_0a_0+n_1a_1 \\
    n_0\cdot s_{01} &= d_0 \\
    n_1\cdot s_{01} &= d_1 \\
    n_{01} &= n_0\cdot n_1 \\
    a_0+n_{01}a_1 &= d_0 \\
    n_{01}a_0+a_1 &= d_1 \\
    a_0 &= (d_0-d_1 n_{01})/(1-n_{01}^2) \\
    a_1 &= (d_1-d_0 n_{01})/(1-n_{01}^2)
    \f}

    \see MeshC::CET
  */
  bool MeshC::CETplane(int sides, double margin)
  {
    if (state_ != State_noT)
      return false; /* Cannot TODO: Add convex enclosure? */

    if (M_->type() != Mesh::Mtype_plane) {
      MESHC_LOG("Mesh type mismatch: "
		<< M_->type() << " should be "
		<< "Mesh::Mtype_plane" << endl);
      return false;
    }

    int nV = (int)M_->nV();

    if (nV<1) {
      MESHC_LOG("Need at least one vertex to calculate enclosure." << endl);
      return false;
    }

    if (sides<3)
      sides = 3;

    /* Calculate tight enclosure. */
    MESHC_LOG("Calculate tight enclosure." << endl);

    int i;

    /* Construct interior boundary normals. */
    std::vector<Point> n(sides); /* Normal vectors. */
    double th;
    for (i=0;i<sides;i++) {
      th = 2.0*M_PI*double(i)/double(sides);
      n[i][0] = -std::sin(th);
      n[i][1] = std::cos(th);
      n[i][2] = 0.0;
    }
    
    /* Initialise enclosure. */
    std::vector<double> d(sides); /* Distances from origin for boundary. */
    for (i=0;i<sides;i++) {
      d[i] = Vec::scalar(n[i],M_->S(0));
    }

    double dist;
    for (int v=1;v<nV;v++) {
      for (i=0;i<sides;i++) {
	dist = Vec::scalar(n[i],M_->S(v));
	if (dist < d[i])
	  d[i] = dist;
      }
    }

    /* Calculate margin */    
    if (margin<0.0) {
      double diameter(0.0);
      double diam;
      if ((sides%2) == 0) { /* Each side has an opposite. */
	for (i=0;i<sides/2;i++) {
	  diam = -d[i]-d[(i+sides/2)%sides];
	  if (diam>diameter)
	    diameter = diam;
	}
	margin = -diameter*margin;
      } else {
	MESHC_LOG("Calculate margin." << endl);
	for (i=0;i<sides/2;i++) {
	  diam = -d[i]-d[(i+sides/2)%sides];
	  if (diam>diameter)
	    diameter = diam;
	  diam = -d[i]-d[(i+sides/2+1)%sides];
	  if (diam>diameter)
	    diameter = diam;
	}
	margin = -diameter*margin;
      }
    }
    
    MESHC_LOG("margin = " << margin << endl);

    MESHC_LOG("Add margin." << endl);
    for (i=0;i<sides;i++) {
      d[i] -= margin;
    }

    MESHC_LOG("Calculate enclosure boundary." << endl);

    std::vector<Point> S(sides);
    double a0, a1, n01;
    int j;
    for (i=0;i<sides;i++) {
      j = (i+1)%sides;
      n01 = Vec::scalar(n[i],n[j]);
      a0 = (d[i]-d[j]*n01);
      a1 = (d[j]-d[i]*n01);
      n01 = 1-n01*n01;
      Vec::scale(S[j],n[i],a0/n01);
      Vec::accum(S[j],n[j],a1/n01);
    }

    /* Add enclosure triangles. */
    MESHC_LOG("Add enclosure triangles." << endl);
    std::vector<Int3> TV(sides-2);
    for (i=0;i<sides-2;i++) {
      TV[i][0] = nV+(0);
      TV[i][1] = nV+(i+1);
      TV[i][2] = nV+((i+2)%sides);
    }

    M_->S_append(Matrix3double(sides,S.data()));
    M_->TV_append(Matrix3int(sides-2,TV.data()));

    MESHC_LOG("CET finished" << endl << *this);
#ifndef FMESHER_NO_X
    M_->redrawX11("CET finished");
#endif
    
    state_ = State_CET;
    return true;
  }

  bool MeshC::CET(int sides, double margin)
  {
    if (state_ != State_noT)
      return false; /* Cannot add enclosure to an existing triangulation. */

    switch (M_->type()) {
    case Mesh::Mtype_plane:
      return CETplane(sides,margin);
      break;
    case Mesh::Mtype_sphere:
      return CETsphere(sides,margin);
      break;
    default:
      MESHC_LOG("Only planar enclosures implemented yet." << endl);
      NOT_IMPLEMENTED;
      return false;
    }
  }








  bool MeshC::DT(const vertexListT& v_set)
  {
    if (is_pruned_) 
      return false; /* ERROR, cannot safely insert nodes into a pruned
		       triangulation. Call insertNode directly if known to
		       be visible/reachable from a given edge.  */

    if (state_ < State_DT)
      if (!prepareDT()) /* Make sure we have a DT. */
	return false;

    int v;
    vertexListT::const_iterator v_iter;
    Dart dh;

    dh = Dart();
    for (v_iter = v_set.begin(); v_iter != v_set.end(); v_iter++) {
      v = *v_iter;
      if (dh.isnull()) dh = Dart(*M_,0);
      dh = insertNode(v,dh); /* Start looking where the previous
				point was found. */
      if (dh.isnull()) {
	MESHC_LOG("DT: Failed to insert node " << v << endl << *this);
      }
    }
      
    MESHC_LOG("DT finished" << endl << *this);
#ifndef FMESHER_NO_X
    M_->redrawX11("DT finished");
#endif

    return true;
  }


  bool MeshC::prepareDT()
  {
    if (state_ < State_CET)
      if (!CET(4)) /* Construct an basic CET. */
	return false;

    if (state_<State_DT) {
      /* We need to make sure the triangulation is a DT. */
      triangleSetT t_set;
      for (int t=0;t<(int)M_->nT();t++)
	t_set.insert(t);
      if (LOP(t_set))
	state_ = State_DT;
    }
    return (state_>=State_DT);
  }


  bool MeshC::prepareCDT()
  {
    if (!prepareDT()) return false; /* Make sure we have a DT. */
    if (state_>=State_CDT)
      return true; /* Nothing to do. Data structures already active. */

    int vi;
    Dart d;
    for (int t=0;t<(int)M_->nT();t++) {
      const Int3& tt = M_->TT()[t];
      for (vi=0;vi<3;vi++)
	if (tt[vi]<0) {
	  d = Dart(*M_,t,1,(vi+1)%3);
	  boundary_.insert(d,MCQsegm::meta_type());
	}
    }

    state_ = State_CDT;
    return true;
  }

  bool MeshC::prepareRCDT(double skinny_limit,
			  double big_limit,
			  const double* big_limits,
			  size_t nQL,
			  int max_n0,
			  int max_n1)
  {
    if (!prepareCDT()) return false; /* Make sure we have a CDT. */

    skinny_.clear();
    big_.clear();
    skinny_.setQ(skinny_limit);
    big_.setQ(big_limit,big_limits,nQL);

    for (int t=0;t<(int)M_->nT();t++) {
      skinny_.insert(Dart(*M_,t));
      big_.insert(Dart(*M_,t));
    }

    max_n0_ = max_n0;
    max_n1_ = max_n1;

    state_ = State_RCDT;
    return true;
  }



  bool MeshC::CDTBoundary(const constrListT& constr)
  {
    if (!prepareCDT()) return false;

    constr_boundary_ = constrListT(constr.begin(),constr.end());
    
    return buildCDT();
  }

  bool MeshC::CDTInterior(const constrListT& constr)
  {
    if (!prepareCDT()) return false;

    constr_interior_ = constrListT(constr.begin(),constr.end());
    
    return buildCDT();
  }

  bool MeshC::CDT(const constrListT& boundary, const constrListT& interior)
  {
    if (!prepareCDT()) return false;

    constr_boundary_ = constrListT(boundary.begin(),boundary.end());
    constr_interior_ = constrListT(interior.begin(),interior.end());
    
    return buildCDT();
  }






  bool MeshC::LOP(MCQswapableD& swapable)
  {
    MESHC_LOG("LOP swapable: "
	      << swapable.countQ() << "/" << swapable.count() << endl);
    /* Swap edges, until none are swapable. */
    Dart dh;
    while (!swapable.emptyQ()) {
      dh = swapable.beginQ()->d_; /* d_ may get erased in swapEdge! */
      swapEdge(dh,swapable);
      MESHC_LOG("LOP swapable: "
		<< swapable.countQ() << "/" << swapable.count() << endl);
    }

    MESHC_LOG("LOP finished" << endl << *this);

    return true;
  }

  bool MeshC::LOP(const triangleSetT& t_set)
  {
    /* Locate interior edges */
    Dart dh, dh2;
    MCQswapableD swapable(this);
    for (triangleSetT::const_iterator ci=t_set.begin();
	 ci != t_set.end(); ci++) {
      dh = Dart(*M_,(*ci));
      for (int vi=0; vi<3; vi++) {
	dh2 = dh;
	dh2.orbit1();
	if ((dh.t() != dh2.t()) /* Only add if not on boundary */
	    && (t_set.find(dh2.t()) != t_set.end())
	    /* Only add if the neighbouring triangle is also in the set. */
	    && ((state_<State_CDT)
		|| ((!boundary_.segm(dh))
		    && (!interior_.segm(dh))))) /* Don't add CDT segments. */
	  swapable.insert(dh); /* MCQswapableD takes care of duplicates. */
	dh.orbit2();
      }
    }

    return LOP(swapable);
  }



  typedef std::list<IntPair> BoundaryList;
  
  void prevnext(BoundaryList::iterator& prev,
		BoundaryList::iterator& curr,
		BoundaryList::iterator& next)
  {
    for (curr--; curr->second == 0; curr--) {};
    for (curr++; curr->second == 0; curr++) {};
    next = curr;
    for (next++; next->second == 0; next++) {};
    prev = curr;
    for (prev--; prev->second == 0; prev--) {};
  }

  void prev(BoundaryList::iterator& prev,
	    BoundaryList::iterator& curr,
	    BoundaryList::iterator& next)
  {
    for (curr--; curr->second == 0; curr--) {};
    next = curr;
    for (next++; next->second == 0; next++) {};
    prev = curr;
    for (prev--; prev->second == 0; prev--) {};
  }

  void next(BoundaryList::iterator& prev,
	    BoundaryList::iterator& curr,
	    BoundaryList::iterator& next)
  {
    for (curr++; curr->second == 0; curr++) {};
    next = curr;
    for (next++; next->second == 0; next++) {};
    prev = curr;
    for (prev--; prev->second == 0; prev--) {};
  }

  void prevnext(BoundaryList::reverse_iterator& prev,
		BoundaryList::reverse_iterator& curr,
		BoundaryList::reverse_iterator& next)
  {
    for (curr--; curr->second == 0; curr--) {};
    for (curr++; curr->second == 0; curr++) {};
    next = curr;
    for (next++; next->second == 0; next++) {};
    prev = curr;
    for (prev--; prev->second == 0; prev--) {};
  }

  void prev(BoundaryList::reverse_iterator& prev,
	    BoundaryList::reverse_iterator& curr,
	    BoundaryList::reverse_iterator& next)
  {
    for (curr--; curr->second == 0; curr--) {};
    next = curr;
    for (next++; next->second == 0; next++) {};
    prev = curr;
    for (prev--; prev->second == 0; prev--) {};
  }

  void next(BoundaryList::reverse_iterator& prev,
	    BoundaryList::reverse_iterator& curr,
	    BoundaryList::reverse_iterator& next)
  {
    for (curr++; curr->second == 0; curr++) {};
    next = curr;
    for (next++; next->second == 0; next++) {};
    prev = curr;
    for (prev--; prev->second == 0; prev--) {};
  }



  /*! Alg. 6.2+6.1 */
  Dart MeshC::CDTInsertSegment(const DartPair& dp,
			       const DartList& trace,
			       triangleSetT& triangles,
			       const bool is_boundary,
			       const constrMetaT& meta)
  {
    if (!prepareCDT()) return Dart();

    Dart dc;

    MESHC_LOG(dp);
    MESHC_LOG(trace);

    MESHC_LOG("Segment crosses " << trace.size()
	      << " edges." << endl);

    Dart dh;
    Dart d0(dp.first);
    Dart d1(dp.second);
    BoundaryList boundary0;
    BoundaryList boundary1;

    int v0(d0.v());
    int v1(d1.v());

    if (trace.size()<=1) {
      /* Only one crossing edge, swap directly. */
      MESHC_LOG("Short trace (" << trace.size() << ").  Swapping directly."
		<< endl);
      
      Dart ds = swapEdge(*trace.begin());
      if (ds.v() == v1) ds.orbit1();
      if (ds.v() != v0) ds = Dart();
      
      MESHC_LOG("Segment dart " << ds << endl);

      if (!ds.isnull()) {      
	if (is_boundary)
	  boundary_.insert(ds, meta);
	else
	  interior_.insert(ds, meta);
      }
      
      return ds;
    }

    for (DartList::const_iterator i(trace.begin());
	 i != trace.end(); i++) {
      dh = *i;
      triangles.insert(dh.t());

#ifndef FMESHER_NO_X
      if (M_->useX11())
	M_->drawX11triangle(dh.t(),true);
#endif

      if (v0 == dh.v()) {
	boundary0.rbegin()->second++;
      } else {
	v0 = dh.v();
	boundary0.push_back(IntPair(v0,1));
      }
      if (v1 == dh.vo()) {
	boundary1.begin()->second++;
      } else {
	v1 = dh.vo();
	boundary1.push_front(IntPair(v1,1));
      }
    }
    dh.alpha2();
    triangles.insert(dh.t());

#ifndef FMESHER_NO_X
    if (M_->useX11())
      M_->drawX11triangle(dh.t(),true);
#endif

    v0 = d0.v();
    v1 = d1.v();

    MESHC_LOG(M_->S(v0) << endl);
    MESHC_LOG(M_->S(v1) << endl);

    boundary0.push_front(IntPair(v0,-1));
    boundary0.push_back(IntPair(v1,-1));
    boundary1.push_back(IntPair(v0,-1));
    boundary1.push_front(IntPair(v1,-1));

    MESHC_LOG("T:  " << triangles << endl);

    BoundaryList::iterator i0, i0_prev, i0_next;
    BoundaryList::reverse_iterator i1, i1_prev, i1_next;

#define	CDTMSG(msg)						\
    MESHC_LOG(msg << endl);				\
    MESHC_LOG("B1: " << boundary1 << endl);		\
    MESHC_LOG("B0: " << boundary0 << endl);		\
    MESHC_LOG("I1: "						\
	      << *i1_prev << *i1 << *i1_next << endl);	\
    MESHC_LOG("I0: "						\
	      << *i0_prev << *i0 << *i0_next << endl);	\
    MESHC_LOG("d0: " << d0 << endl);			\
    MESHC_LOG("vd0: " << vd0 << endl);			\
    MESHC_LOG("dh: " << dh << endl);

    /* Initialise */
    i0 = boundary0.begin(); i0_prev = i0; i0_next = i0;
    i1 = boundary1.rbegin(); i1_prev = i1; i1_next = i1;
    next(i0_prev,i0,i0_next);
    next(i1_prev,i1,i1_next);

    /* Go to first edge */
    Dart vd0(d0); /* First dart in edge bundle from current vertex. */
    vd0.orbit2();
    while (i1->first != vd0.vo()) next(i1_prev,i1,i1_next);

    dh = vd0;

    //    CDTMSG("");
    while (true) {
      bool swapable = true;
      int v10 = vd0.vo(); /* The first opposite vertex. */
      while (swapable) {
	CDTMSG("");
	/* Find swapable edge */
	dh = vd0;
	for (swapable = dh.isSwapable();
	     (dh.vo() != i0_next->first) && (!swapable);
	     swapable = dh.isSwapable()) {
	  dh.orbit0rev();
	  CDTMSG("Looking for swapable");
	  MESHC_LOG("swapable=" << swapable << endl);
	};
	swapable = swapable && (dh.vo() != i0_next->first);
	if (swapable) {
	  /* Swap the edge */
	  bool vd0affected = (dh==vd0); /* vd0 must be updated? */

	  while (i1->first != dh.vo()) next(i1_prev,i1,i1_next);
	  
	  CDTMSG("Before edge swap");

	  i0->second--;
	  i1->second--;
	  dh = swapEdge(dh);
	  MESHC_LOG("dh: " << dh << endl);

	  if (dh.vo() != v0) {
	    if (i0_prev->first == dh.vo()) i0_prev->second++;
	    if (i1_prev->first == dh.vo()) i1_prev->second++;
	  }
	  if (dh.v() != v1) {
	    if (i0_next->first == dh.v()) i0_next->second++;
	    if (i1_next->first == dh.v()) i1_next->second++;
	  }

	  if (dh.vo() == v0) {
	    /* Triangle is linked to v0, needs to regenerate d0. */
	    d0 = dh;
	    d0.orbit2();
	    MESHC_LOG("d0: " << d0 << endl);
	  }

	  CDTMSG("Have swapped, and adjusted counters");

	  if (dh.v() == i0_next->first) {
	    if (dh.vo() == i0_prev->first) {
	      /* Case: i0_next --> i0_prev */
	      /* I0 has been eliminated */
	      if (dh.vo() == v0) {
		/* i0_prev == v0 */
		d0.orbit0();
		MESHC_LOG("d0: " << d0 << endl);
	      }
	      if (i0_prev->first != v0) i0_prev->second--;
	      if (i0_next->first != v1) i0_next->second--;
	    } else if (dh.vo() == i1_prev->first) {
	      /* Case: i0_next --> i1_prev */
	      /* Nothing to do. */
	    }
	  } else if (dh.v() == i1_next->first) {
	    if (dh.vo() == i1_prev->first) {
	      /* Case: i1_next --> i1_prev */
	      /* I1 has been eliminated */
	      if (i1_prev->first != v0) i1_prev->second--;
	      if (i1_next->first != v1) i1_next->second--;
	    } else if (dh.vo() == i1_prev->first) {
	      /* Case: i1_next --> i0_prev */
	      /* Nothing to do. */
	    }
	  }

	  CDTMSG("Have swapped, and possibly eliminated vertices");

	  dh.orbit2rev();
	  if (vd0affected) {
	    vd0 = dh;
	    v10 = vd0.vo();
	  }
	  swapable = (vd0.vo() != i0_next->first);
	  
	  if (swapable && (i0->second>=0)) {
	    prevnext(i1_prev,i1,i1_next);
	    while ((i1_prev->second >= 0) && (i1->first != v10))
	      prev(i1_prev,i1,i1_next);
	  }

	  CDTMSG("Have swapped, and tried to update location");
	  MESHC_LOG((swapable ?
		    "(Restart this vertex)" :
		    "(Leave vertex)") << endl);

      	}
      }

      CDTMSG("Done with vertex, trying to move on.");

      if (d0.vo() == v1) /* Segment has been inserted */
	break; /* This is the only way out of the loop! */

      if (i0->second==0) {
	/* The vertex was eliminated, start from beginning. */
	i0 = boundary0.begin(); i0_prev = i0; i0_next = i0;
	i1 = boundary1.rbegin(); i1_prev = i1; i1_next = i1;
	next(i0_prev,i0,i0_next);
	next(i1_prev,i1,i1_next);
	
	/* Go to first edge */
	vd0 = d0; /* First dart in edge bundle from current vertex. */
	vd0.orbit2();

	CDTMSG("Trying to find opposite vertex.");

	while (i1->first != vd0.vo()) next(i1_prev,i1,i1_next);

	CDTMSG("Vertex eliminated, start from beginning:");
      } else { /* The vertex was not eliminated, go to next. */
	next(i0_prev,i0,i0_next);
	dh.orbit2(); /* First dart in edge bundle from current vertex. */
	vd0 = dh;
	while (i1->first != vd0.vo()) next(i1_prev,i1,i1_next);
	CDTMSG("Vertex not eliminated, go to next:");
      }

    }

    MESHC_LOG("Segment inserted:" << endl);
    MESHC_LOG("B1: " << boundary1 << endl);
    MESHC_LOG("B0: " << boundary0 << endl);
    MESHC_LOG("d0: " << d0 << endl);

    dc = d0;

    MESHC_LOG("Segment dart " << dc << endl);

    if (!dc.isnull()) {      
      if (is_boundary)
	boundary_.insert(dc, meta);
      else
	interior_.insert(dc, meta);
    }

    return dc;
  }


  int MeshC::CDTSplitSegment(const DartPair& dp, const DartList& trace)
  {
    MESHC_LOG(dp);
    MESHC_LOG("Edge trace:" << endl << trace);

    Dart dh;
    Dart d0(dp.first);
    Dart d1(dp.second);

    int v0(d0.v());
    int v1(d1.v());
    bool split(false);
    double delta;
    int dhv1;
    for (DartList::const_iterator i(trace.begin());
	 i != trace.end(); i++) {
      dh = *i;
      MESHC_LOG("Testing edge for interference: " << dh << endl)

      MESHC_LOG("Testing vertex for interference: " << dh.v() << endl);
      delta = M_->inLeftHalfspace(M_->S(v0),M_->S(v1),M_->S(dh.v()));
      split = ((delta >= -MESH_EPSILON) & (delta <= MESH_EPSILON));

      if (!split) {
	/* Test the other edge vertex. */ 
	dh.orbit2();
	dhv1 = dh.v();
	MESHC_LOG("Testing vertex for interference: " << dhv1 << endl);
	delta = M_->inLeftHalfspace(M_->S(v0),M_->S(v1),M_->S(dhv1));
	split = ((delta >= -MESH_EPSILON) & (delta <= MESH_EPSILON));

	if (!split) {
	  dh.orbit2rev();
	  /* Test for interfering segment. */ 
	  split = isSegment(dh);
	  if (split) {
	    MESHC_LOG("Segment interference detected." << endl);

	    Point c;
	    double beta = M_->edgeIntersection(M_->S(v0), M_->S(v1),
					       M_->S(dh.v()), M_->S(dhv1),
					       c);
	    MESHC_LOG("Splitting segment at " << c << endl);

	    int v = addVertex(c);
	    if ((state_>=State_RCDT) && big_.usingQv())
	      big_.setQv(v,std::exp(std::log(big_.getQv(dh.v()))*(1.-beta)+
				    std::log(big_.getQv(dhv1))*beta));
	    dh = splitEdgeDelaunay(dh, v);
	  }
	}
      }

      if (split) {
	/* Vertex is on the new segment line: split. */
	MESHC_LOG("Splitting new segment at vertex " << dh.v() << endl)
	  return dh.v();
      }

      MESHC_LOG("No interference from edge detected.")
    }

    MESHC_LOG("No interference from trace detected.")

    return -1;
  }


  Dart MeshC::CDTInsertSegment(const int v0, const int v1,
			       triangleSetT& triangles,
			       const bool is_boundary,
			       const constrMetaT& meta)
  {
    if (!prepareCDT()) return Dart();
    MESHC_LOG("Inserting segment ("
	      << v0 << "," << v1 << ")" << endl);
    if (v0 == v1) return Dart();

    DartList trace;
    Dart dh(M_->locate_vertex(Dart(),v0));
    if (dh.isnull()) {
      MESHC_LOG("Originating vertex not found "
		<< v0 << endl);
      return Dart();
    }
    DartPair dhp(M_->trace_path(dh,M_->S(v1),v1,&trace));
    if (dhp.second.isnull()) {
      MESHC_LOG("Endpoint vertex not found ("
		<< v0 << "," << v1 << ") "
		<< M_->S(v0) << ", " << M_->S(v0) << endl);
      return Dart();
    }

    Dart dh0(dhp.first);
    Dart dh1(dhp.second);

    /* Already an edge? */
    if (dh0.t() == dh1.t()) {
      MESHC_LOG("Segment already an edge. Darts: "
		<< dhp);
      dh0.alpha0();
      if (v1 == dh0.v()) {
	dh0.alpha0();
	if (is_boundary)
	  boundary_.insert(dh0, meta);
	else
	  interior_.insert(dh0, meta);
	return dh0;
      } else {
	dh1.orbit1();
	if (is_boundary)
	  boundary_.insert(dh1, meta);
	else
	  interior_.insert(dh1, meta);
	return dh1;
      }
    }

    /* Check if the new segment needs to be split at a vertex or
       crossing segment */
    int v2 = CDTSplitSegment(dhp,trace);
    Dart dreturn;
    if (v2 >= 0) {
      /* Recursively insert the split segment pieces, returning a dart
	 originating at the first point. */
      MESHC_LOG("Split segment at vertex " << v2 << endl);
      MESHC_LOG("InsertSegment (" << v2 << "," << v1 << ")" << endl);
      dreturn = CDTInsertSegment(v2, v1, triangles, is_boundary, meta);
      MESHC_LOG("Second part inserted: #1 = " << dreturn << endl);
      MESHC_LOG("Target segment = (" << v2 << "," << v1 << ")" << endl);
      Dart dtmp(dreturn);
      dtmp.orbit1();
      MESHC_LOG("Actual segment = (" << dreturn.v() << "," << dtmp.v() << ")" << endl);
      LOP(triangles);
      triangles.clear();

      MESHC_LOG("InsertSegment (" << v0 << "," << v2 << ")" << endl);
      dreturn = CDTInsertSegment(v0, v2, triangles, is_boundary, meta);
      MESHC_LOG("First part inserted: #1 = " << dreturn << endl);
      MESHC_LOG("Target segment = (" << v0 << "," << v2 << ")" << endl);
      dtmp = dreturn;
      dtmp.orbit1();
      MESHC_LOG("Actual segment = (" << dreturn.v() << "," << dtmp.v() << ")" << endl);
      return dreturn;
    } else {
      /* No crossing segments or points-on-line, so insert the new
	 segment. */
      MESHC_LOG("No splitting required." << endl);
      dreturn = CDTInsertSegment(dhp, trace, triangles, is_boundary, meta);
      MESHC_LOG("One segment inserted: " << dreturn);
      return dreturn;
    }
  }








  Dart MeshC::CDTSegment(const bool boundary,
			 const constrT& constraint)
  {
    if (!prepareCDT()) return Dart();
    const int& v0 = constraint.first.first;
    const int& v1 = constraint.first.second;
    const constrMetaT& meta = constraint.second;

    //      MESHC_LOG_("Trying to add segment: "
    //		 << ci->first.first << "," << ci->first.second << endl);

    if (M_->useVT()) { /* Can check if the vertices are present, and
			  try to add them otherwise. */
      MESHC_LOG("CDT: Checking vertex " << v0 << endl);
      Dart dh = Dart(*M_,0);
      if (M_->VT(v0)==-1) {
	dh = insertNode(v0,dh);
	if (dh.isnull()) {
	  MESHC_LOG("CDT: Failed to insert node " << v0 << endl << *this);
	  return dh;
	}	
      }
      MESHC_LOG("CDT: Checked" << endl);
      MESHC_LOG("CDT: Checking vertex " << v1 << endl);
      if (M_->VT(v1)==-1) {
	dh = insertNode(v1,dh);
	if (dh.isnull()) {
	  MESHC_LOG("CDT: Failed to insert node " << v1 << endl << *this);
	  return dh;
	}
      }
      MESHC_LOG("CDT: Checked" << endl);
    }

    triangleSetT triangles;
    Dart ds(CDTInsertSegment(v0,v1,triangles,boundary,meta));
    if (ds.isnull()) {
      MESHC_LOG((boundary ? "Boundary" : "Interior")
		<< " segment not inserted ("
		<< v0 << "," << v1 << ")" << endl);
      return ds;
    }
    LOP(triangles);
    MESHC_LOG((boundary ? "Boundary" : "Interior")
	      << " segment inserted "
	      << ds << endl);
    return ds;
  }


  bool MeshC::buildCDT()
  {
    if (!prepareCDT()) return false;

    /* Make sure to set useVT true, so that we can easily detect
       missing vertices and try to add them. */
    bool M_useVT = M_->useVT();
    M_->useVT(true);

    constrListT::iterator ci_next;
    for (constrListT::iterator ci = constr_boundary_.begin();
	 ci != constr_boundary_.end(); ) {
      MESHC_LOG("Trying to add boundary segment: "
		 << ci->first.first << "," << ci->first.second <<
		 " group=" << ci->second << endl);
      if (!CDTSegment(true,*ci).isnull()) {
	MESHC_LOG("Success." << endl);
	ci_next = ci;
	ci_next++;
	ci = constr_boundary_.erase(ci);
	ci = ci_next;
      } else {
	MESHC_LOG("Failure." << endl);
	ci++;
      }
    }
    for (constrListT::iterator ci = constr_interior_.begin();
	 ci != constr_interior_.end(); ) {
      MESHC_LOG("Trying to add interior segment: "
		 << ci->first.first << "," << ci->first.second <<
		 " group=" << ci->second << endl);
      if (!CDTSegment(false,*ci).isnull()) {
	MESHC_LOG("Success." << endl);
	ci_next = ci;
	ci_next++;
	ci = constr_interior_.erase(ci);
	ci = ci_next;
      } else {
	MESHC_LOG("Failure." << endl);
	ci++;
      }
    }

    M_->useVT(M_useVT);

    MESHC_LOG("Boundary segments after CDT:" << endl << boundary_);
    MESHC_LOG("Interior segments after CDT:" << endl << interior_);

    MESHC_LOG("CDT finished" << endl << *this);
#ifndef FMESHER_NO_X
    M_->redrawX11("CDT finished");
#endif

    return (constr_boundary_.empty() && constr_interior_.empty());
  }









  bool MeshC::buildRCDTlookahead(MCQsegm* segm, const Point& c)
  {
    MESHC_LOG("Checking for potentially encroached segments at ("
	      << c[0] << ',' << c[1] << ',' << c[2] << ")" << endl);
    for (MCQ::const_iterator ci = segm->begin();
	 ci != segm->end(); ci++) {
      Dart dhc(ci->first);
      double encr = M_->edgeEncroached(dhc,c);
      if (encr>0.0) {
	MESHC_LOG("Potentially encroached segment: "
		  << dhc << " "
		  << encr << endl);
	bisectEdgeDelaunay(dhc);
	// if (!bisectEdgeDelaunay(dhc).isnull())
	//  xtmpl_press_ret("Potentialy encroached segment split.");
	//else
	//  xtmpl_press_ret("Failed to split potentialy encroached segment.");
	return false;
      }
    }
    return true;
  }



  bool MeshC::buildRCDT()
  {
    if (state_<State_RCDT)
      return false; /* ERROR: RCDT not initialised. */

    MESHC_LOG("Encroached boundary segments before RCDT:" << endl
	      << boundary_);
    MESHC_LOG("Encroached interior segments before RCDT:" << endl
	      << interior_);
    MESHC_LOG("Skinny triangles before RCDT:" << endl << skinny_);
    MESHC_LOG("Big triangles before RCDT:" << endl << big_);

    Dart dh;

    while (!(boundary_.emptyQ() && interior_.emptyQ() &&
	     skinny_.emptyQ() && big_.emptyQ())) {

      MESHC_LOG("RCDT: (Bo,In,Sk,Bi) = ("
		<< boundary_.countQ() << ","
		<< interior_.countQ() << ","
		<< skinny_.countQ() << ","
		<< big_.countQ() << ")" << endl);

      dh = boundary_.quality();
      if (!dh.isnull()) {
	MESHC_LOG("Encroached boundary segment: "
		  << dh << " "
		  << big_.quality(dh) << endl);
	bisectEdgeDelaunay(dh);
	//if (!bisectEdgeDelaunay(dh).isnull())
	//  xtmpl_press_ret("Boundary segment has been split");
	//else
	//  xtmpl_press_ret("Boundary segment split failed");
	continue;
      }
      
      dh = interior_.quality();
      if (!dh.isnull()) {
	MESHC_LOG("Encroached interior segment: "
		  << dh << " "
		  << big_.quality(dh) << endl);
	bisectEdgeDelaunay(dh);
	//if (!bisectEdgeDelaunay(dh).isnull())
	//  xtmpl_press_ret("Interior segment has been split");
	//else
	//  xtmpl_press_ret("Interior segment split failed");
	continue;
      }

      if ((max_n0_ >= 0) & (max_n0_ <= int(M_->nV()))) {
	MESHC_LOG("Max vertex count reached: max_n0 = "
		  << max_n0_ << " <= nV = "
		  << M_->nV() << endl);
	skinny_.clear();
	big_.clear();
	continue;
      }

      //      xtmpl_press_ret("No segments need splitting.");

      dh = skinny_.quality();
      if (!dh.isnull()) {
	MESHC_LOG("Skinny triangle: "
		  << dh << " "
		  << skinny_.quality(dh) << endl);
	Point c;
	calcSteinerPoint(dh,c);
	if ((!buildRCDTlookahead(&boundary_,c)) ||
	    (!buildRCDTlookahead(&interior_,c)))
	  continue;
	if (insertNode(addVertex(c),dh).isnull()) {
	  MESHC_LOG("Skinny triangle elimination failed" << endl);
	  M_->removeLastVertex(); /* Failed to add to graph, so delete it. */
	  //	  skinny_.erase(dh); /* Hope we don't fall into
	  //	                        infinite loop; we do... */
	  skinny_.clear(); /* Stop caring about the skinny triangles. */
	  //	  return false;
	  continue;
	}
	continue;
      }
      
      if ((max_n1_ >= 0) & (max_n1_ <= int(M_->nV()))) {
	MESHC_LOG("Max vertex count reached: max_n1 = "
		  << max_n1_ << " <= nV = "
		  << M_->nV() << endl);
	big_.clear();
	continue;
      }

      dh = big_.quality();
      if (!dh.isnull()) {
	MESHC_LOG("Big triangle: "
		  << dh << " "
		  << big_.quality(dh) << endl);
	Point c;
	calcSteinerPoint(dh,c);
	if ((!buildRCDTlookahead(&boundary_,c)) ||
	    (!buildRCDTlookahead(&interior_,c)))
	  continue;
	if (insertNode(addVertex(c),dh).isnull()) {
	  MESHC_LOG("Big triangle elimination failed failed" << endl);
	  M_->removeLastVertex(); /* Failed to add to graph, so delete it. */
	  skinny_.erase(dh); /* Hope we don't fall into infinite loop. */
	  //	  return false;
	  continue;
	}
	continue;
      }
      
    }

    MESHC_LOG("RCDT finished" << endl << *this);
#ifndef FMESHER_NO_X
    M_->redrawX11("RCDT finished");
#endif

    return true;
  }

  bool MeshC::RCDT(double angle_limit,
		   double big_limit,
		   const double* big_limits,
		   size_t nQL,
		   int max_n0,
		   int max_n1)
  {
    if (!prepareRCDT(1./std::sin(M_PI/180.*angle_limit)/2.,
		     big_limit,big_limits,nQL,max_n0,max_n1)) return false;
    return buildRCDT();
  }


  /*!
    Flood-fill algorithm for exterior triangle removal:
    \verbatim
    1. For each boundary segment
    2.   Add a possible exterior triangle to the "exterior set" Ext.
    3.   Unlink the edge.
    4. While Ext is not empty, let t in Ext
    5.   Add any triangle linked to t to Ext.
    6.   Remove t.
    \endverbatim

    Step 1. is delicate, since boundary segments will be added and
    removed in step 3.  Solved by finding the previous segment again,
    and then moving to the next.

    TODO: Allow optional vertex reordering.
   */
  bool MeshC::PruneExterior()
  {
    if (state_ < State_CDT) {
      /* Since there are no constraints at this state, no exterior
	 needs to be pruned, but go to state State_CDT anyway. */
      prepareCDT();
      is_pruned_ = true;
      return true;
    }
    is_pruned_ = true;
    
    Dart d0, dh;
    triangleSetT ext;

    /* Unlink the exterior. */
    for (MCQsegm::const_iterator boundary_segm_i = boundary_.begin();
	 boundary_segm_i != boundary_.end();
	 boundary_segm_i = ++boundary_.find(d0)) {
      d0 = boundary_segm_i->first;
      if (!d0.onBoundary()) {
	dh = d0;
	dh.orbit1();
	if (!boundary_.found(dh))
	  ext.insert(dh.t());
	unlinkEdge(d0);
      }
    }

    /* Make sure useVT is off: */
    bool M_useVT = M_->useVT();
    M_->useVT(false);

    /* Remove exterior triangles. */
    int t_relocated;
    triangleSetT::iterator ext_j;
    for (triangleSetT::iterator ext_i = ext.begin();
	 ext_i != ext.end();
	 ext_i = ext.begin()) {
      dh = Dart(*M_,*ext_i);
      if (!dh.onBoundary()) ext.insert(dh.tadj());
      dh.orbit2();
      if (!dh.onBoundary()) ext.insert(dh.tadj());
      dh.orbit2();
      if (!dh.onBoundary()) ext.insert(dh.tadj());
      t_relocated = removeTriangle(dh);
      if ((ext_j = ext.find(t_relocated)) != ext.end())
	ext.erase(ext_j);
      else
	ext.erase(ext_i);
    }

    { /* Remove unused points. */
      /* Make sure useVT is on: */
      M_->useVT(true);
      int v = M_->nV()-1;
      while ((v>=0) && (M_->VT(v)==-1)) {
	M_->removeLastVertex();
	--v;
      }
    }

    /* Restore useVT-state: */
    M_->useVT(M_useVT);

    MESHC_LOG("Boundary segments after pruning: "
	      << boundary_);

    MESHC_LOG("PruneExterior finished" << endl << *this);
#ifndef FMESHER_NO_X
    M_->redrawX11("PruneExterior finished");
#endif

    return true;
  }






  int MeshC::addVertex(const Point& s)
  {
    int nVorig = (int)M_->nV();
    M_->S_append(s);
    if ((state_>=State_RCDT) && big_.usingQv()) {
      big_.setQv(nVorig,big_.getQ());
    }
    return M_->nV()-1;
  }

  int MeshC::addVertices(const Matrix3double& S)
  {
    size_t nVorig = M_->nV();
    M_->S_append(S);
    if ((state_>=State_RCDT) && big_.usingQv()) {
      for (size_t v=nVorig; v < nVorig+S.rows(); v++)
	big_.setQv(v,big_.getQ());
    }
    return M_->nV()-S.rows();
  }


  Dart MeshC::swapEdge(const Dart& d, MCQswapable& swapable)
  {
    if (!swapable.swapable(d)) {
      /* Not allowed to swap. */
      MESHC_LOG("Not allowed to swap dart " << endl
		<< d << endl);
      return d;
    }
    
    /* Collect swapable data */
    bool edge_list[4];
    Dart dh(d);
    swapable.erase(dh);
    dh.orbit2rev();
    if ((edge_list[1] = swapable.found(dh))) swapable.erase(dh);
    dh.orbit2rev();
    if ((edge_list[2] = swapable.found(dh))) swapable.erase(dh);
    dh.orbit0().orbit2rev();
    if ((edge_list[3] = swapable.found(dh))) swapable.erase(dh);
    dh.orbit2rev();
    if ((edge_list[0] = swapable.found(dh))) swapable.erase(dh);

    // MESHC_LOG("LOP edge list: ("
    // 	      << edge_list[0] << ","
    // 	      << edge_list[1] << ","
    // 	      << edge_list[2] << ","
    // 	      << edge_list[3] << ")" << endl);
    
    Dart dnew(swapEdge(d));
    if (dh == dnew) {
      /* ERROR: this should not happen. */
      MESHC_LOG("Edge swap appears to have failed!" << endl);
      return dnew;
    }

    /* Reassemble swapable data */
    dh = dnew;
    swapable.insert(dh);
    dh.orbit2();
    if (edge_list[1]) swapable.insert(dh);
    dh.orbit2();
    if (edge_list[0]) swapable.insert(dh);
    dh.orbit2().orbit0rev();
    if (edge_list[3]) swapable.insert(dh);
    dh.orbit2();
    if (edge_list[2]) swapable.insert(dh);

    //    MESHC_LOG("Edge swapped:" << endl);

    return dnew;
  }

  Dart MeshC::swapEdge(const Dart& d)
  {
    if (state_ < State_CDT) {
      return M_->swapEdge(d);
    }

    if (boundary_.segm(d) || interior_.segm(d)) {
      /* ERROR: Not allowed to swap. */
      MESHC_LOG("ERROR: Not allowed to swap dart "
		<< d << endl);
      return d;
    }

    /* Collect CDT data */
    bool segm_b[4] = {false,false,false,false};
    bool segm_i[4] = {false,false,false,false};
    constrMetaT meta_b[4] = {constrMetaT(),
			     constrMetaT(),
			     constrMetaT(),
			     constrMetaT()};
    constrMetaT meta_i[4] = {constrMetaT(),
			     constrMetaT(),
			     constrMetaT(),
			     constrMetaT()};
    Dart dh(d);
    if (state_>=State_CDT) {
      dh.orbit2rev();
      if ((segm_b[1] = boundary_.found(dh))) meta_b[1] = boundary_.erase(dh);
      if ((segm_i[1] = interior_.found(dh))) meta_i[1] = interior_.erase(dh);
      dh.orbit2rev();
      if ((segm_b[2] = boundary_.found(dh))) meta_b[2] = boundary_.erase(dh);
      if ((segm_i[2] = interior_.found(dh))) meta_i[2] = interior_.erase(dh);
      dh.orbit0().orbit2rev();
      if ((segm_b[3] = boundary_.found(dh))) meta_b[3] = boundary_.erase(dh);
      if ((segm_i[3] = interior_.found(dh))) meta_i[3] = interior_.erase(dh);
      dh.orbit2rev();
      if ((segm_b[0] = boundary_.found(dh))) meta_b[0] = boundary_.erase(dh);
      if ((segm_i[0] = interior_.found(dh))) meta_i[0] = interior_.erase(dh);
    }

    if (state_>=State_RCDT) {
      /* Collect RCDT data */
      dh = d;
      skinny_.erase(dh);
      big_.erase(dh);
      dh.orbit1();
      skinny_.erase(dh);
      big_.erase(dh);
    }

    Dart dnew(M_->swapEdge(d));

    if (state_>=State_CDT) {
      /* Reassemble CDT data */
      dh = dnew;
      dh.orbit2();
      boundary_.update(dh); if (segm_b[1]) boundary_.insert(dh,meta_b[1]);
      interior_.update(dh); if (segm_i[1]) interior_.insert(dh,meta_i[1]);
      dh.orbit2();
      boundary_.update(dh); if (segm_b[0]) boundary_.insert(dh,meta_b[0]);
      interior_.update(dh); if (segm_i[0]) interior_.insert(dh,meta_i[0]);
      dh.orbit2().orbit0rev();
      boundary_.update(dh); if (segm_b[3]) boundary_.insert(dh,meta_b[3]);
      interior_.update(dh); if (segm_i[3]) interior_.insert(dh,meta_i[3]);
      dh.orbit2();
      boundary_.update(dh); if (segm_b[2]) boundary_.insert(dh,meta_b[2]);
      interior_.update(dh); if (segm_i[2]) interior_.insert(dh,meta_i[2]);
    }
    
    if (state_>=State_RCDT) {
      /* Reassemble RCDT data */
      dh = dnew;
      skinny_.insert(dh);
      big_.insert(dh);
      dh.orbit1();
      skinny_.insert(dh);
      big_.insert(dh);
    }

    //    MESHC_LOG("Edge swapped, boundary segments:" << endl
    //	      << boundary_);

    return dnew;
  }

  Dart MeshC::splitEdge(const Dart& d, int v)
  {
    if (state_ < State_CDT) {
      return M_->splitEdge(d,v);
    }

    /* Collect CDT data */
    Dart dh(d);
    bool segm_b[6];
    bool segm_i[6];
    constrMetaT meta_b[6];
    constrMetaT meta_i[6];
    if (state_>=State_CDT) {
    for (int i=0;i<3;i++) {
      if ((segm_b[i] = boundary_.found(dh))) meta_b[i] = boundary_.erase(dh);
      if ((segm_i[i] = interior_.found(dh))) meta_i[i] = interior_.erase(dh);
      dh.orbit2();
    }
    if (!dh.onBoundary()) {
      dh.orbit1();
      for (int i=3;i<6;i++) {
	if ((segm_b[i] = boundary_.found(dh))) meta_b[i] = boundary_.erase(dh);
	if ((segm_i[i] = interior_.found(dh))) meta_i[i] = interior_.erase(dh);
	dh.orbit2();
      }
    }
    }

    if (state_>=State_RCDT) {
      /* Collect RCDT data */
      dh = d;
      skinny_.erase(dh);
      big_.erase(dh);
      if (!dh.onBoundary()) {
	dh.orbit1();
	skinny_.erase(dh);
	big_.erase(dh);
      }
    }

    Dart dnew(M_->splitEdge(d,v));

    if (state_>=State_CDT) {
      /* Reassemble CDT data */
      dh = dnew;
      boundary_.update(dh); if (segm_b[0]) boundary_.insert(dh,meta_b[0]);
      interior_.update(dh); if (segm_i[0]) interior_.insert(dh,meta_i[0]);
      dh.orbit2();
      boundary_.update(dh); if (segm_b[1]) boundary_.insert(dh,meta_b[1]);
      interior_.update(dh); if (segm_i[1]) interior_.insert(dh,meta_i[1]);
      dh.orbit2().orbit0rev();
      boundary_.update(dh); if (segm_b[2]) boundary_.insert(dh,meta_b[2]);
      interior_.update(dh); if (segm_i[2]) interior_.insert(dh,meta_i[2]);
      dh.orbit2();
      boundary_.update(dh); if (segm_b[0]) boundary_.insert(dh,meta_b[0]);
      interior_.update(dh); if (segm_i[0]) interior_.insert(dh,meta_i[0]);
      if (!dh.onBoundary()) {
	dh.orbit1();
	boundary_.update(dh); if (segm_b[3]) boundary_.insert(dh,meta_b[3]);
	interior_.update(dh); if (segm_i[3]) interior_.insert(dh,meta_i[3]);
	dh.orbit2();
	boundary_.update(dh); if (segm_b[4]) boundary_.insert(dh,meta_b[4]);
	interior_.update(dh); if (segm_i[4]) interior_.insert(dh,meta_i[4]);
	dh.orbit2().orbit0rev();
	boundary_.update(dh); if (segm_b[5]) boundary_.insert(dh,meta_b[5]);
	interior_.update(dh); if (segm_i[5]) interior_.insert(dh,meta_i[5]);
	dh.orbit2();
	boundary_.update(dh); if (segm_b[3]) boundary_.insert(dh,meta_b[3]);
	interior_.update(dh); if (segm_i[3]) interior_.insert(dh,meta_i[3]);
      }
    }

    if (state_>=State_RCDT) {
      /* Reassemble RCDT data */
      dh = dnew;
      skinny_.insert(dh);
      big_.insert(dh);
      dh.orbit0();
      skinny_.insert(dh);
      big_.insert(dh);
      if (!dnew.onBoundary()) {
	dh.orbit0();
	skinny_.insert(dh);
	big_.insert(dh);
	dh.orbit0();
	skinny_.insert(dh);
	big_.insert(dh);
      }
    }

    //    MESHC_LOG("Edge split, boundary segments:" << endl
    //	      << boundary_);
    //    MESHC_LOG("Edge split, interior segments:" << endl
    //	      << boundary_);

    return dnew;
  }

  Dart MeshC::splitTriangle(const Dart& d, int v)
  {
    if (state_ < State_CDT) {
      return M_->splitTriangle(d,v);
    }

    /* Collect CDT data */
    Dart dh(d);
    bool segm_b[3];
    bool segm_i[3];
    constrMetaT meta_b[3];
    constrMetaT meta_i[3];
    if (state_>=State_CDT) {
      for (int i=0;i<3;i++) {
	if ((segm_b[i] = boundary_.found(dh))) meta_b[i] = boundary_.erase(dh);
	if ((segm_i[i] = interior_.found(dh))) meta_i[i] = interior_.erase(dh);
	dh.orbit2();
      }
    }

    if (state_>=State_RCDT) {
      /* Collect RCDT data */
      skinny_.erase(d);
      big_.erase(d);
    }

    Dart dnew(M_->splitTriangle(d,v));

    /* Reassebmle CDT data */
    if (state_>=State_CDT) {
      dh = dnew;
      for (int i=0;i<3;i++) {
	dh.orbit2();
	boundary_.update(dh); if (segm_b[i]) boundary_.insert(dh,meta_b[i]);
	interior_.update(dh); if (segm_i[i]) interior_.insert(dh,meta_i[i]);
	dh.orbit2rev().orbit0();
      }
    }

    if (state_>=State_RCDT) {
      /* Reassemble RCDT data */
      dh = dnew;
      skinny_.insert(dh);
      big_.insert(dh);
      dh.orbit0();
      skinny_.insert(dh);
      big_.insert(dh);
      dh.orbit0();
      skinny_.insert(dh);
      big_.insert(dh);
    }

    //    MESHC_LOG("Triangle split, boundary segments:" << endl
    //	      << boundary_);

    return dnew;
  }




  /*!  */
  void MeshC::unlinkEdge(Dart& d)
  {
    if (state_<State_CDT) {
      M_->unlinkEdge(d);
      return;
    }

    Dart dh(d);
    bool onboundary = d.onBoundary();
    if (!onboundary) {
      dh.orbit1();
      if (interior_.found(dh)) interior_.erase(dh);
    }
    if (interior_.found(d)) interior_.erase(d);
    
    M_->unlinkEdge(d);
        
    if (!onboundary) {
      boundary_.insert(dh,boundary_.erase(dh));
    }
    boundary_.insert(d,boundary_.erase(d));

    return;
  }

  /*!  */
  int MeshC::removeTriangle(Dart& d)
  {
    if (state_ < State_CDT) {
      return M_->removeTriangle(d.t());
    }

    Dart dh(d);
    interior_.erase(dh);
    boundary_.erase(dh);
    if (!dh.onBoundary()) {
      dh.orbit1();
      boundary_.insert(dh,interior_.erase(dh));
      dh.orbit1();
    }
    dh.orbit2();
    interior_.erase(dh);
    boundary_.erase(dh);
    if (!dh.onBoundary()) {
      dh.orbit1();
      boundary_.insert(dh,interior_.erase(dh));
      dh.orbit1();
    }
    dh.orbit2();
    interior_.erase(dh);
    boundary_.erase(dh);
    if (!dh.onBoundary()) {
      dh.orbit1();
      boundary_.insert(dh,interior_.erase(dh));
      dh.orbit1();
    }

    int t_removed = d.t();
    int t_relocated = M_->removeTriangle(t_removed);
    
    dh = Dart(*M_,t_removed,1,0);
    Dart dh_old = Dart(*M_,t_relocated,1,0);
    if (boundary_.found(dh_old)) {
      boundary_.insert(dh,boundary_.erase(dh_old));
    }
    if (interior_.found(dh_old)) {
      interior_.insert(dh,interior_.erase(dh_old));
    }
    dh.orbit2(); dh_old.orbit2();
    if (boundary_.found(dh_old)) {
      boundary_.insert(dh,boundary_.erase(dh_old));
    }
    if (interior_.found(dh_old)) {
      interior_.insert(dh,interior_.erase(dh_old));
    }
    dh.orbit2(); dh_old.orbit2();
    if (boundary_.found(dh_old)) {
      boundary_.insert(dh,boundary_.erase(dh_old));
    }
    if (interior_.found(dh_old)) {
      interior_.insert(dh,interior_.erase(dh_old));
    }
    
    return t_relocated;
  }







  bool MeshC::isSegment(const Dart& d) const
  {
    if (state_<State_CDT) /* No segments */
      return false;

    return (boundary_.segm(d) || interior_.segm(d));
  }







  std::ostream& operator<<(std::ostream& output, const MCQ& Q)
  {
    if (Q.empty()) return output;
    output << "N,n = " << Q.count() << "," << Q.countQ() << endl;
    for (MCQ::map_type::const_iterator qi = Q.darts_.begin();
	 qi != Q.darts_.end(); qi++) {
      output << ' ' << qi->first
	     << ' ' << std::scientific << qi->second
	     << ' ' << Q.foundQ(qi->first)
	     << endl;
    }
    return output;
  }

  std::ostream& operator<<(std::ostream& output, const DartPair& dp)
  {
    output << "d0=(" << dp.first << ") d1=("
	   << dp.second << ")" << endl;
    return output;
  }

  std::ostream& operator<<(std::ostream& output, const DartList& ds)
  {
    output << "n = " << ds.size() << endl;
    if (ds.empty()) return output;
    for (DartList::const_iterator di = ds.begin();
	 di != ds.end(); di++) {
      output << ' ' << *di << endl;
    }
    return output;
  }

  std::ostream& operator<<(std::ostream& output, const std::set<int>& il)
  {
    output << "(n = " << il.size() << ")";
    if (il.empty()) return output;
    for (std::set<int>::const_iterator i = il.begin();
	 i != il.end(); i++) {
      output << ' ' << *i;
    }
    return output;
  }

  std::ostream& operator<<(std::ostream& output, const std::list<int>& il)
  {
    output << "(n = " << il.size() << ")";
    if (il.empty()) return output;
    for (std::list<int>::const_iterator i = il.begin();
	 i != il.end(); i++) {
      output << ' ' << *i;
    }
    return output;
  }

  std::ostream& operator<<(std::ostream& output,
			   const IntPair& p)
  {
    output << "(" << p.first << "," << p.second << ")";
    return output;
  }

  std::ostream& operator<<(std::ostream& output,
			   const std::list<IntPair>& il)
  {
    output << "(n = " << il.size() << ")";
    if (il.empty()) return output;
    for (std::list<IntPair>::const_iterator i = il.begin();
	 i != il.end(); i++) {
      output << " " << *i;
    }
    return output;
  }


  std::ostream& operator<<(std::ostream& output,
			   const MCQsegm::meta_map_type& il)
  {
    output << "(n = " << il.size() << ")" << endl;
    if (il.empty()) return output;
    for (MCQsegm::meta_map_type::const_iterator qi = il.begin();
	 qi != il.end(); qi++) {
      output <<' ' << qi->first
	     << ' ' << qi->second
	     << endl;
    }
    return output;
  }




  std::ostream& operator<<(std::ostream& output,
			   const MCQsegm& segm)
  {
    output << "Segments:\t" << segm.count();
    if (segm.countQ() > 0)
      output << "(" << segm.countQ() << " encroached)";
    output << endl;

    output << "Darts+quality:" << endl << segm.darts_ << endl;
    output << "Metadata:" << endl << segm.meta_;

    return output;
  }


  std::ostream& operator<<(std::ostream& output,
			   const MeshC& MC)
  {
    output << *MC.M_;
    output << "Construction state:\t"
	   << (MC.state_==MeshC::State_noT ? "No triangles" :
	       (MC.state_==MeshC::State_CET ? "CET (Convex enclosure triangulation)" :
		(MC.state_==MeshC::State_DT ? "DT (Delaunay triangulation)" :
		 (MC.state_==MeshC::State_CDT ? "CDT (Constrained DT)" :
		  (MC.state_==MeshC::State_RCDT ? "RCDT (Refined CDT)" : ""
		   )))))
	   << (MC.is_pruned_ ? ", exterior pruned" : "") << endl;
    if (MC.state_>=MeshC::State_CDT) {
      output << "Boundary " << MC.boundary_;
      output << "Interior " << MC.interior_;
    }
    if (MC.state_>=MeshC::State_RCDT) {
      output << "Skinny triangles:\t" << MC.skinny_.countQ() << endl;
      output << "Big triangles:\t\t" << MC.big_.countQ() << endl;
    }
    return output;
  }


} /* namespace fmesh */
