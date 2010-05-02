#include <cstddef>
#include <cstring>
#include <set>
#include <map>
#include <sstream>
#include <cmath>

#include "predicates.h"
#include "meshc.h"

#define WHEREAMI __FILE__ << "(" << __LINE__ << ")\t"

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

  const double MCQ::quality(const Dart& d) const
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







  MCQtri::MCQtri(MeshC* MC, bool only_quality, double quality_limit)
    : MCQ(MC, only_quality), quality_limit_(quality_limit)
  {
    setQ(quality_limit);
  }

  void MCQtri::setQ(double quality_limit)
  {
    quality_limit_ = quality_limit;
    if (!empty()) {
      /* TODO: Implement updating the sets. */
      NOT_IMPLEMENTED;
    }
  }

  double MCQtri::calcQ(const Dart& d) const {
    double quality_lim_ = quality_limit_; /* TODO: min_vi ql[t][vi] */
    return (calcQtri(d) - quality_lim_);
  };

  double MCQskinny::calcQtri(const Dart& d) const
  {
    double quality_ = MC_->skinnyQuality(d.t());
    return quality_;
  }
  
  double MCQbig::calcQtri(const Dart& d) const
  {
    double quality_ = MC_->bigQuality(d.t());
    return quality_;
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

  double MCQswapable::calcQ(const Dart& d) const
  {
    if (MC_->M_->circumcircleOK(d))
      return -1.0; /* Not swapable. */
    else
      return 1.0; /* Swapable. */
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

  const double MCQswapable::quality(const Dart& d) const
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
    
    double encr = M_->edgeEncroached(d,M_->S()[dh.v()]);

    dh.orbit2rev();
    std::cout << WHEREAMI << "encroachedQ("
	      << d.v() << "," << dh.v()
	      << ") = " << encr << std::endl;

    return encr;
  }

  double MeshC::skinnyQuality(int t) const
  {
    if ((t<0) || (t>=(int)M_->nT())) return 0.0;

    double skinny = (M_->triangleCircumcircleRadius(t) / 
		     M_->triangleShortestEdge(t));

    //    std::cout << WHEREAMI << "skinnyQ(" << t << ") = " << skinny << std::endl;

    return skinny;
  }

  double MeshC::bigQuality(int t) const
  {
    return M_->triangleLongestEdge(t);
    //return M_->triangleArea(t);
  }





  /*! Alg 4.3 */
  bool MeshC::recSwapDelaunay(const Dart& d0)
  {
    Dart d1, d2;

    std::cout << WHEREAMI << "Trying to swap " << d0 << std::endl;

    if (d0.isnull() or d0.onBoundary())
      return true; /* OK. Not allowed to swap. */
    if (isSegment(d0))
      return true ; /* OK. Not allowed to swap. */
    if (M_->circumcircleOK(d0))
      return true; /* OK. Need not swap. */

    std::cout << WHEREAMI << "Swap " << d0 << std::endl;

    /* Get opposing darts. */
    d1 = d0;
    d1.alpha1();
    if (d1.onBoundary()) d1 = Dart(); else d1.alpha2();
    d2 = d0;
    d2.orbit2rev().alpha1(); 
    if (d2.onBoundary()) d2 = Dart(); else d2.alpha2();
    
    //    std::cout << WHEREAMI << "TVpre  = " << std::endl << M_->TVO();
    swapEdge(d0);
    //    std::cout << WHEREAMI << "TVpost = " << std::endl << M_->TVO();
    //    std::cout << WHEREAMI << "TTpost = " << std::endl << M_->TTO();

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

    //    std::cout << WHEREAMI << "TV = " << std::endl << M_->TVO();
    std::cout << WHEREAMI << "Split triangle " << td << " with vertex " << v << std::endl;
    d = splitTriangle(td,v);
    
    if (!d0.isnull()) recSwapDelaunay(d0);
    if (!d1.isnull()) recSwapDelaunay(d1);
    if (!d2.isnull()) recSwapDelaunay(d2);

    //    std::cout << WHEREAMI << "TV = " << std::endl << M_->TVO();
    
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

    //    std::cout << WHEREAMI << "TV = " << std::endl << M_->TVO();
    std::cout << WHEREAMI << "Split edge " << ed << " with vertex " << v << std::endl;
    d = splitEdge(ed,v);
    
    if (!d0.isnull()) recSwapDelaunay(d0);
    if (!d1.isnull()) recSwapDelaunay(d1);
    if (!d2.isnull()) recSwapDelaunay(d2);
    if (!d3.isnull()) recSwapDelaunay(d3);

    //    std::cout << WHEREAMI << "TV = " << std::endl << M_->TVO();
    
    return d;
  }

  Dart MeshC::bisectEdgeDelaunay(const Dart& d)
  {
    Dart dh(d);
    int v0(dh.v());
    dh.alpha0();
    int v1(dh.v());

    double s[3];
    Vec::sum(s,M_->S()[v0],M_->S()[v1]);
    Vec::rescale(s,0.5);

    switch (M_->type()) {
    case Mesh::Mtype_manifold:
      /* TODO: Implement. */
      NOT_IMPLEMENTED;
      /* break; For now, fall through to Mtype_plane behaviour! */
    case Mesh::Mtype_plane:
      /* Nothing to do! */
      break;
    case Mesh::Mtype_sphere:
      Vec::rescale(s,1./Vec::length(s));
      break;
    }

    return splitEdgeDelaunay(d,addVertices(&s,1));
  }



  bool MeshC::killTriangle(const Dart& d)
  {
    Point c;
    M_->triangleCircumcenter(d.t(),c);

    std::cout << WHEREAMI << "Center: ("
	      << c[0] << ","
	      << c[1] << ","
	      << c[2] << ")" << std::endl;

    return insertNode(addVertices(&c,1),d);
  }



  /*! Alg 9.3 */
  bool MeshC::insertNode(int v, const Dart& ed)
  {
    Dart td;

    std::cout << WHEREAMI << "Locating node " << v
	      << " " << M_->S()[v] << std::endl;

    td = M_->locatePoint(ed,M_->S()[v]);
    if (td.isnull()) { return false; }; /* ERROR, not found! */
    td = Dart(*M_,td.t());
    Point bary;
    M_->barycentric(td,M_->S()[v],bary);
    size_t pattern(size_t(bary[0]>MESH_EPSILON)*1+
		   size_t(bary[1]>MESH_EPSILON)*2+
		   size_t(bary[2]>MESH_EPSILON)*4);
    std::cout << WHEREAMI << "Triangle dart " << td
	      << " bary=" << bary
	      << " pattern=" << pattern << std::endl;
    switch (pattern) {
    case 7: // +++
      splitTriangleDelaunay(td,v);
      break;
    case 6: // -++ Split e0
      td.orbit2();
      std::cout << WHEREAMI << "Edge dart " << td << std::endl;
      splitEdgeDelaunay(td,v);
      break;
    case 5: // +-+ Split e1
      td.orbit2rev();
      std::cout << WHEREAMI << "Edge dart " << td << std::endl;
      splitEdgeDelaunay(td,v);
      break;
    case 3: // ++- Split e2
      std::cout << WHEREAMI << "Edge dart " << td << std::endl;
      splitEdgeDelaunay(td,v);
      break;
    case 1: // +-- Close to node 0, not allowed
    case 2: // -+- Close to node 1, not allowed
    case 4: // --+ Close to node 2, not allowed
      return false;
      break;
    case 0: // --- Close to all nodes, should not happen!
      return false;
      break;
    }

    return true;
  }

  bool MeshC::DT(const vertexListT& v_set)
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

    if (state_>=State_CDT)
      std::cout << WHEREAMI << "Boundary segments before DT:" << std::endl << boundary_;

    int v;
    vertexListT::const_iterator v_iter;
    Dart td, d, d0, d1, d2;

    for (v_iter = v_set.begin(); v_iter != v_set.end(); v_iter++) {
      v = *v_iter;
      insertNode(v,Dart(*M_,0)); /* TODO: More clever starting edge? */

      if (state_>=State_CDT)
	std::cout << WHEREAMI << "Boundary segments after DT:"
		  << std::endl << boundary_;
    }
      
    M_->redrawX11("DT finished");

    state_ = State_DT;
    return true;
  }


  bool MeshC::prepareDT()
  {
    if (state_<State_DT) {
      /* We need to build a DT first. */
      triangleSetT t_set;
      for (int t=0;t<(int)M_->nT();t++)
	t_set.insert(t);
      if (LOP(t_set))
	state_ = State_DT;
    }
    return (state_>=State_DT) && (!is_pruned_);
  }


  bool MeshC::prepareCDT()
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

  bool MeshC::prepareRCDT(double skinny_limit, double big_limit)
  {
    if (!prepareCDT()) return false; /* Make sure we have a CDT. */

    skinny_.clear();
    big_.clear();
    skinny_.setQ(skinny_limit);
    big_.setQ(big_limit);

    for (int t=0;t<(int)M_->nT();t++) {
      skinny_.insert(Dart(*M_,t));
      big_.insert(Dart(*M_,t));
    }

    state_ = State_RCDT;
    return true;
  }



  bool MeshC::CDTBoundary(const constrListT& constr)
  {
    if (!prepareCDT()) return false;

    constr_boundary_ = constrListT(constr.begin(),constr.end());
    
    return buildCDT();
  };

  bool MeshC::CDTInterior(const constrListT& constr)
  {
    if (!prepareCDT()) return false;

    constr_interior_ = constrListT(constr.begin(),constr.end());
    
    return buildCDT();
  };






  bool MeshC::LOP(MCQswapable& swapable)
  {
    std::cout << WHEREAMI << "LOP swapable: "
	      << swapable.countQ() << "/" << swapable.count() << std::endl;
    /* Swap edges, until none are swapable. */
    while (!swapable.emptyQ()) {
      swapEdgeLOP(swapable.beginQ()->d_,swapable);
      std::cout << WHEREAMI << "LOP swapable: "
		<< swapable.countQ() << "/" << swapable.count() << std::endl;
    }

    M_->redrawX11("LOP finished");

    return true;
  }

  bool MeshC::LOP(const triangleSetT& t_set)
  {
    /* Locate interior edges */
    Dart dh, dh2;
    MCQswapable swapable(this);
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
		|| (!boundary_.segm(dh))
		|| (!interior_.segm(dh)))) /* Don't add CDT segments. */
	  swapable.insert(dh); /* MCQswapable takes care of duplicates. */
	dh.orbit2();
      }
    }

    return LOP(swapable);
  }


  bool MeshC::buildCDT()
  {
    if (!prepareCDT()) return false;

    /* TODO: Implement. */

    constrListT::iterator ci_next;
    for (constrListT::iterator ci = constr_boundary_.begin();
	 ci != constr_boundary_.end(); ) {
      NOT_IMPLEMENTED;
      if (true) {
	ci_next = ci;
	ci_next++;
	ci = constr_boundary_.erase(ci);
	ci = ci_next;
      } else
	ci++;
    }
    for (constrListT::iterator ci = constr_interior_.begin();
	 ci != constr_interior_.end(); ) {
      NOT_IMPLEMENTED;
      if (true) {
	ci_next = ci;
	ci_next++;
	ci = constr_interior_.erase(ci);
	ci = ci_next;
      } else
	ci++;
    }

    std::cout << WHEREAMI << "Boundary segments after CDT:" << std::endl << boundary_;
    std::cout << WHEREAMI << "Interior segments after CDT:" << std::endl << interior_;

    return (constr_boundary_.empty() && constr_interior_.empty());
  };









  bool MeshC::buildRCDTlookahead(MCQsegm* segm, const Point& c)
  {
    std::cout << WHEREAMI << "Checking for potentially encroached segments at ("
	      << c[0] << ',' << c[1] << ',' << c[2] << ")" << std::endl;
    for (MCQ::const_iterator ci = segm->begin();
	 ci != segm->end(); ci++) {
      Dart dhc(ci->first);
      double encr = M_->edgeEncroached(dhc,c);
      if (encr>0.0) {
	std::cout << WHEREAMI << "Potentially encroached segment: "
		  << dhc << " "
		  << encr << std::endl;
	bisectEdgeDelaunay(dhc);
	return false;
      }
    }
    return true;
  }



  bool MeshC::buildRCDT()
  {
    if (state_<State_RCDT)
      return false; /* ERROR: RCDT not initialised. */

    std::cout << WHEREAMI << "Encroached boundary segments before RCDT:" << std::endl
	      << boundary_;
    std::cout << WHEREAMI << "Encroached interior segments before RCDT:" << std::endl
	      << interior_;
    std::cout << WHEREAMI << "Skinny triangles before RCDT:" << std::endl << skinny_;
    std::cout << WHEREAMI << "Big triangles before RCDT:" << std::endl << big_;

    Dart dh;

    int loop = 0;
    while (!(boundary_.emptyQ() && interior_.emptyQ() &&
	     skinny_.emptyQ() && big_.emptyQ())) {
      /* Temporary failsafe exit: */
      /*
      loop++;
      if (loop>50000) return false;
      */

      std::cout << WHEREAMI << "RCDT(" << loop << "): (Bo,In,Sk,Bi) = ("
		<< boundary_.countQ() << ","
		<< interior_.countQ() << ","
		<< skinny_.countQ() << ","
		<< big_.countQ() << ")" << std::endl;

      dh = boundary_.quality();
      if (!dh.isnull()) {
	std::cout << WHEREAMI << "Encroached boundary segment: "
		  << dh << " "
		  << big_.quality(dh) << std::endl;
	bisectEdgeDelaunay(dh);
	continue;
      }
      
      dh = interior_.quality();
      if (!dh.isnull()) {
	std::cout << WHEREAMI << "Encroached interior segment: "
		  << dh << " "
		  << big_.quality(dh) << std::endl;
	bisectEdgeDelaunay(dh);
	continue;
      }

      dh = skinny_.quality();
      if (!dh.isnull()) {
	std::cout << WHEREAMI << "Skinny triangle: "
		  << dh << " "
		  << skinny_.quality(dh) << std::endl;
	Point c;
	M_->triangleCircumcenter(dh.t(),c);
	if ((!buildRCDTlookahead(&boundary_,c)) ||
	    (!buildRCDTlookahead(&interior_,c)))
	  continue;
	killTriangle(dh);
	continue;
      }
      
      dh = big_.quality();
      if (!dh.isnull()) {
	std::cout << WHEREAMI << "Big triangle: "
		  << dh << " "
		  << big_.quality(dh) << std::endl;
	Point c;
	M_->triangleCircumcenter(dh.t(),c);
	if ((!buildRCDTlookahead(&boundary_,c)) ||
	    (!buildRCDTlookahead(&interior_,c)))
	  continue;
	killTriangle(dh);
	continue;
      }
      
    }

    M_->redrawX11("RCDT finished");

    return true;
  };

  bool MeshC::RCDT(double skinny_limit, double big_limit)
  {
    if (!prepareRCDT(skinny_limit,big_limit)) return false;
    return buildRCDT();
  };


  bool MeshC::PruneExterior()
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








  Dart MeshC::swapEdgeLOP(const Dart& d, MCQswapable& swapable)
  {
    if (!swapable.swapable(d)) {
      /* Not allowed to swap. */
      std::cout << WHEREAMI << "LOP: Not allowed to swap dart " << std::endl
		<< d << std::endl;
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

    // std::cout << WHEREAMI << "LOP edge list: ("
    // 	      << edge_list[0] << ","
    // 	      << edge_list[1] << ","
    // 	      << edge_list[2] << ","
    // 	      << edge_list[3] << ")" << std::endl;
    
    Dart dnew(swapEdge(d));
    if (dh == dnew) {
      /* ERROR: this should not happen. */
      std::cout << WHEREAMI << "Edge swap appears to have failed!" << std::endl;
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

    //    std::cout << WHEREAMI << "Edge swapped:" << std::endl;

    return dnew;
  }

  Dart MeshC::swapEdge(const Dart& d)
  {
    if (state_ < State_CDT) {
      return M_->swapEdge(d);
    }

    if (boundary_.segm(d) || interior_.segm(d)) {
      /* ERROR: Not allowed to swap. */
      std::cout << WHEREAMI << "ERROR: Not allowed to swap dart " << std::endl
		<< d << std::endl;
      return d;
    }

    /* Collect CDT data */
    bool segm_b[4];
    bool segm_i[4];
    Dart dh(d);
    dh.orbit2rev();
    if ((segm_b[1] = boundary_.found(dh))) boundary_.erase(dh);
    if ((segm_i[1] = interior_.found(dh))) interior_.erase(dh);
    dh.orbit2rev();
    if ((segm_b[2] = boundary_.found(dh))) boundary_.erase(dh);
    if ((segm_i[2] = interior_.found(dh))) interior_.erase(dh);
    dh.orbit0().orbit2rev();
    if ((segm_b[3] = boundary_.found(dh))) boundary_.erase(dh);
    if ((segm_i[3] = interior_.found(dh))) interior_.erase(dh);
    dh.orbit2rev();
    if ((segm_b[0] = boundary_.found(dh))) boundary_.erase(dh);
    if ((segm_i[0] = interior_.found(dh))) interior_.erase(dh);

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

    /* Reassemble CDT data */
    dh = dnew;
    dh.orbit2();
    if (segm_b[1]) boundary_.insert(dh);
    if (segm_i[1]) interior_.insert(dh);
    dh.orbit2();
    if (segm_b[0]) boundary_.insert(dh);
    if (segm_i[0]) interior_.insert(dh);
    dh.orbit2().orbit0rev();
    if (segm_b[3]) boundary_.insert(dh);
    if (segm_i[3]) interior_.insert(dh);
    dh.orbit2();
    if (segm_b[2]) boundary_.insert(dh);
    if (segm_i[2]) interior_.insert(dh);

    if (state_>=State_RCDT) {
      /* Reassemble RCDT data */
      dh = dnew;
      skinny_.insert(dh);
      big_.insert(dh);
      dh.orbit1();
      skinny_.insert(dh);
      big_.insert(dh);
    }

    //    std::cout << WHEREAMI << "Edge swapped, boundary segments:" << std::endl
    //	      << boundary_;

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
    for (int i=0;i<3;i++) {
      if ((segm_b[i] = boundary_.found(dh))) boundary_.erase(dh);
      if ((segm_i[i] = interior_.found(dh))) interior_.erase(dh);
      dh.orbit2();
    }
    if (!dh.onBoundary()) {
      dh.orbit1();
      for (int i=3;i<6;i++) {
	if ((segm_b[i] = boundary_.found(dh))) boundary_.erase(dh);
	if ((segm_i[i] = interior_.found(dh))) interior_.erase(dh);
	dh.orbit2();
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

    /* Reassemble CDT data */
    dh = dnew;
    if (segm_b[0]) boundary_.insert(dh);
    if (segm_i[0]) interior_.insert(dh);
    dh.orbit2();
    if (segm_b[1]) boundary_.insert(dh);
    if (segm_i[1]) interior_.insert(dh);
    dh.orbit2().orbit0rev();
    if (segm_b[2]) boundary_.insert(dh);
    if (segm_i[2]) interior_.insert(dh);
    dh.orbit2();
    if (segm_b[0]) boundary_.insert(dh);
    if (segm_i[0]) interior_.insert(dh);
    if (!dh.onBoundary()) {
      dh.orbit1();
      if (segm_b[3]) boundary_.insert(dh);
      if (segm_i[3]) interior_.insert(dh);
      dh.orbit2();
      if (segm_b[4]) boundary_.insert(dh);
      if (segm_i[4]) interior_.insert(dh);
      dh.orbit2().orbit0rev();
      if (segm_b[5]) boundary_.insert(dh);
      if (segm_i[5]) interior_.insert(dh);
      dh.orbit2();
      if (segm_b[3]) boundary_.insert(dh);
      if (segm_i[3]) interior_.insert(dh);
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

    //    std::cout << WHEREAMI << "Edge split, boundary segments:" << std::endl
    //	      << boundary_;

    //    xtmpl_press_ret("Edge has been split");

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
    for (int i=0;i<3;i++) {
      if ((segm_b[i] = boundary_.found(dh))) boundary_.erase(dh);
      if ((segm_i[i] = interior_.found(dh))) interior_.erase(dh);
      dh.orbit2();
    }

    if (state_>=State_RCDT) {
      /* Collect RCDT data */
      skinny_.erase(d);
      big_.erase(d);
    }

    Dart dnew(M_->splitTriangle(d,v));

    /* Reassebmle CDT data */
    dh = dnew;
    for (int i=0;i<3;i++) {
      dh.orbit2();
      if (segm_b[i]) boundary_.insert(dh);
      if (segm_i[i]) interior_.insert(dh);
      dh.orbit2rev().orbit0();
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

    //    std::cout << WHEREAMI << "Triangle split, boundary segments:" << std::endl
    //	      << boundary_;

    return dnew;
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
    output << "N,n = " << Q.count() << "," << Q.countQ() << std::endl;
    for (MCQ::map_type::const_iterator qi = Q.darts_.begin();
	 qi != Q.darts_.end(); qi++) {
      output << ' ' << qi->first
	     << ' ' << std::scientific << qi->second
	     << ' ' << Q.foundQ(qi->first)
	     << std::endl;
    }
    return output;
  }

































  
  


} /* namespace fmesh */
