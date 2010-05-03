#ifndef _FMESH_MESHC_
#define _FMESH_MESHC_ 1

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

#include "mesh.h"

#ifndef NOT_IMPLEMENTED
#define NOT_IMPLEMENTED (std::cout				\
			 << __FILE__ << "(" << __LINE__ << ") "	\
			 << "NOT IMPLEMENTED: "			\
			 << __PRETTY_FUNCTION__ << std::endl);
#endif

namespace fmesh {

  class MeshC;
  
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
  public:
    typedef std::map<Dart,double> map_type;
    typedef map_type::value_type map_key_type;
    typedef std::set<MCQdv> set_type;
    typedef map_type::const_iterator const_iterator;
    typedef set_type::const_iterator const_iteratorQ;
  protected:
    MeshC* MC_;
    map_type darts_; /*!< Darts, mapped to quality */
    set_type darts_quality_;
    /*!< Set of "bad quality" (calcQ>0.0) segment darts */
    bool only_quality_;
    /*!< If true, only store darts that are of bad quality */
  public:
    MCQ(MeshC* MC, bool only_quality) : MC_(MC),
					darts_(), darts_quality_(),
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
    const_iterator begin() { return darts_.begin(); };
    const_iterator end() { return darts_.end(); };
    const_iteratorQ beginQ() { return darts_quality_.begin(); };
    const_iteratorQ endQ() { return darts_quality_.end(); };

    friend std::ostream& operator<<(std::ostream& output, const MCQ& Q);
  };

  class MCQtri : public MCQ {
  protected:
    double quality_limit_;
    /*!< Larger values are included in the quality set */
    virtual double calcQtri(const Dart& d) const = 0;
  public:
    MCQtri(MeshC* MC, bool only_quality, double quality_limit);
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
    MCQskinny(MeshC* MC) : MCQtri(MC,true,1.42) {};
    virtual double calcQtri(const Dart& d) const;
  };

  class MCQbig : public MCQtri {
  private:
  public:
    MCQbig(MeshC* MC) : MCQtri(MC,true,1.0) {};
    virtual double calcQtri(const Dart& d) const;
  };

  class MCQsegm : public MCQ {
  private:
    double encroached_limit_;
    /*!< Larger values are included in the quality set */
  public:
    MCQsegm(MeshC* MC) : MCQ(MC,false), encroached_limit_(10*MESH_EPSILON) {};
    double calcQ(const Dart& d) const;
    bool segm(const Dart& d) const; /*!< true if d or d.orbit1() is found */
    void update(const Dart& d); /*!< Update quality, don't insert new */
  };

  class MCQswapable : public MCQ {
  private:
    /* No private data. */
  public:
    MCQswapable(MeshC* MC) : MCQ(MC,false) {};
    bool found(const Dart& d) const;
    bool foundQ(const Dart& d) const;
    const double quality(const Dart& d) const;
    void insert(const Dart& d); /*!< Insert dart if not existing. */
    void erase(const Dart& d); /*!< Remove dart if existing. */
    double calcQ(const Dart& d) const;
    bool swapable(const Dart& d) const; /*!< true if d or d.orbit1() is found */
  };


  /*!
    \brief Class for constructing Delaunay triangulations
  */
  class MeshC {
    friend class MCQswapable;
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
    Dart bisectEdgeDelaunay(const Dart& d);
    bool killTriangle(const Dart& d);
    bool insertNode(int v, const Dart& ed);

    bool isSegment(const Dart& d) const;
    bool buildRCDTlookahead(MCQsegm* segm, const Point& c);

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
      \brief Insert a constraint edge into a CDT.
    */
    Dart CDTinsert(const int v0, const int v1);
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
    MeshC() : M_(NULL), boundary_(this), interior_(this),
	      skinny_(this), big_(this), state_(State_noT),
	      is_pruned_(false) {};
    MeshC(Mesh* M, bool with_conv_hull)
      : M_(M), boundary_(this), interior_(this),
	skinny_(this), big_(this), state_(State_noT),
	is_pruned_(false) {
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

    Dart swapEdgeLOP(const Dart& d, MCQswapable& swapable);
    Dart swapEdge(const Dart& d);
    Dart splitEdge(const Dart& d, int v);
    Dart splitTriangle(const Dart& d, int v);

    double encroachedQuality(const Dart& d) const;
    double skinnyQuality(int t) const;
    double bigQuality(int t) const;

    /*!
      \brief Local Optimisation Procedure (LOP)

      Perform LOP to make the input triangulation Delaunay.
      
      \param swapable The triangulation part to be LOPed, as a set of
      swappable darts.
     */
    bool LOP(MCQswapable& swapable);
    /*!
      \brief Local Optimisation Procedure (LOP)

      Perform LOP to make the input triangulation Delaunay.
      
      \param t_set The triangulation part to be LOPed, as a set
      triangle indices.
     */
    bool LOP(const triangleSetT& t_set);
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

  std::ostream& operator<<(std::ostream& output, const DartOrderedSet& ds);

} /* namespace fmesh */

#endif
