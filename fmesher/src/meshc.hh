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

#include "mesh.hh"

#ifndef NOT_IMPLEMENTED
#define NOT_IMPLEMENTED (std::cout				\
			 << __FILE__ << "(" << __LINE__ << ") "	\
			 << "NOT IMPLEMENTED: "			\
			 << __PRETTY_FUNCTION__ << std::endl);
#endif

namespace fmesh {

  class MeshC;

  typedef int constrMetaT;
  class constrT : public std::pair<IntPair,constrMetaT> {
  public:
    constrT(int v0, int v1) {
      first = IntPair(v0,v1);
      second = constrMetaT();
    }
    constrT(int v0, int v1, constrMetaT meta) {
      first = IntPair(v0,v1);
      second = meta;
    }
  };
  typedef std::list<constrT> constrListT;

  
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
    double quality(const Dart& d) const;
    Dart quality() const;
    const_iterator find(const Dart& d) const { return darts_.find(d); };
    const_iterator begin() const { return darts_.begin(); };
    const_iterator end() const { return darts_.end(); };
    const_iteratorQ beginQ() const { return darts_quality_.begin(); };
    const_iteratorQ endQ() const { return darts_quality_.end(); };

    friend std::ostream& operator<<(std::ostream& output, const MCQ& Q);
  };

  class MCQtri : public MCQ {
  protected:
    double quality_limit_;
    double* quality_limits_;
    size_t quality_limits_cap_;
    /*!< Larger values are included in the quality set */
    virtual double calcQtri(const Dart& d) const = 0;
  public:
    MCQtri(MeshC* MC, bool only_quality,
	   double quality_limit, const double* quality_limits = NULL,
	   size_t nQL = 0);
    ~MCQtri() { if (quality_limits_) delete[] quality_limits_; }

    void setQ(double quality_limit, const double* quality_limits = NULL,
	      size_t nQL = 0);
    void setQv(int v, double quality_limit);
    bool usingQv() const { return (quality_limits_ != NULL); };
    double getQv(int v) const {
      if (quality_limits_ == NULL)
	return quality_limit_;
      else {
	return quality_limits_[v];
      }
    }
    double getQ() const { return quality_limit_; };
    double getQ(int t) const;
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
    double quality(const Dart& d) const {
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
  public:
    typedef constrMetaT meta_type;
    typedef std::map<Dart,meta_type> meta_map_type;
    typedef meta_map_type::value_type meta_map_key_type;
    typedef meta_map_type::const_iterator const_iteratorMeta;
  private:
    meta_map_type meta_; /*!< Darts, mapped to metadata */
    double encroached_limit_;
    /*!< Larger values are included in the quality set */
  public:
    MCQsegm(MeshC* MC) : MCQ(MC,false),
			 meta_(),
			 encroached_limit_(10*MESH_EPSILON) {};
    double calcQ(const Dart& d) const;
    bool segm(const Dart& d) const; /*!< true if d or d.orbit1() is found */
    void update(const Dart& d);
    /*!< Update quality, keep metadata, don't insert new */

    meta_type meta(const Dart& d) const;

    void clear();
    void insert(const Dart& d, const meta_type& meta);
    /*!< Insert dart if not existing, with metadata. */
    meta_type erase(const Dart& d);
    /*!< Remove dart if existing, return metadata. */

    friend std::ostream& operator<<(std::ostream& output, const MCQsegm& segm);
 };

  /*!
    \brief Holds darts, marked swapable or not.
    Treats darts as bi-directional.
   */
  class MCQswapable : public MCQ {
  private:
    /* No private data. */
  public:
    MCQswapable(MeshC* MC) : MCQ(MC,false) {};
    bool found(const Dart& d) const;
    bool foundQ(const Dart& d) const;
    double quality(const Dart& d) const;
    void insert(const Dart& d); /*!< Insert dart if not existing. */
    void erase(const Dart& d); /*!< Remove dart if existing. */
    virtual double calcQ(const Dart& d) const;
    bool swapable(const Dart& d) const; /*!< true if d or d.orbit1() is foundQ */
  };

  /*!
    \brief Holds darts, marked Delaunay swapable or not.
    Treats darts as bi-directional.
   */
  class MCQswapableD : public MCQswapable {
  private:
    /* No private data. */
  public:
    MCQswapableD(MeshC* MC) : MCQswapable(MC) {};
    virtual double calcQ(const Dart& d) const;
  };


  /*!
    \brief Class for constructing Delaunay triangulations
  */
  class MeshC {
    friend class MCQtri;
    friend class MCQswapable;
    friend class MCQswapableD;
  public:
    enum State {State_noT=0, /*!< No triangulation present */
		State_CET, /*!< Convex enclosure triangulation */
		State_DT, /*!< Delaunay triangulation */
		State_CDT, /*!< Constrained DT,
			     segment data structures active. */
		State_RCDT /*!< Refined CDT, triangle quality
                                    data structures active. */
    }; /*!< The current triangulation and data structure state. */
    enum Option {Option_null=0, /*!< No options */
		 Option_offcenter_steiner=1 /*!< Use offcenter steiner points */
    }; /*!< Currently active options. */
  private:
    Mesh *M_;
    /* CDT Constraint and segment data structures: */
    constrListT constr_boundary_; /*!< Boundary edge
				    constraints not yet
				    added as segments. */
    constrListT constr_interior_; /*!< Interior edge
				    constraints not yet
				    added as segments. */
    MCQsegm boundary_; /*!< Boundary segments */
    MCQsegm interior_; /*!< Interior segments. */
    /* RCDT triangle quality data structures: */
    MCQskinny skinny_; /*!< Skinny triangles. */
    MCQbig big_; /*!< Big triangles. */
    double* big_limits_; /*!< Big triangle limits. */
    int max_n0_; /*!< Target number of vertices, overriding skinny triangles. */
    int max_n1_; /*!< Target number of vertices, overriding big triangles. */
    /* State variables: */
    State state_; /*!< The current MeshC::State */
    bool is_pruned_; /*!< True if the mesh is a pruned mesh. */
    unsigned int options_;

    bool recSwapDelaunay(const Dart& d0);
    Dart splitTriangleDelaunay(const Dart& td, int v);
    Dart splitEdgeDelaunay(const Dart& ed, int v);
    Dart bisectEdgeDelaunay(const Dart& d);
    void calcSteinerPoint(const Dart& d, Point& c);
    Dart insertNode(int v, const Dart& ed);

    bool isSegment(const Dart& d) const;
    bool buildRCDTlookahead(MCQsegm* segm, const Point& c);

    /*!
      \brief Make a DT from a CET, calling LOP.
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
    Dart CDTInsertSegment(const DartPair& dp, const DartList& trace,
			  triangleSetT& triangles,
			  const bool is_boundary,
			  const constrMetaT& meta);
    /*!
      \brief Detect if a new constraint edge needs to be split,
      and add a new vertex to the CDT if needed.

      Returns the index of the splitting vertex.
      If the return value is negative, splitting was not done.
      If the return value is non-negative, the "trace" list invalidated.
    */
    int CDTSplitSegment(const DartPair& dp, const DartList& trace);
    /*!
      \brief Insert a constraint edge into a CDT, splitting if needed.
    */
    Dart CDTInsertSegment(const int v0, const int v1,
			  triangleSetT& triangles,
			  const bool is_boundary,
			  const constrMetaT& meta);
    /*!
      \brief Build a CDT from constraint edge lists. Called by prepareCDT,
      CDTBoundary and CDTInterior.
    */
    bool buildCDT();
    /*!
      \brief Initialise the RCDT data structures.
    */
    bool prepareRCDT(double skinny_limit,
		     double big_limit,
		     const double* big_limits = NULL,
		     size_t nQL = 0,
		     int max_n0 = -1,
		     int max_n1 = -1);
    /*!
      \brief Build a RCDT.
    */
    bool buildRCDT();

  public:
    MeshC() : M_(NULL), boundary_(this), interior_(this),
	      skinny_(this), big_(this), big_limits_(NULL),
	      max_n0_(-1), max_n1_(-1),
	      state_(State_noT), is_pruned_(false),
	      options_(Option_null) {};
    MeshC(Mesh* M)
      : M_(M), boundary_(this), interior_(this),
	skinny_(this), big_(this), big_limits_(NULL),
	max_n0_(-1), max_n1_(-1),
	state_(State_noT), is_pruned_(false),
	options_(Option_null) {
      if (M_->nT()>0)
	state_ = State_CET;
    };

    unsigned int getOptions() const { return options_; };
    unsigned int setOptions(unsigned int options) {
      return (options_ = options);
    }

    /*!
      \brief Append a vertex

      Return index of the the added point.
    */
    int addVertex(const Point& s);
    /*!
      \brief Append vertices

      Return index of the first of the added points.
    */
    int addVertices(const Matrix3double& S);

    /*! Swap an edge, keeping track of extra swapable darts information. */
    Dart swapEdge(const Dart& d, MCQswapable& swapable);
    /*! Swap an edge. */
    Dart swapEdge(const Dart& d);
    /*! Split an edge in two. */
    Dart splitEdge(const Dart& d, int v);
    /*! Split a triangle in three. */
    Dart splitTriangle(const Dart& d, int v);

    /*! Unlink an edge, keeping track of possible segments. */
    void unlinkEdge(Dart& d);
    /*! Remove a triangle, keeping track of possible segments.
      Returns index of relocated triangle. */
    int removeTriangle(Dart& d);

    double encroachedQuality(const Dart& d) const;
    double skinnyQuality(int t) const;
    double bigQuality(int t) const;

    /*!
      Append segments to segm, with group metadata in segmgrp.

      \param boundary indicates if boundary or interior segments
      should be appended. Call once with true and once with false
      to extract all segments.
      \param segm Where to append the segments.  Set to NULL if only the
      number of segments is to be returned. 
      \param segmgrp Where to append the group metadata for each segment.
      If NULL, the group metadata is discarded.

      \return The number of appended segments (if segm!=NULL) or
              the number of segments in seg (if segm==NULL)
     */
    int segments(bool boundary,
		 Matrix<int>* segm = NULL,
		 Matrix<int>* segmgrp = NULL) const;

    /*!
      \brief Build a convex enclosure triangulation (CET).
    */
    bool CETsphere(int sides, double margin=-0.05);
    /*!
      \brief Build a convex enclosure triangulation (CET).
    */
    bool CETplane(int sides, double margin=-0.05);
    /*!
      \brief Build a convex enclosure triangulation (CET).

      \param sides The number of sides for the enclosure
      \param margin The absolute margin for the enclosure.  If
      negative, the margin is set to -margin multiplied by the
      approximate diameter.
    */
    bool CET(int sides, double margin=-0.05);
    /*!
      \brief Local Optimisation Procedure (LOP)

      Perform LOP to make the input triangulation Delaunay.
      
      \param swapable The triangulation part to be LOPed, as a set of
      swappable darts.
     */
    bool LOP(MCQswapableD& swapable);
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
      \brief Insert a single segment into a constrained
      Delaunay triangulation (CDT).
    */
    Dart CDTSegment(const bool boundary, const constrT& constraint);
    /*!
      \brief Insert a single segment into a constrained
      Delaunay triangulation (CDT).
    */
    Dart CDTSegment(const bool boundary, const int v0, const int v1,
		    const constrMetaT& meta = constrMetaT()) {
      return CDTSegment(boundary,constrT(v0,v1,meta));
    }
    /*!
      \brief Build constrained Delaunay triangulation (CDT)
    */
    bool CDT(const constrListT& boundary, const constrListT& interior);
    /*!
      \brief Remove exterior triangles from a CDT
    */
    bool PruneExterior();
    /*!
      \brief Refine a CDT
    */
    bool RCDT(double angle_limit,
	      double big_limit,
	      const double* big_limits = NULL,
	      size_t nQL = 0,
	      int max_n0 = -1,
	      int max_n1 = -1);

    friend std::ostream& operator<<(std::ostream& output, const MeshC& MC);
  };

  std::ostream& operator<<(std::ostream& output, const MCQsegm& segm);
  std::ostream& operator<<(std::ostream& output, const MeshC& MC);
  std::ostream& operator<<(std::ostream& output, const IntPair& p);
  std::ostream& operator<<(std::ostream& output, const DartPair& dp);
  std::ostream& operator<<(std::ostream& output, const DartList& ds);
  std::ostream& operator<<(std::ostream& output, const std::set<int>& il);
  std::ostream& operator<<(std::ostream& output, const std::list<int>& il);
  std::ostream& operator<<(std::ostream& output, const std::list<IntPair>& il);

} /* namespace fmesh */

#endif
