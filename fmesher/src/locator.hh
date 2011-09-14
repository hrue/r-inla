#ifndef _FMESH_LOCATOR_
#define _FMESH_LOCATOR_ 1

#include <cstddef>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <vector>
#include <set>
#include <map>
#include <list>
#include <string>
#include <cmath>

#include "vector.hh"
#include "mesh.hh"
#include "trees.hh"

#ifndef WHEREAMI
#define WHEREAMI __FILE__ << "(" << __LINE__ << ")\t"
#endif

#ifndef LOG_
#define LOG_(msg) std::cout << WHEREAMI << msg;
#endif

#ifndef NOT_IMPLEMENTED
#define NOT_IMPLEMENTED (std::cout					\
			 << __FILE__ << "(" << __LINE__ << ")\t"	\
			 << "NOT IMPLEMENTED: "				\
			 << __PRETTY_FUNCTION__ << std::endl);
#endif

namespace fmesh {

  typedef std::pair<double,double> bbox_side_type;
  typedef std::vector<bbox_side_type> bbox_side_list_type;
  typedef std::vector<bbox_side_list_type> bbox_type;

  template <class T>
  class BBoxLocator {

    /*! Container for the bbox search tree */
    class search_tree {
      typedef IntervalTree<T> I_type;
      typedef SegmentTree<T, SegmentSet<double> > S_type;
      typedef SegmentTree<T, I_type > SI_type;
      typedef SegmentTree<T, S_type > SS_type;
      typedef SegmentTree<T, SI_type > SSI_type;
      typedef SegmentTree<T, SS_type > SSS_type;

      int ndim_;
      bool use_interval_tree_;
      I_type* I_;
      S_type* S_;
      SI_type* SI_;
      SS_type* SS_;
      SSI_type* SSI_;
      SSS_type* SSS_;

    private:

      template <class TT>
      void init(TT** tree, const bbox_type::iterator& bbox) {
	(*tree) = new TT(bbox);
      }

      template <class TT>
      void add_segment(TT* tree, int start, int end) {
	(*tree).add_segment(start,end);
      }

      template <class TT>
      void build_tree(TT* tree) {
	(*tree).build_tree();
      }

      template <class TT>
      std::ostream& print(TT* tree, std::ostream& output) const {
	output << *tree;
	return output;
      };

    public:
      search_tree(int ndim, bool use_interval_tree=true) : ndim_(ndim), use_interval_tree_(use_interval_tree), I_(NULL), SI_(NULL), SSI_(NULL), S_(NULL), SS_(NULL), SSS_(NULL) {};

      ~search_tree();

      void init(const bbox_type::iterator& bbox);
      
      std::ostream& print(std::ostream& output);
    };

    /*! Container for the bbox search result iterator */
    class search_iterator {
      /*
      I_type::search_iterator* I_;
      S_type::search_iterator* S_;
      SI_type::search_iterator* SI_;
      SS_type::search_iterator* SS_;
      SSI_type::search_iterator* SSI_;
      SSS_type::search_iterator* SSS_;
      */
    public:
      search_iterator() {};
    };
    
  private:
    int ndim_; /*! The number of search dimensions, at most 3 */
    search_tree search_tree_; /*! The search tree container */

  public:
    BBoxLocator(int ndim,
		bool use_interval_tree = true) : 
      ndim_(ndim),
      search_tree_(ndim,use_interval_tree)
    {
    };

    void init(const bbox_type::iterator& bbox) {
      search_tree_.init(bbox);
    };

    ~BBoxLocator() {
    };
 
    std::ostream& print(std::ostream& output);
 
    template <class TT> friend
    std::ostream& operator<<(std::ostream& output, BBoxLocator<TT>& segm);

  };




  class TriangleLocator {

  private:
    const Mesh* mesh_; /*! The mesh to be searched */
    std::vector<int> dim_; /*! The order of the search dimensions, at most 3 */
    bbox_type bbox_; /*! Bounding boxes */
    BBoxLocator<double> bbox_locator_; /*! The bbox searcher object */

  public:
    TriangleLocator(const Mesh* mesh,
		    const std::vector<int>& dimensions,
		    bool use_interval_tree = true) : 
      mesh_(mesh),
      dim_(dimensions),
      bbox_(),
      bbox_locator_(dimensions.size(),use_interval_tree)
    {
      bbox_.resize(dim_.size());
      if (mesh_) {
	for (int i=0; i<dim_.size(); i++) {
	  bbox_[i].resize(mesh_->nT());
	}
	
	/* Build boxes: */
	int d;
	Point mini;
	Point maxi;
	std::pair<double, double> range;
	for (int t=0; t<mesh_->nT(); t++) {
	    mesh_->triangleBoundingBox(t,mini,maxi);
	    for (int di=0; di<dim_.size(); di++) {
	      d = dim_[di];
	      range.first = mini[d];
	      range.second = maxi[d];
	      bbox_[di][t] = range;
	    }
	}
      }
      
      bbox_locator_.init(bbox_.begin());
    };

    ~TriangleLocator() {
    };


    std::ostream& print(std::ostream& output);

    friend 
    std::ostream& operator<<(std::ostream& output, TriangleLocator& locator);

  };



} /* namespace fmesh */

#include "locator.tcc"

#endif
