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
  public:
    /*! Container for the bbox search tree */
    class Search_tree_type {
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
      }

    public:
      Search_tree_type(int ndim, bool use_interval_tree=true) : ndim_(ndim), use_interval_tree_(use_interval_tree), I_(NULL), S_(NULL), SI_(NULL), SS_(NULL), SSI_(NULL), SSS_(NULL) {};
      
      ~Search_tree_type();
      
      void init(const bbox_type::iterator& bbox);
      
      std::ostream& print(std::ostream& output);
      
      
      /*! Container for the bbox search result iterator */
      class Iterator : public std::iterator< 
	std::forward_iterator_tag, T, int, const T*, const T& >
      {
	bool is_null_;
	const Search_tree_type* search_tree_;
	typename I_type::search_iterator I_;
	typename S_type::search_iterator S_;
	typename SI_type::search_iterator SI_;
	typename SS_type::search_iterator SS_;
	typename SSI_type::search_iterator SSI_;
	typename SSS_type::search_iterator SSS_;
	std::vector<T> loc_;
	
	template <class TreeType>
	void init(TreeType* t, typename TreeType::search_iterator* i) {
	  typename std::vector<T>::const_iterator loc_i = loc_.begin();
	  (*i) = t->search(loc_i);
	  is_null_ = (*i).is_null();
	}
	
	template <class TreeType_iter>
	void next(TreeType_iter* i) {
	  ++(*i);
	  is_null_ = (*i).is_null();
	}
	
      public:
	Iterator(const Search_tree_type* search_tree, const std::vector<T>& loc);
 	Iterator();
	~Iterator();
	
	bool is_null() const {
	  return is_null_;
	}
	
	/*	
		bool operator==(const Iterator& b) const {
		NOT_IMPLEMENTED;
		return (is_null() == b.is_null());
		};
	*/
	/*
	  bool operator!=(const Iterator& b) const {
	  return (!(*this == b));
	  };
	*/
	
	int operator*() const;
	Iterator& operator++();
	
      };
      
      
      
      
    }; // Search_tree_type
    typedef typename Search_tree_type::Iterator search_iterator;

    
  private:
    int ndim_; /*! The number of search dimensions, at most 3 */
    Search_tree_type search_tree_; /*! The search tree container */

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
 
    search_iterator search_begin(const std::vector<double>& loc) const {
      return search_iterator(&search_tree_, loc);
    };
    search_iterator search_end() const {
      return search_iterator();
    };

    std::ostream& print(std::ostream& output);
 
    template <class TT> friend
    std::ostream& operator<<(std::ostream& output, BBoxLocator<TT>& segm);

  };




  class TriangleLocator {

  public:    
    typedef BBoxLocator<double> bbox_locator_type;
    
  private:
    const Mesh* mesh_; /*! The mesh to be searched */
    std::vector<int> dim_; /*! The order of the search dimensions, at most 3 */
    bbox_type bbox_; /*! Bounding boxes */
    bbox_locator_type bbox_locator_; /*! The bbox searcher object */
    
  public:
    TriangleLocator(const Mesh* mesh,
		    const std::vector<int>& dimensions,
		    bool use_interval_tree = true);

    ~TriangleLocator();

    int locate(const Point& s) const;

    std::ostream& print(std::ostream& output);

    friend 
    std::ostream& operator<<(std::ostream& output, TriangleLocator& locator);

  };



} /* namespace fmesh */

#include "locator.tcc"

#endif
