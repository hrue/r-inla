#ifndef _FMESH_TREES_
#define _FMESH_TREES_ 1

#include <cstddef>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <vector>
#include <set>
#include <map>
//#include <list>
#include <string>
//#include <cmath>

#ifndef WHEREAMI
#define WHEREAMI __FILE__ << "(" << __LINE__ << ")\t"
#endif

#ifndef LOG_
#define LOG_(msg) std::cout << WHEREAMI << msg;
#endif

#ifndef NOT_IMPLEMENTED
#define NOT_IMPLEMENTED (std::cout					\
			 << WHEREAMI	\
			 << "NOT IMPLEMENTED: "				\
			 << __PRETTY_FUNCTION__ << std::endl);
#endif

namespace fmesh {


  /*! Static Balanced Binary Tree */
  template < class ValueType >
  class SBBTree {
  public:
    typedef ValueType value_type;
    template < class RefValueType, class TreeRefType > class Iterator;
    typedef Iterator< ValueType, SBBTree< ValueType > > iterator;
    typedef Iterator< const ValueType, const SBBTree< ValueType > > const_iterator;
    
  protected:
    const int n_; /*!< The number of nodes */
    std::vector< ValueType > storage_; /*!< Linear storage of the tree nodes */

  public:
    SBBTree(int n) : n_(n), storage_(n) {};
    ~SBBTree() {
    };

    int size() const { return n_; };

    iterator root() {
      return iterator(this);
    };
    iterator begin() {
      iterator i(this);
      return i.leftmost();
    };
    iterator end() {
      return iterator(this,-1);
    };

    const_iterator const_root() const {
      return const_iterator(this);
    };
    const_iterator const_begin() const {
      return const_iterator(this).leftmost();
    };
    const_iterator const_end() const {
      return const_iterator(this,-1);
    };


    
    template < class RefValueType, class TreeRefType >
    class Iterator : public std::iterator<
      std::bidirectional_iterator_tag, ValueType, int,
      RefValueType*, RefValueType& > {
      
      typedef Iterator< RefValueType, TreeRefType > self_type;
      typedef SBBTree< ValueType > TreeT;
      
    protected:
      TreeRefType* tree_;
      int current_;
      
      int parent_idx() const {
	int ret;
	if (current_>0)
	  ret = (current_-1)/2;
	else
	  ret = -1;
	return ret;
      }
      int left_idx() const {
	int ret;
	if (current_>=0) {
	  ret = 2*current_+1;
	  if (ret >= tree_->size())
	    ret = -1;
	} else
	  ret = -1;
	return ret;
      }
      int right_idx() const {
	int ret;
	if (current_>=0) {
	  ret = 2*current_+2;
	  if (ret >= tree_->size())
	    ret = -1;
	} else
	  ret = -1;
	return ret;
      }
    public:
      Iterator(TreeRefType* tree, int idx=0) : tree_(tree), current_(idx) {
	if (current_>=tree_->size())
	  current_ = -1;
      };
      Iterator(const Iterator& copy) :
	tree_(copy.tree_), current_(copy.current_) { };
      Iterator& operator=(const Iterator& copy) {
	tree_ = copy.tree_;
	current_ = copy.current_;
	return *this;
      };
      
      operator SBBTree<ValueType>::const_iterator() const {
	return SBBTree<ValueType>::const_iterator(tree_, current_);
      };

      int current() const { return current_; };
      
      typename self_type::reference operator*() const {
	return tree_->storage_[current_];
      };
      typename self_type::pointer operator->() const {
	return &(tree_->storage_[current_]);
      };
      bool operator==(const self_type& i) const {
	return (current_==i.current_);
      };
      bool operator!=(const self_type& i) const {
	return (current_!=i.current_);
      };
      
      /*! -1=left, 0=root or null, 1=right */
      int node_classification() const {
	if (current_<=0)
	  return 0;
	else if ((current_ % 2) == 0)
	  return 1;
	else
	  return -1;
      };
      bool is_null() const {
	return (current_<0);
      };
      bool is_root() const {
	return (current_==0);
      };
      bool is_leaf() const {
	return ((left_idx()<0) && (right_idx()<0));
      };
      bool is_left() const {
	return (current_>0) && ((current_ % 2) == 1);
      };
      bool is_right() const {
	return (current_>0) && ((current_ % 2) == 0);
      };
      self_type parent() {
	return self_type(tree_,parent_idx());
      };
      self_type left() {
	return self_type(tree_,left_idx());
      };
      self_type right() {
	return self_type(tree_,right_idx());
      };
      self_type leftmost();
      self_type rightmost();
      self_type next();
      self_type prev();
      self_type& operator--();
      self_type& operator++();
      
    }; // SBBTree::Iterator
    

  }; // SBBTree




  template < class ValueType >
  class SegmentSet {
  public:
    typedef ValueType value_type;
    typedef std::pair<ValueType,ValueType> segment_type;
    typedef std::vector<segment_type> segment_vector_type;
    typedef std::vector<segment_vector_type> multi_segment_type;
    typedef std::set<int> set_type;

  protected:
    typename multi_segment_type::iterator multi_segment_iter_;
    set_type data_;
  public:
    SegmentSet(const typename multi_segment_type::iterator& segm_iter) : multi_segment_iter_(segm_iter), data_() {
    };
    ~SegmentSet() {
    };

    void add_segment(int segm_idx) {
      data_.insert(segm_idx);
    };
    void add_segment(int start_idx, int end_idx) {
      for (int i=start_idx; i<end_idx; ++i)
	add_segment(i);
    };
    void build_tree() {
      /* Nothing to be done. */
    };

    template<class VT> friend 
    std::ostream& operator<<(std::ostream& output, SegmentSet<VT>& segm);
  }; // SegmentSet


  template < class ValueType >
  class OrderedSegmentSet {
  public:
    typedef ValueType value_type;
    typedef std::pair<ValueType,ValueType> segment_type;
    typedef std::vector<segment_type> segment_vector_type;
    typedef std::vector<segment_vector_type> multi_segment_type;
    typedef std::pair<ValueType, int> mapping_type;
    typedef std::multimap<ValueType, int > map_type;

  protected:
    typename multi_segment_type::iterator multi_segment_iter_;
    map_type L_data_;
    map_type R_data_;
  public:
    OrderedSegmentSet(const typename multi_segment_type::iterator& segm_iter) : multi_segment_iter_(segm_iter), L_data_(), R_data_() {
    };
    ~OrderedSegmentSet() {
    };
    
    void add_segment(int segm_idx) {
      const segment_type& segm = (*multi_segment_iter_)[segm_idx];
      L_data_.insert(mapping_type(segm.first, segm_idx));
      R_data_.insert(mapping_type(segm.second, segm_idx));
    };
    void add_segment(int start_idx, int end_idx) {
      for (int i=start_idx; i<end_idx; ++i)
	add_segment(i);
    };
    void build_tree() {
      /* Nothing to be done. */
    };
    
    template<class VT> friend 
    std::ostream& operator<<(std::ostream& output, OrderedSegmentSet<VT>& segm);
  }; // OrderedSegmentSet




  template < class ValueType >
  class IntervalTree {
  public:
    typedef IntervalTree<ValueType> self_type;
    typedef ValueType value_type;
    typedef std::pair<ValueType,ValueType> segment_type;
    typedef std::vector<segment_type> segment_vector_type;
    typedef std::vector<segment_vector_type> multi_segment_type;
    typedef std::pair<ValueType, int> data_pair_type;
    typedef OrderedSegmentSet<ValueType> data_type;

    class node_type {
      friend class IntervalTree<ValueType>;
    protected:
      value_type mid_;
      data_type* data_;
    public:
      node_type() : data_(NULL) {};
      ~node_type() {
	if (data_) {
	  delete data_;
	  data_ = NULL;
	}
      };

      void activate_data(const typename multi_segment_type::iterator& segm_iter) {
	if (!data_) {
	  data_ = new data_type(segm_iter);
	}
      };
    };
  private:
    typedef SBBTree<node_type> tree_type;
    typedef std::set<value_type> breakpoints_type;
    typedef std::vector<int> segment_list_type;

    typename multi_segment_type::iterator multi_segment_iter_;
    segment_list_type segments_;
    breakpoints_type breakpoints_;
    tree_type* tree_;

    void distribute_breakpoints(typename tree_type::iterator i,
				typename breakpoints_type::const_iterator& breakpoint) {
      if (i.current()<0)
	return;
      if (i.is_leaf()) {
	(*i).mid_ = *breakpoint;
	{
	  typename breakpoints_type::const_iterator tmp = breakpoint;
	  ++tmp;
	  if (tmp !=breakpoints_.end())
	    breakpoint = tmp;
	}
      } else {
	distribute_breakpoints(i.left(),breakpoint);
	(*i).mid_ = *breakpoint;
	{
	  typename breakpoints_type::const_iterator tmp = breakpoint;
	  ++tmp;
	  if (tmp !=breakpoints_.end())
	    breakpoint = tmp;
	}
	distribute_breakpoints(i.right(),breakpoint);
      }
    };

    void distribute_segment(typename tree_type::iterator i,
			    int segm_idx) {
      if (i.current()<0)
	return;
      const segment_type& segm = (*multi_segment_iter_)[segm_idx];
      if ((segm.first <= (*i).mid_) && (segm.second >= (*i).mid_)) {
	/* Segment covers the midpoint */
	(*i).activate_data(multi_segment_iter_);
	(*i).data_->add_segment(segm_idx);
      } else if (segm.second < (*i).mid_) {
	/* Segment completely to the left of the midpoint */
	distribute_segment(i.left(), segm_idx);
      } else if (segm.first > (*i).mid_) {
	/* Segment completely to the right of the midpoint */
	distribute_segment(i.right(), segm_idx);
      }
    };
    void distribute_segments() {
      for (typename segment_list_type::const_iterator si = segments_.begin();
	   si != segments_.end(); ++si) {
	distribute_segment(tree_->root(), (*si));
      }
    }

  public:
    IntervalTree(const typename multi_segment_type::iterator& segm_iter) : multi_segment_iter_(segm_iter), tree_(NULL) {
    };
    ~IntervalTree() {
      if (tree_) {
	delete tree_;
	tree_ = NULL;
      }
    };

    void add_segment(int segm_idx) {
      const segment_type& segm = (*multi_segment_iter_)[segm_idx];
      segments_.insert(segments_.end(), segm_idx);
      breakpoints_.insert(segm.first);
      breakpoints_.insert(segm.second);
    };
    void add_segment(int start_idx, int end_idx) {
      for (int i=start_idx; i<end_idx; ++i)
	add_segment(i);
    };
    void build_tree() {
      if (tree_) { delete tree_; tree_ = NULL; }
      if (breakpoints_.size()==0) {
	return;
      }
      tree_ = new tree_type(breakpoints_.size());
      typename breakpoints_type::const_iterator bi = breakpoints_.begin();
      distribute_breakpoints(tree_->root(),bi);
      distribute_segments();
    };



    std::ostream& print_subtree(std::ostream& output,
				typename tree_type::iterator i,
				std::string prefix) {
      if (i.current()<0)
	return output;
      //      output << "IT " << prefix << i.current() << " = (" << (*i).mid_ << ")" << std::endl;
      if ((*i).data_) {
	output << *(*i).data_;
      }
      if (!i.is_leaf()) {
	print_subtree(output, i.left(), prefix+":");
	print_subtree(output, i.right(), prefix+":");
      }
      return output;
    };
    
    template<class VT> friend 
    std::ostream& operator<<(std::ostream& output, IntervalTree<VT>& segm);
  }; // IntervalTree





  template < class ValueType, class SubTreeType >
  class SegmentTree {
  public:
    typedef ValueType value_type;
    typedef SegmentTree<ValueType, SubTreeType> self_type;
    typedef std::pair<ValueType,ValueType> segment_type;
    typedef std::vector<segment_type> segment_vector_type;
    typedef std::vector<segment_vector_type> multi_segment_type;

  public:
    class node_type {
    public:
      value_type left_;
      value_type right_;
      SubTreeType* data_;
    public:
      node_type() : data_(NULL) {};
      ~node_type() {
	if (data_) {
	  delete data_;
	  data_ = NULL;
	}
      };

      void activate_data(const typename multi_segment_type::iterator& segm_iter) {
	if (!data_) {
	  typename multi_segment_type::iterator i = segm_iter;
	  ++i;
	  data_ = new SubTreeType(i);
	}
      };
    };

  private:
    typedef SBBTree<node_type> tree_type;
    typedef std::set<value_type> breakpoints_type;
    typedef std::vector<int> segment_list_type;

    typename multi_segment_type::iterator multi_segment_iter_;
    segment_list_type segments_;
    breakpoints_type breakpoints_;
    tree_type* tree_;

    void distribute_breakpoints(typename tree_type::iterator i,
				typename breakpoints_type::const_iterator& breakpoint) {
      if (i.current()<0)
	return;
      (*i).left_ = *breakpoint;
      if (i.is_leaf()) {
	typename breakpoints_type::const_iterator tmp = breakpoint;
	++tmp;
	if (tmp !=breakpoints_.end())
	  breakpoint = tmp;
      } else {
	distribute_breakpoints(i.left(),breakpoint);
	distribute_breakpoints(i.right(),breakpoint);
      }
      (*i).right_ = *breakpoint;
    };

    bool distribute_segment(typename tree_type::iterator i,
			    int segm_idx) {
      if (i.current()<0)
	return false;
      const segment_type& segm = (*multi_segment_iter_)[segm_idx];
      //      LOG_("I=(" << (*i).left_ << "," << (*i).right_ << ")" <<
      //	   " S=(" << segm.first << "," << segm.second << ")" << std::endl);
      if ((segm.first <= (*i).left_) && (segm.second >= (*i).right_)) {
	/* Segment completely covers the interval */
	(*i).activate_data(multi_segment_iter_);
	(*i).data_->add_segment(segm_idx);
	return true;
      } else if ((segm.first <= (*i).right_) && (segm.second >= (*i).left_)) {
	/* Segment at least partially covers the interval */
	bool left_ok = distribute_segment(i.left(), segm_idx);
	bool right_ok = distribute_segment(i.right(), segm_idx);
	if (!(left_ok || right_ok)) {
	  (*i).activate_data(multi_segment_iter_);
	  (*i).data_->add_segment(segm_idx);
	} 
	return true;
      }
    };
    void distribute_segments() {
      for (typename segment_list_type::const_iterator si = segments_.begin();
	   si != segments_.end(); ++si) {
	distribute_segment(tree_->root(), (*si));
      }
    }

    void build_subtrees(typename tree_type::iterator i) {
      if (i != tree_->end()) {
	if ((*i).data_)
	  (*i).data_->build_tree();
	build_subtrees(i.left());
	build_subtrees(i.right());
      }
    };

  public:
    SegmentTree(const typename multi_segment_type::iterator& segm_iter) : multi_segment_iter_(segm_iter), tree_(NULL) {
    };
    ~SegmentTree() {
      if (tree_) {
	delete tree_;
	tree_ = NULL;
      }
    };

    void add_segment(int segm_idx) {
      const segment_type& segm = (*multi_segment_iter_)[segm_idx];
      segments_.insert(segments_.end(), segm_idx);
      breakpoints_.insert(segm.first);
      breakpoints_.insert(segm.second);
    };
    void add_segment(int start_idx, int end_idx) {
      for (int i=start_idx; i<end_idx; ++i)
	add_segment(i);
    };

    void build_tree(void) {
      if (tree_) { delete tree_; tree_ = NULL; }
      if (breakpoints_.size()==0)
	return;
      if (breakpoints_.size()==1) {
	tree_ = new tree_type(1);
      } else
	tree_ = new tree_type(breakpoints_.size()*2-3);
      typename breakpoints_type::const_iterator bi = breakpoints_.begin();
      distribute_breakpoints(tree_->root(),bi);
      distribute_segments();
      build_subtrees(tree_->root());
    };



    std::ostream& print_subtree(std::ostream& output,
				typename tree_type::iterator i,
				std::string prefix) {
      if (i.current()<0)
	return output;
      //      output << "ST " << prefix << i.current() << " = (" << (*i).left_ << "," << (*i).right_ << ")" << std::endl;
      if ((*i).data_) {
	output << *(*i).data_;
      }
      if (!i.is_leaf()) {
	print_subtree(output, i.left(), prefix+":");
	print_subtree(output, i.right(), prefix+":");
      }
      return output;
    };

    template<class VT, class STT> friend 
    std::ostream& operator<<(std::ostream& output, SegmentTree<VT,STT>& segm);
  }; // SegmentTree






  template<class VT>
  std::ostream& operator<<(std::ostream& output,
			   SegmentSet<VT>& segm)
  {
    output << "SegmentSet";
    output << "(" << segm.data_.size() << ")";
    if (segm.data_.size()>0) {
      output << "  ( ";
      for (typename SegmentSet<VT>::set_type::iterator i = segm.data_.begin(); i != segm.data_.end(); ++i) {
	output << (*i) << " ";
      }
      output << ")" << std::endl;
    } else
      output << std::endl;
    return output;
  };

  template<class VT>
  std::ostream& operator<<(std::ostream& output,
			   OrderedSegmentSet<VT>& segm)
  {
    output << "OrderedSegmentSet";
    output << "(" << segm.L_data_.size() << ")";
    if (segm.L_data_.size()>0) {
      output << " L=( ";
      for (typename OrderedSegmentSet<VT>::map_type::iterator i = segm.L_data_.begin(); i != segm.L_data_.end(); ++i) {
	output << (*i).second << " ";
      }
      output << ") R=( ";
      for (typename OrderedSegmentSet<VT>::map_type::iterator i = segm.R_data_.begin(); i != segm.R_data_.end(); ++i) {
	output << (*i).second << " ";
      }
      output << ")" << std::endl;
    } else
      output << std::endl;
    return output;
  };

  template<class VT>
  std::ostream& operator<<(std::ostream& output,
			   IntervalTree<VT>& segm)
  {
    output << "IntervalTree";
    output << "(" << segm.breakpoints_.size() << ")" << std::endl;
    if (segm.tree_) {
      segm.print_subtree(output, segm.tree_->root(), "");
    }
    return output;
  };

  template<class VT, class STT>
  std::ostream& operator<<(std::ostream& output,
			   SegmentTree<VT,STT>& segm)
  {
    output << "SegmentTree";
    output << "(" << segm.breakpoints_.size() << ")" << std::endl;
    if (segm.tree_) {
      segm.print_subtree(output, segm.tree_->root(), "");
    }
    return output;
  };









} /* namespace fmesh */

#include "trees.tcc"

#endif
