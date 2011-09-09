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

#ifndef NOT_IMPLEMENTED
#define NOT_IMPLEMENTED (std::cout					\
			 << __FILE__ << "(" << __LINE__ << ")\t"	\
			 << "NOT IMPLEMENTED: "				\
			 << __PRETTY_FUNCTION__ << std::endl);
#endif

namespace fmesh {

  template < class TreeValueType, class ValueType >
  class SBBTreeIterator;

  /*! Static Balanced Binary Tree */
  template < class ValueType >
  class SBBTree {
  public:
    typedef ValueType value_type;
    typedef SBBTreeIterator< ValueType, ValueType > iterator;
    typedef SBBTreeIterator< ValueType, ValueType const > const_iterator;
    friend class SBBTreeIterator< ValueType, ValueType >;
    friend class SBBTreeIterator< ValueType, ValueType const >;
    
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
      return iterator(*this).leftmost();
    };
    iterator end() {
      return iterator(*this,-1);
    };

    const_iterator root() const {
      return const_iterator(*this);
    };
    const_iterator begin() const {
      return const_iterator(*this).leftmost();
    };
    const_iterator end() const {
      return const_iterator(*this,-1);
    };

  };



  template < class TreeValueType, class ValueType >
  class SBBTreeIterator : public std::iterator< std::bidirectional_iterator_tag, ValueType, int > {

    typedef SBBTreeIterator< TreeValueType, ValueType > self_type;
    typedef SBBTree< TreeValueType > TreeT;
    friend class SBBTree< TreeValueType >;

  protected:
    TreeT& tree_;
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
	if (ret >= tree_.size())
	  ret = -1;
      } else
	ret = -1;
      return ret;
    }
    int right_idx() const {
      int ret;
      if (current_>=0) {
	ret = 2*current_+2;
	if (ret >= tree_.size())
	  ret = -1;
      } else
	ret = -1;
      return ret;
    }
  public:
    SBBTreeIterator(TreeT& tree, int idx=0) : tree_(tree), current_(idx) {
      if (current_>=tree_.size())
	current_ = -1;
    };

    operator SBBTreeIterator< TreeValueType, TreeValueType const >() const {
      return SBBTreeIterator< TreeValueType, TreeValueType const >(tree_, current_);
    };

    typename self_type::reference operator*() const { return tree_.storage_[current_]; };
    typename self_type::pointer operator->() const { return &(tree_.storage_[current_]); };

    /*! -1=left, 0=root or null, 1=right */
    int node_type() const {
      if (current_<=0)
	return 0;
      else if ((current_ % 2) == 0)
	return 1;
      else
	return -1;
    };
    bool is_root() const {
      return (current_==0);
    };
    bool is_left() const {
      return (current_>0) && ((current_ % 2) == 1);
    };
    bool is_right() const {
      return (current_>0) && ((current_ % 2) == 0);
    };
    self_type& parent() {
      current_ = parent_idx();
      return *this;
    };
    self_type& left() {
      current_ = left_idx();
      return *this;
    };
    self_type& right() {
      current_ = right_idx();
      return *this;
    };
    self_type& leftmost() {
      while (left_idx()>=0)
	left();
      return *this;
    };
    self_type& rightmost() {
      while (right_idx()>=0)
	right();
      return *this;
    };
    self_type& next() {
      if (right_idx()>=0)
	return right().leftmost();
      else {
	while (is_right())
	  parent();
	return parent();
      }
    };
    self_type& prev() {
      if (left_idx()>=0)
	return left().rightmost();
      else {
	while (is_left())
	  parent();
	return parent();
      }
    };
  };


  void test(void) {
    SBBTree<int> tree(10);
    SBBTree<int>::iterator a_tree_iterator = tree.begin();
    SBBTree<int>::const_iterator a_tree_const_iterator = tree.begin();
    SBBTree<int>::iterator a_tree_iterator_(tree.begin());
    SBBTree<int>::const_iterator a_tree_const_iterator_(tree.begin());
  }



} /* namespace fmesh */

//#include "trees.tcc"

#endif
