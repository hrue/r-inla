#ifndef _FMESH_TREES_T_
#define _FMESH_TREES_T_ 1

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

  
  template < class ValueType >
  template < class RefValueType, class TreeRefType >
  SBBTree<ValueType>::Iterator<RefValueType,TreeRefType>
  SBBTree<ValueType>::Iterator<RefValueType,TreeRefType>::leftmost() {
    self_type i(*this);
    while (i.left_idx()>=0)
      i = i.left();
    return i;
  };
  
  
  template < class ValueType >
  template < class RefValueType, class TreeRefType >
  SBBTree<ValueType>::Iterator<RefValueType,TreeRefType>
  SBBTree<ValueType>::Iterator<RefValueType,TreeRefType>::rightmost() {
    self_type i(*this);
    while (i.right_idx()>=0)
      i = i.right();
    return i;
  };
  template < class ValueType >
  template < class RefValueType, class TreeRefType >
  SBBTree<ValueType>::Iterator<RefValueType,TreeRefType>
  SBBTree<ValueType>::Iterator<RefValueType,TreeRefType>::next() {
    self_type i(*this);
    if (i.right_idx()>=0)
      return i.right().leftmost();
    else {
      while (i.is_right())
	i = i.parent();
      return i.parent();
    }
  };
  template < class ValueType >
  template < class RefValueType, class TreeRefType >
  SBBTree<ValueType>::Iterator<RefValueType,TreeRefType>
  SBBTree<ValueType>::Iterator<RefValueType,TreeRefType>::prev() {
    self_type i(*this);
    if (i.left_idx()>=0)
      return i.left().rightmost();
    else {
      while (i.is_left())
	i = i.parent();
      return i.parent();
    }
  };
  
  template < class ValueType >
  template < class RefValueType, class TreeRefType >
  SBBTree<ValueType>::Iterator<RefValueType,TreeRefType>&
  SBBTree<ValueType>::Iterator<RefValueType,TreeRefType>::operator--() {
    if (current_<0) {
      *this = tree_->root().rightmost();
    } else
      *this = prev();
    return *this;
  };
  template < class ValueType >
  template < class RefValueType, class TreeRefType >
  SBBTree<ValueType>::Iterator<RefValueType,TreeRefType>&
  SBBTree<ValueType>::Iterator<RefValueType,TreeRefType>::operator++() {
    *this = next();
    return *this;
  };
  
  
} /* namespace fmesh */

#endif
