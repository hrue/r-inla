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
  
  



  template < class T >
  typename SegmentSet<T>::search_iterator& SegmentSet<T>::search_iterator::search() {
    if (!this->is_null()) {
      for ( ; i_ != this->C_->data_.end() ; ++i_) {
	if (((*this->C_->multi_segment_iter_)[*i_].first <= this->loc_) &&
	    (this->loc_ <= (*this->C_->multi_segment_iter_)[*i_].second)) {
	  break;
	}
      }
      this->is_null_ = (i_ == this->C_->data_.end());
    }
    return *this;
  };

  template < class T >
  template <class map_type, class Compare >
  typename OrderedSegmentSet<T>::template search_iterator<map_type,Compare>& OrderedSegmentSet<T>::search_iterator<map_type,Compare>::search() {
    if (!this->is_null()) {
      for ( ; i_ != the_end_ ; ++i_) {
	if (!Compare()(this->loc_, (*i_).first)) {
	  break;
	}
      }
      this->is_null_ = (i_ == the_end_);
    }
    return *this;
  };





} /* namespace fmesh */

#endif
