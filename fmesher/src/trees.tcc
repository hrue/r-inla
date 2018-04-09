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
  }
  
  
  template < class ValueType >
  template < class RefValueType, class TreeRefType >
  SBBTree<ValueType>::Iterator<RefValueType,TreeRefType>
  SBBTree<ValueType>::Iterator<RefValueType,TreeRefType>::rightmost() {
    self_type i(*this);
    while (i.right_idx()>=0)
      i = i.right();
    return i;
  }
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
  }
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
  }
  
  template < class ValueType >
  template < class RefValueType, class TreeRefType >
  SBBTree<ValueType>::Iterator<RefValueType,TreeRefType>&
  SBBTree<ValueType>::Iterator<RefValueType,TreeRefType>::operator--() {
    if (current_<0) {
      *this = tree_->root().rightmost();
    } else
      *this = prev();
    return *this;
  }
  template < class ValueType >
  template < class RefValueType, class TreeRefType >
  SBBTree<ValueType>::Iterator<RefValueType,TreeRefType>&
  SBBTree<ValueType>::Iterator<RefValueType,TreeRefType>::operator++() {
    *this = next();
    return *this;
  }
  
  



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
  }

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
  }






  template <class T>
  void IntervalTree<T>::distribute_breakpoints(typename tree_type::iterator i,
					       typename breakpoints_type::const_iterator& breakpoint)
  {
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
  }
  
  template <class T>
  void IntervalTree<T>::distribute_segment(typename tree_type::iterator i,
					   int segm_idx)
  {
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
  }
  
  template <class T>
  void IntervalTree<T>::distribute_segments()
  {
    for (typename segment_list_type::const_iterator si = segments_.begin();
	 si != segments_.end(); ++si) {
      distribute_segment(tree_->root(), (*si));
    }
  }

  template <class T>
  void IntervalTree<T>::add_segment(int segm_idx)
  {
    const segment_type& segm = (*multi_segment_iter_)[segm_idx];
    segments_.insert(segments_.end(), segm_idx);
    breakpoints_.insert(segm.first);
    breakpoints_.insert(segm.second);
  }

  template <class T>
  void IntervalTree<T>::add_segment(int start_idx, int end_idx)
  {
    for (int i=start_idx; i<end_idx; ++i)
      add_segment(i);
  }
  template <class T>
  void IntervalTree<T>::build_tree()
  {
    if (tree_) { delete tree_; tree_ = NULL; }
    if (breakpoints_.size()==0) {
      return;
    }
    tree_ = new tree_type(breakpoints_.size());
    typename breakpoints_type::const_iterator bi = breakpoints_.begin();
    distribute_breakpoints(tree_->root(),bi);
    distribute_segments();
  }


  template < class T >
  typename IntervalTree<T>::search_iterator& IntervalTree<T>::search_iterator::search() {
    this->is_null_ = (i_ == this->C_->tree_->end());
    if (!this->is_null()) {
      if (search_mode_ == 0) {
	if (this->loc_ <= (*i_).mid_) {
	  search_mode_ = -1;
	  if (!(*i_).data_) {
	    i_ = i_.left();
	    search_mode_ = 0;
	    search();
	    this->is_null_ = (i_ == this->C_->tree_->end());
	    return *this;
	  }
	  L_i_ = (*i_).data_->L_search(this->loc_i_);
	  if (L_i_.is_null()) {
	    i_ = i_.left();
	    search_mode_ = 0;
	    search();
	    this->is_null_ = (i_ == this->C_->tree_->end());
	    return *this;
	  }
	} else {
	  search_mode_ = +1;
	  if (!(*i_).data_) {
	    i_ = i_.right();
	    search_mode_ = 0;
	    search();
	    this->is_null_ = (i_ == this->C_->tree_->end());
	    return *this;
	  }
	  R_i_ = (*i_).data_->R_search(this->loc_i_);
	  if (R_i_.is_null()) {
	    i_ = i_.right();
	    search_mode_ = 0;
	    search();
	    this->is_null_ = (i_ == this->C_->tree_->end());
	    return *this;
	  }
	}
      } else if (search_mode_ < 0) {
	LOG_("Should not be reached." << std::endl);
	NOT_IMPLEMENTED;
      } else { /* (search_mode_ > 0) */
	LOG_("Should not be reached." << std::endl);
	NOT_IMPLEMENTED;
      }

      this->is_null_ = (i_ == this->C_->tree_->end());
    }
    return *this;
  }
  
  template < class T >
  typename IntervalTree<T>::search_iterator& IntervalTree<T>::search_iterator::operator++() {
    if (!this->is_null()) {
      if (search_mode_<0) { /* Searching left */
	if (!L_i_.is_null()) {
	  ++L_i_;
	}
	if (L_i_.is_null()) {
	  i_ = i_.left();
	  search_mode_ = 0;
	  this->is_null_ = (i_ == this->C_->tree_->end());
	} else
	  return *this;
      } else if (search_mode_>0) { /* Searching right */
	if (!R_i_.is_null()) {
	  ++R_i_;
	}
	if (R_i_.is_null()) {
	  i_ = i_.right();
	  search_mode_ = 0;
	  this->is_null_ = (i_ == this->C_->tree_->end());
	} else
	  return *this;
      }
      search();
    }
    return *this;
  }


  template < class T, class SubTreeType >
  typename SegmentTree<T,SubTreeType>::search_iterator& SegmentTree<T,SubTreeType>::search_iterator::search() {
    this->is_null_ = (i_ == this->C_->tree_->end());
    if (!this->is_null()) {
      sub_i_ = typename SubTreeType::search_iterator();
      if ((*i_).data_) {
	sub_i_ = (*i_).data_->search(this->loc_next_i_);
      }
      if (sub_i_.is_null()) {
	if (this->loc_ <= (*i_).mid_)
	  i_ = i_.left();
	else
	  i_ = i_.right();
	search();
	this->is_null_ = (i_ == this->C_->tree_->end());
	return *this;
      }
      this->is_null_ = (i_ == this->C_->tree_->end());
    }
    return *this;
  }
  
  template < class T, class SubTreeType >
  typename SegmentTree<T,SubTreeType>::search_iterator& SegmentTree<T,SubTreeType>::search_iterator::operator++() {
     if (!this->is_null()) {
       if (!sub_i_.is_null()) {
	 ++sub_i_;
       }
       if (sub_i_.is_null()) {
	 if (this->loc_ <= (*i_).mid_)
	   i_ = i_.left();
	 else
	   i_ = i_.right();
	 search();
	 this->is_null_ = (i_ == this->C_->tree_->end());
	 return *this;
       }
     }
    return *this;
  }


} /* namespace fmesh */

#endif
