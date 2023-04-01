#ifndef _FMESH_LOCATOR_T_
#define _FMESH_LOCATOR_T_ 1

#include "fmesher_debuglog.h"

namespace fmesh {

template <class T>
BBoxLocator<T>::Search_tree_type::Iterator::Iterator()
    : is_null_(true), search_tree_(NULL), loc_() {
  /* Nothing to do. */
}

template <class T>
BBoxLocator<T>::Search_tree_type::Iterator::Iterator(
    const Search_tree_type *search_tree, const std::vector<T> &loc)
    : is_null_(false), search_tree_(search_tree), loc_(loc) {
  if (search_tree_->use_interval_tree_)
    switch (search_tree_->ndim_) {
    case 1:
      init(search_tree->I_, &I_);
      break;
    case 2:
      init(search_tree->SI_, &SI_);
      break;
    case 3:
      init(search_tree->SSI_, &SSI_);
      break;
    }
  else
    switch (search_tree_->ndim_) {
    case 1:
      init(search_tree->S_, &S_);
      break;
    case 2:
      init(search_tree->SS_, &SS_);
      break;
    case 3:
      init(search_tree->SSS_, &SSS_);
      break;
    }
}

template <class T> BBoxLocator<T>::Search_tree_type::Iterator::~Iterator() {}

template <class T>
int BBoxLocator<T>::Search_tree_type::Iterator::operator*() const {
  if (search_tree_->use_interval_tree_) {
    switch (search_tree_->ndim_) {
    case 1:
      return *I_;
      break;
    case 2:
      return *SI_;
      break;
    case 3:
      return *SSI_;
      break;
    }
  } else {
    switch (search_tree_->ndim_) {
    case 1:
      return *S_;
      break;
    case 2:
      return *SS_;
      break;
    case 3:
      return *SSS_;
      break;
    }
  }
  FMLOG_("Error: Invalid internal search_tree structure.");
  return (-1);
}

template <class T>
typename BBoxLocator<T>::Search_tree_type::Iterator &
BBoxLocator<T>::Search_tree_type::Iterator::operator++() {
  if (search_tree_->use_interval_tree_)
    switch (search_tree_->ndim_) {
    case 1:
      next(&I_);
      break;
    case 2:
      next(&SI_);
      break;
    case 3:
      next(&SSI_);
      break;
    }
  else
    switch (search_tree_->ndim_) {
    case 1:
      next(&S_);
      break;
    case 2:
      next(&SS_);
      break;
    case 3:
      next(&SSS_);
      break;
    }
  return *this;
}

template <class T> BBoxLocator<T>::Search_tree_type::~Search_tree_type() {
  if (I_)
    delete I_;
  if (SI_)
    delete SI_;
  if (SSI_)
    delete SSI_;
  if (S_)
    delete S_;
  if (SS_)
    delete SS_;
  if (SSS_)
    delete SSS_;
}

template <class T>
void BBoxLocator<T>::Search_tree_type::init(const bbox_type::iterator &bbox) {
  if (use_interval_tree_)
    switch (ndim_) {
    case 1:
      init(&I_, bbox);
      break;
    case 2:
      init(&SI_, bbox);
      break;
    case 3:
      init(&SSI_, bbox);
      break;
    }
  else
    switch (ndim_) {
    case 1:
      init(&S_, bbox);
      break;
    case 2:
      init(&SS_, bbox);
      break;
    case 3:
      init(&SSS_, bbox);
      break;
    }
  if ((*bbox).size() > 0) {
    if (use_interval_tree_)
      switch (ndim_) {
      case 1:
        add_segment(I_, 0, (*bbox).size());
        break;
      case 2:
        add_segment(SI_, 0, (*bbox).size());
        break;
      case 3:
        add_segment(SSI_, 0, (*bbox).size());
        break;
      }
    else
      switch (ndim_) {
      case 1:
        add_segment(S_, 0, (*bbox).size());
        break;
      case 2:
        add_segment(SS_, 0, (*bbox).size());
        break;
      case 3:
        add_segment(SSS_, 0, (*bbox).size());
        break;
      }
  }
  if ((*bbox).size() > 0) {
    if (use_interval_tree_)
      switch (ndim_) {
      case 1:
        build_tree(I_);
        break;
      case 2:
        build_tree(SI_);
        break;
      case 3:
        build_tree(SSI_);
        break;
      }
    else
      switch (ndim_) {
      case 1:
        build_tree(S_);
        break;
      case 2:
        build_tree(SS_);
        break;
      case 3:
        build_tree(SSS_);
        break;
      }
  }
}

template <class T>
std::ostream &BBoxLocator<T>::Search_tree_type::print(std::ostream &output) {
  if (use_interval_tree_)
    switch (ndim_) {
    case 1:
      return print(I_, output);
      break;
    case 2:
      return print(SI_, output);
      break;
    case 3:
      return print(SSI_, output);
      break;
    }
  else
    switch (ndim_) {
    case 1:
      return print(S_, output);
      break;
    case 2:
      return print(SS_, output);
      break;
    case 3:
      return print(SSS_, output);
      break;
    }
  return output;
}

template <class T> std::ostream &BBoxLocator<T>::print(std::ostream &output) {
  return search_tree_.print(output);
}

template <class T>
std::ostream &operator<<(std::ostream &output, BBoxLocator<T> &bbl) {
  return bbl.print(output);
}

} /* namespace fmesh */

#endif
