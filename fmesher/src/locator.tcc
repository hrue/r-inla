#ifndef _FMESH_LOCATOR_T_
#define _FMESH_LOCATOR_T_ 1

#ifndef WHEREAMI
#define WHEREAMI __FILE__ << "(" << __LINE__ << ")\t"
#endif

#ifndef LOG_
#define LOG_(msg) std::cout << WHEREAMI << msg;
#endif

#ifndef NOT_IMPLEMENTED
#define NOT_IMPLEMENTED (std::cout                                      \
                         << WHEREAMI    \
                         << "NOT IMPLEMENTED: "                         \
                         << __PRETTY_FUNCTION__ << std::endl);
#endif

namespace fmesh {

  template<class T>
  BBoxLocator<T>::search_tree::~search_tree() {
    if (I_) delete I_;
    if (SI_) delete SI_;
    if (SSI_) delete SSI_;
    if (S_) delete S_;
    if (SS_) delete SS_;
    if (SSS_) delete SSS_;
  }

  template<class T>
  void BBoxLocator<T>::search_tree::init(const bbox_type::iterator& bbox) {
    if (use_interval_tree_)
      switch (ndim_) {
      case 1: init(&I_, bbox); break;
      case 2: init(&SI_, bbox); break;
      case 3: init(&SSI_, bbox); break;
      }
    else
      switch (ndim_) {
      case 1: init(&S_, bbox); break;
      case 2: init(&SS_, bbox); break;
      case 3: init(&SSS_, bbox); break;
      }
    if ((*bbox).size()>0)
      if (use_interval_tree_)
	switch (ndim_) {
	case 1: add_segment(I_, 0,(*bbox).size()); break;
	case 2: add_segment(SI_, 0,(*bbox).size()); break;
	case 3: add_segment(SSI_, 0,(*bbox).size()); break;
	}
      else
	switch (ndim_) {
	case 1: add_segment(S_, 0,(*bbox).size()); break;
	case 2: add_segment(SS_, 0,(*bbox).size()); break;
	case 3: add_segment(SSS_, 0,(*bbox).size()); break;
	}
    if ((*bbox).size()>0)
      if (use_interval_tree_)
	switch (ndim_) {
	case 1: build_tree(I_); break;
	case 2: build_tree(SI_); break;
	case 3: build_tree(SSI_); break;
	}
      else
	switch (ndim_) {
	case 1: build_tree(S_); break;
	case 2: build_tree(SS_); break;
	case 3: build_tree(SSS_); break;
	}
  }
  

  template <class T>
  std::ostream& BBoxLocator<T>::search_tree::print(std::ostream& output)
  {
    if (use_interval_tree_)
      switch (ndim_) {
      case 1: return print(I_, output); break;
      case 2: return print(SI_, output); break;
      case 3: return print(SSI_, output); break;
      }
    else
      switch (ndim_) {
      case 1: return print(S_, output); break;
      case 2: return print(SS_, output); break;
      case 3: return print(SSS_, output); break;
      }
    return output;
  }

  template <class T>
  std::ostream& BBoxLocator<T>::print(std::ostream& output)
  {
    return search_tree_.print(output);
  }

  template <class T>
  std::ostream& operator<<(std::ostream& output, BBoxLocator<T>& bbl)
  {
    return bbl.print(output);
  }




} /* namespace fmesh */

#endif
