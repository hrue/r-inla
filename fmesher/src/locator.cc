#include <cstddef>
#include <cstring>
#include <set>
#include <map>
#include <sstream>
#include <cmath>

#include "locator.h"

#ifndef LOG_
#define LOG_(msg) FMLOG_(msg)
#endif
#ifdef DEBUG
#define LOG(msg) LOG_(msg)
#else
#define LOG(msg)
#endif


namespace fmesh {

  TriangleLocator::TriangleLocator(const Mesh* mesh,
				   const std::vector<int>& dimensions,
				   bool use_interval_tree) :
    mesh_(mesh),
    dim_(dimensions),
    bbox_(),
    bbox_locator_(dimensions.size(),use_interval_tree)
  {
    bbox_.resize(dim_.size());
    if (mesh_) {
      for (size_t i=0; i<dim_.size(); ++i) {
	bbox_[i].resize(mesh_->nT());
      }

      /* Build boxes: */
      int d;
      Point mini;
      Point maxi;
      std::pair<double, double> range;
      for (size_t t=0; t<mesh_->nT(); ++t) {
	mesh_->triangleBoundingBox(t,mini,maxi);
	for (size_t di=0; di<dim_.size(); ++di) {
	  d = dim_[di];
	  range.first = mini[d];
	  range.second = maxi[d];
	  bbox_[di][t] = range;
	}
      }
    }

    bbox_locator_.init(bbox_.begin());
  }

  TriangleLocator::~TriangleLocator()
  {
    /* Nothing to do. */
  }


  int TriangleLocator::locate(const Point& s) const {
    LOG("Looking for s=" << s << std::endl)
    std::vector<double> loc(dim_.size());
    for (size_t di=0; di<dim_.size(); ++di) {
      loc[di] = s[dim_[di]];
    }
    Dart d;
    for (bbox_locator_type::search_iterator si =
	   bbox_locator_.search_begin(loc);
	 !si.is_null();
	 ++si) {
      LOG("Starting at "<< *si << std::endl)
      d = mesh_->locate_point(Dart(*mesh_,(*si)),s);
      LOG("Resulting dart " << d << std::endl)
      if (!d.isnull()) {
	  Point b;
	  mesh_->barycentric(Dart(*mesh_,d.t()),s,b);
	  LOG("Barycentric coordinates " << b << std::endl)
	  if ((b[0] >= -10.0*MESH_EPSILON) &&
	      (b[1] >= -10.0*MESH_EPSILON) &&
	      (b[2] >= -10.0*MESH_EPSILON))
	    return (d.t());
	  else {
	    LOG("Mesh::locate_point reported incorrect finding." << std::endl);
	  }
      }
    }
    LOG("Point not found, s=" << s << std::endl)
    return -1;
  }

  std::ostream& TriangleLocator::print(std::ostream& output)
  {
    return bbox_locator_.print(output);
  }

  std::ostream& operator<<(std::ostream& output, TriangleLocator& locator)
  {
    return locator.print(output);
  }


} /* namespace fmesh */
