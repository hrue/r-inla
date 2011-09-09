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

#ifndef NOT_IMPLEMENTED
#define NOT_IMPLEMENTED (std::cout					\
			 << __FILE__ << "(" << __LINE__ << ")\t"	\
			 << "NOT IMPLEMENTED: "				\
			 << __PRETTY_FUNCTION__ << std::endl);
#endif

namespace fmesh {

  template <class T>
  class TriangleLocator {
    //    friend
    //    std::ostream& operator<< <> (std::ostream& output,
    //				 const Matrix<T>& M);
    
  protected:
    const Mesh* mesh_; /*! The mesh to be searched */
    int dimensions_; /*! How many spatial dimensions to use  */
    int use_interval_tree_; /*! Use interval tree for dim=0 ? */
    Matrix<double> bbox_mini_; /*! Bounding box corner */
    Matrix<double> bbox_maxi_; /*! Bounding box corner */

  public:
    TriangleLocator() : mesh_(NULL), dimensions_(0),
			use_interval_tree_(false),
			bbox_mini_(0), bbox_maxi_(0) {};
    TriangleLocator(const Mesh* mesh,
		    int dimensions,
		    bool use_interval_tree) : mesh_(mesh),
					      dimensions_(dimensions),
					      use_interval_tree_(use_interval_tree),
					      bbox_mini_(0), bbox_maxi_(0) {
      bbox_mini_.cols(3).rows(mesh_->nT());
      bbox_maxi_.cols(3).rows(mesh_->nT());
    };

    ~TriangleLocator() {
    };

  };



} /* namespace fmesh */

//#include "locator.tcc"

#endif
