#ifndef _FMESH_MESHER_HELPERS_
#define _FMESH_MESHER_HELPERS_ 1

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "fmesher_debuglog.h"
#ifndef FMESHER_WITH_R
#include "cmdline.h"
#endif
#include "basis.h"
#include "ioutils.h"
#include "locator.h"
#include "mesh.h"
#include "meshc.h"

using std::endl;
using std::ifstream;
using std::ios;
using std::ofstream;
using std::string;

using fmesh::constrListT;
using fmesh::constrMetaT;
using fmesh::constrT;
using fmesh::Dart;
using fmesh::DartList;
using fmesh::DartPair;
using fmesh::Int3;
using fmesh::Int3Raw;
using fmesh::IOHelper;
using fmesh::IOHelperM;
using fmesh::IOHelperSM;
using fmesh::Matrix;
using fmesh::Matrix3double;
using fmesh::MatrixC;
using fmesh::Mesh;
using fmesh::MeshC;
using fmesh::Point;
using fmesh::PointRaw;
using fmesh::SparseMatrix;
using fmesh::TriangleLocator;
using fmesh::Vector3;
using fmesh::vertexListT;

template <class T>
void print_M(string filename, const Matrix<T> &M,
             fmesh::IOMatrixtype matrixt = fmesh::IOMatrixtype_general) {
  M.save(filename, matrixt);
}

template <class T>
void print_SM(string filename, const SparseMatrix<T> &M,
              fmesh::IOMatrixtype matrixt = fmesh::IOMatrixtype_general) {
  M.save(filename, matrixt);
}

template <class T>
void print_M_old(string filename, const Matrix<T> &M,
                 fmesh::IOMatrixtype matrixt = fmesh::IOMatrixtype_general) {
  M.save_ascii_2009(filename, matrixt);
}

template <class T>
void print_SM_old(string filename, const SparseMatrix<T> &M,
                  fmesh::IOMatrixtype matrixt = fmesh::IOMatrixtype_general) {
  M.save_ascii_2009(filename, matrixt);
}

template <typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> v) {
  int last = v.size() - 1;
  out << "[";
  for (int i = 0; i < last; i++)
    out << v[i] << ", ";
  out << v[last] << "]";
  return out;
}

void map_points_to_mesh(const Mesh &M, const Matrix<double> &points,
                        Matrix<int> &point2T, Matrix<double> &point2bary);
void map_points_to_mesh_convex(const Mesh &M, const Matrix<double> &points,
                               Matrix<int> &point2T,
                               Matrix<double> &point2bary);
void filter_locations_slow(Matrix<double> &S, Matrix<int> &idx, double cutoff);

class NNLocator {
  std::multimap<double, size_t> search_map_;
  Matrix<double> const *S_;
  int _dim;

public:
  NNLocator(Matrix<double> *S, int dim) : search_map_(), S_(S), _dim(dim) {}

public:
  double distance2(double const *point, int v) {
    double diff;
    double dist = 0.0;
    for (int d = 0; d < _dim; ++d) {
      diff = point[d] - (*S_)[v][d];
      dist += diff * diff;
    }
    return dist;
  };
  double distance(double const *point, int v) {
    return std::sqrt(distance2(point, v));
  };

  typedef std::multimap<double, size_t>::value_type value_type;
  typedef std::multimap<double, size_t>::iterator iterator;
  typedef std::multimap<double, size_t>::const_iterator const_iterator;
  typedef std::multimap<double, size_t>::reverse_iterator reverse_iterator;
  typedef std::multimap<double, size_t>::const_reverse_iterator
      const_reverse_iterator;

  iterator insert(int idx) {
    return search_map_.insert(value_type((*S_)[idx][0], idx));
  };
  iterator insert(iterator position, int idx) {
    return search_map_.insert(position, value_type((*S_)[idx][0], idx));
  };
  iterator begin() { return search_map_.begin(); };
  iterator end() { return search_map_.end(); };

  // Find nearest neighbour, optionally with distance <= bound
  iterator find_nn_bounded(double const *point, bool have_bound,
                           double distance2_bound) {
    iterator iter, start;
    double dist;
    bool found = false; // true if we've found at least one neighbour.
    double shortest_dist = -1.0;
    iterator found_iter(search_map_.end()); // pointer to the closest
    // found neighbour

    size_t const size = search_map_.size();
    if (size == 0) {
      return search_map_.end();
    } else if (size == 1) {
      found_iter = search_map_.begin();
      if (have_bound &&
          distance2(point, found_iter->second) > distance2_bound) {
        return search_map_.end();
      }
      return found_iter;
    }

    start = search_map_.lower_bound(point[0]);
    // Handle boundary case:
    // If max < point, then lower=end and upper=end, so we need to
    // skip the forward part and only run backward.
    bool forward = (start != found_iter);
    bool outside_bound = false;

    iter = start;
    while ((forward && iter != search_map_.end()) ||
           (!forward && iter != search_map_.begin())) {
      if (!forward) {
        --iter;
      }
      if (found || have_bound) {
        // Check upper bound first
        dist = iter->first - point[0];
        if ((forward && dist > 0.0) || (!forward && dist < 0.0)) {
          dist = dist * dist;
          if ((found && (dist >= shortest_dist)) ||
              (have_bound && (dist > distance2_bound))) {
            outside_bound = true;
          }
        }
      }
      if (!outside_bound) {
        dist = distance2(point, iter->second);
        FMLOG("distance2 = " << dist << endl);
        FMLOG("found = " << found << ", shortest_dist = " << shortest_dist
                         << endl);
        FMLOG("have_bound = " << have_bound << ", distance2_bound = "
                              << distance2_bound << endl);
        if ((!found || (dist < shortest_dist)) &&
            (!have_bound || (dist <= distance2_bound))) {
          found = true;
          found_iter = iter;
          shortest_dist = dist;
          FMLOG("shortest updated" << endl);
        } else {
          FMLOG("no action" << endl);
        }
      }
      if (forward) {
        ++iter;
      }
      if ((iter == search_map_.end()) || (forward && outside_bound)) {
        outside_bound = false;
        forward = false;
        iter = start;
        FMLOG("reverse" << endl);
      } else if (!forward && outside_bound) {
        break;
      }
    }

    FMLOG("Finished. found = " << found << endl);
    return found_iter;
  };
  iterator operator()(double const *point) {
    return find_nn_bounded(point, false, 0.0);
  };
  iterator operator()(double const *point, double cutoff) {
    return find_nn_bounded(point, true, cutoff * cutoff);
  };
};

void filter_locations(Matrix<double> &S, Matrix<int> &idx, double cutoff);

void invalidate_unused_vertex_indices(const Mesh &M, Matrix<int> &idx);
void remap_vertex_indices(const Matrix<int> &idx, Matrix<int> &matrix);
void remap_vertex_indices(const Matrix<int> &idx, constrListT &segm);

void prepare_cdt_input(const Matrix<int> &segm0, const Matrix<int> &segmgrp,
                       constrListT &cdt_segm);

void split_line_segments_on_triangles(
    const Mesh &M, const Matrix<double> &loc0, const Matrix<int> &idx0,
    Matrix<double> &loc1, Matrix<int> &idx1, Matrix<int> &triangle1,
    Matrix<double> &bary1, Matrix<double> &bary2, Matrix<int> &origin1);

#endif
