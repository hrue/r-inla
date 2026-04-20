#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "fmesher_helpers.h"

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

void map_points_to_mesh(const Mesh &M, const Matrix<double> &points,
                        Matrix<int> &point2T, Matrix<double> &point2bary) {
  int t;
  Point s;
  Point b;
  int the_dimensions[] = {0, 1};
  std::vector<int> dimensions(
      the_dimensions, the_dimensions + sizeof(the_dimensions) / sizeof(int));
  TriangleLocator locator(&M, dimensions, true);

  for (size_t i = 0; i < points.rows(); i++) {
    s[0] = points[i][0];
    s[1] = points[i][1];
    s[2] = points[i][2];
    t = locator.locate(s);
    if (t >= 0) {                      /* Point located. */
      M.barycentric(Dart(M, t), s, b); /* Coordinates relative to
    canonical vertex
    ordering. */
      point2T(i, 0) = t;
      point2bary(i, 0) = b[0];
      point2bary(i, 1) = b[1];
      point2bary(i, 2) = b[2];
    } else { /* Point not found. */
      point2T(i, 0) = -1;
    }
  }
}

void map_points_to_mesh_convex(const Mesh &M, const Matrix<double> &points,
                               Matrix<int> &point2T,
                               Matrix<double> &point2bary) {
  Dart d0(M);
  Dart d;
  Point s;
  Point b;
  for (size_t i = 0; i < points.rows(); i++) {
    s[0] = points[i][0];
    s[1] = points[i][1];
    s[2] = points[i][2];
    d = M.locate_point(Dart(M), s);
    if (!d.isnull()) {                     /* Point located. */
      M.barycentric(Dart(M, d.t()), s, b); /* Coordinates relative to
    canonical vertex
    ordering. */
      point2T(i, 0) = d.t();
      point2bary(i, 0) = b[0];
      point2bary(i, 1) = b[1];
      point2bary(i, 2) = b[2];

      d0 = d; /* Bet on the next point being close. */
    } else {  /* Point not found. */
      point2T(i, 0) = -1;
    }
  }
}

void filter_locations_slow(Matrix<double> &S, Matrix<int> &idx, double cutoff) {
  size_t dim = S.cols();
  size_t idx_next = 0;
  typedef std::list<std::pair<int, Point>> excludedT;
  excludedT excluded;

  FMLOG("Filtering locations." << endl);

  /* Extract "unique" points. */
  double dist;
  Point s = Point(0.0, 0.0, 0.0);
  Point diff = Point(0.0, 0.0, 0.0);
  for (size_t v = 0; v < S.rows(); v++) {
    bool was_excluded = false;
    for (size_t d = 0; d < dim; d++)
      s[d] = S[v][d];
    for (size_t v_try = 0; v_try < idx_next; v_try++) {
      for (size_t d = 0; d < dim; d++)
        diff[d] = S[v_try][d] - s[d];
      if (diff.length() <= cutoff) {
        was_excluded = true;
        excluded.push_back(excludedT::value_type(v, s));
        idx(v, 0) = v_try;
        break;
      }
    }
    if (!was_excluded) {
      for (size_t d = 0; d < dim; d++)
        S(idx_next, d) = s[d];
      idx(v, 0) = idx_next;
      idx_next++;
    }
  }

  FMLOG("All vertices handled." << endl);

  /* Remove excess storage. */
  S.rows(idx_next);

  FMLOG("Excess storage removed." << endl);
  FMLOG("Identifying nearest points." << endl);

  /* Identify nearest nodes for excluded locations. */
  for (excludedT::const_iterator i = excluded.begin(); i != excluded.end();
       i++) {
    size_t v = (*i).first;
    for (size_t d = 0; d < dim; d++)
      s[d] = (*i).second[d];
    double nearest_dist = -1.0;
    int nearest_idx = -1;
    for (size_t v_try = 0; v_try < S.rows(); v_try++) {
      for (size_t d = 0; d < dim; d++)
        diff[d] = S[v_try][d] - s[d];
      dist = diff.length();
      if ((nearest_idx < 0) || (dist < nearest_dist)) {
        nearest_idx = v_try;
        nearest_dist = dist;
      }
    }
    if (idx(v, 0) != nearest_idx) {
      FMLOG("Excluded vertex " << v << " remapped from " << idx(v, 0) << " to "
                               << nearest_idx << "." << endl);
    }
    idx(v, 0) = nearest_idx;
  }

  FMLOG("Done identifying nearest points." << endl);
}

void filter_locations(Matrix<double> &S, Matrix<int> &idx, double cutoff) {
  int const dim = S.cols();
  int const Nv = S.rows();
  NNLocator nnl(&S, dim);
  int incl_next = 0;
  int excl_next = Nv - 1;
  std::vector<int> remap(Nv); // New node ordering; included first,
  // then excluded.

  FMLOG("Filtering locations." << endl);

  for (size_t v = 0; v < size_t(Nv); v++) {
    remap[v] = -1;
  }

  NNLocator::iterator nniter;

  FMLOG("Identify 'unique' points." << endl);
  for (size_t v = 0; v < size_t(Nv); v++) {
    nniter = nnl(S[v], cutoff);
    if (nniter != nnl.end()) {
      // Exclude node
      remap[excl_next] = v;
      idx(v, 0) = excl_next;
      --excl_next;
    } else {
      // Include node
      nnl.insert(nniter, v); // Hint position nniter
      remap[incl_next] = v;
      idx(v, 0) = incl_next;
      ++incl_next;
    }
  }
  FMLOG("All vertices handled." << endl);

  FMLOG("Identifying nearest points for excluded locations." << endl);
  for (size_t v = Nv; v > size_t(incl_next);) {
    --v;
    nniter = nnl(S[remap[v]]);
    if (nniter == nnl.end()) {
      FMLOG_("Internal error: No nearest neighbour found." << endl);
    }
    idx(remap[v], 0) = idx(nniter->second, 0);
    FMLOG("Excluded vertex " << remap[v] << " remapped to " << idx(remap[v], 0)
                             << "." << endl);
  }
  FMLOG("Done identifying nearest points." << endl);

  FMLOG("Compactify storage from " << Nv << " to " << incl_next << endl);
  for (size_t v = 0; v < size_t(incl_next); ++v) {
    // In-place overwrite allowed since remapping of included nodes
    // was order preserving.  If no filtering done, no need to copy data.
    if (v != size_t(remap[v])) {
      for (size_t d = 0; d < size_t(dim); d++)
        S(v, d) = S[remap[v]][d];
    }
  }
  FMLOG("Compactify storage done." << endl);

  /* Remove excess storage. */
  FMLOG("Remove excess storage." << endl);
  S.rows(incl_next);

  FMLOG("Excess storage removed." << endl);
}

void invalidate_unused_vertex_indices(const Mesh &M, Matrix<int> &idx) {
  for (size_t v = 0; v < idx.rows(); v++) {
    if ((idx(v, 0) >= 0) &&
        ((idx(v, 0) >= int(M.nV())) || (M.VT(idx(v, 0)) == -1))) {
      idx(v, 0) = -1;
    }
  }
}

void remap_vertex_indices(const Matrix<int> &idx, Matrix<int> &matrix) {
  FMLOG("Remapping vertex indices for an index matrix." << endl);
  FMLOG("Index size: " << idx.rows() << ", " << idx.cols() << endl);
  FMLOG("Matrix size: " << matrix.rows() << ", " << matrix.cols() << endl);
  for (size_t i = 0; i < matrix.rows(); i++) {
    for (size_t j = 0; j < matrix.cols(); j++) {
      matrix(i, j) = idx[matrix[i][j]][0];
    }
  }
  FMLOG("Done." << endl);
}

void remap_vertex_indices(const Matrix<int> &idx, constrListT &segm) {
  FMLOG("Remapping vertex indices constraint segments." << endl);
  FMLOG("Index size: " << idx.rows() << ", " << idx.cols() << endl);
  FMLOG("Segment size: " << segm.size() << endl);
  for (constrListT::iterator i = segm.begin(); i != segm.end(); i++) {
    (*i).first.first = idx[(*i).first.first][0];
    (*i).first.second = idx[(*i).first.second][0];
  }
  FMLOG("Done." << endl);
}

void prepare_cdt_input(const Matrix<int> &segm0, const Matrix<int> &segmgrp,
                       constrListT &cdt_segm) {
  int grp = 0; // Init with default group.

  if (segm0.cols() == 1) {
    int v0 = -1;
    int v1 = -1;
    for (size_t i = 0; i < segm0.rows(); i++) {
      v0 = v1;
      v1 = segm0[i][0];
      if (i < segmgrp.rows()) // Update group index, if available
        grp = segmgrp[i][0];
      if ((v0 >= 0) && (v1 >= 0)) {
        cdt_segm.push_back(constrT(v0, v1, grp));
      }
    }
  } else if (segm0.cols() == 2) {
    int v0 = -1;
    int v1 = -1;
    for (size_t i = 0; i < segm0.rows(); i++) {
      v0 = segm0[i][0];
      v1 = segm0[i][1];
      if (i < segmgrp.rows()) // Update group index, if available
        grp = segmgrp[i][0];
      if ((v0 >= 0) && (v1 >= 0)) {
        cdt_segm.push_back(constrT(v0, v1, grp));
      }
    }
  }
}

/*
 loc0: nloc0-by-3
 idx0: nidx0-by-2
 loc1: nloc1-by-2
 idx1: nidx1-by-2
 triangle1: nidx1-by-1
 bary1: nidx1-by-3, idx1[i,0] coordinates within triangle1[i]
 bary2: nidx1-by-3, idx1[i,1] coordinates within triangle1[i]
 origin1: nidx1-by-1
 */
void split_line_segments_on_triangles(
    const Mesh &M, const Matrix<double> &loc0, const Matrix<int> &idx0,
    Matrix<double> &loc1, Matrix<int> &idx1, Matrix<int> &triangle1,
    Matrix<double> &bary1, Matrix<double> &bary2, Matrix<int> &origin1) {
  FMLOG("Mesh M=" << M << endl);
  FMLOG("Mesh M.TT=" << M.TT() << endl);

  FMLOG("Split line segments into subsegments on triangles." << endl);
  FMLOG("Point size: " << loc0.rows() << ", " << loc0.cols() << endl);
  FMLOG("Index size: " << idx0.rows() << ", " << idx0.cols() << endl);
  Matrix<int> *loc_in_tri = new Matrix<int>(loc0.rows(), 1);
  Matrix<double> *bary_in_tri = new Matrix<double>(loc0.rows(), 3);
  DartList dart_trace;
  map_points_to_mesh(M, loc0, *loc_in_tri, *bary_in_tri);

  /* Initialize output structures. */
  loc1.rows(0);
  idx1.rows(0);
  triangle1.rows(0);
  bary1.rows(0);
  bary2.rows(0);
  origin1.rows(0);
  int i_loc_curr = -1;
  int i_idx_curr = -1;

  FMLOG("Number of lines to split: " << idx0.rows() << std::endl);
  for (size_t i = 0; i < idx0.rows(); ++i) {
    FMLOG("Split line nr " << i << ": (" << idx0[i][0] << ", " << idx0[i][1]
                           << ")" << std::endl);
    Dart d(M, (*loc_in_tri)[idx0[i][0]][0]);
    Point s0(loc0[idx0[i][0]]);
    Point s1(loc0[idx0[i][1]]);
    FMLOG("Tracing path between points" << std::endl);
    FMLOG("s0=" << s0 << std::endl);
    FMLOG("s1=" << s1 << std::endl);
    FMLOG("Starting dart " << d << endl);

    dart_trace.clear();
    DartPair endpoints(M.trace_path(s0, s1, d, &dart_trace));
    FMLOG("Trace:" << endl << dart_trace << std::endl)

    Point b1;
    Point b2;
    Point s_curr(s0);
    Point s_next(s0); /* Initialise the first sub-segment */

    /* Add the first point */
    ++i_loc_curr;
    for (size_t di = 0; di < 3; di++) {
      loc1(i_loc_curr, di) = s_next[di];
    }

    /* Middle sub-segments */
    for (DartList::const_iterator dti(dart_trace.begin());
         dti != dart_trace.end(); ++dti) {
      FMLOG("Making middle subsegment, split on" << endl
                                                 << " " << *dti << std::endl);
      s_curr = s_next;
      FMLOG("Line to split:" << endl
                             << " " << s_curr << endl
                             << " " << s1 << endl);
      FMLOG("Edge to split on:" << endl
                                << " " << M.S((*dti).v()) << endl
                                << " " << M.S((*dti).vo()) << endl);
      M.edgeIntersection(s_curr, s1, M.S((*dti).v()), M.S((*dti).vo()), s_next);
      FMLOG("Split result = " << s_next << endl);
      M.barycentric(*dti, s_curr, b1);
      M.barycentric(*dti, s_next, b2);
      //
      ++i_idx_curr;
      idx1(i_idx_curr, 0) = i_loc_curr;
      ++i_loc_curr;
      idx1(i_idx_curr, 1) = i_loc_curr;
      for (size_t di = 0; di < 3; di++) {
        loc1(i_loc_curr, di) = s_next[di];
        bary1(i_idx_curr, di) = b1[di];
        bary2(i_idx_curr, di) = b2[di];
      }
      origin1(i_idx_curr, 0) = i;
      triangle1(i_idx_curr, 0) = (*dti).t();
    }
    /* Final sub-segment, or both points in the same triangle */
    if (!endpoints.second.isnull()) {
      FMLOG("Making final subsegment." << std::endl);
      s_curr = s_next;
      s_next = s1;
      if (dart_trace.size() == 0) {
        FMLOG("Staying in initial triangle." << std::endl);
        M.barycentric(endpoints.first, s_curr, b1);
        M.barycentric(endpoints.first, s_next, b2);
      } else {
        FMLOG("Moving to final triangle." << std::endl);
        M.barycentric(endpoints.second, s_curr, b1);
        M.barycentric(endpoints.second, s_next, b2);
      }
      //
      ++i_idx_curr;
      idx1(i_idx_curr, 0) = i_loc_curr;
      ++i_loc_curr;
      idx1(i_idx_curr, 1) = i_loc_curr;
      for (size_t di = 0; di < 3; di++) {
        loc1(i_loc_curr, di) = s_next[di];
        bary1(i_idx_curr, di) = b1[di];
        bary2(i_idx_curr, di) = b2[di];
      }
      origin1(i_idx_curr, 0) = i;
      triangle1(i_idx_curr, 0) = endpoints.second.t();
    }
  }
  delete bary_in_tri;
  delete loc_in_tri;
  FMLOG("Done." << endl);
}
