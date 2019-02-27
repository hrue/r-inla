#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "fmesher.hh"

#ifndef WHEREAMI
#define WHEREAMI __FILE__ << "(" << __LINE__ << ")\t"
#endif

#ifndef LOG_
#define LOG_(msg) cout << WHEREAMI << msg;
#endif
#ifndef LOG
#ifdef DEBUG
#define LOG(msg) LOG_(msg)
#else
#define LOG(msg)
#endif
#endif

using std::ios;
using std::ifstream;
using std::ofstream;
using std::string;
using std::cin;
using std::cout;
using std::endl;

using fmesh::Dart;
using fmesh::DartPair;
using fmesh::DartList;
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
using fmesh::Vector3;
using fmesh::constrMetaT;
using fmesh::constrT;
using fmesh::constrListT;
using fmesh::vertexListT;
using fmesh::TriangleLocator;

const bool useVT = true;
const bool useTTi = true;
bool useX11 = false;
const bool useX11text = false;
double x11_delay_factor = 1.0;
double x11_zoom[4];

MatrixC matrices;


template <class T>
void print_M(string filename,
	     const Matrix<T>& M,
	     fmesh::IOMatrixtype matrixt = fmesh::IOMatrixtype_general)
{
  M.save(filename,matrixt);
}

template <class T>
void print_SM(string filename,
	      const SparseMatrix<T>& M,
	      fmesh::IOMatrixtype matrixt = fmesh::IOMatrixtype_general)
{
  M.save(filename,matrixt);
}

template <class T>
void print_M_old(string filename,
		 const Matrix<T>& M,
		 fmesh::IOMatrixtype matrixt = fmesh::IOMatrixtype_general)
{
  M.save_ascii_2009(filename,matrixt);
}

template <class T>
void print_SM_old(string filename,
		  const SparseMatrix<T>& M,
		  fmesh::IOMatrixtype matrixt = fmesh::IOMatrixtype_general)
{
  M.save_ascii_2009(filename,matrixt);
}




void map_points_to_mesh(const Mesh& M,
			const Matrix<double>& points,
			Matrix<int>& point2T,
			Matrix<double>& point2bary)
{
  int t;
  Point s;
  Point b;
  int the_dimensions[] = {0,1};
  std::vector<int> dimensions(the_dimensions,
			      the_dimensions +
			      sizeof(the_dimensions) / sizeof(int) );
  TriangleLocator locator(&M, dimensions, true);

  for (size_t i=0; i<points.rows(); i++) {
    s[0] = points[i][0];
    s[1] = points[i][1];
    s[2] = points[i][2];
    t = locator.locate(s);
    if (t>=0) { /* Point located. */
      M.barycentric(Dart(M,t),s,b); /* Coordinates relative to
					   canonical vertex
					   ordering. */
      point2T(i,0) = t;
      point2bary(i,0) = b[0];
      point2bary(i,1) = b[1];
      point2bary(i,2) = b[2];
    } else { /* Point not found. */
      point2T(i,0) = -1;
    }
  }
}


void map_points_to_mesh_convex(const Mesh& M,
			       const Matrix<double>& points,
			       Matrix<int>& point2T,
			       Matrix<double>& point2bary)
{
  Dart d0(M);
  Dart d;
  Point s;
  Point b;
  for (size_t i=0; i<points.rows(); i++) {
    s[0] = points[i][0];
    s[1] = points[i][1];
    s[2] = points[i][2];
    d = M.locate_point(Dart(M),s);
    if (!d.isnull()) { /* Point located. */
      M.barycentric(Dart(M,d.t()),s,b); /* Coordinates relative to
					   canonical vertex
					   ordering. */
      point2T(i,0) = d.t();
      point2bary(i,0) = b[0];
      point2bary(i,1) = b[1];
      point2bary(i,2) = b[2];

      d0 = d; /* Bet on the next point being close. */
    } else { /* Point not found. */
      point2T(i,0) = -1;
    }
  }
}


void filter_locations_slow(Matrix<double>& S,
			   Matrix<int>& idx,
			   double cutoff)
{
  size_t dim = S.cols();
  size_t idx_next = 0;
  typedef std::list< std::pair<int, Point> > excludedT;
  excludedT excluded;

  LOG("Filtering locations." << endl);

  /* Extract "unique" points. */
  double dist;
  Point s = Point(0.0, 0.0, 0.0);
  Point diff = Point(0.0, 0.0, 0.0);
  for (size_t v=0; v<S.rows(); v++) {
    bool was_excluded = false;
    for (size_t d=0; d<dim; d++)
      s[d] = S[v][d];
    for (size_t v_try=0; v_try<idx_next; v_try++) {
      for (size_t d=0; d<dim; d++)
	diff[d] = S[v_try][d]-s[d];
      if (diff.length() <= cutoff) {
	was_excluded = true;
	excluded.push_back(excludedT::value_type(v,s));
	idx(v,0) = v_try;
	break;
      }
    }
    if (!was_excluded) {
      for (size_t d=0; d<dim; d++)
	S(idx_next,d) = s[d];
      idx(v,0) = idx_next;
      idx_next++;
    }
  }

  LOG("All vertices handled." << endl);

  /* Remove excess storage. */
  S.rows(idx_next);

  LOG("Excess storage removed." << endl);
  LOG("Identifying nearest points." << endl);

  /* Identify nearest nodes for excluded locations. */
  for (excludedT::const_iterator i=excluded.begin();
       i != excluded.end();
       i++) {
    size_t v = (*i).first;
    for (size_t d=0; d<dim; d++)
      s[d] = (*i).second[d];
    double nearest_dist = -1.0;
    int nearest_idx = -1;
    for (size_t v_try=0; v_try<S.rows(); v_try++) {
      for (size_t d=0; d<dim; d++)
	diff[d] = S[v_try][d]-s[d];
      dist = diff.length();
      if ((nearest_idx<0) || (dist<nearest_dist)) {
	nearest_idx = v_try;
	nearest_dist = dist;
      }
    }
    if (idx(v,0) != nearest_idx) {
      LOG("Excluded vertex "
	   << v << " remapped from "
	   << idx(v,0) << " to "
	   << nearest_idx << "."
	   << endl);
    }
    idx(v,0) = nearest_idx;
  }

  LOG("Done identifying nearest points." << endl);
}




template<typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T> v) {
    int last = v.size() - 1;
    out << "[";
    for(int i = 0; i < last; i++)
        out << v[i] << ", ";
    out << v[last] << "]";
    return out;
}



class NNLocator
{
  std::multimap<double, size_t> search_map_;
  Matrix<double> const * S_;
  int _dim;
public:
  NNLocator(Matrix<double> * S,
	    int dim) :
    search_map_(), S_(S), _dim(dim) {}
public:
  double distance2(double const * point, int v) {
    double diff;
    double dist = 0.0;
    for (int d=0; d < _dim; ++d) {
      diff = point[d]-(*S_)[v][d];
      dist += diff*diff;
    }
    return dist;
  };
  double distance(double const * point, int v) {
    return std::sqrt(distance2(point, v));
  };

  typedef std::multimap<double, size_t>::value_type value_type;
  typedef std::multimap<double, size_t>::iterator iterator;
  typedef std::multimap<double, size_t>::const_iterator const_iterator;
  typedef std::multimap<double, size_t>::reverse_iterator reverse_iterator;
  typedef std::multimap<double, size_t>::const_reverse_iterator const_reverse_iterator;

  iterator insert(int idx) {
    return search_map_.insert(value_type((*S_)[idx][0], idx));
  };
  iterator insert(iterator position, int idx) {
    return search_map_.insert(position, value_type((*S_)[idx][0], idx));
  };
  iterator begin() {
    return search_map_.begin();
  };
  iterator end() {
    return search_map_.end();
  };

  // Find nearest neighbour, optionally with distance <= bound
  iterator find_nn_bounded(double const * point,
			   bool have_bound,
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
	if ((forward && dist > 0.0) ||
	    (!forward && dist < 0.0)) {
	  dist = dist*dist;
	  if ((found && (dist >= shortest_dist)) ||
	      (have_bound && (dist > distance2_bound))) {
	    outside_bound = true;
	  }
	}
      }
      if (!outside_bound) {
	dist = distance2(point, iter->second);
	LOG("distance2 = " << dist << endl);
	LOG("found = " << found <<
	    ", shortest_dist = " << shortest_dist << endl);
	LOG("have_bound = " << have_bound <<
	    ", distance2_bound = " << distance2_bound << endl);
	if ((!found || (dist < shortest_dist)) &&
	    (!have_bound || (dist <= distance2_bound))) {
	  found = true;
	  found_iter = iter;
	  shortest_dist = dist;
	  LOG("shortest updated" << endl);
	} else {
	  LOG("no action" << endl);
	}
      }
      if (forward) {
	++iter;
      }
      if ((iter == search_map_.end()) || (forward && outside_bound)) {
	outside_bound = false;
	forward = false;
	iter = start;
	LOG("reverse" << endl);
      } else if (!forward && outside_bound) {
	break;
      }
    }

    LOG("Finished. found = " << found << endl);
    return found_iter;
  };
  iterator operator() (double const * point) {
    return find_nn_bounded(point, false, 0.0);
  };
  iterator operator() (double const * point,
		       double cutoff) {
    return find_nn_bounded(point, true, cutoff*cutoff);
  };
};


void filter_locations(Matrix<double>& S,
		      Matrix<int>& idx,
		      double cutoff)
{
  int const dim = S.cols();
  int const Nv = S.rows();
  NNLocator nnl(&S, dim);
  int incl_next = 0;
  int excl_next = Nv-1;
  std::vector<int> remap(Nv); // New node ordering; included first,
			      // then excluded.

  LOG("Filtering locations." << endl);

  for (size_t v=0; v < size_t(Nv); v++) {
    remap[v] = -1;
  }

  NNLocator::iterator nniter;

  LOG("Identify 'unique' points." << endl);
  for (size_t v=0; v < size_t(Nv); v++) {
    nniter = nnl(S[v], cutoff);
    if (nniter != nnl.end()) {
      // Exclude node
      remap[excl_next] = v;
      idx(v,0) = excl_next;
      --excl_next;
    } else {
      // Include node
      nnl.insert(nniter, v); // Hint position nniter
      remap[incl_next] = v;
      idx(v,0) = incl_next;
      ++incl_next;
    }
  }
  LOG("All vertices handled." << endl);

  LOG("Identifying nearest points for excluded locations." << endl);
  for (size_t v=Nv; v > size_t(incl_next);) {
    --v;
    nniter = nnl(S[remap[v]]);
    if (nniter == nnl.end()) {
      cout << "Internal error: No nearest neighbour found." << endl;
    }
    idx(remap[v],0) = idx(nniter->second,0);
    LOG("Excluded vertex "
    	<< remap[v] << " remapped to "
    	<< idx(remap[v],0) << "."
    	<< endl);
  }
  LOG("Done identifying nearest points." << endl);

  LOG("Compactify storage from " << Nv << " to " << incl_next << endl);
  for (size_t v=0; v < size_t(incl_next); ++v) {
    // In-place overwrite allowed since remapping of included nodes
    // was order preserving.  If no filtering done, no need to copy data.
    if (v != size_t(remap[v])) {
      for (size_t d=0; d < size_t(dim); d++)
	S(v,d) = S[remap[v]][d];
    }
  }
  LOG("Compactify storage done." << endl);

  /* Remove excess storage. */
  LOG("Remove excess storage." << endl);
  S.rows(incl_next);

  LOG("Excess storage removed." << endl);
}


void invalidate_unused_vertex_indices(const Mesh& M,
				      Matrix<int>& idx)
{
  for (size_t v=0; v<idx.rows(); v++) {
    if ((idx(v,0)>=0) &&
	((idx(v,0) >= int(M.nV())) || (M.VT(idx(v,0))==-1))) {
      idx(v,0) = -1;
    }
  }
}

void remap_vertex_indices(const Matrix<int>& idx, Matrix<int>& matrix)
{
  LOG("Remapping vertex indices for an index matrix." << endl);
  LOG("Index size: " << idx.rows() << ", " << idx.cols() << endl);
  LOG("Matrix size: " << matrix.rows() << ", " << matrix.cols() << endl);
  for (size_t i=0; i<matrix.rows(); i++) {
    for (size_t j=0; j<matrix.cols(); j++) {
      matrix(i,j) = idx[matrix[i][j]][0];
    }
  }
  LOG("Done." << endl);
}

void remap_vertex_indices(const Matrix<int>& idx, constrListT& segm)
{
  LOG("Remapping vertex indices constraint segments." << endl);
  LOG("Index size: " << idx.rows() << ", " << idx.cols() << endl);
  LOG("Segment size: " << segm.size() << endl);
  for (constrListT::iterator i=segm.begin(); i != segm.end(); i++) {
    (*i).first.first = idx[(*i).first.first][0];
    (*i).first.second = idx[(*i).first.second][0];
  }
  LOG("Done." << endl);
}



void prepare_cdt_input(const Matrix<int>& segm0,
		       const Matrix<int>& segmgrp,
		       constrListT& cdt_segm)
{
  int grp = 0; // Init with default group.

  if (segm0.cols()==1) {
    int v0 = -1;
    int v1 = -1;
    for (size_t i=0; i < segm0.rows(); i++) {
      v0 = v1;
      v1 = segm0[i][0];
      if (i<segmgrp.rows()) // Update group index, if available
	grp = segmgrp[i][0];
      if ((v0>=0) && (v1>=0)) {
	cdt_segm.push_back(constrT(v0,v1,grp));
      }
    }
  } else if (segm0.cols()==2) {
    int v0 = -1;
    int v1 = -1;
    for (size_t i=0; i < segm0.rows(); i++) {
      v0 = segm0[i][0];
      v1 = segm0[i][1];
      if (i<segmgrp.rows()) // Update group index, if available
	grp = segmgrp[i][0];
      if ((v0>=0) && (v1>=0)) {
	cdt_segm.push_back(constrT(v0,v1,grp));
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
void split_line_segments_on_triangles(const Mesh& M,
				      const Matrix<double>& loc0,
				      const Matrix<int>& idx0,
				      Matrix<double>& loc1,
				      Matrix<int>& idx1,
				      Matrix<int>& triangle1,
				      Matrix<double>& bary1,
				      Matrix<double>& bary2,
				      Matrix<int>& origin1)
{
  LOG("Mesh M=" << M << endl);
  LOG("Mesh M.TT=" << M.TT() << endl);
  
  LOG("Split line segments into subsegments on triangles." << endl);
  LOG("Point size: " << loc0.rows() << ", " << loc0.cols() << endl);
  LOG("Index size: " << idx0.rows() << ", " << idx0.cols() << endl);
  Matrix<int>* loc_in_tri = new Matrix<int>(loc0.rows(), 1);
  Matrix<double>* bary_in_tri = new Matrix<double>(loc0.rows(), 3);
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

  LOG("Number of lines to split: " << idx0.rows() << std::endl);
  for (size_t i=0; i < idx0.rows(); ++i) {
    LOG("Split line nr " << i << ": (" <<
	 idx0[i][0] << ", " <<
	 idx0[i][1] << ")" << std::endl);
    Dart d(M, (*loc_in_tri)[idx0[i][0]][0]);
    Point s0(loc0[idx0[i][0]]);
    Point s1(loc0[idx0[i][1]]);
    LOG("Tracing path between points" << std::endl);
    LOG("s0=" << s0 << std::endl);
    LOG("s1=" << s1 << std::endl);
    LOG("Starting dart " << d << endl);

    dart_trace.clear();
    DartPair endpoints(M.trace_path(s0, s1, d, &dart_trace));
    LOG("Trace:" << endl << dart_trace << std::endl)

    Point b1;
    Point b2;
    Point s_curr(s0);
    Point s_next(s0); /* Initialise the first sub-segment */

    /* Add the first point */
    ++i_loc_curr;
    for (size_t di=0; di<3; di++) {
      loc1(i_loc_curr,di) = s_next[di];
    }

    /* Middle sub-segments */
    for (DartList::const_iterator dti(dart_trace.begin());
	 dti != dart_trace.end();
	 ++dti) {
      LOG("Making middle subsegment, split on" << endl <<
	   " " << *dti << std::endl);
      s_curr = s_next;
      LOG("Line to split:" << endl << " " <<
	  s_curr << endl << " " << s1 << endl);
      LOG("Edge to split on:" << endl << " " <<
	  M.S((*dti).v()) << endl << " " <<
	  M.S((*dti).vo()) << endl);
      M.edgeIntersection(s_curr, s1,
			 M.S((*dti).v()), M.S((*dti).vo()),
			 s_next);
      LOG("Split result = " << s_next << endl);
      M.barycentric(*dti, s_curr, b1);
      M.barycentric(*dti, s_next, b2);
      //
      ++i_idx_curr;
      idx1(i_idx_curr, 0) = i_loc_curr;
      ++i_loc_curr;
      idx1(i_idx_curr, 1) = i_loc_curr;
      for (size_t di=0; di<3; di++) {
	loc1(i_loc_curr, di) = s_next[di];
	bary1(i_idx_curr, di) = b1[di];
	bary2(i_idx_curr, di) = b2[di];
      }
      origin1(i_idx_curr, 0) = i;
      triangle1(i_idx_curr, 0) = (*dti).t();
    }
    /* Final sub-segment, or both points in the same triangle */
    if (!endpoints.second.isnull()) {
      LOG("Making final subsegment." << std::endl);
      s_curr = s_next;
      s_next = s1;
      if (dart_trace.size() == 0) {
	LOG("Staying in initial triangle." << std::endl);
	M.barycentric(endpoints.first, s_curr, b1);
	M.barycentric(endpoints.first, s_next, b2);
      } else {
	LOG("Moving to final triangle." << std::endl);
	M.barycentric(endpoints.second, s_curr, b1);
	M.barycentric(endpoints.second, s_next, b2);
      }
      //
      ++i_idx_curr;
      idx1(i_idx_curr, 0) = i_loc_curr;
      ++i_loc_curr;
      idx1(i_idx_curr, 1) = i_loc_curr;
      for (size_t di=0; di<3; di++) {
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
  LOG("Done." << endl);
}






int main(int argc, char* argv[])
{
  
  gengetopt_args_info args_info;
  struct cmdline_params params;
  
  LOG("checkpoint 1." << std::endl);
  
  cmdline_init(&args_info);
  cmdline_params_init(&params);
  
  LOG("checkpoint 2." << std::endl);
    
  /* call the command line parser */
  if (cmdline_ext(argc, argv, &args_info, &params) != 0) {
    cmdline_free(&args_info);
    return 1;
  }

  LOG("checkpoint 3." << std::endl);

  /* Read an optional config file, but don't override given options */
  if (args_info.config_given) {
    params.initialize = 0;
    params.override = 0;
    /* Call the config file parser */
    if (cmdline_config_file(args_info.config_arg, &args_info, &params) != 0) {
      cmdline_free(&args_info);
      return 1;
    }
  }

  LOG("checkpoint 4." << std::endl);

  if (args_info.dump_config_given)
    cmdline_dump(stdout,&args_info);
  
  std::vector<string> input_s0_names;
  string input_tv0_name = "-";

  LOG("checkpoint 5." << std::endl);
  
  if (args_info.input_given>0)
    input_s0_names.push_back(string(args_info.input_arg[0]));
  if (args_info.input_given>1)
    input_tv0_name = string(args_info.input_arg[1]);
  for (size_t i=2; i < args_info.input_given; i++) {
    input_s0_names.push_back(string(args_info.input_arg[i]));
  }

  LOG("checkpoint 6." << std::endl);

  std::vector<string> quality_names;
  for (size_t i=0; i < args_info.quality_given; i++) {
    quality_names.push_back(string(args_info.quality_arg[i]));
  }

  LOG("checkpoint 7." << std::endl);

  std::vector<string> boundary_names;
  std::vector<string> boundarygrp_names;
  for (size_t i=0; i < args_info.boundary_given; i++) {
    boundary_names.push_back(string(args_info.boundary_arg[i]));
    if (i<args_info.boundarygrp_given)
      boundarygrp_names.push_back(string(args_info.boundarygrp_arg[i]));
    else
      boundarygrp_names.push_back(boundary_names[i]+"grp");
  }

  LOG("checkpoint 8." << std::endl);

  std::vector<string> interior_names;
  std::vector<string> interiorgrp_names;
  for (size_t i=0; i < args_info.interior_given; i++) {
    interior_names.push_back(string(args_info.interior_arg[i]));
    if (i<args_info.interiorgrp_given)
      interiorgrp_names.push_back(string(args_info.interiorgrp_arg[i]));
    else
      interiorgrp_names.push_back(interior_names[i]+"grp");
  }

  LOG("checkpoint 9." << std::endl);

  std::vector<string> aniso_names;
  for (size_t i=0; i < args_info.aniso_given; i++) {
    aniso_names.push_back(string(args_info.aniso_arg[i]));
  }

  std::vector<string> splitlines_names;
  for (size_t i=0; i < args_info.splitlines_given; i++) {
    splitlines_names.push_back(string(args_info.splitlines_arg[i]));
  }

  double cutoff = 1.0e-12;
  if (args_info.cutoff_given>0)
    cutoff = args_info.cutoff_arg;
  double sphere_tolerance = 1.0e-7;
  if (args_info.spheretolerance_given>0)
    sphere_tolerance = args_info.spheretolerance_arg;
  
  int cet_sides = 8;
  double cet_margin = -0.1;
  if (args_info.cet_given>0)
    cet_sides = args_info.cet_arg[0];
  if (args_info.cet_given>1)
    cet_margin = args_info.cet_arg[1];
  
  double rcdt_min_angle = 21;
  double rcdt_big_limit_auto_default = -1.0;
  Matrix<double> rcdt_big_limit_defaults;
  rcdt_big_limit_defaults(0,0) = -0.5;
  for (size_t i=1; (i+2 < args_info.input_given); i++) {
    rcdt_big_limit_defaults(i,0) = -0.5;
  }
  if ((args_info.rcdt_given>0) && (args_info.rcdt_arg[0] != 0))
    rcdt_min_angle = args_info.rcdt_arg[0];
  if (args_info.rcdt_given>1)
    rcdt_big_limit_auto_default = args_info.rcdt_arg[1];
  for (size_t i=0; (i+2 < args_info.rcdt_given); i++) {
    rcdt_big_limit_defaults(i,0) = args_info.rcdt_arg[i+2];
  }

  int rcdt_max_n0 = -1;
  int rcdt_max_n1 = -1;
  if (args_info.max_n0_given > 0)
    rcdt_max_n0 = args_info.max_n0_arg;
  if (args_info.max_n1_given > 0)
    rcdt_max_n1 = args_info.max_n1_arg;

  useX11 = (args_info.x11_given>0) && (args_info.x11_arg>=0);
  x11_delay_factor = args_info.x11_arg;
  if (args_info.x11_zoom_given==4) {
    for (size_t i=0; i<4; i++)
      x11_zoom[i] = args_info.x11_zoom_arg[i];
  } else if (args_info.x11_zoom_given==3) {
    double xmid = args_info.x11_zoom_arg[0];
    double ymid = args_info.x11_zoom_arg[1];
    double margin = args_info.x11_zoom_arg[2];
    x11_zoom[0] = xmid-margin;
    x11_zoom[1] = xmid+margin;
    x11_zoom[2] = ymid-margin;
    x11_zoom[3] = ymid+margin;
  }

  /*
  cout << "CET given:\t" << args_info.cet_given << endl;
  cout << "RCDT given:\t"
       << args_info.rcdt_given << " "
       << args_info.rcdt_min << " "
       << args_info.rcdt_max << " "
       << args_info.rcdt_arg[0] << " "
       << endl;
  if (args_info.boundary_given) {
    cout << "Boundary given:\t"
	 << args_info.boundary_given << " "
	 << args_info.boundary_arg << " "
	 << &(args_info.boundary_arg[0]) << " "
      //	 << string(args_info.boundary_arg[0]) << " "
	 << endl;
  }
  cout << "X11 given:\t"
       << args_info.x11_given << " "
       << args_info.x11_arg << " "
       << endl;

  cout << "CET sides:\t" << cet_sides << endl;
  cout << "CET margin:\t" << cet_margin << endl;
  cout << "RCDT mininmum angle:\t" << rcdt_min_angle << endl;
  cout << "RCDT maximum edge length:\t" << rcdt_big_limit << endl;
  cout << "RCDT maximum edge lengths:\t" << rcdt_big_limits << endl;
  cout << "RCDT maximum n0:\t" << rcdt_max_n0 << endl;
  cout << "RCDT maximum n1:\t" << rcdt_max_n1 << endl;
  cout << "X11 delay factor:\t" << x11_delay_factor << endl;
  */

  LOG("IOprefix init." << std::endl);

  string iprefix("-");
  string oprefix("-");
  if (args_info.inputs_num>0) {
    iprefix = string(args_info.inputs[0]);
  }
  if (iprefix != "-") {
    oprefix = iprefix;
  }
  if (args_info.inputs_num>1) {
    oprefix = string(args_info.inputs[1]);
  }
  if ((iprefix=="-") && (args_info.ic_given==0)) {
    /* Nowhere to read general input!
       May be OK, if the input is available in raw format. */
    if (args_info.ir_given==0) {
      /* No input available.
	 May be OK. */
    }
  }
  if ((oprefix=="-") && (args_info.oc_given==0)) {
    /* Nowhere to place output!
       OK; might just want to see the algorithm at work with --x11. */
  }

  LOG("matrix IO init." << std::endl);

  matrices.io(((args_info.io_arg == io_arg_ba) ||
	       (args_info.io_arg == io_arg_bb)),
	      ((args_info.io_arg == io_arg_ab) ||
	       (args_info.io_arg == io_arg_bb)));
  matrices.input_prefix(iprefix);
  matrices.output_prefix(oprefix);
  for (size_t i=0; i<args_info.ic_given; ++i) {
    matrices.input_file(string(args_info.ic_arg[i]));
  }
  if (args_info.oc_given>0) {
    matrices.output_file(string(args_info.oc_arg));
  }
  for (size_t i=0; i+2<args_info.ir_given; i=i+3) {
    matrices.input_raw(string(args_info.ir_arg[i]),
		       string(args_info.ir_arg[i+1]),
		       string(args_info.ir_arg[i+2]));
  }

  LOG("matrix IO inited." << std::endl);

  for (size_t i=0; i<input_s0_names.size(); i++) {
    if (!matrices.load(input_s0_names[i]).active) {
      cout << "Matrix "+input_s0_names[i]+" not found." << endl;
    }
  }
  LOG("s0 input read." << std::endl);
  if ((args_info.globe_given>0) && (args_info.globe_arg>0)) {
    input_s0_names.push_back(string(".globe"));
    matrices.attach(".globe",
		    (Matrix<double>*)fmesh::make_globe_points(args_info.globe_arg),
		    true);
    LOG("globe points added." << std::endl);
  }

  for (size_t i=0; i<quality_names.size(); i++) {
    if (quality_names[i] != "-")
      if (!matrices.load(quality_names[i]).active) {
	cout << "Matrix "+quality_names[i]+" not found." << endl;
	quality_names[i] = "-";
      }
  }
  LOG("quality input read." << std::endl);

  LOG("iS0" << std::endl);

  if (input_s0_names.size()==0) {
    input_s0_names.push_back(string("s0"));
    matrices.attach(string("s0"), new Matrix<double>(3), true);
  }
  Matrix<double>& iS0 = matrices.DD(input_s0_names[0]);
  Matrix<double>* Quality0_ = new Matrix<double>();
  Matrix<double>& Quality0 = *Quality0_;

  /* Join the location matrices */ 
  for (size_t i=0; i < input_s0_names.size(); i++) {
    Matrix<double>& S0_extra = matrices.DD(input_s0_names[i]);
    if (i>0) /* i=0 is already taken care of above. */
      iS0.append(S0_extra);
    if ((i < quality_names.size()) && (quality_names[i] != "-")) {
      int rows = S0_extra.rows();
      Matrix<double>& quality_extra = matrices.DD(quality_names[i]);
      if (i<rcdt_big_limit_defaults.rows())
	for (int r=quality_extra.rows(); r<rows; r++)
	  quality_extra(r,0) = rcdt_big_limit_defaults[i][0];
      else
	for (int r=quality_extra.rows(); r<rows; r++)
	  quality_extra(r,0) = rcdt_big_limit_defaults[0][0];
      quality_extra.rows(rows); /* Make sure we have the right number
				   of rows */
      Quality0.append(quality_extra);
    } else if (i<size_t(rcdt_big_limit_defaults.rows())) {
      int rows = S0_extra.rows();
      Matrix<double> quality_extra(rows,1);
      for (int r=0; r<rows; r++)
	quality_extra(r,0) = rcdt_big_limit_defaults[i][0];
      Quality0.append(quality_extra);
    } else {
      int rows = S0_extra.rows();
      Matrix<double> quality_extra(rows,1);
      for (int r=0; r<rows; r++)
	quality_extra(r,0) = rcdt_big_limit_defaults[0][0];
      Quality0.append(quality_extra);
    }
  }

  /* OK to overwrite any old quality0 */
  matrices.attach(string("quality0"),Quality0_,true);


  LOG("TV0" << std::endl);

  Matrix<int>* TV0 = NULL;
  if (input_tv0_name != "-") {
    if (!matrices.load(input_tv0_name).active) {
      cout << "Matrix "+input_tv0_name+" not found." << endl;
    } else {
      TV0 = &(matrices.DI(input_tv0_name));
    }
  }


  for (size_t i=0; i<boundary_names.size(); i++) {
    if (!matrices.load(boundary_names[i]).active) {
      cout << "Matrix "+boundary_names[i]+" not found." << endl;
    }
    if (!matrices.load(boundarygrp_names[i]).active) {
      // cout << "Matrix "+boundarygrp_names[i]+" not found. Creating." << endl;
      matrices.attach(boundarygrp_names[i],new Matrix<int>(1),true);
      matrices.DI(boundarygrp_names[i])(0,0) = i+1;
    }
  }

  for (size_t i=0; i<interior_names.size(); i++) {
    if (!matrices.load(interior_names[i]).active) {
      cout << "Matrix "+interior_names[i]+" not found." << endl;
    }
    if (!matrices.load(interiorgrp_names[i]).active) {
      cout << "Matrix "+interiorgrp_names[i]+" not found. Creating." << endl;
      matrices.attach(interiorgrp_names[i],new Matrix<int>(1),true);
      matrices.DI(interiorgrp_names[i])(0,0) = i+1;
    }
  }

  for (size_t i=0; i<aniso_names.size(); i++) {
    if (!matrices.load(aniso_names[i]).active) {
      cout << "Matrix "+aniso_names[i]+" not found." << endl;
    }
  }

  for (size_t i=0; i<splitlines_names.size(); i++) {
    if (!matrices.load(splitlines_names[i]).active) {
      cout << "Matrix "+splitlines_names[i]+" not found." << endl;
    }
  }






  constrListT cdt_boundary;
  for (size_t i=0; i<boundary_names.size(); i++) {
    Matrix<int>& boundary0 = matrices.DI(boundary_names[i]);
    Matrix<int>& boundarygrp = matrices.DI(boundarygrp_names[i]);
    prepare_cdt_input(boundary0,boundarygrp,cdt_boundary);
  }

  constrListT cdt_interior;
  for (size_t i=0; i<interior_names.size(); i++) {
    Matrix<int>& interior0 = matrices.DI(interior_names[i]);
    Matrix<int>& interiorgrp = matrices.DI(interiorgrp_names[i]);
    prepare_cdt_input(interior0,interiorgrp,cdt_interior);
  }


  /* Prepare to filter out points at distance not greater than 'cutoff' */
  matrices.attach("idx", new Matrix<int>(iS0.rows(),1), true);
  matrices.output("idx");
  Matrix<int>& idx = matrices.DI("idx").clear();

  /* Only filter if "smorg" option inactive */
  if (args_info.smorg_given==0) {
    filter_locations(iS0, idx, cutoff);

    /* Remap vertex inputreferences */
    if (TV0)
      remap_vertex_indices(idx, *TV0);
    remap_vertex_indices(idx, cdt_boundary);
    remap_vertex_indices(idx, cdt_interior);
  } else {
    for (size_t v=0; v < iS0.rows(); v++) {
      idx(v,0) = v;
    }
  }


  Mesh M(Mesh::Mtype_plane,0,useVT,useTTi);

  LOG("checkpoint 10." << std::endl);

  bool issphere = false;
  bool isflat = false;
  if ((iS0.rows()>0) && (iS0.cols()<2)) {
    /* 1D data. Not implemented */
    LOG_("1D data not implemented." << std::endl);
    return 0;
  } else if (iS0.rows()>0) {
    Matrix3double S0(iS0); /* Make sure we have a Nx3 matrix. */
    M.S_append(S0);
    size_t nV = iS0.rows();

    isflat = (std::fabs(M.S(0)[2]) < 1.0e-10);
    double radius = M.S(0).length();
    issphere = true;
    for (size_t i=1; i<M.nV(); i++) {
      isflat = (isflat && (std::fabs(M.S(i)[2]) < 1.0e-10));
      issphere = (issphere && (std::fabs(M.S(i).length()-radius) < sphere_tolerance));
    }
    if (!isflat) {
      if (issphere) {
	M.type(Mesh::Mtype_sphere);
      } else {
	M.type(Mesh::Mtype_manifold);
      }
    }

#ifndef FMESHER_NO_X
    M.setX11VBigLimit(M.nV());
#endif
      
    if (TV0) {
      M.TV_set(*TV0);
    }
    
    matrices.attach(string("s"),&M.S(),false);
    matrices.attach("tv",&M.TV(),false);
    matrices.output("s").output("tv");
    
    Point mini(M.S(0));
    Point maxi(M.S(0));
    for (size_t v=1; v<M.nV(); v++)
      for (size_t i=0; i<3; i++) {
	mini[i] = (M.S(v)[i] < mini[i] ? M.S(v)[i] : mini[i]);
	maxi[i] = (M.S(v)[i] > maxi[i] ? M.S(v)[i] : maxi[i]);
      }
    Point sz;
    fmesh::Vec::diff(sz,maxi,mini);
    /*
    double diam;
    diam = (sz[1] < sz[0]
	    ? (sz[2] < sz[0] ? sz[0] : sz[2])
	    : (sz[2] < sz[1] ? sz[1] : sz[2]));
    */
    
#ifndef FMESHER_NO_X
    if (useX11) {
      if (issphere) {
	if (args_info.x11_zoom_given==0) {
	  x11_zoom[0] = -1.1;
	  x11_zoom[1] = 1.1;
	  x11_zoom[2] = -1.1;
	  x11_zoom[3] = 1.1;
	}
	M.useX11(true,useX11text,500,500,
		 x11_zoom[0],x11_zoom[1],x11_zoom[2],x11_zoom[3]);
      } else {
	double w0 = maxi[0]-mini[0];
	double w1 = maxi[1]-mini[1];
	double w = (w0 > w1 ? w0 : w1);
	if (args_info.x11_zoom_given==0) {
	  x11_zoom[0] = mini[0]-w*0.2;
	  x11_zoom[1] = maxi[0]+w*0.2;
	  x11_zoom[2] = mini[1]-w*0.2;
	  x11_zoom[3] = maxi[1]+w*0.2;
	}
	M.useX11(true,useX11text,500,500,
		 x11_zoom[0],x11_zoom[1],x11_zoom[2],x11_zoom[3]);
      }
      M.setX11delay(x11_delay_factor/M.nV());
    }
#endif
    
    if (args_info.smorg_given>0) {
      LOG("Calculating smorg output." << std::endl)
      MeshC MC(&M);
      MC.setOptions(MC.getOptions()|MeshC::Option_offcenter_steiner);

      /* Calculate and collect output. */
      
      matrices.attach("segm.bnd.idx",new Matrix<int>(2),
		      true,fmesh::IOMatrixtype_general);
      matrices.attach("segm.bnd.grp",new Matrix<int>(1),
		      true,fmesh::IOMatrixtype_general);
      MC.segments(true,
		  &matrices.DI("segm.bnd.idx"),
		  &matrices.DI("segm.bnd.grp"));
      
      matrices.output("segm.bnd.idx").output("segm.bnd.grp");
      
      matrices.attach("segm.int.idx",new Matrix<int>(2),
		      true,fmesh::IOMatrixtype_general);
      matrices.attach("segm.int.grp",new Matrix<int>(1),
		      true,fmesh::IOMatrixtype_general);
      MC.segments(false,
		  &matrices.DI("segm.int.idx"),
		  &matrices.DI("segm.int.grp"));
      
      matrices.output("segm.int.idx").output("segm.int.grp");
      
    } else { /* Not smorg.  Build mesh. */

      MeshC MC(&M);
      MC.setOptions(MC.getOptions()|MeshC::Option_offcenter_steiner);

      if (!isflat && !issphere) {
	if (M.nT()==0) {
	  LOG_("Points not in the plane or on a sphere, and triangulation empty."  
	       << std::endl);
	}
	/* Remove everything outside the boundary segments, if any. */
	MC.PruneExterior();
	invalidate_unused_vertex_indices(M, idx);
	/* Nothing more to do here.  Cannot refine non R2/S2 meshes. */
      } else {
	/* If we don't already have a triangulation, we must create one. */
	if (M.nT()==0) {
	  MC.CET(cet_sides,cet_margin);
	}
	
	/* It is more robust to add the constraints before the rest of the
	   nodes are added.  This allows points to fall onto constraint
	   segments, subdividing them as needed. */
	if (cdt_boundary.size()>0)
	  MC.CDTBoundary(cdt_boundary);
	if (cdt_interior.size()>0)
	  MC.CDTInterior(cdt_interior);
	
	/* Add the rest of the nodes. */
	vertexListT vertices;
	for (size_t v=0;v<nV;v++)
	  vertices.push_back(v);
	MC.DT(vertices);
	
	/* Remove everything outside the boundary segments, if any. */
	MC.PruneExterior();
	invalidate_unused_vertex_indices(M, idx);
	
	if (args_info.rcdt_given) {
	  /* Calculate the RCDT: */
	  MC.RCDT(rcdt_min_angle,rcdt_big_limit_auto_default,
		  Quality0.raw(),Quality0.rows(),
		  rcdt_max_n0, rcdt_max_n1);
	  LOG(MC << endl);
	}
	/* Done constructing the triangulation. */
      }
      
      /* Calculate and collect output. */
      
      matrices.attach("segm.bnd.idx",new Matrix<int>(2),
		      true,fmesh::IOMatrixtype_general);
      matrices.attach("segm.bnd.grp",new Matrix<int>(1),
		      true,fmesh::IOMatrixtype_general);
      MC.segments(true,
		  &matrices.DI("segm.bnd.idx"),
		  &matrices.DI("segm.bnd.grp"));
      
      matrices.output("segm.bnd.idx").output("segm.bnd.grp");
      
      matrices.attach("segm.int.idx",new Matrix<int>(2),
		      true,fmesh::IOMatrixtype_general);
      matrices.attach("segm.int.grp",new Matrix<int>(1),
		      true,fmesh::IOMatrixtype_general);
      MC.segments(false,
		  &matrices.DI("segm.int.idx"),
		  &matrices.DI("segm.int.grp"));
      
      matrices.output("segm.int.idx").output("segm.int.grp");
    }
    
    matrices.attach("tt",&M.TT(),false);
    M.useVT(true);
    matrices.attach("vt",&M.VT(),false);
    M.useTTi(true);
    matrices.attach("tti",&M.TTi(),false);
    matrices.attach("vv",new SparseMatrix<int>(M.VV()),
		    true,fmesh::IOMatrixtype_symmetric);
    
    matrices.output("tt").output("tti").output("vt").output("vv");

  }

  LOG("Manifold output." << std::endl)
  /* Output the manifold type. */
  matrices.attach("manifold", new Matrix<int>(1),
		  true, fmesh::IOMatrixtype_general);
  Matrix<int>& manifold = matrices.DI("manifold");
  manifold(0,0) = M.type();
  matrices.output("manifold");

  if (issphere) {
    LOG("issphere output." << std::endl)
    int sph0_order_max = args_info.sph0_arg;
    int sph_order_max = args_info.sph_arg;
    
    if (sph0_order_max >= 0) {
      LOG("sph0 output." << std::endl)
      matrices.attach(string("sph0"),
		      new Matrix<double>(spherical_harmonics(M.S(),
							     sph0_order_max,
							     true)),
		      true);
      matrices.matrixtype("sph0",fmesh::IOMatrixtype_general);
      matrices.output("sph0");
    }

    if (sph_order_max >= 0) {
      LOG("sph output." << std::endl)
      matrices.attach(string("sph"),
		      new Matrix<double>(spherical_harmonics(M.S(),
							     sph_order_max,
							     false)),
		      true);
      matrices.matrixtype("sph",fmesh::IOMatrixtype_general);
      matrices.output("sph");
    }

    if (args_info.bspline_given>0) {
      LOG("bspline output." << std::endl)
      int bspline_n = 2;
      int bspline_degree = 1;
      bool bspline_uniform_knot_angles = true;
      if (args_info.bspline_given>0)
	bspline_n = args_info.bspline_arg[0];
      if (args_info.bspline_given>1)
	bspline_degree = args_info.bspline_arg[1];
      if (args_info.bspline_given>2)
	bspline_uniform_knot_angles = (args_info.bspline_arg[2]>0);

      matrices.attach(string("bspline"),
		      new Matrix<double>(spherical_bsplines(M.S(),
							    bspline_n,
							    bspline_degree,
							    bspline_uniform_knot_angles)),
		      true);
      matrices.matrixtype("bspline",fmesh::IOMatrixtype_general);
      matrices.output("bspline");
    }
  }

  if (args_info.points2mesh_given>0) {
    LOG("points2mesh output." << std::endl);
    if (!isflat && !issphere) {
      LOG_("Cannot calculate points2mesh mapping for non R2/S2 manifolds"
	   << std::endl);
    } else {
      string points2mesh_name(args_info.points2mesh_arg);
      if (!matrices.load(points2mesh_name).active) {
	cout << "Matrix "+points2mesh_name+" not found." << std::endl;
      }
      Matrix<double>& points2mesh = matrices.DD(points2mesh_name);
      size_t points_n = points2mesh.rows();
      Matrix<int>& points2mesh_t =
	matrices.attach(string("p2m.t"),
			new Matrix<int>(points_n,1),
			true);
      Matrix<double>& points2mesh_b =
	matrices.attach(string("p2m.b"),
			new Matrix<double>(points_n,3),
			true);
      matrices.matrixtype("p2m.t",fmesh::IOMatrixtype_general);
      matrices.matrixtype("p2m.b",fmesh::IOMatrixtype_general);
      matrices.output("p2m.t").output("p2m.b");
      
      map_points_to_mesh(M,points2mesh,points2mesh_t,points2mesh_b);
    }
  }
    
  int fem_order_max = args_info.fem_arg;
  if (fem_order_max>=0) {
    LOG("fem output." << std::endl)
    SparseMatrix<double>& C0 = matrices.SD("c0").clear();
    SparseMatrix<double>& C1 = matrices.SD("c1").clear();
    SparseMatrix<double>& B1 = matrices.SD("b1").clear();
    SparseMatrix<double>& G  = matrices.SD("g1").clear();
    SparseMatrix<double>& K  = matrices.SD("k1").clear();
    /* K1=G1-B1, K2=K1*inv(C0)*K1, ... */
    Matrix<double>& Tareas  = matrices.DD("ta").clear();
    
    M.calcQblocks(C0,C1,G,B1,Tareas);

    matrices.attach(string("va"),new Matrix<double>(diag(C0)),true);
    
    K = G-B1;
    
    matrices.matrixtype("c0",fmesh::IOMatrixtype_diagonal);
    matrices.matrixtype("c1",fmesh::IOMatrixtype_symmetric);
    matrices.matrixtype("b1",fmesh::IOMatrixtype_general);
    matrices.matrixtype("g1",fmesh::IOMatrixtype_symmetric);
    matrices.matrixtype("k1",fmesh::IOMatrixtype_general);
    matrices.output("c0");
    matrices.output("c1");
    matrices.output("b1");
    matrices.output("g1");
    matrices.output("k1");
    matrices.output("va");
    matrices.output("ta");
    
    SparseMatrix<double> C0inv = inverse(C0,true);
    SparseMatrix<double> tmp = G*C0inv;
    SparseMatrix<double>* a;
    SparseMatrix<double>* b = &G;
    for (size_t i=1; int(i) < fem_order_max; i++) {
      std::stringstream ss;
      ss << i+1;
      std::string Gname = "g"+ss.str();
      a = b;
      b = &(matrices.SD(Gname).clear());
      *b = tmp*(*a);
      matrices.matrixtype(Gname,fmesh::IOMatrixtype_symmetric);
      matrices.output(Gname);
    }
    tmp = C0inv*K;
    b = &K;
    for (size_t i=1; int(i) < fem_order_max; i++) {
      std::stringstream ss;
      ss << i+1;
      std::string Kname = "k"+ss.str();
      a = b;
      b = &(matrices.SD(Kname).clear());
      *b = (*a)*tmp;
      matrices.matrixtype(Kname,fmesh::IOMatrixtype_general);
      matrices.output(Kname);
    }

    if (aniso_names.size()>0) {
      SparseMatrix<double>& Gani  = matrices.SD("g1aniso").clear();
      M.calcQblocksAni(Gani,
		       matrices.DD(aniso_names[0]),
		       matrices.DD(aniso_names[1]));
      matrices.output("g1aniso");

      SparseMatrix<double> tmp = Gani*C0inv;
      SparseMatrix<double>* a;
      SparseMatrix<double>* b = &Gani;
      for (size_t i=1; int(i) < fem_order_max; i++) {
	std::stringstream ss;
	ss << i+1;
	std::string Gname = "g"+ss.str()+"aniso";
	a = b;
	b = &(matrices.SD(Gname).clear());
	*b = tmp*(*a);
	matrices.matrixtype(Gname,fmesh::IOMatrixtype_symmetric);
	matrices.output(Gname);
      }
    }
    
  }

  if (args_info.grad_given>0) {
    SparseMatrix<double>* D[3];
    M.calcGradientMatrices(D);
    matrices.attach("dx", D[0], true);
    matrices.attach("dy", D[1], true);
    matrices.attach("dz", D[2], true);
    matrices.matrixtype("dx",fmesh::IOMatrixtype_general);
    matrices.matrixtype("dy",fmesh::IOMatrixtype_general);
    matrices.matrixtype("dz",fmesh::IOMatrixtype_general);
    matrices.output("dx").output("dy").output("dz");
  }

  
  if (splitlines_names.size()>0) {
    Matrix<double>& splitlocinput_raw = matrices.DD(splitlines_names[0]);
    /* Make sure we have a Nx3 matrix. */
    Matrix3double splitlocinput(splitlocinput_raw);
    Matrix<double>* splitloc1 = new Matrix<double>(3);
    Matrix<int>* splitidx1 = new Matrix<int>(1);
    Matrix<int>* splittriangle1 = new Matrix<int>(1);
    Matrix<double>* splitbary1 = new Matrix<double>(3);
    Matrix<double>* splitbary2 = new Matrix<double>(3);
    Matrix<int>* splitorigin1 = new Matrix<int>(1);

    split_line_segments_on_triangles(M,
				     splitlocinput,
				     matrices.DI(splitlines_names[1]),
				     *splitloc1,
				     *splitidx1,
				     *splittriangle1,
				     *splitbary1,
				     *splitbary2,
				     *splitorigin1);

    /* Now it's ok to overwrite potential input split* matrices. */
    matrices.attach("split.loc", splitloc1, true);
    matrices.attach("split.idx", splitidx1, true);
    matrices.attach("split.t", splittriangle1, true);
    matrices.attach("split.b1", splitbary1, true);
    matrices.attach("split.b2", splitbary2, true);
    matrices.attach("split.origin", splitorigin1, true);
    matrices.output("split.loc").output("split.idx");
    matrices.output("split.b1").output("split.b2");
    matrices.output("split.t").output("split.origin");
  }

  
  for (size_t i=0; i<args_info.collect_given; i++) {
    string matrix_name = string(args_info.collect_arg[i]);
    if (!(matrix_name=="-") & !(matrix_name=="--")) {
      if (!matrices.activate(matrix_name)) {
	if (!matrices.load(matrix_name).active) {
	  cout << "Matrix "+matrix_name+" not found." << endl;
	} else {
	  cout << "Matrix "+matrix_name+" activated." << endl;
	}
      } else {
	cout << "Matrix "+matrix_name+" active." << endl;
      }
    }
    matrices.output(matrix_name);
  }

  matrices.save();

  cmdline_free(&args_info);

  return 0;
}
