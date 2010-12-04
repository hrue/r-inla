#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "fmesher.hh"

using std::ios;
using std::ifstream;
using std::ofstream;
using std::string;
using std::cin;
using std::cout;
using std::endl;

using fmesh::Dart;
using fmesh::Int3;
using fmesh::Int3Raw;
using fmesh::IOHelper;
using fmesh::IOHelperM;
using fmesh::IOHelperSM;
using fmesh::Matrix;
using fmesh::Matrix3;
using fmesh::Matrix3int;
using fmesh::Matrix3double;
using fmesh::MatrixC;
using fmesh::Mesh;
using fmesh::MeshC;
using fmesh::Point;
using fmesh::PointRaw;
using fmesh::SparseMatrix;
using fmesh::Vector3;

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
  ofstream O;
  O.open(filename.c_str(), ios::out | ios::binary);
  IOHelperM<T> ioh;
  ioh.cD(&M).matrixtype(matrixt);
  ioh.binary().OH(O).OD(O);
  O.close();
}

template <class T>
void print_SM(string filename,
	      const SparseMatrix<T>& M,
	      fmesh::IOMatrixtype matrixt = fmesh::IOMatrixtype_general)
{
  ofstream O;
  O.open(filename.c_str(), ios::out | ios::binary);
  IOHelperSM<T> ioh;
  ioh.cD(&M).matrixtype(matrixt);
  ioh.binary().OH(O).OD(O);
  O.close();
}

template <class T>
void print_M_old(string filename,
		 const Matrix<T>& M,
		 bool header = true,
		 fmesh::IOMatrixtype matrixt = fmesh::IOMatrixtype_general)
{
  ofstream O;
  O.open(filename.c_str(), ios::out);
  IOHelperM<T> ioh;
  ioh.cD(&M).ascii().matrixtype(matrixt);
  if (header) ioh.OH_2009(O);
  ioh.OD_2009(O);
  O.close();
}

template <class T>
void print_SM_old(string filename,
		  const SparseMatrix<T>& M,
		  bool header = false,
		  fmesh::IOMatrixtype matrixt = fmesh::IOMatrixtype_general)
{
  ofstream O;
  O.open(filename.c_str(), ios::out);
  IOHelperSM<T> ioh;
  ioh.cD(&M).ascii().matrixtype(matrixt);
  if (header) ioh.OH_2009(O);
  ioh.OD_2009(O);
  O.close();
}




void map_points_to_mesh(const Mesh& M,
			const Matrix<double>& points,
			Matrix<int>& point2T,
			Matrix<double>& point2bary)
{
  Dart d0(M);
  Dart d;
  Point s;
  Point b;
  for (int i=0; i<points.rows(); i++) {
    s[0] = points[i][0];
    s[1] = points[i][1];
    s[2] = points[i][2];
    //    std::cout << i << endl;
    //    std::cout << s << endl;
    d = M.locate_point(Dart(M),s);
    if (!d.isnull()) { /* Point located. */
      //      std::cout << d << endl;
      M.barycentric(d,s,b);
      //      std::cout << b << endl;
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



int main(int argc, char* argv[])
{
  
  gengetopt_args_info args_info;
  struct cmdline_params params;
  
  cmdline_init(&args_info);
  cmdline_params_init(&params);
     
  /* call the command line parser */
  if (cmdline_ext(argc, argv, &args_info, &params) != 0) {
    cmdline_free(&args_info);
    return 1;
  }

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

  if (args_info.dump_config_given)
    cmdline_dump(stdout,&args_info);
  
  std::vector<string> input_s0_names;
  string input_tv0_name = "-";
  
  input_s0_names.push_back(string(args_info.input_arg[0]));
  if (args_info.input_given>1)
    input_tv0_name = string(args_info.input_arg[1]);
  for (int i=1; i<int(args_info.input_given)-1; i++) {
    input_s0_names.push_back(string(args_info.input_arg[i+1]));
  }

  std::vector<string> quality_names;
  for (int i=0; i<int(args_info.quality_given); i++) {
    quality_names.push_back(string(args_info.quality_arg[i]));
  }


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
  for (int i=1; (i<int(args_info.input_given)-2); i++) {
    rcdt_big_limit_defaults(i,0) = -0.5;
  }
  if ((args_info.rcdt_given>0) && (args_info.rcdt_arg[0] != 0))
    rcdt_min_angle = args_info.rcdt_arg[0];
  if (args_info.rcdt_given>1)
    rcdt_big_limit_auto_default = args_info.rcdt_arg[1];
  for (int i=0; (i < int(args_info.rcdt_given)-2); i++) {
    rcdt_big_limit_defaults(i,0) = args_info.rcdt_arg[i+2];
  }

  useX11 = (args_info.x11_given>0) && (args_info.x11_arg>=0);
  x11_delay_factor = args_info.x11_arg;
  if (args_info.x11_zoom_given==4) {
    for (int i=0; i<4; i++)
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
  cout << "X11 delay factor:\t" << x11_delay_factor << endl;
  */

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

  matrices.io(((args_info.io_arg == io_arg_ba) ||
	       (args_info.io_arg == io_arg_bb)),
	      ((args_info.io_arg == io_arg_ab) ||
	       (args_info.io_arg == io_arg_bb)));
  matrices.input_prefix(iprefix);
  matrices.output_prefix(oprefix);
  for (int i=0; i<(int)args_info.ic_given; ++i) {
    matrices.input_file(string(args_info.ic_arg[i]));
  }
  if (args_info.oc_given>0) {
    matrices.output_file(string(args_info.oc_arg));
  }
  for (int i=0; i+2<(int)args_info.ir_given; i=i+3) {
    matrices.input_raw(string(args_info.ir_arg[i]),
		       string(args_info.ir_arg[i+1]),
		       string(args_info.ir_arg[i+2]));
  }

  for (int i=0; i<input_s0_names.size(); i++) {
    if (!matrices.load(input_s0_names[i]).active) {
      cout << "Matrix "+input_s0_names[i]+" not found." << endl;
    }
  }
  for (int i=0; i<quality_names.size(); i++) {
    if (quality_names[i] != "-")
      if (!matrices.load(quality_names[i]).active) {
	cout << "Matrix "+quality_names[i]+" not found." << endl;
	quality_names[i] = "-";
      }
  }

  Matrix<double>& iS0 = matrices.DD(input_s0_names[0]);
  Matrix<double>* Quality0_ = new Matrix<double>();
  Matrix<double>& Quality0 = *Quality0_;

  /* Join the location matrices */ 
  for (int i=0; i < input_s0_names.size(); i++) {
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
    } else if (i<rcdt_big_limit_defaults.rows()) {
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


  Matrix<int>* TV0 = NULL;
  if (input_tv0_name != "-") {
    if (!matrices.load(input_tv0_name).active) {
      cout << "Matrix "+input_tv0_name+" not found." << endl;
    } else {
      TV0 = &(matrices.DI(input_tv0_name));
    }
  }


  fmesh::constrListT cdt_boundary;
  if (args_info.boundary_given) {
    string b_name = string(args_info.boundary_arg[0]);
    if (!matrices.load(b_name).active) {
      cout << "Matrix "+b_name+" not found." << endl;
    }
    Matrix<int>& boundary0 = matrices.DI(b_name);
    if (boundary0.cols()==1) {
      int v0 = -1;
      int v1 = -1;
      for (int i=0; i < boundary0.rows(); i++) {
	v0 = v1;
	v1 = boundary0[i][0];
	if ((v0>=0) && (v1>=0))
	  cdt_boundary.push_back(fmesh::constrT(v0,v1));
      }
    }
  }

  fmesh::constrListT cdt_interior;
  if (args_info.interior_given) {
    string b_name = string(args_info.interior_arg[0]);
    if (!matrices.load(b_name).active) {
      cout << "Matrix "+b_name+" not found." << endl;
    }
    Matrix<int>& interior0 = matrices.DI(b_name);
    if (interior0.cols()==1) {
      int v0 = -1;
      int v1 = -1;
      for (int i=0; i < interior0.rows(); i++) {
	v0 = v1;
	v1 = interior0[i][0];
	if ((v0>=0) && (v1>=0))
	  cdt_interior.push_back(fmesh::constrT(v0,v1));
      }
    }
  }

  Mesh M(Mesh::Mtype_plane,0,useVT,useTTi);

  int nV = iS0.rows();
  bool issphere = false;
  if ((nV>0) && (iS0.cols()<2)) {
    /* 1D data. Not implemented */
    return 0;
  } else if (nV>0) {
    Matrix3double S0(iS0); /* Make sure we have a Nx3 matrix. */

    double radius = S0[0].length();
    issphere = true;
    bool isflat = (std::abs(S0[0][2]) < 1.0e-10);
    for (int i=1; i<nV; i++) {
      isflat = (isflat && (std::abs(S0[i][2]) < 1.0e-10));
      issphere = (issphere && (std::abs(S0[i].length()-radius) < 1.0e-10));
    }
    if (!isflat) {
      if (issphere) {
	M.type(Mesh::Mtype_sphere);
      } else {
	M.type(Mesh::Mtype_manifold);
      }
    }

    M.S_set(S0);
    M.setX11VBigLimit(nV);
    
    if (TV0) {
      M.TV_set(*TV0);
    }

    matrices.attach(string("s"),&M.S(),false);
    matrices.attach("tv",&M.TV(),false);
    matrices.output("s").output("tv");

    if (isflat || issphere) {
    
    Point mini(S0(0));
    Point maxi(S0(0));
    for (int v=1; v<nV; v++)
      for (int i=0; i<3; i++) {
	mini[i] = (S0(v)[i] < mini[i] ? S0(v)[i] : mini[i]);
	maxi[i] = (S0(v)[i] > maxi[i] ? S0(v)[i] : maxi[i]);
      }
    Point sz;
    fmesh::Vec::diff(sz,maxi,mini);
    double diam;
    diam = (sz[1] < sz[0]
	    ? (sz[2] < sz[0] ? sz[0] : sz[2])
	    : (sz[2] < sz[1] ? sz[1] : sz[2]));
    
    for (int i=0; i<Quality0.rows(); i++) {
      if (Quality0[i][0]<0.0) {
	/* Rudimentary relative biglimit construction: */
	Quality0(i,0) = diam/std::sqrt(nV)*(-Quality0[i][0]);
      }
    }
    
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
	double w1 = maxi[1]-mini[0];
	if (args_info.x11_zoom_given==0) {
	  x11_zoom[0] = mini[0]-w0*0.2;
	  x11_zoom[1] = maxi[0]+w0*0.2;
	  x11_zoom[2] = mini[1]-w1*0.2;
	  x11_zoom[3] = maxi[1]+w1*0.2;
	}
	M.useX11(true,useX11text,500,500,
		 x11_zoom[0],x11_zoom[1],x11_zoom[2],x11_zoom[3]);
      }
      M.setX11delay(x11_delay_factor/M.nV());
    }
    
    MeshC MC(&M);
    MC.setOptions(MC.getOptions()|MeshC::Option_offcenter_steiner);

    if (!TV0)
      MC.CET(cet_sides,cet_margin);
    
    /* TODO: Check that this is ok even when some or all points are
       already in the triangulation. */
    fmesh::vertexListT vertices;
    for (int v=0;v<nV;v++)
      vertices.push_back(v);
    MC.DT(vertices);

    if (cdt_boundary.size()>0) {
      MC.CDTBoundary(cdt_boundary);
    }
    if (cdt_interior.size()>0)
      MC.CDTInterior(cdt_interior);
    MC.PruneExterior();
    
    if (args_info.rcdt_given) {
      /* Calculate the RCDT: */
      MC.RCDT(rcdt_min_angle,rcdt_big_limit_auto_default,
	      Quality0.raw(),Quality0.rows());
    }

    }
    
    if (false) {
    if (oprefix != "-") {
      print_M_old(oprefix+"S.dat",M.S());
      print_M_old(oprefix+"FV.dat",M.TV(),false);
    }
    }

    matrices.attach("tt",&M.TT(),false);
    M.useVT(true);
    matrices.attach("vt",&M.VT(),false);
    M.useTTi(true);
    matrices.attach("tti",&M.TT(),false);
    matrices.attach("vv",new SparseMatrix<int>(M.VV()),
		    true,fmesh::IOMatrixtype_symmetric);
    
    matrices.output("tt").output("tti").output("vt").output("vv");
    
  }

  if (issphere) {
    int sph0_order_max = args_info.sph0_arg;
    int sph_order_max = args_info.sph_arg;
    
    if (sph0_order_max >= 0) {
      matrices.attach(string("sph0"),
		      new Matrix<double>(spherical_harmonics(M.S(),
							     sph0_order_max,
							     true)),
		      true);
      matrices.matrixtype("sph0",fmesh::IOMatrixtype_general);
      matrices.output("sph0");
    }

    if (sph_order_max >= 0) {
      matrices.attach(string("sph"),
		      new Matrix<double>(spherical_harmonics(M.S(),
							     sph_order_max,
							     false)),
		      true);
      matrices.matrixtype("sph",fmesh::IOMatrixtype_general);
      matrices.output("sph");
    }

    if (args_info.bspline_given>0) {
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
    string points2mesh_name(args_info.points2mesh_arg);
    if (!matrices.load(points2mesh_name).active) {
      cout << "Matrix "+points2mesh_name+" not found." << endl;
    }
    Matrix<double>& points2mesh = matrices.DD(points2mesh_name);
    int points_n = points2mesh.rows();
    Matrix<int>& points2mesh_t =
      matrices.attach(string("points2mesh.t"),
		      new Matrix<int>(points_n,1),
		      true);
    Matrix<double>& points2mesh_b =
      matrices.attach(string("points2mesh.b"),
		      new Matrix<double>(points_n,3),
		      true);
    matrices.matrixtype("points2mesh.t",fmesh::IOMatrixtype_general);
    matrices.matrixtype("points2mesh.b",fmesh::IOMatrixtype_general);
    matrices.output("points2mesh.t").output("points2mesh.b");
    
    map_points_to_mesh(M,points2mesh,points2mesh_t,points2mesh_b);
  }
    
  int fem_order_max = args_info.fem_arg;
  if (fem_order_max>0) {
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
    matrices.matrixtype("k1",fmesh::IOMatrixtype_symmetric);
    matrices.output("c0");
    matrices.output("c1");
    matrices.output("b1");
    matrices.output("g1");
    matrices.output("k1");
    matrices.output("va");
    matrices.output("ta");
    
    if (false) {
    if (oprefix != "-") {
      print_SM_old(oprefix+"C.dat",C0,true,fmesh::IOMatrixtype_diagonal);
      print_SM_old(oprefix+"G.dat",G,false,fmesh::IOMatrixtype_symmetric);
      print_SM_old(oprefix+"K.dat",K,false,fmesh::IOMatrixtype_symmetric);
    }
    }
    
    SparseMatrix<double> C0inv = inverse(C0,true);
    SparseMatrix<double> tmp = G*C0inv;
    SparseMatrix<double>* a;
    SparseMatrix<double>* b = &G;
    for (int i=1; i<fem_order_max; i++) {
      std::stringstream ss;
      ss << i+1;
      std::string Gname = "g"+ss.str();
      a = b;
      b = &(matrices.SD(Gname).clear());
      *b = tmp*(*a);
      matrices.matrixtype(Gname,fmesh::IOMatrixtype_symmetric);
      matrices.output(Gname);
      
      if (false) {
      if (oprefix != "-") {
	Gname = "G"+ss.str();
	print_SM_old(oprefix+Gname+".dat",*b,
		     false,fmesh::IOMatrixtype_symmetric);
      }
      }
    }
    tmp = C0inv*K;
    b = &K;
    for (int i=1; i<fem_order_max; i++) {
      std::stringstream ss;
      ss << i+1;
      std::string Kname = "k"+ss.str();
      a = b;
      b = &(matrices.SD(Kname).clear());
      *b = (*a)*tmp;
      matrices.matrixtype(Kname,fmesh::IOMatrixtype_symmetric);
      matrices.output(Kname);
      
      if (false) {
      if (oprefix != "-") {
	Kname = "K"+ss.str();
	print_SM_old(oprefix+Kname+".dat",*b,
		     false,fmesh::IOMatrixtype_symmetric);
      }
      }
    }
    
  }
  

  for (int i=0; i<(int)args_info.collect_given; i++) {
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
