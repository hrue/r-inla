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
using fmesh::Mesh;
using fmesh::MeshC;
using fmesh::Point;
using fmesh::PointRaw;
using fmesh::SparseMatrix;
using fmesh::Vector3;

const bool useVT = true;
const bool useTTi = true;
bool useX11 = true;
const bool useX11text = false;
double x11_delay_factor = 1.0;


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
  ioh.ascii().OH(cout << filename << "\t: ");
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
  ioh.ascii().OH(cout << filename << "\t: ");
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





int main(int argc, char* argv[])
{
  
  gengetopt_args_info args_info;
  struct cmdline_params params;
  
  /* initialize the parameters structure */
  cmdline_params_init(&params);
     
  /* call the command line parser */
  if (cmdline(argc, argv, &args_info) != 0) {
    cmdline_free(&args_info);
    return 1;
  }

  cout << "Config given:\t" << args_info.config_given << endl;
    
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

  int cet_sides = 8;
  double cet_margin = -0.05;
  if ((args_info.cet_given>0) && (args_info.cet_arg[0] != 0))
    cet_sides = args_info.cet_arg[0];
  if (args_info.cet_given>1)
    cet_margin = args_info.cet_arg[1];

  double rcdt_min_angle = 21;
  double rcdt_big_limit = -1.0;
  double rcdt_big_limits = -0.5;
  if ((args_info.rcdt_given>0) && (args_info.rcdt_arg[0] != 0))
    rcdt_min_angle = args_info.rcdt_arg[0];
  if (args_info.rcdt_given>1)
    rcdt_big_limit = args_info.rcdt_arg[1];
  if (args_info.rcdt_given>2)
    rcdt_big_limits = args_info.rcdt_arg[2];

  useX11 = (args_info.x11_given>0);
  x11_delay_factor = args_info.x11_arg;

  cout << "CET given:\t" << args_info.cet_given << endl;
  cout << "RCDT given:\t" << args_info.rcdt_given << endl;
  cout << "X11 given:\t" << args_info.x11_given << endl;

  cout << "CET sides:\t" << cet_sides << endl;
  cout << "CET margin:\t" << cet_margin << endl;
  cout << "RCDT mininmum angle:\t" << rcdt_min_angle << endl;
  cout << "RCDT maximum edge length:\t" << rcdt_big_limit << endl;
  cout << "RCDT maximum edge lengths:\t" << rcdt_big_limits << endl;
  cout << "X11 delay factor:\t" << x11_delay_factor << endl;



  string iprefix("-");
  string oprefix("-");
  if (args_info.inputs_num>0)
    iprefix = string(args_info.inputs[0]);
  if (iprefix=="-")
    oprefix = "fmesher.output.";
  else
    oprefix = iprefix;
  if (args_info.inputs_num>1)
    oprefix = string(args_info.inputs[1]);


  cout << "Input:\t" << iprefix << endl;
  cout << "Output:\t" << oprefix << endl;

  Matrix<double> iS0;
  if (iprefix=="-") {
    IOHelperM<double>().D(&iS0).binary().IH(cin).ID(cin).ascii().OH(cout << "S0: ");
  } else {
    ifstream I((iprefix+"s0").c_str(), ios::in | ios::binary);
    IOHelperM<double>().D(&iS0).binary().IH(I).ID(I).ascii().OH(cout << "S0: ");
    I.close();
  }

  /* Check the input. */
  if (iS0.cols()<2) {
    /* 1D data. Not implemented */
    return 0;
  }
  Matrix3double S0(iS0); /* Make sure we have a Nx3 matrix. */
  int nV = S0.rows();
  if (nV<0) {
    /* No data */
    return 0;
  }
  double radius = S0[0].length();
  bool issphere = true;
  bool isflat = (std::abs(S0[0][2]) < 1.0e-10);
  for (int i=1; i<nV; i++) {
    isflat = (isflat && (std::abs(S0[i][2]) < 1.0e-10));
    issphere = (issphere && (std::abs(S0[i].length()-radius) < 1.0e-10));
  }
  Mesh M(Mesh::Mtype_plane,0,useVT,useTTi);
  if (!isflat) {
    if (issphere) {
      M.type(Mesh::Mtype_sphere);
    } else {
      M.type(Mesh::Mtype_manifold);
    }
  }
  M.S_set(S0);
  M.setX11VBigLimit(nV);

  cout << "Initial mesh:" << endl << M;

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

  double biglim[nV];
  if (rcdt_big_limits==0.0) {
    /* Read limits from a file */
    NOT_IMPLEMENTED;
    return 0;
  } else if (rcdt_big_limits>0.0) {
    for (int v=0;v<nV;v++)
      biglim[v] = rcdt_big_limits;
  } else {
    /* Rudimentary default biglimit construction: */
    for (int v=0;v<nV;v++)
      biglim[v] = diam/std::sqrt(nV)*(-rcdt_big_limits);
  }

  if (useX11) {
    if (issphere) {
      M.useX11(true,useX11text,500,500,
	       -1.1,1.1,
	       -1.1,1.1);
    } else {
      double w0 = maxi[0]-mini[0];
      double w1 = maxi[1]-mini[0];
      M.useX11(true,useX11text,500,500,
	       mini[0]-w0*0.2,maxi[0]+w0*0.2,
	       mini[1]-w1*0.2,maxi[1]+w1*0.2);
    }
    M.setX11delay(x11_delay_factor/M.nV());
  }

  MeshC MC(&M);

  fmesh::vertexListT vertices;
  for (int v=0;v<nV;v++)
    vertices.push_back(v);
  MC.CET(cet_sides,cet_margin);
  MC.DT(vertices);

  if (args_info.rcdt_given) {
    /* Calculate the RCDT: */
    if (rcdt_big_limit<0.0)
      rcdt_big_limit = -rcdt_big_limit*diam;
    MC.RCDT(rcdt_min_angle,rcdt_big_limit,biglim,nV);
  }

  cout << "Final mesh:" << endl << M;

  print_M(oprefix+"s",M.S());
  print_M(oprefix+"tv",M.TV());
  print_M(oprefix+"tt",M.TT());
  M.useTTi(true);
  print_M(oprefix+"tti",M.TTi());
  print_SM(oprefix+"vv",M.VV());

  print_M_old(oprefix+"S.dat",M.S());
  print_M_old(oprefix+"FV.dat",M.TV(),false);

  {
    SparseMatrix<double> C0;
    SparseMatrix<double> C1;
    SparseMatrix<double> B1;
    SparseMatrix<double> G;
    SparseMatrix<double> K; /* K1=G1-B1, K2=K1*inv(C0)*K1, ... */
    M.calcQblocks(C0,C1,G,B1);

    K = G-B1;

    print_SM(oprefix+"c0",C0,fmesh::IOMatrixtype_diagonal);
    print_SM(oprefix+"c1",C1,fmesh::IOMatrixtype_symmetric);
    print_SM(oprefix+"b1",B1,fmesh::IOMatrixtype_symmetric);
    print_SM(oprefix+"g1",G,fmesh::IOMatrixtype_symmetric);
    print_SM(oprefix+"k1",K,fmesh::IOMatrixtype_symmetric);

    print_SM_old(oprefix+"C.dat",C0,true,fmesh::IOMatrixtype_diagonal);
    print_SM_old(oprefix+"G.dat",G,false,fmesh::IOMatrixtype_symmetric);
    print_SM_old(oprefix+"K.dat",K,false,fmesh::IOMatrixtype_symmetric);

    SparseMatrix<double> C0inv = inverse(C0,true);
    SparseMatrix<double> tmp = G*C0inv;
    for (int i=0; i<3; i++) {
      G = tmp*G;
      std::stringstream ss;
      ss << i+2;
      print_SM(oprefix+"g"+ss.str(),G,fmesh::IOMatrixtype_symmetric);
      print_SM_old(oprefix+"G"+ss.str()+".dat",G,
		   false,fmesh::IOMatrixtype_symmetric);
    }
    tmp = C0inv*K;
    for (int i=0; i<3; i++) {
      K = K*tmp;
      std::stringstream ss;
      ss << i+2;
      print_SM(oprefix+"k"+ss.str(),K,fmesh::IOMatrixtype_symmetric);
      print_SM_old(oprefix+"K"+ss.str()+".dat",K,
		   false,fmesh::IOMatrixtype_symmetric);
    }


  }

  cmdline_free(&args_info);

  return 0;
}
