#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "fmesher.h"

using std::ios;
using std::ifstream;
using std::ofstream;
using std::string;
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
const bool useX11 = true;
const bool useX11text = false;
const double x11_delay_factor = 1.0;


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





int main(int argc, const char* argv[])
{
  string prefix("fmesher.io.");
  if (argc>=2)
    prefix = string(argv[1]);

  Matrix<double> iS0;
  ifstream I((prefix+"s0").c_str(), ios::in | ios::binary);
  IOHelperM<double>().D(&iS0).binary().IH(I).ID(I).ascii().OH(cout << "S0: ");
  I.close();

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

  /* Rudimentary default biglimit construction: */
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
  for (int v=0;v<nV;v++)
    biglim[v] = diam/std::sqrt(nV)*0.5;

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
  MC.CET(8,-0.1);
  MC.DT(vertices);

  /* Calculate the RCDT: */
  MC.RCDT(21,diam*2.0,biglim,nV);

  cout << "Final mesh:" << endl << M;

  print_M(prefix+"s",M.S());
  print_M(prefix+"tv",M.TV());
  print_M(prefix+"tt",M.TT());
  M.useTTi(true);
  print_M(prefix+"tti",M.TTi());
  print_SM(prefix+"vv",M.VV());

  print_M_old(prefix+"S.dat",M.S());
  print_M_old(prefix+"FV.dat",M.TV(),false);

  {
    SparseMatrix<double> C0;
    SparseMatrix<double> C1;
    SparseMatrix<double> B1;
    SparseMatrix<double> G;
    SparseMatrix<double> K; /* K1=G1-B1, K2=K1*inv(C0)*K1, ... */
    M.calcQblocks(C0,C1,G,B1);

    K = G-B1;

    print_SM(prefix+"c0",C0,fmesh::IOMatrixtype_diagonal);
    print_SM(prefix+"c1",C1,fmesh::IOMatrixtype_symmetric);
    print_SM(prefix+"b1",B1,fmesh::IOMatrixtype_symmetric);
    print_SM(prefix+"g1",G,fmesh::IOMatrixtype_symmetric);
    print_SM(prefix+"k1",K,fmesh::IOMatrixtype_symmetric);

    print_SM_old(prefix+"C.dat",C0,true,fmesh::IOMatrixtype_diagonal);
    print_SM_old(prefix+"G.dat",G,false,fmesh::IOMatrixtype_symmetric);
    print_SM_old(prefix+"K.dat",K,false,fmesh::IOMatrixtype_symmetric);

    SparseMatrix<double> C0inv = inverse(C0,true);
    SparseMatrix<double> tmp = G*C0inv;
    for (int i=0; i<3; i++) {
      G = tmp*G;
      std::stringstream ss;
      ss << i+2;
      print_SM(prefix+"g"+ss.str(),G,fmesh::IOMatrixtype_symmetric);
      print_SM_old(prefix+"G"+ss.str()+".dat",G,
		   false,fmesh::IOMatrixtype_symmetric);
    }
    tmp = C0inv*K;
    for (int i=0; i<3; i++) {
      K = K*tmp;
      std::stringstream ss;
      ss << i+2;
      print_SM(prefix+"k"+ss.str(),K,fmesh::IOMatrixtype_symmetric);
      print_SM_old(prefix+"K"+ss.str()+".dat",K,
		   false,fmesh::IOMatrixtype_symmetric);
    }


  }

  return 0;
}
