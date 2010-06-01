#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <fstream>

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
    double w0 = maxi[0]-mini[0];
    double w1 = maxi[1]-mini[0];
    M.useX11(true,useX11text,500,500,
	     mini[0]-w0*0.2,maxi[0]+w0*0.2,
	     mini[1]-w1*0.2,maxi[1]+w1*0.2);
  }

  MeshC MC(&M);

  fmesh::vertexListT vertices;
  for (int v=0;v<nV;v++)
    vertices.push_back(v);
  MC.DT(vertices);

  /* Calculate the RCDT: */
  MC.RCDT(21,diam*2.0,biglim,nV);

  cout << "Final mesh:" << endl << M;

  ofstream O;
  O.open((prefix+"s").c_str(), ios::out | ios::binary);
  IOHelperM<double>().cD(&M.S()).binary().OH(O).OD(O).ascii().OH(cout << "S : ");
  O.close();

  O.open((prefix+"tv").c_str(), ios::out | ios::binary);
  IOHelperM<int>().cD(&M.TV()).binary().OH(O).OD(O).ascii().OH(cout << "TV: ");
  O.close();

  {
    SparseMatrix<double> C0;
    SparseMatrix<double> C1;
    SparseMatrix<double> G1;
    SparseMatrix<double> B1;
    M.calcQblocks(C0,C1,G1,B1);
    
    O.open((prefix+"c0").c_str(), ios::out | ios::binary);
    IOHelperSM<double>().cD(&C0).diagonal().binary().OH(O).OD(O).ascii().OH(cout << "C0: ");
    O.close();

    O.open((prefix+"c1").c_str(), ios::out | ios::binary);
    IOHelperSM<double>().cD(&C1).symmetric().binary().OH(O).OD(O).ascii().OH(cout << "C1: ");
    O.close();
    
    O.open((prefix+"g1").c_str(), ios::out | ios::binary);
    IOHelperSM<double>().cD(&G1).symmetric().binary().OH(O).OD(O).ascii().OH(cout << "G1: ");
    O.close();

    O.open((prefix+"b1").c_str(), ios::out | ios::binary);
    IOHelperSM<double>().cD(&B1).symmetric().binary().OH(O).OD(O).ascii().OH(cout << "B1: ");
    O.close();
  }

  return 0;
}
