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
  IOHelperM<double>().D(&iS0).binary().IH(I).ID(I).ascii().OH(cout);
  I.close();
  Matrix3double S0(iS0); /* Make sure we have a Nx3 matrix. */
  int nV = S0.rows();

  Mesh M(Mesh::Mtype_plane,0,useVT,useTTi);
  M.S_set(S0);
  M.setX11VBigLimit(nV);
  if (useX11)
    M.useX11(true,useX11text,500,500);

  MeshC MC(&M);

  fmesh::vertexListT vertices;
  for (int v=0;v<nV;v++)
    vertices.push_back(v);
  MC.DT(vertices);

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

  /* Calculate the RCDT: */
  MC.RCDT(21,diam*2.0,biglim,nV);

  ofstream O;
  O.open((prefix+"s").c_str(), ios::out | ios::binary);
  IOHelperM<double>().cD(&M.S()).binary().OH(O).OD(O).ascii().OH(cout);
  O.close();

  O.open((prefix+"tv").c_str(), ios::out | ios::binary);
  IOHelperM<int>().cD(&M.TV()).binary().OH(O).OD(O).ascii().OH(cout);
  O.close();

  {
    SparseMatrix<double> C;
    SparseMatrix<double> C0;
    SparseMatrix<double> G1;
    M.calcQblocks(C,C0,G1);
    
    O.open((prefix+"c").c_str(), ios::out | ios::binary);
    IOHelperSM<double>().cD(&C).symmetric().binary().OH(O).OD(O).ascii().OH(cout);
    O.close();

    O.open((prefix+"c0").c_str(), ios::out | ios::binary);
    IOHelperSM<double>().cD(&C0).diagonal().binary().OH(O).OD(O).ascii().OH(cout);
    O.close();
    
    O.open((prefix+"g1").c_str(), ios::out | ios::binary);
    IOHelperSM<double>().cD(&G1).symmetric().binary().OH(O).OD(O).ascii().OH(cout);
    O.close();
  }

  return 0;
}
