#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "fmesher.h"
#include "predicates.h"

using std::cout;
using std::endl;
using std::ios;
using std::vector;

using fmesh::Dart;
using fmesh::Int3;
using fmesh::Int3Raw;
using fmesh::IOHelper;
using fmesh::IOHelperM;
using fmesh::IOHelperSM;
using fmesh::Matrix;
using fmesh::Matrix1;
using fmesh::Matrix1int;
using fmesh::Matrix3;
using fmesh::Matrix3double;
using fmesh::Matrix3int;
using fmesh::Mesh;
using fmesh::MeshC;
using fmesh::Point;
using fmesh::PointRaw;
using fmesh::Vector3;

const bool useVT = true;
const bool useTTi = true;
const bool useX11 = true;
const bool useX11text = false;
const int maxiter = 1;

int mesh_test() {
  PointRaw S[4] = {{0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.}, {1., 1., 0.}};
  Int3Raw TV[2] = {{0, 1, 2}, {3, 2, 1}};

  Mesh M(Mesh::Mtype_plane, 0, useVT, useTTi);

  M.S_set(Matrix3double(4, S));
  M.TV_set(Matrix3int(2, TV));

  Mesh M2 = M;

  Dart d = Dart(M, 0, 1, 0);
  cout << d.onBoundary();
  cout << d.orbit2().onBoundary();
  cout << M;
  d = M.swapEdge(Dart(M, 0, 1, 1));
  d = M.swapEdge(Dart(M, 0, 1, 1));
  cout << M;
  PointRaw s = {0.2, 0.2, 0.0};
  M.S_append(Matrix3double(1, &s));
  cout << M;
  d = M.splitTriangle(Dart(M, 1, 1, 0), M.nV() - 1);
  cout << M;
  cout << d << endl;

  /*
cout << "d  :" << d << "\n";
cout << "a0 :" << d.alpha0() << "\n";
cout << "a0 :" << d.alpha0() << "\n";
cout << "a1 :" << d.alpha1() << "\n";
cout << "a1 :" << d.alpha1() << "\n";
cout << "a2 :" << d.alpha2() << "\n";
cout << "o2 :" << d.orbit2() << "\n";
cout << "o2 :" << d.orbit2() << "\n";
cout << "a1 :" << d.alpha1() << "\n";
cout << "a2 :" << d.alpha2() << "\n";
cout << "a0 :" << d.alpha0() << "\n";
cout << "a2 :" << d.alpha2() << "\n";
cout << "o1 :" << d.orbit1() << "\n";
Dart d2(d);
cout << "?? :" << (d2 == d.orbit1()) << "\n";
*/

  return 0;
}

int CDT_test() {
  int n = 5;
  PointRaw S[5] = {{0.3, 0.5, 0},
                   {0.6, 0.6, 0},
                   {0.7, 0.5, 0},
                   {0.7, 0.4, 0},
                   {0.5, 0.4, 0}};

  Mesh M(Mesh::Mtype_plane, 0, useVT, useTTi);

  M.S_set(Matrix3double(n, S));

#ifdef FMESHER_WITH_X
  M.setX11VBigLimit(n);
  if (useX11)
    M.useX11(true, useX11text, 500, 500);
#endif

  MeshC MC(&M);
  MC.CETplane(8, 0.1);

  fmesh::vertexListT vertices;
  MC.DT(vertices);
  for (int v = 0; v < n; v++)
    vertices.push_back(v);
  MC.DT(vertices);

  fmesh::constrListT cinp;
  cinp.push_back(fmesh::constrT(0, 1));
  cinp.push_back(fmesh::constrT(0, 2));
  cinp.push_back(fmesh::constrT(0, 3));
  cinp.push_back(fmesh::constrT(0, 4));
  cinp.push_back(fmesh::constrT(1, 2));
  cinp.push_back(fmesh::constrT(2, 3));
  cinp.push_back(fmesh::constrT(3, 4));
  MC.CDTInterior(cinp);

  cinp.clear();
  cinp.push_back(fmesh::constrT(0, 3, 1));
  cinp.push_back(fmesh::constrT(3, 4, 2));
  cinp.push_back(fmesh::constrT(4, 0, 3));
  MC.CDTBoundary(cinp);

  MC.PruneExterior();

  MC.RCDT(25, 100);
  MC.RCDT(25, 0.05);

  cout << MC;

  return 0;
}

int DT2D_test() {
  int n = 14;
  double S[14][3] = {{0.5, 0.5, 0},   {0.5, 0.6, 0},   {0.3, 0.2, 0},
                     {0.3, 0.6, 0},   {0.5, 0.3, 0},   {0.6, 0.7, 0},
                     {0.7, 0.3, 0},   {0.2, 0.8, 0},   {0.9, 0.5, 0},
                     {0.1, 0.1, 0},   {0.05, 0.05, 0}, {0.95, 0.05, 0},
                     {0.95, 0.95, 0}, {0.05, 0.95, 0}};
  double Sb[4][3] = {{0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.}, {1., 1., 0.}};
  int TVb[2][3] = {{0, 1, 2}, {3, 2, 1}};
  Mesh M(Mesh::Mtype_plane, 0, useVT, useTTi);
  int t, vi, v;

  M.S_set(Matrix3double(n, S));

#ifdef FMESHER_WITH_X
  M.setX11VBigLimit(n);
  if (useX11)
    M.useX11(true, useX11text, 500, 500);
#endif

  M.S_append(Matrix3double(4, Sb));
  for (t = 0; t < 2; t++)
    for (vi = 0; vi < 3; vi++)
      TVb[t][vi] += n;
  M.TV_set(Matrix3int(2, TVb));

  MeshC MC(&M);
  fmesh::vertexListT vertices;
  for (v = 0; v < n; v++)
    vertices.push_back(v);

  MC.DT(vertices);

  fmesh::constrListT cinp;
  cinp.push_back(fmesh::constrT(10, 11));
  cinp.push_back(fmesh::constrT(11, 12));
  cinp.push_back(fmesh::constrT(12, 13));
  cinp.push_back(fmesh::constrT(13, 10));
  MC.CDTBoundary(cinp);

  // cinp.clear();
  // cinp.push_back(fmesh::constrT(10,12));
  // MC.CDTInterior(cinp);

  MC.PruneExterior();

  MC.RCDT(25, 100);
  MC.RCDT(25, 0.05);

  return 0;
}

int DT2D_test2() {
  int n = 25;
  double S[25][3] = {{0.1, 0.1, 0}, {0.3, 0.1, 0}, {0.7, 0.1, 0}, {0.9, 0.1, 0},
                     {0.1, 0.3, 0}, {0.3, 0.3, 0}, {0.7, 0.3, 0}, {0.9, 0.3, 0},
                     {0.1, 0.5, 0}, {0.3, 0.5, 0}, {0.7, 0.5, 0}, {0.9, 0.5, 0},
                     {0.1, 0.7, 0}, {0.3, 0.7, 0}, {0.7, 0.7, 0}, {0.9, 0.7, 0},
                     {0.1, 0.9, 0}, {0.3, 0.9, 0}, {0.7, 0.9, 0}, {0.9, 0.9, 0},
                     {0.5, 0.1, 0}, {0.5, 0.3, 0}, {0.5, 0.5, 0}, {0.5, 0.7, 0},
                     {0.5, 0.9, 0}};
  double Sb[4][3] = {{0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.}, {1., 1., 0.}};
  int TVb[2][3] = {{0, 1, 2}, {3, 2, 1}};
  Mesh M(Mesh::Mtype_plane, 0, useVT, useTTi);
  int t, vi, v;

  /*
  for (int v=0;v<n;v++) {
    S[v][0] = S[v][0]+double(v)/1000.0;
    S[v][1] = S[v][1]+double(v)/1001.0;
  }
  */

  M.S_set(Matrix3double(n, S));

#ifdef FMESHER_WITH_X
  M.setX11VBigLimit(n);
  if (useX11)
    M.useX11(true, useX11text, 500, 500);
#endif

  M.S_append(Matrix3double(4, Sb));
  for (t = 0; t < 2; t++)
    for (vi = 0; vi < 3; vi++)
      TVb[t][vi] += n;
  M.TV_set(Matrix3int(2, TVb));

  MeshC MC(&M);
  fmesh::vertexListT vertices;
  for (v = 0; v < n; v++)
    vertices.push_back(v);
  MC.DT(vertices);

  MC.CDT(fmesh::constrListT(), fmesh::constrListT());

  MC.PruneExterior();

  MC.RCDT(25, 100);
  MC.RCDT(25, 0.05);

  return 0;
}

int DT2D_test3() /* Random points */
{
  int n = 200;
  PointRaw S[200];
  double Sb[4][3] = {{0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.}, {1., 1., 0.}};
  int TVb[2][3] = {{0, 1, 2}, {3, 2, 1}};
  Mesh M(Mesh::Mtype_plane, 0, useVT, useTTi);
  int t, vi, v;

  for (v = 0; v < n; v++) {
    S[v][0] = double(std::rand()) / RAND_MAX * 0.9 + 0.05;
    S[v][1] = double(std::rand()) / RAND_MAX * 0.9 + 0.05;
    S[v][2] = 0.0;
  }

  std::ofstream F("ex.random.io.s0", ios::out | ios::binary);
  Matrix3double FM(n, S);
  fmesh::IOHelperM<double>().cD(&FM).binary().OH(F).OD(F);
  F.close();

  M.S_set(Matrix3double(n, S));

#ifdef FMESHER_WITH_X
  M.setX11VBigLimit(n);
  if (useX11)
    M.useX11(true, useX11text, 500, 500);
    //      M.useX11(true,useX11text,500,500,0.7,0.9,0.4,0.85);
#endif

  M.S_append(Matrix3double(4, Sb));
  for (t = 0; t < 2; t++)
    for (vi = 0; vi < 3; vi++)
      TVb[t][vi] += n;
  M.TV_set(Matrix3int(2, TVb));

  MeshC MC(&M);
  fmesh::vertexListT vertices;
  for (v = 0; v < n; v++)
    vertices.push_back(v);

  MC.DT(vertices);

  /*
  fmesh::triangleSetT triangles;
  for (t=0;t<(int)M.nT();t++)
    triangles.insert(t);
  MC.LOP(triangles);
  */

  //  MC.CDTSegment(false,0,1);
  MC.CDTSegment(true, 0, 1); // Make an interior cut in the manifold.
  MC.CDTSegment(true, 1, 0);

  MC.PruneExterior();

  double *biglim = new double[n];
  for (v = 0; v < n; v++)
    biglim[v] = 0.5;

  MC.RCDT(30, 1.0, biglim, n);
  cout << MC;

  delete[] biglim;

  cout << MC;

  return 0;
}

int DTsphere_test() {
  int n = 12;
  double S[12][3] = {{0.2, 0.2, 0.8},  {-0.1, 0.1, 0.1}, {-0.2, 0.5, -0.10},
                     {0.3, -0.2, 0.2}, {0.4, 0.6, 0.2},  {-0.5, 0.3, -0.3},
                     {-0.6, 0.7, 0.3}, {0.2, 0.3, 2.0},  {-0.8, -0.8, -0.4},
                     {-0.9, 0.5, 0.5}, {-0.1, 0.1, 1.0}, {0.0, 7.0, 7.0}};
  double Sb[4][3] = {
      {1., 0., -0.5}, {-0.7, 0.7, -0.5}, {-0.7, -0.7, -0.5}, {0., 0., 1.}};
  int TVb[4][3] = {{2, 1, 0}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}};
  Mesh M(Mesh::Mtype_sphere, 0, useVT, useTTi);
  int t, vi, v, i;
  double l;

  for (v = 0; v < n; v++) {
    l = std::sqrt(S[v][0] * S[v][0] + S[v][1] * S[v][1] + S[v][2] * S[v][2]);
    for (i = 0; i < 3; i++)
      S[v][i] = S[v][i] / l;
  }
  for (v = 0; v < 4; v++) {
    l = std::sqrt(Sb[v][0] * Sb[v][0] + Sb[v][1] * Sb[v][1] +
                  Sb[v][2] * Sb[v][2]);
    for (i = 0; i < 3; i++)
      Sb[v][i] = Sb[v][i] / l;
  }

  {
    std::ofstream F("ex.sphere1.io.s0", ios::out | ios::binary);
    Matrix3double FM(n, S);
    fmesh::IOHelperM<double>().cD(&FM).binary().OH(F).OD(F);
    F.close();
  }

  M.S_set(Matrix3double(n, S));

#ifdef FMESHER_WITH_X
  M.setX11VBigLimit(n);
  if (useX11)
    //   M.useX11(true,useX11text,500,500,-0.1,0.1,0.6,0.8);
    M.useX11(true, useX11text, 500, 500, -1.05, 1.05, -1.05, 1.05);
#endif

  /*
  M.S_append(Matrix3double(4,Sb));
  for (t=0;t<4;t++)
    for (vi=0;vi<3;vi++)
      TVb[t][vi] += n;
  M.TV_set(Matrix3int(4,TVb));
  */
  MeshC MC(&M);

  fmesh::vertexListT vertices;

  vertices.clear();
  for (v = 0; v < n; v++) {
    //    if ((v!=11) && (v!=12))
    vertices.push_back(v);
  }
  MC.DT(vertices);

  //  vertices.clear();
  // vertices.push_back(11);
  // vertices.push_back(12);
  // MC.DT(vertices);
  MC.CDTSegment(false, 7, 11);
  MC.CDTSegment(false, 10, 11);

  MC.PruneExterior();

  MC.RCDT(20, M_PI / 20.0);

  return 0;
}

int DTsphere_test2() {
  int n = 4;
  //  double S[4][3] = {{0.8,0.3,0.3},
  //		    {0.8,-0.4,0.2},
  //		    {0.8,0.4,-0.2},
  //		    {0.8,-0.4,-0.2}};
  // double S[4][3] = {{0.3,0.8,0.3},
  //		    {0.2,0.8,-0.4},
  //		    {-0.2,0.8,0.4},
  //		    {-0.2,0.8,-0.4}};
  double S[4][3] = {
      {0.4, 0.2, 0.8}, {0.4, -0.2, 0.8}, {-0.4, 0.2, 0.8}, {-0.4, -0.2, -0.3}};

  Mesh M(Mesh::Mtype_sphere, 0, useVT, useTTi);
  int t, vi, v, i;
  double l;

  for (v = 0; v < n; v++) {
    l = std::sqrt(S[v][0] * S[v][0] + S[v][1] * S[v][1] + S[v][2] * S[v][2]);
    for (i = 0; i < 3; i++)
      S[v][i] = S[v][i] / l;
  }

  {
    std::ofstream F("ex.sphere2.io.s0", ios::out | ios::binary);
    Matrix3double FM(n, S);
    fmesh::IOHelperM<double>().cD(&FM).binary().OH(F).OD(F);
    F.close();
  }

  M.S_set(Matrix3double(n, S));

#ifdef FMESHER_WITH_X
  M.setX11VBigLimit(n);
  if (useX11)
    M.useX11(true, useX11text, 500, 500, -1.05, 1.05, -1.05, 1.05);
#endif

  MeshC MC(&M);

  MC.CET(8, 0.2);

  fmesh::vertexListT vertices;
  /*
  vertices.push_back(4);
  vertices.push_back(8);
  MC.DT(vertices);
  MC.CDTSegment(false,4,8);
  */

  vertices.clear();
  for (v = 0; v < n; v++) {
    //    if ((v!=4) && (v!=8))
    vertices.push_back(v);
  }
  MC.DT(vertices);

  MC.PruneExterior();

  MC.RCDT(30, 100);
  MC.RCDT(25, M_PI / 20.0);

  cout << MC;

  return 0;
}

int koala_test() {
  int n = 15;
  PointRaw S[15] = {
      {585.6, 3553.4, 0.},  {537.6, 3497.6, 0.},  {1078.8, 3825.8, 0.},
      {2414.4, 3303.2, 0.}, {1604.4, 3957.2, 0.}, {581.4, 3515.6, 0.},
      {677.4, 3312.2, 0.},  {483.6, 3516.8, 0.},  {1531.2, 3923.6, 0.},
      {548.4, 3573.2, 0.},  {1647.6, 3907.4, 0.}, {1713.6, 3915.8, 0.},
      {2238, 1925, 0.},     {2085, 2151.2, 0.},   {1570.8, 2448.2, 0.}};

  int nb = 45;
  PointRaw Sb[45] = {
      {2125, 4011.1, 0.},   {4.5, 3975.8, 0.},
      {3192.8, 16.9, 0.},   {3248.2, 67.3, 0.},
      {3424.5, 11.9, 0.},   {3434.5, 16.9, 0.},
      {3565.5, 157.9, 0.},  {3948.3, 299, 0.},
      {4114.5, 293.9, 0.},  {4200.1, 263.7, 0.},
      {4336.1, 606.2, 0.},  {4356.3, 777.5, 0.},
      {4346.2, 842.9, 0.},  {4079.3, 913.5, 0.},
      {3998.7, 1276.1, 0.}, {3746.8, 1386.9, 0.},
      {3741.8, 1336.5, 0.}, {3429.5, 1215.7, 0.},
      {3142.4, 1240.8, 0.}, {3076.9, 1633.7, 0.},
      {3157.5, 1729.4, 0.}, {3535.3, 1920.8, 0.},
      {3585.7, 1920.8, 0.}, {3348.9, 2102.1, 0.},
      {3263.3, 2233.1, 0.}, {3182.7, 2253.2, 0.},
      {3001.4, 2238.1, 0.}, {2966.1, 2348.9, 0.},
      {3001.4, 2525.2, 0.}, {2905.7, 2631, 0.},
      {2905.7, 2751.9, 0.}, {2835.2, 2787.1, 0.},
      {2810, 2751.9, 0.},   {2815, 2686.4, 0.},
      {2719.3, 2666.3, 0.}, {2689.1, 2701.5, 0.},
      {2618.6, 3008.8, 0.}, {2658.9, 3134.7, 0.},
      {2699.2, 3149.8, 0.}, {2679, 3240.5, 0.},
      {2603.5, 3356.3, 0.}, {2517.8, 3391.6, 0.},
      {2507.8, 3436.9, 0.}, {2482.6, 3487.3, 0.},
      {2225.7, 3784.4, 0.}
      //,{2125,4011.1,0.}
  };
  int ne = 4;
  PointRaw Se[4] = {
      {0., 0., 0.}, {5000., 0., 0.}, {0., 5000., 0.}, {5000., 5000., 0.}};
  int TVe[2][3] = {{0, 1, 2}, {3, 2, 1}};
  Mesh M(Mesh::Mtype_plane, 0, useVT, useTTi);
  int t, vi, v;

#ifdef FMESHER_WITH_X
  M.setX11VBigLimit(n);
  if (useX11)
    M.useX11(true, useX11text, 500, 500, -440, 4400 + 440, -400 - 240,
             4000 + 400 + 240);
#endif

  M.S_set(Matrix3double(n, S));
  M.S_append(Matrix3double(nb, Sb));

  /*
  M.S_append(Matrix3double(ne,Se));
  for (t=0;t<2;t++)
    for (vi=0;vi<3;vi++)
      TVe[t][vi] += n+nb;
  M.TV_set(Matrix3int(2,TVe));
  MeshC MC(&M);
  */

  MeshC MC(&M);
  MC.CETplane(20, -0.05);

  fmesh::vertexListT vertices;
  for (v = 0; v < nb; v++)
    vertices.push_back(n + v);
  MC.DT(vertices);
  cout << MC;

  {
    Matrix1<int> bnd;
    for (v = 0; v < nb; v++)
      bnd(v) = n + v;
    bnd(nb) = n;

    {
      std::ofstream F("ex.koala.io.s0", ios::out | ios::binary);
      Matrix3double FM(n, S);
      FM.append(Matrix3double(nb, Sb));
      fmesh::IOHelperM<double>().cD(&FM).binary().OH(F).OD(F);
      F.close();
    }

    {
      std::ofstream F("ex.koala.io.bnd0", ios::out | ios::binary);
      fmesh::IOHelperM<int>().cD(&bnd).binary().OH(F).OD(F);
      F.close();
    }
  }

  fmesh::constrListT cinp;
  for (v = 0; v < nb - 1; v++)
    cinp.push_back(fmesh::constrT(n + v, n + v + 1));
  cinp.push_back(fmesh::constrT(n + nb - 1, n + 0));
  M.useVT(true);
  MC.CDTBoundary(cinp);
  M.useVT(useVT);
  cout << MC;

  vertices.clear();
  for (v = 0; v < n; v++)
    vertices.push_back(v);
  MC.DT(vertices);
  cout << MC;

  MC.PruneExterior();
  cout << MC;

  //  MC.setOptions(MC.getOptions()|MeshC::Option_offcenter_steiner);

  double *biglim = new double[n];
  for (v = 0; v < n; v++)
    biglim[v] = 100;

  MC.RCDT(25, 2000, biglim, n);
  cout << MC;

  delete[] biglim;

  return 0;
}

int iohelper_test() {
  std::string prefix("iotest.");
  std::ofstream OA;
  std::ofstream OB;

  {
    fmesh::Matrix3double M(4);
    M(0) = Point(0.0, 0.1, 0.2);
    M(1) = Point(1.0, 1.1, 1.2);
    M(2) = Point(2.0, 2.1, 2.2);
    M(3) = Point(3.0, 3.1, 3.2);

    IOHelperM<double> ioh;
    ioh.D(&M);

    OA.open((prefix + "DDGR.a").c_str(), ios::out | ios::binary);
    OB.open((prefix + "DDGR").c_str(), ios::out | ios::binary);
    ioh.rowmajor();
    ioh.binary().OH(OB).OD(OB).ascii().OH(OA).OD(OA);
    OA.close();
    OA.open((prefix + "DDGC.a").c_str(), ios::out | ios::binary);
    OB.close();
    OB.open((prefix + "DDGC").c_str(), ios::out | ios::binary);
    ioh.colmajor();
    ioh.binary().OH(OB).OD(OB).ascii().OH(OA).OD(OA);
    OA.close();
    OB.close();
  }

  {
    fmesh::Matrix3int M(4);
    M(0) = Int3(00, 01, 02);
    M(1) = Int3(10, 11, 12);
    M(2) = Int3(20, 21, 22);
    M(3) = Int3(30, 31, 32);

    IOHelperM<int> ioh;
    ioh.D(&M);

    OA.open((prefix + "DIGR.a").c_str(), ios::out | ios::binary);
    OB.open((prefix + "DIGR").c_str(), ios::out | ios::binary);
    ioh.rowmajor();
    ioh.binary().OH(OB).OD(OB).ascii().OH(OA).OD(OA);
    OA.close();
    OA.open((prefix + "DIGC.a").c_str(), ios::out | ios::binary);
    OB.close();
    OB.open((prefix + "DIGC").c_str(), ios::out | ios::binary);
    ioh.colmajor();
    ioh.binary().OH(OB).OD(OB).ascii().OH(OA).OD(OA);
    OA.close();
    OB.close();
  }

  {
    fmesh::SparseMatrix<double> M;
    M.cols(3).rows(4);
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 3; j++)
        M(i, j) = (double)i + (double)j / 10.0;

    IOHelperSM<double> ioh;
    ioh.D(&M);

    OA.open((prefix + "SDGR.a").c_str(), ios::out | ios::binary);
    OB.open((prefix + "SDGR").c_str(), ios::out | ios::binary);
    ioh.rowmajor();
    ioh.binary().OH(OB).OD(OB).ascii().OH(OA).OD(OA);
    OA.close();
    OA.open((prefix + "SDGC.a").c_str(), ios::out | ios::binary);
    OB.close();
    OB.open((prefix + "SDGC").c_str(), ios::out | ios::binary);
    ioh.colmajor();
    ioh.binary().OH(OB).OD(OB).ascii().OH(OA).OD(OA);
    OA.close();
    OB.close();

    OA.open((prefix + "SDSR.a").c_str(), ios::out | ios::binary);
    OB.open((prefix + "SDSR").c_str(), ios::out | ios::binary);
    ioh.symmetric();
    ioh.rowmajor();
    ioh.binary().OH(OB).OD(OB).ascii().OH(OA).OD(OA);
    OA.close();
    OA.open((prefix + "SDSC.a").c_str(), ios::out | ios::binary);
    OB.close();
    OB.open((prefix + "SDSC").c_str(), ios::out | ios::binary);
    ioh.colmajor();
    ioh.binary().OH(OB).OD(OB).ascii().OH(OA).OD(OA);
    OA.close();
    OB.close();

    OA.open((prefix + "SDDR.a").c_str(), ios::out | ios::binary);
    OB.open((prefix + "SDDR").c_str(), ios::out | ios::binary);
    ioh.diagonal();
    ioh.rowmajor();
    ioh.binary().OH(OB).OD(OB).ascii().OH(OA).OD(OA);
    OA.close();
    OA.open((prefix + "SDDC.a").c_str(), ios::out | ios::binary);
    OB.close();
    OB.open((prefix + "SDDC").c_str(), ios::out | ios::binary);
    ioh.colmajor();
    ioh.binary().OH(OB).OD(OB).ascii().OH(OA).OD(OA);
    OA.close();
    OB.close();
  }

  return 0;
}

void make_globe_test() {
  fmesh::Mesh M;
  M.type(fmesh::Mesh::Mtype_sphere);
#ifdef FMESHER_WITH_X
  if (useX11)
    M.useX11(true, useX11text, 500, 500, -1.05, 1.05, -1.05, 1.05);
#endif
  M.make_globe(30, 1.0);
#ifdef FMESHER_WITH_X
  if (useX11)
    M.useX11(true, useX11text, 500, 500, -1.05, 1.05, -1.05, 1.05);
#endif
  std::cout << M;
}

int main() {
  //  make_globe_test();
  //  return 0;
  //  iohelper_test();
  for (int i = 0; i < maxiter; i++) {
    CDT_test();
    DTsphere_test2();
    koala_test();
    DT2D_test();
    DT2D_test2();
    DT2D_test3();
    DTsphere_test();
  }

  return 0;
}
