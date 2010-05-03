#include <cstddef>
#include <vector>
#include <iostream>
#include <cmath>
#include <cstdlib>

#include "predicates.h"
#include "fmesher.h"

using std::vector;
using std::cout;
using std::endl;

using fmesh::Point;
using fmesh::Mesh;
using fmesh::Dart;
using fmesh::MeshC;


bool useX11 = true;
bool useX11text = false;
int maxiter = 1;

int predicates_test()
{
  vector< Point* > S1;
  double S2[5][3];
  Point* p;

  double a[3] = {1.0, 2.0, 3.0};
  double b[3] = {1.0, 2.0, 3.0};
  double c[3] = {1.0, 2.0, 3.0};
  p = &a;

  S1.push_back(&a);
  S1.push_back(&b);
  S1.push_back(&c);

  fmesh::predicates::orient2d(a,b,c);
  fmesh::predicates::orient2d(S2[0],S2[1],c);
  fmesh::predicates::orient2d(*S1[0],*S1[1],c);

  return 0;
}

int mesh_test()
{
  double S[4][3] = {{0.,0.,0.},
		    {1.,0.,0.},
		    {0.,1.,0.},
		    {1.,1.,0.}};
  int TV[2][3] = {{0,1,2},
		  {3,2,1}};
  // -std=c++0x
  //  S[0] = {1.,2.,3.};
  //  S[1] = {1.,2.,3.};
  //std::sqrt(2.0);

  Mesh M(Mesh::Mtype_plane,0,true,false);

  M.useTTi(false);
  M.S_set(S,4);
  M.TV_set(TV,2);

  Mesh M2 = M;
  M2.useTTi(false);

  Dart d = Dart(M,0,1,0);
  cout << d.onBoundary();
  cout << d.orbit2().onBoundary();
  cout << M;
  d = M.swapEdge(Dart(M,0,1,1));
  d = M.swapEdge(Dart(M,0,1,1));
  cout << M;
  Point s = {0.2,0.2,0.0};
  M.S_append(&s,1);
  cout << M;
  d = M.splitTriangle(Dart(M,1,1,0),M.nV()-1);
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


int LOP_test()
{
  int n = 0;
  double S[14][3] = {{0.5,0.5,0},
		     {0.5,0.6,0},
		     {0.3,0.2,0},
		     {0.3,0.6,0},
		     {0.5,0.3,0},
		     {0.6,0.7,0},
		     {0.7,0.3,0},
		     {0.2,0.8,0},
		     {0.9,0.5,0},
		     {0.1,0.1,0},
		     {0.05,0.05,0},
		     {0.95,0.05,0},
		     {0.95,0.95,0},
		     {0.05,0.95,0}};
  double Sb[6][3] = {{0.,0.,0.},
		     {1.,0.,0.},
		     {0.,1.,0.},
		     {1.0,1.,0.},
		     {0.6,0.7,0.},
		     {0.7,0.6,0.}};
  int TVb[6][3] = {{0,1,2},
		   {1,5,4},
		   {1,4,2},
		   {3,2,4},
		   {3,4,5},
		   {3,5,1}};
  Mesh M(Mesh::Mtype_plane,0,true,false);
  int t,vi,v;

  M.S_set(S,n);

  M.setX11VBigLimit(n);
  if (useX11)
    M.useX11(true,useX11text,500,500);

  M.S_append(Sb,6);
  for (t=0;t<6;t++)
    for (vi=0;vi<3;vi++)
      TVb[t][vi] += n; 
  M.TV_set(TVb,6);

  MeshC MC(&M,true);

  fmesh::triangleSetT triangles;
  for (t=0;t<(int)M.nT();t++)
    triangles.insert(t);
  MC.LOP(triangles);

  return 0;
}


int CDT_test()
{
  int n = 3;
  double S[3][3] = {{0.3,0.6,0},
		    {0.25,0.2,0},
		     {0.85,0.75,0}};
  double Sb[4][3] = {{0.,0.,0.},
		     {1.,0.,0.},
		     {0.,1.,0.},
		     {1.,1.,0.}};
  int TVb[2][3] = {{0,1,2},
		   {3,2,1}};
  Mesh M(Mesh::Mtype_plane,0,true,false);
  int t,vi,v;

  M.S_set(S,n);

  M.setX11VBigLimit(n);
  if (useX11)
    M.useX11(true,useX11text,500,500);

  M.S_append(Sb,4);
  for (t=0;t<2;t++)
    for (vi=0;vi<3;vi++)
      TVb[t][vi] += n; 
  M.TV_set(TVb,2);

  MeshC MC(&M,true);

  fmesh::vertexListT vertices;
  for (v=0;v<n;v++)
    vertices.push_back(v);
  MC.DT(vertices);

  cout << M;

  fmesh::constrListT cinp;
  cinp.push_back(fmesh::constrT(1,2));
  cinp.push_back(fmesh::constrT(3,6));
  MC.CDTInterior(cinp);

  MC.RCDT(1.415,100);
  MC.RCDT(1.415,0.05);

  return 0;
}


int DT2D_test()
{
  int n = 14;
  double S[14][3] = {{0.5,0.5,0},
		     {0.5,0.6,0},
		     {0.3,0.2,0},
		     {0.3,0.6,0},
		     {0.5,0.3,0},
		     {0.6,0.7,0},
		     {0.7,0.3,0},
		     {0.2,0.8,0},
		     {0.9,0.5,0},
		     {0.1,0.1,0},
		     {0.05,0.05,0},
		     {0.95,0.05,0},
		     {0.95,0.95,0},
		     {0.05,0.95,0}};
  double Sb[4][3] = {{0.,0.,0.},
		     {1.,0.,0.},
		     {0.,1.,0.},
		     {1.,1.,0.}};
  int TVb[2][3] = {{0,1,2},
		   {3,2,1}};
  Mesh M(Mesh::Mtype_plane,0,true,false);
  int t,vi,v;

  M.S_set(S,n);

  M.setX11VBigLimit(n);
  if (useX11)
    M.useX11(true,useX11text,500,500);

  M.S_append(Sb,4);
  for (t=0;t<2;t++)
    for (vi=0;vi<3;vi++)
      TVb[t][vi] += n; 
  M.TV_set(TVb,2);

  MeshC MC(&M,true);
  fmesh::vertexListT vertices;
    for (v=0;v<n;v++)
    vertices.push_back(v);

  MC.DT(vertices);

  cout << M;

  // fmesh::constrListT cinp;
  // cinp.push_back(fmesh::constrT(10,11));
  // cinp.push_back(fmesh::constrT(11,12));
  // cinp.push_back(fmesh::constrT(12,13));
  // cinp.push_back(fmesh::constrT(13,10));
  // MC.CDTBoundary(cinp);

  // cinp.clear();
  // cinp.push_back(fmesh::constrT(10,12));
  // MC.CDTInterior(cinp);

  MC.RCDT(1.415,100);
  MC.RCDT(1.415,0.05);

  return 0;
}


int DT2D_test2()
{
  int n = 25;
  double S[25][3] = {{0.1,0.1,0},
		     {0.3,0.1,0},
		     {0.7,0.1,0},
		     {0.9,0.1,0},
		     {0.1,0.3,0},
		     {0.3,0.3,0},
		     {0.7,0.3,0},
		     {0.9,0.3,0},
		     {0.1,0.5,0},
		     {0.3,0.5,0},
		     {0.7,0.5,0},
		     {0.9,0.5,0},
		     {0.1,0.7,0},
		     {0.3,0.7,0},
		     {0.7,0.7,0},
		     {0.9,0.7,0},
		     {0.1,0.9,0},
		     {0.3,0.9,0},
		     {0.7,0.9,0},
		     {0.9,0.9,0},
		     {0.5,0.1,0},
		     {0.5,0.3,0},
		     {0.5,0.5,0},
		     {0.5,0.7,0},
		     {0.5,0.9,0}};
  double Sb[4][3] = {{0.,0.,0.},
		     {1.,0.,0.},
		     {0.,1.,0.},
		     {1.,1.,0.}};
  int TVb[2][3] = {{0,1,2},
		   {3,2,1}};
  Mesh M(Mesh::Mtype_plane,0,true,false);
  int t,vi,v;

  M.S_set(S,n);

  M.setX11VBigLimit(n);
  if (useX11)
    M.useX11(true,useX11text,500,500);

  M.S_append(Sb,4);
  for (t=0;t<2;t++)
    for (vi=0;vi<3;vi++)
      TVb[t][vi] += n; 
  M.TV_set(TVb,2);

  MeshC MC(&M,true);
  fmesh::vertexListT vertices;
    for (v=0;v<n;v++)
    vertices.push_back(v);
  MC.DT(vertices);

  cout << M;

  MC.CDT(fmesh::constrListT());

  MC.RCDT(1.415,100);
  MC.RCDT(1.415,0.05);

  return 0;
}



int DT2D_test3() /* Random points */
{
  int n = 400;
  fmesh::Point S[400];
  double Sb[4][3] = {{0.,0.,0.},
		     {1.,0.,0.},
		     {0.,1.,0.},
		     {1.,1.,0.}};
  int TVb[2][3] = {{0,1,2},
		   {3,2,1}};
  Mesh M(Mesh::Mtype_plane,0,true,false);
  int t,vi,v;

  for (v=0;v<n;v++) {
    S[v][0] = double(std::rand())/RAND_MAX*0.9+0.05;
    S[v][1] = double(std::rand())/RAND_MAX*0.9+0.05;
    S[v][2] = 0.0;
  }

  M.S_set(S,n);

  M.setX11VBigLimit(n);
  if (useX11)
    M.useX11(true,useX11text,500,500);

  M.S_append(Sb,4);
  for (t=0;t<2;t++)
    for (vi=0;vi<3;vi++)
      TVb[t][vi] += n; 
  M.TV_set(TVb,2);

  MeshC MC(&M,true);
  fmesh::vertexListT vertices;
  for (v=0;v<n;v++)
    vertices.push_back(v);

  MC.DT(vertices);

  // fmesh::triangleSetT triangles;
  // for (t=0;t<(int)M.nT();t++)
  //   triangles.insert(t);
  // MC.LOP(triangles);

  MC.RCDT(1.415,100);
  MC.RCDT(1.415,0.05);

  return 0;
}




int DTsphere_test()
{
  int n = 10;
  double S[11][3] = {{0.2,0.2,0.8},
		     {-0.1,0.1,0.1},
		     {-0.2,0.5,-0.10},
		     {0.3,-0.2,0.2},
		     {0.4,0.6,0.2},
		     {-0.5,0.3,-0.3},
		     {-0.6,0.7,0.3},
		     {0.2,0.3,2.0},
		     {-0.8,-0.8,-0.4},
		     {-0.9,0.5,0.5},
		     {0.2,0.5,0.5}};
  double Sb[4][3] = {{1.,0.,-0.5},
		     {-0.7,0.7,-0.5},
		     {-0.7,-0.7,-0.5},
		     {0.,0.,1.}};
  int TVb[4][3] = {{2,1,0},
		   {0,1,3},
		   {1,2,3},
		   {2,0,3}};
  Mesh M(Mesh::Mtype_sphere,0,true,false);
  int t,vi,v,i;
  double l;

  for (v=0;v<n;v++) {
    l = std::sqrt(S[v][0]*S[v][0]+S[v][1]*S[v][1]+S[v][2]*S[v][2]);
    for (i=0;i<3;i++)
      S[v][i] = S[v][i]/l;
  }
  for (v=0;v<4;v++) {
    l = std::sqrt(Sb[v][0]*Sb[v][0]+Sb[v][1]*Sb[v][1]+Sb[v][2]*Sb[v][2]);
    for (i=0;i<3;i++)
      Sb[v][i] = Sb[v][i]/l;
  }

  M.S_set(S,n);

  M.setX11VBigLimit(n);
  if (useX11)
    M.useX11(true,useX11text,500,500,-1.05,1.05,-1.05,1.05);

  M.S_append(Sb,4);
  for (t=0;t<4;t++)
    for (vi=0;vi<3;vi++)
      TVb[t][vi] += n; 
  M.TV_set(TVb,4);

  cout << M;

  MeshC MC(&M,true);
  fmesh::vertexListT vertices;
  for (v=0;v<n;v++)
    vertices.push_back(v);
  MC.DT(vertices);

  cout << M;

  MC.RCDT(1.415,100);
  MC.RCDT(1.415,0.15);

  return 0;
}


int main()
{
  for (int i=0;i<maxiter;i++) {
    DT2D_test();
    DT2D_test2();
    DT2D_test3();
    CDT_test();
    DTsphere_test();
  }

  return 0;
}
