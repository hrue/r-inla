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

bool useTV = true;
bool useTTi = true;
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

  Mesh M(Mesh::Mtype_plane,0,useTV,useTTi);

  M.S_set(S,4);
  M.TV_set(TV,2);

  Mesh M2 = M;

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


int CDT_test()
{
  int n = 7;
  Point S[7] = {{0.3,0.6,0},
		{0.35,0.1,0},
		{0.85,0.7,0},
		{0.7,0.1,0},
		{0.8,0.3,0},
		{0.9,0.3,0},
		{0.95,0.3,0}};
  Point Sb[4] = {{0.,0.,0.},
		 {1.,0.,0.},
		 {0.,1.,0.},
		 {0.99,1.,0.}};
  int TVb[2][3] = {{0,1,2},
		   {3,2,1}};
  Mesh M(Mesh::Mtype_plane,0,useTV,useTTi);
  int t,vi;

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
  vertices.push_back(0);
  vertices.push_back(1);
  vertices.push_back(2);
  vertices.push_back(3);
  vertices.push_back(4);
  vertices.push_back(5);
  vertices.push_back(6);
  MC.DT(vertices);

  MC.RCDT(1.415,100);

  fmesh::constrListT cinp;
  //  cinp.push_back(fmesh::constrT(1,2));
  cinp.push_back(fmesh::constrT(9,6));
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
  Mesh M(Mesh::Mtype_plane,0,useTV,useTTi);
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

  fmesh::constrListT cinp;
  cinp.push_back(fmesh::constrT(10,11));
  cinp.push_back(fmesh::constrT(11,12));
  cinp.push_back(fmesh::constrT(12,13));
  cinp.push_back(fmesh::constrT(13,10));
  MC.CDTBoundary(cinp);

  // cinp.clear();
  // cinp.push_back(fmesh::constrT(10,12));
  // MC.CDTInterior(cinp);

  MC.PruneExterior();

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
  Mesh M(Mesh::Mtype_plane,0,useTV,useTTi);
  int t,vi,v;

  /*
  for (int v=0;v<n;v++) {
    S[v][0] = S[v][0]+double(v)/1000.0;
    S[v][1] = S[v][1]+double(v)/1001.0;
  }
  */

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

  MC.CDT(fmesh::constrListT(),fmesh::constrListT());

  MC.PruneExterior();

  MC.RCDT(1.415,100);
  MC.RCDT(1.415,0.05);

  return 0;
}



int DT2D_test3() /* Random points */
{
  int n = 200;
  fmesh::Point S[200];
  double Sb[4][3] = {{0.,0.,0.},
		     {1.,0.,0.},
		     {0.,1.,0.},
		     {1.,1.,0.}};
  int TVb[2][3] = {{0,1,2},
		   {3,2,1}};
  Mesh M(Mesh::Mtype_plane,0,useTV,useTTi);
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
  //      M.useX11(true,useX11text,500,500,0.7,0.9,0.4,0.85);

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

  /*  
  fmesh::triangleSetT triangles;
  for (t=0;t<(int)M.nT();t++)
    triangles.insert(t);
  MC.LOP(triangles);
  */

  MC.CDTSegment(false,0,1);

  MC.PruneExterior();

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
  Mesh M(Mesh::Mtype_sphere,0,useTV,useTTi);
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

  MeshC MC(&M,true);

  fmesh::vertexListT vertices;
  vertices.push_back(4);
  vertices.push_back(8);
  MC.DT(vertices);
  MC.CDTSegment(false,4,8);

  vertices.clear();
  for (v=0;v<n;v++) {
    if ((v!=4) && (v!=8))
    vertices.push_back(v);
  }
  MC.DT(vertices);

  MC.PruneExterior();

  MC.RCDT(1.415,100);
  MC.RCDT(1.415,M_PI/20.0);

  return 0;
}


int koala_test()
{
  int n = 15;
  Point S[15] = {
    {585.6,3553.4,0.},
    {537.6,3497.6,0.},
    {1078.8,3825.8,0.},
    {2414.4,3303.2,0.},
    {1604.4,3957.2,0.},
    {581.4,3515.6,0.},
    {677.4,3312.2,0.},
    {483.6,3516.8,0.},
    {1531.2,3923.6,0.},
    {548.4,3573.2,0.},
    {1647.6,3907.4,0.},
    {1713.6,3915.8,0.},
    {2238,1925,0.},
    {2085,2151.2,0.},
    {1570.8,2448.2,0.}
  };

  int nb = 45;
  Point Sb[45] = {{
2125,4011.1,0.},{
4.5,3975.8,0.},{
3192.8,16.9,0.},{
3248.2,67.3,0.},{
3424.5,11.9,0.},{
3434.5,16.9,0.},{
3565.5,157.9,0.},{
3948.3,299,0.},{
4114.5,293.9,0.},{
4200.1,263.7,0.},{
4336.1,606.2,0.},{
4356.3,777.5,0.},{
4346.2,842.9,0.},{
4079.3,913.5,0.},{
3998.7,1276.1,0.},{
3746.8,1386.9,0.},{
3741.8,1336.5,0.},{
3429.5,1215.7,0.},{
3142.4,1240.8,0.},{
3076.9,1633.7,0.},{
3157.5,1729.4,0.},{
3535.3,1920.8,0.},{
3585.7,1920.8,0.},{
3348.9,2102.1,0.},{
3263.3,2233.1,0.},{
3182.7,2253.2,0.},{
3001.4,2238.1,0.},{
2966.1,2348.9,0.},{
3001.4,2525.2,0.},{
2905.7,2631,0.},{
2905.7,2751.9,0.},{
2835.2,2787.1,0.},{
2810,2751.9,0.},{
2815,2686.4,0.},{
2719.3,2666.3,0.},{
2689.1,2701.5,0.},{
2618.6,3008.8,0.},{
2658.9,3134.7,0.},{
2699.2,3149.8,0.},{
2679,3240.5,0.},{
2603.5,3356.3,0.},{
2517.8,3391.6,0.},{
2507.8,3436.9,0.},{
2482.6,3487.3,0.},{
2225.7,3784.4,0.}
		 //,{2125,4011.1,0.}
};
  int ne = 4;
  Point Se[4] = {{0.,0.,0.},
		 {5000.,0.,0.},
		 {0.,5000.,0.},
		 {5000.,5000.,0.}};
  int TVe[2][3] = {{0,1,2},
		   {3,2,1}};
  Mesh M(Mesh::Mtype_plane,0,useTV,useTTi);
  int t,vi,v;

  M.setX11VBigLimit(n);
  if (useX11)
    M.useX11(true,useX11text,500,500,-440,4400+440,-400-240,4000+400+240);

  M.S_set(S,n);
  M.S_append(Sb,nb);


  /*
  M.S_append(Se,ne);
  for (t=0;t<2;t++)
    for (vi=0;vi<3;vi++)
      TVe[t][vi] += n+nb;
  M.TV_set(TVe,2);
  MeshC MC(&M,true);
  */

  MeshC MC(&M,false);
  MC.CET(20,-0.05);

  fmesh::vertexListT vertices;
  for (v=0;v<nb;v++)
    vertices.push_back(n+v);
  MC.DT(vertices);

  fmesh::constrListT cinp;
  for (v=0;v<nb-1;v++)
    cinp.push_back(fmesh::constrT(n+v,n+v+1));
  cinp.push_back(fmesh::constrT(n+nb-1,n+0));
  MC.CDTBoundary(cinp);

  vertices.clear();
  for (v=0;v<n;v++)
    vertices.push_back(v);
  MC.DT(vertices);

  MC.PruneExterior();

  MC.RCDT(1.415,10000);

  MC.RCDT(1.415,100);

  return 0;
}





int main()
{
  for (int i=0;i<maxiter;i++) {
    koala_test();
    DT2D_test();
    DT2D_test2();
    DT2D_test3();
    CDT_test();
    DTsphere_test();
  }

  return 0;
}
