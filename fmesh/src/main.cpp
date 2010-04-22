#include <cstddef>
#include <vector>
#include <iostream>
#include <cmath>

#include "predicates.h"
#include "mesh.h"

using std::vector;
using std::cout;
using std::endl;

using fmesh::Point;
using fmesh::Mesh;
using fmesh::Dart;

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
  double S[4][3] = {0.,0.,0.,
		    1.,0.,0.,
		    0.,1.,0.,
		    1.,1.,0.};
  int TV[2][3] = {0,1,2,
		  3,2,1};
  // -std=c++0x
  //  S[0] = {1.,2.,3.};
  //  S[1] = {1.,2.,3.};
  //std::sqrt(2.0);

  Mesh M(0,false);

  M.useTTi(true);
  M.S_set(S,4);
  M.TV_set(TV,2);

  Mesh M2 = M;
  M2.useTTi(false);

  cout << M;
  cout << M2;

  /*
  Dart d(M);
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


int main()
{
  predicates_test();
  mesh_test();

  return 0;
}
