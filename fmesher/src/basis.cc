#include <cstddef>
#include <cstring>
#include <set>
#include <map>
#include <sstream>
#include <cmath>
#include <vector>
#ifndef NO_SPHERICAL_HARMONICS
#include "gsl/gsl_sf_legendre.h"
#endif

#include "vector.hh"
#include "ioutils.hh"

#define WHEREAMI __FILE__ << "(" << __LINE__ << ")\t"

#ifdef DEBUG
//#define _LOG(msg) std::cout << WHEREAMI << msg;
#define _LOG(msg)
#else
#define _LOG(msg)
#endif

#define M_2_SQRT_PI 3.5449077018110320546


using std::ios;
using std::cout;
using std::cin;
using std::endl;

namespace fmesh {

  int sph_basis_n(int kmax, bool rot_sym)
  {
    if (kmax >= 0) {
      if (rot_sym)
	return (kmax + 1);
      else
	return (kmax + 1) * (kmax + 1);
    } else {
      return 0;
    }
  }

  Matrix<double> spherical_harmonics(const Matrix3<double>& S,
				     size_t max_order,
				     bool rotationally_symmetric)
  {
    Matrix<double> sph(sph_basis_n(max_order,rotationally_symmetric));

#ifndef NO_SPHERICAL_HARMONICS
    size_t i, k, m;
    size_t GSL_res_n = gsl_sf_legendre_array_n(max_order);
    double* GSL_res_array = new double[GSL_res_n];

    if (rotationally_symmetric) {    
      for (i=0; i<S.rows(); i++) {
	gsl_sf_legendre_array(GSL_SF_LEGENDRE_SPHARM, max_order, S[i][2],
			      GSL_res_array);
	for (k = 0; k <= max_order; k++) {
	  sph(i,k) =
	    M_2_SQRT_PI * GSL_res_array[gsl_sf_legendre_array_index(k,0)];
	}
      }
    } else {
        double phi, scaling_sin, scaling_cos;

	std::vector<size_t> Idxs2(max_order + 1);
	// 0, 2, 6, 12, 20, 30, 42, 56, 72, 90, 110, 132, 156, 182, 210, 240
	for (k = 0; k <= max_order; k++) {
	  Idxs2[k] = k*(k+1);
	}

	for (i=0; i<S.rows(); i++) {
	  phi = atan2(S[i][1], S[i][0]);
	  
	  gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_SPHARM, max_order, S[i][2],
				  -1, GSL_res_array);
	  for (k = 0; k <= max_order; k++) {
	    sph(i,Idxs2[k]) =
	      M_2_SQRT_PI * GSL_res_array[gsl_sf_legendre_array_index(k,0)];
	  }
	  for (m = 1; m <= max_order; m++) {
	    scaling_sin = M_2_SQRT_PI * M_SQRT2 * sin(-(m * phi));
	    scaling_cos = M_2_SQRT_PI * M_SQRT2 * cos(m * phi);
	    for (k = m; k <= max_order; k++) {
	      sph(i,Idxs2[k] - m) =
		scaling_sin * GSL_res_array[gsl_sf_legendre_array_index(k,m)];
	      sph(i,Idxs2[k] + m) =
		scaling_cos * GSL_res_array[gsl_sf_legendre_array_index(k,m)];
	    }
	  }
	}

   }

    delete[] GSL_res_array;
#endif

    return sph;
  }


  Matrix<double> spherical_bsplines(const Matrix3<double>& S,
				    size_t n_basis,
				    size_t degree,
				    bool uniform_knot_angle_spacing)
  {
    Matrix<double> basis(n_basis);
    std::vector<double> knots(n_basis+degree+1);
    double s,s1,s2;
    std::vector< Matrix<double> > control(n_basis);
    std::vector< Matrix<double> > control_work(degree+1);
    size_t interval;

    for (size_t i=0; i<=degree; i++) {
      knots[i] = -1.0;
    }
    for (size_t i=degree+1; i<n_basis; i++) {
      knots[i] = (double(i-degree)/double(n_basis-degree))*2.0-1.0;
      if (uniform_knot_angle_spacing) {
	knots[i] = sin(knots[i]*M_PI/2.0);
      }
    }
    for (size_t i=n_basis; i<=n_basis+degree; i++) {
      knots[i] = 1.0;
    }

    for (size_t i=0; i<n_basis; i++) {
      control[i] = Matrix<double>(n_basis);
      control[i](0,i) = 1.0;
    }

    _LOG("degree\t" << degree << endl);
    _LOG("n_basis\t" << n_basis << endl);
    _LOG("n_basis+degree+1\t" << n_basis+degree+1 << endl);

    for (size_t coord_idx=0; coord_idx<S.rows(); coord_idx++) {
      s = S[coord_idx][2];

      _LOG("step 1, coord_idx\t" << coord_idx << endl);
      interval = degree;
      while ((interval+1<n_basis) & (s>=knots[interval+1]))
	interval++;
      
      _LOG("step 2" << endl);
      for (size_t i=0; i<=degree; i++)
	control_work[i] = control[i+interval-degree];
      
      _LOG("step 3" << endl);
      for (size_t k=1; k<=degree; k++)
	for (size_t i=degree; i>=k; i--) {
	  s1 = (knots[i+interval-k+1] - s)/(knots[i+interval-k+1]-knots[i+interval-degree]);
	  s2 = 1.0-s1;

	  for (size_t j=0; j<n_basis; j++)
	    control_work[i](0,j) = (s1 * control_work[i-1](0,j) +
				    s2 * control_work[i](0,j));
	}
      
      for (size_t j=0; j<n_basis; j++)
	basis(coord_idx,j) = control_work[degree](0,j);
    }

    return basis;
}



} /* namespace fmesh */
